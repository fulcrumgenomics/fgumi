//! `SpillWrite` — final step of the block-parallel spill-write split
//! (`SpillGather` → `SpillBlockCompress` → `SpillWrite`).
//!
//! `SpillWrite` (`Serial + Affinity::Writer`, or `Detached` on the standalone
//! sort) receives the compressed [`SpillBlockEvent`]s in dense `ordinal` order
//! (the framework's `ByItemOrdinal` reorder feeds it like `WriteBgzfFile`),
//! demultiplexes `Block`s back to per-`file_id` spill files, and emits the
//! existing [`SortPhase1Event`] so `SortSpillDecompress` / `SortMerge` are
//! unchanged.
//!
//! Because `SortBuffer` (Serial) emits spill chunks one-at-a-time in `seq` order,
//! `SpillGather` (Serial) fans them in order, and the reorder preserves that
//! order, each file's blocks arrive **contiguously** — so `SpillWrite` only ever
//! holds **one** spill file open for writing at a time (`current`). It opens the
//! file on the first block (writing the codec magic), appends each compressed
//! block, and on `is_last_in_file` writes the codec trailer and closes it.
//!
//! # Runs stay closed until the merge, and are bounded by `--max-temp-files`
//!
//! A closed run is recorded in a `RunStack` as a path, **not** opened as a merge
//! slot. The merge cannot start before `AllAnnounced` anyway, and opening a slot
//! per run as it closed held one descriptor plus up to a slot's worth of
//! read-ahead per run for the rest of the sort — unbounded in the number of runs
//! and outside `--max-memory`. On `AllAnnounced` the surviving runs are opened and
//! announced (`SpillReady`, in `file_id` order) followed by `AllAnnounced` with
//! the surviving run count. The trade-off is that `SortSpillDecompress` no
//! longer pre-fills each run's slot queue (up to `PHASE2_DECOMP_CAP` blocks)
//! while spilling continues, so the merge starts on cold slots.
//!
//! When the stack's policy asks for it (see `RunStack`), a contiguous range of
//! runs is merged into one by a consolidation job: `SpillWrite` emits
//! `ConsolidationStart`, the read planner opens the inputs and plans their reads,
//! `ConsolidateRuns` merges them and `WriteConsolidatedRun` writes the merged
//! run, reporting completion through the [`ConsolidationRegistry`], which
//! unparks this step's thread. On completion the merged run replaces its inputs,
//! the inputs are unlinked and `ConsolidationRetire` follows. While a job is in
//! flight no input is taken, so upstream back-pressures while the live-run count
//! is being brought back under the limit.
//!
//! This step owns the `TmpDirAllocator` (Serial ⇒ the pick is uncontended) and
//! the RAII temp-dir handles, which live as long as the step so the spill files
//! outlive their reads.
//!
//! # Run indexes live in the run files until `AllAnnounced`
//!
//! For a sort that can merge by key range, each run's index is built as its
//! blocks are written, at the live stride (folded block by block, never
//! dense), then written as a footer of its own spill file when the run closes
//! and dropped ([`RunIndex::write_footer`]); a consolidated run's index is
//! built the same way by `WriteConsolidatedRun`, at its job's stride. When
//! the final run was spilled, `announce_runs` reads the footers back, after
//! the phase-1 arena is gone, and brings them to the final stride; with a
//! resident final run the merge is serial and the footers are written but
//! never read. Through phase 1 the only run indexes held are the open run's
//! and an in-flight consolidation's, each at its stride (an entry per
//! `stride` frames of that one run), outside `--max-memory`, which bounds the
//! records and refs exactly as in the serial sort.

use std::collections::VecDeque;
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use fgumi_bam_io::ProgressTracker;
use fgumi_sort::{RunIndex, SpillCodec, SpillKeyKind, TmpDirAllocator, spill_magic, spill_trailer};
use log::{info, warn};
use parking_lot::Mutex;
use tempfile::TempDir;

use crate::sort::consolidation::ConsolidationRegistry;
use crate::sort::protocol::{ConsolidationInput, SortPhase1Event, SpillBlockEvent};
use crate::sort::run_stack::{RunStack, SpillRun};
use crate::sort::spill_block_compress::SpillBufPool;
use fgumi_pipeline_core::{
    HeldRetry, Unpushed, external_wait,
    held::HeldSlot,
    outputs::Single,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{
        Affinity, CounterSpec, DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile,
    },
};

/// Counter slot index: bytes written to the spill file this call.
const SPILL_BYTES_WRITTEN: usize = 0;

/// Spill blocks written per `try_run`, whichever of this and
/// [`WRITE_BYTES_PER_CALL`] is reached first. One block per call made the
/// serial writer the phase-1 spill bound: ~12 µs of dispatch, reorder pop and
/// wake per ~9 KiB block capped it near 750 MB/s while ingest waited for the
/// arena. A batch stops at a file's last block, so run closing and
/// consolidation still happen between calls.
const WRITE_BLOCKS_PER_CALL: usize = 64;

/// Compressed bytes written per `try_run` (see [`WRITE_BLOCKS_PER_CALL`]):
/// ~10 ms of page-cache writes, short enough that the detached writer still
/// reports progress promptly.
const WRITE_BYTES_PER_CALL: u64 = 8 * 1024 * 1024;

/// Spill-run accounting shared with the end-of-run sort summary.
#[derive(Debug, Default)]
pub struct SpillRunStats {
    runs_written: AtomicU64,
    consolidations: AtomicU64,
    consolidation_nanos: AtomicU64,
    merge_sources: AtomicU64,
}

impl SpillRunStats {
    /// Spill runs written by the spill phase (before any consolidation).
    #[must_use]
    pub fn runs_written(&self) -> u64 {
        self.runs_written.load(Ordering::Relaxed)
    }

    /// Consolidations performed to honor `--max-temp-files`.
    #[must_use]
    pub fn consolidations(&self) -> u64 {
        self.consolidations.load(Ordering::Relaxed)
    }

    /// Wall time of the consolidation jobs (start to completion), in seconds.
    #[must_use]
    #[allow(clippy::cast_precision_loss)] // display value; > 2^52 ns is ~52 days
    pub fn consolidation_secs(&self) -> f64 {
        self.consolidation_nanos.load(Ordering::Relaxed) as f64 / 1e9
    }

    /// Wall time of the consolidation jobs, in nanoseconds.
    #[must_use]
    pub fn consolidation_nanos(&self) -> u64 {
        self.consolidation_nanos.load(Ordering::Relaxed)
    }

    /// Spill runs handed to the final merge (after consolidation).
    #[must_use]
    pub fn merge_sources(&self) -> u64 {
        self.merge_sources.load(Ordering::Relaxed)
    }

    /// Count one spill run written.
    pub fn record_run_written(&self) {
        self.runs_written.fetch_add(1, Ordering::Relaxed);
    }

    /// Count one consolidation job that took `wall` from start to completion.
    pub fn record_consolidation(&self, wall: Duration) {
        self.consolidations.fetch_add(1, Ordering::Relaxed);
        self.consolidation_nanos
            .fetch_add(u64::try_from(wall.as_nanos()).unwrap_or(u64::MAX), Ordering::Relaxed);
    }

    /// Record how many spill runs the final merge reads.
    pub fn set_merge_sources(&self, sources: u64) {
        self.merge_sources.store(sources, Ordering::Relaxed);
    }
}

/// The one spill file currently being written (open from its first block until
/// its `is_last_in_file` block).
struct OpenSpill {
    file_id: u32,
    key_kind: SpillKeyKind,
    path: PathBuf,
    writer: BufWriter<File>,
    /// Bytes handed to `writer` so far, the codec magic included: the file
    /// offset the next block lands at.
    bytes_written: u64,
    /// The file's frames so far, while every block written carried its
    /// frames; `None` for good once one did not (an index covers every frame
    /// of its file or none).
    index: Option<RunIndex>,
}

/// A consolidation job in flight.
pub(crate) struct ActiveJob {
    job_id: u32,
    /// Range of `runs` being merged. No run is pushed while a job is in flight,
    /// so the range stays valid until the job completes.
    range: std::ops::Range<usize>,
    output: PathBuf,
    started: Instant,
}

/// Fold one written block's frames into its run's `index`, the rule both run
/// writers (`SpillWrite` and `WriteConsolidatedRun`) follow: the first
/// block's frames start the index, at `stride` frames per entry; a later
/// block's extend it at the index's own stride ([`RunIndex::append_folded`],
/// so the run's index is never held denser than that); a block without
/// frames ends it (`None` for good: an index covers every frame or none), and
/// so does an entry whose sums overflow (logged; the merge is then serial). A
/// block's frames are indexed from its own first byte, so they are rebased by
/// `base`, the file position the block lands at.
pub(crate) fn fold_block_frames(
    index: &mut Option<RunIndex>,
    first_block: bool,
    frames: Option<RunIndex>,
    base: u64,
    stride: u32,
) {
    let folded = match (first_block, index.take(), frames) {
        (true, _, Some(frames)) => {
            let mut open = RunIndex::new_strided(frames.key_width(), stride.max(1));
            open.rebase(base);
            open.append_folded(&frames, base).map(|()| open)
        }
        (false, Some(mut open), Some(frames)) => open.append_folded(&frames, base).map(|()| open),
        (_, _, None) | (false, None, Some(_)) => return,
    };
    *index =
        folded.inspect_err(|e| log::warn!("Run index dropped: {e}; the merge will be serial")).ok();
}

/// Sparsify `index` to `stride` (a power of two at or above its own);
/// `None` (logged) when its groups overflow ([`RunIndex::sparsify`]), and the
/// run goes unindexed, so the merge falls back to serial.
pub(crate) fn sparsify_to(mut index: RunIndex, stride: u32) -> Option<RunIndex> {
    match index.sparsify(stride / index.stride().max(1)) {
        Ok(()) => Some(index),
        Err(e) => {
            log::warn!("Run index dropped: {e}; the merge will be serial");
            None
        }
    }
}

/// Write `index` as the footer of a run file whose frames `writer` has just
/// written (see [`RunIndex::write_footer`]). Returns whether a footer was
/// written: `false` for no index.
///
/// # Errors
///
/// Any write error: a run whose footer failed is not a run the merge can
/// read, so the sort fails rather than silently dropping the index.
pub(crate) fn write_index_footer<W: Write>(
    writer: &mut W,
    index: Option<&RunIndex>,
) -> io::Result<bool> {
    let Some(index) = index else { return Ok(false) };
    index.write_footer(writer)?;
    Ok(true)
}

/// `Serial + Affinity::Writer` step that writes per-`file_id` spill files from
/// the compressed block stream, bounds the live runs by consolidating them, and
/// emits `SortPhase1Event`s.
pub struct SpillWrite {
    /// Shared temp-directory allocator (free-space-aware round-robin). `Serial`,
    /// so the lock is effectively uncontended (one writer worker).
    alloc: Arc<Mutex<TmpDirAllocator>>,
    /// Spill codec for chunk files (bgzf or zstd).
    codec: SpillCodec,
    /// The currently-open spill file, if any.
    current: Option<OpenSpill>,
    /// Closed runs awaiting the merge, and the consolidation policy over them.
    runs: RunStack,
    /// Key type of every run so far; all runs of one sort share it.
    key_kind: Option<SpillKeyKind>,
    /// The consolidation job in flight, if any.
    pub(crate) active: Option<ActiveJob>,
    /// Consolidation jobs started: the next job id, and the merged file's name.
    jobs_started: u32,
    /// Where consolidation jobs report completion; `None` leaves consolidation
    /// off (no tail downstream runs a job).
    registry: Option<Arc<ConsolidationRegistry>>,
    /// Spill runs this step has closed; checked against upstream's announced
    /// count. Kept apart from `stats`, which a caller may share.
    runs_closed: u64,
    /// `true` once `AllAnnounced` has been processed; any later event is a
    /// protocol violation, and a drained input without it strands the runs.
    announced: bool,
    /// Events produced but not yet pushed downstream.
    outbox: VecDeque<SortPhase1Event>,
    stats: Arc<SpillRunStats>,
    held: HeldSlot<Unpushed<SortPhase1Event>>,
    output_byte_limit: u64,
    /// Compressed spill bytes written, logged every 256 MiB under `RUST_LOG=info`
    /// so the spill-write rate over wall time is visible alongside ingest.
    spill_progress: ProgressTracker,
    /// RAII temp-dir handles, held for the step's lifetime so spill files survive
    /// while being read by `SortMerge`.
    #[allow(dead_code)]
    temp_dirs: Arc<Vec<TempDir>>,
    /// When `true`, advertise `StepKind::Detached` so the framework drives this
    /// writer on its own dedicated thread (off the pool) instead of as a
    /// pool-scheduled `Serial + Affinity::Writer` step. Set only on the
    /// standalone-sort spill path via [`Self::with_detached`] — the exact
    /// Phase-1 analogue of Lever 2's detached terminal writer — so the single
    /// serial write stream stops consuming a compute worker that could be
    /// compressing. Every other chain leaves it `false`.
    detached: bool,
    /// Where each written block's buffer is returned for reuse by
    /// `SpillBlockCompress` (`None`: dropped).
    buf_pool: Option<Arc<SpillBufPool>>,
    /// `true` once the phase-1 front stopped indexing (a block stamped
    /// stride 0): every run's index was dropped, later blocks' frames are
    /// discarded and later runs are not indexed.
    indexing_off: bool,
    /// Frames per index entry of every live run's index (a power of two):
    /// the largest stride the front has stamped so far.
    index_stride: u32,
}

impl SpillWrite {
    /// Return each written block's buffer to `pool` (shared with
    /// `SpillBlockCompress`).
    #[must_use]
    pub fn with_buf_pool(mut self, pool: Arc<SpillBufPool>) -> Self {
        self.buf_pool = Some(pool);
        self
    }

    /// Build a `SpillWrite`. `alloc` names spill files across the configured temp
    /// dirs; `codec` selects the on-disk format; `temp_dirs` holds the RAII
    /// handles alive for the step's lifetime. `output_byte_limit` byte-bounds the
    /// forwarded-event output queue.
    ///
    /// Consolidation is off until [`Self::with_consolidation`] is applied.
    #[must_use]
    pub fn new(
        alloc: Arc<Mutex<TmpDirAllocator>>,
        codec: SpillCodec,
        output_byte_limit: u64,
        temp_dirs: Arc<Vec<TempDir>>,
    ) -> Self {
        Self {
            alloc,
            codec,
            current: None,
            runs: RunStack::new(0),
            key_kind: None,
            active: None,
            jobs_started: 0,
            registry: None,
            runs_closed: 0,
            announced: false,
            outbox: VecDeque::new(),
            stats: Arc::new(SpillRunStats::default()),
            held: HeldSlot::new(),
            output_byte_limit,
            spill_progress: ProgressTracker::new("Spill bytes written")
                .with_interval(256 * 1024 * 1024),
            temp_dirs,
            detached: false,
            buf_pool: None,
            indexing_off: false,
            index_stride: 1,
        }
    }

    /// `(indexing on, stride)` of the phase-1 index policy (tests).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn index_policy_for_test(&self) -> (bool, u32) {
        (!self.indexing_off, self.index_stride)
    }

    /// `(stride, entries)` of the open run's index, the only run index this
    /// step holds before `AllAnnounced` (tests): a closed run's lives in its
    /// file's footer.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn open_index_shape_for_test(&self) -> Option<(u32, usize)> {
        let index = self.current.as_ref()?.index.as_ref()?;
        Some((index.stride(), index.len()))
    }

    /// Follow the index stride the phase-1 front stamped on a block or on
    /// `AllAnnounced` (which carries the decision made for a resident final
    /// run, which no block carries). The front owns the index policy; this
    /// step keeps none of its own. `0` drops the open run's index, leaves the
    /// closed runs' footers unread and stops indexing for good; a stride
    /// above the live one becomes the stride closing runs are written at and
    /// every footer is brought to at `AllAnnounced`; a lower one changes
    /// nothing, as an index cannot be made denser again.
    fn follow_run_stride(&mut self, index_stride: u32) {
        if self.indexing_off {
            return;
        }
        if index_stride == 0 {
            self.indexing_off = true;
            let mut had_index = self.runs.runs().iter().any(|run| run.index_footer);
            if let Some(open) = self.current.as_mut() {
                had_index |= open.index.take().is_some();
            }
            let msg = format!(
                "Run indexes dropped after {} spill runs: the front found no range fit",
                self.runs_closed
            );
            if had_index {
                info!("{msg}");
            } else {
                log::debug!("{msg}");
            }
            return;
        }
        self.index_stride = self.index_stride.max(index_stride);
    }

    /// Keep at most `max_temp_files` spill runs live, consolidating when the
    /// limit is reached through the chain's consolidation tail, which reports
    /// each job's completion to `registry`. A limit below 2 leaves consolidation
    /// off ([`consolidation_can_trigger`](super::consolidation_can_trigger)).
    #[must_use]
    pub fn with_consolidation(
        mut self,
        max_temp_files: usize,
        registry: Arc<ConsolidationRegistry>,
    ) -> Self {
        self.runs = RunStack::new(max_temp_files);
        self.registry = Some(registry);
        self
    }

    /// Live runs (tests).
    #[cfg(test)]
    pub(crate) fn live_runs_for_test(&self) -> usize {
        self.runs.runs().len()
    }

    /// Share spill-run accounting with the caller (the sort summary).
    #[must_use]
    pub fn with_stats(mut self, stats: Arc<SpillRunStats>) -> Self {
        self.stats = stats;
        self
    }

    /// Run this spill writer on its own dedicated `StepKind::Detached` thread
    /// instead of as a pool-scheduled `Serial + Affinity::Writer` step. Used
    /// ONLY on the standalone-sort spill path (the Phase-1 analogue of Lever 2's
    /// detached terminal writer): it frees a pool worker for the
    /// compression-bound `SpillBlockCompress` work, matching feat-runall's dedicated
    /// spill-I/O thread — but as a single persistent thread for the whole run,
    /// not one per spill chunk.
    ///
    /// The `try_run` body and the bytes it writes are unchanged: the
    /// dedicated-thread driver pops blocks in the same `ByItemOrdinal`
    /// reorder-stage-ordered sequence, so each spill file's blocks still arrive
    /// contiguously (the one-open-file-at-a-time invariant holds) and every
    /// spill file is byte-identical to the pool-scheduled writer's output.
    /// Affinity is ignored for `Detached`.
    #[must_use]
    pub fn with_detached(mut self) -> Self {
        self.detached = true;
        self
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        !matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld)
    }

    /// Allocate a spill path for `file_id` (named by the logical spill index so
    /// the merge tie-break is independent of write order) and create the file,
    /// writing the codec magic prologue.
    fn open_file(&self, file_id: u32, key_kind: SpillKeyKind) -> io::Result<OpenSpill> {
        let path = self.next_path(&format!("chunk_{file_id:04}.keyed"))?;
        // `create_new` fails closed on a duplicate/stale path: a reused `file_id`
        // (or a leftover file) must surface as an error rather than truncate an
        // existing spill and silently corrupt merge input.
        let writer = external_wait(|| {
            let file = OpenOptions::new().write(true).create_new(true).open(&path)?;
            let mut writer = BufWriter::with_capacity(256 * 1024, file);
            writer.write_all(spill_magic(self.codec))?;
            Ok::<_, io::Error>(writer)
        })?;
        let bytes_written = spill_magic(self.codec).len() as u64;
        Ok(OpenSpill { file_id, key_kind, path, writer, bytes_written, index: None })
    }

    /// A path named `name` in the next temp dir the allocator picks.
    fn next_path(&self, name: &str) -> io::Result<PathBuf> {
        let base = self.alloc.lock().next().map_err(|e| {
            io::Error::other(format!("SpillWrite: temp-dir allocation failed: {e:#}"))
        })?;
        Ok(base.join(name))
    }

    /// Process one input event, performing any disk writes and queueing the
    /// `SortPhase1Event`s it produces in `outbox`. A closed run may start a
    /// consolidation job (see [`Self::poll_registry`]). `StepCtx`-free for unit
    /// testing.
    ///
    /// # Errors
    ///
    /// Propagates file-create / write / slot-open errors. Also
    /// errors if a block arrives for a different `file_id` than the open file
    /// while one is open without an intervening `is_last_in_file` — a
    /// framework-ordering invariant violation that must fail loud rather than
    /// corrupt a spill — or if runs of one sort disagree on their key type.
    pub(crate) fn process_event(&mut self, event: SpillBlockEvent) -> io::Result<()> {
        if self.announced {
            // The runs were handed downstream at `AllAnnounced`; anything after
            // it would be written to a stack nobody reads.
            return Err(io::Error::other(
                "SpillWrite: input arrived after AllAnnounced (the run set was already announced)",
            ));
        }
        match event {
            SpillBlockEvent::Range { .. } => Err(io::Error::other(
                "SpillWrite received an unframed Range; SpillBlockCompress frames every range",
            )),
            SpillBlockEvent::Block {
                file_id,
                key_kind,
                is_last_in_file,
                records_ingested_so_far,
                bytes,
                frames,
                index_stride,
                ..
            } => {
                self.follow_run_stride(index_stride);
                // Open the file on its first block; otherwise the open file must
                // match (blocks for one file are contiguous in the ordinal stream).
                // Spill file I/O blocks on the filesystem — minutes on a busy
                // parallel filesystem — which is an external wait, not a
                // wedge (#1068).
                let first_block = self.current.is_none();
                if first_block {
                    self.current = Some(self.open_file(file_id, key_kind)?);
                }
                let open = self.current.as_mut().expect("open file set above");
                if open.file_id != file_id {
                    return Err(io::Error::other(format!(
                        "SpillWrite: block for file_id {file_id} arrived while file_id {} \
                         was still open (blocks must be contiguous per file)",
                        open.file_id
                    )));
                }
                if open.key_kind != key_kind {
                    return Err(io::Error::other(format!(
                        "SpillWrite: spill file_id {file_id} mixes key types {:?} and \
                         {key_kind:?}",
                        open.key_kind
                    )));
                }
                // The first block's frames start the file's index; a later
                // block's extend it, and a block without frames ends it. A
                // block's frames are indexed from its own first byte, so they
                // are rebased by the file position the block lands at.
                let base = open.bytes_written;
                // Once the front stopped indexing, frames are dropped as they
                // arrive.
                let frames = frames.filter(|_| !self.indexing_off);
                fold_block_frames(&mut open.index, first_block, frames, base, self.index_stride);
                external_wait(|| open.writer.write_all(&bytes))?;
                open.bytes_written += bytes.len() as u64;
                self.spill_progress.log_if_needed(bytes.len() as u64);
                if let Some(pool) = &self.buf_pool {
                    pool.give(bytes);
                }

                if is_last_in_file {
                    self.close_run(records_ingested_so_far)?;
                }
                Ok(())
            }
            SpillBlockEvent::Residual { chunk, records_ingested_so_far, .. } => {
                self.ensure_no_open_file("residual")?;
                // Wrap in a fresh, uniquely-owned `Arc`: the chunk is only ever
                // moved (never cloned) onward, so `SortMerge`'s `Arc::try_unwrap`
                // invariant holds.
                self.outbox.push_back(SortPhase1Event::MemoryChunk {
                    chunk: Arc::new(chunk),
                    records_ingested_so_far,
                });
                Ok(())
            }
            SpillBlockEvent::AllAnnounced {
                slot_count,
                memory_chunk_count,
                total_records,
                index_stride,
                ..
            } => {
                self.ensure_no_open_file("AllAnnounced")?;
                // The front's last decision: a stride made for a resident
                // final run reaches this step only here.
                self.follow_run_stride(index_stride);
                self.announce_runs(slot_count, memory_chunk_count, total_records)
            }
        }
    }

    /// Finish the open spill file, record it as a live run, and start a
    /// consolidation if the policy now calls for one.
    /// Process `first` and then events from `next` as one batch, returning the
    /// spill bytes written (only `Block` carries any). A batch holds at most
    /// [`WRITE_BLOCKS_PER_CALL`] events and stops once it has written
    /// [`WRITE_BYTES_PER_CALL`] bytes, and it ends after any event that is not a
    /// mid-file `Block`: a file's last block (which closes the run and may start
    /// a consolidation), a residual or `AllAnnounced` (which queue output). It
    /// also ends as soon as output is queued or a consolidation is in flight, so
    /// `try_run` handles those before writing more.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::process_event`] errors.
    fn write_batch(
        &mut self,
        first: SpillBlockEvent,
        mut next: impl FnMut() -> Option<SpillBlockEvent>,
    ) -> io::Result<u64> {
        let mut spill_bytes: u64 = 0;
        let mut events = 0usize;
        let mut event = first;
        loop {
            let (len, continues) = match &event {
                SpillBlockEvent::Block { bytes, is_last_in_file, .. } => {
                    (bytes.len() as u64, !*is_last_in_file)
                }
                SpillBlockEvent::Range { .. }
                | SpillBlockEvent::Residual { .. }
                | SpillBlockEvent::AllAnnounced { .. } => (0, false),
            };
            self.process_event(event)?;
            spill_bytes += len;
            events += 1;
            if !continues
                || events >= WRITE_BLOCKS_PER_CALL
                || spill_bytes >= WRITE_BYTES_PER_CALL
                || !self.outbox.is_empty()
                || self.active.is_some()
            {
                return Ok(spill_bytes);
            }
            match next() {
                Some(following) => event = following,
                None => return Ok(spill_bytes),
            }
        }
    }

    fn close_run(&mut self, records_ingested_so_far: u64) -> io::Result<()> {
        let OpenSpill { file_id, key_kind, path, mut writer, bytes_written, index } =
            self.current.take().expect("open file present");
        // The run's index, at the live stride (smaller footer, less to read
        // back), goes into the file after the frames and is dropped here: no
        // closed run's index is resident through phase 1. `bytes` is where the
        // frames end, what the merge and a consolidation read.
        let index = index.filter(|_| !self.indexing_off).and_then(|ix| {
            debug_assert!(spill_trailer(self.codec).is_empty(), "only zstd runs are indexed");
            sparsify_to(ix, self.index_stride)
        });
        let bytes = bytes_written + spill_trailer(self.codec).len() as u64;
        // Dropping the writer closes the file, which can block on a parallel
        // filesystem (NFS, Lustre) as long as the writes, so it is bracketed too.
        let index_footer = external_wait(|| {
            writer.write_all(spill_trailer(self.codec))?;
            let index_footer = write_index_footer(&mut writer, index.as_ref())?;
            writer.flush()?;
            drop(writer);
            Ok::<_, io::Error>(index_footer)
        })
        .map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("SpillWrite: writing spill run {}: {e}", path.display()),
            )
        })?;
        drop(index);
        match self.key_kind {
            None => self.key_kind = Some(key_kind),
            Some(kind) if kind == key_kind => {}
            Some(kind) => {
                return Err(io::Error::other(format!(
                    "SpillWrite: spill file_id {file_id} is keyed {key_kind:?} but earlier runs \
                     are keyed {kind:?}"
                )));
            }
        }
        self.runs_closed += 1;
        self.stats.record_run_written();
        self.runs.push(SpillRun { file_id, path, bytes, records_ingested_so_far, index_footer });
        self.start_consolidation_if_needed()
    }

    /// Start the consolidation job the run stack asks for, if any, none is in
    /// flight and a consolidation tail is wired: record its output with the
    /// registry and queue its `ConsolidationStart`.
    fn start_consolidation_if_needed(&mut self) -> io::Result<()> {
        if self.active.is_some() {
            return Ok(());
        }
        let Some(registry) = self.registry.clone() else { return Ok(()) };
        let Some(range) = self.runs.next_merge() else { return Ok(()) };
        let key_kind = self.key_kind.expect("a live run implies a known key kind");
        let inputs: Vec<ConsolidationInput> = self.runs.runs()[range.clone()]
            .iter()
            .map(|r| ConsolidationInput { file_id: r.file_id, path: r.path.clone(), len: r.bytes })
            .collect();
        let input_bytes: u64 = inputs.iter().map(|i| i.len).sum();
        let job_id = self.jobs_started;
        let output = self.next_path(&format!("merged_{job_id:04}.keyed"))?;
        self.jobs_started += 1;
        info!(
            "Consolidating {} spill runs ({input_bytes} bytes) into 1; {} live runs",
            inputs.len(),
            self.runs.runs().len()
        );
        registry.begin(job_id, output.clone());
        let index_stride = if self.indexing_off { 0 } else { self.index_stride };
        self.outbox.push_back(SortPhase1Event::ConsolidationStart {
            job_id,
            inputs,
            key_kind,
            index_stride,
        });
        self.active = Some(ActiveJob { job_id, range, output, started: Instant::now() });
        Ok(())
    }

    /// Consume the in-flight job's completion, if it has arrived: replace its
    /// inputs with the merged run, unlink them, queue `ConsolidationRetire` and
    /// start the next job the policy asks for. Returns whether a completion was
    /// consumed. With none, registers this thread with the registry and checks
    /// once more (a completion between the first check and the registration
    /// would otherwise wait on this thread's timer), so a `false` lets the
    /// caller report `NoProgress` and park: the completion unparks it.
    ///
    /// # Errors
    ///
    /// A completion for a job other than the one in flight (a protocol
    /// violation), or a failure to start the next job.
    pub(crate) fn poll_registry(&mut self) -> io::Result<bool> {
        let Some(registry) = self.registry.clone() else { return Ok(false) };
        if self.active.is_none() {
            return Ok(false);
        }
        let completed = registry.take_completed().or_else(|| {
            #[cfg(test)]
            registry.run_hook_before_register_for_test();
            registry.register_waiter();
            registry.take_completed()
        });
        let Some((job_id, merged_bytes, index_footer)) = completed else { return Ok(false) };
        let ActiveJob { job_id: active_id, range, output, started } =
            self.active.take().expect("active job checked above");
        if job_id != active_id {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "SpillWrite: consolidation job {job_id} completed while job {active_id} was \
                     in flight"
                ),
            ));
        }
        // The merged run's index is in its footer, read back with the other
        // runs' at `AllAnnounced`.
        let absorbed = self.runs.complete_merge(range, output, merged_bytes, index_footer);
        for run in &absorbed {
            if let Err(e) = std::fs::remove_file(&run.path) {
                // The temp dir removes it at exit regardless; only disk headroom
                // is at stake until then.
                warn!("SpillWrite: could not remove consolidated run {}: {e}", run.path.display());
            }
        }
        self.outbox.push_back(SortPhase1Event::ConsolidationRetire { job_id });
        self.stats.record_consolidation(started.elapsed());
        self.start_consolidation_if_needed()?;
        Ok(true)
    }

    /// On `AllAnnounced`: open a merge slot for every surviving run (its
    /// frames ending where they were written), attach the index read back
    /// from its footer (sparsified to the final stride), and queue its
    /// `SpillReady` (in `file_id` order), then `AllAnnounced` with the
    /// surviving run count.
    ///
    /// The footers are read back only when the final run was spilled
    /// (`memory_chunk_count == 0`): the phase-1 arena is retired then, so the
    /// indexes are resident only from here on, in the merge's budget. With a
    /// resident final run the merge is serial whatever the indexes say
    /// ([`SerialReason::ResidualInMemory`](crate::sort::SerialReason::ResidualInMemory)),
    /// so the footers stay unread rather than sit uncharged beside the
    /// resident arena.
    ///
    /// # Errors
    ///
    /// Errors if a slot cannot be opened, if the upstream run count disagrees
    /// with the runs this step closed (a lost or duplicated run), if a run's
    /// frames do not end where they were written or its footer is missing,
    /// unexpected or unreadable (a damaged run file, which must not be read as
    /// records), or if a footer's key width is not the runs' key kind's.
    fn announce_runs(
        &mut self,
        runs_announced: u32,
        memory_chunk_count: u32,
        total_records: u64,
    ) -> io::Result<()> {
        debug_assert!(self.active.is_none(), "no input is taken while a job is in flight");
        let runs_written = self.runs_closed;
        if u64::from(runs_announced) != runs_written {
            return Err(io::Error::other(format!(
                "SpillWrite: upstream announced {runs_announced} spill runs but {runs_written} \
                 were written"
            )));
        }
        self.announced = true;
        let runs = std::mem::replace(&mut self.runs, RunStack::new(0)).into_runs();
        let slot_count = u32::try_from(runs.len()).expect("live runs bounded by the run count");
        self.stats.set_merge_sources(runs.len() as u64);
        let read_back = memory_chunk_count == 0;
        for SpillRun { file_id, path, bytes, records_ingested_so_far, index_footer } in runs {
            let frames_end = fgumi_sort::FramesEnd::At(bytes);
            let slot = external_wait(|| fgumi_sort::open_spill_slot(&path, file_id, frames_end))
                .map_err(|e| {
                    let kind =
                        e.downcast_ref::<io::Error>().map_or(io::ErrorKind::Other, io::Error::kind);
                    io::Error::new(
                        kind,
                        format!("SpillWrite: failed to open spill slot {}: {e:#}", path.display()),
                    )
                })?;
            if slot.has_index_footer() != index_footer {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "SpillWrite: spill run {} was written with {} index footer, but reads \
                         with {} footer",
                        path.display(),
                        if index_footer { "an" } else { "no" },
                        if slot.has_index_footer() { "an" } else { "no" },
                    ),
                ));
            }
            if read_back && let Some(index) = self.read_back_index(&slot, &path)? {
                slot.attach_run_index(Arc::new(index));
            }
            // Every run closed with its blocks' key kind, which `close_run`
            // recorded.
            slot.attach_key_kind(self.key_kind.expect("a closed spill run set the key kind"));
            self.outbox.push_back(SortPhase1Event::SpillReady {
                slot,
                path,
                records_ingested_so_far,
            });
        }
        self.outbox.push_back(SortPhase1Event::AllAnnounced {
            slot_count,
            memory_chunk_count,
            total_records,
        });
        Ok(())
    }

    /// The index of `slot`'s run, read back from its footer and brought to the
    /// final stride, or `None` once indexing is off or when the run has no
    /// footer (a run whose index could not be sparsified; the merge then
    /// reports the slot unindexed and runs serially).
    ///
    /// # Errors
    ///
    /// A footer that cannot be read or decoded, or whose key width is not
    /// the runs' key kind's (`InvalidData`).
    fn read_back_index(
        &self,
        slot: &fgumi_sort::SortMergeSlot,
        path: &std::path::Path,
    ) -> io::Result<Option<RunIndex>> {
        if self.indexing_off {
            return Ok(None);
        }
        let index = external_wait(|| slot.read_run_index()).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("SpillWrite: reading the run index of {}: {e}", path.display()),
            )
        })?;
        let kind = self.key_kind.expect("a closed spill run set the key kind");
        if let Some(ix) = &index
            && Some(ix.key_width()) != fgumi_sort::spill_key_width(kind)
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "SpillWrite: the run index of {} has {}-byte keys, but its runs are keyed \
                     {kind:?}",
                    path.display(),
                    ix.key_width()
                ),
            ));
        }
        Ok(index.and_then(|ix| sparsify_to(ix, self.index_stride)).map(|mut ix| {
            ix.shrink_to_fit();
            ix
        }))
    }

    /// Error if input drained with spill runs that were never announced. Runs are
    /// held until `AllAnnounced`, so without it they would never reach the merge
    /// and the sort would finish "successfully" without their records.
    fn ensure_runs_announced(&self) -> io::Result<()> {
        if !self.announced && self.runs_closed > 0 {
            return Err(io::Error::other(format!(
                "SpillWrite: input drained without AllAnnounced; {} spill runs were never \
                 announced to the merge",
                self.runs_closed
            )));
        }
        Ok(())
    }

    /// Error if a spill file is still open. A `Residual` / `AllAnnounced` event,
    /// or end-of-stream, while `current` holds an unterminated file means a spill
    /// lost its `is_last_in_file` block (a `SpillGather` framing bug) — failing
    /// loud avoids dropping a spill or publishing `AllAnnounced` before its
    /// `SpillReady`.
    fn ensure_no_open_file(&self, at: &str) -> io::Result<()> {
        if let Some(open) = &self.current {
            return Err(io::Error::other(format!(
                "SpillWrite: {at} arrived while spill file_id {} was still open \
                 (missing is_last_in_file block)",
                open.file_id
            )));
        }
        Ok(())
    }
}

impl Step for SpillWrite {
    type Input = SpillBlockEvent;
    type Outputs = Single<SortPhase1Event>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SpillWrite",
            // Detached (own thread) on the standalone-sort spill path; otherwise
            // the default pool-scheduled Serial + sticky writer. `sticky` is
            // irrelevant for Detached (it never enters a worker's worklist).
            kind: if self.detached { StepKind::Detached } else { StepKind::Serial },
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn detached_group(&self) -> DetachedGroup {
        // When detached (standalone-sort spill path), share the sort's I/O
        // writer driver thread with the terminal `WriteBgzfFile` — phase-1 spill
        // and phase-2 output writes are temporally disjoint (true N+2). Consulted
        // only when the step is Detached.
        DetachedGroup::Shared(crate::sort::SORT_IO_GROUP)
    }

    fn affinity(&self) -> Affinity {
        // Ignored for `Detached` (no pool worker drives it); kept for the
        // default Serial path where it pins the writer to the last worker.
        Affinity::Writer
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] = &[CounterSpec::new("spill_bytes_written", "bytes")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        // Produced events go out first, one per call.
        if let Some(out) = self.outbox.pop_front() {
            if let Err(unpushed) = ctx.outputs.push(out) {
                self.held.put(unpushed);
            }
            return Ok(StepOutcome::Progress);
        }

        // A job in flight blocks new input: the live-run count is at its limit
        // until the job completes, so upstream waits. Its completion unparks this
        // thread (`poll_registry`), whether it runs on a pool worker or a driver.
        if self.active.is_some() {
            return Ok(if self.poll_registry()? {
                StepOutcome::Progress
            } else {
                StepOutcome::NoProgress
            });
        }

        if let Some(event) = ctx.input.pop() {
            // A batch of blocks per call (see `write_batch`): one block per call
            // made the per-dispatch cost, not the disk, the phase-1 spill bound.
            let spill_bytes = self.write_batch(event, || ctx.input.pop())?;
            ctx.counters.add(SPILL_BYTES_WRITTEN, spill_bytes);
            return Ok(StepOutcome::Progress);
        }

        if ctx.input.is_drained() {
            // End-of-stream with a file still open means the final spill never
            // got its `is_last_in_file` block — fail loud rather than leave a
            // truncated, unterminated spill on disk.
            self.ensure_no_open_file("input drained")?;
            self.ensure_runs_announced()?;
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

#[cfg(test)]
mod tests;
