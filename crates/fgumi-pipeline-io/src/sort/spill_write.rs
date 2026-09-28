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
//! runs is merged into one with a [`RunMergerDyn`]. The merge is driven
//! **cooperatively** — a bounded batch of records per `try_run`, returning
//! `Progress` — because a `Detached` step must never block inside `try_run`, and a
//! step that stops reporting progress for long enough is cancelled as wedged.
//! While a merge is in flight no input is taken, so upstream back-pressures while
//! the live-run count is being brought back under the limit.
//!
//! This step owns the `TmpDirAllocator` (Serial ⇒ the pick is uncontended) and
//! the RAII temp-dir handles, matching the retired `CompressSpill`'s lifetime
//! contract.

use std::collections::VecDeque;
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use fgumi_bam_io::ProgressTracker;
use fgumi_bgzf::BGZF_MAX_BLOCK_SIZE;
use fgumi_sort::{
    RunMergeProgress, RunMergeSpec, RunMergerDyn, SpillCodec, SpillKeyKind, TmpDirAllocator,
    new_run_merger, spill_magic, spill_trailer,
};
use log::{info, warn};
use parking_lot::Mutex;
use tempfile::TempDir;

use crate::sort::protocol::{SortPhase1Event, SpillBlockEvent};
use crate::sort::run_stack::{RunStack, SpillRun};
use fgumi_pipeline_core::{
    HeldRetry, Unpushed,
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

/// Counter slot index: records rewritten by consolidation this call.
const CONSOLIDATION_RECORDS: usize = 1;

/// Records merged per `try_run` while a consolidation is in flight. Large enough
/// that per-call overhead is negligible for short reads.
const CONSOLIDATION_RECORDS_PER_CALL: usize = 16 * 1024;

/// Record bytes merged per `try_run`, whichever of the two budgets is reached
/// first. Long reads make a record budget alone unbounded in time (16 Ki records
/// of 100 kb reads is ~1.6 GB to decompress and recompress), and a step that
/// reports no progress for long enough is cancelled as wedged.
const CONSOLIDATION_BYTES_PER_CALL: u64 = 8 * 1024 * 1024;

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

    /// Time spent merging runs during consolidation, in seconds.
    #[must_use]
    #[allow(clippy::cast_precision_loss)] // display value; > 2^52 ns is ~52 days
    pub fn consolidation_secs(&self) -> f64 {
        self.consolidation_nanos.load(Ordering::Relaxed) as f64 / 1e9
    }

    /// Time spent merging runs during consolidation, in nanoseconds.
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

    /// Count one consolidation that spent `busy` merging.
    pub fn record_consolidation(&self, busy: Duration) {
        self.consolidations.fetch_add(1, Ordering::Relaxed);
        self.consolidation_nanos
            .fetch_add(u64::try_from(busy.as_nanos()).unwrap_or(u64::MAX), Ordering::Relaxed);
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
}

/// A consolidation in flight.
struct ActiveMerge {
    /// Range of `runs` being merged. No run is pushed while a merge is in flight,
    /// so the range stays valid until the merge completes.
    range: std::ops::Range<usize>,
    merger: Box<dyn RunMergerDyn>,
    output: PathBuf,
    /// Time spent inside the merger's `step` calls (not the time between them).
    busy: Duration,
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
    /// Compression level consolidated runs are written at (the spill level).
    compression: u32,
    /// The currently-open spill file, if any.
    current: Option<OpenSpill>,
    /// Closed runs awaiting the merge, and the consolidation policy over them.
    runs: RunStack,
    /// Key type of every run so far; all runs of one sort share it.
    key_kind: Option<SpillKeyKind>,
    /// The consolidation in flight, if any.
    active: Option<ActiveMerge>,
    /// Consolidations started, for naming merged files.
    merges_started: usize,
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
    /// while being read by `SortMerge`. Matches `CompressSpill`'s lifetime.
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
}

impl SpillWrite {
    /// Build a `SpillWrite`. `alloc` names spill files across the configured temp
    /// dirs; `codec` selects the on-disk format; `temp_dirs` holds the RAII
    /// handles alive for the step's lifetime. `output_byte_limit` byte-bounds the
    /// forwarded-event output queue.
    ///
    /// Consolidation is off until [`Self::with_max_temp_files`] is applied.
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
            compression: 1,
            current: None,
            runs: RunStack::new(0),
            key_kind: None,
            active: None,
            merges_started: 0,
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
        }
    }

    /// Keep at most `max_temp_files` spill runs live, consolidating runs written
    /// at `compression` (the spill compression level) when the limit is reached.
    /// A limit below 2 leaves consolidation off.
    #[must_use]
    pub fn with_max_temp_files(mut self, max_temp_files: usize, compression: u32) -> Self {
        self.runs = RunStack::new(max_temp_files);
        self.compression = compression;
        self
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
        let file = OpenOptions::new().write(true).create_new(true).open(&path)?;
        let mut writer = BufWriter::with_capacity(256 * 1024, file);
        writer.write_all(spill_magic(self.codec))?;
        Ok(OpenSpill { file_id, key_kind, path, writer })
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
    /// consolidation (see [`Self::advance_consolidation`]). `StepCtx`-free for
    /// unit testing.
    ///
    /// # Errors
    ///
    /// Propagates file-create / write / slot-open / merger-open errors. Also
    /// errors if a block arrives for a different `file_id` than the open file
    /// while one is open without an intervening `is_last_in_file` — a
    /// framework-ordering invariant violation that must fail loud rather than
    /// corrupt a spill — or if runs of one sort disagree on their key type.
    fn process_event(&mut self, event: SpillBlockEvent) -> io::Result<()> {
        if self.announced {
            // The runs were handed downstream at `AllAnnounced`; anything after
            // it would be written to a stack nobody reads.
            return Err(io::Error::other(
                "SpillWrite: input arrived after AllAnnounced (the run set was already announced)",
            ));
        }
        match event {
            SpillBlockEvent::Block {
                file_id,
                key_kind,
                is_last_in_file,
                records_ingested_so_far,
                bytes,
                ..
            } => {
                // Open the file on its first block; otherwise the open file must
                // match (blocks for one file are contiguous in the ordinal stream).
                if self.current.is_none() {
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
                open.writer.write_all(&bytes)?;
                self.spill_progress.log_if_needed(bytes.len() as u64);

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
                slot_count, memory_chunk_count, total_records, ..
            } => {
                self.ensure_no_open_file("AllAnnounced")?;
                self.announce_runs(slot_count, memory_chunk_count, total_records)
            }
        }
    }

    /// Finish the open spill file, record it as a live run, and start a
    /// consolidation if the policy now calls for one.
    fn close_run(&mut self, records_ingested_so_far: u64) -> io::Result<()> {
        let OpenSpill { file_id, key_kind, path, mut writer } =
            self.current.take().expect("open file present");
        writer.write_all(spill_trailer(self.codec))?;
        writer.flush()?;
        drop(writer);
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
        let bytes = std::fs::metadata(&path)?.len();
        self.runs.push(SpillRun { file_id, path, bytes, records_ingested_so_far });
        self.start_consolidation_if_needed()
    }

    /// Start the consolidation the run stack asks for, if any and none is in
    /// flight.
    fn start_consolidation_if_needed(&mut self) -> io::Result<()> {
        if self.active.is_some() {
            return Ok(());
        }
        let Some(range) = self.runs.next_merge() else { return Ok(()) };
        let key_kind = self.key_kind.expect("a live run implies a known key kind");
        let inputs: Vec<PathBuf> =
            self.runs.runs()[range.clone()].iter().map(|r| r.path.clone()).collect();
        let input_bytes: u64 = self.runs.runs()[range.clone()].iter().map(|r| r.bytes).sum();
        let output = self.next_path(&format!("merged_{:04}.keyed", self.merges_started))?;
        self.merges_started += 1;
        info!(
            "Consolidating {} spill runs ({input_bytes} bytes) into 1; {} live runs",
            inputs.len(),
            self.runs.runs().len()
        );
        let merger = new_run_merger(
            key_kind,
            &RunMergeSpec {
                inputs: &inputs,
                output: &output,
                codec: self.codec,
                compression: self.compression,
                block_size: BGZF_MAX_BLOCK_SIZE,
            },
        )?;
        self.active = Some(ActiveMerge { range, merger, output, busy: Duration::ZERO });
        Ok(())
    }

    /// Merge the next bounded batch of the in-flight consolidation. When it
    /// completes, replace its inputs with the merged run, delete the inputs, and
    /// start the next consolidation the policy asks for. Returns the records
    /// merged by this call; `0` when no consolidation is in flight.
    ///
    /// # Errors
    ///
    /// Propagates merge (read / decompress / compress / write) errors.
    fn advance_consolidation(&mut self, max_records: usize, max_bytes: u64) -> io::Result<u64> {
        let Some(active) = self.active.as_mut() else { return Ok(0) };
        let before_records = active.merger.records_written();
        let before = Instant::now();
        let progress = active.merger.step(max_records, max_bytes)?;
        active.busy += before.elapsed();
        let merged = active.merger.records_written() - before_records;
        if progress == RunMergeProgress::Working {
            return Ok(merged);
        }
        let ActiveMerge { range, output, busy, .. } =
            self.active.take().expect("active merge present");
        let merged_bytes = std::fs::metadata(&output)?.len();
        let absorbed = self.runs.complete_merge(range, output, merged_bytes);
        for run in &absorbed {
            if let Err(e) = std::fs::remove_file(&run.path) {
                // The temp dir removes it at exit regardless; only disk headroom
                // is at stake until then.
                warn!("SpillWrite: could not remove consolidated run {}: {e}", run.path.display());
            }
        }
        self.stats.record_consolidation(busy);
        self.start_consolidation_if_needed()?;
        Ok(merged)
    }

    /// On `AllAnnounced`: open a merge slot for every surviving run and queue its
    /// `SpillReady` (in `file_id` order), then `AllAnnounced` with the surviving
    /// run count.
    ///
    /// # Errors
    ///
    /// Errors if a slot cannot be opened, or if the upstream run count disagrees
    /// with the runs this step closed (a lost or duplicated run).
    fn announce_runs(
        &mut self,
        runs_announced: u32,
        memory_chunk_count: u32,
        total_records: u64,
    ) -> io::Result<()> {
        debug_assert!(self.active.is_none(), "no input is taken while a merge is in flight");
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
        for SpillRun { file_id, path, records_ingested_so_far, .. } in runs {
            let slot = fgumi_sort::open_spill_slot(&path, file_id).map_err(|e| {
                io::Error::other(format!(
                    "SpillWrite: failed to open spill slot {}: {e:#}",
                    path.display()
                ))
            })?;
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
        const SPECS: &[CounterSpec] = &[
            CounterSpec::new("spill_bytes_written", "bytes"),
            CounterSpec::new("consolidation_records", "records"),
        ];
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

        // An in-flight consolidation takes priority over new input: the live-run
        // count is at its limit until it finishes, so upstream waits (the same
        // back-pressure `RawExternalSorter`'s synchronous consolidation applies),
        // while each call still makes and reports progress.
        if self.active.is_some() {
            let merged = self.advance_consolidation(
                CONSOLIDATION_RECORDS_PER_CALL,
                CONSOLIDATION_BYTES_PER_CALL,
            )?;
            ctx.counters.add(CONSOLIDATION_RECORDS, merged);
            return Ok(StepOutcome::Progress);
        }

        if let Some(event) = ctx.input.pop() {
            // One event per `try_run`, so this is already the batch total for
            // this call. Only `Block` carries spill bytes; read its length
            // before `process_event` consumes the event.
            let spill_bytes = match &event {
                SpillBlockEvent::Block { bytes, .. } => bytes.len() as u64,
                SpillBlockEvent::Residual { .. } | SpillBlockEvent::AllAnnounced { .. } => 0,
            };
            self.process_event(event)?;
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
