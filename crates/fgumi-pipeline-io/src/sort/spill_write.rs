//! `SpillWrite` — final step of the block-parallel spill-write split
//! (`SpillGather` → `SpillCompress` → `SpillWrite`).
//!
//! `SpillWrite` (`Serial + Affinity::Writer`) receives the compressed
//! [`SpillBlockEvent`]s in dense `ordinal` order (the framework's `ByItemOrdinal`
//! reorder feeds it like `WriteBgzfFile`), demultiplexes `Block`s back to
//! per-`file_id` spill files, and emits the existing [`SortPhase1Event`] so
//! `SortSpillDecompress` / `SortMerge` are unchanged.
//!
//! Because `SortBuffer` (Serial) emits spill chunks one-at-a-time in `seq` order,
//! `SpillGather` (Serial) fans them in order, and the reorder preserves that
//! order, each file's blocks arrive **contiguously** — so `SpillWrite` only ever
//! holds **one** open spill file at a time (`current`). It opens the file on the
//! first block (writing the codec magic), appends each compressed block, and on
//! `is_last_in_file` writes the codec trailer, opens the merge slot, and emits
//! `SpillReady`. This step owns the `TmpDirAllocator` (Serial ⇒ the pick is
//! uncontended) and the RAII temp-dir handles, matching the retired
//! `CompressSpill`'s lifetime contract.

use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use fgumi_bam_io::ProgressTracker;
use fgumi_sort::{SpillCodec, TmpDirAllocator, spill_magic, spill_trailer};
use parking_lot::Mutex;
use tempfile::TempDir;

use crate::sort::protocol::{SortPhase1Event, SpillBlockEvent};
use fgumi_pipeline_core::{
    HeldRetry, Unpushed,
    held::HeldSlot,
    outputs::Single,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Affinity, DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// The one spill file currently being written (open from its first block until
/// its `is_last_in_file` block).
struct OpenSpill {
    file_id: u32,
    path: PathBuf,
    writer: BufWriter<File>,
}

/// `Serial + Affinity::Writer` step that writes per-`file_id` spill files from
/// the compressed block stream and emits `SortPhase1Event`s.
pub struct SpillWrite {
    /// Shared temp-directory allocator (free-space-aware round-robin). `Serial`,
    /// so the lock is effectively uncontended (one writer worker).
    alloc: Arc<Mutex<TmpDirAllocator>>,
    /// Spill codec for chunk files (bgzf or zstd).
    codec: SpillCodec,
    /// The currently-open spill file, if any.
    current: Option<OpenSpill>,
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
            held: HeldSlot::new(),
            output_byte_limit,
            spill_progress: ProgressTracker::new("Spill bytes written")
                .with_interval(256 * 1024 * 1024),
            temp_dirs,
            detached: false,
        }
    }

    /// Run this spill writer on its own dedicated `StepKind::Detached` thread
    /// instead of as a pool-scheduled `Serial + Affinity::Writer` step. Used
    /// ONLY on the standalone-sort spill path (the Phase-1 analogue of Lever 2's
    /// detached terminal writer): it frees a pool worker for the
    /// compression-bound `SpillCompress` work, matching feat-runall's dedicated
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
    fn open_file(&self, file_id: u32) -> io::Result<OpenSpill> {
        let base = self.alloc.lock().next().map_err(|e| {
            io::Error::other(format!("SpillWrite: temp-dir allocation failed: {e:#}"))
        })?;
        let path = base.join(format!("chunk_{file_id:04}.keyed"));
        // `create_new` fails closed on a duplicate/stale path: a reused `file_id`
        // (or a leftover file) must surface as an error rather than truncate an
        // existing spill and silently corrupt merge input.
        let file = OpenOptions::new().write(true).create_new(true).open(&path)?;
        let mut writer = BufWriter::with_capacity(256 * 1024, file);
        writer.write_all(spill_magic(self.codec))?;
        Ok(OpenSpill { file_id, path, writer })
    }

    /// Process one input event, performing any disk writes and returning the
    /// `SortPhase1Event` to emit (a `Block` only emits on `is_last_in_file`).
    /// `StepCtx`-free for unit testing.
    ///
    /// # Errors
    ///
    /// Propagates file-create / write / slot-open errors. Also errors if a block
    /// arrives for a different `file_id` than the open file while one is open
    /// without an intervening `is_last_in_file` — a framework-ordering invariant
    /// violation that must fail loud rather than corrupt a spill.
    fn process_event(&mut self, event: SpillBlockEvent) -> io::Result<Option<SortPhase1Event>> {
        match event {
            SpillBlockEvent::Block {
                file_id,
                is_last_in_file,
                records_ingested_so_far,
                bytes,
                ..
            } => {
                // Open the file on its first block; otherwise the open file must
                // match (blocks for one file are contiguous in the ordinal stream).
                if self.current.is_none() {
                    self.current = Some(self.open_file(file_id)?);
                }
                let open = self.current.as_mut().expect("open file set above");
                if open.file_id != file_id {
                    return Err(io::Error::other(format!(
                        "SpillWrite: block for file_id {file_id} arrived while file_id {} \
                         was still open (blocks must be contiguous per file)",
                        open.file_id
                    )));
                }
                open.writer.write_all(&bytes)?;
                self.spill_progress.log_if_needed(bytes.len() as u64);

                if is_last_in_file {
                    let OpenSpill { file_id, path, mut writer } =
                        self.current.take().expect("open file present");
                    writer.write_all(spill_trailer(self.codec))?;
                    writer.flush()?;
                    drop(writer); // close the fd before opening the read slot
                    let slot = fgumi_sort::open_spill_slot(&path, file_id).map_err(|e| {
                        io::Error::other(format!(
                            "SpillWrite: failed to open spill slot {}: {e:#}",
                            path.display()
                        ))
                    })?;
                    Ok(Some(SortPhase1Event::SpillReady { slot, path, records_ingested_so_far }))
                } else {
                    Ok(None)
                }
            }
            SpillBlockEvent::Residual { chunk, records_ingested_so_far, .. } => {
                self.ensure_no_open_file("residual")?;
                // Wrap in a fresh, uniquely-owned `Arc`: the chunk is only ever
                // moved (never cloned) onward, so `SortMerge`'s `Arc::try_unwrap`
                // invariant holds.
                Ok(Some(SortPhase1Event::MemoryChunk {
                    chunk: Arc::new(chunk),
                    records_ingested_so_far,
                }))
            }
            SpillBlockEvent::AllAnnounced {
                slot_count, memory_chunk_count, total_records, ..
            } => {
                self.ensure_no_open_file("AllAnnounced")?;
                Ok(Some(SortPhase1Event::AllAnnounced {
                    slot_count,
                    memory_chunk_count,
                    total_records,
                }))
            }
        }
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

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        if let Some(event) = ctx.input.pop() {
            if let Some(out) = self.process_event(event)? {
                if let Err(unpushed) = ctx.outputs.push(out) {
                    self.held.put(unpushed);
                }
            }
            return Ok(StepOutcome::Progress);
        }

        if ctx.input.is_drained() {
            // End-of-stream with a file still open means the final spill never
            // got its `is_last_in_file` block — fail loud rather than leave a
            // truncated, unterminated spill on disk.
            self.ensure_no_open_file("input drained")?;
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

#[cfg(test)]
mod tests;
