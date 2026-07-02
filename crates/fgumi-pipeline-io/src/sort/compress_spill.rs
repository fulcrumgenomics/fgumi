//! `CompressSpill` — the second step of the P6 Phase-1 split (`SortBuffer` →
//! `CompressSpill` → `SortSpillDecompress` → `SortMerge`).
//!
//! `SortBuffer` (Serial) emits already-sorted chunks; `CompressSpill`
//! (`Parallel`) compresses each spill chunk to disk **inline on its framework
//! worker** (via [`fgumi_sort::write_sorted_chunk_inmem`], or the chunk's own
//! `write_spill` for the template-arena variant, retiring the private
//! `SortWorkerPool` compress path) and forwards the result as the existing
//! [`SortPhase1Event`], so `SortSpillDecompress` / `SortMerge` are unchanged.
//!
//! # Why a `Parallel` step is safe here
//!
//! `SortMerge` collects setup events and gates on **counts** (`slot_count` /
//! `memory_chunk_count` from `AllAnnounced`), not on event arrival order, so
//! multiple `CompressSpill` workers may emit `SpillReady` / `MemoryChunk` events
//! in any order. The one ordering-sensitive concern — the `LoserTree` tie-break
//! for equal sort keys — is handled by stamping each spill slot's `file_id` with
//! the chunk's **logical** spill index (`SortChunkEvent::Spill::seq`, assigned by
//! `SortBuffer`), so the tie-break is independent of which worker writes first.

use std::io;
use std::path::Path;
use std::sync::Arc;

use fgumi_sort::{SpillCodec, TmpDirAllocator};
use parking_lot::Mutex;
use tempfile::TempDir;

use crate::sort::protocol::{MemoryChunkErased, SortChunkEvent, SortPhase1Event};
use fgumi_pipeline_core::{
    HeldRetry, Unpushed,
    held::HeldSlot,
    outputs::Single,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// `Parallel` step that compresses sorted spill chunks to disk and forwards
/// residual chunks, emitting [`SortPhase1Event`]s to `SortSpillDecompress`.
///
/// Not to be confused with the similarly-named
/// [`SpillCompress`](super::SpillCompress): this `CompressSpill` is the
/// **composite compress-and-write-to-disk** step of the coarser
/// `SortBuffer → CompressSpill → SortSpillDecompress → SortMerge` chain, whereas
/// `SpillCompress` is the **pure block-compression** middle step of the finer
/// `SpillGather → SpillCompress → SpillWrite` split (where the disk write is a
/// separate `SpillWrite` step).
pub struct CompressSpill {
    /// Shared temp-directory allocator (free-space-aware round-robin). Behind a
    /// `Mutex` because the step is `Parallel`; the lock is held only for the
    /// brief base-directory pick, never across the (expensive) compress+write.
    alloc: Arc<Mutex<TmpDirAllocator>>,
    /// Spill codec for chunk files (bgzf or zstd).
    codec: SpillCodec,
    /// Temp-file compression level (`0` = uncompressed bgzf; zstd level for zstd).
    compression: u32,
    /// At most one not-yet-pushable output event, parked on downstream
    /// backpressure. This is the framework's standard backpressure idiom (there
    /// is no peek-before-pop), identical to the sibling Parallel step
    /// [`SortSpillDecompress`](super::SortSpillDecompress), and like it this one
    /// held event sits **outside** `output_byte_limit`: total retained memory is
    /// `queue_bytes + (≤1 event) × workers`. For a `SpillReady` the held event is
    /// tiny (an `Arc<SortMergeSlot>` + path); only a `MemoryChunk` (residual)
    /// holds records, and the residual must transit the pipeline occupying
    /// ~`memory_limit` regardless of whether it sits in the queue or this slot, so
    /// the held slot adds no peak beyond what the residual already costs.
    held: HeldSlot<Unpushed<SortPhase1Event>>,
    output_byte_limit: u64,
    /// RAII temp-dir handles, shared across `Parallel` clones. Held for the
    /// step's lifetime so spill files survive while being written; on the last
    /// clone's drop (after the step finishes — i.e. every spill is written and
    /// every slot has an open fd) the dirs are removed. `SortMerge` then reads
    /// each slot via its already-open fd, matching the legacy `SortAndSpill`
    /// unlink-after-emit lifetime (correct on the Unix targets). Empty in tests
    /// that hold their own `TempDir`. Held purely for RAII (`Arc`-cloned to each
    /// worker), never read for its value.
    #[allow(dead_code)]
    temp_dirs: Arc<Vec<TempDir>>,
}

impl CompressSpill {
    /// Build a `CompressSpill` step.
    ///
    /// `alloc` names spill files across the configured temp dirs; `codec` /
    /// `compression` select the on-disk spill format. `temp_dirs` holds the RAII
    /// handles for those dirs alive for the step's lifetime. `output_byte_limit`
    /// byte-bounds the forwarded-event output queue (its `MemoryChunk` variant
    /// retains sorted records, so the queue must budget on bytes, not count).
    #[must_use]
    pub fn new(
        alloc: Arc<Mutex<TmpDirAllocator>>,
        codec: SpillCodec,
        compression: u32,
        output_byte_limit: u64,
        temp_dirs: Arc<Vec<TempDir>>,
    ) -> Self {
        Self { alloc, codec, compression, held: HeldSlot::new(), output_byte_limit, temp_dirs }
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        // `true` once the slot is clear (was empty, or the held event flushed);
        // `false` while it's still held under backpressure. Uses the canonical
        // re-hold helper so the put-back-on-reject invariant lives in one place.
        !matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld)
    }

    /// Allocate a spill path for chunk `seq`. Names the file by the logical spill
    /// index so the path is unique without a shared counter, then draws a base
    /// directory from the shared allocator (the only locked section).
    fn spill_path(&self, seq: u32) -> io::Result<std::path::PathBuf> {
        let base = self.alloc.lock().next().map_err(|e| {
            io::Error::other(format!("CompressSpill: temp-dir allocation failed: {e:#}"))
        })?;
        Ok(base.join(format!("chunk_{seq:04}.keyed")))
    }

    /// Compress one input event into the forwarded [`SortPhase1Event`]. The
    /// `Spill` arm does the file write inline; the others are passthroughs. This
    /// is `StepCtx`-free so it is unit-testable on synthetic chunks.
    fn compress_event(&self, event: SortChunkEvent) -> io::Result<SortPhase1Event> {
        match event {
            SortChunkEvent::Spill { seq, chunk, records_ingested_so_far } => {
                let path = self.spill_path(seq)?;
                write_chunk(&chunk, &path, self.codec, self.compression)?;
                let slot = fgumi_sort::open_spill_slot(&path, seq).map_err(|e| {
                    io::Error::other(format!(
                        "CompressSpill: failed to open spill slot {}: {e:#}",
                        path.display()
                    ))
                })?;
                Ok(SortPhase1Event::SpillReady { slot, path, records_ingested_so_far })
            }
            SortChunkEvent::Residual { chunk, records_ingested_so_far } => {
                // Wrap in a fresh, uniquely-owned `Arc`: the chunk is only ever
                // moved (never cloned) onward, so `SortMerge`'s `Arc::try_unwrap`
                // invariant holds.
                Ok(SortPhase1Event::MemoryChunk { chunk: Arc::new(chunk), records_ingested_so_far })
            }
            SortChunkEvent::AllAnnounced { slot_count, memory_chunk_count, total_records } => {
                Ok(SortPhase1Event::AllAnnounced { slot_count, memory_chunk_count, total_records })
            }
        }
    }
}

/// Dispatch a sorted [`MemoryChunkErased`] to [`fgumi_sort::write_sorted_chunk_inmem`]
/// (or the chunk's own `write_spill`, for the template-arena variant) for the
/// concrete key variant.
fn write_chunk(
    chunk: &MemoryChunkErased,
    path: &Path,
    codec: SpillCodec,
    compression: u32,
) -> io::Result<()> {
    let result = match chunk {
        MemoryChunkErased::Coordinate(c) => {
            fgumi_sort::write_sorted_chunk_inmem(path, codec, compression, c)
        }
        MemoryChunkErased::QuerynameLex(c) => {
            fgumi_sort::write_sorted_chunk_inmem(path, codec, compression, c)
        }
        MemoryChunkErased::QuerynameNatural(c) => {
            fgumi_sort::write_sorted_chunk_inmem(path, codec, compression, c)
        }
        MemoryChunkErased::TemplateCoordinate(c) => c.write_spill(path, codec, compression),
    };
    result.map_err(|e| {
        io::Error::other(format!("CompressSpill: chunk write to {} failed: {e:#}", path.display()))
    })
}

impl Clone for CompressSpill {
    fn clone(&self) -> Self {
        Self {
            alloc: Arc::clone(&self.alloc),
            codec: self.codec,
            compression: self.compression,
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            temp_dirs: Arc::clone(&self.temp_dirs),
        }
    }
}

impl Step for CompressSpill {
    type Input = SortChunkEvent;
    type Outputs = Single<SortPhase1Event>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "CompressSpill",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain held output first (at most one event is ever held, since we
        //    only pop a new input after `flush_held` clears the slot).
        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        // 2. Pop one input event, compress/forward it, hold on a full output.
        if let Some(event) = ctx.input.pop() {
            let forwarded = self.compress_event(event)?;
            if let Err(unpushed) = ctx.outputs.push(forwarded) {
                self.held.put(unpushed);
            }
            return Ok(StepOutcome::Progress);
        }

        // 3. No input available.
        if ctx.input.is_drained() {
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

#[cfg(test)]
mod tests;
