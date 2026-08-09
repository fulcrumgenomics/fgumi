//! `SpillBlockCompress` — middle step of the block-parallel spill-write split
//! (`SpillGather` → `SpillBlockCompress` → `SpillWrite`).
//!
//! `SpillBlockCompress` (`Parallel + ByItemOrdinal`) compresses each raw
//! [`SpillBlockEvent::Block`] from `SpillGather` into a self-contained
//! compressed unit (framed BGZF block(s) for bgzf, or a `[u32 len][zstd frame]`
//! for zstd) via the shared [`SpillBlockCompressor`] kernel, replacing the
//! `bytes` in place. `Residual` and `AllAnnounced` pass straight through. The
//! `ordinal` is preserved on every event, so the framework's `ByItemOrdinal`
//! reorder hands `SpillWrite` a dense, in-order stream regardless of which worker
//! compressed which block — exactly the output path's `BgzfCompress` idiom.
//!
//! Each `Parallel` worker holds its own [`SpillBlockCompressor`], built lazily on
//! the first block (the zstd compressor's construction is fallible, so it is
//! surfaced through `try_run`'s `io::Result` rather than the infallible
//! `new_worker_copy`).

use std::io;

use fgumi_sort::{SpillBlockCompressor, SpillCodec};

use crate::sort::protocol::SpillBlockEvent;
use fgumi_pipeline_core::{
    HeldRetry, Unpushed,
    held::HeldSlot,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// `Parallel + ByItemOrdinal` block compressor for the spill-write split.
///
/// Not to be confused with the similarly-named
/// [`CompressSpill`](super::CompressSpill): this `SpillBlockCompress` is the
/// **pure block-compression** middle step of the finer
/// `SpillGather → SpillBlockCompress → SpillWrite` split (the disk write is the
/// separate `SpillWrite` step), whereas `CompressSpill` is the **composite
/// compress-and-write-to-disk** step of the coarser
/// `SortBuffer → CompressSpill → SortSpillDecompress → SortMerge` chain.
pub struct SpillBlockCompress {
    codec: SpillCodec,
    compression: u32,
    /// Per-worker compressor, built lazily on the first block. `None` until then
    /// (and on fresh `new_worker_copy` clones).
    compressor: Option<SpillBlockCompressor>,
    held: HeldSlot<Unpushed<SpillBlockEvent>>,
    output_byte_limit: u64,
}

impl SpillBlockCompress {
    /// Build a `SpillBlockCompress` for `codec` at `compression`. `output_byte_limit`
    /// byte-bounds the compressed-block output queue.
    #[must_use]
    pub fn new(codec: SpillCodec, compression: u32, output_byte_limit: u64) -> Self {
        Self { codec, compression, compressor: None, held: HeldSlot::new(), output_byte_limit }
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        !matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld)
    }

    /// Compress one event's payload (the `Block` arm) or pass it through. The
    /// `ordinal` and routing fields are preserved. `StepCtx`-free for unit tests.
    ///
    /// # Errors
    ///
    /// Propagates compressor-init or compression errors.
    fn compress_event(&mut self, event: SpillBlockEvent) -> io::Result<SpillBlockEvent> {
        match event {
            SpillBlockEvent::Block {
                ordinal,
                file_id,
                is_last_in_file,
                records_ingested_so_far,
                bytes,
            } => {
                if self.compressor.is_none() {
                    self.compressor =
                        Some(SpillBlockCompressor::new(self.codec, self.compression)?);
                }
                let compressor = self.compressor.as_mut().expect("compressor built above");
                let compressed = compressor.compress_block(&bytes)?;
                Ok(SpillBlockEvent::Block {
                    ordinal,
                    file_id,
                    is_last_in_file,
                    records_ingested_so_far,
                    bytes: compressed,
                })
            }
            // Passthrough variants carry no compressible payload.
            other @ (SpillBlockEvent::Residual { .. } | SpillBlockEvent::AllAnnounced { .. }) => {
                Ok(other)
            }
        }
    }
}

impl Clone for SpillBlockCompress {
    fn clone(&self) -> Self {
        // Fresh per-worker compressor + held slot; shared config copied.
        Self {
            codec: self.codec,
            compression: self.compression,
            compressor: None,
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
        }
    }
}

impl Step for SpillBlockCompress {
    type Input = SpillBlockEvent;
    type Outputs = OrderedBytesSingle<SpillBlockEvent>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SpillBlockCompress",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        if let Some(event) = ctx.input.pop() {
            let forwarded = self.compress_event(event)?;
            if let Err(unpushed) = ctx.outputs.push(forwarded) {
                self.held.put(unpushed);
            }
            return Ok(StepOutcome::Progress);
        }

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
