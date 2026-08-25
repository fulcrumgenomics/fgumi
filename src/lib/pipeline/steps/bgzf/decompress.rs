//! `BgzfDecompress` mid-step. `Parallel + ByItemOrdinal`. Decompresses
//! incoming `BgzfBlock`s into `DecompressedBlock`s using `libdeflater`'s
//! decompressor (one per worker via `Clone`).

use std::io;

use fgumi_bgzf::reader::decompress_block_slice_into_opts;
use libdeflater::Decompressor;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::{BgzfBlock, DecompressedBlock};

/// Per-worker decompression scratch capacity. Sized to the upper end
/// of typical BGZF block decompressed sizes (256 KiB ≈ 4 blocks × 64 KiB).
/// Mirrors legacy `pipeline/bam.rs:1560` `DECOMPRESSION_BUFFER_CAPACITY`.
///
/// Why a fixed size: when each `try_run` calls `Vec::with_capacity(N)` with
/// a *variable* `N` (= `block.uncompressed_size`), every allocation hits a
/// different mimalloc size class and defeats the thread-local cache. A
/// fixed capacity always pulls from the same size class so the freed
/// buffer that just left for downstream gets reused on the next call.
const DECOMPRESS_SCRATCH_CAPACITY: usize = 256 * 1024;

/// `Parallel + ByItemOrdinal` decompressor. Each worker holds its own
/// `libdeflater::Decompressor` (`Clone` constructs a fresh one) and a
/// fixed-capacity output scratch buffer that's emitted via `mem::replace`
/// to avoid variable-size allocations on the hot path.
pub struct BgzfDecompress {
    decompressor: Decompressor,
    /// Per-worker output scratch. `try_run` decompresses into this
    /// buffer, then swaps it out with a fresh fixed-capacity Vec via
    /// `mem::replace`. Keeps allocations on a single mimalloc size
    /// class so the thread-local cache can recycle them. Mirrors legacy
    /// `WorkerState::decompression_buffer`.
    output_scratch: Vec<u8>,
    held: HeldSlot<Unpushed<DecompressedBlock>>,
    output_byte_limit: u64,
    /// Whether to verify each block's CRC32 against its footer. The BAM source
    /// preamble threads the command's `--check-crc`/`--no-check-crc` policy (via
    /// `BamIoOptions::effective_check_crc`) here so the chain reproduces the
    /// non-chain path's CRC behavior; the decompressed-size check always runs.
    verify_crc: bool,
}

impl BgzfDecompress {
    /// Construct a decompressor that verifies CRC32 (the pre-existing default).
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self::new_with_crc(output_byte_limit, true)
    }

    /// Construct a decompressor with an explicit CRC-verification policy.
    ///
    /// The chain's BAM source preamble uses this to honor the command's
    /// `--check-crc` / `--no-check-crc` policy, matching the non-chain path.
    #[must_use]
    pub fn new_with_crc(output_byte_limit: u64, verify_crc: bool) -> Self {
        Self {
            decompressor: Decompressor::new(),
            output_scratch: Vec::with_capacity(DECOMPRESS_SCRATCH_CAPACITY),
            held: HeldSlot::new(),
            output_byte_limit,
            verify_crc,
        }
    }
}

impl Clone for BgzfDecompress {
    fn clone(&self) -> Self {
        Self {
            decompressor: Decompressor::new(),
            output_scratch: Vec::with_capacity(DECOMPRESS_SCRATCH_CAPACITY),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            verify_crc: self.verify_crc,
        }
    }
}

impl Step for BgzfDecompress {
    type Input = BgzfBlock;
    type Outputs = OrderedBytesSingle<DecompressedBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "BgzfDecompress",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    // `Contention` (not `NoProgress`) keeps the worker
                    // alive for retry — `NoProgress` would let the
                    // framework mark this worker `Skip` if input is also
                    // drained, silently dropping the held item.
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(block) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };

        // Decompress into the per-worker scratch buffer. After the
        // decompressor fills it, `mem::replace` swaps it out (the filled
        // buffer goes into the emitted block; a fresh fixed-capacity
        // buffer takes its place). Single mimalloc size class keeps the
        // thread-local cache hot.
        decompress_block_slice_into_opts(
            &block.bytes,
            &mut self.decompressor,
            &mut self.output_scratch,
            self.verify_crc,
        )?;
        let bytes = std::mem::replace(
            &mut self.output_scratch,
            Vec::with_capacity(DECOMPRESS_SCRATCH_CAPACITY),
        );

        let decompressed = DecompressedBlock { batch_serial: block.batch_serial, bytes };

        match ctx.outputs.push(decompressed) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_advertises_parallel_byordinal() {
        let s = BgzfDecompress::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "BgzfDecompress");
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn clone_constructs_fresh_decompressor() {
        let s = BgzfDecompress::new(1024);
        let _cloned = s.clone();
    }
}
