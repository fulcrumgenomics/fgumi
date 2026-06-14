//! `BgzfCompress` mid-step. `Parallel + ByItemOrdinal`. Compresses each
//! `DecompressedBlock` into one logical `BgzfBlock` whose `bytes` field
//! contains the concatenation of however many physical 64-KiB BGZF blocks
//! `libdeflater` produces internally.
//!
//! ## 1:1 input-to-output mapping (Phase 3 design)
//!
//! Each `try_run` consumes exactly one input and emits exactly one output
//! with `batch_serial = input.batch_serial`. This preserves the
//! consecutive-ordinal invariant required by the downstream
//! `ReorderStage` without any sub-serial encoding.
//!
//! Internally, `InlineBgzfCompressor` may produce multiple physical BGZF
//! blocks when the input exceeds `BGZF_MAX_BLOCK_SIZE = 64 KiB`. We
//! concatenate those physical blocks into a single `BgzfBlock`'s bytes —
//! a valid BGZF stream is just a concatenation of independent blocks, so
//! the downstream writer can emit them verbatim.

use std::io;

use fgumi_bgzf::InlineBgzfCompressor;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::{BgzfBlock, DecompressedBlock};

/// Per-worker BGZF compressed-output scratch capacity. The compressed
/// output of `BgzfCompress` is the concatenation of physical 64 KiB BGZF
/// blocks; a single batch typically produces 1-4 such blocks. Sized to
/// the typical case so the same mimalloc size class is hit every call.
/// Matches legacy `bam.rs:1561` `SERIALIZATION_BUFFER_CAPACITY` × 4.
const COMPRESS_SCRATCH_CAPACITY: usize = 256 * 1024;

/// `Parallel + ByItemOrdinal` BGZF compressor.
pub struct BgzfCompress {
    /// Per-worker compressor (each worker holds its own clone).
    compressor: InlineBgzfCompressor,
    compression_level: u32,
    /// Per-worker output scratch buffer; `mem::replace`d on each emit so
    /// the freshly-allocated replacement is always the same size class.
    /// See `BgzfDecompress::output_scratch` for rationale.
    output_scratch: Vec<u8>,
    held: HeldSlot<Unpushed<BgzfBlock>>,
    output_byte_limit: u64,
}

impl BgzfCompress {
    #[must_use]
    pub fn new(compression_level: u32, output_byte_limit: u64) -> Self {
        Self {
            compressor: InlineBgzfCompressor::new(compression_level),
            compression_level,
            output_scratch: Vec::with_capacity(COMPRESS_SCRATCH_CAPACITY),
            held: HeldSlot::new(),
            output_byte_limit,
        }
    }
}

impl Clone for BgzfCompress {
    fn clone(&self) -> Self {
        Self {
            compressor: InlineBgzfCompressor::new(self.compression_level),
            compression_level: self.compression_level,
            output_scratch: Vec::with_capacity(COMPRESS_SCRATCH_CAPACITY),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
        }
    }
}

impl Step for BgzfCompress {
    type Input = DecompressedBlock;
    type Outputs = OrderedBytesSingle<BgzfBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "BgzfCompress",
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
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(parent) = ctx.input.pop() else {
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
        let DecompressedBlock { batch_serial, bytes } = parent;
        let uncompressed_size = u32::try_from(bytes.len()).unwrap_or(u32::MAX);

        // Compress and harvest physical BGZF blocks.
        self.compressor.write_all(&bytes)?;
        self.compressor.flush()?;
        let children = self.compressor.take_blocks();
        // INVARIANT: `child.serial` here is `InlineBgzfCompressor`'s
        // monotonic per-worker counter; it has no relation to the
        // pipeline-wide `batch_serial`. We deliberately drop it (read
        // only `child.data`) and inherit the parent's `batch_serial` for
        // the emitted `BgzfBlock`. Future code must NOT propagate
        // `child.serial` into the framework's ordinal stream.

        // Concatenate the physical block bytes into the per-worker
        // scratch, then mem::replace it out — single mimalloc size class
        // every call.
        for child in children {
            self.output_scratch.extend_from_slice(&child.data);
        }
        let bytes = std::mem::replace(
            &mut self.output_scratch,
            Vec::with_capacity(COMPRESS_SCRATCH_CAPACITY),
        );

        let out = BgzfBlock { batch_serial, bytes, uncompressed_size };
        match ctx.outputs.push(out) {
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
        let s = BgzfCompress::new(1, 1024);
        let p = s.profile();
        assert_eq!(p.name, "BgzfCompress");
        assert_eq!(p.kind, StepKind::Parallel);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn clone_constructs_fresh_compressor() {
        let s = BgzfCompress::new(1, 1024);
        let _cloned = s.clone();
    }

    #[test]
    fn small_input_produces_single_concatenated_output() {
        let mut s = BgzfCompress::new(1, 1024);
        s.compressor.write_all(b"hello bgzf").unwrap();
        s.compressor.flush().unwrap();
        let blocks = s.compressor.take_blocks();
        assert!(!blocks.is_empty(), "expected at least one BGZF block");
        // Each block is a valid BGZF block (starts with gzip magic 0x1f 0x8b).
        for b in &blocks {
            assert_eq!(&b.data[..2], &[0x1f, 0x8b], "BGZF magic");
        }
    }
}
