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

impl BgzfCompress {
    /// Drain the compressor's physical BGZF blocks and assemble the single
    /// output `BgzfBlock` byte buffer.
    ///
    /// INVARIANT: each block's `serial` is `InlineBgzfCompressor`'s monotonic
    /// per-worker counter; it has no relation to the pipeline-wide
    /// `batch_serial`. We deliberately drop it (read only `child.data`) and
    /// inherit the parent's `batch_serial` for the emitted `BgzfBlock`. Future
    /// code must NOT propagate `child.serial` into the framework's ordinal stream.
    ///
    /// The output is the concatenation of the physical block bytes (a valid BGZF
    /// stream is just independent blocks back to back). When there is exactly one
    /// physical block — the common case, since a batch usually fits in one 64 KiB
    /// block — its buffer is moved out directly (no scratch allocation and no
    /// concatenating memcpy); this is byte-identical to concatenating a single
    /// block. For the multi-block case the bytes are concatenated into the
    /// per-worker scratch (then `mem::replace`d out so the replacement stays in a
    /// fixed mimalloc size class), and each drained child buffer is handed back to
    /// the compressor's pool so the next compression reuses it.
    fn assemble_output_bytes(&mut self) -> Vec<u8> {
        let mut children = self.compressor.take_blocks();
        match children.len() {
            // No physical blocks — reachable from the all-clean rejects path,
            // where an empty `DecompressedBlock` produces no BGZF output. Return
            // an empty `Vec` rather than `mem::replace`-ing out the scratch:
            // after a prior multi-block batch the scratch retains a
            // `COMPRESS_SCRATCH_CAPACITY` allocation, and handing that out as an
            // "empty" output would inflate queue memory/backpressure on the
            // ordinal-preserving fast path.
            0 => Vec::new(),
            1 => children.pop().expect("len == 1").data,
            _ => {
                for child in children {
                    self.output_scratch.extend_from_slice(&child.data);
                    self.compressor.recycle_buffer(child.data);
                }
                std::mem::replace(
                    &mut self.output_scratch,
                    Vec::with_capacity(COMPRESS_SCRATCH_CAPACITY),
                )
            }
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
        let bytes = self.assemble_output_bytes();

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

    /// Reference assembly: the straightforward concatenation the move/recycle
    /// fast path must reproduce byte-for-byte.
    fn reference_concat(level: u32, input: &[u8]) -> Vec<u8> {
        let mut c = InlineBgzfCompressor::new(level);
        c.write_all(input).unwrap();
        c.flush().unwrap();
        let mut out = Vec::new();
        for block in c.take_blocks() {
            out.extend_from_slice(&block.data);
        }
        out
    }

    /// Drive `assemble_output_bytes` on a fresh step over `input`.
    fn assemble(level: u32, input: &[u8]) -> Vec<u8> {
        let mut s = BgzfCompress::new(level, 1024);
        s.compressor.write_all(input).unwrap();
        s.compressor.flush().unwrap();
        s.assemble_output_bytes()
    }

    #[test]
    fn single_block_fast_path_is_byte_identical() {
        // Small input → one physical block → move fast path.
        let input = b"the quick brown fox jumps over the lazy dog";
        let mut s = BgzfCompress::new(1, 1024);
        s.compressor.write_all(input).unwrap();
        s.compressor.flush().unwrap();
        assert_eq!(s.compressor.take_blocks().len(), 1, "input should fit one block");

        let assembled = assemble(1, input);
        assert_eq!(assembled, reference_concat(1, input));
        assert_eq!(&assembled[..2], &[0x1f, 0x8b], "valid BGZF magic");
    }

    #[test]
    fn multi_block_concat_path_is_byte_identical() {
        // Input larger than one 64 KiB block forces the multi-block concat +
        // recycle path; output must still match the plain concatenation.
        let mut input = Vec::with_capacity(200 * 1024);
        for i in 0..(200 * 1024u32) {
            input.push((i % 251) as u8);
        }
        let mut s = BgzfCompress::new(1, 1024);
        s.compressor.write_all(&input).unwrap();
        s.compressor.flush().unwrap();
        assert!(s.compressor.take_blocks().len() > 1, "input should span >1 block");

        let assembled = assemble(1, &input);
        assert_eq!(assembled, reference_concat(1, &input));
    }

    #[test]
    fn recycle_buffer_repopulates_pool_for_reuse() {
        // After the multi-block path recycles buffers, a subsequent compression
        // reuses them, and output remains byte-identical to a fresh compressor.
        let big: Vec<u8> = (0..(200 * 1024u32)).map(|i| (i % 251) as u8).collect();
        let mut s = BgzfCompress::new(1, 1024);

        s.compressor.write_all(&big).unwrap();
        s.compressor.flush().unwrap();
        let first = s.assemble_output_bytes();
        assert_eq!(first, reference_concat(1, &big));

        // Second block through the same (now pool-primed) compressor.
        s.compressor.write_all(&big).unwrap();
        s.compressor.flush().unwrap();
        let second = s.assemble_output_bytes();
        assert_eq!(second, reference_concat(1, &big), "recycled buffers must not corrupt output");
    }

    #[test]
    fn empty_batch_returns_unallocated_vec() {
        // An all-clean rejects batch produces no physical BGZF blocks. The
        // assembled output must be a truly empty `Vec`, not the per-worker
        // scratch — otherwise an "empty" output carries a retained
        // `COMPRESS_SCRATCH_CAPACITY` allocation and inflates queue memory.
        let mut s = BgzfCompress::new(1, 1024);

        // Fresh compressor, nothing written → no blocks.
        let fresh = s.assemble_output_bytes();
        assert!(fresh.is_empty());
        assert_eq!(fresh.capacity(), 0, "empty output must not carry a scratch allocation");

        // After a real multi-block batch (which routes through the scratch and
        // replaces it with a `COMPRESS_SCRATCH_CAPACITY` Vec), a following empty
        // batch must still return an unallocated Vec.
        let big: Vec<u8> = (0..(200 * 1024u32)).map(|i| (i % 251) as u8).collect();
        s.compressor.write_all(&big).unwrap();
        s.compressor.flush().unwrap();
        let multi = s.assemble_output_bytes();
        assert!(!multi.is_empty(), "multi-block batch should produce output");

        let after = s.assemble_output_bytes(); // no new input → empty batch
        assert!(after.is_empty());
        assert_eq!(
            after.capacity(),
            0,
            "empty batch after a multi-block batch must not retain scratch"
        );
    }
}
