//! BGZF compression pipeline stage.
//!
//! This stage consumes a flat `Vec<u8>` of raw BAM record bytes (without `block_size`
//! prefixes) produced by [`FilterStage`] and compresses them into one or more BGZF blocks
//! using [`InlineBgzfCompressor`].
//!
//! Each call to `process` produces a [`CompressedBatch`] containing all BGZF blocks
//! generated for the input batch. The blocks are ordered by their internal serial counter
//! (see [`InlineBgzfCompressor`]) and may be written directly to the output file once
//! the `ReorderBuffer` has reassembled pipeline batches in sequence.
//!
//! ## Thread safety
//!
//! `CompressStage` is shared across worker threads (`Send + Sync`). Each worker thread
//! maintains its own [`InlineBgzfCompressor`] in thread-local storage so that compression
//! state — including buffer allocations — is not shared or synchronised.

use std::cell::RefCell;
use std::io;

use anyhow::Result;
use fgumi_bgzf::writer::{CompressedBlock, InlineBgzfCompressor};

use crate::stage::{PipelineStage, with_thread_local};

// ============================================================================
// Output type
// ============================================================================

/// One or more BGZF blocks produced by compressing a single pipeline batch.
///
/// Blocks are in the order they were produced and should be written sequentially
/// to the output file.
pub struct CompressedBatch {
    /// BGZF-compressed blocks ready for writing to the output BAM file.
    pub blocks: Vec<CompressedBlock>,
}

// ============================================================================
// Thread-local compressor storage
// ============================================================================

thread_local! {
    /// Per-worker [`InlineBgzfCompressor`].  Initialised on first use.
    static THREAD_COMPRESSOR: RefCell<Option<InlineBgzfCompressor>> = const { RefCell::new(None) };
}

// ============================================================================
// CompressStage
// ============================================================================

/// Pipeline stage that BGZF-compresses filtered raw BAM bytes into [`CompressedBatch`]es.
///
/// Each worker thread owns an [`InlineBgzfCompressor`] stored in thread-local storage.
/// The compressor accumulates data until its internal 64 KB BGZF block fills, then
/// flushes the remaining bytes when `process` completes.
pub struct CompressStage {
    /// BGZF compression level (1 = fastest, 12 = smallest).
    compression_level: u32,
}

impl CompressStage {
    /// Create a new `CompressStage`.
    ///
    /// # Arguments
    ///
    /// * `compression_level` - BGZF compression level in the range `[1, 12]`.
    ///   Values outside this range are clamped by [`InlineBgzfCompressor`].
    ///   A value of 1 is fastest; 6 is a good balance of speed and size.
    #[must_use]
    pub fn new(compression_level: u32) -> Self {
        Self { compression_level }
    }
}

impl PipelineStage for CompressStage {
    /// Input: flat raw BAM bytes (no `block_size` prefixes) from [`FilterStage`].
    type Input = Vec<u8>;
    /// Output: one or more BGZF-compressed blocks.
    type Output = CompressedBatch;

    /// Compress `input` bytes into BGZF blocks using the calling thread's compressor.
    ///
    /// A new [`InlineBgzfCompressor`] is created on the first call within each worker
    /// thread and reused for subsequent calls. The compressor is flushed at the end of
    /// each call so that all bytes for this batch are contained in the returned blocks.
    ///
    /// # Errors
    ///
    /// Returns an error if BGZF compression fails (e.g., libdeflate returns an error).
    fn process(&self, input: Self::Input) -> Result<Self::Output> {
        let compression_level = self.compression_level;
        with_thread_local(
            &THREAD_COMPRESSOR,
            || InlineBgzfCompressor::new(compression_level),
            |compressor| {
                compressor
                    .write_all(&input)
                    .map_err(|e| anyhow::anyhow!("BGZF write_all failed: {e}"))?;
                compressor
                    .flush()
                    .map_err(|e: io::Error| anyhow::anyhow!("BGZF flush failed: {e}"))?;

                let blocks = compressor.take_blocks();
                Ok(CompressedBatch { blocks })
            },
        )
    }

    /// Estimate memory usage as the sum of all compressed block sizes.
    fn output_memory_estimate(&self, output: &Self::Output) -> usize {
        output.blocks.iter().map(|b| b.data.len()).sum()
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Verify that compressing a small non-empty input produces at least one BGZF block
    /// with the correct magic bytes.
    #[test]
    fn test_compress_stage_small_input_produces_bgzf_block() {
        let stage = CompressStage::new(6);
        let input: Vec<u8> = b"Hello, BGZF pipeline!".to_vec();

        let output = stage.process(input).expect("compression should succeed");
        assert!(!output.blocks.is_empty(), "should produce at least one block");

        let first_block = &output.blocks[0];
        assert!(!first_block.data.is_empty(), "block data should be non-empty");
        // BGZF magic: gzip magic 0x1f 0x8b
        assert_eq!(&first_block.data[0..2], &[0x1f, 0x8b], "should start with gzip magic");
    }

    /// Verify that compressing an empty input produces no blocks.
    #[test]
    fn test_compress_stage_empty_input_no_blocks() {
        let stage = CompressStage::new(6);
        let input: Vec<u8> = Vec::new();

        let output = stage.process(input).expect("empty input should succeed");
        assert!(output.blocks.is_empty(), "empty input should produce no blocks");
    }

    /// Verify that `output_memory_estimate` returns the sum of block data lengths.
    #[test]
    fn test_compress_stage_memory_estimate() {
        let stage = CompressStage::new(6);
        let input: Vec<u8> = b"test data for memory estimate".to_vec();

        let output = stage.process(input).expect("should succeed");
        let expected: usize = output.blocks.iter().map(|b| b.data.len()).sum();
        assert_eq!(stage.output_memory_estimate(&output), expected);
    }

    /// Verify that compressing more than 64 KB of data produces multiple blocks.
    #[test]
    fn test_compress_stage_large_input_multiple_blocks() {
        let stage = CompressStage::new(1); // fastest compression for speed
        // 65280 is the BGZF_MAX_BLOCK_SIZE; write slightly more to force 2 blocks.
        let input: Vec<u8> = vec![b'X'; 65_380];

        let output = stage.process(input).expect("large input should succeed");
        assert!(output.blocks.len() >= 2, "large input should produce multiple blocks");
    }
}
