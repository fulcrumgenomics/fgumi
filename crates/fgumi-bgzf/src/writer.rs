//! BGZF compression utilities for BAM output.
//!
//! This module provides inline BGZF compression for use with the unified pipeline.
//! Each worker thread compresses data inline using `InlineBgzfCompressor`, producing
//! `CompressedBlock` instances that can be written directly to output.
//!
//! Uses libdeflate (via the `bgzf` crate) for high-performance BGZF compression.

use bgzf::{CompressionLevel, Compressor as BgzfCompressor};
use std::io;

// ============================================================================
// Constants
// ============================================================================

/// Maximum uncompressed size for a BGZF block (64KB - header/footer overhead).
pub const BGZF_MAX_BLOCK_SIZE: usize = bgzf::BGZF_BLOCK_SIZE;

// ============================================================================
// Block types
// ============================================================================

/// A compressed BGZF block ready for writing.
#[derive(Debug, Clone)]
pub struct CompressedBlock {
    /// Serial number for ordering.
    pub serial: u64,
    /// Complete BGZF block (header + compressed data + footer).
    pub data: Vec<u8>,
}

// ============================================================================
// Inline BGZF Compressor (for unified pipeline)
// ============================================================================

/// Per-worker state for inline BGZF compression.
///
/// This struct is designed to be used within worker threads of the unified
/// pipeline, where each worker compresses data inline rather than sending
/// to a separate compression thread pool.
///
/// # Usage
///
/// ```ignore
/// let mut compressor = InlineBgzfCompressor::new(6);
///
/// // Write encoded BAM data
/// compressor.write_all(&encoded_record)?;
///
/// // When done with a batch, flush and take the compressed blocks
/// compressor.flush()?;
/// let blocks = compressor.take_blocks();
/// ```
pub struct InlineBgzfCompressor {
    /// Buffer accumulating uncompressed data (up to 64KB).
    buffer: Vec<u8>,
    /// bgzf crate compressor (reused for efficiency).
    compressor: BgzfCompressor,
    /// Completed compressed blocks ready to return.
    completed_blocks: Vec<CompressedBlock>,
    /// Pool of reusable compression buffers.
    buffer_pool: Vec<Vec<u8>>,
    /// Next serial number for block ordering.
    next_serial: u64,
}

impl InlineBgzfCompressor {
    /// Create a new inline compressor with the specified compression level.
    ///
    /// # Arguments
    ///
    /// * `compression_level` - Compression level (0-12, higher = smaller but slower).
    ///   Level 0 writes uncompressed (stored) BGZF blocks, level 1 is fastest DEFLATE,
    ///   level 6 is a good balance.
    ///
    /// # Panics
    ///
    /// Panics if `compression_level` is outside `0..=12`.
    #[must_use]
    pub fn new(compression_level: u32) -> Self {
        let level = u8::try_from(compression_level)
            .ok()
            .filter(|&l| l <= 12)
            .unwrap_or_else(|| panic!("compression level {compression_level} is outside 0..=12"));
        let compression_level_obj = CompressionLevel::new(level)
            .unwrap_or_else(|_| panic!("bgzf rejected compression level {level}"));
        let compressor = BgzfCompressor::new(compression_level_obj);
        Self {
            buffer: Vec::with_capacity(BGZF_MAX_BLOCK_SIZE),
            compressor,
            completed_blocks: Vec::new(),
            buffer_pool: Vec::new(),
            next_serial: 0,
        }
    }

    /// Get mutable access to the internal buffer for direct writes.
    ///
    /// This enables zero-copy serialization by allowing callers to write
    /// directly into the compression buffer, avoiding an intermediate copy.
    ///
    /// After writing to the buffer, call `maybe_compress()` to check if
    /// the buffer is full and needs compression.
    #[inline]
    pub fn buffer_mut(&mut self) -> &mut Vec<u8> {
        &mut self.buffer
    }

    /// Get the current buffer length.
    #[inline]
    #[must_use]
    pub fn buffer_len(&self) -> usize {
        self.buffer.len()
    }

    /// Compress if the buffer is full (>= 64KB).
    ///
    /// Call this after writing to `buffer_mut()` to ensure data is
    /// compressed when the buffer reaches the BGZF block size limit.
    ///
    /// Returns `Ok(true)` if compression occurred, `Ok(false)` otherwise.
    ///
    /// # Errors
    ///
    /// Returns an error if BGZF compression fails.
    #[inline]
    pub fn maybe_compress(&mut self) -> io::Result<bool> {
        if self.buffer.len() >= BGZF_MAX_BLOCK_SIZE {
            self.compress_current_buffer()?;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Write data to the compressor, compressing when the buffer fills.
    ///
    /// Data is buffered until it reaches ~64KB, then compressed into
    /// a BGZF block and added to the completed blocks list.
    ///
    /// # Errors
    ///
    /// Returns an error if BGZF compression fails.
    pub fn write_all(&mut self, data: &[u8]) -> io::Result<()> {
        let mut offset = 0;

        while offset < data.len() {
            let remaining_in_buffer = BGZF_MAX_BLOCK_SIZE - self.buffer.len();
            let to_copy = remaining_in_buffer.min(data.len() - offset);

            self.buffer.extend_from_slice(&data[offset..offset + to_copy]);
            offset += to_copy;

            // If buffer is full, compress it
            if self.buffer.len() >= BGZF_MAX_BLOCK_SIZE {
                self.compress_current_buffer()?;
            }
        }

        Ok(())
    }

    /// Flush any remaining buffered data.
    ///
    /// This compresses any data remaining in the buffer, even if it's
    /// less than 64KB. Call this before `take_blocks()` to ensure
    /// all data is compressed.
    ///
    /// # Errors
    ///
    /// Returns an error if BGZF compression fails.
    pub fn flush(&mut self) -> io::Result<()> {
        if !self.buffer.is_empty() {
            self.compress_current_buffer()?;
        }
        Ok(())
    }

    /// Take all completed compressed blocks.
    ///
    /// Returns the blocks and clears the internal list. The blocks
    /// can then be written to the output file.
    pub fn take_blocks(&mut self) -> Vec<CompressedBlock> {
        std::mem::take(&mut self.completed_blocks)
    }

    /// Return a drained block buffer to the internal pool for reuse by a later
    /// [`compress_current_buffer`](Self::compress_current_buffer).
    ///
    /// Consumers that drive the compressor with [`write_all`](Self::write_all) +
    /// [`flush`](Self::flush) + [`take_blocks`](Self::take_blocks) (rather than
    /// [`write_blocks_to`](Self::write_blocks_to), which recycles automatically)
    /// otherwise leave `buffer_pool` empty, so every block compression allocates
    /// a fresh output `Vec`. Such a consumer can hand back any block `Vec` it is
    /// done with via this method to restore the recycling.
    ///
    /// The buffer is cleared before being pooled. The pool is bounded to a small
    /// number of buffers so a bursty consumer cannot grow it without limit; once
    /// full, the handed-back buffer is simply dropped.
    pub fn recycle_buffer(&mut self, mut buffer: Vec<u8>) {
        /// Cap on pooled buffers. A single compressor processes one block at a
        /// time, so a handful of recycled buffers fully covers steady state.
        const MAX_POOLED_BUFFERS: usize = 4;
        if self.buffer_pool.len() < MAX_POOLED_BUFFERS {
            buffer.clear();
            self.buffer_pool.push(buffer);
        }
    }

    /// Write all completed compressed blocks directly to output and recycle buffers.
    ///
    /// This is efficient for single-threaded use as it writes blocks directly
    /// and recycles their buffers to the pool for reuse, avoiding repeated allocations.
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the output fails.
    pub fn write_blocks_to<W: io::Write + ?Sized>(&mut self, output: &mut W) -> io::Result<()> {
        // Drain into a temporary so we can call `recycle_buffer` (which borrows
        // `self` mutably) for each block without holding a borrow on
        // `self.completed_blocks`. On a write error, restore the unwritten block
        // and the remaining tail to `completed_blocks` so the queue is never
        // silently emptied — a caller can retry or surface the partial state.
        let mut remaining = std::mem::take(&mut self.completed_blocks).into_iter();
        while let Some(block) = remaining.next() {
            if let Err(e) = output.write_all(&block.data) {
                self.completed_blocks.push(block);
                self.completed_blocks.extend(remaining);
                return Err(e);
            }
            // Route the drained buffer through the capped recycle path so the
            // pool stays bounded by MAX_POOLED_BUFFERS, just like the
            // steady-state recycle path.
            self.recycle_buffer(block.data);
        }
        Ok(())
    }

    /// Compress the current buffer and add to completed blocks.
    fn compress_current_buffer(&mut self) -> io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }

        // Get buffer from pool or allocate new
        let mut compressed_data = self.buffer_pool.pop().unwrap_or_default();
        compressed_data.clear();

        // Compress using bgzf crate's Compressor
        self.compressor
            .compress(&self.buffer, &mut compressed_data)
            .map_err(|e| io::Error::other(format!("BGZF compression failed: {e}")))?;

        let serial = self.next_serial;
        self.next_serial += 1;
        self.completed_blocks.push(CompressedBlock { serial, data: compressed_data });

        // Reset buffer for next block
        self.buffer.clear();

        Ok(())
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reader::{
        BGZF_EOF, BGZF_FOOTER_SIZE as READER_FOOTER_SIZE, BGZF_HEADER_SIZE as READER_HEADER_SIZE,
    };

    #[test]
    fn test_bgzf_constants() {
        assert_eq!(BGZF_MAX_BLOCK_SIZE, 65280);
        assert_eq!(READER_HEADER_SIZE, 18);
        assert_eq!(READER_FOOTER_SIZE, 8);
        assert_eq!(BGZF_EOF.len(), 28);
    }

    #[test]
    fn test_compress_small() {
        let mut compressor = InlineBgzfCompressor::new(6);

        compressor.write_all(b"Hello, BGZF!").expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");

        let blocks = compressor.take_blocks();
        assert_eq!(blocks.len(), 1);

        let compressed = &blocks[0];
        assert!(!compressed.data.is_empty());
        // Check BGZF magic
        assert_eq!(&compressed.data[0..2], &[0x1f, 0x8b]);
        // Check BC subfield ID
        assert_eq!(&compressed.data[12..14], b"BC");
    }

    #[test]
    fn test_compress_max_size() {
        let mut compressor = InlineBgzfCompressor::new(6);

        // Create a max-size block with compressible data
        let data = vec![b'A'; BGZF_MAX_BLOCK_SIZE];
        compressor.write_all(&data).expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");

        let blocks = compressor.take_blocks();
        assert_eq!(blocks.len(), 1);

        let compressed = &blocks[0];
        // Max compressed BGZF block is 65536 bytes
        assert!(compressed.data.len() <= 65536);
    }

    #[test]
    fn test_compress_multiple_blocks() {
        let mut compressor = InlineBgzfCompressor::new(6);

        // Write more than one block's worth of data
        let data = vec![b'X'; BGZF_MAX_BLOCK_SIZE + 100];
        compressor.write_all(&data).expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");

        let blocks = compressor.take_blocks();
        // Should produce 2 blocks: one full, one with remaining 100 bytes
        assert_eq!(blocks.len(), 2);
    }

    #[test]
    fn test_write_blocks_to_single() {
        let mut compressor = InlineBgzfCompressor::new(6);

        compressor.write_all(b"Test data").expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");

        let mut output = Vec::new();
        compressor.write_blocks_to(&mut output).expect("writing blocks to output should succeed");

        // Should have written something
        assert!(!output.is_empty());
        // Check BGZF magic
        assert_eq!(&output[0..2], &[0x1f, 0x8b]);
        // Completed blocks should be cleared
        assert!(compressor.take_blocks().is_empty());
    }

    #[test]
    fn test_write_blocks_to_multiple() {
        let mut compressor = InlineBgzfCompressor::new(6);

        // Write more than one block's worth
        let data = vec![b'Y'; BGZF_MAX_BLOCK_SIZE + 50];
        compressor.write_all(&data).expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");

        let mut output = Vec::new();
        compressor.write_blocks_to(&mut output).expect("writing blocks to output should succeed");

        // Should have written two blocks
        assert!(!output.is_empty());
        // Completed blocks should be cleared
        assert!(compressor.take_blocks().is_empty());

        // Count BGZF block headers (0x1f 0x8b magic)
        let block_count = output.windows(2).filter(|w| w == &[0x1f, 0x8b]).count();
        assert_eq!(block_count, 2);
    }

    #[test]
    fn test_write_blocks_to_respects_pool_cap() {
        // Draining many blocks through write_blocks_to must not grow buffer_pool
        // beyond the cap enforced by recycle_buffer (MAX_POOLED_BUFFERS).
        let mut compressor = InlineBgzfCompressor::new(6);

        // Produce many full blocks so there are far more drained buffers than
        // the pool cap. Each full block's worth of data yields one block.
        let data = vec![b'Z'; BGZF_MAX_BLOCK_SIZE * 20];
        compressor.write_all(&data).expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");

        let mut output = Vec::new();
        compressor.write_blocks_to(&mut output).expect("writing blocks to output should succeed");

        // With 20 drained blocks against a cap of 4, the pool must end *exactly*
        // at the cap. A `<= 4` check alone would also pass if `write_blocks_to`
        // stopped recycling and just dropped the drained buffers (the pool would
        // stay at 0), so the equality assertion pins both postconditions at once:
        // bounded growth (not > 4) AND actual reuse (not < 4 / lost recycling).
        assert_eq!(
            compressor.buffer_pool.len(),
            4,
            "buffer_pool should be filled to the cap by recycling, got {}",
            compressor.buffer_pool.len()
        );
    }

    #[test]
    fn test_write_blocks_to_preserves_unwritten_blocks_on_error() {
        // A writer that accepts the first `write` and then fails. Each block is
        // emitted with a single `write_all` (the writer returns the full length
        // each time), so the failure lands partway through the block queue.
        struct FailAfterFirst {
            writes: usize,
        }
        impl io::Write for FailAfterFirst {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                self.writes += 1;
                if self.writes > 1 {
                    return Err(io::Error::other("simulated write failure"));
                }
                Ok(buf.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }

        let mut compressor = InlineBgzfCompressor::new(6);
        // Produce several full blocks so the failure has unwritten blocks to drop.
        let data = vec![b'Q'; BGZF_MAX_BLOCK_SIZE * 3];
        compressor.write_all(&data).expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");
        let serials: Vec<u64> =
            compressor.completed_blocks.iter().map(|block| block.serial).collect();
        assert!(serials.len() >= 3, "expected multiple blocks, got {}", serials.len());

        let mut output = FailAfterFirst { writes: 0 };
        let err =
            compressor.write_blocks_to(&mut output).expect_err("write must surface the error");
        assert_eq!(err.kind(), io::ErrorKind::Other);

        // Only the first block (successfully written) was consumed; the block
        // that failed plus every block after it must survive in the queue, in
        // their original order, so a retry/recovery path can re-emit them.
        let remaining: Vec<u64> =
            compressor.completed_blocks.iter().map(|block| block.serial).collect();
        assert_eq!(remaining, serials[1..].to_vec());
    }

    #[test]
    fn test_write_blocks_to_equivalence() {
        // Test that write_blocks_to produces same output as take_blocks
        let test_data = b"Equivalence test data for blocks";

        // First, use take_blocks
        let mut compressor1 = InlineBgzfCompressor::new(6);
        compressor1.write_all(test_data).expect("writing test data should succeed");
        compressor1.flush().expect("flushing compressor should succeed");
        let blocks = compressor1.take_blocks();
        let mut output1 = Vec::new();
        for block in &blocks {
            output1.extend_from_slice(&block.data);
        }

        // Then, use write_blocks_to
        let mut compressor2 = InlineBgzfCompressor::new(6);
        compressor2.write_all(test_data).expect("writing test data should succeed");
        compressor2.flush().expect("flushing compressor should succeed");
        let mut output2 = Vec::new();
        compressor2.write_blocks_to(&mut output2).expect("writing blocks to output should succeed");

        // Both should produce identical output
        assert_eq!(output1, output2);
    }

    #[test]
    fn test_recycle_buffer_reuses_pool_without_corrupting_output() {
        let data = b"recycle pool roundtrip data";

        // Reference output from a fresh compressor.
        let mut reference = InlineBgzfCompressor::new(6);
        reference.write_all(data).expect("write");
        reference.flush().expect("flush");
        let expected: Vec<u8> = reference.take_blocks().into_iter().flat_map(|b| b.data).collect();

        // Compressor that hands its drained block buffer back to the pool, then
        // compresses again — the recycled buffer must not corrupt the output.
        let mut compressor = InlineBgzfCompressor::new(6);
        compressor.write_all(data).expect("write");
        compressor.flush().expect("flush");
        for block in compressor.take_blocks() {
            compressor.recycle_buffer(block.data);
        }
        compressor.write_all(data).expect("write");
        compressor.flush().expect("flush");
        let second: Vec<u8> = compressor.take_blocks().into_iter().flat_map(|b| b.data).collect();
        assert_eq!(second, expected);
    }

    #[test]
    fn test_recycle_buffer_is_bounded() {
        let mut compressor = InlineBgzfCompressor::new(6);
        // Hand back far more buffers than the cap; the pool must stay bounded.
        for _ in 0..100 {
            compressor.recycle_buffer(vec![0u8; 64]);
        }
        assert!(compressor.buffer_pool.len() <= 4, "pool must be capped");
    }

    #[test]
    fn test_write_blocks_to_empty() {
        let mut compressor = InlineBgzfCompressor::new(6);

        // Don't write anything, just flush
        compressor.flush().expect("flushing compressor should succeed");

        let mut output = Vec::new();
        compressor.write_blocks_to(&mut output).expect("writing blocks to output should succeed");

        // Should be empty since no data was written
        assert!(output.is_empty());
    }

    #[test]
    fn test_compression_level_zero_is_uncompressed() {
        // Highly compressible data: 32 KiB of zeros. At any real compression level
        // (>= 1) this collapses to a tiny block. At level 0 libdeflate emits stored
        // deflate blocks, so the BGZF payload is at least as large as the input.
        // If level 0 is silently clamped up to level 1, the size assertion fails.
        let data = vec![0u8; 32 * 1024];

        let mut compressor_l0 = InlineBgzfCompressor::new(0);
        compressor_l0.write_all(&data).expect("level 0 write should succeed");
        compressor_l0.flush().expect("level 0 flush should succeed");
        let blocks_l0 = compressor_l0.take_blocks();
        assert_eq!(blocks_l0.len(), 1, "expected a single BGZF block");
        let size_l0 = blocks_l0[0].data.len();

        let mut compressor_l1 = InlineBgzfCompressor::new(1);
        compressor_l1.write_all(&data).expect("level 1 write should succeed");
        compressor_l1.flush().expect("level 1 flush should succeed");
        let blocks_l1 = compressor_l1.take_blocks();
        let size_l1 = blocks_l1[0].data.len();

        assert!(
            size_l1 < 1024,
            "sanity: level 1 must compress 32 KiB of zeros to < 1 KiB (got {size_l1} bytes)"
        );
        assert!(
            size_l0 >= data.len(),
            "level 0 must be uncompressed (>= {} bytes), got {size_l0} bytes",
            data.len()
        );

        let mut decompressor = libdeflater::Decompressor::new();
        let mut out = Vec::new();
        crate::reader::decompress_block_slice_into(&blocks_l0[0].data, &mut decompressor, &mut out)
            .expect("decompression of level-0 block should succeed");
        assert_eq!(out, data, "level-0 BGZF block must round-trip to the input");
    }

    #[test]
    fn test_serial_increments() {
        let mut compressor = InlineBgzfCompressor::new(6);

        // Write enough data to produce multiple blocks
        let data = vec![b'X'; BGZF_MAX_BLOCK_SIZE + 100];
        compressor.write_all(&data).expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");

        let blocks = compressor.take_blocks();
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].serial, 0);
        assert_eq!(blocks[1].serial, 1);

        // Write another batch — serial should continue from where it left off
        compressor.write_all(b"more data").expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");

        let blocks2 = compressor.take_blocks();
        assert_eq!(blocks2.len(), 1);
        assert_eq!(blocks2[0].serial, 2);
    }

    #[test]
    #[should_panic(expected = "outside 0..=12")]
    fn test_compress_level_out_of_range_panics() {
        let _ = InlineBgzfCompressor::new(13);
    }
}
