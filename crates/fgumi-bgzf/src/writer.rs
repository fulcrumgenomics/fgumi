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
const BGZF_MAX_BLOCK_SIZE: usize = bgzf::BGZF_BLOCK_SIZE;

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
}

impl InlineBgzfCompressor {
    /// Create a new inline compressor with the specified compression level.
    ///
    /// # Arguments
    ///
    /// * `compression_level` - Compression level (1-12, higher = smaller but slower).
    ///   Level 1 is fastest compression, level 6 is a good balance.
    ///
    /// # Panics
    ///
    /// Panics if compression level 6 is rejected by the bgzf library (should never happen).
    #[must_use]
    pub fn new(compression_level: u32) -> Self {
        let level = u8::try_from(compression_level.clamp(1, 12))
            .expect("value in [1, 12] always fits in u8");
        let compression_level_obj =
            CompressionLevel::new(level).unwrap_or_else(|_| CompressionLevel::new(6).unwrap());
        let compressor = BgzfCompressor::new(compression_level_obj);
        Self {
            buffer: Vec::with_capacity(BGZF_MAX_BLOCK_SIZE),
            compressor,
            completed_blocks: Vec::new(),
            buffer_pool: Vec::new(),
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

            self.buffer
                .extend_from_slice(&data[offset..offset + to_copy]);
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

    /// Write all completed compressed blocks directly to output and recycle buffers.
    ///
    /// This is efficient for single-threaded use as it writes blocks directly
    /// and recycles their buffers to the pool for reuse, avoiding repeated allocations.
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the output fails.
    pub fn write_blocks_to<W: io::Write + ?Sized>(&mut self, output: &mut W) -> io::Result<()> {
        for block in self.completed_blocks.drain(..) {
            output.write_all(&block.data)?;
            // Recycle the buffer for reuse
            let mut buf = block.data;
            buf.clear();
            self.buffer_pool.push(buf);
        }
        Ok(())
    }

    /// Recycle block buffers back to the pool for reuse.
    ///
    /// Call this after writing blocks to return their buffers for reuse,
    /// reducing allocation overhead in the single-threaded pipeline.
    pub fn recycle_blocks(&mut self, blocks: Vec<CompressedBlock>) {
        for block in blocks {
            let mut buf = block.data;
            buf.clear();
            self.buffer_pool.push(buf);
        }
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

        self.completed_blocks.push(CompressedBlock {
            serial: 0,
            data: compressed_data,
        });

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
    use crate::reader::{BGZF_FOOTER_SIZE as READER_FOOTER_SIZE, BGZF_HEADER_SIZE as READER_HEADER_SIZE, BGZF_EOF};

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

        compressor.write_all(b"Hello, BGZF!").unwrap();
        compressor.flush().unwrap();

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
        compressor.write_all(&data).unwrap();
        compressor.flush().unwrap();

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
        compressor.write_all(&data).unwrap();
        compressor.flush().unwrap();

        let blocks = compressor.take_blocks();
        // Should produce 2 blocks: one full, one with remaining 100 bytes
        assert_eq!(blocks.len(), 2);
    }

    #[test]
    fn test_write_blocks_to_single() {
        let mut compressor = InlineBgzfCompressor::new(6);

        compressor.write_all(b"Test data").unwrap();
        compressor.flush().unwrap();

        let mut output = Vec::new();
        compressor.write_blocks_to(&mut output).unwrap();

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
        compressor.write_all(&data).unwrap();
        compressor.flush().unwrap();

        let mut output = Vec::new();
        compressor.write_blocks_to(&mut output).unwrap();

        // Should have written two blocks
        assert!(!output.is_empty());
        // Completed blocks should be cleared
        assert!(compressor.take_blocks().is_empty());

        // Count BGZF block headers (0x1f 0x8b magic)
        let block_count = output.windows(2).filter(|w| w == &[0x1f, 0x8b]).count();
        assert_eq!(block_count, 2);
    }

    #[test]
    fn test_write_blocks_to_equivalence() {
        // Test that write_blocks_to produces same output as take_blocks
        let test_data = b"Equivalence test data for blocks";

        // First, use take_blocks
        let mut compressor1 = InlineBgzfCompressor::new(6);
        compressor1.write_all(test_data).unwrap();
        compressor1.flush().unwrap();
        let blocks = compressor1.take_blocks();
        let mut output1 = Vec::new();
        for block in &blocks {
            output1.extend_from_slice(&block.data);
        }

        // Then, use write_blocks_to
        let mut compressor2 = InlineBgzfCompressor::new(6);
        compressor2.write_all(test_data).unwrap();
        compressor2.flush().unwrap();
        let mut output2 = Vec::new();
        compressor2.write_blocks_to(&mut output2).unwrap();

        // Both should produce identical output
        assert_eq!(output1, output2);
    }

    #[test]
    fn test_write_blocks_to_empty() {
        let mut compressor = InlineBgzfCompressor::new(6);

        // Don't write anything, just flush
        compressor.flush().unwrap();

        let mut output = Vec::new();
        compressor.write_blocks_to(&mut output).unwrap();

        // Should be empty since no data was written
        assert!(output.is_empty());
    }
}
