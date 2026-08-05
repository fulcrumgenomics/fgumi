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

/// Cap on the number of compression buffers [`InlineBgzfCompressor`] pools for
/// reuse.
///
/// The pool only has to cover the *compress* side, which takes one buffer at a
/// time: `compress_current_buffer` pops exactly one per block. The queue on the
/// other side can be arbitrarily long — `write_all` appends a block per 64 KiB
/// written — so a drain of many blocks hands back more buffers than the
/// compressor will ever want at once. The surplus is dropped rather than parked
/// for the compressor's lifetime, which is what keeps the pool's memory a
/// function of this constant instead of of how much the caller buffered before
/// draining.
pub const MAX_POOLED_BUFFERS: usize = 4;

/// Largest buffer [`InlineBgzfCompressor::recycle_buffer`] will pool.
///
/// A block buffer is sized by `Compressor::compress` to
/// `header + deflate_compress_bound(input) + footer`. That bound exceeds the
/// input for incompressible data, so a legitimate buffer can run somewhat over
/// [`BGZF_MAX_BLOCK_SIZE`]; double it so no real block buffer is ever refused,
/// while still rejecting a wildly oversized `Vec` a caller hands back.
pub const MAX_POOLED_BUFFER_BYTES: usize = BGZF_MAX_BLOCK_SIZE * 2;

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
    /// block compression.
    ///
    /// Consumers that drive the compressor with [`write_all`](Self::write_all) +
    /// [`flush`](Self::flush) + [`take_blocks`](Self::take_blocks) (rather than
    /// [`write_blocks_to`](Self::write_blocks_to), which recycles automatically)
    /// otherwise leave `buffer_pool` empty, so every block compression allocates
    /// a fresh output `Vec`. Such a consumer can hand back any block `Vec` it is
    /// done with via this method to restore the recycling.
    ///
    /// The buffer is cleared — but not shrunk, since its capacity is the whole
    /// point — before being pooled. A buffer is only pooled if the pool has
    /// room ([`MAX_POOLED_BUFFERS`]) *and* the buffer is a plausible block
    /// buffer ([`MAX_POOLED_BUFFER_BYTES`]); otherwise it is dropped. Bounding
    /// the count alone would let one oversized `Vec` handed in by a caller sit
    /// in the pool for the compressor's lifetime, so the pool's memory is
    /// bounded on both axes.
    ///
    /// Not to be confused with the pipeline's `WorkerCoreState::recycle_buffer`,
    /// which pools the *uncompressed* worker buffers under a different cap.
    /// This one pools the compressor's own compressed-block output buffers.
    pub fn recycle_buffer(&mut self, mut buffer: Vec<u8>) {
        if self.buffer_pool.len() < MAX_POOLED_BUFFERS
            && buffer.capacity() <= MAX_POOLED_BUFFER_BYTES
        {
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
    /// Returns an error if writing to the output fails. The block that failed
    /// to write, and every block after it, is left in the queue, so they are
    /// available to **inspect** rather than being silently discarded.
    ///
    /// They are *not* safe to replay. [`io::Write::write_all`] loops over
    /// [`write`](io::Write::write), advancing past each `Ok(n)`, so it can
    /// commit part of the failing block to `output` before an error surfaces —
    /// and the whole block is re-queued, not the unwritten remainder. Writing
    /// the queue again would repeat those bytes inside a gzip member and
    /// produce a stream no BGZF reader can decode. Treat `output` as
    /// unrecoverable after this returns an error: the queued blocks are for
    /// diagnostics, or for writing somewhere new.
    pub fn write_blocks_to<W: io::Write + ?Sized>(&mut self, output: &mut W) -> io::Result<()> {
        // Nothing queued is the common call (a flush with no pending blocks).
        // Returning here keeps `completed_blocks`'s allocation, which the
        // `mem::take` below would otherwise hand away and force a realloc on
        // the next compression.
        if self.completed_blocks.is_empty() {
            return Ok(());
        }
        // Drain into a temporary so `recycle_buffer` (which borrows `self`
        // mutably) can be called per block without holding a borrow on
        // `self.completed_blocks`.
        let mut remaining = std::mem::take(&mut self.completed_blocks).into_iter();
        while let Some(block) = remaining.next() {
            if let Err(e) = output.write_all(&block.data) {
                self.completed_blocks.push(block);
                self.completed_blocks.extend(remaining);
                return Err(e);
            }
            // Route the drained buffer through the capped recycle path so the
            // pool stays bounded, just like the steady-state recycle path.
            self.recycle_buffer(block.data);
        }
        Ok(())
    }

    /// Compress the current buffer and add to completed blocks.
    fn compress_current_buffer(&mut self) -> io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }

        // Get buffer from pool or allocate new. No clear is needed: the only
        // thing this buffer is handed to is `Compressor::compress`, which
        // resizes it from zero before writing (bgzf 0.4's `resize_uninit`
        // begins with `Vec::clear`). Capacity is what the pool is for, and that
        // survives.
        let mut compressed_data = self.buffer_pool.pop().unwrap_or_default();

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

    // ── Buffer recycling ────────────────────────────────────────────────────

    /// Compress `payload` into one block and hand the compressor back.
    fn compressor_with_one_block(payload: &[u8]) -> (InlineBgzfCompressor, Vec<CompressedBlock>) {
        let mut compressor = InlineBgzfCompressor::new(6);
        compressor.write_all(payload).expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");
        let blocks = compressor.take_blocks();
        (compressor, blocks)
    }

    /// A `take_blocks` consumer keeps ownership of every block buffer, so the
    /// pool is left empty and each subsequent compression allocates afresh.
    /// Handing a drained buffer back must restore the recycling.
    #[test]
    fn test_recycle_buffer_repopulates_the_pool() {
        let (mut compressor, blocks) = compressor_with_one_block(b"payload for one block");
        assert!(
            compressor.buffer_pool.is_empty(),
            "take_blocks should leave the pool empty, that is the gap recycle_buffer closes"
        );

        let buffer = blocks.into_iter().next().expect("one block").data;
        let recycled_capacity = buffer.capacity();
        // The identity of the allocation, not just of the `Vec`. Everything
        // below is about this exact heap block being handed back to the next
        // compression rather than freed and replaced.
        let recycled_ptr = buffer.as_ptr();
        compressor.recycle_buffer(buffer);
        assert_eq!(compressor.buffer_pool.len(), 1);
        // The allocation is what is being reused, so the capacity has to survive
        // pooling -- a `shrink_to_fit` alongside the `clear` would satisfy every
        // other assertion here while deleting the point of the method.
        assert_eq!(
            compressor.buffer_pool[0].capacity(),
            recycled_capacity,
            "pooling must retain the buffer's allocation, not just the buffer"
        );

        // The next compression must consume the pooled buffer rather than
        // allocate. Asserting the pool merely drained would also pass for a
        // pop-and-discard implementation, so pin the produced block to the
        // recycled allocation instead.
        compressor.write_all(b"a second block").expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");
        assert!(compressor.buffer_pool.is_empty(), "the pooled buffer should have been taken");
        let next = compressor.take_blocks();
        assert_eq!(
            next[0].data.as_ptr(),
            recycled_ptr,
            "the new block should be built in the recycled allocation, not a fresh one"
        );
    }

    /// The pool is bounded on capacity as well as count: a caller handing back
    /// an outsized `Vec` must not park that allocation for the compressor's
    /// lifetime. A real block buffer (~64 KiB) is well under the limit, so this
    /// rejects only what it is meant to.
    #[test]
    fn test_recycle_buffer_refuses_an_oversized_buffer() {
        let (mut compressor, blocks) = compressor_with_one_block(b"payload for one block");
        let real = blocks.into_iter().next().expect("one block").data;
        assert!(
            real.capacity() <= MAX_POOLED_BUFFER_BYTES,
            "a real block buffer ({} bytes) must still be poolable",
            real.capacity()
        );

        compressor.recycle_buffer(Vec::with_capacity(MAX_POOLED_BUFFER_BYTES + 1));
        assert!(compressor.buffer_pool.is_empty(), "an oversized buffer should not be pooled");

        compressor.recycle_buffer(real);
        assert_eq!(compressor.buffer_pool.len(), 1, "a real block buffer should be pooled");
    }

    /// `write_blocks_to` is the other path into the pool, and the line this
    /// change rewrote. It must recycle what it drains, and honour the same cap
    /// as the explicit entry point -- a drain of many blocks cannot grow the
    /// pool past it.
    #[test]
    fn test_write_blocks_to_recycles_up_to_the_cap() {
        // More blocks than the cap, so both halves of the behaviour are visible.
        let payload = vec![b'x'; BGZF_MAX_BLOCK_SIZE * (MAX_POOLED_BUFFERS + 2)];
        let mut compressor = InlineBgzfCompressor::new(6);
        compressor.write_all(&payload).expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");
        let queued = compressor.completed_blocks.len();
        assert!(queued > MAX_POOLED_BUFFERS, "need more than {MAX_POOLED_BUFFERS} blocks");
        assert!(compressor.buffer_pool.is_empty(), "pool starts empty");

        let mut sink = Vec::new();
        compressor.write_blocks_to(&mut sink).expect("writing blocks should succeed");

        assert!(compressor.completed_blocks.is_empty(), "every block should have been written");
        assert_eq!(
            compressor.buffer_pool.len(),
            MAX_POOLED_BUFFERS,
            "draining {queued} blocks should fill the pool to the cap and drop the surplus"
        );
    }

    /// A consumer that hands back more buffers than the compressor can use
    /// must not grow the pool without limit; surplus buffers are dropped.
    #[test]
    fn test_recycle_buffer_bounds_the_pool() {
        let (mut compressor, _blocks) = compressor_with_one_block(b"payload for one block");
        for _ in 0..32 {
            compressor.recycle_buffer(vec![0u8; 1024]);
        }
        assert_eq!(
            compressor.buffer_pool.len(),
            MAX_POOLED_BUFFERS,
            "pool should stop growing at the cap"
        );
    }

    /// A sink that accepts `accept` bytes total and then fails, so a
    /// multi-block `write_blocks_to` fails partway through.
    struct ShortSink {
        accept: usize,
    }

    impl io::Write for ShortSink {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if buf.len() > self.accept {
                return Err(io::Error::new(io::ErrorKind::WriteZero, "sink is full"));
            }
            self.accept -= buf.len();
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// On a write failure, the unwritten blocks must stay queued. Draining them
    /// into the void would lose compressed output with no way for the caller to
    /// notice, since the error says nothing about how far the write got.
    #[test]
    fn test_write_blocks_to_keeps_unwritten_blocks_on_error() {
        // Two full blocks, then a sink that has room for the first one only.
        let payload = vec![b'x'; BGZF_MAX_BLOCK_SIZE * 2];
        let mut compressor = InlineBgzfCompressor::new(6);
        compressor.write_all(&payload).expect("writing data should succeed");
        compressor.flush().expect("flushing compressor should succeed");
        let first_block_len = compressor.completed_blocks[0].data.len();
        let queued = compressor.completed_blocks.len();
        assert!(queued >= 2, "expected at least two blocks, got {queued}");

        let mut sink = ShortSink { accept: first_block_len };
        let err = compressor
            .write_blocks_to(&mut sink)
            .expect_err("the sink should reject the second block");
        assert_eq!(err.kind(), io::ErrorKind::WriteZero);
        // Identity, not just the count: a queue of the right length holding the
        // wrong blocks (say, the tail restored but the failing block dropped)
        // loses exactly the output this behaviour exists to preserve. Serials
        // are contiguous from 0, so the retained queue must start at 1 and run
        // to the end.
        let retained: Vec<u64> = compressor.completed_blocks.iter().map(|b| b.serial).collect();
        let expected: Vec<u64> = (1..u64::try_from(queued).expect("block count fits")).collect();
        assert_eq!(retained, expected, "the failing block and its tail should still be queued");
    }
}
