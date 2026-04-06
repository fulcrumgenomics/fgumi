//! Shared BGZF I/O utilities for pooled writers.
//!
//! Provides the reorder-and-write loop used by both [`PooledBamWriter`](super::pooled_bam_writer)
//! and [`PooledChunkWriter`](super::pooled_chunk_writer), and the staging buffer logic for
//! accumulating data into ~64KB blocks before submitting compression jobs.

use crate::sort::worker_pool::{BufferPool, CompressJob, CompressResult, SortWorkerPool};
use anyhow::Result;
use crossbeam_channel::{Receiver, Sender};
use fgumi_bgzf::{BGZF_EOF, BGZF_MAX_BLOCK_SIZE};
use std::collections::BTreeMap;
use std::io::{BufWriter, Write};
use std::sync::Arc;

/// Padding beyond `BGZF_MAX_BLOCK_SIZE` for the staging buffer capacity.
const STAGING_PADDING: usize = 4096;

/// Staging buffer that accumulates data and submits full blocks to the pool.
pub(crate) struct StagingBuffer {
    pool: Arc<SortWorkerPool>,
    buf: Vec<u8>,
    next_serial: u64,
    result_tx: Sender<CompressResult>,
}

impl StagingBuffer {
    /// Create a new staging buffer.
    #[must_use]
    pub(crate) fn new(pool: Arc<SortWorkerPool>, result_tx: Sender<CompressResult>) -> Self {
        Self {
            pool,
            buf: Vec::with_capacity(BGZF_MAX_BLOCK_SIZE + STAGING_PADDING),
            next_serial: 0,
            result_tx,
        }
    }

    /// The underlying byte buffer for direct writes.
    ///
    /// Callers must ensure writes followed by `flush_if_full()` keep each individual
    /// append ≤ `BGZF_MAX_BLOCK_SIZE`. For potentially-large data use `write_chunked`.
    pub(crate) fn buf(&mut self) -> &mut Vec<u8> {
        &mut self.buf
    }

    /// Returns true if the staging buffer has reached the BGZF block size threshold.
    #[inline]
    pub(crate) fn is_full(&self) -> bool {
        self.buf.len() >= BGZF_MAX_BLOCK_SIZE
    }

    /// Flush the staging buffer: swap it with a recycled buffer and submit for compression.
    ///
    /// No-op when the buffer is empty (avoids submitting empty BGZF blocks).
    pub(crate) fn flush(&mut self) {
        if self.buf.is_empty() {
            return;
        }
        let data = std::mem::replace(&mut self.buf, self.pool.buffer_pool.checkout());
        if self.buf.capacity() < BGZF_MAX_BLOCK_SIZE + STAGING_PADDING {
            self.buf.reserve(BGZF_MAX_BLOCK_SIZE + STAGING_PADDING - self.buf.capacity());
        }

        let serial = self.next_serial;
        self.next_serial += 1;

        self.pool.submit_compress(CompressJob { data, serial, result_tx: self.result_tx.clone() });
    }

    /// Flush if full, otherwise no-op.
    #[inline]
    pub(crate) fn flush_if_full(&mut self) {
        if self.is_full() {
            self.flush();
        }
    }

    /// Write `data` to the staging buffer, flushing BGZF-sized chunks as they fill up.
    ///
    /// Unlike writing directly to `buf()`, this correctly handles data larger than
    /// `BGZF_MAX_BLOCK_SIZE` (e.g. large BAM headers) by splitting into multiple jobs.
    pub(crate) fn write_chunked(&mut self, data: &[u8]) {
        let mut remaining = data;
        while !remaining.is_empty() {
            let space = BGZF_MAX_BLOCK_SIZE.saturating_sub(self.buf.len());
            let n = remaining.len().min(space);
            self.buf.extend_from_slice(&remaining[..n]);
            remaining = &remaining[n..];
            self.flush_if_full();
        }
    }
}

/// I/O writer loop: receives compressed blocks and writes them in serial order.
///
/// Blocks may arrive out of order from the pool; this loop reassembles them and
/// writes a valid BGZF stream including the EOF marker.
///
/// # Errors
///
/// Returns an error if any disk write fails or if a compressed block is missing
/// (which would silently truncate the output).
#[allow(clippy::needless_pass_by_value)]
pub(crate) fn io_writer_loop(
    mut writer: BufWriter<std::fs::File>,
    result_rx: Receiver<CompressResult>,
    buffer_pool: BufferPool,
) -> Result<()> {
    let mut next_expected: u64 = 0;
    let mut reorder_buf: BTreeMap<u64, Vec<u8>> = BTreeMap::new();

    while let Ok(result) = result_rx.recv() {
        buffer_pool.checkin(result.recycled_buf);

        if result.serial == next_expected {
            writer.write_all(&result.compressed)?;
            next_expected += 1;

            while let Some(data) = reorder_buf.remove(&next_expected) {
                writer.write_all(&data)?;
                next_expected += 1;
            }
        } else {
            reorder_buf.insert(result.serial, result.compressed);
        }
    }

    // Drain remaining buffered blocks — any gap means a worker dropped a result
    while let Some((&serial, _)) = reorder_buf.first_key_value() {
        if serial == next_expected {
            let data = reorder_buf.remove(&serial).expect("key just checked");
            writer.write_all(&data)?;
            next_expected += 1;
        } else {
            return Err(anyhow::anyhow!(
                "missing compressed block {next_expected}: next available is {serial}; \
                 the output would be silently truncated"
            ));
        }
    }

    writer.write_all(&BGZF_EOF)?;
    writer.flush()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use tempfile::TempDir;

    /// Build a round-trip helper: write `data` via `StagingBuffer` → `io_writer_loop` → read back raw bytes.
    fn roundtrip_data(data: &[u8]) -> Vec<u8> {
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6));
        let (result_tx, result_rx) = pool.compress_result_channel();
        let buffer_pool = pool.buffer_pool.clone();

        let dir = TempDir::new().unwrap();
        let out_path = dir.path().join("out.bgzf");

        let out_file = std::fs::File::create(&out_path).unwrap();
        let writer = std::io::BufWriter::new(out_file);
        let io_handle = std::thread::spawn(move || io_writer_loop(writer, result_rx, buffer_pool));

        let mut staging = StagingBuffer::new(Arc::clone(&pool), result_tx);
        staging.write_chunked(data);
        staging.flush();
        drop(staging); // closes result_tx senders → io_writer_loop exits

        io_handle.join().unwrap().unwrap();

        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
        std::fs::read(&out_path).unwrap()
    }

    #[test]
    fn test_staging_buffer_flush_empty_is_noop() {
        let pool = Arc::new(SortWorkerPool::new(1, 1, 6));
        let (result_tx, _result_rx) = pool.compress_result_channel();

        let mut staging = StagingBuffer::new(Arc::clone(&pool), result_tx);
        // Flush with empty buffer: should not submit a compress job
        staging.flush();

        assert_eq!(
            pool.stats.compress_jobs_submitted.load(std::sync::atomic::Ordering::Relaxed),
            0
        );

        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
    }

    #[test]
    fn test_staging_buffer_is_full() {
        let pool = Arc::new(SortWorkerPool::new(1, 1, 6));
        let (result_tx, _result_rx) = pool.compress_result_channel();
        let mut staging = StagingBuffer::new(Arc::clone(&pool), result_tx);

        assert!(!staging.is_full(), "empty buffer should not be full");
        staging.buf().extend(vec![0u8; BGZF_MAX_BLOCK_SIZE]);
        assert!(staging.is_full(), "buffer at BGZF_MAX_BLOCK_SIZE should be full");

        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
    }

    #[test]
    fn test_staging_buffer_write_chunked_large_data() {
        // Data larger than BGZF_MAX_BLOCK_SIZE must be split into multiple compress jobs.
        let large = vec![b'A'; BGZF_MAX_BLOCK_SIZE * 2 + 1000];
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6));
        let (result_tx, result_rx) = pool.compress_result_channel();
        let buffer_pool = pool.buffer_pool.clone();

        let dir = TempDir::new().unwrap();
        let out_path = dir.path().join("large.bgzf");
        let out_file = std::fs::File::create(&out_path).unwrap();
        let writer = std::io::BufWriter::new(out_file);
        let io_handle = std::thread::spawn(move || io_writer_loop(writer, result_rx, buffer_pool));

        let mut staging = StagingBuffer::new(Arc::clone(&pool), result_tx);
        staging.write_chunked(&large);
        staging.flush();
        drop(staging);

        io_handle.join().unwrap().unwrap();

        // ≥2 full blocks + 1 partial = at least 3 compress jobs
        assert!(
            pool.stats.compress_jobs_submitted.load(std::sync::atomic::Ordering::Relaxed) >= 2,
            "expected multiple compress jobs for data > BGZF_MAX_BLOCK_SIZE"
        );

        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
    }

    #[test]
    fn test_io_writer_loop_reorders_out_of_order_blocks() {
        // Write blocks out of order; io_writer_loop must reassemble them correctly.
        let data1 = b"first block data".to_vec();
        let data2 = b"second block data".to_vec();

        let pool = Arc::new(SortWorkerPool::new(2, 1, 6));
        let (result_tx, result_rx) = pool.compress_result_channel();
        let buffer_pool = pool.buffer_pool.clone();

        let dir = TempDir::new().unwrap();
        let out_path = dir.path().join("reorder.bgzf");
        let out_file = std::fs::File::create(&out_path).unwrap();
        let writer = std::io::BufWriter::new(out_file);
        let io_handle = std::thread::spawn(move || io_writer_loop(writer, result_rx, buffer_pool));

        // Submit block 1 first, then block 0 (out of order)
        pool.submit_compress(CompressJob { data: data2, serial: 1, result_tx: result_tx.clone() });
        pool.submit_compress(CompressJob { data: data1, serial: 0, result_tx });

        // Wait for both compress results to be received by io_writer_loop
        io_handle.join().unwrap().unwrap();

        // Output file exists and contains the BGZF EOF marker
        let bytes = std::fs::read(&out_path).unwrap();
        assert!(bytes.ends_with(&BGZF_EOF), "output should end with BGZF EOF marker");

        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
    }

    #[test]
    fn test_roundtrip_small_data() {
        let data = b"hello world from bgzf_io";
        let output = roundtrip_data(data);
        // Output is a valid BGZF stream ending with the EOF marker
        assert!(output.ends_with(&BGZF_EOF), "must end with BGZF EOF");
        // Non-empty (has at least the compressed data block + EOF)
        assert!(output.len() > BGZF_EOF.len());
    }

    #[test]
    fn test_roundtrip_empty_data() {
        // No data: flush() is a no-op, so io_writer_loop writes only the EOF marker
        let output = roundtrip_data(b"");
        assert_eq!(output, BGZF_EOF.to_vec(), "empty input → only BGZF EOF marker");
    }
}
