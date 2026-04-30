//! Parallel chunk writer using `SortWorkerPool` for BGZF compression.
//!
//! `PooledChunkWriter` replaces single-threaded `GenericKeyedChunkWriter` during spill,
//! distributing BGZF block compression across the shared worker pool. This targets the
//! #1 bottleneck: spill write is 62-80% of total sort wall time.
//!
//! # Architecture
//!
//! ```text
//! Main thread                   Worker pool            I/O thread
//! ─────────────                 ───────────            ──────────
//! write_record() ──►            (compress)   ──►       write to disk
//!   buffer ~64KB    submit job   parallel      result   sequential
//!   submit job      ──────────► workers        ──────► reorder buf
//!                                              ──────► write ordered
//! ```
//!
//! The reorder buffer is bounded (capacity = `num_workers * 4`) via a `PermitPool`:
//! `StagingBuffer::flush()` acquires a permit before each submit; `io_writer_loop`
//! releases it after each block is handed to the OS write buffer (`write_all` on a
//! `BufWriter`), blocking the main thread when the budget is exhausted rather than
//! accumulating blocks in an unbounded reorder buffer.

use crate::sort::bgzf_io::{StagingBuffer, io_writer_loop};
use crate::sort::keys::RawSortKey;
use crate::sort::worker_pool::{CompressResult, PermitPool, SortWorkerPool};
use anyhow::Result;
use crossbeam_channel::bounded;
use fgumi_bgzf::BGZF_MAX_BLOCK_SIZE;
use std::io::BufWriter;
use std::marker::PhantomData;
use std::path::Path;
use std::sync::Arc;
use std::thread::{self, JoinHandle};

/// A chunk writer that uses `SortWorkerPool` for parallel BGZF compression.
///
/// Records are buffered into ~64KB staging blocks, then submitted to the pool
/// for compression. An I/O thread receives compressed blocks and writes them
/// in serial order using a bounded reorder buffer.
pub struct PooledChunkWriter<K: RawSortKey> {
    /// `None` only after `start_finish()` transfers ownership to `SpillWriteHandle`.
    staging: Option<StagingBuffer>,
    /// Reusable scratch buffer for key serialization (non-embedded keys only).
    key_buf: Vec<u8>,
    io_handle: Option<JoinHandle<Result<()>>>,
    _phantom: PhantomData<K>,
}

impl<K: RawSortKey> PooledChunkWriter<K> {
    /// Create a new pooled chunk writer.
    ///
    /// Opens the output file and spawns an I/O writer thread that receives
    /// compressed blocks and writes them in serial order.
    ///
    /// # Errors
    ///
    /// Returns an error if the output file cannot be created.
    pub fn new(pool: Arc<SortWorkerPool>, path: &Path) -> Result<Self> {
        let file = std::fs::File::create(path)?;
        let writer = BufWriter::with_capacity(256 * 1024, file);

        let reorder_capacity = pool.num_workers() * 4;
        let (result_tx, result_rx) = bounded::<CompressResult>(reorder_capacity);
        let buffer_pool = pool.buffer_pool.clone();
        let permit_pool = Arc::new(PermitPool::new(reorder_capacity));

        let pp = Arc::clone(&permit_pool);
        let io_handle = thread::spawn(move || io_writer_loop(writer, result_rx, buffer_pool, pp));

        Ok(Self {
            staging: Some(StagingBuffer::new(pool, result_tx, permit_pool)),
            key_buf: Vec::new(),
            io_handle: Some(io_handle),
            _phantom: PhantomData,
        })
    }

    /// Write a keyed record to the chunk file.
    ///
    /// Buffers the record and its key into a staging area. When the staging
    /// area reaches ~64KB, it's submitted to the pool for compression.
    ///
    /// # Errors
    ///
    /// Returns an error if key serialization fails.
    ///
    /// # Panics
    ///
    /// Panics if called after [`start_finish`](Self::start_finish) has been called.
    #[allow(clippy::cast_possible_truncation)]
    pub fn write_record(&mut self, key: &K, record: &[u8]) -> Result<()> {
        let staging = self.staging.as_mut().expect("write_record called after start_finish");
        if K::EMBEDDED_IN_RECORD {
            // Fast path: key is part of the record bytes, no extra serialization.
            // Budget: 4-byte length prefix + record bytes.
            let needed = 4 + record.len();
            if staging.buf().len() + needed > BGZF_MAX_BLOCK_SIZE {
                staging.flush()?;
            }
            staging.buf().extend_from_slice(&(record.len() as u32).to_le_bytes());
            if record.len() > BGZF_MAX_BLOCK_SIZE.saturating_sub(4) {
                staging.write_chunked(record)?;
            } else {
                staging.buf().extend_from_slice(record);
                staging.flush_if_full()?;
            }
        } else {
            // Non-embedded key: serialize key into a reusable scratch buffer so we
            // know its exact size before the pre-flush check.
            self.key_buf.clear();
            key.write_to(&mut self.key_buf)?;
            let needed = self.key_buf.len() + 4 + record.len();
            // No size limit check: records larger than one BGZF block are handled
            // by write_chunked(), which splits them across multiple blocks. The
            // reader uses streaming read_exact() that transparently spans blocks.
            if staging.buf().len() + needed > BGZF_MAX_BLOCK_SIZE {
                staging.flush()?;
            }
            staging.buf().extend_from_slice(&self.key_buf);
            staging.buf().extend_from_slice(&(record.len() as u32).to_le_bytes());
            staging.write_chunked(record)?;
        }
        Ok(())
    }

    /// Finish writing: flush remaining data, wait for I/O thread.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing or the I/O thread encountered an error.
    pub fn finish(self) -> Result<()> {
        self.start_finish()?.wait()
    }

    /// Flush remaining data and signal the I/O thread, but don't wait for it.
    ///
    /// Returns a [`SpillWriteHandle`] that can be waited on later. This allows
    /// the caller to overlap the I/O thread's disk writes with other work
    /// (e.g., reading the next batch of records).
    ///
    /// # Errors
    ///
    /// Returns an error if flushing the staging buffer fails.
    pub fn start_finish(mut self) -> Result<SpillWriteHandle> {
        if let Some(mut staging) = self.staging.take() {
            if !staging.buf().is_empty() {
                staging.flush()?;
            }
            drop(staging); // closes result_tx → I/O thread exits after draining
        }
        Ok(SpillWriteHandle::new(self.io_handle.take()))
    }
}

impl<K: RawSortKey> Drop for PooledChunkWriter<K> {
    fn drop(&mut self) {
        if self.io_handle.is_some() {
            // Writer dropped before finish()/start_finish() (e.g. early error return).
            // Drop staging first — this closes result_tx, signaling the I/O thread to
            // drain and exit. Then join the thread to avoid silently detaching it.
            drop(self.staging.take());
            if let Some(handle) = self.io_handle.take() {
                match handle.join() {
                    Ok(Ok(())) => {}
                    Ok(Err(e)) => log::error!("PooledChunkWriter: I/O writer thread error: {e}"),
                    Err(_) => log::error!("PooledChunkWriter: I/O writer thread panicked"),
                }
            }
        }
    }
}

/// Handle for a spill write that is finishing in the background.
///
/// Created by [`PooledChunkWriter::start_finish`]. The I/O thread continues
/// writing compressed blocks to disk. Call [`wait`](SpillWriteHandle::wait)
/// to block until all data is written and the file is closed.
///
/// If dropped without calling `wait`, the `Drop` impl joins the thread and
/// logs any error rather than silently detaching it.
#[must_use = "call wait() to propagate write errors; dropping silently logs them"]
pub struct SpillWriteHandle {
    io_handle: Option<JoinHandle<Result<()>>>,
}

impl SpillWriteHandle {
    /// Create a new handle wrapping an I/O thread join handle.
    pub(crate) fn new(io_handle: Option<JoinHandle<Result<()>>>) -> Self {
        Self { io_handle }
    }

    /// Wait for the background I/O thread to finish writing all blocks.
    ///
    /// # Errors
    ///
    /// Returns an error if the I/O thread panicked or encountered a write error.
    pub fn wait(mut self) -> Result<()> {
        if let Some(handle) = self.io_handle.take() {
            handle.join().map_err(|_| anyhow::anyhow!("I/O writer thread panicked"))??;
        }
        Ok(())
    }
}

impl Drop for SpillWriteHandle {
    fn drop(&mut self) {
        if let Some(handle) = self.io_handle.take() {
            match handle.join() {
                Ok(Ok(())) => {}
                Ok(Err(e)) => log::error!("SpillWriteHandle: I/O writer thread error: {e}"),
                Err(_) => log::error!("SpillWriteHandle: I/O writer thread panicked"),
            }
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sort::inline_buffer::TemplateKey;
    use crate::sort::raw::GenericKeyedChunkReader;
    use tempfile::TempDir;

    /// Create a test `TemplateKey` with distinct values for roundtrip verification.
    #[allow(clippy::cast_possible_truncation)]
    fn make_key(i: u64) -> TemplateKey {
        TemplateKey::new(
            i as i32,   // tid1
            i as i32,   // pos1
            false,      // neg1
            i32::MAX,   // tid2
            i32::MAX,   // pos2
            false,      // neg2
            0,          // cb_hash
            0,          // library
            (0, false), // mi
            i,          // name_hash
            false,      // is_upper
        )
    }

    #[test]
    #[allow(clippy::cast_possible_truncation)]
    fn test_pooled_writer_roundtrip() {
        let dir = TempDir::new().unwrap();
        let chunk_path = dir.path().join("test_chunk.keyed");
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6));

        let records: Vec<(TemplateKey, Vec<u8>)> = (0..100)
            .map(|i| {
                let key = make_key(i);
                let record = vec![(i % 256) as u8; 200 + (i as usize % 50)];
                (key, record)
            })
            .collect();

        {
            let mut writer = PooledChunkWriter::<TemplateKey>::new(Arc::clone(&pool), &chunk_path)
                .expect("create writer");

            for (key, record) in &records {
                writer.write_record(key, record).expect("write record");
            }
            writer.finish().expect("finish writer");
        }

        let mut reader =
            GenericKeyedChunkReader::<TemplateKey>::open(&chunk_path, None).expect("open reader");

        let mut buf = Vec::new();
        let mut read_records = Vec::new();
        while let Some(key) = reader.next_record(&mut buf).expect("read record") {
            read_records.push((key, buf.clone()));
        }

        assert_eq!(records.len(), read_records.len(), "record count mismatch");
        for (i, ((expected_key, expected_data), (actual_key, actual_data))) in
            records.iter().zip(read_records.iter()).enumerate()
        {
            assert_eq!(*expected_key, *actual_key, "key mismatch at {i}");
            assert_eq!(expected_data, actual_data, "data mismatch at {i}");
        }

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }

    #[test]
    fn test_pooled_writer_empty() {
        let dir = TempDir::new().unwrap();
        let chunk_path = dir.path().join("empty_chunk.keyed");
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6));

        {
            let writer = PooledChunkWriter::<TemplateKey>::new(Arc::clone(&pool), &chunk_path)
                .expect("create writer");
            writer.finish().expect("finish empty writer");
        }

        assert!(chunk_path.exists());
        let metadata = std::fs::metadata(&chunk_path).expect("stat file");
        assert!(metadata.len() > 0, "file should not be empty (has EOF marker)");

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }

    #[test]
    #[allow(clippy::cast_possible_truncation)]
    fn test_pooled_writer_large_records() {
        let dir = TempDir::new().unwrap();
        let chunk_path = dir.path().join("large_chunk.keyed");
        let pool = Arc::new(SortWorkerPool::new(4, 1, 6));

        let records: Vec<(TemplateKey, Vec<u8>)> = (0..500)
            .map(|i| {
                let key = make_key(i);
                let record = vec![(i % 256) as u8; 1000];
                (key, record)
            })
            .collect();

        {
            let mut writer = PooledChunkWriter::<TemplateKey>::new(Arc::clone(&pool), &chunk_path)
                .expect("create writer");

            for (key, record) in &records {
                writer.write_record(key, record).expect("write record");
            }
            writer.finish().expect("finish writer");
        }

        let mut reader =
            GenericKeyedChunkReader::<TemplateKey>::open(&chunk_path, None).expect("open reader");

        let mut buf = Vec::new();
        let mut count = 0;
        while let Some(key) = reader.next_record(&mut buf).expect("read record") {
            assert_eq!(key, records[count].0, "key mismatch at {count}");
            assert_eq!(buf, records[count].1, "data mismatch at {count}");
            count += 1;
        }
        assert_eq!(count, records.len());

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }

    #[test]
    #[allow(clippy::cast_possible_truncation)]
    fn test_start_finish_and_wait() {
        // `start_finish()` returns a handle while the I/O thread runs in the background.
        // `handle.wait()` must join it and surface any errors.
        let dir = TempDir::new().unwrap();
        let chunk_path = dir.path().join("pipelined_chunk.keyed");
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6));

        let records: Vec<(TemplateKey, Vec<u8>)> =
            (0..50).map(|i| (make_key(i), vec![(i % 256) as u8; 100])).collect();

        let handle = {
            let mut writer = PooledChunkWriter::<TemplateKey>::new(Arc::clone(&pool), &chunk_path)
                .expect("create writer");
            for (key, record) in &records {
                writer.write_record(key, record).expect("write record");
            }
            writer.start_finish().expect("start_finish")
        };

        // I/O is completing in background; wait for it
        handle.wait().expect("wait should succeed");

        // Verify all records are readable back
        let mut reader =
            GenericKeyedChunkReader::<TemplateKey>::open(&chunk_path, None).expect("open reader");
        let mut buf = Vec::new();
        let mut count = 0;
        while let Some(key) = reader.next_record(&mut buf).expect("read record") {
            assert_eq!(key, records[count].0, "key mismatch at {count}");
            count += 1;
        }
        assert_eq!(count, records.len());

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }

    #[test]
    fn test_spill_write_handle_drop_without_wait() {
        // Dropping a `SpillWriteHandle` without calling `wait()` must not panic —
        // the `Drop` impl joins the thread and logs any error.
        let dir = TempDir::new().unwrap();
        let chunk_path = dir.path().join("dropped_chunk.keyed");
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6));

        let handle = {
            let mut writer = PooledChunkWriter::<TemplateKey>::new(Arc::clone(&pool), &chunk_path)
                .expect("create writer");
            writer.write_record(&make_key(0), &[1, 2, 3]).expect("write");
            writer.start_finish().expect("start_finish")
        };

        // Drop handle without calling wait() — Drop impl joins thread silently
        drop(handle);

        // File must exist and be non-empty (the I/O thread completed via Drop)
        assert!(chunk_path.exists());
        assert!(std::fs::metadata(&chunk_path).unwrap().len() > 0);

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }

    #[test]
    fn test_drop_before_finish() {
        // Dropping `PooledChunkWriter` without calling finish() must not panic or
        // deadlock — the Drop impl signals the I/O thread and joins it.
        let dir = TempDir::new().unwrap();
        let chunk_path = dir.path().join("dropped_writer.keyed");
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6));

        {
            let mut writer = PooledChunkWriter::<TemplateKey>::new(Arc::clone(&pool), &chunk_path)
                .expect("create writer");
            writer.write_record(&make_key(0), &[1, 2, 3]).expect("write");
            // Drop without calling finish() — exercises the Drop impl
        }

        // File must exist (I/O thread completed via Drop)
        assert!(chunk_path.exists());

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }
}
