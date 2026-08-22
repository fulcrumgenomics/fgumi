//! Parallel chunk writer using `SortWorkerPool` for spill compression.
//!
//! `PooledChunkWriter` replaces single-threaded `GenericKeyedChunkWriter` during spill,
//! distributing per-block compression (BGZF or zstd, picked at construction time)
//! across the shared worker pool. This targets the #1 bottleneck: spill write is
//! 62-80% of total sort wall time.
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

use crate::bgzf_io::{StagingBuffer, io_writer_loop};
use crate::codec::{SpillCodec, ZSPILL_MAGIC};
use crate::keys::RawSortKey;
use crate::worker_pool::{CompressResult, CompressTarget, PermitPool, SortWorkerPool};
use anyhow::Result;
use crossbeam_channel::bounded;
use fgumi_bgzf::BGZF_EOF;
use std::io::BufWriter;
use std::marker::PhantomData;
use std::path::Path;
use std::sync::Arc;
use std::thread::{self, JoinHandle};

/// A chunk writer that uses `SortWorkerPool` for parallel spill compression.
///
/// Records are buffered into ~64KB staging blocks, then submitted to the pool
/// for compression. An I/O thread receives compressed blocks and writes them
/// in serial order using a bounded reorder buffer.
pub struct PooledChunkWriter<K: RawSortKey> {
    /// `None` only after `start_finish()` transfers ownership to `SpillWriteHandle`.
    staging: Option<StagingBuffer>,
    /// Reusable scratch buffer for key serialization (non-embedded keys only).
    key_buf: Vec<u8>,
    /// Bytes a frame may hold, resolved from the codec once.
    ///
    /// This writer pre-flushes so a record never straddles a frame boundary --
    /// which is what lets the merge borrow most records in place instead of
    /// reassembling them. That budget must be the size the staging buffer will
    /// actually flush at: when this was `BGZF_MAX_BLOCK_SIZE` outright, it
    /// pre-flushed at 64 KiB no matter what the staging buffer was configured
    /// for, so raising the frame size changed the block count not at all
    /// (measured: 5,368,249 blocks at both 64 KiB and 256 KiB).
    frame_bytes: usize,
    io_handle: Option<JoinHandle<Result<()>>>,
    _phantom: PhantomData<K>,
}

/// Remove the BGZF end-of-stream marker from the end of a spill run.
///
/// Called before appending another chunk, so the finished run carries exactly one
/// terminator. `read_raw_blocks` would in fact skip an interior marker
/// (`fgumi-bgzf`), so this is not what makes appending correct — it keeps the run
/// well-formed for anything else that reads it, htslib included, and avoids the
/// edge where a read batch of nothing but markers looks like end-of-file.
///
/// Verifies the tail before truncating: silently lopping 28 bytes off a file that
/// does not end the way we expect would corrupt a record instead of removing a
/// marker.
fn truncate_bgzf_terminator(path: &Path) -> Result<()> {
    use std::io::{Read, Seek, SeekFrom};

    let mut file = std::fs::OpenOptions::new().read(true).write(true).open(path)?;
    let len = file.metadata()?.len();
    let marker_len = BGZF_EOF.len() as u64;
    anyhow::ensure!(
        len >= marker_len,
        "cannot append to {}: file is {len} bytes, shorter than a BGZF terminator",
        path.display()
    );

    // Seek from the end by the marker's length. `BGZF_EOF` is a 28-byte
    // constant, so this conversion cannot wrap in practice; do it explicitly
    // rather than silencing the lint.
    let marker_offset = i64::try_from(marker_len).expect("BGZF terminator length fits in i64");
    file.seek(SeekFrom::End(-marker_offset))?;
    let mut tail = vec![0u8; BGZF_EOF.len()];
    file.read_exact(&mut tail)?;
    anyhow::ensure!(
        tail == BGZF_EOF,
        "cannot append to {}: it does not end with a BGZF terminator",
        path.display()
    );

    file.set_len(len - marker_len)?;
    Ok(())
}

impl<K: RawSortKey> PooledChunkWriter<K> {
    /// Create a new pooled chunk writer with an explicit spill codec.
    ///
    /// Opens the output file and spawns an I/O writer thread that receives
    /// compressed blocks and writes them in serial order. Zstd-format spill
    /// files start with the four-byte `ZSPILL_MAGIC` so the reader can detect
    /// the format without consulting external state.
    ///
    /// # Errors
    ///
    /// Returns an error if the output file cannot be created.
    pub fn new(pool: Arc<SortWorkerPool>, path: &Path, codec: SpillCodec) -> Result<Self> {
        Self::open(pool, path, codec, false)
    }

    /// Create a run, or append another sorted chunk to an existing one.
    ///
    /// When `appending`, the caller must have established that every key in this
    /// chunk sorts at or after the run's last key; this type does not and cannot
    /// check that.
    ///
    /// Two things differ when appending. The file is opened for append rather than
    /// truncated, and the zstd format magic is written only when a run is created —
    /// it identifies the file, not each chunk within it. For BGZF the trailing
    /// end-of-stream marker left by the previous chunk is removed first, so a
    /// finished run carries exactly one terminator, at its end.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened, or if a BGZF run being
    /// appended to does not end with the terminator this writer would have
    /// written.
    pub fn open(
        pool: Arc<SortWorkerPool>,
        path: &Path,
        codec: SpillCodec,
        appending: bool,
    ) -> Result<Self> {
        let file = if appending {
            if matches!(codec, SpillCodec::Bgzf) {
                truncate_bgzf_terminator(path)?;
            }
            std::fs::OpenOptions::new().append(true).open(path)?
        } else {
            std::fs::File::create(path)?
        };
        let mut writer = BufWriter::with_capacity(256 * 1024, file);
        if matches!(codec, SpillCodec::Zstd) && !appending {
            use std::io::Write;
            writer.write_all(&ZSPILL_MAGIC)?;
        }

        let reorder_capacity = pool.num_workers() * 4;
        let (result_tx, result_rx) = bounded::<CompressResult>(reorder_capacity);
        let buffer_pool = pool.buffer_pool.clone();
        let permit_pool = Arc::new(PermitPool::new(reorder_capacity));

        let pp = Arc::clone(&permit_pool);
        let io_handle =
            thread::spawn(move || io_writer_loop(writer, result_rx, buffer_pool, pp, codec, None));

        Ok(Self {
            // `CompressTarget::Spill`: every block this writer submits is a
            // Phase 1 spill chunk and must be compressed at `temp_compression`.
            staging: Some(StagingBuffer::new(
                pool,
                result_tx,
                permit_pool,
                codec,
                CompressTarget::Spill,
            )),
            key_buf: Vec::new(),
            frame_bytes: crate::bgzf_io::spill_frame_bytes(codec),
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
            let frame = self.frame_bytes;
            if staging.buf().len() + needed > frame {
                staging.flush()?;
            }
            staging.buf().extend_from_slice(&(record.len() as u32).to_le_bytes());
            if record.len() > frame.saturating_sub(4) {
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
            // No size limit check: records larger than one frame are handled by
            // write_chunked(), which splits them across multiple blocks. The
            // reader uses streaming read_exact() that transparently spans blocks.
            if staging.buf().len() + needed > self.frame_bytes {
                staging.flush()?;
            }
            staging.buf().extend_from_slice(&self.key_buf);
            staging.buf().extend_from_slice(&(record.len() as u32).to_le_bytes());
            staging.write_chunked(record)?;
        }
        Ok(())
    }

    /// The frame budget this writer pre-flushes against.
    ///
    /// Exposed so it can be pinned to [`crate::bgzf_io::spill_frame_bytes`]
    /// rather than trusted to match it.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn frame_bytes(&self) -> usize {
        self.frame_bytes
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
    use crate::external::GenericKeyedChunkReader;
    use crate::inline::TemplateKey;
    use tempfile::TempDir;

    /// A finished run carries exactly one BGZF terminator, at its end, however
    /// many chunks were appended to it.
    ///
    /// Tested here rather than through a full sort, because a sort deletes its
    /// spill files on success — the property is only observable at this level.
    #[test]
    fn test_appending_leaves_one_bgzf_terminator() {
        let dir = TempDir::new().expect("tempdir");
        let path = dir.path().join("run.keyed");
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, SpillCodec::Bgzf, false));

        // Two chunks, the second appended to the first.
        for (chunk, appending) in [(0u64, false), (1u64, true)] {
            let mut writer = PooledChunkWriter::<TemplateKey>::open(
                Arc::clone(&pool),
                &path,
                SpillCodec::Bgzf,
                appending,
            )
            .expect("open run");
            writer.write_record(&make_key(chunk), b"record-bytes").expect("write");
            writer.finish().expect("finish");
        }

        let bytes = std::fs::read(&path).expect("read run");
        assert!(bytes.ends_with(&BGZF_EOF), "a finished run must end with the terminator");
        let interior = &bytes[..bytes.len() - BGZF_EOF.len()];
        assert!(
            !interior.windows(BGZF_EOF.len()).any(|w| w == BGZF_EOF),
            "an appended run must not carry an interior terminator"
        );

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }

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
    fn test_pooled_writer_roundtrip_zstd() {
        // Mirror of `test_pooled_writer_roundtrip` but exercising the zstd
        // spill path: the file must start with `ZSPILL_MAGIC` and round-trip
        // identical records back through `GenericKeyedChunkReader`, which
        // auto-detects the magic and routes to `ZspillStreamReader`.
        let dir = TempDir::new().unwrap();
        let chunk_path = dir.path().join("test_chunk_zstd.keyed");
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Zstd, false));

        let records: Vec<(TemplateKey, Vec<u8>)> = (0..100)
            .map(|i| {
                let key = make_key(i);
                let record = vec![(i % 256) as u8; 200 + (i as usize % 50)];
                (key, record)
            })
            .collect();

        {
            let mut writer = PooledChunkWriter::<TemplateKey>::new(
                Arc::clone(&pool),
                &chunk_path,
                SpillCodec::Zstd,
            )
            .expect("create writer");

            for (key, record) in &records {
                writer.write_record(key, record).expect("write record");
            }
            writer.finish().expect("finish writer");
        }

        // The first four bytes must be the ZSPILL magic so readers (and
        // future debugging tools) can route the file via the zstd path.
        let bytes = std::fs::read(&chunk_path).expect("read chunk file");
        assert!(bytes.len() >= ZSPILL_MAGIC.len(), "zstd spill file too short to hold magic");
        assert_eq!(
            &bytes[..ZSPILL_MAGIC.len()],
            &ZSPILL_MAGIC[..],
            "zstd spill file missing ZSPILL_MAGIC prefix"
        );

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
    #[allow(clippy::cast_possible_truncation)]
    fn test_pooled_writer_roundtrip() {
        let dir = TempDir::new().unwrap();
        let chunk_path = dir.path().join("test_chunk.keyed");
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf, false));

        let records: Vec<(TemplateKey, Vec<u8>)> = (0..100)
            .map(|i| {
                let key = make_key(i);
                let record = vec![(i % 256) as u8; 200 + (i as usize % 50)];
                (key, record)
            })
            .collect();

        {
            let mut writer = PooledChunkWriter::<TemplateKey>::new(
                Arc::clone(&pool),
                &chunk_path,
                SpillCodec::Bgzf,
            )
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
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf, false));

        {
            let writer = PooledChunkWriter::<TemplateKey>::new(
                Arc::clone(&pool),
                &chunk_path,
                SpillCodec::Bgzf,
            )
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
        let pool = Arc::new(SortWorkerPool::new(4, 1, 6, crate::codec::SpillCodec::Bgzf, false));

        let records: Vec<(TemplateKey, Vec<u8>)> = (0..500)
            .map(|i| {
                let key = make_key(i);
                let record = vec![(i % 256) as u8; 1000];
                (key, record)
            })
            .collect();

        {
            let mut writer = PooledChunkWriter::<TemplateKey>::new(
                Arc::clone(&pool),
                &chunk_path,
                SpillCodec::Bgzf,
            )
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
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf, false));

        let records: Vec<(TemplateKey, Vec<u8>)> =
            (0..50).map(|i| (make_key(i), vec![(i % 256) as u8; 100])).collect();

        let handle = {
            let mut writer = PooledChunkWriter::<TemplateKey>::new(
                Arc::clone(&pool),
                &chunk_path,
                SpillCodec::Bgzf,
            )
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
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf, false));

        let handle = {
            let mut writer = PooledChunkWriter::<TemplateKey>::new(
                Arc::clone(&pool),
                &chunk_path,
                SpillCodec::Bgzf,
            )
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
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf, false));

        {
            let mut writer = PooledChunkWriter::<TemplateKey>::new(
                Arc::clone(&pool),
                &chunk_path,
                SpillCodec::Bgzf,
            )
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
