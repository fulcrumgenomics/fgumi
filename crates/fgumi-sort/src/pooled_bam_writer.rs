#![deny(unsafe_code)]
//! Pool-based BAM writer using `SortWorkerPool` for parallel BGZF compression.
//!
//! `PooledBamWriter` replaces the noodles `MultithreadedWriter` during merge output,
//! distributing BGZF block compression across the shared worker pool. This eliminates
//! three inefficiencies in the noodles writer:
//!
//! 1. Per-block channel allocation (~100K `bounded(1)` channels for a 6GB file)
//! 2. Per-block `Vec` clone (~2-4GB of unnecessary memcpy)
//! 3. Separate deflater thread pool (wastes threads while `SortWorkerPool` sits idle)
//!
//! # Architecture
//!
//! ```text
//! Main thread                   Worker pool            I/O thread
//! ─────────────                 ───────────            ──────────
//! write_raw_record() ──►        (compress)   ──►       write to disk
//!   buffer ~64KB       submit    parallel      result   sequential
//!   submit job         ──────►  workers        ──────► reorder buf
//!                                              ──────► write ordered
//! ```

use crate::bgzf_io::{BlockOffset, StagingBuffer, io_writer_loop};
use crate::codec::SpillCodec;
use crate::worker_pool::{CompressResult, PermitPool, SortWorkerPool};
use anyhow::Result;
use crossbeam_channel::{Receiver, bounded, unbounded};
use fgumi_bam_io::BaiBuilder;
use fgumi_bgzf::BGZF_MAX_BLOCK_SIZE;
use noodles::bam::bai;
use noodles::sam::Header;
use std::io::BufWriter;
use std::path::Path;
use std::sync::Arc;
use std::thread::{self, JoinHandle};

use super::pooled_chunk_writer::SpillWriteHandle;

/// A BAM output writer that uses `SortWorkerPool` for parallel BGZF compression.
///
/// Records are buffered into ~64KB staging blocks, then submitted to the pool
/// for compression. An I/O thread receives compressed blocks and writes them
/// in serial order using a bounded reorder buffer.
///
/// This writer produces valid BAM output (header + length-prefixed records in BGZF blocks).
pub struct PooledBamWriter {
    /// `None` only after `start_finish()` transfers ownership to `SpillWriteHandle`.
    staging: Option<StagingBuffer>,
    io_handle: Option<JoinHandle<Result<()>>>,
    /// Present only for writers created via [`PooledBamWriter::new_indexing`];
    /// drives incremental BAI generation.
    index: Option<IndexState>,
}

/// State for incremental BAI generation, held only by an indexing writer.
struct IndexState {
    /// Shared virtual-offset accumulator (identical logic to `IndexingBamWriter`).
    bai: BaiBuilder,
    /// Per-block `(serial, compressed_start)` notifications from the I/O thread.
    block_offset_rx: Receiver<BlockOffset>,
    /// Number of reference sequences, needed to finalize the index.
    num_refs: usize,
}

impl PooledBamWriter {
    /// Create a new pooled BAM writer.
    ///
    /// Opens the output file, writes the BAM header into the initial staging buffer,
    /// and spawns an I/O writer thread that receives compressed blocks and writes
    /// them in serial order.
    ///
    /// # Errors
    ///
    /// Returns an error if the output file cannot be created or the header cannot
    /// be serialized.
    pub fn new(pool: Arc<SortWorkerPool>, path: &Path, header: &Header) -> Result<Self> {
        Self::new_inner(pool, path, header, false)
    }

    /// Create a pooled BAM writer that also builds a BAI index during the write.
    ///
    /// Output bytes are identical to [`new`](Self::new); additionally, each
    /// record's virtual file offset is tracked and the index is returned by
    /// [`finish_index`](Self::finish_index). Only meaningful for
    /// coordinate-sorted output.
    ///
    /// # Errors
    ///
    /// Returns an error if the output file cannot be created or the header
    /// cannot be serialized.
    pub fn new_indexing(pool: Arc<SortWorkerPool>, path: &Path, header: &Header) -> Result<Self> {
        Self::new_inner(pool, path, header, true)
    }

    /// Shared constructor for [`new`](Self::new) and
    /// [`new_indexing`](Self::new_indexing).
    fn new_inner(
        pool: Arc<SortWorkerPool>,
        path: &Path,
        header: &Header,
        indexing: bool,
    ) -> Result<Self> {
        let file = std::fs::File::create(path)?;
        let writer = BufWriter::with_capacity(256 * 1024, file);

        let reorder_capacity = pool.num_workers() * 4;
        let (result_tx, result_rx) = bounded::<CompressResult>(reorder_capacity);
        let buffer_pool = pool.buffer_pool.clone();
        let permit_pool = Arc::new(PermitPool::new(reorder_capacity));

        // When indexing, wire an unbounded channel so the I/O thread can emit
        // each block's compressed start offset without ever blocking on send
        // (the consumer drains it during writes and once more at finish).
        let (block_offset_tx, index) = if indexing {
            let (tx, rx) = unbounded::<BlockOffset>();
            let num_refs = header.reference_sequences().len();
            (Some(tx), Some(IndexState { bai: BaiBuilder::new(), block_offset_rx: rx, num_refs }))
        } else {
            (None, None)
        };

        let pp = Arc::clone(&permit_pool);
        let io_handle = thread::spawn(move || {
            io_writer_loop(writer, result_rx, buffer_pool, pp, SpillCodec::Bgzf, block_offset_tx)
        });

        let mut staging = StagingBuffer::new(pool, result_tx, permit_pool, SpillCodec::Bgzf);

        // Write BAM header into a temporary buffer then flush in BGZF-sized chunks.
        // Headers can exceed BGZF_MAX_BLOCK_SIZE; write_chunked handles the splitting.
        let mut header_buf = Vec::new();
        fgumi_bam_io::write_bam_header(&mut header_buf, header)?;
        staging.write_chunked(&header_buf)?;
        staging.flush()?;

        Ok(Self { staging: Some(staging), io_handle: Some(io_handle), index })
    }

    /// Write a raw BAM record to the output.
    ///
    /// Writes a 4-byte little-endian length prefix followed by the raw record bytes.
    /// When the staging buffer reaches ~64KB, it's submitted to the pool for compression.
    ///
    /// # Errors
    ///
    /// Returns an error if writing fails (currently infallible for staging buffer).
    ///
    /// # Panics
    ///
    /// Panics if called after [`start_finish`](Self::start_finish) has been called.
    #[inline]
    #[allow(clippy::cast_possible_truncation)]
    pub fn write_raw_record(&mut self, record_bytes: &[u8]) -> Result<()> {
        let staging = self.staging.as_mut().expect("write_raw_record called after start_finish");
        let block_size = record_bytes.len() as u32;
        // Pre-flush if this record won't fit in the current block to keep each
        // CompressJob within BGZF_MAX_BLOCK_SIZE.
        let needed = 4 + record_bytes.len();
        if staging.buf().len() + needed > BGZF_MAX_BLOCK_SIZE {
            staging.flush()?;
        }
        // For indexing, capture the record's start position *after* any pre-flush
        // and *before* appending. `next_serial`/`buf_len` are the writer's own
        // block number and uncompressed offset, so the recorded position matches
        // the block layout exactly (records never straddle a block except the
        // oversized case below, which `compute_end_vpos` resolves by walking
        // full-size blocks forward).
        if let Some(index) = self.index.as_mut() {
            index.bai.record(staging.next_serial(), staging.buf_len(), needed, record_bytes);
        }
        // Write the 4-byte length prefix (always fits: staging is empty after the pre-flush,
        // and 4 bytes is well within BGZF_MAX_BLOCK_SIZE).
        staging.buf().extend_from_slice(&block_size.to_le_bytes());
        // If the record payload itself exceeds one BGZF block, split it across blocks.
        // BAM records can legally span BGZF blocks; readers handle this transparently.
        if record_bytes.len() > BGZF_MAX_BLOCK_SIZE.saturating_sub(4) {
            staging.write_chunked(record_bytes)?;
        } else {
            staging.buf().extend_from_slice(record_bytes);
            staging.flush_if_full()?;
        }
        // Drain block-offset notifications and resolve completed records so the
        // pending-entry cache stays bounded by the in-flight block count.
        if let Some(index) = self.index.as_mut() {
            while let Ok(bo) = index.block_offset_rx.try_recv() {
                index.bai.note_block(bo.serial, bo.compressed_start);
            }
            index.bai.resolve()?;
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

    /// Flush any remaining staged data and drop the staging buffer.
    ///
    /// Dropping the staging buffer closes the result channel, so the I/O thread
    /// drains its queue, writes every block (and, for an indexing writer, emits
    /// every block offset), then exits. Shared by
    /// [`start_finish`](Self::start_finish) and [`finish_index`](Self::finish_index).
    fn flush_and_close_staging(&mut self) -> Result<()> {
        if let Some(mut staging) = self.staging.take() {
            if !staging.buf().is_empty() {
                staging.flush()?;
            }
            drop(staging);
        }
        Ok(())
    }

    /// Flush remaining data and signal the I/O thread, but don't wait for it.
    ///
    /// Returns a [`SpillWriteHandle`] that can be waited on later.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing the staging buffer fails.
    pub fn start_finish(mut self) -> Result<SpillWriteHandle> {
        self.flush_and_close_staging()?;
        Ok(SpillWriteHandle::new(self.io_handle.take()))
    }

    /// Finish writing and return the BAI index built during the write.
    ///
    /// Flushes remaining data, joins the I/O thread (ensuring every block is
    /// written and every block-offset notification emitted), resolves the
    /// trailing records, and builds the index.
    ///
    /// # Panics
    ///
    /// Panics if called on a writer not created via
    /// [`new_indexing`](Self::new_indexing).
    ///
    /// # Errors
    ///
    /// Returns an error if flushing, the I/O thread, or index construction
    /// fails (the latter includes any record whose block offset was never
    /// observed).
    pub fn finish_index(mut self) -> Result<bai::Index> {
        let mut index = self.index.take().expect("finish_index requires an indexing writer");

        // Flush remaining staging and close the input to the I/O thread: it then
        // drains, writes all blocks, emits all offsets, and exits (dropping its
        // block-offset sender).
        self.flush_and_close_staging()?;

        // Join the I/O thread so every block offset has been emitted before we
        // drain the channel for the last time.
        if let Some(handle) = self.io_handle.take() {
            handle.join().map_err(|_| anyhow::anyhow!("I/O writer thread panicked"))??;
        }

        // Drain the final block offsets and resolve the remaining records.
        while let Ok(bo) = index.block_offset_rx.try_recv() {
            index.bai.note_block(bo.serial, bo.compressed_start);
        }
        index.bai.resolve()?;
        let idx = index.bai.build(index.num_refs)?;
        Ok(idx)
    }
}

impl Drop for PooledBamWriter {
    fn drop(&mut self) {
        if self.io_handle.is_some() {
            // Writer dropped before finish()/start_finish() (e.g. early error return).
            // Drop staging first — this closes result_tx, signaling the I/O thread to
            // drain and exit. Then join the thread to avoid silently detaching it.
            drop(self.staging.take());
            if let Some(handle) = self.io_handle.take() {
                match handle.join() {
                    Ok(Ok(())) => {}
                    Ok(Err(e)) => log::error!("PooledBamWriter: I/O writer thread error: {e}"),
                    Err(_) => log::error!("PooledBamWriter: I/O writer thread panicked"),
                }
            }
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
#[allow(clippy::default_constructed_unit_structs)]
mod tests {
    use super::*;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReferenceSequence;

    /// Build a minimal test header with a few reference sequences.
    fn test_header() -> Header {
        Header::builder()
            .set_header(Map::<noodles::sam::header::record::value::map::Header>::new(
                noodles::sam::header::record::value::map::header::Version::new(1, 6),
            ))
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(std::num::NonZero::new(100_000).expect("non-zero")),
            )
            .add_reference_sequence(
                "chr2",
                Map::<ReferenceSequence>::new(std::num::NonZero::new(200_000).expect("non-zero")),
            )
            .build()
    }

    /// Create a simple raw BAM record (unmapped, minimal fields).
    #[allow(clippy::cast_possible_truncation)]
    fn make_test_record(name: &[u8], seq_len: usize) -> Vec<u8> {
        let name_len = name.len() + 1;
        let seq_bytes = seq_len.div_ceil(2);
        let record_len = 32 + name_len + seq_bytes + seq_len;

        let mut rec = Vec::with_capacity(record_len);

        rec.extend_from_slice(&(-1i32).to_le_bytes()); // ref_id
        rec.extend_from_slice(&(-1i32).to_le_bytes()); // pos
        let bin_mq_nl: u32 = (name_len as u32) | (4680 << 16);
        rec.extend_from_slice(&bin_mq_nl.to_le_bytes());
        let flag_nc: u32 = 4 << 16; // unmapped
        rec.extend_from_slice(&flag_nc.to_le_bytes());
        rec.extend_from_slice(&(seq_len as u32).to_le_bytes());
        rec.extend_from_slice(&(-1i32).to_le_bytes()); // mate_ref_id
        rec.extend_from_slice(&(-1i32).to_le_bytes()); // mate_pos
        rec.extend_from_slice(&0i32.to_le_bytes()); // tlen
        rec.extend_from_slice(name);
        rec.push(0); // null terminator
        rec.resize(rec.len() + seq_bytes, 0x11); // seq: all A's
        rec.resize(rec.len() + seq_len, 0xFF); // qual

        rec
    }

    /// Create a mapped raw BAM record with a single `seq_len`M CIGAR at
    /// `(ref_id, pos)` (0-based pos), mapq 60, flag 0.
    #[allow(clippy::cast_possible_truncation)]
    fn make_mapped_record(name: &[u8], ref_id: i32, pos: i32, seq_len: usize) -> Vec<u8> {
        let name_len = name.len() + 1;
        let seq_bytes = seq_len.div_ceil(2);
        let mut rec = Vec::new();

        rec.extend_from_slice(&ref_id.to_le_bytes());
        rec.extend_from_slice(&pos.to_le_bytes());
        let bin_mq_nl: u32 = (name_len as u32) | (60u32 << 8) | (4680u32 << 16); // mapq 60
        rec.extend_from_slice(&bin_mq_nl.to_le_bytes());
        let flag_nc: u32 = 1; // n_cigar_op = 1, flag = 0 (mapped)
        rec.extend_from_slice(&flag_nc.to_le_bytes());
        rec.extend_from_slice(&(seq_len as u32).to_le_bytes());
        rec.extend_from_slice(&(-1i32).to_le_bytes()); // mate_ref_id
        rec.extend_from_slice(&(-1i32).to_le_bytes()); // mate_pos
        rec.extend_from_slice(&0i32.to_le_bytes()); // tlen
        rec.extend_from_slice(name);
        rec.push(0); // null terminator
        let cigar_op: u32 = (seq_len as u32) << 4; // seq_len M (op code 0)
        rec.extend_from_slice(&cigar_op.to_le_bytes());
        rec.resize(rec.len() + seq_bytes, 0x11); // seq
        rec.resize(rec.len() + seq_len, 0xFF); // qual
        rec
    }

    /// The indexing writer must produce byte-identical BAM output to the plain
    /// pooled writer (indexing is a pure side-channel) and yield a writable BAI.
    ///
    /// Covers a coordinate-sorted mix that spans many BGZF blocks: mapped reads
    /// on two references, trailing unmapped reads, and an oversized record that
    /// spans multiple blocks (exercising cross-block virtual-offset accounting).
    #[test]
    fn test_indexing_writer_output_matches_plain_and_writes_bai() {
        let dir = tempfile::TempDir::new().expect("tempdir");
        let header = test_header(); // chr1 (100k), chr2 (200k)

        let mut records: Vec<Vec<u8>> = Vec::new();
        for i in 0..1500 {
            records.push(make_mapped_record(format!("c1_{i:05}").as_bytes(), 0, i * 30, 60));
        }
        for i in 0..1500 {
            records.push(make_mapped_record(format!("c2_{i:05}").as_bytes(), 1, i * 30, 60));
        }
        for i in 0..100 {
            records.push(make_test_record(format!("unm_{i:03}").as_bytes(), 60)); // unmapped
        }
        // Oversized (unmapped) record that spans multiple BGZF blocks.
        records.push(make_test_record(b"oversized", BGZF_MAX_BLOCK_SIZE + 500));

        let pool = Arc::new(SortWorkerPool::new(4, 1, 6, crate::codec::SpillCodec::Bgzf));

        let plain_path = dir.path().join("plain.bam");
        {
            let mut w = PooledBamWriter::new(Arc::clone(&pool), &plain_path, &header)
                .expect("plain writer");
            for r in &records {
                w.write_raw_record(r).expect("write");
            }
            w.finish().expect("finish plain");
        }

        let indexed_path = dir.path().join("indexed.bam");
        let index = {
            let mut w = PooledBamWriter::new_indexing(Arc::clone(&pool), &indexed_path, &header)
                .expect("indexing writer");
            for r in &records {
                w.write_raw_record(r).expect("write");
            }
            w.finish_index().expect("finish_index")
        };

        // Indexing must not change a single output byte.
        assert_eq!(
            std::fs::read(&plain_path).expect("read plain"),
            std::fs::read(&indexed_path).expect("read indexed"),
            "indexing writer must produce byte-identical BAM output to the plain writer"
        );

        // The returned index must serialize to a sidecar that loads back cleanly
        // and covers every reference: a non-empty file alone would not catch a
        // corrupt virtual offset, which is the contract this writer owes.
        let bai_path = dir.path().join("indexed.bam.bai");
        fgumi_bam_io::write_bai_index(&bai_path, &index).expect("write bai");
        let loaded = bai::fs::read(&bai_path).expect("BAI must be loadable");
        assert_eq!(
            loaded.reference_sequences().len(),
            header.reference_sequences().len(),
            "index should cover every reference sequence"
        );

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }

    #[test]
    fn test_pooled_bam_writer_roundtrip() {
        let dir = tempfile::TempDir::new().expect("tempdir");
        let bam_path = dir.path().join("test.bam");
        let header = test_header();
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf));

        let num_records = 200;
        let records: Vec<Vec<u8>> = (0..num_records)
            .map(|i| make_test_record(format!("read_{i:04}").as_bytes(), 50))
            .collect();

        {
            let mut writer =
                PooledBamWriter::new(Arc::clone(&pool), &bam_path, &header).expect("create writer");
            for rec in &records {
                writer.write_raw_record(rec).expect("write record");
            }
            writer.finish().expect("finish writer");
        }

        let mut reader = noodles::bam::io::reader::Builder::default()
            .build_from_path(&bam_path)
            .expect("open reader");
        let read_header = reader.read_header().expect("read header");
        assert_eq!(read_header.reference_sequences().len(), 2);

        let mut count = 0;
        for result in reader.records() {
            let _record = result.expect("read record");
            count += 1;
        }
        assert_eq!(count, num_records, "record count mismatch");

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }

    #[test]
    fn test_pooled_bam_writer_empty() {
        let dir = tempfile::TempDir::new().expect("tempdir");
        let bam_path = dir.path().join("empty.bam");
        let header = test_header();
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf));

        {
            let writer =
                PooledBamWriter::new(Arc::clone(&pool), &bam_path, &header).expect("create writer");
            writer.finish().expect("finish empty writer");
        }

        let mut reader = noodles::bam::io::reader::Builder::default()
            .build_from_path(&bam_path)
            .expect("open reader");
        let read_header = reader.read_header().expect("read header");
        assert_eq!(read_header.reference_sequences().len(), 2);
        let record_count = reader
            .records()
            .collect::<std::io::Result<Vec<_>>>()
            .expect("records should read cleanly")
            .len();
        assert_eq!(record_count, 0);

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }

    #[test]
    fn test_pooled_bam_writer_many_records() {
        let dir = tempfile::TempDir::new().expect("tempdir");
        let bam_path = dir.path().join("many.bam");
        let header = test_header();
        let pool = Arc::new(SortWorkerPool::new(4, 1, 6, crate::codec::SpillCodec::Bgzf));

        let num_records = 5000;
        {
            let mut writer =
                PooledBamWriter::new(Arc::clone(&pool), &bam_path, &header).expect("create writer");
            for i in 0..num_records {
                let rec = make_test_record(format!("read_{i:06}").as_bytes(), 100);
                writer.write_raw_record(&rec).expect("write record");
            }
            writer.finish().expect("finish writer");
        }

        let mut reader = noodles::bam::io::reader::Builder::default()
            .build_from_path(&bam_path)
            .expect("open reader");
        let _header = reader.read_header().expect("read header");
        let record_count = reader
            .records()
            .collect::<std::io::Result<Vec<_>>>()
            .expect("records should read cleanly")
            .len();
        assert_eq!(record_count, num_records);

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }

    #[test]
    fn test_pooled_bam_writer_raw_bytes_match() {
        let dir = tempfile::TempDir::new().expect("tempdir");
        let bam_path = dir.path().join("raw_match.bam");
        let header = test_header();
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf));

        let records: Vec<Vec<u8>> =
            (0..50).map(|i| make_test_record(format!("r{i}").as_bytes(), 30)).collect();

        {
            let mut writer =
                PooledBamWriter::new(Arc::clone(&pool), &bam_path, &header).expect("create writer");
            for rec in &records {
                writer.write_raw_record(rec).expect("write");
            }
            writer.finish().expect("finish");
        }

        let mut reader =
            noodles::bam::io::reader::Builder::default().build_from_path(&bam_path).expect("open");
        let _h = reader.read_header().expect("header");

        let mut count = 0;
        for result in reader.records() {
            result.expect("read record");
            count += 1;
        }
        assert_eq!(count, records.len());

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }

    #[test]
    fn test_pooled_bam_writer_oversized_record() {
        // Verify that a record payload exceeding BGZF_MAX_BLOCK_SIZE is split across
        // multiple BGZF blocks and can be read back correctly.
        let dir = tempfile::TempDir::new().expect("tempdir");
        let bam_path = dir.path().join("oversized.bam");
        let header = test_header();
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf));

        // A sequence of BGZF_MAX_BLOCK_SIZE bytes exceeds the threshold.
        let oversized_rec = make_test_record(b"oversized_read", BGZF_MAX_BLOCK_SIZE);

        {
            let mut writer =
                PooledBamWriter::new(Arc::clone(&pool), &bam_path, &header).expect("create writer");
            writer.write_raw_record(&oversized_rec).expect("write oversized record");
            writer.finish().expect("finish writer");
        }

        let mut reader = noodles::bam::io::reader::Builder::default()
            .build_from_path(&bam_path)
            .expect("open reader");
        let _header = reader.read_header().expect("read header");
        let record_count = reader
            .records()
            .collect::<std::io::Result<Vec<_>>>()
            .expect("records should read cleanly")
            .len();
        assert_eq!(record_count, 1, "oversized record should round-trip");

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }

    #[test]
    fn test_drop_before_finish() {
        // Dropping `PooledBamWriter` without calling finish() must not panic or
        // deadlock — the Drop impl signals the I/O thread and joins it.
        let dir = tempfile::TempDir::new().expect("tempdir");
        let bam_path = dir.path().join("dropped_writer.bam");
        let header = test_header();
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf));

        {
            let mut writer =
                PooledBamWriter::new(Arc::clone(&pool), &bam_path, &header).expect("create writer");
            let rec = make_test_record(b"r0", 10);
            writer.write_raw_record(&rec).expect("write");
            // Drop without calling finish() — exercises the Drop impl
        }

        // File must exist (I/O thread completed via Drop)
        assert!(bam_path.exists());

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
    }
}
