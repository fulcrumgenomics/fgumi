//! Shared BGZF I/O utilities for pooled writers.
//!
//! Provides the reorder-and-write loop used by both [`PooledBamWriter`](super::pooled_bam_writer)
//! and [`PooledChunkWriter`](super::pooled_chunk_writer), and the staging buffer logic for
//! accumulating data into ~64KB blocks before submitting compression jobs.

use crate::codec::SpillCodec;
use crate::worker_pool::{
    BufferPool, CompressJob, CompressResult, CompressTarget, PermitPool, SortWorkerPool,
};
use anyhow::Result;
use crossbeam_channel::{Receiver, Sender};
use fgumi_bgzf::{BGZF_EOF, BGZF_MAX_BLOCK_SIZE};
use std::collections::BTreeMap;
use std::io::{BufWriter, Write};
use std::sync::Arc;

/// Padding beyond `BGZF_MAX_BLOCK_SIZE` for the staging buffer capacity.
const STAGING_PADDING: usize = 4096;

/// Per-block position notification emitted by the I/O writer loop when index
/// generation is enabled.
///
/// `serial` is the block's serial number (issued by the writer's
/// [`PermitPool`] at flush time, and equal to its block number since one compress job produces
/// exactly one BGZF block). `compressed_start` is the cumulative number of
/// compressed bytes written to the file *before* this block — i.e. its on-disk
/// byte offset. The pooled indexing writer uses these to resolve BAI virtual
/// offsets.
#[derive(Debug, Clone, Copy)]
pub(crate) struct BlockOffset {
    pub serial: u64,
    pub compressed_start: u64,
}

/// Staging buffer that accumulates data and submits full blocks to the pool.
pub(crate) struct StagingBuffer {
    pool: Arc<SortWorkerPool>,
    buf: Vec<u8>,
    result_tx: Sender<CompressResult>,
    permit_pool: Arc<PermitPool>,
    codec: SpillCodec,
    /// Stamped onto every job this buffer submits, so the worker that pops it
    /// compresses at the level this writer asked for regardless of the pool's
    /// phase at that moment.
    target: CompressTarget,
}

impl StagingBuffer {
    /// Seconds this buffer's producer spent blocked waiting for an output
    /// permit, and the number of waits. See [`PermitPool::blocked`].
    pub(crate) fn write_backpressure(&self) -> (f64, u64) {
        self.permit_pool.blocked()
    }

    /// Create a new staging buffer.
    #[must_use]
    pub(crate) fn new(
        pool: Arc<SortWorkerPool>,
        result_tx: Sender<CompressResult>,
        permit_pool: Arc<PermitPool>,
        codec: SpillCodec,
        target: CompressTarget,
    ) -> Self {
        Self {
            pool,
            buf: Vec::with_capacity(BGZF_MAX_BLOCK_SIZE + STAGING_PADDING),
            result_tx,
            permit_pool,
            codec,
            target,
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

    /// Current uncompressed length of the pending (not-yet-flushed) block.
    ///
    /// This is the uncompressed offset at which the next appended byte will
    /// land in the current block — used by the indexing writer to record where
    /// a record starts.
    #[inline]
    pub(crate) fn buf_len(&self) -> usize {
        self.buf.len()
    }

    /// Serial number the pending block will be assigned when flushed.
    ///
    /// Because one compress job produces exactly one BGZF block, this is also
    /// the block number the indexing writer pairs with [`BlockOffset`].
    #[inline]
    pub(crate) fn next_serial(&self) -> u64 {
        self.permit_pool.submitted()
    }

    /// Flush the staging buffer: swap it with a recycled buffer and submit for compression.
    ///
    /// Acquires a permit from the pool before submitting, blocking if the reorder
    /// budget is exhausted. This bounds the number of in-flight compressed blocks to
    /// the pool capacity, preventing unbounded reorder buffer growth.
    ///
    /// No-op when the buffer is empty (avoids submitting empty BGZF blocks).
    ///
    /// # Errors
    ///
    /// Returns an error if the permit pool has been closed (I/O writer exited).
    pub(crate) fn flush(&mut self) -> anyhow::Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        self.permit_pool.acquire()?;

        let data = std::mem::replace(&mut self.buf, self.pool.buffer_pool.checkout());
        if self.buf.capacity() < BGZF_MAX_BLOCK_SIZE + STAGING_PADDING {
            self.buf.reserve(BGZF_MAX_BLOCK_SIZE + STAGING_PADDING - self.buf.capacity());
        }

        // The serial comes from the permit pool, which counts it as submitted,
        // so the I/O writer can tell a lost final block from a clean end.
        let serial = self.permit_pool.issue_serial();
        self.pool.submit_compress(CompressJob {
            data,
            serial,
            result_tx: self.result_tx.clone(),
            codec: self.codec,
            target: self.target,
        });
        Ok(())
    }

    /// Flush if full, otherwise no-op.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`flush`](Self::flush).
    #[inline]
    pub(crate) fn flush_if_full(&mut self) -> anyhow::Result<()> {
        if self.is_full() { self.flush() } else { Ok(()) }
    }

    /// Write `data` to the staging buffer, flushing BGZF-sized chunks as they fill up.
    ///
    /// Unlike writing directly to `buf()`, this correctly handles data larger than
    /// `BGZF_MAX_BLOCK_SIZE` (e.g. large BAM headers) by splitting into multiple jobs.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`flush`](Self::flush).
    pub(crate) fn write_chunked(&mut self, data: &[u8]) -> anyhow::Result<()> {
        let mut remaining = data;
        while !remaining.is_empty() {
            let space = BGZF_MAX_BLOCK_SIZE.saturating_sub(self.buf.len());
            let n = remaining.len().min(space);
            self.buf.extend_from_slice(&remaining[..n]);
            remaining = &remaining[n..];
            self.flush_if_full()?;
        }
        Ok(())
    }
}

/// Write one output block in serial order: emit its compressed start offset
/// (when indexing is enabled), write the bytes, advance the running compressed
/// offset, and release one reorder permit.
///
/// `compressed_start` is captured *before* the write, so it is the on-disk byte
/// offset at which this block begins. The BGZF EOF marker is written separately
/// and is intentionally never passed here (no record references it).
fn write_block_in_order<W: Write>(
    writer: &mut BufWriter<W>,
    serial: u64,
    data: &[u8],
    compressed_offset: &mut u64,
    block_offset_tx: Option<&Sender<BlockOffset>>,
    permit_pool: &Arc<PermitPool>,
) -> Result<()> {
    if let Some(tx) = block_offset_tx {
        // Best-effort: the indexing consumer keeps the receiver alive until finish.
        let _ = tx.send(BlockOffset { serial, compressed_start: *compressed_offset });
    }
    writer.write_all(data)?;
    *compressed_offset += data.len() as u64;
    permit_pool.release();
    Ok(())
}

/// I/O writer loop: receives compressed blocks and writes them in serial order.
///
/// Uses a `BTreeMap` as a reorder buffer. When the next expected serial arrives,
/// writes it immediately. Out-of-order blocks are buffered until their turn.
/// Releases one permit to `permit_pool` after each block is written out,
/// unblocking the corresponding `StagingBuffer::flush()` call and bounding the
/// number of in-flight compressed blocks to the pool capacity.
/// Writes BGZF EOF marker and flushes when all blocks are received.
///
/// Generic over the sink rather than fixed to `File`: spill chunks are always
/// files, but the sort's *output* may be stdout, which reaches here as the
/// boxed writer `open_output_writer` hands back.
///
/// When `block_offset_tx` is `Some`, each written block's `(serial,
/// compressed_start)` is emitted on it (in strict block order) for BAI virtual
/// offset resolution; when `None`, this is a no-op with zero overhead.
///
/// # Errors
///
/// Returns an error if any write fails, if a compressed block is missing
/// (which would silently truncate the output) — including the final one,
/// detected by comparing the blocks written with [`PermitPool::submitted`] —
/// or if a block arrives with a serial `permit_pool` did not issue.
#[allow(clippy::needless_pass_by_value)]
pub(crate) fn io_writer_loop<W: Write>(
    mut writer: BufWriter<W>,
    result_rx: Receiver<CompressResult>,
    buffer_pool: BufferPool,
    permit_pool: Arc<PermitPool>,
    codec: SpillCodec,
    block_offset_tx: Option<Sender<BlockOffset>>,
) -> Result<()> {
    let result = io_writer_loop_inner(
        &mut writer,
        &result_rx,
        &buffer_pool,
        &permit_pool,
        codec,
        block_offset_tx.as_ref(),
    );
    if result.is_err() {
        // Unblock any producers waiting on acquire() so they don't park forever.
        permit_pool.close();
    }
    result
}

fn io_writer_loop_inner<W: Write>(
    writer: &mut BufWriter<W>,
    result_rx: &Receiver<CompressResult>,
    buffer_pool: &BufferPool,
    permit_pool: &Arc<PermitPool>,
    codec: SpillCodec,
    block_offset_tx: Option<&Sender<BlockOffset>>,
) -> Result<()> {
    let mut next_expected: u64 = 0;
    let mut reorder_buf: BTreeMap<u64, Vec<u8>> = BTreeMap::new();
    let mut compressed_offset: u64 = 0;
    let tx = block_offset_tx;

    while let Ok(result) = result_rx.recv() {
        buffer_pool.checkin(result.recycled_buf);

        // Serials are issued (and counted) by the permit pool before the job is
        // queued, so one at or past the count was never counted, and the
        // lost-final-block check below would be blind to it.
        let issued = permit_pool.submitted();
        if result.serial >= issued {
            return Err(anyhow::anyhow!(
                "compressed block {} was not issued by the permit pool ({issued} issued)",
                result.serial
            ));
        }

        if result.serial == next_expected {
            write_block_in_order(
                writer,
                next_expected,
                &result.compressed,
                &mut compressed_offset,
                tx,
                permit_pool,
            )?;
            next_expected += 1;

            while let Some(data) = reorder_buf.remove(&next_expected) {
                write_block_in_order(
                    writer,
                    next_expected,
                    &data,
                    &mut compressed_offset,
                    tx,
                    permit_pool,
                )?;
                next_expected += 1;
            }
        } else {
            reorder_buf.insert(result.serial, result.compressed);
            // Permit held: released when this block is written in the cascade above.
        }
    }

    // Drain remaining buffered blocks — any gap means a worker dropped a result.
    while let Some((&serial, _)) = reorder_buf.first_key_value() {
        if serial == next_expected {
            let data = reorder_buf.remove(&serial).expect("key just checked");
            write_block_in_order(
                writer,
                next_expected,
                &data,
                &mut compressed_offset,
                tx,
                permit_pool,
            )?;
            next_expected += 1;
        } else {
            return Err(anyhow::anyhow!(
                "missing compressed block {next_expected}: next available is {serial}; \
                 the output would be silently truncated"
            ));
        }
    }

    // A lost final block leaves no gap above: its job was abandoned (pool shut
    // down while it was queued or compressing, or its worker panicked), its
    // result sender dropped, and the input simply closed. Only the count of
    // issued serials can tell that apart from a clean end of stream.
    let submitted = permit_pool.submitted();
    if next_expected < submitted {
        return Err(anyhow::anyhow!(
            "missing compressed block {next_expected} of {submitted} submitted; \
             the output would be silently truncated"
        ));
    }

    if matches!(codec, SpillCodec::Bgzf) {
        writer.write_all(&BGZF_EOF)?;
    }
    writer.flush()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::sync::Arc;
    use tempfile::TempDir;

    fn make_permit_pool(pool: &Arc<SortWorkerPool>) -> Arc<PermitPool> {
        Arc::new(PermitPool::new(pool.num_workers() * 4))
    }

    /// Build a round-trip helper: write `data` via `StagingBuffer` → `io_writer_loop` → read back raw bytes.
    ///
    /// `io_writer_loop` is the unit under test, so the output is the raw stream
    /// produced by the loop for the given codec — for zstd that means the
    /// sequence of `[u32 LE frame-len][zstd frame]` blocks the worker emits,
    /// without the `ZSPILL_MAGIC` prefix (which a real chunk writer would write
    /// before invoking the loop).
    fn roundtrip_data(data: &[u8], codec: SpillCodec) -> Vec<u8> {
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, codec));
        let (result_tx, result_rx) = pool.compress_result_channel();
        let buffer_pool = pool.buffer_pool.clone();
        let permit_pool = make_permit_pool(&pool);

        let dir = TempDir::new().unwrap();
        let out_path = dir.path().join("out.spill");

        let out_file = std::fs::File::create(&out_path).unwrap();
        let writer = std::io::BufWriter::new(out_file);
        let pp = Arc::clone(&permit_pool);
        let io_handle = std::thread::spawn(move || {
            io_writer_loop(writer, result_rx, buffer_pool, pp, codec, None)
        });

        let mut staging = StagingBuffer::new(
            Arc::clone(&pool),
            result_tx,
            permit_pool,
            codec,
            CompressTarget::Spill,
        );
        staging.write_chunked(data).unwrap();
        staging.flush().unwrap();
        drop(staging); // closes result_tx senders → io_writer_loop exits

        io_handle.join().unwrap().unwrap();

        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
        std::fs::read(&out_path).unwrap()
    }

    #[rstest]
    #[case(SpillCodec::Bgzf)]
    #[case(SpillCodec::Zstd)]
    fn test_staging_buffer_flush_empty_is_noop(#[case] codec: SpillCodec) {
        let pool = Arc::new(SortWorkerPool::new(1, 1, 6, codec));
        let (result_tx, _result_rx) = pool.compress_result_channel();
        let permit_pool = make_permit_pool(&pool);

        let mut staging = StagingBuffer::new(
            Arc::clone(&pool),
            result_tx,
            permit_pool,
            codec,
            CompressTarget::Spill,
        );
        // Flush with empty buffer: should not submit a compress job
        staging.flush().unwrap();

        assert_eq!(
            pool.stats.compress_jobs_submitted.load(std::sync::atomic::Ordering::Relaxed),
            0
        );
        assert_eq!(staging.next_serial(), 0, "an empty flush issues no serial");

        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
    }

    #[rstest]
    #[case(SpillCodec::Bgzf)]
    #[case(SpillCodec::Zstd)]
    fn test_staging_buffer_is_full(#[case] codec: SpillCodec) {
        let pool = Arc::new(SortWorkerPool::new(1, 1, 6, codec));
        let (result_tx, _result_rx) = pool.compress_result_channel();
        let permit_pool = make_permit_pool(&pool);
        let mut staging = StagingBuffer::new(
            Arc::clone(&pool),
            result_tx,
            permit_pool,
            codec,
            CompressTarget::Spill,
        );

        assert!(!staging.is_full(), "empty buffer should not be full");
        staging.buf().extend(vec![0u8; BGZF_MAX_BLOCK_SIZE]);
        assert!(staging.is_full(), "buffer at BGZF_MAX_BLOCK_SIZE should be full");

        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
    }

    #[rstest]
    #[case(SpillCodec::Bgzf)]
    #[case(SpillCodec::Zstd)]
    fn test_staging_buffer_write_chunked_large_data(#[case] codec: SpillCodec) {
        // Data larger than BGZF_MAX_BLOCK_SIZE must be split into multiple compress jobs.
        let large = vec![b'A'; BGZF_MAX_BLOCK_SIZE * 2 + 1000];
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, codec));
        let (result_tx, result_rx) = pool.compress_result_channel();
        let buffer_pool = pool.buffer_pool.clone();
        let permit_pool = make_permit_pool(&pool);

        let dir = TempDir::new().unwrap();
        let out_path = dir.path().join("large.spill");
        let out_file = std::fs::File::create(&out_path).unwrap();
        let writer = std::io::BufWriter::new(out_file);
        let pp = Arc::clone(&permit_pool);
        let io_handle = std::thread::spawn(move || {
            io_writer_loop(writer, result_rx, buffer_pool, pp, codec, None)
        });

        let mut staging = StagingBuffer::new(
            Arc::clone(&pool),
            result_tx,
            permit_pool,
            codec,
            CompressTarget::Spill,
        );
        staging.write_chunked(&large).unwrap();
        staging.flush().unwrap();
        // Every submitted job was issued a serial, so the writer can detect a
        // lost tail.
        assert_eq!(
            staging.next_serial(),
            pool.stats.compress_jobs_submitted.load(std::sync::atomic::Ordering::Relaxed),
        );
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

    #[rstest]
    #[case(SpillCodec::Bgzf)]
    #[case(SpillCodec::Zstd)]
    fn test_io_writer_loop_reorders_out_of_order_blocks(#[case] codec: SpillCodec) {
        // Write blocks out of order; io_writer_loop must reassemble them correctly.
        let data1 = b"first block data".to_vec();
        let data2 = b"second block data".to_vec();

        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, codec));
        let (result_tx, result_rx) = pool.compress_result_channel();
        let buffer_pool = pool.buffer_pool.clone();
        let permit_pool = Arc::new(PermitPool::new(4));

        let dir = TempDir::new().unwrap();
        let out_path = dir.path().join("reorder.spill");
        let out_file = std::fs::File::create(&out_path).unwrap();
        let writer = std::io::BufWriter::new(out_file);
        let pp = Arc::clone(&permit_pool);
        let io_handle = std::thread::spawn(move || {
            io_writer_loop(writer, result_rx, buffer_pool, pp, codec, None)
        });

        // Submit block 1 first, then block 0 (out of order).
        // Each needs a pre-acquired permit since they bypass StagingBuffer::flush().
        // Serials are issued in order (0, 1); submit 1 before 0.
        let serial0 = permit_pool.issue_serial();
        let serial1 = permit_pool.issue_serial();
        permit_pool.acquire().unwrap();
        pool.submit_compress(CompressJob {
            data: data2,
            serial: serial1,
            result_tx: result_tx.clone(),
            codec,
            target: CompressTarget::Spill,
        });
        permit_pool.acquire().unwrap();
        pool.submit_compress(CompressJob {
            data: data1,
            serial: serial0,
            result_tx,
            codec,
            target: CompressTarget::Spill,
        });

        // Wait for both compress results to be received by io_writer_loop
        io_handle.join().unwrap().unwrap();

        // Bgzf appends an EOF marker, zstd does not.
        let bytes = std::fs::read(&out_path).unwrap();
        match codec {
            SpillCodec::Bgzf => {
                assert!(bytes.ends_with(&BGZF_EOF), "bgzf output should end with BGZF EOF marker");
            }
            SpillCodec::Zstd => {
                assert!(!bytes.is_empty(), "zstd output should contain the two compressed frames");
                assert!(
                    !bytes.ends_with(&BGZF_EOF),
                    "zstd output must not append the BGZF EOF marker"
                );
            }
        }

        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
    }

    #[rstest]
    #[case(SpillCodec::Bgzf)]
    #[case(SpillCodec::Zstd)]
    fn test_roundtrip_small_data(#[case] codec: SpillCodec) {
        let data = b"hello world from bgzf_io";
        let output = roundtrip_data(data, codec);
        match codec {
            SpillCodec::Bgzf => {
                // Valid BGZF stream ending with the EOF marker, plus a compressed data block.
                assert!(output.ends_with(&BGZF_EOF), "bgzf output must end with BGZF EOF");
                assert!(output.len() > BGZF_EOF.len());
            }
            SpillCodec::Zstd => {
                // Zstd output is `[u32 LE frame-len][zstd frame]`; no EOF marker is appended.
                assert!(!output.is_empty(), "zstd output must contain the compressed frame");
                assert!(!output.ends_with(&BGZF_EOF), "zstd output must not append BGZF EOF");
            }
        }
    }

    #[rstest]
    #[case(SpillCodec::Bgzf)]
    #[case(SpillCodec::Zstd)]
    fn test_roundtrip_empty_data(#[case] codec: SpillCodec) {
        // No data: flush() is a no-op, so the loop sees zero compress results.
        // Bgzf still writes the EOF marker; zstd writes nothing.
        let output = roundtrip_data(b"", codec);
        match codec {
            SpillCodec::Bgzf => {
                assert_eq!(output, BGZF_EOF.to_vec(), "empty bgzf input → only BGZF EOF marker");
            }
            SpillCodec::Zstd => {
                assert!(output.is_empty(), "empty zstd input → empty output (no EOF marker)");
            }
        }
    }

    /// Delivers `serials` to an I/O writer loop over an in-memory buffer whose
    /// permit pool issued `issued` serials, returning the loop's result and the
    /// bytes it wrote.
    fn run_writer_loop(
        issued: u64,
        serials: &[u64],
        codec: SpillCodec,
    ) -> (Result<()>, Vec<u8>, Arc<PermitPool>) {
        let (result_tx, result_rx) = crossbeam_channel::bounded::<CompressResult>(8);
        let permit_pool = Arc::new(PermitPool::new(4));
        for _ in 0..issued {
            permit_pool.issue_serial();
        }
        for &serial in serials {
            result_tx
                .send(CompressResult {
                    serial,
                    compressed: format!("block {serial}").into_bytes(),
                    recycled_buf: Vec::new(),
                })
                .unwrap();
        }
        drop(result_tx);

        let mut writer = BufWriter::new(Vec::new());
        let result = io_writer_loop_inner(
            &mut writer,
            &result_rx,
            &BufferPool::new(4),
            &permit_pool,
            codec,
            None,
        );
        (result, writer.into_inner().unwrap(), permit_pool)
    }

    /// A final block that was submitted but never delivered (its job abandoned
    /// by a pool shutdown, or its worker panicking mid-compress) leaves no gap
    /// in the serials the writer receives. The write must still fail rather
    /// than stamp an EOF marker onto a truncated stream.
    #[rstest]
    #[case(SpillCodec::Bgzf)]
    #[case(SpillCodec::Zstd)]
    fn test_io_writer_loop_fails_when_final_block_never_arrives(#[case] codec: SpillCodec) {
        let (result, bytes, _) = run_writer_loop(2, &[0], codec);
        let err = result.expect_err("a lost final block must fail the write");
        assert!(
            err.to_string().contains("missing compressed block 1 of 2 submitted"),
            "unexpected error: {err}"
        );
        assert!(!bytes.ends_with(&BGZF_EOF), "no EOF marker on a truncated stream");
    }

    /// The public loop closes the permit pool on that failure, so a producer
    /// parked on a permit is released.
    #[test]
    fn test_io_writer_loop_closes_permit_pool_on_lost_final_block() {
        let (result_tx, result_rx) = crossbeam_channel::bounded::<CompressResult>(1);
        let permit_pool = Arc::new(PermitPool::new(4));
        permit_pool.issue_serial();
        drop(result_tx);
        let err = io_writer_loop(
            BufWriter::new(Vec::new()),
            result_rx,
            BufferPool::new(4),
            Arc::clone(&permit_pool),
            SpillCodec::Bgzf,
            None,
        )
        .expect_err("a lost final block must fail the write");
        assert!(err.to_string().contains("missing compressed block 0 of 1"), "{err}");
        assert!(permit_pool.acquire().is_err(), "a failed writer must close the permit pool");
    }

    /// Every issued block delivered: the submitted-count check passes.
    #[test]
    fn test_io_writer_loop_accepts_all_submitted_blocks() {
        let (result, bytes, _) = run_writer_loop(2, &[1, 0], SpillCodec::Bgzf);
        result.expect("all submitted blocks arrived");
        assert_eq!(bytes, [b"block 0".as_slice(), b"block 1", &BGZF_EOF].concat());
    }

    /// A block whose serial the permit pool never issued means a producer
    /// submitted without counting it, which would blind the lost-final-block
    /// check; the writer must refuse it rather than trust the count.
    #[test]
    fn test_io_writer_loop_rejects_serial_not_issued_by_permit_pool() {
        let (result, _, _) = run_writer_loop(1, &[0, 1], SpillCodec::Bgzf);
        let err = result.expect_err("an uncounted block must fail the write");
        assert!(
            err.to_string().contains("compressed block 1 was not issued by the permit pool"),
            "unexpected error: {err}"
        );
        let (result, _, _) = run_writer_loop(0, &[0], SpillCodec::Bgzf);
        result.expect_err("a producer that never counts must fail on its first block");
    }

    /// The staging buffer's serials come from its permit pool, so the pool's
    /// count always matches the blocks submitted.
    #[rstest]
    #[case(SpillCodec::Bgzf)]
    #[case(SpillCodec::Zstd)]
    fn test_staging_buffer_serials_come_from_permit_pool(#[case] codec: SpillCodec) {
        let pool = Arc::new(SortWorkerPool::new(1, 1, 6, codec));
        let (result_tx, result_rx) = pool.compress_result_channel();
        let permit_pool = make_permit_pool(&pool);
        let mut staging = StagingBuffer::new(
            Arc::clone(&pool),
            result_tx,
            Arc::clone(&permit_pool),
            codec,
            CompressTarget::Spill,
        );
        for expected in 0..3 {
            assert_eq!(staging.next_serial(), expected);
            staging.buf().extend_from_slice(b"data");
            staging.flush().unwrap();
        }
        assert_eq!(permit_pool.submitted(), 3);
        drop(staging);
        let mut serials: Vec<u64> = result_rx.iter().map(|r| r.serial).collect();
        serials.sort_unstable();
        assert_eq!(serials, [0, 1, 2]);

        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
    }
}
