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

/// Padding beyond the frame size for the staging buffer capacity.
const STAGING_PADDING: usize = 4096;

/// Uncompressed bytes per zstd **spill** frame.
///
/// BGZF mandates blocks of at most 64 KiB; zstd has no such limit, and spill
/// files are a private temporary format that nothing outside this crate reads.
/// Pinning zstd frames to the BGZF ceiling was therefore vestigial, and it set
/// the merge's block count -- 5,368,249 blocks on a measured 89-way merge, of
/// which the consumer takes roughly one per round trip. Fewer, larger frames cut
/// round trips proportionally.
///
/// Raising this costs memory in two places, both worth watching: the per-worker
/// zstd decompress buffer ([`crate::worker_pool::zstd_decomp_cap`]) and the
/// per-file reorder buffer, which holds up to `PHASE2_DECOMP_CAP` *uncompressed*
/// frames for each of K files. Read-ahead itself does not scale with it, because
/// the refill allowance is sized in bytes rather than blocks.
pub(crate) const SPILL_FRAME_BYTES: usize = BGZF_MAX_BLOCK_SIZE;

/// Bytes a staging buffer accumulates before it submits a compress job.
///
/// Keyed on the codec so BAM output stays inside the BGZF format: output always
/// uses [`SpillCodec::Bgzf`], so it cannot pick up a spill-only frame size by
/// accident.
#[must_use]
pub(crate) const fn spill_frame_bytes(codec: SpillCodec) -> usize {
    match codec {
        SpillCodec::Bgzf => BGZF_MAX_BLOCK_SIZE,
        SpillCodec::Zstd => SPILL_FRAME_BYTES,
    }
}

/// Per-block position notification emitted by the I/O writer loop when index
/// generation is enabled.
///
/// `serial` is the block's serial number (assigned by [`StagingBuffer`] at
/// flush time, and equal to its block number since one compress job produces
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
    next_serial: u64,
    result_tx: Sender<CompressResult>,
    permit_pool: Arc<PermitPool>,
    codec: SpillCodec,
    /// Stamped onto every job this buffer submits, so the worker that pops it
    /// compresses at the level this writer asked for regardless of the pool's
    /// phase at that moment.
    target: CompressTarget,
}

impl StagingBuffer {
    /// The permit pool that carries this writer's histograms.
    ///
    /// Exposed so a caller can retain the [`Arc`] and harvest
    /// [`PermitPool::writer_stats`] *after* the output drain. The writer's own
    /// `finish` consumes the staging to run that drain, so the pool is the only
    /// handle that outlives it. See that method for why the drain must be
    /// included.
    pub(crate) fn permit_pool(&self) -> &Arc<PermitPool> {
        &self.permit_pool
    }

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
            buf: Vec::with_capacity(spill_frame_bytes(codec) + STAGING_PADDING),
            next_serial: 0,
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
        self.buf.len() >= spill_frame_bytes(self.codec)
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
        self.next_serial
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
        let want = spill_frame_bytes(self.codec) + STAGING_PADDING;
        if self.buf.capacity() < want {
            self.buf.reserve(want - self.buf.capacity());
        }

        let serial = self.next_serial;
        self.next_serial += 1;

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
            let space = spill_frame_bytes(self.codec).saturating_sub(self.buf.len());
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
    permit_pool.write_dur.time(|| writer.write_all(data))?;
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
/// Returns an error if any write fails or if a compressed block is missing
/// (which would silently truncate the output).
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
    // Value carries the arrival instant so the wait for an earlier serial is
    // measurable; a block written straight through never enters the map and so
    // records no wait, which is the correct reading.
    let mut reorder_buf: BTreeMap<u64, (Vec<u8>, std::time::Instant)> = BTreeMap::new();
    let mut compressed_offset: u64 = 0;
    let tx = block_offset_tx;

    while let Ok(result) = result_rx.recv() {
        buffer_pool.checkin(result.recycled_buf);

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

            while let Some((data, arrived)) = reorder_buf.remove(&next_expected) {
                permit_pool.write_reorder_wait.record(crate::merge_trace::elapsed_nanos(arrived));
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
            reorder_buf.insert(result.serial, (result.compressed, std::time::Instant::now()));
            // Permit held: released when this block is written in the cascade above.
        }
        // Sampled on every arrival, in-order or not, so the depth reflects the
        // queue the writer is actually carrying rather than only its bad moments.
        // A block count rides a duration histogram, so `record_count` scales it
        // into the microsecond lane; recording the raw count would bucket every
        // depth below `BLOCKS_TO_NANOS` to zero and the report would read zero.
        permit_pool.write_reorder_depth.record_count(reorder_buf.len() as u64);
    }

    // Drain remaining buffered blocks — any gap means a worker dropped a result.
    while let Some((&serial, _)) = reorder_buf.first_key_value() {
        if serial == next_expected {
            let (data, arrived) = reorder_buf.remove(&serial).expect("key just checked");
            permit_pool.write_reorder_wait.record(crate::merge_trace::elapsed_nanos(arrived));
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
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, codec, false));
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
        let pool = Arc::new(SortWorkerPool::new(1, 1, 6, codec, false));
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

        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
    }

    #[rstest]
    #[case(SpillCodec::Bgzf)]
    #[case(SpillCodec::Zstd)]
    fn test_staging_buffer_is_full(#[case] codec: SpillCodec) {
        let pool = Arc::new(SortWorkerPool::new(1, 1, 6, codec, false));
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
        staging.buf().extend(vec![0u8; spill_frame_bytes(codec)]);
        assert!(staging.is_full(), "buffer at the codec's frame size should be full");

        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
    }

    /// BGZF's 64 KiB block ceiling is a *format* requirement; zstd has none.
    ///
    /// Spill files are a private temporary format, so pinning their zstd frames
    /// to the BGZF limit is vestigial — and it sets the block count, which is
    /// what the merge pays per round trip (measured: 5,368,249 blocks on an
    /// 89-way merge, ~1 consumed per consumer round trip). BAM output must stay
    /// at the BGZF size regardless, and it always uses the BGZF codec, so keying
    /// the threshold on the codec keeps output correct by construction.
    #[test]
    fn test_only_zstd_spill_frames_escape_the_bgzf_block_ceiling() {
        assert_eq!(
            spill_frame_bytes(SpillCodec::Bgzf),
            BGZF_MAX_BLOCK_SIZE,
            "BGZF blocks are capped by the format and must not grow"
        );
        assert!(
            spill_frame_bytes(SpillCodec::Zstd) >= BGZF_MAX_BLOCK_SIZE,
            "zstd frames must never be smaller than the BGZF block they replace"
        );
        assert!(
            crate::worker_pool::zstd_decomp_cap() >= spill_frame_bytes(SpillCodec::Zstd),
            "the decompress buffer must hold the largest frame the writer can emit, or every \
             frame at the new size fails to decompress"
        );
        assert!(
            crate::worker_pool::MAX_ZSTD_FRAME_BYTES >= spill_frame_bytes(SpillCodec::Zstd),
            "the read-side length guard must admit a frame the writer can emit"
        );
    }

    /// Both zstd frame caps must admit the largest frame the writer can emit.
    ///
    /// There are two, on two different read paths: the pool's, used by the merge,
    /// and `zspill_stream`'s, used by consolidation. They were previously held
    /// equal by a comment saying "kept in sync", which is not enforcement -- and
    /// the consolidation path is exercised only by the spill-heavy
    /// configurations, so a mismatch would pass the standard matrix and fail
    /// exactly where spill volume is largest.
    #[test]
    fn test_frame_caps_admit_the_largest_frame_the_writer_emits() {
        let frame = spill_frame_bytes(SpillCodec::Zstd);
        assert!(
            crate::worker_pool::zstd_decomp_cap() >= frame,
            "the merge path's decompress cap is below the writer's frame size"
        );
        assert!(
            crate::zspill_stream::frame_decomp_cap() >= frame,
            "the consolidation path's decompress cap is below the writer's frame size"
        );
    }

    /// The spill writer's pre-flush budget must *be* the staging buffer's frame
    /// size, not a copy of it.
    ///
    /// `PooledChunkWriter` pre-flushes so a record never straddles a frame
    /// boundary -- which is what lets the merge borrow most records in place. That
    /// budget was `BGZF_MAX_BLOCK_SIZE` outright, so it flushed at 64 KiB no
    /// matter what the staging buffer was configured for, and raising the frame
    /// size changed the block count *not at all*: 5,368,249 blocks measured at
    /// both 64 KiB and 256 KiB. The sweep looked like a 4% regression rather than
    /// a no-op, so nothing about the wall time revealed that the knob was inert.
    ///
    /// Third duplicated threshold found this way, after the two zstd decompress
    /// caps. Derive, then pin.
    #[rstest]
    #[case(SpillCodec::Bgzf)]
    #[case(SpillCodec::Zstd)]
    fn test_spill_writer_pre_flushes_at_the_staging_frame_size(#[case] codec: SpillCodec) {
        let pool = Arc::new(SortWorkerPool::new(1, 1, 6, codec, false));
        let dir = tempfile::tempdir().expect("tempdir");
        let writer =
            crate::pooled_chunk_writer::PooledChunkWriter::<crate::keys::RawCoordinateKey>::new(
                Arc::clone(&pool),
                &dir.path().join("c.keyed"),
                codec,
            )
            .expect("writer");
        assert_eq!(
            writer.frame_bytes(),
            spill_frame_bytes(codec),
            "the writer's pre-flush budget must equal the frame size the staging buffer              flushes at, or the frame size has no effect on the block count"
        );
        drop(writer);
        if let Ok(p) = Arc::try_unwrap(pool) {
            p.shutdown();
        }
    }

    /// At the default frame size the derived cap must equal the 256 KiB constant
    /// it replaced.
    ///
    /// This buffer is per-worker scratch touched once per decompressed frame, so
    /// its size is a cache parameter. Deriving it with a *fixed* 4 MiB of slack
    /// instead of a proportional 4x inflated it 16x at the default frame size and
    /// cost 25% of merge wall (249.5s against a 199.2s baseline) with peak RSS
    /// essentially unchanged -- a regression invisible to any memory check.
    #[test]
    fn test_default_frame_size_preserves_the_original_decompress_cap() {
        assert_eq!(SPILL_FRAME_BYTES, BGZF_MAX_BLOCK_SIZE, "guard: default frame size");
        assert_eq!(
            crate::worker_pool::zstd_decomp_cap(),
            256 * 1024,
            "the derived cap must reproduce the tuned constant at the default frame size"
        );
    }

    #[rstest]
    #[case(SpillCodec::Bgzf)]
    #[case(SpillCodec::Zstd)]
    fn test_staging_buffer_write_chunked_large_data(#[case] codec: SpillCodec) {
        // Data larger than BGZF_MAX_BLOCK_SIZE must be split into multiple compress jobs.
        let large = vec![b'A'; BGZF_MAX_BLOCK_SIZE * 2 + 1000];
        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, codec, false));
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

        let pool = Arc::new(SortWorkerPool::new(2, 1, 6, codec, false));
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
        permit_pool.acquire().unwrap();
        pool.submit_compress(CompressJob {
            data: data2,
            serial: 1,
            result_tx: result_tx.clone(),
            codec,
            target: CompressTarget::Spill,
        });
        permit_pool.acquire().unwrap();
        pool.submit_compress(CompressJob {
            data: data1,
            serial: 0,
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
}
