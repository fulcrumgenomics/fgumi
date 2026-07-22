#![deny(unsafe_code)]
//! Shared worker pool for all sort pipeline work (N+2 thread model).
//!
//! This module provides `SortWorkerPool`, a fixed-size thread pool where each worker
//! owns reusable compression and decompression state. Workers perform ALL CPU/IO work:
//! block reading, decompression, and compression across all sort phases.
//!
//! # Thread Model (N+2)
//!
//! - **N pool workers**: Do ALL CPU/IO work (block reading, decompression, compression)
//! - **+1 I/O writer thread**: `PooledBamWriter` / `PooledChunkWriter` (sequential disk write)
//! - **+1 main thread**: Phase 1 orchestration + Phase 2 merge loop
//!
//! # Design
//!
//! - **Phase-aware scheduling**: Workers check the current phase to pick eligible steps.
//!   Phase 1: `DecompressInput` > `ReadInputBlocks` > `CompressSpill`.
//!   Phase 2: `Phase2FileWork` (read+decompress the next block of any file that has
//!   room in its reorder buffer) > `CompressOutput`.
//! - **Per-worker state**: Each worker owns an `InlineBgzfCompressor`, a
//!   `libdeflater::Decompressor`, a `zstd::bulk::Compressor`, and a
//!   `zstd::bulk::Decompressor` (the latter pair are used only when the
//!   selected spill codec is `Zstd`). Per-worker contexts avoid cross-thread
//!   synchronization on the compression hot path.
//! - **Per-file work-stealing (Phase 2)**: Each spill file has its own `Phase2FileState`
//!   with a `Mutex`-guarded reader, a raw-block FIFO, and a decompressed reorder buffer.
//!   Workers scan the shared snapshot, pick a file with work to do, and advance it
//!   one block at a time — any worker can make progress on any file.
//! - **Held-item pattern**: Workers never block on queue push. If output is full, they
//!   hold the result locally and advance it first on the next iteration.
//! - **Buffer recycling**: A shared buffer pool (`crossbeam` channel) recycles
//!   `Vec<u8>` buffers to avoid per-job heap allocations.
//! - **Bounded backpressure**: All channels are bounded to prevent unbounded memory
//!   growth when producers outpace consumers.

use crate::codec::{SpillCodec, ZSPILL_MAGIC};
use crossbeam_channel::{Receiver, Sender, bounded};
use crossbeam_queue::ArrayQueue;
use fgumi_bgzf::reader::read_raw_blocks;
use fgumi_bgzf::writer::InlineBgzfCompressor;
use fgumi_bgzf::{RawBgzfBlock, decompress_block};
use log::debug;
use std::collections::VecDeque;
use std::fmt::Write as FmtWrite;
use std::io::{BufReader, Read};
use std::io::{Seek, SeekFrom};
use std::sync::atomic::{AtomicBool, AtomicU8, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::Instant;
use zstd::bulk::{Compressor as ZstdCompressor, Decompressor as ZstdDecompressor};

use fgumi_bam_io::ReorderBuffer;

/// Read up to `n` length-prefixed zstd frames from `reader`.
///
/// On-disk format for a zstd spill file is the four-byte file magic followed
/// by a sequence of `[u32 LE compressed-len][zstd frame bytes]` records. The
/// magic is consumed by `set_phase2_files`, so `reader` is positioned at the
/// first length prefix on entry.
///
/// Returns the frames read. Stops cleanly at a frame boundary on EOF and
/// returns an error if the file is truncated inside a length prefix or frame
/// body, or if a length prefix exceeds `MAX_ZSTD_FRAME_BYTES`.
pub(crate) fn read_raw_zstd_frames<R: std::io::Read + ?Sized>(
    reader: &mut R,
    n: usize,
) -> std::io::Result<Vec<Vec<u8>>> {
    let mut out: Vec<Vec<u8>> = Vec::with_capacity(n);
    for _ in 0..n {
        match read_length_prefix(reader)? {
            None => break,
            Some(frame_len) => {
                let mut frame = vec![0u8; frame_len];
                reader.read_exact(&mut frame)?;
                out.push(frame);
            }
        }
    }
    Ok(out)
}

/// Read a `u32 LE` length prefix and validate it against `MAX_ZSTD_FRAME_BYTES`.
///
/// Returns `Ok(None)` only when the reader is at a clean frame boundary (no
/// bytes consumed before EOF). A 1–3-byte partial prefix surfaces as
/// `UnexpectedEof` so truncation is not silently treated as a clean stop.
pub(crate) fn read_length_prefix<R: std::io::Read + ?Sized>(
    reader: &mut R,
) -> std::io::Result<Option<usize>> {
    use std::io::ErrorKind;
    let mut first = [0u8; 1];
    match reader.read(&mut first)? {
        0 => return Ok(None),
        1 => {}
        n => {
            return Err(std::io::Error::other(format!(
                "Read returned {n} bytes for a 1-byte buffer",
            )));
        }
    }
    let mut rest = [0u8; 3];
    reader.read_exact(&mut rest)?;
    let len_buf = [first[0], rest[0], rest[1], rest[2]];
    let frame_len = u32::from_le_bytes(len_buf) as usize;
    if frame_len > MAX_ZSTD_FRAME_BYTES {
        return Err(std::io::Error::new(
            ErrorKind::InvalidData,
            format!(
                "zstd spill frame length {frame_len} exceeds MAX_ZSTD_FRAME_BYTES ({MAX_ZSTD_FRAME_BYTES}): file likely corrupted",
            ),
        ));
    }
    Ok(Some(frame_len))
}

/// Cap on uncompressed size of a zstd spill frame. Production frames are
/// bounded by the staging buffer (`BGZF_MAX_BLOCK_SIZE` + padding ~= 68 KB);
/// this leaves slack but stays small enough that per-frame allocations don't
/// dominate the merge phase when there are many tens of thousands of frames.
const ZSTD_FRAME_DECOMP_CAP: usize = 256 * 1024;

/// Hard cap on the `u32 LE` length prefix of any zstd spill frame. Frames are
/// produced one per ~64 KiB of input by `handle_compress_job`; even
/// pathological expansion can't reach this. Beyond it, we treat the value as
/// corruption rather than allocate gigabytes.
pub(crate) const MAX_ZSTD_FRAME_BYTES: usize = 2 * 1024 * 1024;

/// Maximum zstd compression level recognized by the `zstd` crate.
const ZSTD_MAX_CLEVEL: u32 = 22;

// ============================================================================
// Job and Result Types
// ============================================================================

/// A compression job: compress uncompressed data into one compressed block.
///
/// The codec determines the output framing:
/// - `SpillCodec::Bgzf` produces one BGZF block (header + deflate + footer).
/// - `SpillCodec::Zstd` produces `[u32 LE frame-len][zstd frame]`.
pub struct CompressJob {
    /// Uncompressed data to compress.
    pub data: Vec<u8>,
    /// Serial number for ordering output blocks.
    pub serial: u64,
    /// Channel to send the compressed result back.
    pub result_tx: Sender<CompressResult>,
    /// Codec to use when compressing.
    pub codec: SpillCodec,
}

/// Result of a compression job.
pub struct CompressResult {
    /// Serial number matching the input job.
    pub serial: u64,
    /// Compressed BGZF block data.
    pub compressed: Vec<u8>,
    /// The original uncompressed buffer, cleared for reuse.
    pub recycled_buf: Vec<u8>,
}

// ============================================================================
// Pool Instrumentation
// ============================================================================

/// Tracks pool activity for diagnostics.
#[derive(Debug, Default)]
pub(crate) struct PoolStats {
    pub(crate) compress_jobs_submitted: AtomicU64,
}

impl PoolStats {
    pub fn log_summary(&self) {
        let compress = self.compress_jobs_submitted.load(Ordering::Relaxed);
        debug!("  Pool stats: {compress} compress jobs");
    }
}

// ============================================================================
// Sort Pipeline Steps
// ============================================================================

/// Sort-specific work steps that pool workers can perform.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum SortStep {
    /// Read raw BGZF blocks from input BAM during Phase 1.
    ReadInputBlocks = 0,
    /// Decompress input BGZF blocks during Phase 1.
    DecompressInput = 1,
    /// Compress data for spill files during Phase 1.
    CompressSpill = 2,
    /// Read+decompress one unit of work for some Phase 2 spill file (work-stealing).
    Phase2FileWork = 3,
    /// Compress data for merge output during Phase 2.
    CompressOutput = 4,
}

impl SortStep {
    /// Number of distinct sort steps.
    pub const COUNT: usize = 5;

    /// Short label for display.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::ReadInputBlocks => "RdInp",
            Self::DecompressInput => "DecInp",
            Self::CompressSpill => "CmpSpl",
            Self::Phase2FileWork => "P2File",
            Self::CompressOutput => "CmpOut",
        }
    }
}

/// Maximum number of worker threads for per-thread stats arrays.
const SORT_MAX_THREADS: usize = 32;

// ============================================================================
// Sort Pipeline Stats
// ============================================================================

/// Comprehensive instrumentation for sort pipeline worker pool.
///
/// Modeled on `PipelineStats` from the unified pipeline but tailored to sort's
/// step set. All fields are `AtomicU64` for lock-free updates from workers.
pub(crate) struct SortPipelineStats {
    // Per-step timing (nanoseconds) — indexed by SortStep as usize
    pub step_ns: [AtomicU64; SortStep::COUNT],
    // Per-step success counts
    pub step_count: [AtomicU64; SortStep::COUNT],

    // Per-thread work distribution: [thread_id][step_index]
    /// How many times each thread completed each step.
    pub per_thread_step_counts: Box<[[AtomicU64; SortStep::COUNT]; SORT_MAX_THREADS]>,
    /// Per-thread idle time in nanoseconds (time in backoff/yield).
    pub per_thread_idle_ns: Box<[AtomicU64; SORT_MAX_THREADS]>,

    /// Number of worker threads (for display bounds).
    pub num_threads: usize,
}

impl SortPipelineStats {
    /// Create a new stats collector for the given number of workers.
    #[must_use]
    pub fn new(num_threads: usize) -> Self {
        Self {
            step_ns: std::array::from_fn(|_| AtomicU64::new(0)),
            step_count: std::array::from_fn(|_| AtomicU64::new(0)),
            per_thread_step_counts: new_sort_2d_array(),
            per_thread_idle_ns: new_sort_1d_array(),
            num_threads,
        }
    }

    /// Record a completed step with timing for a given worker thread.
    pub fn record_step(&self, thread_id: usize, step: SortStep, elapsed_ns: u64) {
        let step_idx = step as usize;
        self.step_ns[step_idx].fetch_add(elapsed_ns, Ordering::Relaxed);
        self.step_count[step_idx].fetch_add(1, Ordering::Relaxed);
        if thread_id < SORT_MAX_THREADS {
            self.per_thread_step_counts[thread_id][step_idx].fetch_add(1, Ordering::Relaxed);
        }
    }

    /// Record idle/backoff time for a worker thread.
    pub fn record_idle(&self, thread_id: usize, idle_ns: u64) {
        if thread_id < SORT_MAX_THREADS {
            self.per_thread_idle_ns[thread_id].fetch_add(idle_ns, Ordering::Relaxed);
        }
    }

    /// Log a comprehensive summary of pipeline statistics.
    #[allow(clippy::cast_precision_loss)]
    pub fn log_summary(&self) {
        let mut s = String::with_capacity(1024);
        writeln!(s, "=== Sort Pipeline Stats ===").expect("write to String");

        // Per-step summary
        let all_steps = [
            SortStep::ReadInputBlocks,
            SortStep::DecompressInput,
            SortStep::CompressSpill,
            SortStep::Phase2FileWork,
            SortStep::CompressOutput,
        ];

        for &step in &all_steps {
            let idx = step as usize;
            let count = self.step_count[idx].load(Ordering::Relaxed);
            if count > 0 {
                let ns = self.step_ns[idx].load(Ordering::Relaxed);
                let secs = ns as f64 / 1_000_000_000.0;
                writeln!(s, "  {:<22} {count:>8} jobs, {secs:>6.1}s total", format!("{step:?}"))
                    .expect("write");
            }
        }

        // Per-thread work distribution
        let nt = self.num_threads.min(SORT_MAX_THREADS);
        if nt > 0 {
            writeln!(s).expect("write");
            writeln!(s, "  Per-Thread Work Distribution:").expect("write");

            // Header
            write!(s, "    Thread").expect("write");
            for &step in &all_steps {
                write!(s, " {:>8}", step.label()).expect("write");
            }
            writeln!(s, "   Idle ms").expect("write");

            // Per-thread rows
            for tid in 0..nt {
                write!(s, "    T{tid:<5}").expect("write");
                for step_idx in 0..SortStep::COUNT {
                    let count = self.per_thread_step_counts[tid][step_idx].load(Ordering::Relaxed);
                    write!(s, " {count:>8}").expect("write");
                }
                let idle_ns = self.per_thread_idle_ns[tid].load(Ordering::Relaxed);
                writeln!(s, " {:>10.1}", idle_ns as f64 / 1_000_000.0).expect("write");
            }

            // Total row
            write!(s, "    Total ").expect("write");
            for step_idx in 0..SortStep::COUNT {
                let mut total = 0u64;
                for tid in 0..nt {
                    total += self.per_thread_step_counts[tid][step_idx].load(Ordering::Relaxed);
                }
                write!(s, " {total:>8}").expect("write");
            }
            let total_idle: u64 =
                (0..nt).map(|tid| self.per_thread_idle_ns[tid].load(Ordering::Relaxed)).sum();
            writeln!(s, " {:>10.1}", total_idle as f64 / 1_000_000.0).expect("write");

            // Utilization
            let total_work_ns: u64 =
                (0..SortStep::COUNT).map(|i| self.step_ns[i].load(Ordering::Relaxed)).sum();
            let work_ms = total_work_ns as f64 / 1_000_000.0;
            let idle_ms = total_idle as f64 / 1_000_000.0;
            let total_ms = work_ms + idle_ms;
            if total_ms > 0.0 {
                let utilization = (work_ms / total_ms) * 100.0;
                writeln!(s).expect("write");
                writeln!(s, "  Thread Utilization: {utilization:.1}% (work={work_ms:.1}ms idle={idle_ms:.1}ms)")
                    .expect("write");
            }
        }

        // Log as a single multiline message
        for line in s.trim_end().lines() {
            debug!("{line}");
        }
    }
}

/// Helper to create a boxed 2D array of `AtomicU64` for sort stats.
#[allow(clippy::unnecessary_box_returns)]
fn new_sort_2d_array() -> Box<[[AtomicU64; SortStep::COUNT]; SORT_MAX_THREADS]> {
    let v: Vec<[AtomicU64; SortStep::COUNT]> =
        (0..SORT_MAX_THREADS).map(|_| std::array::from_fn(|_| AtomicU64::new(0))).collect();
    v.into_boxed_slice().try_into().expect("Vec length matches SORT_MAX_THREADS")
}

/// Helper to create a boxed 1D array of `AtomicU64` for sort stats.
#[allow(clippy::unnecessary_box_returns)]
fn new_sort_1d_array() -> Box<[AtomicU64; SORT_MAX_THREADS]> {
    let v: Vec<AtomicU64> = (0..SORT_MAX_THREADS).map(|_| AtomicU64::new(0)).collect();
    v.into_boxed_slice().try_into().expect("Vec length matches SORT_MAX_THREADS")
}

// ============================================================================
// Buffer Pool
// ============================================================================

/// Recycling pool for `Vec<u8>` buffers.
///
/// Uses a bounded crossbeam channel: producers return used buffers,
/// consumers check out buffers (falling back to new allocation if empty).
pub struct BufferPool {
    tx: Sender<Vec<u8>>,
    rx: Receiver<Vec<u8>>,
}

impl BufferPool {
    /// Create a buffer pool with the given capacity.
    #[must_use]
    pub fn new(capacity: usize) -> Self {
        let (tx, rx) = bounded(capacity);
        Self { tx, rx }
    }

    /// Get a buffer from the pool, or allocate a new one if the pool is empty.
    #[must_use]
    pub fn checkout(&self) -> Vec<u8> {
        self.rx.try_recv().unwrap_or_default()
    }

    /// Returns the number of buffers currently available in the pool.
    #[must_use]
    pub fn len(&self) -> usize {
        self.rx.len()
    }

    /// Returns true if no buffers are currently available in the pool.
    #[must_use]
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.rx.is_empty()
    }

    /// Return a buffer to the pool for reuse.
    /// If the pool is full, the buffer is dropped.
    pub fn checkin(&self, mut buf: Vec<u8>) {
        buf.clear();
        let _ = self.tx.try_send(buf);
    }
}

impl Clone for BufferPool {
    fn clone(&self) -> Self {
        Self { tx: self.tx.clone(), rx: self.rx.clone() }
    }
}

// ============================================================================
// Permit Pool
// ============================================================================

/// Bounded semaphore for controlling in-flight compressed blocks.
///
/// Pre-filled with `capacity` permits. `StagingBuffer::flush()` acquires a
/// permit (blocking) before submitting each compress job; `io_writer_loop`
/// releases a permit after each `write_all`. At most `capacity` compressed
/// blocks exist anywhere in the pipeline simultaneously, bounding the reorder
/// buffer to `capacity × BGZF_MAX_BLOCK_SIZE` bytes.
pub(crate) struct PermitPool {
    tx: std::sync::Mutex<Option<Sender<()>>>,
    rx: Receiver<()>,
    /// Set by [`PermitPool::close`] so `acquire` fails immediately instead of
    /// draining permits that are still buffered in the channel.
    ///
    /// Dropping the `Sender` alone is not enough: `recv` keeps succeeding until
    /// the buffered permits run out, so after the I/O writer errors, producers
    /// would each take one more permit and compress another block -- up to
    /// `num_workers * 4` blocks of work thrown away -- before the disconnect
    /// finally surfaced. This never produced wrong output or a hang; it just
    /// wasted work on the failure path.
    closed: AtomicBool,
}

impl PermitPool {
    /// Create a permit pool pre-filled with `capacity` permits.
    pub(crate) fn new(capacity: usize) -> Self {
        let (tx, rx) = bounded(capacity);
        for _ in 0..capacity {
            tx.try_send(()).expect("fresh channel has capacity for initial permits");
        }
        Self { tx: std::sync::Mutex::new(Some(tx)), rx, closed: AtomicBool::new(false) }
    }

    /// Acquire a permit, blocking until one is available.
    ///
    /// Returns an error if the pool has been closed (I/O writer exited with an
    /// error). The closed flag is checked both before *and* after taking a
    /// permit: the pre-check fails fast when the pool is already closed, and the
    /// post-`recv` re-check closes the race where [`close`](Self::close) sets the
    /// flag (and drops the sender) while this producer is parked, yet `recv`
    /// still hands back a permit that was buffered before the close. Without the
    /// re-check that producer would compress another block no writer will ever
    /// service. The permit received in that race is simply dropped — a closed
    /// pool never returns permits to the channel.
    pub(crate) fn acquire(&self) -> anyhow::Result<()> {
        if self.closed.load(Ordering::Acquire) {
            anyhow::bail!("permit pool closed: I/O writer thread exited");
        }
        self.rx
            .recv()
            .map_err(|_| anyhow::anyhow!("permit pool closed: I/O writer thread exited"))?;
        if self.closed.load(Ordering::Acquire) {
            anyhow::bail!("permit pool closed: I/O writer thread exited");
        }
        Ok(())
    }

    /// Release a permit back to the pool after a block has been written to disk.
    ///
    /// A release after `close` is a no-op: the sender is already gone.
    pub(crate) fn release(&self) {
        if let Ok(guard) = self.tx.lock()
            && let Some(tx) = guard.as_ref()
        {
            let _ = tx.try_send(());
        }
    }

    /// Close the pool, unblocking any threads waiting on `acquire()`.
    ///
    /// Drops the sending half of the channel so that `rx.recv()` in `acquire()`
    /// returns `Err`, which is mapped to an `anyhow` error. Called by
    /// `io_writer_loop` on write error to prevent producers from parking forever.
    pub(crate) fn close(&self) {
        // Set the flag before dropping the sender so any producer that observes
        // the disconnect also observes the closed state.
        self.closed.store(true, Ordering::Release);
        if let Ok(mut guard) = self.tx.lock() {
            guard.take(); // drops the Sender, closing the channel
        }
    }
}

/// Number of raw BGZF blocks to read in each I/O batch from input file.
const INPUT_READ_BATCH_SIZE: usize = 16;

// ============================================================================
// Phase 2 per-file state (work-stealing across files)
// ============================================================================

/// Number of raw BGZF blocks to keep read-ahead per spill file.
///
/// Bounds disk read-ahead memory: K files × `PHASE2_RAW_CAP` × ~64 KB.
pub(crate) const PHASE2_RAW_CAP: usize = 8;

/// Number of decompressed blocks the per-file reorder buffer may hold before
/// workers stop decompressing more for that file.
///
/// Bounds in-flight decompressed memory: K files × `PHASE2_DECOMP_CAP` ×
/// payload size. Payload size is ~64 KB for BGZF (one block) and up to
/// `ZSTD_FRAME_DECOMP_CAP` (256 KB) for zstd frames. This is a soft cap — the
/// "always accept the next-expected serial" rule lets it transiently exceed by
/// up to ~`num_workers` blocks per file.
pub(crate) const PHASE2_DECOMP_CAP: usize = 8;

/// Number of raw blocks to read from disk per `ReadRawBlocks` call.
pub(crate) const PHASE2_READ_BATCH: usize = 4;

/// Reader state for a single spill file. Locked when reading from disk.
pub(crate) struct Phase2Reader {
    pub(crate) inner: BufReader<std::fs::File>,
    pub(crate) next_serial: u64,
    pub(crate) eof: bool,
}

/// Per-spill-file state shared between all pool workers and the main thread.
///
/// Phase 2 uses work stealing across files: any worker can grab work from any
/// file. The locks here are deliberately fine-grained so different workers can
/// be reading, decompressing, and the main thread can be popping records all
/// concurrently as long as they touch different sub-states.
pub(crate) struct Phase2FileState {
    /// Disk reader. Held only while popping bytes from disk.
    pub(crate) reader: Mutex<Phase2Reader>,
    /// Lock-free copy of `reader.eof`. Set immediately after any code path
    /// that sets `reader_guard.eof = true`, so `is_drained` can fast-path
    /// without touching the reader mutex (called per decompressed block on
    /// the hot Phase 2 path).
    pub(crate) reader_eof: AtomicBool,
    /// Codec used to compress this file. Detected at open time from magic.
    pub(crate) codec: SpillCodec,
    /// Raw compressed blocks read from disk, in serial order. For BGZF, each
    /// entry is the raw block bytes (header + deflate + footer). For zstd,
    /// each entry is one complete zstd frame's bytes.
    pub(crate) raw_blocks: Mutex<VecDeque<(u64, Vec<u8>)>>,
    /// Decompressed blocks reordered by serial. Main thread pops the next-in-order
    /// block here when its parser exhausts the current buffer.
    pub(crate) decompressed: Mutex<ReorderBuffer<Vec<u8>>>,
    /// Number of raw blocks currently being decompressed (popped from
    /// `raw_blocks` but not yet inserted into `decompressed`). Used by
    /// `is_drained` to avoid a race where the consumer exits while a worker
    /// is mid-decompress.
    pub(crate) decomp_in_flight: AtomicUsize,
}

impl Phase2FileState {
    pub(crate) fn new(reader: BufReader<std::fs::File>, codec: SpillCodec) -> Self {
        Self {
            reader: Mutex::new(Phase2Reader { inner: reader, next_serial: 0, eof: false }),
            reader_eof: AtomicBool::new(false),
            codec,
            raw_blocks: Mutex::new(VecDeque::with_capacity(PHASE2_RAW_CAP)),
            decompressed: Mutex::new(ReorderBuffer::new()),
            decomp_in_flight: AtomicUsize::new(0),
        }
    }

    /// Mark the disk reader as having reached EOF. Updates both the
    /// reader-internal flag and the lock-free atomic copy. Must be called
    /// while holding the reader `Mutex` (pass the guard to prove it).
    pub(crate) fn mark_reader_eof(&self, reader_guard: &mut Phase2Reader) {
        reader_guard.eof = true;
        self.reader_eof.store(true, Ordering::Release);
    }

    /// Gather probe statistics for this file: `(pending_blocks, pending_bytes, active)`.
    ///
    /// `pending_blocks` counts raw + decompressed blocks in flight.
    /// `pending_bytes` sums the byte length of decompressed blocks.
    /// `active` is true if the disk reader has not yet reached EOF.
    pub(crate) fn probe_stats(&self) -> (u64, u64, bool) {
        let raw_len =
            self.raw_blocks.lock().expect("phase2 raw_blocks mutex poisoned").len() as u64;
        let decomp_guard = self.decompressed.lock().expect("phase2 decompressed mutex poisoned");
        let decomp_len = decomp_guard.len() as u64;
        let decomp_bytes: u64 = decomp_guard.iter().map(|buf| buf.len() as u64).sum();
        drop(decomp_guard);
        let active = !self.reader_eof.load(Ordering::Relaxed);
        (raw_len + decomp_len, decomp_bytes, active)
    }

    /// Returns true when this file has produced all its data: disk reader at
    /// EOF, no raw blocks waiting, no decompressed blocks waiting, and no
    /// decompression in progress.
    ///
    /// Fast path: if the disk reader has not yet reached EOF, returns `false`
    /// without acquiring any mutex. This is the overwhelmingly common case on
    /// the Phase 2 hot path — every successful decompression calls this to
    /// decide whether to wake the consumer.
    pub(crate) fn is_drained(&self) -> bool {
        if !self.reader_eof.load(Ordering::Acquire) {
            return false;
        }
        let raw_empty =
            self.raw_blocks.lock().expect("phase2 raw_blocks mutex poisoned").is_empty();
        if !raw_empty {
            return false;
        }
        if self.decomp_in_flight.load(Ordering::Acquire) > 0 {
            return false;
        }
        self.decompressed.lock().expect("phase2 decompressed mutex poisoned").is_empty()
    }
}

/// Phase constants for the sort pipeline.
pub mod phase {
    /// Pool is shut down.
    pub const SHUTDOWN: u8 = 0;
    /// Phase 1: Reading input, sorting, spilling.
    pub const PHASE1: u8 = 1;
    /// Phase 2: Merge reading from spill files + compressing output.
    pub const PHASE2: u8 = 2;
    /// Legacy mode: channel-based compress/decompress only (no phase-aware scheduling).
    /// This is the default state when the pool is first created.
    pub const LEGACY: u8 = 255;
}

// ============================================================================
// SortWorkerPool
// ============================================================================

/// Shared pool of N worker threads for ALL sort pipeline work.
///
/// Workers perform all CPU/IO work: block reading, decompression, and compression.
/// The pool is created once per sort invocation and reused across all phases.
///
/// # Thread Model (N+2)
///
/// - **N pool workers**: Do ALL CPU/IO work (block reading, decompression, compression)
/// - **+1 I/O writer thread**: `PooledBamWriter` / `PooledChunkWriter` (sequential disk write)
/// - **+1 main thread**: Phase 1 orchestration + Phase 2 merge loop
///
/// # Phase-Aware Scheduling
///
/// Workers check the current phase to determine eligible steps:
/// - **Phase 1**: `DecompressInput` > `ReadInputBlocks` > `CompressSpill`
/// - **Phase 2**: `Phase2FileWork` (per-file work stealing: read next raw block
///   batch from any file with FIFO room, OR decompress the next queued raw
///   block for any file whose reorder buffer has capacity) > `CompressOutput`
pub struct SortWorkerPool {
    // Shared pipeline state (visible to workers and main thread)
    shared: Arc<SharedPipelineState>,

    /// Worker thread handles. `None` after shutdown (taken by `Drop` or `shutdown`).
    workers: Option<Vec<JoinHandle<()>>>,
    pub(crate) stats: PoolStats,
    pub(crate) pipeline_stats: Arc<SortPipelineStats>,
    pub buffer_pool: BufferPool,
    num_workers: usize,
    pub(crate) spill_codec: SpillCodec,
}

/// Shared state visible to all workers and the main thread.
///
/// Uses `ArrayQueue` for all inter-step data queues (lock-free, non-blocking
/// `push()`/`pop()` only). The compress result channel stays as `crossbeam_channel`
/// because the I/O writer thread needs blocking `recv()`.
pub(crate) struct SharedPipelineState {
    /// Current phase: 0=shutdown, 1=Phase1, 2=Phase2, 255=Legacy.
    pub(crate) phase: AtomicU8,

    /// Number of workers permitted to be active in the current phase. Workers
    /// with `worker_id >= active_worker_limit` idle (backoff) instead of taking
    /// work. Lets the sort driver run Phase 1 (accumulate/sort/spill) on fewer
    /// threads than Phase 2 (merge/write) — e.g. to cede cores to an upstream
    /// producer in a pipeline. Defaults to `num_workers` (all active), so the
    /// behavior is unchanged unless the driver calls `set_active_workers`.
    /// Capped workers re-check this within `MAX_BACKOFF_US`, so raising the
    /// limit reactivates them without an explicit wake.
    pub(crate) active_worker_limit: AtomicUsize,

    // --- Phase 1 queues ---
    /// Input BAM file (Mutex because only one worker reads at a time).
    pub(crate) input_file: std::sync::Mutex<Option<Box<dyn Read + Send>>>,
    /// Input EOF flag — set by the worker that reads the last block.
    pub(crate) input_eof: AtomicBool,
    /// Set when an I/O error occurs reading the input file.
    ///
    /// `PooledInputStream::read()` checks this and returns `io::Error` so the
    /// error surfaces to the caller rather than appearing as a truncated stream.
    pub(crate) input_read_error: Arc<AtomicBool>,
    /// Set when a worker fails to decompress a BGZF block (Phase 1 or Phase 2).
    ///
    /// Workers set this and wake the main thread rather than panicking, so the
    /// main thread can surface the error instead of parking forever.
    pub(crate) decompression_error: Arc<AtomicBool>,
    /// Set when a worker encounters an I/O error reading a chunk file (Phase 2).
    ///
    /// Without this flag, a chunk read error silently marks the source as EOF and
    /// the merge loop produces a truncated output BAM.  The main thread checks this
    /// in `poll_decompressed_blocks` and surfaces the error instead.
    pub(crate) chunk_read_error: Arc<AtomicBool>,
    /// Set when a worker thread panics unexpectedly.
    ///
    /// `do_shutdown` checks join results and sets this flag so the main thread
    /// does not park forever waiting for work that will never arrive.
    pub(crate) worker_panicked: Arc<AtomicBool>,
    /// Next serial for input block reading (atomic increment for ordering).
    input_read_serial: AtomicU64,
    /// Raw input blocks: `ReadInputBlocks` → `DecompressInput`.
    pub(crate) raw_input_blocks: Arc<ArrayQueue<(u64, RawBgzfBlock)>>,
    /// Decompressed input blocks: `DecompressInput` → main thread.
    pub(crate) decompressed_input: Arc<ArrayQueue<(u64, Vec<u8>)>>,
    /// Set by the last worker to decompress the final input block.
    pub(crate) decompressed_input_done: Arc<AtomicBool>,
    /// Count of input blocks successfully pushed to `decompressed_input`.
    ///
    /// Unlike `input_blocks_decompressed` (which counted blocks entering decompression),
    /// this only increments when a block reaches the queue. `decompressed_input_done`
    /// is set only when this equals `input_read_serial`, preventing the race where a
    /// worker with a held (not-yet-queued) block causes premature EOF signalling.
    ///
    /// This comparison is only sound because `try_read_input_blocks` reserves a
    /// batch's whole serial range from `input_read_serial` *while still holding the
    /// input-file lock*. Reserving after releasing the lock would let a worker's
    /// already-read blocks go unaccounted for while a sibling reads the EOF batch and
    /// sets `input_eof`, so the two counters could compare equal with reads still in
    /// flight and the merge would finalize early.
    input_blocks_queued: AtomicU64,

    // --- Phase 2 per-file state (work-stealing) ---
    /// Per-spill-file state, indexed by merge source id (`0..K`). Built before
    /// `set_phase(PHASE2)` and cleared when phase 2 ends.
    ///
    /// `RwLock` so the main thread can swap the vector across phase 2
    /// invocations. Workers and the main thread only hold the read guard long
    /// enough to `Arc::clone` the inner `Vec` (see `phase2_files_snapshot`),
    /// so writers never block on a long-held reader.
    pub(crate) phase2_files: std::sync::RwLock<Arc<Vec<Phase2FileState>>>,
    /// All chunk files have reached disk EOF (raw queues may still hold blocks).
    pub(crate) all_chunks_eof: Arc<AtomicBool>,
    /// Number of files whose disk reader has hit EOF (when == K, set `all_chunks_eof`).
    sources_at_eof: AtomicU64,
    /// Total number of chunk sources (set before Phase 2 starts).
    pub(crate) total_sources: AtomicU64,

    // --- Compress queue (shared Phase 1 + Phase 2) ---
    /// Compress jobs: main thread → workers (`ArrayQueue`, non-blocking push).
    pub(crate) compress_queue: Arc<ArrayQueue<CompressJob>>,

    /// Number of workers (for `low_water` threshold in backpressure).
    num_workers: usize,

    /// Main thread handle for `park()`/`unpark()` notification.
    /// Workers call `unpark()` after inserting into `decompressed_input`
    /// (Phase 1) or into a per-file `Phase2FileState.decompressed` reorder
    /// buffer (Phase 2) so the main thread wakes immediately instead of
    /// spin-yielding.
    main_thread_handle: std::thread::Thread,
}

impl SharedPipelineState {
    fn new(num_workers: usize, main_thread_handle: std::thread::Thread) -> Self {
        let data_queue_cap = num_workers * 8;
        let compress_queue_cap = num_workers * 4;

        Self {
            phase: AtomicU8::new(phase::LEGACY),
            active_worker_limit: AtomicUsize::new(num_workers),

            input_file: std::sync::Mutex::new(None),
            input_eof: AtomicBool::new(false),
            input_read_error: Arc::new(AtomicBool::new(false)),
            decompression_error: Arc::new(AtomicBool::new(false)),
            chunk_read_error: Arc::new(AtomicBool::new(false)),
            worker_panicked: Arc::new(AtomicBool::new(false)),
            input_read_serial: AtomicU64::new(0),
            raw_input_blocks: Arc::new(ArrayQueue::new(data_queue_cap)),
            decompressed_input: Arc::new(ArrayQueue::new(data_queue_cap)),
            decompressed_input_done: Arc::new(AtomicBool::new(false)),
            input_blocks_queued: AtomicU64::new(0),

            phase2_files: std::sync::RwLock::new(Arc::new(Vec::new())),
            all_chunks_eof: Arc::new(AtomicBool::new(false)),
            sources_at_eof: AtomicU64::new(0),
            total_sources: AtomicU64::new(0),

            compress_queue: Arc::new(ArrayQueue::new(compress_queue_cap)),

            num_workers,
            main_thread_handle,
        }
    }

    /// Snapshot the current Phase 2 file vector. Cheap (just clones the `Arc`).
    pub(crate) fn phase2_files_snapshot(&self) -> Arc<Vec<Phase2FileState>> {
        Arc::clone(&self.phase2_files.read().expect("phase2_files rwlock poisoned"))
    }

    /// Snapshot current queue depths for backpressure-driven scheduling.
    fn get_backpressure(&self) -> SortBackpressureState {
        let current_phase = self.phase.load(Ordering::Acquire);
        let low_water = self.num_workers;

        SortBackpressureState {
            decompressed_input_low: self.decompressed_input.len() < low_water,
            input_eof: self.input_eof.load(Ordering::Acquire),
            decompressed_input_done: self.decompressed_input_done.load(Ordering::Acquire),

            compress_has_items: !self.compress_queue.is_empty(),
            phase: current_phase,
        }
    }
}

/// Per-worker mutable state — no sharing, no locks.
///
/// Every step output has a held-item slot. If an `ArrayQueue::push()` fails
/// because the downstream queue is full, the item is stored here. On the next
/// loop iteration, `try_advance_all_held` drains held items before any new work
/// is attempted. This prevents deadlock: a worker holding an output item will
/// try other steps (e.g., compress) instead of blocking.
struct SortWorkerState {
    worker_id: usize,
    /// Compressor used for Phase 1 spill writes (temp compression level).
    compressor: InlineBgzfCompressor,
    /// Compressor used for Phase 2 merge output (output compression level).
    output_compressor: InlineBgzfCompressor,
    /// Zstd compressor reused across spill frames when `SpillCodec::Zstd`.
    zstd_compressor: ZstdCompressor<'static>,
    /// Zstd decompressor reused across Phase 2 frames when `SpillCodec::Zstd`.
    zstd_decompressor: ZstdDecompressor<'static>,
    /// Scratch buffer reused across zstd frame decompressions to avoid
    /// allocating a fresh Vec for every frame on the merge hot path.
    zstd_decompress_buf: Vec<u8>,
    decompressor: libdeflater::Decompressor,
    /// Phase 2 file scan cursor — starts at `worker_id` and advances on success
    /// for cache locality and reduced lock contention. Workers no longer own a
    /// fixed subset of files; any worker can do work on any file.
    phase2_file_cursor: usize,

    // Held items (one per step output) — see plan §Worker State
    /// Held raw input blocks from `ReadInputBlocks` (couldn't push to `raw_input_blocks` queue).
    held_raw_input_blocks: Vec<(u64, RawBgzfBlock)>,
    /// Held decompressed input block from `DecompressInput`.
    held_decompressed_input: Option<(u64, Vec<u8>)>,
    // Compress output goes directly to result_tx channel (I/O thread) — no held item.
    /// Backoff microseconds for idle spinning.
    backoff_us: u64,
    /// Monotonic counter incremented on each idle sleep; mixed with `worker_id` to
    /// produce per-worker jitter so all workers don't wake simultaneously.
    idle_iter: u64,
}

impl SortWorkerState {
    /// Returns true if this worker is holding any items that need advancement.
    /// CRITICAL: Workers must not exit while holding items — they would be lost.
    fn has_any_held_items(&self) -> bool {
        !self.held_raw_input_blocks.is_empty() || self.held_decompressed_input.is_some()
    }
}

// ============================================================================
// Backpressure State and Priority Selection
// ============================================================================

/// Snapshot of queue depths for backpressure-driven scheduling.
///
/// Sampled once per worker loop iteration via `SharedPipelineState::get_backpressure()`.
/// All checks use `ArrayQueue::len()` which is O(1) and lock-free.
#[allow(clippy::struct_excessive_bools)]
struct SortBackpressureState {
    // Phase 1
    decompressed_input_low: bool,
    input_eof: bool,
    decompressed_input_done: bool,

    // Shared
    compress_has_items: bool,
    phase: u8,
}

/// Backpressure-driven priority selection — the sort pipeline's equivalent
/// of `BalancedChaseDrain.build_priorities()`.
///
/// Returns a static slice of steps ordered by priority. The scheduler naturally
/// adapts to all 7 sub-phases without explicit phase tracking because:
/// - During 1A (reading): compress queue empty, decompressed low → read/decompress
/// - During 1C (sort): decompressed full → skip decompress, do compress if available
/// - During 1D (spill): compress queue fills → prioritize compress
/// - During 1E (overlap): both compress and decompress needed → split by queue depths
/// - During Phase 2: both compress (output) and decompress (chunks) → split
fn get_sort_priorities(bp: &SortBackpressureState) -> &'static [SortStep] {
    match bp.phase {
        phase::PHASE1 => {
            if bp.input_eof && !bp.compress_has_items && bp.decompressed_input_done {
                // Input fully decompressed and no compress work — nothing productive to do
                &[]
            } else if bp.compress_has_items && !bp.decompressed_input_low {
                // Spill compression is the bottleneck (13.7s at t4). Drain compress
                // while decompressed blocks are plentiful for the main thread.
                &[SortStep::CompressSpill, SortStep::DecompressInput, SortStep::ReadInputBlocks]
            } else {
                // Default/starving: feed the main thread first, compress if available
                &[SortStep::DecompressInput, SortStep::ReadInputBlocks, SortStep::CompressSpill]
            }
        }
        phase::PHASE2 => {
            // Phase 2 file work and output compression are independent: each
            // worker grabs whatever has work. We never gate on `all_chunks_eof`
            // — even after disk reads finish, decompression and parser drain
            // continue until all per-file reorder buffers empty.
            if bp.compress_has_items {
                // Drain output compression while we can; it's the writer-side bottleneck.
                &[SortStep::CompressOutput, SortStep::Phase2FileWork]
            } else {
                &[SortStep::Phase2FileWork, SortStep::CompressOutput]
            }
        }
        // Legacy/transition: compress only (drain any remaining jobs)
        _ => &[SortStep::CompressSpill],
    }
}

// ============================================================================
// Backoff (ported from unified pipeline base.rs)
// ============================================================================

/// Minimum backoff duration in microseconds.
const MIN_BACKOFF_US: u64 = 10;
/// Maximum backoff duration in microseconds (1ms).
const MAX_BACKOFF_US: u64 = 1000;

/// Sleep for the given backoff duration with ±25% jitter.
///
/// Uses `yield_now()` at minimum backoff to avoid sleep syscall overhead.
/// `worker_id` and `iter` are mixed into the seed so concurrent workers
/// do not produce identical sleep durations and thundering-herd on wakeup.
fn sleep_with_jitter(backoff_us: u64, worker_id: usize, iter: u64) {
    if backoff_us <= MIN_BACKOFF_US {
        std::thread::yield_now();
    } else {
        let jitter_range = backoff_us / 4;
        // Cheap deterministic seed — no syscall, differs per worker and iteration.
        let jitter_seed = (worker_id as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15).wrapping_add(iter);
        let jitter_offset = (jitter_seed % (jitter_range * 2)).saturating_sub(jitter_range);
        let actual_us = backoff_us.saturating_add(jitter_offset).max(MIN_BACKOFF_US);
        std::thread::sleep(std::time::Duration::from_micros(actual_us));
    }
}

// ============================================================================
// Held-Item Helpers (ported from unified pipeline)
// ============================================================================

/// Try to advance a held item to its output `ArrayQueue`.
/// Returns true if the held slot is now empty (item pushed or was already None).
fn try_advance_held<T>(queue: &ArrayQueue<T>, held: &mut Option<T>) -> bool {
    if let Some(item) = held.take() {
        match queue.push(item) {
            Ok(()) => true,
            Err(item) => {
                *held = Some(item);
                false
            }
        }
    } else {
        true // nothing held
    }
}

/// Try to advance a batch of held items to an `ArrayQueue`.
/// Returns true if all items were pushed (held vec is now empty).
fn try_advance_held_batch<T>(queue: &ArrayQueue<T>, held: &mut Vec<T>) -> bool {
    if held.is_empty() {
        return true;
    }
    let batch = std::mem::take(held);
    let mut iter = batch.into_iter();
    for item in iter.by_ref() {
        match queue.push(item) {
            Ok(()) => {}
            Err(item) => {
                held.push(item);
                break;
            }
        }
    }
    held.extend(iter); // keep remaining
    held.is_empty()
}

/// Result of attempting a work step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum StepResult {
    /// Step completed successfully.
    Success,
    /// Step has input but output queue is full.
    OutputFull,
    /// Step has no input available.
    InputEmpty,
}

/// Publishes a worker panic to `worker_panicked` and wakes the main thread if
/// the worker loop unwinds.
///
/// `SortWorkerPool::do_shutdown` also sets this flag when a join reports a
/// panic, but that only happens after Phase 2 has finished. A worker that
/// panics mid-run -- especially while holding a reader mutex, which leaves
/// `eof` unset and every later `try_lock` returning `Poisoned` -- would
/// otherwise leave the merge loop spinning forever. Publishing from the worker
/// itself turns that silent hang into a surfaced error.
///
/// Call [`Self::disarm`] on the normal exit path so a clean shutdown does not
/// report a panic.
struct WorkerPanicGuard {
    shared: Arc<SharedPipelineState>,
    /// `true` until [`disarm`](Self::disarm) runs. `Drop` reports a panic only
    /// while armed, so the clean exit path can defuse the guard without leaking.
    armed: bool,
}

impl WorkerPanicGuard {
    /// Arm a guard that reports a worker panic on drop unless disarmed first.
    fn new(shared: Arc<SharedPipelineState>) -> Self {
        Self { shared, armed: true }
    }

    /// Defuse the guard after the worker loop returns normally.
    ///
    /// This clears the armed flag and lets the guard drop normally, releasing
    /// its `Arc<SharedPipelineState>`. Using `mem::forget` here would skip the
    /// destructor and leak that `Arc`, keeping the shared pipeline state (and
    /// its queued buffers / phase-2 state) alive for the whole process after
    /// every clean worker exit.
    fn disarm(mut self) {
        self.armed = false;
    }
}

impl Drop for WorkerPanicGuard {
    fn drop(&mut self) {
        if !self.armed {
            // Clean exit: `disarm` already ran, so there is nothing to report.
            return;
        }
        // Reached only when the worker loop unwound.
        self.shared.worker_panicked.store(true, Ordering::Release);
        self.shared.main_thread_handle.unpark();
    }
}

impl SortWorkerPool {
    /// Create a new worker pool with `num_workers` threads.
    ///
    /// Each worker owns its own compressor (×2: spill and output) and decompressor.
    /// Workers are phase-aware and perform all CPU/IO work.
    ///
    /// - `temp_compression`: BGZF level for Phase 1 spill writes (typically 1 for speed).
    /// - `output_compression`: BGZF level for Phase 2 merge output (typically 6 for size).
    /// - `spill_codec`: codec used for spill chunks (BGZF or Zstd). Output is always BGZF.
    #[must_use]
    pub fn new(
        num_workers: usize,
        temp_compression: u32,
        output_compression: u32,
        spill_codec: SpillCodec,
    ) -> Self {
        let buffer_pool = BufferPool::new(num_workers * 4);
        let stats = PoolStats::default();
        let pipeline_stats = Arc::new(SortPipelineStats::new(num_workers));
        let main_thread_handle = std::thread::current();
        let shared = Arc::new(SharedPipelineState::new(num_workers, main_thread_handle));

        let workers: Vec<JoinHandle<()>> = (0..num_workers)
            .map(|worker_id| {
                let pstats = Arc::clone(&pipeline_stats);
                let shared = Arc::clone(&shared);

                thread::spawn(move || {
                    let zstd_level =
                        i32::try_from(temp_compression.clamp(1, ZSTD_MAX_CLEVEL)).expect("clamped");
                    let zstd_decompress_buf = if matches!(spill_codec, SpillCodec::Zstd) {
                        vec![0u8; ZSTD_FRAME_DECOMP_CAP]
                    } else {
                        Vec::new()
                    };
                    let mut worker = SortWorkerState {
                        worker_id,
                        compressor: InlineBgzfCompressor::new(temp_compression),
                        output_compressor: InlineBgzfCompressor::new(output_compression),
                        zstd_compressor: ZstdCompressor::new(zstd_level)
                            .expect("zstd compressor init"),
                        zstd_decompressor: ZstdDecompressor::new().expect("zstd decompressor init"),
                        zstd_decompress_buf,
                        decompressor: libdeflater::Decompressor::new(),
                        phase2_file_cursor: worker_id,
                        held_raw_input_blocks: Vec::new(),
                        held_decompressed_input: None,
                        backoff_us: MIN_BACKOFF_US,
                        idle_iter: 0,
                    };

                    // Publish a panic as soon as the worker unwinds, not at
                    // join time. `do_shutdown` also sets this flag, but only
                    // after the merge has finished; a worker that panics
                    // mid-run while holding a reader mutex would otherwise
                    // leave `eof` unset and the merge loop spinning, turning a
                    // diagnosable panic into a silent 100%-CPU hang.
                    let panic_guard = WorkerPanicGuard::new(Arc::clone(&shared));
                    Self::worker_loop(&shared, &mut worker, &pstats);
                    panic_guard.disarm();
                })
            })
            .collect();

        Self {
            shared,
            workers: Some(workers),
            stats,
            pipeline_stats,
            buffer_pool,
            num_workers,
            spill_codec,
        }
    }

    // ========================================================================
    // Worker Loop — modeled on generic_worker_loop (base.rs:4379)
    // ========================================================================

    /// The main worker loop — phase-aware, non-blocking, with held-item pattern.
    ///
    /// Follows the unified pipeline's `generic_worker_loop` pattern exactly:
    /// 1. Check shutdown
    /// 2. Check completion (safe exit requires no held items)
    /// 3. Try to advance ALL held items first (deadlock prevention)
    /// 4. Get priorities based on backpressure (queue depths)
    /// 5. Try owned exclusive step first (prevents starvation)
    /// 6. Try each priority step (break on first success)
    /// 7. Backoff with jitter
    fn worker_loop(
        shared: &SharedPipelineState,
        worker: &mut SortWorkerState,
        pstats: &SortPipelineStats,
    ) {
        loop {
            // 1. Check shutdown
            let current_phase = shared.phase.load(Ordering::Acquire);
            if current_phase == phase::SHUTDOWN {
                break;
            }

            // 1b. Honor the per-phase active-worker cap: workers above the limit
            //     take no NEW work, but must still drain any items they already
            //     hold — held items are per-worker, so a capped worker that
            //     froze with held work would strand that output (no other worker
            //     can advance it). So advance held items first, then idle without
            //     acquiring fresh work. Capped workers re-check the limit within
            //     MAX_BACKOFF_US, so a later raise reactivates them. Default
            //     limit == num_workers ⇒ no worker is ever capped (no-op).
            if worker.worker_id >= shared.active_worker_limit.load(Ordering::Acquire) {
                if Self::try_advance_all_held(shared, worker) {
                    worker.backoff_us = MIN_BACKOFF_US;
                } else {
                    sleep_with_jitter(worker.backoff_us, worker.worker_id, worker.idle_iter);
                    worker.idle_iter = worker.idle_iter.wrapping_add(1);
                    worker.backoff_us = (worker.backoff_us * 2).min(MAX_BACKOFF_US);
                }
                continue;
            }

            // 2. Check phase completion — wait for next phase, only exit on SHUTDOWN.
            //    Workers must survive across Phase 1 → Phase 2 transitions.
            if Self::is_phase_complete(shared, current_phase) && !worker.has_any_held_items() {
                sleep_with_jitter(worker.backoff_us, worker.worker_id, worker.idle_iter);
                worker.idle_iter = worker.idle_iter.wrapping_add(1);
                worker.backoff_us = (worker.backoff_us * 2).min(MAX_BACKOFF_US);
                continue;
            }

            let mut did_work = false;

            // 3. Try to advance ALL held items first (deadlock prevention)
            did_work |= Self::try_advance_all_held(shared, worker);

            // 4. Get backpressure state and resolve priorities (done inline in step 6)
            let owned_step = Self::exclusive_step_for(worker.worker_id, shared, current_phase);

            // 5. Try owned exclusive step first (prevents starvation)
            if !did_work
                && let Some(step) = owned_step
                && Self::is_step_eligible(step, shared, worker, current_phase)
            {
                let t0 = Instant::now();
                let result = Self::execute_step(shared, worker, step);
                if result == StepResult::Success {
                    pstats.record_step(worker.worker_id, step, Self::nanos_u64(t0.elapsed()));
                    did_work = true;
                }
            }

            // 6. Try each priority step (break on first success or OutputFull)
            //    Only compute backpressure/priorities when needed (skip if work already found)
            if !did_work {
                let bp = shared.get_backpressure();
                let priorities = get_sort_priorities(&bp);

                for &step in priorities {
                    if !Self::is_step_eligible(step, shared, worker, current_phase) {
                        continue;
                    }
                    // Skip exclusive steps this worker doesn't own
                    if Self::is_exclusive_step(step)
                        && !Self::can_attempt_exclusive(owned_step, step, shared)
                    {
                        continue;
                    }
                    // Skip the owned step (already tried above)
                    if owned_step == Some(step) {
                        continue;
                    }

                    let t0 = Instant::now();
                    match Self::execute_step(shared, worker, step) {
                        StepResult::Success => {
                            pstats.record_step(
                                worker.worker_id,
                                step,
                                Self::nanos_u64(t0.elapsed()),
                            );
                            did_work = true;
                            break; // Restart priority evaluation
                        }
                        StepResult::OutputFull => break, // Try downstream via held-item advancement
                        StepResult::InputEmpty => {}     // Try next step
                    }
                }
            }

            // 7. Backoff with jitter (ported from unified pipeline)
            if did_work {
                worker.backoff_us = MIN_BACKOFF_US;
            } else {
                let idle_start = Instant::now();
                sleep_with_jitter(worker.backoff_us, worker.worker_id, worker.idle_iter);
                worker.idle_iter = worker.idle_iter.wrapping_add(1);
                worker.backoff_us = (worker.backoff_us * 2).min(MAX_BACKOFF_US);
                pstats.record_idle(worker.worker_id, Self::nanos_u64(idle_start.elapsed()));
            }
        }
    }

    /// Check if the current phase is "complete" (no more work to do).
    ///
    /// This does NOT mean the worker should exit — it must also have no held items.
    fn is_phase_complete(shared: &SharedPipelineState, current_phase: u8) -> bool {
        match current_phase {
            phase::PHASE1 => {
                shared.decompressed_input_done.load(Ordering::Acquire)
                    && shared.compress_queue.is_empty()
            }
            phase::PHASE2 => {
                if !shared.all_chunks_eof.load(Ordering::Acquire) {
                    return false;
                }
                if !shared.compress_queue.is_empty() {
                    return false;
                }
                let files = shared.phase2_files_snapshot();
                files.iter().all(Phase2FileState::is_drained)
            }
            // Legacy mode never "completes" — it runs until phase changes
            _ => false,
        }
    }

    // ========================================================================
    // Exclusive Step Ownership (plan §Exclusive Step Ownership)
    // ========================================================================

    /// Return the exclusive step this worker owns, if any.
    ///
    /// For `num_workers >= 2`: Worker 0 owns `ReadInputBlocks` (Phase 1) so
    /// only one worker at a time can hold the input file lock. For
    /// `num_workers == 1`: the single worker does everything (returns `None`,
    /// no ownership restrictions). Phase 2 has no exclusive step — work
    /// stealing across files handles contention via per-file mutexes.
    fn exclusive_step_for(
        worker_id: usize,
        shared: &SharedPipelineState,
        current_phase: u8,
    ) -> Option<SortStep> {
        if shared.num_workers < 2 {
            return None; // Single worker does everything
        }
        if worker_id != 0 {
            return None; // Only worker 0 owns read steps
        }
        match current_phase {
            phase::PHASE1 => Some(SortStep::ReadInputBlocks),
            _ => None,
        }
    }

    /// Whether a step is exclusive (requires ownership).
    ///
    /// Only `ReadInputBlocks` is exclusive — it reads from a shared input file
    /// protected by a Mutex.
    fn is_exclusive_step(step: SortStep) -> bool {
        matches!(step, SortStep::ReadInputBlocks)
    }

    /// Whether this worker can attempt an exclusive step it doesn't own.
    ///
    /// Non-owner workers can attempt exclusive steps only if `num_workers == 1`
    /// (single worker mode). `owned_step` is pre-computed by the caller to avoid
    /// recomputing `exclusive_step_for` on every step in the priority loop.
    fn can_attempt_exclusive(
        owned_step: Option<SortStep>,
        step: SortStep,
        shared: &SharedPipelineState,
    ) -> bool {
        // Owner can always attempt
        if owned_step == Some(step) {
            return true;
        }
        // Single worker mode: no restrictions
        shared.num_workers < 2
    }

    // ========================================================================
    // Step Eligibility and Dispatch
    // ========================================================================

    /// Check whether a step is eligible to run in the current state.
    fn is_step_eligible(
        step: SortStep,
        shared: &SharedPipelineState,
        worker: &SortWorkerState,
        current_phase: u8,
    ) -> bool {
        match step {
            SortStep::ReadInputBlocks => {
                current_phase == phase::PHASE1
                    && !shared.input_eof.load(Ordering::Acquire)
                    && worker.held_raw_input_blocks.is_empty()
            }
            SortStep::DecompressInput => {
                current_phase == phase::PHASE1
                    && worker.held_decompressed_input.is_none()
                    && (!shared.raw_input_blocks.is_empty()
                        || (shared.input_eof.load(Ordering::Acquire)
                            && !shared.decompressed_input_done.load(Ordering::Acquire)))
            }
            SortStep::CompressSpill | SortStep::CompressOutput => !shared.compress_queue.is_empty(),
            SortStep::Phase2FileWork => current_phase == phase::PHASE2,
        }
    }

    /// Dispatch to the appropriate step function.
    fn execute_step(
        shared: &SharedPipelineState,
        worker: &mut SortWorkerState,
        step: SortStep,
    ) -> StepResult {
        match step {
            SortStep::ReadInputBlocks => Self::try_read_input_blocks(shared, worker),
            SortStep::DecompressInput => Self::try_decompress_input(shared, worker),
            // Bind the compressor to the dispatched step, not shared.phase, to avoid
            // the race where a worker pops a spill job then set_phase(PHASE2) fires
            // before the compressor is chosen, causing spill data to be compressed
            // at the output level.
            SortStep::CompressSpill => {
                Self::try_compress(shared, &mut worker.compressor, &mut worker.zstd_compressor)
            }
            SortStep::CompressOutput => Self::try_compress(
                shared,
                &mut worker.output_compressor,
                &mut worker.zstd_compressor,
            ),
            SortStep::Phase2FileWork => Self::try_phase2_file_work(shared, worker),
        }
    }

    // ========================================================================
    // Held-item advancement (deadlock prevention)
    // ========================================================================

    /// Try to push all held items to their output `ArrayQueue`s.
    /// Returns true if any held item was successfully advanced.
    fn try_advance_all_held(shared: &SharedPipelineState, worker: &mut SortWorkerState) -> bool {
        let mut advanced = false;

        // Raw input blocks (batch)
        if !worker.held_raw_input_blocks.is_empty() {
            let before = worker.held_raw_input_blocks.len();
            try_advance_held_batch(&shared.raw_input_blocks, &mut worker.held_raw_input_blocks);
            if worker.held_raw_input_blocks.len() < before {
                advanced = true;
            }
        }

        // Decompressed input (single) — unpark main only when the push succeeds so the
        // main thread drains only when there is actually new data available.
        if worker.held_decompressed_input.is_some() {
            let pushed =
                try_advance_held(&shared.decompressed_input, &mut worker.held_decompressed_input);
            if pushed {
                shared.main_thread_handle.unpark();
                advanced = true;
                // This may have been the last block. Increment `input_blocks_queued` now
                // that the block is actually in the queue (not held), then re-check the
                // done condition. This prevents premature EOF when multiple workers hold
                // blocks simultaneously — decompressed_input_done is only set once all
                // blocks have been queued, not just decompressed.
                let queued = shared.input_blocks_queued.fetch_add(1, Ordering::AcqRel) + 1;
                let total = shared.input_read_serial.load(Ordering::Acquire);
                if shared.input_eof.load(Ordering::Acquire)
                    && shared.raw_input_blocks.is_empty()
                    && queued >= total
                    && !shared.decompressed_input_done.load(Ordering::Acquire)
                {
                    shared.decompressed_input_done.store(true, Ordering::Release);
                    shared.main_thread_handle.unpark();
                }
            }
        }

        advanced
    }

    // ========================================================================
    // Phase 1 Steps
    // ========================================================================

    /// `ReadInputBlocks`: read a batch of raw BGZF blocks from the input BAM file.
    ///
    /// Uses `try_lock` to avoid blocking — only one worker reads at a time.
    /// Blocks that can't be pushed to the `ArrayQueue` are stored in `held_raw_input_blocks`.
    fn try_read_input_blocks(
        shared: &SharedPipelineState,
        worker: &mut SortWorkerState,
    ) -> StepResult {
        if shared.input_eof.load(Ordering::Acquire) {
            return StepResult::InputEmpty;
        }

        // Try to acquire the input file lock (non-blocking)
        let Ok(mut guard) = shared.input_file.try_lock() else {
            return StepResult::InputEmpty; // Another worker is reading
        };

        let Some(reader) = guard.as_mut() else {
            return StepResult::InputEmpty; // No input file set
        };

        // Read a batch of raw BGZF blocks
        let blocks = match read_raw_blocks(reader.as_mut(), INPUT_READ_BATCH_SIZE) {
            Ok(b) => b,
            Err(e) => {
                log::error!("I/O error reading input BAM: {e}");
                shared.input_read_error.store(true, Ordering::Release);
                shared.input_eof.store(true, Ordering::Release);
                shared.main_thread_handle.unpark();
                return StepResult::InputEmpty;
            }
        };

        if blocks.is_empty() {
            shared.input_eof.store(true, Ordering::Release);
            return StepResult::InputEmpty;
        }

        // Reserve this batch's serial range while STILL HOLDING the input lock.
        //
        // Assigning serials after releasing the lock is racy: a worker preempted
        // between `drop(guard)` and the per-block `fetch_add` leaves its already-read
        // blocks unaccounted for in `input_read_serial`, while another worker acquires
        // the lock, reads the empty EOF batch, and sets `input_eof`. The done-check
        // (`input_blocks_queued == input_read_serial && input_eof`) then fires early,
        // the merge finalizes, and this worker's trailing blocks are pushed after
        // termination — silently truncating the sorted output. Reserving under the lock
        // guarantees `input_read_serial` accounts for every block that has been read
        // before `input_eof` can be observed as set.
        // `AcqRel` rather than `Relaxed`, matching the other coordination counters
        // (`input_blocks_queued`, `sources_at_eof`, `decomp_in_flight`): it makes the
        // reservation's visibility explicit at this site instead of resting on an implicit
        // happens-before chain through the input-file mutex, which a future refactor of the
        // locking or EOF signalling could break with no local signal here.
        let batch_len = blocks.len() as u64;
        let base_serial = shared.input_read_serial.fetch_add(batch_len, Ordering::AcqRel);

        // Drop the lock before pushing to queue
        drop(guard);

        dispatch_reserved_blocks(
            base_serial,
            blocks,
            &shared.raw_input_blocks,
            &mut worker.held_raw_input_blocks,
        );

        StepResult::Success
    }

    /// `DecompressInput`: decompress a raw input BGZF block.
    fn try_decompress_input(
        shared: &SharedPipelineState,
        worker: &mut SortWorkerState,
    ) -> StepResult {
        // Don't take new work if we're holding an item
        if worker.held_decompressed_input.is_some() {
            return StepResult::OutputFull;
        }

        let Some((serial, block)) = shared.raw_input_blocks.pop() else {
            // No blocks to decompress. Check if all blocks are done — the EOF
            // condition can only be detected here (not after a successful pop)
            // when the last block was already queued before input_eof was set.
            // Use `input_blocks_queued` (not `input_blocks_decompressed`) to avoid
            // a race where workers with held blocks trigger premature EOF.
            if shared.input_eof.load(Ordering::Acquire)
                && !shared.decompressed_input_done.load(Ordering::Acquire)
            {
                let queued = shared.input_blocks_queued.load(Ordering::Acquire);
                let total = shared.input_read_serial.load(Ordering::Acquire);
                if queued >= total {
                    shared.decompressed_input_done.store(true, Ordering::Release);
                    shared.main_thread_handle.unpark();
                }
            }
            return StepResult::InputEmpty;
        };

        let data = match decompress_block(&block, &mut worker.decompressor) {
            Ok(d) => d,
            Err(e) => {
                log::error!("BGZF decompression error (input block serial {serial}): {e}");
                shared.decompression_error.store(true, Ordering::Release);
                shared.main_thread_handle.unpark();
                return StepResult::InputEmpty;
            }
        };

        let input_eof = shared.input_eof.load(Ordering::Acquire);
        let raw_empty = shared.raw_input_blocks.is_empty();

        // Try to push to decompressed_input ArrayQueue
        let pushed = match shared.decompressed_input.push((serial, data)) {
            Ok(()) => {
                shared.main_thread_handle.unpark();
                true
            }
            Err(item) => {
                worker.held_decompressed_input = Some(item);
                // Queue full — wake main thread to drain so we can push next time
                shared.main_thread_handle.unpark();
                false
            }
        };

        // Only increment `input_blocks_queued` (and check done) when the block is
        // actually in the queue. If held, the count stays low — preventing the race
        // where another worker sees count==total and fires decompressed_input_done
        // prematurely while this block is still in held_decompressed_input.
        // The done check runs again in try_advance_all_held once the block is pushed.
        if pushed {
            let queued = shared.input_blocks_queued.fetch_add(1, Ordering::AcqRel) + 1;
            let total = shared.input_read_serial.load(Ordering::Acquire);
            if input_eof && raw_empty && queued >= total {
                shared.decompressed_input_done.store(true, Ordering::Release);
                shared.main_thread_handle.unpark();
            }
        }

        StepResult::Success
    }

    // ========================================================================
    // Phase 2 Step: work-stealing across spill files
    // ========================================================================

    /// `Phase2FileWork`: do one unit of work for some Phase 2 spill file.
    ///
    /// This is the unified Phase 2 work step. Each call attempts to find a file
    /// where we can productively make progress and does ONE unit of work for it
    /// (either decompress one block or read a batch of raw blocks). On success,
    /// the worker's per-file cursor advances so the next call rotates through
    /// the file set fairly.
    ///
    /// # Per-file state
    ///
    /// Each [`Phase2FileState`] has three independently-locked sub-states:
    /// - `reader`: the disk reader. Held while pulling raw bytes from disk.
    /// - `raw_blocks`: FIFO of raw BGZF blocks waiting to be decompressed.
    /// - `decompressed`: per-file [`ReorderBuffer`] of decompressed blocks the
    ///   main thread will pop in serial order.
    ///
    /// Workers always prefer decompression to disk reads (decompressed blocks
    /// directly feed the merge loop). They use `try_lock` everywhere so a
    /// blocked file never starves the others.
    ///
    /// # Deadlock-free admission
    ///
    /// Workers refuse to pull a new raw block when the per-file decompressed
    /// reorder buffer is at `PHASE2_DECOMP_CAP`, EXCEPT in the case where the
    /// raw FIFO's head matches the buffer's `next_seq` and the buffer is stuck
    /// (`!can_pop`). That is the gap-filler the main thread is waiting for, and
    /// failing to admit it would deadlock the merge.
    #[allow(clippy::too_many_lines)]
    fn try_phase2_file_work(
        shared: &SharedPipelineState,
        worker: &mut SortWorkerState,
    ) -> StepResult {
        let files = shared.phase2_files_snapshot();
        let n = files.len();
        if n == 0 {
            return StepResult::InputEmpty;
        }

        for offset in 0..n {
            let i = (worker.phase2_file_cursor + offset) % n;
            let file = &files[i];

            // -- Try decompression first (highest-value work) ----------------
            // `try_pop_raw_for_decompress` increments `decomp_in_flight` before
            // returning, so the consumer's `is_drained` check sees this work as
            // outstanding even after the raw FIFO becomes empty. We MUST
            // decrement after inserting into the reorder buffer (or on the
            // error path) to keep the counter balanced.
            let popped = Self::try_pop_raw_for_decompress(file);
            if let Some((serial, raw_bytes)) = popped {
                let data = match file.codec {
                    SpillCodec::Bgzf => {
                        let raw_block = RawBgzfBlock { data: raw_bytes };
                        match decompress_block(&raw_block, &mut worker.decompressor) {
                            Ok(d) => d,
                            Err(e) => {
                                log::error!(
                                    "BGZF decompression error (chunk source {i} serial {serial}): {e}"
                                );
                                shared.decompression_error.store(true, Ordering::Release);
                                file.decomp_in_flight.fetch_sub(1, Ordering::AcqRel);
                                shared.main_thread_handle.unpark();
                                worker.phase2_file_cursor = (i + 1) % n;
                                return StepResult::Success;
                            }
                        }
                    }
                    SpillCodec::Zstd => {
                        // Allocate the scratch buffer lazily so BGZF-only sorts
                        // don't pay 256 KiB × num_workers of dead memory.
                        if worker.zstd_decompress_buf.len() < ZSTD_FRAME_DECOMP_CAP {
                            worker.zstd_decompress_buf.resize(ZSTD_FRAME_DECOMP_CAP, 0);
                        }
                        match worker
                            .zstd_decompressor
                            .decompress_to_buffer(&raw_bytes, &mut worker.zstd_decompress_buf)
                        {
                            Ok(n) => {
                                // Copy the `n` decompressed bytes (≤ one
                                // staging-buffer's worth, typically ~65 KB)
                                // into a fresh Vec for the consumer. The
                                // scratch buffer keeps its 256 KiB capacity
                                // so the next frame on this worker reuses it
                                // without reallocating.
                                worker.zstd_decompress_buf[..n].to_vec()
                            }
                            Err(e) => {
                                log::error!(
                                    "zstd decompression error (chunk source {i} serial {serial}): {e}"
                                );
                                shared.decompression_error.store(true, Ordering::Release);
                                file.decomp_in_flight.fetch_sub(1, Ordering::AcqRel);
                                shared.main_thread_handle.unpark();
                                worker.phase2_file_cursor = (i + 1) % n;
                                return StepResult::Success;
                            }
                        }
                    }
                };
                let now_poppable = {
                    let mut dec_guard =
                        file.decompressed.lock().expect("phase2 decompressed mutex poisoned");
                    dec_guard.insert(serial, data);
                    dec_guard.can_pop()
                };
                // Decrement AFTER the insert is published. The unpark below
                // wakes the consumer in case it has been parked waiting on
                // this specific file (now_poppable=true) or is in the
                // is_drained path waiting for in_flight to reach zero.
                file.decomp_in_flight.fetch_sub(1, Ordering::AcqRel);
                if now_poppable || file.is_drained() {
                    // Wake the consumer either because new data is available
                    // or because the last in-flight decompression for this
                    // file just completed and the file is now fully drained.
                    shared.main_thread_handle.unpark();
                }
                worker.phase2_file_cursor = (i + 1) % n;
                return StepResult::Success;
            }

            // -- Try reading raw blocks from disk ----------------------------
            // Skip if disk reader is contended OR already at EOF.
            let Ok(mut reader_guard) = file.reader.try_lock() else {
                continue; // another worker is reading this file
            };
            if reader_guard.eof {
                continue;
            }

            // Bound disk read-ahead per file: don't keep pulling if the raw
            // FIFO is already full. Use try_lock so a momentarily contended
            // raw FIFO doesn't block the reader.
            let raw_full = match file.raw_blocks.try_lock() {
                Ok(g) => g.len() >= PHASE2_RAW_CAP,
                Err(_) => true,
            };
            if raw_full {
                continue;
            }

            let raw_bytes: Vec<Vec<u8>> = match file.codec {
                SpillCodec::Bgzf => {
                    match read_raw_blocks(&mut reader_guard.inner, PHASE2_READ_BATCH) {
                        Ok(blocks) => blocks.into_iter().map(|b| b.data).collect(),
                        Err(e) => {
                            log::error!("I/O error reading chunk file (source {i}): {e}");
                            shared.chunk_read_error.store(true, Ordering::Release);
                            file.mark_reader_eof(&mut reader_guard);
                            drop(reader_guard);
                            shared.main_thread_handle.unpark();
                            Self::maybe_mark_all_eof(shared);
                            worker.phase2_file_cursor = (i + 1) % n;
                            return StepResult::Success;
                        }
                    }
                }
                SpillCodec::Zstd => {
                    match read_raw_zstd_frames(&mut reader_guard.inner, PHASE2_READ_BATCH) {
                        Ok(frames) => frames,
                        Err(e) => {
                            log::error!("I/O error reading zstd chunk file (source {i}): {e}");
                            shared.chunk_read_error.store(true, Ordering::Release);
                            file.mark_reader_eof(&mut reader_guard);
                            drop(reader_guard);
                            shared.main_thread_handle.unpark();
                            Self::maybe_mark_all_eof(shared);
                            worker.phase2_file_cursor = (i + 1) % n;
                            return StepResult::Success;
                        }
                    }
                }
            };

            if raw_bytes.is_empty() {
                file.mark_reader_eof(&mut reader_guard);
                drop(reader_guard);
                shared.main_thread_handle.unpark();
                Self::maybe_mark_all_eof(shared);
                worker.phase2_file_cursor = (i + 1) % n;
                return StepResult::Success;
            }

            // Acquire the raw_blocks lock BEFORE releasing the reader lock and
            // assigning serials. This is critical for FIFO order: if we dropped
            // the reader lock before pushing, two workers could each read a
            // batch and bump `next_serial`, then race on the raw_blocks push,
            // landing higher serials in front of lower ones. The merge
            // consumer's gap-filler admission rule cannot recover from that
            // and would deadlock. Lock order `reader → raw_blocks` is the only
            // nested-lock path in this function.
            let mut raw_guard = file.raw_blocks.lock().expect("phase2 raw_blocks mutex poisoned");
            let start_serial = reader_guard.next_serial;
            reader_guard.next_serial += raw_bytes.len() as u64;
            for (idx, b) in raw_bytes.into_iter().enumerate() {
                raw_guard.push_back((start_serial + idx as u64, b));
            }
            drop(raw_guard);
            drop(reader_guard);

            worker.phase2_file_cursor = (i + 1) % n;
            return StepResult::Success;
        }

        StepResult::InputEmpty
    }

    /// Try to pop a raw block from `file` for decompression, applying
    /// deadlock-free admission control against the file's reorder buffer.
    ///
    /// Returns `Some((serial, raw_block))` if a block was popped, `None` otherwise.
    /// `None` is returned when:
    /// - either lock is contended (`try_lock` failed),
    /// - the raw FIFO is empty,
    /// - or the reorder buffer is at cap and the head raw block isn't a gap-filler.
    ///
    /// On success, `decomp_in_flight` is incremented so the consumer's
    /// `is_drained()` check correctly reflects the in-progress decompression.
    /// The caller is responsible for the matching decrement after inserting
    /// (or on the decompression-error path).
    fn try_pop_raw_for_decompress(file: &Phase2FileState) -> Option<(u64, Vec<u8>)> {
        let mut raw_guard = file.raw_blocks.try_lock().ok()?;
        let head_serial = raw_guard.front().map(|(s, _)| *s)?;

        // Cheap admission check using the per-file reorder buffer.
        // Two cases admit: (1) under cap (normal), (2) reorder buffer is stuck
        // and this serial is the gap-filler. Otherwise: backpressure.
        let admit = {
            let dec_guard = file.decompressed.try_lock().ok()?;
            dec_guard.len() < PHASE2_DECOMP_CAP
                || (!dec_guard.can_pop() && head_serial == dec_guard.next_seq())
        };
        if !admit {
            return None;
        }

        // Reserve the in-flight slot under the raw_blocks lock so the consumer
        // can never observe (raw_empty && in_flight==0 && decompressed_empty)
        // while a worker is still in the middle of decompressing this block.
        let popped = raw_guard.pop_front();
        if popped.is_some() {
            file.decomp_in_flight.fetch_add(1, Ordering::AcqRel);
        }
        popped
    }

    // ========================================================================
    // Compress Step (shared by Phase 1 + Phase 2)
    // ========================================================================

    /// Try to pick up a compress job from the `ArrayQueue` (non-blocking).
    ///
    /// The compressor is passed in by the caller (dispatched from `execute_step`)
    /// so the choice is bound to the scheduled `SortStep`, not `shared.phase`.
    /// This avoids the race where a worker pops a Phase-1 spill job and then
    /// `set_phase(PHASE2)` fires before the compressor is selected.
    fn try_compress(
        shared: &SharedPipelineState,
        bgzf_compressor: &mut InlineBgzfCompressor,
        zstd_compressor: &mut ZstdCompressor<'static>,
    ) -> StepResult {
        let Some(job) = shared.compress_queue.pop() else {
            return StepResult::InputEmpty;
        };
        Self::handle_compress_job(shared, job, bgzf_compressor, zstd_compressor);
        StepResult::Success
    }

    // ========================================================================
    // Helper: EOF tracking for Phase 2
    // ========================================================================

    /// Check if all sources have reached EOF and set the global flag if so.
    fn maybe_mark_all_eof(shared: &SharedPipelineState) {
        let eof_count = shared.sources_at_eof.fetch_add(1, Ordering::AcqRel) + 1;
        let total = shared.total_sources.load(Ordering::Acquire);
        if total > 0 && eof_count >= total {
            shared.all_chunks_eof.store(true, Ordering::Release);
            shared.main_thread_handle.unpark();
        }
    }

    // ========================================================================
    // Public API
    // ========================================================================

    /// Number of worker threads in the pool.
    pub fn num_workers(&self) -> usize {
        self.num_workers
    }

    /// Codec used to compress Phase 1 spill chunks.
    #[must_use]
    pub fn spill_codec(&self) -> SpillCodec {
        self.spill_codec
    }

    /// Phase 1 input pipeline queue depths: `(raw_input_blocks, decompressed_input, buffer_pool)`.
    pub(crate) fn phase1_queue_depths(&self) -> (usize, usize, usize) {
        (
            self.shared.raw_input_blocks.len(),
            self.shared.decompressed_input.len(),
            self.buffer_pool.len(),
        )
    }

    /// Get a clone of the decompressed input `ArrayQueue` for `PooledInputStream`.
    pub(crate) fn decompressed_input_queue(
        &self,
    ) -> Arc<crossbeam_queue::ArrayQueue<(u64, Vec<u8>)>> {
        Arc::clone(&self.shared.decompressed_input)
    }

    /// Get a clone of the decompressed input done flag for `PooledInputStream`.
    pub(crate) fn decompressed_input_done_flag(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.shared.decompressed_input_done)
    }

    /// Get the input read error flag for `PooledInputStream` error surfacing.
    pub(crate) fn input_read_error_flag(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.shared.input_read_error)
    }

    /// Get the chunk read error flag for chunk consumer error surfacing.
    pub(crate) fn chunk_read_error_flag(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.shared.chunk_read_error)
    }

    /// Get the worker-panicked flag for chunk consumer error surfacing.
    pub(crate) fn worker_panicked_flag(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.shared.worker_panicked)
    }

    /// Get the decompression error flag for `PooledInputStream` and chunk consumer error surfacing.
    pub(crate) fn decompress_error_flag(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.shared.decompression_error)
    }

    /// Snapshot the Phase 2 per-file state vector for the merge consumer.
    ///
    /// Returns the same `Arc` workers see — the consumer reads from per-file
    /// reorder buffers via this snapshot.
    pub(crate) fn phase2_files(&self) -> Arc<Vec<Phase2FileState>> {
        self.shared.phase2_files_snapshot()
    }

    /// Set the current pipeline phase.
    pub fn set_phase(&self, new_phase: u8) {
        self.shared.phase.store(new_phase, Ordering::Release);
    }

    /// Cap the number of workers active in the current phase to `n` (clamped to
    /// `[1, num_workers]`). Workers with `worker_id >= n` idle until the limit is
    /// raised. Used to run Phase 1 on fewer threads than Phase 2; raising the
    /// limit reactivates idled workers within `MAX_BACKOFF_US` (no explicit wake
    /// needed). Defaults to `num_workers` (all active) when never called.
    pub fn set_active_workers(&self, n: usize) {
        let n = n.clamp(1, self.num_workers);
        self.shared.active_worker_limit.store(n, Ordering::Release);
    }

    /// Hand the pool over to Phase 2 once ingest is done, widening it to
    /// `active_workers`.
    ///
    /// Both halves belong to the same transition and must happen together.
    /// Widening alone was not enough: the phase is what makes workers schedule
    /// [`SortStep::CompressOutput`], and that step is what selects
    /// `output_compressor` over the Phase 1 spill compressor. A sort that stayed in
    /// Phase 1 through its output write therefore compressed the output BAM at
    /// `temp_compression`, silently discarding the caller's `output_compression`.
    ///
    /// Call this once ingest has finished — after that point every remaining byte
    /// the pool touches is merge input or output. Callers must not leave spill
    /// compression outstanding: Phase 1 spill jobs still queued when the phase flips
    /// would be picked up by `CompressOutput` and written at the output level.
    pub fn begin_phase2(&self, active_workers: usize) {
        self.set_active_workers(active_workers);
        self.set_phase(phase::PHASE2);
    }

    /// Set the input file for Phase 1 reading.
    ///
    /// Must be called before `set_phase(PHASE1)`.
    ///
    /// # Panics
    ///
    /// Panics if the input file mutex is poisoned.
    pub fn set_input_file(&self, reader: Box<dyn Read + Send>) {
        *self.shared.input_file.lock().expect("input_file mutex should not be poisoned") =
            Some(reader);
    }

    /// Build the Phase 2 per-file state vector and publish it to all workers.
    ///
    /// Workers do not own files — they cooperatively scan all files and steal
    /// work via `try_lock`. Must be called before `set_phase(PHASE2)`.
    ///
    /// # Errors
    ///
    /// Returns an error if a chunk file cannot be opened.
    ///
    /// # Panics
    ///
    /// Panics if the `phase2_files` rwlock is poisoned.
    pub fn set_phase2_files(&self, files: &[std::path::PathBuf]) -> anyhow::Result<()> {
        let total_sources = files.len();
        self.shared.total_sources.store(total_sources as u64, Ordering::Release);

        // Reset EOF state
        self.shared.all_chunks_eof.store(false, Ordering::Release);
        self.shared.sources_at_eof.store(0, Ordering::Release);

        let mut states: Vec<Phase2FileState> = Vec::with_capacity(total_sources);
        for path in files {
            let mut file = std::fs::File::open(path).map_err(|e| {
                anyhow::anyhow!("Failed to open chunk file {}: {e}", path.display())
            })?;
            // Detect codec by peeking at the first 4 bytes. BGZF starts with
            // 0x1f 0x8b; zstd-spill starts with "ZSP1". On zstd, consume the
            // four magic bytes so the reader is positioned at the first frame.
            //
            // Use `read_exact_or_eof` so a short `read()` (legal for `Read`)
            // can't truncate `ZSPILL_MAGIC` and silently misroute a zstd spill
            // through the BGZF fallback. Clean EOF preserves the existing
            // BGZF-fallback behavior (empty file → BGZF reader will error).
            let mut magic = [0u8; 4];
            let read_n =
                if crate::external::read_exact_or_eof(&mut file, &mut magic).map_err(|e| {
                    anyhow::anyhow!("Failed to read chunk magic {}: {e}", path.display())
                })? {
                    4
                } else {
                    0
                };
            let codec = SpillCodec::from_magic(&magic[..read_n]).unwrap_or(SpillCodec::Bgzf);
            // Zstd consumes the magic itself; BGZF wants the file rewound to
            // byte 0 since its decoder reads the gzip header.
            let body_start = match codec {
                SpillCodec::Bgzf => 0,
                SpillCodec::Zstd => ZSPILL_MAGIC.len() as u64,
            };
            file.seek(SeekFrom::Start(body_start)).map_err(|e| {
                anyhow::anyhow!("Failed to seek chunk file {}: {e}", path.display())
            })?;
            let reader = BufReader::with_capacity(2 * 1024 * 1024, file);
            states.push(Phase2FileState::new(reader, codec));
        }

        let mut guard = self.shared.phase2_files.write().expect("phase2_files rwlock poisoned");
        *guard = Arc::new(states);
        Ok(())
    }

    /// Clear the Phase 2 file vector. Call this after Phase 2 finishes (and
    /// before any subsequent Phase 1) so the file descriptors are released.
    ///
    /// # Panics
    ///
    /// Panics if the `phase2_files` rwlock is poisoned.
    pub fn clear_phase2_files(&self) {
        let mut guard = self.shared.phase2_files.write().expect("phase2_files rwlock poisoned");
        *guard = Arc::new(Vec::new());
    }

    /// Submit a compression job to the pool (non-blocking, spin-yield on full).
    ///
    /// The main thread calls this during spill writes and merge output. If the
    /// compress `ArrayQueue` is full, spins briefly with `yield_now()` — acceptable
    /// because the main thread has no other productive work during spill writes.
    pub fn submit_compress(&self, job: CompressJob) {
        self.stats.compress_jobs_submitted.fetch_add(1, Ordering::Relaxed);
        let mut job = job;
        loop {
            if self.shared.phase.load(Ordering::Acquire) == phase::SHUTDOWN {
                return; // Workers have exited; no one will pop the queue
            }
            match self.shared.compress_queue.push(job) {
                Ok(()) => return,
                Err(returned) => {
                    job = returned;
                    std::thread::yield_now();
                }
            }
        }
    }

    /// Create a new result channel pair for compress results.
    ///
    /// The result channel stays as `crossbeam_channel::bounded()` because the
    /// I/O writer thread needs blocking `recv()`.
    #[allow(dead_code)]
    pub fn compress_result_channel(&self) -> (Sender<CompressResult>, Receiver<CompressResult>) {
        bounded(self.num_workers * 2)
    }

    /// Shut down the pool, waiting for all workers to finish.
    ///
    /// Logs pipeline statistics before joining workers. After this call the pool is fully
    /// stopped. It is also safe to simply drop the pool — `Drop` performs the same cleanup
    /// (minus the debug logging) if `shutdown` was not called explicitly.
    pub fn shutdown(mut self) {
        if log::log_enabled!(log::Level::Debug) {
            self.stats.log_summary();
            self.pipeline_stats.log_summary();
        }
        self.do_shutdown();
    }

    /// Internal shutdown: signal workers and join them. Safe to call multiple times
    /// (idempotent via `Option::take`). Called by both `shutdown` and `Drop`.
    fn do_shutdown(&mut self) {
        self.shared.phase.store(phase::SHUTDOWN, Ordering::Release);
        if let Some(workers) = self.workers.take() {
            for w in workers {
                if w.join().is_err() {
                    // Worker panicked — set flag and wake main thread so it doesn't
                    // park forever waiting for work that will never arrive.
                    self.shared.worker_panicked.store(true, Ordering::Release);
                    self.shared.main_thread_handle.unpark();
                }
            }
        }
    }

    // ========================================================================
    // Job handlers
    // ========================================================================

    /// Convert `Duration::as_nanos()` (u128) to u64 nanoseconds for stats.
    #[allow(clippy::cast_possible_truncation)]
    fn nanos_u64(d: std::time::Duration) -> u64 {
        d.as_nanos() as u64
    }

    /// Handle a compress job on a worker thread.
    fn handle_compress_job(
        shared: &SharedPipelineState,
        job: CompressJob,
        bgzf_compressor: &mut InlineBgzfCompressor,
        zstd_compressor: &mut ZstdCompressor<'static>,
    ) {
        let compressed: Vec<u8> = match job.codec {
            SpillCodec::Bgzf => {
                bgzf_compressor
                    .write_all(&job.data)
                    .expect("BGZF compression write should not fail for valid data");
                bgzf_compressor.flush().expect("BGZF compression flush should not fail");
                let blocks = bgzf_compressor.take_blocks();
                let mut out = Vec::with_capacity(blocks.iter().map(|b| b.data.len()).sum());
                for block in &blocks {
                    out.extend_from_slice(&block.data);
                }
                out
            }
            SpillCodec::Zstd => {
                // One self-contained zstd frame per job, length-prefixed so the
                // reader can split frames without scanning the stream.
                let frame =
                    zstd_compressor.compress(&job.data).expect("zstd compression should not fail");
                let frame_len = u32::try_from(frame.len())
                    .expect("zstd frame larger than 4 GiB cannot fit in a u32 length prefix");
                let mut out = Vec::with_capacity(4 + frame.len());
                out.extend_from_slice(&frame_len.to_le_bytes());
                out.extend_from_slice(&frame);
                out
            }
        };

        let mut recycled = job.data;
        recycled.clear();

        let serial = job.serial;
        let mut result = CompressResult { serial, compressed, recycled_buf: recycled };

        // Use try_send in a yield loop rather than blocking send() so workers
        // remain responsive to SHUTDOWN during the result phase. A blocking send()
        // on a full bounded channel would prevent do_shutdown() from joining the
        // worker if the writer stopped draining before dropping its receiver.
        loop {
            match job.result_tx.try_send(result) {
                Ok(()) => return,
                Err(crossbeam_channel::TrySendError::Disconnected(_)) => {
                    log::warn!(
                        "compress result discarded (serial {serial}): I/O writer thread disconnected"
                    );
                    return;
                }
                Err(crossbeam_channel::TrySendError::Full(r)) => {
                    if shared.phase.load(Ordering::Acquire) == phase::SHUTDOWN {
                        return; // Abandon on shutdown to unblock do_shutdown() join
                    }
                    result = r;
                    std::thread::yield_now();
                }
            }
        }
    }
}

impl Drop for SortWorkerPool {
    /// Ensures workers are joined even if `shutdown` was not called explicitly.
    ///
    /// This prevents thread leaks on early `?` exits in the sort pipeline. Statistics
    /// are not logged here (only in `shutdown`).
    fn drop(&mut self) {
        self.do_shutdown();
    }
}

/// Push a batch of already-serial-reserved blocks to the shared input queue,
/// holding any the queue cannot accept.
///
/// Serials are pre-reserved under the input-file lock by `try_read_input_blocks`
/// (see `SharedPipelineState::input_blocks_queued`), so this assigns them
/// sequentially from `base_serial`. Every block is accounted for: it is either
/// queued or moved to `held`. Once the queue rejects a push, this block and all
/// remaining ones are held, preserving serial order within the batch.
fn dispatch_reserved_blocks(
    base_serial: u64,
    blocks: Vec<RawBgzfBlock>,
    queue: &ArrayQueue<(u64, RawBgzfBlock)>,
    held: &mut Vec<(u64, RawBgzfBlock)>,
) {
    let mut next_serial = base_serial;
    let mut blocks_iter = blocks.into_iter();
    for block in blocks_iter.by_ref() {
        let serial = next_serial;
        next_serial += 1;
        match queue.push((serial, block)) {
            Ok(()) => {}
            Err((serial, block)) => {
                held.push((serial, block));
                break;
            }
        }
    }
    // Hold any remaining blocks we didn't attempt to push.
    for block in blocks_iter {
        let serial = next_serial;
        next_serial += 1;
        held.push((serial, block));
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    /// Build `n` distinguishable placeholder blocks.
    fn dummy_blocks(n: usize) -> Vec<RawBgzfBlock> {
        (0..n).map(|i| RawBgzfBlock { data: vec![u8::try_from(i % 256).unwrap()] }).collect()
    }

    /// Every reserved serial must be accounted for exactly once, whether the block
    /// lands in the queue or is held — a serial that is reserved but never attached
    /// to a block would leave `input_read_serial` permanently ahead of
    /// `input_blocks_queued` and stall the done-check. When the queue fills
    /// mid-batch, the rejected block and every block after it are held in ascending
    /// serial order.
    #[rstest]
    // Queue has room for the whole batch: everything is queued, nothing held.
    #[case::all_fit(8, 100, 5, vec![100, 101, 102, 103, 104], vec![])]
    // Queue fills after two blocks: the third is rejected and the 4th/5th follow it.
    #[case::queue_fills_mid_batch(2, 10, 5, vec![10, 11], vec![12, 13, 14])]
    // An empty batch must not consume serials or touch the queue.
    #[case::empty_batch(4, 7, 0, vec![], vec![])]
    fn test_dispatch_reserved_blocks_accounts_for_every_serial(
        #[case] capacity: usize,
        #[case] base_serial: u64,
        #[case] num_blocks: usize,
        #[case] expected_queued: Vec<u64>,
        #[case] expected_held: Vec<u64>,
    ) {
        let queue: ArrayQueue<(u64, RawBgzfBlock)> = ArrayQueue::new(capacity);
        let mut held = Vec::new();

        dispatch_reserved_blocks(base_serial, dummy_blocks(num_blocks), &queue, &mut held);

        let mut queued: Vec<u64> = std::iter::from_fn(|| queue.pop()).map(|(s, _)| s).collect();
        let held_serials: Vec<u64> = held.iter().map(|(s, _)| *s).collect();

        queued.sort_unstable();
        assert_eq!(queued, expected_queued, "queued serials");
        assert_eq!(held_serials, expected_held, "held blocks keep ascending serial order");

        // The union must be exactly the reserved range — no gaps, no duplicates.
        queued.extend(&held_serials);
        queued.sort_unstable();
        let reserved: Vec<u64> = (base_serial..base_serial + num_blocks as u64).collect();
        assert_eq!(queued, reserved, "reserved range fully accounted for");
    }

    // ── PermitPool fail-closed ───────────────────────────────────────────────

    #[test]
    fn permit_pool_acquire_fails_immediately_after_close() {
        // Buffered permits remain in the channel after close. Without the
        // closed flag, `recv` would hand them out and each producer would
        // compress another block before noticing the writer is gone.
        let pool = PermitPool::new(4);
        pool.acquire().expect("permit available before close");

        pool.close();

        for attempt in 0..4 {
            let err = pool
                .acquire()
                .expect_err("acquire must fail once closed, even with buffered permits");
            assert!(
                err.to_string().contains("permit pool closed"),
                "attempt {attempt} gave an unexpected error: {err}",
            );
        }
    }

    #[test]
    fn permit_pool_release_after_close_is_a_noop() {
        let pool = PermitPool::new(2);
        pool.close();
        pool.release(); // must not panic or resurrect the pool
        assert!(pool.acquire().is_err(), "a release after close must not make permits available");
    }

    #[test]
    fn permit_pool_hands_out_permits_up_to_capacity() {
        let pool = PermitPool::new(2);
        pool.acquire().expect("first permit");
        pool.acquire().expect("second permit");
        pool.release();
        pool.acquire().expect("released permit is reusable");
    }

    #[test]
    fn permit_pool_close_unblocks_parked_acquirers_with_error() {
        // Drain every permit so the acquirers below must park in `recv()`.
        let pool = Arc::new(PermitPool::new(2));
        pool.acquire().expect("first permit");
        pool.acquire().expect("second permit");

        // `close()` must wake every parked producer with an error and never hand
        // out a permit: whichever side of the race a producer lands on, the
        // pre-check or the disconnected `recv` returns the closure error. The
        // join also proves no producer hangs after the writer is gone.
        let handles: Vec<_> = (0..4)
            .map(|_| {
                let pool = Arc::clone(&pool);
                std::thread::spawn(move || pool.acquire())
            })
            .collect();

        pool.close();

        for handle in handles {
            let result = handle.join().expect("acquirer thread must not panic");
            let err = result.expect_err("a parked acquirer must fail once the pool is closed");
            assert!(
                err.to_string().contains("permit pool closed"),
                "unexpected error from a closed pool: {err}",
            );
        }
    }

    // ── Worker panic surfaces mid-run ────────────────────────────────────────

    /// A worker that panics must publish `worker_panicked` immediately, not at
    /// join time. `do_shutdown` sets the same flag, but only after Phase 2; a
    /// mid-run panic reported that late leaves the merge loop spinning while
    /// the main thread waits for work that will never arrive.
    ///
    /// The wait is bounded so a regression fails fast instead of hanging CI.
    #[test]
    fn worker_panic_sets_flag_before_shutdown() {
        let shared = Arc::new(SharedPipelineState::new(1, std::thread::current()));
        let flag = Arc::clone(&shared.worker_panicked);
        assert!(!flag.load(Ordering::Acquire), "flag starts clear");

        let guard_shared = Arc::clone(&shared);
        let handle = std::thread::spawn(move || {
            let _panic_guard = WorkerPanicGuard::new(guard_shared);
            panic!("simulated worker panic");
        });

        // Bounded wait: the flag must be set by the time the panicking thread
        // has unwound, which is well before any pool shutdown would run.
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(5);
        while !flag.load(Ordering::Acquire) && std::time::Instant::now() < deadline {
            std::thread::sleep(std::time::Duration::from_millis(5));
        }

        assert!(
            flag.load(Ordering::Acquire),
            "worker_panicked must be set by the unwinding worker, not deferred to join",
        );
        assert!(handle.join().is_err(), "the worker thread should have panicked");
    }

    /// The guard must not report a panic on the normal exit path.
    #[test]
    fn disarmed_worker_panic_guard_does_not_set_flag() {
        let shared = Arc::new(SharedPipelineState::new(1, std::thread::current()));
        let flag = Arc::clone(&shared.worker_panicked);

        let guard = WorkerPanicGuard::new(Arc::clone(&shared));
        guard.disarm();

        assert!(!flag.load(Ordering::Acquire), "a clean worker exit must not report a panic");
    }

    /// `disarm` must let the guard drop normally so its `Arc<SharedPipelineState>`
    /// is released. The old `mem::forget` path leaked one `Arc` per clean worker
    /// exit, pinning the shared pipeline state (and its buffers) for the process.
    #[test]
    fn disarmed_worker_panic_guard_releases_its_shared_arc() {
        let shared = Arc::new(SharedPipelineState::new(1, std::thread::current()));
        assert_eq!(Arc::strong_count(&shared), 1, "only the test holds the Arc");

        let guard = WorkerPanicGuard::new(Arc::clone(&shared));
        assert_eq!(Arc::strong_count(&shared), 2, "the guard holds a second reference");

        guard.disarm();
        assert_eq!(
            Arc::strong_count(&shared),
            1,
            "disarm must drop the guard's Arc rather than leaking it via mem::forget",
        );
    }

    #[test]
    fn test_buffer_pool_checkout_empty() {
        let pool = BufferPool::new(4);
        let buf = pool.checkout();
        assert!(buf.is_empty());
    }

    #[test]
    fn test_buffer_pool_recycle() {
        let pool = BufferPool::new(4);
        let mut buf = Vec::with_capacity(1024);
        buf.extend_from_slice(b"hello");
        pool.checkin(buf);

        let recycled = pool.checkout();
        // Buffer should be cleared but retain capacity
        assert!(recycled.is_empty());
        assert!(recycled.capacity() >= 1024);
    }

    #[test]
    fn test_pool_stats_log_summary() {
        let stats = PoolStats::default();
        stats.compress_jobs_submitted.fetch_add(42, Ordering::Relaxed);
        // Verify log_summary executes without panicking (logging may be a no-op in tests).
        stats.log_summary();
        assert_eq!(stats.compress_jobs_submitted.load(Ordering::Relaxed), 42);
    }

    #[test]
    fn test_pool_compress_roundtrip() {
        let pool = SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf);
        let (result_tx, result_rx) = pool.compress_result_channel();

        // Submit a compress job
        let data = vec![b'A'; 1000];
        pool.submit_compress(CompressJob { data, serial: 0, result_tx, codec: SpillCodec::Bgzf });

        // Wait for result
        let result = result_rx.recv().expect("should receive compress result");
        assert_eq!(result.serial, 0);
        assert!(!result.compressed.is_empty());
        // Compressed data should start with BGZF magic
        assert_eq!(&result.compressed[0..2], &[0x1f, 0x8b]);
        // Recycled buffer should be empty
        assert!(result.recycled_buf.is_empty());

        pool.shutdown();
    }

    #[test]
    fn test_pool_many_jobs() {
        let pool = SortWorkerPool::new(4, 1, 6, crate::codec::SpillCodec::Bgzf);
        let (result_tx, result_rx) = pool.compress_result_channel();

        let num_jobs = 100usize;

        // Submit in a separate thread to avoid deadlock: the result channel
        // is bounded, so workers block on send when it's full, which blocks
        // the compress channel, which blocks submit.
        let submit_tx = result_tx.clone();
        let submit_handle = std::thread::spawn(move || {
            for i in 0..num_jobs {
                let data = vec![b'X'; 500 + i];
                pool.submit_compress(CompressJob {
                    data,
                    serial: i as u64,
                    result_tx: submit_tx.clone(),
                    codec: SpillCodec::Bgzf,
                });
            }
            drop(submit_tx);
            pool
        });
        drop(result_tx);

        let mut received = 0;
        while let Ok(_result) = result_rx.recv() {
            received += 1;
        }
        assert_eq!(received, num_jobs);

        let pool = submit_handle.join().expect("submit thread should not panic");
        pool.shutdown();
    }

    #[test]
    fn set_active_workers_clamps() {
        let pool = SortWorkerPool::new(4, 1, 6, crate::codec::SpillCodec::Bgzf);
        // Clamp to [1, num_workers].
        pool.set_active_workers(0);
        assert_eq!(pool.shared.active_worker_limit.load(Ordering::Acquire), 1);
        pool.set_active_workers(100);
        assert_eq!(pool.shared.active_worker_limit.load(Ordering::Acquire), 4);
        pool.set_active_workers(2);
        assert_eq!(pool.shared.active_worker_limit.load(Ordering::Acquire), 2);
        pool.shutdown();
    }

    #[test]
    fn active_worker_limit_caps_then_reactivates() {
        // Reactivation must be observable on the SAME pool: a fresh pool can never
        // prove that raising the cap re-enables workers idled by an earlier cap.
        // With a cap of 2 over 6 workers, only workers 0..2 may run compress steps;
        // workers >= 2 stay idle. Raising the cap to 6 and running a second batch
        // must bring those workers online — and because per-worker counts are
        // cumulative and they did nothing under the cap, any work they show after
        // the second batch was done exclusively in that batch.
        let cs = SortStep::CompressSpill as usize;

        // Run one batch on `pool` (moved in, returned out so the next batch can
        // reuse it) and report the cumulative per-worker CompressSpill counts.
        // `cumulative_jobs` is the total submitted across all batches so far —
        // the counters are cumulative, so it is what they must settle at.
        let run_batch = |pool: SortWorkerPool,
                         num_jobs: usize,
                         cumulative_jobs: u64|
         -> (SortWorkerPool, Vec<u64>) {
            let (result_tx, result_rx) = pool.compress_result_channel();
            let submit_tx = result_tx.clone();
            let submit_handle = std::thread::spawn(move || {
                for i in 0..num_jobs {
                    pool.submit_compress(CompressJob {
                        data: vec![b'X'; 4096],
                        serial: i as u64,
                        result_tx: submit_tx.clone(),
                        codec: SpillCodec::Bgzf,
                    });
                }
                drop(submit_tx);
                pool
            });
            drop(result_tx);
            let mut received = 0;
            while result_rx.recv().is_ok() {
                received += 1;
            }
            assert_eq!(received, num_jobs);
            let pool = submit_handle.join().expect("submit thread should not panic");

            // A worker sends the job's result — and drops the job, closing its
            // sender — inside `execute_step`, and only records the step counter
            // after `execute_step` returns. So draining the result channel does
            // NOT imply the counters are up to date: the last worker can still
            // be between the send and the increment. Wait for the counters to
            // settle instead of reading a value that is still in flight.
            let read_counts = || -> Vec<u64> {
                (0..6)
                    .map(|t| {
                        pool.pipeline_stats.per_thread_step_counts[t][cs].load(Ordering::Relaxed)
                    })
                    .collect()
            };
            let deadline = Instant::now() + Duration::from_secs(30);
            let counts = loop {
                let counts = read_counts();
                if counts.iter().sum::<u64>() >= cumulative_jobs {
                    break counts;
                }
                assert!(
                    Instant::now() < deadline,
                    "timed out waiting for step counters to reach {cumulative_jobs}; \
                     per-worker counts {counts:?}"
                );
                std::thread::sleep(Duration::from_micros(50));
            };

            (pool, counts)
        };

        let pool = SortWorkerPool::new(6, 1, 6, crate::codec::SpillCodec::Bgzf);

        // Batch 1: capped at 2 — only workers 0..2 may run.
        pool.set_active_workers(2);
        let (pool, after_capped) = run_batch(pool, 300, 300);
        assert!(
            after_capped.iter().filter(|&&c| c > 0).count() <= 2,
            "cap=2 must bound active workers; per-worker counts {after_capped:?}"
        );
        assert!(
            after_capped[2..].iter().all(|&c| c == 0),
            "workers >= cap must do no work; per-worker counts {after_capped:?}"
        );
        assert_eq!(after_capped.iter().sum::<u64>(), 300, "all batch-1 jobs accounted for");

        // Batch 2: raise the cap on the SAME pool and run again.
        pool.set_active_workers(6);
        let (pool, after_full) = run_batch(pool, 300, 600);
        // Workers >= the old cap did nothing in batch 1, so any cumulative work
        // they now show was processed solely in batch 2 — i.e. reactivation worked.
        let reactivated_work: u64 = after_full[2..].iter().sum();
        assert!(
            reactivated_work > 0,
            "raising the cap on the same pool must reactivate workers >= old cap; \
             after_capped {after_capped:?}, after_full {after_full:?}"
        );
        assert_eq!(
            after_full.iter().sum::<u64>(),
            600,
            "all jobs across both batches accounted for; per-worker counts {after_full:?}"
        );
        pool.shutdown();
    }

    #[test]
    fn test_pool_stats() {
        let pool = SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf);
        let (c_tx, c_rx) = pool.compress_result_channel();

        // Submit one compress job
        pool.submit_compress(CompressJob {
            data: vec![b'A'; 100],
            serial: 0,
            result_tx: c_tx,
            codec: SpillCodec::Bgzf,
        });
        let _ = c_rx.recv();

        assert_eq!(pool.stats.compress_jobs_submitted.load(Ordering::Relaxed), 1);

        pool.shutdown();
    }

    #[test]
    fn test_pipeline_stats_record_step_and_idle() {
        let stats = SortPipelineStats::new(2);

        stats.record_step(0, SortStep::ReadInputBlocks, 1_000_000);
        stats.record_step(0, SortStep::ReadInputBlocks, 500_000);
        stats.record_step(1, SortStep::DecompressInput, 2_000_000);
        stats.record_step(0, SortStep::CompressSpill, 300_000);
        stats.record_idle(0, 100_000);
        stats.record_idle(1, 200_000);

        let read_idx = SortStep::ReadInputBlocks as usize;
        let decomp_idx = SortStep::DecompressInput as usize;
        let compress_idx = SortStep::CompressSpill as usize;

        assert_eq!(stats.step_count[read_idx].load(Ordering::Relaxed), 2);
        assert_eq!(stats.step_ns[read_idx].load(Ordering::Relaxed), 1_500_000);
        assert_eq!(stats.step_count[decomp_idx].load(Ordering::Relaxed), 1);
        assert_eq!(stats.step_count[compress_idx].load(Ordering::Relaxed), 1);
        assert_eq!(stats.per_thread_step_counts[0][read_idx].load(Ordering::Relaxed), 2);
        assert_eq!(stats.per_thread_step_counts[1][decomp_idx].load(Ordering::Relaxed), 1);
        assert_eq!(stats.per_thread_idle_ns[0].load(Ordering::Relaxed), 100_000);
        assert_eq!(stats.per_thread_idle_ns[1].load(Ordering::Relaxed), 200_000);
    }

    #[test]
    fn test_pipeline_stats_log_summary_does_not_panic() {
        let stats = SortPipelineStats::new(4);
        stats.record_step(0, SortStep::ReadInputBlocks, 1_000_000_000);
        stats.record_step(1, SortStep::CompressSpill, 500_000_000);
        stats.record_idle(0, 10_000_000);
        // Verify log_summary doesn't panic (output goes to log, not captured in tests)
        stats.log_summary();
    }

    #[test]
    fn test_buffer_pool_full_drops_excess() {
        let pool = BufferPool::new(2);

        // Put 3 items in a capacity-2 pool; checkin clears len but preserves capacity.
        pool.checkin(Vec::with_capacity(256));
        pool.checkin(Vec::with_capacity(512));
        pool.checkin(Vec::with_capacity(1024)); // silently dropped (pool full)

        // Drain the 2 pooled items — both are len=0 but retain non-zero capacity.
        let a = pool.checkout();
        let b = pool.checkout();
        assert!(
            a.capacity() > 0 || b.capacity() > 0,
            "at least one pooled buffer should retain allocated capacity"
        );

        // Third checkout: pool exhausted → fresh Vec::default() with zero capacity.
        let fresh = pool.checkout();
        assert_eq!(fresh.len(), 0);
        assert_eq!(fresh.capacity(), 0, "fresh allocation has no pre-allocated capacity");
    }

    #[test]
    fn test_sort_priorities_phase1_default_feeds_main_thread() {
        let bp = SortBackpressureState {
            decompressed_input_low: true,
            input_eof: false,
            decompressed_input_done: false,
            compress_has_items: false,
            phase: phase::PHASE1,
        };
        let priorities = get_sort_priorities(&bp);
        assert_eq!(priorities[0], SortStep::DecompressInput);
    }

    #[test]
    fn test_sort_priorities_phase1_compress_backpressure() {
        let bp = SortBackpressureState {
            decompressed_input_low: false,
            input_eof: false,
            decompressed_input_done: false,
            compress_has_items: true,
            phase: phase::PHASE1,
        };
        let priorities = get_sort_priorities(&bp);
        assert_eq!(priorities[0], SortStep::CompressSpill);
    }

    #[test]
    fn test_sort_priorities_phase1_all_done_returns_empty() {
        let bp = SortBackpressureState {
            decompressed_input_low: false,
            input_eof: true,
            decompressed_input_done: true,
            compress_has_items: false,
            phase: phase::PHASE1,
        };
        assert!(get_sort_priorities(&bp).is_empty());
    }

    #[test]
    fn test_sort_priorities_phase2_default_feeds_merge_loop() {
        let bp = SortBackpressureState {
            decompressed_input_low: false,
            input_eof: false,
            decompressed_input_done: false,
            compress_has_items: false,
            phase: phase::PHASE2,
        };
        let priorities = get_sort_priorities(&bp);
        assert_eq!(priorities[0], SortStep::Phase2FileWork);
    }

    #[test]
    fn test_sort_priorities_phase2_compress_backpressure() {
        let bp = SortBackpressureState {
            decompressed_input_low: false,
            input_eof: false,
            decompressed_input_done: false,
            compress_has_items: true,
            phase: phase::PHASE2,
        };
        let priorities = get_sort_priorities(&bp);
        assert_eq!(priorities[0], SortStep::CompressOutput);
    }

    #[test]
    fn test_sort_priorities_phase2_after_eof_still_drains_files() {
        // Even after all_chunks_eof, file work continues until per-file
        // reorder buffers drain — `is_phase_complete` is the actual exit gate.
        let bp = SortBackpressureState {
            decompressed_input_low: false,
            input_eof: false,
            decompressed_input_done: false,
            compress_has_items: false,
            phase: phase::PHASE2,
        };
        let priorities = get_sort_priorities(&bp);
        assert_eq!(priorities[0], SortStep::Phase2FileWork);
    }

    #[test]
    fn test_sort_priorities_legacy_returns_compress_only() {
        let bp = SortBackpressureState {
            decompressed_input_low: false,
            input_eof: false,
            decompressed_input_done: false,
            compress_has_items: true,
            phase: phase::LEGACY,
        };
        let priorities = get_sort_priorities(&bp);
        assert_eq!(priorities.len(), 1);
        assert_eq!(priorities[0], SortStep::CompressSpill);
    }

    #[test]
    fn test_worker_pool_num_workers() {
        let pool = SortWorkerPool::new(3, 1, 6, crate::codec::SpillCodec::Bgzf);
        assert_eq!(pool.num_workers(), 3);
        pool.shutdown();
    }

    // ========================================================================
    // Phase2FileState admission and drain tests
    // ========================================================================

    /// Build an empty `Phase2FileState` for unit-testing the admission rule.
    /// The reader is backed by a temporary empty file — we don't exercise it
    /// here; we only manipulate `raw_blocks` and `decompressed` directly.
    fn empty_phase2_file() -> Phase2FileState {
        let tmp = tempfile::tempfile().expect("failed to create tempfile");
        let reader = BufReader::with_capacity(1024, tmp);
        Phase2FileState::new(reader, SpillCodec::Bgzf)
    }

    /// Build a tiny placeholder raw block whose contents we never decode.
    fn dummy_raw_block(byte: u8) -> Vec<u8> {
        vec![byte; 8]
    }

    #[test]
    fn test_admission_under_cap_admits() {
        let file = empty_phase2_file();
        file.raw_blocks.lock().expect("raw lock").push_back((0, dummy_raw_block(0)));
        let popped = SortWorkerPool::try_pop_raw_for_decompress(&file);
        assert!(popped.is_some(), "under cap with empty reorder buffer should admit");
        assert_eq!(file.decomp_in_flight.load(Ordering::Acquire), 1);
        assert!(file.raw_blocks.lock().expect("raw lock").is_empty());
    }

    #[test]
    fn test_admission_at_cap_poppable_rejects() {
        let file = empty_phase2_file();
        // Fill the reorder buffer to PHASE2_DECOMP_CAP starting at serial 0,
        // so the buffer is poppable (next_seq = 0 is present) AND at cap.
        {
            let mut dec = file.decompressed.lock().expect("dec lock");
            for s in 0..PHASE2_DECOMP_CAP as u64 {
                dec.insert(s, vec![0u8; 4]);
            }
            assert_eq!(dec.len(), PHASE2_DECOMP_CAP);
            assert!(dec.can_pop());
        }
        // Head raw is the next serial after the buffer's contents — not the
        // gap-filler, so admission should reject (consumer is supposed to
        // drain the buffer first).
        file.raw_blocks
            .lock()
            .expect("raw lock")
            .push_back((PHASE2_DECOMP_CAP as u64, dummy_raw_block(1)));
        let popped = SortWorkerPool::try_pop_raw_for_decompress(&file);
        assert!(popped.is_none(), "at cap and poppable should reject (apply backpressure)");
        assert_eq!(file.decomp_in_flight.load(Ordering::Acquire), 0);
        assert_eq!(file.raw_blocks.lock().expect("raw lock").len(), 1);
    }

    #[test]
    fn test_admission_at_cap_stuck_admits_gap_filler() {
        let file = empty_phase2_file();
        // Fill the reorder buffer with serials 1..=PHASE2_DECOMP_CAP, leaving
        // serial 0 as the gap. Buffer is at cap and !can_pop.
        {
            let mut dec = file.decompressed.lock().expect("dec lock");
            for s in 1..=PHASE2_DECOMP_CAP as u64 {
                dec.insert(s, vec![0u8; 4]);
            }
            assert_eq!(dec.len(), PHASE2_DECOMP_CAP);
            assert!(!dec.can_pop(), "buffer should be stuck waiting for serial 0");
        }
        // The head raw is serial 0 — the gap-filler. Admission must take it
        // even though we're at cap, otherwise the consumer deadlocks.
        file.raw_blocks.lock().expect("raw lock").push_back((0, dummy_raw_block(0)));
        let popped = SortWorkerPool::try_pop_raw_for_decompress(&file);
        assert!(popped.is_some(), "at cap and stuck should admit gap-filler at next_seq");
        assert_eq!(file.decomp_in_flight.load(Ordering::Acquire), 1);
    }

    #[test]
    fn test_admission_at_cap_stuck_wrong_head_rejects() {
        let file = empty_phase2_file();
        // Same setup as gap-filler test, but the head raw is NOT the
        // gap-filler — admission must reject and the merge must rely on
        // another file's progress to make this one drainable.
        {
            let mut dec = file.decompressed.lock().expect("dec lock");
            for s in 1..=PHASE2_DECOMP_CAP as u64 {
                dec.insert(s, vec![0u8; 4]);
            }
        }
        file.raw_blocks
            .lock()
            .expect("raw lock")
            .push_back((PHASE2_DECOMP_CAP as u64 + 1, dummy_raw_block(2)));
        let popped = SortWorkerPool::try_pop_raw_for_decompress(&file);
        assert!(popped.is_none(), "at cap, stuck, but head != next_seq should reject");
        assert_eq!(file.decomp_in_flight.load(Ordering::Acquire), 0);
    }

    #[test]
    fn test_admission_empty_raw_returns_none() {
        let file = empty_phase2_file();
        // Empty raw FIFO — try_pop must return None without touching in_flight.
        let popped = SortWorkerPool::try_pop_raw_for_decompress(&file);
        assert!(popped.is_none());
        assert_eq!(file.decomp_in_flight.load(Ordering::Acquire), 0);
    }

    #[test]
    fn test_is_drained_respects_in_flight_counter() {
        let file = empty_phase2_file();
        // Mark reader as EOF and ensure both queues are empty.
        file.mark_reader_eof(&mut file.reader.lock().expect("reader lock"));
        assert!(file.is_drained(), "reader_eof + empty queues + no in-flight should be drained");

        // Simulate a worker mid-decompression: in_flight > 0 must hide drain.
        file.decomp_in_flight.fetch_add(1, Ordering::AcqRel);
        assert!(!file.is_drained(), "in-flight decompression must keep is_drained=false");

        // Decrementing brings us back to drained.
        file.decomp_in_flight.fetch_sub(1, Ordering::AcqRel);
        assert!(file.is_drained());
    }

    #[test]
    fn test_is_drained_blocks_on_pending_raw() {
        let file = empty_phase2_file();
        file.mark_reader_eof(&mut file.reader.lock().expect("reader lock"));
        file.raw_blocks.lock().expect("raw lock").push_back((0, dummy_raw_block(0)));
        assert!(!file.is_drained(), "raw blocks pending must keep is_drained=false");
    }

    #[test]
    fn test_is_drained_blocks_on_pending_decompressed() {
        let file = empty_phase2_file();
        file.mark_reader_eof(&mut file.reader.lock().expect("reader lock"));
        file.decompressed.lock().expect("dec lock").insert(0, vec![1, 2, 3]);
        assert!(!file.is_drained(), "decompressed blocks pending must keep is_drained=false");
    }

    // ========================================================================
    // Zstd-framing parser tests
    //
    // `read_length_prefix` and `read_raw_zstd_frames` are the workhorses of the
    // ZSP1 reader path; both are reused by `ZspillStreamReader`. These tests
    // pin the boundary behaviour (clean EOF vs. truncation, oversized length
    // cap) that the file-format invariants rely on.
    // ========================================================================

    use rstest::rstest;
    use std::io::Cursor;

    #[test]
    fn test_read_length_prefix_clean_eof_returns_none() {
        // No bytes available at all → clean stream end.
        let mut reader = Cursor::new(Vec::<u8>::new());
        let got = read_length_prefix(&mut reader).expect("clean EOF must not be an error");
        assert!(got.is_none(), "empty reader must yield Ok(None), got {got:?}");
    }

    #[rstest]
    #[case(1)]
    #[case(2)]
    #[case(3)]
    fn test_read_length_prefix_partial_prefix_is_unexpected_eof(#[case] partial: usize) {
        // 1–3 bytes after a frame boundary is a truncated length prefix, not a
        // clean EOF; the parser must surface UnexpectedEof so corruption is not
        // silently swallowed.
        let bytes = vec![0xABu8; partial];
        let err = read_length_prefix(&mut Cursor::new(bytes))
            .expect_err("partial prefix must error, partial={partial}");
        assert_eq!(err.kind(), std::io::ErrorKind::UnexpectedEof, "partial={partial}, got {err:?}");
    }

    #[test]
    fn test_read_length_prefix_oversized_is_invalid_data() {
        let bogus_len: u32 =
            u32::try_from(MAX_ZSTD_FRAME_BYTES).expect("MAX_ZSTD_FRAME_BYTES fits u32") + 1;
        let buf = bogus_len.to_le_bytes().to_vec();
        let err =
            read_length_prefix(&mut Cursor::new(buf)).expect_err("oversized length must error");
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    }

    #[test]
    fn test_read_length_prefix_valid() {
        let valid_len: u32 = 1234;
        let buf = valid_len.to_le_bytes().to_vec();
        let got = read_length_prefix(&mut Cursor::new(buf)).expect("valid prefix");
        assert_eq!(got, Some(valid_len as usize));
    }

    #[test]
    fn test_read_raw_zstd_frames_clean_eof_returns_empty() {
        let mut reader = Cursor::new(Vec::<u8>::new());
        let frames = read_raw_zstd_frames(&mut reader, 4).expect("clean EOF should be Ok");
        assert!(frames.is_empty(), "no frames in an empty stream");
    }

    #[test]
    fn test_read_raw_zstd_frames_truncated_body_is_unexpected_eof() {
        // Length prefix promises N body bytes, but only N/2 are present.
        let body = b"hello there, this is a frame body".to_vec();
        let frame_len = u32::try_from(body.len()).expect("fits");
        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(&frame_len.to_le_bytes());
        buf.extend_from_slice(&body[..body.len() / 2]);

        let err =
            read_raw_zstd_frames(&mut Cursor::new(buf), 1).expect_err("truncated body must error");
        assert_eq!(err.kind(), std::io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn test_read_raw_zstd_frames_reads_multiple_frames() {
        // Two opaque "frames" — the parser does not decompress, so any bytes
        // serve as a stand-in for a real zstd frame.
        let first_frame = b"frame zero".to_vec();
        let second_frame = b"frame one (slightly longer)".to_vec();
        let mut buf: Vec<u8> = Vec::new();
        for f in [&first_frame, &second_frame] {
            let len = u32::try_from(f.len()).expect("fits");
            buf.extend_from_slice(&len.to_le_bytes());
            buf.extend_from_slice(f);
        }
        let frames = read_raw_zstd_frames(&mut Cursor::new(buf), 8).expect("two frames");
        assert_eq!(frames, vec![first_frame, second_frame]);
    }

    // ========================================================================
    // set_phase2_files codec-detection tests
    //
    // After `set_phase2_files`, each file's `Phase2FileState.codec` reflects
    // the codec implied by its 4-byte magic, and the underlying reader is
    // positioned past the magic for zstd (which consumes it) and at byte 0
    // for BGZF (whose decoder reads the gzip header itself).
    // ========================================================================

    /// Snapshot the file position by locking the per-file reader. `BufReader`'s
    /// `stream_position` accounts for any buffered bytes — for a fresh
    /// `BufReader` whose buffer hasn't been filled this equals the underlying
    /// `File`'s seek position.
    fn phase2_file_position(state: &Phase2FileState) -> u64 {
        let mut guard = state.reader.lock().expect("reader lock");
        guard.inner.stream_position().expect("stream_position")
    }

    #[test]
    fn test_set_phase2_files_detects_zstd_magic_and_seeks_past_it() {
        use std::io::Write;
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("zstd-spill.zsp");
        let mut file = std::fs::File::create(&path).expect("create");
        file.write_all(&ZSPILL_MAGIC).expect("write magic");
        file.write_all(&[0xAA, 0xBB, 0xCC]).expect("write body");
        drop(file);

        let pool = SortWorkerPool::new(1, 1, 6, SpillCodec::Bgzf);
        pool.set_phase2_files(std::slice::from_ref(&path)).expect("set_phase2_files");
        let files = pool.phase2_files();
        assert_eq!(files.len(), 1);
        assert_eq!(files[0].codec, SpillCodec::Zstd, "ZSPILL_MAGIC must select zstd codec");
        assert_eq!(
            phase2_file_position(&files[0]),
            ZSPILL_MAGIC.len() as u64,
            "zstd reader must be positioned past the 4-byte magic"
        );
        pool.shutdown();
    }

    #[test]
    fn test_set_phase2_files_detects_bgzf_magic_and_keeps_position_zero() {
        use std::io::Write;
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("bgzf-spill.bgzf");
        let mut file = std::fs::File::create(&path).expect("create");
        // BGZF/gzip magic followed by junk; from_magic only inspects the first
        // two bytes and the decoder is not invoked here.
        file.write_all(&[0x1f, 0x8b, 0x00, 0x00, 0x55, 0x66]).expect("write magic");
        drop(file);

        let pool = SortWorkerPool::new(1, 1, 6, SpillCodec::Zstd);
        pool.set_phase2_files(std::slice::from_ref(&path)).expect("set_phase2_files");
        let files = pool.phase2_files();
        assert_eq!(files.len(), 1);
        assert_eq!(files[0].codec, SpillCodec::Bgzf, "BGZF magic must select bgzf codec");
        assert_eq!(
            phase2_file_position(&files[0]),
            0,
            "bgzf reader must be rewound to byte 0 so the decoder sees the header"
        );
        pool.shutdown();
    }

    #[test]
    fn test_set_phase2_files_empty_file_falls_back_to_bgzf() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("empty.spill");
        std::fs::File::create(&path).expect("create empty");

        let pool = SortWorkerPool::new(1, 1, 6, SpillCodec::Zstd);
        pool.set_phase2_files(std::slice::from_ref(&path)).expect("set_phase2_files");
        let files = pool.phase2_files();
        assert_eq!(files.len(), 1);
        assert_eq!(
            files[0].codec,
            SpillCodec::Bgzf,
            "no magic must fall back to bgzf (the legacy decoder reports the truncation)"
        );
        assert_eq!(phase2_file_position(&files[0]), 0);
        pool.shutdown();
    }
}
