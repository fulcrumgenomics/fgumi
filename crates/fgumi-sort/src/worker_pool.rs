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
//!   Phase 1: `DecompressInput` > `ReadInputBlocks` > `Compress`.
//!   Phase 2: `Phase2FileWork` (read+decompress the next block of any file that has
//!   room in its reorder buffer) > `Compress`.
//!   The phase orders the steps; it never decides which *compressor* a block goes
//!   through — that is carried per-job as a `CompressTarget`.
//! - **Per-worker state**: Each worker owns two `InlineBgzfCompressor`s (one at
//!   `temp_compression` for spill blocks, one at `output_compression` for output
//!   blocks, selected by the job's `CompressTarget`), a
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
use crate::merge_stalls::{
    Phase2ScanReport, Phase2ScanTally, PopSkip, ReadSkip, WakeLatencyReport, combine_skip,
};
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
use std::sync::{Arc, Mutex, MutexGuard};
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
/// produced one per ~64 KiB of input by `compress_job`; even
/// pathological expansion can't reach this. Beyond it, we treat the value as
/// corruption rather than allocate gigabytes.
pub(crate) const MAX_ZSTD_FRAME_BYTES: usize = 2 * 1024 * 1024;

/// Maximum zstd compression level recognized by the `zstd` crate.
const ZSTD_MAX_CLEVEL: u32 = 22;

// ============================================================================
// Job and Result Types
// ============================================================================

/// Which of the pool's two compressors a queued block must go through, and so
/// which compression level it is written at.
///
/// One `compress_queue` carries both spill and output blocks, so a worker
/// popping a job has to know which kind it has. That is recorded here, at submit
/// time, by the only code that knows it for certain — the staging buffer of the
/// writer that produced the block. It is deliberately *not* inferred from the
/// pool's current phase: the phase is shared mutable state that advances
/// independently of what is already sitting in the queue, so inferring from it
/// silently mis-compresses every block that outlives the transition (an output
/// block popped in Phase 1 or `LEGACY` at `temp_compression`, a Phase 1 spill
/// block popped after `begin_phase2` at `output_compression`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompressTarget {
    /// A Phase 1 spill chunk — compress at `temp_compression`.
    Spill,
    /// The sort's output BAM — compress at `output_compression`.
    Output,
}

/// A compression job: compress uncompressed data into one compressed block.
///
/// The codec determines the output framing:
/// - `SpillCodec::Bgzf` produces one BGZF block (header + deflate + footer).
/// - `SpillCodec::Zstd` produces `[u32 LE frame-len][zstd frame]`.
///
/// The target determines the compression *level* — see [`CompressTarget`].
pub struct CompressJob {
    /// Uncompressed data to compress.
    pub data: Vec<u8>,
    /// Serial number for ordering output blocks.
    pub serial: u64,
    /// Channel to send the compressed result back.
    pub result_tx: Sender<CompressResult>,
    /// Codec to use when compressing.
    pub codec: SpillCodec,
    /// Which compressor (and therefore which level) this block belongs to.
    pub target: CompressTarget,
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
    /// Compress one queued block — spill or output.
    ///
    /// There is a single `compress_queue`, so this is a single step. Which
    /// compressor the popped block goes through is carried on the job as a
    /// [`CompressTarget`], never derived from the step or the phase.
    Compress = 2,
    /// Read+decompress one unit of work for some Phase 2 spill file (work-stealing).
    Phase2FileWork = 3,
}

impl SortStep {
    /// Number of distinct sort steps.
    pub const COUNT: usize = 4;

    /// Short label for display.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::ReadInputBlocks => "RdInp",
            Self::DecompressInput => "DecInp",
            Self::Compress => "Cmprs",
            Self::Phase2FileWork => "P2File",
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
            SortStep::Compress,
            SortStep::Phase2FileWork,
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
/// Bounds disk read-ahead memory for every file except the one at the drain
/// frontier, which is allowed [`PHASE2_STARVING_RAW_CAP`] instead: K files ×
/// `PHASE2_RAW_CAP` × ~64 KB, plus that single deeper file. Still a function of
/// config rather than input size, because the deep allowance is scoped to one
/// file at a time.
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

/// Read-ahead depth for the single merge source the consumer is draining while
/// it has nothing buffered. See the read path in
/// [`SortWorkerPool::try_phase2_file_work`].
///
/// Applies to one file at a time, so the extra read-ahead memory is this many
/// compressed blocks — a few MB — not K times that.
pub(crate) const PHASE2_STARVING_RAW_CAP: usize = 64;

/// Blocks to read per call for that same file. Large enough that one read is a
/// sequential run rather than a 256 KB pinprick, since only one worker may hold
/// a given file's reader mutex and that read is therefore the whole refill rate.
pub(crate) const PHASE2_STARVING_READ_BATCH: usize = 32;

/// One compressed block in a file's raw FIFO, with when it got there.
///
/// The timestamp is what turns "the pool did not decompress this yet" into a
/// measurement: `enqueued_nanos` to the moment a worker claims it is pure
/// scheduling latency, work that was available and unclaimed.
pub(crate) struct RawEntry {
    /// Position in the file's block sequence.
    pub(crate) serial: u64,
    /// Raw block bytes: a whole BGZF block, or a whole zstd frame.
    pub(crate) bytes: Vec<u8>,
    /// When this entry was pushed, in pool-epoch nanoseconds.
    pub(crate) enqueued_nanos: u64,
}

/// A decompressed block in a file's reorder buffer, with when it was published.
///
/// `inserted_nanos` to the moment the consumer takes it is how long the block
/// spent as genuine lookahead. Near-zero across the run means the reorder
/// buffer is a pass-through and its cap is not the constraint, however full the
/// other files look.
pub(crate) struct TimedBlock {
    /// Decompressed payload.
    pub(crate) data: Vec<u8>,
    /// When this block was inserted, in pool-epoch nanoseconds.
    pub(crate) inserted_nanos: u64,
}

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
    pub(crate) raw_blocks: Mutex<VecDeque<RawEntry>>,
    /// Decompressed blocks reordered by serial. Main thread pops the next-in-order
    /// block here when its parser exhausts the current buffer.
    pub(crate) decompressed: Mutex<ReorderBuffer<TimedBlock>>,
    /// Number of raw blocks currently being decompressed (popped from
    /// `raw_blocks` but not yet inserted into `decompressed`). Used by
    /// `is_drained` to avoid a race where the consumer exits while a worker
    /// is mid-decompress.
    pub(crate) decomp_in_flight: AtomicUsize,

    /// Lock-free mirror of `raw_blocks.len()`, maintained under that mutex.
    ///
    /// Reading a depth should not require taking the lock that a worker might
    /// be holding to *change* it: the park path needs both depths on every
    /// park, and acquiring two mutexes per park would both cost more than the
    /// measurement is worth and perturb the contention it is trying to measure.
    pub(crate) raw_len: AtomicUsize,
    /// Lock-free mirror of `decompressed.len()`, maintained under that mutex.
    pub(crate) decomp_len: AtomicUsize,
    /// Pool-epoch nanoseconds at which the reorder buffer last drained to
    /// nothing, or `0` if it currently holds something.
    ///
    /// The start of a refill cycle: the interval from here to the next claim,
    /// and to the next insert, is what says whether a hungry file waits to be
    /// noticed or waits to be served.
    pub(crate) emptied_at_nanos: AtomicU64,
    /// The [`EmptyCause`](crate::merge_trace::EmptyCause) of the current empty,
    /// as a `u8`.
    ///
    /// Stored on the file rather than only tallied, because `read_lag` is
    /// defined over the `Dry` empties alone: on a `RawReady` or
    /// `Decompressing` empty there was already something in the pipeline, so
    /// the next disk read was never what the consumer was waiting for, and
    /// folding those cycles in makes a scheduling delay look like a storage
    /// one. Meaningless while `emptied_at_nanos` is `0`.
    pub(crate) emptied_cause: AtomicU8,
    /// Whether some worker has already claimed a block since the buffer last
    /// emptied, so `claim_lag` records the *first* response, not every one.
    pub(crate) claimed_since_empty: AtomicBool,
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
            raw_len: AtomicUsize::new(0),
            decomp_len: AtomicUsize::new(0),
            emptied_at_nanos: AtomicU64::new(0),
            emptied_cause: AtomicU8::new(crate::merge_trace::EmptyCause::Dry.as_u8()),
            claimed_since_empty: AtomicBool::new(false),
        }
    }

    /// The file's depths without taking either mutex: `(raw, decompressed,
    /// in-flight)`.
    ///
    /// A momentarily inconsistent triple is fine here and a lock is not: these
    /// feed diagnostics sampled at a park, never a correctness decision.
    pub(crate) fn depths(&self) -> (usize, usize, usize) {
        (
            self.raw_len.load(Ordering::Relaxed),
            self.decomp_len.load(Ordering::Relaxed),
            self.decomp_in_flight.load(Ordering::Relaxed),
        )
    }

    /// Whether the consumer would block on this file right now: its reorder
    /// buffer is empty and nothing is on its way into it.
    ///
    /// Two relaxed atomic loads, no mutex — cheap enough to consult on every
    /// worker scan. Raw blocks may still be queued; that is precisely a file
    /// that needs a worker's attention, not one that has it.
    pub(crate) fn is_starving(&self) -> bool {
        self.decomp_len.load(Ordering::Relaxed) == 0
            && self.decomp_in_flight.load(Ordering::Relaxed) == 0
    }

    /// Note that the reorder buffer just drained to nothing at `now`, for
    /// `cause`.
    ///
    /// Opens a refill cycle. Called by the consumer, which is the only thing
    /// that removes from the buffer, so there is no race to lose here. The
    /// cause is stored as well as tallied so `read_lag` can be restricted to
    /// the empties it is defined over -- see [`Self::emptied_while_dry`].
    pub(crate) fn mark_emptied(&self, now: u64, cause: crate::merge_trace::EmptyCause) {
        self.emptied_cause.store(cause.as_u8(), Ordering::Relaxed);
        self.emptied_at_nanos.store(now, Ordering::Relaxed);
        self.claimed_since_empty.store(false, Ordering::Relaxed);
    }

    /// Whether the current refill cycle opened with the file holding nothing
    /// anywhere -- no queued raw block, no decompression in flight.
    ///
    /// Only these cycles belong in `read_lag`: on the others a block was
    /// already in the pipeline when the buffer drained, so the next disk read
    /// completing was not what the consumer waited on.
    pub(crate) fn emptied_while_dry(&self) -> bool {
        matches!(
            crate::merge_trace::EmptyCause::from_u8(self.emptied_cause.load(Ordering::Relaxed)),
            crate::merge_trace::EmptyCause::Dry
        )
    }

    /// Pool-epoch nanoseconds at which this file's buffer emptied, if it is
    /// still empty.
    pub(crate) fn emptied_at(&self) -> Option<u64> {
        match self.emptied_at_nanos.load(Ordering::Relaxed) {
            0 => None,
            t => Some(t),
        }
    }

    /// Claim the "first response since the buffer emptied" flag.
    ///
    /// Returns `true` for exactly one caller per refill cycle, so several
    /// workers piling onto the same hungry file record one claim latency
    /// between them rather than one each.
    pub(crate) fn take_first_claim_since_empty(&self) -> bool {
        !self.claimed_since_empty.swap(true, Ordering::Relaxed)
    }

    /// Close the refill cycle: the buffer holds something again.
    pub(crate) fn mark_refilled(&self) {
        self.emptied_at_nanos.store(0, Ordering::Relaxed);
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
        let decomp_bytes: u64 = decomp_guard.iter().map(|b| b.data.len() as u64).sum();
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
/// - **Phase 1**: `DecompressInput` > `ReadInputBlocks` > `Compress`
/// - **Phase 2**: `Phase2FileWork` (per-file work stealing: read next raw block
///   batch from any file with FIFO room, OR decompress the next queued raw
///   block for any file whose reorder buffer has capacity) > `Compress`
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
    /// Busy-time counters for the four components of the k-way merge. See
    /// [`crate::merge_phases`] -- these overlap and do not partition wall time.
    pub(crate) merge_phases: crate::merge_phases::MergePhaseCounters,

    /// Why Phase 2 file scans found no work. Complements the idle totals in
    /// [`SortPipelineStats`], which say how long workers idled but not what for.
    pub(crate) phase2_scans: crate::merge_stalls::Phase2ScanStats,

    /// How deep workers were sleeping when they found work, which bounds how
    /// much of their idle time is discovery latency rather than absent work.
    pub(crate) wake_latency: crate::merge_stalls::WakeLatencyStats,

    /// Start of the pool's monotonic clock. Every `*_nanos` field elsewhere is
    /// an offset from here, so timestamps fit in a `u64` atomic instead of
    /// needing a lock around an `Instant`.
    pub(crate) epoch: Instant,
    /// Per-stage dwell times for spill blocks -- see [`crate::merge_trace`].
    pub(crate) block_lifecycle: crate::merge_trace::BlockLifecycleStats,
    /// How long a file that ran dry takes to be refilled, and by which stage.
    pub(crate) refill: crate::merge_trace::RefillStats,
    /// What the merge consumer's access pattern looks like and what it costs.
    pub(crate) consumer_trace: crate::merge_trace::ConsumerTraceStats,
    /// How long a worker's *fruitless* scan of every spill file takes -- the
    /// scans that found no work, which are the ones on the path by which a
    /// hungry file gets noticed. A scan that finds work returns early and is
    /// not timed.
    pub(crate) fruitless_scan: crate::merge_trace::DurationHistogram,

    /// `wake_latency` as it stood when the pool entered Phase 2, so the merge's
    /// share can be reported without Phase 1's sleeps folded in.
    pub(crate) wake_latency_at_merge_start: std::sync::OnceLock<WakeLatencyReport>,

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

    /// Worker thread handles, indexed by `worker_id`, so the consumer can wake
    /// an idle worker instead of waiting out its backoff.
    ///
    /// Each worker publishes its own handle once, at the top of its loop; a
    /// slot that is still empty simply cannot be woken yet, which is harmless
    /// because a worker that has not reached its loop is not sleeping either.
    worker_threads: Vec<std::sync::OnceLock<std::thread::Thread>>,
    /// Rotates the target of [`Self::wake_one_worker`] so repeated wakes spread
    /// across the pool rather than always hitting worker 0.
    wake_cursor: AtomicUsize,

    /// Read batches taken at the deep frontier allowance, and the blocks they
    /// returned; and the same for batches taken at the uniform allowance.
    ///
    /// Blocks-per-batch is the number that says whether the deep path is
    /// actually engaging. Without it, a gate that silently almost never fires
    /// looks exactly like a change that did not help — which is precisely the
    /// mistake this counter exists to prevent.
    pub(crate) deep_read_batches: AtomicU64,
    pub(crate) deep_read_blocks: AtomicU64,
    pub(crate) shallow_read_batches: AtomicU64,
    pub(crate) shallow_read_blocks: AtomicU64,

    /// Lowest merge source index that has not yet delivered all its records.
    ///
    /// Phase 1 chunks its input sequentially, so an input that is already in
    /// the requested order spills coordinate-*disjoint* runs: source 0 holds
    /// the earliest records, source 1 the next, and the merge drains them one
    /// at a time rather than interleaving. This is the frontier of that drain,
    /// and [`SortWorkerPool::try_phase2_file_work`] uses it to look where the
    /// work actually is. Only ever moves forward, and only a scan *hint* — the
    /// scan still wraps over every file, so a stale value costs nothing but a
    /// few microseconds of walking.
    pub(crate) phase2_lowest_active: AtomicUsize,
}

impl SharedPipelineState {
    fn new(num_workers: usize, main_thread_handle: std::thread::Thread) -> Self {
        let data_queue_cap = num_workers * 8;
        let compress_queue_cap = num_workers * 4;

        Self {
            merge_phases: crate::merge_phases::MergePhaseCounters::default(),
            phase2_scans: crate::merge_stalls::Phase2ScanStats::default(),
            wake_latency: crate::merge_stalls::WakeLatencyStats::default(),
            wake_latency_at_merge_start: std::sync::OnceLock::new(),
            epoch: Instant::now(),
            block_lifecycle: crate::merge_trace::BlockLifecycleStats::default(),
            refill: crate::merge_trace::RefillStats::default(),
            consumer_trace: crate::merge_trace::ConsumerTraceStats::default(),
            fruitless_scan: crate::merge_trace::DurationHistogram::default(),
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
            worker_threads: (0..num_workers).map(|_| std::sync::OnceLock::new()).collect(),
            wake_cursor: AtomicUsize::new(0),
            deep_read_batches: AtomicU64::new(0),
            deep_read_blocks: AtomicU64::new(0),
            shallow_read_batches: AtomicU64::new(0),
            shallow_read_blocks: AtomicU64::new(0),
            phase2_lowest_active: AtomicUsize::new(0),
        }
    }

    /// Publish a worker's thread handle so it can be woken by name.
    fn register_worker_thread(&self, worker_id: usize) {
        if let Some(slot) = self.worker_threads.get(worker_id) {
            let _ = slot.set(std::thread::current());
        }
    }

    /// Wake one idle worker, rotating which one.
    ///
    /// The consumer calls this the instant a reorder buffer drains. Without it
    /// the wake path runs one way — workers to main thread — so a file that
    /// runs dry waits out an idle worker's backoff (up to [`MAX_BACKOFF_US`])
    /// before anyone even looks at it. That wait is pure latency: the work is
    /// available and unclaimed.
    ///
    /// One worker, not all: the refill it needs to start is a disk read, which
    /// is serialized by the file's reader mutex anyway, and a woken worker
    /// resets its own backoff to [`MIN_BACKOFF_US`] on success and so stays hot
    /// for the blocks that follow. Waking the whole pool would spend N fruitless
    /// scans to get the one read that matters.
    pub(crate) fn wake_one_worker(&self) {
        if self.worker_threads.is_empty() {
            return;
        }
        let idx = Self::wake_target(
            self.wake_cursor.fetch_add(1, Ordering::Relaxed),
            self.active_worker_limit.load(Ordering::Acquire),
            self.worker_threads.len(),
        );
        if let Some(handle) = self.worker_threads[idx].get() {
            handle.unpark();
        }
    }

    /// Which worker slot the next wake should target.
    ///
    /// Rotates over the *active* workers, not the pool width. A pool sized to a
    /// wide Phase 1 and capped to a narrower Phase 2 leaves workers with
    /// `worker_id >= active_worker_limit` idle by policy, and unparking one of
    /// those buys nothing: it re-enters its idle path without ever looking at
    /// the starving file, while the workers that could refill it stay parked
    /// for the rest of their backoff — exactly the latency this wake exists to
    /// remove.
    ///
    /// Pure so the selection is testable without standing up a pool, following
    /// [`classify_scan`](crate::merge_stalls::classify_scan). Clamped to at
    /// least 1 so the modulo is always defined, however the limit was set.
    fn wake_target(cursor: usize, active_limit: usize, pool_width: usize) -> usize {
        cursor % active_limit.clamp(1, pool_width.max(1))
    }

    /// Advance the drain frontier past every source that is now fully drained.
    ///
    /// Called by the merge consumer when a source reports drained. Walks rather
    /// than jumping to `source_id + 1` so the frontier stays truthful when
    /// sources finish out of order — otherwise a single early finisher would
    /// park the hint past files that are still live, and the hint would point at
    /// the one place there is guaranteed to be no work.
    pub(crate) fn advance_phase2_frontier(&self, files: &[Phase2FileState]) {
        let mut frontier = self.phase2_lowest_active.load(Ordering::Relaxed);
        while frontier < files.len() && files[frontier].is_drained() {
            frontier += 1;
        }
        self.phase2_lowest_active.store(frontier, Ordering::Relaxed);
    }

    /// Nanoseconds since the pool's epoch, the unit every `*_nanos` field uses.
    pub(crate) fn now_nanos(&self) -> u64 {
        crate::merge_trace::elapsed_nanos(self.epoch)
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
    /// The wait the previous loop iteration took, if it waited. Taken (and
    /// cleared) by the next iteration that finds work, which is what makes that
    /// wait "productive" — see [`crate::merge_stalls::WakeLatencyStats`].
    last_wait: Option<PendingWait>,
}

impl SortWorkerState {
    /// Returns true if this worker is holding any items that need advancement.
    /// CRITICAL: Workers must not exit while holding items — they would be lost.
    fn has_any_held_items(&self) -> bool {
        !self.held_raw_input_blocks.is_empty() || self.held_decompressed_input.is_some()
    }
}

/// A wait whose cost is not yet known, and the phase it ran in.
///
/// The duration is the *observed* one rather than the backoff level that was
/// requested: waits are `park_timeout`, so a wake from
/// [`SharedPipelineState::wake_one_worker`] ends one early, and the gap between
/// requested and actual is precisely what that wake exists to create.
///
/// The phase travels with it because the wait and the step that redeems it are
/// separate loop iterations and the phase can change in between — see
/// [`productive_wait_nanos`].
#[derive(Debug, Clone, Copy)]
struct PendingWait {
    /// Pipeline phase the wait ran in (see [`phase`]).
    phase: u8,
    /// Observed duration of the wait, in nanoseconds.
    nanos: u64,
}

/// The wait to charge as productive, for a step that succeeded in `current_phase`.
///
/// `None` when nothing was pending, or when the wait ran in a different phase.
/// [`crate::merge_stalls::WakeLatencyStats`] runs for the pool's whole life and
/// the merge's share is taken by subtracting a snapshot pinned at the Phase 1 ->
/// Phase 2 boundary, so a Phase 1 wait recorded after that snapshot lands inside
/// the merge's window and reports discovery lag for a merge that had not started
/// when the worker went to sleep. Nor can it be credited to Phase 1, whose
/// snapshot is already taken. A wait that straddles the boundary is therefore
/// dropped rather than misattributed to whichever side happens to be reporting.
///
/// Pure so the rule is testable without running the pool through a phase
/// transition, following [`jittered_wait_micros`].
fn productive_wait_nanos(pending: Option<PendingWait>, current_phase: u8) -> Option<u64> {
    pending.filter(|wait| wait.phase == current_phase).map(|wait| wait.nanos)
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
                &[SortStep::Compress, SortStep::DecompressInput, SortStep::ReadInputBlocks]
            } else {
                // Default/starving: feed the main thread first, compress if available
                &[SortStep::DecompressInput, SortStep::ReadInputBlocks, SortStep::Compress]
            }
        }
        phase::PHASE2 => {
            // Phase 2 file work and output compression are independent: each
            // worker grabs whatever has work. We never gate on `all_chunks_eof`
            // — even after disk reads finish, decompression and parser drain
            // continue until all per-file reorder buffers empty.
            if bp.compress_has_items {
                // Drain output compression while we can; it's the writer-side bottleneck.
                &[SortStep::Compress, SortStep::Phase2FileWork]
            } else {
                &[SortStep::Phase2FileWork, SortStep::Compress]
            }
        }
        // Legacy/transition: compress only (drain any remaining jobs, of either
        // kind — each carries its own `CompressTarget`).
        _ => &[SortStep::Compress],
    }
}

// ============================================================================
// Backoff (ported from unified pipeline base.rs)
// ============================================================================

/// Minimum backoff duration in microseconds.
pub(crate) const MIN_BACKOFF_US: u64 = 10;
/// Maximum backoff duration in microseconds (1ms).
///
/// This bounds how late an idle worker notices work that appeared while it was
/// backing off. [`SharedPipelineState::wake_one_worker`] cuts that wait short
/// for the case that measurably mattered — a merge source running dry — and
/// [`crate::merge_stalls`] measures how often the full wait is still paid.
pub(crate) const MAX_BACKOFF_US: u64 = 1000;

/// Wait out the given backoff with ±25% jitter, or until someone unparks us.
///
/// Uses `yield_now()` at minimum backoff to avoid a syscall on the hot path;
/// above that it parks with a timeout rather than sleeping, so a wake from
/// [`SharedPipelineState::wake_one_worker`] returns immediately instead of
/// waiting out the full duration. `worker_id` and `iter` are mixed into the seed
/// so concurrent workers do not produce identical durations and thundering-herd
/// on wakeup.
///
/// A pending unpark token makes the next `park_timeout` return at once. That is
/// a spurious early wake, which costs one scan and is harmless — the caller
/// loops.
fn idle_wait_with_jitter(backoff_us: u64, worker_id: usize, iter: u64) {
    if backoff_us <= MIN_BACKOFF_US {
        std::thread::yield_now();
    } else {
        let actual_us = jittered_wait_micros(backoff_us, worker_id, iter);
        std::thread::park_timeout(std::time::Duration::from_micros(actual_us));
    }
}

/// `backoff_us` spread over ±25%, deterministically per `(worker_id, iter)`.
///
/// Pure so the spread is testable without parking a thread, following
/// [`classify_scan`](crate::merge_stalls::classify_scan).
///
/// The two directions are computed separately because the offset is unsigned.
/// Subtracting the half-range from `seed % (range * 2)` in `u64` — the obvious
/// spelling — saturates every negative offset to zero, which silently turns
/// ±25% into 0..+25%: no worker ever wakes *earlier* than nominal, so the pool
/// keeps exactly the synchronized wake-ups the jitter exists to break up.
///
/// Callers must have already excluded `backoff_us <= MIN_BACKOFF_US`, which is
/// what keeps `jitter_range` non-zero and the modulo defined.
fn jittered_wait_micros(backoff_us: u64, worker_id: usize, iter: u64) -> u64 {
    debug_assert!(backoff_us > MIN_BACKOFF_US, "caller must handle the yield_now case");
    let jitter_range = (backoff_us / 4).max(1);
    // Cheap deterministic seed — no syscall, differs per worker and iteration.
    let jitter_seed = (worker_id as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15).wrapping_add(iter);
    let jitter_offset = jitter_seed % (jitter_range * 2);
    let actual_us = if jitter_offset < jitter_range {
        backoff_us.saturating_sub(jitter_range - jitter_offset)
    } else {
        backoff_us.saturating_add(jitter_offset - jitter_range)
    };
    actual_us.max(MIN_BACKOFF_US)
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
                        last_wait: None,
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
        // Publish this thread so the consumer can wake it out of a backoff.
        shared.register_worker_thread(worker.worker_id);
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
                    idle_wait_with_jitter(worker.backoff_us, worker.worker_id, worker.idle_iter);
                    worker.idle_iter = worker.idle_iter.wrapping_add(1);
                    worker.backoff_us = (worker.backoff_us * 2).min(MAX_BACKOFF_US);
                }
                // A capped worker is idle by policy, not because work was slow
                // to appear, so its wait must not count as discovery latency.
                worker.last_wait = None;
                continue;
            }

            // 2. Check phase completion — wait for next phase, only exit on SHUTDOWN.
            //    Workers must survive across Phase 1 → Phase 2 transitions.
            if Self::is_phase_complete(shared, current_phase) && !worker.has_any_held_items() {
                idle_wait_with_jitter(worker.backoff_us, worker.worker_id, worker.idle_iter);
                worker.idle_iter = worker.idle_iter.wrapping_add(1);
                worker.backoff_us = (worker.backoff_us * 2).min(MAX_BACKOFF_US);
                // Waiting for the next phase, not for work — see above.
                worker.last_wait = None;
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
                // The wait that preceded this step was slept *through* work
                // becoming available, so the work may have arrived at any point
                // during it. That is the part of idle time the pipeline actually
                // pays for, and it is invisible in the idle total, which counts
                // the productive and unproductive waits alike. `wake_one_worker`
                // shortens the wait when the consumer knows work is coming; what
                // is recorded here is what the wait actually cost.
                //
                // `take` unconditionally: a wait from another phase is dropped,
                // not carried forward to the next step in this one.
                if let Some(waited_ns) =
                    productive_wait_nanos(worker.last_wait.take(), current_phase)
                {
                    shared.wake_latency.record_productive_sleep(waited_ns);
                }
                worker.backoff_us = MIN_BACKOFF_US;
            } else {
                let slept_us = worker.backoff_us;
                let idle_start = Instant::now();
                idle_wait_with_jitter(slept_us, worker.worker_id, worker.idle_iter);
                worker.idle_iter = worker.idle_iter.wrapping_add(1);
                worker.backoff_us = (worker.backoff_us * 2).min(MAX_BACKOFF_US);
                let idle_ns = Self::nanos_u64(idle_start.elapsed());
                pstats.record_idle(worker.worker_id, idle_ns);
                shared.wake_latency.record_sleep(slept_us, idle_ns);
                worker.last_wait = Some(PendingWait { phase: current_phase, nanos: idle_ns });
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
            SortStep::Compress => !shared.compress_queue.is_empty(),
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
            SortStep::Compress => Self::try_compress(shared, worker),
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

        // Why each file was passed over. Lives on the stack and is thrown away
        // the moment the scan finds work, so a productive scan pays a handful of
        // increments; only the fall-through at the bottom — the path that is
        // about to sleep anyway — publishes it. See [`crate::merge_stalls`].
        let mut tally = Phase2ScanTally::default();
        let scan_start = shared.now_nanos();

        // Where to start looking. Normally the worker's own round-robin cursor,
        // which spreads N workers over the file set and is what keeps an
        // interleaved merge — where every file is hot — saturated.
        //
        // But when the input was already in the requested order, the spill runs
        // are disjoint and only one file is ever being consumed. Round-robin
        // then spends most of a scan on files sitting at `PHASE2_DECOMP_CAP`
        // with nothing to give, while the one file the merge is blocked on waits
        // to be reached. So if the drain frontier is starving, start there: it
        // is the file the consumer is either blocked on or about to be.
        //
        // Gated on `is_starving` rather than applied unconditionally, because a
        // frontier that is merely *active* is the common case in an interleaved
        // merge, and sending every worker to file 0 there would undo the spread
        // that makes that case fast.
        let frontier = shared.phase2_lowest_active.load(Ordering::Relaxed);
        let frontier_starving = frontier < n && files[frontier].is_starving();
        let start = if frontier_starving { frontier } else { worker.phase2_file_cursor };

        for offset in 0..n {
            let i = (start + offset) % n;
            let file = &files[i];

            // -- Try decompression first (highest-value work) ----------------
            // `try_pop_raw_for_decompress` increments `decomp_in_flight` before
            // returning, so the consumer's `is_drained` check sees this work as
            // outstanding even after the raw FIFO becomes empty. We MUST
            // decrement after inserting into the reorder buffer (or on the
            // error path) to keep the counter balanced.
            // Diverges on success -- `decompress_and_publish` always returns --
            // so the binding is why the decompress half declined.
            let pop_skip = match Self::try_pop_raw_for_decompress(file) {
                Err(skip) => skip,
                Ok(entry) => {
                    worker.phase2_file_cursor = (i + 1) % n;
                    Self::note_claim(shared, file, &entry);
                    return Self::decompress_and_publish(shared, worker, file, i, entry);
                }
            };

            // -- Try reading raw blocks from disk ----------------------------
            // Skip if disk reader is contended OR already at EOF.
            let Ok(mut reader_guard) = file.reader.try_lock() else {
                // Another worker is inside a disk read for this file. That is
                // I/O in flight, not lock contention, and conflating the two
                // would blame the mutex for the disk.
                tally.note(combine_skip(pop_skip, ReadSkip::ReaderBusy));
                continue;
            };
            if reader_guard.eof {
                tally.note(combine_skip(pop_skip, ReadSkip::ReaderEof));
                continue;
            }

            // The file the merge is draining gets a deeper allowance than the
            // uniform per-file one.
            //
            // Uniform read-ahead is what makes the disjoint-run case slow. Only
            // one worker may read a given file (its reader mutex), so at
            // `PHASE2_READ_BATCH` = 4 blocks the hot file is refilled 4 blocks at
            // a time by one thread, and the per-batch read latency is amortized
            // over only those 4. Reading a larger run turns many small
            // serialized reads into one sequential one, without needing more
            // readers.
            //
            // Keyed on being the frontier alone, NOT on the file also looking
            // starved. An earlier version required `is_starving()` here, which
            // needs `decomp_in_flight == 0` -- but by the time a worker reaches
            // the read path some other worker has usually already begun
            // decompressing, so the file no longer looked starved and the deep
            // path fired for only 6% of batches (5.6 blocks per read against the
            // intended 32). Whether a worker happens to be mid-decompress says
            // nothing about how far ahead the file should read.
            //
            // Still scoped to one file: at merge start every file is empty, and
            // letting all K deepen at once would multiply read-ahead memory by K.
            let (raw_cap, read_batch) = if i == frontier {
                (PHASE2_STARVING_RAW_CAP, PHASE2_STARVING_READ_BATCH)
            } else {
                (PHASE2_RAW_CAP, PHASE2_READ_BATCH)
            };

            // Bound disk read-ahead per file: don't keep pulling if the raw
            // FIFO is already full. Use try_lock so a momentarily contended
            // raw FIFO doesn't block the reader.
            let read_skip = match file.raw_blocks.try_lock() {
                Ok(g) if g.len() >= raw_cap => Some(ReadSkip::RawFull),
                Ok(_) => None,
                // Treated as full, as it always has been -- but recorded as
                // contention, because a FIFO whose depth we could not read is
                // not evidence that read-ahead depth is the constraint.
                Err(_) => Some(ReadSkip::RawLockContended),
            };
            if let Some(read_skip) = read_skip {
                tally.note(combine_skip(pop_skip, read_skip));
                continue;
            }

            // Both codecs read into the same `Vec<Vec<u8>>`, so the read is
            // dispatched on the codec but the failure is handled once. The two
            // arms previously carried identical error handling differing only
            // in one word of the message, which had to be kept in step by hand.
            let read = match file.codec {
                SpillCodec::Bgzf => shared
                    .merge_phases
                    .read
                    .time(|| read_raw_blocks(&mut reader_guard.inner, read_batch))
                    .map(|blocks| blocks.into_iter().map(|b| b.data).collect()),
                SpillCodec::Zstd => shared
                    .merge_phases
                    .read
                    .time(|| read_raw_zstd_frames(&mut reader_guard.inner, read_batch)),
            };
            let raw_bytes: Vec<Vec<u8>> = match read {
                Ok(bytes) => bytes,
                Err(e) => {
                    log::error!("I/O error reading {} chunk file (source {i}): {e}", file.codec);
                    shared.chunk_read_error.store(true, Ordering::Release);
                    return Self::retire_phase2_source(shared, worker, file, reader_guard, i, n);
                }
            };

            // Record which allowance this batch was taken at, and what it
            // returned, so "the deep path did not help" can be told apart from
            // "the deep path did not run".
            let (batches, blocks) = if i == frontier {
                (&shared.deep_read_batches, &shared.deep_read_blocks)
            } else {
                (&shared.shallow_read_batches, &shared.shallow_read_blocks)
            };
            batches.fetch_add(1, Ordering::Relaxed);
            blocks.fetch_add(raw_bytes.len() as u64, Ordering::Relaxed);

            if raw_bytes.is_empty() {
                return Self::retire_phase2_source(shared, worker, file, reader_guard, i, n);
            }

            // Acquire the raw_blocks lock BEFORE releasing the reader lock and
            // assigning serials. This is critical for FIFO order: if we dropped
            // the reader lock before pushing, two workers could each read a
            // batch and bump `next_serial`, then race on the raw_blocks push,
            // landing higher serials in front of lower ones. The merge
            // consumer's gap-filler admission rule cannot recover from that
            // and would deadlock. Lock order `reader → raw_blocks` is the only
            // nested-lock path in this function.
            let enqueued_nanos = shared.now_nanos();
            let mut raw_guard = file.raw_blocks.lock().expect("phase2 raw_blocks mutex poisoned");
            let start_serial = reader_guard.next_serial;
            reader_guard.next_serial += raw_bytes.len() as u64;
            for (idx, bytes) in raw_bytes.into_iter().enumerate() {
                raw_guard.push_back(RawEntry {
                    serial: start_serial + idx as u64,
                    bytes,
                    enqueued_nanos,
                });
            }
            file.raw_len.store(raw_guard.len(), Ordering::Relaxed);
            drop(raw_guard);
            drop(reader_guard);

            // If this file was dry when its buffer emptied, the read is what
            // the consumer has actually been waiting for -- the pool could not
            // have decompressed anything sooner. Separating that from the
            // claim/insert lags is what keeps a storage problem from reading as
            // a scheduling one.
            if let Some(emptied) = file.emptied_at()
                && file.emptied_while_dry()
            {
                shared.refill.read_lag.record(enqueued_nanos.saturating_sub(emptied));
            }

            worker.phase2_file_cursor = (i + 1) % n;
            return StepResult::Success;
        }

        // Every file declined. This worker is about to sleep, so publishing the
        // reasons here costs nothing it was not already about to spend, and it
        // turns an entry in the idle total into an explanation of that entry.
        shared.phase2_scans.record_fruitless_scan(tally);
        // A scan touches every spill file, so at 86 files and ~2M fruitless
        // scans this is not obviously negligible -- and it sits directly on the
        // path by which a hungry file gets noticed. Measure it rather than
        // assume.
        shared.fruitless_scan.record(shared.now_nanos().saturating_sub(scan_start));
        StepResult::InputEmpty
    }

    /// Record that a worker claimed `entry`, closing the scheduling half of a
    /// refill cycle if this is the first claim since the file ran dry.
    fn note_claim(shared: &SharedPipelineState, file: &Phase2FileState, entry: &RawEntry) {
        let now = shared.now_nanos();
        shared.block_lifecycle.raw_dwell.record(now.saturating_sub(entry.enqueued_nanos));
        if let Some(emptied) = file.emptied_at()
            && file.take_first_claim_since_empty()
        {
            shared.refill.claim_lag.record(now.saturating_sub(emptied));
        }
    }

    /// Decompress one raw block and publish it into its file's reorder buffer.
    ///
    /// Split out of [`Self::try_phase2_file_work`] so the scan loop reads as the
    /// scan it is. Always returns [`StepResult::Success`]: a decompression
    /// failure is still work done — it sets the error flag and wakes the
    /// consumer, which must surface the failure rather than spin.
    ///
    /// The caller must have advanced `worker.phase2_file_cursor` already, and
    /// must have obtained `entry` from [`Self::try_pop_raw_for_decompress`],
    /// which reserved the matching `decomp_in_flight` slot this function is
    /// responsible for releasing.
    fn decompress_and_publish(
        shared: &SharedPipelineState,
        worker: &mut SortWorkerState,
        file: &Phase2FileState,
        source_idx: usize,
        entry: RawEntry,
    ) -> StepResult {
        let RawEntry { serial, bytes: raw_bytes, .. } = entry;
        let data = match file.codec {
            SpillCodec::Bgzf => {
                let raw_block = RawBgzfBlock { data: raw_bytes };
                match shared
                    .merge_phases
                    .decompress
                    .time(|| decompress_block(&raw_block, &mut worker.decompressor))
                {
                    Ok(d) => d,
                    Err(e) => {
                        log::error!(
                            "BGZF decompression error (chunk source {source_idx} serial {serial}): {e}"
                        );
                        shared.decompression_error.store(true, Ordering::Release);
                        file.decomp_in_flight.fetch_sub(1, Ordering::AcqRel);
                        shared.main_thread_handle.unpark();
                        return StepResult::Success;
                    }
                }
            }
            SpillCodec::Zstd => {
                // Allocate the scratch buffer lazily so BGZF-only sorts don't
                // pay 256 KiB × num_workers of dead memory.
                if worker.zstd_decompress_buf.len() < ZSTD_FRAME_DECOMP_CAP {
                    worker.zstd_decompress_buf.resize(ZSTD_FRAME_DECOMP_CAP, 0);
                }
                match shared.merge_phases.decompress.time(|| {
                    worker
                        .zstd_decompressor
                        .decompress_to_buffer(&raw_bytes, &mut worker.zstd_decompress_buf)
                }) {
                    // Copy the `n` decompressed bytes (≤ one staging-buffer's
                    // worth, typically ~65 KB) into a fresh Vec for the
                    // consumer. The scratch buffer keeps its 256 KiB capacity so
                    // the next frame on this worker reuses it without
                    // reallocating.
                    Ok(n) => worker.zstd_decompress_buf[..n].to_vec(),
                    Err(e) => {
                        log::error!(
                            "zstd decompression error (chunk source {source_idx} serial {serial}): {e}"
                        );
                        shared.decompression_error.store(true, Ordering::Release);
                        file.decomp_in_flight.fetch_sub(1, Ordering::AcqRel);
                        shared.main_thread_handle.unpark();
                        return StepResult::Success;
                    }
                }
            }
        };

        let inserted_nanos = shared.now_nanos();
        let now_poppable = {
            let mut dec_guard =
                file.decompressed.lock().expect("phase2 decompressed mutex poisoned");
            dec_guard.insert(serial, TimedBlock { data, inserted_nanos });
            file.decomp_len.store(dec_guard.len(), Ordering::Relaxed);
            dec_guard.can_pop()
        };
        // Close the refill cycle this insert may have ended. Recorded before
        // the wake below so the latency excludes however long the consumer then
        // takes to notice -- that delay is the consumer's, and it is measured
        // separately as park time.
        if let Some(emptied) = file.emptied_at() {
            shared.refill.insert_lag.record(inserted_nanos.saturating_sub(emptied));
            file.mark_refilled();
        }
        // Decrement AFTER the insert is published. The unpark below wakes the
        // consumer in case it has been parked waiting on this specific file
        // (now_poppable=true) or is in the is_drained path waiting for
        // in_flight to reach zero.
        file.decomp_in_flight.fetch_sub(1, Ordering::AcqRel);
        if now_poppable || file.is_drained() {
            // Wake the consumer either because new data is available or because
            // the last in-flight decompression for this file just completed and
            // the file is now fully drained.
            shared.main_thread_handle.unpark();
        }
        StepResult::Success
    }

    /// Try to pop a raw block from `file` for decompression, applying
    /// deadlock-free admission control against the file's reorder buffer.
    ///
    /// Returns `Ok(RawEntry)` if a block was popped, and otherwise
    /// the [`PopSkip`] saying why not: either lock was contended (`try_lock`
    /// failed), the raw FIFO was empty, or the reorder buffer was at cap and the
    /// head raw block was not a gap-filler. The caller records that reason —
    /// "the pool declined this file" and "the pool had nothing to take" look
    /// identical in an idle total and mean opposite things.
    ///
    /// On success, `decomp_in_flight` is incremented so the consumer's
    /// `is_drained()` check correctly reflects the in-progress decompression.
    /// The caller is responsible for the matching decrement after inserting
    /// (or on the decompression-error path).
    fn try_pop_raw_for_decompress(
        file: &Phase2FileState,
    ) -> std::result::Result<RawEntry, PopSkip> {
        let Ok(mut raw_guard) = file.raw_blocks.try_lock() else {
            return Err(PopSkip::RawLockContended);
        };
        let Some(head_serial) = raw_guard.front().map(|e| e.serial) else {
            return Err(PopSkip::RawEmpty);
        };

        // Cheap admission check using the per-file reorder buffer.
        // Two cases admit: (1) under cap (normal), (2) reorder buffer is stuck
        // and this serial is the gap-filler. Otherwise: backpressure.
        let admit = {
            let Ok(dec_guard) = file.decompressed.try_lock() else {
                return Err(PopSkip::DecompLockContended);
            };
            dec_guard.len() < PHASE2_DECOMP_CAP
                || (!dec_guard.can_pop() && head_serial == dec_guard.next_seq())
        };
        if !admit {
            return Err(PopSkip::DecompCapped);
        }

        // Reserve the in-flight slot under the raw_blocks lock so the consumer
        // can never observe (raw_empty && in_flight==0 && decompressed_empty)
        // while a worker is still in the middle of decompressing this block.
        // The FIFO was non-empty above and is still held, so this cannot fail.
        let popped = raw_guard.pop_front().expect("raw FIFO emptied under its own lock");
        file.raw_len.store(raw_guard.len(), Ordering::Relaxed);
        file.decomp_in_flight.fetch_add(1, Ordering::AcqRel);
        Ok(popped)
    }

    // ========================================================================
    // Compress Step (shared by Phase 1 + Phase 2)
    // ========================================================================

    /// Try to pick up a compress job from the `ArrayQueue` (non-blocking).
    ///
    /// The compressor is selected from the job's own [`CompressTarget`], which
    /// the producing staging buffer stamped on it at submit time. Nothing about
    /// this choice depends on shared mutable state, so a block is compressed at
    /// the level its writer asked for no matter how long it waited in the queue
    /// or which phase the pool is in when a worker gets to it.
    ///
    /// This is the third and final form of that selection. It was originally
    /// read from `shared.phase` at pop time; that was narrowed to the dispatched
    /// `SortStep` to close the window where `set_phase` fired between the pop
    /// and the choice — but a step is still chosen from the phase, so blocks
    /// that outlived a phase transition were still mis-compressed.
    fn try_compress(shared: &SharedPipelineState, worker: &mut SortWorkerState) -> StepResult {
        let Some(job) = shared.compress_queue.pop() else {
            return StepResult::InputEmpty;
        };
        // Destructured so the BGZF compressor and the zstd one are borrowed as
        // disjoint fields rather than through `worker` twice.
        let SortWorkerState { compressor, output_compressor, zstd_compressor, .. } = worker;
        let target = job.target;
        let bgzf_compressor = match job.target {
            CompressTarget::Spill => compressor,
            CompressTarget::Output => output_compressor,
        };
        // Attributed by target, not lumped: this one step compresses both Phase 1
        // spill blocks and merge output blocks, so a single counter would charge
        // spill compression to the merge.
        let counter = match target {
            CompressTarget::Spill => &shared.merge_phases.spill_compress,
            CompressTarget::Output => &shared.merge_phases.output_compress,
        };
        // Only the compression is timed; the handoff below can block on a full
        // result channel, which is writer backpressure rather than CPU work.
        let compressed =
            counter.time(|| Self::compress_job(&job, bgzf_compressor, zstd_compressor));
        Self::deliver_compress_result(shared, job, compressed);
        StepResult::Success
    }

    // ========================================================================
    // Helper: EOF tracking for Phase 2
    // ========================================================================

    /// Retire a phase-2 source: mark its reader EOF, wake the consumer, and
    /// advance this worker's cursor past it.
    ///
    /// Shared by the two ways a source stops producing — a clean end of file
    /// and an unreadable chunk — which differ only in whether the caller
    /// latches `chunk_read_error` before calling. Both must retire the reader
    /// rather than leave it live: a reader still marked readable is retried by
    /// every worker forever, and a consumer parked on this source never wakes.
    ///
    /// Takes the guard by value so the reader lock is provably released before
    /// the consumer is unparked; waking it while still holding the lock would
    /// have it contend for a lock this thread has not dropped yet.
    ///
    /// Returns [`StepResult::Success`] because the *step* completed. A read
    /// failure reaches the consumer through `chunk_read_error`, not through
    /// this result.
    fn retire_phase2_source(
        shared: &SharedPipelineState,
        worker: &mut SortWorkerState,
        file: &Phase2FileState,
        mut reader_guard: MutexGuard<'_, Phase2Reader>,
        source: usize,
        num_sources: usize,
    ) -> StepResult {
        file.mark_reader_eof(&mut reader_guard);
        drop(reader_guard);
        shared.main_thread_handle.unpark();
        Self::maybe_mark_all_eof(shared);
        worker.phase2_file_cursor = (source + 1) % num_sources;
        StepResult::Success
    }

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

    /// Number of workers currently eligible to run work, as last set by
    /// [`Self::set_active_workers`].
    ///
    /// Not the same as [`Self::num_workers`]: the pool is sized to the wider of
    /// the two phases and then capped per phase, so on a run whose Phase 1 is
    /// wider than its Phase 2 the pool holds threads that cannot take Phase 2
    /// work. Anything dividing by a thread count — utilization above all —
    /// wants this number, or it charges the merge for workers that were never
    /// allowed to help and reports a saturated pool as idle.
    pub fn active_workers(&self) -> usize {
        self.shared.active_worker_limit.load(Ordering::Acquire)
    }

    /// Busy-time breakdown of worker-side merge work.
    ///
    /// Overlapping sums, not a partition of merge wall time -- see
    /// [`crate::merge_phases`].
    pub(crate) fn merge_phase_breakdown(&self) -> crate::merge_phases::MergePhaseBreakdown {
        self.shared.merge_phases.snapshot()
    }

    /// The pool's shared state, for the merge consumer's own instrumentation.
    ///
    /// The consumer runs on the main thread and is not a pool worker, but it is
    /// one end of the same pipeline: it needs the shared clock so its
    /// timestamps are comparable with the workers', and the trace counters so
    /// both halves of a handoff are recorded against one another rather than
    /// separately.
    pub(crate) fn shared_state(&self) -> Arc<SharedPipelineState> {
        Arc::clone(&self.shared)
    }

    /// A spill block's full journey through the merge -- see
    /// [`crate::merge_trace`].
    ///
    /// The two working stages come from the merge phase counters the workers
    /// already record into, and the two dwell stages from `block_lifecycle`.
    /// This is the one place that holds both, which is why it is the only place
    /// a complete report can be assembled.
    pub(crate) fn block_lifecycle_report(&self) -> crate::merge_trace::BlockLifecycleReport {
        let phases = &self.shared.merge_phases;
        crate::merge_trace::BlockLifecycleReport {
            read_batch: phases.read.histogram(),
            decompress: phases.decompress.histogram(),
            ..self.shared.block_lifecycle.snapshot()
        }
    }

    /// How long files that ran dry took to be refilled.
    pub(crate) fn refill_report(&self) -> crate::merge_trace::RefillReport {
        self.shared.refill.snapshot()
    }

    /// The merge consumer's access pattern and what it cost.
    pub(crate) fn consumer_trace_report(&self) -> crate::merge_trace::ConsumerTraceReport {
        self.shared.consumer_trace.snapshot()
    }

    /// How long a worker's fruitless scan of every spill file takes.
    pub(crate) fn fruitless_scan_report(&self) -> crate::merge_trace::HistogramReport {
        self.shared.fruitless_scan.snapshot()
    }

    /// Why Phase 2 file scans found no work -- see [`crate::merge_stalls`].
    pub(crate) fn phase2_scan_report(&self) -> Phase2ScanReport {
        self.shared.phase2_scans.snapshot()
    }

    /// How deep workers were sleeping when they found work, during the merge.
    ///
    /// Scoped to Phase 2 by subtracting the snapshot taken at the phase
    /// boundary; falls back to the running total only when the pool never
    /// entered Phase 2, in which case there is no merge to misattribute to.
    pub(crate) fn wake_latency_report(&self) -> WakeLatencyReport {
        let now = self.shared.wake_latency.snapshot();
        match self.shared.wake_latency_at_merge_start.get() {
            Some(&baseline) => now.since(baseline),
            None => now,
        }
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

    /// Phase 2 read batches split by allowance:
    /// `(frontier_batches, frontier_blocks, other_batches, other_blocks)`.
    pub(crate) fn read_batch_split(&self) -> (u64, u64, u64, u64) {
        (
            self.shared.deep_read_batches.load(Ordering::Relaxed),
            self.shared.deep_read_blocks.load(Ordering::Relaxed),
            self.shared.shallow_read_batches.load(Ordering::Relaxed),
            self.shared.shallow_read_blocks.load(Ordering::Relaxed),
        )
    }

    /// Set the current pipeline phase.
    pub fn set_phase(&self, new_phase: u8) {
        // Pin the wake-latency counters as they stand entering the merge. They
        // run for the pool's whole life, so without a baseline to subtract, a
        // merge-scoped report would include every Phase 1 sleep -- see
        // `WakeLatencyReport::since`. `OnceLock` so a workflow that re-enters
        // PHASE2 keeps the first boundary rather than silently rebasing.
        if new_phase == phase::PHASE2 {
            let _ =
                self.shared.wake_latency_at_merge_start.set(self.shared.wake_latency.snapshot());
        }
        self.shared.phase.store(new_phase, Ordering::Release);
    }

    /// The current pipeline phase (see [`phase`]), the read counterpart to
    /// [`set_phase`](Self::set_phase).
    ///
    /// The phase orders the steps workers schedule, and nothing else — notably
    /// *not* which compressor a block goes through, which `try_compress` takes
    /// from the job's own [`CompressTarget`]. It is the observable the Phase 2
    /// teardown tests assert on. Workers read the phase constantly, but no
    /// caller *outside* the pool needs it, so this accessor exists for those
    /// tests.
    #[cfg(test)]
    pub(crate) fn current_phase(&self) -> u8 {
        self.shared.phase.load(Ordering::Acquire)
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
    /// Both halves belong to the same transition and must happen together: the
    /// phase is what makes workers schedule [`SortStep::Phase2FileWork`], and
    /// widening is what gives them the threads to do it on.
    ///
    /// Call this once ingest has finished — after that point every remaining byte
    /// the pool touches is merge input or output.
    ///
    /// The phase deliberately says nothing about *compression levels*. Each
    /// queued block carries its own [`CompressTarget`], so spill jobs still
    /// outstanding when the phase flips are compressed at `temp_compression`
    /// regardless, and output blocks are compressed at `output_compression` even
    /// if the pool has already left Phase 2. Callers used to owe this method a
    /// drained spill queue for exactly that reason; they no longer do.
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
        // The drain frontier indexes the file vector being replaced, so a pool
        // that merges twice would carry the prior set's frontier into the new
        // one: at or past the old source count it disables frontier
        // prioritization outright, and below it points at whichever file now
        // happens to sit at that index. `advance_phase2_frontier` only walks
        // forward, so neither is self-correcting.
        self.shared.phase2_lowest_active.store(0, Ordering::Release);

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

    /// Compress one job's payload on a worker thread.
    ///
    /// Split from [`Self::deliver_compress_result`] so the merge-phase counters
    /// can time compression alone. Delivery spins in a `try_send` retry loop
    /// whenever the writer is behind, and timing the two together charges that
    /// backpressure to `spill_compress` / `output_compress` — inflating worker
    /// busy time, and with it the utilization the merge verdict reads, so an
    /// I/O-bound merge reports as CPU-bound. That distinction is the entire
    /// point of the breakdown, so the wait is measured nowhere rather than in
    /// the wrong place.
    fn compress_job(
        job: &CompressJob,
        bgzf_compressor: &mut InlineBgzfCompressor,
        zstd_compressor: &mut ZstdCompressor<'static>,
    ) -> Vec<u8> {
        match job.codec {
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
                // Unlike BGZF, there is one zstd compressor per worker, fixed at
                // `temp_compression` — zstd is a spill-only codec (the output BAM
                // is always BGZF, so `PooledBamWriter` always submits
                // `SpillCodec::Bgzf`). An Output job reaching here would be
                // compressed at the temp level, which is the silent wrong-level
                // output this whole selection mechanism exists to prevent — so
                // assert unconditionally rather than only in debug builds. The
                // check is one compare on a `Copy` enum per spill block, inside
                // the zstd arm, so the BGZF path never evaluates it.
                assert_eq!(
                    job.target,
                    CompressTarget::Spill,
                    "zstd is spill-only; its per-worker compressor is fixed at temp_compression"
                );
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
        }
    }

    /// Hand a finished compress job to the I/O writer, recycling its input
    /// buffer.
    ///
    /// Deliberately untimed — see [`Self::compress_job`].
    fn deliver_compress_result(
        shared: &SharedPipelineState,
        job: CompressJob,
        compressed: Vec<u8>,
    ) {
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

    /// A job's [`CompressTarget`] decides its compression level, in every phase.
    ///
    /// This is the invariant that replaced deriving the compressor from the
    /// pool's phase at pop time. That derivation produced the same bug twice:
    /// output blocks popped outside Phase 2 were written at `temp_compression`,
    /// and Phase 1 spill blocks still queued when `begin_phase2` fired were
    /// written at `output_compression`. Both are impossible if the level travels
    /// with the block, so the phase is swept across all three values here and
    /// must not change either result.
    ///
    /// `temp_compression = 0` and `output_compression = 9` make the two levels
    /// separable by size on compressible data: a `Spill` job must stay at the
    /// (larger) stored size and an `Output` job must come back compressed.
    #[rstest]
    #[case::legacy(phase::LEGACY)]
    #[case::phase1(phase::PHASE1)]
    #[case::phase2(phase::PHASE2)]
    fn test_compress_target_decides_level_regardless_of_phase(#[case] phase_under_test: u8) {
        // Compressible: 8 KiB of one byte shrinks to almost nothing at level 9.
        let data = vec![b'A'; 8192];
        let pool = SortWorkerPool::new(2, 0, 9, crate::codec::SpillCodec::Bgzf);
        pool.set_phase(phase_under_test);

        let compressed_len = |target: CompressTarget| {
            let (result_tx, result_rx) = pool.compress_result_channel();
            pool.submit_compress(CompressJob {
                data: data.clone(),
                serial: 0,
                result_tx,
                codec: SpillCodec::Bgzf,
                target,
            });
            result_rx.recv().expect("compress result").compressed.len()
        };

        let spill_len = compressed_len(CompressTarget::Spill);
        let output_len = compressed_len(CompressTarget::Output);

        assert!(
            spill_len > data.len(),
            "a Spill job must be compressed at temp_compression=0 (stored), so it must not shrink \
             below the {} input bytes; got {spill_len} in phase {phase_under_test}",
            data.len()
        );
        assert!(
            output_len < data.len() / 4,
            "an Output job must be compressed at output_compression=9, so 8 KiB of one repeated \
             byte must shrink far below a quarter of its size; got {output_len} in phase \
             {phase_under_test}"
        );

        pool.shutdown();
    }

    /// zstd is spill-only, and `compress_job` asserts that unconditionally.
    ///
    /// The assert is the backstop for the one target/codec combination the
    /// per-job `CompressTarget` cannot express correctly: there is a single zstd
    /// compressor per worker, fixed at `temp_compression`, so an `Output` job
    /// arriving on the zstd arm would be written at the *spill* level — exactly
    /// the silent wrong-level output this mechanism exists to prevent. Without a
    /// test, a future output-zstd feature (or a mis-stamped call site) would only
    /// trip it in production.
    ///
    /// `compress_job` is called directly, on the test thread, rather than
    /// through `submit_compress`: a panic on a pool worker unwinds only that
    /// worker — its guard publishes `worker_panicked` and the unwind stops at the
    /// thread boundary — so `#[should_panic]` would never observe it, and the
    /// test would need a sleep just to reach the panic. Calling the handler
    /// inline puts the assert on the test thread and makes it deterministic.
    #[test]
    #[should_panic(expected = "zstd is spill-only")]
    fn test_zstd_output_job_trips_the_spill_only_assert() {
        let mut bgzf_compressor = InlineBgzfCompressor::new(0);
        let mut zstd_compressor = ZstdCompressor::new(1).expect("zstd compressor init");
        let (result_tx, _result_rx) = bounded(1);

        SortWorkerPool::compress_job(
            &CompressJob {
                data: vec![b'A'; 64],
                serial: 0,
                result_tx,
                codec: SpillCodec::Zstd,
                target: CompressTarget::Output,
            },
            &mut bgzf_compressor,
            &mut zstd_compressor,
        );
    }

    #[test]
    fn test_pool_compress_roundtrip() {
        let pool = SortWorkerPool::new(2, 1, 6, crate::codec::SpillCodec::Bgzf);
        let (result_tx, result_rx) = pool.compress_result_channel();

        // Submit a compress job
        let data = vec![b'A'; 1000];
        pool.submit_compress(CompressJob {
            data,
            serial: 0,
            result_tx,
            codec: SpillCodec::Bgzf,
            target: CompressTarget::Spill,
        });

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
                    target: CompressTarget::Spill,
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
        // A fresh pool is fully active, so anything dividing by the active count
        // before a cap is applied still gets the pool width.
        assert_eq!(pool.active_workers(), 4, "a pool starts with every worker active");
        // Clamp to [1, num_workers].
        pool.set_active_workers(0);
        assert_eq!(pool.active_workers(), 1);
        pool.set_active_workers(100);
        assert_eq!(pool.active_workers(), 4);
        pool.set_active_workers(2);
        assert_eq!(pool.active_workers(), 2);
        pool.shutdown();
    }

    /// The merge breakdown divides worker busy time by the *active* count, so a
    /// pool sized to a wide Phase 1 and capped to a narrower Phase 2 must report
    /// the narrower number — reporting the pool width there would understate
    /// utilization and flip the merge verdict from CPU-bound to I/O-bound.
    #[test]
    fn active_workers_reports_the_phase2_cap_not_the_pool_width() {
        let pool = SortWorkerPool::new(8, 1, 6, crate::codec::SpillCodec::Bgzf);
        pool.set_active_workers(8);
        pool.begin_phase2(3);
        assert_eq!(pool.num_workers(), 8, "the pool is still eight threads wide");
        assert_eq!(pool.active_workers(), 3, "but only three may take Phase 2 work");
        pool.shutdown();
    }

    /// A wait that straddles the Phase 1 -> Phase 2 boundary belongs to
    /// neither report.
    ///
    /// `wake_latency` runs for the pool's whole life and the merge's share is
    /// taken by subtracting a snapshot pinned at the boundary, so a Phase 1 wait
    /// recorded after that snapshot lands inside the merge's window — booking
    /// discovery lag against a merge that had not started when the worker went
    /// to sleep. The wait cannot be credited to Phase 1 either, since its
    /// snapshot is already taken, so it is dropped.
    #[rstest]
    #[case::same_phase_wait_is_charged(
        Some(PendingWait { phase: phase::PHASE2, nanos: 5_000 }), phase::PHASE2, Some(5_000)
    )]
    #[case::phase1_wait_is_not_charged_to_phase2(
        Some(PendingWait { phase: phase::PHASE1, nanos: 5_000 }), phase::PHASE2, None
    )]
    #[case::phase2_wait_is_not_charged_to_phase1(
        Some(PendingWait { phase: phase::PHASE2, nanos: 5_000 }), phase::PHASE1, None
    )]
    #[case::a_phase1_wait_within_phase1_is_charged(
        Some(PendingWait { phase: phase::PHASE1, nanos: 7_000 }), phase::PHASE1, Some(7_000)
    )]
    #[case::no_pending_wait(None, phase::PHASE2, None)]
    fn productive_wait_is_scoped_to_the_phase_it_ran_in(
        #[case] pending: Option<PendingWait>,
        #[case] current_phase: u8,
        #[case] expected: Option<u64>,
    ) {
        assert_eq!(productive_wait_nanos(pending, current_phase), expected);
    }

    /// The jitter must spread waits in *both* directions.
    ///
    /// Its whole job is to desynchronize workers that backed off together, and a
    /// one-sided spread cannot: if every wait is nominal-or-longer, the pool
    /// keeps the alignment the jitter was added to break. The unsigned-subtract
    /// spelling that produces exactly that failure is the one this pins against.
    #[rstest]
    fn jitter_spreads_waits_both_sides_of_the_backoff(
        #[values(20u64, 40, 80, 160, 320, 640, 1000)] backoff_us: u64,
    ) {
        let range = backoff_us / 4;
        let waits: Vec<u64> =
            (0..512).map(|iter| jittered_wait_micros(backoff_us, 3, iter)).collect();

        let min = *waits.iter().min().expect("non-empty");
        let max = *waits.iter().max().expect("non-empty");
        assert!(
            min >= backoff_us.saturating_sub(range).max(MIN_BACKOFF_US),
            "waits must stay within -25%: min {min} for backoff {backoff_us}"
        );
        assert!(max <= backoff_us + range, "waits must stay within +25%: max {max}");
        assert!(
            min < backoff_us,
            "some wait must be SHORTER than nominal (backoff {backoff_us}, min {min}) -- \
             a saturating unsigned subtract clamps these away and leaves 0..+25%"
        );
        assert!(max > backoff_us, "and some wait must be longer (max {max})");
    }

    /// The consumer wakes a worker the instant a merge source runs dry, so the
    /// woken worker has to be one that may actually take Phase 2 work. Rotating
    /// over the pool width instead sends most wakes to workers idled by the
    /// Phase 2 cap, which re-park without looking at the starving file — and
    /// the workers that could have refilled it wait out their full backoff.
    #[rstest]
    // A pool 8 wide capped to 3 must never select a worker above the cap, at
    // any point in the rotation.
    #[case::capped_pool_rotates_only_over_active_workers(0, 3, 8, 0)]
    #[case::capped_pool_wraps_at_the_cap(2, 3, 8, 2)]
    #[case::capped_pool_never_selects_a_capped_worker(3, 3, 8, 0)]
    #[case::capped_pool_keeps_rotating_past_the_wrap(7, 3, 8, 1)]
    // Uncapped, the cap equals the width and the rotation is unchanged.
    #[case::uncapped_pool_rotates_over_the_whole_width(7, 8, 8, 7)]
    #[case::uncapped_pool_wraps_at_the_width(8, 8, 8, 0)]
    // A limit above the width cannot select a slot that does not exist.
    #[case::limit_above_the_width_is_clamped(5, 99, 4, 1)]
    // A zero limit must not divide by zero; slot 0 is always a real worker.
    #[case::zero_limit_falls_back_to_one_worker(7, 0, 8, 0)]
    fn wake_target_rotates_only_over_workers_that_can_take_phase2_work(
        #[case] cursor: usize,
        #[case] active_limit: usize,
        #[case] pool_width: usize,
        #[case] expected: usize,
    ) {
        let idx = SharedPipelineState::wake_target(cursor, active_limit, pool_width);
        assert_eq!(idx, expected);
        assert!(idx < pool_width, "the slot must exist in the pool");
        if active_limit > 0 {
            assert!(idx < active_limit, "a capped worker cannot refill the starving source");
        }
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
        let cs = SortStep::Compress as usize;

        // Run one batch on `pool` (moved in, returned out so the next batch can
        // reuse it) and report the cumulative per-worker Compress counts.
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
                        target: CompressTarget::Spill,
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
            target: CompressTarget::Spill,
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
        stats.record_step(0, SortStep::Compress, 300_000);
        stats.record_idle(0, 100_000);
        stats.record_idle(1, 200_000);

        let read_idx = SortStep::ReadInputBlocks as usize;
        let decomp_idx = SortStep::DecompressInput as usize;
        let compress_idx = SortStep::Compress as usize;

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
        stats.record_step(1, SortStep::Compress, 500_000_000);
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
        assert_eq!(priorities[0], SortStep::Compress);
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
        assert_eq!(priorities[0], SortStep::Compress);
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
        assert_eq!(priorities[0], SortStep::Compress);
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

    /// Build a tiny placeholder raw entry whose contents we never decode.
    fn raw_entry(serial: u64, byte: u8) -> RawEntry {
        RawEntry { serial, bytes: vec![byte; 8], enqueued_nanos: 0 }
    }

    /// A placeholder decompressed block; the insert timestamp is irrelevant to
    /// admission control, which is what these tests exercise.
    fn timed_block(len: usize) -> TimedBlock {
        TimedBlock { data: vec![0u8; len], inserted_nanos: 0 }
    }

    #[test]
    fn test_admission_under_cap_admits() {
        let file = empty_phase2_file();
        file.raw_blocks.lock().expect("raw lock").push_back(raw_entry(0, 0));
        let popped = SortWorkerPool::try_pop_raw_for_decompress(&file);
        assert!(popped.is_ok(), "under cap with empty reorder buffer should admit");
        assert_eq!(file.decomp_in_flight.load(Ordering::Acquire), 1);
        assert!(file.raw_blocks.lock().expect("raw lock").is_empty());
    }

    /// Seed the FIFO through the deque *and* the mirror, the way the read path
    /// leaves them, then pop and assert both moved together.
    ///
    /// Seeding matters: the other admission tests push straight into the deque
    /// and leave `raw_len` at zero, so asserting the mirror is zero after a pop
    /// would pass even with the store deleted. Starting from a non-zero mirror is
    /// what gives the assertion teeth.
    ///
    /// The mirror is what `EmptyCause::classify` and the consumer's park census
    /// read instead of taking the lock, so a dropped store does not corrupt data
    /// -- it silently reclassifies every `Dry` empty as `RawReady`, sending the
    /// next investigation after a scheduling bug that is not there.
    #[rstest]
    #[case::pop_leaves_the_fifo_empty(1)]
    #[case::pop_leaves_the_rest_behind(3)]
    fn test_raw_len_mirror_tracks_the_deque_across_a_pop(#[case] seeded: usize) {
        let file = empty_phase2_file();
        {
            let mut raw = file.raw_blocks.lock().expect("raw lock");
            for serial in 0..seeded as u64 {
                raw.push_back(raw_entry(serial, 0));
            }
            file.raw_len.store(raw.len(), Ordering::Relaxed);
        }
        assert_eq!(file.depths().0, seeded, "the mirror starts where the read path leaves it");

        SortWorkerPool::try_pop_raw_for_decompress(&file).expect("under cap should admit");

        let deque_len = file.raw_blocks.lock().expect("raw lock").len();
        assert_eq!(deque_len, seeded - 1, "exactly one entry leaves the FIFO");
        assert_eq!(
            file.depths().0,
            deque_len,
            "the lock-free raw_len mirror must track the deque; EmptyCause::classify reads it"
        );
    }

    #[test]
    fn test_admission_at_cap_poppable_rejects() {
        let file = empty_phase2_file();
        // Fill the reorder buffer to PHASE2_DECOMP_CAP starting at serial 0,
        // so the buffer is poppable (next_seq = 0 is present) AND at cap.
        {
            let mut dec = file.decompressed.lock().expect("dec lock");
            for s in 0..PHASE2_DECOMP_CAP as u64 {
                dec.insert(s, timed_block(4));
            }
            assert_eq!(dec.len(), PHASE2_DECOMP_CAP);
            assert!(dec.can_pop());
        }
        // Head raw is the next serial after the buffer's contents — not the
        // gap-filler, so admission should reject (consumer is supposed to
        // drain the buffer first).
        file.raw_blocks.lock().expect("raw lock").push_back(raw_entry(PHASE2_DECOMP_CAP as u64, 1));
        let popped = SortWorkerPool::try_pop_raw_for_decompress(&file);
        assert_eq!(
            popped.err(),
            Some(PopSkip::DecompCapped),
            "rejecting at cap must be reported as backpressure, not as an empty file"
        );
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
                dec.insert(s, timed_block(4));
            }
            assert_eq!(dec.len(), PHASE2_DECOMP_CAP);
            assert!(!dec.can_pop(), "buffer should be stuck waiting for serial 0");
        }
        // The head raw is serial 0 — the gap-filler. Admission must take it
        // even though we're at cap, otherwise the consumer deadlocks.
        file.raw_blocks.lock().expect("raw lock").push_back(raw_entry(0, 0));
        let popped = SortWorkerPool::try_pop_raw_for_decompress(&file);
        assert!(popped.is_ok(), "at cap and stuck should admit gap-filler at next_seq");
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
                dec.insert(s, timed_block(4));
            }
        }
        file.raw_blocks
            .lock()
            .expect("raw lock")
            .push_back(raw_entry(PHASE2_DECOMP_CAP as u64 + 1, 2));
        let popped = SortWorkerPool::try_pop_raw_for_decompress(&file);
        assert_eq!(popped.err(), Some(PopSkip::DecompCapped));
        assert_eq!(file.decomp_in_flight.load(Ordering::Acquire), 0);
    }

    /// `read_lag` is defined over the empties where the file held nothing
    /// anywhere, so the cause must be readable back off the file. Recording it
    /// for a `RawReady` or `Decompressing` empty would time a read that was
    /// never on the critical path and make a scheduling delay read as a storage
    /// one.
    #[rstest]
    #[case::nothing_anywhere(0, 0, true)]
    #[case::raw_block_already_queued(1, 0, false)]
    #[case::decompression_already_in_flight(0, 1, false)]
    #[case::both(1, 1, false)]
    fn test_emptied_while_dry_tracks_the_cause_the_cycle_opened_with(
        #[case] raw_len: usize,
        #[case] in_flight: usize,
        #[case] expected_dry: bool,
    ) {
        let file = empty_phase2_file();
        file.mark_emptied(42, crate::merge_trace::EmptyCause::classify(raw_len, in_flight));
        assert_eq!(file.emptied_at(), Some(42));
        assert_eq!(file.emptied_while_dry(), expected_dry);
    }

    /// `claim_lag` is one sample per refill cycle, not one per worker that
    /// piles onto the same hungry file. Only the flag's take-once behaviour
    /// enforces that; a plain load here would inflate `claim_lag`'s count and
    /// shift `claim_share`, the most actionable number in the trace report.
    #[test]
    fn test_first_claim_since_empty_fires_once_per_refill_cycle() {
        let file = empty_phase2_file();
        file.mark_emptied(10, crate::merge_trace::EmptyCause::Dry);
        assert!(file.take_first_claim_since_empty(), "the first responder records");
        assert!(!file.take_first_claim_since_empty(), "a peer on the same cycle does not");

        file.mark_refilled();
        assert_eq!(file.emptied_at(), None, "the cycle is closed");

        file.mark_emptied(20, crate::merge_trace::EmptyCause::RawReady);
        assert!(file.take_first_claim_since_empty(), "a new cycle re-arms the flag");
        assert!(!file.take_first_claim_since_empty(), "and only re-arms it once");
    }

    #[test]
    fn test_admission_empty_raw_reports_raw_empty() {
        let file = empty_phase2_file();
        // Empty raw FIFO — try_pop must report RawEmpty without touching
        // in_flight.
        let popped = SortWorkerPool::try_pop_raw_for_decompress(&file);
        assert_eq!(
            popped.err(),
            Some(PopSkip::RawEmpty),
            "an empty FIFO must not be conflated with a file the pool declined"
        );
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
        file.raw_blocks.lock().expect("raw lock").push_back(raw_entry(0, 0));
        assert!(!file.is_drained(), "raw blocks pending must keep is_drained=false");
    }

    #[test]
    fn test_is_drained_blocks_on_pending_decompressed() {
        let file = empty_phase2_file();
        file.mark_reader_eof(&mut file.reader.lock().expect("reader lock"));
        file.decompressed.lock().expect("dec lock").insert(0, timed_block(3));
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

    /// A pool that merges twice must start the second merge's frontier at its
    /// own first source.
    ///
    /// The frontier is an index into the file vector `set_phase2_files`
    /// replaces, so carrying it across would leave the second merge pointing at
    /// a file the number no longer describes — and `advance_phase2_frontier`
    /// only walks forward, so nothing would ever pull it back.
    #[test]
    fn test_set_phase2_files_restarts_the_drain_frontier() {
        use std::io::Write;
        let dir = tempfile::tempdir().expect("tempdir");
        let spill = |name: &str| {
            let path = dir.path().join(name);
            let mut file = std::fs::File::create(&path).expect("create");
            file.write_all(&[0x1f, 0x8b, 0x00, 0x00]).expect("write magic");
            path
        };

        let pool = SortWorkerPool::new(1, 1, 6, SpillCodec::Bgzf);

        // First merge: one source, drained to completion.
        pool.set_phase2_files(std::slice::from_ref(&spill("first.spill")))
            .expect("set_phase2_files");
        let first = pool.phase2_files();
        first[0].reader_eof.store(true, Ordering::Release);
        assert!(first[0].is_drained(), "the lone source must read as fully consumed");
        pool.shared.advance_phase2_frontier(&first);
        assert_eq!(
            pool.shared.phase2_lowest_active.load(Ordering::Relaxed),
            1,
            "a completed merge parks the frontier past its last source"
        );

        // Second merge: a fresh set, whose file 0 is live.
        pool.set_phase2_files(&[spill("second-a.spill"), spill("second-b.spill")])
            .expect("set_phase2_files");
        assert_eq!(
            pool.shared.phase2_lowest_active.load(Ordering::Relaxed),
            0,
            "the new source set must be prioritized from its own first file"
        );

        pool.shutdown();
    }

    /// Repeated wakes must rotate rather than always hitting worker 0:
    /// re-waking one worker while the others sleep is the starvation the wake
    /// path exists to avoid. Registered handles are all this thread, so the
    /// observable contract here is the cursor advancing once per wake.
    #[test]
    fn test_wake_one_worker_rotates_across_the_pool() {
        let shared = SharedPipelineState::new(4, std::thread::current());
        for expected in 1..=8 {
            shared.wake_one_worker();
            assert_eq!(
                shared.wake_cursor.load(Ordering::Relaxed),
                expected,
                "each wake must advance the cursor exactly once"
            );
        }
    }

    /// The rotation test above only proves the cursor moves; every slot it
    /// rotates through is empty, so it passes with the `unpark` dispatch
    /// deleted. This one registers a real worker and asserts the wake is
    /// actually delivered, which is the behaviour the starving-source path
    /// depends on.
    ///
    /// No sleep is needed to order the two threads: `unpark` on a thread that
    /// has not parked yet leaves a token, so the next `park` returns
    /// immediately. Failure surfaces as the `recv_timeout` below rather than as
    /// a hang.
    #[test]
    fn test_wake_one_worker_unparks_a_registered_worker() {
        let shared = Arc::new(SharedPipelineState::new(1, std::thread::current()));
        let (registered_tx, registered_rx) = bounded::<()>(1);
        let (woken_tx, woken_rx) = bounded::<()>(1);

        let worker_shared = Arc::clone(&shared);
        let worker = std::thread::spawn(move || {
            worker_shared.register_worker_thread(0);
            registered_tx.send(()).expect("main thread is still waiting");
            std::thread::park();
            woken_tx.send(()).expect("main thread is still waiting");
        });

        registered_rx.recv().expect("the worker must publish its handle");
        shared.wake_one_worker();
        woken_rx
            .recv_timeout(Duration::from_secs(10))
            .expect("wake_one_worker must unpark the worker registered in slot 0");
        worker.join().expect("worker thread must not panic");
    }
}
