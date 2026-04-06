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
//!   Phase 2: `DecompressChunks` > `ReadChunkBlocks` > `CompressOutput`.
//! - **Per-worker state**: Each worker owns an `InlineBgzfCompressor` and a
//!   `libdeflater::Decompressor`, avoiding cross-thread synchronization.
//! - **Worker-owns-files**: In Phase 2, each worker exclusively owns a subset of spill
//!   files (worker i owns files i, i+N, i+2N...). No locks needed for file access.
//! - **Held-item pattern**: Workers never block on queue push. If output is full, they
//!   hold the result locally and advance it first on the next iteration.
//! - **Buffer recycling**: A shared buffer pool (`crossbeam` channel) recycles
//!   `Vec<u8>` buffers to avoid per-job heap allocations.
//! - **Bounded backpressure**: All channels are bounded to prevent unbounded memory
//!   growth when producers outpace consumers.

use crossbeam_channel::{Receiver, Sender, bounded};
use crossbeam_queue::ArrayQueue;
use fgumi_bgzf::reader::read_raw_blocks;
use fgumi_bgzf::writer::InlineBgzfCompressor;
use fgumi_bgzf::{RawBgzfBlock, decompress_block};
use log::info;
use std::fmt::Write as FmtWrite;
use std::io::{BufReader, Read};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU8, AtomicU64, Ordering};
use std::thread::{self, JoinHandle};
use std::time::Instant;

// ============================================================================
// Job and Result Types
// ============================================================================

/// A compression job: compress uncompressed data into a BGZF block.
pub struct CompressJob {
    /// Uncompressed data to compress.
    pub data: Vec<u8>,
    /// Serial number for ordering output blocks.
    pub serial: u64,
    /// Channel to send the compressed result back.
    pub result_tx: Sender<CompressResult>,
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
    pub(crate) compress_jobs: AtomicU64,
}

impl PoolStats {
    pub fn log_summary(&self) {
        let compress = self.compress_jobs.load(Ordering::Relaxed);
        info!("  Pool stats: {compress} compress jobs");
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
    /// Read raw BGZF blocks from spill files during Phase 2.
    ReadChunkBlocks = 3,
    /// Decompress spill file blocks during Phase 2.
    DecompressChunks = 4,
    /// Compress data for merge output during Phase 2.
    CompressOutput = 5,
}

impl SortStep {
    /// Number of distinct sort steps.
    pub const COUNT: usize = 6;

    /// Short label for display.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::ReadInputBlocks => "RdInp",
            Self::DecompressInput => "DecInp",
            Self::CompressSpill => "CmpSpl",
            Self::ReadChunkBlocks => "RdChk",
            Self::DecompressChunks => "DecChk",
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
            SortStep::ReadChunkBlocks,
            SortStep::DecompressChunks,
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
            info!("{line}");
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
// Phase 2 Data Types (chunk reading)
// ============================================================================

/// Number of raw BGZF blocks to read in each I/O batch per file.
const CHUNK_READ_BATCH_SIZE: usize = 16;

/// Number of raw BGZF blocks to read in each I/O batch from input file.
const INPUT_READ_BATCH_SIZE: usize = 16;

/// A raw BGZF block tagged with its source (spill file) ID and serial number.
pub struct TaggedRawBlock {
    /// The raw compressed BGZF block.
    pub block: RawBgzfBlock,
    /// Which spill file this block came from.
    pub source_id: usize,
    /// Per-source serial number for reordering after decompression.
    pub serial: u64,
}

/// A decompressed block tagged with its source ID and serial number.
pub struct TaggedDecompressedBlock {
    /// Decompressed data.
    pub data: Vec<u8>,
    /// Which spill file this block came from.
    pub source_id: usize,
    /// Per-source serial number for reordering.
    pub serial: u64,
}

/// Per-spill-file state for the worker-owns-files pattern.
///
/// Each worker exclusively owns a set of `ChunkFileState` instances.
/// Worker `i` owns files `i, i+N, i+2N, ...` where `N` is the number of workers.
/// No locks needed — each worker has exclusive access to its files.
pub struct ChunkFileState {
    /// Buffered reader for this spill file.
    reader: BufReader<std::fs::File>,
    /// Source ID for tagging blocks.
    source_id: usize,
    /// Next serial number to assign to blocks from this file.
    next_serial: u64,
    /// Whether this file has reached EOF.
    eof: bool,
}

impl ChunkFileState {
    /// Create a new chunk file state.
    fn new(reader: BufReader<std::fs::File>, source_id: usize) -> Self {
        Self { reader, source_id, next_serial: 0, eof: false }
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
/// - **Phase 2**: `DecompressChunks` > `ReadChunkBlocks` > `CompressOutput`
pub struct SortWorkerPool {
    // Shared pipeline state (visible to workers and main thread)
    shared: Arc<SharedPipelineState>,

    /// Worker thread handles. `None` after shutdown (taken by `Drop` or `shutdown`).
    workers: Option<Vec<JoinHandle<()>>>,
    pub(crate) stats: PoolStats,
    pub(crate) pipeline_stats: Arc<SortPipelineStats>,
    pub buffer_pool: BufferPool,
    num_workers: usize,
}

/// Shared state visible to all workers and the main thread.
///
/// Uses `ArrayQueue` for all inter-step data queues (lock-free, non-blocking
/// `push()`/`pop()` only). The compress result channel stays as `crossbeam_channel`
/// because the I/O writer thread needs blocking `recv()`.
pub(crate) struct SharedPipelineState {
    /// Current phase: 0=shutdown, 1=Phase1, 2=Phase2, 255=Legacy.
    pub(crate) phase: AtomicU8,

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
    input_blocks_queued: AtomicU64,

    // --- Phase 2 queues ---
    /// Raw chunk blocks: `ReadChunkBlocks` → `DecompressChunks`.
    pub(crate) raw_chunk_blocks: Arc<ArrayQueue<TaggedRawBlock>>,
    /// Decompressed chunk blocks: `DecompressChunks` → main thread.
    pub(crate) decompressed_chunks: Arc<ArrayQueue<TaggedDecompressedBlock>>,
    /// All chunk files have reached EOF.
    pub(crate) all_chunks_eof: Arc<AtomicBool>,
    /// Number of sources that have reached EOF (when == `total_sources`, set `all_chunks_eof`).
    sources_at_eof: AtomicU64,
    /// Total number of chunk sources (set before Phase 2 starts).
    pub(crate) total_sources: AtomicU64,

    // --- Compress queue (shared Phase 1 + Phase 2) ---
    /// Compress jobs: main thread → workers (`ArrayQueue`, non-blocking push).
    pub(crate) compress_queue: Arc<ArrayQueue<CompressJob>>,

    // --- Per-worker chunk file distribution (Phase 2) ---
    /// One sender per worker; main thread sends `Vec<ChunkFileState>` before Phase 2.
    chunk_file_senders: Vec<std::sync::Mutex<Option<Sender<Vec<ChunkFileState>>>>>,
    /// One receiver per worker; worker receives its assigned files.
    chunk_file_receivers: Vec<Receiver<Vec<ChunkFileState>>>,

    /// Number of workers (for `low_water` threshold in backpressure).
    num_workers: usize,

    /// Main thread handle for `park()`/`unpark()` notification.
    /// Workers call `unpark()` after pushing to `decompressed_input` or `decompressed_chunks`
    /// so the main thread wakes immediately instead of spin-yielding.
    main_thread_handle: std::thread::Thread,
}

impl SharedPipelineState {
    fn new(num_workers: usize, main_thread_handle: std::thread::Thread) -> Self {
        let data_queue_cap = num_workers * 8;
        let compress_queue_cap = num_workers * 4;

        // Per-worker channels for chunk file distribution (one-shot, keep as channel)
        let mut senders = Vec::with_capacity(num_workers);
        let mut receivers = Vec::with_capacity(num_workers);
        for _ in 0..num_workers {
            let (tx, rx) = bounded::<Vec<ChunkFileState>>(1);
            senders.push(std::sync::Mutex::new(Some(tx)));
            receivers.push(rx);
        }

        Self {
            phase: AtomicU8::new(phase::LEGACY),

            input_file: std::sync::Mutex::new(None),
            input_eof: AtomicBool::new(false),
            input_read_error: Arc::new(AtomicBool::new(false)),
            input_read_serial: AtomicU64::new(0),
            raw_input_blocks: Arc::new(ArrayQueue::new(data_queue_cap)),
            decompressed_input: Arc::new(ArrayQueue::new(data_queue_cap)),
            decompressed_input_done: Arc::new(AtomicBool::new(false)),
            input_blocks_queued: AtomicU64::new(0),

            raw_chunk_blocks: Arc::new(ArrayQueue::new(data_queue_cap)),
            decompressed_chunks: Arc::new(ArrayQueue::new(data_queue_cap)),
            all_chunks_eof: Arc::new(AtomicBool::new(false)),
            sources_at_eof: AtomicU64::new(0),
            total_sources: AtomicU64::new(0),

            compress_queue: Arc::new(ArrayQueue::new(compress_queue_cap)),

            chunk_file_senders: senders,
            chunk_file_receivers: receivers,

            num_workers,
            main_thread_handle,
        }
    }

    /// Snapshot current queue depths for backpressure-driven scheduling.
    fn get_backpressure(&self) -> SortBackpressureState {
        let current_phase = self.phase.load(Ordering::Acquire);
        let low_water = self.num_workers;

        SortBackpressureState {
            decompressed_input_low: self.decompressed_input.len() < low_water,
            input_eof: self.input_eof.load(Ordering::Acquire),
            decompressed_input_done: self.decompressed_input_done.load(Ordering::Acquire),

            decompressed_chunks_low: self.decompressed_chunks.len() < low_water,
            all_chunks_eof: self.all_chunks_eof.load(Ordering::Acquire),

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
    decompressor: libdeflater::Decompressor,
    /// Worker's assigned chunk files for Phase 2 (worker i owns files i, i+N, i+2N...).
    chunk_files: Vec<ChunkFileState>,

    // Held items (one per step output) — see plan §Worker State
    /// Held raw input blocks from `ReadInputBlocks` (couldn't push to `raw_input_blocks` queue).
    held_raw_input_blocks: Vec<(u64, RawBgzfBlock)>,
    /// Held decompressed input block from `DecompressInput`.
    held_decompressed_input: Option<(u64, Vec<u8>)>,
    /// Held raw chunk blocks from `ReadChunkBlocks` (couldn't push to `raw_chunk_blocks` queue).
    held_raw_chunk_blocks: Vec<TaggedRawBlock>,
    /// Held decompressed chunk block from `DecompressChunks`.
    held_decompressed_chunk: Option<TaggedDecompressedBlock>,
    // Compress output goes directly to result_tx channel (I/O thread) — no held item.
    /// Backoff microseconds for idle spinning.
    backoff_us: u64,
}

impl SortWorkerState {
    /// Returns true if this worker is holding any items that need advancement.
    /// CRITICAL: Workers must not exit while holding items — they would be lost.
    fn has_any_held_items(&self) -> bool {
        !self.held_raw_input_blocks.is_empty()
            || self.held_decompressed_input.is_some()
            || !self.held_raw_chunk_blocks.is_empty()
            || self.held_decompressed_chunk.is_some()
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

    // Phase 2
    decompressed_chunks_low: bool,
    all_chunks_eof: bool,

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
            if bp.all_chunks_eof && !bp.compress_has_items {
                // All chunks at EOF and no compress work — nothing productive to do
                &[]
            } else if bp.compress_has_items && !bp.decompressed_chunks_low {
                // Output compression backpressure (15.4s at t4). Drain it.
                &[SortStep::CompressOutput, SortStep::DecompressChunks, SortStep::ReadChunkBlocks]
            } else {
                // Default: feed the merge loop, compress output when available
                &[SortStep::DecompressChunks, SortStep::ReadChunkBlocks, SortStep::CompressOutput]
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
/// Uses `yield_now()` at minimum backoff to avoid sleep syscall overhead.
fn sleep_with_jitter(backoff_us: u64) {
    if backoff_us <= MIN_BACKOFF_US {
        std::thread::yield_now();
    } else {
        let jitter_range = backoff_us / 4;
        let jitter_seed = u64::from(
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.subsec_nanos())
                .unwrap_or(0),
        );
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

impl SortWorkerPool {
    /// Create a new worker pool with `num_workers` threads.
    ///
    /// Each worker owns its own compressor (×2: spill and output) and decompressor.
    /// Workers are phase-aware and perform all CPU/IO work.
    ///
    /// - `temp_compression`: BGZF level for Phase 1 spill writes (typically 1 for speed).
    /// - `output_compression`: BGZF level for Phase 2 merge output (typically 6 for size).
    #[must_use]
    pub fn new(num_workers: usize, temp_compression: u32, output_compression: u32) -> Self {
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
                    let mut worker = SortWorkerState {
                        worker_id,
                        compressor: InlineBgzfCompressor::new(temp_compression),
                        output_compressor: InlineBgzfCompressor::new(output_compression),
                        decompressor: libdeflater::Decompressor::new(),
                        chunk_files: Vec::new(),
                        held_raw_input_blocks: Vec::new(),
                        held_decompressed_input: None,
                        held_raw_chunk_blocks: Vec::new(),
                        held_decompressed_chunk: None,
                        backoff_us: MIN_BACKOFF_US,
                    };

                    Self::worker_loop(&shared, &mut worker, &pstats);
                })
            })
            .collect();

        Self { shared, workers: Some(workers), stats, pipeline_stats, buffer_pool, num_workers }
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

            // 2. Check phase completion — wait for next phase, only exit on SHUTDOWN.
            //    Workers must survive across Phase 1 → Phase 2 transitions.
            if Self::is_phase_complete(shared, worker, current_phase)
                && !worker.has_any_held_items()
            {
                sleep_with_jitter(worker.backoff_us);
                worker.backoff_us = (worker.backoff_us * 2).min(MAX_BACKOFF_US);
                continue;
            }

            let mut did_work = false;

            // 3. Try to advance ALL held items first (deadlock prevention)
            did_work |= Self::try_advance_all_held(shared, worker);

            // Phase 2: receive chunk files on first iteration
            if current_phase == phase::PHASE2
                && worker.chunk_files.is_empty()
                && shared.total_sources.load(Ordering::Acquire) > 0
            {
                if let Ok(files) = shared.chunk_file_receivers[worker.worker_id].try_recv() {
                    worker.chunk_files = files;
                }
            }

            // 4. Get backpressure state and resolve priorities (done inline in step 6)
            let owned_step = Self::exclusive_step_for(worker.worker_id, shared, current_phase);

            // 5. Try owned exclusive step first (prevents starvation)
            if !did_work {
                if let Some(step) = owned_step {
                    if Self::is_step_eligible(step, shared, worker, current_phase) {
                        let t0 = Instant::now();
                        let result = Self::execute_step(shared, worker, step);
                        if result == StepResult::Success {
                            pstats.record_step(
                                worker.worker_id,
                                step,
                                Self::nanos_u64(t0.elapsed()),
                            );
                            did_work = true;
                        }
                    }
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
                        && !Self::can_attempt_exclusive(worker, step, shared, current_phase)
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
                sleep_with_jitter(worker.backoff_us);
                worker.backoff_us = (worker.backoff_us * 2).min(MAX_BACKOFF_US);
                pstats.record_idle(worker.worker_id, Self::nanos_u64(idle_start.elapsed()));
            }
        }
    }

    /// Check if the current phase is "complete" (no more work to do).
    ///
    /// This does NOT mean the worker should exit — it must also have no held items.
    fn is_phase_complete(
        shared: &SharedPipelineState,
        worker: &SortWorkerState,
        current_phase: u8,
    ) -> bool {
        match current_phase {
            phase::PHASE1 => {
                shared.decompressed_input_done.load(Ordering::Acquire)
                    && shared.compress_queue.is_empty()
            }
            phase::PHASE2 => {
                shared.all_chunks_eof.load(Ordering::Acquire)
                    && shared.raw_chunk_blocks.is_empty()
                    && shared.decompressed_chunks.is_empty()
                    && shared.compress_queue.is_empty()
                    && worker.chunk_files.iter().all(|f| f.eof)
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
    /// For `num_workers >= 2`: Worker 0 owns `ReadInputBlocks` (Phase 1) and
    /// `ReadChunkBlocks` (Phase 2). For `num_workers == 1`: single worker does
    /// everything (returns `None`, no ownership restrictions).
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
            // Phase 2: each worker reads its own files, no exclusive ownership needed
            _ => None,
        }
    }

    /// Whether a step is exclusive (requires ownership).
    ///
    /// Only `ReadInputBlocks` is exclusive — it reads from a shared input file
    /// protected by a Mutex. `ReadChunkBlocks` is NOT exclusive because each
    /// worker reads from its own assigned chunk files (no contention).
    fn is_exclusive_step(step: SortStep) -> bool {
        matches!(step, SortStep::ReadInputBlocks)
    }

    /// Whether this worker can attempt an exclusive step it doesn't own.
    ///
    /// Non-owner workers can attempt exclusive steps only if `num_workers == 1`
    /// (single worker mode).
    fn can_attempt_exclusive(
        worker: &SortWorkerState,
        step: SortStep,
        shared: &SharedPipelineState,
        current_phase: u8,
    ) -> bool {
        // Owner can always attempt
        if Self::exclusive_step_for(worker.worker_id, shared, current_phase) == Some(step) {
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
            SortStep::ReadChunkBlocks => {
                current_phase == phase::PHASE2
                    && worker.held_raw_chunk_blocks.is_empty()
                    && !worker.chunk_files.is_empty()
                    && !worker.chunk_files.iter().all(|f| f.eof)
            }
            SortStep::DecompressChunks => {
                current_phase == phase::PHASE2
                    && worker.held_decompressed_chunk.is_none()
                    && !shared.raw_chunk_blocks.is_empty()
            }
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
            SortStep::CompressSpill | SortStep::CompressOutput => {
                Self::try_compress(shared, worker)
            }
            SortStep::ReadChunkBlocks => Self::try_read_chunk_blocks(shared, worker),
            SortStep::DecompressChunks => Self::try_decompress_chunks(shared, worker),
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

        // Decompressed input (single) — always unpark main (drain full queue or deliver data)
        if worker.held_decompressed_input.is_some() {
            let pushed =
                try_advance_held(&shared.decompressed_input, &mut worker.held_decompressed_input);
            shared.main_thread_handle.unpark();
            if pushed {
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

        // Raw chunk blocks (batch)
        if !worker.held_raw_chunk_blocks.is_empty() {
            let before = worker.held_raw_chunk_blocks.len();
            try_advance_held_batch(&shared.raw_chunk_blocks, &mut worker.held_raw_chunk_blocks);
            if worker.held_raw_chunk_blocks.len() < before {
                advanced = true;
            }
        }

        // Decompressed chunks (single) — always unpark main
        if worker.held_decompressed_chunk.is_some() {
            advanced |=
                try_advance_held(&shared.decompressed_chunks, &mut worker.held_decompressed_chunk);
            shared.main_thread_handle.unpark();
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

        // Drop the lock before pushing to queue
        drop(guard);

        // Assign serial numbers and push to raw_input_blocks ArrayQueue.
        // Once the queue is full, hold this block and all remaining blocks.
        let mut blocks_iter = blocks.into_iter();
        for block in blocks_iter.by_ref() {
            let serial = shared.input_read_serial.fetch_add(1, Ordering::Relaxed);
            match shared.raw_input_blocks.push((serial, block)) {
                Ok(()) => {}
                Err((serial, block)) => {
                    worker.held_raw_input_blocks.push((serial, block));
                    break;
                }
            }
        }
        // Hold any remaining blocks we didn't attempt to push
        for block in blocks_iter {
            let serial = shared.input_read_serial.fetch_add(1, Ordering::Relaxed);
            worker.held_raw_input_blocks.push((serial, block));
        }

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

        let data = decompress_block(&block, &mut worker.decompressor)
            .expect("BGZF decompression should not fail for valid input blocks");

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
    // Phase 2 Steps
    // ========================================================================

    /// `ReadChunkBlocks`: read raw BGZF blocks from this worker's assigned chunk files.
    ///
    /// Worker-owns-files pattern: no locks needed. Worker i exclusively owns its files.
    fn try_read_chunk_blocks(
        shared: &SharedPipelineState,
        worker: &mut SortWorkerState,
    ) -> StepResult {
        if worker.chunk_files.is_empty() || shared.all_chunks_eof.load(Ordering::Acquire) {
            return StepResult::InputEmpty;
        }

        let mut read_any = false;

        // Round-robin through this worker's assigned files
        for file_state in &mut worker.chunk_files {
            if file_state.eof {
                continue;
            }

            let blocks = match read_raw_blocks(&mut file_state.reader, CHUNK_READ_BATCH_SIZE) {
                Ok(b) => b,
                Err(e) => {
                    log::error!(
                        "I/O error reading chunk file (source {}): {e}",
                        file_state.source_id
                    );
                    // Mark as EOF so the merge loop surfaces the truncation via
                    // the missing-serial check in io_writer_loop rather than hanging.
                    file_state.eof = true;
                    Self::maybe_mark_all_eof(shared);
                    continue;
                }
            };

            if blocks.is_empty() {
                file_state.eof = true;
                Self::maybe_mark_all_eof(shared);
                continue;
            }

            let mut blocks_iter = blocks.into_iter();
            let source_id = file_state.source_id;
            let mut backed_up = false;
            for block in blocks_iter.by_ref() {
                let serial = file_state.next_serial;
                file_state.next_serial += 1;
                let tagged = TaggedRawBlock { block, source_id, serial };
                match shared.raw_chunk_blocks.push(tagged) {
                    Ok(()) => {}
                    Err(tagged) => {
                        worker.held_raw_chunk_blocks.push(tagged);
                        backed_up = true;
                        break;
                    }
                }
            }
            // Hold any remaining blocks we didn't attempt to push
            if backed_up {
                for block in blocks_iter {
                    let serial = file_state.next_serial;
                    file_state.next_serial += 1;
                    worker.held_raw_chunk_blocks.push(TaggedRawBlock { block, source_id, serial });
                }
            }
            read_any = true;
            break; // Read from one file per step to be fair
        }

        if read_any { StepResult::Success } else { StepResult::InputEmpty }
    }

    /// `DecompressChunks`: decompress a raw chunk block.
    fn try_decompress_chunks(
        shared: &SharedPipelineState,
        worker: &mut SortWorkerState,
    ) -> StepResult {
        // Don't take new work if we're holding an item
        if worker.held_decompressed_chunk.is_some() {
            return StepResult::OutputFull;
        }

        let Some(tagged) = shared.raw_chunk_blocks.pop() else {
            return StepResult::InputEmpty;
        };

        let data = decompress_block(&tagged.block, &mut worker.decompressor)
            .expect("BGZF decompression should not fail for valid chunk blocks");

        let result =
            TaggedDecompressedBlock { data, source_id: tagged.source_id, serial: tagged.serial };

        match shared.decompressed_chunks.push(result) {
            Ok(()) => {
                shared.main_thread_handle.unpark();
            }
            Err(item) => {
                worker.held_decompressed_chunk = Some(item);
                // Queue full — wake main thread to drain so we can push next time
                shared.main_thread_handle.unpark();
            }
        }

        StepResult::Success
    }

    // ========================================================================
    // Compress Step (shared by Phase 1 + Phase 2)
    // ========================================================================

    /// Try to pick up a compress job from the `ArrayQueue` (non-blocking).
    fn try_compress(shared: &SharedPipelineState, worker: &mut SortWorkerState) -> StepResult {
        let Some(job) = shared.compress_queue.pop() else {
            return StepResult::InputEmpty;
        };

        // Use phase-appropriate compressor: Phase 1 spill uses temp compression level,
        // Phase 2 merge output uses output compression level.
        let compressor = if shared.phase.load(Ordering::Acquire) == phase::PHASE2 {
            &mut worker.output_compressor
        } else {
            &mut worker.compressor
        };
        Self::handle_compress_job(job, compressor);
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

    /// Get a clone of the decompressed chunks `ArrayQueue` for Phase 2 merge.
    pub(crate) fn decompressed_chunks_queue(
        &self,
    ) -> Arc<crossbeam_queue::ArrayQueue<TaggedDecompressedBlock>> {
        Arc::clone(&self.shared.decompressed_chunks)
    }

    /// Get a clone of the `all_chunks_eof` flag for Phase 2 merge.
    pub(crate) fn all_chunks_eof_flag(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.shared.all_chunks_eof)
    }

    /// Get a clone of the raw (pre-decompression) chunk blocks queue for Phase 2 merge.
    ///
    /// Used by `MainThreadChunkConsumer` to verify all blocks have been decompressed
    /// before declaring EOF — `all_chunks_eof` is set when raw reads finish, but blocks
    /// may still be in this queue awaiting decompression.
    pub(crate) fn raw_chunk_blocks_queue(&self) -> Arc<ArrayQueue<TaggedRawBlock>> {
        Arc::clone(&self.shared.raw_chunk_blocks)
    }

    /// Set the current pipeline phase.
    pub fn set_phase(&self, new_phase: u8) {
        self.shared.phase.store(new_phase, Ordering::Release);
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

    /// Distribute chunk files to workers for Phase 2.
    ///
    /// Worker `i` gets files `i, i+N, i+2N, ...` where `N` is `num_workers`.
    /// Must be called before `set_phase(PHASE2)`.
    ///
    /// # Errors
    ///
    /// Returns an error if a chunk file cannot be opened.
    ///
    /// # Panics
    ///
    /// Panics if a sender mutex is poisoned.
    pub fn distribute_chunk_files(
        &self,
        files: Vec<(std::path::PathBuf, usize)>,
    ) -> anyhow::Result<()> {
        let n = self.num_workers;
        let total_sources = files.len();
        self.shared.total_sources.store(total_sources as u64, Ordering::Release);

        // Reset EOF state
        self.shared.all_chunks_eof.store(false, Ordering::Release);
        self.shared.sources_at_eof.store(0, Ordering::Release);

        // Assign files to workers: worker i gets files i, i+N, i+2N...
        let mut per_worker: Vec<Vec<ChunkFileState>> = (0..n).map(|_| Vec::new()).collect();

        for (idx, (path, source_id)) in files.into_iter().enumerate() {
            let worker_idx = idx % n;
            let file = std::fs::File::open(&path).map_err(|e| {
                anyhow::anyhow!("Failed to open chunk file {}: {e}", path.display())
            })?;
            let reader = BufReader::with_capacity(2 * 1024 * 1024, file);
            per_worker[worker_idx].push(ChunkFileState::new(reader, source_id));
        }

        // Send to each worker via their dedicated channel
        for (worker_idx, files) in per_worker.into_iter().enumerate() {
            let guard = self.shared.chunk_file_senders[worker_idx]
                .lock()
                .expect("chunk_file_sender mutex should not be poisoned");
            if let Some(tx) = guard.as_ref() {
                let _ = tx.send(files);
            }
        }
        Ok(())
    }

    /// Submit a compression job to the pool (non-blocking, spin-yield on full).
    ///
    /// The main thread calls this during spill writes and merge output. If the
    /// compress `ArrayQueue` is full, spins briefly with `yield_now()` — acceptable
    /// because the main thread has no other productive work during spill writes.
    pub fn submit_compress(&self, job: CompressJob) {
        self.stats.compress_jobs.fetch_add(1, Ordering::Relaxed);
        let mut job = job;
        loop {
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
        // Drop chunk file senders so workers blocked on receiving them unblock.
        for sender in &self.shared.chunk_file_senders {
            if let Ok(mut guard) = sender.lock() {
                *guard = None;
            }
        }
        if let Some(workers) = self.workers.take() {
            for w in workers {
                let _ = w.join();
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
    fn handle_compress_job(job: CompressJob, compressor: &mut InlineBgzfCompressor) {
        compressor
            .write_all(&job.data)
            .expect("BGZF compression write should not fail for valid data");
        compressor.flush().expect("BGZF compression flush should not fail");

        let blocks = compressor.take_blocks();
        let total_len: usize = blocks.iter().map(|b| b.data.len()).sum();
        let mut compressed = Vec::with_capacity(total_len);
        for block in &blocks {
            compressed.extend_from_slice(&block.data);
        }

        let mut recycled = job.data;
        recycled.clear();

        let serial = job.serial;
        if job
            .result_tx
            .send(CompressResult { serial, compressed, recycled_buf: recycled })
            .is_err()
        {
            log::warn!(
                "compress result discarded (serial {serial}): I/O writer thread disconnected"
            );
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

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

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
        stats.compress_jobs.fetch_add(42, Ordering::Relaxed);
        // Verify log_summary executes without panicking (logging may be a no-op in tests).
        stats.log_summary();
        assert_eq!(stats.compress_jobs.load(Ordering::Relaxed), 42);
    }

    #[test]
    fn test_pool_compress_roundtrip() {
        let pool = SortWorkerPool::new(2, 1, 6);
        let (result_tx, result_rx) = pool.compress_result_channel();

        // Submit a compress job
        let data = vec![b'A'; 1000];
        pool.submit_compress(CompressJob { data, serial: 0, result_tx });

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
        let pool = SortWorkerPool::new(4, 1, 6);
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
    fn test_pool_stats() {
        let pool = SortWorkerPool::new(2, 1, 6);
        let (c_tx, c_rx) = pool.compress_result_channel();

        // Submit one compress job
        pool.submit_compress(CompressJob { data: vec![b'A'; 100], serial: 0, result_tx: c_tx });
        let _ = c_rx.recv();

        assert_eq!(pool.stats.compress_jobs.load(Ordering::Relaxed), 1);

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
            decompressed_chunks_low: false,
            all_chunks_eof: false,
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
            decompressed_chunks_low: false,
            all_chunks_eof: false,
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
            decompressed_chunks_low: false,
            all_chunks_eof: false,
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
            decompressed_chunks_low: true,
            all_chunks_eof: false,
            compress_has_items: false,
            phase: phase::PHASE2,
        };
        let priorities = get_sort_priorities(&bp);
        assert_eq!(priorities[0], SortStep::DecompressChunks);
    }

    #[test]
    fn test_sort_priorities_phase2_compress_backpressure() {
        let bp = SortBackpressureState {
            decompressed_input_low: false,
            input_eof: false,
            decompressed_input_done: false,
            decompressed_chunks_low: false,
            all_chunks_eof: false,
            compress_has_items: true,
            phase: phase::PHASE2,
        };
        let priorities = get_sort_priorities(&bp);
        assert_eq!(priorities[0], SortStep::CompressOutput);
    }

    #[test]
    fn test_sort_priorities_phase2_all_done_returns_empty() {
        let bp = SortBackpressureState {
            decompressed_input_low: false,
            input_eof: false,
            decompressed_input_done: false,
            decompressed_chunks_low: false,
            all_chunks_eof: true,
            compress_has_items: false,
            phase: phase::PHASE2,
        };
        assert!(get_sort_priorities(&bp).is_empty());
    }

    #[test]
    fn test_sort_priorities_legacy_returns_compress_only() {
        let bp = SortBackpressureState {
            decompressed_input_low: false,
            input_eof: false,
            decompressed_input_done: false,
            decompressed_chunks_low: false,
            all_chunks_eof: false,
            compress_has_items: true,
            phase: phase::LEGACY,
        };
        let priorities = get_sort_priorities(&bp);
        assert_eq!(priorities.len(), 1);
        assert_eq!(priorities[0], SortStep::CompressSpill);
    }

    #[test]
    fn test_worker_pool_num_workers() {
        let pool = SortWorkerPool::new(3, 1, 6);
        assert_eq!(pool.num_workers(), 3);
        pool.shutdown();
    }
}
