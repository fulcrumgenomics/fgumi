//! Core infrastructure, traits, and shared types for the unified pipeline.
//!
//! This module contains:
//! - Batch buffer types (`RawBlockBatch`, `CompressedBlockBatch`, etc.)
//! - BGZF compression types
//! - `PipelineConfig` and `PipelineStats`
//! - `OutputPipelineState` trait for shared step functions
//! - `HasCompressor` trait
//! - `BatchWeight` trait for template-based batching
//! - Shared step functions (`shared_try_step_compress`, `shared_try_step_write`)

use crossbeam_queue::ArrayQueue;
use log::info;
use parking_lot::Mutex;
use std::io::{self, Read, Write};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU8, AtomicU64, Ordering};
use std::thread;
use std::time::{Duration, Instant};

use crate::progress::ProgressTracker;

use super::deadlock::{DeadlockState, QueueSnapshot, check_deadlock_and_restore};

use crate::bgzf_reader::{RawBgzfBlock, decompress_block_into, read_raw_blocks};
use crate::read_info::LibraryIndex;
use crate::reorder_buffer::ReorderBuffer;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::data::field::Tag;

use super::scheduler::{BackpressureState, Scheduler, SchedulerStrategy, create_scheduler};

// ============================================================================
// Memory Breakdown Structure
// ============================================================================

/// Memory breakdown for comprehensive debugging
#[cfg(feature = "memory-debug")]
#[derive(Debug, Clone)]
pub struct MemoryBreakdown {
    /// System RSS in GB
    pub system_rss_gb: f64,
    /// Total tracked memory in GB
    pub tracked_total_gb: f64,
    /// Untracked memory (allocator overhead, etc.) in GB
    pub untracked_gb: f64,

    // Queue memory in MB/GB
    pub q1_mb: f64,
    pub q2_mb: f64,
    pub q3_mb: f64,
    pub q4_gb: f64,
    pub q5_gb: f64,
    pub q6_mb: f64,
    pub q7_mb: f64,

    // Processing memory in GB/MB
    pub position_groups_gb: f64,
    pub templates_gb: f64,
    pub reorder_buffers_mb: f64,
    pub grouper_mb: f64,
    pub worker_local_mb: f64,

    // Infrastructure memory in MB
    pub decompressors_mb: f64,
    pub compressors_mb: f64,
    pub worker_buffers_mb: f64,
    pub io_buffers_mb: f64,
    pub thread_stacks_mb: f64,
    pub queue_capacity_mb: f64,
    /// Total infrastructure memory in GB (sum of above)
    pub infrastructure_gb: f64,
}

/// Comprehensive memory tracking fields, grouped into a single struct to keep
/// `PipelineStats` readable. All fields are `AtomicU64` for lock-free updates.
#[cfg(feature = "memory-debug")]
#[derive(Debug)]
pub struct MemoryDebugStats {
    // Queue-specific memory tracking
    /// Memory held in Q1 (raw BGZF blocks)
    pub q1_memory_bytes: AtomicU64,
    /// Memory held in Q2 (decompressed blocks)
    pub q2_memory_bytes: AtomicU64,
    /// Memory held in Q3 (decoded records)
    pub q3_memory_bytes: AtomicU64,
    /// Memory held in Q4 (position groups) - likely the big one
    pub q4_memory_bytes: AtomicU64,
    /// Memory held in Q5 (processed groups)
    pub q5_memory_bytes: AtomicU64,
    /// Memory held in Q6 (serialized data)
    pub q6_memory_bytes: AtomicU64,
    /// Memory held in Q7 (compressed output blocks)
    pub q7_memory_bytes: AtomicU64,

    // Processing memory tracking
    /// Memory held in position groups during processing
    pub position_group_processing_bytes: AtomicU64,
    /// Memory held in templates during operations
    pub template_processing_bytes: AtomicU64,
    /// Memory held in reorder buffers
    pub reorder_buffer_bytes: AtomicU64,
    /// Memory held by grouper state
    pub grouper_memory_bytes: AtomicU64,
    /// Memory held by worker-local allocations
    pub worker_local_memory_bytes: AtomicU64,

    // Infrastructure memory (known constants set once at pipeline startup)
    /// Memory for per-thread decompressor instances (libdeflater)
    pub decompressor_memory_bytes: AtomicU64,
    /// Memory for per-thread compressor instances (InlineBgzfCompressor)
    pub compressor_memory_bytes: AtomicU64,
    /// Memory for per-thread serialization + decompression buffers
    pub worker_buffer_memory_bytes: AtomicU64,
    /// Memory for I/O buffers (BufReader + BufWriter)
    pub io_buffer_memory_bytes: AtomicU64,
    /// Memory for thread stacks (2MB per thread default)
    pub thread_stack_memory_bytes: AtomicU64,
    /// Memory for ArrayQueue pre-allocation overhead
    pub queue_capacity_memory_bytes: AtomicU64,

    // System memory tracking
    /// Actual system RSS (from /proc/self/status or sysinfo)
    pub system_rss_bytes: AtomicU64,
}

#[cfg(feature = "memory-debug")]
impl MemoryDebugStats {
    /// Create a new memory debug stats collector with all counters at zero.
    #[must_use]
    pub fn new() -> Self {
        Self {
            q1_memory_bytes: AtomicU64::new(0),
            q2_memory_bytes: AtomicU64::new(0),
            q3_memory_bytes: AtomicU64::new(0),
            q4_memory_bytes: AtomicU64::new(0),
            q5_memory_bytes: AtomicU64::new(0),
            q6_memory_bytes: AtomicU64::new(0),
            q7_memory_bytes: AtomicU64::new(0),
            position_group_processing_bytes: AtomicU64::new(0),
            template_processing_bytes: AtomicU64::new(0),
            reorder_buffer_bytes: AtomicU64::new(0),
            grouper_memory_bytes: AtomicU64::new(0),
            worker_local_memory_bytes: AtomicU64::new(0),
            decompressor_memory_bytes: AtomicU64::new(0),
            compressor_memory_bytes: AtomicU64::new(0),
            worker_buffer_memory_bytes: AtomicU64::new(0),
            io_buffer_memory_bytes: AtomicU64::new(0),
            thread_stack_memory_bytes: AtomicU64::new(0),
            queue_capacity_memory_bytes: AtomicU64::new(0),
            system_rss_bytes: AtomicU64::new(0),
        }
    }
}

// ============================================================================
// Thread-Local Memory Tracking
// ============================================================================

#[cfg(feature = "memory-debug")]
/// Sentinel value meaning "no thread ID assigned yet".
const THREAD_ID_UNSET: usize = usize::MAX;

#[cfg(feature = "memory-debug")]
thread_local! {
    static THREAD_ID: std::cell::Cell<usize> = std::cell::Cell::new(THREAD_ID_UNSET);
}

#[cfg(feature = "memory-debug")]
/// Get or assign a thread ID for memory tracking.
/// Returns IDs in range 0..MAX_THREADS. IDs wrap if more than MAX_THREADS
/// unique threads call this function (per-thread data will be shared).
pub fn get_or_assign_thread_id() -> usize {
    THREAD_ID.with(|id| {
        let current = id.get();
        if current == THREAD_ID_UNSET {
            use std::sync::atomic::{AtomicUsize, Ordering};
            static THREAD_COUNTER: AtomicUsize = AtomicUsize::new(0);

            let new_id = THREAD_COUNTER.fetch_add(1, Ordering::Relaxed) % MAX_THREADS;
            id.set(new_id);
            new_id
        } else {
            current
        }
    })
}

// ============================================================================
// System Memory Utilities
// ============================================================================

#[cfg(feature = "memory-debug")]
/// Get current process RSS from /proc/self/status (Linux) or sysinfo (cross-platform)
pub fn get_process_rss_bytes() -> Option<u64> {
    // Try Linux-style /proc/self/status first (most accurate)
    if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
        return status
            .lines()
            .find(|line| line.starts_with("VmRSS:"))?
            .split_whitespace()
            .nth(1)?
            .parse::<u64>()
            .ok()
            .map(|kb| kb * 1024); // Convert KB to bytes
    }

    // Cross-platform fallback using sysinfo (cached instance)
    use std::sync::Mutex;
    use sysinfo::{ProcessRefreshKind, RefreshKind, System};

    static RSS_SYSTEM: std::sync::OnceLock<Mutex<System>> = std::sync::OnceLock::new();

    let sys = RSS_SYSTEM.get_or_init(|| {
        Mutex::new(System::new_with_specifics(
            RefreshKind::nothing().with_processes(ProcessRefreshKind::nothing().with_memory()),
        ))
    });

    let mut sys_guard = sys.lock().ok()?;
    sys_guard.refresh_processes_specifics(
        sysinfo::ProcessesToUpdate::All,
        false,
        ProcessRefreshKind::nothing().with_memory(),
    );

    // Get current process RSS
    let pid = sysinfo::get_current_pid().ok()?;
    let process = sys_guard.process(pid)?;
    Some(process.memory()) // sysinfo returns bytes
}

#[cfg(feature = "memory-debug")]
/// Log comprehensive memory statistics with accuracy tracking
pub fn log_comprehensive_memory_stats(stats: &PipelineStats) {
    // Update RSS
    if let Some(rss) = get_process_rss_bytes() {
        stats.update_system_rss(rss);
    }

    // TODO: Queue memory integration deferred - requires pipeline reference to read actual queues
    // For now, queue memory is tracked through estimates only

    let breakdown = stats.get_memory_breakdown();

    // Main memory line with RSS vs tracked accuracy
    if breakdown.system_rss_gb > 0.0 {
        let pct = (breakdown.tracked_total_gb / breakdown.system_rss_gb * 100.0) as u32;
        log::info!(
            "MEMORY: RSS={:.1}GB Tracked={:.1}GB ({}%) | Queue: Q1:{:.0}MB Q2:{:.0}MB Q3:{:.0}MB Q4:{:.1}GB Q5:{:.1}GB Q6:{:.0}MB Q7:{:.0}MB | Proc: Pos={:.1}GB Tmpl={:.1}GB | Infra={:.0}MB",
            breakdown.system_rss_gb,
            breakdown.tracked_total_gb,
            pct,
            breakdown.q1_mb,
            breakdown.q2_mb,
            breakdown.q3_mb,
            breakdown.q4_gb,
            breakdown.q5_gb,
            breakdown.q6_mb,
            breakdown.q7_mb,
            breakdown.position_groups_gb,
            breakdown.templates_gb,
            breakdown.infrastructure_gb * 1e3, // Convert GB to MB for display
        );
    } else {
        log::info!(
            "MEMORY: Tracked={:.1}GB | Queue: Q1:{:.0}MB Q2:{:.0}MB Q3:{:.0}MB Q4:{:.1}GB Q5:{:.1}GB Q6:{:.0}MB Q7:{:.0}MB | Proc: Pos={:.1}GB Tmpl={:.1}GB | Infra={:.0}MB",
            breakdown.tracked_total_gb,
            breakdown.q1_mb,
            breakdown.q2_mb,
            breakdown.q3_mb,
            breakdown.q4_gb,
            breakdown.q5_gb,
            breakdown.q6_mb,
            breakdown.q7_mb,
            breakdown.position_groups_gb,
            breakdown.templates_gb,
            breakdown.infrastructure_gb * 1e3,
        );
    }

    // Untracked memory details (only if RSS is available)
    if breakdown.system_rss_gb > 0.0 {
        let untracked_pct = ((breakdown.untracked_gb / breakdown.system_rss_gb) * 100.0) as u32;
        if breakdown.untracked_gb > 1.0 {
            log::info!(
                "   Untracked: {:.1}GB ({}%) = allocator fragmentation + noodles internals",
                breakdown.untracked_gb,
                untracked_pct,
            );
        }
    }
}

#[cfg(feature = "memory-debug")]
/// Start memory monitoring thread that reports at the given interval.
pub fn start_memory_monitor(
    stats: Arc<PipelineStats>,
    shutdown_signal: Arc<AtomicBool>,
    report_interval_secs: u64,
) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        let mut last_report = Instant::now();
        let report_interval = Duration::from_secs(report_interval_secs);

        let mut last_rss: u64 = 0;
        let mut peak_rss: u64 = 0;
        let mut stats_printed = false;
        while !shutdown_signal.load(Ordering::Relaxed) {
            if last_report.elapsed() >= report_interval {
                log_comprehensive_memory_stats(&stats);
                let current_rss = stats.memory.system_rss_bytes.load(Ordering::Relaxed);
                // Print mimalloc stats when RSS starts declining (just past peak)
                if !stats_printed
                    && current_rss > 0
                    && last_rss > 0
                    && current_rss < last_rss
                    && peak_rss > 4_000_000_000
                {
                    log::info!("=== MIMALLOC STATS AT PEAK (no mi_collect) ===");
                    // SAFETY: `mi_stats_print_out(None, null_mut())` prints mimalloc allocator
                    // statistics to stderr using the default output handler. mimalloc internally
                    // synchronizes stats collection, making this safe to call concurrently with
                    // allocation/deallocation on other threads.
                    unsafe {
                        libmimalloc_sys::mi_stats_print_out(None, std::ptr::null_mut());
                    }
                    stats_printed = true;
                }
                if current_rss > peak_rss {
                    peak_rss = current_rss;
                }
                last_rss = current_rss;
                last_report = Instant::now();
            }
            thread::sleep(Duration::from_millis(100));
        }
    })
}

// ============================================================================
// Batch Weight Trait (for template-based batching)
// ============================================================================

/// Trait for groups that can report their "weight" for batching purposes.
///
/// The weight is typically the number of templates in the group, allowing
/// the pipeline to batch groups based on total templates rather than group count.
/// This provides more consistent batch sizes across datasets with varying
/// templates-per-group ratios.
///
/// # Example
///
/// For a position group with 10 templates, `batch_weight()` returns 10.
/// The pipeline accumulates groups until the total weight reaches a threshold
/// (e.g., 500 templates), then flushes the batch.
pub trait BatchWeight {
    /// Returns the weight of this group for batching purposes.
    /// For position groups, this is typically the number of templates.
    fn batch_weight(&self) -> usize;
}

// ============================================================================
// Memory Estimate Trait (for memory-based queue limits)
// ============================================================================

/// Trait for types that can estimate their heap memory usage.
///
/// This is used by memory-bounded queues to track how much memory is
/// consumed by items in the queue. The estimate should include heap
/// allocations (Vec contents, String data, etc.) but not the struct itself.
///
/// # Example
///
/// For a batch containing a `Vec<u8>` with 1000 bytes, `estimate_heap_size()`
/// returns approximately 1000 (the Vec's heap allocation).
pub trait MemoryEstimate {
    /// Returns an estimate of the heap memory used by this value in bytes.
    ///
    /// This should include:
    /// - Vec/String heap allocations (capacity, not just len)
    /// - Nested struct heap allocations
    ///
    /// This should NOT include:
    /// - The size of the struct itself (stack size)
    /// - Shared/reference-counted data (counted once, not per-reference)
    fn estimate_heap_size(&self) -> usize;
}

// ============================================================================
// Memory Tracking for Queues
// ============================================================================

/// Tracks memory usage across pipeline queues.
///
/// This tracker uses atomic operations to enable lock-free memory accounting.
/// It's designed to provide backpressure when queue memory exceeds a limit,
/// preventing unbounded memory growth in threaded pipelines.
///
/// # Usage
///
/// ```ignore
/// let tracker = MemoryTracker::new(1_000_000_000); // 1 GB limit
///
/// // Before pushing to a queue
/// let item_size = item.estimate_heap_size();
/// if tracker.try_add(item_size) {
///     queue.push(item);
/// } else {
///     // Queue is at memory limit - apply backpressure
/// }
///
/// // After popping from a queue
/// let item = queue.pop();
/// let item_size = item.estimate_heap_size();
/// tracker.remove(item_size);
/// ```
#[derive(Debug)]
#[allow(clippy::struct_field_names)]
pub struct MemoryTracker {
    /// Current memory usage in bytes.
    current_bytes: AtomicU64,
    /// Peak memory usage in bytes (for stats).
    peak_bytes: AtomicU64,
    /// Maximum allowed memory in bytes. 0 = no limit.
    limit_bytes: u64,
}

impl MemoryTracker {
    /// Create a new memory tracker with the specified limit.
    ///
    /// # Arguments
    /// * `limit_bytes` - Maximum allowed memory in bytes. Use 0 for no limit.
    #[must_use]
    pub fn new(limit_bytes: u64) -> Self {
        Self { current_bytes: AtomicU64::new(0), peak_bytes: AtomicU64::new(0), limit_bytes }
    }

    /// Create a new memory tracker with no limit.
    #[must_use]
    pub fn unlimited() -> Self {
        Self::new(0)
    }

    /// Try to add bytes to the tracker.
    ///
    /// Returns `true` if the memory was added, `false` if we're already at or
    /// above the limit. The key behavior is: if we're currently under the limit,
    /// any single addition succeeds (even if it would exceed the limit). We only
    /// reject additions when we're already at or above capacity.
    ///
    /// Note: When the limit is 0 (no limit), this always returns `true`.
    pub fn try_add(&self, bytes: usize) -> bool {
        if self.limit_bytes == 0 {
            // No limit - just track without checking
            let new_total =
                self.current_bytes.fetch_add(bytes as u64, Ordering::Relaxed) + bytes as u64;
            self.update_peak(new_total);
            return true;
        }

        // Use compare-and-swap loop to atomically check and add
        loop {
            let current = self.current_bytes.load(Ordering::Relaxed);

            // If we're already at or above the limit, reject the addition
            if current >= self.limit_bytes {
                return false;
            }

            // We're under the limit, so allow this addition (even if it exceeds)
            let new_total = current + bytes as u64;

            // Try to update atomically
            if self
                .current_bytes
                .compare_exchange_weak(current, new_total, Ordering::Relaxed, Ordering::Relaxed)
                .is_ok()
            {
                self.update_peak(new_total);
                return true;
            }
            // CAS failed - another thread modified, retry
        }
    }

    /// Update peak memory if `new_total` is higher.
    #[inline]
    fn update_peak(&self, new_total: u64) {
        let mut current_peak = self.peak_bytes.load(Ordering::Relaxed);
        while new_total > current_peak {
            match self.peak_bytes.compare_exchange_weak(
                current_peak,
                new_total,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(actual) => current_peak = actual,
            }
        }
    }

    /// Remove bytes from the tracker (call when item is consumed from queue).
    /// Uses saturating subtraction to prevent underflow from estimation mismatches.
    pub fn remove(&self, bytes: usize) {
        let bytes = bytes as u64;
        let mut current = self.current_bytes.load(Ordering::Relaxed);
        loop {
            let new_val = current.saturating_sub(bytes);
            match self.current_bytes.compare_exchange_weak(
                current,
                new_val,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(actual) => current = actual,
            }
        }
    }

    /// Add bytes unconditionally (for tracking after a successful push).
    /// Unlike `try_add`, this always adds and doesn't check the limit.
    pub fn add(&self, bytes: usize) {
        let new_total =
            self.current_bytes.fetch_add(bytes as u64, Ordering::Relaxed) + bytes as u64;
        self.update_peak(new_total);
    }

    /// Get the current memory usage in bytes.
    #[must_use]
    pub fn current(&self) -> u64 {
        self.current_bytes.load(Ordering::Relaxed)
    }

    /// Get the memory limit in bytes. Returns 0 if no limit.
    #[must_use]
    pub fn limit(&self) -> u64 {
        self.limit_bytes
    }

    /// Get the peak memory usage in bytes.
    #[must_use]
    pub fn peak(&self) -> u64 {
        self.peak_bytes.load(Ordering::Relaxed)
    }

    /// Check if we're at or above the backpressure threshold.
    ///
    /// For optimal performance, we apply backpressure at a capped threshold
    /// (512MB) regardless of the configured limit. This keeps the pipeline
    /// lean with good cache behavior. The configured limit acts as a safety
    /// cap but doesn't delay backpressure.
    ///
    /// Even with no limit (`limit_bytes` == 0), we still apply backpressure
    /// at 512MB to maintain optimal throughput.
    #[must_use]
    pub fn is_at_limit(&self) -> bool {
        // Always apply backpressure at 512MB for optimal cache behavior.
        // If a lower limit is configured, use that instead.
        let backpressure_threshold = if self.limit_bytes == 0 {
            BACKPRESSURE_THRESHOLD_BYTES
        } else {
            self.limit_bytes.min(BACKPRESSURE_THRESHOLD_BYTES)
        };
        self.current() >= backpressure_threshold
    }

    /// Check if we've drained below the low-water mark.
    ///
    /// This is used for hysteresis to prevent thrashing: we enter drain mode
    /// when at the backpressure threshold, but only exit when we've drained
    /// to half that threshold.
    ///
    /// The threshold is 256MB (half of the 512MB backpressure cap), or half
    /// of a lower configured limit if one is set.
    #[must_use]
    pub fn is_below_drain_threshold(&self) -> bool {
        // Drain threshold is half of the backpressure threshold
        let backpressure_threshold = if self.limit_bytes == 0 {
            BACKPRESSURE_THRESHOLD_BYTES
        } else {
            self.limit_bytes.min(BACKPRESSURE_THRESHOLD_BYTES)
        };
        self.current() < backpressure_threshold / 2
    }
}

impl Default for MemoryTracker {
    fn default() -> Self {
        Self::unlimited()
    }
}

// ============================================================================
// Step Result Enum
// ============================================================================

/// Result of attempting a pipeline step.
///
/// This enum enables non-blocking step functions that can report why they
/// couldn't make progress, allowing the scheduler to make informed decisions
/// about which step to try next.
///
/// # Deadlock Prevention
///
/// By returning immediately with `OutputFull` or `InputEmpty` instead of
/// blocking, workers can try other steps (especially Write to drain queues),
/// preventing the deadlock that occurs when all threads block on push.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepResult {
    /// Successfully processed and advanced an item.
    Success,
    /// Output queue is full - should try downstream steps to drain.
    OutputFull,
    /// Input queue is empty - should try upstream steps to fill.
    InputEmpty,
}

impl StepResult {
    /// Returns true if the step made progress.
    #[inline]
    #[must_use]
    pub fn is_success(self) -> bool {
        matches!(self, StepResult::Success)
    }
}

/// Progress logging interval for the 7-step pipeline.
pub const PROGRESS_LOG_INTERVAL: u64 = 1_000_000;

/// Memory threshold for scheduler backpressure signaling (512 MB).
///
/// This constant controls when the pipeline signals memory pressure to the scheduler,
/// which responds by prioritizing downstream steps (Process, Serialize, Compress, Write)
/// to drain buffered data before allowing more input.
///
/// # Architecture
///
/// Memory is tracked at the reorder buffer BEFORE the Group step:
/// - **BAM pipeline**: `q3_reorder_heap_bytes` (after Decode, before Group)
/// - **FASTQ pipeline**: `q2_reorder_heap_bytes` (after Decompress, before Group)
///
/// This placement is critical: we track memory before the exclusive Group step rather
/// than after it. Tracking after Group creates a coordination problem where memory is
/// released from the pre-Group buffer before knowing if the post-Group queue can accept
/// it, potentially leaving data in an untracked intermediate buffer.
///
/// # Threshold Behavior
///
/// - When tracked memory >= `BACKPRESSURE_THRESHOLD_BYTES`, `is_memory_high()` returns true
/// - The scheduler enters "drain mode" and prioritizes output steps
/// - When tracked memory < `BACKPRESSURE_THRESHOLD_BYTES / 2`, `is_memory_drained()` returns true
/// - The scheduler exits drain mode (hysteresis prevents thrashing)
///
/// # Relationship to Queue Memory Limit
///
/// This threshold is independent of `queue_memory_limit` (default 4GB), which controls
/// the hard limit for `can_decompress_proceed()` / `can_decode_proceed()`. The backpressure
/// threshold is deliberately lower to allow gradual slowdown before hitting hard limits.
///
/// The effective threshold is `min(queue_memory_limit, BACKPRESSURE_THRESHOLD_BYTES)`,
/// so if a smaller memory limit is configured, backpressure activates at that limit.
pub const BACKPRESSURE_THRESHOLD_BYTES: u64 = 512 * 1024 * 1024; // 512 MB

/// Q5 (processed queue) backpressure threshold.
///
/// This is set lower than the Q3 threshold (256 MB vs 512 MB) because items in Q5
/// are typically larger (e.g., `SimplexProcessedBatch` with `RecordBuf` vectors).
/// When Q5 memory exceeds this threshold, the Process step pauses to let
/// downstream steps (Serialize, Compress, Write) catch up.
pub const Q5_BACKPRESSURE_THRESHOLD_BYTES: u64 = 256 * 1024 * 1024; // 256 MB

// ============================================================================
// Input-Half Reorder Buffer State (shared between BAM and FASTQ pipelines)
// ============================================================================

/// Atomic state for tracking a reorder buffer's memory and sequence position.
///
/// This struct encapsulates the atomic counters needed for lock-free backpressure
/// decisions on reorder buffers. Both BAM and FASTQ pipelines use this pattern
/// for their pre-Group reorder buffers (Q3 in BAM, Q2 in FASTQ).
///
/// # Usage Pattern
///
/// ```ignore
/// // Producer (Decompress/Decode step):
/// if state.can_proceed(serial, memory_limit) {
///     // Push to reorder buffer and update heap bytes
///     state.add_heap_bytes(item_size);
/// }
///
/// // Consumer (Group step):
/// if let Some(item) = reorder_buffer.try_pop_next() {
///     state.update_next_seq(new_seq);
///     state.sub_heap_bytes(item_size);
/// }
/// ```
#[derive(Debug)]
pub struct ReorderBufferState {
    /// Next sequence number the reorder buffer needs to make progress.
    /// Updated by the consumer (Group step) after popping from the buffer.
    pub next_seq: AtomicU64,
    /// Current heap bytes tracked in the reorder buffer.
    /// Updated by producer (add) and consumer (sub).
    pub heap_bytes: AtomicU64,
    /// Memory limit for backpressure. 0 means use default threshold.
    memory_limit: u64,
}

impl ReorderBufferState {
    /// Create a new reorder buffer state with the given memory limit.
    ///
    /// # Arguments
    /// * `memory_limit` - Memory limit in bytes. Use 0 for default threshold.
    #[must_use]
    pub fn new(memory_limit: u64) -> Self {
        Self { next_seq: AtomicU64::new(0), heap_bytes: AtomicU64::new(0), memory_limit }
    }

    /// Check if a producer can proceed with the given serial number.
    ///
    /// Always allows the serial that the consumer needs (`next_seq`) to proceed,
    /// even if over memory limit. This ensures the consumer can always make progress.
    /// For other serials, applies backpressure if over 50% of the limit.
    #[must_use]
    pub fn can_proceed(&self, serial: u64) -> bool {
        let limit = self.effective_limit();
        let next_seq = self.next_seq.load(Ordering::Acquire);
        let heap_bytes = self.heap_bytes.load(Ordering::Acquire);

        // Always accept the serial the reorder buffer needs to make progress.
        // This fills gaps and allows the consumer to pop items, reducing memory.
        if serial == next_seq {
            return true;
        }

        // Apply backpressure if over 50% of limit.
        let effective_limit = limit / 2;
        heap_bytes < effective_limit
    }

    /// Check if memory is at the backpressure threshold.
    ///
    /// Returns true when tracked memory >= threshold, signaling that the
    /// scheduler should enter drain mode and prioritize output steps.
    #[must_use]
    pub fn is_memory_high(&self) -> bool {
        let threshold = self.effective_limit();
        self.heap_bytes.load(Ordering::Acquire) >= threshold
    }

    /// Check if memory has drained below the low-water mark.
    ///
    /// Provides hysteresis to prevent thrashing: enter drain mode at backpressure
    /// threshold, only exit when drained to half that threshold.
    #[must_use]
    pub fn is_memory_drained(&self) -> bool {
        let threshold = self.effective_limit();
        self.heap_bytes.load(Ordering::Acquire) < threshold / 2
    }

    /// Get the effective memory limit (respecting default threshold).
    #[inline]
    #[must_use]
    fn effective_limit(&self) -> u64 {
        if self.memory_limit == 0 {
            BACKPRESSURE_THRESHOLD_BYTES
        } else {
            self.memory_limit.min(BACKPRESSURE_THRESHOLD_BYTES)
        }
    }

    /// Add heap bytes to the tracker (after pushing to reorder buffer).
    #[inline]
    pub fn add_heap_bytes(&self, bytes: u64) {
        self.heap_bytes.fetch_add(bytes, Ordering::AcqRel);
    }

    /// Subtract heap bytes from the tracker (after popping from reorder buffer).
    #[inline]
    pub fn sub_heap_bytes(&self, bytes: u64) {
        self.heap_bytes.fetch_sub(bytes, Ordering::AcqRel);
    }

    /// Update the next sequence number (after consumer advances).
    #[inline]
    pub fn update_next_seq(&self, new_seq: u64) {
        self.next_seq.store(new_seq, Ordering::Release);
    }

    /// Get the current next sequence number.
    #[inline]
    #[must_use]
    pub fn get_next_seq(&self) -> u64 {
        self.next_seq.load(Ordering::Acquire)
    }

    /// Get the current heap bytes.
    #[inline]
    #[must_use]
    pub fn get_heap_bytes(&self) -> u64 {
        self.heap_bytes.load(Ordering::Acquire)
    }

    /// Get the configured memory limit.
    #[inline]
    #[must_use]
    pub fn get_memory_limit(&self) -> u64 {
        self.memory_limit
    }
}

// ============================================================================
// BGZF Batch Buffer Types
// ============================================================================

use crate::bgzf_writer::{CompressedBlock, InlineBgzfCompressor};

/// Batch of raw BGZF blocks (input side).
///
/// This is the input buffer type for the unified pipeline when processing
/// BAM files. It holds multiple raw (compressed) BGZF blocks that will be
/// decompressed and processed by workers.
#[derive(Default)]
pub struct RawBlockBatch {
    /// The raw BGZF blocks in this batch.
    pub blocks: Vec<RawBgzfBlock>,
}

impl RawBlockBatch {
    /// Create a new empty batch.
    #[must_use]
    pub fn new() -> Self {
        Self { blocks: Vec::new() }
    }

    /// Create a batch with pre-allocated capacity.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self { blocks: Vec::with_capacity(capacity) }
    }

    /// Number of blocks in this batch.
    #[must_use]
    pub fn len(&self) -> usize {
        self.blocks.len()
    }

    /// Returns true if the batch contains no blocks.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }

    /// Total uncompressed size of all blocks.
    #[must_use]
    pub fn total_uncompressed_size(&self) -> usize {
        self.blocks.iter().map(RawBgzfBlock::uncompressed_size).sum()
    }

    /// Total compressed size of all blocks (raw bytes read from file).
    #[must_use]
    pub fn total_compressed_size(&self) -> usize {
        self.blocks.iter().map(RawBgzfBlock::len).sum()
    }

    /// Clear all blocks from this batch.
    pub fn clear(&mut self) {
        self.blocks.clear();
    }
}

impl MemoryEstimate for RawBlockBatch {
    fn estimate_heap_size(&self) -> usize {
        // Vec overhead + sum of block data sizes
        self.blocks.iter().map(|b| b.data.capacity()).sum::<usize>()
            + self.blocks.capacity() * std::mem::size_of::<RawBgzfBlock>()
    }
}

/// Batch of compressed BGZF blocks (output side).
///
/// This is the output buffer type for the unified pipeline. It holds
/// compressed blocks ready to be written to the output file.
#[derive(Default)]
pub struct CompressedBlockBatch {
    /// The compressed blocks in this batch.
    pub blocks: Vec<CompressedBlock>,
    /// Number of records/templates in this batch (for progress logging).
    pub record_count: u64,
}

impl CompressedBlockBatch {
    /// Create a new empty batch.
    #[must_use]
    pub fn new() -> Self {
        Self { blocks: Vec::new(), record_count: 0 }
    }

    /// Create a batch with pre-allocated capacity.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self { blocks: Vec::with_capacity(capacity), record_count: 0 }
    }

    /// Number of blocks in this batch.
    #[must_use]
    pub fn len(&self) -> usize {
        self.blocks.len()
    }

    /// Returns true if the batch contains no blocks.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }

    /// Total size of all compressed blocks.
    #[must_use]
    pub fn total_size(&self) -> usize {
        self.blocks.iter().map(|b| b.data.len()).sum()
    }

    /// Clear all blocks from this batch.
    pub fn clear(&mut self) {
        self.blocks.clear();
        self.record_count = 0;
    }
}

impl MemoryEstimate for CompressedBlockBatch {
    fn estimate_heap_size(&self) -> usize {
        // Vec overhead + sum of block data sizes
        self.blocks.iter().map(|b| b.data.capacity()).sum::<usize>()
            + self.blocks.capacity() * std::mem::size_of::<CompressedBlock>()
    }
}

/// Configuration for BGZF batch processing.
#[derive(Debug, Clone)]
pub struct BgzfBatchConfig {
    /// Number of raw blocks to read per batch.
    pub blocks_per_batch: usize,
    /// Compression level for output (0-12, default 6).
    pub compression_level: u32,
}

impl Default for BgzfBatchConfig {
    fn default() -> Self {
        Self { blocks_per_batch: 16, compression_level: 6 }
    }
}

impl BgzfBatchConfig {
    /// Create a new configuration with the given blocks per batch.
    #[must_use]
    pub fn new(blocks_per_batch: usize) -> Self {
        Self { blocks_per_batch, ..Default::default() }
    }

    /// Set the compression level.
    #[must_use]
    pub fn with_compression_level(mut self, level: u32) -> Self {
        self.compression_level = level;
        self
    }
}

/// Read a batch of raw BGZF blocks into a buffer.
///
/// This is the read function for the unified pipeline when processing BAM files.
/// Returns `Ok(true)` if blocks were read, `Ok(false)` at EOF.
///
/// # Errors
///
/// Returns an I/O error if reading from the underlying reader fails.
pub fn read_raw_block_batch(
    reader: &mut dyn Read,
    buffer: &mut RawBlockBatch,
    blocks_per_batch: usize,
) -> io::Result<bool> {
    buffer.clear();
    let blocks = read_raw_blocks(reader, blocks_per_batch)?;
    if blocks.is_empty() {
        return Ok(false);
    }
    buffer.blocks = blocks;
    Ok(true)
}

/// Write a batch of compressed blocks to output.
///
/// This is the write function for the unified pipeline.
///
/// # Errors
///
/// Returns an I/O error if writing to the underlying writer fails.
pub fn write_compressed_batch(
    writer: &mut dyn Write,
    batch: &CompressedBlockBatch,
) -> io::Result<()> {
    for block in &batch.blocks {
        writer.write_all(&block.data)?;
    }
    Ok(())
}

/// Per-worker state for BGZF processing.
///
/// Each worker thread has its own decompressor and compressor instances
/// to avoid synchronization overhead.
pub struct BgzfWorkerState {
    /// Decompressor for input blocks.
    pub decompressor: libdeflater::Decompressor,
    /// Compressor for output blocks.
    pub compressor: InlineBgzfCompressor,
}

impl BgzfWorkerState {
    /// Create new worker state with the given compression level.
    #[must_use]
    pub fn new(compression_level: u32) -> Self {
        Self {
            decompressor: libdeflater::Decompressor::new(),
            compressor: InlineBgzfCompressor::new(compression_level),
        }
    }

    /// Decompress a batch of raw blocks into uncompressed data.
    ///
    /// Returns the concatenated uncompressed data from all blocks.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if decompression fails.
    pub fn decompress_batch(&mut self, batch: &RawBlockBatch) -> io::Result<Vec<u8>> {
        let total_size = batch.total_uncompressed_size();
        let mut result = Vec::with_capacity(total_size);

        for block in &batch.blocks {
            decompress_block_into(block, &mut self.decompressor, &mut result)?;
        }

        Ok(result)
    }
}

// ============================================================================
// 9-Step Pipeline Infrastructure
// ============================================================================
//
// This section implements the 9-step unified pipeline where any thread
// can execute any step, but some steps are exclusive (only one thread at a time).
//
// Steps:
// 1. Read (exclusive) - Read raw BGZF blocks from input file
// 2. Decompress (parallel) - Decompress blocks using libdeflater
// 3. FindBoundaries (exclusive, FAST) - Find record boundaries in decompressed data
// 4. Decode (parallel) - Decode BAM records at known boundaries
// 5. Group (exclusive) - Group decoded records
// 6. Process (parallel) - Command-specific processing
// 7. Serialize (parallel) - Encode output records to BAM bytes
// 8. Compress (parallel) - Compress to BGZF blocks
// 9. Write (exclusive) - Write blocks to output file
//
// Queues: Q1(1→2), Q2(2→3, reorder), Q2b(3→4), Q3(4→5, reorder), Q4(5→6), Q5(6→7), Q6(7→8), Q7(8→9, reorder)

/// The 9 pipeline steps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PipelineStep {
    /// Step 1: Read raw BGZF blocks from input (exclusive - file mutex)
    Read,
    /// Step 2: Decompress blocks (parallel)
    Decompress,
    /// Step 3: Find record boundaries in decompressed data (exclusive, FAST ~0.1μs/block)
    FindBoundaries,
    /// Step 4: Decode BAM records at known boundaries (parallel)
    Decode,
    /// Step 5: Group decoded records (exclusive - grouper state)
    Group,
    /// Step 6: Process groups - command-specific work (parallel)
    Process,
    /// Step 7: Serialize output records to bytes (parallel)
    Serialize,
    /// Step 8: Compress bytes to BGZF blocks (parallel)
    Compress,
    /// Step 9: Write blocks to output file (exclusive - file mutex)
    Write,
}

impl PipelineStep {
    /// Returns true if this step requires exclusive access (only one thread at a time).
    #[must_use]
    pub const fn is_exclusive(&self) -> bool {
        matches!(self, Self::Read | Self::FindBoundaries | Self::Group | Self::Write)
    }

    /// Returns all steps in pipeline order.
    #[must_use]
    pub const fn all() -> [PipelineStep; 9] {
        [
            PipelineStep::Read,
            PipelineStep::Decompress,
            PipelineStep::FindBoundaries,
            PipelineStep::Decode,
            PipelineStep::Group,
            PipelineStep::Process,
            PipelineStep::Serialize,
            PipelineStep::Compress,
            PipelineStep::Write,
        ]
    }

    /// Convert from step index to `PipelineStep`.
    ///
    /// # Panics
    ///
    /// Panics if `idx` is greater than 8.
    #[must_use]
    pub const fn from_index(idx: usize) -> PipelineStep {
        match idx {
            0 => PipelineStep::Read,
            1 => PipelineStep::Decompress,
            2 => PipelineStep::FindBoundaries,
            3 => PipelineStep::Decode,
            4 => PipelineStep::Group,
            5 => PipelineStep::Process,
            6 => PipelineStep::Serialize,
            7 => PipelineStep::Compress,
            8 => PipelineStep::Write,
            _ => panic!("PipelineStep::from_index: invalid index (must be 0..=8)"),
        }
    }

    /// Get the 0-based index of this step (Read=0, ..., Write=8).
    #[must_use]
    pub const fn index(&self) -> usize {
        match self {
            PipelineStep::Read => 0,
            PipelineStep::Decompress => 1,
            PipelineStep::FindBoundaries => 2,
            PipelineStep::Decode => 3,
            PipelineStep::Group => 4,
            PipelineStep::Process => 5,
            PipelineStep::Serialize => 6,
            PipelineStep::Compress => 7,
            PipelineStep::Write => 8,
        }
    }

    /// Get short name for display.
    #[must_use]
    pub const fn short_name(&self) -> &'static str {
        match self {
            PipelineStep::Read => "Rd",
            PipelineStep::Decompress => "Dc",
            PipelineStep::FindBoundaries => "Fb",
            PipelineStep::Decode => "De",
            PipelineStep::Group => "Gr",
            PipelineStep::Process => "Pr",
            PipelineStep::Serialize => "Se",
            PipelineStep::Compress => "Co",
            PipelineStep::Write => "Wr",
        }
    }
}

// ============================================================================
// ActiveSteps - Configurable set of active pipeline steps
// ============================================================================

/// Tracks which pipeline steps are active for a given pipeline configuration.
///
/// Some steps may be skipped depending on input type and mode:
/// - Gzip+synchronized: skips `Decompress`, `FindBoundaries`, `Group`
/// - BGZF+synchronized: skips `Group`
/// - BAM (all active): all 9 steps
#[derive(Debug, Clone)]
pub struct ActiveSteps {
    /// Active steps in pipeline order.
    steps: Vec<PipelineStep>,
    /// Fast lookup: `active[step.index()]` is true iff step is active.
    active: [bool; 9],
}

impl ActiveSteps {
    /// Create from a list of active steps (must be in pipeline order, unique).
    ///
    /// # Panics
    ///
    /// Panics if steps are not in ascending pipeline order or contain duplicates.
    #[must_use]
    pub fn new(steps: &[PipelineStep]) -> Self {
        assert!(
            steps.windows(2).all(|w| w[0].index() < w[1].index()),
            "ActiveSteps must be unique and in pipeline order"
        );
        let mut active = [false; 9];
        for &step in steps {
            active[step.index()] = true;
        }
        Self { steps: steps.to_vec(), active }
    }

    /// All 9 steps active (BAM pipeline, default).
    #[must_use]
    pub fn all() -> Self {
        Self::new(&PipelineStep::all())
    }

    /// Check if a step is active.
    #[must_use]
    pub fn is_active(&self, step: PipelineStep) -> bool {
        self.active[step.index()]
    }

    /// Get the active steps in pipeline order.
    #[must_use]
    pub fn steps(&self) -> &[PipelineStep] {
        &self.steps
    }

    /// Get only the active exclusive steps in pipeline order.
    #[must_use]
    pub fn exclusive_steps(&self) -> Vec<PipelineStep> {
        self.steps.iter().copied().filter(PipelineStep::is_exclusive).collect()
    }

    /// Number of active steps.
    #[must_use]
    pub fn len(&self) -> usize {
        self.steps.len()
    }

    /// Returns true if no steps are active.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.steps.is_empty()
    }

    /// Filter a priority buffer in-place, keeping only active steps.
    /// Returns the number of active steps remaining.
    pub fn filter_in_place(&self, buffer: &mut [PipelineStep; 9]) -> usize {
        let mut write = 0;
        for read in 0..9 {
            if self.active[buffer[read].index()] {
                buffer[write] = buffer[read];
                write += 1;
            }
        }
        write
    }
}

// ============================================================================
// GroupKey - Pre-computed grouping key for fast comparison in Group step
// ============================================================================

/// Pre-computed grouping key for fast comparison in Group step.
///
/// All fields are integers/hashes for O(1) comparison. This is computed during
/// the parallel Decode step so the serial Group step only does integer comparisons.
///
/// For paired-end reads, positions are normalized so the lower position comes first.
/// For single-end reads, the mate fields use `UNKNOWN_*` sentinel values.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GroupKey {
    // Position info (normalized: lower position first)
    /// Reference sequence index for position 1 (lower).
    pub ref_id1: i32,
    /// Unclipped 5' position for position 1.
    pub pos1: i32,
    /// Strand for position 1 (0=forward, 1=reverse).
    pub strand1: u8,
    /// Reference sequence index for position 2 (higher or mate).
    pub ref_id2: i32,
    /// Unclipped 5' position for position 2.
    pub pos2: i32,
    /// Strand for position 2.
    pub strand2: u8,

    // Grouping metadata
    /// Library index (pre-computed from RG tag via header lookup).
    pub library_idx: u16,
    /// Hash of cell barcode (0 if none).
    pub cell_hash: u64,

    // For name-based grouping within position groups
    /// Hash of QNAME for fast name comparison.
    pub name_hash: u64,
}

impl GroupKey {
    /// Sentinel value for unknown reference ID (unpaired reads).
    pub const UNKNOWN_REF: i32 = i32::MAX;
    /// Sentinel value for unknown position (unpaired reads).
    pub const UNKNOWN_POS: i32 = i32::MAX;
    /// Sentinel value for unknown strand (unpaired reads).
    pub const UNKNOWN_STRAND: u8 = u8::MAX;

    /// Create a `GroupKey` for a paired-end read with mate info.
    ///
    /// Positions are automatically normalized so the lower position comes first.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn paired(
        ref_id: i32,
        pos: i32,
        strand: u8,
        mate_ref_id: i32,
        mate_pos: i32,
        mate_strand: u8,
        library_idx: u16,
        cell_hash: u64,
        name_hash: u64,
    ) -> Self {
        // Normalize: put lower position first (matching ReadInfo behavior)
        let (ref_id1, pos1, strand1, ref_id2, pos2, strand2) =
            if (ref_id, pos, strand) <= (mate_ref_id, mate_pos, mate_strand) {
                (ref_id, pos, strand, mate_ref_id, mate_pos, mate_strand)
            } else {
                (mate_ref_id, mate_pos, mate_strand, ref_id, pos, strand)
            };

        Self { ref_id1, pos1, strand1, ref_id2, pos2, strand2, library_idx, cell_hash, name_hash }
    }

    /// Create a `GroupKey` for a single-end/unpaired read.
    #[must_use]
    pub fn single(
        ref_id: i32,
        pos: i32,
        strand: u8,
        library_idx: u16,
        cell_hash: u64,
        name_hash: u64,
    ) -> Self {
        Self {
            ref_id1: ref_id,
            pos1: pos,
            strand1: strand,
            ref_id2: Self::UNKNOWN_REF,
            pos2: Self::UNKNOWN_POS,
            strand2: Self::UNKNOWN_STRAND,
            library_idx,
            cell_hash,
            name_hash,
        }
    }

    /// Returns the position-only key for grouping by genomic position.
    ///
    /// This is used by `RecordPositionGrouper` to determine if records belong to
    /// the same position group (ignoring name).
    #[must_use]
    pub fn position_key(&self) -> (i32, i32, u8, i32, i32, u8, u16, u64) {
        (
            self.ref_id1,
            self.pos1,
            self.strand1,
            self.ref_id2,
            self.pos2,
            self.strand2,
            self.library_idx,
            self.cell_hash,
        )
    }
}

impl PartialOrd for GroupKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for GroupKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.position_key()
            .cmp(&other.position_key())
            .then_with(|| self.name_hash.cmp(&other.name_hash))
    }
}

impl Default for GroupKey {
    fn default() -> Self {
        Self {
            ref_id1: Self::UNKNOWN_REF,
            pos1: Self::UNKNOWN_POS,
            strand1: Self::UNKNOWN_STRAND,
            ref_id2: Self::UNKNOWN_REF,
            pos2: Self::UNKNOWN_POS,
            strand2: Self::UNKNOWN_STRAND,
            library_idx: 0,
            cell_hash: 0,
            name_hash: 0,
        }
    }
}

// ============================================================================
// DecodedRecord - Record with pre-computed grouping key
// ============================================================================

/// A decoded BAM record with its pre-computed grouping key.
///
/// This is the output of the Decode step and input to the Group step.
/// The key is computed during the parallel Decode step so that the
/// serial Group step only needs to do fast integer comparisons.
///
/// Uses an enum for the record data to avoid carrying both `RecordBuf` and
/// `Vec<u8>` — saving ~24 bytes per record in parsed mode and ~200 bytes
/// in raw mode.
#[derive(Debug)]
pub struct DecodedRecord {
    /// Pre-computed grouping key.
    pub key: GroupKey,
    /// The record data — either a parsed `RecordBuf` or raw bytes.
    pub(crate) data: DecodedRecordData,
}

/// Record data: either a parsed noodles `RecordBuf` or raw BAM bytes.
#[derive(Debug)]
pub enum DecodedRecordData {
    Parsed(RecordBuf),
    Raw(Vec<u8>),
}

impl DecodedRecord {
    /// Create a new decoded record with its grouping key.
    #[must_use]
    pub fn new(record: RecordBuf, key: GroupKey) -> Self {
        Self { key, data: DecodedRecordData::Parsed(record) }
    }

    /// Create a decoded record from raw bytes, skipping noodles decode.
    #[must_use]
    pub fn from_raw_bytes(raw: Vec<u8>, key: GroupKey) -> Self {
        Self { key, data: DecodedRecordData::Raw(raw) }
    }

    /// Returns the raw bytes if this is a raw-mode record.
    #[must_use]
    pub fn raw_bytes(&self) -> Option<&[u8]> {
        match &self.data {
            DecodedRecordData::Raw(v) => Some(v),
            DecodedRecordData::Parsed(_) => None,
        }
    }

    /// Takes the raw bytes out if this is a raw-mode record.
    #[must_use]
    pub fn into_raw_bytes(self) -> Option<Vec<u8>> {
        match self.data {
            DecodedRecordData::Raw(v) => Some(v),
            DecodedRecordData::Parsed(_) => None,
        }
    }

    /// Returns a reference to the `RecordBuf` if this is a parsed-mode record.
    #[must_use]
    pub fn record(&self) -> Option<&RecordBuf> {
        match &self.data {
            DecodedRecordData::Parsed(r) => Some(r),
            DecodedRecordData::Raw(_) => None,
        }
    }

    /// Takes the `RecordBuf` out if this is a parsed-mode record.
    #[must_use]
    pub fn into_record(self) -> Option<RecordBuf> {
        match self.data {
            DecodedRecordData::Parsed(r) => Some(r),
            DecodedRecordData::Raw(_) => None,
        }
    }
}

impl MemoryEstimate for DecodedRecord {
    fn estimate_heap_size(&self) -> usize {
        match &self.data {
            DecodedRecordData::Parsed(record) => {
                crate::template::estimate_record_buf_heap_size(record)
            }
            DecodedRecordData::Raw(raw) => raw.capacity(),
        }
    }
}

impl MemoryEstimate for Vec<DecodedRecord> {
    fn estimate_heap_size(&self) -> usize {
        self.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
            + self.capacity() * std::mem::size_of::<DecodedRecord>()
    }
}

impl MemoryEstimate for RecordBuf {
    fn estimate_heap_size(&self) -> usize {
        // Delegate to the detailed estimation in template.rs (single source of truth)
        crate::template::estimate_record_buf_heap_size(self)
    }
}

impl MemoryEstimate for Vec<RecordBuf> {
    fn estimate_heap_size(&self) -> usize {
        self.iter().map(MemoryEstimate::estimate_heap_size).sum()
    }
}

impl MemoryEstimate for Vec<u8> {
    fn estimate_heap_size(&self) -> usize {
        self.capacity()
    }
}

// ============================================================================
// GroupKeyConfig - Configuration for computing `GroupKey` during Decode
// ============================================================================

/// Configuration for computing `GroupKey` during the Decode step.
///
/// When this is provided to the pipeline, the Decode step will compute
/// full `GroupKey` values for each record. This moves expensive computations
/// (CIGAR parsing, tag extraction) from the serial Group step to the parallel
/// Decode step.
#[derive(Debug, Clone)]
pub struct GroupKeyConfig {
    /// Library index for fast RG → library lookup.
    pub library_index: Arc<LibraryIndex>,
    /// Tag used for cell barcode extraction. None skips cell extraction.
    pub cell_tag: Option<Tag>,
    /// When true, skip noodles decode and work with raw BAM bytes.
    pub raw_byte_mode: bool,
}

impl GroupKeyConfig {
    /// Create a new `GroupKeyConfig`.
    #[must_use]
    pub fn new(library_index: LibraryIndex, cell_tag: Tag) -> Self {
        Self {
            library_index: Arc::new(library_index),
            cell_tag: Some(cell_tag),
            raw_byte_mode: false,
        }
    }

    /// Create a `GroupKeyConfig` for raw-byte mode.
    #[must_use]
    pub fn new_raw(library_index: LibraryIndex, cell_tag: Tag) -> Self {
        Self {
            library_index: Arc::new(library_index),
            cell_tag: Some(cell_tag),
            raw_byte_mode: true,
        }
    }

    /// Create a `GroupKeyConfig` for raw-byte mode without cell barcode extraction.
    #[must_use]
    pub fn new_raw_no_cell(library_index: LibraryIndex) -> Self {
        Self { library_index: Arc::new(library_index), cell_tag: None, raw_byte_mode: true }
    }
}

impl Default for GroupKeyConfig {
    fn default() -> Self {
        Self {
            library_index: Arc::new(LibraryIndex::default()),
            cell_tag: Some(Tag::from([b'C', b'B'])), // Default cell barcode tag
            raw_byte_mode: false,
        }
    }
}

/// Decompressed data batch (output of Step 2, input to Step 3).
#[derive(Default)]
pub struct DecompressedBatch {
    /// Concatenated decompressed bytes from multiple BGZF blocks.
    pub data: Vec<u8>,
}

impl DecompressedBatch {
    /// Create a new empty batch.
    #[must_use]
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    /// Create a batch with pre-allocated capacity.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self { data: Vec::with_capacity(capacity) }
    }

    /// Returns true if the batch contains no data.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Clear all data from this batch.
    pub fn clear(&mut self) {
        self.data.clear();
    }
}

impl MemoryEstimate for DecompressedBatch {
    fn estimate_heap_size(&self) -> usize {
        self.data.capacity()
    }
}

/// Serialized data batch (output of Step 5, input to Step 6).
#[derive(Default)]
pub struct SerializedBatch {
    /// Encoded BAM bytes ready for compression.
    pub data: Vec<u8>,
    /// Number of records/templates in this batch (for progress logging).
    pub record_count: u64,
}

impl SerializedBatch {
    /// Create a new empty batch.
    #[must_use]
    pub fn new() -> Self {
        Self { data: Vec::new(), record_count: 0 }
    }

    /// Returns true if the batch contains no data.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Clear all data from this batch.
    pub fn clear(&mut self) {
        self.data.clear();
        self.record_count = 0;
    }
}

impl MemoryEstimate for SerializedBatch {
    fn estimate_heap_size(&self) -> usize {
        self.data.capacity()
    }
}

// ============================================================================
// OutputPipelineQueues - Shared output-half state for both pipelines
// ============================================================================

/// Shared output-half state for both BAM and FASTQ pipelines.
///
/// Contains all queues and state from Group step onwards:
/// Group → Process → Serialize → Compress → Write
///
/// # Type Parameters
/// - `G`: Group/template type (input to Process step)
/// - `P`: Processed type (output of Process step)
pub struct OutputPipelineQueues<G, P: MemoryEstimate> {
    // ========== Queue: Group → Process ==========
    /// Batches of groups/templates waiting to be processed.
    pub groups: ArrayQueue<(u64, Vec<G>)>,
    /// Current heap bytes in groups queue (stats/reporting only, not used for backpressure).
    /// Mutations are gated behind the `memory-debug` feature; reads are zero without the feature.
    pub groups_heap_bytes: AtomicU64,

    // ========== Queue: Process → Serialize ==========
    /// Batches of processed data waiting for serialization.
    pub processed: ArrayQueue<(u64, Vec<P>)>,
    /// Current heap bytes in processed queue.
    pub processed_heap_bytes: AtomicU64,

    // ========== Queue: Serialize → Compress ==========
    /// Serialized bytes waiting for compression.
    pub serialized: ArrayQueue<(u64, SerializedBatch)>,
    /// Current heap bytes in serialized queue.
    pub serialized_heap_bytes: AtomicU64,

    // ========== Queue: Compress → Write (with reorder) ==========
    /// Compressed blocks waiting to be written.
    pub compressed: ArrayQueue<(u64, CompressedBlockBatch)>,
    /// Current heap bytes in compressed queue.
    pub compressed_heap_bytes: AtomicU64,
    /// Reorder buffer to ensure blocks are written in order.
    pub write_reorder: Mutex<ReorderBuffer<CompressedBlockBatch>>,

    // ========== Output ==========
    /// Output file, mutex-protected for exclusive access.
    pub output: Mutex<Option<Box<dyn Write + Send>>>,

    // ========== Completion and Error Tracking ==========
    /// Flag indicating an error occurred.
    pub error_flag: AtomicBool,
    /// Storage for the first error.
    pub error: Mutex<Option<io::Error>>,
    /// Count of items written (groups for BAM, templates for FASTQ).
    pub items_written: AtomicU64,
    /// Flag indicating pipeline is draining (input complete, finishing up).
    pub draining: AtomicBool,
    /// Progress tracker for logging.
    pub progress: ProgressTracker,

    // ========== Performance Statistics ==========
    /// Optional performance statistics collector.
    pub stats: Option<Arc<PipelineStats>>,
}

impl<G: Send, P: Send + MemoryEstimate> OutputPipelineQueues<G, P> {
    /// Create new output pipeline queues.
    ///
    /// # Arguments
    /// - `queue_capacity`: Capacity for all `ArrayQueue`s
    /// - `output`: Output writer (will be wrapped in Mutex)
    /// - `stats`: Optional performance statistics collector
    /// - `progress_name`: Name for the progress tracker (e.g., "Processed records")
    #[must_use]
    pub fn new(
        queue_capacity: usize,
        output: Box<dyn Write + Send>,
        stats: Option<Arc<PipelineStats>>,
        progress_name: &str,
    ) -> Self {
        Self {
            groups: ArrayQueue::new(queue_capacity),
            groups_heap_bytes: AtomicU64::new(0),
            processed: ArrayQueue::new(queue_capacity),
            processed_heap_bytes: AtomicU64::new(0),
            serialized: ArrayQueue::new(queue_capacity),
            serialized_heap_bytes: AtomicU64::new(0),
            compressed: ArrayQueue::new(queue_capacity),
            compressed_heap_bytes: AtomicU64::new(0),
            write_reorder: Mutex::new(ReorderBuffer::new()),
            output: Mutex::new(Some(output)),
            error_flag: AtomicBool::new(false),
            error: Mutex::new(None),
            items_written: AtomicU64::new(0),
            draining: AtomicBool::new(false),
            progress: ProgressTracker::new(progress_name).with_interval(PROGRESS_LOG_INTERVAL),
            stats,
        }
    }

    // ========== Error Handling ==========

    /// Record an error and signal threads to stop.
    pub fn set_error(&self, error: io::Error) {
        self.error_flag.store(true, Ordering::SeqCst);
        let mut guard = self.error.lock();
        if guard.is_none() {
            *guard = Some(error);
        }
    }

    /// Check if an error has occurred.
    #[must_use]
    pub fn has_error(&self) -> bool {
        self.error_flag.load(Ordering::Relaxed)
    }

    /// Take the stored error.
    pub fn take_error(&self) -> Option<io::Error> {
        self.error.lock().take()
    }

    // ========== Memory Backpressure ==========

    /// Check if processed queue memory is at the backpressure threshold.
    #[must_use]
    pub fn is_processed_memory_high(&self) -> bool {
        self.processed_heap_bytes.load(Ordering::Acquire) >= Q5_BACKPRESSURE_THRESHOLD_BYTES
    }

    /// Check if pipeline is in drain mode (bypasses memory backpressure).
    #[must_use]
    pub fn is_draining(&self) -> bool {
        self.draining.load(Ordering::Relaxed)
    }

    /// Set the draining flag.
    pub fn set_draining(&self, value: bool) {
        self.draining.store(value, Ordering::Relaxed);
    }

    /// Check backpressure for Process step (queue full OR memory high, unless draining).
    #[must_use]
    pub fn should_apply_process_backpressure(&self) -> bool {
        self.processed.is_full() || (!self.is_draining() && self.is_processed_memory_high())
    }

    // ========== Queue Depths ==========

    /// Get current lengths of all output queues.
    #[must_use]
    pub fn queue_depths(&self) -> OutputQueueDepths {
        OutputQueueDepths {
            groups: self.groups.len(),
            processed: self.processed.len(),
            serialized: self.serialized.len(),
            compressed: self.compressed.len(),
        }
    }

    /// Check if all output queues are empty.
    #[must_use]
    pub fn are_queues_empty(&self) -> bool {
        self.groups.is_empty()
            && self.processed.is_empty()
            && self.serialized.is_empty()
            && self.compressed.is_empty()
            && self.write_reorder.lock().is_empty()
    }
}

/// Snapshot of output queue depths.
#[derive(Debug, Clone, Copy)]
pub struct OutputQueueDepths {
    pub groups: usize,
    pub processed: usize,
    pub serialized: usize,
    pub compressed: usize,
}

// ============================================================================
// WorkerCoreState - Shared worker state for both pipelines
// ============================================================================

/// Minimum backoff duration in microseconds.
pub const MIN_BACKOFF_US: u64 = 10;
/// Maximum backoff duration in microseconds (1ms).
pub const MAX_BACKOFF_US: u64 = 1000;
/// Default serialization buffer capacity (256KB).
pub const SERIALIZATION_BUFFER_CAPACITY: usize = 256 * 1024;

/// Core state shared by all worker threads in both BAM and FASTQ pipelines.
///
/// This struct contains the common fields that every pipeline worker needs:
/// compressor, scheduler, serialization buffer, and adaptive backoff state.
///
/// Pipeline-specific worker states (like held items for different queue stages)
/// can embed this struct and add their own fields.
///
/// # Example
///
/// ```ignore
/// pub struct MyPipelineWorkerState<P: Send> {
///     pub core: WorkerCoreState,
///     pub held_processed: Option<(u64, Vec<P>, usize)>,
///     // ... pipeline-specific held items
/// }
/// ```
pub struct WorkerCoreState {
    /// Compressor for output BGZF blocks.
    pub compressor: InlineBgzfCompressor,
    /// Scheduler for step selection (strategy-based).
    pub scheduler: Box<dyn Scheduler>,
    /// Reusable buffer for serialization.
    /// Swapped out each batch to avoid per-batch allocation.
    pub serialization_buffer: Vec<u8>,
    /// Current backoff duration in microseconds (for adaptive backoff).
    pub backoff_us: u64,
}

impl WorkerCoreState {
    /// Create new worker core state.
    ///
    /// # Arguments
    /// * `compression_level` - BGZF compression level (0-12)
    /// * `thread_id` - Thread index (0-based)
    /// * `num_threads` - Total number of worker threads
    /// * `scheduler_strategy` - Strategy for step scheduling
    /// * `active_steps` - Which pipeline steps are active
    #[must_use]
    pub fn new(
        compression_level: u32,
        thread_id: usize,
        num_threads: usize,
        scheduler_strategy: SchedulerStrategy,
        active_steps: ActiveSteps,
    ) -> Self {
        Self {
            compressor: InlineBgzfCompressor::new(compression_level),
            scheduler: create_scheduler(scheduler_strategy, thread_id, num_threads, active_steps),
            serialization_buffer: Vec::with_capacity(SERIALIZATION_BUFFER_CAPACITY),
            backoff_us: MIN_BACKOFF_US,
        }
    }

    /// Reset backoff to minimum (after successful work).
    #[inline]
    pub fn reset_backoff(&mut self) {
        self.backoff_us = MIN_BACKOFF_US;
    }

    /// Increase backoff exponentially (after no work done).
    #[inline]
    pub fn increase_backoff(&mut self) {
        self.backoff_us = (self.backoff_us * 2).min(MAX_BACKOFF_US);
    }

    /// Sleep for the current backoff duration with jitter.
    ///
    /// Uses `yield_now()` for minimum backoff to avoid sleep syscall overhead.
    /// Adds jitter (±25%) to prevent thundering herd synchronization.
    #[inline]
    pub fn sleep_backoff(&self) {
        if self.backoff_us <= MIN_BACKOFF_US {
            // At minimum backoff, yield is cheaper than sleep
            std::thread::yield_now();
        } else {
            // Add jitter to prevent thundering herd (±25%)
            // Use thread ID hash + time for simple pseudo-random jitter
            let jitter_range = self.backoff_us / 4;
            let jitter_seed = u64::from(
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map(|d| d.subsec_nanos())
                    .unwrap_or(0),
            );
            // Simple hash: mix in some bits to vary across threads
            let jitter_offset = (jitter_seed % (jitter_range * 2)).saturating_sub(jitter_range);
            let actual_us = self.backoff_us.saturating_add(jitter_offset).max(MIN_BACKOFF_US);
            std::thread::sleep(std::time::Duration::from_micros(actual_us));
        }
    }
}

impl HasCompressor for WorkerCoreState {
    fn compressor_mut(&mut self) -> &mut InlineBgzfCompressor {
        &mut self.compressor
    }
}

// ============================================================================
// Pipeline Validation - Post-shutdown data loss detection
// ============================================================================

/// Error returned when pipeline validation detects data loss or inconsistency.
///
/// This error provides detailed diagnostics about what went wrong, including
/// which queues still contain data, which counters don't match, and how much
/// memory is leaked in tracking.
#[derive(Debug, Clone)]
pub struct PipelineValidationError {
    /// Names of queues that are not empty at completion.
    pub non_empty_queues: Vec<String>,
    /// Descriptions of counter mismatches.
    pub counter_mismatches: Vec<String>,
    /// Total heap bytes still tracked (should be 0 at completion).
    pub leaked_heap_bytes: u64,
}

impl std::fmt::Display for PipelineValidationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Pipeline validation failed - potential data loss detected:")?;
        if !self.non_empty_queues.is_empty() {
            writeln!(f, "  Non-empty queues: {}", self.non_empty_queues.join(", "))?;
        }
        if !self.counter_mismatches.is_empty() {
            writeln!(f, "  Counter mismatches:")?;
            for mismatch in &self.counter_mismatches {
                writeln!(f, "    - {mismatch}")?;
            }
        }
        if self.leaked_heap_bytes > 0 {
            writeln!(f, "  Leaked heap bytes: {}", self.leaked_heap_bytes)?;
        }
        Ok(())
    }
}

impl std::error::Error for PipelineValidationError {}

// ============================================================================
// PipelineLifecycle Trait - Common patterns for pipeline state management
// ============================================================================

/// Trait for common pipeline lifecycle operations.
///
/// This trait captures the patterns shared between BAM and FASTQ pipelines:
/// - Completion checking
/// - Error handling
/// - Drain mode management
/// - Output queue access
///
/// Implementing this trait allows sharing monitor thread logic, worker loop
/// exit conditions, and pipeline finalization code.
pub trait PipelineLifecycle {
    /// Check if the pipeline has completed all work.
    ///
    /// Returns true when:
    /// 1. Input is fully consumed (`read_done`, `group_done`)
    /// 2. All queues are empty
    /// 3. All reorder buffers are empty
    fn is_complete(&self) -> bool;

    /// Check if an error has occurred.
    fn has_error(&self) -> bool;

    /// Take the stored error (if any).
    fn take_error(&self) -> Option<io::Error>;

    /// Set an error on the pipeline.
    fn set_error(&self, error: io::Error);

    /// Check if the pipeline is in drain mode.
    ///
    /// When draining, memory-based backpressure is bypassed to prevent
    /// deadlock during pipeline completion.
    fn is_draining(&self) -> bool;

    /// Set the draining flag.
    fn set_draining(&self, value: bool);

    /// Get reference to optional pipeline statistics.
    fn stats(&self) -> Option<&PipelineStats>;

    /// Get reference to the progress tracker.
    fn progress(&self) -> &ProgressTracker;

    /// Get the number of items written.
    fn items_written(&self) -> u64;

    /// Flush the output writer and finalize.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if flushing fails.
    fn flush_output(&self) -> io::Result<()>;

    /// Validate pipeline completion to detect data loss.
    ///
    /// This method should be called after the pipeline has completed to verify:
    /// 1. All queues are empty (no data stuck in transit)
    /// 2. All batch counters match (no batches lost between stages)
    /// 3. Internal buffers are empty (grouper state, boundary leftovers)
    ///
    /// Note: Heap byte tracking is reported in `PipelineValidationError::leaked_heap_bytes`
    /// but is currently advisory only (implementations set it to 0) because estimation
    /// can be imprecise. Only queue emptiness and counter checks cause validation failure.
    ///
    /// Returns `Ok(())` if validation passes, or `Err(PipelineValidationError)`
    /// with detailed diagnostics if any issues are detected.
    ///
    /// # Errors
    ///
    /// Returns `PipelineValidationError` if validation detects data loss or inconsistency.
    fn validate_completion(&self) -> Result<(), PipelineValidationError>;
}

/// Monitor thread exit condition - returns true when pipeline should stop.
///
/// This helper function encapsulates the common exit condition check for
/// monitor threads in both pipelines.
#[inline]
pub fn should_monitor_exit<P: PipelineLifecycle>(pipeline: &P) -> bool {
    pipeline.is_complete() || pipeline.has_error()
}

/// Trait for states that support monitor thread operations.
///
/// This trait extends `PipelineLifecycle` with deadlock detection capabilities,
/// allowing shared monitor thread logic between BAM and FASTQ pipelines.
pub trait MonitorableState: PipelineLifecycle {
    /// Get reference to the deadlock state.
    fn deadlock_state(&self) -> &DeadlockState;

    /// Build a queue snapshot for deadlock detection.
    ///
    /// This method captures the current state of all queues and reorder buffers
    /// for deadlock analysis. The implementation is pipeline-specific since
    /// BAM and FASTQ have different queue structures.
    fn build_queue_snapshot(&self) -> QueueSnapshot;
}

/// Run monitor loop until pipeline completes or errors.
///
/// This function provides the common monitor thread pattern:
/// 1. Sleep for the sample interval
/// 2. Check exit condition (complete or error)
/// 3. Call the pipeline-specific sampling callback
/// 4. Periodically perform deadlock detection
///
/// # Arguments
/// - `state`: Arc-wrapped pipeline state implementing `MonitorableState`
/// - `sample_interval_ms`: Milliseconds between samples (typically 100)
/// - `deadlock_check_samples`: Number of samples between deadlock checks (10 = ~1 second)
/// - `on_sample`: Callback for pipeline-specific per-sample work (stats collection, logging)
///
/// # Example
/// ```ignore
/// run_monitor_loop(&state, 100, 10, |s| {
///     if let Some(ref stats) = s.stats() {
///         stats.add_queue_sample(...);
///     }
/// });
/// ```
pub fn run_monitor_loop<S, F>(
    state: &Arc<S>,
    sample_interval_ms: u64,
    deadlock_check_samples: u32,
    on_sample: F,
) where
    S: MonitorableState,
    F: Fn(&S),
{
    let mut deadlock_counter = 0u32;
    loop {
        thread::sleep(Duration::from_millis(sample_interval_ms));

        if should_monitor_exit(state.as_ref()) {
            break;
        }

        // Call pipeline-specific sampling function
        on_sample(state.as_ref());

        // Deadlock check at specified interval
        if state.deadlock_state().is_enabled() {
            deadlock_counter += 1;
            if deadlock_counter >= deadlock_check_samples {
                deadlock_counter = 0;
                let snapshot = state.build_queue_snapshot();
                check_deadlock_and_restore(state.deadlock_state(), &snapshot);
            }
        }
    }
}

// ============================================================================
// Panic Handling Helpers
// ============================================================================

/// Extract a human-readable message from a panic payload.
///
/// This is used to provide useful error messages when worker threads panic.
/// It handles the common cases of &str and String panic payloads, with a
/// fallback for other types.
#[must_use]
#[allow(clippy::needless_pass_by_value)]
pub fn extract_panic_message(panic_info: Box<dyn std::any::Any + Send>) -> String {
    if let Some(s) = panic_info.downcast_ref::<&str>() {
        (*s).to_string()
    } else if let Some(s) = panic_info.downcast_ref::<String>() {
        s.clone()
    } else {
        "Unknown panic".to_string()
    }
}

/// Handle a worker thread panic by setting an error on the pipeline state.
///
/// This combines panic message extraction with error setting for cleaner
/// worker thread spawning code.
pub fn handle_worker_panic<S: PipelineLifecycle>(
    state: &S,
    thread_id: usize,
    panic_info: Box<dyn std::any::Any + Send>,
) {
    let msg = extract_panic_message(panic_info);
    state.set_error(io::Error::other(format!("Worker thread {thread_id} panicked: {msg}")));
}

// ============================================================================
// Pipeline Finalization Helpers
// ============================================================================

/// Join all worker threads, waiting for completion.
///
/// Returns an error if any thread panicked without setting an error on the state.
///
/// # Errors
///
/// Returns an I/O error if any worker thread panicked.
pub fn join_worker_threads(handles: Vec<thread::JoinHandle<()>>) -> io::Result<()> {
    for handle in handles {
        handle.join().map_err(|_| io::Error::other("Worker thread panicked"))?;
    }
    Ok(())
}

/// Join the monitor thread if present.
pub fn join_monitor_thread(handle: Option<thread::JoinHandle<()>>) {
    if let Some(h) = handle {
        let _ = h.join();
    }
}

/// Finalize a pipeline after all threads have completed.
///
/// This function:
/// 1. Checks for errors that occurred during processing
/// 2. Flushes the output writer
/// 3. Logs pipeline statistics if enabled
/// 4. Returns the count of items written
///
/// # Arguments
/// - `state`: The pipeline state implementing `PipelineLifecycle`
///
/// # Returns
/// - `Ok(items_written)` on success
/// - `Err(error)` if an error occurred during processing
///
/// # Errors
///
/// Returns an I/O error if an error occurred during processing or output flush fails.
pub fn finalize_pipeline<S: PipelineLifecycle>(state: &S) -> io::Result<u64> {
    // Check for errors
    if let Some(error) = state.take_error() {
        return Err(error);
    }

    // Validate pipeline completion to detect data loss
    state.validate_completion().map_err(io::Error::other)?;

    // Flush output
    state.flush_output()?;

    // Log pipeline statistics if enabled
    if let Some(stats) = state.stats() {
        stats.log_summary();
    }

    Ok(state.items_written())
}

// ============================================================================
// Trait Implementations for OutputPipelineQueues
// ============================================================================

impl<G: Send + 'static, P: Send + MemoryEstimate + 'static> ProcessPipelineState<G, P>
    for OutputPipelineQueues<G, P>
{
    fn process_input_pop(&self) -> Option<(u64, Vec<G>)> {
        self.groups.pop()
    }

    fn process_output_is_full(&self) -> bool {
        self.processed.is_full()
    }

    fn process_output_push(&self, item: (u64, Vec<P>)) -> Result<(), (u64, Vec<P>)> {
        let heap_size: usize = item.1.iter().map(MemoryEstimate::estimate_heap_size).sum();
        let result = self.processed.push(item);
        if result.is_ok() {
            self.processed_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
        }
        result
    }

    fn has_error(&self) -> bool {
        self.error_flag.load(Ordering::Acquire)
    }

    fn set_error(&self, error: io::Error) {
        OutputPipelineQueues::set_error(self, error);
    }

    fn should_apply_process_backpressure(&self) -> bool {
        OutputPipelineQueues::should_apply_process_backpressure(self)
    }

    fn is_draining(&self) -> bool {
        OutputPipelineQueues::is_draining(self)
    }
}

impl<G: Send + 'static, P: Send + MemoryEstimate + 'static> SerializePipelineState<P>
    for OutputPipelineQueues<G, P>
{
    fn serialize_input_pop(&self) -> Option<(u64, Vec<P>)> {
        let result = self.processed.pop();
        if let Some((_, ref batch)) = result {
            let heap_size: usize = batch.iter().map(MemoryEstimate::estimate_heap_size).sum();
            self.processed_heap_bytes.fetch_sub(heap_size as u64, Ordering::AcqRel);
        }
        result
    }

    fn serialize_output_is_full(&self) -> bool {
        self.serialized.is_full()
    }

    fn serialize_output_push(
        &self,
        item: (u64, SerializedBatch),
    ) -> Result<(), (u64, SerializedBatch)> {
        let heap_size = item.1.estimate_heap_size();
        let result = self.serialized.push(item);
        if result.is_ok() {
            self.serialized_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
        }
        result
    }

    fn has_error(&self) -> bool {
        self.error_flag.load(Ordering::Acquire)
    }

    fn set_error(&self, error: io::Error) {
        OutputPipelineQueues::set_error(self, error);
    }

    fn record_serialized_bytes(&self, bytes: u64) {
        if let Some(ref stats) = self.stats {
            stats.serialized_bytes.fetch_add(bytes, Ordering::Relaxed);
        }
    }
}

impl<G: Send + 'static, P: Send + MemoryEstimate + 'static> WritePipelineState
    for OutputPipelineQueues<G, P>
{
    fn write_input_queue(&self) -> &ArrayQueue<(u64, CompressedBlockBatch)> {
        &self.compressed
    }

    fn write_reorder_buffer(&self) -> &Mutex<ReorderBuffer<CompressedBlockBatch>> {
        &self.write_reorder
    }

    fn write_output(&self) -> &Mutex<Option<Box<dyn Write + Send>>> {
        &self.output
    }

    fn has_error(&self) -> bool {
        self.error_flag.load(Ordering::Acquire)
    }

    fn set_error(&self, error: io::Error) {
        OutputPipelineQueues::set_error(self, error);
    }

    fn record_written(&self, count: u64) {
        self.items_written.fetch_add(count, Ordering::Release);
    }

    fn stats(&self) -> Option<&PipelineStats> {
        self.stats.as_deref()
    }
}

impl<G: Send + 'static, P: Send + MemoryEstimate + 'static> OutputPipelineState
    for OutputPipelineQueues<G, P>
{
    type Processed = P;

    fn has_error(&self) -> bool {
        self.error_flag.load(Ordering::Acquire)
    }

    fn set_error(&self, error: io::Error) {
        OutputPipelineQueues::set_error(self, error);
    }

    fn q5_pop(&self) -> Option<(u64, SerializedBatch)> {
        self.serialized.pop()
    }

    fn q5_push(&self, item: (u64, SerializedBatch)) -> Result<(), (u64, SerializedBatch)> {
        self.serialized.push(item)
    }

    fn q5_is_full(&self) -> bool {
        self.serialized.is_full()
    }

    fn q5_track_pop(&self, heap_size: u64) {
        self.serialized_heap_bytes.fetch_sub(heap_size, Ordering::AcqRel);
    }

    fn q6_pop(&self) -> Option<(u64, CompressedBlockBatch)> {
        self.compressed.pop()
    }

    fn q6_push(
        &self,
        item: (u64, CompressedBlockBatch),
    ) -> Result<(), (u64, CompressedBlockBatch)> {
        self.compressed.push(item)
    }

    fn q6_is_full(&self) -> bool {
        self.compressed.is_full()
    }

    fn q6_reorder_insert(&self, serial: u64, batch: CompressedBlockBatch) {
        self.write_reorder.lock().insert(serial, batch);
    }

    fn q6_reorder_try_pop_next(&self) -> Option<CompressedBlockBatch> {
        self.write_reorder.lock().try_pop_next()
    }

    fn output_try_lock(
        &self,
    ) -> Option<parking_lot::MutexGuard<'_, Option<Box<dyn Write + Send>>>> {
        self.output.try_lock()
    }

    fn increment_written(&self) -> u64 {
        self.items_written.fetch_add(1, Ordering::Release)
    }

    fn record_compressed_bytes_out(&self, bytes: u64) {
        if let Some(ref stats) = self.stats {
            stats.compressed_bytes_out.fetch_add(bytes, Ordering::Relaxed);
        }
    }

    fn record_q6_pop_progress(&self) {
        // Deadlock tracking handled by caller if needed
    }

    fn record_q7_push_progress(&self) {
        // Deadlock tracking handled by caller if needed
    }

    fn stats(&self) -> Option<&PipelineStats> {
        self.stats.as_deref()
    }
}

/// Unit type always has zero heap size (used in tests).
impl MemoryEstimate for () {
    fn estimate_heap_size(&self) -> usize {
        0
    }
}

/// Configuration for the BAM pipeline.
#[derive(Debug, Clone)]
pub struct PipelineConfig {
    /// Number of threads in the pool.
    pub num_threads: usize,
    /// Capacity for each queue.
    pub queue_capacity: usize,
    /// Low water mark - prioritize reading when Q1 is below this.
    pub input_low_water: usize,
    /// High water mark - prioritize writing when Q6 is above this.
    pub output_high_water: usize,
    /// Compression level for BGZF output (0-12).
    pub compression_level: u32,
    /// Number of raw BGZF blocks to read per batch.
    pub blocks_per_read_batch: usize,
    /// Whether to collect pipeline performance statistics.
    pub collect_stats: bool,
    /// Number of groups to batch before processing (1 = no batching).
    ///
    /// This is critical for performance when work units are small. Each work unit
    /// goes through compress → write individually, and the compress step flushes
    /// after each unit, creating BGZF blocks. Small work units (~100-500 bytes)
    /// create tiny blocks instead of optimal 64KB blocks.
    ///
    /// - **`batch_size=1`**: For commands with naturally large work units
    ///   (e.g., `group` with `RawPositionGroups` containing many records).
    /// - **`batch_size>1`**: For commands with small work units
    ///   (e.g., `filter` with single records). Use ~400 for ~60KB batches.
    pub batch_size: usize,
    /// Target templates per batch for weight-based batching.
    ///
    /// When set to non-zero, groups are batched by total template count rather
    /// than group count. This provides consistent batch sizes across datasets
    /// with varying templates-per-group ratios.
    ///
    /// - **`target_templates_per_batch=0`**: Use `batch_size` for group-count batching.
    /// - **`target_templates_per_batch>0`**: Accumulate groups until total templates >= this value.
    ///
    /// Recommended value: 500 templates per batch for good performance.
    pub target_templates_per_batch: usize,
    /// Whether the BAM header has already been read from the input stream.
    ///
    /// When true, the boundary finder will skip header parsing since the header
    /// has already been consumed by the caller (e.g., for streaming from stdin).
    /// This enables streaming input support where the header must be read once
    /// and passed separately.
    pub header_already_read: bool,
    /// Scheduler strategy for thread work assignment.
    ///
    /// - `FixedPriority` (default): Thread roles are fixed (reader, writer, workers).
    /// - `ChaseBottleneck`: Threads dynamically follow work through the pipeline.
    pub scheduler_strategy: SchedulerStrategy,
    /// Memory limit for queue contents in bytes. 0 = no limit.
    ///
    /// When set, the pipeline will apply backpressure (refuse to push to queues)
    /// when total queue memory exceeds this limit. This prevents unbounded memory
    /// growth in threaded pipelines processing large datasets.
    ///
    /// Recommended values:
    /// - 0: No limit (default, uses count-based queue capacity only)
    /// - `500_000_000`: 500 MB limit
    /// - `1_000_000_000`: 1 GB limit
    /// - `2_000_000_000`: 2 GB limit
    pub queue_memory_limit: u64,
    /// Deadlock detection timeout in seconds (default 10, 0 = disabled).
    ///
    /// When no progress is made for this duration, a warning is logged with
    /// diagnostic info (queue depths, memory usage, per-queue timestamps).
    pub deadlock_timeout_secs: u64,
    /// Whether automatic deadlock recovery is enabled (default false).
    ///
    /// When enabled, uses progressive doubling: 2x -> 4x -> 8x -> unbind,
    /// then restores toward original limits after 30s of sustained progress.
    pub deadlock_recover_enabled: bool,
    /// Shared statistics instance for external memory monitoring.
    /// When provided, the pipeline will use this instead of creating its own.
    pub shared_stats: Option<Arc<PipelineStats>>,
}

impl PipelineConfig {
    /// Create a new configuration with the specified thread count and compression level.
    #[must_use]
    pub fn new(num_threads: usize, compression_level: u32) -> Self {
        Self {
            num_threads: num_threads.max(1),
            queue_capacity: 64,
            input_low_water: 8,
            output_high_water: 32,
            compression_level,
            blocks_per_read_batch: 16,
            collect_stats: false,
            batch_size: 1,
            target_templates_per_batch: 0, // 0 = use batch_size instead
            header_already_read: false,
            scheduler_strategy: SchedulerStrategy::default(),
            queue_memory_limit: 0,           // No limit by default
            deadlock_timeout_secs: 10,       // Default 10 second timeout
            deadlock_recover_enabled: false, // Detection only by default
            shared_stats: None,              // No shared stats by default
        }
    }

    /// Set the compression level.
    #[must_use]
    pub fn with_compression_level(mut self, level: u32) -> Self {
        self.compression_level = level;
        self
    }

    /// Set blocks per read batch.
    #[must_use]
    pub fn with_blocks_per_batch(mut self, blocks: usize) -> Self {
        self.blocks_per_read_batch = blocks;
        self
    }

    /// Enable or disable performance statistics collection.
    #[must_use]
    pub fn with_stats(mut self, collect: bool) -> Self {
        self.collect_stats = collect;
        self
    }

    /// Set custom statistics instance for memory debugging.
    /// This allows external monitoring to access the same stats used by the pipeline.
    #[must_use]
    pub fn with_shared_stats(mut self, stats: Arc<PipelineStats>) -> Self {
        self.collect_stats = true; // Enable stats collection
        self.shared_stats = Some(stats);
        self
    }

    /// Set the batch size for grouping work units before processing.
    #[must_use]
    pub fn with_batch_size(mut self, size: usize) -> Self {
        self.batch_size = size.max(1);
        self
    }

    /// Set the target templates per batch for weight-based batching.
    ///
    /// When set to non-zero, groups are batched by total template count rather
    /// than group count. Set to 0 to use `batch_size` instead.
    #[must_use]
    pub fn with_target_templates_per_batch(mut self, count: usize) -> Self {
        self.target_templates_per_batch = count;
        self
    }

    /// Create a configuration auto-tuned for the given thread count.
    ///
    /// This adjusts queue capacity and batch sizes based on the number of threads
    /// to optimize throughput and reduce contention.
    #[must_use]
    pub fn auto_tuned(num_threads: usize, compression_level: u32) -> Self {
        let num_threads = num_threads.max(1);

        // Scale queue capacity with thread count (min 64, max 256)
        let queue_capacity = (num_threads * 16).clamp(64, 256);

        // Low water mark: trigger read refill when below thread count
        let input_low_water = num_threads.max(4);

        // High water mark: trigger write when output is above 4x threads
        let output_high_water = (num_threads * 4).max(32);

        // Larger batches for more threads (amortize syscall overhead)
        let blocks_per_read_batch = if num_threads >= 8 { 32 } else { 16 };

        Self {
            num_threads,
            queue_capacity,
            input_low_water,
            output_high_water,
            compression_level,
            blocks_per_read_batch,
            collect_stats: false,
            batch_size: 1,
            // Template-based batching: accumulate groups until ~500 templates
            // This provides consistent batch sizes regardless of templates-per-group ratio
            target_templates_per_batch: 500,
            header_already_read: false,
            scheduler_strategy: SchedulerStrategy::default(),
            queue_memory_limit: 0,           // No limit by default
            deadlock_timeout_secs: 10,       // Default 10 second timeout
            deadlock_recover_enabled: false, // Detection only by default
            shared_stats: None,              // No shared stats by default
        }
    }

    /// Set the scheduler strategy.
    #[must_use]
    pub fn with_scheduler_strategy(mut self, strategy: SchedulerStrategy) -> Self {
        self.scheduler_strategy = strategy;
        self
    }

    /// Set the queue memory limit.
    ///
    /// When set to a non-zero value, the pipeline will apply backpressure
    /// (refuse to push to queues) when total queue memory exceeds this limit.
    ///
    /// # Arguments
    /// * `limit_bytes` - Maximum allowed queue memory in bytes. Use 0 for no limit.
    #[must_use]
    pub fn with_queue_memory_limit(mut self, limit_bytes: u64) -> Self {
        self.queue_memory_limit = limit_bytes;
        self
    }

    /// Set the deadlock detection timeout.
    ///
    /// # Arguments
    /// * `timeout_secs` - Timeout in seconds (0 = disabled, default 10)
    #[must_use]
    pub fn with_deadlock_timeout(mut self, timeout_secs: u64) -> Self {
        self.deadlock_timeout_secs = timeout_secs;
        self
    }

    /// Enable or disable deadlock recovery.
    ///
    /// When enabled, the pipeline will double memory limits when a deadlock
    /// is detected, and restore toward original limits after sustained progress.
    #[must_use]
    pub fn with_deadlock_recovery(mut self, enabled: bool) -> Self {
        self.deadlock_recover_enabled = enabled;
        self
    }
}

// ============================================================================
// Pipeline Statistics
// ============================================================================

/// Maximum number of threads supported for per-thread statistics.
pub const MAX_THREADS: usize = 32;

/// Number of pipeline steps for array sizing.
const NUM_STEPS: usize = 9;

/// A snapshot of queue sizes at a point in time.
#[derive(Debug, Clone)]
pub struct QueueSample {
    /// Time since pipeline start in milliseconds.
    pub time_ms: u64,
    /// Size of each queue: `[Q1, Q2, Q2b, Q3, Q4, Q5, Q6, Q7]`.
    pub queue_sizes: [usize; 8],
    /// Size of reorder buffers: `[Q2_reorder, Q3_reorder, Q7_reorder]`.
    pub reorder_sizes: [usize; 3],
    /// Estimated memory usage per queue in bytes: `[Q1, Q2, Q2b, Q3, Q4, Q5, Q6, Q7]`.
    /// Zero means not measured.
    pub queue_memory_bytes: [u64; 8],
    /// Estimated memory in reorder buffers in bytes: `[Q2_reorder, Q3_reorder, Q7_reorder]`.
    pub reorder_memory_bytes: [u64; 3],
    /// Current step each thread is on (0-8, or 255 for idle).
    pub thread_steps: Vec<u8>,
}

/// Statistics collected during pipeline execution.
///
/// These metrics help identify bottlenecks and optimization opportunities.
/// All counters are atomic to allow lock-free updates from multiple threads.
#[derive(Debug)]
pub struct PipelineStats {
    // Per-step timing (nanoseconds)
    /// Total time spent in Step 1: Read
    pub step_read_ns: AtomicU64,
    /// Total time spent in Step 2: Decompress
    pub step_decompress_ns: AtomicU64,
    /// Total time spent in Step 3: `FindBoundaries`
    pub step_find_boundaries_ns: AtomicU64,
    /// Total time spent in Step 4: Decode
    pub step_decode_ns: AtomicU64,
    /// Total time spent in Step 5: Group
    pub step_group_ns: AtomicU64,
    /// Total time spent in Step 6: Process
    pub step_process_ns: AtomicU64,
    /// Total time spent in Step 7: Serialize
    pub step_serialize_ns: AtomicU64,
    /// Total time spent in Step 8: Compress
    pub step_compress_ns: AtomicU64,
    /// Total time spent in Step 9: Write
    pub step_write_ns: AtomicU64,

    // Per-step success counts
    /// Number of successful Read operations
    pub step_read_count: AtomicU64,
    /// Number of successful Decompress operations
    pub step_decompress_count: AtomicU64,
    /// Number of successful `FindBoundaries` operations
    pub step_find_boundaries_count: AtomicU64,
    /// Number of successful Decode operations
    pub step_decode_count: AtomicU64,
    /// Number of successful Group operations
    pub step_group_count: AtomicU64,
    /// Number of successful Process operations
    pub step_process_count: AtomicU64,
    /// Number of successful Serialize operations
    pub step_serialize_count: AtomicU64,
    /// Number of successful Compress operations
    pub step_compress_count: AtomicU64,
    /// Number of successful Write operations
    pub step_write_count: AtomicU64,

    // Contention metrics
    /// Number of failed `try_lock` attempts on input file
    pub read_contention: AtomicU64,
    /// Number of failed `try_lock` attempts on boundary state
    pub boundary_contention: AtomicU64,
    /// Number of failed `try_lock` attempts on group state
    pub group_contention: AtomicU64,
    /// Number of failed `try_lock` attempts on output file
    pub write_contention: AtomicU64,

    // Queue metrics
    /// Number of times Q1 was empty when trying to decompress
    pub q1_empty: AtomicU64,
    /// Number of times Q2 was empty when trying to find boundaries
    pub q2_empty: AtomicU64,
    /// Number of times Q2b was empty when trying to decode
    pub q2b_empty: AtomicU64,
    /// Number of times Q3 was empty when trying to group
    pub q3_empty: AtomicU64,
    /// Number of times Q4 was empty when trying to process
    pub q4_empty: AtomicU64,
    /// Number of times Q5 was empty when trying to serialize
    pub q5_empty: AtomicU64,
    /// Number of times Q6 was empty when trying to compress
    pub q6_empty: AtomicU64,
    /// Number of times Q7 was empty when trying to write
    pub q7_empty: AtomicU64,

    // Idle tracking
    /// Number of `yield_now()` calls (no work available)
    pub idle_yields: AtomicU64,

    // ========================================================================
    // NEW: Per-thread statistics for bottleneck analysis
    // ========================================================================
    /// Per-thread step counts (successes): `[thread_id][step_index]`
    /// Tracks how many times each thread successfully executed each step.
    pub per_thread_step_counts: Box<[[AtomicU64; NUM_STEPS]; MAX_THREADS]>,

    /// Per-thread step attempts: `[thread_id][step_index]`
    /// Tracks how many times each thread attempted each step (success + failure).
    pub per_thread_step_attempts: Box<[[AtomicU64; NUM_STEPS]; MAX_THREADS]>,

    /// Per-thread idle time in nanoseconds.
    /// Time spent in `yield_now()` when no work was available.
    pub per_thread_idle_ns: Box<[AtomicU64; MAX_THREADS]>,

    /// Current step each thread is working on (for activity snapshots).
    /// Value is step index (0-8) or 255 for idle.
    pub per_thread_current_step: Box<[AtomicU8; MAX_THREADS]>,

    // Queue wait time (time spent waiting for exclusive step locks)
    /// Total time spent waiting for `FindBoundaries` lock (nanoseconds).
    pub boundary_wait_ns: AtomicU64,
    /// Total time spent waiting for Group lock (nanoseconds).
    pub group_wait_ns: AtomicU64,
    /// Total time spent waiting for Write lock (nanoseconds).
    pub write_wait_ns: AtomicU64,

    // Batch size tracking
    /// Sum of all batch sizes (for computing average).
    pub batch_size_sum: AtomicU64,
    /// Number of batches processed.
    pub batch_count: AtomicU64,
    /// Minimum batch size seen.
    pub batch_size_min: AtomicU64,
    /// Maximum batch size seen.
    pub batch_size_max: AtomicU64,

    /// Number of threads configured (for display).
    pub num_threads: AtomicU64,

    // ========================================================================
    // Throughput metrics (bytes/records processed per step)
    // ========================================================================
    /// Total bytes read from input file.
    pub bytes_read: AtomicU64,
    /// Total bytes written to output file.
    pub bytes_written: AtomicU64,
    /// Compressed bytes input to decompress step.
    pub compressed_bytes_in: AtomicU64,
    /// Decompressed bytes output from decompress step.
    pub decompressed_bytes: AtomicU64,
    /// Serialized bytes input to compress step.
    pub serialized_bytes: AtomicU64,
    /// Compressed bytes output from compress step.
    pub compressed_bytes_out: AtomicU64,
    /// Number of records decoded.
    pub records_decoded: AtomicU64,
    /// Number of groups produced by grouper.
    pub groups_produced: AtomicU64,

    // ========================================================================
    // Queue monitoring samples (periodic snapshots of queue sizes)
    // ========================================================================
    /// Collected queue size samples for timeline analysis.
    /// Protected by Mutex since samples are collected periodically from monitor thread.
    pub queue_samples: Mutex<Vec<QueueSample>>,

    // ========================================================================
    // Memory limiting statistics
    // ========================================================================
    /// Number of times memory drain mode was activated.
    pub memory_drain_activations: AtomicU64,
    /// Number of times Group step was rejected due to memory limit.
    pub group_memory_rejects: AtomicU64,
    /// Peak memory usage in the reorder buffer (bytes).
    pub peak_memory_bytes: AtomicU64,

    // ========================================================================
    // COMPREHENSIVE MEMORY TRACKING SYSTEM (behind memory-debug feature)
    // ========================================================================
    /// Comprehensive memory tracking stats, only present with `memory-debug` feature.
    #[cfg(feature = "memory-debug")]
    pub memory: MemoryDebugStats,
}

/// Helper to create a boxed array of `AtomicU64` initialized to zero.
/// We use Box to avoid stack overflow for large arrays.
#[allow(clippy::unnecessary_box_returns)]
fn new_atomic_array<const N: usize>() -> Box<[AtomicU64; N]> {
    // Create a Vec and convert to boxed slice, then try_into boxed array
    let v: Vec<AtomicU64> = (0..N).map(|_| AtomicU64::new(0)).collect();
    v.into_boxed_slice().try_into().unwrap()
}

/// Helper to create a boxed 2D array of `AtomicU64` initialized to zero.
/// We use Box to avoid stack overflow for large arrays.
#[allow(clippy::unnecessary_box_returns)]
fn new_atomic_2d_array<const R: usize, const C: usize>() -> Box<[[AtomicU64; C]; R]> {
    let v: Vec<[AtomicU64; C]> =
        (0..R).map(|_| std::array::from_fn(|_| AtomicU64::new(0))).collect();
    v.into_boxed_slice().try_into().unwrap()
}

/// Helper to create a boxed array of `AtomicU8` initialized to a value.
#[allow(clippy::unnecessary_box_returns)]
fn new_atomic_u8_array<const N: usize>(init: u8) -> Box<[AtomicU8; N]> {
    let v: Vec<AtomicU8> = (0..N).map(|_| AtomicU8::new(init)).collect();
    v.into_boxed_slice().try_into().unwrap()
}

impl Default for PipelineStats {
    fn default() -> Self {
        Self::new()
    }
}

impl PipelineStats {
    /// Create a new stats collector.
    #[must_use]
    pub fn new() -> Self {
        Self {
            step_read_ns: AtomicU64::new(0),
            step_decompress_ns: AtomicU64::new(0),
            step_find_boundaries_ns: AtomicU64::new(0),
            step_decode_ns: AtomicU64::new(0),
            step_group_ns: AtomicU64::new(0),
            step_process_ns: AtomicU64::new(0),
            step_serialize_ns: AtomicU64::new(0),
            step_compress_ns: AtomicU64::new(0),
            step_write_ns: AtomicU64::new(0),
            step_read_count: AtomicU64::new(0),
            step_decompress_count: AtomicU64::new(0),
            step_find_boundaries_count: AtomicU64::new(0),
            step_decode_count: AtomicU64::new(0),
            step_group_count: AtomicU64::new(0),
            step_process_count: AtomicU64::new(0),
            step_serialize_count: AtomicU64::new(0),
            step_compress_count: AtomicU64::new(0),
            step_write_count: AtomicU64::new(0),
            read_contention: AtomicU64::new(0),
            boundary_contention: AtomicU64::new(0),
            group_contention: AtomicU64::new(0),
            write_contention: AtomicU64::new(0),
            q1_empty: AtomicU64::new(0),
            q2_empty: AtomicU64::new(0),
            q2b_empty: AtomicU64::new(0),
            q3_empty: AtomicU64::new(0),
            q4_empty: AtomicU64::new(0),
            q5_empty: AtomicU64::new(0),
            q6_empty: AtomicU64::new(0),
            q7_empty: AtomicU64::new(0),
            idle_yields: AtomicU64::new(0),
            // New per-thread stats
            per_thread_step_counts: new_atomic_2d_array::<MAX_THREADS, NUM_STEPS>(),
            per_thread_step_attempts: new_atomic_2d_array::<MAX_THREADS, NUM_STEPS>(),
            per_thread_idle_ns: new_atomic_array::<MAX_THREADS>(),
            per_thread_current_step: new_atomic_u8_array::<MAX_THREADS>(255), // 255 = idle
            boundary_wait_ns: AtomicU64::new(0),
            group_wait_ns: AtomicU64::new(0),
            write_wait_ns: AtomicU64::new(0),
            batch_size_sum: AtomicU64::new(0),
            batch_count: AtomicU64::new(0),
            batch_size_min: AtomicU64::new(u64::MAX),
            batch_size_max: AtomicU64::new(0),
            num_threads: AtomicU64::new(0),
            // Throughput metrics
            bytes_read: AtomicU64::new(0),
            bytes_written: AtomicU64::new(0),
            compressed_bytes_in: AtomicU64::new(0),
            decompressed_bytes: AtomicU64::new(0),
            serialized_bytes: AtomicU64::new(0),
            compressed_bytes_out: AtomicU64::new(0),
            records_decoded: AtomicU64::new(0),
            groups_produced: AtomicU64::new(0),
            // Queue monitoring
            queue_samples: Mutex::new(Vec::new()),
            // Memory limiting stats
            memory_drain_activations: AtomicU64::new(0),
            group_memory_rejects: AtomicU64::new(0),
            peak_memory_bytes: AtomicU64::new(0),

            #[cfg(feature = "memory-debug")]
            memory: MemoryDebugStats::new(),
        }
    }

    /// Set the number of threads for display purposes.
    pub fn set_num_threads(&self, n: usize) {
        self.num_threads.store(n as u64, Ordering::Relaxed);
    }

    /// Record time and success for a step (without per-thread tracking).
    #[inline]
    pub fn record_step(&self, step: PipelineStep, elapsed_ns: u64) {
        self.record_step_for_thread(step, elapsed_ns, None);
    }

    /// Record time and success for a step with per-thread tracking.
    #[inline]
    pub fn record_step_for_thread(
        &self,
        step: PipelineStep,
        elapsed_ns: u64,
        thread_id: Option<usize>,
    ) {
        // Record global stats
        match step {
            PipelineStep::Read => {
                self.step_read_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                self.step_read_count.fetch_add(1, Ordering::Relaxed);
            }
            PipelineStep::Decompress => {
                self.step_decompress_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                self.step_decompress_count.fetch_add(1, Ordering::Relaxed);
            }
            PipelineStep::FindBoundaries => {
                self.step_find_boundaries_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                self.step_find_boundaries_count.fetch_add(1, Ordering::Relaxed);
            }
            PipelineStep::Decode => {
                self.step_decode_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                self.step_decode_count.fetch_add(1, Ordering::Relaxed);
            }
            PipelineStep::Group => {
                self.step_group_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                self.step_group_count.fetch_add(1, Ordering::Relaxed);
            }
            PipelineStep::Process => {
                self.step_process_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                self.step_process_count.fetch_add(1, Ordering::Relaxed);
            }
            PipelineStep::Serialize => {
                self.step_serialize_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                self.step_serialize_count.fetch_add(1, Ordering::Relaxed);
            }
            PipelineStep::Compress => {
                self.step_compress_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                self.step_compress_count.fetch_add(1, Ordering::Relaxed);
            }
            PipelineStep::Write => {
                self.step_write_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                self.step_write_count.fetch_add(1, Ordering::Relaxed);
            }
        }

        // Record per-thread stats if thread_id provided
        if let Some(tid) = thread_id {
            if tid < MAX_THREADS {
                let step_idx = step as usize;
                self.per_thread_step_counts[tid][step_idx].fetch_add(1, Ordering::Relaxed);
            }
        }
    }

    /// Record lock contention for a mutex.
    #[inline]
    pub fn record_contention(&self, step: PipelineStep) {
        match step {
            PipelineStep::Read => {
                self.read_contention.fetch_add(1, Ordering::Relaxed);
            }
            PipelineStep::FindBoundaries => {
                self.boundary_contention.fetch_add(1, Ordering::Relaxed);
            }
            PipelineStep::Group => {
                self.group_contention.fetch_add(1, Ordering::Relaxed);
            }
            PipelineStep::Write => {
                self.write_contention.fetch_add(1, Ordering::Relaxed);
            }
            _ => {}
        }
    }

    /// Record time spent waiting for a lock (exclusive step contention time).
    #[inline]
    pub fn record_wait_time(&self, step: PipelineStep, wait_ns: u64) {
        match step {
            PipelineStep::FindBoundaries => {
                self.boundary_wait_ns.fetch_add(wait_ns, Ordering::Relaxed);
            }
            PipelineStep::Group => {
                self.group_wait_ns.fetch_add(wait_ns, Ordering::Relaxed);
            }
            PipelineStep::Write => {
                self.write_wait_ns.fetch_add(wait_ns, Ordering::Relaxed);
            }
            _ => {}
        }
    }

    /// Record empty queue poll.
    ///
    /// Queue numbers: 1=raw, 2=decompressed, 25=boundaries (Q2b), 3=decoded,
    /// 4=groups, 5=processed, 6=serialized, 7=compressed
    #[inline]
    pub fn record_queue_empty(&self, queue_num: usize) {
        match queue_num {
            1 => self.q1_empty.fetch_add(1, Ordering::Relaxed),
            2 => self.q2_empty.fetch_add(1, Ordering::Relaxed),
            25 => self.q2b_empty.fetch_add(1, Ordering::Relaxed), // Q2b (boundaries)
            3 => self.q3_empty.fetch_add(1, Ordering::Relaxed),
            4 => self.q4_empty.fetch_add(1, Ordering::Relaxed),
            5 => self.q5_empty.fetch_add(1, Ordering::Relaxed),
            6 => self.q6_empty.fetch_add(1, Ordering::Relaxed),
            7 => self.q7_empty.fetch_add(1, Ordering::Relaxed),
            _ => 0,
        };
    }

    /// Record an idle yield (without timing).
    #[inline]
    pub fn record_idle(&self) {
        self.idle_yields.fetch_add(1, Ordering::Relaxed);
    }

    /// Record idle time for a specific thread.
    #[inline]
    pub fn record_idle_for_thread(&self, thread_id: usize, idle_ns: u64) {
        self.idle_yields.fetch_add(1, Ordering::Relaxed);
        if thread_id < MAX_THREADS {
            self.per_thread_idle_ns[thread_id].fetch_add(idle_ns, Ordering::Relaxed);
        }
    }

    /// Record a step attempt for a specific thread (called before attempting).
    /// This tracks total attempts regardless of success/failure.
    #[inline]
    #[allow(clippy::cast_possible_truncation)]
    pub fn record_step_attempt(&self, thread_id: usize, step: PipelineStep) {
        if thread_id < MAX_THREADS {
            let step_idx = step as usize;
            self.per_thread_step_attempts[thread_id][step_idx].fetch_add(1, Ordering::Relaxed);
            // Also update current step
            self.per_thread_current_step[thread_id].store(step_idx as u8, Ordering::Relaxed);
        }
    }

    /// Set the current step a thread is working on.
    #[inline]
    pub fn set_current_step(&self, thread_id: usize, step: PipelineStep) {
        if thread_id < MAX_THREADS {
            self.per_thread_current_step[thread_id].store(step as u8, Ordering::Relaxed);
        }
    }

    /// Clear the current step (thread is idle or between steps).
    #[inline]
    pub fn clear_current_step(&self, thread_id: usize) {
        if thread_id < MAX_THREADS {
            self.per_thread_current_step[thread_id].store(255, Ordering::Relaxed);
        }
    }

    /// Get current step for all threads (for activity snapshot).
    #[allow(clippy::cast_possible_truncation)]
    pub fn get_thread_activity(&self, num_threads: usize) -> Vec<Option<PipelineStep>> {
        (0..num_threads.min(MAX_THREADS))
            .map(|tid| {
                let step_idx = self.per_thread_current_step[tid].load(Ordering::Relaxed);
                if step_idx < NUM_STEPS as u8 {
                    Some(PipelineStep::from_index(step_idx as usize))
                } else {
                    None // idle
                }
            })
            .collect()
    }

    /// Record a batch size for batch size distribution tracking.
    #[inline]
    pub fn record_batch_size(&self, size: usize) {
        let size = size as u64;
        self.batch_size_sum.fetch_add(size, Ordering::Relaxed);
        self.batch_count.fetch_add(1, Ordering::Relaxed);

        // Update min (using compare-exchange loop)
        let mut current_min = self.batch_size_min.load(Ordering::Relaxed);
        while size < current_min {
            match self.batch_size_min.compare_exchange_weak(
                current_min,
                size,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(actual) => current_min = actual,
            }
        }

        // Update max
        let mut current_max = self.batch_size_max.load(Ordering::Relaxed);
        while size > current_max {
            match self.batch_size_max.compare_exchange_weak(
                current_max,
                size,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(actual) => current_max = actual,
            }
        }
    }

    /// Add a queue size sample for periodic monitoring.
    /// Called from a background monitor thread every ~100ms.
    pub fn add_queue_sample(&self, sample: QueueSample) {
        self.queue_samples.lock().push(sample);
    }

    /// Get the collected queue samples for analysis.
    pub fn get_queue_samples(&self) -> Vec<QueueSample> {
        self.queue_samples.lock().clone()
    }

    /// Record a memory drain mode activation.
    #[inline]
    pub fn record_memory_drain_activation(&self) {
        self.memory_drain_activations.fetch_add(1, Ordering::Relaxed);
    }

    /// Record a Group step rejection due to memory limit.
    #[inline]
    pub fn record_group_memory_reject(&self) {
        self.group_memory_rejects.fetch_add(1, Ordering::Relaxed);
    }

    /// Update peak memory usage if current is higher.
    #[inline]
    pub fn record_memory_usage(&self, bytes: u64) {
        let mut current_peak = self.peak_memory_bytes.load(Ordering::Relaxed);
        while bytes > current_peak {
            match self.peak_memory_bytes.compare_exchange_weak(
                current_peak,
                bytes,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(actual) => current_peak = actual,
            }
        }
    }

    /// Format statistics as a human-readable summary.
    #[allow(clippy::similar_names)] // Intentional: xxx_ns vs xxx_ms for nanoseconds vs milliseconds
    #[allow(
        clippy::too_many_lines,
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    pub fn format_summary(&self) -> String {
        use std::fmt::Write;

        let mut s = String::new();
        writeln!(s, "Pipeline Statistics:").unwrap();
        writeln!(s).unwrap();

        // Helper to format step stats
        #[allow(clippy::uninlined_format_args)]
        let format_step = |name: &str, ns: u64, count: u64| -> String {
            if count == 0 {
                format!("  {:<20} {:>10} ops, {:>12}", name, 0, "-")
            } else {
                let total_ms = ns as f64 / 1_000_000.0;
                let avg_us = (ns as f64 / count as f64) / 1_000.0;
                format!(
                    "  {:<20} {:>10} ops, {:>10.1}ms total, {:>8.1}µs avg",
                    name, count, total_ms, avg_us
                )
            }
        };

        writeln!(s, "Step Timing:").unwrap();
        writeln!(
            s,
            "{}",
            format_step(
                "Read",
                self.step_read_ns.load(Ordering::Relaxed),
                self.step_read_count.load(Ordering::Relaxed)
            )
        )
        .unwrap();
        writeln!(
            s,
            "{}",
            format_step(
                "Decompress",
                self.step_decompress_ns.load(Ordering::Relaxed),
                self.step_decompress_count.load(Ordering::Relaxed)
            )
        )
        .unwrap();
        writeln!(
            s,
            "{}",
            format_step(
                "FindBoundaries",
                self.step_find_boundaries_ns.load(Ordering::Relaxed),
                self.step_find_boundaries_count.load(Ordering::Relaxed)
            )
        )
        .unwrap();
        writeln!(
            s,
            "{}",
            format_step(
                "Decode",
                self.step_decode_ns.load(Ordering::Relaxed),
                self.step_decode_count.load(Ordering::Relaxed)
            )
        )
        .unwrap();
        writeln!(
            s,
            "{}",
            format_step(
                "Group",
                self.step_group_ns.load(Ordering::Relaxed),
                self.step_group_count.load(Ordering::Relaxed)
            )
        )
        .unwrap();
        writeln!(
            s,
            "{}",
            format_step(
                "Process",
                self.step_process_ns.load(Ordering::Relaxed),
                self.step_process_count.load(Ordering::Relaxed)
            )
        )
        .unwrap();
        writeln!(
            s,
            "{}",
            format_step(
                "Serialize",
                self.step_serialize_ns.load(Ordering::Relaxed),
                self.step_serialize_count.load(Ordering::Relaxed)
            )
        )
        .unwrap();
        writeln!(
            s,
            "{}",
            format_step(
                "Compress",
                self.step_compress_ns.load(Ordering::Relaxed),
                self.step_compress_count.load(Ordering::Relaxed)
            )
        )
        .unwrap();
        writeln!(
            s,
            "{}",
            format_step(
                "Write",
                self.step_write_ns.load(Ordering::Relaxed),
                self.step_write_count.load(Ordering::Relaxed)
            )
        )
        .unwrap();

        writeln!(s).unwrap();
        writeln!(s, "Contention:").unwrap();
        writeln!(
            s,
            "  Read lock:     {:>10} failed attempts",
            self.read_contention.load(Ordering::Relaxed)
        )
        .unwrap();
        writeln!(
            s,
            "  Boundary lock: {:>10} failed attempts",
            self.boundary_contention.load(Ordering::Relaxed)
        )
        .unwrap();
        writeln!(
            s,
            "  Group lock:    {:>10} failed attempts",
            self.group_contention.load(Ordering::Relaxed)
        )
        .unwrap();
        writeln!(
            s,
            "  Write lock:    {:>10} failed attempts",
            self.write_contention.load(Ordering::Relaxed)
        )
        .unwrap();

        writeln!(s).unwrap();
        writeln!(s, "Queue Empty Polls:").unwrap();
        writeln!(s, "  Q1 (raw):        {:>10}", self.q1_empty.load(Ordering::Relaxed)).unwrap();
        writeln!(s, "  Q2 (decomp):     {:>10}", self.q2_empty.load(Ordering::Relaxed)).unwrap();
        writeln!(s, "  Q2b (boundary):  {:>10}", self.q2b_empty.load(Ordering::Relaxed)).unwrap();
        writeln!(s, "  Q3 (decoded):    {:>10}", self.q3_empty.load(Ordering::Relaxed)).unwrap();
        writeln!(s, "  Q4 (groups):     {:>10}", self.q4_empty.load(Ordering::Relaxed)).unwrap();
        writeln!(s, "  Q5 (processed):  {:>10}", self.q5_empty.load(Ordering::Relaxed)).unwrap();
        writeln!(s, "  Q6 (serialized): {:>10}", self.q6_empty.load(Ordering::Relaxed)).unwrap();
        writeln!(s, "  Q7 (compressed): {:>10}", self.q7_empty.load(Ordering::Relaxed)).unwrap();

        writeln!(s).unwrap();
        writeln!(s, "Idle Yields: {:>10}", self.idle_yields.load(Ordering::Relaxed)).unwrap();

        // NEW: Wait time for exclusive steps
        let boundary_wait = self.boundary_wait_ns.load(Ordering::Relaxed);
        let group_wait = self.group_wait_ns.load(Ordering::Relaxed);
        let write_wait = self.write_wait_ns.load(Ordering::Relaxed);
        if boundary_wait > 0 || group_wait > 0 || write_wait > 0 {
            writeln!(s).unwrap();
            writeln!(s, "Lock Wait Time:").unwrap();
            writeln!(s, "  Boundary lock: {:>10.1}ms", boundary_wait as f64 / 1_000_000.0).unwrap();
            writeln!(s, "  Group lock:    {:>10.1}ms", group_wait as f64 / 1_000_000.0).unwrap();
            writeln!(s, "  Write lock:    {:>10.1}ms", write_wait as f64 / 1_000_000.0).unwrap();
        }

        // NEW: Batch size statistics
        let batch_count = self.batch_count.load(Ordering::Relaxed);
        if batch_count > 0 {
            let batch_sum = self.batch_size_sum.load(Ordering::Relaxed);
            let batch_min = self.batch_size_min.load(Ordering::Relaxed);
            let batch_max = self.batch_size_max.load(Ordering::Relaxed);
            let batch_avg = batch_sum as f64 / batch_count as f64;

            writeln!(s).unwrap();
            writeln!(s, "Batch Size (records per batch):").unwrap();
            writeln!(s, "  Count:   {batch_count:>10}").unwrap();
            writeln!(s, "  Min:     {batch_min:>10}").unwrap();
            writeln!(s, "  Max:     {batch_max:>10}").unwrap();
            writeln!(s, "  Average: {batch_avg:>10.1}").unwrap();
        }

        // NEW: Per-thread work distribution
        let num_threads = self.num_threads.load(Ordering::Relaxed) as usize;
        if num_threads > 0 {
            writeln!(s).unwrap();
            writeln!(s, "Per-Thread Work Distribution:").unwrap();

            // Step names for the header
            let step_names = ["Rd", "Dc", "Fb", "De", "Gr", "Pr", "Se", "Co", "Wr"];

            // Header
            write!(s, "  Thread ").unwrap();
            for name in &step_names {
                write!(s, " {name:>6}").unwrap();
            }
            writeln!(s, "    Idle ms").unwrap();

            // Per-thread rows
            for tid in 0..num_threads.min(MAX_THREADS) {
                write!(s, "  T{tid:<5} ").unwrap();
                for step_idx in 0..NUM_STEPS {
                    let count = self.per_thread_step_counts[tid][step_idx].load(Ordering::Relaxed);
                    write!(s, " {count:>6}").unwrap();
                }
                let idle_ns = self.per_thread_idle_ns[tid].load(Ordering::Relaxed);
                writeln!(s, " {:>10.1}", idle_ns as f64 / 1_000_000.0).unwrap();
            }

            // Total row
            write!(s, "  Total  ").unwrap();
            for step_idx in 0..NUM_STEPS {
                let mut total = 0u64;
                for tid in 0..num_threads.min(MAX_THREADS) {
                    total += self.per_thread_step_counts[tid][step_idx].load(Ordering::Relaxed);
                }
                write!(s, " {total:>6}").unwrap();
            }
            let total_idle: u64 = (0..num_threads.min(MAX_THREADS))
                .map(|tid| self.per_thread_idle_ns[tid].load(Ordering::Relaxed))
                .sum();
            writeln!(s, " {:>10.1}", total_idle as f64 / 1_000_000.0).unwrap();

            // Per-thread attempt statistics (shows success rate)
            writeln!(s).unwrap();
            writeln!(s, "Per-Thread Attempt Success Rate:").unwrap();

            // Header
            write!(s, "  Thread ").unwrap();
            for name in &step_names {
                write!(s, " {name:>6}").unwrap();
            }
            writeln!(s, "   Total%").unwrap();

            // Per-thread rows with success rates
            for tid in 0..num_threads.min(MAX_THREADS) {
                write!(s, "  T{tid:<5} ").unwrap();
                let mut thread_attempts = 0u64;
                let mut thread_successes = 0u64;
                for step_idx in 0..NUM_STEPS {
                    let attempts =
                        self.per_thread_step_attempts[tid][step_idx].load(Ordering::Relaxed);
                    let successes =
                        self.per_thread_step_counts[tid][step_idx].load(Ordering::Relaxed);
                    thread_attempts += attempts;
                    thread_successes += successes;
                    if attempts == 0 {
                        write!(s, "   -  ").unwrap();
                    } else {
                        let rate = (successes as f64 / attempts as f64) * 100.0;
                        write!(s, " {rate:>5.0}%").unwrap();
                    }
                }
                if thread_attempts == 0 {
                    writeln!(s, "      -").unwrap();
                } else {
                    let total_rate = (thread_successes as f64 / thread_attempts as f64) * 100.0;
                    writeln!(s, "  {total_rate:>5.1}%").unwrap();
                }
            }
        }

        // Thread utilization summary
        let total_work_ns = self.step_read_ns.load(Ordering::Relaxed)
            + self.step_decompress_ns.load(Ordering::Relaxed)
            + self.step_find_boundaries_ns.load(Ordering::Relaxed)
            + self.step_decode_ns.load(Ordering::Relaxed)
            + self.step_group_ns.load(Ordering::Relaxed)
            + self.step_process_ns.load(Ordering::Relaxed)
            + self.step_serialize_ns.load(Ordering::Relaxed)
            + self.step_compress_ns.load(Ordering::Relaxed)
            + self.step_write_ns.load(Ordering::Relaxed);

        let total_idle_ns: u64 = (0..num_threads.min(MAX_THREADS))
            .map(|tid| self.per_thread_idle_ns[tid].load(Ordering::Relaxed))
            .sum();

        let total_contention = self.read_contention.load(Ordering::Relaxed)
            + self.boundary_contention.load(Ordering::Relaxed)
            + self.group_contention.load(Ordering::Relaxed)
            + self.write_contention.load(Ordering::Relaxed);

        if total_work_ns > 0 {
            writeln!(s).unwrap();
            writeln!(s, "Thread Utilization:").unwrap();

            let work_ms = total_work_ns as f64 / 1_000_000.0;
            let idle_ms = total_idle_ns as f64 / 1_000_000.0;
            let total_thread_ms = work_ms + idle_ms;

            if total_thread_ms > 0.0 {
                let utilization = (work_ms / total_thread_ms) * 100.0;
                writeln!(s, "  Work time:       {work_ms:>10.1}ms").unwrap();
                writeln!(s, "  Idle time:       {idle_ms:>10.1}ms").unwrap();
                writeln!(s, "  Utilization:     {utilization:>10.1}%").unwrap();
                writeln!(s, "  Contention attempts: {total_contention:>7}").unwrap();
            }
        }

        // ========== Throughput Metrics ==========
        let bytes_read = self.bytes_read.load(Ordering::Relaxed);
        let bytes_written = self.bytes_written.load(Ordering::Relaxed);
        let compressed_bytes_in = self.compressed_bytes_in.load(Ordering::Relaxed);
        let decompressed_bytes = self.decompressed_bytes.load(Ordering::Relaxed);
        let serialized_bytes = self.serialized_bytes.load(Ordering::Relaxed);
        let compressed_bytes_out = self.compressed_bytes_out.load(Ordering::Relaxed);
        let records_decoded = self.records_decoded.load(Ordering::Relaxed);
        let groups_produced = self.groups_produced.load(Ordering::Relaxed);

        // Only show throughput section if we have data
        if bytes_read > 0 || bytes_written > 0 {
            writeln!(s).unwrap();
            writeln!(s, "Throughput:").unwrap();

            // Helper function to format bytes nicely
            let format_bytes = |bytes: u64| -> String {
                if bytes >= 1_000_000_000 {
                    format!("{:.2} GB", bytes as f64 / 1_000_000_000.0)
                } else if bytes >= 1_000_000 {
                    format!("{:.1} MB", bytes as f64 / 1_000_000.0)
                } else if bytes >= 1_000 {
                    format!("{:.1} KB", bytes as f64 / 1_000.0)
                } else {
                    format!("{bytes} B")
                }
            };

            // Helper to format count with K/M suffix
            let format_count = |count: u64| -> String {
                if count >= 1_000_000 {
                    format!("{:.2}M", count as f64 / 1_000_000.0)
                } else if count >= 1_000 {
                    format!("{:.1}K", count as f64 / 1_000.0)
                } else {
                    format!("{count}")
                }
            };

            // Read throughput
            let read_ns = self.step_read_ns.load(Ordering::Relaxed);
            if bytes_read > 0 && read_ns > 0 {
                let read_ms = read_ns as f64 / 1_000_000.0;
                let read_mb_s = (bytes_read as f64 / 1_000_000.0) / (read_ms / 1000.0);
                writeln!(
                    s,
                    "  Read:        {:>10} in {:>8.1}ms = {:>8.1} MB/s",
                    format_bytes(bytes_read),
                    read_ms,
                    read_mb_s
                )
                .unwrap();
            }

            // Decompress throughput (input → output with expansion ratio)
            let decompress_ns = self.step_decompress_ns.load(Ordering::Relaxed);
            if compressed_bytes_in > 0 && decompressed_bytes > 0 && decompress_ns > 0 {
                let decompress_ms = decompress_ns as f64 / 1_000_000.0;
                let in_mb_s = (compressed_bytes_in as f64 / 1_000_000.0) / (decompress_ms / 1000.0);
                let out_mb_s = (decompressed_bytes as f64 / 1_000_000.0) / (decompress_ms / 1000.0);
                let expansion = decompressed_bytes as f64 / compressed_bytes_in as f64;
                writeln!(
                    s,
                    "  Decompress:  {:>10} → {:>10} = {:>6.1} → {:>6.1} MB/s ({:.2}x expansion)",
                    format_bytes(compressed_bytes_in),
                    format_bytes(decompressed_bytes),
                    in_mb_s,
                    out_mb_s,
                    expansion
                )
                .unwrap();
            }

            // Decode throughput (records per second)
            let decode_ns = self.step_decode_ns.load(Ordering::Relaxed);
            if records_decoded > 0 && decode_ns > 0 {
                let decode_ms = decode_ns as f64 / 1_000_000.0;
                let records_per_s = records_decoded as f64 / (decode_ms / 1000.0);
                writeln!(
                    s,
                    "  Decode:      {:>10} records     = {:>8} records/s",
                    format_count(records_decoded),
                    format_count(records_per_s as u64)
                )
                .unwrap();
            }

            // Group throughput (records in → groups out)
            let group_ns = self.step_group_ns.load(Ordering::Relaxed);
            if records_decoded > 0 && groups_produced > 0 && group_ns > 0 {
                let group_ms = group_ns as f64 / 1_000_000.0;
                let records_in_per_s = records_decoded as f64 / (group_ms / 1000.0);
                let groups_out_per_s = groups_produced as f64 / (group_ms / 1000.0);
                writeln!(
                    s,
                    "  Group:       {:>10} → {:>10} = {:>6} records/s in, {:>6} groups/s out",
                    format_count(records_decoded),
                    format_count(groups_produced),
                    format_count(records_in_per_s as u64),
                    format_count(groups_out_per_s as u64)
                )
                .unwrap();
            }

            // Process throughput (groups per second)
            let process_ns = self.step_process_ns.load(Ordering::Relaxed);
            if groups_produced > 0 && process_ns > 0 {
                let process_ms = process_ns as f64 / 1_000_000.0;
                let groups_per_s = groups_produced as f64 / (process_ms / 1000.0);
                writeln!(
                    s,
                    "  Process:     {:>10} groups      = {:>8} groups/s",
                    format_count(groups_produced),
                    format_count(groups_per_s as u64)
                )
                .unwrap();
            }

            // Serialize throughput
            let serialize_ns = self.step_serialize_ns.load(Ordering::Relaxed);
            if serialized_bytes > 0 && serialize_ns > 0 {
                let serialize_ms = serialize_ns as f64 / 1_000_000.0;
                let mb_per_s = (serialized_bytes as f64 / 1_000_000.0) / (serialize_ms / 1000.0);
                writeln!(
                    s,
                    "  Serialize:   {:>10}             = {:>8.1} MB/s",
                    format_bytes(serialized_bytes),
                    mb_per_s
                )
                .unwrap();
            }

            // Compress throughput (input → output with compression ratio)
            let compress_ns = self.step_compress_ns.load(Ordering::Relaxed);
            if serialized_bytes > 0 && compressed_bytes_out > 0 && compress_ns > 0 {
                let compress_ms = compress_ns as f64 / 1_000_000.0;
                let in_mb_s = (serialized_bytes as f64 / 1_000_000.0) / (compress_ms / 1000.0);
                let out_mb_s = (compressed_bytes_out as f64 / 1_000_000.0) / (compress_ms / 1000.0);
                let compression = serialized_bytes as f64 / compressed_bytes_out as f64;
                writeln!(
                    s,
                    "  Compress:    {:>10} → {:>10} = {:>6.1} → {:>6.1} MB/s ({:.2}x compression)",
                    format_bytes(serialized_bytes),
                    format_bytes(compressed_bytes_out),
                    in_mb_s,
                    out_mb_s,
                    compression
                )
                .unwrap();
            }

            // Write throughput
            let write_ns = self.step_write_ns.load(Ordering::Relaxed);
            if bytes_written > 0 && write_ns > 0 {
                let write_ms = write_ns as f64 / 1_000_000.0;
                let write_mb_s = (bytes_written as f64 / 1_000_000.0) / (write_ms / 1000.0);
                writeln!(
                    s,
                    "  Write:       {:>10} in {:>8.1}ms = {:>8.1} MB/s",
                    format_bytes(bytes_written),
                    write_ms,
                    write_mb_s
                )
                .unwrap();
            }
        }

        // Queue sample summary if we have any
        let samples = self.queue_samples.lock();
        if !samples.is_empty() {
            writeln!(s).unwrap();
            writeln!(s, "Queue Size Timeline ({} samples at ~100ms intervals):", samples.len())
                .unwrap();
            writeln!(
                s,
                "  Time   Q1   Q2  Q2b   Q3   Q4   Q5   Q6   Q7 | R2  R3  R7 |  R3_MB  Threads"
            )
            .unwrap();

            // Show all samples
            for sample in samples.iter() {
                let r3_mb = sample.reorder_memory_bytes[1] as f64 / 1_048_576.0;
                write!(
                    s,
                    "  {:>4}  {:>3}  {:>3}  {:>3}  {:>3}  {:>3}  {:>3}  {:>3}  {:>3} | {:>3} {:>3} {:>3} | {:>6.1}  ",
                    sample.time_ms,
                    sample.queue_sizes[0],
                    sample.queue_sizes[1],
                    sample.queue_sizes[2],
                    sample.queue_sizes[3],
                    sample.queue_sizes[4],
                    sample.queue_sizes[5],
                    sample.queue_sizes[6],
                    sample.queue_sizes[7],
                    sample.reorder_sizes[0],
                    sample.reorder_sizes[1],
                    sample.reorder_sizes[2],
                    r3_mb,
                )
                .unwrap();
                // Show thread activity as compact string
                for &step_idx in &sample.thread_steps {
                    if step_idx < NUM_STEPS as u8 {
                        let short = match step_idx {
                            0 => "R",
                            1 => "D",
                            2 => "F",
                            3 => "d",
                            4 => "G",
                            5 => "P",
                            6 => "S",
                            7 => "C",
                            8 => "W",
                            _ => "?",
                        };
                        write!(s, "{short}").unwrap();
                    } else {
                        write!(s, ".").unwrap();
                    }
                }
                writeln!(s).unwrap();
            }

            // Summary of peak reorder buffer usage
            let peak_r3_items = samples.iter().map(|s| s.reorder_sizes[1]).max().unwrap_or(0);
            let peak_r3_bytes =
                samples.iter().map(|s| s.reorder_memory_bytes[1]).max().unwrap_or(0);
            let peak_r3_mb = peak_r3_bytes as f64 / 1_048_576.0;
            writeln!(s).unwrap();
            writeln!(s, "Peak Q3 Reorder Buffer: {peak_r3_items} items, {peak_r3_mb:.1} MB")
                .unwrap();
        }

        // Memory limiting statistics
        let group_rejects = self.group_memory_rejects.load(Ordering::Relaxed);
        let peak_memory = self.peak_memory_bytes.load(Ordering::Relaxed);

        if group_rejects > 0 || peak_memory > 0 {
            writeln!(s).unwrap();
            writeln!(s, "Memory Limiting:").unwrap();
            if group_rejects > 0 {
                writeln!(s, "  Group rejects (memory): {group_rejects:>10}").unwrap();
            }
            if peak_memory > 0 {
                let peak_mb = peak_memory as f64 / 1_048_576.0;
                writeln!(s, "  Peak memory usage:      {peak_mb:>10.1} MB").unwrap();
            }
        }

        s
    }

    /// Log statistics using the log crate at info level.
    pub fn log_summary(&self) {
        for line in self.format_summary().lines() {
            info!("{line}");
        }
    }
}

// ============================================================================
// Comprehensive Memory Tracking Methods (behind memory-debug feature)
// ============================================================================

#[cfg(feature = "memory-debug")]
impl PipelineStats {
    /// Track queue memory addition (queue memory is separate from processing memory)
    pub fn track_queue_memory_add(&self, queue_name: &str, size: usize) {
        let m = &self.memory;
        let size_u64 = size as u64;
        match queue_name {
            "q1" => m.q1_memory_bytes.fetch_add(size_u64, Ordering::Relaxed),
            "q2" => m.q2_memory_bytes.fetch_add(size_u64, Ordering::Relaxed),
            "q3" => m.q3_memory_bytes.fetch_add(size_u64, Ordering::Relaxed),
            "q4" => m.q4_memory_bytes.fetch_add(size_u64, Ordering::Relaxed),
            "q5" => m.q5_memory_bytes.fetch_add(size_u64, Ordering::Relaxed),
            "q6" => m.q6_memory_bytes.fetch_add(size_u64, Ordering::Relaxed),
            "q7" => m.q7_memory_bytes.fetch_add(size_u64, Ordering::Relaxed),
            _ => 0,
        };
    }

    /// Track queue memory removal (queue memory is separate from processing memory).
    /// Uses saturating subtraction to prevent u64 underflow from estimation mismatches.
    pub fn track_queue_memory_remove(&self, queue_name: &str, size: usize) {
        let m = &self.memory;
        let counter = match queue_name {
            "q1" => &m.q1_memory_bytes,
            "q2" => &m.q2_memory_bytes,
            "q3" => &m.q3_memory_bytes,
            "q4" => &m.q4_memory_bytes,
            "q5" => &m.q5_memory_bytes,
            "q6" => &m.q6_memory_bytes,
            "q7" => &m.q7_memory_bytes,
            _ => return,
        };
        let size_u64 = size as u64;
        let mut current = counter.load(Ordering::Relaxed);
        loop {
            let new_val = current.saturating_sub(size_u64);
            match counter.compare_exchange_weak(
                current,
                new_val,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(actual) => current = actual,
            }
        }
    }

    /// Track position group processing memory (separate from queue memory).
    /// Uses saturating subtraction to prevent u64 underflow from estimation mismatches.
    pub fn track_position_group_memory(&self, size: usize, is_allocation: bool) {
        let counter = &self.memory.position_group_processing_bytes;
        let size_u64 = size as u64;
        if is_allocation {
            counter.fetch_add(size_u64, Ordering::Relaxed);
        } else {
            let mut current = counter.load(Ordering::Relaxed);
            loop {
                let new_val = current.saturating_sub(size_u64);
                match counter.compare_exchange_weak(
                    current,
                    new_val,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => break,
                    Err(actual) => current = actual,
                }
            }
        }
    }

    /// Track template processing memory (separate from queue memory).
    /// Uses saturating subtraction to prevent u64 underflow from estimation mismatches.
    pub fn track_template_memory(&self, size: usize, is_allocation: bool) {
        let counter = &self.memory.template_processing_bytes;
        let size_u64 = size as u64;
        if is_allocation {
            counter.fetch_add(size_u64, Ordering::Relaxed);
        } else {
            let mut current = counter.load(Ordering::Relaxed);
            loop {
                let new_val = current.saturating_sub(size_u64);
                match counter.compare_exchange_weak(
                    current,
                    new_val,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => break,
                    Err(actual) => current = actual,
                }
            }
        }
    }

    /// Update system RSS.
    pub fn update_system_rss(&self, rss_bytes: u64) {
        self.memory.system_rss_bytes.store(rss_bytes, Ordering::Relaxed);
    }

    /// Set infrastructure memory estimates (call once at pipeline startup).
    pub fn set_infrastructure_memory(&self, num_threads: usize, queue_capacity: usize) {
        let m = &self.memory;
        m.decompressor_memory_bytes.store(num_threads as u64 * 32 * 1024, Ordering::Relaxed);
        m.compressor_memory_bytes.store(num_threads as u64 * 280 * 1024, Ordering::Relaxed);
        m.worker_buffer_memory_bytes.store(num_threads as u64 * 512 * 1024, Ordering::Relaxed);
        m.io_buffer_memory_bytes.store(16u64 * 1024 * 1024, Ordering::Relaxed);
        m.thread_stack_memory_bytes
            .store((num_threads as u64 + 1) * 2 * 1024 * 1024, Ordering::Relaxed);
        m.queue_capacity_memory_bytes.store(7u64 * queue_capacity as u64 * 128, Ordering::Relaxed);
    }

    /// Update queue memory stats from actual queues (call periodically during monitoring)
    pub fn update_queue_memory_from_external(&self, queue_stats: &[(&str, u64)]) {
        let m = &self.memory;
        for (queue_name, current_bytes) in queue_stats {
            match *queue_name {
                "q1" => m.q1_memory_bytes.store(*current_bytes, Ordering::Relaxed),
                "q2" => m.q2_memory_bytes.store(*current_bytes, Ordering::Relaxed),
                "q3" => m.q3_memory_bytes.store(*current_bytes, Ordering::Relaxed),
                "q4" => m.q4_memory_bytes.store(*current_bytes, Ordering::Relaxed),
                "q5" => m.q5_memory_bytes.store(*current_bytes, Ordering::Relaxed),
                "q6" => m.q6_memory_bytes.store(*current_bytes, Ordering::Relaxed),
                "q7" => m.q7_memory_bytes.store(*current_bytes, Ordering::Relaxed),
                _ => {}
            }
        }
    }

    /// Get current memory breakdown
    pub fn get_memory_breakdown(&self) -> MemoryBreakdown {
        let m = &self.memory;

        // Load each counter once to avoid divergence under contention
        let q1 = m.q1_memory_bytes.load(Ordering::Relaxed);
        let q2 = m.q2_memory_bytes.load(Ordering::Relaxed);
        let q3 = m.q3_memory_bytes.load(Ordering::Relaxed);
        let q4 = m.q4_memory_bytes.load(Ordering::Relaxed);
        let q5 = m.q5_memory_bytes.load(Ordering::Relaxed);
        let q6 = m.q6_memory_bytes.load(Ordering::Relaxed);
        let q7 = m.q7_memory_bytes.load(Ordering::Relaxed);
        let queue_total = q1 + q2 + q3 + q4 + q5 + q6 + q7;

        let pos_groups = m.position_group_processing_bytes.load(Ordering::Relaxed);
        let templates = m.template_processing_bytes.load(Ordering::Relaxed);
        let reorder = m.reorder_buffer_bytes.load(Ordering::Relaxed);
        let grouper = m.grouper_memory_bytes.load(Ordering::Relaxed);
        let worker_local = m.worker_local_memory_bytes.load(Ordering::Relaxed);
        let processing_total = pos_groups + templates + reorder + grouper + worker_local;

        let infra_decompressors = m.decompressor_memory_bytes.load(Ordering::Relaxed);
        let infra_compressors = m.compressor_memory_bytes.load(Ordering::Relaxed);
        let infra_buffers = m.worker_buffer_memory_bytes.load(Ordering::Relaxed);
        let infra_io = m.io_buffer_memory_bytes.load(Ordering::Relaxed);
        let infra_stacks = m.thread_stack_memory_bytes.load(Ordering::Relaxed);
        let infra_queues = m.queue_capacity_memory_bytes.load(Ordering::Relaxed);
        let infra_total = infra_decompressors
            + infra_compressors
            + infra_buffers
            + infra_io
            + infra_stacks
            + infra_queues;

        let tracked_total = queue_total + processing_total + infra_total;
        let system_rss = m.system_rss_bytes.load(Ordering::Relaxed);
        let untracked = system_rss.saturating_sub(tracked_total);

        MemoryBreakdown {
            system_rss_gb: system_rss as f64 / 1e9,
            tracked_total_gb: tracked_total as f64 / 1e9,
            untracked_gb: untracked as f64 / 1e9,

            q1_mb: q1 as f64 / 1e6,
            q2_mb: q2 as f64 / 1e6,
            q3_mb: q3 as f64 / 1e6,
            q4_gb: q4 as f64 / 1e9,
            q5_gb: q5 as f64 / 1e9,
            q6_mb: q6 as f64 / 1e6,
            q7_mb: q7 as f64 / 1e6,

            position_groups_gb: pos_groups as f64 / 1e9,
            templates_gb: templates as f64 / 1e9,
            reorder_buffers_mb: reorder as f64 / 1e6,
            grouper_mb: grouper as f64 / 1e6,
            worker_local_mb: worker_local as f64 / 1e6,

            decompressors_mb: infra_decompressors as f64 / 1e6,
            compressors_mb: infra_compressors as f64 / 1e6,
            worker_buffers_mb: infra_buffers as f64 / 1e6,
            io_buffers_mb: infra_io as f64 / 1e6,
            thread_stacks_mb: infra_stacks as f64 / 1e6,
            queue_capacity_mb: infra_queues as f64 / 1e6,
            infrastructure_gb: infra_total as f64 / 1e9,
        }
    }
}

// ============================================================================
// Shared Traits for Steps 4-7 (used by both BAM and FASTQ pipelines)
// ============================================================================

/// Trait for accessing shared pipeline state needed by steps 5-7.
///
/// This allows the serialize, compress, and write steps to be generic over
/// both BAM and FASTQ pipelines, reducing code duplication.
pub trait OutputPipelineState: Send + Sync {
    /// The type of processed items (output of process step).
    type Processed: Send;

    // Error handling
    fn has_error(&self) -> bool;
    fn set_error(&self, error: io::Error);

    // Queue 5 access (Serialize → Compress)
    fn q5_pop(&self) -> Option<(u64, SerializedBatch)>;
    /// # Errors
    ///
    /// Returns the item if the queue is full.
    fn q5_push(&self, item: (u64, SerializedBatch)) -> Result<(), (u64, SerializedBatch)>;
    fn q5_is_full(&self) -> bool;
    /// Track memory released when popping from Q5.
    fn q5_track_pop(&self, _heap_size: u64) {}

    // Queue 6 access (Compress → Write)
    fn q6_pop(&self) -> Option<(u64, CompressedBlockBatch)>;
    /// # Errors
    ///
    /// Returns the item if the queue is full.
    fn q6_push(&self, item: (u64, CompressedBlockBatch))
    -> Result<(), (u64, CompressedBlockBatch)>;
    fn q6_is_full(&self) -> bool;
    /// Track memory released when popping from Q6.
    fn q6_track_pop(&self, _heap_size: u64) {}

    // Reorder buffer for Q6
    fn q6_reorder_insert(&self, serial: u64, batch: CompressedBlockBatch);
    fn q6_reorder_try_pop_next(&self) -> Option<CompressedBlockBatch>;

    // Output access
    fn output_try_lock(&self)
    -> Option<parking_lot::MutexGuard<'_, Option<Box<dyn Write + Send>>>>;

    // Completion tracking
    fn increment_written(&self) -> u64;

    // Throughput metrics (optional, with default no-op implementation)
    /// Record compressed bytes output from compress step.
    fn record_compressed_bytes_out(&self, _bytes: u64) {}

    // Deadlock detection progress tracking (optional, with default no-op)
    /// Record progress when popping from serialized queue (Q6).
    fn record_q6_pop_progress(&self) {}
    /// Record progress when pushing to compressed queue (Q7).
    fn record_q7_push_progress(&self) {}

    // Stats access (optional, with default no-op implementation)
    /// Get reference to pipeline stats for recording metrics.
    fn stats(&self) -> Option<&PipelineStats> {
        None
    }
}

/// Trait for types that have a BGZF compressor.
pub trait HasCompressor {
    fn compressor_mut(&mut self) -> &mut InlineBgzfCompressor;
}

/// Trait for worker states that may hold items between iterations.
///
/// When a worker tries to push to a full queue, it "holds" the item and returns
/// immediately (non-blocking). The next iteration tries to advance the held item
/// before doing new work.
///
/// **CRITICAL**: Worker loops MUST check `has_any_held_items()` before exiting!
/// If a worker exits while holding items, that data is lost. The correct exit
/// condition is:
///
/// ```ignore
/// if state.is_complete() && !worker.has_any_held_items() {
///     break;
/// }
/// ```
///
/// This trait ensures both BAM and FASTQ pipelines use the same exit logic,
/// preventing the race condition where data in held items is lost on exit.
pub trait WorkerStateCommon {
    /// Check if this worker is holding any items that haven't been pushed yet.
    ///
    /// Returns `true` if any `held_*` field contains data. Used to prevent
    /// premature exit from the worker loop - workers must not exit while
    /// holding data or it will be lost.
    fn has_any_held_items(&self) -> bool;

    /// Clear all held items (used during error recovery or cleanup).
    ///
    /// **WARNING**: This discards data! Only use when the pipeline is in an
    /// error state and data loss is acceptable.
    fn clear_held_items(&mut self);
}

/// Trait for workers that have a `WorkerCoreState`.
///
/// This enables shared worker loop functions to access scheduler and backoff state.
pub trait HasWorkerCore {
    fn core(&self) -> &WorkerCoreState;
    fn core_mut(&mut self) -> &mut WorkerCoreState;
}

/// Handle worker backoff with optional idle time tracking.
///
/// This consolidates the common backoff pattern used in both BAM and FASTQ worker loops:
/// - If work was done, reset backoff
/// - If no work, mark thread as idle, sleep and record idle time (if stats enabled), then increase backoff
#[inline]
#[allow(clippy::cast_possible_truncation)]
pub fn handle_worker_backoff<W: HasWorkerCore>(
    worker: &mut W,
    stats: Option<&PipelineStats>,
    did_work: bool,
) {
    if did_work {
        worker.core_mut().reset_backoff();
    } else {
        if let Some(stats) = stats {
            let tid = worker.core().scheduler.thread_id();
            // Mark thread as idle before sleeping so activity snapshots are accurate
            stats.clear_current_step(tid);
            let idle_start = Instant::now();
            worker.core_mut().sleep_backoff();
            stats.record_idle_for_thread(tid, idle_start.elapsed().as_nanos() as u64);
        } else {
            worker.core_mut().sleep_backoff();
        }
        worker.core_mut().increase_backoff();
    }
}

// ============================================================================
// Generic Worker Loop (Consolidated Implementation)
// ============================================================================

/// Trait for pipeline step execution context.
///
/// This trait abstracts over the differences between BAM and FASTQ pipelines,
/// allowing a single `generic_worker_loop` implementation to handle both.
/// Each pipeline provides a context struct that implements this trait.
pub trait StepContext {
    /// The worker type for this pipeline (e.g., `WorkerState<P>` or `FastqWorkerState<P>`)
    type Worker: WorkerStateCommon + HasWorkerCore;

    /// Execute a pipeline step, returning `(success, was_contention)`.
    fn execute_step(&self, worker: &mut Self::Worker, step: PipelineStep) -> (bool, bool);

    /// Compute backpressure state for the scheduler.
    fn get_backpressure(&self, worker: &Self::Worker) -> BackpressureState;

    /// Check and set drain mode if appropriate (called each iteration).
    fn check_drain_mode(&self);

    /// Check if the pipeline has encountered an error.
    fn has_error(&self) -> bool;

    /// Check if the pipeline is complete.
    fn is_complete(&self) -> bool;

    /// Get stats for recording (if enabled).
    fn stats(&self) -> Option<&PipelineStats>;

    /// Whether this context should skip the Read step.
    /// Returns `true` for non-reader worker threads.
    fn skip_read(&self) -> bool;

    /// Whether to check completion at end of loop iteration (original BAM behavior).
    /// If false (default), checks at start (original FASTQ behavior).
    fn check_completion_at_end(&self) -> bool {
        false
    }

    // -------------------------------------------------------------------------
    // Sticky Read Methods (consolidated)
    // -------------------------------------------------------------------------

    /// Fast check before entering sticky read loop.
    /// BAM uses this to skip the loop entirely when `read_done` is true.
    /// Default returns false (no sticky read).
    fn should_attempt_sticky_read(&self) -> bool {
        false
    }

    /// Condition for continuing the sticky read loop.
    /// Called before each read attempt.
    /// - BAM: `!error && !read_done && q1.len() < capacity`
    /// - FASTQ: `true` (relies on `execute_read_step` returning false)
    fn sticky_read_should_continue(&self) -> bool {
        false
    }

    /// Execute the read step. Returns true if read succeeded.
    fn execute_read_step(&self, _worker: &mut Self::Worker) -> bool {
        false
    }

    /// Get drain mode flag for exclusive step relaxation.
    /// Default is false.
    fn is_drain_mode(&self) -> bool {
        false
    }

    /// Check if worker should attempt an exclusive step.
    /// Default implementation always returns true (no ownership checks).
    fn should_attempt_step(
        &self,
        _worker: &Self::Worker,
        _step: PipelineStep,
        _drain_mode: bool,
    ) -> bool {
        true
    }

    /// Get the exclusive step owned by this worker (if any).
    /// Used by BAM to prioritize owned steps before the normal priority loop.
    /// Default returns None (no exclusive step ownership).
    fn exclusive_step_owned(&self, _worker: &Self::Worker) -> Option<PipelineStep> {
        None
    }
}

/// Generic worker loop implementation used by both BAM and FASTQ pipelines.
///
/// This consolidates the duplicated worker loop logic into a single function.
/// The `StepContext` trait provides pipeline-specific behavior.
#[allow(clippy::too_many_lines, clippy::cast_possible_truncation)]
pub fn generic_worker_loop<C: StepContext>(ctx: &C, worker: &mut C::Worker) {
    let collect_stats = ctx.stats().is_some();
    let check_completion_at_end = ctx.check_completion_at_end();

    loop {
        // Check for errors
        if ctx.has_error() {
            break;
        }

        // Check for completion at start (FASTQ behavior, default)
        // CRITICAL: Don't exit while holding items - they would be lost!
        if !check_completion_at_end && ctx.is_complete() && !worker.has_any_held_items() {
            break;
        }

        let mut did_work = false;

        // Sticky read for reader threads (fills Q1 before doing other work)
        // Uses outer guard to skip when not applicable (BAM optimization)
        if ctx.should_attempt_sticky_read() {
            while ctx.sticky_read_should_continue() {
                // Record attempt before executing
                if let Some(stats) = ctx.stats() {
                    stats.record_step_attempt(
                        worker.core().scheduler.thread_id(),
                        PipelineStep::Read,
                    );
                }

                let success = if collect_stats {
                    let start = Instant::now();
                    let success = ctx.execute_read_step(worker);
                    if success {
                        if let Some(stats) = ctx.stats() {
                            stats.record_step_for_thread(
                                PipelineStep::Read,
                                start.elapsed().as_nanos() as u64,
                                Some(worker.core().scheduler.thread_id()),
                            );
                        }
                    }
                    success
                } else {
                    ctx.execute_read_step(worker)
                };

                if success {
                    did_work = true;
                } else {
                    break;
                }
            }
        }

        // Check/set drain mode
        ctx.check_drain_mode();

        // Get step priorities from scheduler
        let backpressure = ctx.get_backpressure(worker);
        let priorities_slice = worker.core_mut().scheduler.get_priorities(backpressure);
        let priority_count = priorities_slice.len().min(9);
        let mut priorities = [PipelineStep::Read; 9];
        priorities[..priority_count].copy_from_slice(&priorities_slice[..priority_count]);

        let drain_mode = ctx.is_drain_mode();

        // BAM optimization: try owned exclusive step first (before priority loop)
        // This prevents starvation: since only this thread can do the step, prioritize it
        let owned_step = ctx.exclusive_step_owned(worker);
        if let Some(step) = owned_step {
            if step != PipelineStep::Read && !ctx.has_error() {
                // Record attempt
                if let Some(stats) = ctx.stats() {
                    stats.record_step_attempt(worker.core().scheduler.thread_id(), step);
                }

                let (success, elapsed_ns, was_contention) = if collect_stats {
                    let start = Instant::now();
                    let (success, was_contention) = ctx.execute_step(worker, step);
                    (success, start.elapsed().as_nanos() as u64, was_contention)
                } else {
                    let (success, was_contention) = ctx.execute_step(worker, step);
                    (success, 0, was_contention)
                };

                worker.core_mut().scheduler.record_outcome(step, success, was_contention);

                if success {
                    if let Some(stats) = ctx.stats() {
                        stats.record_step_for_thread(
                            step,
                            elapsed_ns,
                            Some(worker.core().scheduler.thread_id()),
                        );
                    }
                    did_work = true;
                }
            }
        }

        // Execute steps in priority order (if owned step didn't succeed)
        if !did_work {
            for &step in &priorities[..priority_count] {
                if ctx.has_error() {
                    break;
                }

                // Skip Read step for non-reader workers
                if ctx.skip_read() && step == PipelineStep::Read {
                    continue;
                }

                // Skip the owned step (already tried above)
                if Some(step) == owned_step {
                    continue;
                }

                // Skip exclusive steps this worker doesn't own (BAM optimization)
                if !ctx.should_attempt_step(worker, step, drain_mode) {
                    continue;
                }

                // Record attempt
                if let Some(stats) = ctx.stats() {
                    stats.record_step_attempt(worker.core().scheduler.thread_id(), step);
                }

                // Execute with timing
                let (success, elapsed_ns, was_contention) = if collect_stats {
                    let start = Instant::now();
                    let (success, was_contention) = ctx.execute_step(worker, step);
                    (success, start.elapsed().as_nanos() as u64, was_contention)
                } else {
                    let (success, was_contention) = ctx.execute_step(worker, step);
                    (success, 0, was_contention)
                };

                // Record outcome for adaptive scheduler
                worker.core_mut().scheduler.record_outcome(step, success, was_contention);

                if success {
                    if let Some(stats) = ctx.stats() {
                        stats.record_step_for_thread(
                            step,
                            elapsed_ns,
                            Some(worker.core().scheduler.thread_id()),
                        );
                    }
                    did_work = true;
                    break; // Restart priority evaluation
                }
            }
        }

        // Check for completion at end (original BAM behavior)
        // CRITICAL: Don't exit while holding items - they would be lost!
        if check_completion_at_end && ctx.is_complete() && !worker.has_any_held_items() {
            break;
        }

        // Backoff handling
        handle_worker_backoff(worker, ctx.stats(), did_work);
    }
}

/// Trait for workers that can hold compressed batches when output queue is full.
///
/// This trait enables non-blocking compress steps. When the output queue is full,
/// instead of blocking, the worker holds the compressed batch and returns
/// `StepResult::OutputFull`. The next time the compress step runs, it first
/// tries to advance the held batch before processing new work.
pub trait HasHeldCompressed {
    /// Get mutable reference to the held compressed batch.
    /// Includes `heap_size` for memory tracking.
    fn held_compressed_mut(&mut self) -> &mut Option<(u64, CompressedBlockBatch, usize)>;
}

/// Trait for workers that hold boundary batches when output queue is full.
///
/// IMPORTANT: This pattern must be kept in sync between BAM and FASTQ pipelines.
/// See: bam.rs `try_step_find_boundaries()` and fastq.rs `fastq_try_step_find_boundaries()`
///
/// The pattern:
/// 1. Check/advance held item first (priority 1)
/// 2. Acquire ordering lock (BAM: `boundary_state`, FASTQ: `boundary_lock`)
/// 3. Brief lock for reorder buffer insert/pop
/// 4. Do boundary work (under ordering lock to ensure sequential processing)
/// 5. Push result or hold if output queue full
///
/// Note: FASTQ requires strict ordering due to per-stream leftover state, so it
/// uses a separate `boundary_lock`. BAM uses `boundary_state` which serves both
/// as the ordering lock and contains the boundary-finding state.
pub trait HasHeldBoundaries<B> {
    fn held_boundaries_mut(&mut self) -> &mut Option<(u64, B)>;
}

// ============================================================================
// Shared Step Functions (used by both BAM and FASTQ pipelines)
// ============================================================================

/// Shared Step: Compress serialized data to BGZF blocks (non-blocking).
///
/// This step is parallel - multiple threads can compress concurrently.
/// It takes data from Q5, compresses it using the worker's compressor,
/// and pushes the resulting blocks to Q6.
///
/// # Non-Blocking Behavior
///
/// Instead of blocking when Q6 is full, this function:
/// 1. First tries to advance any held compressed batch from a previous attempt
/// 2. If held batch can't be advanced, returns `OutputFull` immediately
/// 3. Otherwise, processes new work and tries to push the result
/// 4. If push fails, holds the result for the next attempt
///
/// This prevents deadlock by ensuring workers can always return and try
/// other steps (especially Write to drain queues).
pub fn shared_try_step_compress<S, W>(state: &S, worker: &mut W) -> StepResult
where
    S: OutputPipelineState,
    W: HasCompressor + HasHeldCompressed,
{
    // =========================================================================
    // Priority 1: Try to advance any held compressed batch
    // =========================================================================
    if let Some((serial, held, _heap_size)) = worker.held_compressed_mut().take() {
        match state.q6_push((serial, held)) {
            Ok(()) => {
                // Successfully advanced held item, continue to process more
                state.record_q7_push_progress();
            }
            Err((serial, held)) => {
                // Still can't push - put it back and signal output full
                let heap_size = held.estimate_heap_size();
                *worker.held_compressed_mut() = Some((serial, held, heap_size));
                return StepResult::OutputFull;
            }
        }
    }

    // =========================================================================
    // Priority 2: Check if output queue has space (soft check)
    // =========================================================================
    if state.q6_is_full() {
        return StepResult::OutputFull;
    }

    // =========================================================================
    // Priority 3: Pop from Q5
    // =========================================================================
    let Some((serial, serialized)) = state.q5_pop() else {
        if let Some(stats) = state.stats() {
            stats.record_queue_empty(6);
        }
        return StepResult::InputEmpty;
    };
    state.record_q6_pop_progress();

    // Track memory released from Q5
    let q5_heap_size = serialized.estimate_heap_size() as u64;
    state.q5_track_pop(q5_heap_size);

    // =========================================================================
    // Priority 4: Compress the data
    // =========================================================================
    let compressor = worker.compressor_mut();

    if let Err(e) = compressor.write_all(&serialized.data) {
        state.set_error(e);
        return StepResult::InputEmpty;
    }
    if let Err(e) = compressor.flush() {
        state.set_error(e);
        return StepResult::InputEmpty;
    }

    // Take the compressed blocks, preserving the record count
    let blocks = compressor.take_blocks();

    // Record compressed bytes for throughput metrics
    let compressed_bytes: u64 = blocks.iter().map(|b| b.data.len() as u64).sum();
    state.record_compressed_bytes_out(compressed_bytes);

    let batch = CompressedBlockBatch { blocks, record_count: serialized.record_count };

    // =========================================================================
    // Priority 5: Try to push result
    // =========================================================================
    match state.q6_push((serial, batch)) {
        Ok(()) => {
            state.record_q7_push_progress();
            StepResult::Success
        }
        Err((serial, batch)) => {
            // Output full - hold the result for next attempt
            let heap_size = batch.estimate_heap_size();
            *worker.held_compressed_mut() = Some((serial, batch, heap_size));
            StepResult::OutputFull
        }
    }
}

/// Shared Step: Write compressed blocks to output.
///
/// This step is exclusive - only one thread at a time can write.
/// It drains Q6 into a reorder buffer, then writes batches in order.
///
/// Note: Currently unused as both pipelines have their own write functions
/// with specific features (BAM: stats/progress logging, FASTQ: similar).
/// Kept for future unification in Phase 2b.
#[allow(dead_code)]
fn shared_try_step_write<S: OutputPipelineState>(state: &S) -> bool {
    // Try to acquire exclusive access to output file FIRST
    let Some(mut guard) = state.output_try_lock() else {
        return false; // Another thread is writing
    };

    let Some(ref mut writer) = *guard else {
        return false; // File already closed
    };

    // Drain Q6 into reorder buffer
    while let Some((serial, batch)) = state.q6_pop() {
        let q7_heap = batch.estimate_heap_size() as u64;
        state.q6_track_pop(q7_heap);
        state.q6_reorder_insert(serial, batch);
    }

    // Write all ready batches in order
    let mut wrote_any = false;
    while let Some(batch) = state.q6_reorder_try_pop_next() {
        // Write all blocks in the batch
        for block in &batch.blocks {
            if let Err(e) = writer.write_all(&block.data) {
                state.set_error(e);
                return false;
            }
        }
        state.increment_written();
        wrote_any = true;
    }

    wrote_any
}

// ============================================================================
// Held-Item Pattern Helpers
// ============================================================================

/// Try to advance a held item to a queue.
///
/// Returns true if the item was successfully pushed (or nothing was held).
/// Returns false if the queue is still full and the item is still held.
///
/// This is a helper for the non-blocking held-item pattern used to prevent deadlock.
#[inline]
pub fn try_advance_held<T>(queue: &ArrayQueue<(u64, T)>, held: &mut Option<(u64, T)>) -> bool {
    if let Some((serial, item)) = held.take() {
        match queue.push((serial, item)) {
            Ok(()) => true,
            Err(returned) => {
                *held = Some(returned);
                false
            }
        }
    } else {
        true // Nothing held
    }
}

/// Try to push an item to a queue, holding it if the queue is full.
///
/// Returns `Success` if pushed, `OutputFull` if held for later.
#[inline]
pub fn try_push_or_hold<T>(
    queue: &ArrayQueue<(u64, T)>,
    serial: u64,
    item: T,
    held: &mut Option<(u64, T)>,
) -> StepResult {
    match queue.push((serial, item)) {
        Ok(()) => StepResult::Success,
        Err(returned) => {
            *held = Some(returned);
            StepResult::OutputFull
        }
    }
}

// ============================================================================
// Process Step Traits (Step 4)
// ============================================================================

/// Trait for pipeline states that support the Process step.
///
/// This trait abstracts over the queue access patterns for the process step,
/// allowing `shared_try_step_process` to work with both BAM and FASTQ pipelines.
pub trait ProcessPipelineState<G, P>: Send + Sync {
    /// Pop a batch from the input queue (groups/templates).
    fn process_input_pop(&self) -> Option<(u64, Vec<G>)>;

    /// Check if output queue is full.
    fn process_output_is_full(&self) -> bool;

    /// Push processed results to output queue.
    ///
    /// # Errors
    ///
    /// Returns the item if the queue is full.
    fn process_output_push(&self, item: (u64, Vec<P>)) -> Result<(), (u64, Vec<P>)>;

    /// Check if an error has occurred.
    fn has_error(&self) -> bool;

    /// Set an error.
    fn set_error(&self, e: io::Error);

    /// Check if backpressure should be applied before processing new work.
    /// Returns true if queue is full OR memory is high (unless draining).
    /// Default: just checks queue capacity (backwards compatible).
    fn should_apply_process_backpressure(&self) -> bool {
        self.process_output_is_full()
    }

    /// Check if pipeline is in draining mode (completing after input exhausted).
    /// When draining, memory backpressure is bypassed to prevent deadlock.
    /// Default: false (backwards compatible).
    fn is_draining(&self) -> bool {
        false
    }
}

/// Trait for workers that can hold processed batches.
///
/// This enables non-blocking process steps - when the output queue is full,
/// the worker holds the result and returns `OutputFull` instead of blocking.
pub trait HasHeldProcessed<P> {
    /// Get mutable reference to held processed batch.
    /// Includes `heap_size` for memory tracking.
    fn held_processed_mut(&mut self) -> &mut Option<(u64, Vec<P>, usize)>;
}

/// Shared Process step - pops batch, applies `process_fn`, pushes results.
///
/// # Non-Blocking Behavior
///
/// 1. First tries to advance any held batch from a previous attempt
/// 2. If held batch can't be advanced, returns `OutputFull`
/// 3. Otherwise processes new work and tries to push
/// 4. If push fails, holds the result for next attempt
#[inline]
pub fn shared_try_step_process<S, W, G, P, F>(
    state: &S,
    worker: &mut W,
    process_fn: F,
) -> StepResult
where
    S: ProcessPipelineState<G, P>,
    W: HasHeldProcessed<P>,
    P: MemoryEstimate,
    F: Fn(G) -> io::Result<P>,
{
    // Priority 1: Advance held item
    let held = worker.held_processed_mut();
    if let Some((serial, items, _heap_size)) = held.take() {
        match state.process_output_push((serial, items)) {
            Ok(()) => {
                // Successfully pushed - memory tracking handled by trait impl
            }
            Err((serial, items)) => {
                // Re-calculate heap_size for accurate memory tracking
                let heap_size: usize = items.iter().map(MemoryEstimate::estimate_heap_size).sum();
                *held = Some((serial, items, heap_size));
                return StepResult::OutputFull;
            }
        }
    }

    // Priority 2: Check errors
    if state.has_error() {
        return StepResult::InputEmpty;
    }

    // Priority 3: Check backpressure (queue capacity AND memory)
    if state.should_apply_process_backpressure() {
        return StepResult::OutputFull;
    }

    // Priority 4: Pop input
    let Some((serial, batch)) = state.process_input_pop() else {
        return StepResult::InputEmpty;
    };

    // Priority 5: Process items
    let mut results = Vec::with_capacity(batch.len());
    for item in batch {
        match process_fn(item) {
            Ok(processed) => results.push(processed),
            Err(e) => {
                state.set_error(e);
                return StepResult::InputEmpty;
            }
        }
    }

    // Priority 6: Push result
    match state.process_output_push((serial, results)) {
        Ok(()) => StepResult::Success,
        Err((serial, results)) => {
            // Calculate heap_size for accurate memory tracking
            let heap_size: usize = results.iter().map(MemoryEstimate::estimate_heap_size).sum();
            *worker.held_processed_mut() = Some((serial, results, heap_size));
            StepResult::OutputFull
        }
    }
}

// ============================================================================
// Serialize Step Traits (Step 5)
// ============================================================================

/// Trait for pipeline states that support the Serialize step.
///
/// This trait abstracts over the queue access patterns for the serialize step,
/// allowing `shared_try_step_serialize` to work with both BAM and FASTQ pipelines.
pub trait SerializePipelineState<P>: Send + Sync {
    /// Pop a batch from the input queue (processed items).
    fn serialize_input_pop(&self) -> Option<(u64, Vec<P>)>;

    /// Check if output queue is full.
    fn serialize_output_is_full(&self) -> bool;

    /// Push serialized batch to output queue.
    ///
    /// # Errors
    ///
    /// Returns the item if the queue is full.
    fn serialize_output_push(
        &self,
        item: (u64, SerializedBatch),
    ) -> Result<(), (u64, SerializedBatch)>;

    /// Check if an error has occurred.
    fn has_error(&self) -> bool;

    /// Set an error.
    fn set_error(&self, e: io::Error);

    /// Record serialized bytes for throughput metrics (optional).
    fn record_serialized_bytes(&self, _bytes: u64) {}

    /// Record serialized record count for completion tracking (optional).
    fn record_serialized_records(&self, _count: u64) {}
}

/// Trait for workers that can hold serialized batches.
pub trait HasHeldSerialized {
    /// Get mutable reference to held serialized batch.
    /// Includes `heap_size` for memory tracking.
    fn held_serialized_mut(&mut self) -> &mut Option<(u64, SerializedBatch, usize)>;

    /// Get mutable reference to serialization buffer (for reuse).
    fn serialization_buffer_mut(&mut self) -> &mut Vec<u8>;

    /// Get the capacity for the serialization buffer.
    /// BAM uses 64KB, FASTQ uses 256KB.
    fn serialization_buffer_capacity(&self) -> usize;
}

/// Shared Serialize step - pops batch, serializes items, pushes concatenated result.
///
/// Uses a buffer-based serialize function signature for efficiency: the `serialize_fn`
/// writes directly into the provided buffer, avoiding intermediate allocations.
///
/// # Non-Blocking Behavior
///
/// Same pattern as `shared_try_step_process` - holds result if output full.
#[inline]
pub fn shared_try_step_serialize<S, W, P, F>(
    state: &S,
    worker: &mut W,
    mut serialize_fn: F,
) -> StepResult
where
    S: SerializePipelineState<P>,
    W: HasHeldSerialized,
    F: FnMut(P, &mut Vec<u8>) -> io::Result<u64>,
{
    // Priority 1: Advance held item
    if let Some((serial, held_batch, _heap_size)) = worker.held_serialized_mut().take() {
        match state.serialize_output_push((serial, held_batch)) {
            Ok(()) => {
                // Note: Memory tracking would be added here when needed
            }
            Err((serial, held_batch)) => {
                let heap_size = held_batch.estimate_heap_size();
                *worker.held_serialized_mut() = Some((serial, held_batch, heap_size));
                return StepResult::OutputFull;
            }
        }
    }

    // Priority 2: Check errors
    if state.has_error() {
        return StepResult::InputEmpty;
    }

    // Priority 3: Check output space
    if state.serialize_output_is_full() {
        return StepResult::OutputFull;
    }

    // Priority 4: Pop input
    let Some((serial, batch)) = state.serialize_input_pop() else {
        return StepResult::InputEmpty;
    };

    // Get capacity first (before mutable borrow of buffer)
    let capacity = worker.serialization_buffer_capacity();

    // Priority 5: Serialize directly into buffer (no intermediate allocation)
    let buffer = worker.serialization_buffer_mut();
    buffer.clear();
    let mut total_records: u64 = 0;

    for item in batch {
        match serialize_fn(item, buffer) {
            Ok(record_count) => {
                total_records += record_count;
            }
            Err(e) => {
                state.set_error(e);
                return StepResult::InputEmpty;
            }
        }
    }

    // Swap buffer (avoid allocation) - use worker's configured capacity
    let data = std::mem::replace(buffer, Vec::with_capacity(capacity));
    state.record_serialized_bytes(data.len() as u64);
    state.record_serialized_records(total_records);

    let result_batch = SerializedBatch { data, record_count: total_records };

    // Priority 6: Push result
    match state.serialize_output_push((serial, result_batch)) {
        Ok(()) => StepResult::Success,
        Err((serial, result_batch)) => {
            let heap_size = result_batch.estimate_heap_size();
            *worker.held_serialized_mut() = Some((serial, result_batch, heap_size));
            StepResult::OutputFull
        }
    }
}

// ============================================================================
// Write Step Traits (Step 7)
// ============================================================================

/// Trait for pipeline states that support the Write step.
///
/// This trait abstracts over the write step's needs for queue access,
/// reorder buffer, and output file.
pub trait WritePipelineState: Send + Sync {
    /// Get reference to the input queue for write step.
    fn write_input_queue(&self) -> &ArrayQueue<(u64, CompressedBlockBatch)>;

    /// Get reference to the reorder buffer for maintaining output order.
    fn write_reorder_buffer(&self) -> &Mutex<ReorderBuffer<CompressedBlockBatch>>;

    /// Get reference to the output writer.
    fn write_output(&self) -> &Mutex<Option<Box<dyn Write + Send>>>;

    /// Check if an error has occurred.
    fn has_error(&self) -> bool;

    /// Set an error.
    fn set_error(&self, e: io::Error);

    /// Record that records were written.
    fn record_written(&self, count: u64);

    /// Get reference to stats for recording contention (optional).
    fn stats(&self) -> Option<&PipelineStats>;
}

/// Shared Write step - drains queue to reorder buffer, writes in order.
///
/// Returns `Success` if any data was written, `InputEmpty` otherwise.
/// This step is exclusive - uses `try_lock` to avoid blocking.
pub fn shared_try_step_write_new<S: WritePipelineState>(state: &S) -> StepResult {
    if state.has_error() {
        return StepResult::InputEmpty;
    }

    // Try to acquire output lock (non-blocking)
    let Some(mut output_guard) = state.write_output().try_lock() else {
        if let Some(stats) = state.stats() {
            stats.record_contention(PipelineStep::Write);
        }
        return StepResult::InputEmpty;
    };

    let Some(ref mut output) = *output_guard else {
        return StepResult::InputEmpty;
    };

    let mut wrote_any = false;
    {
        let mut reorder = state.write_reorder_buffer().lock();
        let queue = state.write_input_queue();

        // Drain queue into reorder buffer
        while let Some((serial, batch)) = queue.pop() {
            reorder.insert(serial, batch);
        }

        // Write in-order batches
        while let Some(batch) = reorder.try_pop_next() {
            for block in &batch.blocks {
                if let Err(e) = output.write_all(&block.data) {
                    state.set_error(e);
                    return StepResult::InputEmpty;
                }
            }
            state.record_written(batch.record_count);
            wrote_any = true;
        }
    }

    if wrote_any { StepResult::Success } else { StepResult::InputEmpty }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stats_record_step_timing() {
        let stats = PipelineStats::new();

        // Record some step timing
        stats.record_step(PipelineStep::Decompress, 1_000_000); // 1ms
        stats.record_step(PipelineStep::Decompress, 2_000_000); // 2ms

        assert_eq!(stats.step_decompress_ns.load(Ordering::Relaxed), 3_000_000);
        assert_eq!(stats.step_decompress_count.load(Ordering::Relaxed), 2);
    }

    #[test]
    fn test_stats_record_step_for_thread() {
        let stats = PipelineStats::new();

        stats.record_step_for_thread(PipelineStep::Read, 500_000, Some(0));
        stats.record_step_for_thread(PipelineStep::Read, 500_000, Some(1));

        assert_eq!(stats.step_read_ns.load(Ordering::Relaxed), 1_000_000);
        assert_eq!(stats.step_read_count.load(Ordering::Relaxed), 2);
        assert_eq!(stats.per_thread_step_counts[0][0].load(Ordering::Relaxed), 1);
        assert_eq!(stats.per_thread_step_counts[1][0].load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_stats_record_queue_empty() {
        let stats = PipelineStats::new();

        stats.record_queue_empty(1);
        stats.record_queue_empty(1);
        stats.record_queue_empty(2);
        stats.record_queue_empty(25); // Q2b
        stats.record_queue_empty(7);

        assert_eq!(stats.q1_empty.load(Ordering::Relaxed), 2);
        assert_eq!(stats.q2_empty.load(Ordering::Relaxed), 1);
        assert_eq!(stats.q2b_empty.load(Ordering::Relaxed), 1);
        assert_eq!(stats.q7_empty.load(Ordering::Relaxed), 1);
        assert_eq!(stats.q3_empty.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_stats_record_idle_for_thread() {
        let stats = PipelineStats::new();

        stats.record_idle_for_thread(0, 100_000);
        stats.record_idle_for_thread(0, 200_000);
        stats.record_idle_for_thread(1, 50_000);

        assert_eq!(stats.idle_yields.load(Ordering::Relaxed), 3);
        assert_eq!(stats.per_thread_idle_ns[0].load(Ordering::Relaxed), 300_000);
        assert_eq!(stats.per_thread_idle_ns[1].load(Ordering::Relaxed), 50_000);
    }

    // ========================================================================
    // MemoryTracker tests
    // ========================================================================

    #[test]
    fn test_memory_tracker_new() {
        let tracker = MemoryTracker::new(1000);
        assert_eq!(tracker.current(), 0);
        assert_eq!(tracker.peak(), 0);
        assert_eq!(tracker.limit(), 1000);
    }

    #[test]
    fn test_memory_tracker_unlimited() {
        let tracker = MemoryTracker::unlimited();
        assert_eq!(tracker.limit(), 0);
        // Unlimited tracker should always succeed
        assert!(tracker.try_add(1_000_000));
        assert!(tracker.try_add(1_000_000_000));
        assert_eq!(tracker.current(), 1_001_000_000);
    }

    #[test]
    fn test_memory_tracker_try_add_under_limit() {
        let tracker = MemoryTracker::new(1000);
        assert!(tracker.try_add(500));
        assert_eq!(tracker.current(), 500);
    }

    #[test]
    fn test_memory_tracker_try_add_at_limit() {
        let tracker = MemoryTracker::new(1000);
        assert!(tracker.try_add(1000));
        // Now at the limit, next add should be rejected
        assert!(!tracker.try_add(1));
    }

    #[test]
    fn test_memory_tracker_try_add_single_exceeds() {
        // Key behavior: if currently under limit, a single addition that
        // would exceed the limit still succeeds.
        let tracker = MemoryTracker::new(1000);
        assert!(tracker.try_add(500)); // under limit
        assert!(tracker.try_add(600)); // would exceed, but we're under limit so it succeeds
        assert_eq!(tracker.current(), 1100);
        // Now over limit, next add should be rejected
        assert!(!tracker.try_add(1));
    }

    #[test]
    fn test_memory_tracker_remove_saturating() {
        let tracker = MemoryTracker::new(1000);
        tracker.try_add(100);
        // Remove more than current -> saturates to 0
        tracker.remove(200);
        assert_eq!(tracker.current(), 0);
    }

    #[test]
    fn test_memory_tracker_peak_tracking() {
        let tracker = MemoryTracker::new(0); // unlimited
        tracker.try_add(100);
        tracker.try_add(200);
        assert_eq!(tracker.peak(), 300);
        tracker.remove(250);
        assert_eq!(tracker.current(), 50);
        // Peak should still reflect the high-water mark
        assert_eq!(tracker.peak(), 300);
    }

    #[test]
    fn test_memory_tracker_is_at_limit() {
        // Backpressure threshold is min(limit, 512MB).
        // With a limit of 1000, backpressure threshold = 1000.
        let tracker = MemoryTracker::new(1000);
        assert!(!tracker.is_at_limit());
        tracker.try_add(999);
        assert!(!tracker.is_at_limit());
        tracker.try_add(1);
        assert!(tracker.is_at_limit());
    }

    #[test]
    fn test_memory_tracker_drain_threshold() {
        // Drain threshold is half of backpressure threshold.
        // With limit=1000, backpressure=1000, drain=500.
        let tracker = MemoryTracker::new(1000);
        tracker.try_add(1000);
        assert!(!tracker.is_below_drain_threshold()); // at 1000, threshold is 500
        tracker.remove(501);
        assert!(tracker.is_below_drain_threshold()); // at 499, below 500
    }

    #[test]
    fn test_memory_tracker_default_is_unlimited() {
        let tracker = MemoryTracker::default();
        assert_eq!(tracker.limit(), 0);
        // Should behave like unlimited
        assert!(tracker.try_add(1_000_000));
    }

    // ========================================================================
    // ReorderBufferState tests
    // ========================================================================

    #[test]
    fn test_reorder_buffer_state_new() {
        let state = ReorderBufferState::new(1000);
        assert_eq!(state.get_next_seq(), 0);
        assert_eq!(state.get_heap_bytes(), 0);
        assert_eq!(state.get_memory_limit(), 1000);
    }

    #[test]
    fn test_reorder_buffer_state_can_proceed_next_seq() {
        // Always allows the next_seq serial, even if over memory limit
        let state = ReorderBufferState::new(100);
        state.add_heap_bytes(10_000); // way over limit
        // next_seq is 0, so serial 0 should always proceed
        assert!(state.can_proceed(0));
    }

    #[test]
    fn test_reorder_buffer_state_can_proceed_over_limit() {
        // Blocks non-next serials when heap_bytes >= effective_limit / 2
        let state = ReorderBufferState::new(1000);
        // effective_limit = min(1000, 512MB) = 1000
        // backpressure at 50% = 500
        state.add_heap_bytes(500);
        // Serial 1 is not next_seq (0), and heap_bytes >= 500, so should block
        assert!(!state.can_proceed(1));
        // Serial 0 is next_seq, should still proceed
        assert!(state.can_proceed(0));
    }

    #[test]
    fn test_reorder_buffer_state_is_memory_high() {
        let state = ReorderBufferState::new(1000);
        assert!(!state.is_memory_high());
        state.add_heap_bytes(1000);
        assert!(state.is_memory_high());
    }

    #[test]
    fn test_reorder_buffer_state_is_memory_drained() {
        // Hysteresis: enter drain at threshold, exit at half threshold
        let state = ReorderBufferState::new(1000);
        state.add_heap_bytes(1000);
        assert!(!state.is_memory_drained()); // at 1000, threshold/2 = 500
        state.sub_heap_bytes(501);
        assert!(state.is_memory_drained()); // at 499, below 500
    }

    #[test]
    fn test_reorder_buffer_state_effective_limit() {
        // 0 uses BACKPRESSURE_THRESHOLD_BYTES, non-zero uses min
        let state_zero = ReorderBufferState::new(0);
        // is_memory_high checks heap_bytes >= effective_limit
        // With 0 limit, effective_limit = BACKPRESSURE_THRESHOLD_BYTES (512MB)
        assert!(!state_zero.is_memory_high()); // 0 < 512MB

        let state_small = ReorderBufferState::new(100);
        state_small.add_heap_bytes(100);
        assert!(state_small.is_memory_high()); // 100 >= min(100, 512MB) = 100
    }

    #[test]
    fn test_reorder_buffer_state_add_sub_heap_bytes() {
        let state = ReorderBufferState::new(0);
        state.add_heap_bytes(100);
        assert_eq!(state.get_heap_bytes(), 100);
        state.add_heap_bytes(50);
        assert_eq!(state.get_heap_bytes(), 150);
        state.sub_heap_bytes(30);
        assert_eq!(state.get_heap_bytes(), 120);
    }

    // ========================================================================
    // GroupKey tests
    // ========================================================================

    #[test]
    fn test_group_key_single() {
        let key = GroupKey::single(1, 100, 0, 5, 0, 42);
        assert_eq!(key.ref_id1, 1);
        assert_eq!(key.pos1, 100);
        assert_eq!(key.strand1, 0);
        assert_eq!(key.ref_id2, GroupKey::UNKNOWN_REF);
        assert_eq!(key.pos2, GroupKey::UNKNOWN_POS);
        assert_eq!(key.strand2, GroupKey::UNKNOWN_STRAND);
        assert_eq!(key.library_idx, 5);
        assert_eq!(key.cell_hash, 0);
        assert_eq!(key.name_hash, 42);
    }

    #[test]
    fn test_group_key_paired() {
        // When (ref_id, pos, strand) <= (mate_ref_id, mate_pos, mate_strand),
        // the positions stay as given.
        let key = GroupKey::paired(1, 100, 0, 2, 200, 1, 3, 0, 99);
        assert_eq!(key.ref_id1, 1);
        assert_eq!(key.pos1, 100);
        assert_eq!(key.strand1, 0);
        assert_eq!(key.ref_id2, 2);
        assert_eq!(key.pos2, 200);
        assert_eq!(key.strand2, 1);
    }

    #[test]
    fn test_group_key_paired_swap() {
        // When (ref_id, pos, strand) > (mate_ref_id, mate_pos, mate_strand),
        // the positions are swapped so lower comes first.
        let key = GroupKey::paired(5, 500, 1, 1, 100, 0, 3, 0, 99);
        assert_eq!(key.ref_id1, 1);
        assert_eq!(key.pos1, 100);
        assert_eq!(key.strand1, 0);
        assert_eq!(key.ref_id2, 5);
        assert_eq!(key.pos2, 500);
        assert_eq!(key.strand2, 1);
    }

    #[test]
    fn test_group_key_position_key() {
        let key = GroupKey::single(1, 100, 0, 5, 7, 42);
        let pk = key.position_key();
        // position_key returns tuple without name_hash
        assert_eq!(
            pk,
            (
                1,
                100,
                0,
                GroupKey::UNKNOWN_REF,
                GroupKey::UNKNOWN_POS,
                GroupKey::UNKNOWN_STRAND,
                5,
                7
            )
        );
    }

    #[test]
    fn test_group_key_ord_by_position() {
        let key_a = GroupKey::single(1, 100, 0, 0, 0, 0);
        let key_b = GroupKey::single(2, 50, 0, 0, 0, 0);
        // key_a has lower ref_id1, so it should come first
        assert!(key_a < key_b);
    }

    #[test]
    fn test_group_key_ord_tiebreak_name_hash() {
        let key_a = GroupKey::single(1, 100, 0, 0, 0, 10);
        let key_b = GroupKey::single(1, 100, 0, 0, 0, 20);
        // Same position, name_hash breaks tie
        assert!(key_a < key_b);
    }

    #[test]
    fn test_group_key_default() {
        let key = GroupKey::default();
        assert_eq!(key.ref_id1, GroupKey::UNKNOWN_REF);
        assert_eq!(key.pos1, GroupKey::UNKNOWN_POS);
        assert_eq!(key.strand1, GroupKey::UNKNOWN_STRAND);
        assert_eq!(key.ref_id2, GroupKey::UNKNOWN_REF);
        assert_eq!(key.pos2, GroupKey::UNKNOWN_POS);
        assert_eq!(key.strand2, GroupKey::UNKNOWN_STRAND);
        assert_eq!(key.library_idx, 0);
        assert_eq!(key.cell_hash, 0);
        assert_eq!(key.name_hash, 0);
    }

    #[test]
    fn test_group_key_eq() {
        let key_a = GroupKey::single(1, 100, 0, 5, 0, 42);
        let key_b = GroupKey::single(1, 100, 0, 5, 0, 42);
        assert_eq!(key_a, key_b);
    }

    #[test]
    fn test_group_key_paired_same_position() {
        // Same ref_id and pos, different strand normalizes correctly
        let key = GroupKey::paired(1, 100, 1, 1, 100, 0, 0, 0, 0);
        // (1,100,0) < (1,100,1) so mate comes first
        assert_eq!(key.ref_id1, 1);
        assert_eq!(key.pos1, 100);
        assert_eq!(key.strand1, 0);
        assert_eq!(key.strand2, 1);
    }

    #[test]
    fn test_group_key_hash() {
        use std::collections::HashSet;
        let key_a = GroupKey::single(1, 100, 0, 0, 0, 42);
        let key_b = GroupKey::single(1, 100, 0, 0, 0, 42);
        let key_c = GroupKey::single(2, 200, 1, 0, 0, 99);
        let mut set = HashSet::new();
        set.insert(key_a);
        assert!(set.contains(&key_b));
        set.insert(key_c);
        assert_eq!(set.len(), 2);
    }

    // ========================================================================
    // PipelineStep tests
    // ========================================================================

    #[test]
    fn test_pipeline_step_is_exclusive() {
        assert!(PipelineStep::Read.is_exclusive());
        assert!(!PipelineStep::Decompress.is_exclusive());
        assert!(PipelineStep::FindBoundaries.is_exclusive());
        assert!(!PipelineStep::Decode.is_exclusive());
        assert!(PipelineStep::Group.is_exclusive());
        assert!(!PipelineStep::Process.is_exclusive());
        assert!(!PipelineStep::Serialize.is_exclusive());
        assert!(!PipelineStep::Compress.is_exclusive());
        assert!(PipelineStep::Write.is_exclusive());
    }

    #[test]
    fn test_pipeline_step_all() {
        let all = PipelineStep::all();
        assert_eq!(all.len(), 9);
        assert_eq!(all[0], PipelineStep::Read);
        assert_eq!(all[1], PipelineStep::Decompress);
        assert_eq!(all[2], PipelineStep::FindBoundaries);
        assert_eq!(all[3], PipelineStep::Decode);
        assert_eq!(all[4], PipelineStep::Group);
        assert_eq!(all[5], PipelineStep::Process);
        assert_eq!(all[6], PipelineStep::Serialize);
        assert_eq!(all[7], PipelineStep::Compress);
        assert_eq!(all[8], PipelineStep::Write);
    }

    #[test]
    fn test_pipeline_step_from_index() {
        assert_eq!(PipelineStep::from_index(0), PipelineStep::Read);
        assert_eq!(PipelineStep::from_index(1), PipelineStep::Decompress);
        assert_eq!(PipelineStep::from_index(2), PipelineStep::FindBoundaries);
        assert_eq!(PipelineStep::from_index(3), PipelineStep::Decode);
        assert_eq!(PipelineStep::from_index(4), PipelineStep::Group);
        assert_eq!(PipelineStep::from_index(5), PipelineStep::Process);
        assert_eq!(PipelineStep::from_index(6), PipelineStep::Serialize);
        assert_eq!(PipelineStep::from_index(7), PipelineStep::Compress);
        assert_eq!(PipelineStep::from_index(8), PipelineStep::Write);
    }

    #[test]
    fn test_pipeline_step_short_name() {
        assert_eq!(PipelineStep::Read.short_name(), "Rd");
        assert_eq!(PipelineStep::Decompress.short_name(), "Dc");
        assert_eq!(PipelineStep::FindBoundaries.short_name(), "Fb");
        assert_eq!(PipelineStep::Decode.short_name(), "De");
        assert_eq!(PipelineStep::Group.short_name(), "Gr");
        assert_eq!(PipelineStep::Process.short_name(), "Pr");
        assert_eq!(PipelineStep::Serialize.short_name(), "Se");
        assert_eq!(PipelineStep::Compress.short_name(), "Co");
        assert_eq!(PipelineStep::Write.short_name(), "Wr");
        // Verify each is exactly 2 chars
        for step in PipelineStep::all() {
            assert_eq!(step.short_name().len(), 2);
        }
    }

    // ========================================================================
    // StepResult tests
    // ========================================================================

    #[test]
    fn test_step_result_is_success() {
        assert!(StepResult::Success.is_success());
        assert!(!StepResult::OutputFull.is_success());
        assert!(!StepResult::InputEmpty.is_success());
    }

    #[test]
    fn test_step_result_variants() {
        // All three variants exist and are distinct
        let s = StepResult::Success;
        let o = StepResult::OutputFull;
        let i = StepResult::InputEmpty;
        assert_ne!(s, o);
        assert_ne!(s, i);
        assert_ne!(o, i);
    }

    // ========================================================================
    // Batch type tests
    // ========================================================================

    #[test]
    fn test_raw_block_batch_new_empty() {
        let batch = RawBlockBatch::new();
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);
    }

    #[test]
    fn test_raw_block_batch_with_capacity() {
        let batch = RawBlockBatch::with_capacity(32);
        assert!(batch.is_empty());
        assert!(batch.blocks.capacity() >= 32);
    }

    #[test]
    fn test_compressed_block_batch_new() {
        let batch = CompressedBlockBatch::new();
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);
        assert_eq!(batch.record_count, 0);
    }

    #[test]
    fn test_compressed_block_batch_clear() {
        let mut batch = CompressedBlockBatch::new();
        batch.record_count = 42;
        batch.clear();
        assert!(batch.is_empty());
        assert_eq!(batch.record_count, 0);
    }

    #[test]
    fn test_bgzf_batch_config_default() {
        let config = BgzfBatchConfig::default();
        assert_eq!(config.blocks_per_batch, 16);
        assert_eq!(config.compression_level, 6);
    }

    #[test]
    fn test_bgzf_batch_config_new() {
        let config = BgzfBatchConfig::new(64);
        assert_eq!(config.blocks_per_batch, 64);
        // compression_level should still be default
        assert_eq!(config.compression_level, 6);
    }

    #[test]
    fn test_decompressed_batch_new_empty() {
        let batch = DecompressedBatch::new();
        assert!(batch.is_empty());
        assert!(batch.data.is_empty());
    }

    #[test]
    fn test_serialized_batch_clear() {
        let mut batch = SerializedBatch::new();
        batch.data.extend_from_slice(&[1, 2, 3]);
        batch.record_count = 10;
        batch.clear();
        assert!(batch.is_empty());
        assert_eq!(batch.record_count, 0);
    }

    // ========================================================================
    // PipelineConfig tests
    // ========================================================================

    #[test]
    fn test_pipeline_config_new_defaults() {
        let config = PipelineConfig::new(4, 6);
        assert_eq!(config.num_threads, 4);
        assert_eq!(config.compression_level, 6);
        assert_eq!(config.queue_capacity, 64);
        assert_eq!(config.batch_size, 1);
        assert_eq!(config.queue_memory_limit, 0);
        assert!(!config.collect_stats);
    }

    #[test]
    fn test_pipeline_config_builder_chain() {
        let config = PipelineConfig::new(4, 6)
            .with_compression_level(9)
            .with_batch_size(100)
            .with_stats(true)
            .with_queue_memory_limit(1_000_000);
        assert_eq!(config.compression_level, 9);
        assert_eq!(config.batch_size, 100);
        assert!(config.collect_stats);
        assert_eq!(config.queue_memory_limit, 1_000_000);
    }

    #[test]
    fn test_pipeline_config_auto_tuned_1_thread() {
        let config = PipelineConfig::auto_tuned(1, 6);
        assert_eq!(config.num_threads, 1);
        // queue_capacity = (1*16).clamp(64, 256) = 64
        assert_eq!(config.queue_capacity, 64);
    }

    #[test]
    fn test_pipeline_config_auto_tuned_8_threads() {
        let config = PipelineConfig::auto_tuned(8, 6);
        assert_eq!(config.num_threads, 8);
        // queue_capacity = (8*16).clamp(64, 256) = 128
        assert_eq!(config.queue_capacity, 128);
        // blocks_per_read_batch = 32 for >= 8 threads
        assert_eq!(config.blocks_per_read_batch, 32);
    }

    #[test]
    fn test_pipeline_config_auto_tuned_32_threads() {
        let config = PipelineConfig::auto_tuned(32, 6);
        assert_eq!(config.num_threads, 32);
        // queue_capacity = (32*16).clamp(64, 256) = 256 (capped)
        assert_eq!(config.queue_capacity, 256);
    }

    #[test]
    fn test_pipeline_config_with_compression_level() {
        let config = PipelineConfig::new(4, 6).with_compression_level(12);
        assert_eq!(config.compression_level, 12);
    }

    #[test]
    fn test_pipeline_config_with_batch_size_min_1() {
        let config = PipelineConfig::new(4, 6).with_batch_size(0);
        // batch_size of 0 gets clamped to 1
        assert_eq!(config.batch_size, 1);
    }

    #[test]
    fn test_pipeline_config_with_queue_memory_limit() {
        let config = PipelineConfig::new(4, 6).with_queue_memory_limit(500_000_000);
        assert_eq!(config.queue_memory_limit, 500_000_000);
    }

    // ========================================================================
    // PipelineValidationError tests
    // ========================================================================

    #[test]
    fn test_pipeline_validation_error_display_empty() {
        let err = PipelineValidationError {
            non_empty_queues: vec![],
            counter_mismatches: vec![],
            leaked_heap_bytes: 0,
        };
        let display = format!("{err}");
        // Even empty error prints the header
        assert!(display.contains("Pipeline validation failed"));
    }

    #[test]
    fn test_pipeline_validation_error_display_full() {
        let err = PipelineValidationError {
            non_empty_queues: vec!["Q1".to_string(), "Q2".to_string()],
            counter_mismatches: vec!["read_count != write_count".to_string()],
            leaked_heap_bytes: 1024,
        };
        let display = format!("{err}");
        assert!(display.contains("Pipeline validation failed"));
        assert!(display.contains("Q1"));
        assert!(display.contains("Q2"));
        assert!(display.contains("read_count != write_count"));
        assert!(display.contains("1024"));
    }

    // ========================================================================
    // extract_panic_message tests
    // ========================================================================

    #[test]
    fn test_extract_panic_message_str() {
        let payload: Box<dyn std::any::Any + Send> = Box::new("something went wrong");
        let msg = extract_panic_message(payload);
        assert_eq!(msg, "something went wrong");
    }

    #[test]
    fn test_extract_panic_message_string() {
        let payload: Box<dyn std::any::Any + Send> = Box::new(String::from("an error occurred"));
        let msg = extract_panic_message(payload);
        assert_eq!(msg, "an error occurred");
    }

    #[test]
    fn test_extract_panic_message_other() {
        let payload: Box<dyn std::any::Any + Send> = Box::new(42_i32);
        let msg = extract_panic_message(payload);
        assert_eq!(msg, "Unknown panic");
    }

    // ========================================================================
    // WorkerCoreState tests
    // ========================================================================

    #[test]
    fn test_worker_core_state_initial_values() {
        use super::super::scheduler::SchedulerStrategy;
        let state = WorkerCoreState::new(6, 0, 4, SchedulerStrategy::default(), ActiveSteps::all());
        assert_eq!(state.backoff_us, MIN_BACKOFF_US);
    }

    #[test]
    fn test_worker_core_state_reset_backoff() {
        use super::super::scheduler::SchedulerStrategy;
        let mut state =
            WorkerCoreState::new(6, 0, 4, SchedulerStrategy::default(), ActiveSteps::all());
        state.increase_backoff();
        assert!(state.backoff_us > MIN_BACKOFF_US);
        state.reset_backoff();
        assert_eq!(state.backoff_us, MIN_BACKOFF_US);
    }

    #[test]
    fn test_worker_core_state_increase_backoff() {
        use super::super::scheduler::SchedulerStrategy;
        let mut state =
            WorkerCoreState::new(6, 0, 4, SchedulerStrategy::default(), ActiveSteps::all());
        assert_eq!(state.backoff_us, MIN_BACKOFF_US); // 10
        state.increase_backoff();
        assert_eq!(state.backoff_us, MIN_BACKOFF_US * 2); // 20
        state.increase_backoff();
        assert_eq!(state.backoff_us, MIN_BACKOFF_US * 4); // 40
        // Keep increasing until we hit the cap
        for _ in 0..20 {
            state.increase_backoff();
        }
        assert_eq!(state.backoff_us, MAX_BACKOFF_US);
    }

    // ========================================================================
    // OutputPipelineQueues tests
    // ========================================================================

    struct TestProcessed {
        size: usize,
    }

    impl MemoryEstimate for TestProcessed {
        fn estimate_heap_size(&self) -> usize {
            self.size
        }
    }

    #[test]
    fn test_output_queues_new() {
        let output: Box<dyn std::io::Write + Send> = Box::new(Vec::<u8>::new());
        let queues: OutputPipelineQueues<(), TestProcessed> =
            OutputPipelineQueues::new(16, output, None, "test");
        assert!(queues.groups.is_empty());
        assert!(queues.processed.is_empty());
        assert!(queues.serialized.is_empty());
        assert!(queues.compressed.is_empty());
    }

    #[test]
    fn test_output_queues_set_take_error() {
        let output: Box<dyn std::io::Write + Send> = Box::new(Vec::<u8>::new());
        let queues: OutputPipelineQueues<(), TestProcessed> =
            OutputPipelineQueues::new(16, output, None, "test");
        assert!(!queues.has_error());
        queues.set_error(io::Error::other("test error"));
        assert!(queues.has_error());
        let err = queues.take_error();
        assert!(err.is_some());
        assert_eq!(err.unwrap().to_string(), "test error");
    }

    #[test]
    fn test_output_queues_draining() {
        let output: Box<dyn std::io::Write + Send> = Box::new(Vec::<u8>::new());
        let queues: OutputPipelineQueues<(), TestProcessed> =
            OutputPipelineQueues::new(16, output, None, "test");
        assert!(!queues.is_draining());
        queues.set_draining(true);
        assert!(queues.is_draining());
    }

    #[test]
    fn test_output_queues_queue_depths_empty() {
        let output: Box<dyn std::io::Write + Send> = Box::new(Vec::<u8>::new());
        let queues: OutputPipelineQueues<(), TestProcessed> =
            OutputPipelineQueues::new(16, output, None, "test");
        let depths = queues.queue_depths();
        assert_eq!(depths.groups, 0);
        assert_eq!(depths.processed, 0);
        assert_eq!(depths.serialized, 0);
        assert_eq!(depths.compressed, 0);
    }

    #[test]
    fn test_output_queues_are_queues_empty() {
        let output: Box<dyn std::io::Write + Send> = Box::new(Vec::<u8>::new());
        let queues: OutputPipelineQueues<(), TestProcessed> =
            OutputPipelineQueues::new(16, output, None, "test");
        assert!(queues.are_queues_empty());
    }

    // ========================================================================
    // MemoryEstimate impl tests
    // ========================================================================

    #[test]
    fn test_memory_estimate_unit() {
        let unit = ();
        assert_eq!(unit.estimate_heap_size(), 0);
    }

    #[test]
    fn test_decoded_record_record_accessor() {
        let rec = RecordBuf::default();
        let parsed = DecodedRecord::new(rec, GroupKey::default());
        assert!(parsed.record().is_some());
        assert!(parsed.raw_bytes().is_none());

        let raw = DecodedRecord::from_raw_bytes(vec![0u8; 32], GroupKey::default());
        assert!(raw.record().is_none());
        assert!(raw.raw_bytes().is_some());
    }

    #[test]
    fn test_memory_estimate_serialized_batch() {
        let mut batch = SerializedBatch::new();
        batch.data.reserve(1024);
        assert!(batch.estimate_heap_size() >= 1024);
    }

    #[test]
    fn test_memory_estimate_decompressed_batch() {
        let mut batch = DecompressedBatch::new();
        batch.data.reserve(2048);
        assert!(batch.estimate_heap_size() >= 2048);
    }
}
