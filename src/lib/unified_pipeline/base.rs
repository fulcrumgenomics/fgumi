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
    pub fn remove(&self, bytes: usize) {
        self.current_bytes.fetch_sub(bytes as u64, Ordering::Relaxed);
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
        const BACKPRESSURE_CAP: u64 = 512 * 1024 * 1024;
        let backpressure_threshold = if self.limit_bytes == 0 {
            BACKPRESSURE_CAP
        } else {
            self.limit_bytes.min(BACKPRESSURE_CAP)
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
        const BACKPRESSURE_CAP: u64 = 512 * 1024 * 1024;
        let backpressure_threshold = if self.limit_bytes == 0 {
            BACKPRESSURE_CAP
        } else {
            self.limit_bytes.min(BACKPRESSURE_CAP)
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
pub fn write_compressed_batch(
    writer: &mut dyn Write,
    batch: CompressedBlockBatch,
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
            _ => PipelineStep::Write,
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
    /// This is used by `PositionGrouper` to determine if records belong to
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
#[derive(Debug)]
pub struct DecodedRecord {
    /// The decoded BAM record.
    pub record: RecordBuf,
    /// Pre-computed grouping key.
    pub key: GroupKey,
}

impl DecodedRecord {
    /// Create a new decoded record with its grouping key.
    #[must_use]
    pub fn new(record: RecordBuf, key: GroupKey) -> Self {
        Self { record, key }
    }
}

impl MemoryEstimate for DecodedRecord {
    fn estimate_heap_size(&self) -> usize {
        // RecordBuf contains:
        // - name: Option<BString> (Vec<u8>)
        // - sequence: Vec<u8>
        // - quality_scores: Vec<u8>
        // - cigar: Vec<Op> (Op is 4 bytes)
        // - data: Vec<(Tag, Value)> - variable size tags
        //
        // O(1) estimation based on accessible fields:
        let name_size = self.record.name().map_or(0, |n| n.len());
        let seq_len = self.record.sequence().len();
        let qual_len = self.record.quality_scores().len();
        let cigar_ops = self.record.cigar().as_ref().len();

        // Core size (O(1) - no iteration)
        let core_size = name_size + seq_len + qual_len + (cigar_ops * 4);

        // Estimate auxiliary data as ~25% of core size (typical for BAM records)
        // This avoids O(n) iteration over data fields while staying reasonably accurate
        let data_estimate = core_size / 4;

        // Total with 64 bytes overhead for struct/Vec allocations
        core_size + data_estimate + 64
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
        // RecordBuf contains:
        // - name: Option<BString> (Vec<u8>)
        // - sequence: Sequence (Vec<u8>)
        // - quality_scores: QualityScores (Vec<u8>)
        // - cigar: Cigar (Vec<Op>)
        // - data: Data (IndexMap of tags to values)
        //
        // O(1) estimation based on accessible fields:
        let name_size = self.name().map_or(0, |n| n.len());
        let seq_len = self.sequence().len();
        let qual_len = self.quality_scores().len();
        let cigar_ops = self.cigar().as_ref().len();

        // Core size (O(1) - no iteration)
        let core_size = name_size + seq_len + qual_len + (cigar_ops * 4);

        // Estimate auxiliary data as ~25% of core size (typical for BAM records)
        // This avoids O(n) iteration over data fields while staying reasonably accurate
        let data_estimate = core_size / 4;

        // Total with 64 bytes overhead for struct/Vec allocations
        core_size + data_estimate + 64
    }
}

impl MemoryEstimate for Vec<RecordBuf> {
    fn estimate_heap_size(&self) -> usize {
        self.iter().map(MemoryEstimate::estimate_heap_size).sum()
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
    /// Tag used for cell barcode extraction.
    pub cell_tag: Tag,
}

impl GroupKeyConfig {
    /// Create a new `GroupKeyConfig`.
    #[must_use]
    pub fn new(library_index: LibraryIndex, cell_tag: Tag) -> Self {
        Self { library_index: Arc::new(library_index), cell_tag }
    }
}

impl Default for GroupKeyConfig {
    fn default() -> Self {
        Self {
            library_index: Arc::new(LibraryIndex::default()),
            cell_tag: Tag::from([b'C', b'B']), // Default cell barcode tag
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
    /// Current heap bytes in groups queue.
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
    pub stats: Option<PipelineStats>,
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
        stats: Option<PipelineStats>,
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
    #[must_use]
    pub fn new(
        compression_level: u32,
        thread_id: usize,
        num_threads: usize,
        scheduler_strategy: SchedulerStrategy,
    ) -> Self {
        Self {
            compressor: InlineBgzfCompressor::new(compression_level),
            scheduler: create_scheduler(scheduler_strategy, thread_id, num_threads),
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
    fn flush_output(&self) -> io::Result<()>;
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
    state.set_error(io::Error::other(format!("Worker thread {} panicked: {}", thread_id, msg)));
}

// ============================================================================
// Pipeline Finalization Helpers
// ============================================================================

/// Join all worker threads, waiting for completion.
///
/// Returns an error if any thread panicked without setting an error on the state.
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
pub fn finalize_pipeline<S: PipelineLifecycle>(state: &S) -> io::Result<u64> {
    // Check for errors
    if let Some(error) = state.take_error() {
        return Err(error);
    }

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
        self.stats.as_ref()
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
        let heap_size = item.1.estimate_heap_size();
        let result = self.compressed.push(item);
        if result.is_ok() {
            self.compressed_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
        }
        result
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
        self.stats.as_ref()
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
    ///   (e.g., `group` with `PositionGroups` containing many records).
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
            writeln!(s, "  Count:   {:>10}", batch_count).unwrap();
            writeln!(s, "  Min:     {:>10}", batch_min).unwrap();
            writeln!(s, "  Max:     {:>10}", batch_max).unwrap();
            writeln!(s, "  Average: {:>10.1}", batch_avg).unwrap();
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
                write!(s, " {:>6}", name).unwrap();
            }
            writeln!(s, "    Idle ms").unwrap();

            // Per-thread rows
            for tid in 0..num_threads.min(MAX_THREADS) {
                write!(s, "  T{:<5} ", tid).unwrap();
                for step_idx in 0..NUM_STEPS {
                    let count = self.per_thread_step_counts[tid][step_idx].load(Ordering::Relaxed);
                    write!(s, " {:>6}", count).unwrap();
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
                write!(s, " {:>6}", total).unwrap();
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
                write!(s, " {:>6}", name).unwrap();
            }
            writeln!(s, "   Total%").unwrap();

            // Per-thread rows with success rates
            for tid in 0..num_threads.min(MAX_THREADS) {
                write!(s, "  T{:<5} ", tid).unwrap();
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
                        write!(s, " {:>5.0}%", rate).unwrap();
                    }
                }
                if thread_attempts == 0 {
                    writeln!(s, "      -").unwrap();
                } else {
                    let total_rate = (thread_successes as f64 / thread_attempts as f64) * 100.0;
                    writeln!(s, "  {:>5.1}%", total_rate).unwrap();
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
                writeln!(s, "  Work time:       {:>10.1}ms", work_ms).unwrap();
                writeln!(s, "  Idle time:       {:>10.1}ms", idle_ms).unwrap();
                writeln!(s, "  Utilization:     {:>10.1}%", utilization).unwrap();
                writeln!(s, "  Contention attempts: {:>7}", total_contention).unwrap();
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
                    format!("{} B", bytes)
                }
            };

            // Helper to format count with K/M suffix
            let format_count = |count: u64| -> String {
                if count >= 1_000_000 {
                    format!("{:.2}M", count as f64 / 1_000_000.0)
                } else if count >= 1_000 {
                    format!("{:.1}K", count as f64 / 1_000.0)
                } else {
                    format!("{}", count)
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
                        write!(s, "{}", short).unwrap();
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
            writeln!(s, "Peak Q3 Reorder Buffer: {} items, {:.1} MB", peak_r3_items, peak_r3_mb)
                .unwrap();
        }

        // Memory limiting statistics
        let group_rejects = self.group_memory_rejects.load(Ordering::Relaxed);
        let peak_memory = self.peak_memory_bytes.load(Ordering::Relaxed);

        if group_rejects > 0 || peak_memory > 0 {
            writeln!(s).unwrap();
            writeln!(s, "Memory Limiting:").unwrap();
            if group_rejects > 0 {
                writeln!(s, "  Group rejects (memory): {:>10}", group_rejects).unwrap();
            }
            if peak_memory > 0 {
                let peak_mb = peak_memory as f64 / 1_048_576.0;
                writeln!(s, "  Peak memory usage:      {:>10.1} MB", peak_mb).unwrap();
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
    fn q5_push(&self, item: (u64, SerializedBatch)) -> Result<(), (u64, SerializedBatch)>;
    fn q5_is_full(&self) -> bool;
    /// Track memory released when popping from Q5.
    fn q5_track_pop(&self, _heap_size: u64) {}

    // Queue 6 access (Compress → Write)
    fn q6_pop(&self) -> Option<(u64, CompressedBlockBatch)>;
    fn q6_push(&self, item: (u64, CompressedBlockBatch))
    -> Result<(), (u64, CompressedBlockBatch)>;
    fn q6_is_full(&self) -> bool;

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
                let heap_size: usize = items.iter().map(|i| i.estimate_heap_size()).sum();
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
            let heap_size: usize = results.iter().map(|i| i.estimate_heap_size()).sum();
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
}
