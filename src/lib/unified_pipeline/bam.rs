//! BAM pipeline: Read → Decompress → `FindBoundaries` → Decode → Group → Process → Serialize → Compress → Write.
//!
//! This module implements the BAM-specific pipeline for processing BAM files
//! with grouping operations like `group` and `codec`.

use crossbeam_queue::ArrayQueue;
use noodles::bam::{self};
use noodles::sam::{Header, alignment::RecordBuf};
use parking_lot::Mutex;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::thread;
use std::time::{Duration, Instant};

use crate::bgzf_reader::{BGZF_EOF, decompress_block_into, read_raw_blocks};
use crate::bgzf_writer::InlineBgzfCompressor;
use crate::progress::ProgressTracker;
use crate::reorder_buffer::ReorderBuffer;

use super::base::{
    BatchWeight, CompressedBlockBatch, DecodedRecord, DecompressedBatch, GroupKeyConfig,
    HasCompressor, HasHeldBoundaries, HasHeldCompressed, HasHeldProcessed, HasHeldSerialized,
    HasWorkerCore, MemoryEstimate, MonitorableState, OutputPipelineQueues, OutputPipelineState,
    PROGRESS_LOG_INTERVAL, PipelineConfig, PipelineLifecycle, PipelineStats, PipelineStep,
    ProcessPipelineState, QueueSample, RawBlockBatch, ReorderBufferState, SerializePipelineState,
    SerializedBatch, StepContext, WorkerCoreState, WorkerStateCommon, WritePipelineState,
    finalize_pipeline, generic_worker_loop, handle_worker_panic, join_monitor_thread,
    join_worker_threads, shared_try_step_compress,
};
use super::deadlock::{DeadlockConfig, DeadlockState, QueueSnapshot, check_deadlock_and_restore};
use super::scheduler::{BackpressureState, SchedulerStrategy};
use crate::read_info::{LibraryIndex, compute_group_key};

/// Buffer size for buffered I/O (8 MB).
/// This reduces syscalls by batching reads/writes into larger chunks.
const IO_BUFFER_SIZE: usize = 8 * 1024 * 1024;

/// Default target templates per batch for template-based batching.
/// Groups are accumulated until the total template count reaches this threshold.
/// This provides consistent batch sizes regardless of templates-per-group variation.
pub const DEFAULT_TARGET_TEMPLATES_PER_BATCH: usize = 500;

// ============================================================================
// Boundary Finding Types (for 8-step pipeline)
// ============================================================================

/// Output of `FindBoundaries` step: buffer + record offsets for parallel decoding.
///
/// This struct enables parallel BAM record decoding by pre-computing where
/// each record starts in the decompressed data. The actual parsing/decoding
/// can then be parallelized across multiple threads.
#[derive(Debug, Clone)]
pub struct BoundaryBatch {
    /// The decompressed bytes (with leftover prepended, suffix removed).
    pub buffer: Vec<u8>,
    /// Byte offsets where each record starts (offsets into buffer).
    /// Length = `num_records` + 1 (last entry is `buffer.len()` for easy slicing).
    pub offsets: Vec<usize>,
}

/// State for the `FindBoundaries` step (sequential).
///
/// This state maintains leftover bytes from incomplete records that span
/// across BGZF block boundaries. The boundary finding is very fast (~0.1μs
/// per block) since it only reads 4-byte integers without decoding records.
///
/// Uses a reusable work buffer to minimize allocations on the hot path.
pub struct BoundaryState {
    /// Leftover bytes from previous block (incomplete record at end).
    leftover: Vec<u8>,
    /// Reusable working buffer to avoid per-call allocations.
    work_buffer: Vec<u8>,
    /// Whether the BAM header has been skipped.
    header_skipped: bool,
}

impl BoundaryState {
    /// Create a new boundary state.
    #[must_use]
    pub fn new() -> Self {
        Self { leftover: Vec::new(), work_buffer: Vec::new(), header_skipped: false }
    }

    /// Create a new boundary state that doesn't skip the header.
    /// Use this when the input stream is already positioned past the header.
    #[must_use]
    pub fn new_no_header() -> Self {
        Self { leftover: Vec::new(), work_buffer: Vec::new(), header_skipped: true }
    }

    /// Parse BAM header and return the number of bytes consumed.
    /// Returns None if more data is needed.
    fn parse_header_size(data: &[u8]) -> Option<usize> {
        // BAM header structure:
        // - magic: 4 bytes ("BAM\1")
        // - l_text: 4 bytes (header text length)
        // - text: l_text bytes
        // - n_ref: 4 bytes (number of references)
        // - for each reference:
        //   - l_name: 4 bytes
        //   - name: l_name bytes
        //   - l_ref: 4 bytes

        if data.len() < 8 {
            return None;
        }

        // Check magic
        if &data[0..4] != b"BAM\x01" {
            // Not a valid BAM file, but let's not error here
            // Just return 0 so records start immediately
            return Some(0);
        }

        let l_text = u32::from_le_bytes([data[4], data[5], data[6], data[7]]) as usize;
        let mut offset = 8 + l_text;

        if data.len() < offset + 4 {
            return None;
        }

        let n_ref = u32::from_le_bytes([
            data[offset],
            data[offset + 1],
            data[offset + 2],
            data[offset + 3],
        ]) as usize;
        offset += 4;

        // Parse each reference
        for _ in 0..n_ref {
            if data.len() < offset + 4 {
                return None;
            }
            let l_name = u32::from_le_bytes([
                data[offset],
                data[offset + 1],
                data[offset + 2],
                data[offset + 3],
            ]) as usize;
            offset += 4 + l_name + 4; // l_name + name + l_ref

            if data.len() < offset {
                return None;
            }
        }

        Some(offset)
    }

    /// Find record boundaries in decompressed data.
    ///
    /// This is FAST (~0.1μs per block) because it only scans 4-byte integers
    /// to find where records start - no actual record decoding is performed.
    ///
    /// # Arguments
    ///
    /// * `decompressed` - Decompressed bytes from one or more BGZF blocks
    ///
    /// # Returns
    ///
    /// A `BoundaryBatch` containing the complete records and their offsets.
    /// Any incomplete record at the end is saved as leftover for the next call.
    pub fn find_boundaries(&mut self, decompressed: &[u8]) -> io::Result<BoundaryBatch> {
        // Step 1: Combine leftover with new data into reusable work_buffer
        // This avoids allocating a new Vec on every call
        self.work_buffer.clear();
        if !self.leftover.is_empty() {
            self.work_buffer.append(&mut self.leftover);
        }
        self.work_buffer.extend_from_slice(decompressed);

        // Step 2: Skip header if not already done
        let mut cursor = 0usize;
        if !self.header_skipped {
            if let Some(header_size) = Self::parse_header_size(&self.work_buffer) {
                cursor = header_size;
                self.header_skipped = true;
            } else {
                // Not enough data to parse header, save as leftover and return empty batch
                std::mem::swap(&mut self.leftover, &mut self.work_buffer);
                return Ok(BoundaryBatch { buffer: Vec::new(), offsets: vec![0] });
            }
        }

        // Step 3: Scan for record boundaries (FAST - just read integers)
        let start_cursor = cursor;
        let mut offsets = vec![0usize]; // First offset is 0 (relative to start of records)

        while cursor + 4 <= self.work_buffer.len() {
            let block_size = u32::from_le_bytes([
                self.work_buffer[cursor],
                self.work_buffer[cursor + 1],
                self.work_buffer[cursor + 2],
                self.work_buffer[cursor + 3],
            ]) as usize;

            let record_end = cursor + 4 + block_size;
            if record_end > self.work_buffer.len() {
                break; // Incomplete record - becomes leftover
            }

            cursor = record_end;
            // Offset is relative to start of records (after header)
            offsets.push(cursor - start_cursor);
        }

        // Step 4: Save leftover for next block (reuse allocation)
        // Split work_buffer: [0..start_cursor | start_cursor..cursor | cursor..]
        //                     header (discard) | records (output)    | leftover
        self.leftover.clear();
        self.leftover.extend_from_slice(&self.work_buffer[cursor..]);

        // Extract the records buffer - this allocation is unavoidable as we return ownership
        let buffer = self.work_buffer[start_cursor..cursor].to_vec();

        // Validate: verify each record's block_size matches the offset difference
        #[cfg(debug_assertions)]
        for i in 0..offsets.len().saturating_sub(1) {
            let start = offsets[i];
            let end = offsets[i + 1];
            if end > start + 4 {
                let stored = u32::from_le_bytes([
                    buffer[start],
                    buffer[start + 1],
                    buffer[start + 2],
                    buffer[start + 3],
                ]) as usize;
                let expected = end - start - 4;
                debug_assert_eq!(
                    stored, expected,
                    "find_boundaries: block_size mismatch at record {i}: stored={stored}, expected={expected}"
                );
            }
        }

        Ok(BoundaryBatch { buffer, offsets })
    }

    /// Call at EOF to get any remaining leftover.
    ///
    /// This validates that any remaining bytes form complete records.
    /// If there are incomplete bytes at EOF, an error is returned.
    pub fn finish(&mut self) -> io::Result<Option<BoundaryBatch>> {
        if self.leftover.is_empty() {
            return Ok(None);
        }

        // Try to parse remaining leftover
        let mut offsets = vec![0usize];
        let mut cursor = 0usize;

        while cursor + 4 <= self.leftover.len() {
            let block_size = u32::from_le_bytes([
                self.leftover[cursor],
                self.leftover[cursor + 1],
                self.leftover[cursor + 2],
                self.leftover[cursor + 3],
            ]) as usize;

            let record_end = cursor + 4 + block_size;
            if record_end > self.leftover.len() {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!(
                        "Incomplete BAM record at EOF: need {} bytes, have {}",
                        record_end - cursor,
                        self.leftover.len() - cursor
                    ),
                ));
            }

            cursor = record_end;
            offsets.push(cursor);
        }

        if cursor == 0 {
            return Ok(None);
        }

        Ok(Some(BoundaryBatch { buffer: std::mem::take(&mut self.leftover), offsets }))
    }
}

impl Default for BoundaryState {
    fn default() -> Self {
        Self::new()
    }
}

/// Decode BAM records from a boundary batch (parallel step).
///
/// This function takes pre-computed record boundaries and decodes the actual
/// BAM records. Since the boundaries are known, this can be called in parallel
/// on different batches.
///
/// # Arguments
///
/// * `batch` - A `BoundaryBatch` with record offsets
/// * `group_key_config` - Config for computing `GroupKey` (library index and cell tag)
///
/// # Returns
///
/// A vector of decoded `DecodedRecord` instances (record + pre-computed `GroupKey`).
pub fn decode_records(
    batch: &BoundaryBatch,
    group_key_config: &GroupKeyConfig,
) -> io::Result<Vec<DecodedRecord>> {
    use noodles::bam::record::codec::decode;

    let num_records = batch.offsets.len().saturating_sub(1);
    let mut records = Vec::with_capacity(num_records);

    for i in 0..num_records {
        let start = batch.offsets[i];
        let end = batch.offsets[i + 1];

        // Validate: end > start + 4 (need at least block_size prefix)
        if end <= start + 4 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Invalid record bounds: start={start}, end={end}, record_index={i}, \
                     num_records={num_records}, buffer_len={}",
                    batch.buffer.len()
                ),
            ));
        }

        // Validate: block_size in buffer matches offset difference
        let stored_block_size = u32::from_le_bytes([
            batch.buffer[start],
            batch.buffer[start + 1],
            batch.buffer[start + 2],
            batch.buffer[start + 3],
        ]) as usize;
        let expected_block_size = end - start - 4;
        if stored_block_size != expected_block_size {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Block size mismatch: stored={stored_block_size}, expected={expected_block_size}, \
                     record_index={i}, start={start}, end={end}, buffer_len={}",
                    batch.buffer.len()
                ),
            ));
        }

        // Skip the 4-byte block_size prefix
        let record_data = &batch.buffer[start + 4..end];

        let mut record = RecordBuf::default();
        decode(&mut &record_data[..], &mut record)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let key =
            compute_group_key(&record, &group_key_config.library_index, group_key_config.cell_tag);
        records.push(DecodedRecord::new(record, key));
    }

    Ok(records)
}

// ============================================================================
// BAM Pipeline State
// ============================================================================

/// Shared state for the BAM pipeline.
///
/// Generic parameters:
/// - `G`: Group type (output of `Group` step, input to Process step)
/// - `P`: Processed type (output of Process step, input to Serialize step).
///   Must implement `MemoryEstimate` for queue memory tracking.
pub struct BamPipelineState<G, P: MemoryEstimate> {
    /// Pipeline configuration.
    pub config: PipelineConfig,

    // ========== Step 1: Read ==========
    /// Input file, mutex-protected for exclusive access.
    pub input_file: Mutex<Option<Box<dyn Read + Send>>>,
    /// Flag indicating EOF has been reached.
    pub read_done: AtomicBool,
    /// Next serial number to assign when reading.
    pub next_read_serial: AtomicU64,

    // ========== Queue 1: Read → Decompress ==========
    /// Raw BGZF blocks waiting to be decompressed.
    pub q1_raw_blocks: ArrayQueue<(u64, RawBlockBatch)>,

    // ========== Step 2: Decompress (parallel) ==========
    // No state needed - each thread has its own Decompressor
    /// Flag indicating all decompression is done.
    pub decompress_done: AtomicBool,
    /// Count of batches that have completed decompression (for completion tracking).
    pub batches_decompressed: AtomicU64,

    // ========== Queue 2: Decompress → FindBoundaries (with reorder) ==========
    /// Decompressed data waiting for boundary finding.
    pub q2_decompressed: ArrayQueue<(u64, DecompressedBatch)>,
    /// Reorder buffer to ensure step 3 receives data in order.
    pub q2_reorder: Mutex<ReorderBuffer<DecompressedBatch>>,

    // ========== Q2 Reorder Buffer Atomic State (for lock-free admission control) ==========
    /// Atomic state for Q2 reorder buffer (`next_seq` and `heap_bytes`).
    /// Used by Decompress and `FindBoundaries` steps for memory backpressure.
    pub q2_reorder_state: ReorderBufferState,

    // ========== Step 3: FindBoundaries (exclusive, FAST) ==========
    /// Boundary finding state, mutex-protected for exclusive access.
    pub boundary_state: Mutex<BoundaryState>,
    /// Flag indicating all boundaries have been found.
    pub boundary_done: AtomicBool,
    /// Next serial number for boundary batches.
    pub next_boundary_serial: AtomicU64,
    /// Count of batches that have completed boundary finding (for completion tracking).
    pub batches_boundary_found: AtomicU64,
    /// Count of batches that `FindBoundaries` has processed (popped from `q2_reorder`).
    /// Used for completion tracking - `FindBoundaries` only finishes when
    /// `batches_boundary_processed == batches_decompressed`.
    pub batches_boundary_processed: AtomicU64,

    // ========== Queue 2b: FindBoundaries → Decode ==========
    /// Boundary batches waiting to be decoded.
    pub q2b_boundaries: ArrayQueue<(u64, BoundaryBatch)>,

    // ========== Step 4: Decode (parallel) ==========
    /// Flag indicating all decoding is done.
    pub decode_done: AtomicBool,
    /// Count of batches that have completed decoding (for completion tracking).
    pub batches_decoded: AtomicU64,
    /// Configuration for computing `GroupKey` during decode.
    pub group_key_config: GroupKeyConfig,

    // ========== Queue 3: Decode → Group (with reorder) ==========
    /// Decoded records waiting to be grouped.
    pub q3_decoded: ArrayQueue<(u64, Vec<DecodedRecord>)>,
    /// Reorder buffer to ensure step 5 receives records in order.
    pub q3_reorder: Mutex<ReorderBuffer<Vec<DecodedRecord>>>,

    // ========== Q3 Reorder Buffer Atomic State (for lock-free admission control) ==========
    /// Atomic state for Q3 reorder buffer (`next_seq` and `heap_bytes`).
    /// Used by Decode and Group steps for memory backpressure.
    pub q3_reorder_state: ReorderBufferState,
    /// Whether the reorder buffer can currently pop (mirrors `q3_reorder.can_pop()`).
    /// Updated by Group step after inserting/popping from reorder buffer.
    pub q3_reorder_can_pop: AtomicBool,

    // ========== Step 5: Group (exclusive) ==========
    /// Flag indicating all grouping is done.
    pub group_done: AtomicBool,
    /// Next serial number for output groups.
    pub next_group_serial: AtomicU64,
    /// Count of batches that have been processed by Group step (popped from `q3_reorder`).
    /// Used for completion tracking - Group only finishes when `batches_grouped == batches_boundary_found`.
    pub batches_grouped: AtomicU64,

    // ========== Output-Half State (Group → Process → Serialize → Compress → Write) ==========
    /// Shared output pipeline queues and state.
    pub output: OutputPipelineQueues<G, P>,

    // ========== Deadlock Detection ==========
    /// Deadlock detection and recovery state.
    pub deadlock_state: DeadlockState,
}

impl<G: Send, P: Send + MemoryEstimate> BamPipelineState<G, P> {
    /// Create a new pipeline state.
    #[must_use]
    pub fn new(
        config: PipelineConfig,
        input: Box<dyn Read + Send>,
        output: Box<dyn Write + Send>,
        group_key_config: GroupKeyConfig,
    ) -> Self {
        let cap = config.queue_capacity;
        let memory_limit = config.queue_memory_limit;
        let stats = if config.collect_stats { Some(PipelineStats::new()) } else { None };
        // Create boundary state based on whether header was already read
        let boundary_state = if config.header_already_read {
            BoundaryState::new_no_header()
        } else {
            BoundaryState::new()
        };
        // Create deadlock state from config
        let deadlock_config =
            DeadlockConfig::new(config.deadlock_timeout_secs, config.deadlock_recover_enabled);
        let deadlock_state = DeadlockState::new(&deadlock_config, memory_limit);
        Self {
            config,
            // Step 1: Read
            input_file: Mutex::new(Some(input)),
            read_done: AtomicBool::new(false),
            next_read_serial: AtomicU64::new(0),
            // Q1: Read → Decompress
            q1_raw_blocks: ArrayQueue::new(cap),
            // Step 2: Decompress
            decompress_done: AtomicBool::new(false),
            batches_decompressed: AtomicU64::new(0),
            // Q2: Decompress → FindBoundaries (with reorder)
            q2_decompressed: ArrayQueue::new(cap),
            q2_reorder: Mutex::new(ReorderBuffer::new()),
            // Q2 reorder buffer atomic state (for lock-free admission control)
            // Note: Q2 uses the same memory_limit as Q3 for backpressure
            q2_reorder_state: ReorderBufferState::new(memory_limit),
            // Step 3: FindBoundaries
            boundary_state: Mutex::new(boundary_state),
            boundary_done: AtomicBool::new(false),
            next_boundary_serial: AtomicU64::new(0),
            batches_boundary_found: AtomicU64::new(0),
            batches_boundary_processed: AtomicU64::new(0),
            // Q2b: FindBoundaries → Decode
            q2b_boundaries: ArrayQueue::new(cap),
            // Step 4: Decode
            decode_done: AtomicBool::new(false),
            batches_decoded: AtomicU64::new(0),
            group_key_config,
            // Q3: Decode → Group (with reorder)
            q3_decoded: ArrayQueue::new(cap),
            q3_reorder: Mutex::new(ReorderBuffer::new()),
            // Q3 reorder buffer atomic state (for lock-free admission control)
            q3_reorder_state: ReorderBufferState::new(memory_limit),
            q3_reorder_can_pop: AtomicBool::new(false),
            // Step 5: Group
            group_done: AtomicBool::new(false),
            next_group_serial: AtomicU64::new(0),
            batches_grouped: AtomicU64::new(0),
            // Output-half state (Group → Process → Serialize → Compress → Write)
            output: OutputPipelineQueues::new(cap, output, stats, "Processed records"),
            // Deadlock detection
            deadlock_state,
        }
    }

    /// Record an error and signal threads to stop.
    pub fn set_error(&self, error: io::Error) {
        self.output.set_error(error);
    }

    /// Check if an error has occurred.
    #[must_use]
    pub fn has_error(&self) -> bool {
        self.output.has_error()
    }

    /// Take the stored error.
    pub fn take_error(&self) -> Option<io::Error> {
        self.output.take_error()
    }

    /// Check if the pipeline is complete.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        // First check atomic flags - all stages must be done
        if !self.read_done.load(Ordering::Acquire) || !self.group_done.load(Ordering::Acquire) {
            return false;
        }

        // Check input-half queues
        if !self.q1_raw_blocks.is_empty()
            || !self.q2_decompressed.is_empty()
            || !self.q2b_boundaries.is_empty()
            || !self.q3_decoded.is_empty()
        {
            return false;
        }

        // Check input-half reorder buffers
        let q2_empty = self.q2_reorder.lock().is_empty();
        let q3_empty = self.q3_reorder.lock().is_empty();
        if !q2_empty || !q3_empty {
            return false;
        }

        // Delegate output-half check
        self.output.are_queues_empty()
    }

    /// Get queue lengths for priority scheduling.
    #[must_use]
    pub fn queue_depths(&self) -> QueueDepths {
        let output_depths = self.output.queue_depths();
        QueueDepths {
            q1: self.q1_raw_blocks.len(),
            q2: self.q2_decompressed.len(),
            q2b: self.q2b_boundaries.len(),
            q3: self.q3_decoded.len(),
            q4: output_depths.groups,
            q5: output_depths.processed,
            q6: output_depths.serialized,
            q7: output_depths.compressed,
        }
    }

    /// Check if Decompress step can proceed with pushing a batch to Q2.
    ///
    /// This implements memory-based backpressure on the Q2 reorder buffer to prevent
    /// unbounded memory growth when `FindBoundaries` (exclusive step) falls behind.
    ///
    /// # Deadlock Prevention
    ///
    /// Always allows the serial that `FindBoundaries` needs (`next_seq`) to proceed,
    /// even if over memory limit. This ensures `FindBoundaries` can always make progress.
    #[must_use]
    pub fn can_decompress_proceed(&self, serial: u64) -> bool {
        // Delegate to Q2 reorder state's can_proceed method
        self.q2_reorder_state.can_proceed(serial)
    }

    /// Check if Decode step can proceed with pushing decoded records to Q3.
    ///
    /// This implements memory-based backpressure on the Q3 reorder buffer to prevent
    /// unbounded memory growth when Group (exclusive step) falls behind.
    ///
    /// # Deadlock Prevention
    ///
    /// Always allows the serial that Group needs (`next_seq`) to proceed, even if
    /// over memory limit. This ensures Group can always make progress.
    #[must_use]
    pub fn can_decode_proceed(&self, serial: u64) -> bool {
        // Delegate to Q3 reorder state's can_proceed method
        self.q3_reorder_state.can_proceed(serial)
    }

    /// Check if memory is at the backpressure threshold.
    ///
    /// Uses Q3 reorder buffer tracking (before Group step) to signal memory pressure
    /// to the scheduler. See [`BACKPRESSURE_THRESHOLD_BYTES`] for architecture details.
    #[must_use]
    pub fn is_memory_high(&self) -> bool {
        self.q3_reorder_state.is_memory_high()
    }

    /// Check if memory has drained below the low-water mark.
    ///
    /// Provides hysteresis to prevent thrashing: enter drain mode at backpressure
    /// threshold, only exit when drained to half that threshold.
    #[must_use]
    pub fn is_memory_drained(&self) -> bool {
        self.q3_reorder_state.is_memory_drained()
    }

    /// Check if Q5 (processed queue) memory is at the backpressure threshold.
    ///
    /// When true, the Process step should pause to let downstream steps
    /// (Serialize, Compress, Write) drain the queue. This prevents unbounded
    /// memory growth when processing is faster than serialization.
    #[must_use]
    pub fn is_q5_memory_high(&self) -> bool {
        self.output.is_processed_memory_high()
    }

    /// Check if pipeline is in drain mode (bypasses memory backpressure).
    ///
    /// When draining, memory-based backpressure is bypassed to prevent deadlock
    /// during pipeline completion. Queue-full checks still apply.
    #[must_use]
    pub fn is_draining(&self) -> bool {
        self.output.is_draining()
    }

    /// Get optional reference to pipeline statistics.
    #[must_use]
    pub fn stats(&self) -> Option<&PipelineStats> {
        self.output.stats.as_ref()
    }

    /// Get optional reference to progress tracker.
    #[must_use]
    pub fn progress(&self) -> &ProgressTracker {
        &self.output.progress
    }

    /// Get items written count.
    #[must_use]
    pub fn items_written(&self) -> u64 {
        self.output.items_written.load(Ordering::Relaxed)
    }

    /// Set the draining flag.
    pub fn set_draining(&self, value: bool) {
        self.output.set_draining(value);
    }

    /// Flush the output writer and finalize.
    pub fn flush_output(&self) -> io::Result<()> {
        if let Some(mut writer) = self.output.output.lock().take() {
            writer.flush()?;
        }
        Ok(())
    }
}

// ============================================================================
// PipelineLifecycle Trait Implementation
// ============================================================================

impl<G: Send + 'static, P: Send + MemoryEstimate + 'static> PipelineLifecycle
    for BamPipelineState<G, P>
{
    fn is_complete(&self) -> bool {
        BamPipelineState::is_complete(self)
    }

    fn has_error(&self) -> bool {
        BamPipelineState::has_error(self)
    }

    fn take_error(&self) -> Option<io::Error> {
        BamPipelineState::take_error(self)
    }

    fn set_error(&self, error: io::Error) {
        BamPipelineState::set_error(self, error);
    }

    fn is_draining(&self) -> bool {
        BamPipelineState::is_draining(self)
    }

    fn set_draining(&self, value: bool) {
        BamPipelineState::set_draining(self, value);
    }

    fn stats(&self) -> Option<&PipelineStats> {
        BamPipelineState::stats(self)
    }

    fn progress(&self) -> &ProgressTracker {
        BamPipelineState::progress(self)
    }

    fn items_written(&self) -> u64 {
        BamPipelineState::items_written(self)
    }

    fn flush_output(&self) -> io::Result<()> {
        BamPipelineState::flush_output(self)
    }
}

// ============================================================================
// MonitorableState Trait Implementation
// ============================================================================

impl<G: Send + 'static, P: Send + MemoryEstimate + 'static> MonitorableState
    for BamPipelineState<G, P>
{
    fn deadlock_state(&self) -> &DeadlockState {
        &self.deadlock_state
    }

    fn build_queue_snapshot(&self) -> QueueSnapshot {
        // Collect reorder buffer memory (requires locks)
        let q2_reorder_mem = {
            let reorder = self.q2_reorder.lock();
            reorder.total_heap_size() as u64
        };
        let q3_reorder_mem = {
            let reorder = self.q3_reorder.lock();
            reorder.total_heap_size() as u64
        };

        QueueSnapshot {
            q1_len: self.q1_raw_blocks.len(),
            q2_len: self.q2_decompressed.len(),
            q2b_len: self.q2b_boundaries.len(),
            q3_len: self.q3_decoded.len(),
            q4_len: self.output.groups.len(),
            q5_len: self.output.processed.len(),
            q6_len: self.output.serialized.len(),
            q7_len: self.output.compressed.len(),
            q2_reorder_mem,
            q3_reorder_mem,
            memory_limit: self.deadlock_state.get_memory_limit(),
            read_done: self.read_done.load(Ordering::Relaxed),
            group_done: self.group_done.load(Ordering::Relaxed),
            draining: self.output.draining.load(Ordering::Relaxed),
        }
    }
}

impl<G: Send + 'static, P: Send + MemoryEstimate + 'static> OutputPipelineState
    for BamPipelineState<G, P>
{
    type Processed = P;

    fn has_error(&self) -> bool {
        self.output.has_error()
    }

    fn set_error(&self, error: io::Error) {
        self.output.set_error(error);
    }

    fn q5_pop(&self) -> Option<(u64, SerializedBatch)> {
        self.output.serialized.pop()
    }

    fn q5_push(&self, item: (u64, SerializedBatch)) -> Result<(), (u64, SerializedBatch)> {
        self.output.serialized.push(item)
    }

    fn q5_is_full(&self) -> bool {
        self.output.serialized.is_full()
    }

    fn q5_track_pop(&self, heap_size: u64) {
        self.output.serialized_heap_bytes.fetch_sub(heap_size, Ordering::AcqRel);
    }

    fn q6_pop(&self) -> Option<(u64, CompressedBlockBatch)> {
        self.output.compressed.pop()
    }

    fn q6_push(
        &self,
        item: (u64, CompressedBlockBatch),
    ) -> Result<(), (u64, CompressedBlockBatch)> {
        self.output.compressed.push(item)
    }

    fn q6_is_full(&self) -> bool {
        self.output.compressed.is_full()
    }

    fn q6_reorder_insert(&self, serial: u64, batch: CompressedBlockBatch) {
        self.output.write_reorder.lock().insert(serial, batch);
    }

    fn q6_reorder_try_pop_next(&self) -> Option<CompressedBlockBatch> {
        self.output.write_reorder.lock().try_pop_next()
    }

    fn output_try_lock(
        &self,
    ) -> Option<parking_lot::MutexGuard<'_, Option<Box<dyn Write + Send>>>> {
        self.output.output.try_lock()
    }

    fn increment_written(&self) -> u64 {
        self.output.items_written.fetch_add(1, Ordering::Release)
    }

    fn record_compressed_bytes_out(&self, bytes: u64) {
        if let Some(ref stats) = self.output.stats {
            stats.compressed_bytes_out.fetch_add(bytes, Ordering::Relaxed);
        }
    }

    fn record_q6_pop_progress(&self) {
        self.deadlock_state.record_q6_pop();
    }

    fn record_q7_push_progress(&self) {
        self.deadlock_state.record_q7_push();
    }

    fn stats(&self) -> Option<&PipelineStats> {
        self.output.stats.as_ref()
    }
}

// ============================================================================
// New Shared Traits (Phase 2 - Pipeline Consolidation)
// ============================================================================

impl<G: Send + 'static, P: Send + MemoryEstimate + 'static> ProcessPipelineState<G, P>
    for BamPipelineState<G, P>
{
    fn process_input_pop(&self) -> Option<(u64, Vec<G>)> {
        let result = self.output.groups.pop();
        if result.is_some() {
            self.deadlock_state.record_q4_pop();
        }
        result
    }

    fn process_output_is_full(&self) -> bool {
        self.output.processed.is_full()
    }

    fn process_output_push(&self, item: (u64, Vec<P>)) -> Result<(), (u64, Vec<P>)> {
        // Calculate heap size before push for memory tracking
        let heap_size: usize = item.1.iter().map(MemoryEstimate::estimate_heap_size).sum();
        let result = self.output.processed.push(item);
        if result.is_ok() {
            self.output.processed_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
            self.deadlock_state.record_q5_push();
        }
        result
    }

    fn has_error(&self) -> bool {
        self.output.has_error()
    }

    fn set_error(&self, error: io::Error) {
        self.output.set_error(error);
    }

    fn should_apply_process_backpressure(&self) -> bool {
        self.output.should_apply_process_backpressure()
    }

    fn is_draining(&self) -> bool {
        self.output.is_draining()
    }
}

impl<G: Send + 'static, P: Send + MemoryEstimate + 'static> SerializePipelineState<P>
    for BamPipelineState<G, P>
{
    fn serialize_input_pop(&self) -> Option<(u64, Vec<P>)> {
        let result = self.output.processed.pop();
        if let Some((_, ref batch)) = result {
            // Track memory being removed from processed queue
            let heap_size: usize = batch.iter().map(MemoryEstimate::estimate_heap_size).sum();
            self.output.processed_heap_bytes.fetch_sub(heap_size as u64, Ordering::AcqRel);
            self.deadlock_state.record_q5_pop();
        }
        result
    }

    fn serialize_output_is_full(&self) -> bool {
        self.output.serialized.is_full()
    }

    fn serialize_output_push(
        &self,
        item: (u64, SerializedBatch),
    ) -> Result<(), (u64, SerializedBatch)> {
        let heap_size = item.1.estimate_heap_size();
        let result = self.output.serialized.push(item);
        if result.is_ok() {
            self.output.serialized_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
            self.deadlock_state.record_q6_push();
        }
        result
    }

    fn has_error(&self) -> bool {
        self.output.has_error()
    }

    fn set_error(&self, error: io::Error) {
        self.output.set_error(error);
    }

    fn record_serialized_bytes(&self, bytes: u64) {
        if let Some(ref stats) = self.output.stats {
            stats.serialized_bytes.fetch_add(bytes, Ordering::Relaxed);
        }
    }
}

impl<G: Send + 'static, P: Send + MemoryEstimate + 'static> WritePipelineState
    for BamPipelineState<G, P>
{
    fn write_input_queue(&self) -> &ArrayQueue<(u64, CompressedBlockBatch)> {
        &self.output.compressed
    }

    fn write_reorder_buffer(&self) -> &Mutex<ReorderBuffer<CompressedBlockBatch>> {
        &self.output.write_reorder
    }

    fn write_output(&self) -> &Mutex<Option<Box<dyn Write + Send>>> {
        &self.output.output
    }

    fn has_error(&self) -> bool {
        self.output.has_error()
    }

    fn set_error(&self, error: io::Error) {
        self.output.set_error(error);
    }

    fn record_written(&self, count: u64) {
        self.output.items_written.fetch_add(count, Ordering::Release);
    }

    fn stats(&self) -> Option<&PipelineStats> {
        self.output.stats.as_ref()
    }
}

/// Snapshot of queue depths for priority scheduling.
#[derive(Debug, Clone, Copy)]
pub struct QueueDepths {
    pub q1: usize,
    pub q2: usize,
    pub q2b: usize,
    pub q3: usize,
    pub q4: usize,
    pub q5: usize,
    pub q6: usize,
    pub q7: usize,
}

impl QueueDepths {
    /// Check if a step has input available in its input queue.
    /// Returns true if the step might have work, false if input queue is definitely empty.
    #[inline]
    #[must_use]
    pub fn has_input_for_step(&self, step: PipelineStep) -> bool {
        match step {
            PipelineStep::Read => true, // Read always can try (reads from file, not queue)
            PipelineStep::Decompress => self.q1 > 0,
            PipelineStep::FindBoundaries => self.q2 > 0,
            PipelineStep::Decode => self.q2b > 0,
            PipelineStep::Group => self.q3 > 0,
            PipelineStep::Process => self.q4 > 0,
            PipelineStep::Serialize => self.q5 > 0,
            PipelineStep::Compress => self.q6 > 0,
            PipelineStep::Write => self.q7 > 0,
        }
    }
}

// ============================================================================
// Grouper Trait (for Step 5: Group)
// ============================================================================

/// Trait for command-specific record grouping logic.
///
/// The Grouper receives already-decoded BAM records and groups them according
/// to command-specific rules. Since records are pre-decoded, grouping is very
/// fast - just comparing record names, positions, or tags.
///
/// Different commands use different grouping strategies:
/// - `group`: Groups by genomic position
/// - `simplex/duplex/codec`: Groups by MI tag
/// - `filter/clip/correct`: No grouping (each record is its own "group")
pub trait Grouper: Send {
    /// The type of group produced by this grouper.
    type Group: Send;

    /// Add decoded records to the grouper.
    ///
    /// Records are guaranteed to be in order (from template-coordinate sorted BAM).
    /// The grouper maintains partial groups waiting for more records.
    ///
    /// Each `DecodedRecord` contains the record plus a pre-computed `GroupKey`
    /// for fast comparison (position, name hash, library, etc.).
    ///
    /// Returns completed groups (may be empty if more records are needed).
    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>>;

    /// Signal that no more input will arrive (EOF).
    ///
    /// Returns any remaining partial group.
    fn finish(&mut self) -> io::Result<Option<Self::Group>>;

    /// Returns true if the grouper has a partial group.
    fn has_pending(&self) -> bool;
}

/// State for the exclusive `Group` step.
///
/// This is held under a mutex and accessed by whichever thread
/// is currently executing the `Group` step.
pub struct GroupState<G: Send> {
    /// The grouper instance that performs grouping of decoded records.
    pub grouper: Box<dyn Grouper<Group = G> + Send>,
    /// Flag indicating EOF has been signaled to the grouper.
    finished: bool,
    /// Groups waiting to be pushed to Q4 (backpressure buffer).
    pending_groups: VecDeque<G>,
    /// Accumulated weight (total templates) of pending groups.
    /// Used for template-based batching.
    pending_weight: usize,
}

impl<G: Send> GroupState<G> {
    /// Create a new state with the given grouper.
    #[must_use]
    pub fn new(grouper: Box<dyn Grouper<Group = G> + Send>) -> Self {
        Self { grouper, finished: false, pending_groups: VecDeque::new(), pending_weight: 0 }
    }

    /// Check if there are pending groups waiting to be pushed to Q4.
    #[must_use]
    pub fn has_pending_output(&self) -> bool {
        !self.pending_groups.is_empty()
    }

    /// Process decoded records and return completed groups.
    pub fn process(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<G>> {
        self.grouper.add_records(records)
    }

    /// Signal EOF and get any remaining group.
    pub fn finish(&mut self) -> io::Result<Option<G>> {
        if self.finished {
            return Ok(None);
        }
        self.finished = true;
        self.grouper.finish()
    }

    /// Check if EOF has been signaled.
    #[must_use]
    pub fn is_finished(&self) -> bool {
        self.finished
    }

    /// Check if grouper has pending data.
    #[must_use]
    pub fn has_pending(&self) -> bool {
        self.grouper.has_pending()
    }
}

// ============================================================================
// Step Functions Trait (for BAM pipeline)
// ============================================================================

/// Functions provided by the command for each pipeline step.
///
/// Generic parameters:
/// - `G`: Group type (output of `DeserializeGroup`, input to Process)
/// - `P`: Processed type (output of Process, input to Serialize)
#[allow(clippy::type_complexity)]
pub struct PipelineFunctions<G: Send, P: Send> {
    /// Step 5: Process a group. Called in parallel by multiple threads.
    /// Returns `io::Result` for proper error propagation.
    pub process_fn: Box<dyn Fn(G) -> io::Result<P> + Send + Sync>,

    /// Step 6: Serialize processed output to a provided buffer.
    /// Appends serialized BAM bytes to the buffer and returns the record count.
    /// This enables buffer reuse in single-threaded mode to avoid allocations.
    pub serialize_fn: Box<dyn Fn(P, &mut Vec<u8>) -> io::Result<u64> + Send + Sync>,
}

impl<G: Send, P: Send> PipelineFunctions<G, P> {
    /// Create new step functions.
    pub fn new<ProcessFn, SerializeFn>(process_fn: ProcessFn, serialize_fn: SerializeFn) -> Self
    where
        ProcessFn: Fn(G) -> io::Result<P> + Send + Sync + 'static,
        SerializeFn: Fn(P, &mut Vec<u8>) -> io::Result<u64> + Send + Sync + 'static,
    {
        Self { process_fn: Box::new(process_fn), serialize_fn: Box::new(serialize_fn) }
    }
}

// ============================================================================
// Per-Thread Worker State
// ============================================================================

/// Default buffer capacities based on `SingleThreadedBuffers` patterns.
const DECOMPRESSION_BUFFER_CAPACITY: usize = 256 * 1024; // 256KB (4 blocks × 64KB)
const SERIALIZATION_BUFFER_CAPACITY: usize = 64 * 1024; // 64KB (typical group size)

/// Per-thread state for parallel steps.
///
/// Each worker thread has its own compressor, decompressor, scheduler,
/// reusable buffers, and held items for non-blocking pipeline execution.
///
/// # Held Items for Deadlock Prevention
///
/// Each held_* field stores an item that couldn't be pushed to the next queue.
/// Instead of blocking forever (which causes deadlock when all threads block),
/// workers hold the item and return immediately, allowing them to try other
/// steps (especially Write to drain the pipeline).
pub struct WorkerState<P: Send> {
    /// Core worker state (compressor, scheduler, serialization buffer, backoff).
    pub core: WorkerCoreState,
    /// Decompressor for step 2.
    pub decompressor: libdeflater::Decompressor,
    /// Reusable buffer for decompression (Step 2).
    /// Swapped out each batch to avoid per-batch allocation.
    pub decompression_buffer: Vec<u8>,

    // ==================== Held Items for Non-Blocking Steps ====================
    /// Held raw blocks from Read step (couldn't push to `q1_raw_blocks`).
    pub held_raw: Option<(u64, RawBlockBatch)>,
    /// Held decompressed data from Decompress step (couldn't push to `q2_decompressed`).
    /// Includes `heap_size` for memory tracking.
    pub held_decompressed: Option<(u64, DecompressedBatch, usize)>,
    /// Held boundaries from `FindBoundaries` step (couldn't push to `q2b_boundaries`).
    pub held_boundaries: Option<(u64, BoundaryBatch)>,
    /// Held decoded records from Decode step (couldn't push to `q3_decoded`).
    /// Includes `heap_size` for memory tracking - held items have their memory
    /// released from tracking and must re-reserve when retrying.
    pub held_decoded: Option<(u64, Vec<DecodedRecord>, usize)>,
    /// Held processed results from Process step (couldn't push to `q5_processed`).
    /// Includes `heap_size` for memory tracking.
    pub held_processed: Option<(u64, Vec<P>, usize)>,
    /// Held serialized batch from Serialize step (couldn't push to `q6_serialized`).
    /// Includes `heap_size` for memory tracking.
    pub held_serialized: Option<(u64, SerializedBatch, usize)>,
    /// Held compressed batch from Compress step (couldn't push to `q7_compressed`).
    /// Includes `heap_size` for memory tracking.
    pub held_compressed: Option<(u64, CompressedBlockBatch, usize)>,
}

impl<P: Send> WorkerState<P> {
    /// Create new worker state with the given compression level, thread info, and scheduler strategy.
    #[must_use]
    pub fn new(
        compression_level: u32,
        thread_id: usize,
        num_threads: usize,
        scheduler_strategy: SchedulerStrategy,
    ) -> Self {
        Self {
            core: WorkerCoreState::new(
                compression_level,
                thread_id,
                num_threads,
                scheduler_strategy,
            ),
            decompressor: libdeflater::Decompressor::new(),
            decompression_buffer: Vec::with_capacity(DECOMPRESSION_BUFFER_CAPACITY),
            // Initialize all held items to None
            held_raw: None,
            held_decompressed: None,
            held_boundaries: None,
            held_decoded: None,
            held_processed: None,
            held_serialized: None,
            held_compressed: None,
        }
    }

    /// Returns true if any held item fields are Some.
    ///
    /// Used to check if a worker still has pending work before completion.
    #[inline]
    #[must_use]
    pub fn has_any_held_items(&self) -> bool {
        self.held_raw.is_some()
            || self.held_decompressed.is_some()
            || self.held_boundaries.is_some()
            || self.held_decoded.is_some()
            || self.held_processed.is_some()
            || self.held_serialized.is_some()
            || self.held_compressed.is_some()
    }

    /// Clear all held items (for cleanup/error handling).
    pub fn clear_held_items(&mut self) {
        self.held_raw = None;
        self.held_decompressed = None;
        self.held_boundaries = None;
        self.held_decoded = None;
        self.held_processed = None;
        self.held_serialized = None;
        self.held_compressed = None;
    }
}

impl<P: Send> HasCompressor for WorkerState<P> {
    fn compressor_mut(&mut self) -> &mut InlineBgzfCompressor {
        &mut self.core.compressor
    }
}

impl<P: Send> HasHeldCompressed for WorkerState<P> {
    fn held_compressed_mut(&mut self) -> &mut Option<(u64, CompressedBlockBatch, usize)> {
        &mut self.held_compressed
    }
}

impl<P: Send> HasHeldBoundaries<BoundaryBatch> for WorkerState<P> {
    fn held_boundaries_mut(&mut self) -> &mut Option<(u64, BoundaryBatch)> {
        &mut self.held_boundaries
    }
}

impl<P: Send> HasHeldProcessed<P> for WorkerState<P> {
    fn held_processed_mut(&mut self) -> &mut Option<(u64, Vec<P>, usize)> {
        &mut self.held_processed
    }
}

impl<P: Send> HasHeldSerialized for WorkerState<P> {
    fn held_serialized_mut(&mut self) -> &mut Option<(u64, SerializedBatch, usize)> {
        &mut self.held_serialized
    }

    fn serialization_buffer_mut(&mut self) -> &mut Vec<u8> {
        &mut self.core.serialization_buffer
    }

    fn serialization_buffer_capacity(&self) -> usize {
        SERIALIZATION_BUFFER_CAPACITY // 64KB for BAM
    }
}

impl<P: Send> WorkerStateCommon for WorkerState<P> {
    fn has_any_held_items(&self) -> bool {
        WorkerState::has_any_held_items(self)
    }

    fn clear_held_items(&mut self) {
        WorkerState::clear_held_items(self);
    }
}

impl<P: Send> HasWorkerCore for WorkerState<P> {
    fn core(&self) -> &WorkerCoreState {
        &self.core
    }

    fn core_mut(&mut self) -> &mut WorkerCoreState {
        &mut self.core
    }
}

// ============================================================================
// Step Execution Functions
// ============================================================================

/// Try to execute Step 1: Read raw BGZF blocks.
///
/// Returns true if work was done.
///
/// # Non-Blocking Design
///
/// Uses the held-item pattern to prevent deadlock. If the output queue is full,
/// the batch is stored in `worker.held_raw` and the function returns immediately.
/// This allows the worker to try other steps (especially Write) to drain the pipeline.
fn try_step_read<G: Send, P: Send + MemoryEstimate>(
    state: &BamPipelineState<G, P>,
    worker: &mut WorkerState<P>,
) -> bool {
    // =========================================================================
    // Priority 1: Try to advance any held raw batch first
    // =========================================================================
    if let Some((serial, held)) = worker.held_raw.take() {
        match state.q1_raw_blocks.push((serial, held)) {
            Ok(()) => {
                // Successfully advanced held item, continue to read more
                state.deadlock_state.record_q1_push();
            }
            Err((serial, held)) => {
                // Still can't push - put it back and signal output full
                worker.held_raw = Some((serial, held));
                return false;
            }
        }
    }

    // =========================================================================
    // Priority 2: Skip if reading is done
    // =========================================================================
    if state.read_done.load(Ordering::Relaxed) {
        return false;
    }

    // =========================================================================
    // Priority 3: Check if output queue has space (soft check)
    // =========================================================================
    if state.q1_raw_blocks.len() >= state.config.queue_capacity {
        return false;
    }

    // =========================================================================
    // Priority 4: Try to acquire exclusive access to input file
    // =========================================================================
    let Some(mut guard) = state.input_file.try_lock() else {
        // Record contention for diagnostics
        if let Some(stats) = state.stats() {
            stats.record_contention(PipelineStep::Read);
        }
        return false; // Another thread is reading
    };

    let Some(ref mut reader) = *guard else {
        return false; // File already closed
    };

    // =========================================================================
    // Priority 5: Read a batch of raw BGZF blocks
    // =========================================================================
    // Read FIRST, then assign serial - ensures we don't waste serial numbers on EOF
    match read_raw_blocks(reader.as_mut(), state.config.blocks_per_read_batch) {
        Ok(blocks) if blocks.is_empty() => {
            // EOF - no data read, don't increment serial
            state.read_done.store(true, Ordering::SeqCst);
            false
        }
        Ok(blocks) => {
            // Data was read, now assign serial number
            let serial = state.next_read_serial.fetch_add(1, Ordering::SeqCst);
            let batch = RawBlockBatch { blocks };

            // Record bytes read for throughput metrics
            if let Some(stats) = state.stats() {
                stats.bytes_read.fetch_add(batch.total_compressed_size() as u64, Ordering::Relaxed);
            }

            // =========================================================================
            // Priority 6: Try to push result (non-blocking)
            // =========================================================================
            match state.q1_raw_blocks.push((serial, batch)) {
                Ok(()) => {
                    state.deadlock_state.record_q1_push();
                    true
                }
                Err((serial, batch)) => {
                    // Output full - hold the result for next attempt
                    worker.held_raw = Some((serial, batch));
                    false
                }
            }
        }
        Err(e) => {
            state.set_error(e);
            false
        }
    }
}

/// Try to execute Step 2: Decompress blocks.
///
/// This step is parallel - multiple threads can decompress concurrently.
///
/// # Non-Blocking Design
///
/// Uses the held-item pattern to prevent deadlock. If the output queue is full,
/// the batch is stored in `worker.held_decompressed` and the function returns immediately.
///
/// # Memory-Based Backpressure
///
/// Before pushing to `q2_decompressed`, checks if the Q2 reorder buffer would accept
/// this serial number. This prevents unbounded memory growth in the reorder buffer
/// while avoiding deadlock by always accepting the serial that `FindBoundaries` needs.
fn try_step_decompress<G: Send, P: Send + MemoryEstimate>(
    state: &BamPipelineState<G, P>,
    worker: &mut WorkerState<P>,
) -> bool {
    // =========================================================================
    // Priority 1: Try to advance any held decompressed batch first
    // =========================================================================
    // Held items had their heap_size released when held. Re-reserve before checking.
    if let Some((serial, held, heap_size)) = worker.held_decompressed.take() {
        // Re-reserve memory for this held batch
        state.q2_reorder_state.add_heap_bytes(heap_size as u64);

        // Check memory-based backpressure before trying to push
        if !state.can_decompress_proceed(serial) {
            // Reorder buffer is over limit - release reservation and hold again
            state.q2_reorder_state.sub_heap_bytes(heap_size as u64);
            worker.held_decompressed = Some((serial, held, heap_size));
            return false;
        }
        match state.q2_decompressed.push((serial, held)) {
            Ok(()) => {
                // Successfully pushed - keep reservation (released by FindBoundaries when popping)
                state.batches_decompressed.fetch_add(1, Ordering::Release);
                state.deadlock_state.record_q2_push();
            }
            Err((serial, held)) => {
                // Can't push - release reservation and hold
                state.q2_reorder_state.sub_heap_bytes(heap_size as u64);
                worker.held_decompressed = Some((serial, held, heap_size));
                return false;
            }
        }
    }

    // =========================================================================
    // Priority 2: Check if output queue has space (soft check)
    // =========================================================================
    if state.q2_decompressed.is_full() {
        return false;
    }

    // =========================================================================
    // Priority 3: Pop input and process
    // =========================================================================
    let Some((serial, raw_batch)) = state.q1_raw_blocks.pop() else {
        if let Some(stats) = state.stats() {
            stats.record_queue_empty(1);
        }
        return false;
    };
    state.deadlock_state.record_q1_pop();

    // Prepare worker's buffer: clear and reserve capacity
    worker.decompression_buffer.clear();
    let expected_size = raw_batch.total_uncompressed_size();
    worker.decompression_buffer.reserve(expected_size);

    // Decompress directly into worker's buffer (no intermediate allocations)
    for block in &raw_batch.blocks {
        if let Err(e) =
            decompress_block_into(block, &mut worker.decompressor, &mut worker.decompression_buffer)
        {
            state.set_error(e);
            return false;
        }
    }

    // Swap buffer into batch, replace with fresh pre-allocated buffer
    let decompressed = std::mem::replace(
        &mut worker.decompression_buffer,
        Vec::with_capacity(DECOMPRESSION_BUFFER_CAPACITY),
    );

    // Record compression metrics
    if let Some(stats) = state.stats() {
        stats
            .compressed_bytes_in
            .fetch_add(raw_batch.total_compressed_size() as u64, Ordering::Relaxed);
        stats.decompressed_bytes.fetch_add(decompressed.len() as u64, Ordering::Relaxed);
    }

    // =========================================================================
    // Priority 4: Calculate and reserve memory BEFORE checking backpressure
    // =========================================================================
    let batch = DecompressedBatch { data: decompressed };
    let heap_size = batch.estimate_heap_size();
    state.q2_reorder_state.add_heap_bytes(heap_size as u64);

    // =========================================================================
    // Priority 5: Check memory backpressure and try to push result
    // =========================================================================
    if !state.can_decompress_proceed(serial) {
        // Over limit - release reservation and hold
        state.q2_reorder_state.sub_heap_bytes(heap_size as u64);
        worker.held_decompressed = Some((serial, batch, heap_size));
        return false;
    }

    match state.q2_decompressed.push((serial, batch)) {
        Ok(()) => {
            // Successfully pushed - keep reservation (released by FindBoundaries when popping)
            state.batches_decompressed.fetch_add(1, Ordering::Release);
            state.deadlock_state.record_q2_push();
            true
        }
        Err((serial, batch)) => {
            // Output full - release reservation and hold
            state.q2_reorder_state.sub_heap_bytes(heap_size as u64);
            worker.held_decompressed = Some((serial, batch, heap_size));
            false
        }
    }
}

/// Try to execute Step 3: Find record boundaries in decompressed data.
///
/// SYNC WITH: fastq.rs `fastq_try_step_find_boundaries()`
/// Both implementations use the "held boundaries" pattern for parallelism.
/// See base.rs `HasHeldBoundaries` trait for pattern documentation.
///
/// This step is exclusive but FAST (~0.1μs per block) - only scans 4-byte integers.
/// Processes multiple batches per lock acquisition to reduce contention.
/// Result of attempting an exclusive step.
/// Returns (`did_work`, `was_contention`).
///
/// # Non-Blocking Design
///
/// Uses the held-item pattern to prevent deadlock. If the output queue is full,
/// the batch is stored in `worker.held_boundaries` and the function returns immediately.
fn try_step_find_boundaries<G: Send, P: Send + MemoryEstimate>(
    state: &BamPipelineState<G, P>,
    worker: &mut WorkerState<P>,
) -> (bool, bool) {
    // =========================================================================
    // Priority 1: Try to advance any held boundary batch first
    // =========================================================================
    let mut did_work = false;
    if let Some((serial, held)) = worker.held_boundaries.take() {
        match state.q2b_boundaries.push((serial, held)) {
            Ok(()) => {
                // Successfully advanced held item, increment completion counter
                state.batches_boundary_found.fetch_add(1, Ordering::Release);
                state.deadlock_state.record_q2b_push();
                did_work = true;
            }
            Err((serial, held)) => {
                // Still can't push - put it back and signal backpressure
                worker.held_boundaries = Some((serial, held));
                return (false, false); // Backpressure, not contention
            }
        }
    }

    // =========================================================================
    // Priority 2: Check if output queue has space
    // =========================================================================
    if state.q2b_boundaries.is_full() {
        return (false, false); // Backpressure, not contention
    }

    // =========================================================================
    // Priority 3: Try to acquire exclusive access to boundary state
    // =========================================================================
    let Some(mut boundary_guard) = state.boundary_state.try_lock() else {
        if let Some(stats) = state.stats() {
            stats.record_contention(PipelineStep::FindBoundaries);
        }
        return (did_work, true); // Contention (but may have advanced held item)
    };

    // Process multiple batches per lock acquisition to reduce contention
    // Max batches per acquisition (tune for balance between throughput and latency)
    const MAX_BATCHES_PER_LOCK: usize = 8;

    for _ in 0..MAX_BATCHES_PER_LOCK {
        // Check if output queue still has space
        if state.q2b_boundaries.is_full() {
            break;
        }

        // Drain Q2 into reorder buffer AND get next in-order batch
        let batch_with_size = {
            let mut reorder = state.q2_reorder.lock();

            // Insert all pending decompressed batches into reorder buffer.
            // Memory was already reserved by Decompress - just insert for ordering.
            while let Some((serial, batch)) = state.q2_decompressed.pop() {
                state.deadlock_state.record_q2_pop();
                let heap_size = batch.estimate_heap_size();
                reorder.insert_with_size(serial, batch, heap_size);
            }

            // Try to pop the next in-order batch
            let result = reorder.try_pop_next_with_size();

            // Update next_seq atomic for Decompress's backpressure check
            state.q2_reorder_state.update_next_seq(reorder.next_seq());

            result
        };

        // Release memory from atomic tracker when popping from reorder buffer
        let Some((batch, heap_size)) = batch_with_size else {
            if !did_work {
                if let Some(stats) = state.stats() {
                    stats.record_queue_empty(2);
                }
            }
            break; // No more data available
        };
        state.q2_reorder_state.sub_heap_bytes(heap_size as u64);

        // Track that we've processed a batch from q2 (for completion tracking)
        state.batches_boundary_processed.fetch_add(1, Ordering::Release);

        // Find boundaries in the decompressed data
        match boundary_guard.find_boundaries(&batch.data) {
            Ok(boundary_batch) => {
                // Only push if there are records
                if boundary_batch.offsets.len() > 1 {
                    // Record batch size for statistics
                    let num_records = boundary_batch.offsets.len() - 1;
                    if let Some(stats) = state.stats() {
                        stats.record_batch_size(num_records);
                    }

                    let serial = state.next_boundary_serial.fetch_add(1, Ordering::SeqCst);
                    // Try non-blocking push
                    match state.q2b_boundaries.push((serial, boundary_batch)) {
                        Ok(()) => {
                            // Successfully pushed, increment completion counter
                            state.batches_boundary_found.fetch_add(1, Ordering::Release);
                            state.deadlock_state.record_q2b_push();
                        }
                        Err((serial, boundary_batch)) => {
                            // Output full - hold the result and stop processing
                            worker.held_boundaries = Some((serial, boundary_batch));
                            return (true, false); // Did work (processed data), will retry push later
                        }
                    }
                }
                did_work = true;
            }
            Err(e) => {
                state.set_error(e);
                return (false, false);
            }
        }
    }

    if did_work {
        return (true, false); // Success, no contention (we held the lock)
    }

    // No batches processed - check if we should finish
    // Completion check: Only finish when THIS step has processed all input batches.
    // We check batches_boundary_processed == total_read directly, which implies that
    // Decompress has also finished (since we can't process more than was decompressed).
    // This prevents a race where FindBoundaries sets boundary_done while data is still
    // in q2_decompressed waiting to be processed.
    let read_done = state.read_done.load(Ordering::Acquire);
    let total_read = state.next_read_serial.load(Ordering::Acquire);
    let batches_boundary_processed = state.batches_boundary_processed.load(Ordering::Acquire);

    if read_done
        && batches_boundary_processed == total_read
        && !state.boundary_done.load(Ordering::Acquire)
    {
        // All input processed - finish and emit any remaining boundaries
        match boundary_guard.finish() {
            Ok(Some(final_batch)) => {
                if final_batch.offsets.len() > 1 {
                    let serial = state.next_boundary_serial.fetch_add(1, Ordering::SeqCst);
                    // Try non-blocking push for final batch
                    match state.q2b_boundaries.push((serial, final_batch)) {
                        Ok(()) => {
                            // Successfully pushed, increment completion counter
                            state.batches_boundary_found.fetch_add(1, Ordering::Release);
                            state.deadlock_state.record_q2b_push();
                        }
                        Err((serial, final_batch)) => {
                            // Hold final batch for next attempt
                            worker.held_boundaries = Some((serial, final_batch));
                            return (true, false);
                        }
                    }
                }
                state.boundary_done.store(true, Ordering::SeqCst);
                (true, false)
            }
            Ok(None) => {
                state.boundary_done.store(true, Ordering::SeqCst);
                (false, false)
            }
            Err(e) => {
                state.set_error(e);
                (false, false)
            }
        }
    } else {
        (false, false) // No work available, no contention
    }
}

/// Try to execute Step 4: Decode BAM records at known boundaries.
///
/// This step is parallel - multiple threads can decode concurrently.
///
/// # Non-Blocking Design
///
/// Uses the held-item pattern to prevent deadlock. If the output queue is full,
/// the batch is stored in `worker.held_decoded` and the function returns immediately.
///
/// # Memory-Based Backpressure
///
/// Before pushing to `q3_decoded`, checks if the Q3 reorder buffer would accept
/// this serial number. This prevents unbounded memory growth in the reorder buffer
/// while avoiding deadlock by always accepting the serial that Group is waiting for.
fn try_step_decode<G: Send, P: Send + MemoryEstimate>(
    state: &BamPipelineState<G, P>,
    worker: &mut WorkerState<P>,
) -> bool {
    // =========================================================================
    // Priority 1: Try to advance any held decoded batch first
    // =========================================================================
    // Held items had their heap_size released when held. Re-reserve before checking.
    if let Some((serial, held, heap_size)) = worker.held_decoded.take() {
        // Re-reserve memory for this held batch
        state.q3_reorder_state.add_heap_bytes(heap_size as u64);

        // Check memory-based backpressure before trying to push
        if !state.can_decode_proceed(serial) {
            // Reorder buffer is over limit - release reservation and hold again
            state.q3_reorder_state.sub_heap_bytes(heap_size as u64);
            worker.held_decoded = Some((serial, held, heap_size));
            return false;
        }
        match state.q3_decoded.push((serial, held)) {
            Ok(()) => {
                // Successfully pushed - keep reservation (released by Group when popping)
                state.batches_decoded.fetch_add(1, Ordering::Release);
                state.deadlock_state.record_q3_push();
            }
            Err((serial, held)) => {
                // Can't push - release reservation and hold
                state.q3_reorder_state.sub_heap_bytes(heap_size as u64);
                worker.held_decoded = Some((serial, held, heap_size));
                return false;
            }
        }
    }

    // =========================================================================
    // Priority 2: Check if output queue has space
    // =========================================================================
    if state.q3_decoded.is_full() {
        return false;
    }

    // =========================================================================
    // Priority 3: Pop input
    // =========================================================================
    let Some((serial, boundary_batch)) = state.q2b_boundaries.pop() else {
        if let Some(stats) = state.stats() {
            stats.record_queue_empty(25); // Q2b (boundaries queue)
        }
        // Check if boundary finding is done and queue is empty
        if state.boundary_done.load(Ordering::SeqCst) && state.q2b_boundaries.is_empty() {
            state.decode_done.store(true, Ordering::SeqCst);
        } else if let Some(stats) = state.stats() {
            // Q2b is extension of Q2
            stats.record_queue_empty(2);
        }
        return false;
    };
    state.deadlock_state.record_q2b_pop();

    // =========================================================================
    // Priority 4: Decode records and compute GroupKey
    // =========================================================================
    match decode_records(&boundary_batch, &state.group_key_config) {
        Ok(records) => {
            // Record decoded count for throughput metrics
            if let Some(stats) = state.stats() {
                stats.records_decoded.fetch_add(records.len() as u64, Ordering::Relaxed);
            }

            // Calculate and reserve memory BEFORE checking backpressure
            let heap_size = records.estimate_heap_size();
            state.q3_reorder_state.add_heap_bytes(heap_size as u64);

            // =========================================================================
            // Priority 5: Check memory backpressure and try to push result
            // =========================================================================
            if !state.can_decode_proceed(serial) {
                // Over limit - release reservation and hold
                state.q3_reorder_state.sub_heap_bytes(heap_size as u64);
                worker.held_decoded = Some((serial, records, heap_size));
                return false;
            }

            match state.q3_decoded.push((serial, records)) {
                Ok(()) => {
                    // Successfully pushed - keep reservation (released by Group when popping)
                    state.batches_decoded.fetch_add(1, Ordering::Release);
                    state.deadlock_state.record_q3_push();
                    true
                }
                Err((serial, records)) => {
                    // Output full - release reservation and hold
                    state.q3_reorder_state.sub_heap_bytes(heap_size as u64);
                    worker.held_decoded = Some((serial, records, heap_size));
                    false
                }
            }
        }
        Err(e) => {
            state.set_error(e);
            false
        }
    }
}

/// Try to execute Step 5: Group decoded records.
///
/// This step is exclusive - only one thread at a time.
/// Groups are accumulated and pushed to Q4 as batches for compression efficiency.
///
/// Batching mode depends on configuration:
/// - `target_templates_per_batch > 0`: Weight-based batching using `BatchWeight::batch_weight()`
/// - `target_templates_per_batch == 0`: Count-based batching using `batch_size`
fn try_step_group<G: Send + BatchWeight + MemoryEstimate + 'static, P: Send + MemoryEstimate>(
    state: &BamPipelineState<G, P>,
    group_state: &Mutex<GroupState<G>>,
) -> (bool, bool) {
    // Try to acquire exclusive access to group state
    let Some(mut guard) = group_state.try_lock() else {
        if let Some(stats) = state.stats() {
            stats.record_contention(PipelineStep::Group);
        }
        return (false, true); // Contention!
    };

    let batch_size = state.config.batch_size;
    let target_weight = state.config.target_templates_per_batch;
    let use_weight_batching = target_weight > 0;

    // Helper to check if we should flush the pending batch
    let should_flush = |pending_len: usize, pending_weight: usize| -> bool {
        if use_weight_batching {
            pending_weight >= target_weight
        } else {
            pending_len >= batch_size
        }
    };

    // Helper to push a batch of groups to Q4
    // Returns Ok(()) on success, Err(batch) on failure so caller can restore
    // Note: We don't track Q4 memory here - only Q3 (decoded records) is tracked.
    // Mixing different estimation methods (O(1) for Q3, estimate_heap_size for Q4)
    // on the same tracker causes imbalance and deadlock.
    let push_batch = |groups: Vec<G>, state: &BamPipelineState<G, P>| -> Result<(), Vec<G>> {
        if state.output.groups.is_full() {
            return Err(groups);
        }

        let serial = state.next_group_serial.fetch_add(1, Ordering::SeqCst);
        state
            .output
            .groups
            .push((serial, groups))
            .unwrap_or_else(|_| panic!("groups push failed after is_full check"));
        state.deadlock_state.record_q4_push();
        Ok(())
    };

    // Helper to flush all pending groups and reset weight
    let flush_all = |guard: &mut GroupState<G>, state: &BamPipelineState<G, P>| -> Option<bool> {
        if guard.pending_groups.is_empty() {
            return Some(true);
        }
        if state.output.groups.is_full() {
            return None; // Backpressure
        }
        let batch: Vec<G> = guard.pending_groups.drain(..).collect();
        guard.pending_weight = 0;
        match push_batch(batch, state) {
            Ok(()) => Some(true),
            Err(batch) => {
                // Restore the batch on failure (race condition)
                for group in batch.into_iter().rev() {
                    guard.pending_weight += group.batch_weight();
                    guard.pending_groups.push_front(group);
                }
                None
            }
        }
    };

    // First, try to drain pending_groups to Q4 as batches
    while should_flush(guard.pending_groups.len(), guard.pending_weight) {
        if state.output.groups.is_full() {
            return (false, false); // Queue backpressure, not contention
        }
        // For weight-based batching, flush all pending groups at once
        // For count-based batching, drain exactly batch_size
        if use_weight_batching {
            let batch: Vec<G> = guard.pending_groups.drain(..).collect();
            guard.pending_weight = 0;
            if let Err(batch) = push_batch(batch, state) {
                // Restore the batch on race condition failure
                for group in batch.into_iter().rev() {
                    guard.pending_weight += group.batch_weight();
                    guard.pending_groups.push_front(group);
                }
                return (false, false);
            }
            break; // Only one flush per check for weight-based
        }
        let batch: Vec<G> = guard.pending_groups.drain(..batch_size).collect();
        if let Err(batch) = push_batch(batch, state) {
            // Restore the batch on race condition failure
            for group in batch.into_iter().rev() {
                guard.pending_weight += group.batch_weight();
                guard.pending_groups.push_front(group);
            }
            return (false, false); // Backpressure
        }
    }

    // If finish() was called and all pending groups have been drained
    if guard.is_finished() && !state.group_done.load(Ordering::SeqCst) {
        // CRITICAL: Must retry until flush succeeds to avoid data loss
        let mut retries = 0;
        const MAX_RETRIES: u32 = 10_000;
        loop {
            if flush_all(&mut guard, state).is_some() {
                state.group_done.store(true, Ordering::SeqCst);
                return (true, false); // Success
            }
            retries += 1;
            if retries > MAX_RETRIES {
                state.set_error(std::io::Error::other(format!(
                    "Failed to flush final groups after {} retries",
                    MAX_RETRIES
                )));
                return (false, false);
            }
            std::thread::yield_now();
        }
    }

    // Process multiple record batches per lock acquisition to reduce contention
    const MAX_BATCHES_PER_LOCK: usize = 8;
    const MAX_PENDING_DRAIN: usize = 16;
    let mut did_work = false;

    // Reusable buffer for pre-draining (allocated once per try_step_group call)
    let mut pending: Vec<(u64, Vec<DecodedRecord>, usize)> = Vec::with_capacity(MAX_PENDING_DRAIN);

    for _ in 0..MAX_BATCHES_PER_LOCK {
        // Pre-drain q3_decoded BEFORE taking the reorder lock (lock-free operations)
        // This reduces critical section time by moving ArrayQueue ops outside the lock
        pending.clear();
        while pending.len() < MAX_PENDING_DRAIN {
            if let Some((serial, batch)) = state.q3_decoded.pop() {
                state.deadlock_state.record_q3_pop();
                let heap_size = batch.estimate_heap_size();
                pending.push((serial, batch, heap_size));
            } else {
                break;
            }
        }

        // Now take the reorder lock and insert all pending batches
        let records = {
            let mut reorder = state.q3_reorder.lock();

            // Insert all pre-drained batches into reorder buffer
            for (serial, batch, heap_size) in pending.drain(..) {
                reorder.insert_with_size(serial, batch, heap_size);
            }

            // Try to pop the next in-order batch
            let result = reorder.try_pop_next_with_size();

            // Update can_pop and next_seq atomics for Decode's backpressure check
            state.q3_reorder_state.update_next_seq(reorder.next_seq());
            state.q3_reorder_can_pop.store(reorder.can_pop(), Ordering::Release);

            result
        };

        // Release memory from atomic tracker when popping from reorder buffer
        let Some((records, heap_size)) = records else {
            if !did_work {
                if let Some(stats) = state.stats() {
                    stats.record_queue_empty(3);
                }
            }
            break; // No more data available
        };
        state.q3_reorder_state.sub_heap_bytes(heap_size as u64);

        // Track that we've processed a batch from q3 (for completion tracking)
        state.batches_grouped.fetch_add(1, Ordering::Release);

        // Process the decoded records
        match guard.process(records) {
            Ok(groups) => {
                // Record groups produced for throughput metrics
                if let Some(stats) = state.stats() {
                    stats.groups_produced.fetch_add(groups.len() as u64, Ordering::Relaxed);
                }

                // Add groups and track weight
                for group in groups {
                    guard.pending_weight += group.batch_weight();
                    guard.pending_groups.push_back(group);
                }

                // Push batches when threshold is reached
                while should_flush(guard.pending_groups.len(), guard.pending_weight) {
                    if state.output.groups.is_full() {
                        return (true, false); // Did work, stopping due to backpressure
                    }
                    if use_weight_batching {
                        let batch: Vec<G> = guard.pending_groups.drain(..).collect();
                        guard.pending_weight = 0;
                        if let Err(batch) = push_batch(batch, state) {
                            // Restore the batch on race condition failure
                            for group in batch.into_iter().rev() {
                                guard.pending_weight += group.batch_weight();
                                guard.pending_groups.push_front(group);
                            }
                            return (true, false); // Memory limit race - backpressure
                        }
                        break;
                    }
                    let batch: Vec<G> = guard.pending_groups.drain(..batch_size).collect();
                    if let Err(batch) = push_batch(batch, state) {
                        // Restore the batch on race condition failure
                        for group in batch.into_iter().rev() {
                            guard.pending_weight += group.batch_weight();
                            guard.pending_groups.push_front(group);
                        }
                        return (true, false); // Memory limit race - backpressure
                    }
                }
                did_work = true;
            }
            Err(e) => {
                state.set_error(e);
                return (false, false);
            }
        }
    }

    if did_work {
        return (true, false); // Success
    }

    // No records processed - check if we should finish
    if guard.is_finished() {
        return (false, false); // Already finished
    }

    // Use atomic counters for completion check (not reorder buffer next_seq)
    // This avoids TOCTOU races because counters are incremented atomically after push.
    // CRITICAL: We use batches_grouped (batches processed by Group), NOT batches_decoded
    // (batches pushed to q3). Using batches_decoded caused a race where Group would finish
    // before actually processing all the data in q3.
    let boundary_done = state.boundary_done.load(Ordering::Acquire);
    let total_boundary_batches = state.batches_boundary_found.load(Ordering::Acquire);
    let batches_grouped = state.batches_grouped.load(Ordering::Acquire);

    if boundary_done && batches_grouped == total_boundary_batches {
        // All input processed - finish and emit any remaining group
        match guard.finish() {
            Ok(Some(group)) => {
                guard.pending_weight += group.batch_weight();
                guard.pending_groups.push_back(group);

                // Flush remaining batches
                while should_flush(guard.pending_groups.len(), guard.pending_weight) {
                    if state.output.groups.is_full() {
                        return (true, false); // Did work, backpressure
                    }
                    if use_weight_batching {
                        let batch: Vec<G> = guard.pending_groups.drain(..).collect();
                        guard.pending_weight = 0;
                        if let Err(batch) = push_batch(batch, state) {
                            // Restore the batch on race condition failure
                            for group in batch.into_iter().rev() {
                                guard.pending_weight += group.batch_weight();
                                guard.pending_groups.push_front(group);
                            }
                            return (true, false); // backpressure
                        }
                        break;
                    }
                    let batch: Vec<G> = guard.pending_groups.drain(..batch_size).collect();
                    if let Err(batch) = push_batch(batch, state) {
                        // Restore the batch on race condition failure
                        for group in batch.into_iter().rev() {
                            guard.pending_weight += group.batch_weight();
                            guard.pending_groups.push_front(group);
                        }
                        return (true, false); // backpressure
                    }
                }

                // Flush any remaining
                // CRITICAL: Must retry until flush succeeds to avoid data loss
                {
                    let mut retries = 0;
                    const MAX_RETRIES: u32 = 10_000;
                    loop {
                        if flush_all(&mut guard, state).is_some() {
                            state.group_done.store(true, Ordering::SeqCst);
                            break;
                        }
                        retries += 1;
                        if retries > MAX_RETRIES {
                            state.set_error(std::io::Error::other(format!(
                                "Failed to flush final groups after {} retries",
                                MAX_RETRIES
                            )));
                            return (false, false);
                        }
                        std::thread::yield_now();
                    }
                }
                (true, false) // Success
            }
            Ok(None) => {
                // CRITICAL: Must retry until flush succeeds to avoid data loss
                let mut retries = 0;
                const MAX_RETRIES: u32 = 10_000;
                loop {
                    if flush_all(&mut guard, state).is_some() {
                        state.group_done.store(true, Ordering::SeqCst);
                        break;
                    }
                    retries += 1;
                    if retries > MAX_RETRIES {
                        state.set_error(std::io::Error::other(format!(
                            "Failed to flush final groups after {} retries",
                            MAX_RETRIES
                        )));
                        return (false, false);
                    }
                    std::thread::yield_now();
                }
                (false, false) // Finished, no new work
            }
            Err(e) => {
                state.set_error(e);
                (false, false)
            }
        }
    } else {
        (false, false) // No work available
    }
}

/// Try to execute Step 6: Process groups.
///
/// This step is parallel - multiple threads can process concurrently.
/// Receives a batch of groups from Q4, processes each, and pushes
/// the batch of processed results to Q5.
///
/// # Non-Blocking Design
///
/// Uses the held-item pattern to prevent deadlock. If the output queue is full,
/// the batch is stored in `worker.held_processed` and the function returns immediately.
///
/// # Memory Backpressure
///
/// This function checks both queue capacity AND memory pressure before doing new work.
/// The backpressure check happens AFTER advancing any held item, which matches the
/// baseline behavior and prevents memory spikes at high thread counts.
fn try_step_process<G: Send + MemoryEstimate + 'static, P: Send + MemoryEstimate + 'static>(
    state: &BamPipelineState<G, P>,
    fns: &PipelineFunctions<G, P>,
    worker: &mut WorkerState<P>,
) -> bool {
    // =========================================================================
    // Priority 1: Try to advance any held processed batch first
    // =========================================================================
    if let Some((serial, held, heap_size)) = worker.held_processed.take() {
        match state.output.processed.push((serial, held)) {
            Ok(()) => {
                // Successfully advanced held item - memory tracking handled by trait impl
                state.output.processed_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
                state.deadlock_state.record_q5_push();
            }
            Err((serial, held)) => {
                // Still can't push - put it back and signal output full
                worker.held_processed = Some((serial, held, heap_size));
                return false;
            }
        }
    }

    // =========================================================================
    // Priority 2: Check if output queue has space (count and memory)
    // When draining, bypass memory backpressure to prevent deadlock
    // =========================================================================
    if state.output.processed.is_full() || (!state.is_draining() && state.is_q5_memory_high()) {
        return false;
    }

    // =========================================================================
    // Priority 3: Pop and process batches
    // Always drain multiple batches when work is available for better throughput.
    // Q5 memory backpressure above prevents unbounded growth.
    // =========================================================================
    const MAX_BATCHES: usize = 8;
    let mut did_work = false;

    for _ in 0..MAX_BATCHES {
        // Check output space (count and memory) before each batch
        // When draining, bypass memory backpressure to prevent deadlock
        if state.output.processed.is_full() || (!state.is_draining() && state.is_q5_memory_high()) {
            break;
        }

        let Some((serial, batch)) = state.output.groups.pop() else {
            if let Some(stats) = state.stats() {
                stats.record_queue_empty(4);
            }
            break;
        };
        state.deadlock_state.record_q4_pop();
        // Process each group in the batch
        let mut results: Vec<P> = Vec::with_capacity(batch.len());
        for group in batch {
            match (fns.process_fn)(group) {
                Ok(processed) => results.push(processed),
                Err(e) => {
                    state.set_error(e);
                    return false;
                }
            }
        }

        // Calculate heap size for memory tracking
        let heap_size: usize = results.iter().map(MemoryEstimate::estimate_heap_size).sum();

        // Try to push result
        match state.output.processed.push((serial, results)) {
            Ok(()) => {
                state.output.processed_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
                state.deadlock_state.record_q5_push();
                did_work = true;
            }
            Err((serial, results)) => {
                // Output full - hold the result for next attempt
                worker.held_processed = Some((serial, results, heap_size));
                break;
            }
        }
    }

    did_work
}

/// Try to execute Step 7: Serialize records.
///
/// This step is parallel - multiple threads can serialize concurrently.
/// Receives a batch of processed items from Q5, serializes each, and
/// concatenates all serialized data into a single `SerializedBatch`.
///
/// # Non-Blocking Design
///
/// Uses the held-item pattern to prevent deadlock. If the output queue is full,
/// the batch is stored in `worker.held_serialized` and the function returns immediately.
fn try_step_serialize<G: Send + 'static, P: Send + MemoryEstimate + 'static>(
    state: &BamPipelineState<G, P>,
    fns: &PipelineFunctions<G, P>,
    worker: &mut WorkerState<P>,
) -> bool {
    // =========================================================================
    // Priority 1: Try to advance any held serialized batch first
    // =========================================================================
    if let Some((serial, held, heap_size)) = worker.held_serialized.take() {
        match state.output.serialized.push((serial, held)) {
            Ok(()) => {
                // Successfully advanced held item
                state.output.serialized_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
                state.deadlock_state.record_q6_push();
            }
            Err((serial, held)) => {
                // Still can't push - put it back and signal output full
                worker.held_serialized = Some((serial, held, heap_size));
                return false;
            }
        }
    }

    // =========================================================================
    // Priority 2: Check if output queue has space
    // =========================================================================
    if state.output.serialized.is_full() {
        return false;
    }

    // =========================================================================
    // Priority 3: Pop input
    // =========================================================================
    let Some((serial, batch)) = state.output.processed.pop() else {
        if let Some(stats) = state.stats() {
            stats.record_queue_empty(5);
        }
        return false;
    };
    state.deadlock_state.record_q5_pop();

    // Track memory being removed from Q5
    let q5_heap_size: usize = batch.iter().map(MemoryEstimate::estimate_heap_size).sum();
    state.output.processed_heap_bytes.fetch_sub(q5_heap_size as u64, Ordering::AcqRel);

    // =========================================================================
    // Priority 4: Serialize all items
    // =========================================================================
    // Prepare worker's serialization buffer
    worker.core.serialization_buffer.clear();

    // Serialize all items into worker's buffer
    let mut total_record_count: u64 = 0;
    for item in batch {
        match (fns.serialize_fn)(item, &mut worker.core.serialization_buffer) {
            Ok(record_count) => {
                total_record_count += record_count;
            }
            Err(e) => {
                state.set_error(e);
                return false;
            }
        }
    }

    // Swap buffer into batch, replace with fresh pre-allocated buffer
    let combined_data = std::mem::replace(
        &mut worker.core.serialization_buffer,
        Vec::with_capacity(SERIALIZATION_BUFFER_CAPACITY),
    );

    // Record serialized bytes for throughput metrics
    if let Some(stats) = state.stats() {
        stats.serialized_bytes.fetch_add(combined_data.len() as u64, Ordering::Relaxed);
    }

    // =========================================================================
    // Priority 5: Try to push result (non-blocking)
    // =========================================================================
    let batch = SerializedBatch { data: combined_data, record_count: total_record_count };
    let heap_size = batch.estimate_heap_size();
    match state.output.serialized.push((serial, batch)) {
        Ok(()) => {
            state.output.serialized_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
            state.deadlock_state.record_q6_push();
            true
        }
        Err((serial, batch)) => {
            // Output full - hold the result for next attempt
            worker.held_serialized = Some((serial, batch, heap_size));
            false
        }
    }
}

/// Try to execute Step 8: Compress to BGZF blocks.
///
/// This step is parallel - multiple threads can compress concurrently.
/// Delegates to the shared implementation which uses the held-item pattern.
fn try_step_compress<G: Send + 'static, P: Send + MemoryEstimate + 'static>(
    state: &BamPipelineState<G, P>,
    worker: &mut WorkerState<P>,
) -> bool {
    shared_try_step_compress(state, worker).is_success()
}

/// Try to execute Step 9: Write blocks to output.
///
/// This step is exclusive - only one thread at a time.
fn try_step_write<G: Send, P: Send + MemoryEstimate>(
    state: &BamPipelineState<G, P>,
) -> (bool, bool) {
    // Try to acquire exclusive access to output file FIRST
    // This avoids wasting time draining if we can't get the lock
    let Some(mut guard) = state.output.output.try_lock() else {
        // Record contention for diagnostics
        if let Some(stats) = state.stats() {
            stats.record_contention(PipelineStep::Write);
        }
        return (false, true); // Contention!
    };

    let Some(ref mut writer) = *guard else {
        return (false, false); // File already closed, not contention
    };

    // Drain Q7 into reorder buffer AND write all ready batches in single lock scope
    let mut wrote_any = false;
    let q7_truly_empty;
    {
        let mut reorder = state.output.write_reorder.lock();

        // Drain Q7 into reorder buffer
        while let Some((serial, batch)) = state.output.compressed.pop() {
            state.deadlock_state.record_q7_pop();
            reorder.insert(serial, batch);
        }

        // Write all ready batches
        while let Some(batch) = reorder.try_pop_next() {
            // Write all blocks in the batch
            let mut batch_bytes: u64 = 0;
            for block in &batch.blocks {
                if let Err(e) = writer.write_all(&block.data) {
                    state.set_error(e);
                    return (false, false); // Error, not contention
                }
                batch_bytes += block.data.len() as u64;
            }

            // Record bytes written for throughput metrics
            if let Some(stats) = state.stats() {
                stats.bytes_written.fetch_add(batch_bytes, Ordering::Relaxed);
            }

            // Use actual record count from the batch
            let records_in_batch = batch.record_count;
            state.output.items_written.fetch_add(records_in_batch, Ordering::Relaxed);
            state.output.progress.log_if_needed(records_in_batch);
            wrote_any = true;
        }

        // Check if truly empty (queue drained and reorder buffer has no pending items)
        q7_truly_empty = reorder.is_empty();
    }

    // Record queue empty only if both Q7 queue AND reorder buffer are empty
    // (not when items are waiting out-of-order in the reorder buffer)
    if !wrote_any && q7_truly_empty {
        if let Some(stats) = state.stats() {
            stats.record_queue_empty(7);
        }
    }

    (wrote_any, false) // Success or no work, no contention (we held lock)
}

// ============================================================================
// Step Context (for consolidated generic_worker_loop)
// ============================================================================

/// Context for BAM pipeline step execution.
///
/// This struct holds references to all the state needed to execute pipeline steps,
/// and implements `StepContext` to work with `generic_worker_loop`.
pub struct BamStepContext<'a, G: Send, P: Send + MemoryEstimate> {
    pub state: &'a BamPipelineState<G, P>,
    pub group_state: &'a Mutex<GroupState<G>>,
    pub fns: &'a PipelineFunctions<G, P>,
    pub is_reader: bool,
}

impl<G, P> StepContext for BamStepContext<'_, G, P>
where
    G: Send + BatchWeight + MemoryEstimate + 'static,
    P: Send + MemoryEstimate + 'static,
{
    type Worker = WorkerState<P>;

    fn execute_step(&self, worker: &mut Self::Worker, step: PipelineStep) -> (bool, bool) {
        execute_step(self.state, self.group_state, self.fns, worker, step)
    }

    fn get_backpressure(&self, _worker: &Self::Worker) -> BackpressureState {
        let depths = self.state.queue_depths();
        let read_done = self.state.read_done.load(Ordering::Relaxed);
        BackpressureState {
            output_high: depths.q7 > self.state.config.output_high_water,
            input_low: depths.q1 < self.state.config.input_low_water,
            read_done,
            memory_high: !self.state.is_draining() && self.state.is_memory_high(),
            memory_drained: self.state.is_memory_drained(),
        }
    }

    fn check_drain_mode(&self) {
        let read_done = self.state.read_done.load(Ordering::Relaxed);
        if read_done && self.state.q1_raw_blocks.is_empty() {
            self.state.output.draining.store(true, Ordering::Relaxed);
        }
    }

    fn has_error(&self) -> bool {
        self.state.has_error()
    }

    fn is_complete(&self) -> bool {
        self.state.is_complete()
    }

    fn stats(&self) -> Option<&PipelineStats> {
        self.state.stats()
    }

    fn skip_read(&self) -> bool {
        // Always skip Read in priority loop:
        // - Readers handle reading via sticky read before the priority loop
        // - Workers don't read at all
        true
    }

    fn check_completion_at_end(&self) -> bool {
        true // Original BAM behavior: check completion at end of loop
    }

    fn should_attempt_sticky_read(&self) -> bool {
        // Outer guard: skip entirely when read_done (original BAM optimization)
        self.is_reader && !self.state.read_done.load(Ordering::Relaxed)
    }

    fn sticky_read_should_continue(&self) -> bool {
        // Full condition checked each iteration
        !self.state.has_error()
            && !self.state.read_done.load(Ordering::Relaxed)
            && self.state.q1_raw_blocks.len() < self.state.config.queue_capacity
    }

    fn execute_read_step(&self, worker: &mut Self::Worker) -> bool {
        try_step_read(self.state, worker)
    }

    fn is_drain_mode(&self) -> bool {
        let read_done = self.state.read_done.load(Ordering::Relaxed);
        let group_done = self.state.group_done.load(Ordering::Relaxed);
        read_done && group_done
    }

    fn should_attempt_step(
        &self,
        worker: &Self::Worker,
        step: PipelineStep,
        drain_mode: bool,
    ) -> bool {
        worker.core.scheduler.should_attempt_step_with_drain(step, drain_mode)
    }

    fn exclusive_step_owned(&self, worker: &Self::Worker) -> Option<PipelineStep> {
        if self.is_reader {
            // Reader thread doesn't use the "try owned first" pattern
            // (it has sticky read instead)
            None
        } else {
            worker.core.scheduler.exclusive_step_owned()
        }
    }
}

/// Execute a single pipeline step, returning (success, `was_contention`).
///
/// `was_contention` indicates if failure was due to lock contention (true)
/// or due to queue being full/empty (false). Only contention failures should
/// be recorded for Thompson Sampling updates.
fn execute_step<
    G: Send + BatchWeight + MemoryEstimate + 'static,
    P: Send + MemoryEstimate + 'static,
>(
    state: &BamPipelineState<G, P>,
    group_state: &Mutex<GroupState<G>>,
    fns: &PipelineFunctions<G, P>,
    worker: &mut WorkerState<P>,
    step: PipelineStep,
) -> (bool, bool) {
    match step {
        PipelineStep::Read => (false, false), // Never read from worker threads
        PipelineStep::Decompress => (try_step_decompress(state, worker), false),
        PipelineStep::FindBoundaries => try_step_find_boundaries(state, worker),
        PipelineStep::Decode => (try_step_decode(state, worker), false),
        PipelineStep::Group => try_step_group(state, group_state),
        PipelineStep::Process => (try_step_process(state, fns, worker), false),
        PipelineStep::Serialize => (try_step_serialize(state, fns, worker), false),
        PipelineStep::Compress => (try_step_compress(state, worker), false),
        PipelineStep::Write => try_step_write(state),
    }
}

// ============================================================================
// Single-Threaded Fast Path
// ============================================================================

/// Pre-allocated buffers for single-threaded pipeline execution.
///
/// These buffers are created once and reused each iteration to avoid
/// per-iteration allocation overhead.
struct SingleThreadedBuffers {
    /// Buffer for concatenated decompressed BGZF block data.
    /// Cleared and reused each iteration.
    decompressed: Vec<u8>,
    /// Buffer for serialized BAM record data.
    /// Cleared and reused each group to avoid allocation.
    serialized: Vec<u8>,
}

impl SingleThreadedBuffers {
    /// Create new buffers with reasonable initial capacity.
    fn new() -> Self {
        Self {
            // 4 blocks * 64KB max uncompressed = 256KB typical
            decompressed: Vec::with_capacity(256 * 1024),
            // Typical group serializes to ~64KB
            serialized: Vec::with_capacity(64 * 1024),
        }
    }
}

/// Run the BAM pipeline in single-threaded mode.
///
/// This avoids the overhead of thread spawning, queues, and atomic
/// operations when only one thread is requested. Significantly faster
/// for small inputs or when parallelization overhead exceeds the benefit.
fn run_bam_pipeline_single_threaded<G, P>(
    config: &PipelineConfig,
    mut input: Box<dyn Read + Send>,
    mut output: Box<dyn Write + Send>,
    mut grouper: Box<dyn Grouper<Group = G> + Send>,
    fns: PipelineFunctions<G, P>,
    group_key_config: GroupKeyConfig,
) -> io::Result<u64>
where
    G: Send + 'static,
    P: Send + MemoryEstimate + 'static,
{
    // Step 1+2: Reader and decompressor
    let mut decompressor = libdeflater::Decompressor::new();

    // Step 3: Boundary finder state (skip header if already read)
    let mut boundary_state = if config.header_already_read {
        BoundaryState::new_no_header()
    } else {
        BoundaryState::new()
    };

    // Step 8: Compressor
    let mut compressor = InlineBgzfCompressor::new(config.compression_level);

    // Pre-allocated reusable buffers
    let mut buffers = SingleThreadedBuffers::new();

    // Progress tracking
    let progress = ProgressTracker::new("Processed records").with_interval(PROGRESS_LOG_INTERVAL);

    // Main loop: read -> decompress -> find_boundaries -> decode -> group -> process -> serialize -> compress -> write
    loop {
        // Step 1: Read a batch of raw BGZF blocks
        let blocks = read_raw_blocks(input.as_mut(), 4)?; // Read 4 blocks at a time
        if blocks.is_empty() {
            break; // EOF
        }

        // Clear decompression buffer for reuse (keeps capacity)
        buffers.decompressed.clear();

        // Step 2: Decompress all blocks into reusable buffer
        let expected_size: usize =
            blocks.iter().map(super::super::bgzf_reader::RawBgzfBlock::uncompressed_size).sum();
        buffers.decompressed.reserve(expected_size);

        for block in &blocks {
            decompress_block_into(block, &mut decompressor, &mut buffers.decompressed)?;
        }

        // Step 3: Find record boundaries
        let boundary_batch = boundary_state.find_boundaries(&buffers.decompressed)?;

        // Step 4: Decode records (only if there are any)
        if boundary_batch.offsets.len() > 1 {
            let decoded = decode_records(&boundary_batch, &group_key_config)?;

            // Step 5: Feed decoded records to grouper
            let groups = grouper.add_records(decoded)?;

            // Process each group through steps 6-9
            for group in groups {
                // Step 6: Process
                let processed = (fns.process_fn)(group)?;

                // Step 7: Serialize (reuse buffer)
                buffers.serialized.clear();
                let record_count = (fns.serialize_fn)(processed, &mut buffers.serialized)?;

                // Step 8: Compress (only when buffer reaches 64KB)
                compressor.write_all(&buffers.serialized)?;
                compressor.maybe_compress()?;

                // Step 9: Write any completed blocks to output
                compressor.write_blocks_to(output.as_mut())?;

                progress.log_if_needed(record_count);
            }
        }
    }

    // Handle any remaining bytes from boundary finding
    if let Some(final_batch) = boundary_state.finish()? {
        if final_batch.offsets.len() > 1 {
            let decoded = decode_records(&final_batch, &group_key_config)?;
            let groups = grouper.add_records(decoded)?;

            for group in groups {
                let processed = (fns.process_fn)(group)?;
                buffers.serialized.clear();
                let record_count = (fns.serialize_fn)(processed, &mut buffers.serialized)?;
                compressor.write_all(&buffers.serialized)?;
                compressor.maybe_compress()?;
                compressor.write_blocks_to(output.as_mut())?;
                progress.log_if_needed(record_count);
            }
        }
    }

    // Finish grouper - process any remaining partial group
    if let Some(final_group) = grouper.finish()? {
        // Step 6: Process
        let processed = (fns.process_fn)(final_group)?;

        // Step 7: Serialize (reuse buffer)
        buffers.serialized.clear();
        let record_count = (fns.serialize_fn)(processed, &mut buffers.serialized)?;

        // Step 8: Compress (only when buffer reaches 64KB)
        compressor.write_all(&buffers.serialized)?;
        compressor.maybe_compress()?;

        // Step 9: Write any completed blocks to output
        compressor.write_blocks_to(output.as_mut())?;

        progress.log_if_needed(record_count);
    }

    // Flush any remaining data in compression buffer
    compressor.flush()?;
    compressor.write_blocks_to(output.as_mut())?;

    // Flush output
    output.flush()?;

    Ok(progress.count())
}

// ============================================================================
// Public Run Function
// ============================================================================

/// Run the BAM pipeline.
///
/// # Type Parameters
///
/// - `G`: Group type produced by the grouper
/// - `P`: Processed type produced by the process function
///
/// # Arguments
///
/// - `config`: Pipeline configuration
/// - `input`: Input reader (e.g., BAM file)
/// - `output`: Output writer (e.g., BAM file)
/// - `grouper`: The grouper that groups decoded records
/// - `fns`: Step functions for processing and serialization
///
/// # Returns
///
/// Number of groups processed, or an error.
pub fn run_bam_pipeline<G, P>(
    config: PipelineConfig,
    input: Box<dyn Read + Send>,
    output: Box<dyn Write + Send>,
    grouper: Box<dyn Grouper<Group = G> + Send>,
    fns: PipelineFunctions<G, P>,
    group_key_config: GroupKeyConfig,
) -> io::Result<u64>
where
    G: Send + BatchWeight + MemoryEstimate + 'static,
    P: Send + MemoryEstimate + 'static,
{
    let num_threads = config.num_threads;
    let compression_level = config.compression_level;
    let scheduler_strategy = config.scheduler_strategy;

    // ============================================================
    // Single-threaded fast path: avoid thread/queue overhead
    // ============================================================
    if num_threads == 1 {
        return run_bam_pipeline_single_threaded(
            &config,
            input,
            output,
            grouper,
            fns,
            group_key_config,
        );
    }

    let state = Arc::new(BamPipelineState::<G, P>::new(config, input, output, group_key_config));

    // Set num_threads for per-thread stats display
    if let Some(stats) = state.stats() {
        stats.set_num_threads(num_threads);
    }

    let group_state = Arc::new(Mutex::new(GroupState::new(grouper)));
    let fns = Arc::new(fns);

    // Spawn worker threads
    // Thread 0 is the sticky reader, threads 1..N-1 are workers only
    let handles: Vec<_> = (0..num_threads)
        .map(|thread_id| {
            let state = Arc::clone(&state);
            let group_state = Arc::clone(&group_state);
            let fns = Arc::clone(&fns);

            thread::spawn(move || {
                // Wrap worker logic in catch_unwind to handle panics gracefully
                let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    let mut worker = WorkerState::new(
                        compression_level,
                        thread_id,
                        num_threads,
                        scheduler_strategy,
                    );
                    let ctx = BamStepContext {
                        state: &state,
                        group_state: &group_state,
                        fns: &fns,
                        is_reader: thread_id == 0,
                    };
                    generic_worker_loop(&ctx, &mut worker);
                }));

                // If a panic occurred, set the error flag so other threads exit
                if let Err(panic_info) = result {
                    handle_worker_panic(&*state, thread_id, panic_info);
                }
            })
        })
        .collect();

    // Spawn queue monitor thread if stats or deadlock detection are enabled
    let monitor_handle = if state.stats().is_some() || state.deadlock_state.is_enabled() {
        let state_clone = Arc::clone(&state);
        Some(thread::spawn(move || {
            let start_time = Instant::now();
            let mut deadlock_check_counter = 0u32;
            loop {
                // Sleep 100ms between samples
                thread::sleep(Duration::from_millis(100));

                // Exit if pipeline is done or has error
                if state_clone.is_complete() || state_clone.has_error() {
                    break;
                }

                // Collect queue sizes (needed for both stats and deadlock detection)
                let queue_sizes = [
                    state_clone.q1_raw_blocks.len(),
                    state_clone.q2_decompressed.len(),
                    state_clone.q2b_boundaries.len(),
                    state_clone.q3_decoded.len(),
                    state_clone.output.groups.len(),
                    state_clone.output.processed.len(),
                    state_clone.output.serialized.len(),
                    state_clone.output.compressed.len(),
                ];

                // Collect reorder buffer sizes and memory (need locks)
                let (q2_reorder_len, q2_reorder_mem) = {
                    let reorder = state_clone.q2_reorder.lock();
                    (reorder.len(), reorder.total_heap_size() as u64)
                };
                let (q3_reorder_len, q3_reorder_mem) = {
                    let reorder = state_clone.q3_reorder.lock();
                    (reorder.len(), reorder.total_heap_size() as u64)
                };
                let (q7_reorder_len, q7_reorder_mem) = {
                    let reorder = state_clone.output.write_reorder.lock();
                    (reorder.len(), reorder.total_heap_size() as u64)
                };
                let reorder_sizes = [q2_reorder_len, q3_reorder_len, q7_reorder_len];
                let reorder_memory_bytes = [q2_reorder_mem, q3_reorder_mem, q7_reorder_mem];

                // Queue memory from AtomicU64 counters (Q4-Q7)
                let q4_mem = state_clone.output.groups_heap_bytes.load(Ordering::Relaxed);
                let q5_mem = state_clone.output.processed_heap_bytes.load(Ordering::Relaxed);
                let q6_mem = state_clone.output.serialized_heap_bytes.load(Ordering::Relaxed);
                let q7_mem = state_clone.output.compressed_heap_bytes.load(Ordering::Relaxed);
                let queue_memory_bytes = [0, 0, 0, 0, q4_mem, q5_mem, q6_mem, q7_mem];

                // Collect thread activity
                let thread_steps: Vec<u8> = if let Some(stats) = state_clone.stats() {
                    let num_threads = stats.num_threads.load(Ordering::Relaxed) as usize;
                    (0..num_threads)
                        .map(|tid| stats.per_thread_current_step[tid].load(Ordering::Relaxed))
                        .collect()
                } else {
                    Vec::new()
                };

                // Record sample and track peak memory
                if let Some(stats) = state_clone.stats() {
                    // Track peak memory from all queues (reorder buffers + ArrayQueues)
                    let total_mem = q2_reorder_mem
                        + q3_reorder_mem
                        + q7_reorder_mem
                        + q4_mem
                        + q5_mem
                        + q6_mem
                        + q7_mem;
                    stats.record_memory_usage(total_mem);

                    stats.add_queue_sample(QueueSample {
                        time_ms: start_time.elapsed().as_millis() as u64,
                        queue_sizes,
                        reorder_sizes,
                        queue_memory_bytes,
                        reorder_memory_bytes,
                        thread_steps,
                    });
                }

                // Check for deadlock every ~1 second (10 iterations * 100ms)
                if state_clone.deadlock_state.is_enabled() {
                    deadlock_check_counter += 1;
                    if deadlock_check_counter >= 10 {
                        deadlock_check_counter = 0;
                        let snapshot = state_clone.build_queue_snapshot();
                        check_deadlock_and_restore(&state_clone.deadlock_state, &snapshot);
                    }
                }
            }
        }))
    } else {
        None
    };

    // Wait for all threads to complete
    join_worker_threads(handles)?;
    join_monitor_thread(monitor_handle);

    // Finalize: check errors, flush output, log stats
    finalize_pipeline(&*state)
}

// ============================================================================
// BAM Pipeline Helpers
// ============================================================================

// Thread-local buffer for serializing BAM records.
// Reusing this buffer across calls avoids repeated allocations.
thread_local! {
    static SERIALIZE_RECORD_BUFFER: std::cell::RefCell<Vec<u8>> =
        std::cell::RefCell::new(Vec::with_capacity(512));
}

/// Serialize a batch of BAM records to bytes.
///
/// This function encodes multiple BAM records into a single byte buffer,
/// ready for BGZF compression. Uses thread-local buffer to avoid allocations.
/// Serialize BAM records into a provided buffer.
///
/// Appends serialized BAM bytes to the provided buffer and returns the record count.
/// This variant is more efficient when the caller wants to reuse a buffer.
pub fn serialize_bam_records_into(
    records: &[RecordBuf],
    header: &Header,
    output: &mut Vec<u8>,
) -> io::Result<u64> {
    use noodles::bam::record::codec::encode_record_buf;

    log::trace!("serialize_bam_records_into: {} records", records.len());

    // Pre-allocate output buffer based on batch size estimate
    // Typical record: ~300-500 bytes, plus 4-byte block_size prefix
    // Use 400 bytes as a reasonable average
    let estimated_batch_size = records.len() * 400;
    output.reserve(estimated_batch_size);

    SERIALIZE_RECORD_BUFFER.with(|buf| {
        let mut record_data = buf.borrow_mut();

        for (i, record) in records.iter().enumerate() {
            // Clear and reuse the buffer for each record
            record_data.clear();
            if let Err(e) = encode_record_buf(&mut record_data, header, record) {
                log::error!(
                    "serialize_bam_records_into: failed to encode record {}: {:?}, name={:?}, seq_len={}, qual_len={}",
                    i,
                    e,
                    record.name(),
                    record.sequence().len(),
                    record.quality_scores().len(),
                );
                return Err(e);
            }

            // Add `block_size` prefix + record data
            let block_size = record_data.len() as u32;
            output.extend_from_slice(&block_size.to_le_bytes());
            output.extend_from_slice(&record_data);
        }

        Ok(records.len() as u64)
    })
}

/// Serialize BAM records to a new `SerializedBatch`.
///
/// This is a convenience wrapper around `serialize_bam_records_into` that allocates
/// a new buffer. Use `serialize_bam_records_into` for better performance when
/// buffer reuse is possible.
pub fn serialize_bam_records(
    records: &[RecordBuf],
    header: &Header,
) -> io::Result<SerializedBatch> {
    let mut data = Vec::with_capacity(records.len() * 256);
    let record_count = serialize_bam_records_into(records, header, &mut data)?;
    Ok(SerializedBatch { data, record_count })
}

/// Serialize a single BAM record to bytes.
///
/// This produces raw BAM record bytes (`block_size` prefix + record data),
/// suitable for BGZF compression in the pipeline. Uses thread-local buffer.
pub fn serialize_bam_record(record: &RecordBuf, header: &Header) -> io::Result<SerializedBatch> {
    let mut data = Vec::with_capacity(256);
    let record_count = serialize_bam_record_into(record, header, &mut data)?;
    Ok(SerializedBatch { data, record_count })
}

/// Serialize a single BAM record to bytes, appending to the provided buffer.
///
/// This produces raw BAM record bytes (`block_size` prefix + record data),
/// suitable for BGZF compression in the pipeline. Uses thread-local buffer for encoding.
///
/// Returns the number of records serialized (always 1 for a single record).
pub fn serialize_bam_record_into(
    record: &RecordBuf,
    header: &Header,
    output: &mut Vec<u8>,
) -> io::Result<u64> {
    use noodles::bam::record::codec::encode_record_buf;

    SERIALIZE_RECORD_BUFFER.with(|buf| {
        let mut record_data = buf.borrow_mut();
        record_data.clear();
        encode_record_buf(&mut record_data, header, record)?;

        // Append `block_size` prefix + record data to output
        let block_size = record_data.len() as u32;
        output.extend_from_slice(&block_size.to_le_bytes());
        output.extend_from_slice(&record_data);

        Ok(1)
    })
}

/// Serialize BAM records directly to a BGZF compressor (zero-copy).
///
/// This writes records directly to the compressor's internal buffer, avoiding
/// the intermediate serialization buffer copy. Records are compressed into
/// BGZF blocks as the buffer fills.
///
/// Returns the number of records serialized.
pub fn serialize_bam_records_to_compressor(
    records: &[RecordBuf],
    header: &Header,
    compressor: &mut crate::bgzf_writer::InlineBgzfCompressor,
) -> io::Result<u64> {
    use noodles::bam::record::codec::encode_record_buf;

    log::trace!("serialize_bam_records_to_compressor: {} records", records.len());

    SERIALIZE_RECORD_BUFFER.with(|buf| {
        let mut record_data = buf.borrow_mut();

        for (i, record) in records.iter().enumerate() {
            // Clear and reuse the buffer for each record
            record_data.clear();
            if let Err(e) = encode_record_buf(&mut record_data, header, record) {
                log::error!(
                    "serialize_bam_records_to_compressor: failed to encode record {}: {:?}, name={:?}, seq_len={}, qual_len={}",
                    i,
                    e,
                    record.name(),
                    record.sequence().len(),
                    record.quality_scores().len(),
                );
                return Err(e);
            }

            // Write block_size prefix + record data directly to compressor buffer
            let block_size = record_data.len() as u32;
            let buffer = compressor.buffer_mut();
            buffer.extend_from_slice(&block_size.to_le_bytes());
            buffer.extend_from_slice(&record_data);

            // Compress if buffer is full
            compressor.maybe_compress()?;
        }

        Ok(records.len() as u64)
    })
}

/// Configuration for running a BAM file through the pipeline.
#[derive(Debug, Clone)]
pub struct BamPipelineConfig {
    /// Base pipeline configuration.
    pub pipeline: PipelineConfig,
    /// Compression level for output (0-12).
    pub compression_level: u32,
    /// Configuration for computing `GroupKey` during decode.
    /// If None, a default config will be built from the header.
    pub group_key_config: Option<GroupKeyConfig>,
}

impl BamPipelineConfig {
    /// Create a new BAM pipeline configuration.
    #[must_use]
    pub fn new(num_threads: usize, compression_level: u32) -> Self {
        Self {
            pipeline: PipelineConfig::new(num_threads, compression_level),
            compression_level,
            group_key_config: None,
        }
    }

    /// Create a configuration auto-tuned for the given thread count.
    ///
    /// This adjusts queue capacity and batch sizes based on the number of threads
    /// to optimize throughput and reduce contention.
    #[must_use]
    pub fn auto_tuned(num_threads: usize, compression_level: u32) -> Self {
        Self {
            pipeline: PipelineConfig::auto_tuned(num_threads, compression_level),
            compression_level,
            group_key_config: None,
        }
    }

    /// Set the compression level.
    #[must_use]
    pub fn with_compression_level(mut self, level: u32) -> Self {
        self.compression_level = level;
        self.pipeline.compression_level = level;
        self
    }

    /// Set the `GroupKey` configuration for position-based grouping.
    #[must_use]
    pub fn with_group_key_config(mut self, config: GroupKeyConfig) -> Self {
        self.group_key_config = Some(config);
        self
    }
}

/// Run a BAM file through the pipeline with a grouper factory.
///
/// This is a convenience function that handles BAM header I/O and
/// sets up the pipeline correctly for BAM processing.
///
/// # Type Parameters
///
/// - `G`: Group type produced by the grouper (e.g., `RecordBuf`, `PositionGroup`, `MiGroup`)
/// - `P`: Processed type produced by the process function (e.g., `Vec<RecordBuf>`)
///
/// # Arguments
///
/// - `config`: Pipeline configuration
/// - `input_path`: Path to input BAM file
/// - `output_path`: Path to output BAM file
/// - `grouper_fn`: Function that creates a grouper given the header
/// - `process_fn`: Function to process each group
/// - `serialize_fn`: Function to serialize processed output (receives header reference and output buffer)
///
/// # Returns
///
/// Number of groups processed, or an error.
pub fn run_bam_pipeline_with_grouper<G, P, GrouperFn, ProcessFn, SerializeFn>(
    config: BamPipelineConfig,
    input_path: &Path,
    output_path: &Path,
    grouper_fn: GrouperFn,
    process_fn: ProcessFn,
    serialize_fn: SerializeFn,
) -> io::Result<u64>
where
    G: Send + BatchWeight + MemoryEstimate + 'static,
    P: Send + MemoryEstimate + 'static,
    GrouperFn: FnOnce(&Header) -> Box<dyn Grouper<Group = G> + Send>,
    ProcessFn: Fn(G) -> io::Result<P> + Send + Sync + 'static,
    SerializeFn: Fn(P, &Header, &mut Vec<u8>) -> io::Result<u64> + Send + Sync + 'static,
{
    // First pass: Read header only (we need it for grouper creation and output writing)
    let header = {
        let input_file = File::open(input_path)
            .map_err(|e| io::Error::new(e.kind(), format!("Failed to open input: {e}")))?;
        let mut bam_reader = bam::io::Reader::new(input_file);
        bam_reader.read_header().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Failed to read BAM header: {e}"))
        })?
    };

    // Create output BAM and write header
    let output_file = File::create(output_path)
        .map_err(|e| io::Error::new(e.kind(), format!("Failed to create output: {e}")))?;

    // Write BAM header using BGZF compression
    let mut header_writer = bam::io::Writer::new(output_file);
    header_writer
        .write_header(&header)
        .map_err(|e| io::Error::other(format!("Failed to write BAM header: {e}")))?;

    // Finish the BGZF writer and get the underlying file handle for the pipeline.
    // We need to:
    // 1. Get the BGZF writer from the BAM writer
    // 2. Flush/finish the BGZF stream (writes any pending data)
    // 3. Get the underlying file handle
    // This ensures the pipeline writes raw BGZF blocks directly to the file,
    // not through another BGZF compression layer.
    let mut bgzf_writer = header_writer.into_inner();
    bgzf_writer
        .try_finish()
        .map_err(|e| io::Error::other(format!("Failed to finish BGZF header: {e}")))?;
    let output = bgzf_writer.into_inner();

    // Re-open input file for pipeline (raw file reader for BGZF block reading)
    // Wrap in BufReader to reduce syscalls
    let input = File::open(input_path)
        .map_err(|e| io::Error::new(e.kind(), format!("Failed to re-open input: {e}")))?;
    let input = BufReader::with_capacity(IO_BUFFER_SIZE, input);

    // Wrap output in BufWriter to reduce syscalls
    let output = BufWriter::with_capacity(IO_BUFFER_SIZE, output);

    // Build GroupKeyConfig from header if not provided
    let group_key_config = config.group_key_config.unwrap_or_else(|| {
        use noodles::sam::alignment::record::data::field::Tag;
        let library_index = LibraryIndex::from_header(&header);
        let cell_tag = Tag::from([b'C', b'B']); // Default cell tag
        GroupKeyConfig::new(library_index, cell_tag)
    });

    // Create the grouper
    // Note: Header skipping is now handled by BoundaryState in the pipeline
    let grouper = grouper_fn(&header);

    // Create step functions with header captured
    let header_clone = header.clone();
    let fns = PipelineFunctions::new(process_fn, move |p: P, buf: &mut Vec<u8>| {
        serialize_fn(p, &header_clone, buf)
    });

    // Run the pipeline
    let result = run_bam_pipeline(
        config.pipeline,
        Box::new(input),
        Box::new(output),
        grouper,
        fns,
        group_key_config,
    );

    // After pipeline completes, write BGZF EOF block to finalize the BAM file
    if result.is_ok() {
        use std::io::Write as _;

        // Append EOF block to output file
        let mut output_file = std::fs::OpenOptions::new()
            .append(true)
            .open(output_path)
            .map_err(|e| io::Error::new(e.kind(), format!("Failed to open output for EOF: {e}")))?;

        output_file
            .write_all(&BGZF_EOF)
            .map_err(|e| io::Error::new(e.kind(), format!("Failed to write BGZF EOF: {e}")))?;
    }

    result
}

/// Run a BAM file through the pipeline with a custom output header.
///
/// This variant allows specifying a different output header than the input header,
/// useful for commands like `simplex` that produce unmapped consensus reads.
///
/// # Type Parameters
///
/// - `G`: Group type produced by the grouper (e.g., `RecordBuf`, `PositionGroup`, `MiGroup`)
/// - `P`: Processed type produced by the process function (e.g., `Vec<RecordBuf>`)
///
/// # Arguments
///
/// - `config`: Pipeline configuration
/// - `input_path`: Path to input BAM file
/// - `output_path`: Path to output BAM file
/// - `output_header`: Custom header to write to output file (and use for serialization)
/// - `grouper_fn`: Function that creates a grouper given the input header
/// - `process_fn`: Function to process each group
/// - `serialize_fn`: Function to serialize processed output (receives output header reference)
///
/// # Returns
///
/// Number of groups processed, or an error.
pub fn run_bam_pipeline_with_header<G, P, GrouperFn, ProcessFn, SerializeFn>(
    config: BamPipelineConfig,
    input_path: &Path,
    output_path: &Path,
    output_header: Header,
    grouper_fn: GrouperFn,
    process_fn: ProcessFn,
    serialize_fn: SerializeFn,
) -> io::Result<u64>
where
    G: Send + BatchWeight + MemoryEstimate + 'static,
    P: Send + MemoryEstimate + 'static,
    GrouperFn: FnOnce(&Header) -> Box<dyn Grouper<Group = G> + Send>,
    ProcessFn: Fn(G) -> io::Result<P> + Send + Sync + 'static,
    SerializeFn: Fn(P, &Header, &mut Vec<u8>) -> io::Result<u64> + Send + Sync + 'static,
{
    // First pass: Read input header only (for grouper creation)
    let input_header = {
        let input_file = File::open(input_path)
            .map_err(|e| io::Error::new(e.kind(), format!("Failed to open input: {e}")))?;
        let mut bam_reader = bam::io::Reader::new(input_file);
        bam_reader.read_header().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Failed to read BAM header: {e}"))
        })?
    };

    // Create output BAM and write the custom output header
    let output_file = File::create(output_path)
        .map_err(|e| io::Error::new(e.kind(), format!("Failed to create output: {e}")))?;

    // Write BAM header using BGZF compression
    let mut header_writer = bam::io::Writer::new(output_file);
    header_writer
        .write_header(&output_header)
        .map_err(|e| io::Error::other(format!("Failed to write BAM header: {e}")))?;

    // Finish the BGZF writer and get the underlying file handle for the pipeline.
    // We need to:
    // 1. Get the BGZF writer from the BAM writer
    // 2. Flush/finish the BGZF stream (writes any pending data + EOF)
    // 3. Get the underlying file handle
    let mut bgzf_writer = header_writer.into_inner();
    bgzf_writer
        .try_finish()
        .map_err(|e| io::Error::other(format!("Failed to finish BGZF header: {e}")))?;
    let output = bgzf_writer.into_inner();

    // Re-open input file for pipeline (raw file reader for BGZF block reading)
    // Wrap in BufReader to reduce syscalls
    let input = File::open(input_path)
        .map_err(|e| io::Error::new(e.kind(), format!("Failed to re-open input: {e}")))?;
    let input = BufReader::with_capacity(IO_BUFFER_SIZE, input);

    // Wrap output in BufWriter to reduce syscalls
    let output = BufWriter::with_capacity(IO_BUFFER_SIZE, output);

    // Build GroupKeyConfig from input header if not provided
    let group_key_config = config.group_key_config.unwrap_or_else(|| {
        use noodles::sam::alignment::record::data::field::Tag;
        let library_index = LibraryIndex::from_header(&input_header);
        let cell_tag = Tag::from([b'C', b'B']); // Default cell tag
        GroupKeyConfig::new(library_index, cell_tag)
    });

    // Create the grouper using INPUT header
    // Note: Header skipping is now handled by BoundaryState in the pipeline
    let grouper = grouper_fn(&input_header);

    // Create step functions with OUTPUT header captured (for serialization)
    let fns = PipelineFunctions::new(process_fn, move |p: P, buf: &mut Vec<u8>| {
        serialize_fn(p, &output_header, buf)
    });

    // Run the pipeline
    let result = run_bam_pipeline(
        config.pipeline,
        Box::new(input),
        Box::new(output),
        grouper,
        fns,
        group_key_config,
    );

    // After pipeline completes, write BGZF EOF block to finalize the BAM file
    if result.is_ok() {
        use std::io::Write as _;

        // Append EOF block to output file
        let mut output_file = std::fs::OpenOptions::new()
            .append(true)
            .open(output_path)
            .map_err(|e| io::Error::new(e.kind(), format!("Failed to open output for EOF: {e}")))?;

        output_file
            .write_all(&BGZF_EOF)
            .map_err(|e| io::Error::new(e.kind(), format!("Failed to write BGZF EOF: {e}")))?;
    }

    result
}

// ============================================================================
// Reader-based Pipeline Functions (for streaming support)
// ============================================================================

/// Run a BAM pipeline from an already-opened reader.
///
/// This variant accepts a pre-opened reader and header, enabling streaming from
/// stdin or other non-seekable sources.
///
/// # Type Parameters
///
/// - `G`: Group type produced by the grouper (e.g., `RecordBuf`, `PositionGroup`, `MiGroup`)
/// - `P`: Processed type produced by the process function (e.g., `Vec<RecordBuf>`)
/// - `R`: Reader type that implements `Read + Send`
///
/// # Arguments
///
/// - `config`: Pipeline configuration
/// - `input`: Pre-opened input reader (header already read)
/// - `input_header`: Header that was read from the input (used for grouping)
/// - `output_path`: Path to output BAM file
/// - `output_header`: Optional custom header for output file and serialization.
///   If `None`, uses `input_header` for both.
/// - `grouper_fn`: Function that creates a grouper given the input header
/// - `process_fn`: Function to process each group
/// - `serialize_fn`: Function to serialize processed output (receives output header reference)
///
/// # Returns
///
/// Number of groups processed, or an error.
#[allow(clippy::too_many_arguments)]
pub fn run_bam_pipeline_from_reader<G, P, R, GrouperFn, ProcessFn, SerializeFn>(
    config: BamPipelineConfig,
    input: R,
    input_header: Header,
    output_path: &Path,
    output_header: Option<Header>,
    grouper_fn: GrouperFn,
    process_fn: ProcessFn,
    serialize_fn: SerializeFn,
) -> io::Result<u64>
where
    G: Send + BatchWeight + MemoryEstimate + 'static,
    P: Send + MemoryEstimate + 'static,
    R: Read + Send + 'static,
    GrouperFn: FnOnce(&Header) -> Box<dyn Grouper<Group = G> + Send>,
    ProcessFn: Fn(G) -> io::Result<P> + Send + Sync + 'static,
    SerializeFn: Fn(P, &Header, &mut Vec<u8>) -> io::Result<u64> + Send + Sync + 'static,
{
    // Use output_header if provided, otherwise clone input_header
    let output_header = output_header.unwrap_or_else(|| input_header.clone());

    // Create output BAM and write the output header
    let output_file = File::create(output_path)
        .map_err(|e| io::Error::new(e.kind(), format!("Failed to create output: {e}")))?;

    // Write BAM header using BGZF compression
    let mut header_writer = bam::io::Writer::new(output_file);
    header_writer
        .write_header(&output_header)
        .map_err(|e| io::Error::other(format!("Failed to write BAM header: {e}")))?;

    // Finish the BGZF writer and get the underlying file handle for the pipeline.
    let mut bgzf_writer = header_writer.into_inner();
    bgzf_writer
        .try_finish()
        .map_err(|e| io::Error::other(format!("Failed to finish BGZF header: {e}")))?;
    let output = bgzf_writer.into_inner();

    // Wrap output in BufWriter to reduce syscalls
    let output = BufWriter::with_capacity(IO_BUFFER_SIZE, output);

    // Build GroupKeyConfig from input header if not provided
    let group_key_config = config.group_key_config.unwrap_or_else(|| {
        use noodles::sam::alignment::record::data::field::Tag;
        let library_index = LibraryIndex::from_header(&input_header);
        let cell_tag = Tag::from([b'C', b'B']); // Default cell tag
        GroupKeyConfig::new(library_index, cell_tag)
    });

    // Create the grouper using INPUT header
    let grouper = grouper_fn(&input_header);

    // Create step functions with OUTPUT header captured (for serialization)
    let fns = PipelineFunctions::new(process_fn, move |p: P, buf: &mut Vec<u8>| {
        serialize_fn(p, &output_header, buf)
    });

    // Run the pipeline with the already-opened reader.
    // NOTE: The input stream starts at position 0 (including header bytes), so the pipeline
    // must still skip the header. We don't set header_already_read since the bytes are present.
    let pipeline_config = config.pipeline;

    let result = run_bam_pipeline(
        pipeline_config,
        Box::new(input),
        Box::new(output),
        grouper,
        fns,
        group_key_config,
    );

    // After pipeline completes, write BGZF EOF block to finalize the BAM file
    if result.is_ok() {
        use std::io::Write as _;

        // Append EOF block to output file
        let mut output_file = std::fs::OpenOptions::new()
            .append(true)
            .open(output_path)
            .map_err(|e| io::Error::new(e.kind(), format!("Failed to open output for EOF: {e}")))?;

        output_file
            .write_all(&BGZF_EOF)
            .map_err(|e| io::Error::new(e.kind(), format!("Failed to write BGZF EOF: {e}")))?;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_info::LibraryIndex;

    /// Create a minimal `BamPipelineState` for testing memory backpressure.
    fn create_test_state(memory_limit: u64) -> BamPipelineState<(), ()> {
        let config = PipelineConfig::new(2, 6).with_queue_memory_limit(memory_limit);
        let input: Box<dyn Read + Send> = Box::new(std::io::empty());
        let output: Box<dyn Write + Send> = Box::new(std::io::sink());
        // Create minimal GroupKeyConfig for testing
        let header = Header::default();
        let library_index = LibraryIndex::from_header(&header);
        let cell_tag: [u8; 2] = [b'C', b'B'];
        let group_key_config = GroupKeyConfig::new(library_index, cell_tag.into());
        BamPipelineState::new(config, input, output, group_key_config)
    }

    #[test]
    fn test_can_decompress_proceed_no_limit() {
        let state = create_test_state(0); // No limit
        // Should always proceed when no limit
        assert!(state.can_decompress_proceed(0));
        assert!(state.can_decompress_proceed(100));
    }

    #[test]
    fn test_can_decompress_proceed_under_limit() {
        let state = create_test_state(1024 * 1024); // 1MB limit
        // Under 50% of limit, should proceed
        state.q2_reorder_state.heap_bytes.store(100_000, Ordering::SeqCst);
        assert!(state.can_decompress_proceed(5));
    }

    #[test]
    fn test_can_decompress_proceed_over_limit_but_needed_serial() {
        let state = create_test_state(1024 * 1024); // 1MB limit
        // Over 50% of limit
        state.q2_reorder_state.heap_bytes.store(600_000, Ordering::SeqCst);
        state.q2_reorder_state.next_seq.store(5, Ordering::SeqCst);
        // Should still proceed for the needed serial (deadlock prevention)
        assert!(state.can_decompress_proceed(5));
        // But not for other serials
        assert!(!state.can_decompress_proceed(6));
        assert!(!state.can_decompress_proceed(10));
    }

    #[test]
    fn test_can_decompress_proceed_over_limit() {
        let state = create_test_state(1024 * 1024); // 1MB limit
        // Over 50% of limit
        state.q2_reorder_state.heap_bytes.store(600_000, Ordering::SeqCst);
        state.q2_reorder_state.next_seq.store(0, Ordering::SeqCst);
        // Should not proceed for non-needed serials
        assert!(!state.can_decompress_proceed(5));
    }

    #[test]
    fn test_can_decode_proceed_no_limit() {
        let state = create_test_state(0); // No limit (uses default 512MB threshold)
        // Under default threshold, should proceed
        assert!(state.can_decode_proceed(0));
        assert!(state.can_decode_proceed(100));
    }

    #[test]
    fn test_can_decode_proceed_under_limit() {
        let state = create_test_state(1024 * 1024); // 1MB limit
        // Under 50% of limit, should proceed
        state.q3_reorder_state.heap_bytes.store(100_000, Ordering::SeqCst);
        assert!(state.can_decode_proceed(5));
    }

    #[test]
    fn test_can_decode_proceed_over_limit_but_needed_serial() {
        let state = create_test_state(1024 * 1024); // 1MB limit
        // Over 50% of limit
        state.q3_reorder_state.heap_bytes.store(600_000, Ordering::SeqCst);
        state.q3_reorder_state.next_seq.store(5, Ordering::SeqCst);
        // Should still proceed for the needed serial (deadlock prevention)
        assert!(state.can_decode_proceed(5));
        // But not for other serials
        assert!(!state.can_decode_proceed(6));
    }

    #[test]
    fn test_is_memory_high_threshold() {
        let state = create_test_state(1024 * 1024 * 1024); // 1GB limit (uses 512MB cap)
        // Under 512MB threshold
        state.q3_reorder_state.heap_bytes.store(500 * 1024 * 1024, Ordering::SeqCst);
        assert!(!state.is_memory_high());
        // At 512MB threshold
        state.q3_reorder_state.heap_bytes.store(512 * 1024 * 1024, Ordering::SeqCst);
        assert!(state.is_memory_high());
    }

    #[test]
    fn test_is_memory_drained_threshold() {
        let state = create_test_state(1024 * 1024 * 1024); // 1GB limit (uses 512MB cap)
        // Under 256MB (half of 512MB), should be drained
        state.q3_reorder_state.heap_bytes.store(200 * 1024 * 1024, Ordering::SeqCst);
        assert!(state.is_memory_drained());
        // At 256MB, not drained
        state.q3_reorder_state.heap_bytes.store(256 * 1024 * 1024, Ordering::SeqCst);
        assert!(!state.is_memory_drained());
    }
}
