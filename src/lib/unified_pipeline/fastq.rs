//! Unified FASTQ pipeline with configurable active steps.
//!
//! This module implements the FASTQ-specific pipeline for the `extract` command,
//! which reads from multiple synchronized FASTQ streams (R1, R2, I1, I2) and
//! groups reads by template name.
//!
//! A single unified pipeline handles both gzip and BGZF inputs. The set of active
//! pipeline steps is configured via [`ActiveSteps`] based on input type and mode:
//!
//! - **Gzip + synchronized** (e.g. `extract`): Read → Decode → Process → Serialize →
//!   Compress → Write. Uses [`MultiRecordCountReader`] to produce pre-bounded batches,
//!   skipping Decompress, `FindBoundaries`, and Group entirely.
//! - **BGZF + synchronized**: Read → Decompress → `FindBoundaries` → Decode → Process →
//!   Serialize → Compress → Write. Decompress handles BGZF blocks, `FindBoundaries`
//!   locates record boundaries, Decode creates templates directly (bypassing Group).
//! - **Non-synchronized** (gzip or BGZF): Includes Group step for template grouping.
//!   Optionally includes `FindBoundaries` + Decode when `use_parallel_parse` is enabled.
//!
//! Inactive steps are excluded from scheduling — the conditionals are outside the
//! worker loop, not inside. Queues for skipped steps are allocated at minimum
//! capacity and never accessed.

use crossbeam_queue::ArrayQueue;
use noodles::sam::Header;
use parking_lot::Mutex;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::thread;

use crate::bgzf_reader::{
    BGZF_EOF, BGZF_FOOTER_SIZE, BGZF_HEADER_SIZE, decompress_block_slice_into, read_raw_blocks,
};
use crate::bgzf_writer::InlineBgzfCompressor;
use crate::grouper::{FastqGrouper, FastqRecord, FastqTemplate};
use crate::progress::ProgressTracker;
use crate::reorder_buffer::ReorderBuffer;
use libdeflater::Decompressor;

use super::base::{
    ActiveSteps, CompressedBlockBatch, HasCompressor, HasHeldBoundaries, HasHeldCompressed,
    HasHeldProcessed, HasHeldSerialized, HasWorkerCore, MemoryEstimate, MonitorableState,
    OutputPipelineQueues, OutputPipelineState, PipelineLifecycle, PipelineStats, PipelineStep,
    PipelineValidationError, ProcessPipelineState, ReorderBufferState, SerializePipelineState,
    SerializedBatch, StepContext, WorkerCoreState, WorkerStateCommon, WritePipelineState,
    finalize_pipeline, generic_worker_loop, handle_worker_panic, join_monitor_thread,
    join_worker_threads, run_monitor_loop, shared_try_step_compress,
};
use super::deadlock::{DeadlockConfig, DeadlockState, QueueSnapshot};
use super::scheduler::{BackpressureState, SchedulerStrategy};

// ============================================================================
// Multi-Stream Reader for FASTQ Input
// ============================================================================
//
// This module provides infrastructure for reading from multiple synchronized
// FASTQ input streams (R1, R2, I1, I2, etc.) in the unified pipeline.

// ============================================================================
// FASTQ Multi-Stream Types for Unified 7-Step Pipeline
// ============================================================================

/// Chunk from one FASTQ stream (may be compressed or decompressed).
///
/// Used in the unified 7-step pipeline for FASTQ → BAM conversion.
/// - For BGZF inputs: `is_decompressed = false`, Step 2 decompresses
/// - For Gzip/Plain inputs: `is_decompressed = true`, Step 2 passes through
#[derive(Debug)]
pub struct FastqStreamChunk {
    /// Which stream this chunk came from (0=R1, 1=R2, etc.)
    pub stream_idx: usize,
    /// Raw or decompressed bytes
    pub data: Vec<u8>,
    /// True if already decompressed (Gzip/Plain path)
    pub is_decompressed: bool,
}

/// Batch of chunks from multiple FASTQ streams (output of Step 1).
///
/// Each batch contains one chunk per stream (when available) and a serial
/// number for maintaining order through the pipeline.
#[derive(Debug, Default)]
pub struct FastqReadBatch {
    /// Chunks from each stream in this batch
    pub chunks: Vec<FastqStreamChunk>,
    /// Serial number for ordering
    pub serial: u64,
}

impl FastqReadBatch {
    /// Estimate heap memory usage of this batch.
    #[must_use]
    pub fn estimate_heap_size(&self) -> usize {
        // Vec capacity for chunks
        let chunks_capacity = self.chunks.capacity() * std::mem::size_of::<FastqStreamChunk>();
        // Data in each chunk
        let data_size: usize = self.chunks.iter().map(|c| c.data.capacity()).sum();
        chunks_capacity + data_size
    }
}

// ============================================================================
// Types for Parallel Parse Pipeline (moved from formats/fastq.rs)
// ============================================================================

/// Decompressed batch - data ready for boundary finding.
#[derive(Debug, Clone)]
pub struct FastqDecompressedBatch {
    /// Decompressed chunks from each stream.
    pub chunks: Vec<FastqDecompressedChunk>,
    /// Serial number for ordering.
    pub serial: u64,
}

/// A chunk of decompressed data from a single stream.
#[derive(Debug, Clone)]
pub struct FastqDecompressedChunk {
    /// Which stream this came from.
    pub stream_idx: usize,
    /// Decompressed data.
    pub data: Vec<u8>,
}

impl MemoryEstimate for FastqDecompressedBatch {
    fn estimate_heap_size(&self) -> usize {
        self.chunks.iter().map(|c| c.data.capacity()).sum()
    }
}

/// Batch with record boundaries identified.
/// Contains byte ranges for complete records only.
#[derive(Debug, Clone)]
pub struct FastqBoundaryBatch {
    /// Per-stream data with only complete records.
    pub streams: Vec<FastqStreamBoundaries>,
    /// Serial number for ordering.
    pub serial: u64,
}

/// Record boundaries for a single stream.
#[derive(Debug, Clone)]
pub struct FastqStreamBoundaries {
    /// Which stream this is for.
    pub stream_idx: usize,
    /// The data buffer (complete records only, leftover removed).
    pub data: Vec<u8>,
    /// Byte offsets where each record starts (including position 0).
    /// The end of record N is `offsets[N+1]` or `data.len()` for the last.
    pub offsets: Vec<usize>,
}

impl MemoryEstimate for FastqBoundaryBatch {
    fn estimate_heap_size(&self) -> usize {
        self.streams
            .iter()
            .map(|s| s.data.capacity() + s.offsets.capacity() * std::mem::size_of::<usize>())
            .sum()
    }
}

/// Parsed records for a single stream, carrying the stream identity.
#[derive(Debug, Clone)]
pub struct FastqParsedStream {
    /// Which stream this is for (0=R1, 1=R2, etc.).
    pub stream_idx: usize,
    /// Parsed records from this stream.
    pub records: Vec<FastqRecord>,
}

/// Parsed records batch - ready for grouping.
#[derive(Debug, Clone)]
pub struct FastqParsedBatch {
    /// Parsed records per stream, each carrying its `stream_idx`.
    pub streams: Vec<FastqParsedStream>,
    /// Serial number for ordering.
    pub serial: u64,
}

impl MemoryEstimate for FastqParsedBatch {
    fn estimate_heap_size(&self) -> usize {
        self.streams
            .iter()
            .map(|stream| {
                stream
                    .records
                    .iter()
                    .map(|r| r.name.capacity() + r.sequence.capacity() + r.quality.capacity())
                    .sum::<usize>()
            })
            .sum()
    }
}

/// State for finding FASTQ record boundaries across chunks.
///
/// Per-stream state for finding FASTQ record boundaries.
#[derive(Debug, Clone, Default)]
pub struct StreamBoundaryState {
    /// Leftover bytes (incomplete record from previous chunk).
    pub leftover: Vec<u8>,
    /// Reusable work buffer to reduce allocations.
    pub work_buffer: Vec<u8>,
}

/// Since chunks may split in the middle of a record, we need to:
/// 1. Save leftover bytes from incomplete records
/// 2. Prepend them to the next chunk's data
///
/// This struct uses per-stream locks to allow parallel processing of
/// different streams. For synchronized FASTQs, this eliminates lock
/// contention since each stream can be processed independently.
#[derive(Debug, Default)]
pub struct FastqBoundaryState {
    /// Per-stream state, each with its own lock for parallel access.
    pub stream_states: Vec<parking_lot::Mutex<StreamBoundaryState>>,
}

impl FastqBoundaryState {
    /// Create state for the given number of streams.
    #[must_use]
    pub fn new(num_streams: usize) -> Self {
        Self {
            stream_states: (0..num_streams)
                .map(|_| parking_lot::Mutex::new(StreamBoundaryState::default()))
                .collect(),
        }
    }
}

/// FASTQ format operations for the parallel parse pipeline.
///
/// Provides static methods for finding boundaries and parsing records.
pub struct FastqFormat;

impl FastqFormat {
    /// Find record boundaries in decompressed data.
    ///
    /// Each stream is processed independently with its own lock, allowing
    /// parallel processing of different streams. Within a stream, chunks
    /// must be processed in order (due to leftover handling).
    ///
    /// Uses reusable work buffers per stream to minimize allocations.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if boundary finding encounters invalid data.
    ///
    /// # Panics
    ///
    /// Panics if the boundary state is not initialized for all streams in the batch.
    pub fn find_boundaries(
        state: &FastqBoundaryState,
        batch: FastqDecompressedBatch,
    ) -> io::Result<FastqBoundaryBatch> {
        let max_stream = batch.chunks.iter().map(|c| c.stream_idx).max().unwrap_or(0);
        // Streams must be pre-allocated since we take &self (not &mut self)
        assert!(
            state.stream_states.len() > max_stream,
            "FastqBoundaryState not initialized for stream {max_stream}"
        );

        // Track which streams have chunks in this batch
        let mut streams_with_chunks = vec![false; state.stream_states.len()];
        let mut streams = Vec::with_capacity(state.stream_states.len());

        for chunk in batch.chunks {
            let stream_idx = chunk.stream_idx;
            streams_with_chunks[stream_idx] = true;

            // Lock only this stream's state - other streams can be processed in parallel
            let mut stream_state = state.stream_states[stream_idx].lock();

            // Use reusable work buffer to combine leftover with new data
            stream_state.work_buffer.clear();
            // Move leftover to work buffer (avoids double borrow)
            let leftover = std::mem::take(&mut stream_state.leftover);
            if !leftover.is_empty() {
                stream_state.work_buffer.extend_from_slice(&leftover);
            }
            stream_state.work_buffer.extend_from_slice(&chunk.data);
            // chunk.data is consumed here, freeing its allocation

            // Find complete FASTQ records
            let (data, offsets, leftover_start) =
                find_fastq_boundaries_inplace(&stream_state.work_buffer);

            // Save leftover for next chunk
            stream_state.leftover = stream_state.work_buffer[leftover_start..].to_vec();

            streams.push(FastqStreamBoundaries { stream_idx, data, offsets });
            // stream_state lock is dropped here, allowing other threads to access this stream
        }

        // For synchronized mode with multiple streams, process leftover from streams
        // that had no chunk in this batch. This handles the case where one stream's
        // reader reaches EOF before the other, but still has complete records in leftover.
        if state.stream_states.len() > 1 {
            for (stream_idx, &had_chunk) in streams_with_chunks.iter().enumerate() {
                if had_chunk {
                    continue; // Already processed above
                }

                let mut stream_state = state.stream_states[stream_idx].lock();
                if stream_state.leftover.is_empty() {
                    continue; // No leftover to process
                }

                // Process leftover as if it were a chunk
                stream_state.work_buffer.clear();
                let leftover = std::mem::take(&mut stream_state.leftover);
                stream_state.work_buffer.extend_from_slice(&leftover);

                let (data, offsets, leftover_start) =
                    find_fastq_boundaries_inplace(&stream_state.work_buffer);

                stream_state.leftover = stream_state.work_buffer[leftover_start..].to_vec();

                streams.push(FastqStreamBoundaries { stream_idx, data, offsets });
            }
        }

        // Align record counts for synchronized mode.
        // When multiple streams are processed together (e.g., R1/R2 paired FASTQs),
        // byte-aligned chunks may contain different numbers of complete records.
        // Move excess records back to leftover for the next batch.
        if streams.len() > 1 {
            // Find minimum record count across all streams
            // (offsets includes position 0, so record_count = offsets.len() - 1)
            let min_records =
                streams.iter().map(|s| s.offsets.len().saturating_sub(1)).min().unwrap_or(0);

            // Move excess records back to leftover for each stream
            for stream in &mut streams {
                let record_count = stream.offsets.len().saturating_sub(1);
                if record_count > min_records {
                    let excess_start = stream.offsets[min_records];

                    // Re-acquire lock for this stream
                    let mut stream_state = state.stream_states[stream.stream_idx].lock();

                    // Prepend excess bytes to existing leftover
                    // (existing leftover has incomplete record from end of this chunk)
                    let excess_bytes = stream.data[excess_start..].to_vec();
                    let incomplete_leftover = std::mem::take(&mut stream_state.leftover);
                    stream_state.leftover = excess_bytes;
                    stream_state.leftover.extend(incomplete_leftover);

                    // Truncate stream data and offsets
                    stream.data.truncate(excess_start);
                    stream.offsets.truncate(min_records + 1);
                }
            }
        }

        Ok(FastqBoundaryBatch { streams, serial: batch.serial })
    }

    /// Parse records from boundary data.
    ///
    /// This step is parallel - it's the key optimization for FASTQ.
    /// Given the byte offsets from `find_boundaries`, this constructs
    /// the actual record objects.
    /// # Errors
    ///
    /// Returns an I/O error if parsing any record fails.
    pub fn parse_records(batch: FastqBoundaryBatch) -> io::Result<FastqParsedBatch> {
        // This is the KEY PARALLEL STEP!
        // Parse records from boundary information.
        // Since each record is independent, this can be parallelized.
        //
        // We preserve `stream_idx` from the boundary batch so that downstream
        // consumers can correctly identify which stream each set of records
        // belongs to, regardless of the order in the Vec.

        let streams = batch
            .streams
            .into_iter()
            .map(|stream| {
                let records = parse_fastq_records_from_boundaries(&stream.data, &stream.offsets)?;
                Ok(FastqParsedStream { stream_idx: stream.stream_idx, records })
            })
            .collect::<io::Result<Vec<_>>>()?;

        Ok(FastqParsedBatch { streams, serial: batch.serial })
    }
}

// ============================================================================
// Boundary Finding Helper Functions
// ============================================================================

/// Find FASTQ record boundaries in data (allocation-reducing version).
///
/// Returns (`complete_data`, offsets, `leftover_start`).
/// - `complete_data`: Bytes containing only complete records (owned)
/// - offsets: Start positions of each record (including 0)
/// - `leftover_start`: Index where leftover begins in the original data
///
/// The caller should extract `data[leftover_start..]` for the leftover bytes.
fn find_fastq_boundaries_inplace(data: &[u8]) -> (Vec<u8>, Vec<usize>, usize) {
    if data.is_empty() {
        return (Vec::new(), vec![0], 0);
    }

    let mut offsets = vec![0];
    let mut pos = 0;

    // Scan for complete FASTQ records (4 lines each)
    while pos < data.len() {
        // Find end of this record (4 newlines)
        if let Some(record_end) = find_fastq_record_end(&data[pos..]) {
            pos += record_end;
            offsets.push(pos);
        } else {
            // Incomplete record - everything from pos onwards is leftover
            break;
        }
    }

    // Only allocate for the complete data (unavoidable - we return ownership)
    let complete_data = data[..pos].to_vec();

    (complete_data, offsets, pos)
}

/// Find the end of a complete FASTQ record starting at the given position.
///
/// Returns the number of bytes consumed (position after the last newline),
/// or None if the record is incomplete.
///
/// FASTQ format:
/// ```text
/// @name
/// ACGT...
/// +
/// IIII...
/// ```
fn find_fastq_record_end(data: &[u8]) -> Option<usize> {
    if data.is_empty() || data[0] != b'@' {
        return None;
    }

    // Find 4 newlines using SIMD-accelerated memchr
    let mut pos = 0;
    for _ in 0..4 {
        match memchr::memchr(b'\n', &data[pos..]) {
            Some(nl) => pos += nl + 1,
            None => return None,
        }
    }
    Some(pos)
}

/// Parse FASTQ records from boundary data.
///
/// This function takes pre-computed boundaries and constructs `FastqRecord` objects.
/// It's designed to be called in parallel for different chunks.
fn parse_fastq_records_from_boundaries(
    data: &[u8],
    offsets: &[usize],
) -> io::Result<Vec<FastqRecord>> {
    if offsets.len() <= 1 {
        return Ok(Vec::new());
    }

    let num_records = offsets.len() - 1;
    let mut records = Vec::with_capacity(num_records);

    for i in 0..num_records {
        let start = offsets[i];
        let end = offsets[i + 1];

        if start >= end || start >= data.len() {
            continue;
        }

        let record_data = &data[start..end.min(data.len())];
        let record = parse_single_fastq_record(record_data)?;
        records.push(record);
    }

    Ok(records)
}

/// Parse a single FASTQ record from bytes.
///
/// Expected format:
/// ```text
/// @name
/// ACGT...
/// +
/// IIII...
/// ```
fn parse_single_fastq_record(data: &[u8]) -> io::Result<FastqRecord> {
    // Line 1: @name
    if data.is_empty() || data[0] != b'@' {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "FASTQ record must start with @"));
    }

    let mut lines = data.split(|&b| b == b'\n');

    // Line 1: name (skip @)
    let name_line = lines
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing name line"))?;
    let name = name_line[1..].to_vec(); // Skip @

    // Line 2: sequence
    let sequence = lines
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing sequence line"))?
        .to_vec();

    // Line 3: + line (skip it)
    let plus_line =
        lines.next().ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing + line"))?;
    if plus_line.is_empty() || plus_line[0] != b'+' {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Third line must start with +"));
    }

    // Line 4: quality
    let quality = lines
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing quality line"))?
        .to_vec();

    // Validate lengths match
    if sequence.len() != quality.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Sequence length ({}) != quality length ({})", sequence.len(), quality.len()),
        ));
    }

    Ok(FastqRecord { name, sequence, quality })
}

/// Configuration for FASTQ 7-step pipeline.
///
/// This configuration controls the unified pipeline that works for both
/// BGZF-compressed and Gzip/Plain FASTQ inputs.
#[derive(Debug, Clone)]
#[allow(clippy::struct_excessive_bools)]
pub struct FastqPipelineConfig {
    /// Number of worker threads
    pub num_threads: usize,
    /// Capacity of inter-step queues
    pub queue_capacity: usize,
    /// BGZF compression level for output (0-12)
    pub compression_level: u32,
    /// Whether to collect pipeline statistics
    pub collect_stats: bool,
    /// True if inputs are BGZF (need decompression in Step 2)
    pub inputs_are_bgzf: bool,
    /// Number of templates to batch before processing.
    /// ~400 templates × ~150 bytes = ~60KB per batch (optimal for BGZF).
    pub batch_size: usize,
    /// Scheduler strategy for thread work assignment.
    pub scheduler_strategy: SchedulerStrategy,
    /// Memory limit for template queue in bytes (0 = no limit, but backpressure
    /// is still applied at 512MB for optimal performance).
    pub queue_memory_limit: u64,
    /// Deadlock detection timeout in seconds (0 = disabled).
    pub deadlock_timeout_secs: u64,
    /// Whether automatic deadlock recovery is enabled.
    pub deadlock_recover_enabled: bool,
    /// Whether to use parallel Parse step (the key t8 scaling optimization).
    /// When true, parsing is done in parallel before the Group step.
    /// When false, parsing is done under the Group step's lock (original behavior).
    pub use_parallel_parse: bool,
    /// Whether input FASTQs are synchronized (records at same position match).
    /// When true, skips name validation in Group step (done in Process instead).
    /// This eliminates lock contention for synchronized FASTQs.
    pub synchronized: bool,
    /// Shared statistics instance for external memory monitoring.
    /// When provided, the pipeline will use this instead of creating its own.
    pub shared_stats: Option<Arc<PipelineStats>>,
}

impl FastqPipelineConfig {
    /// Create a new configuration with the given thread count.
    ///
    /// # Arguments
    /// * `num_threads` - Number of worker threads (minimum 1)
    /// * `inputs_are_bgzf` - True if inputs are BGZF-compressed
    /// * `compression_level` - BGZF compression level (1-12)
    #[must_use]
    pub fn new(num_threads: usize, inputs_are_bgzf: bool, compression_level: u32) -> Self {
        let threads = num_threads.max(1);
        // Scale queue capacity with threads to limit memory usage
        // 128 per thread, capped at 1024 (original default for 8+ threads)
        let queue_capacity = (threads * 128).min(1024);
        Self {
            num_threads: threads,
            queue_capacity,
            compression_level,
            collect_stats: false,
            inputs_are_bgzf,
            batch_size: 400, // ~400 templates × ~150 bytes = ~60KB per batch
            scheduler_strategy: SchedulerStrategy::default(),
            queue_memory_limit: 4 * 1024 * 1024 * 1024, // 4GB default
            deadlock_timeout_secs: 10,                  // Default 10s
            deadlock_recover_enabled: false,
            use_parallel_parse: false, // Disabled - overhead hurts t4 performance
            synchronized: false,       // Not synchronized by default
            shared_stats: None,        // No shared stats by default
        }
    }

    /// Enable or disable statistics collection.
    #[must_use]
    pub fn with_stats(mut self, enabled: bool) -> Self {
        self.collect_stats = enabled;
        self
    }

    /// Set the compression level for output.
    #[must_use]
    pub fn with_compression_level(mut self, level: u32) -> Self {
        self.compression_level = level.min(12);
        self
    }

    /// Set the queue capacity.
    #[must_use]
    pub fn with_queue_capacity(mut self, capacity: usize) -> Self {
        self.queue_capacity = capacity.max(8);
        self
    }

    /// Set the batch size for template batching.
    #[must_use]
    pub fn with_batch_size(mut self, size: usize) -> Self {
        self.batch_size = size.max(1);
        self
    }

    /// Set the scheduler strategy.
    #[must_use]
    pub fn with_scheduler_strategy(mut self, strategy: SchedulerStrategy) -> Self {
        self.scheduler_strategy = strategy;
        self
    }

    /// Set the memory limit for template queue in bytes.
    /// Use 0 to disable the limit (backpressure still applied at 512MB).
    #[must_use]
    pub fn with_queue_memory_limit(mut self, limit_bytes: u64) -> Self {
        self.queue_memory_limit = limit_bytes;
        self
    }

    /// Set the deadlock detection timeout in seconds.
    /// Use 0 to disable deadlock detection.
    #[must_use]
    pub fn with_deadlock_timeout(mut self, timeout_secs: u64) -> Self {
        self.deadlock_timeout_secs = timeout_secs;
        self
    }

    /// Enable or disable automatic deadlock recovery.
    #[must_use]
    pub fn with_deadlock_recovery(mut self, enabled: bool) -> Self {
        self.deadlock_recover_enabled = enabled;
        self
    }

    /// Enable synchronized mode for FASTQs known to be in sync.
    ///
    /// In synchronized mode, records are zipped by position without name validation
    /// in the Group step. Name validation is done in the Process step instead.
    /// This eliminates lock contention for synchronized FASTQs.
    ///
    /// Enabling synchronized mode automatically enables parallel parse, as the
    /// synchronized optimization requires the parallel parse pipeline.
    ///
    /// Note: Synchronized mode is automatically disabled for single-threaded execution
    /// because the optimization relies on multiple threads to avoid deadlock. With one
    /// thread, there is no lock contention to optimize away.
    #[must_use]
    pub fn with_synchronized(mut self, enabled: bool) -> Self {
        // Synchronized mode uses spin-wait which requires multiple threads to drain
        // queues. With a single thread, this causes deadlock. Additionally, the
        // optimization eliminates lock contention which doesn't exist with one thread.
        self.synchronized = enabled && self.num_threads > 1;
        // Synchronized mode requires parallel parse to work
        if self.synchronized {
            self.use_parallel_parse = true;
        }
        self
    }

    /// Compute the active pipeline steps based on input type and mode.
    ///
    /// - Gzip+synchronized: Read, Decode, Process, Serialize, Compress, Write
    /// - BGZF+synchronized: Read, Decompress, `FindBoundaries`, Decode, Process,
    ///   Serialize, Compress, Write
    /// - BGZF (not synchronized): depends on `use_parallel_parse`
    #[must_use]
    pub fn active_steps(&self) -> ActiveSteps {
        // Synchronized mode requires parallel parse; enforce the invariant here
        // so callers that construct config with pub fields directly get correct steps.
        let use_parallel_parse = self.use_parallel_parse || self.synchronized;

        let mut steps = vec![PipelineStep::Read];

        // RecordCount reader (gzip+synchronized) bypasses q1_raw entirely,
        // so Decompress is not needed. All other paths need it.
        let is_record_count = !self.inputs_are_bgzf && self.synchronized;
        if !is_record_count {
            steps.push(PipelineStep::Decompress);
        }

        // `FindBoundaries` scans decompressed data for record boundaries.
        // Only needed for BGZF when parallel parse or synchronized mode is active.
        if self.inputs_are_bgzf && use_parallel_parse {
            steps.push(PipelineStep::FindBoundaries);
        }

        // Decode (Parse) parses FASTQ records from boundary data.
        // Only active when parallel parse is enabled (or forced by synchronized mode).
        if use_parallel_parse {
            steps.push(PipelineStep::Decode);
        }

        // Group collects parsed records into templates.
        // Skipped in synchronized mode (Decode creates templates directly).
        if !self.synchronized {
            steps.push(PipelineStep::Group);
        }

        steps.extend([
            PipelineStep::Process,
            PipelineStep::Serialize,
            PipelineStep::Compress,
            PipelineStep::Write,
        ]);
        ActiveSteps::new(&steps)
    }
}

// ============================================================================
// Phase 2: MultiBgzfBlockReader - Step 1 for BGZF Path
// ============================================================================

/// Reads raw BGZF blocks from multiple FASTQ files.
///
/// Outputs `FastqReadBatch` with `is_decompressed=false`, indicating that
/// Step 2 should decompress the data.
///
/// This is Step 1 of the 7-step pipeline for BGZF-compressed inputs.
pub struct MultiBgzfBlockReader {
    /// Buffered readers for each input file
    readers: Vec<BufReader<File>>,
    /// Next serial number for batches
    next_serial: u64,
    /// EOF flags for each stream
    eof_flags: Vec<bool>,
    /// Number of BGZF blocks to read per stream per batch
    blocks_per_batch: usize,
}

impl MultiBgzfBlockReader {
    /// Create a new reader for the given FASTQ file paths.
    ///
    /// # Arguments
    /// * `paths` - Paths to BGZF-compressed FASTQ files (R1, R2, etc.)
    ///
    /// # Errors
    /// Returns an error if any file cannot be opened.
    pub fn new(paths: &[PathBuf]) -> io::Result<Self> {
        let readers = paths
            .iter()
            .map(|p| {
                let file = File::open(p)?;
                Ok(BufReader::with_capacity(256 * 1024, file))
            })
            .collect::<io::Result<Vec<_>>>()?;
        let num_streams = readers.len();
        Ok(Self {
            readers,
            next_serial: 0,
            eof_flags: vec![false; num_streams],
            blocks_per_batch: 4, // Read 4 BGZF blocks per stream per batch
        })
    }

    /// Get the number of input streams.
    #[must_use]
    pub fn num_streams(&self) -> usize {
        self.readers.len()
    }

    /// Check if all streams have reached EOF.
    #[must_use]
    pub fn is_done(&self) -> bool {
        self.eof_flags.iter().all(|&eof| eof)
    }

    /// Read next batch of raw BGZF blocks from all streams.
    ///
    /// Returns `Ok(None)` when all streams have reached EOF.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if reading BGZF blocks fails.
    pub fn read_next_batch(&mut self) -> io::Result<Option<FastqReadBatch>> {
        if self.is_done() {
            return Ok(None);
        }

        let mut batch = FastqReadBatch {
            chunks: Vec::with_capacity(self.readers.len()),
            serial: self.next_serial,
        };
        self.next_serial += 1;

        for (stream_idx, reader) in self.readers.iter_mut().enumerate() {
            if self.eof_flags[stream_idx] {
                continue;
            }

            // Read multiple BGZF blocks for this stream
            let blocks = read_raw_blocks(reader, self.blocks_per_batch)?;

            if blocks.is_empty() {
                self.eof_flags[stream_idx] = true;
                continue;
            }

            // Concatenate raw block data
            let total_size: usize = blocks.iter().map(|b| b.data.len()).sum();
            let mut raw_data = Vec::with_capacity(total_size);
            for block in blocks {
                raw_data.extend_from_slice(&block.data);
            }

            batch.chunks.push(FastqStreamChunk {
                stream_idx,
                data: raw_data,
                is_decompressed: false, // BGZF: needs decompression in Step 2
            });
        }

        if batch.chunks.is_empty() { Ok(None) } else { Ok(Some(batch)) }
    }

    /// Set the number of BGZF blocks to read per stream per batch.
    pub fn set_blocks_per_batch(&mut self, blocks: usize) {
        self.blocks_per_batch = blocks.max(1);
    }
}

// ============================================================================
// Phase 3: MultiDecompressedReader - Step 1 for Gzip/Plain Path
// ============================================================================

/// Reads already-decompressed chunks from multiple FASTQ files.
///
/// Outputs `FastqReadBatch` with `is_decompressed=true`, indicating that
/// Step 2 should pass through the data without decompression.
///
/// This is Step 1 of the 7-step pipeline for Gzip or Plain inputs,
/// where decompression is handled by fgoxide upfront.
pub struct MultiDecompressedReader<R: BufRead + Send> {
    /// Readers for each input stream
    readers: Vec<R>,
    /// Next serial number for batches
    next_serial: u64,
    /// EOF flags for each stream
    eof_flags: Vec<bool>,
    /// Size of chunks to read from each stream
    chunk_size: usize,
}

impl<R: BufRead + Send> MultiDecompressedReader<R> {
    /// Create a new reader from pre-opened decompressed readers.
    ///
    /// # Arguments
    /// * `readers` - Already-opened readers (e.g., from fgoxide)
    #[must_use]
    pub fn new(readers: Vec<R>) -> Self {
        let num_streams = readers.len();
        Self {
            readers,
            next_serial: 0,
            eof_flags: vec![false; num_streams],
            chunk_size: 64 * 1024, // 64KB chunks
        }
    }

    /// Get the number of input streams.
    #[must_use]
    pub fn num_streams(&self) -> usize {
        self.readers.len()
    }

    /// Check if all streams have reached EOF.
    #[must_use]
    pub fn is_done(&self) -> bool {
        self.eof_flags.iter().all(|&eof| eof)
    }

    /// Read next batch of decompressed chunks from all streams.
    ///
    /// Returns `Ok(None)` when all streams have reached EOF.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if reading from any stream fails.
    pub fn read_next_batch(&mut self) -> io::Result<Option<FastqReadBatch>> {
        log::trace!(
            "MultiDecompressedReader::read_next_batch: is_done={}, eof_flags={:?}",
            self.is_done(),
            self.eof_flags
        );
        if self.is_done() {
            return Ok(None);
        }

        let mut batch = FastqReadBatch {
            chunks: Vec::with_capacity(self.readers.len()),
            serial: self.next_serial,
        };
        self.next_serial += 1;

        for (stream_idx, reader) in self.readers.iter_mut().enumerate() {
            if self.eof_flags[stream_idx] {
                continue;
            }

            let mut data = vec![0u8; self.chunk_size];
            let mut total_read = 0;

            // Fill the buffer
            while total_read < self.chunk_size {
                match reader.read(&mut data[total_read..])? {
                    0 => {
                        self.eof_flags[stream_idx] = true;
                        break;
                    }
                    n => total_read += n,
                }
            }

            if total_read > 0 {
                data.truncate(total_read);
                batch.chunks.push(FastqStreamChunk {
                    stream_idx,
                    data,
                    is_decompressed: true, // Gzip/Plain: already decompressed
                });
            }
        }

        log::trace!(
            "MultiDecompressedReader::read_next_batch: chunks={}, serial={}",
            batch.chunks.len(),
            batch.serial
        );
        if batch.chunks.is_empty() { Ok(None) } else { Ok(Some(batch)) }
    }

    /// Set the chunk size for reading.
    pub fn set_chunk_size(&mut self, size: usize) {
        self.chunk_size = size.max(4096);
    }
}

// ============================================================================
// MultiRecordCountReader - Record-count-based reader for gzip FASTQs
// ============================================================================

/// Default number of FASTQ records per stream per batch.
/// ~200 records × ~300 bytes ≈ 60KB, comparable to the 64KB chunk size.
const DEFAULT_RECORDS_PER_BATCH: usize = 200;

/// Reads a fixed number of complete FASTQ records from each stream per batch.
///
/// Unlike `MultiDecompressedReader` which reads fixed byte-size chunks,
/// this reader guarantees equal record counts across all streams in each batch.
/// This eliminates the O(N²) regression caused by `find_boundaries()` having to
/// re-scan growing leftover buffers when R1/R2 records have different lengths.
///
/// Produces `FastqBoundaryBatch` directly, skipping the `Decompress` and
/// `FindBoundaries` steps entirely.
pub struct MultiRecordCountReader<R: BufRead + Send> {
    /// Readers for each input stream.
    readers: Vec<R>,
    /// Next serial number for batches.
    next_serial: u64,
    /// EOF flags for each stream.
    eof_flags: Vec<bool>,
    /// Number of FASTQ records to read per stream per batch.
    records_per_batch: usize,
}

impl<R: BufRead + Send> MultiRecordCountReader<R> {
    /// Create a new record-count reader from pre-opened decompressed readers.
    #[must_use]
    pub fn new(readers: Vec<R>) -> Self {
        let num_streams = readers.len();
        Self {
            readers,
            next_serial: 0,
            eof_flags: vec![false; num_streams],
            records_per_batch: DEFAULT_RECORDS_PER_BATCH,
        }
    }

    /// Get the number of input streams.
    #[must_use]
    pub fn num_streams(&self) -> usize {
        self.readers.len()
    }

    /// Check if all streams have reached EOF.
    #[must_use]
    pub fn is_done(&self) -> bool {
        self.eof_flags.iter().all(|&eof| eof)
    }

    /// Set the number of records per batch.
    pub fn set_records_per_batch(&mut self, n: usize) {
        self.records_per_batch = n.max(1);
    }

    /// Read the next batch of FASTQ records from all streams.
    ///
    /// Reads exactly `records_per_batch` complete FASTQ records from each
    /// non-EOF stream. At EOF, aligns to the minimum record count across
    /// streams so templates can be formed.
    ///
    /// Returns `Ok(None)` when all streams have reached EOF.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if reading from any stream fails.
    pub fn read_next_batch(&mut self) -> io::Result<Option<FastqBoundaryBatch>> {
        if self.is_done() {
            return Ok(None);
        }

        #[allow(clippy::type_complexity)]
        let mut stream_results: Vec<Option<(Vec<u8>, Vec<usize>, bool)>> =
            Vec::with_capacity(self.readers.len());

        for (stream_idx, reader) in self.readers.iter_mut().enumerate() {
            if self.eof_flags[stream_idx] {
                stream_results.push(None);
                continue;
            }

            let (data, offsets, at_eof) = read_n_fastq_records(reader, self.records_per_batch)?;

            if at_eof {
                self.eof_flags[stream_idx] = true;
            }

            if offsets.len() <= 1 {
                // No complete records read (empty or EOF)
                stream_results.push(None);
            } else {
                stream_results.push(Some((data, offsets, at_eof)));
            }
        }

        // Find min record count across active streams for alignment at EOF
        let min_records = stream_results
            .iter()
            .filter_map(|r| r.as_ref().map(|(_, offsets, _)| offsets.len() - 1))
            .min();

        let Some(min_records) = min_records else {
            // All streams exhausted
            return Ok(None);
        };

        // If any stream is exhausted while others still have data, stop to avoid
        // producing batches with missing streams (which would drop mates).
        if stream_results.iter().any(Option::is_none) {
            for flag in &mut self.eof_flags {
                *flag = true;
            }
            return Ok(None);
        }

        if min_records == 0 {
            return Ok(None);
        }

        // Build FastqBoundaryBatch, truncating to min_records at EOF boundaries
        let mut streams = Vec::with_capacity(self.readers.len());
        for (stream_idx, result) in stream_results.into_iter().enumerate() {
            let Some((data, offsets, _at_eof)) = result else {
                continue;
            };

            let record_count = offsets.len() - 1;
            let (final_data, final_offsets) = if record_count > min_records {
                // Truncate to min_records. Excess records are permanently discarded
                // (reader already advanced). This only occurs at EOF, where excess
                // records have no mate and cannot form valid templates.
                let truncated_offsets: Vec<usize> = offsets[..=min_records].to_vec();
                let truncated_data = data[..offsets[min_records]].to_vec();
                (truncated_data, truncated_offsets)
            } else {
                (data, offsets)
            };

            streams.push(FastqStreamBoundaries {
                stream_idx,
                data: final_data,
                offsets: final_offsets,
            });
        }

        if streams.is_empty() {
            return Ok(None);
        }

        let serial = self.next_serial;
        self.next_serial += 1;

        Ok(Some(FastqBoundaryBatch { streams, serial }))
    }
}

/// Read exactly `n` complete FASTQ records from a buffered reader.
///
/// Uses `BufRead::read_until` (which delegates to `memchr` internally) for
/// SIMD-accelerated newline scanning. Each FASTQ record
/// consists of exactly 4 lines (name, sequence, plus, quality).
///
/// Returns `(data, offsets, at_eof)` where:
/// - `data`: concatenated bytes of all complete records
/// - `offsets`: byte offset of each record start (includes 0, so len = records + 1)
/// - `at_eof`: true if the reader reached EOF
fn read_n_fastq_records<R: BufRead>(
    reader: &mut R,
    n: usize,
) -> io::Result<(Vec<u8>, Vec<usize>, bool)> {
    // Pre-allocate: ~300 bytes per record is typical
    let mut data = Vec::with_capacity(n * 300);
    let mut offsets = Vec::with_capacity(n + 1);
    offsets.push(0);
    let mut at_eof = false;

    for _ in 0..n {
        let record_start = data.len();

        // Read 4 lines per FASTQ record
        let mut lines_read = 0;
        for _ in 0..4 {
            let before = data.len();
            let bytes_read = reader.read_until(b'\n', &mut data)?;
            if bytes_read == 0 {
                // EOF mid-record: truncate back to record start
                data.truncate(record_start);
                at_eof = true;
                break;
            }
            lines_read += 1;

            // If the line doesn't end with \n (EOF without trailing newline),
            // add one for consistency
            if data[data.len() - 1] != b'\n' {
                data.push(b'\n');
                at_eof = true;
            }

            // Validate first line starts with '@' (name line)
            if lines_read == 1 && data[before] != b'@' {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "Expected FASTQ record to start with '@', got '{}'",
                        data[before] as char
                    ),
                ));
            }
        }

        if lines_read < 4 {
            // Incomplete record at EOF
            break;
        }

        offsets.push(data.len());

        if at_eof {
            break;
        }
    }

    Ok((data, offsets, at_eof))
}

// ============================================================================
// Phase 4: Step 2 - Decompress or Passthrough
// ============================================================================

/// Decompress a batch of FASTQ chunks (Step 2 of the 7-step pipeline).
///
/// This function handles both BGZF and Gzip/Plain paths:
/// - If `chunk.is_decompressed == true`: passthrough (no-op)
/// - If `chunk.is_decompressed == false`: decompress BGZF blocks
///
/// # Arguments
/// * `batch` - The batch of chunks to process
/// * `decompressor` - A reusable libdeflater decompressor
///
/// # Returns
/// A new batch with all chunks marked as decompressed
///
/// # Errors
///
/// Returns an I/O error if BGZF decompression fails.
pub fn decompress_fastq_batch(
    batch: FastqReadBatch,
    decompressor: &mut Decompressor,
) -> io::Result<FastqReadBatch> {
    let mut output =
        FastqReadBatch { chunks: Vec::with_capacity(batch.chunks.len()), serial: batch.serial };

    for chunk in batch.chunks {
        let decompressed_data = if chunk.is_decompressed {
            // Gzip/Plain path: passthrough
            chunk.data
        } else {
            // BGZF path: decompress all blocks in the chunk
            decompress_bgzf_chunk(&chunk.data, decompressor)?
        };

        output.chunks.push(FastqStreamChunk {
            stream_idx: chunk.stream_idx,
            data: decompressed_data,
            is_decompressed: true,
        });
    }

    Ok(output)
}

/// Estimate total uncompressed size by summing ISIZE fields from all BGZF blocks.
///
/// The ISIZE field is the last 4 bytes of each BGZF block and contains the
/// uncompressed size of that block. By summing these, we can pre-allocate the
/// result buffer to avoid reallocations during decompression.
fn estimate_uncompressed_size(raw_data: &[u8]) -> usize {
    let mut total = 0;
    let mut offset = 0;

    while offset + BGZF_HEADER_SIZE + BGZF_FOOTER_SIZE <= raw_data.len() {
        let bsize = u16::from_le_bytes([raw_data[offset + 16], raw_data[offset + 17]]) as usize;
        let block_size = bsize + 1;

        if offset + block_size > raw_data.len() {
            break;
        }

        // ISIZE is last 4 bytes of block (uncompressed size mod 2^32)
        let isize_offset = offset + block_size - 4;
        if isize_offset + 4 <= raw_data.len() {
            let isize = u32::from_le_bytes([
                raw_data[isize_offset],
                raw_data[isize_offset + 1],
                raw_data[isize_offset + 2],
                raw_data[isize_offset + 3],
            ]) as usize;
            total += isize;
        }
        offset += block_size;
    }
    total
}

/// Decompress concatenated BGZF blocks from a single chunk.
///
/// The chunk contains multiple BGZF blocks concatenated together.
/// This function parses each block and decompresses it.
fn decompress_bgzf_chunk(raw_data: &[u8], decompressor: &mut Decompressor) -> io::Result<Vec<u8>> {
    // Pre-allocate based on BGZF ISIZE fields to avoid reallocations
    let estimated_size = estimate_uncompressed_size(raw_data);
    let mut result = Vec::with_capacity(estimated_size);
    let mut offset = 0;

    while offset < raw_data.len() {
        // Need at least header size to parse block
        if offset + BGZF_HEADER_SIZE > raw_data.len() {
            break; // Incomplete block at end
        }

        // Parse BSIZE from header (bytes 16-17, little-endian)
        let bsize = u16::from_le_bytes([raw_data[offset + 16], raw_data[offset + 17]]) as usize;
        let block_size = bsize + 1;

        // Validate block size
        if offset + block_size > raw_data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "BGZF block extends past chunk: offset={}, block_size={}, chunk_len={}",
                    offset,
                    block_size,
                    raw_data.len()
                ),
            ));
        }

        // Decompress directly from slice - no allocation for RawBgzfBlock
        decompress_block_slice_into(
            &raw_data[offset..offset + block_size],
            decompressor,
            &mut result,
        )?;

        offset += block_size;
    }

    Ok(result)
}

// ============================================================================
// Phase 5: FastqMultiStreamGrouper - Step 3 Grouping
// ============================================================================

/// Grouper wrapper for multi-stream FASTQ that receives `FastqReadBatch`.
///
/// This wraps the existing `FastqGrouper` from grouper.rs and provides
/// a simpler interface for the 7-step pipeline.
///
/// Note: For synchronized FASTQs, the grouper is bypassed entirely -
/// templates are created directly in the Parse step. This struct is only
/// used for non-synchronized mode.
pub struct FastqMultiStreamGrouper {
    /// The underlying grouper
    inner: FastqGrouper,
}

impl FastqMultiStreamGrouper {
    /// Create a new grouper for the specified number of input streams.
    #[must_use]
    pub fn new(num_streams: usize) -> Self {
        Self { inner: FastqGrouper::new(num_streams) }
    }

    /// Add a decompressed batch and return all complete templates.
    ///
    /// This feeds each chunk in the batch to the appropriate stream
    /// and then drains all complete templates.
    ///
    /// # Arguments
    /// * `batch` - A batch with all chunks marked as decompressed
    ///
    /// # Returns
    /// All templates that are complete (have records from all streams)
    ///
    /// # Errors
    ///
    /// Returns an I/O error if parsing or grouping fails.
    pub fn add_batch(&mut self, batch: FastqReadBatch) -> io::Result<Vec<FastqTemplate>> {
        // Feed each chunk to the appropriate stream (parses and accumulates)
        for chunk in batch.chunks {
            debug_assert!(chunk.is_decompressed, "Step 3 expects decompressed data, got raw BGZF");
            self.inner.add_bytes_for_stream(chunk.stream_idx, &chunk.data)?;
        }

        // Drain all complete templates (records matched across all streams)
        self.inner.drain_complete_templates()
    }

    /// Add pre-parsed records from a `FastqParsedBatch` and return all complete templates.
    ///
    /// This is the key method for the parallel Parse optimization - it receives
    /// already-parsed records so no parsing happens under the grouper lock.
    ///
    /// # Arguments
    /// * `streams` - Vector of parsed records per stream
    ///
    /// # Returns
    /// All templates that are complete (have records from all streams)
    ///
    /// # Errors
    ///
    /// Returns an I/O error if grouping fails.
    pub fn add_parsed_batch(
        &mut self,
        streams: Vec<FastqParsedStream>,
    ) -> io::Result<Vec<FastqTemplate>> {
        // Feed to grouper using the explicit stream_idx from each stream,
        // not positional index — this is robust to any ordering.
        for stream in streams {
            self.inner.add_records_for_stream(stream.stream_idx, stream.records)?;
        }

        // Drain all complete templates (records matched across all streams)
        self.inner.drain_complete_templates()
    }

    /// Check if there are pending records or leftover bytes.
    #[must_use]
    pub fn has_pending(&self) -> bool {
        self.inner.has_pending()
    }

    /// Finish processing and return any remaining template.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if the grouper fails to finalize.
    pub fn finish(&mut self) -> io::Result<Option<FastqTemplate>> {
        self.inner.finish()
    }
}

// ============================================================================
// Phase 6 & 7: FASTQ Pipeline State and Entry Point
// ============================================================================

/// Reader type for the FASTQ pipeline.
///
/// This enum abstracts over the two reader types (BGZF and Gzip/Plain).
pub enum FastqReader<R: BufRead + Send> {
    /// BGZF-compressed inputs - reads raw blocks for parallel decompression
    Bgzf(MultiBgzfBlockReader),
    /// Gzip/Plain inputs - reads already-decompressed chunks
    Decompressed(MultiDecompressedReader<R>),
    /// Record-count-based reader - produces `FastqBoundaryBatch` directly.
    /// Used for gzip/plain inputs to skip `Decompress` + `FindBoundaries` steps.
    RecordCount(MultiRecordCountReader<R>),
}

impl<R: BufRead + Send> FastqReader<R> {
    /// Read the next batch (for Bgzf and Decompressed variants).
    ///
    /// # Errors
    ///
    /// Returns an I/O error if reading fails or if called on a `RecordCount` reader.
    pub fn read_next_batch(&mut self) -> io::Result<Option<FastqReadBatch>> {
        match self {
            FastqReader::Bgzf(r) => r.read_next_batch(),
            FastqReader::Decompressed(r) => r.read_next_batch(),
            FastqReader::RecordCount(_) => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "RecordCount reader does not produce FastqReadBatch; use read_next_boundary_batch()",
            )),
        }
    }

    /// Read the next boundary batch (for `RecordCount` variant only).
    ///
    /// # Errors
    ///
    /// Returns an I/O error if reading fails or if called on a non-RecordCount reader.
    pub fn read_next_boundary_batch(&mut self) -> io::Result<Option<FastqBoundaryBatch>> {
        match self {
            FastqReader::RecordCount(r) => r.read_next_batch(),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "read_next_boundary_batch() only works with RecordCount reader",
            )),
        }
    }

    /// Check if all streams have reached EOF.
    #[must_use]
    pub fn is_done(&self) -> bool {
        match self {
            FastqReader::Bgzf(r) => r.is_done(),
            FastqReader::Decompressed(r) => r.is_done(),
            FastqReader::RecordCount(r) => r.is_done(),
        }
    }

    /// Get the number of streams.
    #[must_use]
    pub fn num_streams(&self) -> usize {
        match self {
            FastqReader::Bgzf(r) => r.num_streams(),
            FastqReader::Decompressed(r) => r.num_streams(),
            FastqReader::RecordCount(r) => r.num_streams(),
        }
    }

    /// Returns true if this is a record-count reader.
    #[must_use]
    pub fn is_record_count(&self) -> bool {
        matches!(self, FastqReader::RecordCount(_))
    }
}

/// Default serialization buffer capacity.
/// Sized for typical batch: 400 templates × ~300 bytes = 120KB, plus headroom.
const SERIALIZATION_BUFFER_CAPACITY: usize = 256 * 1024; // 256KB

/// Worker state for FASTQ pipeline threads.
///
/// Generic over `P` (the processed record type) to support holding processed batches.
pub struct FastqWorkerState<P: Send> {
    /// Core worker state (compressor, scheduler, serialization buffer, backoff).
    pub core: WorkerCoreState,
    /// Decompressor for BGZF blocks.
    pub decompressor: Decompressor,

    // =========================================================================
    // Held Items for Non-Blocking Steps
    // =========================================================================
    // These fields allow step functions to return immediately when output queues
    // are full, instead of blocking in push_with_backoff. This prevents deadlock
    // by ensuring threads can always return to try other steps (especially Write).
    /// Held raw batch from Read step (couldn't push to `q1_raw`).
    pub held_raw: Option<(u64, FastqReadBatch)>,
    /// Held decompressed batch from Decompress step (couldn't push to `q2_decompressed`).
    /// Includes `heap_size` for memory tracking.
    pub held_decompressed: Option<(u64, FastqReadBatch, usize)>,
    /// Held boundary batch from `FindBoundaries` step (couldn't push to `q2_5_boundaries`).
    pub held_boundaries: Option<(u64, FastqBoundaryBatch)>,
    /// Held processed batch from Process step (couldn't push to `q4_processed`).
    /// Includes `heap_size` for memory tracking.
    pub held_processed: Option<(u64, Vec<P>, usize)>,
    /// Held serialized batch from Serialize step (couldn't push to `q5_serialized`).
    /// Includes `heap_size` for memory tracking.
    pub held_serialized: Option<(u64, SerializedBatch, usize)>,
    /// Held compressed batch from Compress step (couldn't push to `q6_compressed`).
    /// Includes `heap_size` for memory tracking.
    pub held_compressed: Option<(u64, CompressedBlockBatch, usize)>,
}

impl<P: Send> FastqWorkerState<P> {
    /// Create new worker state with the given compression level, thread info, and scheduler strategy.
    #[must_use]
    pub fn new(
        compression_level: u32,
        thread_id: usize,
        num_threads: usize,
        scheduler_strategy: SchedulerStrategy,
        active_steps: ActiveSteps,
    ) -> Self {
        Self {
            core: WorkerCoreState::new(
                compression_level,
                thread_id,
                num_threads,
                scheduler_strategy,
                active_steps,
            ),
            decompressor: Decompressor::new(),
            held_raw: None,
            held_decompressed: None,
            held_boundaries: None,
            held_processed: None,
            held_serialized: None,
            held_compressed: None,
        }
    }

    /// Returns true if any held item fields are Some.
    ///
    /// Used to ensure all held items are flushed before pipeline completion.
    #[must_use]
    pub fn has_any_held_items(&self) -> bool {
        self.held_raw.is_some()
            || self.held_decompressed.is_some()
            || self.held_boundaries.is_some()
            || self.held_processed.is_some()
            || self.held_serialized.is_some()
            || self.held_compressed.is_some()
    }

    /// Clear all held items (for cleanup/error handling).
    pub fn clear_held_items(&mut self) {
        self.held_raw = None;
        self.held_decompressed = None;
        self.held_boundaries = None;
        self.held_processed = None;
        self.held_serialized = None;
        self.held_compressed = None;
    }
}

impl<P: Send> HasCompressor for FastqWorkerState<P> {
    fn compressor_mut(&mut self) -> &mut InlineBgzfCompressor {
        &mut self.core.compressor
    }
}

impl<P: Send> HasHeldCompressed for FastqWorkerState<P> {
    fn held_compressed_mut(&mut self) -> &mut Option<(u64, CompressedBlockBatch, usize)> {
        &mut self.held_compressed
    }
}

impl<P: Send> HasHeldBoundaries<FastqBoundaryBatch> for FastqWorkerState<P> {
    fn held_boundaries_mut(&mut self) -> &mut Option<(u64, FastqBoundaryBatch)> {
        &mut self.held_boundaries
    }
}

impl<P: Send> HasHeldProcessed<P> for FastqWorkerState<P> {
    fn held_processed_mut(&mut self) -> &mut Option<(u64, Vec<P>, usize)> {
        &mut self.held_processed
    }
}

impl<P: Send> HasHeldSerialized for FastqWorkerState<P> {
    fn held_serialized_mut(&mut self) -> &mut Option<(u64, SerializedBatch, usize)> {
        &mut self.held_serialized
    }

    fn serialization_buffer_mut(&mut self) -> &mut Vec<u8> {
        &mut self.core.serialization_buffer
    }

    fn serialization_buffer_capacity(&self) -> usize {
        SERIALIZATION_BUFFER_CAPACITY // 256KB for FASTQ
    }
}

impl<P: Send> WorkerStateCommon for FastqWorkerState<P> {
    fn has_any_held_items(&self) -> bool {
        FastqWorkerState::has_any_held_items(self)
    }

    fn clear_held_items(&mut self) {
        FastqWorkerState::clear_held_items(self);
    }
}

impl<P: Send> HasWorkerCore for FastqWorkerState<P> {
    fn core(&self) -> &WorkerCoreState {
        &self.core
    }

    fn core_mut(&mut self) -> &mut WorkerCoreState {
        &mut self.core
    }
}

/// Shared state for FASTQ 7-step pipeline.
///
/// Generic parameter P is the processed type (output of `process_fn`).
pub struct FastqPipelineState<R: BufRead + Send, P: Send + MemoryEstimate> {
    /// Pipeline configuration.
    pub config: FastqPipelineConfig,

    // ========== Step 1: Read ==========
    /// The reader (BGZF or Decompressed).
    pub reader: Mutex<Option<FastqReader<R>>>,
    /// Whether the reader is a `RecordCount` reader (produces `FastqBoundaryBatch` directly).
    pub reader_is_record_count: bool,
    /// Flag indicating EOF has been reached.
    pub read_done: AtomicBool,
    /// Count of batches successfully read and pushed to q1.
    pub batches_read: AtomicU64,

    // ========== Queue 1: Read → Decompress ==========
    /// Batches waiting to be decompressed.
    pub q1_raw: ArrayQueue<(u64, FastqReadBatch)>,

    // ========== Queue 2: Decompress → Group (with reorder) ==========
    /// Decompressed batches waiting for grouping.
    pub q2_decompressed: ArrayQueue<(u64, FastqReadBatch)>,
    /// Reorder buffer to ensure grouping receives data in order.
    pub q2_reorder: Mutex<ReorderBuffer<FastqReadBatch>>,

    // ========== Q2 Reorder Buffer Atomic State (for lock-free admission control) ==========
    /// Atomic state for Q2 reorder buffer (`next_seq` and `heap_bytes`).
    /// Used by Decompress and Group steps for memory backpressure.
    pub q2_reorder_state: ReorderBufferState,

    // ========== Queue 2.5: FindBoundaries → Parse (parallel parse only) ==========
    /// State for finding FASTQ record boundaries (only used when `use_parallel_parse=true`).
    /// Uses per-stream internal locks for parallel access.
    pub boundary_state: FastqBoundaryState,
    /// Lock to ensure boundary batches are processed in order.
    /// Required because per-stream leftover handling needs chunks in order.
    pub boundary_lock: Mutex<()>,
    /// Flag indicating boundary finding is complete.
    pub boundaries_done: AtomicBool,
    /// Count of batches that have had boundaries found.
    pub batches_boundaries_found: AtomicU64,
    /// Batches with boundaries found, waiting to be parsed.
    pub q2_5_boundaries: ArrayQueue<(u64, FastqBoundaryBatch)>,
    /// Current heap bytes in Q2.5 boundaries queue.
    pub q2_5_boundaries_heap_bytes: AtomicU64,

    // ========== Queue 2.75: Parse → Group (parallel parse only) ==========
    /// Flag indicating parsing is complete.
    pub parse_done: AtomicBool,
    /// Count of batches that have been parsed.
    pub batches_parsed: AtomicU64,
    /// Parsed batches waiting to be grouped (overflow from reorder).
    pub q2_75_parsed: ArrayQueue<(u64, FastqParsedBatch)>,
    /// Reorder buffer for parsed batches (ensures grouper receives data in order).
    pub q2_75_reorder: Mutex<ReorderBuffer<FastqParsedBatch>>,
    /// Atomic state for Q2.75 reorder buffer (`next_seq` and `heap_bytes`).
    /// Used by Parse and Group steps for memory backpressure.
    pub q2_75_reorder_state: ReorderBufferState,

    // ========== Step 3: Group (exclusive) ==========
    /// The grouper.
    pub grouper: Mutex<FastqMultiStreamGrouper>,
    /// Count of batches that have been fully processed by the grouper.
    pub batches_grouped: AtomicU64,
    /// Flag indicating grouping is complete.
    pub group_done: AtomicBool,
    /// Next serial for templates (batch counter).
    pub next_template_serial: AtomicU64,
    /// Total number of individual templates pushed to Q3 (for debugging).
    pub total_templates_pushed: AtomicU64,
    /// Total number of individual records serialized (for completion check).
    /// This is incremented in the Serialize step and compared to `templates_written`.
    pub total_records_serialized: AtomicU64,

    // ========== Pending Templates ==========
    /// Pending templates being accumulated into a batch.
    pub pending_templates: Mutex<Vec<FastqTemplate>>,

    // ========== Output-Half State (Group → Process → Serialize → Compress → Write) ==========
    /// Shared output pipeline queues and state.
    pub output: OutputPipelineQueues<FastqTemplate, P>,

    // ========== Deadlock Detection ==========
    /// State for deadlock detection and recovery.
    pub deadlock_state: DeadlockState,
}

impl<R: BufRead + Send, P: Send + MemoryEstimate> FastqPipelineState<R, P> {
    /// Create a new pipeline state.
    #[must_use]
    pub fn new(
        config: FastqPipelineConfig,
        reader: FastqReader<R>,
        output: Box<dyn Write + Send>,
    ) -> Self {
        let cap = config.queue_capacity;
        let batch_size = config.batch_size;
        let num_streams = reader.num_streams();
        let is_record_count = reader.is_record_count();
        let stats = if config.collect_stats {
            config.shared_stats.clone().or_else(|| Some(Arc::new(PipelineStats::new())))
        } else {
            None
        };

        // Create deadlock detection config and state
        let deadlock_config =
            DeadlockConfig::new(config.deadlock_timeout_secs, config.deadlock_recover_enabled);
        let memory_limit = config.queue_memory_limit;
        let deadlock_state = DeadlockState::new(&deadlock_config, memory_limit);

        Self {
            config,
            reader: Mutex::new(Some(reader)),
            reader_is_record_count: is_record_count,
            read_done: AtomicBool::new(false),
            batches_read: AtomicU64::new(0),
            q1_raw: ArrayQueue::new(cap),
            q2_decompressed: ArrayQueue::new(cap),
            q2_reorder: Mutex::new(ReorderBuffer::new()),
            // Q2 reorder buffer atomic state (for lock-free admission control)
            q2_reorder_state: ReorderBufferState::new(memory_limit),
            // Parallel Parse step state (Q2.5: FindBoundaries → Parse)
            boundary_state: FastqBoundaryState::new(num_streams),
            boundary_lock: Mutex::new(()),
            boundaries_done: AtomicBool::new(false),
            batches_boundaries_found: AtomicU64::new(0),
            q2_5_boundaries: ArrayQueue::new(cap),
            q2_5_boundaries_heap_bytes: AtomicU64::new(0),
            // Parallel Parse step state (Q2.75: Parse → Group)
            parse_done: AtomicBool::new(false),
            batches_parsed: AtomicU64::new(0),
            q2_75_parsed: ArrayQueue::new(cap),
            q2_75_reorder: Mutex::new(ReorderBuffer::new()),
            q2_75_reorder_state: ReorderBufferState::new(memory_limit),
            // Note: For synchronized mode, the grouper is bypassed entirely -
            // templates are created directly in the Parse step.
            grouper: Mutex::new(FastqMultiStreamGrouper::new(num_streams)),
            batches_grouped: AtomicU64::new(0),
            group_done: AtomicBool::new(false),
            next_template_serial: AtomicU64::new(0),
            total_templates_pushed: AtomicU64::new(0),
            total_records_serialized: AtomicU64::new(0),
            pending_templates: Mutex::new(Vec::with_capacity(batch_size)),
            // Output-half state (Group → Process → Serialize → Compress → Write)
            output: OutputPipelineQueues::new(cap, output, stats, "Processed records"),
            deadlock_state,
        }
    }

    /// Record an error.
    pub fn set_error(&self, error: io::Error) {
        self.output.set_error(error);
    }

    /// Check if an error has occurred.
    pub fn has_error(&self) -> bool {
        self.output.has_error()
    }

    /// Check if Decompress step can proceed with pushing a batch to Q2.
    ///
    /// This implements memory-based backpressure on the Q2 reorder buffer to prevent
    /// unbounded memory growth when Group (exclusive step) falls behind.
    ///
    /// # Deadlock Prevention
    ///
    /// Always allows the serial that Group needs (`next_seq`) to proceed,
    /// even if over memory limit. This ensures Group can always make progress.
    #[must_use]
    pub fn can_decompress_proceed(&self, serial: u64) -> bool {
        // Delegate to Q2 reorder state's can_proceed method
        self.q2_reorder_state.can_proceed(serial)
    }

    /// Check if memory is at the backpressure threshold.
    ///
    /// Uses Q2 reorder buffer tracking (before Group step) to signal memory pressure
    /// to the scheduler. See [`super::base::BACKPRESSURE_THRESHOLD_BYTES`] for architecture details.
    #[must_use]
    pub fn is_memory_high(&self) -> bool {
        self.q2_reorder_state.is_memory_high()
    }

    /// Check if memory has drained below the low-water mark.
    ///
    /// Provides hysteresis to prevent thrashing: enter drain mode at backpressure
    /// threshold, only exit when drained to half that threshold.
    #[must_use]
    pub fn is_memory_drained(&self) -> bool {
        self.q2_reorder_state.is_memory_drained()
    }

    /// Check if Q4 processed queue memory is high (for backpressure in Process step).
    ///
    /// Uses the same threshold as BAM's Q5 processed queue to ensure consistent
    /// memory backpressure behavior across both pipelines.
    #[must_use]
    pub fn is_q4_memory_high(&self) -> bool {
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

    /// Take the stored error.
    pub fn take_error(&self) -> Option<io::Error> {
        self.output.take_error()
    }

    /// Check if the pipeline is complete.
    ///
    /// Uses queue-emptiness check to ensure all data has flowed through
    /// the entire pipeline and been written.
    pub fn is_complete(&self) -> bool {
        // Check flags - both read and group must be done
        let read_done = self.read_done.load(Ordering::Acquire);
        let group_done = self.group_done.load(Ordering::Acquire);
        if !read_done || !group_done {
            return false;
        }

        // When using parallel parse, also check intermediate flags
        if self.config.use_parallel_parse {
            let boundaries_done = self.boundaries_done.load(Ordering::Acquire);
            let parse_done = self.parse_done.load(Ordering::Acquire);
            if !boundaries_done || !parse_done {
                log::trace!(
                    "is_complete: parallel parse flags not done: boundaries_done={boundaries_done}, parse_done={parse_done}"
                );
                return false;
            }
        }

        // Check input-half ArrayQueues are empty (lock-free checks)
        if !self.q1_raw.is_empty() || !self.q2_decompressed.is_empty() {
            return false;
        }

        // When using parallel parse, also check intermediate queues
        if self.config.use_parallel_parse {
            let boundaries_queue_len = self.q2_5_boundaries.len();
            let parsed_queue_len = self.q2_75_parsed.len();
            if boundaries_queue_len > 0 || parsed_queue_len > 0 {
                log::trace!(
                    "is_complete: parallel parse queues not empty: q2_5={boundaries_queue_len}, q2_75={parsed_queue_len}"
                );
                return false;
            }
        }

        // Check input-half reorder buffers and pending templates (requires locks)
        let q2_empty = self.q2_reorder.lock().is_empty();
        let pending_empty = self.pending_templates.lock().is_empty();

        // When using parallel parse, also check the Q2.75 reorder buffer
        if self.config.use_parallel_parse {
            let q2_75_reorder_empty = self.q2_75_reorder.lock().is_empty();
            if !q2_empty || !pending_empty || !q2_75_reorder_empty {
                log::trace!(
                    "is_complete: reorder buffers: q2={}, pending={}, q2_75_reorder={}",
                    !q2_empty,
                    !pending_empty,
                    !q2_75_reorder_empty
                );
                return false;
            }
        } else if !q2_empty || !pending_empty {
            return false;
        }

        // Delegate output-half check
        self.output.are_queues_empty()
    }

    /// Get optional reference to pipeline statistics.
    #[must_use]
    pub fn stats(&self) -> Option<&PipelineStats> {
        self.output.stats.as_deref()
    }

    /// Get reference to progress tracker.
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
    ///
    /// # Errors
    ///
    /// Returns an I/O error if writing the BGZF EOF or flushing fails.
    pub fn flush_output(&self) -> io::Result<()> {
        if let Some(mut writer) = self.output.output.lock().take() {
            // Write BGZF EOF marker before flushing
            writer.write_all(&BGZF_EOF)?;
            writer.flush()?;
        }
        Ok(())
    }

    /// Validate pipeline completion to detect data loss.
    ///
    /// Checks that:
    /// 1. All queues are empty
    /// 2. All batch counters match between stages
    /// 3. Internal buffers are empty (grouper state, boundary leftovers)
    ///
    /// Note: Heap byte tracking is reported but advisory only (set to 0) because
    /// estimation can be imprecise. Only queue/buffer emptiness and counter checks
    /// cause validation failure.
    ///
    /// # Errors
    ///
    /// Returns `PipelineValidationError` with diagnostics if any issues are detected.
    #[allow(clippy::too_many_lines)]
    pub fn validate_completion(&self) -> Result<(), PipelineValidationError> {
        let mut non_empty_queues = Vec::new();
        let mut counter_mismatches = Vec::new();

        // Check input-half queues are empty
        if !self.q1_raw.is_empty() {
            non_empty_queues.push(format!("q1_raw ({})", self.q1_raw.len()));
        }
        if !self.q2_decompressed.is_empty() {
            non_empty_queues.push(format!("q2_decompressed ({})", self.q2_decompressed.len()));
        }

        // Check parallel parse queues if enabled
        if self.config.use_parallel_parse {
            if !self.q2_5_boundaries.is_empty() {
                non_empty_queues.push(format!("q2_5_boundaries ({})", self.q2_5_boundaries.len()));
            }
            if !self.q2_75_parsed.is_empty() {
                non_empty_queues.push(format!("q2_75_parsed ({})", self.q2_75_parsed.len()));
            }
        }

        // Check output-half queues are empty
        if !self.output.groups.is_empty() {
            non_empty_queues.push(format!("q3_templates ({})", self.output.groups.len()));
        }
        if !self.output.processed.is_empty() {
            non_empty_queues.push(format!("q4_processed ({})", self.output.processed.len()));
        }
        if !self.output.serialized.is_empty() {
            non_empty_queues.push(format!("q5_serialized ({})", self.output.serialized.len()));
        }
        if !self.output.compressed.is_empty() {
            non_empty_queues.push(format!("q6_compressed ({})", self.output.compressed.len()));
        }

        // Check reorder buffers are empty
        {
            let q2_reorder = self.q2_reorder.lock();
            if !q2_reorder.is_empty() {
                non_empty_queues.push(format!("q2_reorder ({})", q2_reorder.len()));
            }
        }
        if self.config.use_parallel_parse {
            let q2_75_reorder = self.q2_75_reorder.lock();
            if !q2_75_reorder.is_empty() {
                non_empty_queues.push(format!("q2_75_reorder ({})", q2_75_reorder.len()));
            }
        }
        {
            let write_reorder = self.output.write_reorder.lock();
            if !write_reorder.is_empty() {
                non_empty_queues.push(format!("write_reorder ({})", write_reorder.len()));
            }
        }
        {
            let pending = self.pending_templates.lock();
            if !pending.is_empty() {
                non_empty_queues.push(format!("pending_templates ({})", pending.len()));
            }
        }

        // Check grouper has no pending data
        {
            let grouper = self.grouper.lock();
            if grouper.has_pending() {
                non_empty_queues.push("grouper_pending".to_string());
            }
        }

        // Check boundary state has no leftover bytes (only in parallel parse mode)
        if self.config.use_parallel_parse {
            for (idx, stream_state) in self.boundary_state.stream_states.iter().enumerate() {
                let leftover_len = stream_state.lock().leftover.len();
                if leftover_len > 0 {
                    non_empty_queues.push(format!("boundary_leftover[{idx}] ({leftover_len})"));
                }
            }
        }

        // Check batch counter invariants
        let batches_read = self.batches_read.load(Ordering::Acquire);
        let batches_grouped = self.batches_grouped.load(Ordering::Acquire);

        // In non-parallel mode: batches flow Read -> Group directly
        // In parallel mode: batches flow Read -> FindBoundaries -> Parse -> Group
        if self.config.use_parallel_parse {
            let batches_boundaries_found = self.batches_boundaries_found.load(Ordering::Acquire);
            let batches_parsed = self.batches_parsed.load(Ordering::Acquire);

            if batches_boundaries_found != batches_read {
                counter_mismatches.push(format!(
                    "batches_boundaries_found ({batches_boundaries_found}) != batches_read ({batches_read})"
                ));
            }
            if batches_parsed != batches_boundaries_found {
                counter_mismatches.push(format!(
                    "batches_parsed ({batches_parsed}) != batches_boundaries_found ({batches_boundaries_found})"
                ));
            }
            if batches_grouped != batches_parsed {
                counter_mismatches.push(format!(
                    "batches_grouped ({batches_grouped}) != batches_parsed ({batches_parsed})"
                ));
            }
        } else {
            // Non-parallel mode: batches_grouped should match batches_read
            if batches_grouped != batches_read {
                counter_mismatches.push(format!(
                    "batches_grouped ({batches_grouped}) != batches_read ({batches_read})"
                ));
            }
        }

        // Note: Heap byte tracking can have small imbalances due to estimation errors,
        // so we don't fail validation on heap bytes. The important checks are queues
        // (actual data) and counters (batch flow).
        let leaked_heap_bytes = 0u64;

        // Return error if any issues found
        if !non_empty_queues.is_empty() || !counter_mismatches.is_empty() {
            return Err(PipelineValidationError {
                non_empty_queues,
                counter_mismatches,
                leaked_heap_bytes,
            });
        }

        Ok(())
    }
}

// ============================================================================
// PipelineLifecycle Implementation
// ============================================================================

impl<R: BufRead + Send + 'static, P: Send + MemoryEstimate + 'static> PipelineLifecycle
    for FastqPipelineState<R, P>
{
    fn is_complete(&self) -> bool {
        FastqPipelineState::is_complete(self)
    }

    fn has_error(&self) -> bool {
        FastqPipelineState::has_error(self)
    }

    fn take_error(&self) -> Option<io::Error> {
        FastqPipelineState::take_error(self)
    }

    fn set_error(&self, error: io::Error) {
        self.output.set_error(error);
    }

    fn is_draining(&self) -> bool {
        self.output.is_draining()
    }

    fn set_draining(&self, value: bool) {
        FastqPipelineState::set_draining(self, value);
    }

    fn stats(&self) -> Option<&PipelineStats> {
        FastqPipelineState::stats(self)
    }

    fn progress(&self) -> &ProgressTracker {
        FastqPipelineState::progress(self)
    }

    fn items_written(&self) -> u64 {
        FastqPipelineState::items_written(self)
    }

    fn flush_output(&self) -> io::Result<()> {
        FastqPipelineState::flush_output(self)
    }

    fn validate_completion(&self) -> Result<(), PipelineValidationError> {
        FastqPipelineState::validate_completion(self)
    }
}

// ============================================================================
// MonitorableState Trait Implementation
// ============================================================================

impl<R: BufRead + Send + 'static, P: Send + MemoryEstimate + 'static> MonitorableState
    for FastqPipelineState<R, P>
{
    fn deadlock_state(&self) -> &DeadlockState {
        &self.deadlock_state
    }

    fn build_queue_snapshot(&self) -> QueueSnapshot {
        // Track parallel parse queues when enabled
        let (q2b_len, q3_reorder_mem) = if self.config.use_parallel_parse {
            (
                self.q2_5_boundaries.len() + self.q2_75_parsed.len(),
                self.q2_75_reorder_state.get_heap_bytes(),
            )
        } else {
            (0, 0)
        };

        QueueSnapshot {
            q1_len: self.q1_raw.len(),
            q2_len: self.q2_decompressed.len(),
            q2b_len,
            q3_len: self.output.groups.len(),
            q4_len: self.output.processed.len(),
            q5_len: self.output.serialized.len(),
            q6_len: self.output.compressed.len(),
            q7_len: 0, // Not used in FASTQ (q6_compressed is the write input)
            q2_reorder_mem: self.q2_reorder_state.get_heap_bytes(),
            q3_reorder_mem,
            memory_limit: self.deadlock_state.get_memory_limit(),
            read_done: self.read_done.load(Ordering::Relaxed),
            group_done: self.group_done.load(Ordering::Relaxed),
            draining: self.output.draining.load(Ordering::Relaxed),
        }
    }
}

impl<R: BufRead + Send + 'static, P: Send + MemoryEstimate + 'static> OutputPipelineState
    for FastqPipelineState<R, P>
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
        let heap_size = item.1.estimate_heap_size();
        let result = self.output.compressed.push(item);
        if result.is_ok() {
            self.output.compressed_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
        }
        result
    }

    fn q6_is_full(&self) -> bool {
        self.output.compressed.is_full()
    }

    fn q6_track_pop(&self, heap_size: u64) {
        self.output.compressed_heap_bytes.fetch_sub(heap_size, Ordering::AcqRel);
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
        self.output.stats.as_deref()
    }
}

// ============================================================================
// New Shared Traits (Phase 3 - Pipeline Consolidation)
// ============================================================================

impl<R: BufRead + Send + 'static, P: Send + MemoryEstimate + 'static>
    ProcessPipelineState<FastqTemplate, P> for FastqPipelineState<R, P>
{
    fn process_input_pop(&self) -> Option<(u64, Vec<FastqTemplate>)> {
        self.output.groups.pop()
    }

    fn process_output_is_full(&self) -> bool {
        self.output.processed.is_full()
    }

    fn process_output_push(&self, item: (u64, Vec<P>)) -> Result<(), (u64, Vec<P>)> {
        self.output.processed.push(item)
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

impl<R: BufRead + Send + 'static, P: Send + MemoryEstimate + 'static> SerializePipelineState<P>
    for FastqPipelineState<R, P>
{
    fn serialize_input_pop(&self) -> Option<(u64, Vec<P>)> {
        self.output.processed.pop()
    }

    fn serialize_output_is_full(&self) -> bool {
        self.output.serialized.is_full()
    }

    fn serialize_output_push(
        &self,
        item: (u64, SerializedBatch),
    ) -> Result<(), (u64, SerializedBatch)> {
        self.output.serialized.push(item)
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

impl<R: BufRead + Send + 'static, P: Send + MemoryEstimate + 'static> WritePipelineState
    for FastqPipelineState<R, P>
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
        self.output.stats.as_deref()
    }
}

// ========== Step Functions ==========

/// Try to read the next batch (Step 1).
///
/// For `RecordCount` readers (gzip+synchronized), reads complete FASTQ records and
/// pushes `FastqBoundaryBatch` directly to `q2_5_boundaries`, skipping Decompress
/// and `FindBoundaries` steps entirely.
///
/// For `Bgzf` and `Decompressed` readers, reads raw/decompressed chunks to `q1_raw`.
fn fastq_try_step_read<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
    worker: &mut FastqWorkerState<P>,
) -> bool {
    // Check reader type once: RecordCount path skips q1_raw entirely
    if state.reader_is_record_count {
        return fastq_try_step_read_record_count(state, worker);
    }

    // =========================================================================
    // Standard path (Bgzf / Decompressed): read raw batches to q1_raw
    // =========================================================================

    // Priority 1: Try to advance any held item first
    if let Some((serial, held)) = worker.held_raw.take() {
        match state.q1_raw.push((serial, held)) {
            Ok(()) => {
                state.deadlock_state.record_q1_push();
            }
            Err((serial, held)) => {
                worker.held_raw = Some((serial, held));
                return false;
            }
        }
    }

    // Priority 2: Check for completion/error
    if state.read_done.load(Ordering::Relaxed) {
        return false;
    }
    if state.has_error() {
        return false;
    }

    // Priority 3: Check if output queue has space
    if state.q1_raw.is_full() {
        return false;
    }

    // Priority 4: Acquire exclusive reader lock
    let Some(mut guard) = state.reader.try_lock() else {
        if let Some(stats) = state.stats() {
            stats.record_contention(PipelineStep::Read);
        }
        return false;
    };
    let Some(ref mut reader) = *guard else {
        return false;
    };

    // Priority 5: Read next batch
    match reader.read_next_batch() {
        Ok(Some(batch)) => {
            let serial = batch.serial;
            match state.q1_raw.push((serial, batch)) {
                Ok(()) => {
                    state.batches_read.fetch_add(1, Ordering::Release);
                    state.deadlock_state.record_q1_push();
                    true
                }
                Err((serial, batch)) => {
                    worker.held_raw = Some((serial, batch));
                    state.batches_read.fetch_add(1, Ordering::Release);
                    false
                }
            }
        }
        Ok(None) => {
            state.read_done.store(true, Ordering::SeqCst);
            false
        }
        Err(e) => {
            state.set_error(e);
            false
        }
    }
}

/// Read step for `RecordCount` readers: produces `FastqBoundaryBatch` directly
/// to `q2_5_boundaries`, bypassing Decompress and `FindBoundaries`.
fn fastq_try_step_read_record_count<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
    worker: &mut FastqWorkerState<P>,
) -> bool {
    // Priority 1: Try to advance held boundary batch
    if let Some((serial, held)) = worker.held_boundaries.take() {
        let heap_size = held.estimate_heap_size();
        match state.q2_5_boundaries.push((serial, held)) {
            Ok(()) => {
                state.q2_5_boundaries_heap_bytes.fetch_add(heap_size as u64, Ordering::Relaxed);
                state.deadlock_state.record_q2_5_push();
            }
            Err((serial, held)) => {
                worker.held_boundaries = Some((serial, held));
                return false;
            }
        }
    }

    // Priority 2: Check for completion/error
    if state.read_done.load(Ordering::Relaxed) {
        return false;
    }
    if state.has_error() {
        return false;
    }

    // Priority 3: Check if output queue has space
    if state.q2_5_boundaries.is_full() {
        return false;
    }

    // Priority 4: Acquire exclusive reader lock
    let Some(mut guard) = state.reader.try_lock() else {
        if let Some(stats) = state.stats() {
            stats.record_contention(PipelineStep::Read);
        }
        return false;
    };
    let Some(ref mut reader) = *guard else {
        return false;
    };

    // Priority 5: Read next boundary batch
    match reader.read_next_boundary_batch() {
        Ok(Some(batch)) => {
            let serial = batch.serial;
            let heap_size = batch.estimate_heap_size();
            match state.q2_5_boundaries.push((serial, batch)) {
                Ok(()) => {
                    state.q2_5_boundaries_heap_bytes.fetch_add(heap_size as u64, Ordering::Relaxed);
                    state.deadlock_state.record_q2_5_push();
                    state.batches_read.fetch_add(1, Ordering::Release);
                    // RecordCount reader produces boundaries directly (Read + FindBoundaries combined)
                    state.batches_boundaries_found.fetch_add(1, Ordering::Release);
                    true
                }
                Err((serial, batch)) => {
                    worker.held_boundaries = Some((serial, batch));
                    state.batches_read.fetch_add(1, Ordering::Release);
                    state.batches_boundaries_found.fetch_add(1, Ordering::Release);
                    false
                }
            }
        }
        Ok(None) => {
            // EOF: set read_done and boundaries_done (skip Decompress + FindBoundaries)
            state.read_done.store(true, Ordering::SeqCst);
            state.boundaries_done.store(true, Ordering::SeqCst);
            false
        }
        Err(e) => {
            state.set_error(e);
            false
        }
    }
}

/// Try to decompress a batch (Step 2).
///
/// # Memory-Based Backpressure
///
/// Before pushing to `q2_decompressed`, checks if the Q2 reorder buffer would accept
/// this serial number. This prevents unbounded memory growth in the reorder buffer
/// while avoiding deadlock by always accepting the serial that Group needs.
fn fastq_try_step_decompress<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
    worker: &mut FastqWorkerState<P>,
) -> bool {
    // =========================================================================
    // Priority 1: Try to advance any held item first
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
                // Successfully pushed - keep reservation (released by Group when popping)
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
    // Priority 2: Check for errors
    // =========================================================================
    if state.has_error() {
        return false;
    }

    // =========================================================================
    // Priority 3: Check if output queue has space (soft check)
    // =========================================================================
    if state.q2_decompressed.is_full() {
        return false;
    }

    // =========================================================================
    // Priority 4: Pop from input queue
    // =========================================================================
    let Some((serial, batch)) = state.q1_raw.pop() else {
        if let Some(stats) = state.stats() {
            stats.record_queue_empty(1);
        }
        return false;
    };
    state.deadlock_state.record_q1_pop();

    // =========================================================================
    // Priority 5: Decompress (or passthrough)
    // =========================================================================
    match decompress_fastq_batch(batch, &mut worker.decompressor) {
        Ok(decompressed) => {
            // =========================================================================
            // Priority 6: Calculate and reserve memory BEFORE checking backpressure
            // =========================================================================
            let heap_size = decompressed.estimate_heap_size();
            state.q2_reorder_state.add_heap_bytes(heap_size as u64);

            // Check memory-based backpressure before trying to push
            if !state.can_decompress_proceed(serial) {
                // Reorder buffer is over limit - release reservation and hold
                state.q2_reorder_state.sub_heap_bytes(heap_size as u64);
                worker.held_decompressed = Some((serial, decompressed, heap_size));
                return false;
            }

            // =========================================================================
            // Priority 7: Try to push result (non-blocking)
            // =========================================================================
            match state.q2_decompressed.push((serial, decompressed)) {
                Ok(()) => {
                    state.deadlock_state.record_q2_push();
                    true
                }
                Err((serial, decompressed)) => {
                    // Output full - release reservation and hold
                    state.q2_reorder_state.sub_heap_bytes(heap_size as u64);
                    worker.held_decompressed = Some((serial, decompressed, heap_size));
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

// ============================================================================
// Parallel Parse Step Functions (for use_parallel_parse=true)
// ============================================================================

/// Try to find FASTQ record boundaries (Step 2.5).
///
/// SYNC WITH: bam.rs `try_step_find_boundaries()`
/// Both implementations use the "held boundaries" pattern for parallelism.
/// See base.rs `HasHeldBoundaries` trait for pattern documentation.
///
/// This step scans decompressed data for record boundaries (@ characters at
/// line starts). Uses per-stream internal locks for parallel access.
///
/// Only used when `config.use_parallel_parse` is true.
///
/// Returns (`did_work`, `had_contention`).
fn fastq_try_step_find_boundaries<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
    worker: &mut FastqWorkerState<P>,
) -> (bool, bool) {
    if !state.config.use_parallel_parse {
        return (false, false); // Not using parallel parse
    }
    if state.boundaries_done.load(Ordering::Relaxed) || state.has_error() {
        return (false, false);
    }

    // =========================================================================
    // Priority 1: Try to advance any held boundary batch first
    // =========================================================================
    let mut did_work = false;
    if let Some((serial, held)) = worker.held_boundaries.take() {
        let boundary_heap_size = held.estimate_heap_size();
        match state.q2_5_boundaries.push((serial, held)) {
            Ok(()) => {
                state
                    .q2_5_boundaries_heap_bytes
                    .fetch_add(boundary_heap_size as u64, Ordering::Relaxed);
                state.batches_boundaries_found.fetch_add(1, Ordering::Release);
                state.deadlock_state.record_q2_5_push();
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
    if state.q2_5_boundaries.is_full() {
        return (false, false); // Backpressure, not contention
    }

    // =========================================================================
    // Priority 3: Acquire boundary lock to ensure ordered processing
    // =========================================================================
    // FASTQ per-stream leftover handling requires batches to be processed in order.
    // Use try_lock() so other threads can do different work instead of blocking.
    let Some(_boundary_guard) = state.boundary_lock.try_lock() else {
        return (did_work, true); // Contention on boundary lock
    };

    // =========================================================================
    // Priority 4: Get next batch from reorder buffer (under boundary lock)
    // =========================================================================
    let batch_opt = {
        let mut reorder = state.q2_reorder.lock();

        // Insert pending items from Q2 ArrayQueue to reorder buffer
        while let Some((serial, batch)) = state.q2_decompressed.pop() {
            let heap_size = batch.estimate_heap_size();
            reorder.insert_with_size(serial, batch, heap_size);
            state.deadlock_state.record_q2_pop();
        }

        // Try to pop next in-order batch
        let batch_opt = reorder.try_pop_next();

        // Update Q2 reorder tracking
        state.q2_reorder_state.update_next_seq(reorder.next_seq());

        // Check for completion while we have the lock
        let read_done = state.read_done.load(Ordering::Acquire);
        let batches_read = state.batches_read.load(Ordering::Acquire);
        let batches_boundaries_found = state.batches_boundaries_found.load(Ordering::Acquire);
        let all_processed = batches_read == batches_boundaries_found;
        let reorder_empty = reorder.is_empty();

        if batch_opt.is_none() && read_done && all_processed && reorder_empty {
            state.boundaries_done.store(true, Ordering::Release);
            log::info!(
                "fastq_try_step_find_boundaries: set boundaries_done=true (read={batches_read}, found={batches_boundaries_found})"
            );
        }

        batch_opt
    }; // Reorder lock released here (boundary lock still held)

    let Some(batch) = batch_opt else {
        if let Some(stats) = state.stats() {
            stats.record_queue_empty(2);
        }
        return (did_work, false); // No more data available
    };

    let serial = batch.serial;
    let heap_size = batch.estimate_heap_size();

    // Release memory from atomic tracker when popping from reorder buffer
    state.q2_reorder_state.sub_heap_bytes(heap_size as u64);

    // =========================================================================
    // Priority 5: Find boundaries (under boundary lock for ordering)
    // =========================================================================
    // Convert FastqReadBatch to FastqDecompressedBatch for the format trait
    let decompressed = FastqDecompressedBatch {
        chunks: batch
            .chunks
            .into_iter()
            .map(|c| FastqDecompressedChunk { stream_idx: c.stream_idx, data: c.data })
            .collect(),
        serial,
    };

    // Find boundaries - must be done in order due to per-stream leftover state
    let boundary_batch = match FastqFormat::find_boundaries(&state.boundary_state, decompressed) {
        Ok(batch) => batch,
        Err(e) => {
            state.set_error(e);
            return (true, false); // Did work (processed batch from reorder buffer)
        }
    };

    // =========================================================================
    // Priority 6: Push or hold (boundary lock released after this)
    // =========================================================================
    let boundary_heap_size = boundary_batch.estimate_heap_size();
    match state.q2_5_boundaries.push((serial, boundary_batch)) {
        Ok(()) => {
            state
                .q2_5_boundaries_heap_bytes
                .fetch_add(boundary_heap_size as u64, Ordering::Relaxed);
            state.batches_boundaries_found.fetch_add(1, Ordering::Release);
            state.deadlock_state.record_q2_5_push();
            (true, false) // Did work, no contention
        }
        Err((serial, batch)) => {
            // Output full - hold the result and return
            worker.held_boundaries = Some((serial, batch));
            (true, false) // Did work (processed data), will retry push later
        }
    }
}

/// Try to parse FASTQ records (Step 2.75 - parallel).
///
/// This step takes boundary batches and constructs `FastqRecord` objects.
/// **This is the KEY PARALLEL STEP** that fixes the t8 scaling bottleneck.
///
/// Only used when `config.use_parallel_parse` is true.
///
/// For synchronized FASTQs, this step also creates templates and pushes them
/// directly to Q3, bypassing the Group step entirely. This eliminates all lock
/// contention since each worker can create templates independently.
#[allow(clippy::too_many_lines)]
fn fastq_try_step_parse<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
) -> bool {
    if !state.config.use_parallel_parse {
        return false; // Not using parallel parse
    }
    if state.parse_done.load(Ordering::Relaxed) || state.has_error() {
        return false;
    }

    // CRITICAL: Check output queue capacity BEFORE popping input to prevent data loss
    // For synchronized mode, check Q3 (templates) instead of Q2.75 (parsed)
    if state.config.synchronized {
        if state.output.groups.is_full() {
            return false; // Output queue full
        }
    } else if state.q2_75_parsed.is_full() {
        return false; // Output queue full, wait for Group step to drain it
    }

    // Pop from Q2.5 boundaries queue
    let Some((serial, boundary_batch)) = state.q2_5_boundaries.pop() else {
        if let Some(stats) = state.stats() {
            stats.record_queue_empty(25); // Q2b/Q2.5 (boundaries queue)
        }
        // Check if parsing is complete
        let boundaries_done = state.boundaries_done.load(Ordering::Acquire);
        let all_parsed = state.batches_boundaries_found.load(Ordering::Acquire)
            == state.batches_parsed.load(Ordering::Acquire);

        if boundaries_done && all_parsed && state.q2_5_boundaries.is_empty() {
            state.parse_done.store(true, Ordering::Release);
            // In synchronized mode, Group is skipped — set group_done alongside parse_done
            if state.config.synchronized {
                state.group_done.store(true, Ordering::Release);
                log::trace!(
                    "fastq_try_step_parse: set parse_done=true, group_done=true (synchronized)"
                );
            } else {
                log::trace!("fastq_try_step_parse: set parse_done=true");
            }
        } else if let Some(stats) = state.stats() {
            // Record Q2.5 as extension of Q2
            stats.record_queue_empty(2);
        }
        return false;
    };
    state.deadlock_state.record_q2_5_pop();

    let input_heap_size = boundary_batch.estimate_heap_size();

    // Parse records - THIS IS THE KEY PARALLEL OPERATION
    match FastqFormat::parse_records(boundary_batch) {
        Ok(parsed_batch) => {
            // Only decrement input memory AFTER successful parse
            state.q2_5_boundaries_heap_bytes.fetch_sub(input_heap_size as u64, Ordering::Relaxed);

            // SYNCHRONIZED FAST PATH: Create templates directly, bypass Group step
            if state.config.synchronized {
                // Create templates by zipping records at matching positions
                let templates = match create_templates_from_streams(parsed_batch.streams) {
                    Ok(t) => t,
                    Err(e) => {
                        state.set_error(e);
                        return false;
                    }
                };

                if templates.is_empty() {
                    // Empty batch - just count it as parsed and grouped
                    state.batches_parsed.fetch_add(1, Ordering::Release);
                    state.batches_grouped.fetch_add(1, Ordering::Release);
                    return true;
                }

                let count = templates.len();
                // Use the original batch serial to preserve input order.
                // The downstream reorder buffer will ensure batches are written in order
                // even if they complete out-of-order due to parallel processing.
                let template_serial = serial;

                // Push templates directly to Q3 with spin-wait
                let mut templates_to_push = templates;
                #[cfg(feature = "memory-debug")]
                let q4_heap: u64 =
                    templates_to_push.iter().map(|t| t.estimate_heap_size() as u64).sum();
                loop {
                    match state.output.groups.push((template_serial, templates_to_push)) {
                        Ok(()) => {
                            #[cfg(feature = "memory-debug")]
                            state.output.groups_heap_bytes.fetch_add(q4_heap, Ordering::AcqRel);
                            state.total_templates_pushed.fetch_add(count as u64, Ordering::Release);
                            if let Some(stats) = state.stats() {
                                stats.groups_produced.fetch_add(count as u64, Ordering::Relaxed);
                            }
                            state.deadlock_state.record_q4_push();
                            state.batches_parsed.fetch_add(1, Ordering::Release);
                            state.batches_grouped.fetch_add(1, Ordering::Release);
                            return true;
                        }
                        Err(returned) => {
                            templates_to_push = returned.1;
                            if state.has_error() {
                                return false;
                            }
                            std::hint::spin_loop();
                        }
                    }
                }
            }

            // NON-SYNCHRONIZED PATH: Push to Q2.75 for Group step
            let parsed_heap_size = parsed_batch.estimate_heap_size();

            // Try to insert into Q2.75 reorder buffer first
            if let Some(mut reorder) = state.q2_75_reorder.try_lock() {
                reorder.insert_with_size(serial, parsed_batch, parsed_heap_size);
                state.q2_75_reorder_state.add_heap_bytes(parsed_heap_size as u64);
                state.batches_parsed.fetch_add(1, Ordering::Release);
                true
            } else {
                // Couldn't get reorder lock - push to ArrayQueue with spin-wait
                // CRITICAL: We MUST NOT lose data since we already consumed input
                let mut batch_to_push = parsed_batch;
                loop {
                    match state.q2_75_parsed.push((serial, batch_to_push)) {
                        Ok(()) => {
                            state.batches_parsed.fetch_add(1, Ordering::Release);
                            return true;
                        }
                        Err(returned) => {
                            batch_to_push = returned.1;
                            if state.has_error() {
                                return false;
                            }
                            // Try reorder buffer again
                            if let Some(mut reorder) = state.q2_75_reorder.try_lock() {
                                reorder.insert_with_size(serial, batch_to_push, parsed_heap_size);
                                state.q2_75_reorder_state.add_heap_bytes(parsed_heap_size as u64);
                                state.batches_parsed.fetch_add(1, Ordering::Release);
                                return true;
                            }
                            std::hint::spin_loop();
                        }
                    }
                }
            }
        }
        Err(e) => {
            // Batch already removed from Q2.5; keep heap tracking consistent
            state.q2_5_boundaries_heap_bytes.fetch_sub(input_heap_size as u64, Ordering::Relaxed);
            state.set_error(e);
            false
        }
    }
}

/// Create templates by zipping records from synchronized FASTQ streams.
///
/// For single-end: each record becomes its own template.
/// For paired-end: R1[i] and R2[i] are zipped into a single template.
fn create_templates_from_streams(
    mut streams: Vec<FastqParsedStream>,
) -> io::Result<Vec<FastqTemplate>> {
    let num_streams = streams.len();

    match num_streams {
        0 => Ok(Vec::new()),
        1 => {
            // Single-end: each record becomes its own template
            let records = streams.pop().unwrap().records;
            Ok(records
                .into_iter()
                .map(|r| {
                    let name = r.name.clone();
                    FastqTemplate { name, records: vec![r] }
                })
                .collect())
        }
        2 => {
            // Paired-end: zip R1 and R2 by position.
            // Sort by stream_idx so streams[0] is always R1 and streams[1] is
            // always R2, regardless of the order produced by find_boundaries().
            streams.sort_by_key(|s| s.stream_idx);
            let mut drain = streams.into_iter();
            let r1_records = drain.next().unwrap().records;
            let r2_records = drain.next().unwrap().records;

            // Validate batch sizes match
            if r1_records.len() != r2_records.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "FASTQ batch size mismatch: R1 has {} records, R2 has {} records",
                        r1_records.len(),
                        r2_records.len()
                    ),
                ));
            }

            // Zip directly by position - no name validation here (done in Process step)
            let templates: Vec<FastqTemplate> = r1_records
                .into_iter()
                .zip(r2_records)
                .map(|(r1, r2)| {
                    let name = r1.name.clone();
                    FastqTemplate { name, records: vec![r1, r2] }
                })
                .collect();

            Ok(templates)
        }
        n => Err(io::Error::new(
            io::ErrorKind::Unsupported,
            format!("Synchronized mode not supported for {n} streams (max 2)"),
        )),
    }
}

/// Try to group parsed records (Step 3 with parallel parse - exclusive).
///
/// This version receives pre-parsed records from the Parse step instead of
/// parsing under the lock. This is much faster since parsing is the bottleneck.
///
/// Only used when `config.use_parallel_parse` is true.
///
/// For synchronized FASTQs, this step is a no-op since templates are created
/// directly in the Parse step. It only handles completion detection.
#[allow(clippy::too_many_lines)]
fn fastq_try_step_group_parsed<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
) -> bool {
    if !state.config.use_parallel_parse {
        return false; // Use the regular group step
    }
    if state.group_done.load(Ordering::Relaxed) || state.has_error() {
        return false;
    }

    // SYNCHRONIZED FAST PATH: Templates are created in Parse step, no grouping needed
    if state.config.synchronized {
        // Check if grouping is complete (all parsed batches have been grouped)
        let parse_done = state.parse_done.load(Ordering::Acquire);
        let all_grouped = state.batches_parsed.load(Ordering::Acquire)
            == state.batches_grouped.load(Ordering::Acquire);

        if parse_done && all_grouped {
            state.group_done.store(true, Ordering::Release);
            log::trace!("fastq_try_step_group_parsed: synchronized mode - set group_done=true");
        }
        return false; // No work to do - Parse step handles everything
    }

    // NON-SYNCHRONIZED PATH: Group templates under lock
    // Try to acquire grouper
    let Some(mut grouper_guard) = state.grouper.try_lock() else {
        if let Some(stats) = state.stats() {
            stats.record_contention(PipelineStep::Group);
        }
        return false;
    };

    // Also acquire pending_templates
    let Some(mut pending) = state.pending_templates.try_lock() else {
        return false;
    };

    // Also acquire the Q2.75 reorder buffer
    let Some(mut reorder) = state.q2_75_reorder.try_lock() else {
        return false;
    };

    let batch_size = state.config.batch_size;
    let mut processed_any = false;

    // Process from ArrayQueue first (lock already held, safe to drain)
    while let Some((serial, parsed_batch)) = state.q2_75_parsed.pop() {
        let heap_size = parsed_batch.estimate_heap_size();
        reorder.insert_with_size(serial, parsed_batch, heap_size);
    }

    // Process in order from reorder buffer
    while let Some(parsed_batch) = reorder.try_pop_next() {
        let heap_size = parsed_batch.estimate_heap_size();
        state.q2_75_reorder_state.update_next_seq(reorder.next_seq());
        state.q2_75_reorder_state.sub_heap_bytes(heap_size as u64);

        // Feed pre-parsed records to grouper - NO PARSING UNDER LOCK!
        match grouper_guard.add_parsed_batch(parsed_batch.streams) {
            Ok(templates) => {
                pending.extend(templates);

                // Flush full batches
                while pending.len() >= batch_size {
                    // Check if queue has space before allocating serial
                    if state.output.groups.is_full() {
                        break; // Queue full, keep pending for next iteration
                    }

                    let batch: Vec<FastqTemplate> = pending.drain(..batch_size).collect();
                    let count = batch.len();
                    #[cfg(feature = "memory-debug")]
                    let q4_heap: u64 = batch.iter().map(|t| t.estimate_heap_size() as u64).sum();
                    let serial = state.next_template_serial.fetch_add(1, Ordering::Relaxed);

                    match state.output.groups.push((serial, batch)) {
                        Ok(()) => {
                            #[cfg(feature = "memory-debug")]
                            state.output.groups_heap_bytes.fetch_add(q4_heap, Ordering::AcqRel);
                            state.total_templates_pushed.fetch_add(count as u64, Ordering::Release);
                            if let Some(stats) = state.stats() {
                                stats.groups_produced.fetch_add(count as u64, Ordering::Relaxed);
                            }
                            state.deadlock_state.record_q4_push();
                        }
                        Err(returned) => {
                            // Push failed (race condition) - restore templates and revert serial
                            // Prepend to maintain order
                            let mut restored = returned.1;
                            restored.extend(std::mem::take(&mut *pending));
                            *pending = restored;
                            state.next_template_serial.fetch_sub(1, Ordering::Release);
                            break;
                        }
                    }
                }

                processed_any = true;
                state.batches_grouped.fetch_add(1, Ordering::Release);
            }
            Err(e) => {
                state.set_error(e);
                return false;
            }
        }
    }

    // Check if grouping is complete
    let parse_done = state.parse_done.load(Ordering::Acquire);
    let all_grouped = state.batches_parsed.load(Ordering::Acquire)
        == state.batches_grouped.load(Ordering::Acquire);

    if parse_done && all_grouped && reorder.is_empty() && !grouper_guard.has_pending() {
        // Finish grouper
        match grouper_guard.finish() {
            Ok(Some(template)) => {
                pending.push(template);
            }
            Ok(None) => {}
            Err(e) => {
                state.set_error(e);
                return false;
            }
        }

        // Flush any remaining templates - CRITICAL: spin-wait to prevent data loss
        if !pending.is_empty() {
            let serial = state.next_template_serial.fetch_add(1, Ordering::Relaxed);
            let batch: Vec<FastqTemplate> = std::mem::take(&mut *pending);
            let count = batch.len();
            #[cfg(feature = "memory-debug")]
            let q4_heap: u64 = batch.iter().map(|t| t.estimate_heap_size() as u64).sum();

            let mut batch_to_push = batch;
            loop {
                match state.output.groups.push((serial, batch_to_push)) {
                    Ok(()) => {
                        #[cfg(feature = "memory-debug")]
                        state.output.groups_heap_bytes.fetch_add(q4_heap, Ordering::AcqRel);
                        state.total_templates_pushed.fetch_add(count as u64, Ordering::Release);
                        break;
                    }
                    Err(returned) => {
                        batch_to_push = returned.1;
                        if state.has_error() {
                            // On error, restore to pending so we don't lose data
                            pending.extend(batch_to_push);
                            return false;
                        }
                        std::hint::spin_loop();
                    }
                }
            }
        }

        state.group_done.store(true, Ordering::Release);
        log::trace!("fastq_try_step_group_parsed: set group_done=true");
    }

    processed_any
}

/// Try to group decompressed batches (Step 3 - exclusive).
///
/// This step batches templates together for efficient downstream processing.
/// Templates are accumulated in `pending_templates` until `batch_size` is reached,
/// then pushed to Q3 as a batch.
#[allow(clippy::too_many_lines)]
fn fastq_try_step_group<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
) -> bool {
    const MAX_PENDING_DRAIN: usize = 16;
    const MAX_RETRIES: u32 = 10_000;

    if state.group_done.load(Ordering::Relaxed) || state.has_error() {
        return false;
    }

    // Try to acquire grouper
    let Some(mut grouper_guard) = state.grouper.try_lock() else {
        // Record contention for diagnostics
        if let Some(stats) = state.stats() {
            stats.record_contention(PipelineStep::Group);
        }
        return false;
    };

    // Also acquire pending_templates lock
    let mut pending = state.pending_templates.lock();
    let batch_size = state.config.batch_size;

    // Pre-drain q2_decompressed BEFORE taking the reorder lock (lock-free operations)
    // This reduces critical section time by moving ArrayQueue ops outside the lock
    let mut pre_drained: Vec<(u64, FastqReadBatch, usize)> = Vec::with_capacity(MAX_PENDING_DRAIN);
    while pre_drained.len() < MAX_PENDING_DRAIN {
        if let Some((serial, batch)) = state.q2_decompressed.pop() {
            let heap_size = batch.estimate_heap_size();
            pre_drained.push((serial, batch, heap_size));
            state.deadlock_state.record_q2_pop();
        } else {
            break;
        }
    }

    // Drain reorder buffer for in-order batches
    let mut reorder = state.q2_reorder.lock();

    // Insert pre-drained items into reorder buffer
    for (serial, batch, heap_size) in pre_drained {
        log::trace!("fastq_try_step_group: inserting batch serial={serial} into reorder buffer");
        reorder.insert_with_size(serial, batch, heap_size);
    }

    // Helper to flush pending templates as a batch. Returns true if flushed, false if queue full.
    let flush_pending = |pending: &mut Vec<FastqTemplate>,
                         state: &FastqPipelineState<R, P>|
     -> bool {
        if pending.is_empty() {
            return true;
        }
        if state.output.groups.is_full() {
            return false; // Queue full, keep pending for next iteration
        }
        let batch: Vec<FastqTemplate> = std::mem::take(pending);
        let count = batch.len() as u64;
        #[cfg(feature = "memory-debug")]
        let q4_heap: u64 = batch.iter().map(|t| t.estimate_heap_size() as u64).sum();
        let serial = state.next_template_serial.fetch_add(1, Ordering::Release);
        log::trace!("fastq_try_step_group: pushing batch of {count} templates, serial={serial}");
        // Handle race condition: queue could fill between is_full check and push
        match state.output.groups.push((serial, batch)) {
            Ok(()) => {
                #[cfg(feature = "memory-debug")]
                state.output.groups_heap_bytes.fetch_add(q4_heap, Ordering::AcqRel);
                state.total_templates_pushed.fetch_add(count, Ordering::Release);
                // Record groups produced for throughput metrics
                if let Some(stats) = state.stats() {
                    stats.groups_produced.fetch_add(count, Ordering::Relaxed);
                }
                // Record progress for deadlock detection (Q3 in FASTQ maps to Q4 output in BAM model)
                state.deadlock_state.record_q4_push();
                true
            }
            Err(returned) => {
                // Push failed - restore templates to pending for next iteration
                pending.extend(returned.1);
                // Revert the serial increment
                state.next_template_serial.fetch_sub(1, Ordering::Release);
                false
            }
        }
    };

    // First, try to flush any pending templates from previous iterations
    // Always try to flush if non-empty (not just when >= batch_size) to handle final partial batches
    if !pending.is_empty() && !flush_pending(&mut pending, state) {
        return false; // Q3 still full, let other threads drain it
    }

    // Then, drain in-order batches from reorder buffer
    let mut processed_any = false;
    while let Some((batch, heap_size)) = reorder.try_pop_next_with_size() {
        // Update next_seq atomic for backpressure tracking
        state.q2_reorder_state.update_next_seq(reorder.next_seq());
        // Release the memory tracked for this batch
        state.q2_reorder_state.sub_heap_bytes(heap_size as u64);

        log::trace!("fastq_try_step_group: processing batch with {} chunks", batch.chunks.len());
        // Feed batch to grouper
        log::trace!("fastq_try_step_group: calling add_batch");
        match grouper_guard.add_batch(batch) {
            Ok(templates) => {
                log::trace!("fastq_try_step_group: got {} templates", templates.len());
                // Add ALL templates to pending batch first (so none are lost if flush fails)
                pending.extend(templates);
                // Then try to flush full batches
                while pending.len() >= batch_size {
                    // If flush fails (queue full), stop and retry later
                    if !flush_pending(&mut pending, state) {
                        // Note: We've already processed this batch in the grouper,
                        // so we just need to wait for Q3 to have space.
                        // Templates remain in pending and will be flushed later.
                        state.batches_grouped.fetch_add(1, Ordering::Release);
                        return true; // Did some work, will retry pending flush later
                    }
                }
                // Mark this batch as fully grouped
                state.batches_grouped.fetch_add(1, Ordering::Release);
                processed_any = true;
            }
            Err(e) => {
                log::error!("fastq_try_step_group: add_batch failed: {e:?}");
                state.set_error(e);
                return false;
            }
        }
    }

    // Check if grouping is complete
    // We need to verify:
    // 1. Reading is done (no more batches will be added)
    // 2. All batches that were read have been grouped (batches_grouped == batches_read)
    // 3. The grouper has no pending data
    let read_done = state.read_done.load(Ordering::Acquire);
    let batches_read = state.batches_read.load(Ordering::Acquire);
    let batches_grouped = state.batches_grouped.load(Ordering::Acquire);
    let all_batches_grouped = batches_grouped == batches_read;
    log::trace!(
        "fastq_try_step_group: check complete - read_done={}, batches_read={}, batches_grouped={}, has_pending={}",
        read_done,
        batches_read,
        batches_grouped,
        grouper_guard.has_pending()
    );
    if read_done && all_batches_grouped {
        // Finish grouper
        log::trace!("fastq_try_step_group: finishing grouper");
        match grouper_guard.finish() {
            Ok(Some(template)) => {
                log::trace!("fastq_try_step_group: finish produced 1 template");
                pending.push(template);
            }
            Ok(None) => {
                log::trace!("fastq_try_step_group: finish produced 0 templates");
            }
            Err(e) => {
                state.set_error(e);
                return false;
            }
        }
        // Flush any remaining pending templates as final batch.
        // CRITICAL: Must retry until flush succeeds to avoid data loss.
        // Use blocking loop with yield to allow consumer threads to drain.
        let mut retries: u32 = 0;
        while !pending.is_empty() {
            if flush_pending(&mut pending, state) {
                continue; // Flushed one batch, try next
            }
            retries += 1;
            if retries > MAX_RETRIES {
                state.set_error(io::Error::other(format!(
                    "Failed to flush {} final templates after {} retries",
                    pending.len(),
                    MAX_RETRIES
                )));
                return false;
            }
            // Yield to allow consumer threads to drain Q3
            std::thread::yield_now();
        }
        // All pending templates flushed
        state.group_done.store(true, Ordering::SeqCst);
        log::trace!("fastq_try_step_group: set group_done=true");
    }

    processed_any
}

/// Try to process a batch of templates (Step 4).
/// Always drains multiple batches when available for better throughput.
fn fastq_try_step_process<R: BufRead + Send, P: Send + MemoryEstimate, PF>(
    state: &FastqPipelineState<R, P>,
    process_fn: &PF,
    worker: &mut FastqWorkerState<P>,
) -> bool
where
    PF: Fn(FastqTemplate) -> io::Result<P>,
{
    // =========================================================================
    // Priority 1: Try to advance any held item first
    // =========================================================================
    if let Some((serial, held, heap_size)) = worker.held_processed.take() {
        match state.output.processed.push((serial, held)) {
            Ok(()) => {
                // Successfully advanced held item
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
    // Priority 2: Check for errors
    // =========================================================================
    if state.has_error() {
        return false;
    }

    // =========================================================================
    // Priority 3: Check if output queue has space (count and memory)
    // When draining, bypass memory backpressure to prevent deadlock
    // =========================================================================
    if state.output.processed.is_full() || (!state.is_draining() && state.is_q4_memory_high()) {
        return false;
    }

    // =========================================================================
    // Priority 4: Pop and process batches (multi-batch drain for throughput)
    // Q4 memory backpressure above prevents unbounded growth.
    // =========================================================================
    let max_batches = 8;
    let mut did_work = false;

    for _ in 0..max_batches {
        // Check output space (count and memory) before each batch
        // When draining, bypass memory backpressure to prevent deadlock
        if state.output.processed.is_full() || (!state.is_draining() && state.is_q4_memory_high()) {
            break;
        }

        let Some((serial, batch)) = state.output.groups.pop() else {
            if let Some(stats) = state.stats() {
                stats.record_queue_empty(4);
            }
            break;
        };
        state.deadlock_state.record_q4_pop();

        #[cfg(feature = "memory-debug")]
        {
            let q4_heap: u64 = batch.iter().map(|t| t.estimate_heap_size() as u64).sum();
            state.output.groups_heap_bytes.fetch_sub(q4_heap, Ordering::AcqRel);
        }

        log::trace!(
            "fastq_try_step_process: processing batch of {} templates, serial={}",
            batch.len(),
            serial
        );

        // Process each template in the batch
        let mut results: Vec<P> = Vec::with_capacity(batch.len());
        for template in batch {
            match process_fn(template) {
                Ok(processed) => results.push(processed),
                Err(e) => {
                    log::error!("fastq_try_step_process: error: {e:?}");
                    state.set_error(e);
                    return false;
                }
            }
        }

        log::trace!("fastq_try_step_process: processed {} items successfully", results.len());

        // Calculate heap size for memory tracking
        let heap_size: usize = results.iter().map(MemoryEstimate::estimate_heap_size).sum();

        // Try to push result (non-blocking)
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

/// Try to serialize a batch of processed data (Step 5).
fn fastq_try_step_serialize<R: BufRead + Send, P: Send + MemoryEstimate, SF>(
    state: &FastqPipelineState<R, P>,
    serialize_fn: &SF,
    header: &Header,
    worker: &mut FastqWorkerState<P>,
) -> bool
where
    SF: Fn(P, &Header, &mut Vec<u8>) -> io::Result<u64>,
{
    // =========================================================================
    // Priority 1: Try to advance any held item first
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
    // Priority 2: Check for errors
    // =========================================================================
    if state.has_error() {
        return false;
    }

    // =========================================================================
    // Priority 3: Check if output queue has space (soft check)
    // =========================================================================
    if state.output.serialized.is_full() {
        return false;
    }

    // =========================================================================
    // Priority 4: Pop batch from input queue
    // =========================================================================
    let Some((serial, batch)) = state.output.processed.pop() else {
        if let Some(stats) = state.stats() {
            stats.record_queue_empty(5);
        }
        return false;
    };
    state.deadlock_state.record_q5_pop();

    // Track memory being removed from Q4
    let q4_heap_size: usize = batch.iter().map(MemoryEstimate::estimate_heap_size).sum();
    state.output.processed_heap_bytes.fetch_sub(q4_heap_size as u64, Ordering::AcqRel);

    // =========================================================================
    // Priority 5: Serialize all items
    // =========================================================================
    // Prepare worker's serialization buffer
    worker.core.serialization_buffer.clear();

    // Serialize all items directly into worker's buffer (no intermediate allocation)
    log::trace!(
        "fastq_try_step_serialize: serializing batch of {} items, serial={}",
        batch.len(),
        serial
    );
    let mut total_record_count: u64 = 0;
    for item in batch {
        match serialize_fn(item, header, &mut worker.core.serialization_buffer) {
            Ok(record_count) => {
                total_record_count += record_count;
            }
            Err(e) => {
                log::error!("fastq_try_step_serialize: error: {e:?}");
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

    log::trace!(
        "fastq_try_step_serialize: serialized successfully, total_data_len={}, record_count={}",
        combined_data.len(),
        total_record_count
    );

    // Track total records serialized for completion check
    // This must use Release ordering to synchronize with Acquire in is_complete()
    state.total_records_serialized.fetch_add(total_record_count, Ordering::Release);

    // Record serialized bytes for throughput metrics
    if let Some(stats) = state.stats() {
        stats.serialized_bytes.fetch_add(combined_data.len() as u64, Ordering::Relaxed);
    }

    // =========================================================================
    // Priority 6: Try to push result (non-blocking)
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

/// Try to compress serialized data (Step 6).
/// Delegates to the shared implementation which uses the held-item pattern.
fn fastq_try_step_compress<R: BufRead + Send + 'static, P: Send + MemoryEstimate + 'static>(
    state: &FastqPipelineState<R, P>,
    worker: &mut FastqWorkerState<P>,
) -> bool {
    shared_try_step_compress(state, worker).is_success()
}

/// Try to write compressed data (Step 7 - exclusive).
fn fastq_try_step_write<R: BufRead + Send + 'static, P: Send + MemoryEstimate + 'static>(
    state: &FastqPipelineState<R, P>,
) -> bool {
    if state.has_error() {
        return false;
    }

    // Try to acquire output
    let Some(mut output_guard) = state.output.output.try_lock() else {
        // Record contention for diagnostics
        if let Some(stats) = state.stats() {
            stats.record_contention(PipelineStep::Write);
        }
        return false;
    };
    let Some(ref mut output) = *output_guard else {
        return false;
    };

    // Drain Q6 into reorder buffer AND write all ready batches in single lock scope
    let mut wrote_any = false;
    let q7_truly_empty;
    {
        let mut reorder = state.output.write_reorder.lock();

        // Move items from Q6 to reorder buffer
        while let Some((serial, batch)) = state.output.compressed.pop() {
            let q7_heap = batch.estimate_heap_size() as u64;
            state.q6_track_pop(q7_heap);
            reorder.insert(serial, batch);
            state.deadlock_state.record_q7_pop();
        }

        // Write in-order batches
        while let Some(batch) = reorder.try_pop_next() {
            // Write all blocks in the batch
            let mut batch_bytes: u64 = 0;
            for block in &batch.blocks {
                batch_bytes += block.data.len() as u64;
                if let Err(e) = output.write_all(&block.data) {
                    state.set_error(e);
                    return false;
                }
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

    wrote_any
}

/// Execute a pipeline step for FASTQ, returning `(success, was_contention)`.
/// Steps not applicable to FASTQ (`FindBoundaries`, `Decode`) return `(false, false)`.
fn fastq_execute_step<R: BufRead + Send + 'static, P: Send + MemoryEstimate + 'static, PF, SF>(
    state: &FastqPipelineState<R, P>,
    header: &Header,
    process_fn: &PF,
    serialize_fn: &SF,
    worker: &mut FastqWorkerState<P>,
    step: PipelineStep,
) -> (bool, bool)
where
    PF: Fn(FastqTemplate) -> io::Result<P>,
    SF: Fn(P, &Header, &mut Vec<u8>) -> io::Result<u64>,
{
    match step {
        PipelineStep::Read => (fastq_try_step_read(state, worker), false),
        PipelineStep::Decompress => (fastq_try_step_decompress(state, worker), false),
        PipelineStep::FindBoundaries => {
            // When parallel parse is enabled, FindBoundaries scans for record boundaries
            if state.config.use_parallel_parse {
                // Uses held-boundaries pattern for parallelism (see base.rs HasHeldBoundaries)
                fastq_try_step_find_boundaries(state, worker)
            } else {
                (false, false) // Not applicable when parallel parse is disabled
            }
        }
        PipelineStep::Decode => {
            // When parallel parse is enabled, Decode step is used for parsing FASTQ records
            if state.config.use_parallel_parse {
                let success = fastq_try_step_parse(state);
                // Parse is parallel - no contention tracking
                (success, false)
            } else {
                (false, false) // Not applicable when parallel parse is disabled
            }
        }
        PipelineStep::Group => {
            if state.config.use_parallel_parse {
                // Use the new group step that receives pre-parsed records
                let success = fastq_try_step_group_parsed(state);
                // Group is exclusive - contention if we couldn't get the lock
                (success, !success && !state.q2_75_parsed.is_empty())
            } else {
                // Use the original group step that does parsing under the lock
                let success = fastq_try_step_group(state);
                // Group is exclusive - contention if we couldn't get the lock
                (success, !success && !state.q2_decompressed.is_empty())
            }
        }
        PipelineStep::Process => (fastq_try_step_process(state, process_fn, worker), false),
        PipelineStep::Serialize => {
            (fastq_try_step_serialize(state, serialize_fn, header, worker), false)
        }
        PipelineStep::Compress => (fastq_try_step_compress(state, worker), false),
        PipelineStep::Write => {
            let success = fastq_try_step_write(state);
            // Write is exclusive - contention if we couldn't get the lock
            (success, !success && !state.output.compressed.is_empty())
        }
    }
}

// ============================================================================
// Step Context (for consolidated generic_worker_loop)
// ============================================================================

/// Context for FASTQ pipeline step execution.
///
/// This struct holds references to all the state needed to execute pipeline steps,
/// and implements `StepContext` to work with `generic_worker_loop`.
pub struct FastqStepContext<'a, R: BufRead + Send, P: Send + MemoryEstimate, PF, SF> {
    pub state: &'a FastqPipelineState<R, P>,
    pub header: &'a Header,
    pub process_fn: &'a PF,
    pub serialize_fn: &'a SF,
    pub is_reader: bool,
}

impl<R, P, PF, SF> StepContext for FastqStepContext<'_, R, P, PF, SF>
where
    R: BufRead + Send + 'static,
    P: Send + MemoryEstimate + 'static,
    PF: Fn(FastqTemplate) -> io::Result<P>,
    SF: Fn(P, &Header, &mut Vec<u8>) -> io::Result<u64>,
{
    type Worker = FastqWorkerState<P>;

    fn execute_step(&self, worker: &mut Self::Worker, step: PipelineStep) -> (bool, bool) {
        fastq_execute_step(
            self.state,
            self.header,
            self.process_fn,
            self.serialize_fn,
            worker,
            step,
        )
    }

    fn get_backpressure(&self, _worker: &Self::Worker) -> BackpressureState {
        let cap = self.state.config.queue_capacity;
        let read_done = self.state.read_done.load(Ordering::Relaxed);
        BackpressureState {
            output_high: self.state.output.compressed.len() > cap * 3 / 4,
            input_low: self.state.q1_raw.len() < cap / 4,
            read_done,
            memory_high: !self.state.is_draining() && self.state.is_memory_high(),
            memory_drained: self.state.is_memory_drained(),
        }
    }

    fn check_drain_mode(&self) {
        let read_done = self.state.read_done.load(Ordering::Relaxed);
        if read_done && self.state.q1_raw.is_empty() {
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

    // check_completion_at_end defaults to false (original FASTQ behavior)

    fn should_attempt_sticky_read(&self) -> bool {
        self.is_reader
    }

    fn sticky_read_should_continue(&self) -> bool {
        // FASTQ just loops until read fails - no pre-condition
        true
    }

    fn execute_read_step(&self, worker: &mut Self::Worker) -> bool {
        fastq_try_step_read(self.state, worker)
    }

    fn exclusive_step_owned(&self, worker: &Self::Worker) -> Option<PipelineStep> {
        if self.is_reader {
            // Reader thread doesn't use the "try owned first" pattern
            // (it has sticky read instead)
            None
        } else {
            // Non-reader workers may own an exclusive step
            worker.core.scheduler.exclusive_step_owned()
        }
    }
}

/// Run the unified 7-step pipeline for FASTQ → BAM conversion.
///
/// This pipeline works for BGZF, Gzip, and Plain FASTQ inputs.
///
/// # Type Parameters
/// * `P` - Processed type (output of `process_fn`, input to `serialize_fn`)
/// * `PF` - Process function type
/// * `SF` - Serialize function type
///
/// # Arguments
/// * `config` - Pipeline configuration
/// * `fastq_paths` - Paths to input FASTQ files (for BGZF inputs)
/// * `decompressed_readers` - Pre-opened readers (for Gzip/Plain inputs)
/// * `header` - BAM header for output
/// * `output` - Output writer
/// * `process_fn` - Function to convert `FastqTemplate` → P
/// * `serialize_fn` - Function to convert P → `SerializedBatch`
///
/// # Returns
/// Number of templates written, or an error.
///
/// # Errors
///
/// Returns an I/O error if any pipeline step or file I/O fails.
#[allow(clippy::too_many_lines)]
pub fn run_fastq_pipeline<P, PF, SF>(
    mut config: FastqPipelineConfig,
    fastq_paths: &[PathBuf],
    decompressed_readers: Option<Vec<Box<dyn BufRead + Send>>>,
    header: &Header,
    mut output: Box<dyn Write + Send>,
    process_fn: PF,
    serialize_fn: SF,
) -> io::Result<u64>
where
    P: Send + MemoryEstimate + 'static,
    PF: Fn(FastqTemplate) -> io::Result<P> + Send + Sync + 'static,
    SF: Fn(P, &Header, &mut Vec<u8>) -> io::Result<u64> + Send + Sync + 'static,
{
    log::debug!("run_fastq_pipeline: starting, num_threads={}", config.num_threads);

    // Synchronized mode requires parallel parse: Decode creates templates directly,
    // and FindBoundaries is needed for BGZF inputs. This is a safety net for callers
    // that set config fields directly; active_steps() also enforces this invariant.
    if config.synchronized {
        config.use_parallel_parse = true;
    }

    // Write BAM header first
    log::debug!("run_fastq_pipeline: writing BAM header");
    write_bam_header(&mut output, header)?;
    log::debug!("run_fastq_pipeline: BAM header written successfully");

    // Create reader based on input type
    log::debug!(
        "run_fastq_pipeline: creating reader, decompressed_readers.is_some()={}, inputs_are_bgzf={}, synchronized={}",
        decompressed_readers.is_some(),
        config.inputs_are_bgzf,
        config.synchronized,
    );
    let reader: FastqReader<Box<dyn BufRead + Send>> = if let Some(readers) = decompressed_readers {
        if config.synchronized {
            // Gzip+synchronized: use record-count reader (skips Decompress + FindBoundaries)
            log::debug!(
                "run_fastq_pipeline: using RecordCount path (synchronized), {} readers",
                readers.len()
            );
            FastqReader::RecordCount(MultiRecordCountReader::new(readers))
        } else {
            // Gzip/Plain path: byte-chunk reader
            log::debug!("run_fastq_pipeline: using Decompressed path, {} readers", readers.len());
            FastqReader::Decompressed(MultiDecompressedReader::new(readers))
        }
    } else {
        // BGZF path: block-level reader (unified through same pipeline)
        log::debug!("run_fastq_pipeline: using BGZF path");
        let bgzf_reader = MultiBgzfBlockReader::new(fastq_paths)?;
        FastqReader::Bgzf(bgzf_reader)
    };

    // Create state
    log::debug!("run_fastq_pipeline: creating pipeline state");
    let state = Arc::new(FastqPipelineState::<Box<dyn BufRead + Send>, P>::new(
        config.clone(),
        reader,
        output,
    ));
    log::debug!("run_fastq_pipeline: state created, spawning {} workers", config.num_threads);

    let process_fn = Arc::new(process_fn);
    let serialize_fn = Arc::new(serialize_fn);
    let header = Arc::new(header.clone());
    let num_threads = config.num_threads;
    let compression_level = config.compression_level;
    let scheduler_strategy = config.scheduler_strategy;
    let active_steps = config.active_steps();

    // Spawn workers
    log::debug!("run_fastq_pipeline: spawning {num_threads} worker threads");
    let handles: Vec<_> = (0..num_threads)
        .map(|thread_id| {
            let state = Arc::clone(&state);
            let process_fn = Arc::clone(&process_fn);
            let serialize_fn = Arc::clone(&serialize_fn);
            let header = Arc::clone(&header);
            let active_steps = active_steps.clone();

            thread::spawn(move || {
                // Wrap worker logic in catch_unwind to handle panics gracefully
                let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    log::debug!("Worker thread {thread_id} starting");
                    let mut worker = FastqWorkerState::new(
                        compression_level,
                        thread_id,
                        num_threads,
                        scheduler_strategy,
                        active_steps,
                    );
                    log::debug!("Worker thread {thread_id} created worker state");
                    let ctx = FastqStepContext {
                        state: &state,
                        header: &header,
                        process_fn: &*process_fn,
                        serialize_fn: &*serialize_fn,
                        is_reader: thread_id == 0,
                    };
                    generic_worker_loop(&ctx, &mut worker);
                    log::debug!("Worker thread {thread_id} finished");
                }));

                // If a panic occurred, set the error flag so other threads exit
                if let Err(panic_info) = result {
                    handle_worker_panic(&*state, thread_id, panic_info);
                }
            })
        })
        .collect();
    log::debug!("run_fastq_pipeline: all workers spawned");

    // Spawn monitor thread for deadlock detection
    let monitor_handle = if state.stats().is_some() || state.deadlock_state.is_enabled() {
        let state_clone = Arc::clone(&state);
        Some(thread::spawn(move || {
            // Use shared monitor loop: 100ms sample interval, deadlock check every 10 samples (~1s)
            run_monitor_loop(&state_clone, 100, 10, |s| {
                // Log parallel parse state for debugging (at trace level)
                if s.config.use_parallel_parse && s.deadlock_state.is_enabled() {
                    let bd = s.boundaries_done.load(Ordering::Relaxed);
                    let pd = s.parse_done.load(Ordering::Relaxed);
                    let br = s.batches_read.load(Ordering::Relaxed);
                    let bf = s.batches_boundaries_found.load(Ordering::Relaxed);
                    let bp = s.batches_parsed.load(Ordering::Relaxed);
                    let bg = s.batches_grouped.load(Ordering::Relaxed);
                    log::trace!(
                        "Parallel parse state: boundaries_done={bd}, parse_done={pd}, batches: read={br}, boundaries={bf}, parsed={bp}, grouped={bg}"
                    );
                }
            });
        }))
    } else {
        None
    };

    // Wait for completion
    log::debug!("run_fastq_pipeline: waiting for workers to complete");
    join_worker_threads(handles)?;
    log::debug!("run_fastq_pipeline: all workers joined");
    join_monitor_thread(monitor_handle);

    // Debug: Log pipeline counters
    log::info!(
        "Pipeline counters: batches_read={}, batches_grouped={}, total_templates_pushed={}, total_records_serialized={}, templates_written={}",
        state.batches_read.load(Ordering::Relaxed),
        state.batches_grouped.load(Ordering::Relaxed),
        state.total_templates_pushed.load(Ordering::Relaxed),
        state.total_records_serialized.load(Ordering::Relaxed),
        state.output.items_written.load(Ordering::Relaxed)
    );

    // Finalize: check errors, flush output (with BGZF EOF), log stats
    finalize_pipeline(&*state)
}

/// Write BAM header to output.
fn write_bam_header(writer: &mut dyn Write, header: &Header) -> io::Result<()> {
    use noodles::bam;

    let mut encoder = bam::io::Writer::new(writer);
    encoder.write_header(header)?;
    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::super::bam::*;
    use super::super::base::*;
    use super::*;
    use crate::bgzf_reader::RawBgzfBlock;
    use crate::bgzf_writer::CompressedBlock;
    use PipelineStep::*;
    use rstest::rstest;
    use std::io::Cursor;

    #[test]
    fn test_raw_block_batch() {
        let mut batch = RawBlockBatch::new();
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);

        // Add a mock block
        batch.blocks.push(RawBgzfBlock { data: vec![0u8; 100] });
        assert!(!batch.is_empty());
        assert_eq!(batch.len(), 1);

        // Clear
        batch.clear();
        assert!(batch.is_empty());
    }

    #[test]
    fn test_compressed_block_batch() {
        let mut batch = CompressedBlockBatch::new();
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);

        // Add a block
        batch.blocks.push(CompressedBlock { serial: 0, data: vec![1, 2, 3, 4, 5] });
        assert!(!batch.is_empty());
        assert_eq!(batch.len(), 1);
        assert_eq!(batch.total_size(), 5);

        // Clear
        batch.clear();
        assert!(batch.is_empty());
    }

    #[test]
    fn test_bgzf_batch_config() {
        let config = BgzfBatchConfig::default();
        assert_eq!(config.blocks_per_batch, 16);
        assert_eq!(config.compression_level, 6);

        let config = BgzfBatchConfig::new(32).with_compression_level(9);
        assert_eq!(config.blocks_per_batch, 32);
        assert_eq!(config.compression_level, 9);
    }

    // ========================================================================
    // 7-Step Pipeline Tests
    // ========================================================================

    #[test]
    fn test_pipeline_config() {
        let config = PipelineConfig::new(4, 6);
        assert_eq!(config.num_threads, 4);
        assert_eq!(config.compression_level, 6);

        let config = PipelineConfig::new(8, 6).with_compression_level(9);
        assert_eq!(config.num_threads, 8);
        assert_eq!(config.compression_level, 9);
    }

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
        assert_eq!(all[8], PipelineStep::Write);
    }

    #[test]
    fn test_decompressed_batch() {
        let mut batch = DecompressedBatch::new();
        assert!(batch.is_empty());

        batch.data.extend_from_slice(b"hello");
        assert!(!batch.is_empty());

        batch.clear();
        assert!(batch.is_empty());
    }

    #[test]
    fn test_serialized_batch() {
        let mut batch = SerializedBatch::new();
        assert!(batch.is_empty());

        batch.data.extend_from_slice(b"data");
        assert!(!batch.is_empty());

        batch.clear();
        assert!(batch.is_empty());
    }

    #[test]
    fn test_bam_pipeline_config() {
        let config = BamPipelineConfig::new(4, 6);
        assert_eq!(config.pipeline.num_threads, 4);
        assert_eq!(config.compression_level, 6);

        let config = BamPipelineConfig::new(8, 6).with_compression_level(9);
        assert_eq!(config.compression_level, 9);
        assert_eq!(config.pipeline.compression_level, 9);
    }

    #[test]
    fn test_read_raw_blocks_from_memory() {
        // Create a minimal BGZF block in memory
        // This is a compressed empty block
        let bgzf_empty_block: Vec<u8> = vec![
            0x1f, 0x8b, // gzip magic
            0x08, // compression method
            0x04, // flags (FEXTRA)
            0x00, 0x00, 0x00, 0x00, // mtime
            0x00, // extra flags
            0xff, // OS
            0x06, 0x00, // extra length (6)
            0x42, 0x43, // BC subfield
            0x02, 0x00, // subfield length
            0x1b, 0x00, // block size - 1 (27)
            0x03, 0x00, // deflate empty block
            0x00, 0x00, 0x00, 0x00, // CRC32
            0x00, 0x00, 0x00, 0x00, // uncompressed size
        ];

        let mut reader = Cursor::new(bgzf_empty_block);
        let blocks = read_raw_blocks(&mut reader, 10).unwrap();

        // Should read 0 blocks (EOF blocks are skipped by read_raw_blocks)
        assert_eq!(blocks.len(), 0);
    }

    /// Test that demonstrates data loss when using `let _ = queue.push()` pattern.
    ///
    /// This test uses a mock pipeline to show that with a tiny queue capacity,
    /// the old pattern of discarding push failures leads to data loss and
    /// mismatched counters.
    #[test]
    fn test_discarded_push_causes_data_loss() {
        use std::sync::Arc;
        use std::thread;

        // Simulate a producer-consumer with a tiny queue
        let queue: Arc<ArrayQueue<u64>> = Arc::new(ArrayQueue::new(2));
        let items_pushed = Arc::new(AtomicU64::new(0));
        let items_received = Arc::new(AtomicU64::new(0));

        // Producer: push 100 items, discarding failures (OLD buggy pattern)
        let queue_producer = Arc::clone(&queue);
        let pushed = Arc::clone(&items_pushed);
        let producer = thread::spawn(move || {
            for i in 0..100u64 {
                // OLD BUGGY PATTERN: discard push result
                let _ = queue_producer.push(i);
                pushed.fetch_add(1, Ordering::Relaxed);
            }
        });

        // Consumer: slowly drain the queue
        let queue_consumer = Arc::clone(&queue);
        let received = Arc::clone(&items_received);
        let consumer = thread::spawn(move || {
            // Wait for producer to finish pushing
            thread::sleep(std::time::Duration::from_millis(50));
            while queue_consumer.pop().is_some() {
                received.fetch_add(1, Ordering::Relaxed);
            }
        });

        producer.join().unwrap();
        consumer.join().unwrap();

        let pushed_count = items_pushed.load(Ordering::Relaxed);
        let received_count = items_received.load(Ordering::Relaxed);

        // This demonstrates the bug: we "pushed" 100 items but only received 2
        // because the queue was full and we discarded the failures
        assert_eq!(pushed_count, 100, "Producer thought it pushed 100 items");
        assert!(
            received_count < pushed_count,
            "Data was lost! Pushed {pushed_count} but only received {received_count}",
        );
        // With queue capacity 2, we should only receive 2 items
        assert_eq!(received_count, 2, "Only queue capacity items should be received");
    }

    #[test]
    fn test_multi_stream_grouper_add_parsed_batch() {
        // Test the add_parsed_batch method with pre-parsed records
        let mut grouper = FastqMultiStreamGrouper::new(2);

        // Create pre-parsed records for 2 streams (R1 and R2)
        // Both have the same template "read1"
        let streams = vec![
            FastqParsedStream {
                stream_idx: 0,
                records: vec![FastqRecord {
                    name: b"read1".to_vec(),
                    sequence: b"ACGT".to_vec(),
                    quality: b"IIII".to_vec(),
                }],
            },
            FastqParsedStream {
                stream_idx: 1,
                records: vec![FastqRecord {
                    name: b"read1".to_vec(),
                    sequence: b"GGGG".to_vec(),
                    quality: b"JJJJ".to_vec(),
                }],
            },
        ];

        let templates = grouper.add_parsed_batch(streams).unwrap();

        // Should produce one complete template with both reads
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].records.len(), 2);
        assert_eq!(templates[0].records[0].name, b"read1");
        assert_eq!(templates[0].records[1].name, b"read1");
    }

    #[test]
    fn test_multi_stream_grouper_add_parsed_batch_multiple_templates() {
        let mut grouper = FastqMultiStreamGrouper::new(2);

        // Create pre-parsed records for multiple templates
        let streams = vec![
            FastqParsedStream {
                stream_idx: 0,
                records: vec![
                    FastqRecord {
                        name: b"read1".to_vec(),
                        sequence: b"AAAA".to_vec(),
                        quality: b"IIII".to_vec(),
                    },
                    FastqRecord {
                        name: b"read2".to_vec(),
                        sequence: b"CCCC".to_vec(),
                        quality: b"JJJJ".to_vec(),
                    },
                ],
            },
            FastqParsedStream {
                stream_idx: 1,
                records: vec![
                    FastqRecord {
                        name: b"read1".to_vec(),
                        sequence: b"GGGG".to_vec(),
                        quality: b"KKKK".to_vec(),
                    },
                    FastqRecord {
                        name: b"read2".to_vec(),
                        sequence: b"TTTT".to_vec(),
                        quality: b"LLLL".to_vec(),
                    },
                ],
            },
        ];

        let templates = grouper.add_parsed_batch(streams).unwrap();

        // Should produce two complete templates
        assert_eq!(templates.len(), 2);
    }

    #[test]
    fn test_parallel_parse_config_default_disabled() {
        // Parallel parse is disabled by default because overhead hurts t4 performance.
        // It can be enabled via with_parallel_parse(true) for high-thread scenarios.
        let config = FastqPipelineConfig::new(4, false, 6);
        assert!(!config.use_parallel_parse, "Parallel parse should be disabled by default");
    }

    #[test]
    fn test_parallel_parse_boundary_finding_integration() {
        // Test the full boundary finding -> parse flow
        let boundary_state = FastqBoundaryState::new(2);

        // Create a batch with data for 2 streams
        let batch = FastqDecompressedBatch {
            chunks: vec![
                FastqDecompressedChunk {
                    stream_idx: 0,
                    data: b"@read1\nACGT\n+\nIIII\n@read2\nGGGG\n+\nJJJJ\n".to_vec(),
                },
                FastqDecompressedChunk {
                    stream_idx: 1,
                    data: b"@read1\nTTTT\n+\nKKKK\n@read2\nCCCC\n+\nLLLL\n".to_vec(),
                },
            ],
            serial: 0,
        };

        // Step 1: Find boundaries
        let boundary_batch = FastqFormat::find_boundaries(&boundary_state, batch).unwrap();
        assert_eq!(boundary_batch.streams.len(), 2);
        // Each stream should have 2 complete records (offsets at 0, 19, 38)
        assert_eq!(boundary_batch.streams[0].offsets.len(), 3);
        assert_eq!(boundary_batch.streams[1].offsets.len(), 3);

        // Step 2: Parse records
        let parsed_batch = FastqFormat::parse_records(boundary_batch).unwrap();
        assert_eq!(parsed_batch.streams.len(), 2);
        assert_eq!(parsed_batch.streams[0].stream_idx, 0);
        assert_eq!(parsed_batch.streams[1].stream_idx, 1);
        assert_eq!(parsed_batch.streams[0].records.len(), 2);
        assert_eq!(parsed_batch.streams[1].records.len(), 2);

        // Verify record contents
        assert_eq!(parsed_batch.streams[0].records[0].name, b"read1");
        assert_eq!(parsed_batch.streams[0].records[0].sequence, b"ACGT");
        assert_eq!(parsed_batch.streams[0].records[1].name, b"read2");
        assert_eq!(parsed_batch.streams[1].records[0].name, b"read1");
        assert_eq!(parsed_batch.streams[1].records[0].sequence, b"TTTT");

        // Step 3: Group - use the add_parsed_batch method
        let mut grouper = FastqMultiStreamGrouper::new(2);
        let templates = grouper.add_parsed_batch(parsed_batch.streams).unwrap();

        // Should produce 2 complete templates (read1 and read2)
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].records.len(), 2); // R1 and R2 for read1
        assert_eq!(templates[1].records.len(), 2); // R1 and R2 for read2
    }

    #[test]
    fn test_parallel_parse_records_spanning_chunks() {
        // Test records that span chunk boundaries
        let boundary_state = FastqBoundaryState::new(1);

        // First chunk: one complete record + incomplete record
        let batch1 = FastqDecompressedBatch {
            chunks: vec![FastqDecompressedChunk {
                stream_idx: 0,
                data: b"@read1\nACGT\n+\nIIII\n@read2\nGG".to_vec(),
            }],
            serial: 0,
        };

        let boundary_batch1 = FastqFormat::find_boundaries(&boundary_state, batch1).unwrap();
        assert_eq!(boundary_batch1.streams[0].offsets.len(), 2); // One complete record
        assert!(!boundary_state.stream_states[0].lock().leftover.is_empty()); // Leftover from incomplete

        // Second chunk: completes the record
        let batch2 = FastqDecompressedBatch {
            chunks: vec![FastqDecompressedChunk { stream_idx: 0, data: b"GG\n+\nJJJJ\n".to_vec() }],
            serial: 1,
        };

        let boundary_batch2 = FastqFormat::find_boundaries(&boundary_state, batch2).unwrap();
        // Leftover + new data should form complete record
        assert!(boundary_batch2.streams[0].offsets.len() >= 2);
        assert!(boundary_state.stream_states[0].lock().leftover.is_empty()); // No more leftover

        // Parse both batches
        let parsed1 = FastqFormat::parse_records(boundary_batch1).unwrap();
        let parsed2 = FastqFormat::parse_records(boundary_batch2).unwrap();

        assert_eq!(parsed1.streams[0].records.len(), 1);
        assert_eq!(parsed1.streams[0].records[0].name, b"read1");
        assert_eq!(parsed2.streams[0].records.len(), 1);
        assert_eq!(parsed2.streams[0].records[0].name, b"read2");
    }

    #[test]
    fn test_parallel_parse_thread_safety() {
        use std::thread;

        // Test that concurrent parsing produces consistent results
        // Create multiple boundary batches that can be parsed in parallel
        let num_threads = 4;
        let batches_per_thread = 10;

        let results: Vec<_> = (0..num_threads)
            .map(|thread_id| {
                thread::spawn(move || {
                    let mut records_parsed = 0;
                    for batch_id in 0..batches_per_thread {
                        // Create a boundary batch for this thread/batch
                        let name = format!("read_t{thread_id}_b{batch_id}");
                        let data = format!("@{name}\nACGT\n+\nIIII\n");

                        let boundary_batch = FastqBoundaryBatch {
                            streams: vec![FastqStreamBoundaries {
                                stream_idx: 0,
                                data: data.as_bytes().to_vec(),
                                offsets: vec![0, data.len()],
                            }],
                            serial: (thread_id * batches_per_thread + batch_id) as u64,
                        };

                        // Parse the batch
                        let parsed = FastqFormat::parse_records(boundary_batch).unwrap();
                        assert_eq!(parsed.streams[0].stream_idx, 0);
                        assert_eq!(parsed.streams[0].records.len(), 1);
                        assert_eq!(
                            String::from_utf8_lossy(&parsed.streams[0].records[0].name),
                            name
                        );
                        records_parsed += 1;
                    }
                    records_parsed
                })
            })
            .collect();

        // Verify all threads completed successfully
        let total_parsed: usize =
            results.into_iter().map(|h| h.join().expect("Thread panicked")).sum();

        assert_eq!(total_parsed, num_threads * batches_per_thread, "All records should be parsed");
    }

    #[test]
    fn test_grouper_thread_safe_access() {
        use parking_lot::Mutex;
        use std::sync::Arc;
        use std::thread;

        // Test that the grouper produces consistent results when accessed
        // through a mutex (simulating the pipeline's locking pattern)
        let grouper = Arc::new(Mutex::new(FastqMultiStreamGrouper::new(2)));
        let num_threads = 4;
        let records_per_thread = 10;

        let handles: Vec<_> = (0..num_threads)
            .map(|thread_id| {
                let grouper = Arc::clone(&grouper);
                thread::spawn(move || {
                    for i in 0..records_per_thread {
                        let name = format!("read_t{thread_id}_i{i}");
                        let streams = vec![
                            FastqParsedStream {
                                stream_idx: 0,
                                records: vec![FastqRecord {
                                    name: name.as_bytes().to_vec(),
                                    sequence: b"ACGT".to_vec(),
                                    quality: b"IIII".to_vec(),
                                }],
                            },
                            FastqParsedStream {
                                stream_idx: 1,
                                records: vec![FastqRecord {
                                    name: name.as_bytes().to_vec(),
                                    sequence: b"TTTT".to_vec(),
                                    quality: b"JJJJ".to_vec(),
                                }],
                            },
                        ];

                        let templates = grouper.lock().add_parsed_batch(streams).unwrap();
                        // Each batch should produce exactly one template
                        assert_eq!(templates.len(), 1);
                        assert_eq!(templates[0].records.len(), 2);
                    }
                })
            })
            .collect();

        // Wait for all threads
        for h in handles {
            h.join().expect("Thread panicked");
        }

        // Verify grouper state is clean
        let grouper = grouper.lock();
        assert!(!grouper.has_pending(), "Grouper should have no pending data");
    }

    // ========================================================================
    // Synchronized Stream Record Alignment Tests
    // ========================================================================

    #[test]
    fn test_find_boundaries_aligns_unequal_record_counts() {
        // Test that find_boundaries aligns record counts when streams have different
        // numbers of complete records in their chunks.
        let boundary_state = FastqBoundaryState::new(2);

        // Stream 0 has 3 complete records, Stream 1 has 2 complete records
        let batch = FastqDecompressedBatch {
            chunks: vec![
                FastqDecompressedChunk {
                    stream_idx: 0,
                    data: b"@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n@r3\nACGT\n+\nIIII\n".to_vec(),
                },
                FastqDecompressedChunk {
                    stream_idx: 1,
                    data: b"@r1\nTTTT\n+\nJJJJ\n@r2\nTTTT\n+\nJJJJ\n".to_vec(),
                },
            ],
            serial: 0,
        };

        let boundary_batch = FastqFormat::find_boundaries(&boundary_state, batch).unwrap();

        // Both streams should have exactly 2 records (the minimum)
        assert_eq!(boundary_batch.streams.len(), 2);
        // offsets.len() = record_count + 1 (includes position 0)
        assert_eq!(
            boundary_batch.streams[0].offsets.len(),
            3,
            "Stream 0 should have 2 records (3 offsets)"
        );
        assert_eq!(
            boundary_batch.streams[1].offsets.len(),
            3,
            "Stream 1 should have 2 records (3 offsets)"
        );

        // Stream 0's excess record should be in leftover
        let leftover = &boundary_state.stream_states[0].lock().leftover;
        assert!(!leftover.is_empty(), "Stream 0 should have leftover containing the excess record");
        assert!(leftover.starts_with(b"@r3\n"), "Leftover should contain the third record");

        // Stream 1's leftover should be empty (no excess)
        let leftover1 = &boundary_state.stream_states[1].lock().leftover;
        assert!(leftover1.is_empty(), "Stream 1 should have no leftover");
    }

    #[test]
    fn test_find_boundaries_leftover_persists_to_next_batch() {
        // Test that excess records moved to leftover are correctly processed
        // in the next batch.
        let boundary_state = FastqBoundaryState::new(2);

        // Batch 1: Stream 0 has 3 records, Stream 1 has 2 records
        let batch1 = FastqDecompressedBatch {
            chunks: vec![
                FastqDecompressedChunk {
                    stream_idx: 0,
                    data: b"@r1\nAAAA\n+\nIIII\n@r2\nAAAA\n+\nIIII\n@r3\nAAAA\n+\nIIII\n".to_vec(),
                },
                FastqDecompressedChunk {
                    stream_idx: 1,
                    data: b"@r1\nTTTT\n+\nJJJJ\n@r2\nTTTT\n+\nJJJJ\n".to_vec(),
                },
            ],
            serial: 0,
        };

        let boundary_batch1 = FastqFormat::find_boundaries(&boundary_state, batch1).unwrap();
        assert_eq!(boundary_batch1.streams[0].offsets.len() - 1, 2); // 2 records from stream 0
        assert_eq!(boundary_batch1.streams[1].offsets.len() - 1, 2); // 2 records from stream 1

        // Batch 2: Stream 0 has 1 record, Stream 1 has 2 records
        // Stream 0's leftover (r3) + new record (r4) = 2 records total
        let batch2 = FastqDecompressedBatch {
            chunks: vec![
                FastqDecompressedChunk { stream_idx: 0, data: b"@r4\nAAAA\n+\nIIII\n".to_vec() },
                FastqDecompressedChunk {
                    stream_idx: 1,
                    data: b"@r3\nTTTT\n+\nJJJJ\n@r4\nTTTT\n+\nJJJJ\n".to_vec(),
                },
            ],
            serial: 1,
        };

        let boundary_batch2 = FastqFormat::find_boundaries(&boundary_state, batch2).unwrap();

        // Stream 0: leftover(r3) + new(r4) = 2 records
        // Stream 1: 2 new records (r3, r4)
        // Both should have 2 records
        assert_eq!(
            boundary_batch2.streams[0].offsets.len() - 1,
            2,
            "Stream 0 should have 2 records (leftover + new)"
        );
        assert_eq!(
            boundary_batch2.streams[1].offsets.len() - 1,
            2,
            "Stream 1 should have 2 records"
        );

        // Parse and verify record names
        let parsed = FastqFormat::parse_records(boundary_batch2).unwrap();
        assert_eq!(
            parsed.streams[0].records[0].name, b"r3",
            "First record should be r3 from leftover"
        );
        assert_eq!(parsed.streams[0].records[1].name, b"r4", "Second record should be r4");
    }

    #[test]
    fn test_find_boundaries_processes_leftover_without_new_chunk() {
        // Test that when one stream has no new chunk but has leftover,
        // the leftover is still processed.
        let boundary_state = FastqBoundaryState::new(2);

        // Batch 1: Both streams have data, but stream 0 has an extra record
        let batch1 = FastqDecompressedBatch {
            chunks: vec![
                FastqDecompressedChunk {
                    stream_idx: 0,
                    data: b"@r1\nAAAA\n+\nIIII\n@r2\nAAAA\n+\nIIII\n".to_vec(),
                },
                FastqDecompressedChunk { stream_idx: 1, data: b"@r1\nTTTT\n+\nJJJJ\n".to_vec() },
            ],
            serial: 0,
        };

        let boundary_batch1 = FastqFormat::find_boundaries(&boundary_state, batch1).unwrap();
        // Both aligned to 1 record
        assert_eq!(boundary_batch1.streams[0].offsets.len() - 1, 1);
        assert_eq!(boundary_batch1.streams[1].offsets.len() - 1, 1);

        // Stream 0's leftover should have r2
        assert!(!boundary_state.stream_states[0].lock().leftover.is_empty());

        // Batch 2: Only stream 1 has a new chunk (simulating stream 0 at EOF)
        let batch2 = FastqDecompressedBatch {
            chunks: vec![FastqDecompressedChunk {
                stream_idx: 1,
                data: b"@r2\nTTTT\n+\nJJJJ\n".to_vec(),
            }],
            serial: 1,
        };

        let boundary_batch2 = FastqFormat::find_boundaries(&boundary_state, batch2).unwrap();

        // Both streams should be present (stream 0 from leftover, stream 1 from new chunk)
        assert_eq!(boundary_batch2.streams.len(), 2, "Both streams should be present");

        // Find stream 0 and stream 1 in the result (order may vary)
        let stream0 = boundary_batch2.streams.iter().find(|s| s.stream_idx == 0).unwrap();
        let stream1 = boundary_batch2.streams.iter().find(|s| s.stream_idx == 1).unwrap();

        assert_eq!(stream0.offsets.len() - 1, 1, "Stream 0 should have 1 record from leftover");
        assert_eq!(stream1.offsets.len() - 1, 1, "Stream 1 should have 1 record");

        // Stream 0's leftover should now be empty
        assert!(
            boundary_state.stream_states[0].lock().leftover.is_empty(),
            "Stream 0 leftover should be consumed"
        );
    }

    #[test]
    fn test_find_boundaries_equal_counts_no_alignment_needed() {
        // Test that when both streams have equal record counts, no alignment is needed
        let boundary_state = FastqBoundaryState::new(2);

        let batch = FastqDecompressedBatch {
            chunks: vec![
                FastqDecompressedChunk {
                    stream_idx: 0,
                    data: b"@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n".to_vec(),
                },
                FastqDecompressedChunk {
                    stream_idx: 1,
                    data: b"@r1\nTTTT\n+\nJJJJ\n@r2\nTTTT\n+\nJJJJ\n".to_vec(),
                },
            ],
            serial: 0,
        };

        let boundary_batch = FastqFormat::find_boundaries(&boundary_state, batch).unwrap();

        // Both streams should have 2 records
        assert_eq!(boundary_batch.streams[0].offsets.len() - 1, 2);
        assert_eq!(boundary_batch.streams[1].offsets.len() - 1, 2);

        // No leftover for either stream
        assert!(boundary_state.stream_states[0].lock().leftover.is_empty());
        assert!(boundary_state.stream_states[1].lock().leftover.is_empty());
    }

    #[test]
    fn test_find_boundaries_single_stream_no_alignment() {
        // Test that single-stream mode doesn't apply alignment (no need)
        let boundary_state = FastqBoundaryState::new(1);

        let batch = FastqDecompressedBatch {
            chunks: vec![FastqDecompressedChunk {
                stream_idx: 0,
                data: b"@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n@r3\nACGT\n+\nIIII\n".to_vec(),
            }],
            serial: 0,
        };

        let boundary_batch = FastqFormat::find_boundaries(&boundary_state, batch).unwrap();

        // All 3 records should be present (no alignment needed for single stream)
        assert_eq!(boundary_batch.streams[0].offsets.len() - 1, 3);
        assert!(boundary_state.stream_states[0].lock().leftover.is_empty());
    }

    #[test]
    fn test_find_boundaries_empty_batch_with_leftover() {
        // Test handling when a batch has no chunks but streams have leftover
        let boundary_state = FastqBoundaryState::new(2);

        // First, create leftover by processing a batch with unequal records
        let batch1 = FastqDecompressedBatch {
            chunks: vec![
                FastqDecompressedChunk {
                    stream_idx: 0,
                    data: b"@r1\nAAAA\n+\nIIII\n@r2\nAAAA\n+\nIIII\n".to_vec(),
                },
                FastqDecompressedChunk { stream_idx: 1, data: b"@r1\nTTTT\n+\nJJJJ\n".to_vec() },
            ],
            serial: 0,
        };

        let _ = FastqFormat::find_boundaries(&boundary_state, batch1).unwrap();

        // Verify stream 0 has leftover
        assert!(!boundary_state.stream_states[0].lock().leftover.is_empty());

        // Now process a batch with only stream 1 data
        let batch2 = FastqDecompressedBatch {
            chunks: vec![FastqDecompressedChunk {
                stream_idx: 1,
                data: b"@r2\nTTTT\n+\nJJJJ\n".to_vec(),
            }],
            serial: 1,
        };

        let boundary_batch2 = FastqFormat::find_boundaries(&boundary_state, batch2).unwrap();

        // Both streams should be present
        assert_eq!(boundary_batch2.streams.len(), 2);

        // After processing, stream 0's leftover should be consumed
        assert!(boundary_state.stream_states[0].lock().leftover.is_empty());
    }

    // ========================================================================
    // Stream Identity Preservation Tests (regression for R1/R2 swap bug)
    // ========================================================================

    #[test]
    fn test_parse_records_preserves_stream_idx() {
        // Verify that parse_records carries stream_idx from boundary batch
        // to parsed batch, even when streams are in non-sequential order.
        let boundary_batch = FastqBoundaryBatch {
            streams: vec![
                FastqStreamBoundaries {
                    stream_idx: 1, // R2 first (reversed order)
                    data: b"@read1\nTTTT\n+\nJJJJ\n".to_vec(),
                    offsets: vec![0, 20],
                },
                FastqStreamBoundaries {
                    stream_idx: 0, // R1 second
                    data: b"@read1\nACGT\n+\nIIII\n".to_vec(),
                    offsets: vec![0, 20],
                },
            ],
            serial: 42,
        };

        let parsed = FastqFormat::parse_records(boundary_batch).unwrap();
        assert_eq!(parsed.serial, 42);
        assert_eq!(parsed.streams.len(), 2);
        // stream_idx must be preserved, not assumed from position
        assert_eq!(
            parsed.streams[0].stream_idx, 1,
            "First parsed stream should be R2 (stream_idx=1)"
        );
        assert_eq!(
            parsed.streams[1].stream_idx, 0,
            "Second parsed stream should be R1 (stream_idx=0)"
        );
        assert_eq!(parsed.streams[0].records[0].sequence, b"TTTT");
        assert_eq!(parsed.streams[1].records[0].sequence, b"ACGT");
    }

    #[test]
    fn test_create_templates_from_reversed_streams() {
        // Regression test: when find_boundaries() produces streams in reversed
        // order [R2, R1] near EOF, create_templates_from_streams() must still
        // produce templates with records in the correct R1, R2 order.
        //
        // Before the fix, this would swap R1 and R2 data, causing read
        // structures to be applied to the wrong reads.
        let streams = vec![
            FastqParsedStream {
                stream_idx: 1, // R2 comes first in the Vec (reversed!)
                records: vec![FastqRecord {
                    name: b"read1".to_vec(),
                    sequence: b"TTTT".to_vec(), // R2 sequence
                    quality: b"JJJJ".to_vec(),
                }],
            },
            FastqParsedStream {
                stream_idx: 0, // R1 comes second
                records: vec![FastqRecord {
                    name: b"read1".to_vec(),
                    sequence: b"ACGT".to_vec(), // R1 sequence
                    quality: b"IIII".to_vec(),
                }],
            },
        ];

        let templates = create_templates_from_streams(streams).unwrap();

        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].records.len(), 2);
        // Critical: R1 must be first (records[0]), R2 must be second (records[1])
        assert_eq!(
            templates[0].records[0].sequence, b"ACGT",
            "First record in template must be R1 (ACGT), not R2"
        );
        assert_eq!(
            templates[0].records[1].sequence, b"TTTT",
            "Second record in template must be R2 (TTTT), not R1"
        );
    }

    #[test]
    fn test_create_templates_from_correctly_ordered_streams() {
        // Verify correct behavior when streams are already in order
        let streams = vec![
            FastqParsedStream {
                stream_idx: 0,
                records: vec![FastqRecord {
                    name: b"read1".to_vec(),
                    sequence: b"ACGT".to_vec(),
                    quality: b"IIII".to_vec(),
                }],
            },
            FastqParsedStream {
                stream_idx: 1,
                records: vec![FastqRecord {
                    name: b"read1".to_vec(),
                    sequence: b"TTTT".to_vec(),
                    quality: b"JJJJ".to_vec(),
                }],
            },
        ];

        let templates = create_templates_from_streams(streams).unwrap();

        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].records[0].sequence, b"ACGT", "R1 should be first");
        assert_eq!(templates[0].records[1].sequence, b"TTTT", "R2 should be second");
    }

    #[test]
    fn test_add_parsed_batch_with_reversed_streams() {
        // Regression test: add_parsed_batch must use stream_idx, not positional
        // index, when feeding records to the grouper.
        let mut grouper = FastqMultiStreamGrouper::new(2);

        // Streams in reversed order (R2 first, R1 second)
        let streams = vec![
            FastqParsedStream {
                stream_idx: 1, // R2 first
                records: vec![FastqRecord {
                    name: b"read1".to_vec(),
                    sequence: b"TTTT".to_vec(),
                    quality: b"JJJJ".to_vec(),
                }],
            },
            FastqParsedStream {
                stream_idx: 0, // R1 second
                records: vec![FastqRecord {
                    name: b"read1".to_vec(),
                    sequence: b"ACGT".to_vec(),
                    quality: b"IIII".to_vec(),
                }],
            },
        ];

        let templates = grouper.add_parsed_batch(streams).unwrap();

        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].records.len(), 2);
        // The grouper stores records by stream_idx, so records[0] = stream 0 (R1),
        // records[1] = stream 1 (R2)
        assert_eq!(
            templates[0].records[0].sequence, b"ACGT",
            "R1 record should be first in template"
        );
        assert_eq!(
            templates[0].records[1].sequence, b"TTTT",
            "R2 record should be second in template"
        );
    }

    #[test]
    fn test_end_to_end_reversed_stream_order_at_eof() {
        // End-to-end regression test simulating the EOF boundary condition.
        //
        // This reproduces the exact bug scenario: when one stream reaches EOF
        // before the other, find_boundaries() processes the leftover stream
        // after the non-EOF stream, producing a reversed streams Vec.
        // The full pipeline (find_boundaries -> parse_records ->
        // create_templates_from_streams) must handle this correctly.
        let boundary_state = FastqBoundaryState::new(2);

        // Batch 1: Both streams have data, but stream 0 has an extra record.
        // This creates leftover for stream 0.
        let batch1 = FastqDecompressedBatch {
            chunks: vec![
                FastqDecompressedChunk {
                    stream_idx: 0,
                    data: b"@read1\nACGT\n+\nIIII\n@read2\nGGGG\n+\nJJJJ\n".to_vec(),
                },
                FastqDecompressedChunk { stream_idx: 1, data: b"@read1\nTTTT\n+\nKKKK\n".to_vec() },
            ],
            serial: 0,
        };

        let boundary_batch1 = FastqFormat::find_boundaries(&boundary_state, batch1).unwrap();
        // Both aligned to 1 record; stream 0 has leftover
        assert_eq!(boundary_batch1.streams[0].offsets.len() - 1, 1);
        assert_eq!(boundary_batch1.streams[1].offsets.len() - 1, 1);

        // Batch 2: Only stream 1 has a new chunk (stream 0 at EOF).
        // Stream 0's leftover is processed in the "missing chunk" path,
        // which appends AFTER stream 1 — producing reversed order.
        let batch2 = FastqDecompressedBatch {
            chunks: vec![FastqDecompressedChunk {
                stream_idx: 1,
                data: b"@read2\nCCCC\n+\nLLLL\n".to_vec(),
            }],
            serial: 1,
        };

        let boundary_batch2 = FastqFormat::find_boundaries(&boundary_state, batch2).unwrap();
        assert_eq!(boundary_batch2.streams.len(), 2);

        // The order of streams in boundary_batch2 may be [stream_idx=1, stream_idx=0]
        // (reversed). This is the trigger for the bug.
        let first_stream_idx = boundary_batch2.streams[0].stream_idx;
        let second_stream_idx = boundary_batch2.streams[1].stream_idx;

        // Parse records — stream_idx must be preserved
        let parsed = FastqFormat::parse_records(boundary_batch2).unwrap();
        assert_eq!(parsed.streams[0].stream_idx, first_stream_idx);
        assert_eq!(parsed.streams[1].stream_idx, second_stream_idx);

        // Create templates — must produce correct R1/R2 ordering regardless
        // of the stream order in the Vec
        let templates = create_templates_from_streams(parsed.streams).unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, b"read2");
        assert_eq!(templates[0].records.len(), 2);

        // THE KEY ASSERTION: R1 data (GGGG from leftover) must be records[0],
        // R2 data (CCCC from new chunk) must be records[1].
        assert_eq!(
            templates[0].records[0].sequence, b"GGGG",
            "records[0] must be R1 (stream 0) data, not R2"
        );
        assert_eq!(
            templates[0].records[1].sequence, b"CCCC",
            "records[1] must be R2 (stream 1) data, not R1"
        );
    }

    // ========================================================================
    // MultiRecordCountReader tests
    // ========================================================================

    fn make_fastq_records(records: &[(&str, &str)]) -> Vec<u8> {
        let mut data = Vec::new();
        for (name, seq) in records {
            data.extend_from_slice(format!("@{name}\n").as_bytes());
            data.extend_from_slice(seq.as_bytes());
            data.push(b'\n');
            data.extend_from_slice(b"+\n");
            // Quality = same length as sequence, all 'I'
            data.extend(std::iter::repeat_n(b'I', seq.len()));
            data.push(b'\n');
        }
        data
    }

    #[test]
    fn test_read_n_fastq_records_basic() {
        let data = make_fastq_records(&[("r1", "ACGT"), ("r2", "TGCA"), ("r3", "AAAA")]);
        let mut cursor = Cursor::new(data);

        let (buf, offsets, at_eof) = read_n_fastq_records(&mut cursor, 2).unwrap();

        assert_eq!(offsets.len(), 3); // 2 records + initial 0
        assert!(!at_eof);
        // First record starts at 0
        assert_eq!(offsets[0], 0);
        // Each record: "@rN\nACGT\n+\nIIII\n" = 4+4+2+4+1 = 15 bytes each
        assert!(offsets[1] > 0);
        assert!(offsets[2] > offsets[1]);
        // Verify first record starts with @
        assert_eq!(buf[0], b'@');
    }

    #[test]
    fn test_read_n_fastq_records_at_eof() {
        let data = make_fastq_records(&[("r1", "ACGT")]);
        let mut cursor = Cursor::new(data);

        let (_, offsets, at_eof) = read_n_fastq_records(&mut cursor, 5).unwrap();

        // Only 1 record available, requested 5
        assert_eq!(offsets.len(), 2); // 1 record
        assert!(at_eof);
    }

    #[test]
    fn test_read_n_fastq_records_empty_input() {
        let mut cursor = Cursor::new(Vec::<u8>::new());

        let (_, offsets, at_eof) = read_n_fastq_records(&mut cursor, 5).unwrap();

        assert_eq!(offsets.len(), 1); // Just the initial 0
        assert!(at_eof);
    }

    #[test]
    fn test_multi_record_count_reader_variable_length() {
        // R1 has short reads, R2 has long reads (the regression scenario)
        let r1_data = make_fastq_records(&[("r1_a", "AC"), ("r1_b", "TG"), ("r1_c", "AA")]);
        let r2_data = make_fastq_records(&[
            ("r2_a", "ACGTACGTACGT"),
            ("r2_b", "TGCATGCATGCA"),
            ("r2_c", "AAAACCCCGGGG"),
        ]);

        let r1 = Cursor::new(r1_data);
        let r2 = Cursor::new(r2_data);

        let mut reader = MultiRecordCountReader::new(vec![r1, r2]);
        reader.set_records_per_batch(2);

        let batch = reader.read_next_batch().unwrap().unwrap();

        assert_eq!(batch.streams.len(), 2);
        // Both streams should have exactly 2 records
        let r1_records = batch.streams[0].offsets.len() - 1;
        let r2_records = batch.streams[1].offsets.len() - 1;
        assert_eq!(r1_records, 2, "R1 should have 2 records");
        assert_eq!(r2_records, 2, "R2 should have 2 records");
    }

    #[test]
    fn test_multi_record_count_reader_eof_alignment() {
        // R1 has 3 records, R2 has 2 records
        let r1_data = make_fastq_records(&[("r1_a", "ACGT"), ("r1_b", "TGCA"), ("r1_c", "AAAA")]);
        let r2_data = make_fastq_records(&[("r2_a", "ACGT"), ("r2_b", "TGCA")]);

        let r1 = Cursor::new(r1_data);
        let r2 = Cursor::new(r2_data);

        let mut reader = MultiRecordCountReader::new(vec![r1, r2]);
        reader.set_records_per_batch(5); // Request more than available

        let batch = reader.read_next_batch().unwrap().unwrap();

        assert_eq!(batch.streams.len(), 2);
        // At EOF, should align to minimum (2 records)
        let r1_records = batch.streams[0].offsets.len() - 1;
        let r2_records = batch.streams[1].offsets.len() - 1;
        assert_eq!(r1_records, r2_records, "Record counts should be aligned at EOF");
        assert_eq!(r2_records, 2, "Should have min(3, 2) = 2 records");
    }

    #[test]
    fn test_multi_record_count_reader_single_stream() {
        let data = make_fastq_records(&[("r1", "ACGT"), ("r2", "TGCA")]);
        let cursor = Cursor::new(data);

        let mut reader = MultiRecordCountReader::new(vec![cursor]);
        reader.set_records_per_batch(1);

        let batch1 = reader.read_next_batch().unwrap().unwrap();
        assert_eq!(batch1.streams.len(), 1);
        assert_eq!(batch1.streams[0].offsets.len() - 1, 1);
        assert_eq!(batch1.serial, 0);

        let batch2 = reader.read_next_batch().unwrap().unwrap();
        assert_eq!(batch2.streams[0].offsets.len() - 1, 1);
        assert_eq!(batch2.serial, 1);

        // Third read should be None (EOF)
        assert!(reader.read_next_batch().unwrap().is_none());
        assert!(reader.is_done());
    }

    #[test]
    fn test_multi_record_count_reader_empty_input() {
        let cursor = Cursor::new(Vec::<u8>::new());
        let mut reader = MultiRecordCountReader::new(vec![cursor]);

        assert!(reader.read_next_batch().unwrap().is_none());
        assert!(reader.is_done());
    }

    #[test]
    fn test_multi_record_count_reader_produces_valid_boundaries() {
        let data = make_fastq_records(&[("read1", "ACGTACGT"), ("read2", "TGCATGCA")]);
        let cursor = Cursor::new(data);

        let mut reader = MultiRecordCountReader::new(vec![cursor]);
        reader.set_records_per_batch(2);

        let batch = reader.read_next_batch().unwrap().unwrap();
        let stream = &batch.streams[0];

        // Verify offsets point to valid record boundaries
        for i in 0..stream.offsets.len() - 1 {
            let start = stream.offsets[i];
            let end = if i + 1 < stream.offsets.len() {
                stream.offsets[i + 1]
            } else {
                stream.data.len()
            };
            // Each record should start with '@'
            assert_eq!(stream.data[start], b'@', "Record {i} should start with '@'");
            // Each record should end with '\n'
            assert_eq!(stream.data[end - 1], b'\n', "Record {i} should end with newline");
        }
    }

    // ========================================================================
    // ActiveSteps Configuration Tests
    // ========================================================================

    #[rstest]
    // Gzip+sync: active_steps() enforces synchronized → use_parallel_parse
    #[case::gzip_synchronized(false, true, false,
        &[Read, Decode, Process, Serialize, Compress, Write],
        &[Decompress, FindBoundaries, Group])]
    // BGZF+sync: Decompress+FindBoundaries+Decode active, Group skipped
    #[case::bgzf_synchronized(true, true, false,
        &[Read, Decompress, FindBoundaries, Decode, Process, Serialize, Compress, Write],
        &[Group])]
    // BGZF+non-sync+ppp=false: no FindBoundaries/Decode
    #[case::bgzf_non_synchronized(true, false, false,
        &[Read, Decompress, Group, Process, Serialize, Compress, Write],
        &[FindBoundaries, Decode])]
    // Gzip+non-sync+ppp=false: no FindBoundaries/Decode
    #[case::gzip_non_synchronized(false, false, false,
        &[Read, Decompress, Group, Process, Serialize, Compress, Write],
        &[FindBoundaries, Decode])]
    // BGZF+non-sync+ppp=true: all 9 steps
    #[case::all_steps_with_parallel_parse(true, false, true,
        &[Read, Decompress, FindBoundaries, Decode, Group, Process, Serialize, Compress, Write],
        &[])]
    fn test_active_steps(
        #[case] inputs_are_bgzf: bool,
        #[case] synchronized: bool,
        #[case] use_parallel_parse: bool,
        #[case] expected_active: &[PipelineStep],
        #[case] expected_inactive: &[PipelineStep],
    ) {
        let mut config = FastqPipelineConfig::new(4, inputs_are_bgzf, 1);
        if synchronized {
            config = config.with_synchronized(true);
        }
        config.use_parallel_parse = use_parallel_parse;
        let steps = config.active_steps();

        assert_eq!(steps.steps(), expected_active);
        for &step in expected_inactive {
            assert!(!steps.is_active(step));
        }
    }
}
