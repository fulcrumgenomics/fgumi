//! Unified FASTQ pipeline with per-stream parallel reading.
//!
//! This module implements the FASTQ-specific pipeline for the `extract` command,
//! which reads from multiple synchronized FASTQ streams (R1, R2, I1, I2) and
//! creates templates by zipping records at matching positions.
//!
//! The pipeline always uses 8 active steps for both gzip and BGZF inputs:
//!
//! Read → Decompress → `FindBoundaries` (Pair) → Decode → Process → Serialize → Compress → Write
//!
//! - **Read**: Per-stream parallel reading. Each stream has its own mutex,
//!   allowing multiple threads to read different streams concurrently.
//! - **Decompress**: For BGZF, decompresses raw blocks. For gzip/plain, passthrough.
//! - **`FindBoundaries` (Pair)**: Assembles per-stream chunks by `batch_num` into
//!   multi-stream batches. For BGZF, finds record boundaries. For gzip, uses
//!   pre-computed offsets directly.
//! - **Decode**: Parses FASTQ records and creates templates by zipping records
//!   at matching positions, bypassing the Group step entirely.

use std::collections::BTreeMap;

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
use crate::fastq_parse::FastqRecord;
use crate::grouper::FastqTemplate;
use crate::progress::ProgressTracker;
use crate::reorder_buffer::ReorderBuffer;
use libdeflater::Decompressor;

use super::base::{
    ActiveSteps, CompressedBlockBatch, HasCompressor, HasHeldBoundaries, HasHeldCompressed,
    HasHeldProcessed, HasHeldSerialized, HasRecycledBuffers, HasWorkerCore, MemoryEstimate,
    MonitorableState, OutputPipelineQueues, OutputPipelineState, PipelineLifecycle, PipelineStats,
    PipelineStep, PipelineValidationError, ProcessPipelineState, SerializePipelineState,
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
// Per-Stream Parallel Reading Types
// ============================================================================

/// Reader for a single FASTQ input stream. Each stream has its own lock
/// to allow multiple threads to read different streams concurrently.
pub enum StreamReader<R: BufRead + Send> {
    /// Raw BGZF file — Read step produces raw blocks, Decompress step decompresses.
    Bgzf(BufReader<File>),
    /// Pre-decompressed reader (gzip/plain) — Read step produces record-aligned data.
    Decompressed(R),
}

/// A chunk of data from a single FASTQ stream, produced by the per-stream Read step.
///
/// For BGZF inputs: `data` contains raw concatenated BGZF blocks, `offsets` is `None`.
/// For gzip/plain inputs: `data` contains record-aligned data, `offsets` is `Some`.
#[derive(Debug, Clone)]
pub struct PerStreamChunk {
    /// Which stream this chunk came from (0=R1, 1=R2, etc.)
    pub stream_idx: usize,
    /// Per-stream batch number (monotonically increasing per stream).
    pub batch_num: u64,
    /// Chunk data (raw BGZF or decompressed record-aligned).
    pub data: Vec<u8>,
    /// Record boundary offsets. Present for gzip (record-aligned), absent for BGZF.
    pub offsets: Option<Vec<usize>>,
}

impl PerStreamChunk {
    /// Estimate heap memory usage.
    #[must_use]
    pub fn estimate_heap_size(&self) -> usize {
        self.data.capacity()
            + self.offsets.as_ref().map_or(0, |o| o.capacity() * std::mem::size_of::<usize>())
    }
}

/// State for the Pair step that assembles per-stream chunks into multi-stream batches.
///
/// Accumulates per-stream decompressed chunks by `batch_num`. When all streams
/// have delivered their chunk for a given `batch_num`, the batch is emitted.
///
/// Completion is detected by the Pair step function using count-based tracking:
/// when `read_done && chunks_paired == batches_read`, all chunks have arrived
/// and remaining incomplete batches can be flushed. This follows the same pattern
/// as the BAM pipeline's `batches_boundary_processed == total_read`.
pub(crate) struct PairState {
    /// Per `batch_num`: Vec of `Option<PerStreamChunk>`, one slot per stream.
    pending: BTreeMap<u64, Vec<Option<PerStreamChunk>>>,
    /// Next `batch_num` to emit (for ordered output).
    next_emit: u64,
    /// Number of input streams.
    num_streams: usize,
}

impl PairState {
    fn new(num_streams: usize) -> Self {
        Self { pending: BTreeMap::new(), next_emit: 0, num_streams }
    }

    /// Insert a data chunk into the pending map.
    fn insert(&mut self, chunk: PerStreamChunk) {
        let stream_idx = chunk.stream_idx;
        let batch_num = chunk.batch_num;
        let slots = self.pending.entry(batch_num).or_insert_with(|| vec![None; self.num_streams]);
        slots[stream_idx] = Some(chunk);
    }

    /// Try to pop a complete set of chunks for the next `batch_num`.
    ///
    /// When `all_arrived` is false (normal operation), ALL streams must have
    /// delivered their chunk. When `all_arrived` is true (all Read chunks have
    /// been consumed), incomplete batches are emitted with whatever data is present.
    fn try_pop_complete(&mut self, all_arrived: bool) -> Option<Vec<PerStreamChunk>> {
        let slots = self.pending.get(&self.next_emit)?;
        let complete = if all_arrived {
            slots.iter().any(Option::is_some)
        } else {
            slots.iter().all(Option::is_some)
        };
        if !complete {
            return None;
        }

        let slots = self
            .pending
            .remove(&self.next_emit)
            .expect("next_emit key must exist in pending map after get() succeeded");
        self.next_emit += 1;
        Some(slots.into_iter().flatten().collect())
    }

    fn is_empty(&self) -> bool {
        self.pending.is_empty()
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
        self.chunks.iter().map(|c| c.data.capacity()).sum::<usize>()
            + self.chunks.capacity() * std::mem::size_of::<FastqDecompressedChunk>()
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
            .sum::<usize>()
            + self.streams.capacity() * std::mem::size_of::<FastqStreamBoundaries>()
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
                stream.records.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
                    + stream.records.capacity() * std::mem::size_of::<FastqRecord>()
            })
            .sum::<usize>()
            + self.streams.capacity() * std::mem::size_of::<FastqParsedStream>()
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

/// Find FASTQ record boundaries in data using SIMD-accelerated newline detection.
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

    let offsets = fgumi_simd_fastq::find_record_offsets(data);
    let last_offset = offsets.last().copied().unwrap_or(0);
    let complete_data = data[..last_offset].to_vec();

    (complete_data, offsets, last_offset)
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
    FastqRecord::from_slice(data)
}

// ============================================================================
// BlockParseFast / BlockMerge — Parallel BGZF Boundary Detection
/// Max blocks processed per lock acquisition in `BlockMerge` / `FindBoundaries` steps.
const MAX_BATCHES_PER_LOCK: usize = 8;
// ============================================================================

/// Result of parallel per-chunk parsing in the `BlockParseFast` step.
///
/// Each decompressed BGZF chunk is split into:
/// - `prefix_bytes`: incomplete record fragment at the start (belongs to previous block)
/// - `records`: fully-parsed complete records from the middle
/// - `suffix_bytes`: incomplete record fragment at the end (belongs to next block)
#[derive(Debug)]
pub struct BlockParsed {
    /// Monotonically increasing block index (per stream), used for ordering.
    pub block_idx: u64,
    /// Which stream this block came from (0=R1, 1=R2, etc.).
    pub stream_idx: usize,
    /// Fully parsed complete records from the middle of the block.
    pub records: Vec<FastqRecord>,
    /// Incomplete prefix bytes (start of the block, part of the previous record).
    pub prefix_bytes: Vec<u8>,
    /// Incomplete suffix bytes (end of the block, part of the next record).
    pub suffix_bytes: Vec<u8>,
    // is_eof reserved for future use (EOF signalling between steps).
}

impl MemoryEstimate for BlockParsed {
    fn estimate_heap_size(&self) -> usize {
        self.records.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
            + self.records.capacity() * std::mem::size_of::<FastqRecord>()
            + self.prefix_bytes.capacity()
            + self.suffix_bytes.capacity()
    }
}

/// Mutable state owned by the serial `BlockMerge` step.
///
/// The merge step pairs corresponding R1 and R2 blocks (by `block_idx`), stitches
/// the suffix bytes of block N with the prefix bytes of block N+1 to recover
/// cross-block records, and zips R1/R2 records into `FastqTemplate`s.
pub(crate) struct BlockMergeState {
    /// Received R1 blocks not yet processed, keyed by `block_idx`.
    r1_pending: BTreeMap<u64, BlockParsed>,
    /// Received R2 blocks not yet processed, keyed by `block_idx`.
    r2_pending: BTreeMap<u64, BlockParsed>,
    /// Next R1 block index expected by the merge step.
    r1_next: u64,
    /// Next R2 block index expected by the merge step.
    r2_next: u64,
    /// Trailing suffix bytes from the last R1 block (cross-block record fragment).
    r1_suffix_bytes: Vec<u8>,
    /// Trailing suffix bytes from the last R2 block (cross-block record fragment).
    r2_suffix_bytes: Vec<u8>,
    /// Excess R1 records left over from the previous round (R1 had more than R2).
    r1_surplus: Vec<FastqRecord>,
    /// Excess R2 records left over from the previous round (R2 had more than R1).
    r2_surplus: Vec<FastqRecord>,
    /// Monotonically increasing output serial for the template queue.
    serial_out: u64,
    /// Estimated heap bytes held in `r1_pending` + `r2_pending` maps.
    ///
    /// Used for local backpressure: when this exceeds [`PENDING_BACKPRESSURE_BYTES`]
    /// and the merge step *can* process in-order blocks, we skip draining `q2` so
    /// the queue fills up and naturally throttles `BlockParseFast` workers.
    pending_heap_bytes: u64,
}

impl BlockMergeState {
    fn new() -> Self {
        Self {
            r1_pending: BTreeMap::new(),
            r2_pending: BTreeMap::new(),
            r1_next: 0,
            r2_next: 0,
            r1_suffix_bytes: Vec::new(),
            r2_suffix_bytes: Vec::new(),
            r1_surplus: Vec::new(),
            r2_surplus: Vec::new(),
            serial_out: 0,
            pending_heap_bytes: 0,
        }
    }

    fn is_empty(&self) -> bool {
        let empty = self.r1_pending.is_empty()
            && self.r2_pending.is_empty()
            && self.r1_surplus.is_empty()
            && self.r2_surplus.is_empty()
            && self.r1_suffix_bytes.is_empty()
            && self.r2_suffix_bytes.is_empty();
        debug_assert!(
            !empty || self.pending_heap_bytes == 0,
            "pending_heap_bytes={} but maps are empty",
            self.pending_heap_bytes
        );
        empty
    }
}

/// Detect where the first complete FASTQ record begins in `data`.
///
/// Scans forward collecting newline positions and tests sliding windows of 4
/// consecutive lines against the FASTQ invariants (`@` prefix, `+` separator,
/// equal seq/qual lengths). Returns the byte offset at which the first
/// complete record starts.
///
/// Returns `0` if the first complete record starts at the beginning (including
/// after validating all four FASTQ lines). Returns `data.len()` if no complete
/// record boundary can be detected.
///
/// Note: does NOT use a `data[0] == b'@'` fast-path because FASTQ quality
/// strings can legally contain `@` bytes (Phred+33 quality ≥ 31).
fn detect_prefix_end(data: &[u8]) -> usize {
    if data.is_empty() {
        return 0;
    }

    // Collect the first 8 newline positions.
    let mut newlines = [0usize; 8];
    let mut count = 0;
    for (i, &b) in data.iter().enumerate() {
        if b == b'\n' {
            newlines[count] = i;
            count += 1;
            if count == 8 {
                break;
            }
        }
    }

    // Need at least 4 newlines to form a complete FASTQ record.
    if count < 4 {
        return data.len(); // No complete record found — entire buffer is prefix.
    }

    // Try each sliding window of 4 consecutive newlines.
    for start in 0..count.saturating_sub(3) {
        let line0_start = if start == 0 { 0 } else { newlines[start - 1] + 1 };
        let line0_end = newlines[start]; // end of "@name" line (newline pos)
        let line1_end = newlines[start + 1]; // end of sequence line
        let line2_end = newlines[start + 2]; // end of "+" line
        let line3_end = newlines[start + 3]; // end of quality line

        // Line 0 must start with '@'
        if data[line0_start] != b'@' {
            continue;
        }

        // Line 2 must start with '+'
        let line2_start = line1_end + 1;
        if line2_start >= data.len() || data[line2_start] != b'+' {
            continue;
        }

        // Sequence and quality lengths must match.
        let seq_len = line1_end - (line0_end + 1);
        let qual_len = line3_end - (line2_end + 1);
        if seq_len != qual_len {
            continue;
        }

        // Valid record found at line0_start.
        return line0_start;
    }

    // Could not detect boundary.
    data.len()
}

/// Detect where the last complete FASTQ record ends in `data`.
///
/// Scans backward collecting newline positions and tests the last 4 newlines
/// against FASTQ invariants. Returns the byte offset of the end of the last
/// complete record (i.e., `data[..suffix_start]` contains all complete records
/// and `data[suffix_start..]` is an incomplete trailing fragment).
/// Returns `data.len()` if all data forms complete records.
fn detect_suffix_start(data: &[u8]) -> usize {
    if data.is_empty() {
        return 0;
    }

    // Collect the last 8 newline positions scanning backward.
    let mut newlines = [0usize; 8];
    let mut count = 0;
    let mut i = data.len();
    while i > 0 {
        i -= 1;
        if data[i] == b'\n' {
            // Store in reverse order, then reverse at end.
            newlines[count] = i;
            count += 1;
            if count == 8 {
                break;
            }
        }
    }

    if count < 4 {
        // Not enough newlines to form a complete record.
        return 0;
    }

    // Reverse so newlines are in ascending order.
    newlines[..count].reverse();

    // Try windows from the last 4 newlines backward.
    // We want the latest valid window (rightmost complete record).
    let window_start = count.saturating_sub(4);
    for start in (0..=window_start).rev() {
        if start + 3 >= count {
            continue;
        }
        let line0_start = if start == 0 { 0 } else { newlines[start - 1] + 1 };
        let line0_end = newlines[start];
        let line1_end = newlines[start + 1];
        let line2_end = newlines[start + 2];
        let line3_end = newlines[start + 3];

        if data[line0_start] != b'@' {
            continue;
        }
        let line2_start = line1_end + 1;
        if line2_start >= data.len() || data[line2_start] != b'+' {
            continue;
        }
        let seq_len = line1_end - (line0_end + 1);
        let qual_len = line3_end - (line2_end + 1);
        if seq_len != qual_len {
            continue;
        }

        // Valid: the record ends at line3_end + 1.
        return line3_end + 1;
    }

    0
}

/// Stitch `suffix_bytes` from the previous block with `prefix_bytes` from the
/// current block to recover the cross-block FASTQ record (if any).
///
/// Returns `None` if both slices are empty or the combined data is not a valid record.
/// Any error is propagated as an `io::Error`.
fn stitch_cross_block_record(
    suffix_bytes: &[u8],
    prefix_bytes: &[u8],
) -> io::Result<Option<FastqRecord>> {
    if suffix_bytes.is_empty() && prefix_bytes.is_empty() {
        return Ok(None);
    }
    let mut combined = Vec::with_capacity(suffix_bytes.len() + prefix_bytes.len());
    combined.extend_from_slice(suffix_bytes);
    combined.extend_from_slice(prefix_bytes);
    let record = FastqRecord::from_slice(&combined)?;
    Ok(Some(record))
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
    /// Shared statistics instance for external memory monitoring.
    /// When provided, the pipeline will use this instead of creating its own.
    pub shared_stats: Option<Arc<PipelineStats>>,
    /// Number of FASTQ records per stream per batch for `RecordCount` readers.
    /// Scales with thread count to reduce queue contention at higher parallelism.
    /// Default: 200 (auto-scaled in `new()` based on thread count).
    pub records_per_batch: usize,
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
        // Scale FASTQ batch size: more threads = larger batches to reduce queue ops
        // Cap at 800 (beyond which diminishing returns and memory increases)
        let records_per_batch = (200 * threads.min(4)).min(800);
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
            shared_stats: None, // No shared stats by default
            records_per_batch,
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

    /// Compute the active pipeline steps.
    ///
    /// Always returns 8 steps: Read, Decompress, `FindBoundaries`, Decode,
    /// Process, Serialize, Compress, Write.
    ///
    /// - Decompress: for BGZF decompresses blocks; for gzip/plain passes through.
    /// - `FindBoundaries`: Pair step that assembles per-stream chunks into
    ///   multi-stream batches and finds record boundaries.
    /// - Decode: parses FASTQ records and creates templates directly.
    /// - Group is never active.
    #[must_use]
    pub fn active_steps(&self) -> ActiveSteps {
        ActiveSteps::new(&[
            PipelineStep::Read,
            PipelineStep::Decompress,
            PipelineStep::FindBoundaries,
            PipelineStep::Decode,
            PipelineStep::Process,
            PipelineStep::Serialize,
            PipelineStep::Compress,
            PipelineStep::Write,
        ])
    }
}

// ============================================================================
// BGZF Block Reading Constants
// ============================================================================

/// Default number of BGZF blocks to read per stream per batch.
const DEFAULT_BLOCKS_PER_BATCH: usize = 4;

/// Memory backpressure threshold for `q2_block_parsed` (128 MB).
///
/// When the total heap bytes in the parsed-block queue exceeds this threshold,
/// `BlockParseFast` workers stop pulling new decompressed chunks, giving the
/// serial `BlockMerge` step time to drain. This prevents unbounded RSS growth
/// when workers parse faster than the single-threaded merge can consume.
const Q2_BLOCK_PARSED_BACKPRESSURE_BYTES: u64 = 128 * 1024 * 1024;

/// Memory backpressure threshold for the `BlockMerge` pending maps (256 MB).
///
/// When the total heap bytes sitting in `r1_pending` + `r2_pending` exceeds
/// this threshold **and** the merge step can process the next in-order blocks,
/// we skip draining `q2_block_parsed` so the queue fills up, creating natural
/// backpressure on `BlockParseFast` workers via `ArrayQueue::is_full()`.
///
/// When the next in-order blocks are *not* available, we always drain regardless
/// of this threshold — ensuring the needed blocks can reach the pending maps
/// and preventing deadlock.
const PENDING_BACKPRESSURE_BYTES: u64 = 256 * 1024 * 1024;

// ============================================================================
// Helper functions (used by per-stream pipeline steps)
// ============================================================================

/// Read exactly `n` complete FASTQ records from a buffered reader.
///
/// Uses `BufRead::read_until` for line-by-line reading. Each FASTQ record
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
// Phase 6 & 7: FASTQ Pipeline State and Entry Point
// ============================================================================

/// Align record counts across multiple streams, truncating excess to the minimum.
///
/// When streams have different record counts (e.g., at EOF), excess records
/// are discarded since they have no mate and cannot form valid templates.
fn align_stream_records(
    mut streams: Vec<FastqStreamBoundaries>,
    serial: u64,
) -> FastqBoundaryBatch {
    if streams.len() > 1 {
        let min_records =
            streams.iter().map(|s| s.offsets.len().saturating_sub(1)).min().unwrap_or(0);
        for stream in &mut streams {
            let record_count = stream.offsets.len().saturating_sub(1);
            if record_count > min_records && min_records > 0 {
                let excess_start = stream.offsets[min_records];
                stream.data.truncate(excess_start);
                stream.offsets.truncate(min_records + 1);
            }
        }
    }
    FastqBoundaryBatch { streams, serial }
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
    /// Next stream index to try reading from (round-robin across streams).
    pub next_stream: usize,
    /// Held per-stream chunk from Read step (couldn't push to `q0_chunks`).
    pub held_chunk: Option<(u64, PerStreamChunk)>,
    /// Held decompressed chunk from Decompress step (couldn't push to `q1_decompressed`).
    pub held_decompressed_chunk: Option<(u64, PerStreamChunk)>,
    /// Held boundary batch from `FindBoundaries` step (couldn't push to `q2_5_boundaries`).
    pub held_boundaries: Option<(u64, FastqBoundaryBatch)>,
    /// Held `BlockParsed` item from `BlockParseFast` step (couldn't push to `q2_block_parsed`).
    pub held_block_parsed: Option<BlockParsed>,
    /// Held parsed templates from Parse step (couldn't push to `output.groups`).
    /// Includes template count for metrics tracking on successful push.
    pub held_parsed: Option<(u64, Vec<FastqTemplate>, usize)>,
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
            next_stream: thread_id, // stagger starting stream across workers
            held_chunk: None,
            held_decompressed_chunk: None,
            held_boundaries: None,
            held_block_parsed: None,
            held_parsed: None,
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
        self.held_chunk.is_some()
            || self.held_decompressed_chunk.is_some()
            || self.held_boundaries.is_some()
            || self.held_block_parsed.is_some()
            || self.held_parsed.is_some()
            || self.held_processed.is_some()
            || self.held_serialized.is_some()
            || self.held_compressed.is_some()
    }

    /// Clear all held items (for cleanup/error handling).
    pub fn clear_held_items(&mut self) {
        self.held_chunk = None;
        self.held_decompressed_chunk = None;
        self.held_boundaries = None;
        self.held_block_parsed = None;
        self.held_parsed = None;
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

impl<P: Send> HasRecycledBuffers for FastqWorkerState<P> {
    fn take_or_alloc_buffer(&mut self, capacity: usize) -> Vec<u8> {
        self.core.take_or_alloc_buffer(capacity)
    }

    fn recycle_buffer(&mut self, buf: Vec<u8>) {
        self.core.recycle_buffer(buf);
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

    // ========== Step 1: Per-Stream Read ==========
    /// Per-stream readers, each independently lockable for parallel reading.
    pub readers: Vec<Mutex<Option<StreamReader<R>>>>,
    /// Per-stream monotonic batch counter.
    pub batch_counters: Vec<AtomicU64>,
    /// Per-stream EOF flags.
    pub stream_eof: Vec<AtomicBool>,
    /// Number of input streams.
    pub num_streams: usize,
    /// Number of BGZF blocks to read per stream per batch.
    pub blocks_per_batch: usize,
    /// Number of FASTQ records to read per stream per batch (for gzip/plain).
    pub records_per_batch: usize,
    /// Flag indicating all streams have reached EOF.
    pub read_done: AtomicBool,
    /// Count of per-stream chunks successfully read.
    pub batches_read: AtomicU64,

    // ========== Q0: Read → Decompress ==========
    /// Per-stream chunks waiting to be decompressed.
    pub q0_chunks: ArrayQueue<(u64, PerStreamChunk)>,

    // ========== Q1: Decompress → BlockParseFast (parallel) ==========
    /// Decompressed per-stream chunks waiting for parallel block parsing.
    pub q1_decompressed: ArrayQueue<(u64, PerStreamChunk)>,

    // ========== BlockParseFast → BlockMerge ==========
    /// Count of per-stream chunks consumed by `BlockParseFast` from q1.
    pub chunks_block_parsed: AtomicU64,
    /// Parsed blocks waiting for the serial `BlockMerge` step.
    /// Each element is a `BlockParsed` from one per-stream chunk.
    pub q2_block_parsed: ArrayQueue<BlockParsed>,
    /// Current heap bytes in `q2_block_parsed` (for memory backpressure).
    pub q2_block_parsed_heap_bytes: AtomicU64,
    /// Serial `BlockMerge` step state (locked via `try_lock` for serial execution).
    pub(crate) block_merge_state: Mutex<BlockMergeState>,
    /// Flag indicating all blocks have been merged and templates emitted.
    pub block_merge_done: AtomicBool,
    /// Count of `BlockParsed` items consumed by `BlockMerge`.
    pub blocks_merged: AtomicU64,

    // ========== Legacy: Pair (FindBoundaries) → Decode (gzip path) ==========
    /// Count of per-stream chunks consumed by the Pair step from q1 (gzip path).
    pub chunks_paired: AtomicU64,
    /// Pair state: accumulates per-stream chunks by `batch_num` (gzip path).
    pub(crate) pair_state: Mutex<PairState>,
    /// State for finding FASTQ record boundaries (gzip path).
    pub boundary_state: FastqBoundaryState,
    /// Flag indicating pair assembly / boundary finding is complete (gzip path).
    pub boundaries_done: AtomicBool,
    /// Count of multi-stream batches assembled by Pair step (gzip path).
    pub batches_boundaries_found: AtomicU64,
    /// Batches with boundaries found, waiting to be parsed (gzip path).
    pub q2_5_boundaries: ArrayQueue<(u64, FastqBoundaryBatch)>,
    /// Current heap bytes in Q2.5 boundaries queue (gzip path).
    pub q2_5_boundaries_heap_bytes: AtomicU64,

    // ========== Decode → Process ==========
    /// Flag indicating parsing is complete.
    pub parse_done: AtomicBool,
    /// Count of batches that have been parsed.
    pub batches_parsed: AtomicU64,
    /// Count of batches that have been grouped (same as parsed).
    pub batches_grouped: AtomicU64,
    /// Flag indicating grouping is complete.
    pub group_done: AtomicBool,
    /// Total number of individual templates pushed to Q3 (for debugging).
    pub total_templates_pushed: AtomicU64,
    /// Total number of individual records serialized (for completion check).
    pub total_records_serialized: AtomicU64,

    // ========== Output-Half State (Process → Serialize → Compress → Write) ==========
    /// Shared output pipeline queues and state.
    pub output: OutputPipelineQueues<FastqTemplate, P>,

    // ========== Deadlock Detection ==========
    /// State for deadlock detection and recovery.
    pub deadlock_state: DeadlockState,
}

impl<R: BufRead + Send, P: Send + MemoryEstimate> FastqPipelineState<R, P> {
    /// Create a new pipeline state with per-stream readers.
    #[must_use]
    pub fn new(
        config: FastqPipelineConfig,
        readers: Vec<StreamReader<R>>,
        output: Box<dyn Write + Send>,
    ) -> Self {
        let cap = config.queue_capacity;
        let num_streams = readers.len();
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

        let per_stream_readers: Vec<Mutex<Option<StreamReader<R>>>> =
            readers.into_iter().map(|r| Mutex::new(Some(r))).collect();
        let batch_counters: Vec<AtomicU64> = (0..num_streams).map(|_| AtomicU64::new(0)).collect();
        let stream_eof: Vec<AtomicBool> =
            (0..num_streams).map(|_| AtomicBool::new(false)).collect();

        Self {
            readers: per_stream_readers,
            batch_counters,
            stream_eof,
            num_streams,
            blocks_per_batch: DEFAULT_BLOCKS_PER_BATCH,
            records_per_batch: config.records_per_batch,
            read_done: AtomicBool::new(false),
            batches_read: AtomicU64::new(0),
            q0_chunks: ArrayQueue::new(cap),
            q1_decompressed: ArrayQueue::new(cap),
            chunks_block_parsed: AtomicU64::new(0),
            q2_block_parsed: ArrayQueue::new(cap),
            q2_block_parsed_heap_bytes: AtomicU64::new(0),
            block_merge_state: Mutex::new(BlockMergeState::new()),
            block_merge_done: AtomicBool::new(false),
            blocks_merged: AtomicU64::new(0),
            chunks_paired: AtomicU64::new(0),
            pair_state: Mutex::new(PairState::new(num_streams)),
            boundary_state: FastqBoundaryState::new(num_streams),
            boundaries_done: AtomicBool::new(false),
            batches_boundaries_found: AtomicU64::new(0),
            q2_5_boundaries: ArrayQueue::new(cap),
            q2_5_boundaries_heap_bytes: AtomicU64::new(0),
            parse_done: AtomicBool::new(false),
            batches_parsed: AtomicU64::new(0),
            batches_grouped: AtomicU64::new(0),
            group_done: AtomicBool::new(false),
            total_templates_pushed: AtomicU64::new(0),
            total_records_serialized: AtomicU64::new(0),
            output: OutputPipelineQueues::new(
                cap,
                output,
                stats,
                "Processed records",
                memory_limit,
            ),
            deadlock_state,
            config,
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

    /// Check if `q2_block_parsed` memory exceeds the backpressure threshold.
    ///
    /// Used by `BlockParseFast` to avoid piling up parsed data faster than
    /// the serial `BlockMerge` can drain it.
    #[must_use]
    pub fn is_q2_block_parsed_memory_high(&self) -> bool {
        self.q2_block_parsed_heap_bytes.load(Ordering::Acquire)
            >= Q2_BLOCK_PARSED_BACKPRESSURE_BYTES
    }

    /// Check if Q4 processed queue memory is high (for backpressure in Process step).
    ///
    /// Uses the same threshold as BAM's Q5 processed queue to ensure consistent
    /// memory backpressure behavior across both pipelines.
    #[must_use]
    pub fn is_q4_memory_high(&self) -> bool {
        self.output.is_processed_memory_high()
    }

    /// Check if the pipeline is in drain mode (input exhausted, completing remaining work).
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
            log::trace!("is_complete: read_done={read_done}, group_done={group_done}");
            return false;
        }

        if self.config.inputs_are_bgzf {
            // BGZF path: check block_merge_done and parse_done
            let block_merge_done = self.block_merge_done.load(Ordering::Acquire);
            let parse_done = self.parse_done.load(Ordering::Acquire);
            if !block_merge_done || !parse_done {
                log::debug!(
                    "is_complete: BGZF flags not done: block_merge_done={block_merge_done}, parse_done={parse_done}"
                );
                return false;
            }
            // Check BGZF-specific queues
            if !self.q0_chunks.is_empty()
                || !self.q1_decompressed.is_empty()
                || !self.q2_block_parsed.is_empty()
            {
                return false;
            }
        } else {
            // Gzip/plain path: check boundaries_done and parse_done
            let boundaries_done = self.boundaries_done.load(Ordering::Acquire);
            let parse_done = self.parse_done.load(Ordering::Acquire);
            if !boundaries_done || !parse_done {
                log::debug!(
                    "is_complete: gzip flags not done: boundaries_done={boundaries_done}, parse_done={parse_done}"
                );
                return false;
            }
            // Check input-half ArrayQueues are empty (lock-free checks)
            if !self.q0_chunks.is_empty() || !self.q1_decompressed.is_empty() {
                return false;
            }
            // Check intermediate queues
            if !self.q2_5_boundaries.is_empty() {
                log::trace!(
                    "is_complete: q2_5_boundaries not empty: {}",
                    self.q2_5_boundaries.len()
                );
                return false;
            }
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
        if !self.q0_chunks.is_empty() {
            non_empty_queues.push(format!("q0_chunks ({})", self.q0_chunks.len()));
        }
        if !self.q1_decompressed.is_empty() {
            non_empty_queues.push(format!("q1_decompressed ({})", self.q1_decompressed.len()));
        }

        if self.config.inputs_are_bgzf {
            if !self.q2_block_parsed.is_empty() {
                non_empty_queues.push(format!("q2_block_parsed ({})", self.q2_block_parsed.len()));
            }
            // Check BlockMerge state is empty
            {
                let merge = self.block_merge_state.lock();
                if !merge.is_empty() {
                    non_empty_queues.push("block_merge_state (non-empty)".to_string());
                }
            }
        } else {
            if !self.q2_5_boundaries.is_empty() {
                non_empty_queues.push(format!("q2_5_boundaries ({})", self.q2_5_boundaries.len()));
            }
            // Check pair state is empty
            {
                let pair = self.pair_state.lock();
                if !pair.is_empty() {
                    non_empty_queues.push("pair_state (non-empty)".to_string());
                }
            }
            // Check boundary state has no leftover bytes
            for (idx, stream_state) in self.boundary_state.stream_states.iter().enumerate() {
                let leftover_len = stream_state.lock().leftover.len();
                if leftover_len > 0 {
                    non_empty_queues.push(format!("boundary_leftover[{idx}] ({leftover_len})"));
                }
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

        // Check write reorder buffer is empty
        {
            let write_reorder = self.output.write_reorder.lock();
            if !write_reorder.is_empty() {
                non_empty_queues.push(format!("write_reorder ({})", write_reorder.len()));
            }
        }

        // Check batch counter invariants (gzip path only — BGZF uses different counters)
        if !self.config.inputs_are_bgzf {
            let batches_grouped = self.batches_grouped.load(Ordering::Acquire);
            let batches_boundaries_found = self.batches_boundaries_found.load(Ordering::Acquire);
            let batches_parsed = self.batches_parsed.load(Ordering::Acquire);

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
        let parse_done = self.parse_done.load(Ordering::Relaxed);
        let batches_read = self.batches_read.load(Ordering::Relaxed);
        let (q2b_len, extra_state) = if self.config.inputs_are_bgzf {
            let block_merge_done = self.block_merge_done.load(Ordering::Relaxed);
            let chunks_bp = self.chunks_block_parsed.load(Ordering::Relaxed);
            let blocks_merged = self.blocks_merged.load(Ordering::Relaxed);
            let q2_heap_mb =
                self.q2_block_parsed_heap_bytes.load(Ordering::Relaxed) / (1024 * 1024);
            (
                self.q2_block_parsed.len(),
                Some(format!(
                    "block_merge_done={block_merge_done}, parse_done={parse_done}, batches: read={batches_read} block_parsed={chunks_bp} merged={blocks_merged}, q2_heap={q2_heap_mb}MB"
                )),
            )
        } else {
            let boundaries_done = self.boundaries_done.load(Ordering::Relaxed);
            let chunks_paired = self.chunks_paired.load(Ordering::Relaxed);
            let batches_found = self.batches_boundaries_found.load(Ordering::Relaxed);
            let batches_parsed = self.batches_parsed.load(Ordering::Relaxed);
            (
                self.q2_5_boundaries.len(),
                Some(format!(
                    "boundaries_done={boundaries_done}, parse_done={parse_done}, batches: read={batches_read} paired={chunks_paired} found={batches_found} parsed={batches_parsed}"
                )),
            )
        };
        QueueSnapshot {
            q1_len: self.q0_chunks.len(),
            q2_len: self.q1_decompressed.len(),
            q2b_len,
            q3_len: self.output.groups.len(),
            q4_len: self.output.processed.len(),
            q5_len: self.output.serialized.len(),
            q6_len: self.output.compressed.len(),
            q7_len: 0,         // Not used in FASTQ (q6_compressed is the write input)
            q2_reorder_mem: 0, // No reorder buffer in new per-stream pipeline
            q3_reorder_mem: 0,
            memory_limit: self.deadlock_state.get_memory_limit(),
            read_done: self.read_done.load(Ordering::Relaxed),
            group_done: self.group_done.load(Ordering::Relaxed),
            draining: self.output.draining.load(Ordering::Relaxed),
            extra_state,
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

    fn write_reorder_can_proceed(&self, serial: u64) -> bool {
        self.output.write_reorder_state.can_proceed(serial)
    }

    fn write_reorder_is_memory_high(&self) -> bool {
        self.output.write_reorder_state.is_memory_high()
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

    fn write_reorder_state(&self) -> &super::base::ReorderBufferState {
        &self.output.write_reorder_state
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

/// Try to read from any available per-stream reader (Step 1).
///
/// Multiple worker threads can read different streams concurrently via per-stream
/// `try_lock()`. For BGZF inputs, produces raw BGZF data (no offsets). For gzip/plain
/// inputs, produces record-aligned data with offsets.
#[allow(clippy::too_many_lines)]
fn fastq_try_step_read<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
    worker: &mut FastqWorkerState<P>,
) -> bool {
    // Priority 1: Try to advance held chunk
    if let Some((serial, held)) = worker.held_chunk.take() {
        match state.q0_chunks.push((serial, held)) {
            Ok(()) => {
                state.deadlock_state.record_q1_push();
            }
            Err((serial, held)) => {
                worker.held_chunk = Some((serial, held));
                return false;
            }
        }
    }

    // Priority 2: Check for completion/error
    if state.read_done.load(Ordering::Relaxed) || state.has_error() {
        return false;
    }

    // Priority 3: Check if output queue has space
    if state.q0_chunks.is_full() {
        return false;
    }

    // Priority 4: Try to acquire any stream's reader (round-robin to balance reads)
    let start = worker.next_stream % state.num_streams;
    for i in 0..state.num_streams {
        let stream_idx = (start + i) % state.num_streams;
        if state.stream_eof[stream_idx].load(Ordering::Relaxed) {
            continue;
        }
        let Some(mut guard) = state.readers[stream_idx].try_lock() else {
            continue; // Another thread has this stream
        };
        let Some(ref mut reader) = *guard else {
            continue;
        };

        // Advance round-robin so next call starts from a different stream
        worker.next_stream = stream_idx + 1;

        let batch_num = state.batch_counters[stream_idx].fetch_add(1, Ordering::Relaxed);

        match reader {
            StreamReader::Bgzf(r) => {
                match read_raw_blocks(r, state.blocks_per_batch) {
                    Ok(blocks) if blocks.is_empty() => {
                        // EOF — undo counter, set flags
                        state.batch_counters[stream_idx].fetch_sub(1, Ordering::Relaxed);
                        state.stream_eof[stream_idx].store(true, Ordering::Release);
                        if state.stream_eof.iter().all(|f| f.load(Ordering::Acquire)) {
                            state.read_done.store(true, Ordering::Release);
                        }
                    }
                    Ok(blocks) => {
                        // Concatenate raw block data
                        let total_size: usize = blocks.iter().map(|b| b.data.len()).sum();
                        let mut raw_data = Vec::with_capacity(total_size);
                        for block in blocks {
                            raw_data.extend_from_slice(&block.data);
                        }
                        let serial = state.batches_read.fetch_add(1, Ordering::Release);
                        let chunk =
                            PerStreamChunk { stream_idx, batch_num, data: raw_data, offsets: None };
                        match state.q0_chunks.push((serial, chunk)) {
                            Ok(()) => {
                                state.deadlock_state.record_q1_push();
                                return true;
                            }
                            Err((serial, chunk)) => {
                                worker.held_chunk = Some((serial, chunk));
                                return true;
                            }
                        }
                    }
                    Err(e) => {
                        state.set_error(e);
                        return false;
                    }
                }
            }
            StreamReader::Decompressed(r) => {
                match read_n_fastq_records(r, state.records_per_batch) {
                    Ok((data, offsets, at_eof)) => {
                        if offsets.len() <= 1 {
                            // No complete records — undo counter
                            state.batch_counters[stream_idx].fetch_sub(1, Ordering::Relaxed);
                            if at_eof {
                                state.stream_eof[stream_idx].store(true, Ordering::Release);
                                if state.stream_eof.iter().all(|f| f.load(Ordering::Acquire)) {
                                    state.read_done.store(true, Ordering::Release);
                                }
                            }
                            continue;
                        }
                        // Set EOF flag now if at_eof, even though we have records
                        // to push. This avoids an extra wasted read call and
                        // ensures read_done is set promptly for pipeline shutdown.
                        if at_eof {
                            state.stream_eof[stream_idx].store(true, Ordering::Release);
                            if state.stream_eof.iter().all(|f| f.load(Ordering::Acquire)) {
                                state.read_done.store(true, Ordering::Release);
                            }
                        }
                        let serial = state.batches_read.fetch_add(1, Ordering::Release);
                        let chunk =
                            PerStreamChunk { stream_idx, batch_num, data, offsets: Some(offsets) };
                        match state.q0_chunks.push((serial, chunk)) {
                            Ok(()) => {
                                state.deadlock_state.record_q1_push();
                                return true;
                            }
                            Err((serial, chunk)) => {
                                worker.held_chunk = Some((serial, chunk));
                                return true;
                            }
                        }
                    }
                    Err(e) => {
                        state.set_error(e);
                        return false;
                    }
                }
            }
        }
    }
    false // no stream available
}

/// Try to decompress a per-stream chunk (Step 2).
///
/// For BGZF: decompresses raw blocks. For gzip/plain: passthrough.
/// No completion flag — the Pair step tracks completion via count-based
/// checking (`chunks_paired == batches_read`), following the BAM pipeline pattern.
fn fastq_try_step_decompress<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
    worker: &mut FastqWorkerState<P>,
) -> bool {
    // Priority 1: Try to advance held decompressed chunk
    if let Some((serial, held)) = worker.held_decompressed_chunk.take() {
        match state.q1_decompressed.push((serial, held)) {
            Ok(()) => {
                state.deadlock_state.record_q2_push();
            }
            Err((serial, held)) => {
                worker.held_decompressed_chunk = Some((serial, held));
                return false;
            }
        }
    }

    // Priority 2: Check for errors
    if state.has_error() {
        return false;
    }

    // Priority 3: Check if output queue has space
    if state.q1_decompressed.is_full() {
        return false;
    }

    // Priority 4: Pop from input queue
    let Some((serial, chunk)) = state.q0_chunks.pop() else {
        if let Some(stats) = state.stats() {
            stats.record_queue_empty(1);
        }
        return false;
    };
    state.deadlock_state.record_q1_pop();

    // Priority 5: Decompress or passthrough
    let decompressed = if chunk.offsets.is_some() {
        // Gzip/plain: passthrough — already has record boundaries
        chunk
    } else {
        // BGZF: decompress raw blocks
        match decompress_bgzf_chunk(&chunk.data, &mut worker.decompressor) {
            Ok(decompressed_data) => PerStreamChunk {
                stream_idx: chunk.stream_idx,
                batch_num: chunk.batch_num,
                data: decompressed_data,
                offsets: None, // BGZF: boundary finding done in Pair step
            },
            Err(e) => {
                state.set_error(e);
                return false;
            }
        }
    };

    // Priority 6: Push result
    match state.q1_decompressed.push((serial, decompressed)) {
        Ok(()) => {
            state.deadlock_state.record_q2_push();
            true
        }
        Err((serial, chunk)) => {
            worker.held_decompressed_chunk = Some((serial, chunk));
            true // did work (decompressed)
        }
    }
}

// ============================================================================
// BlockParseFast Step (parallel, BGZF path)
// ============================================================================

/// `BlockParseFast` step: in parallel, parse all complete records from a single
/// decompressed BGZF chunk and emit a `BlockParsed` item to `q2_block_parsed`.
///
/// Each worker independently pops one `PerStreamChunk` from `q1_decompressed`,
/// detects where the chunk's phase boundary lies (using `detect_prefix_end` /
/// `detect_suffix_start`), runs SIMD-accelerated `find_record_offsets` on the
/// middle portion, parses those records, and pushes a `BlockParsed` item.
///
/// The cross-block records (prefix + suffix fragments) are forwarded to the
/// serial `BlockMerge` step rather than parsed here, because they depend on
/// data from the adjacent block.
fn fastq_try_step_block_parse<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
    worker: &mut FastqWorkerState<P>,
) -> bool {
    if state.has_error() {
        return false;
    }

    // Priority 1: Try to push any held BlockParsed item.
    if let Some(held) = worker.held_block_parsed.take() {
        match state.q2_block_parsed.push(held) {
            Ok(()) => {
                // Count here: the item is now in Q2b (was held from a previous call).
                // Heap bytes were already added when the item was created.
                state.chunks_block_parsed.fetch_add(1, Ordering::Release);
                state.deadlock_state.record_q2_5_push();
            }
            Err(held) => {
                worker.held_block_parsed = Some(held);
                return false;
            }
        }
    }

    // Check if BlockMerge is already done (all chunks processed).
    if state.block_merge_done.load(Ordering::Relaxed) {
        return false;
    }

    // Priority 2: Check output queue capacity and memory backpressure.
    if state.q2_block_parsed.is_full() || state.is_q2_block_parsed_memory_high() {
        return false;
    }

    // Priority 3: Pop a decompressed chunk from q1.
    let Some((serial, chunk)) = state.q1_decompressed.pop() else {
        return false;
    };
    state.deadlock_state.record_q2_pop();
    let _ = serial; // serial is not used for ordering in the BGZF path (block_idx is used)

    let stream_idx = chunk.stream_idx;
    let block_idx = chunk.batch_num;

    // For gzip/plain chunks (offsets pre-computed), delegate to the existing logic.
    // BGZF chunks have offsets=None.
    if chunk.offsets.is_some() {
        // Should not happen in BGZF mode, but guard against it.
        state.set_error(io::Error::new(
            io::ErrorKind::InvalidData,
            "BlockParseFast received a gzip/plain chunk in BGZF mode",
        ));
        return false;
    }

    let data = &chunk.data;

    // Detect phase: where does the first complete record start and the last end?
    let prefix_end = detect_prefix_end(data);
    let suffix_start = if prefix_end >= data.len() {
        data.len() // entire buffer is prefix (no complete records)
    } else {
        detect_suffix_start(&data[prefix_end..]) + prefix_end
    };

    let prefix_bytes = data[..prefix_end].to_vec();
    let suffix_bytes = data[suffix_start..].to_vec();
    let middle = &data[prefix_end..suffix_start];

    // Parse complete records from the middle using SIMD-accelerated offset detection.
    let offsets = fgumi_simd_fastq::find_record_offsets(middle);
    let num_records = offsets.len().saturating_sub(1);
    let mut records = Vec::with_capacity(num_records);
    for i in 0..num_records {
        let start = offsets[i];
        let end = offsets[i + 1];
        if start >= end || start >= middle.len() {
            continue;
        }
        match FastqRecord::from_slice(&middle[start..end.min(middle.len())]) {
            Ok(rec) => records.push(rec),
            Err(e) => {
                state.set_error(e);
                return false;
            }
        }
    }

    let block_parsed = BlockParsed { block_idx, stream_idx, records, prefix_bytes, suffix_bytes };

    // Track heap bytes at creation time. BlockMerge subtracts when it pops.
    let heap_bytes = block_parsed.estimate_heap_size() as u64;
    state.q2_block_parsed_heap_bytes.fetch_add(heap_bytes, Ordering::Release);

    match state.q2_block_parsed.push(block_parsed) {
        Ok(()) => {
            // Count here: the item is now in Q2b.
            state.chunks_block_parsed.fetch_add(1, Ordering::Release);
            state.deadlock_state.record_q2_5_push();
            true
        }
        Err(held) => {
            // Don't count yet — item is held. It will be counted when pushed in Priority 1.
            // Heap bytes already tracked above.
            worker.held_block_parsed = Some(held);
            true // did work (parsed the block)
        }
    }
}

// ============================================================================
// BlockMerge Step (serial via try_lock, BGZF path)
// ============================================================================

/// Result of [`drain_exhausted_stream`].
enum DrainResult {
    /// All available blocks drained successfully.
    Ok { did_work: bool, batches_this_call: usize },
    /// Output queue is full; caller must store the held data in `worker.held_parsed`.
    HeldParsed(u64, Vec<FastqTemplate>, usize),
    /// A stitching error occurred.
    Error(io::Error),
}

/// Drain remaining blocks from a single stream when the other stream is exhausted.
///
/// This is the shared implementation for both "R1 exhausted, drain R2" and
/// "R2 exhausted, drain R1" paths. The `drain_r1` flag selects which stream
/// to drain: `true` drains R1 blocks (pairing with R2 surplus), `false` drains
/// R2 blocks (pairing with R1 surplus).
///
/// Returns [`DrainResult::HeldParsed`] if the output queue is full (caller stores
/// it in `worker.held_parsed`), [`DrainResult::Error`] if a stitching error
/// occurs, or [`DrainResult::Ok`] with updated counters.
fn drain_exhausted_stream<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
    merge: &mut BlockMergeState,
    drain_r1: bool,
    mut did_work: bool,
    mut batches_this_call: usize,
) -> DrainResult {
    let (pending, suffix, surplus, other_surplus, next) = if drain_r1 {
        (
            &mut merge.r1_pending,
            &mut merge.r1_suffix_bytes,
            &mut merge.r1_surplus,
            &mut merge.r2_surplus,
            &mut merge.r1_next,
        )
    } else {
        (
            &mut merge.r2_pending,
            &mut merge.r2_suffix_bytes,
            &mut merge.r2_surplus,
            &mut merge.r1_surplus,
            &mut merge.r2_next,
        )
    };

    while batches_this_call < MAX_BATCHES_PER_LOCK {
        let block_next = *next;
        let Some(block) = pending.remove(&block_next) else {
            break;
        };
        merge.pending_heap_bytes =
            merge.pending_heap_bytes.saturating_sub(block.estimate_heap_size() as u64);
        *next += 1;

        let cross = match stitch_cross_block_record(suffix, &block.prefix_bytes) {
            Ok(rec) => rec,
            Err(e) => return DrainResult::Error(e),
        };

        let mut all_this: Vec<FastqRecord> = std::mem::take(surplus);
        if let Some(rec) = cross {
            all_this.push(rec);
        }
        all_this.extend(block.records);
        *suffix = block.suffix_bytes;

        // Pair with any available surplus from the other stream.
        let all_other: Vec<FastqRecord> = std::mem::take(other_surplus);
        let pair_count = all_this.len().min(all_other.len());
        if pair_count > 0 {
            // Ensure R1 is always first in the template regardless of which stream
            // we're draining.
            let (mut r1_vec, mut r2_vec) =
                if drain_r1 { (all_this, all_other) } else { (all_other, all_this) };

            let templates: Vec<FastqTemplate> = r1_vec
                .drain(..pair_count)
                .zip(r2_vec.drain(..pair_count))
                .map(|(r1, r2)| {
                    let name = r1.name().to_vec();
                    FastqTemplate { name, records: vec![r1, r2] }
                })
                .collect();

            // Restore surplus to the correct fields.
            if drain_r1 {
                *surplus = r1_vec;
                *other_surplus = r2_vec;
            } else {
                *surplus = r2_vec;
                *other_surplus = r1_vec;
            }

            let serial = merge.serial_out;
            merge.serial_out += 1;
            let count = templates.len();

            match state.output.groups.push((serial, templates)) {
                Ok(()) => {
                    state.total_templates_pushed.fetch_add(count as u64, Ordering::Release);
                    if let Some(stats) = state.stats() {
                        stats.groups_produced.fetch_add(count as u64, Ordering::Relaxed);
                    }
                    state.deadlock_state.record_q4_push();
                    did_work = true;
                }
                Err((_serial, returned)) => {
                    return DrainResult::HeldParsed(serial, returned, count);
                }
            }
        } else {
            *surplus = all_this;
            *other_surplus = all_other;
        }
        batches_this_call += 1;
    }

    DrainResult::Ok { did_work, batches_this_call }
}

/// `BlockMerge` step: serial step that assembles `BlockParsed` items into
/// `FastqTemplate`s and pushes them to the output groups queue.
///
/// This step uses `try_lock` so only one worker runs it at a time, matching
/// the serial semantics needed for in-order cross-block record stitching and
/// R1/R2 pairing.
///
/// Returns `(did_work, had_contention)`.
#[allow(clippy::too_many_lines)]
fn fastq_try_step_block_merge<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
    worker: &mut FastqWorkerState<P>,
) -> (bool, bool) {
    if state.has_error() {
        return (false, false);
    }

    // Priority 1: Advance held templates first (same pattern as FindBoundaries).
    let mut did_work = false;
    if let Some((serial, held_templates, count)) = worker.held_parsed.take() {
        match state.output.groups.push((serial, held_templates)) {
            Ok(()) => {
                state.total_templates_pushed.fetch_add(count as u64, Ordering::Release);
                if let Some(stats) = state.stats() {
                    stats.groups_produced.fetch_add(count as u64, Ordering::Relaxed);
                }
                state.deadlock_state.record_q4_push();
                // In the BGZF path, we don't use batches_parsed/batches_grouped.
                did_work = true;
            }
            Err(returned) => {
                worker.held_parsed = Some((serial, returned.1, count));
                return (false, false);
            }
        }
    }

    // Check if already done.
    if state.block_merge_done.load(Ordering::Relaxed) {
        return (did_work, false);
    }

    // Check output capacity.
    if state.output.groups.is_full() {
        return (did_work, false);
    }

    // Acquire the merge state lock (try_lock for serial execution).
    let Some(mut merge) = state.block_merge_state.try_lock() else {
        return (did_work, true); // contention
    };

    let num_streams = state.num_streams;

    // Determine whether to drain q2 into the pending maps.
    //
    // When the pending maps are large AND we can already process the next
    // in-order blocks, we SKIP draining.  This lets q2 fill up, creating
    // natural backpressure on BlockParseFast workers via queue fullness —
    // bounding total memory without risking deadlock.
    //
    // When we CANNOT process (next in-order blocks missing), we ALWAYS
    // drain so the needed blocks can reach the pending maps.
    let can_process = if num_streams == 1 {
        merge.r1_pending.contains_key(&merge.r1_next)
    } else {
        merge.r1_pending.contains_key(&merge.r1_next)
            && merge.r2_pending.contains_key(&merge.r2_next)
    };
    let within_limit = merge.pending_heap_bytes < PENDING_BACKPRESSURE_BYTES;

    let mut drained = 0;
    if within_limit || !can_process {
        while let Some(block) = state.q2_block_parsed.pop() {
            let heap_bytes = block.estimate_heap_size() as u64;
            state.q2_block_parsed_heap_bytes.fetch_sub(heap_bytes, Ordering::Release);
            merge.pending_heap_bytes += heap_bytes;
            state.deadlock_state.record_q2_5_pop();
            state.blocks_merged.fetch_add(1, Ordering::Release);
            if block.stream_idx == 0 {
                merge.r1_pending.insert(block.block_idx, block);
            } else {
                merge.r2_pending.insert(block.block_idx, block);
            }
            drained += 1;
        }
    }
    if drained == 0 && merge.r1_pending.is_empty() && merge.r2_pending.is_empty() {
        // Nothing to do. Check for completion.
        let all_chunks_block_parsed = state.read_done.load(Ordering::Acquire)
            && state.chunks_block_parsed.load(Ordering::Acquire)
                == state.batches_read.load(Ordering::Acquire);
        if all_chunks_block_parsed && merge.is_empty() {
            state.block_merge_done.store(true, Ordering::Release);
            state.parse_done.store(true, Ordering::Release);
            state.group_done.store(true, Ordering::Release);
        }
        return (did_work, false);
    }

    let mut batches_this_call = 0;

    if num_streams == 1 {
        // Single-stream (R1 only): process R1 blocks independently.
        while batches_this_call < MAX_BATCHES_PER_LOCK {
            let r1_next = merge.r1_next;
            let Some(r1_block) = merge.r1_pending.remove(&r1_next) else {
                break;
            };
            merge.pending_heap_bytes =
                merge.pending_heap_bytes.saturating_sub(r1_block.estimate_heap_size() as u64);
            merge.r1_next += 1;

            // Stitch cross-block record.
            let cross =
                match stitch_cross_block_record(&merge.r1_suffix_bytes, &r1_block.prefix_bytes) {
                    Ok(rec) => rec,
                    Err(e) => {
                        state.set_error(e);
                        return (true, false);
                    }
                };

            // Build full record list for this block.
            let mut all_records: Vec<FastqRecord> = std::mem::take(&mut merge.r1_surplus);
            if let Some(rec) = cross {
                all_records.push(rec);
            }
            all_records.extend(r1_block.records);
            merge.r1_suffix_bytes = r1_block.suffix_bytes;

            // Emit templates (single-stream: each record is its own template).
            let templates: Vec<FastqTemplate> = all_records
                .into_iter()
                .map(|r| {
                    let name = r.name().to_vec();
                    FastqTemplate { name, records: vec![r] }
                })
                .collect();

            let serial = merge.serial_out;
            merge.serial_out += 1;
            let count = templates.len();

            match state.output.groups.push((serial, templates)) {
                Ok(()) => {
                    state.total_templates_pushed.fetch_add(count as u64, Ordering::Release);
                    if let Some(stats) = state.stats() {
                        stats.groups_produced.fetch_add(count as u64, Ordering::Relaxed);
                    }
                    state.deadlock_state.record_q4_push();
                    did_work = true;
                }
                Err((serial, returned)) => {
                    worker.held_parsed = Some((serial, returned, count));
                    did_work = true;
                    break;
                }
            }
            batches_this_call += 1;
        }
    } else {
        // Paired-end: pair R1 and R2 blocks by block_idx.
        while batches_this_call < MAX_BATCHES_PER_LOCK {
            let r1_next = merge.r1_next;
            let r2_next = merge.r2_next;
            if !merge.r1_pending.contains_key(&r1_next) || !merge.r2_pending.contains_key(&r2_next)
            {
                break;
            }

            let r1_block = merge.r1_pending.remove(&r1_next).expect("just checked");
            let r2_block = merge.r2_pending.remove(&r2_next).expect("just checked");
            merge.pending_heap_bytes = merge.pending_heap_bytes.saturating_sub(
                (r1_block.estimate_heap_size() + r2_block.estimate_heap_size()) as u64,
            );
            merge.r1_next += 1;
            merge.r2_next += 1;

            // Stitch cross-block records.
            let r1_cross =
                match stitch_cross_block_record(&merge.r1_suffix_bytes, &r1_block.prefix_bytes) {
                    Ok(rec) => rec,
                    Err(e) => {
                        state.set_error(e);
                        return (true, false);
                    }
                };
            let r2_cross =
                match stitch_cross_block_record(&merge.r2_suffix_bytes, &r2_block.prefix_bytes) {
                    Ok(rec) => rec,
                    Err(e) => {
                        state.set_error(e);
                        return (true, false);
                    }
                };

            // Build full record lists for this round.
            let mut all_r1: Vec<FastqRecord> = std::mem::take(&mut merge.r1_surplus);
            if let Some(rec) = r1_cross {
                all_r1.push(rec);
            }
            all_r1.extend(r1_block.records);

            let mut all_r2: Vec<FastqRecord> = std::mem::take(&mut merge.r2_surplus);
            if let Some(rec) = r2_cross {
                all_r2.push(rec);
            }
            all_r2.extend(r2_block.records);

            merge.r1_suffix_bytes = r1_block.suffix_bytes;
            merge.r2_suffix_bytes = r2_block.suffix_bytes;

            // Zip min(r1, r2) pairs into templates, moving records (no clone).
            let pair_count = all_r1.len().min(all_r2.len());
            let templates: Vec<FastqTemplate> = all_r1
                .drain(..pair_count)
                .zip(all_r2.drain(..pair_count))
                .map(|(r1, r2)| {
                    let name = r1.name().to_vec();
                    FastqTemplate { name, records: vec![r1, r2] }
                })
                .collect();

            // Save surplus for the next round (drain left the remainder).
            merge.r1_surplus = all_r1;
            merge.r2_surplus = all_r2;

            let serial = merge.serial_out;
            merge.serial_out += 1;
            let count = templates.len();

            match state.output.groups.push((serial, templates)) {
                Ok(()) => {
                    state.total_templates_pushed.fetch_add(count as u64, Ordering::Release);
                    if let Some(stats) = state.stats() {
                        stats.groups_produced.fetch_add(count as u64, Ordering::Relaxed);
                    }
                    state.deadlock_state.record_q4_push();
                    did_work = true;
                }
                Err((serial, returned)) => {
                    worker.held_parsed = Some((serial, returned, count));
                    did_work = true;
                    break;
                }
            }
            batches_this_call += 1;
        }

        // Drain remaining blocks when one stream is exhausted.
        //
        // R1 and R2 BGZF files can have different numbers of blocks (different
        // compressed sizes), so after the paired loop above one stream may have
        // blocks left in its pending map. Once the shorter stream is fully
        // consumed (EOF + all its blocks have been block-parsed and merged), we
        // drain the longer stream's remaining blocks here so the pipeline can
        // complete.
        if worker.held_parsed.is_none() && batches_this_call < MAX_BATCHES_PER_LOCK {
            let r1_total = state.batch_counters[0].load(Ordering::Acquire);
            let r2_total = state.batch_counters[1].load(Ordering::Acquire);
            let r1_exhausted = state.stream_eof[0].load(Ordering::Acquire)
                && merge.r1_next == r1_total
                && !merge.r1_pending.contains_key(&merge.r1_next);
            let r2_exhausted = state.stream_eof[1].load(Ordering::Acquire)
                && merge.r2_next == r2_total
                && !merge.r2_pending.contains_key(&merge.r2_next);

            // Drain extra R2 blocks when R1 is exhausted.
            if r1_exhausted && !merge.r2_pending.is_empty() {
                match drain_exhausted_stream(state, &mut merge, false, did_work, batches_this_call)
                {
                    DrainResult::Ok { did_work: dw, batches_this_call: bc } => {
                        did_work = dw;
                        batches_this_call = bc;
                    }
                    DrainResult::HeldParsed(serial, templates, count) => {
                        worker.held_parsed = Some((serial, templates, count));
                        did_work = true;
                    }
                    DrainResult::Error(e) => {
                        state.set_error(e);
                        return (true, false);
                    }
                }
            }

            // Drain extra R1 blocks when R2 is exhausted.
            if worker.held_parsed.is_none() && r2_exhausted && !merge.r1_pending.is_empty() {
                match drain_exhausted_stream(state, &mut merge, true, did_work, batches_this_call) {
                    DrainResult::Ok { did_work: dw, batches_this_call: _bc } => {
                        did_work = dw;
                    }
                    DrainResult::HeldParsed(serial, templates, count) => {
                        worker.held_parsed = Some((serial, templates, count));
                        did_work = true;
                    }
                    DrainResult::Error(e) => {
                        state.set_error(e);
                        return (true, false);
                    }
                }
            }
        }
    }

    // Check for completion: all chunks processed, queue drained, and merge state empty.
    let all_chunks_block_parsed = state.read_done.load(Ordering::Acquire)
        && state.chunks_block_parsed.load(Ordering::Acquire)
            == state.batches_read.load(Ordering::Acquire);
    if all_chunks_block_parsed
        && state.q2_block_parsed.is_empty()
        && merge.r1_pending.is_empty()
        && merge.r2_pending.is_empty()
        && merge.r1_surplus.is_empty()
        && merge.r2_surplus.is_empty()
        && merge.r1_suffix_bytes.is_empty()
        && merge.r2_suffix_bytes.is_empty()
        && worker.held_parsed.is_none()
    {
        state.block_merge_done.store(true, Ordering::Release);
        state.parse_done.store(true, Ordering::Release);
        state.group_done.store(true, Ordering::Release);
    }

    (did_work, false)
}

// ============================================================================
// Pair Step (FindBoundaries)
// ============================================================================

/// Pair step: assemble per-stream chunks into multi-stream boundary batches.
///
/// Accumulates decompressed per-stream chunks by `batch_num`. When all non-EOF
/// streams have delivered a given `batch_num`, the chunks are assembled and
/// record boundaries are found (for BGZF) or directly used (for gzip).
///
/// Returns (`did_work`, `had_contention`).
fn fastq_try_step_find_boundaries<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
    worker: &mut FastqWorkerState<P>,
) -> (bool, bool) {
    if state.has_error() {
        return (false, false);
    }

    // Priority 1: Try to advance held boundary batch BEFORE checking boundaries_done.
    // This is critical: the last batch may have gone to held_boundaries in the same
    // call that set boundaries_done=true (pair was emptied). If we checked boundaries_done
    // first, the held batch would never be pushed to q2_5, deadlocking the pipeline.
    let mut did_work = false;
    if let Some((serial, held)) = worker.held_boundaries.take() {
        let boundary_heap_size = held.estimate_heap_size();
        match state.q2_5_boundaries.push((serial, held)) {
            Ok(()) => {
                state
                    .q2_5_boundaries_heap_bytes
                    .fetch_add(boundary_heap_size as u64, Ordering::Relaxed);
                // Note: batches_boundaries_found was already incremented at serial
                // assignment time (fetch_add in Priority 5), not here.
                state.deadlock_state.record_q2_5_push();
                did_work = true;
            }
            Err((serial, held)) => {
                worker.held_boundaries = Some((serial, held));
                return (false, false);
            }
        }
    }

    // Now safe to check boundaries_done (held items already handled above).
    if state.boundaries_done.load(Ordering::Relaxed) {
        return (did_work, false);
    }

    // Priority 2: Check if output queue has space
    if state.q2_5_boundaries.is_full() {
        return (did_work, false);
    }

    // Priority 3: Acquire pair state lock
    let Some(mut pair) = state.pair_state.try_lock() else {
        return (did_work, true); // Contention
    };

    // Priority 4: Drain q1_decompressed into pair buffer
    while let Some((_, chunk)) = state.q1_decompressed.pop() {
        state.deadlock_state.record_q2_pop();
        state.chunks_paired.fetch_add(1, Ordering::Release);
        pair.insert(chunk);
    }

    // Check if ALL chunks from Read have arrived at the Pair.
    // Same pattern as BAM: `read_done && chunks_paired == batches_read`.
    let all_arrived = state.read_done.load(Ordering::Acquire)
        && state.chunks_paired.load(Ordering::Acquire)
            == state.batches_read.load(Ordering::Acquire);

    // Priority 5: Try to emit complete batches.
    // Cap batches processed per lock hold to avoid starving other workers (mirrors BAM pipeline).
    let mut batches_this_call = 0;
    while batches_this_call < MAX_BATCHES_PER_LOCK {
        let Some(chunks) = pair.try_pop_complete(all_arrived) else { break };
        // Atomically assign a unique serial. This must be fetch_add (not load)
        // because the held_boundaries path can race: Worker A creates a batch
        // but push fails (goes to held_boundaries without incrementing), then
        // Worker B enters this loop and would get the same serial from load().
        let serial = state.batches_boundaries_found.fetch_add(1, Ordering::Release);

        let boundary_batch = if chunks.iter().all(|c| c.offsets.is_some()) {
            // Gzip path: all chunks already have record boundaries.
            let streams: Vec<FastqStreamBoundaries> = chunks
                .into_iter()
                .map(|c| FastqStreamBoundaries {
                    stream_idx: c.stream_idx,
                    data: c.data,
                    offsets: c.offsets.expect("gzip chunks must have pre-computed offsets"),
                })
                .collect();
            align_stream_records(streams, serial)
        } else {
            // BGZF path: need to find record boundaries in decompressed data.
            let decompressed = FastqDecompressedBatch {
                chunks: chunks
                    .into_iter()
                    .map(|c| FastqDecompressedChunk { stream_idx: c.stream_idx, data: c.data })
                    .collect(),
                serial,
            };
            match FastqFormat::find_boundaries(&state.boundary_state, decompressed) {
                Ok(batch) => batch,
                Err(e) => {
                    state.set_error(e);
                    return (true, false);
                }
            }
        };

        let boundary_heap_size = boundary_batch.estimate_heap_size();
        match state.q2_5_boundaries.push((serial, boundary_batch)) {
            Ok(()) => {
                state
                    .q2_5_boundaries_heap_bytes
                    .fetch_add(boundary_heap_size as u64, Ordering::Relaxed);
                // Note: batches_boundaries_found was already incremented by
                // fetch_add above when the serial was assigned.
                state.deadlock_state.record_q2_5_push();
                did_work = true;
            }
            Err((serial, batch)) => {
                worker.held_boundaries = Some((serial, batch));
                did_work = true;
                break; // output queue full, stop emitting
            }
        }
        batches_this_call += 1;
    }

    // Completion: all chunks paired and all batches emitted.
    if all_arrived && pair.is_empty() {
        state.boundaries_done.store(true, Ordering::Release);
    }

    (did_work, false)
}

/// Try to parse FASTQ records (Step 2.75 - parallel).
///
/// This step takes boundary batches and constructs `FastqRecord` objects,
/// then creates templates directly by zipping records at matching positions.
/// This bypasses the Group step entirely, eliminating all lock contention.
///
/// Uses held-item pattern (like all other pipeline steps): if the output queue
/// (`output.groups`) is full, the parsed templates are stored in
/// `worker.held_parsed` and the function returns immediately. This allows the
/// thread to work on downstream steps (especially Write) to drain backpressure,
/// preventing deadlocks at low thread counts.
///
/// **This is the KEY PARALLEL STEP** that fixes the t8 scaling bottleneck.
fn fastq_try_step_parse<R: BufRead + Send, P: Send + MemoryEstimate>(
    state: &FastqPipelineState<R, P>,
    worker: &mut FastqWorkerState<P>,
) -> bool {
    if state.parse_done.load(Ordering::Relaxed) || state.has_error() {
        return false;
    }

    // Priority 1: Try to advance held parsed templates
    if let Some((serial, held_templates, count)) = worker.held_parsed.take() {
        match state.output.groups.push((serial, held_templates)) {
            Ok(()) => {
                #[cfg(feature = "memory-debug")]
                {
                    let q4_heap: u64 = 0; // already tracked when first parsed
                    state.output.groups_heap_bytes.fetch_add(q4_heap, Ordering::AcqRel);
                }
                state.total_templates_pushed.fetch_add(count as u64, Ordering::Release);
                if let Some(stats) = state.stats() {
                    stats.groups_produced.fetch_add(count as u64, Ordering::Relaxed);
                }
                state.deadlock_state.record_q4_push();
                state.batches_parsed.fetch_add(1, Ordering::Release);
                state.batches_grouped.fetch_add(1, Ordering::Release);
                // Continue to try popping more input below
            }
            Err(returned) => {
                worker.held_parsed = Some((serial, returned.1, count));
                return false; // Output still full — let thread work on downstream steps
            }
        }
    }

    // Priority 2: Check output queue capacity BEFORE popping input to prevent data loss
    if state.output.groups.is_full() {
        return false;
    }

    // Priority 3: Pop from Q2.5 boundaries queue
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
            // Group is skipped — set group_done alongside parse_done
            state.group_done.store(true, Ordering::Release);
            log::trace!("PARSE: set parse_done=true, group_done=true");
        } else if let Some(stats) = state.stats() {
            // Record Q2.5 as extension of Q2
            stats.record_queue_empty(2);
        }
        return false;
    };
    state.deadlock_state.record_q2_5_pop();

    let input_heap_size = boundary_batch.estimate_heap_size();

    // Priority 4: Parse records — THIS IS THE KEY PARALLEL OPERATION
    match FastqFormat::parse_records(boundary_batch) {
        Ok(parsed_batch) => {
            // Only decrement input memory AFTER successful parse
            state.q2_5_boundaries_heap_bytes.fetch_sub(input_heap_size as u64, Ordering::Relaxed);

            // Create templates by zipping records at matching positions
            let templates = match create_templates_from_streams(parsed_batch.streams) {
                Ok(t) => t,
                Err(e) => {
                    state.set_error(e);
                    return false;
                }
            };

            let count = templates.len();

            // Priority 5: Push templates to Q3 using held-item pattern
            match state.output.groups.push((serial, templates)) {
                Ok(()) => {
                    #[cfg(feature = "memory-debug")]
                    {
                        // Heap size tracking for pushed templates is not yet implemented;
                        // queue memory is tracked through estimates only for now.
                        let q4_heap: u64 = 0;
                        state.output.groups_heap_bytes.fetch_add(q4_heap, Ordering::AcqRel);
                    }
                    state.total_templates_pushed.fetch_add(count as u64, Ordering::Release);
                    if let Some(stats) = state.stats() {
                        stats.groups_produced.fetch_add(count as u64, Ordering::Relaxed);
                    }
                    state.deadlock_state.record_q4_push();
                    state.batches_parsed.fetch_add(1, Ordering::Release);
                    state.batches_grouped.fetch_add(1, Ordering::Release);
                    true
                }
                Err(returned) => {
                    // Output queue full — store in held_parsed for retry on next call.
                    // This allows the thread to work on downstream steps (Write, Compress, etc.)
                    // instead of spinning, preventing deadlocks at low thread counts.
                    worker.held_parsed = Some((serial, returned.1, count));
                    true // Did work (parsed the batch)
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
            let records = streams.pop().expect("streams is non-empty in single-end branch").records;
            Ok(records
                .into_iter()
                .map(|r| {
                    let name = r.name().to_vec();
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
            let r1_records = drain.next().expect("sorted streams must have R1 at index 0").records;
            let r2_records = drain.next().expect("sorted streams must have R2 at index 1").records;

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
                    let name = r1.name().to_vec();
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
    // Memory backpressure is always enforced (including during draining) to
    // prevent OOM.  The slot-based is_full() check is sufficient to guarantee
    // forward progress: Serialize drains the processed queue -> slots free ->
    // Process resumes.
    // =========================================================================
    if state.output.processed.is_full() || state.is_q4_memory_high() {
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
        if state.output.processed.is_full() || state.is_q4_memory_high() {
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
    let batch = SerializedBatch {
        data: combined_data,
        record_count: total_record_count,
        secondary_data: None,
    };
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
///
/// Drains Q6 (compressed) into the write reorder buffer, then writes
/// consecutive batches in serial order.
///
/// Returns true if any data was actually written to the output file.
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

    // Drain Q6 into reorder buffer AND write all ready batches in single lock scope.
    let mut wrote_any = false;
    let q7_truly_empty;
    {
        let mut reorder = state.output.write_reorder.lock();

        // Drain Q6 into reorder buffer.
        while let Some((serial, batch)) = state.output.compressed.pop() {
            let heap_size = batch.estimate_heap_size();
            let q7_heap = heap_size as u64;
            state.q6_track_pop(q7_heap);
            reorder.insert_with_size(serial, batch, heap_size);
            state.output.write_reorder_state.add_heap_bytes(q7_heap);
            state.deadlock_state.record_q7_pop();
        }

        // Write in-order batches.
        while let Some((batch, heap_size)) = reorder.try_pop_next_with_size() {
            let mut batch_bytes: u64 = 0;
            for block in &batch.blocks {
                batch_bytes += block.data.len() as u64;
                if let Err(e) = output.write_all(&block.data) {
                    state.set_error(e);
                    return false;
                }
            }
            state.output.write_reorder_state.sub_heap_bytes(heap_size as u64);
            state.output.write_reorder_state.update_next_seq(reorder.next_seq());
            if let Some(stats) = state.stats() {
                stats.bytes_written.fetch_add(batch_bytes, Ordering::Relaxed);
            }
            let records_in_batch = batch.record_count;
            state.output.items_written.fetch_add(records_in_batch, Ordering::Relaxed);
            state.output.progress.log_if_needed(records_in_batch);
            wrote_any = true;
        }

        q7_truly_empty = reorder.is_empty();
    }

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
            if state.config.inputs_are_bgzf {
                // BGZF path: BlockMerge is the serial stitch+pair step.
                fastq_try_step_block_merge(state, worker)
            } else {
                // Gzip/plain path: original serial boundary-finding + pair assembly.
                fastq_try_step_find_boundaries(state, worker)
            }
        }
        PipelineStep::Decode => {
            if state.config.inputs_are_bgzf {
                // BGZF path: BlockParseFast is the parallel per-chunk parse step.
                let success = fastq_try_step_block_parse(state, worker);
                (success, false)
            } else {
                // Gzip/plain path: original parallel parse step.
                let success = fastq_try_step_parse(state, worker);
                (success, false)
            }
        }
        PipelineStep::Group => {
            // Group step is never active — synchronized mode creates templates
            // directly in the Decode/BlockMerge step
            (false, false)
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
            input_low: self.state.q0_chunks.len() < cap / 4,
            read_done,
            memory_high: !self.state.is_draining()
                && self.state.output.write_reorder_state.is_memory_high(),
            memory_drained: self.state.output.write_reorder_state.is_memory_drained(),
        }
    }

    fn check_drain_mode(&self) {
        let read_done = self.state.read_done.load(Ordering::Relaxed);
        if read_done && self.state.q0_chunks.is_empty() {
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
#[allow(clippy::too_many_lines, clippy::needless_pass_by_value)]
pub fn run_fastq_pipeline<P, PF, SF>(
    config: FastqPipelineConfig,
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

    // Write BAM header first
    log::debug!("run_fastq_pipeline: writing BAM header");
    write_bam_header(&mut output, header)?;
    log::debug!("run_fastq_pipeline: BAM header written successfully");

    // Create per-stream readers based on input type
    log::debug!(
        "run_fastq_pipeline: creating readers, decompressed_readers.is_some()={}, inputs_are_bgzf={}",
        decompressed_readers.is_some(),
        config.inputs_are_bgzf,
    );
    let stream_readers: Vec<StreamReader<Box<dyn BufRead + Send>>> =
        if let Some(readers) = decompressed_readers {
            // Gzip/plain: wrap each reader as StreamReader::Decompressed
            let num_readers = readers.len();
            log::debug!("run_fastq_pipeline: using {num_readers} Decompressed readers");
            readers.into_iter().map(StreamReader::Decompressed).collect()
        } else {
            // BGZF: open each file as StreamReader::Bgzf
            log::debug!("run_fastq_pipeline: using {} BGZF readers", fastq_paths.len());
            fastq_paths
                .iter()
                .map(|p| {
                    let file = File::open(p)?;
                    Ok(StreamReader::Bgzf(BufReader::with_capacity(256 * 1024, file)))
                })
                .collect::<io::Result<Vec<_>>>()?
        };

    // Create state
    log::debug!("run_fastq_pipeline: creating pipeline state");
    let state = Arc::new(FastqPipelineState::<Box<dyn BufRead + Send>, P>::new(
        config.clone(),
        stream_readers,
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
                        is_reader: thread_id < state.num_streams,
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
                // Log parse state for debugging (at trace level)
                if s.deadlock_state.is_enabled() {
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
        let blocks = read_raw_blocks(&mut reader, 10).expect("failed to read raw blocks");

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

        producer.join().expect("thread should not panic");
        consumer.join().expect("thread should not panic");

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
    fn test_config_defaults() {
        let config = FastqPipelineConfig::new(4, false, 6);
        assert!(!config.inputs_are_bgzf);
        assert_eq!(config.num_threads, 4);
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
        let boundary_batch =
            FastqFormat::find_boundaries(&boundary_state, batch).expect("find_boundaries failed");
        assert_eq!(boundary_batch.streams.len(), 2);
        // Each stream should have 2 complete records (offsets at 0, 19, 38)
        assert_eq!(boundary_batch.streams[0].offsets.len(), 3);
        assert_eq!(boundary_batch.streams[1].offsets.len(), 3);

        // Step 2: Parse records
        let parsed_batch =
            FastqFormat::parse_records(boundary_batch).expect("parse_records failed");
        assert_eq!(parsed_batch.streams.len(), 2);
        assert_eq!(parsed_batch.streams[0].stream_idx, 0);
        assert_eq!(parsed_batch.streams[1].stream_idx, 1);
        assert_eq!(parsed_batch.streams[0].records.len(), 2);
        assert_eq!(parsed_batch.streams[1].records.len(), 2);

        // Verify record contents
        assert_eq!(parsed_batch.streams[0].records[0].name(), b"read1");
        assert_eq!(parsed_batch.streams[0].records[0].sequence(), b"ACGT");
        assert_eq!(parsed_batch.streams[0].records[1].name(), b"read2");
        assert_eq!(parsed_batch.streams[1].records[0].name(), b"read1");
        assert_eq!(parsed_batch.streams[1].records[0].sequence(), b"TTTT");
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

        let boundary_batch1 =
            FastqFormat::find_boundaries(&boundary_state, batch1).expect("find_boundaries failed");
        assert_eq!(boundary_batch1.streams[0].offsets.len(), 2); // One complete record
        assert!(!boundary_state.stream_states[0].lock().leftover.is_empty()); // Leftover from incomplete

        // Second chunk: completes the record
        let batch2 = FastqDecompressedBatch {
            chunks: vec![FastqDecompressedChunk { stream_idx: 0, data: b"GG\n+\nJJJJ\n".to_vec() }],
            serial: 1,
        };

        let boundary_batch2 =
            FastqFormat::find_boundaries(&boundary_state, batch2).expect("find_boundaries failed");
        // Leftover + new data should form complete record
        assert!(boundary_batch2.streams[0].offsets.len() >= 2);
        assert!(boundary_state.stream_states[0].lock().leftover.is_empty()); // No more leftover

        // Parse both batches
        let parsed1 = FastqFormat::parse_records(boundary_batch1).expect("parse_records failed");
        let parsed2 = FastqFormat::parse_records(boundary_batch2).expect("parse_records failed");

        assert_eq!(parsed1.streams[0].records.len(), 1);
        assert_eq!(parsed1.streams[0].records[0].name(), b"read1");
        assert_eq!(parsed2.streams[0].records.len(), 1);
        assert_eq!(parsed2.streams[0].records[0].name(), b"read2");
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
                        let parsed = FastqFormat::parse_records(boundary_batch)
                            .expect("parse_records failed");
                        assert_eq!(parsed.streams[0].stream_idx, 0);
                        assert_eq!(parsed.streams[0].records.len(), 1);
                        assert_eq!(
                            String::from_utf8_lossy(parsed.streams[0].records[0].name()),
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

        let boundary_batch =
            FastqFormat::find_boundaries(&boundary_state, batch).expect("find_boundaries failed");

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

        let boundary_batch1 =
            FastqFormat::find_boundaries(&boundary_state, batch1).expect("find_boundaries failed");
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

        let boundary_batch2 =
            FastqFormat::find_boundaries(&boundary_state, batch2).expect("find_boundaries failed");

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
        let parsed = FastqFormat::parse_records(boundary_batch2).expect("parse_records failed");
        assert_eq!(
            parsed.streams[0].records[0].name(),
            b"r3",
            "First record should be r3 from leftover"
        );
        assert_eq!(parsed.streams[0].records[1].name(), b"r4", "Second record should be r4");
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

        let boundary_batch1 =
            FastqFormat::find_boundaries(&boundary_state, batch1).expect("find_boundaries failed");
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

        let boundary_batch2 =
            FastqFormat::find_boundaries(&boundary_state, batch2).expect("find_boundaries failed");

        // Both streams should be present (stream 0 from leftover, stream 1 from new chunk)
        assert_eq!(boundary_batch2.streams.len(), 2, "Both streams should be present");

        // Find stream 0 and stream 1 in the result (order may vary)
        let stream0 =
            boundary_batch2.streams.iter().find(|s| s.stream_idx == 0).expect("stream not found");
        let stream1 =
            boundary_batch2.streams.iter().find(|s| s.stream_idx == 1).expect("stream not found");

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

        let boundary_batch =
            FastqFormat::find_boundaries(&boundary_state, batch).expect("find_boundaries failed");

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

        let boundary_batch =
            FastqFormat::find_boundaries(&boundary_state, batch).expect("find_boundaries failed");

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

        let _ =
            FastqFormat::find_boundaries(&boundary_state, batch1).expect("find_boundaries failed");

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

        let boundary_batch2 =
            FastqFormat::find_boundaries(&boundary_state, batch2).expect("find_boundaries failed");

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

        let parsed = FastqFormat::parse_records(boundary_batch).expect("parse_records failed");
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
        assert_eq!(parsed.streams[0].records[0].sequence(), b"TTTT");
        assert_eq!(parsed.streams[1].records[0].sequence(), b"ACGT");
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
                records: vec![
                    FastqRecord::from_slice(b"@read1\nTTTT\n+\nJJJJ\n").unwrap(), // R2 sequence
                ],
            },
            FastqParsedStream {
                stream_idx: 0, // R1 comes second
                records: vec![
                    FastqRecord::from_slice(b"@read1\nACGT\n+\nIIII\n").unwrap(), // R1 sequence
                ],
            },
        ];

        let templates =
            create_templates_from_streams(streams).expect("create templates from streams");

        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].records.len(), 2);
        // Critical: R1 must be first (records[0]), R2 must be second (records[1])
        assert_eq!(
            templates[0].records[0].sequence(),
            b"ACGT",
            "First record in template must be R1 (ACGT), not R2"
        );
        assert_eq!(
            templates[0].records[1].sequence(),
            b"TTTT",
            "Second record in template must be R2 (TTTT), not R1"
        );
    }

    #[test]
    fn test_create_templates_from_correctly_ordered_streams() {
        // Verify correct behavior when streams are already in order
        let streams = vec![
            FastqParsedStream {
                stream_idx: 0,
                records: vec![FastqRecord::from_slice(b"@read1\nACGT\n+\nIIII\n").unwrap()],
            },
            FastqParsedStream {
                stream_idx: 1,
                records: vec![FastqRecord::from_slice(b"@read1\nTTTT\n+\nJJJJ\n").unwrap()],
            },
        ];

        let templates =
            create_templates_from_streams(streams).expect("create templates from streams");

        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].records[0].sequence(), b"ACGT", "R1 should be first");
        assert_eq!(templates[0].records[1].sequence(), b"TTTT", "R2 should be second");
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

        let boundary_batch1 =
            FastqFormat::find_boundaries(&boundary_state, batch1).expect("find_boundaries failed");
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

        let boundary_batch2 =
            FastqFormat::find_boundaries(&boundary_state, batch2).expect("find_boundaries failed");
        assert_eq!(boundary_batch2.streams.len(), 2);

        // The order of streams in boundary_batch2 may be [stream_idx=1, stream_idx=0]
        // (reversed). This is the trigger for the bug.
        let first_stream_idx = boundary_batch2.streams[0].stream_idx;
        let second_stream_idx = boundary_batch2.streams[1].stream_idx;

        // Parse records — stream_idx must be preserved
        let parsed = FastqFormat::parse_records(boundary_batch2).expect("parse_records failed");
        assert_eq!(parsed.streams[0].stream_idx, first_stream_idx);
        assert_eq!(parsed.streams[1].stream_idx, second_stream_idx);

        // Create templates — must produce correct R1/R2 ordering regardless
        // of the stream order in the Vec
        let templates =
            create_templates_from_streams(parsed.streams).expect("create templates from streams");
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, b"read2");
        assert_eq!(templates[0].records.len(), 2);

        // THE KEY ASSERTION: R1 data (GGGG from leftover) must be records[0],
        // R2 data (CCCC from new chunk) must be records[1].
        assert_eq!(
            templates[0].records[0].sequence(),
            b"GGGG",
            "records[0] must be R1 (stream 0) data, not R2"
        );
        assert_eq!(
            templates[0].records[1].sequence(),
            b"CCCC",
            "records[1] must be R2 (stream 1) data, not R1"
        );
    }

    // ========================================================================
    // read_n_fastq_records tests (helper function used by per-stream reader)
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

        let (buf, offsets, at_eof) =
            read_n_fastq_records(&mut cursor, 2).expect("read N FASTQ records");

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

        let (_, offsets, at_eof) =
            read_n_fastq_records(&mut cursor, 5).expect("read N FASTQ records");

        // Only 1 record available, requested 5
        assert_eq!(offsets.len(), 2); // 1 record
        assert!(at_eof);
    }

    #[test]
    fn test_read_n_fastq_records_empty_input() {
        let mut cursor = Cursor::new(Vec::<u8>::new());

        let (_, offsets, at_eof) =
            read_n_fastq_records(&mut cursor, 5).expect("read N FASTQ records");

        assert_eq!(offsets.len(), 1); // Just the initial 0
        assert!(at_eof);
    }

    // ========================================================================
    // Per-Stream Types Tests
    // ========================================================================

    #[test]
    fn test_pair_state_basic() {
        let pair = PairState::new(2);
        assert!(pair.is_empty());
    }

    #[test]
    fn test_pair_state_insert_and_pop() {
        let mut pair = PairState::new(2);

        // Insert both streams for batch 0
        pair.insert(PerStreamChunk {
            stream_idx: 0,
            batch_num: 0,
            data: b"data0".to_vec(),
            offsets: Some(vec![0, 5]),
        });
        assert!(pair.try_pop_complete(false).is_none()); // Not complete yet

        pair.insert(PerStreamChunk {
            stream_idx: 1,
            batch_num: 0,
            data: b"data1".to_vec(),
            offsets: Some(vec![0, 5]),
        });
        let chunks = pair.try_pop_complete(false).expect("try_pop_complete should succeed");
        assert_eq!(chunks.len(), 2);
        assert!(pair.is_empty());
    }

    #[test]
    fn test_pair_state_uneven_streams() {
        // Stream 0 produces 2 batches, stream 1 produces 1 batch.
        // With all_arrived=false, batch 1 won't emit (stream 1 missing).
        // With all_arrived=true, batch 1 emits with just stream 0's data.
        let mut pair = PairState::new(2);

        pair.insert(PerStreamChunk {
            stream_idx: 0,
            batch_num: 0,
            data: b"d00".to_vec(),
            offsets: Some(vec![0, 3]),
        });
        pair.insert(PerStreamChunk {
            stream_idx: 1,
            batch_num: 0,
            data: b"d10".to_vec(),
            offsets: Some(vec![0, 3]),
        });
        pair.insert(PerStreamChunk {
            stream_idx: 0,
            batch_num: 1,
            data: b"d01".to_vec(),
            offsets: Some(vec![0, 3]),
        });

        // Batch 0: complete
        let chunks = pair.try_pop_complete(false).expect("try_pop_complete should succeed");
        assert_eq!(chunks.len(), 2);

        // Batch 1: only stream 0 — not complete without all_arrived
        assert!(pair.try_pop_complete(false).is_none());

        // With all_arrived, batch 1 emits with just stream 0
        let chunks = pair.try_pop_complete(true).expect("try_pop_complete should succeed");
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0].stream_idx, 0);
        assert!(pair.is_empty());
    }

    /// Tests that the count-based completion logic (`chunks_paired == batches_read`)
    /// correctly drives `all_arrived` for `PairState::try_pop_complete`.
    ///
    /// This validates the simplified completion pattern borrowed from the BAM pipeline:
    /// no intermediate "done" flags — just compare a counter at the source (Read) with
    /// a counter at the destination (Pair).
    #[test]
    fn test_pair_state_count_based_completion() {
        let mut pair = PairState::new(2);

        // Simulate: 2 streams, stream 0 produces batches 0..2, stream 1 produces batch 0 only.
        // Total batches_read = 3, chunks arriving at pair counted by chunks_paired.

        // Insert batch 0 from both streams.
        pair.insert(PerStreamChunk {
            stream_idx: 0,
            batch_num: 0,
            data: b"s0b0".to_vec(),
            offsets: Some(vec![0, 4]),
        });
        pair.insert(PerStreamChunk {
            stream_idx: 1,
            batch_num: 0,
            data: b"s1b0".to_vec(),
            offsets: Some(vec![0, 4]),
        });

        // Simulate: batches_read=2 (still reading), chunks_paired=2.
        // all_arrived = read_done(false) && ... → false.
        let all_arrived = false;
        let chunks = pair.try_pop_complete(all_arrived).expect("try_pop_complete should succeed");
        assert_eq!(chunks.len(), 2);

        // Insert batch 1 from stream 0 only (stream 1 hit EOF earlier).
        pair.insert(PerStreamChunk {
            stream_idx: 0,
            batch_num: 1,
            data: b"s0b1".to_vec(),
            offsets: Some(vec![0, 4]),
        });

        // Simulate: batches_read=3, chunks_paired=2 (third chunk just inserted, not yet counted).
        // all_arrived = read_done(true) && 2 == 3 → false.
        let all_arrived = false;
        assert!(pair.try_pop_complete(all_arrived).is_none());

        // Simulate: chunks_paired catches up to 3.
        // all_arrived = read_done(true) && 3 == 3 → true.
        let all_arrived = true;
        let chunks = pair.try_pop_complete(all_arrived).expect("try_pop_complete should succeed");
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0].stream_idx, 0);
        assert!(pair.is_empty());
    }

    #[test]
    fn test_align_stream_records_equal() {
        let streams = vec![
            FastqStreamBoundaries {
                stream_idx: 0,
                data: b"@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n".to_vec(),
                offsets: vec![0, 18, 36],
            },
            FastqStreamBoundaries {
                stream_idx: 1,
                data: b"@r1\nTTTT\n+\nJJJJ\n@r2\nTTTT\n+\nJJJJ\n".to_vec(),
                offsets: vec![0, 18, 36],
            },
        ];
        let batch = align_stream_records(streams, 0);
        assert_eq!(batch.streams[0].offsets.len() - 1, 2);
        assert_eq!(batch.streams[1].offsets.len() - 1, 2);
    }

    #[test]
    fn test_align_stream_records_unequal() {
        let streams = vec![
            FastqStreamBoundaries {
                stream_idx: 0,
                data: b"@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n@r3\nACGT\n+\nIIII\n".to_vec(),
                offsets: vec![0, 18, 36, 54],
            },
            FastqStreamBoundaries {
                stream_idx: 1,
                data: b"@r1\nTTTT\n+\nJJJJ\n@r2\nTTTT\n+\nJJJJ\n".to_vec(),
                offsets: vec![0, 18, 36],
            },
        ];
        let batch = align_stream_records(streams, 0);
        // Both should be truncated to min(3, 2) = 2 records
        assert_eq!(batch.streams[0].offsets.len() - 1, 2);
        assert_eq!(batch.streams[1].offsets.len() - 1, 2);
    }

    // ========================================================================
    // ActiveSteps Configuration Tests
    // ========================================================================

    #[rstest]
    // Both gzip and BGZF: always 8 active steps (Group is always inactive)
    #[case::gzip(false)]
    #[case::bgzf(true)]
    fn test_active_steps(#[case] inputs_are_bgzf: bool) {
        let config = FastqPipelineConfig::new(4, inputs_are_bgzf, 1);
        let steps = config.active_steps();

        assert_eq!(
            steps.steps(),
            &[Read, Decompress, FindBoundaries, Decode, Process, Serialize, Compress, Write]
        );
        assert!(!steps.is_active(Group));
    }

    #[test]
    fn test_fastq_decompressed_batch_memory_estimate() {
        let mut data = Vec::with_capacity(2048);
        data.extend_from_slice(&[0u8; 100]);
        let mut chunks = Vec::with_capacity(4);
        chunks.push(FastqDecompressedChunk { stream_idx: 0, data });

        let batch = FastqDecompressedBatch { chunks, serial: 0 };
        let estimate = batch.estimate_heap_size();

        // Should use data capacity (2048) not len (100)
        assert!(estimate >= 2048, "estimate {estimate} should be >= 2048 (data capacity)");
        // Should include Vec<FastqDecompressedChunk> overhead
        let vec_overhead = 4 * std::mem::size_of::<FastqDecompressedChunk>();
        assert!(
            estimate >= 2048 + vec_overhead,
            "estimate {estimate} should include chunk Vec overhead"
        );
    }

    #[test]
    fn test_fastq_boundary_batch_memory_estimate() {
        let mut data = Vec::with_capacity(1024);
        data.extend_from_slice(&[0u8; 100]);
        let mut offsets = Vec::with_capacity(16);
        offsets.push(0usize);

        let mut streams = Vec::with_capacity(4);
        streams.push(FastqStreamBoundaries { stream_idx: 0, data, offsets });

        let batch = FastqBoundaryBatch { streams, serial: 0 };
        let estimate = batch.estimate_heap_size();

        let expected_min = 1024
            + 16 * std::mem::size_of::<usize>()
            + 4 * std::mem::size_of::<FastqStreamBoundaries>();
        assert!(
            estimate >= expected_min,
            "estimate {estimate} should be >= {expected_min} (capacities + overhead)"
        );
    }

    #[test]
    fn test_fastq_parsed_batch_memory_estimate() {
        use crate::fastq_parse::FastqRecord;

        // Single-allocation record: data = "@read1\nACGT\n+\nIIII\n"
        let record =
            FastqRecord::from_slice(b"@read1\nACGT\n+\nIIII\n").expect("valid FASTQ record");
        let mut records = Vec::with_capacity(8);
        records.push(record);

        let mut streams = Vec::with_capacity(4);
        streams.push(FastqParsedStream { stream_idx: 0, records });

        let batch = FastqParsedBatch { streams, serial: 0 };
        let estimate = batch.estimate_heap_size();

        // Estimate should include at least the data allocation + struct overheads
        let records_overhead = 8 * std::mem::size_of::<FastqRecord>();
        let streams_overhead = 4 * std::mem::size_of::<FastqParsedStream>();
        let expected_min = b"@read1\nACGT\n+\nIIII\n".len() + records_overhead + streams_overhead;
        assert!(
            estimate >= expected_min,
            "estimate {estimate} should be >= {expected_min} (data capacity + overhead)"
        );
    }

    // ========================================================================
    // FastqRecord::from_slice Tests
    // ========================================================================

    #[test]
    fn test_fastq_record_from_slice_normal() {
        use crate::fastq_parse::FastqRecord;
        let data = b"@read1\nACGT\n+\nIIII\n";
        let rec = FastqRecord::from_slice(data).expect("valid FASTQ record");
        assert_eq!(rec.name(), b"read1");
        assert_eq!(rec.sequence(), b"ACGT");
        assert_eq!(rec.quality(), b"IIII");
    }

    #[test]
    fn test_fastq_record_from_slice_at_in_quality() {
        // '@' and '+' characters in quality scores must not confuse the parser.
        use crate::fastq_parse::FastqRecord;
        let data = b"@read1\nACGT\n+\n@+!I\n";
        let rec = FastqRecord::from_slice(data).expect("record with @ in quality");
        assert_eq!(rec.name(), b"read1");
        assert_eq!(rec.sequence(), b"ACGT");
        assert_eq!(rec.quality(), b"@+!I");
    }

    #[test]
    fn test_fastq_record_from_slice_mismatched_lengths() {
        use crate::fastq_parse::FastqRecord;
        // seq=4 bases but qual=3 chars — should error.
        let data = b"@read1\nACGT\n+\nIII\n";
        let result = FastqRecord::from_slice(data);
        assert!(result.is_err(), "mismatched seq/qual lengths must return an error");
    }

    #[test]
    fn test_fastq_record_from_slice_empty_data() {
        use crate::fastq_parse::FastqRecord;
        let result = FastqRecord::from_slice(b"");
        assert!(result.is_err(), "empty data must return an error");
    }

    #[test]
    fn test_fastq_record_from_slice_no_leading_at() {
        use crate::fastq_parse::FastqRecord;
        let data = b"read1\nACGT\n+\nIIII\n";
        let result = FastqRecord::from_slice(data);
        assert!(result.is_err(), "missing @ prefix must return an error");
    }

    #[test]
    fn test_fastq_record_from_slice_no_trailing_newline() {
        // Quality without trailing newline should still parse (trimmed).
        use crate::fastq_parse::FastqRecord;
        let data = b"@read1\nACGT\n+\nIIII";
        let rec = FastqRecord::from_slice(data).expect("record without trailing newline");
        assert_eq!(rec.quality(), b"IIII");
    }

    // ========================================================================
    // detect_prefix_end Tests
    // ========================================================================

    #[test]
    fn test_detect_prefix_end_empty() {
        // Empty data → 0 (no prefix).
        assert_eq!(detect_prefix_end(b""), 0);
    }

    #[test]
    fn test_detect_prefix_end_starts_on_boundary() {
        // Data begins with '@' → record boundary immediately, prefix_end = 0.
        let data = b"@read1\nACGT\n+\nIIII\n";
        assert_eq!(detect_prefix_end(data), 0);
    }

    #[test]
    fn test_detect_prefix_end_mid_record() {
        // Data starts mid-sequence; the first '@' appears after the partial sequence.
        // prefix is "CGT\n+\nIIII\n", then '@read2\n...' begins.
        let suffix = b"@read2\nGGGG\n+\nJJJJ\n";
        let prefix = b"CGT\n+\nIIII\n";
        let mut data = prefix.to_vec();
        data.extend_from_slice(suffix);
        let end = detect_prefix_end(&data);
        // The returned offset should point to the '@' of read2.
        assert_eq!(&data[end..=end], b"@");
        assert!(end > 0, "prefix_end must be > 0 when data starts mid-record");
    }

    #[test]
    fn test_detect_prefix_end_single_record() {
        // Data is exactly one complete record that starts on a boundary.
        let data = b"@r\nA\n+\nI\n";
        assert_eq!(detect_prefix_end(data), 0);
    }

    #[test]
    fn test_detect_prefix_end_insufficient_data() {
        // Fewer than 4 newlines → cannot determine boundary, return data.len().
        let data = b"partial_no_newlines";
        assert_eq!(detect_prefix_end(data), data.len());
    }

    #[test]
    fn test_detect_prefix_end_at_in_quality_is_not_boundary() {
        // A '@' in quality scores must not be treated as a record boundary.
        // Layout: partial quality "!@\n" then a real record "@r2\nAA\n+\nII\n".
        // The prefix is "!@\n" (3 bytes), so prefix_end should be 3.
        //
        // We must construct data that starts mid-quality (after seq+plus lines).
        // Full record: @r1\nAA\n+\n!@\n => but we start mid-quality at '!'
        // So data = "!@\n" + "@r2\nAA\n+\nII\n"
        let data = b"!@\n@r2\nAA\n+\nII\n";
        let end = detect_prefix_end(data);
        // Should skip the '!@\n' prefix and land on '@r2'.
        assert_eq!(&data[end..=end], b"@");
    }

    // ========================================================================
    // detect_suffix_start Tests
    // ========================================================================

    #[test]
    fn test_detect_suffix_start_empty() {
        assert_eq!(detect_suffix_start(b""), 0);
    }

    #[test]
    fn test_detect_suffix_start_ends_on_boundary() {
        // Data ends exactly on a record boundary (trailing '\n') → data.len().
        let data = b"@read1\nACGT\n+\nIIII\n";
        assert_eq!(detect_suffix_start(data), data.len());
    }

    #[test]
    fn test_detect_suffix_start_ends_mid_record() {
        // Data has one complete record followed by an incomplete one.
        // The complete record is "@r1\nACGT\n+\nIIII\n" (18 bytes).
        // Then "@r2\nGG" is incomplete (only 6 bytes with no quality).
        let complete = b"@r1\nACGT\n+\nIIII\n";
        let partial = b"@r2\nGG";
        let mut data = complete.to_vec();
        data.extend_from_slice(partial);
        let suffix_start = detect_suffix_start(&data);
        // Everything from suffix_start onwards is the incomplete fragment.
        assert_eq!(suffix_start, complete.len());
    }

    #[test]
    fn test_detect_suffix_start_single_record_whole_buffer() {
        // A single complete record → suffix_start = data.len() (no suffix).
        let data = b"@r\nA\n+\nI\n";
        assert_eq!(detect_suffix_start(data), data.len());
    }

    #[test]
    fn test_detect_suffix_start_insufficient_data() {
        // Fewer than 4 newlines → cannot confirm a full record, return 0.
        let data = b"@r1\nAC";
        assert_eq!(detect_suffix_start(data), 0);
    }

    #[test]
    fn test_detect_suffix_start_multiple_records() {
        // Two complete records; suffix_start should be data.len().
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nGGGG\n+\nJJJJ\n";
        assert_eq!(detect_suffix_start(data), data.len());
    }

    // ========================================================================
    // detect_prefix_end / detect_suffix_start Round-trip
    // ========================================================================

    #[test]
    fn test_prefix_suffix_round_trip() {
        // Simulate a BGZF chunk that has a partial record at start and end.
        // Layout: "GT\n+\nIIII\n" (suffix of previous record)
        //       + "@r2\nACGT\n+\nJJJJ\n" (complete middle record)
        //       + "@r3\nAA" (prefix of next record)
        let prefix_bytes = b"GT\n+\nIIII\n";
        let middle_bytes = b"@r2\nACGT\n+\nJJJJ\n";
        let suffix_bytes = b"@r3\nAA";
        let mut data = prefix_bytes.to_vec();
        data.extend_from_slice(middle_bytes);
        data.extend_from_slice(suffix_bytes);

        let prefix_end = detect_prefix_end(&data);
        let suffix_start = detect_suffix_start(&data[prefix_end..]) + prefix_end;

        // prefix_end should skip "GT\n+\nIIII\n" (10 bytes) to land on '@r2'.
        assert_eq!(prefix_end, prefix_bytes.len());
        // suffix_start (relative to start of data) should be where '@r3' begins.
        assert_eq!(suffix_start, prefix_bytes.len() + middle_bytes.len());
        // Middle slice should be exactly middle_bytes.
        assert_eq!(&data[prefix_end..suffix_start], middle_bytes);
    }

    // ========================================================================
    // stitch_cross_block_record Tests
    // ========================================================================

    #[test]
    fn test_stitch_cross_block_record_both_empty() {
        let result = stitch_cross_block_record(b"", b"").expect("no error for empty slices");
        assert!(result.is_none(), "both empty → None");
    }

    #[test]
    fn test_stitch_cross_block_record_valid() {
        // The record "@r1\nACGT\n+\nIIII\n" split across two blocks.
        let suffix = b"@r1\nACGT\n+\n";
        let prefix = b"IIII\n";
        let result = stitch_cross_block_record(suffix, prefix)
            .expect("valid cross-block record")
            .expect("record should be Some");
        assert_eq!(result.name(), b"r1");
        assert_eq!(result.sequence(), b"ACGT");
        assert_eq!(result.quality(), b"IIII");
    }

    // ========================================================================
    // BlockParseFast Integration Tests
    // ========================================================================

    /// Verify that `detect_prefix_end` + `detect_suffix_start` + SIMD offsets correctly
    /// split a decompressed chunk into prefix, middle records, and suffix.
    #[test]
    fn test_block_parse_split_prefix_middle_suffix() {
        // Simulate a decompressed chunk starting mid-record and ending mid-record.
        // Previous block ended with: "@r1\nACGT\n+" (partial)
        // This chunk continues:      "IIII\n" (completes r1)
        //                       then "@r2\nGGGG\n+\nJJJJ\n" (complete)
        //                       then "@r3\nTT" (starts r3, incomplete)
        let prefix_frag = b"\nIIII\n"; // completes r1 (name+seq+plus already in prev block)
        let complete_r2 = b"@r2\nGGGG\n+\nJJJJ\n";
        let suffix_frag = b"@r3\nTT";

        // Build mock "previous suffix" + "current block" data.
        // For the block parse step the data is just the decompressed chunk.
        // We arrange it so that prefix = "\nIIII\n" (6 bytes starting with \n, not @).
        let mut data = prefix_frag.to_vec();
        data.extend_from_slice(complete_r2);
        data.extend_from_slice(suffix_frag);

        let prefix_end = detect_prefix_end(&data);
        let suffix_start = if prefix_end >= data.len() {
            data.len()
        } else {
            detect_suffix_start(&data[prefix_end..]) + prefix_end
        };

        let prefix_bytes = &data[..prefix_end];
        let suffix_bytes = &data[suffix_start..];
        let middle = &data[prefix_end..suffix_start];

        // prefix is the partial fragment before the first complete record.
        assert_eq!(prefix_bytes, prefix_frag);
        // suffix is the incomplete trailing fragment.
        assert_eq!(suffix_bytes, suffix_frag);
        // middle is exactly the complete r2 record.
        assert_eq!(middle, complete_r2);

        // Verify SIMD offsets on middle.
        let offsets = fgumi_simd_fastq::find_record_offsets(middle);
        assert_eq!(offsets.len(), 2, "one complete record → two offsets (start + end)");

        // Parse the middle record.
        let rec =
            FastqRecord::from_slice(&middle[offsets[0]..offsets[1]]).expect("valid r2 record");
        assert_eq!(rec.name(), b"r2");
        assert_eq!(rec.sequence(), b"GGGG");
        assert_eq!(rec.quality(), b"JJJJ");
    }

    #[test]
    fn test_block_parse_starts_on_boundary() {
        // Chunk starts exactly on a record boundary → prefix_end = 0.
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nGGGG\n+\nJJJJ\n";
        let prefix_end = detect_prefix_end(data);
        assert_eq!(prefix_end, 0, "no prefix when data starts on @");

        let suffix_start = detect_suffix_start(&data[prefix_end..]) + prefix_end;
        assert_eq!(suffix_start, data.len(), "no suffix when data ends on record boundary");

        let offsets = fgumi_simd_fastq::find_record_offsets(data);
        assert_eq!(offsets.len(), 3, "two complete records → three offsets");
    }

    #[test]
    fn test_block_parse_entire_buffer_is_prefix() {
        // Chunk has no complete records (all data is a partial record fragment).
        let data = b"only_partial_no_newlines";
        let prefix_end = detect_prefix_end(data);
        // Should return data.len() — entire buffer is prefix.
        assert_eq!(prefix_end, data.len());
    }

    // ========================================================================
    // BlockMerge Integration Tests
    // ========================================================================

    /// Helper: build a `BlockParsed` for stream 0 (R1) with given records.
    fn make_block_parsed(
        block_idx: u64,
        stream_idx: usize,
        records: Vec<FastqRecord>,
        prefix_bytes: Vec<u8>,
        suffix_bytes: Vec<u8>,
    ) -> BlockParsed {
        BlockParsed { block_idx, stream_idx, records, prefix_bytes, suffix_bytes }
    }

    fn make_record(name: &str, seq: &str, qual: &str) -> FastqRecord {
        let raw = format!("@{name}\n{seq}\n+\n{qual}\n");
        FastqRecord::from_slice(raw.as_bytes()).unwrap()
    }

    #[test]
    fn test_block_merge_state_single_stream_basic() {
        // Single-stream (R1 only): 2 records across two blocks, no cross-block record.
        let mut state = BlockMergeState::new();
        assert!(state.is_empty());

        let r1 = make_record("r1", "ACGT", "IIII");
        let r2 = make_record("r2", "GGGG", "JJJJ");

        // Block 0: 1 record, no prefix/suffix.
        let block0 = make_block_parsed(0, 0, vec![r1.clone()], vec![], vec![]);
        // Block 1: 1 record, no prefix/suffix.
        let block1 = make_block_parsed(1, 0, vec![r2.clone()], vec![], vec![]);

        state.r1_pending.insert(0, block0);
        state.r1_pending.insert(1, block1);

        // Process block 0.
        let r1_next = state.r1_next;
        let r1_block = state.r1_pending.remove(&r1_next).unwrap();
        state.r1_next += 1;

        let cross0 =
            stitch_cross_block_record(&state.r1_suffix_bytes, &r1_block.prefix_bytes).unwrap();
        assert!(cross0.is_none(), "no cross-block record in block 0");

        let mut all_r1: Vec<FastqRecord> = std::mem::take(&mut state.r1_surplus);
        if let Some(rec) = cross0 {
            all_r1.push(rec);
        }
        all_r1.extend(r1_block.records);
        state.r1_suffix_bytes = r1_block.suffix_bytes;

        assert_eq!(all_r1.len(), 1);
        assert_eq!(all_r1[0].name(), b"r1");
        assert_eq!(all_r1[0].sequence(), b"ACGT");
    }

    #[test]
    fn test_block_merge_cross_block_record() {
        // The record "@r1\nACGT\n+\nIIII\n" is split: "@r1\nACGT\n+\n" in block 0's suffix,
        // "IIII\n" in block 1's prefix.
        let suffix_bytes = b"@r1\nACGT\n+\n".to_vec();
        let prefix_bytes = b"IIII\n".to_vec();

        let cross = stitch_cross_block_record(&suffix_bytes, &prefix_bytes)
            .expect("valid cross-block record")
            .expect("should be Some");

        assert_eq!(cross.name(), b"r1");
        assert_eq!(cross.sequence(), b"ACGT");
        assert_eq!(cross.quality(), b"IIII");
    }

    #[test]
    fn test_block_merge_r1_surplus_carries() {
        // R1 block has 3 records, R2 block has 2 records → 1 R1 surplus.
        let r1_records = [
            make_record("r1", "AAAA", "IIII"),
            make_record("r2", "CCCC", "IIII"),
            make_record("r3", "GGGG", "IIII"),
        ];
        let r2_records = [make_record("r1", "TTTT", "JJJJ"), make_record("r2", "CCCC", "JJJJ")];

        let pair_count = r1_records.len().min(r2_records.len()); // 2
        let surplus: Vec<FastqRecord> = r1_records[pair_count..].to_vec();

        assert_eq!(pair_count, 2);
        assert_eq!(surplus.len(), 1);
        assert_eq!(surplus[0].name(), b"r3");
    }

    #[test]
    fn test_block_merge_r2_surplus_carries() {
        // R2 block has more records than R1 block → R2 surplus.
        let r1_records = [make_record("r1", "AAAA", "IIII")];
        let r2_records = [make_record("r1", "TTTT", "JJJJ"), make_record("r2", "CCCC", "JJJJ")];

        let pair_count = r1_records.len().min(r2_records.len()); // 1
        let r2_surplus: Vec<FastqRecord> = r2_records[pair_count..].to_vec();

        assert_eq!(pair_count, 1);
        assert_eq!(r2_surplus.len(), 1);
        assert_eq!(r2_surplus[0].name(), b"r2");
    }

    #[test]
    fn test_block_merge_state_out_of_order_insertion() {
        // Insert block 1 before block 0; verify r1_next ordering.
        let mut state = BlockMergeState::new();
        let r1 = make_record("r1", "ACGT", "IIII");
        let r2 = make_record("r2", "GGGG", "JJJJ");

        state.r1_pending.insert(1, make_block_parsed(1, 0, vec![r2], vec![], vec![]));
        state.r1_pending.insert(0, make_block_parsed(0, 0, vec![r1], vec![], vec![]));

        // Block 0 must be processed first (BTreeMap ordering by key).
        let first_key = *state.r1_pending.keys().next().unwrap();
        assert_eq!(first_key, 0, "BTreeMap must yield block 0 first");
    }

    #[test]
    fn test_block_merge_state_empty() {
        let state = BlockMergeState::new();
        assert!(state.is_empty());
        assert!(state.r1_pending.is_empty());
        assert!(state.r2_pending.is_empty());
        assert!(state.r1_surplus.is_empty());
        assert!(state.r2_surplus.is_empty());
        assert!(state.r1_suffix_bytes.is_empty());
        assert!(state.r2_suffix_bytes.is_empty());
        assert_eq!(state.serial_out, 0);
    }

    /// Helper to create a minimal FASTQ pipeline state for backpressure tests.
    fn make_fastq_state() -> FastqPipelineState<Cursor<Vec<u8>>, Vec<u8>> {
        let config = FastqPipelineConfig::new(2, false, 6);
        let readers = vec![StreamReader::Decompressed(Cursor::new(Vec::new()))];
        let output: Box<dyn std::io::Write + Send> = Box::new(Vec::<u8>::new());
        FastqPipelineState::new(config, readers, output)
    }

    /// Verify that the FASTQ pipeline's scheduler backpressure reports memory
    /// state from the write reorder buffer, matching the BAM pipeline's behavior.
    ///
    /// When the write reorder buffer is above its memory threshold, the scheduler
    /// must see `memory_high: true` so it can redirect threads to drain Compress→Write.
    /// Without this, threads pile up on upstream steps while nobody does output work.
    #[test]
    fn test_fastq_backpressure_reports_write_reorder_memory() {
        let state = make_fastq_state();
        let header = Header::default();
        let process_fn = |_t: FastqTemplate| -> io::Result<Vec<u8>> { Ok(Vec::new()) };
        let serialize_fn =
            |_p: Vec<u8>, _h: &Header, _buf: &mut Vec<u8>| -> io::Result<u64> { Ok(0) };

        let ctx = FastqStepContext {
            state: &state,
            header: &header,
            process_fn: &process_fn,
            serialize_fn: &serialize_fn,
            is_reader: false,
        };

        let active = ActiveSteps::all();
        let worker = FastqWorkerState::new(6, 0, 2, SchedulerStrategy::default(), active);

        // Initially: memory should not be high
        let bp = ctx.get_backpressure(&worker);
        assert!(!bp.memory_high, "memory_high should be false when write reorder buffer is empty");
        assert!(
            bp.memory_drained,
            "memory_drained should be true when write reorder buffer is empty"
        );

        // Push write reorder buffer above its threshold (512 MB default)
        let threshold = BACKPRESSURE_THRESHOLD_BYTES;
        state.output.write_reorder_state.add_heap_bytes(threshold + 1);

        let bp = ctx.get_backpressure(&worker);
        assert!(
            bp.memory_high,
            "memory_high should be true when write reorder buffer exceeds limit"
        );
        assert!(
            !bp.memory_drained,
            "memory_drained should be false when write reorder buffer exceeds limit"
        );

        // When draining, memory_high should be suppressed (same as BAM pipeline)
        state.output.set_draining(true);
        let bp = ctx.get_backpressure(&worker);
        assert!(!bp.memory_high, "memory_high should be false during drain even with high memory");

        // Drain memory below low-water mark (50% of threshold)
        state.output.write_reorder_state.sub_heap_bytes(threshold + 1);
        state.output.set_draining(false);
        let bp = ctx.get_backpressure(&worker);
        assert!(!bp.memory_high, "memory_high should be false after draining memory");
        assert!(
            bp.memory_drained,
            "memory_drained should be true after draining below low-water mark"
        );
    }
}
