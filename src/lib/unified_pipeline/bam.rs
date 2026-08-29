//! BAM pipeline: Read → Decompress → `FindBoundaries` → Decode → Group → Process → Serialize → Compress → Write.
//!
//! This module implements the BAM-specific pipeline for processing BAM files
//! with grouping operations like `group` and `codec`.

use crossbeam_queue::ArrayQueue;
use noodles::bam::{self};
use noodles::sam::{Header, alignment::RecordBuf};
use parking_lot::Mutex;
use std::collections::VecDeque;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::thread;
use std::time::{Duration, Instant};

use crate::bgzf_reader::{BGZF_EOF, decompress_block_into_opts, read_raw_blocks};
use crate::bgzf_writer::InlineBgzfCompressor;
#[cfg(test)]
use crate::sam::SamTag;
use fgumi_bam_io::ProgressTracker;
use fgumi_bam_io::ReorderBuffer;

use super::base::{
    ActiveSteps, BatchWeight, CompressedBlockBatch, DecodedRecord, DecompressedBatch,
    GroupKeyConfig, HasCompressor, HasHeldBoundaries, HasHeldCompressed, HasHeldProcessed,
    HasHeldSerialized, HasRecycledBuffers, HasWorkerCore, MemoryEstimate, MonitorableState,
    OutputPipelineQueues, OutputPipelineState, PROGRESS_LOG_INTERVAL, PipelineConfig,
    PipelineLifecycle, PipelineStats, PipelineStep, PipelineValidationError, ProcessPipelineState,
    QueueSample, RawBlockBatch, ReorderBufferState, SerializePipelineState, SerializedBatch,
    StepContext, WorkerCoreState, WorkerStateCommon, WritePipelineState,
    finalize_pipeline_with_buffers, generic_worker_loop, handle_worker_panic, join_monitor_thread,
    join_worker_threads, push_charged, refund_queue_bytes, shared_try_step_compress,
};
use super::deadlock::{
    DeadlockAction, DeadlockConfig, DeadlockState, QueueSnapshot, check_deadlock_and_restore,
};
use super::scheduler::{BackpressureState, SchedulerStrategy};
use crate::read_info::LibraryIndex;
use fgumi_raw_bam;
use fgumi_raw_bam::RawRecordView;

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

impl MemoryEstimate for BoundaryBatch {
    fn estimate_heap_size(&self) -> usize {
        self.buffer.capacity() + self.offsets.capacity() * std::mem::size_of::<usize>()
    }
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

    /// Number of bytes currently held from an incomplete record.
    ///
    /// Non-zero at pipeline completion means input was read but never became a
    /// record, so completion validation reports it.
    #[must_use]
    pub fn leftover_len(&self) -> usize {
        self.leftover.len()
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
        if &data[0..4] != fgumi_raw_bam::BAM_MAGIC {
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
    ///
    /// # Errors
    ///
    /// Returns an I/O error if the BAM header is malformed.
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
    ///
    /// # Errors
    ///
    /// Returns an I/O error if there are incomplete BAM records at EOF.
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

        // Fewer than four bytes left: not even a `block_size`, so the record's
        // own length is unknown. The scan loop above cannot see this — it stops
        // at `cursor + 4 > len` — and the tail used to be discarded, silently
        // with `Ok(None)` when no complete record preceded it and by sitting
        // past the final offset when one did.
        if cursor < self.leftover.len() {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "Incomplete BAM record at EOF: {} trailing byte(s) are too few to hold a \
                     record length",
                    self.leftover.len() - cursor
                ),
            ));
        }

        debug_assert!(cursor > 0, "cursor == 0 is caught by the trailing-byte check above");

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
/// * `group_key_config` - `Some` config for computing each record's `GroupKey`
///   (library index, cell tag, UMI tag), or `None` to skip key computation
///   entirely and attach a default placeholder key. `None` is only for groupers
///   that never read the key (e.g. `SingleRawRecordGrouper`, which discards the
///   decoded record). Pairing `None` with a key-consuming grouper mis-groups:
///   key-only groupers (position/MI) collapse every record into one group, and
///   `TemplateGrouper` loses its `name_hash` fast-path (its QNAME byte-compare
///   still keeps it correct, only slower). Every key-consuming command sets
///   `Some` explicitly, so `None` is never reached by accident.
///
/// # Returns
///
/// A vector of decoded `DecodedRecord` instances (record + `GroupKey`).
///
/// # Errors
///
/// Returns an I/O error if any BAM record is malformed.
pub fn decode_records(
    batch: &BoundaryBatch,
    group_key_config: Option<&GroupKeyConfig>,
) -> io::Result<Vec<DecodedRecord>> {
    let num_records = batch.offsets.len().saturating_sub(1);
    let mut records = Vec::with_capacity(num_records);

    let buffer_len = batch.buffer.len();
    for i in 0..num_records {
        let start = batch.offsets[i];
        let end = batch.offsets[i + 1];

        // Validate: ordering, in-bounds, and room for the block_size prefix.
        // Guards the raw indexes into batch.buffer below against malformed input.
        if start >= end || end > buffer_len || end - start <= 4 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Invalid record bounds: start={start}, end={end}, record_index={i}, \
                     num_records={num_records}, buffer_len={buffer_len}"
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
                     record_index={i}, start={start}, end={end}, buffer_len={buffer_len}"
                ),
            ));
        }

        // Skip the 4-byte block_size prefix and copy record data.
        let raw = batch.buffer[start + 4..end].to_vec();
        if raw.len() < 32 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("BAM record too short: len={}", raw.len()),
            ));
        }
        // Validate l_read_name fits within the record to prevent panics in
        // read_name(). CIGAR, sequence, and aux access still rely on the raw
        // bytes being a well-formed BAM record produced by the BAM reader;
        // callers must not feed arbitrary external bytes here.
        let l_rn = raw[8] as usize;
        if raw.len() < 32 + l_rn {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "BAM record truncated: len={}, l_read_name={l_rn} (need >= {})",
                    raw.len(),
                    32 + l_rn
                ),
            ));
        }
        // Skip the per-record key computation when no grouper needs it; the
        // placeholder key is never read in that case (see the `None` contract).
        let (key, umi_position) = match group_key_config {
            Some(cfg) => {
                compute_group_key_from_raw(&raw, &cfg.library_index, cfg.cell_tag, cfg.umi_tag)
            }
            None => (super::base::GroupKey::default(), None),
        };
        let mut decoded = DecodedRecord::from_raw_bytes(raw, key);
        // The UMI position returned above is record-relative (already converted
        // from aux-relative inside `compute_group_key_from_raw`).
        if let Some((value_offset, value_len)) = umi_position {
            decoded.set_cached_umi(value_offset, value_len);
        }
        records.push(decoded);
    }

    Ok(records)
}

/// Compute a `GroupKey` directly from raw BAM bytes, matching `compute_group_key()` exactly.
///
/// Uses 1-based coordinate helpers to produce identical keys to the noodles path.
///
/// When `umi_tag` is `Some`, the same aux-data scan that locates RG, the cell
/// barcode, and MC also locates the UMI tag's value. The returned tuple's
/// second element is the **record-relative** `(value_offset, value_len)` of
/// the UMI tag's value bytes (excluding NUL terminator), suitable for
/// `DecodedRecord::set_cached_umi`. It is `None` when no `umi_tag` was
/// supplied, when the tag is absent or non-Z-typed, or when the record's aux
/// offset cannot be resolved.
///
/// # Panics
///
/// Panics if `raw` is not a validated BAM record payload (as produced by the BAM
/// reader / raw-record pipeline). Callers must not pass arbitrary external bytes;
/// raw-field accessors will panic on malformed or truncated input.
#[must_use]
pub fn compute_group_key_from_raw(
    raw: &[u8],
    library_index: &LibraryIndex,
    cell_tag: Option<noodles::sam::alignment::record::data::field::Tag>,
    umi_tag: Option<[u8; 2]>,
) -> (super::base::GroupKey, Option<(u32, u16)>) {
    use super::base::GroupKey;

    // Extract name hash (match noodles path: empty name → None → hash 0)
    let name = fgumi_raw_bam::read_name(raw);
    let name_hash = if name.is_empty() {
        LibraryIndex::hash_name(None)
    } else {
        LibraryIndex::hash_name(Some(name))
    };

    // Check secondary/supplementary
    let flg = RawRecordView::new(raw).flags();
    let is_secondary = (flg & fgumi_raw_bam::flags::SECONDARY) != 0;
    let is_supplementary = (flg & fgumi_raw_bam::flags::SUPPLEMENTARY) != 0;
    if is_secondary || is_supplementary {
        // A secondary/supplementary read cannot compute its own template
        // coordinate (it lacks its own primary's position). When `fgumi zipper`
        // has stamped the exact coordinate into `tc`, key on it so the read
        // groups into the same position group as its primary — instead of an
        // UNKNOWN position that relies on the read sorting adjacent to its
        // primary. Library/cell are extracted the same way as primaries so the
        // position key matches (they share their template's RG/CB).
        let aux_offset = fgumi_raw_bam::aux_data_offset_from_record(raw).unwrap_or(raw.len());
        let aux_data = &raw[aux_offset..];
        if let Some([tid1, pos1, neg1, tid2, pos2, neg2]) =
            fgumi_raw_bam::read_tc_template_coordinate(aux_data)
        {
            let cell_tag_bytes = cell_tag.map_or([0u8; 2], |t| [t.as_ref()[0], t.as_ref()[1]]);
            let aux_tags =
                fgumi_raw_bam::extract_aux_string_tags(aux_data, cell_tag_bytes, umi_tag);
            let library_idx =
                aux_tags.rg.map_or(0, |rg| library_index.get(LibraryIndex::hash_rg(rg)));
            let cell_hash = aux_tags.cell.map_or(0, |cb| LibraryIndex::hash_cell_barcode(Some(cb)));
            let key = GroupKey::paired(
                tid1,
                pos1,
                u8::from(neg1 != 0),
                tid2,
                pos2,
                u8::from(neg2 != 0),
                library_idx,
                cell_hash,
                name_hash,
            );
            return (key, None);
        }
        return (GroupKey { name_hash, ..GroupKey::default() }, None);
    }

    // Own position (1-based, matching noodles) — zero-allocation CIGAR iteration
    let reverse = (flg & fgumi_raw_bam::flags::REVERSE) != 0;
    let own_pos = fgumi_raw_bam::unclipped_5prime_from_raw_bam(raw);

    // A mapped record with an empty or truncated CIGAR has no computable
    // unclipped 5' position, so `unclipped_5prime_from_raw_bam` returns the
    // `i32::MAX` sentinel. Fall back to the name-only key rather than letting the
    // sentinel reach a position slot — otherwise distinct templates sharing
    // ref/strand/library/cell would collide on `i32::MAX` (`position_key`
    // excludes `name_hash`). Matches the secondary/supplementary fallback above.
    if own_pos == i32::MAX {
        return (GroupKey { name_hash, ..GroupKey::default() }, None);
    }

    let own_ref_id = fgumi_raw_bam::ref_id(raw);
    let strand = u8::from(reverse);

    // Single-pass aux tag extraction (RG, cell barcode, MC, optional UMI).
    // Resolve aux_offset once so we can both slice aux data and convert the
    // returned UMI position from aux-relative to record-relative.
    let aux_offset = fgumi_raw_bam::aux_data_offset_from_record(raw).unwrap_or(raw.len());
    let aux_data = &raw[aux_offset..];
    let cell_tag_bytes = cell_tag.map_or([0u8; 2], |t| [t.as_ref()[0], t.as_ref()[1]]);
    let aux_tags = fgumi_raw_bam::extract_aux_string_tags(aux_data, cell_tag_bytes, umi_tag);

    // Convert UMI position from aux-relative to record-relative.
    let umi_position = aux_tags.umi_position.and_then(|(value_off_in_aux, value_len)| {
        let aux_offset_u32 = u32::try_from(aux_offset).ok()?;
        let value_offset = aux_offset_u32.checked_add(value_off_in_aux)?;
        Some((value_offset, value_len))
    });

    let library_idx = if let Some(rg) = aux_tags.rg {
        let rg_hash = LibraryIndex::hash_rg(rg);
        library_index.get(rg_hash)
    } else {
        0
    };

    let cell_hash =
        if let Some(cb) = aux_tags.cell { LibraryIndex::hash_cell_barcode(Some(cb)) } else { 0 };

    // Check if paired
    let is_paired = (flg & fgumi_raw_bam::flags::PAIRED) != 0;
    if !is_paired {
        return (
            GroupKey::single(own_ref_id, own_pos, strand, library_idx, cell_hash, name_hash),
            umi_position,
        );
    }

    // Mate info — guard against MATE_UNMAPPED (matching noodles path)
    let mate_unmapped = (flg & fgumi_raw_bam::flags::MATE_UNMAPPED) != 0;
    let mate_reverse = (flg & fgumi_raw_bam::flags::MATE_REVERSE) != 0;
    let mate_strand = u8::from(mate_reverse);
    let raw_mate_ref_id = fgumi_raw_bam::mate_ref_id(raw);
    let raw_mate_pos = fgumi_raw_bam::mate_pos(raw);

    // Get mate unclipped 5' position via MC tag (skip if mate is unmapped)
    let mate_pos_result = if mate_unmapped {
        None
    } else {
        aux_tags
            .mc
            .map(|mc| fgumi_raw_bam::mate_unclipped_5prime_1based(raw_mate_pos, mate_reverse, mc))
    };

    let key = match mate_pos_result {
        Some(mp) => GroupKey::paired(
            own_ref_id,
            own_pos,
            strand,
            raw_mate_ref_id,
            mp,
            mate_strand,
            library_idx,
            cell_hash,
            name_hash,
        ),
        None => {
            // No MC tag — fall back to single-end behavior.
            //
            // `RecordPositionGrouper::validate_mc_tag` detects a missing MC tag by
            // testing `GroupKey::has_mate_position()` on the key produced here,
            // instead of re-walking the aux data on its serial step. Keep this the
            // only way an eligible paired record can receive a single-ended key.
            GroupKey::single(own_ref_id, own_pos, strand, library_idx, cell_hash, name_hash)
        }
    };
    (key, umi_position)
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
    ///
    /// Private: push and pop through the `q1_push` / `q1_pop` helpers so
    /// `q1_heap_bytes` stays synchronized. External mutation would let a batch
    /// enter or leave the queue without its paired byte accounting.
    q1_raw_blocks: ArrayQueue<(u64, RawBlockBatch)>,
    /// Heap bytes currently held in Q1 (compressed BGZF blocks).
    ///
    /// Charged against the queue memory budget by
    /// [`BamPipelineState::queue_bytes_in_flight`], so it is maintained
    /// unconditionally rather than only under `memory-debug`. The charge is
    /// `RawBlockBatch::estimate_heap_size` (allocation capacity), not
    /// `total_compressed_size` (logical length) — the budget has to account for
    /// the memory actually held, and `read_raw_blocks` over-allocates the block
    /// `Vec` to `blocks_per_read_batch` whether or not it fills it.
    ///
    /// Private: maintained only by the `q1_push` / `q1_pop` helpers.
    q1_heap_bytes: AtomicU64,

    // ========== Step 2: Decompress (parallel) ==========
    // No state needed - each thread has its own Decompressor
    //
    // There is no `decompress_done` flag: downstream stages track this step through
    // `batches_decompressed`, and `FindBoundaries` infers it from
    // `batches_boundary_processed == next_read_serial` (nothing can be processed that
    // was not decompressed first).
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
    ///
    /// Push and pop through the private `q2b_push` / `q2b_pop` helpers so
    /// `q2b_heap_bytes` stays accurate. Private for the same reason.
    q2b_boundaries: ArrayQueue<(u64, BoundaryBatch)>,
    /// Heap bytes currently held in Q2b (decompressed record buffers).
    ///
    /// Q2b holds whole decompressed BGZF payloads, so it is one of the larger
    /// consumers when the pipeline backs up; it is charged against the queue
    /// memory budget by [`BamPipelineState::queue_bytes_in_flight`].
    ///
    /// Private: maintained only by the `q2b_push` / `q2b_pop` helpers.
    q2b_heap_bytes: AtomicU64,

    // ========== Step 4: Decode (parallel) ==========
    //
    // There is no `decode_done` flag: `Group` tracks this step through
    // `batches_grouped == batches_boundary_found`, which counts batches only once they
    // have been pushed, and so cannot conclude early while one is still held for a
    // retried push.
    /// Count of batches that have completed decoding (for completion tracking).
    pub batches_decoded: AtomicU64,
    /// Configuration for computing `GroupKey` during decode, or `None` to skip
    /// key computation for groupers that never read it (see [`decode_records`]).
    pub group_key_config: Option<GroupKeyConfig>,

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
        group_key_config: Option<GroupKeyConfig>,
    ) -> Self {
        let cap = config.queue_capacity;
        // Q3 decode admission window: bound reorder *skew* to concurrency (a few
        // batches per worker), capped by queue capacity. This is not a byte
        // limit — total memory is bounded upstream by the Read admission gate —
        // see `ReorderBufferState::can_proceed`. Restores decode parallelism
        // lost to fgumi#787's per-serial byte backpressure.
        let q3_decode_window = config.num_threads.saturating_mul(4).min(cap).max(1);
        let memory_limit = config.queue_memory_limit;
        let stats = if config.collect_stats {
            config.shared_stats.clone().or_else(|| Some(Arc::new(PipelineStats::new())))
        } else {
            None
        };
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
            q1_heap_bytes: AtomicU64::new(0),
            // Step 2: Decompress
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
            q2b_heap_bytes: AtomicU64::new(0),
            // Step 4: Decode
            batches_decoded: AtomicU64::new(0),
            group_key_config,
            // Q3: Decode → Group (with reorder)
            q3_decoded: ArrayQueue::new(cap),
            q3_reorder: Mutex::new(ReorderBuffer::new()),
            // Q3 reorder buffer atomic state (for lock-free admission control)
            q3_reorder_state: ReorderBufferState::new(memory_limit)
                .with_window(q3_decode_window as u64),
            q3_reorder_can_pop: AtomicBool::new(false),
            // Step 5: Group
            group_done: AtomicBool::new(false),
            next_group_serial: AtomicU64::new(0),
            batches_grouped: AtomicU64::new(0),
            // Output-half state (Group → Process → Serialize → Compress → Write)
            output: OutputPipelineQueues::new(
                cap,
                output,
                stats,
                "Processed records",
                memory_limit,
            ),
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
    /// Bounds the Q3 reorder buffer by *serial skew* -- how far ahead of the
    /// serial Group needs (`next_seq`) decode may run -- rather than by a byte
    /// threshold. Total in-flight memory is bounded upstream by the Read
    /// admission gate; see [`super::base::ReorderBufferState::can_proceed`].
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
    /// to the scheduler. See [`super::base::BACKPRESSURE_THRESHOLD_BYTES`] for architecture details.
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

    /// Push a boundary batch onto Q2b, charging its heap bytes on success.
    ///
    /// Returns the batch unchanged when Q2b is out of slots, so callers keep
    /// the existing held-item retry behaviour.
    fn q2b_push(&self, serial: u64, batch: BoundaryBatch) -> Result<(), (u64, BoundaryBatch)> {
        // `push_charged` charges before the push so the batch is never visible to
        // a consumer while its bytes are still uncharged: a pop in that window
        // would refund bytes never added and, because `q2b_heap_bytes` gates Read
        // via `queue_bytes_in_flight`, permanently over-state the counter and
        // close the Read gate. It refunds on the out-of-slots path so a rejected
        // push leaves the counter unchanged.
        let heap_size = batch.estimate_heap_size() as u64;
        push_charged(&self.q2b_boundaries, &self.q2b_heap_bytes, heap_size, (serial, batch))
    }

    /// Pop a boundary batch from Q2b, refunding its heap bytes.
    fn q2b_pop(&self) -> Option<(u64, BoundaryBatch)> {
        let (serial, batch) = self.q2b_boundaries.pop()?;
        refund_queue_bytes(&self.q2b_heap_bytes, batch.estimate_heap_size() as u64);
        Some((serial, batch))
    }

    /// Decide whether a boundary batch just popped from Q2b should be decoded
    /// now or handed back to Q2b under per-serial memory backpressure (fgumi#746).
    ///
    /// The Q3 reorder buffer's own [`ReorderBufferState::can_proceed`] always
    /// admits the serial Group needs next (the gap-filler), even when memory is
    /// high, and backpressures only *future* serials that fall outside the
    /// serial-skew window. Returns:
    ///
    /// - `Some((serial, batch))` — decode it now: either it is the gap-filler,
    ///   or it is a future batch whose requeue lost the race to a full Q2b (a
    ///   rare, bounded overshoot — never a lost or reordered-into-loss batch,
    ///   since Q3 reassembles by serial). Decoding it here skips the per-serial
    ///   Q3 memory backstop, but only as defense-in-depth: the Read admission
    ///   gate enforces the same `queue_memory_limit` on `queue_bytes_in_flight`,
    ///   a sum that already includes the Q3 and Q2b heap, so the overshoot
    ///   cannot push total in-flight memory past the global bound.
    /// - `None` — the batch was requeued to Q2b so another worker can still
    ///   reach the gap-filler, avoiding the "all workers hold a non-`next_seq`
    ///   batch and nobody can produce `next_seq`" deadlock.
    fn admit_or_requeue_decode(
        &self,
        serial: u64,
        batch: BoundaryBatch,
    ) -> Option<(u64, BoundaryBatch)> {
        if self.q3_reorder_state.can_proceed(serial) {
            return Some((serial, batch));
        }
        // On requeue-success (`Ok`) return `None` so the caller advances no work;
        // on requeue-failure the batch is handed back to be decoded directly
        // (bounded overshoot). `Result::err` maps exactly that: `Ok` -> `None`,
        // `Err(returned)` -> `Some(returned)`.
        self.q2b_push(serial, batch).err()
    }

    /// Push a raw block batch onto Q1, charging its heap bytes on success.
    ///
    /// Mirrors [`Self::q2b_push`]: the charge precedes the push (see that method
    /// for why a visible-but-uncharged batch permanently over-states the Read
    /// gate), and the out-of-slots path refunds so a rejected push leaves the
    /// counter unchanged. Callers record the deadlock-detector push on `Ok`.
    fn q1_push(&self, serial: u64, batch: RawBlockBatch) -> Result<(), (u64, RawBlockBatch)> {
        let heap_size = batch.estimate_heap_size() as u64;
        push_charged(&self.q1_raw_blocks, &self.q1_heap_bytes, heap_size, (serial, batch))
    }

    /// Pop a raw block batch from Q1, refunding its heap bytes.
    fn q1_pop(&self) -> Option<(u64, RawBlockBatch)> {
        let (serial, batch) = self.q1_raw_blocks.pop()?;
        refund_queue_bytes(&self.q1_heap_bytes, batch.estimate_heap_size() as u64);
        Some((serial, batch))
    }

    /// Heap bytes the pipeline is currently holding in its accounted queues.
    ///
    /// Sums every byte counter the pipeline maintains: Q1 (compressed blocks),
    /// Q2 and Q3 — whose `ReorderBufferState` counters span both the queue and
    /// its reorder buffer — Q2b, Q4, Q5 (which also covers the MI-assign
    /// reorder buffer), Q6, Q7, and the write reorder buffer.
    ///
    /// This is an estimate, not an exact figure: it counts queued batches, not
    /// the per-thread working memory a stage allocates while operating on one,
    /// and it uses `Vec::capacity` rather than the allocator's true block
    /// sizes. It is nonetheless the whole of the data the pipeline parks
    /// between stages, which is what grows without bound when the writer
    /// stalls.
    #[must_use]
    pub fn queue_bytes_in_flight(&self) -> u64 {
        // Saturating sum: every debit that feeds these counters now floors at
        // zero, but the nine loads are not one atomic snapshot, so an in-flight
        // debit/credit interleaving could still transiently read a counter high.
        // Fold with `saturating_add` so the aggregate can never overflow the
        // checked `+` in a debug/coverage build even under that interleaving —
        // the gate is conservative for an instant rather than panicking.
        [
            self.q1_heap_bytes.load(Ordering::Acquire),
            self.q2_reorder_state.get_heap_bytes(),
            self.q2b_heap_bytes.load(Ordering::Acquire),
            self.q3_reorder_state.get_heap_bytes(),
            self.output.groups_heap_bytes.load(Ordering::Acquire),
            self.output.processed_heap_bytes.load(Ordering::Acquire),
            self.output.serialized_heap_bytes.load(Ordering::Acquire),
            self.output.compressed_heap_bytes.load(Ordering::Acquire),
            self.output.write_reorder_state.get_heap_bytes(),
        ]
        .into_iter()
        .fold(0u64, u64::saturating_add)
    }

    /// Whether the Read step may admit another batch under the queue memory
    /// budget.
    ///
    /// Read is the pipeline's only source of new bytes, and every other stage
    /// is bounded by a slot count rather than by bytes — so when the output
    /// device stalls, each stage fills with however many bytes its slots
    /// happen to hold and the configured budget bounds nothing. Gating Read on
    /// [`Self::queue_bytes_in_flight`] is what turns that budget into a real
    /// ceiling (issue #746).
    ///
    /// This cannot deadlock. Nothing downstream waits on Read to make
    /// progress: `read_done` is set only at EOF, so declining to read leaves
    /// every stage free to drain, which lowers the in-flight total and lets
    /// reading resume. Admission is also always granted when nothing is
    /// accounted for in flight, so a single batch larger than the whole budget
    /// still gets through instead of stalling the pipeline forever.
    ///
    /// A `queue_memory_limit` of 0 means "no limit" and disables the gate.
    #[must_use]
    pub fn read_admission_allowed(&self) -> bool {
        let limit = self.config.queue_memory_limit;
        if limit == 0 {
            return true;
        }
        let in_flight = self.queue_bytes_in_flight();
        // Test-only observation hook (`None` in production — a single predictable
        // branch). Publishing the in-flight total from the very check that gates
        // the reader lets backpressure tests wait on the *actual* park condition
        // (`in_flight >= limit`) instead of inferring it from wall-clock
        // quiescence, which false-settles when the reader is descheduled (#809).
        if let Some(probe) = &self.config.queue_bytes_probe {
            probe.store(in_flight, Ordering::Relaxed);
        }
        in_flight == 0 || in_flight < limit
    }

    /// Check if the pipeline is in drain mode (input exhausted, completing remaining work).
    #[must_use]
    pub fn is_draining(&self) -> bool {
        self.output.is_draining()
    }

    /// Get optional reference to pipeline statistics.
    #[must_use]
    pub fn stats(&self) -> Option<&PipelineStats> {
        self.output.stats.as_deref()
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

    /// Flush the output writer, write the BGZF EOF marker, and finalize.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if writing the BGZF EOF or flushing fails.
    pub fn flush_output(&self) -> io::Result<()> {
        if let Some(mut writer) = self.output.output.lock().take() {
            writer.flush()?;
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
    /// 3. The boundary finder holds no partial record
    ///
    /// The Group step's buffers are *not* visible here — the grouper and its
    /// pending groups sit behind their own mutex — so
    /// [`super::base::finalize_pipeline_with_buffers`] folds
    /// [`GroupState::unflushed_buffers`] into this result.
    ///
    /// Note: Heap byte tracking is reported but advisory only (set to 0) because
    /// estimation can be imprecise. Only queue emptiness and counter checks
    /// cause validation failure.
    ///
    /// # Errors
    ///
    /// Returns `PipelineValidationError` with diagnostics if any issues are detected.
    pub fn validate_completion(&self) -> Result<(), PipelineValidationError> {
        let mut non_empty_queues = Vec::new();
        let mut counter_mismatches = Vec::new();

        // Check all input-half queues are empty
        if !self.q1_raw_blocks.is_empty() {
            non_empty_queues.push(format!("q1_raw_blocks ({})", self.q1_raw_blocks.len()));
        }
        if !self.q2_decompressed.is_empty() {
            non_empty_queues.push(format!("q2_decompressed ({})", self.q2_decompressed.len()));
        }
        if !self.q2b_boundaries.is_empty() {
            non_empty_queues.push(format!("q2b_boundaries ({})", self.q2b_boundaries.len()));
        }
        if !self.q3_decoded.is_empty() {
            non_empty_queues.push(format!("q3_decoded ({})", self.q3_decoded.len()));
        }

        // Check output-half queues are empty
        if !self.output.groups.is_empty() {
            non_empty_queues.push(format!("q4_groups ({})", self.output.groups.len()));
        }
        if !self.output.processed.is_empty() {
            non_empty_queues.push(format!("q5_processed ({})", self.output.processed.len()));
        }
        if !self.output.serialized.is_empty() {
            non_empty_queues.push(format!("q6_serialized ({})", self.output.serialized.len()));
        }
        if !self.output.compressed.is_empty() {
            non_empty_queues.push(format!("q7_compressed ({})", self.output.compressed.len()));
        }

        // Check reorder buffers are empty
        {
            let q2_reorder = self.q2_reorder.lock();
            if !q2_reorder.is_empty() {
                non_empty_queues.push(format!("q2_reorder ({})", q2_reorder.len()));
            }
        }
        {
            let q3_reorder = self.q3_reorder.lock();
            if !q3_reorder.is_empty() {
                non_empty_queues.push(format!("q3_reorder ({})", q3_reorder.len()));
            }
        }
        {
            let mi_assign_reorder = self.output.mi_assign_reorder.lock();
            if !mi_assign_reorder.is_empty() {
                non_empty_queues.push(format!("mi_assign_reorder ({})", mi_assign_reorder.len()));
            }
        }
        {
            let write_reorder = self.output.write_reorder.lock();
            if !write_reorder.is_empty() {
                non_empty_queues.push(format!("write_reorder ({})", write_reorder.len()));
            }
        }

        // Check the boundary finder holds no partial record. Leftover bytes at
        // completion are input that never became a record, which is exactly the
        // internal-buffer case `PipelineLifecycle::validate_completion`
        // documents (and which the FASTQ side already reports per stream).
        {
            let leftover_len = self.boundary_state.lock().leftover_len();
            if leftover_len > 0 {
                non_empty_queues.push(format!("boundary_leftover ({leftover_len})"));
            }
        }

        // Check batch counter invariants
        let total_read = self.next_read_serial.load(Ordering::Acquire);
        let batches_decompressed = self.batches_decompressed.load(Ordering::Acquire);
        let batches_boundary_processed = self.batches_boundary_processed.load(Ordering::Acquire);
        let batches_boundary_found = self.batches_boundary_found.load(Ordering::Acquire);
        let batches_decoded = self.batches_decoded.load(Ordering::Acquire);
        let batches_grouped = self.batches_grouped.load(Ordering::Acquire);

        // Each batch flows through: Read -> Decompress -> FindBoundaries -> Decode -> Group
        if batches_decompressed != total_read {
            counter_mismatches.push(format!(
                "batches_decompressed ({batches_decompressed}) != total_read ({total_read})"
            ));
        }
        if batches_boundary_processed != total_read {
            counter_mismatches.push(format!(
                "batches_boundary_processed ({batches_boundary_processed}) != total_read ({total_read})"
            ));
        }
        // Note: batches_boundary_found may differ from total_read because FindBoundaries
        // can split or combine batches. But batches_decoded should match batches_boundary_found.
        if batches_decoded != batches_boundary_found {
            counter_mismatches.push(format!(
                "batches_decoded ({batches_decoded}) != batches_boundary_found ({batches_boundary_found})"
            ));
        }
        if batches_grouped != batches_boundary_found {
            counter_mismatches.push(format!(
                "batches_grouped ({batches_grouped}) != batches_boundary_found ({batches_boundary_found})"
            ));
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

    fn validate_completion(&self) -> Result<(), PipelineValidationError> {
        BamPipelineState::validate_completion(self)
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
        let mi_assign_reorder_len = self.output.mi_assign_reorder.lock().len();

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
            mi_assign_reorder_len,
            memory_limit: self.deadlock_state.get_memory_limit(),
            read_done: self.read_done.load(Ordering::Relaxed),
            group_done: self.group_done.load(Ordering::Relaxed),
            draining: self.output.draining.load(Ordering::Relaxed),
            extra_state: None,
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
        // Saturating (via `refund_queue_bytes`): this counter is summed into the
        // overflow-checked `queue_bytes_in_flight`, so an under-subtraction must
        // floor at zero rather than wrap to `u64::MAX`.
        refund_queue_bytes(&self.output.serialized_heap_bytes, heap_size);
    }

    fn q6_pop(&self) -> Option<(u64, CompressedBlockBatch)> {
        self.output.compressed.pop()
    }

    fn q6_push(
        &self,
        item: (u64, CompressedBlockBatch),
    ) -> Result<(), (u64, CompressedBlockBatch)> {
        let heap_size = item.1.estimate_heap_size();
        // Charge *before* publishing: once the batch is visible a consumer can
        // pop and refund it, so a charge recorded after the push can land after
        // that refund and strand a phantom nonzero counter (see `q6_track_pop`).
        // Refund if the push fails, since nothing was published.
        self.output.compressed_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
        let result = self.output.compressed.push(item);
        if result.is_err() {
            refund_queue_bytes(&self.output.compressed_heap_bytes, heap_size as u64);
        }
        result
    }

    fn q6_is_full(&self) -> bool {
        self.output.compressed.is_full()
    }

    fn q6_track_pop(&self, heap_size: u64) {
        // Saturating (via `refund_queue_bytes`): this counter is summed into the
        // overflow-checked `queue_bytes_in_flight`, so an under-subtraction must
        // floor at zero rather than wrap to `u64::MAX`. `q6_push` now charges
        // before publishing, closing the pop-before-charge race; the floor
        // stays as defense in depth against any residual accounting skew.
        refund_queue_bytes(&self.output.compressed_heap_bytes, heap_size);
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
// New Shared Traits (Phase 2 - Pipeline Consolidation)
// ============================================================================

impl<G: Send + MemoryEstimate + 'static, P: Send + MemoryEstimate + 'static>
    ProcessPipelineState<G, P> for BamPipelineState<G, P>
{
    fn process_input_pop(&self) -> Option<(u64, Vec<G>)> {
        let result = self.output.groups.pop();
        if result.is_some() {
            // Q4 byte accounting is handled by try_step_process (which does its own
            // direct pop). Do NOT refund here to avoid double-subtraction.
            self.deadlock_state.record_q4_pop();
        }
        result
    }

    fn process_output_is_full(&self) -> bool {
        self.output.processed.is_full()
    }

    fn process_output_push(&self, item: (u64, Vec<P>)) -> Result<(), (u64, Vec<P>)> {
        // Charge before publishing (via `push_charged`) so a consumer cannot pop
        // and refund the batch before the charge lands; it refunds on a failed
        // push. Charging after the push strands a phantom counter (see `q6_push`).
        let heap_size: u64 =
            item.1.iter().map(|p| MemoryEstimate::estimate_heap_size(p) as u64).sum();
        let result = push_charged(
            &self.output.processed,
            &self.output.processed_heap_bytes,
            heap_size,
            item,
        );
        if result.is_ok() {
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
            // Track memory being removed from processed queue. Saturating, like
            // every other debit summed into `queue_bytes_in_flight`: a raw
            // `fetch_sub` past zero wraps to `u64::MAX` and slams the `Read` gate
            // shut (or panics the overflow-checked sum in debug/coverage).
            let heap_size: u64 =
                batch.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>() as u64;
            refund_queue_bytes(&self.output.processed_heap_bytes, heap_size);
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
        // Charge before publishing so a consumer cannot pop and refund the batch
        // before the charge lands; refund on a failed push (see `q6_push`).
        self.output.serialized_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
        let result = self.output.serialized.push(item);
        if result.is_ok() {
            self.deadlock_state.record_q6_push();
        } else {
            refund_queue_bytes(&self.output.serialized_heap_bytes, heap_size as u64);
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

// The `Grouper` trait is defined once, in `fgumi-bam-io`, alongside the
// `DecodedRecord` it consumes. Re-exported so existing
// `crate::unified_pipeline::Grouper` paths keep resolving.
pub use fgumi_bam_io::Grouper;

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
    ///
    /// # Errors
    ///
    /// Returns an I/O error if grouping fails.
    pub fn process(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<G>> {
        self.grouper.add_records(records)
    }

    /// Signal EOF and get any remaining group.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if the grouper fails to finalize.
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

    /// Names of this step's buffers that still hold data, for completion
    /// validation.
    ///
    /// Two of the three states [`super::base::PipelineLifecycle::validate_completion`]
    /// documents live here rather than on the pipeline state: groups that were
    /// never pushed to Q4, and records the grouper is still holding after
    /// [`Self::finish`]. The state cannot reach them — they sit behind this
    /// step's own mutex — so the run function folds this list into the
    /// validation result.
    ///
    /// Empty means the Group step is genuinely drained.
    #[must_use]
    pub fn unflushed_buffers(&self) -> Vec<String> {
        let mut buffers = Vec::new();
        if !self.pending_groups.is_empty() {
            buffers.push(format!("group_pending_groups ({})", self.pending_groups.len()));
        }
        if self.grouper.has_pending() {
            buffers.push("grouper (partial group pending)".to_string());
        }
        buffers
    }
}

// ============================================================================
// Step Functions Trait (for BAM pipeline)
// ============================================================================

/// Position of a single processed item within its pipeline batch, published
/// by the BAM pipeline harness when invoking [`PipelineFunctions::mi_assign_fn`].
///
/// `batch_serial` is the monotonically-increasing serial assigned by the
/// Group step (one per `(serial, Vec<P>)` push to Q5), and `idx_in_batch`
/// identifies one item within that batch. Together `(batch_serial,
/// idx_in_batch)` uniquely names every `ProcessedPositionGroup` produced by
/// a run of the pipeline, independent of completion timing.
///
/// `batch_len` is the total number of items in the enclosing batch; the
/// hook can compare `idx_in_batch + 1 == batch_len` to detect the final
/// item if it needs to finalize per-batch state before returning.
///
/// See `docs/design/deterministic-mi-numbering.md` for the design that
/// motivates this struct.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BatchOrdinal {
    /// Group-step serial of the enclosing batch.
    pub batch_serial: u64,
    /// Zero-based index of this item within the batch.
    pub idx_in_batch: usize,
    /// Total number of items in the batch.
    pub batch_len: usize,
}

impl BatchOrdinal {
    /// True when this is the final item in the batch.
    ///
    /// Convenient shorthand for `self.idx_in_batch + 1 == self.batch_len`,
    /// used by `mi_assign_fn` callers that need to run per-batch
    /// finalization on the last item only.
    #[inline]
    #[must_use]
    pub fn is_last(&self) -> bool {
        self.idx_in_batch + 1 == self.batch_len
    }
}

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

    /// Optional secondary serialization (e.g., rejected records).
    /// Called with a borrow of P BEFORE the primary `serialize_fn` consumes it.
    pub secondary_serialize_fn:
        Option<Box<dyn Fn(&P, &mut Vec<u8>) -> io::Result<u64> + Send + Sync>>,

    /// Optional serial-ordered hook that runs after `process_fn` and before
    /// `serialize_fn`, with mutable access to each batched item.
    ///
    /// When set, the BAM pipeline routes batches through a single-threaded
    /// "MI Assign" stage that pops in pipeline-serial order and invokes
    /// this closure once per item, in `(batch_serial, idx_in_batch)`
    /// order, before the batch is handed to `serialize_fn`. The hook is
    /// permitted to mutate the item — `fgumi group` and `fgumi dedup` use
    /// it to rewrite each template's local `MoleculeId` into a global one
    /// via `MoleculeId::with_offset`, so that `serialize_fn` itself stays
    /// synchronization-free.
    ///
    /// When `None` (the default), the pipeline behaves exactly as before:
    /// Process pushes onto Q5 and Serialize pops from Q5 in completion
    /// order. Commands that do not need cumulative ordering should leave
    /// this field unset.
    pub mi_assign_fn: Option<Box<dyn Fn(BatchOrdinal, &mut P) -> io::Result<()> + Send + Sync>>,
}

impl<G: Send, P: Send> PipelineFunctions<G, P> {
    /// Create new step functions.
    pub fn new<ProcessFn, SerializeFn>(process_fn: ProcessFn, serialize_fn: SerializeFn) -> Self
    where
        ProcessFn: Fn(G) -> io::Result<P> + Send + Sync + 'static,
        SerializeFn: Fn(P, &mut Vec<u8>) -> io::Result<u64> + Send + Sync + 'static,
    {
        Self {
            process_fn: Box::new(process_fn),
            serialize_fn: Box::new(serialize_fn),
            secondary_serialize_fn: None,
            mi_assign_fn: None,
        }
    }

    /// Attach a secondary serialize function for dual-output pipelines.
    #[must_use]
    pub fn with_secondary_serialize<F>(mut self, f: F) -> Self
    where
        F: Fn(&P, &mut Vec<u8>) -> io::Result<u64> + Send + Sync + 'static,
    {
        self.secondary_serialize_fn = Some(Box::new(f));
        self
    }

    /// Install a serial-ordered MI-assign hook. See
    /// [`PipelineFunctions::mi_assign_fn`] for the contract.
    #[must_use]
    pub fn with_mi_assign<F>(mut self, f: F) -> Self
    where
        F: Fn(BatchOrdinal, &mut P) -> io::Result<()> + Send + Sync + 'static,
    {
        self.mi_assign_fn = Some(Box::new(f));
        self
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
                ActiveSteps::all(),
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

impl<P: Send> HasRecycledBuffers for WorkerState<P> {
    fn take_or_alloc_buffer(&mut self, capacity: usize) -> Vec<u8> {
        self.core.take_or_alloc_buffer(capacity)
    }

    fn recycle_buffer(&mut self, buf: Vec<u8>) {
        self.core.recycle_buffer(buf);
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
        // `q1_push` charges the batch's heap bytes before the push and refunds
        // on the out-of-slots path; see `q2b_push` for why an uncharged-but-
        // visible batch permanently over-states the Read gate.
        match state.q1_push(serial, held) {
            Ok(()) => {
                // Successfully advanced held item, continue to read more
                state.deadlock_state.record_q1_push();
            }
            Err((serial, held)) => {
                // Still can't push - put it back and signal output full.
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
    // Priority 3b: Check the queue memory budget
    // =========================================================================
    // Every other stage is bounded by a slot count, not by bytes, so this is
    // the only place the configured budget can actually cap the pipeline. See
    // `BamPipelineState::read_admission_allowed` for why declining here cannot
    // deadlock.
    if !state.read_admission_allowed() {
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
            // `q1_push` charges the batch's heap bytes before the push and
            // refunds on the out-of-slots path. See `q2b_push`.
            match state.q1_push(serial, batch) {
                Ok(()) => {
                    state.deadlock_state.record_q1_push();
                    true
                }
                Err((serial, batch)) => {
                    // Output full - hold the result for the next attempt.
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
/// the batch is stored in `worker.held_decompressed` and the function returns
/// immediately. Held batches push unconditionally at Priority 1 (physical
/// capacity only — no memory check). Memory backpressure is applied between
/// P1 and P3 to gate new work when the reorder buffer is large.
fn try_step_decompress<G: Send, P: Send + MemoryEstimate>(
    state: &BamPipelineState<G, P>,
    worker: &mut WorkerState<P>,
) -> bool {
    // =========================================================================
    // Priority 1: Try to advance any held decompressed batch first
    // =========================================================================
    // Push unconditionally — only physical queue capacity can block.
    // Memory backpressure is NOT checked here: the held batch is already in
    // memory (just untracked while held), so pushing it doesn't increase actual
    // usage. Checking can_proceed here caused deadlocks when all workers held
    // non-next_seq batches and nobody could produce next_seq.
    let mut advanced_held = false;
    if let Some((serial, held, heap_size)) = worker.held_decompressed.take() {
        state.q2_reorder_state.add_heap_bytes(heap_size as u64);
        match state.q2_decompressed.push((serial, held)) {
            Ok(()) => {
                state.batches_decompressed.fetch_add(1, Ordering::Release);
                state.deadlock_state.record_q2_push();
                advanced_held = true;
            }
            Err((serial, held)) => {
                state.q2_reorder_state.sub_heap_bytes(heap_size as u64);
                worker.held_decompressed = Some((serial, held, heap_size));
                return false;
            }
        }
    }

    // =========================================================================
    // Priority 2: Check backpressure (physical capacity AND memory)
    // =========================================================================
    // Memory check gates new work when the Q2 reorder buffer is large.
    // This is safe from deadlock: Q1 is FIFO, so `next_seq` for q2_reorder
    // has already been decompressed (it's in q2_decompressed, q2_reorder,
    // or a held slot pushed at P1). FindBoundaries continues draining →
    // memory drops → workers resume. Same pattern as Process step's
    // `should_apply_process_backpressure()`.
    if state.q2_decompressed.is_full() || state.q2_reorder_state.is_memory_high() {
        return advanced_held;
    }

    // =========================================================================
    // Priority 3: Pop input and process
    // =========================================================================
    // `q1_pop` refunds the batch's heap bytes as it leaves the queue.
    let Some((serial, raw_batch)) = state.q1_pop() else {
        if let Some(stats) = state.stats() {
            stats.record_queue_empty(1);
        }
        return advanced_held;
    };
    state.deadlock_state.record_q1_pop();

    // Prepare worker's buffer: clear and reserve capacity
    worker.decompression_buffer.clear();
    let expected_size = raw_batch.total_uncompressed_size();
    worker.decompression_buffer.reserve(expected_size);

    // Decompress directly into worker's buffer (no intermediate allocations)
    for block in &raw_batch.blocks {
        if let Err(e) = decompress_block_into_opts(
            block,
            &mut worker.decompressor,
            &mut worker.decompression_buffer,
            state.config.verify_crc,
        ) {
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
    // Priority 4: Calculate and reserve memory for tracking
    // =========================================================================
    let batch = DecompressedBatch { data: decompressed };
    let heap_size = batch.estimate_heap_size();
    state.q2_reorder_state.add_heap_bytes(heap_size as u64);

    // =========================================================================
    // Priority 5: Try to push result
    // =========================================================================
    match state.q2_decompressed.push((serial, batch)) {
        Ok(()) => {
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

/// Whether `FindBoundaries` has consumed all its input and pushed all its output,
/// and so may set `boundary_done`.
///
/// Two things must hold, and the second is easy to miss.
///
/// **All input consumed.** `batches_boundary_processed == next_read_serial` implies
/// `Decompress` also finished, since no more can be processed than was
/// decompressed, so this also rules out data still sitting in `q2_decompressed`.
///
/// **No output still outstanding.** `batches_boundary_processed` is incremented
/// when a batch is popped from `q2_reorder` — *before* its boundaries are pushed —
/// so a batch whose push failed is already counted as processed while it waits in
/// some worker's `held_boundaries`, in no queue and invisible to the other workers.
/// Setting `boundary_done` there lets `Group` finish on
/// `batches_grouped == batches_boundary_found`, which excludes that batch, and latch
/// its one-way `finished` flag; the batch is then pushed and decoded, and reaches an
/// already-finalized grouper.
///
/// `next_boundary_serial == batches_boundary_found` is what rules that out. A serial
/// is minted only immediately before a push is attempted, and
/// `batches_boundary_found` is incremented only when a push succeeds (including the
/// retry of a held batch, which reuses its serial rather than minting a new one), so
/// the two are equal exactly when nothing is outstanding. Note this cannot be done
/// with `batches_boundary_processed`: a batch with no records is counted as processed
/// but mints no serial and is never pushed, so `found == processed` is not an
/// invariant even when nothing is held.
///
/// `batches_boundary_found` is loaded before `next_boundary_serial` so the comparison
/// errs toward "not complete": a batch outstanding at the first load forces
/// `next_boundary_serial` above the value read, and it only ever grows.
fn boundary_finding_is_complete<G: Send, P: Send + MemoryEstimate>(
    state: &BamPipelineState<G, P>,
) -> bool {
    state.read_done.load(Ordering::Acquire)
        && state.batches_boundary_processed.load(Ordering::Acquire)
            == state.next_read_serial.load(Ordering::Acquire)
        && state.batches_boundary_found.load(Ordering::Acquire)
            == state.next_boundary_serial.load(Ordering::Acquire)
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
#[allow(clippy::too_many_lines)]
fn try_step_find_boundaries<G: Send, P: Send + MemoryEstimate>(
    state: &BamPipelineState<G, P>,
    worker: &mut WorkerState<P>,
) -> (bool, bool) {
    const MAX_BATCHES_PER_LOCK: usize = 8;

    // =========================================================================
    // Priority 1: Try to advance any held boundary batch first
    // =========================================================================
    let mut did_work = false;
    if let Some((serial, held)) = worker.held_boundaries.take() {
        match state.q2b_push(serial, held) {
            Ok(()) => {
                // Successfully advanced held item, increment completion counter
                state.batches_boundary_found.fetch_add(1, Ordering::Release);
                state.deadlock_state.record_q2b_push();
                did_work = true;
            }
            Err((serial, held)) => {
                // Still can't push - put it back and signal backpressure. The serial
                // it already holds keeps it visible to `boundary_finding_is_complete`.
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
            if !did_work && let Some(stats) = state.stats() {
                stats.record_queue_empty(2);
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
                    match state.q2b_push(serial, boundary_batch) {
                        Ok(()) => {
                            // Successfully pushed, increment completion counter
                            state.batches_boundary_found.fetch_add(1, Ordering::Release);
                            state.deadlock_state.record_q2b_push();
                        }
                        Err((serial, boundary_batch)) => {
                            // Output full - hold the result and stop processing. The
                            // serial is already minted and `batches_boundary_found` is
                            // not, which is what keeps this batch visible to
                            // `boundary_finding_is_complete` while it waits.
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

    // No batches processed - check if we should finish. The flag is tested first so
    // the drain-down spin, which reaches here on every pass, does not re-read the
    // counters only to be rejected by an already-set flag.
    if !state.boundary_done.load(Ordering::Acquire) && boundary_finding_is_complete(state) {
        // All input processed - finish and emit any remaining boundaries
        match boundary_guard.finish() {
            Ok(Some(final_batch)) => {
                if final_batch.offsets.len() > 1 {
                    let serial = state.next_boundary_serial.fetch_add(1, Ordering::SeqCst);
                    // Try non-blocking push for final batch
                    match state.q2b_push(serial, final_batch) {
                        Ok(()) => {
                            // Successfully pushed, increment completion counter
                            state.batches_boundary_found.fetch_add(1, Ordering::Release);
                            state.deadlock_state.record_q2b_push();
                        }
                        Err((serial, final_batch)) => {
                            // Hold final batch for next attempt. `boundary_done` is
                            // deliberately not set on this path — the retry re-enters
                            // `finish()` after Priority 1 pushes this batch.
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
/// the batch is stored in `worker.held_decoded` and the function returns
/// immediately. Held batches push unconditionally at Priority 1 (physical
/// capacity only — no memory check). Memory backpressure is applied between
/// P1 and P3 to gate new work when the reorder buffer is large.
fn try_step_decode<G: Send, P: Send + MemoryEstimate>(
    state: &BamPipelineState<G, P>,
    worker: &mut WorkerState<P>,
) -> bool {
    // =========================================================================
    // Priority 1: Try to advance any held decoded batch first
    // =========================================================================
    // Push unconditionally — only physical queue capacity can block.
    // See Decompress Priority 1 for rationale.
    let mut advanced_held = false;
    if let Some((serial, held, heap_size)) = worker.held_decoded.take() {
        state.q3_reorder_state.add_heap_bytes(heap_size as u64);
        match state.q3_decoded.push((serial, held)) {
            Ok(()) => {
                state.batches_decoded.fetch_add(1, Ordering::Release);
                state.deadlock_state.record_q3_push();
                advanced_held = true;
            }
            Err((serial, held)) => {
                state.q3_reorder_state.sub_heap_bytes(heap_size as u64);
                worker.held_decoded = Some((serial, held, heap_size));
                return false;
            }
        }
    }

    // =========================================================================
    // Priority 2: Check physical-capacity backpressure
    // =========================================================================
    // Only the output queue being physically full blocks new work here. Memory
    // backpressure is applied *per serial* after the pop (below), because a
    // blanket `is_memory_high()` guard at this point is not deadlock-safe for
    // Decode the way it is for Decompress: it can stall the very serial the Q3
    // reorder buffer needs next, starving Group (fgumi#746).
    if state.q3_decoded.is_full() {
        return advanced_held;
    }

    // =========================================================================
    // Priority 3: Pop input
    // =========================================================================
    let Some((serial, boundary_batch)) = state.q2b_pop() else {
        if let Some(stats) = state.stats() {
            stats.record_queue_empty(25); // Q2b (boundaries queue)

            // Q2b is an extension of Q2, but an empty Q2b only counts as a stall
            // while `FindBoundaries` might still deliver more. Once it is done and the
            // queue has drained, this is the terminal state, not starvation.
            if !state.boundary_done.load(Ordering::SeqCst) || !state.q2b_boundaries.is_empty() {
                stats.record_queue_empty(2);
            }
        }
        return advanced_held;
    };

    // Per-serial memory backpressure (fgumi#746). A blanket `is_memory_high()`
    // guard before the pop instead stalled the gap-filler: when a slow writer
    // backed the pipeline up, Decode declined to decode the very serial the Q3
    // reorder buffer was waiting on, so it never released to Group and the
    // pipeline wedged. `admit_or_requeue_decode` applies the buffer's own
    // per-serial `can_proceed` after the pop instead — see it for the requeue
    // and bounded-overshoot semantics.
    let Some((serial, boundary_batch)) = state.admit_or_requeue_decode(serial, boundary_batch)
    else {
        return advanced_held;
    };
    state.deadlock_state.record_q2b_pop();

    // =========================================================================
    // Priority 4: Decode records and compute GroupKey
    // =========================================================================
    match decode_records(&boundary_batch, state.group_key_config.as_ref()) {
        Ok(records) => {
            // Record decoded count for throughput metrics
            if let Some(stats) = state.stats() {
                stats.records_decoded.fetch_add(records.len() as u64, Ordering::Relaxed);
            }

            // Calculate and reserve memory for tracking
            let heap_size = records.estimate_heap_size();
            state.q3_reorder_state.add_heap_bytes(heap_size as u64);

            // =========================================================================
            // Priority 5: Try to push result
            // =========================================================================
            match state.q3_decoded.push((serial, records)) {
                Ok(()) => {
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
#[allow(clippy::too_many_lines)]
fn try_step_group<G: Send + BatchWeight + MemoryEstimate + 'static, P: Send + MemoryEstimate>(
    state: &BamPipelineState<G, P>,
    group_state: &Mutex<GroupState<G>>,
) -> (bool, bool) {
    const MAX_BATCHES_PER_LOCK: usize = 8;
    const MAX_PENDING_DRAIN: usize = 16;

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
    // Q4 bytes are charged here against the queue memory budget and refunded by
    // `try_step_process` on pop. Q3 keeps its own separate tracker.
    let push_batch = |groups: Vec<G>, state: &BamPipelineState<G, P>| -> Result<(), Vec<G>> {
        if state.output.groups.is_full() {
            return Err(groups);
        }

        // Charge Q4 against the queue memory budget. This is O(records) in the
        // batch, on the exclusive Group step, but it is a `capacity()` read per
        // record against grouping work that already walked every one of them —
        // and leaving Q4 uncharged is what let a stalled writer park unbounded
        // grouped input here (issue #746).
        let heap_size: u64 = groups.iter().map(|g| g.estimate_heap_size() as u64).sum();
        state.output.groups_heap_bytes.fetch_add(heap_size, Ordering::AcqRel);

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
        // Try to flush remaining groups non-blocking. If Q4 is full, return
        // and let the scheduler retry — the thread can do useful Process work
        // in between to drain Q4. The deadlock detector handles true deadlocks.
        if flush_all(&mut guard, state).is_some() {
            state.group_done.store(true, Ordering::SeqCst);
            return (true, false);
        }
        return (false, false); // Q4 full, scheduler will retry
    }

    // Process multiple record batches per lock acquisition to reduce contention
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
            if !did_work && let Some(stats) = state.stats() {
                stats.record_queue_empty(3);
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

                // Flush any remaining groups non-blocking. If Q4 is full,
                // return and let the scheduler retry via is_finished() check.
                if flush_all(&mut guard, state).is_some() {
                    state.group_done.store(true, Ordering::SeqCst);
                }
                (true, false) // Did work (finished grouper)
            }
            Ok(None) => {
                // No final group — flush any pending. If Q4 is full, return
                // and let the scheduler retry via is_finished() check.
                if flush_all(&mut guard, state).is_some() {
                    state.group_done.store(true, Ordering::SeqCst);
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
    const MAX_BATCHES: usize = 8;

    // =========================================================================
    // Priority 1: Try to advance any held processed batch first
    // =========================================================================
    if let Some((serial, held, heap_size)) = worker.held_processed.take() {
        // Charge before publishing so a consumer cannot pop and refund the batch
        // before the charge lands; refund on a failed push (see `q6_push`).
        state.output.processed_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
        match state.output.processed.push((serial, held)) {
            Ok(()) => {
                state.deadlock_state.record_q5_push();
            }
            Err((serial, held)) => {
                // Still can't push - refund the charge, put it back, signal output full
                refund_queue_bytes(&state.output.processed_heap_bytes, heap_size as u64);
                worker.held_processed = Some((serial, held, heap_size));
                return false;
            }
        }
    }

    // =========================================================================
    // Priority 2: Check if output queue has space (count and memory)
    // Memory backpressure is always enforced (including during draining) to
    // prevent OOM.  The slot-based is_full() check is sufficient to guarantee
    // forward progress: Serialize drains Q5 → slots free → Process resumes.
    // =========================================================================
    if state.output.processed.is_full() || state.is_q5_memory_high() {
        return false;
    }

    // =========================================================================
    // Priority 3: Pop and process batches
    // Always drain multiple batches when work is available for better throughput.
    // Q5 memory backpressure above prevents unbounded growth.
    // =========================================================================
    let mut did_work = false;

    for _ in 0..MAX_BATCHES {
        // Check output space (count and memory) before each batch
        if state.output.processed.is_full() || state.is_q5_memory_high() {
            break;
        }

        let Some((serial, batch)) = state.output.groups.pop() else {
            if let Some(stats) = state.stats() {
                stats.record_queue_empty(4);
            }
            break;
        };
        state.deadlock_state.record_q4_pop();

        // Refund Q4's charge. The batch is unchanged since `push_batch`
        // charged it, so the two estimates agree exactly.
        let q4_heap: u64 = batch.iter().map(|g| g.estimate_heap_size() as u64).sum();
        refund_queue_bytes(&state.output.groups_heap_bytes, q4_heap);

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

        // Charge before publishing so a consumer cannot pop and refund the batch
        // before the charge lands; refund on a failed push (see `q6_push`).
        state.output.processed_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
        match state.output.processed.push((serial, results)) {
            Ok(()) => {
                state.deadlock_state.record_q5_push();
                did_work = true;
            }
            Err((serial, results)) => {
                // Output full - refund the charge and hold the result for next attempt
                refund_queue_bytes(&state.output.processed_heap_bytes, heap_size as u64);
                worker.held_processed = Some((serial, results, heap_size));
                break;
            }
        }
    }

    did_work
}

/// Outcome of an attempt to pop a serial-ordered batch out of the MI Assign
/// reorder zone. Distinguishing `Stalled` from `Empty` lets the caller
/// avoid polluting the queue-emptiness stat with reorder-buffer waits;
/// distinguishing `Errored` from both lets the call site treat hook
/// failures as the terminal abort they actually are.
enum MiAssignPopOutcome<P> {
    /// A batch is ready to serialize, already mutated by the hook.
    Ready { serial: u64, batch: Vec<P> },
    /// Q5 was drained into the reorder buffer, but the next-expected
    /// pipeline-serial-order batch hasn't arrived yet — there is in-flight
    /// work parked in the reorder buffer, so this is *not* a Q5-empty
    /// stall.
    Stalled,
    /// Q5 was empty and the reorder buffer was empty — nothing to do.
    Empty,
    /// The hook returned an error; `state.set_error` has already been
    /// invoked and the worker loop will short-circuit on the next
    /// iteration. Functionally equivalent to `Stalled` at the call
    /// site, but named distinctly so the contract is self-documenting.
    Errored,
}

/// Drain `processed` (Q5) into the per-pipeline MI-assign reorder buffer
/// and try to pop the next-in-pipeline-serial-order batch.
///
/// Called only when `mi_assign_fn` is installed. Holds the
/// `mi_assign_reorder` mutex for the duration of (drain + pop + hook
/// invocation), so the hook always sees `(batch_serial, idx_in_batch)`
/// pairs in monotonic order even though parallel workers race to pop
/// from `processed`.
fn try_pop_mi_assigned<G: Send + 'static, P: Send + MemoryEstimate + 'static>(
    state: &BamPipelineState<G, P>,
    hook: &(dyn Fn(BatchOrdinal, &mut P) -> io::Result<()> + Send + Sync),
) -> MiAssignPopOutcome<P> {
    let mut reorder = state.output.mi_assign_reorder.lock();

    // Drain whatever is currently in Q5 into the reorder buffer. This keeps
    // Q5 drainable (so producers don't get stuck on slot backpressure) and
    // gives the reorder buffer a chance to receive the next-expected serial.
    //
    // Memory accounting: heap bytes stay charged to `processed_heap_bytes`
    // while the batch is parked in the reorder buffer (carried in via
    // `insert_with_size`, refunded on `try_pop_next_with_size`). This keeps
    // `is_q5_memory_high()` honest so producers stay backpressured if many
    // batches stack up waiting for an out-of-order serial.
    while let Some((serial, batch)) = state.output.processed.pop() {
        let heap_size: usize = batch.iter().map(MemoryEstimate::estimate_heap_size).sum();
        state.deadlock_state.record_q5_pop();
        reorder.insert_with_size(serial, batch, heap_size);
    }

    // Pop the next batch in pipeline-serial order if it's available.
    let serial = reorder.next_seq();
    let Some((mut batch, heap_size)) = reorder.try_pop_next_with_size() else {
        // No in-order batch — distinguish a true Q5-empty stall (reorder
        // buffer also empty after the drain) from a reorder-stalled wait
        // (reorder buffer non-empty, just missing an earlier serial).
        return if reorder.is_empty() {
            MiAssignPopOutcome::Empty
        } else {
            MiAssignPopOutcome::Stalled
        };
    };
    // Saturating, like every other debit summed into `queue_bytes_in_flight`.
    refund_queue_bytes(&state.output.processed_heap_bytes, heap_size as u64);

    // Run the hook on every item in the popped batch, in order. We hold the
    // reorder mutex for the whole loop so that two workers cannot interleave
    // hook invocations across batches and race on whatever shared state the
    // hook is advancing (e.g. `group`/`dedup`'s cumulative MI counter).
    let batch_len = batch.len();
    for (idx_in_batch, item) in batch.iter_mut().enumerate() {
        let ordinal = BatchOrdinal { batch_serial: serial, idx_in_batch, batch_len };
        if let Err(e) = (hook)(ordinal, item) {
            state.set_error(e);
            return MiAssignPopOutcome::Errored;
        }
    }

    MiAssignPopOutcome::Ready { serial, batch }
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
        // Charge before publishing so a consumer cannot pop and refund the batch
        // before the charge lands; refund on a failed push (see `q6_push`).
        state.output.serialized_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
        match state.output.serialized.push((serial, held)) {
            Ok(()) => {
                // Successfully advanced held item
                state.deadlock_state.record_q6_push();
            }
            Err((serial, held)) => {
                // Still can't push - refund the charge, put it back, signal output full
                refund_queue_bytes(&state.output.serialized_heap_bytes, heap_size as u64);
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
    // Priority 3: Pop input.
    //
    // When `mi_assign_fn` is installed, route the batch through a serial-
    // ordered "MI Assign" zone so the hook always observes batches (and
    // items within them) in monotonic `(batch_serial, idx_in_batch)` order
    // even though `try_step_serialize` itself runs on parallel workers.
    // The hook may mutate items in place (e.g. fold a cumulative offset
    // into per-template `MoleculeId`s); the subsequent serialize work then
    // proceeds on the mutated batch outside the reorder lock.
    //
    // When `mi_assign_fn` is `None`, the existing fast path applies:
    // pop directly from `processed` (Q5) in completion order. No mutex,
    // no reorder buffer, no behavior change for any other command.
    // =========================================================================
    let (serial, batch) = if let Some(ref hook) = fns.mi_assign_fn {
        match try_pop_mi_assigned(state, hook.as_ref()) {
            MiAssignPopOutcome::Ready { serial, batch } => (serial, batch),
            // Stalled: reorder buffer has batches parked but the
            // next-expected serial isn't here yet — Q5 is *not* empty
            // in any pipeline-meaningful sense, so counting this as
            // Q5-empty would skew the queue-emptiness stat.
            // Errored: hook returned an error and `set_error` is already
            // recorded; the worker loop's `state.has_error()` check tears
            // the pipeline down on the next iteration. Both cases yield
            // here without recording a Q5-empty stat.
            MiAssignPopOutcome::Stalled | MiAssignPopOutcome::Errored => return false,
            MiAssignPopOutcome::Empty => {
                if let Some(stats) = state.stats() {
                    stats.record_queue_empty(5);
                }
                return false;
            }
        }
    } else {
        let Some(item) = state.output.processed.pop() else {
            if let Some(stats) = state.stats() {
                stats.record_queue_empty(5);
            }
            return false;
        };
        state.deadlock_state.record_q5_pop();
        // Track memory being removed from Q5 (only when we go straight from
        // Q5 to Serialize; the MI Assign path subtracts after the reorder
        // buffer hands the batch off to serialize).
        // Saturating, like every other debit summed into `queue_bytes_in_flight`.
        let q5_heap_size: u64 =
            item.1.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>() as u64;
        refund_queue_bytes(&state.output.processed_heap_bytes, q5_heap_size);
        item
    };

    // =========================================================================
    // Priority 4: Serialize all items
    // =========================================================================
    // Prepare worker's serialization buffers
    worker.core.serialization_buffer.clear();
    worker.core.secondary_serialization_buffer.clear();

    // Serialize all items into worker's buffer
    let mut total_record_count: u64 = 0;
    for item in batch {
        // Secondary serialize (borrows item) — must run before primary consumes it
        if let Some(ref secondary_fn) = fns.secondary_serialize_fn
            && let Err(e) = (secondary_fn)(&item, &mut worker.core.secondary_serialization_buffer)
        {
            state.set_error(e);
            return false;
        }
        // Primary serialize (consumes item)
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

    // Build secondary data if any was serialized
    let secondary_data = if worker.core.secondary_serialization_buffer.is_empty() {
        None
    } else {
        Some(std::mem::replace(
            &mut worker.core.secondary_serialization_buffer,
            Vec::with_capacity(SERIALIZATION_BUFFER_CAPACITY),
        ))
    };

    // Record serialized bytes for throughput metrics
    if let Some(stats) = state.stats() {
        stats.serialized_bytes.fetch_add(combined_data.len() as u64, Ordering::Relaxed);
    }

    // =========================================================================
    // Priority 5: Try to push result (non-blocking)
    // =========================================================================
    let batch =
        SerializedBatch { data: combined_data, record_count: total_record_count, secondary_data };
    let heap_size = batch.estimate_heap_size();
    // Charge before publishing so a consumer cannot pop and refund the batch
    // before the charge lands; refund on a failed push (see `q6_push`).
    state.output.serialized_heap_bytes.fetch_add(heap_size as u64, Ordering::AcqRel);
    match state.output.serialized.push((serial, batch)) {
        Ok(()) => {
            state.deadlock_state.record_q6_push();
            true
        }
        Err((serial, batch)) => {
            // Output full - refund the charge and hold the result for next attempt
            refund_queue_bytes(&state.output.serialized_heap_bytes, heap_size as u64);
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
fn try_step_write<G: Send + 'static, P: Send + MemoryEstimate + 'static>(
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

        // Drain Q7 into reorder buffer, tracking heap bytes for admission control.
        while let Some((serial, batch)) = state.output.compressed.pop() {
            let q7_heap = batch.estimate_heap_size() as u64;
            state.q6_track_pop(q7_heap);
            state.deadlock_state.record_q7_pop();
            let heap_size = batch.estimate_heap_size();
            reorder.insert_with_size(serial, batch, heap_size);
            state.output.write_reorder_state.add_heap_bytes(heap_size as u64);
        }

        // Write all ready batches
        while let Some((batch, heap_size)) = reorder.try_pop_next_with_size() {
            // Write all blocks in the batch
            let mut batch_bytes: u64 = 0;
            for block in &batch.blocks {
                if let Err(e) = writer.write_all(&block.data) {
                    state.set_error(e);
                    return (false, false); // Error, not contention
                }
                batch_bytes += block.data.len() as u64;
            }

            // Write secondary data (e.g., rejects) in the same serial order
            if let Some(ref secondary_data) = batch.secondary_data
                && !secondary_data.is_empty()
                && let Some(ref secondary_mutex) = state.output.secondary_output
            {
                let mut sw_guard = secondary_mutex.lock();
                if let Some(ref mut sw) = *sw_guard
                    && let Err(e) = sw.write_raw_bytes(secondary_data)
                {
                    state.set_error(e);
                    return (false, false);
                }
            }

            // Update admission-control state so the scheduler-level
            // write_reorder_is_memory_high signal stays accurate. (Compress
            // does not consume this gate — see WritePipelineState docs.)
            state.output.write_reorder_state.sub_heap_bytes(heap_size as u64);
            state.output.write_reorder_state.update_next_seq(reorder.next_seq());

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
    if !wrote_any
        && q7_truly_empty
        && let Some(stats) = state.stats()
    {
        stats.record_queue_empty(7);
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
            && self.state.read_admission_allowed()
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
    /// Buffer for secondary serialized data (e.g., rejected records).
    /// Only used when a secondary serialize function is set.
    secondary: Vec<u8>,
}

impl SingleThreadedBuffers {
    /// Create new buffers with reasonable initial capacity.
    fn new() -> Self {
        Self {
            // 4 blocks * 64KB max uncompressed = 256KB typical
            decompressed: Vec::with_capacity(256 * 1024),
            // Typical group serializes to ~64KB
            serialized: Vec::with_capacity(64 * 1024),
            secondary: Vec::new(),
        }
    }
}

/// Type alias for an `mi_assign_fn` reference, used by
/// [`run_mi_assign_single_threaded`] to keep the parameter list readable
/// for clippy's `type_complexity` lint.
type MiAssignFnRef<'a, P> =
    Option<&'a (dyn Fn(BatchOrdinal, &mut P) -> io::Result<()> + Send + Sync)>;

/// Optionally invoke the MI Assign hook on a processed group before it is
/// serialized.
///
/// The single-threaded BAM pipeline path is intrinsically in serial order,
/// so the hook receives a strictly increasing `batch_serial` (one per
/// processed group, treated as a one-item batch) with no synchronization
/// needed. Returns the consumed `processed` (the hook may have mutated it
/// in place) so the caller can pass it on to `serialize_fn`.
fn run_mi_assign_single_threaded<P>(
    mi_assign_fn: MiAssignFnRef<'_, P>,
    next_serial: &mut u64,
    mut processed: P,
) -> io::Result<P> {
    if let Some(hook) = mi_assign_fn {
        (hook)(
            BatchOrdinal { batch_serial: *next_serial, idx_in_batch: 0, batch_len: 1 },
            &mut processed,
        )?;
        *next_serial += 1;
    }
    Ok(processed)
}

/// Run the BAM pipeline in single-threaded mode.
///
/// This avoids the overhead of thread spawning, queues, and atomic
/// operations when only one thread is requested. Significantly faster
/// for small inputs or when parallelization overhead exceeds the benefit.
#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn run_bam_pipeline_single_threaded<G, P>(
    config: &PipelineConfig,
    mut input: Box<dyn Read + Send>,
    mut output: Box<dyn Write + Send>,
    mut grouper: Box<dyn Grouper<Group = G> + Send>,
    fns: PipelineFunctions<G, P>,
    group_key_config: Option<GroupKeyConfig>,
    mut secondary_writer: Option<fgumi_bam_io::RawBamWriter>,
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

    // Per-position-group serial counter, incremented every time a processed
    // group is handed off to `serialize_fn`. Only meaningful when
    // `fns.mi_assign_fn` is set; otherwise the harness never reads it. The
    // single-threaded path treats each processed group as a single-item
    // batch (`idx_in_batch=0, batch_len=1`), which is the natural mapping
    // since serialization here is synchronous and one-at-a-time.
    let mut next_serial: u64 = 0;

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
            decompress_block_into_opts(
                block,
                &mut decompressor,
                &mut buffers.decompressed,
                config.verify_crc,
            )?;
        }

        // Step 3: Find record boundaries
        let boundary_batch = boundary_state.find_boundaries(&buffers.decompressed)?;

        // Step 4: Decode records (only if there are any)
        if boundary_batch.offsets.len() > 1 {
            let decoded = decode_records(&boundary_batch, group_key_config.as_ref())?;

            // Step 5: Feed decoded records to grouper
            let groups = grouper.add_records(decoded)?;

            // Process each group through steps 6-9
            for group in groups {
                // Step 6: Process
                let processed = (fns.process_fn)(group)?;

                // Step 7a: Run the optional MI Assign hook so any per-template
                // rewriting it does is reflected in BOTH outputs. The parallel
                // path runs MI Assign before secondary_serialize, so doing it
                // here keeps single-threaded byte-for-byte equivalent.
                let processed = run_mi_assign_single_threaded(
                    fns.mi_assign_fn.as_deref(),
                    &mut next_serial,
                    processed,
                )?;

                // Step 7b: Secondary serialize (borrows processed)
                buffers.secondary.clear();
                if let Some(ref secondary_fn) = fns.secondary_serialize_fn {
                    (secondary_fn)(&processed, &mut buffers.secondary)?;
                }

                // Step 7c: Primary serialize (consumes processed, reuse buffer).
                buffers.serialized.clear();
                let record_count = (fns.serialize_fn)(processed, &mut buffers.serialized)?;

                // Step 8: Compress (only when buffer reaches 64KB)
                compressor.write_all(&buffers.serialized)?;
                compressor.maybe_compress()?;

                // Step 9: Write any completed blocks to output
                compressor.write_blocks_to(output.as_mut())?;

                // Write secondary data
                if !buffers.secondary.is_empty()
                    && let Some(ref mut sw) = secondary_writer
                {
                    sw.write_raw_bytes(&buffers.secondary)?;
                }

                progress.log_if_needed(record_count);
            }
        }
    }

    // Handle any remaining bytes from boundary finding
    if let Some(final_batch) = boundary_state.finish()?
        && final_batch.offsets.len() > 1
    {
        let decoded = decode_records(&final_batch, group_key_config.as_ref())?;
        let groups = grouper.add_records(decoded)?;

        for group in groups {
            let processed = (fns.process_fn)(group)?;

            // Run MI Assign first so secondary_serialize sees the same
            // post-assign state as the parallel path (matches Step 7 above).
            let processed = run_mi_assign_single_threaded(
                fns.mi_assign_fn.as_deref(),
                &mut next_serial,
                processed,
            )?;

            buffers.secondary.clear();
            if let Some(ref secondary_fn) = fns.secondary_serialize_fn {
                (secondary_fn)(&processed, &mut buffers.secondary)?;
            }

            buffers.serialized.clear();
            let record_count = (fns.serialize_fn)(processed, &mut buffers.serialized)?;
            compressor.write_all(&buffers.serialized)?;
            compressor.maybe_compress()?;
            compressor.write_blocks_to(output.as_mut())?;

            if !buffers.secondary.is_empty()
                && let Some(ref mut sw) = secondary_writer
            {
                sw.write_raw_bytes(&buffers.secondary)?;
            }

            progress.log_if_needed(record_count);
        }
    }

    // Finish grouper - process any remaining partial group
    if let Some(final_group) = grouper.finish()? {
        // Step 6: Process
        let processed = (fns.process_fn)(final_group)?;

        // Step 7a: Run MI Assign first so per-template rewriting is reflected
        // in BOTH outputs and matches the parallel path's ordering.
        let processed = run_mi_assign_single_threaded(
            fns.mi_assign_fn.as_deref(),
            &mut next_serial,
            processed,
        )?;

        // Step 7b: Secondary serialize (borrows processed)
        buffers.secondary.clear();
        if let Some(ref secondary_fn) = fns.secondary_serialize_fn {
            (secondary_fn)(&processed, &mut buffers.secondary)?;
        }

        // Step 7c: Primary serialize (consumes processed, reuse buffer).
        buffers.serialized.clear();
        let record_count = (fns.serialize_fn)(processed, &mut buffers.serialized)?;

        // Step 8: Compress (only when buffer reaches 64KB)
        compressor.write_all(&buffers.serialized)?;
        compressor.maybe_compress()?;

        // Step 9: Write any completed blocks to output
        compressor.write_blocks_to(output.as_mut())?;

        // Write secondary data
        if !buffers.secondary.is_empty()
            && let Some(ref mut sw) = secondary_writer
        {
            sw.write_raw_bytes(&buffers.secondary)?;
        }

        progress.log_if_needed(record_count);
    }

    // The threaded path validates completion against the pipeline state; this
    // path has no such state, so the same contract is checked directly here. A
    // grouper that still reports pending records after `finish()` has swallowed
    // them, and finishing the output now would leave a valid, complete-looking
    // BAM missing exactly those records.
    if grouper.has_pending() {
        return Err(io::Error::other(PipelineValidationError {
            non_empty_queues: vec!["grouper (partial group pending after finish)".to_string()],
            counter_mismatches: Vec::new(),
            leaked_heap_bytes: 0,
        }));
    }

    // Flush any remaining data in compression buffer
    compressor.flush()?;
    compressor.write_blocks_to(output.as_mut())?;

    // Flush output and write BGZF EOF marker
    output.flush()?;
    output.write_all(&BGZF_EOF)?;
    output.flush()?;

    // Finalize secondary output writer
    if let Some(writer) = secondary_writer {
        writer.finish().map_err(|e| {
            io::Error::new(e.kind(), format!("Failed to finalize secondary output: {e}"))
        })?;
    }

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
///
/// # Errors
///
/// Returns an I/O error if any pipeline step fails.
#[allow(clippy::too_many_lines, clippy::cast_possible_truncation)]
pub fn run_bam_pipeline<G, P>(
    config: PipelineConfig,
    input: Box<dyn Read + Send>,
    output: Box<dyn Write + Send>,
    grouper: Box<dyn Grouper<Group = G> + Send>,
    fns: PipelineFunctions<G, P>,
    group_key_config: Option<GroupKeyConfig>,
    secondary_writer: Option<fgumi_bam_io::RawBamWriter>,
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
            secondary_writer,
        );
    }

    let mut state = BamPipelineState::<G, P>::new(config, input, output, group_key_config);
    if let Some(sw) = secondary_writer {
        state.output.set_secondary_output(sw);
    }
    let state = Arc::new(state);

    // Set num_threads for stats display
    if let Some(stats) = state.stats() {
        stats.set_num_threads(num_threads);
        #[cfg(feature = "memory-debug")]
        stats.set_infrastructure_memory(num_threads, state.config.queue_capacity);
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

                // Queue memory from AtomicU64 counters
                // Q1: tracked unconditionally (charged against the queue budget)
                let q1_mem = state_clone.q1_heap_bytes.load(Ordering::Relaxed);
                // Q2-Q3: reorder buffer heap_bytes tracked unconditionally (used for backpressure)
                let q2_mem = state_clone.q2_reorder_state.heap_bytes.load(Ordering::Relaxed);
                // Q2b: tracked unconditionally (charged against the queue budget)
                let boundaries_mem = state_clone.q2b_heap_bytes.load(Ordering::Relaxed);
                let q3_mem = state_clone.q3_reorder_state.heap_bytes.load(Ordering::Relaxed);
                // Q4: tracked unconditionally (charged against the queue budget)
                let q4_mem = state_clone.output.groups_heap_bytes.load(Ordering::Relaxed);
                // Q5-Q7: tracked unconditionally via OutputPipelineQueues atomic counters
                let q5_mem = state_clone.output.processed_heap_bytes.load(Ordering::Relaxed);
                let q6_mem = state_clone.output.serialized_heap_bytes.load(Ordering::Relaxed);
                let q7_mem = state_clone.output.compressed_heap_bytes.load(Ordering::Relaxed);
                let queue_memory_bytes =
                    [q1_mem, q2_mem, boundaries_mem, q3_mem, q4_mem, q5_mem, q6_mem, q7_mem];

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
                    let total_mem = q1_mem
                        + q2_mem
                        + boundaries_mem
                        + q3_mem
                        + q7_reorder_mem
                        + q4_mem
                        + q5_mem
                        + q6_mem
                        + q7_mem;
                    stats.record_memory_usage(total_mem);

                    // Update shared PipelineStats with actual queue memory bytes
                    // so the memory debugging monitor can see them
                    #[cfg(feature = "memory-debug")]
                    stats.update_queue_memory_from_external(&[
                        ("q1", q1_mem),
                        ("q2", q2_mem),
                        ("q2b", boundaries_mem),
                        ("q3", q3_mem),
                        ("q4", q4_mem),
                        ("q5", q5_mem),
                        ("q6", q6_mem),
                        ("q7", q7_mem),
                    ]);

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
                        if let DeadlockAction::Detected =
                            check_deadlock_and_restore(&state_clone.deadlock_state, &snapshot)
                        {
                            state_clone.set_error(io::Error::new(
                                io::ErrorKind::TimedOut,
                                "pipeline deadlock detected with recovery disabled; \
                                 use --deadlock-recover to enable automatic recovery",
                            ));
                            break;
                        }
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

    // Finalize: check errors, flush output, log stats. The Group step's own
    // buffers are folded in: the state cannot reach them, and records stranded
    // in the grouper (or groups never pushed to Q4) are data loss that would
    // otherwise be reported as a clean run.
    let result = finalize_pipeline_with_buffers(&*state, || group_state.lock().unflushed_buffers());

    // Finalize secondary output writer (if present)
    if let Some(ref secondary_mutex) = state.output.secondary_output {
        let mut guard = secondary_mutex.lock();
        if let Some(writer) = guard.take()
            && let Err(e) = writer.finish().map_err(|e| {
                io::Error::new(e.kind(), format!("Failed to finalize secondary output: {e}"))
            })
        {
            if result.is_err() {
                log::error!("Secondary output finalization also failed: {e}");
            } else {
                return Err(e);
            }
        }
    }

    result
}

// ============================================================================
// BAM Pipeline Helpers
// ============================================================================

/// Serialize a single BAM record to bytes, appending to the provided buffer.
///
/// This produces raw BAM record bytes (`block_size` prefix + record data),
/// suitable for BGZF compression in the pipeline.
///
/// Returns the number of records serialized (always 1 for a single record).
///
/// Accepts `RecordBuf` because callers produce records via noodles typed operations.
///
/// # Errors
///
/// Returns an I/O error if record encoding fails.
pub fn serialize_bam_record_into(
    record: &RecordBuf,
    header: &Header,
    output: &mut Vec<u8>,
) -> io::Result<u64> {
    use noodles::sam::alignment::io::Write as _;

    // Writer<&mut Vec<u8>> writes the 4-byte little-endian `block_size` prefix
    // followed by the record payload directly into `output`, which matches the
    // exact on-disk BAM framing expected by the pipeline and BGZF compressor.
    // Writing in place avoids the per-call allocations that a temporary
    // `Vec<u8>` scratch would incur on every record in the hot path.
    let mut writer = noodles::bam::io::Writer::from(&mut *output);
    writer.write_alignment_record(header, record).map_err(io::Error::other)?;

    Ok(1)
}

/// Configuration for running a BAM file through the pipeline.
#[derive(Debug, Clone)]
pub struct BamPipelineConfig {
    /// Base pipeline configuration.
    pub pipeline: PipelineConfig,
    /// Compression level for output (0-12).
    pub compression_level: u32,
    /// Configuration for computing each record's `GroupKey` during decode.
    ///
    /// `Some(cfg)` computes the key from `cfg`; `None` skips the per-record
    /// computation entirely and attaches a default placeholder key. `None` is
    /// only valid for groupers that never read the key (e.g.
    /// `SingleRawRecordGrouper`) — a key-consuming grouper paired with `None`
    /// would mis-group (see [`decode_records`]). Every key-consuming command
    /// sets `Some` explicitly, so it is never left unset by accident.
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

/// Open an output writer for pipeline use, supporting stdout via "-" or "/dev/stdout".
///
/// Delegates to [`fgumi_bam_io::open_output_writer`], the single place the `-`
/// convention is honoured for BAM output. This previously re-implemented the
/// dispatch and returned a `LineWriter`-backed [`std::io::stdout`], which tears
/// every BGZF flush at each `0x0a`; `open_output_writer` returns a block-buffered
/// stdout handle instead.
fn open_pipeline_output(output_path: &Path) -> io::Result<Box<dyn Write + Send>> {
    fgumi_bam_io::open_output_writer(output_path).map_err(|e| {
        let kind = e.downcast_ref::<io::Error>().map_or(io::ErrorKind::Other, io::Error::kind);
        io::Error::new(kind, format!("{e:#}"))
    })
}

/// Convert an input-open failure into an `io::Error`, preserving its kind.
///
/// `create_bam_reader_for_pipeline` returns `anyhow::Error`; flattening every
/// cause to `InvalidData` would lose `NotFound`/`PermissionDenied`, which the
/// bare `File::open` this replaced used to surface.
#[cfg(any(test, feature = "test-utils"))]
fn open_input_error(e: &anyhow::Error) -> io::Error {
    let kind =
        e.downcast_ref::<io::Error>().map_or(io::ErrorKind::InvalidData, std::io::Error::kind);
    io::Error::new(kind, format!("{e:#}"))
}

/// Run a BAM file through the pipeline with a grouper factory.
///
/// This is a convenience function that handles BAM header I/O and
/// sets up the pipeline correctly for BAM processing.
///
/// **Test-only.** Behind the `test-utils` feature: every production command
/// drives the pipeline through the [`run_bam_pipeline_from_reader`] family,
/// which takes an already-opened reader and header. This path-taking wrapper
/// exists only so integration tests can go from a fixture path to a finished
/// output in one call, and it is not part of the supported public API.
///
/// # Type Parameters
///
/// - `G`: Group type produced by the grouper (e.g., `RecordBuf`, `RawPositionGroup`, `MiGroup`)
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
///
/// # Errors
///
/// Returns an I/O error if any pipeline step or file I/O fails.
#[cfg(any(test, feature = "test-utils"))]
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
    // Open once, through the normalizing reader: this accepts uncompressed SAM
    // and stdin, and hands back both the header and a stream replaying from byte
    // zero for the pipeline's block reader. Opening the path twice — once for the
    // header, once for the blocks — is what used to make these entry points
    // reject the inputs every command accepts.
    let (input, header) = fgumi_bam_io::create_bam_reader_for_pipeline(input_path)
        .map_err(|e| open_input_error(&e))?;

    // Create output writer (supports stdout via "-" or "/dev/stdout")
    let output_writer = open_pipeline_output(output_path)?;

    // Write BAM header using BGZF compression
    let mut header_writer = bam::io::Writer::new(output_writer);
    header_writer
        .write_header(&header)
        .map_err(|e| io::Error::other(format!("Failed to write BAM header: {e}")))?;

    // Finish the BGZF writer and get the underlying writer for the pipeline.
    // We need to:
    // 1. Get the BGZF writer from the BAM writer
    // 2. Flush/finish the BGZF stream (writes any pending data)
    // 3. Get the underlying writer
    // This ensures the pipeline writes raw BGZF blocks directly to the output,
    // not through another BGZF compression layer.
    let mut bgzf_writer = header_writer.into_inner();
    bgzf_writer
        .try_finish()
        .map_err(|e| io::Error::other(format!("Failed to finish BGZF header: {e}")))?;
    let output = bgzf_writer.into_inner();

    // Wrap in BufReader to reduce syscalls.
    let input = BufReader::with_capacity(IO_BUFFER_SIZE, input);

    let output = BufWriter::with_capacity(IO_BUFFER_SIZE, output);

    // Pass the caller's `Option` straight through: `Some` computes keys during
    // decode, `None` skips it for groupers that never read the key. Every
    // key-consuming command sets `Some` explicitly, so no default is needed here.
    let group_key_config = config.group_key_config;

    // Create the grouper
    // Note: Header skipping is now handled by BoundaryState in the pipeline
    let grouper = grouper_fn(&header);

    // Create step functions with header captured
    let header_clone = header.clone();
    let fns = PipelineFunctions::new(process_fn, move |p: P, buf: &mut Vec<u8>| {
        serialize_fn(p, &header_clone, buf)
    });

    run_bam_pipeline(
        config.pipeline,
        Box::new(input),
        Box::new(output),
        grouper,
        fns,
        group_key_config,
        None,
    )
}

// ============================================================================
// Reader-based Pipeline Functions (for streaming support)
// ============================================================================

/// Shared driver behind [`run_bam_pipeline_from_reader`] and
/// [`run_bam_pipeline_from_reader_with_mi_assign`].
///
/// Resolves the output header, opens the output writer, writes the BAM
/// header, builds the `GroupKeyConfig` and grouper, then defers to the
/// caller-supplied `build_fns` to construct the `PipelineFunctions`
/// (with or without an `mi_assign_fn` installed) before delegating to
/// [`run_bam_pipeline`].
///
/// `input_header` is taken by value so it can be cloned into the output
/// header when the caller passes `output_header: None`; both public
/// wrappers do the same.
#[allow(clippy::needless_pass_by_value)]
fn run_bam_pipeline_from_reader_inner<G, P, R, GrouperFn, BuildFns>(
    config: BamPipelineConfig,
    input: R,
    input_header: Header,
    output_path: &Path,
    output_header: Option<Header>,
    grouper_fn: GrouperFn,
    build_fns: BuildFns,
) -> io::Result<u64>
where
    G: Send + BatchWeight + MemoryEstimate + 'static,
    P: Send + MemoryEstimate + 'static,
    R: Read + Send + 'static,
    GrouperFn: FnOnce(&Header) -> Box<dyn Grouper<Group = G> + Send>,
    BuildFns: FnOnce(Header) -> PipelineFunctions<G, P>,
{
    // Use output_header if provided, otherwise clone input_header
    let output_header = output_header.unwrap_or_else(|| input_header.clone());

    // Create output writer (supports stdout via "-" or "/dev/stdout")
    let output_writer = open_pipeline_output(output_path)?;

    // Write BAM header using BGZF compression
    let mut header_writer = bam::io::Writer::new(output_writer);
    header_writer
        .write_header(&output_header)
        .map_err(|e| io::Error::other(format!("Failed to write BAM header: {e}")))?;

    // Finish the BGZF writer and get the underlying writer for the pipeline.
    let mut bgzf_writer = header_writer.into_inner();
    bgzf_writer
        .try_finish()
        .map_err(|e| io::Error::other(format!("Failed to finish BGZF header: {e}")))?;
    let output = bgzf_writer.into_inner();

    let output = BufWriter::with_capacity(IO_BUFFER_SIZE, output);

    // Pass the caller's `Option` straight through: `Some` computes keys during
    // decode, `None` skips it for groupers that never read the key. Every
    // key-consuming command sets `Some` explicitly, so no default is needed here.
    let group_key_config = config.group_key_config;

    // Create the grouper using INPUT header
    let grouper = grouper_fn(&input_header);

    // Hand the resolved OUTPUT header to the caller so its serialize closure
    // can capture it. This is the only piece that must stay caller-side
    // because the closure type is otherwise unnameable.
    let fns = build_fns(output_header);

    // Wrap in BufReader to reduce syscalls. The caller hands us the unbuffered
    // stream from `create_bam_reader_for_pipeline*`, and the pipeline's block
    // reader (`read_raw_blocks`) issues a one-byte probe, a 17-byte header read
    // and a body read per BGZF block — three syscalls per block straight off a
    // file or pipe without this. The by-path entry point buffers identically.
    let input = BufReader::with_capacity(IO_BUFFER_SIZE, input);

    run_bam_pipeline(
        config.pipeline,
        Box::new(input),
        Box::new(output),
        grouper,
        fns,
        group_key_config,
        None,
    )
}

/// Run a BAM pipeline from an already-opened reader.
///
/// This variant accepts a pre-opened reader and header, enabling streaming from
/// stdin or other non-seekable sources.
///
/// # Type Parameters
///
/// - `G`: Group type produced by the grouper (e.g., `RecordBuf`, `RawPositionGroup`, `MiGroup`)
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
///
/// # Errors
///
/// Returns an I/O error if any pipeline step or file I/O fails.
#[allow(clippy::too_many_arguments, clippy::needless_pass_by_value)]
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
    run_bam_pipeline_from_reader_inner(
        config,
        input,
        input_header,
        output_path,
        output_header,
        grouper_fn,
        |output_header| {
            PipelineFunctions::new(process_fn, move |p: P, buf: &mut Vec<u8>| {
                serialize_fn(p, &output_header, buf)
            })
        },
    )
}

/// Run a BAM pipeline from an already-opened reader, with an additional
/// serial-ordered MI-assign hook.
///
/// Same as [`run_bam_pipeline_from_reader`] but installs an `mi_assign_fn`
/// on the underlying [`PipelineFunctions`]. The hook is invoked once per
/// processed item, in pipeline-serial order
/// (`(batch_serial, idx_in_batch)`), before that item is handed to
/// `serialize_fn` — see
/// [`PipelineFunctions::mi_assign_fn`] for the contract. Used by
/// `fgumi group` and `fgumi dedup` to rewrite per-template `MoleculeId`s
/// from local to global so that `serialize_fn` itself stays
/// synchronization-free.
///
/// # Errors
///
/// Returns an I/O error if any pipeline step or file I/O fails.
#[allow(clippy::too_many_arguments, clippy::needless_pass_by_value)]
pub fn run_bam_pipeline_from_reader_with_mi_assign<
    G,
    P,
    R,
    GrouperFn,
    ProcessFn,
    SerializeFn,
    MiAssignFn,
>(
    config: BamPipelineConfig,
    input: R,
    input_header: Header,
    output_path: &Path,
    output_header: Option<Header>,
    grouper_fn: GrouperFn,
    process_fn: ProcessFn,
    serialize_fn: SerializeFn,
    mi_assign_fn: MiAssignFn,
) -> io::Result<u64>
where
    G: Send + BatchWeight + MemoryEstimate + 'static,
    P: Send + MemoryEstimate + 'static,
    R: Read + Send + 'static,
    GrouperFn: FnOnce(&Header) -> Box<dyn Grouper<Group = G> + Send>,
    ProcessFn: Fn(G) -> io::Result<P> + Send + Sync + 'static,
    SerializeFn: Fn(P, &Header, &mut Vec<u8>) -> io::Result<u64> + Send + Sync + 'static,
    MiAssignFn: Fn(BatchOrdinal, &mut P) -> io::Result<()> + Send + Sync + 'static,
{
    run_bam_pipeline_from_reader_inner(
        config,
        input,
        input_header,
        output_path,
        output_header,
        grouper_fn,
        |output_header| {
            PipelineFunctions::new(process_fn, move |p: P, buf: &mut Vec<u8>| {
                serialize_fn(p, &output_header, buf)
            })
            .with_mi_assign(mi_assign_fn)
        },
    )
}

/// Run a BAM pipeline from an already-opened reader, with a secondary output file.
///
/// Pick the header that the rejects (secondary) BAM should advertise.
///
/// `None` falls back to the resolved primary output header — that matches
/// `filter.rs`'s prior behavior, where rejects share the primary's header.
/// Consensus commands pass `Some(input_header)` so the rejects BAM carries
/// the raw-input header rather than the consensus output header.
#[inline]
fn resolve_secondary_header(
    secondary_output_header: Option<&Header>,
    primary_output_header: &Header,
) -> Header {
    secondary_output_header.cloned().unwrap_or_else(|| primary_output_header.clone())
}

/// This variant routes rejected/secondary records through the pipeline's ordering
/// infrastructure so both primary and secondary output files maintain input order.
///
/// The secondary serialize function is called with a borrow of the processed batch
/// BEFORE the primary serialize function consumes it.
///
/// # Parameters
///
/// - `output_header`: header for the primary output BAM. `None` reuses
///   `input_header`.
/// - `secondary_output_header`: header for the rejects (secondary) BAM. `None`
///   reuses the resolved primary `output_header` — the filter-style case where
///   primary and secondary share a header. Consensus commands pass
///   `Some(input_header)` so the rejects BAM advertises the raw-input header
///   (with the input's RG/PG/contig metadata and sort order) instead of the
///   transformed primary output header.
///
/// # Errors
///
/// Returns an I/O error if any pipeline step or file I/O fails.
#[allow(clippy::too_many_arguments, clippy::needless_pass_by_value)]
pub fn run_bam_pipeline_from_reader_with_secondary<
    G,
    P,
    R,
    GrouperFn,
    ProcessFn,
    SerializeFn,
    SecondaryFn,
>(
    config: BamPipelineConfig,
    input: R,
    input_header: Header,
    output_path: &Path,
    output_header: Option<Header>,
    secondary_output_path: &Path,
    secondary_output_header: Option<Header>,
    grouper_fn: GrouperFn,
    process_fn: ProcessFn,
    serialize_fn: SerializeFn,
    secondary_serialize_fn: SecondaryFn,
) -> io::Result<u64>
where
    G: Send + BatchWeight + MemoryEstimate + 'static,
    P: Send + MemoryEstimate + 'static,
    R: Read + Send + 'static,
    GrouperFn: FnOnce(&Header) -> Box<dyn Grouper<Group = G> + Send>,
    ProcessFn: Fn(G) -> io::Result<P> + Send + Sync + 'static,
    SerializeFn: Fn(P, &Header, &mut Vec<u8>) -> io::Result<u64> + Send + Sync + 'static,
    SecondaryFn: Fn(&P, &mut Vec<u8>) -> io::Result<u64> + Send + Sync + 'static,
{
    // Use output_header if provided, otherwise clone input_header
    let output_header = output_header.unwrap_or_else(|| input_header.clone());

    // Resolve which header the rejects (secondary) BAM should advertise.
    // Defaults to the primary output header to match prior behavior; consensus
    // commands override this with the input header so rejects carry raw-input
    // RG/PG/contig metadata and the input's claimed sort order (rejects are
    // emitted in input order, so the input sort order remains valid).
    let secondary_header =
        resolve_secondary_header(secondary_output_header.as_ref(), &output_header);

    // Create primary output BAM and write the output header
    let output_writer = open_pipeline_output(output_path)?;

    let mut header_writer = bam::io::Writer::new(output_writer);
    header_writer
        .write_header(&output_header)
        .map_err(|e| io::Error::other(format!("Failed to write BAM header: {e}")))?;

    let mut bgzf_writer = header_writer.into_inner();
    bgzf_writer
        .try_finish()
        .map_err(|e| io::Error::other(format!("Failed to finish BGZF header: {e}")))?;
    let output = bgzf_writer.into_inner();

    let output = BufWriter::with_capacity(IO_BUFFER_SIZE, output);

    // Create secondary output (reject BAM) with its own BGZF compression
    let secondary_writer = fgumi_bam_io::create_raw_bam_writer(
        secondary_output_path,
        &secondary_header,
        1, // single-threaded BGZF for secondary
        config.compression_level,
    )
    .map_err(|e| io::Error::other(format!("Failed to create secondary output: {e}")))?;

    // Pass the caller's `Option` straight through: `Some` computes keys during
    // decode, `None` skips it for groupers that never read the key. Every
    // key-consuming command sets `Some` explicitly, so no default is needed here.
    let group_key_config = config.group_key_config;

    let grouper = grouper_fn(&input_header);

    let fns = PipelineFunctions::new(process_fn, move |p: P, buf: &mut Vec<u8>| {
        serialize_fn(p, &output_header, buf)
    })
    .with_secondary_serialize(secondary_serialize_fn);

    let pipeline_config = config.pipeline;

    // Buffered for the same reason as the other from-reader entry point: the
    // caller's stream is unbuffered and `read_raw_blocks` reads block by block.
    let input = BufReader::with_capacity(IO_BUFFER_SIZE, input);

    run_bam_pipeline(
        pipeline_config,
        Box::new(input),
        Box::new(output),
        grouper,
        fns,
        group_key_config,
        Some(secondary_writer),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_info::LibraryIndex;
    use rstest::rstest;

    // ========================================================================
    // BatchOrdinal
    // ========================================================================

    #[rstest]
    // Final-item cases: idx_in_batch == batch_len - 1.
    #[case(0, 0, 1, true)]
    #[case(7, 4, 5, true)]
    // Non-final-item cases: idx_in_batch < batch_len - 1.
    #[case(0, 0, 4, false)]
    #[case(9, 2, 4, false)]
    fn batch_ordinal_is_last(
        #[case] batch_serial: u64,
        #[case] idx_in_batch: usize,
        #[case] batch_len: usize,
        #[case] expected: bool,
    ) {
        assert_eq!(BatchOrdinal { batch_serial, idx_in_batch, batch_len }.is_last(), expected);
    }

    #[test]
    fn pipeline_functions_default_has_no_optional_hooks() {
        // Sanity guard: adding `mi_assign_fn` (or any future optional field)
        // must not change defaults — every existing command relies on
        // `secondary_serialize_fn` and `mi_assign_fn` being `None` unless
        // explicitly installed.
        let fns: PipelineFunctions<(), ()> = PipelineFunctions::new(|_g| Ok(()), |_p, _buf| Ok(0));
        assert!(fns.secondary_serialize_fn.is_none());
        assert!(fns.mi_assign_fn.is_none());
    }

    #[test]
    fn pipeline_functions_with_mi_assign_installs_hook() {
        let fns: PipelineFunctions<(), ()> =
            PipelineFunctions::new(|_g| Ok(()), |_p, _buf| Ok(0)).with_mi_assign(|_ord, _p| Ok(()));
        assert!(fns.mi_assign_fn.is_some());
        // Other optional hooks remain unset.
        assert!(fns.secondary_serialize_fn.is_none());
    }

    /// Create a minimal `BamPipelineState` for testing memory backpressure.
    fn create_test_state(memory_limit: u64) -> BamPipelineState<(), ()> {
        let config = PipelineConfig::new(2, 6).with_queue_memory_limit(memory_limit);
        let input: Box<dyn Read + Send> = Box::new(std::io::empty());
        let output: Box<dyn Write + Send> = Box::new(std::io::sink());
        // Create minimal GroupKeyConfig for testing
        let header = Header::default();
        let library_index = LibraryIndex::from_header(&header);
        let group_key_config = GroupKeyConfig::new(library_index, SamTag::CB.into());
        BamPipelineState::new(config, input, output, Some(group_key_config))
    }

    /// A `queue_memory_limit` of 0 means "no limit", so Read is never gated.
    #[test]
    fn read_admission_is_unconditional_without_a_limit() {
        let state = create_test_state(0);
        state.q1_heap_bytes.store(u64::MAX / 2, Ordering::Release);
        assert!(state.read_admission_allowed());
    }

    /// Read is gated once the accounted queues reach the budget, and reopens
    /// as soon as a stage drains below it.
    #[test]
    fn read_admission_tracks_the_queue_memory_budget() {
        let state = create_test_state(1024);

        state.q2b_heap_bytes.store(512, Ordering::Release);
        assert!(state.read_admission_allowed(), "under budget: reading continues");

        state.q2b_heap_bytes.store(1024, Ordering::Release);
        assert!(!state.read_admission_allowed(), "at budget: reading stops");

        state.q2b_heap_bytes.store(256, Ordering::Release);
        assert!(state.read_admission_allowed(), "drained below budget: reading resumes");
    }

    /// The budget keeps discriminating above 512 MiB (issue #765).
    ///
    /// The per-stage high-water marks are capped at
    /// [`BACKPRESSURE_THRESHOLD_BYTES`](crate::unified_pipeline::BACKPRESSURE_THRESHOLD_BYTES),
    /// and for a long time those marks were the only thing `--max-memory`
    /// reached — which made every budget above 512 MiB behave identically. The
    /// budget's capacity role lives here instead, where nothing caps it: at
    /// 4 GiB in flight a 512 MiB budget must refuse and a 32 GiB budget must
    /// admit. If both answers ever agree again, the knob has gone inert.
    #[test]
    fn read_admission_still_discriminates_above_the_per_stage_marks() {
        let in_flight = 4 * 1024 * 1024 * 1024;

        let at_the_mark =
            create_test_state(crate::unified_pipeline::base::BACKPRESSURE_THRESHOLD_BYTES);
        at_the_mark.q2b_heap_bytes.store(in_flight, Ordering::Release);
        assert!(!at_the_mark.read_admission_allowed(), "512 MiB budget must gate 4 GiB in flight");

        let above_the_mark = create_test_state(32 * 1024 * 1024 * 1024);
        above_the_mark.q2b_heap_bytes.store(in_flight, Ordering::Release);
        assert!(
            above_the_mark.read_admission_allowed(),
            "32 GiB budget must admit 4 GiB in flight; a budget above 512 MiB is not inert"
        );
    }

    /// With nothing accounted for in flight, Read is always admitted — so an
    /// input whose first batch alone exceeds the whole budget still makes
    /// progress instead of wedging the pipeline.
    #[test]
    fn read_admission_always_allows_the_first_batch() {
        let state = create_test_state(1);
        assert_eq!(state.queue_bytes_in_flight(), 0);
        assert!(state.read_admission_allowed());
    }

    /// Q2b's charge is refunded exactly on pop, so the counter returns to zero.
    #[test]
    fn q2b_charge_is_refunded_on_pop() {
        let state = create_test_state(1024);
        assert_eq!(state.q2b_heap_bytes.load(Ordering::Acquire), 0);

        let batch = BoundaryBatch { buffer: vec![0u8; 4096], offsets: vec![0, 4096] };
        let charged = batch.estimate_heap_size() as u64;
        assert!(state.q2b_push(7, batch).is_ok());
        assert_eq!(state.q2b_heap_bytes.load(Ordering::Acquire), charged);
        assert!(charged >= 4096, "a 4 KiB payload must be charged as at least 4 KiB");

        let (serial, _batch) = state.q2b_pop().expect("pushed batch should pop");
        assert_eq!(serial, 7);
        assert_eq!(state.q2b_heap_bytes.load(Ordering::Acquire), 0);
    }

    /// Q5 (serialized) charges its byte counter *before* publishing and refunds
    /// it on pop. Regression for #810: charging after publication let a
    /// consumer's refund race ahead of the producer's charge and strand a
    /// phantom nonzero counter that keeps the `Read` gate shut.
    #[test]
    fn serialize_output_push_charges_before_publish_and_refunds_on_pop() {
        let state = create_test_state(0);
        let batch =
            SerializedBatch { data: vec![0u8; 4096], record_count: 0, secondary_data: None };
        let charge = batch.estimate_heap_size() as u64;
        assert!(charge >= 4096, "a 4 KiB payload must be charged as at least 4 KiB");

        assert!(state.serialize_output_push((3, batch)).is_ok());
        assert_eq!(state.output.serialized_heap_bytes.load(Ordering::Acquire), charge);

        let (serial, _batch) = state.q5_pop().expect("pushed batch should pop");
        assert_eq!(serial, 3);
        state.q5_track_pop(charge);
        assert_eq!(state.output.serialized_heap_bytes.load(Ordering::Acquire), 0);
    }

    /// Q6 (compressed) charges before publishing and refunds a push the queue
    /// rejects, so a full-queue rejection leaves no phantom charge behind — the
    /// same #810 ordering fix as the serialized queue.
    #[test]
    fn q6_push_charges_on_publish_and_refunds_a_rejected_push() {
        let state = create_test_state(0);
        let make = || CompressedBlockBatch {
            blocks: vec![],
            record_count: 0,
            secondary_data: Some(vec![0u8; 4096]),
        };
        let charge = make().estimate_heap_size() as u64;
        assert!(charge >= 4096, "a 4 KiB secondary payload must be charged as at least 4 KiB");

        // Fill the bounded output queue until a push is rejected.
        let mut accepted = 0u64;
        while state.q6_push((accepted, make())).is_ok() {
            accepted += 1;
            assert!(accepted < 100_000, "the bounded queue should fill well before this");
        }
        assert!(accepted > 0, "at least one push must be accepted");

        // The rejected push refunded its charge, so the counter reflects exactly
        // the accepted batches — not the one that bounced off a full queue.
        assert_eq!(
            state.output.compressed_heap_bytes.load(Ordering::Acquire),
            accepted * charge,
            "a rejected push must refund its own charge"
        );

        // Draining through the track-pop refund path returns the counter to zero.
        while state.q6_pop().is_some() {
            state.q6_track_pop(charge);
        }
        assert_eq!(state.output.compressed_heap_bytes.load(Ordering::Acquire), 0);
    }

    /// `serialize_input_pop` pops the processed (Q5-input) queue and returns its
    /// charge through the saturating `refund_queue_bytes` path (#810). The test
    /// state's processed payload is `Vec<()>`, whose heap charge is zero, so this
    /// pins the pop + refund *path* (the debit runs, saturating, and the popped
    /// batch is returned intact) rather than a nonzero byte amount.
    #[test]
    fn serialize_input_pop_refunds_the_processed_queue_charge() {
        let state = create_test_state(0);
        assert!(
            state.process_output_push((7, vec![(), (), ()])).is_ok(),
            "the empty processed queue must accept a push"
        );

        let (serial, batch) = state.serialize_input_pop().expect("pushed batch must pop");
        assert_eq!(serial, 7);
        assert_eq!(batch.len(), 3, "the popped batch must be returned intact");
        assert_eq!(
            state.output.processed_heap_bytes.load(Ordering::Acquire),
            0,
            "serialize_input_pop must leave the processed-queue counter floored at zero"
        );

        // Popping the now-empty queue is a no-op that touches no counter.
        assert!(state.serialize_input_pop().is_none());
    }

    /// `serialize_output_push` charges before publishing and refunds a push the
    /// bounded serialized queue rejects, so a full-queue rejection strands no
    /// phantom charge — the Q6 ordering fix (#810) applied to the serialized
    /// queue's rejected-push branch.
    #[test]
    fn serialize_output_push_refunds_a_rejected_push() {
        let state = create_test_state(0);
        let make =
            || SerializedBatch { data: vec![0u8; 4096], record_count: 0, secondary_data: None };
        let charge = make().estimate_heap_size() as u64;
        assert!(charge >= 4096, "a 4 KiB payload must be charged as at least 4 KiB");

        // Fill the bounded serialized queue until a push is rejected.
        let mut accepted = 0u64;
        while state.serialize_output_push((accepted, make())).is_ok() {
            accepted += 1;
            assert!(accepted < 100_000, "the bounded queue should fill well before this");
        }
        assert!(accepted > 0, "at least one push must be accepted");

        // The rejected push refunded its own charge, so the counter reflects
        // exactly the accepted batches — not the one that bounced off a full queue.
        assert_eq!(
            state.output.serialized_heap_bytes.load(Ordering::Acquire),
            accepted * charge,
            "a rejected serialized push must refund its own charge"
        );
    }

    /// `queue_bytes_in_flight` must include every accounted counter — not just
    /// `q2b_heap_bytes` — and the Read gate must react to their aggregate. A
    /// counter dropped from the sum, or a `q2b_push` that failed to charge,
    /// would slip past a test that only drives `q2b_heap_bytes` with a direct
    /// `store`.
    #[test]
    fn queue_bytes_in_flight_sums_charges_and_gates_on_the_aggregate() {
        // Fixed charges for two accounted counters, plus Q2b charged through
        // its real push path (not a direct store).
        const Q1_CHARGE: u64 = 256;
        const Q4_CHARGE: u64 = 512;
        let batch = BoundaryBatch { buffer: vec![0u8; 4096], offsets: vec![0, 4096] };
        let q2b_charge = batch.estimate_heap_size() as u64;

        // Budget sits above Q1 + Q4 alone but at the full aggregate, so the Q2b
        // push is what tips the gate shut.
        let limit = Q1_CHARGE + Q4_CHARGE + q2b_charge;
        let state = create_test_state(limit);

        state.q1_heap_bytes.store(Q1_CHARGE, Ordering::Release);
        state.output.groups_heap_bytes.store(Q4_CHARGE, Ordering::Release);
        assert_eq!(state.queue_bytes_in_flight(), Q1_CHARGE + Q4_CHARGE);
        assert!(state.read_admission_allowed(), "below budget with only Q1 and Q4 charged");

        assert!(state.q2b_push(3, batch).is_ok());
        assert_eq!(
            state.queue_bytes_in_flight(),
            Q1_CHARGE + Q4_CHARGE + q2b_charge,
            "aggregate must sum Q1, Q4 and the Q2b push charge"
        );
        assert!(!state.read_admission_allowed(), "at the aggregate budget: reading stops");

        // Draining Q2b through its pop path refunds only its own charge, and the
        // gate reopens because the aggregate falls back below the budget.
        state.q2b_pop().expect("pushed batch should pop");
        assert_eq!(state.queue_bytes_in_flight(), Q1_CHARGE + Q4_CHARGE);
        assert!(state.read_admission_allowed(), "Q2b drained: reading resumes");

        // Every remaining summed counter must contribute to the aggregate. Q1,
        // Q4 and Q2b are pinned above; charge the other six — Q2/Q3 reorder, Q5
        // processed, Q6 serialized, Q7 compressed, and the write-reorder state —
        // with distinct powers of two so dropping any single term from
        // `queue_bytes_in_flight` changes the total and fails here. The
        // write-reorder counter is the one the stalled-writer scenario this PR
        // fixes depends on most, so its omission would otherwise go uncaught.
        state.q2_reorder_state.add_heap_bytes(1);
        state.q3_reorder_state.add_heap_bytes(2);
        state.output.processed_heap_bytes.store(4, Ordering::Release);
        state.output.serialized_heap_bytes.store(8, Ordering::Release);
        state.output.compressed_heap_bytes.store(16, Ordering::Release);
        state.output.write_reorder_state.add_heap_bytes(32);
        assert_eq!(
            state.queue_bytes_in_flight(),
            Q1_CHARGE + Q4_CHARGE + 1 + 2 + 4 + 8 + 16 + 32,
            "every accounted counter must contribute to the aggregate"
        );
    }

    /// A `q2b_push` rejected because Q2b is out of slots must refund its charge
    /// and hand the batch back, leaving `q2b_heap_bytes` exactly where it was —
    /// otherwise a run of rejected pushes would permanently over-state the Read
    /// gate.
    #[test]
    fn q2b_push_refunds_and_returns_the_batch_when_out_of_slots() {
        let state = create_test_state(0);

        // Fill Q2b to capacity through the raw queue so its byte counter starts
        // at zero and only the rejected push under test can move it.
        let cap = state.q2b_boundaries.capacity();
        for serial in 0..cap as u64 {
            let filler = BoundaryBatch { buffer: Vec::new(), offsets: vec![0] };
            assert!(state.q2b_boundaries.push((serial, filler)).is_ok());
        }
        assert!(state.q2b_boundaries.is_full());
        assert_eq!(state.q2b_heap_bytes.load(Ordering::Acquire), 0);

        let rejected = BoundaryBatch { buffer: vec![0u8; 4096], offsets: vec![0, 4096] };
        let Err((serial, returned)) = state.q2b_push(99, rejected) else {
            panic!("push into a full Q2b must be rejected");
        };
        assert_eq!(serial, 99, "the rejected serial is handed back unchanged");
        assert_eq!(returned.buffer.len(), 4096, "the rejected batch is handed back unchanged");
        assert_eq!(
            state.q2b_heap_bytes.load(Ordering::Acquire),
            0,
            "a rejected push must leave the Q2b counter unchanged"
        );
    }

    /// The gap-filler serial (`serial == next_seq`) is admitted for decoding
    /// even when the Q3 reorder buffer is over its memory limit — otherwise the
    /// very serial Group is waiting on would be requeued and the pipeline would
    /// wedge (fgumi#746).
    #[test]
    fn admit_or_requeue_decode_admits_the_gap_filler_even_over_the_limit() {
        let state = create_test_state(1024);
        // next_seq = 0 and heap well over the 50% (512-byte) backpressure mark,
        // so `can_proceed` would reject any *future* serial here.
        state.q3_reorder_state.next_seq.store(0, Ordering::SeqCst);
        state.q3_reorder_state.heap_bytes.store(800, Ordering::SeqCst);

        let batch = BoundaryBatch { buffer: vec![0u8; 4096], offsets: vec![0, 4096] };
        let admitted = state.admit_or_requeue_decode(0, batch);
        let Some((serial, returned)) = admitted else {
            panic!("the gap-filler serial must be admitted for decoding");
        };
        assert_eq!(serial, 0, "the gap-filler serial is returned for decoding");
        assert_eq!(returned.buffer.len(), 4096, "the gap-filler batch is returned unchanged");
        assert!(state.q2b_boundaries.is_empty(), "the gap-filler is not requeued to Q2b");
    }

    /// A future serial over the backpressure mark is requeued to Q2b (not
    /// decoded) so another worker can still reach the gap-filler.
    #[test]
    fn admit_or_requeue_decode_requeues_a_future_serial_under_pressure() {
        let state = create_test_state(1024);
        state.q3_reorder_state.next_seq.store(0, Ordering::SeqCst);
        state.q3_reorder_state.heap_bytes.store(800, Ordering::SeqCst);

        let batch = BoundaryBatch { buffer: vec![0u8; 4096], offsets: vec![0, 4096] };
        assert!(
            state.admit_or_requeue_decode(50, batch).is_none(),
            "a future serial under memory pressure is requeued, not decoded"
        );
        assert_eq!(state.q2b_boundaries.len(), 1, "the future serial is handed back to Q2b");
        let (serial, _) = state.q2b_pop().expect("the requeued batch should pop");
        assert_eq!(serial, 50, "the requeued serial is preserved");
    }

    /// When a future serial cannot be requeued because Q2b is full, it is
    /// decoded directly — the bounded-overshoot path. This is safe because Q3
    /// reassembles by serial, so the batch is never lost or reordered.
    #[test]
    fn admit_or_requeue_decode_decodes_directly_when_q2b_is_full() {
        let state = create_test_state(1024);
        state.q3_reorder_state.next_seq.store(0, Ordering::SeqCst);
        state.q3_reorder_state.heap_bytes.store(800, Ordering::SeqCst);

        // Fill Q2b to capacity so the requeue push must fail.
        let cap = state.q2b_boundaries.capacity();
        for filler in 0..cap as u64 {
            let batch = BoundaryBatch { buffer: Vec::new(), offsets: vec![0] };
            assert!(state.q2b_boundaries.push((filler, batch)).is_ok());
        }
        assert!(state.q2b_boundaries.is_full());
        let charged_before = state.q2b_heap_bytes.load(Ordering::Acquire);

        let batch = BoundaryBatch { buffer: vec![0u8; 4096], offsets: vec![0, 4096] };
        let admitted = state.admit_or_requeue_decode(50, batch);
        let Some((serial, returned)) = admitted else {
            panic!("a future serial that cannot be requeued must be decoded directly");
        };
        assert_eq!(serial, 50, "the overshoot serial is returned for decoding");
        assert_eq!(returned.buffer.len(), 4096, "the overshoot batch is returned unchanged");
        assert_eq!(
            state.q2b_heap_bytes.load(Ordering::Acquire),
            charged_before,
            "the failed requeue refunds its charge, leaving the Q2b counter unchanged"
        );
    }

    /// A held raw batch that still cannot push at Priority 1 (Q1 full) must
    /// refund its Q1 charge, be handed back to the worker, and signal
    /// output-full — never leave bytes charged for a batch that is not in Q1.
    #[test]
    fn try_step_read_refunds_the_held_charge_when_q1_stays_full() {
        let state = create_test_state(0);

        // Fill Q1 so the held batch cannot be pushed.
        let cap = state.q1_raw_blocks.capacity();
        for serial in 0..cap as u64 {
            assert!(state.q1_raw_blocks.push((serial, RawBlockBatch::new())).is_ok());
        }
        assert!(state.q1_raw_blocks.is_full());
        assert_eq!(state.q1_heap_bytes.load(Ordering::Acquire), 0);

        // A held batch with real capacity so the (charge, refund) pair is
        // non-zero and a leak would be visible in the counter.
        let held = RawBlockBatch::with_capacity(8);
        assert!(held.estimate_heap_size() > 0, "the held batch must carry a non-zero charge");
        let mut worker = create_test_worker();
        worker.held_raw = Some((7, held));

        assert!(!try_step_read(&state, &mut worker), "a full Q1 signals output-full");
        assert!(worker.held_raw.is_some(), "the held batch is handed back to the worker");
        assert_eq!(worker.held_raw.as_ref().unwrap().0, 7, "the held serial is preserved");
        assert_eq!(
            state.q1_heap_bytes.load(Ordering::Acquire),
            0,
            "the rejected held push must refund its Q1 charge"
        );
    }

    /// With no held batch and Q1 physically open, Read must still decline to
    /// admit new input once the accounted in-flight bytes reach the queue
    /// memory budget — the Priority 3b gate.
    #[test]
    fn try_step_read_declines_when_over_the_queue_memory_budget() {
        let state = create_test_state(100);
        state.q1_heap_bytes.store(200, Ordering::Release);
        assert!(!state.read_admission_allowed(), "precondition: over budget");
        assert!(state.q1_raw_blocks.is_empty(), "precondition: Q1 has physical room");

        let mut worker = create_test_worker();
        assert!(worker.held_raw.is_none(), "no held batch, so Priority 1 is skipped");

        assert!(!try_step_read(&state, &mut worker), "over-budget Read declines to admit input");
        assert!(worker.held_raw.is_none(), "the gate reads nothing, so no batch is held");
        assert!(state.q1_raw_blocks.is_empty(), "the gate reads nothing into Q1");
        assert!(
            !state.read_done.load(Ordering::Acquire),
            "the gate must return before `read_raw_blocks`; a set `read_done` means Read hit EOF \
             instead of being declined"
        );
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
        // With no memory limit, Q3 admission is governed purely by the
        // serial-skew window (W = min(4 * num_threads, queue_capacity) = 8 for
        // this 2-thread test config); memory is bounded upstream by the Read
        // admission gate.
        let state = create_test_state(0);
        // next_seq (0) always proceeds.
        assert!(state.can_decode_proceed(0));
        // Serials within [next_seq, next_seq + W) proceed.
        assert!(state.can_decode_proceed(7));
        // Serials at/beyond the window are held back (bounded skew), even with
        // no memory limit.
        assert!(!state.can_decode_proceed(8));
        assert!(!state.can_decode_proceed(100));
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
        // Heap at the raw memory backstop, next_seq = 5.
        state.q3_reorder_state.heap_bytes.store(1024 * 1024, Ordering::SeqCst);
        state.q3_reorder_state.next_seq.store(5, Ordering::SeqCst);
        // The needed serial (next_seq) always proceeds, even over memory
        // (deadlock prevention).
        assert!(state.can_decode_proceed(5));
        // A future serial within the window is still held back by the
        // giant-batch memory backstop once heap reaches the limit.
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

    // ========================================================================
    // BoundaryState Tests
    // ========================================================================

    /// Frame `payload` as a BAM record: a little-endian `block_size` prefix
    /// followed by the payload bytes. `find_boundaries` only reads that prefix,
    /// so the payload's contents are irrelevant to boundary finding.
    fn framed_record(payload: &[u8]) -> Vec<u8> {
        let mut framed =
            u32::try_from(payload.len()).expect("payload fits in u32").to_le_bytes()[..].to_vec();
        framed.extend_from_slice(payload);
        framed
    }

    /// A record fragment too short to hold even a `block_size` must fail at EOF.
    ///
    /// `finish` scanned with `cursor + 4 <= leftover.len()`, so a 1-3 byte tail
    /// fell out of the loop and was dropped: with no complete record ahead of it
    /// the whole leftover was discarded with `Ok(None)`, and with one it was
    /// left past the final offset where nothing reads it. Either way the run
    /// ended successfully on a truncated file, which its own documentation
    /// already promised was an error.
    #[rstest]
    #[case::fragment_only(0)]
    #[case::record_then_fragment(1)]
    fn test_boundary_state_finish_rejects_a_sub_header_fragment(#[case] complete_records: usize) {
        let mut state = BoundaryState::new_no_header();
        let record = framed_record(&[b'A'; 32]);
        let mut data: Vec<u8> =
            std::iter::repeat_n(record.as_slice(), complete_records).flatten().copied().collect();
        // Three bytes: a truncated `block_size`, so its record length is unknown.
        data.extend_from_slice(&[0x01, 0x02, 0x03]);

        let batch = state.find_boundaries(&data).expect("complete records must be found");
        assert_eq!(
            batch.offsets.len() - 1,
            complete_records,
            "only the complete records may be emitted before EOF"
        );
        assert_eq!(
            batch.buffer,
            data[..data.len() - 3],
            "the complete records must come through byte-for-byte"
        );

        let err = state
            .finish()
            .expect_err("a record fragment at EOF must fail, as the doc comment promises");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        assert!(
            err.to_string().contains("3 trailing byte"),
            "the error must name the number of bytes left unread, got: {err}"
        );
    }

    /// A record split across two calls must be reassembled, not rejected.
    ///
    /// The control for the case above: leftover bytes are the normal state
    /// between blocks, so the EOF check must fire only at EOF and only on bytes
    /// that never became a record.
    #[test]
    fn test_boundary_state_reassembles_a_record_split_across_calls() {
        let mut state = BoundaryState::new_no_header();
        let record = framed_record(&[b'C'; 32]);
        let (head, tail) = record.split_at(20);

        let first = state.find_boundaries(head).expect("a partial record is held, not rejected");
        assert_eq!(first.offsets.len(), 1, "no complete record yet");

        let second = state.find_boundaries(tail).expect("the completing chunk must parse");
        assert_eq!(second.buffer, record, "the reassembled record must be byte-identical");

        assert!(
            state.finish().expect("nothing may remain after the record completed").is_none(),
            "a fully consumed stream leaves no final batch"
        );
    }

    // ========================================================================
    // Pipeline Validation Tests
    // ========================================================================

    /// A grouper that holds records forever, to drive completion validation.
    struct StuckGrouper;

    impl Grouper for StuckGrouper {
        type Group = ();

        fn add_records(&mut self, _records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>> {
            Ok(Vec::new())
        }

        fn finish(&mut self) -> io::Result<Option<Self::Group>> {
            Ok(None)
        }

        fn has_pending(&self) -> bool {
            true
        }
    }

    /// A grouper that is genuinely drained.
    struct DrainedGrouper;

    impl Grouper for DrainedGrouper {
        type Group = ();

        fn add_records(&mut self, _records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>> {
            Ok(Vec::new())
        }

        fn finish(&mut self) -> io::Result<Option<Self::Group>> {
            Ok(None)
        }

        fn has_pending(&self) -> bool {
            false
        }
    }

    /// Leftover bytes in the boundary finder are data loss and must be reported.
    ///
    /// `base.rs` names boundary leftovers as one of the three internal buffers
    /// completion validation exists to catch; the BAM implementation checked
    /// none of them. The FASTQ side already reports its own per-stream
    /// leftovers, so this brings the two into line.
    #[test]
    fn test_validation_detects_boundary_leftover() {
        let state = create_test_state(0);
        state.read_done.store(true, Ordering::SeqCst);
        state.group_done.store(true, Ordering::SeqCst);

        {
            let mut boundary = state.boundary_state.lock();
            *boundary = BoundaryState::new_no_header();
            // A record announcing 32 bytes but delivering 8: held as leftover.
            let mut partial = 32u32.to_le_bytes().to_vec();
            partial.extend_from_slice(&[b'A'; 8]);
            boundary.find_boundaries(&partial).expect("a partial record is held");
        }

        let err = state.validate_completion().expect_err("held bytes must fail validation");
        assert!(
            err.non_empty_queues.iter().any(|s| s.contains("boundary_leftover")),
            "validation must name the boundary leftover, got: {err}"
        );
    }

    /// The Group step's own buffers are the other two states `base.rs` names.
    ///
    /// The pipeline state cannot see them — the grouper and its pending groups
    /// live behind a separate mutex — so they are folded into the same
    /// validation result by the run function.
    #[test]
    fn test_group_state_reports_unflushed_records_and_groups() {
        let mut group_state = GroupState::new(Box::new(StuckGrouper));
        assert_eq!(
            group_state.unflushed_buffers(),
            vec!["grouper (partial group pending)".to_string()],
            "a grouper reporting pending records must be named"
        );

        group_state.pending_groups.push_back(());
        group_state.pending_groups.push_back(());
        assert_eq!(
            group_state.unflushed_buffers(),
            vec![
                "group_pending_groups (2)".to_string(),
                "grouper (partial group pending)".to_string(),
            ],
            "groups never pushed to Q4 must be named alongside the grouper"
        );

        let drained: GroupState<()> = GroupState::new(Box::new(DrainedGrouper));
        assert!(drained.unflushed_buffers().is_empty(), "a drained Group step must report nothing");
    }

    #[test]
    fn test_validation_passes_when_complete() {
        let state = create_test_state(0);
        // Set all done flags
        state.read_done.store(true, Ordering::SeqCst);
        state.group_done.store(true, Ordering::SeqCst);
        // Counters are all 0, queues are all empty
        let result = state.validate_completion();
        assert!(result.is_ok(), "Validation should pass: {result:?}");
    }

    #[test]
    fn test_validation_detects_non_empty_q1() {
        let state = create_test_state(0);
        state.read_done.store(true, Ordering::SeqCst);
        state.group_done.store(true, Ordering::SeqCst);

        // Add item to q1_raw_blocks
        let batch = RawBlockBatch { blocks: vec![] };
        assert!(state.q1_raw_blocks.push((0, batch)).is_ok());

        let result = state.validate_completion();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.non_empty_queues.iter().any(|s| s.contains("q1_raw_blocks")));
    }

    #[test]
    fn test_validation_detects_non_empty_q2() {
        let state = create_test_state(0);
        state.read_done.store(true, Ordering::SeqCst);
        state.group_done.store(true, Ordering::SeqCst);

        // Add item to q2_decompressed
        let batch = DecompressedBatch { data: vec![] };
        assert!(state.q2_decompressed.push((0, batch)).is_ok());

        let result = state.validate_completion();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.non_empty_queues.iter().any(|s| s.contains("q2_decompressed")));
    }

    /// `FindBoundaries` may only set `boundary_done` once every input batch has been
    /// consumed *and* every boundary batch it minted has actually been pushed.
    ///
    /// `batch_held` is the regression: a batch whose push failed is already counted in
    /// `batches_boundary_processed` while it waits in a worker's `held_boundaries`. See
    /// [`boundary_finding_is_complete`] for why that orphans the batch.
    #[rstest]
    #[case::all_pushed(true, 5, 5, true)]
    #[case::batch_held(true, 5, 4, false)]
    #[case::input_outstanding(true, 3, 5, false)]
    #[case::reader_still_going(false, 5, 5, false)]
    fn test_boundary_finding_completion_requires_no_held_batch(
        #[case] read_done: bool,
        #[case] batches_processed: u64,
        #[case] batches_pushed: u64,
        #[case] expected_complete: bool,
    ) {
        let state = create_test_state(0);
        state.read_done.store(read_done, Ordering::SeqCst);
        state.next_read_serial.store(5, Ordering::SeqCst);
        state.batches_boundary_processed.store(batches_processed, Ordering::SeqCst);
        // Five batches minted a serial; `batches_pushed` of them reached the queue, so
        // a shortfall is a batch held for a retried push.
        state.next_boundary_serial.store(5, Ordering::SeqCst);
        state.batches_boundary_found.store(batches_pushed, Ordering::SeqCst);

        assert_eq!(
            boundary_finding_is_complete(&state),
            expected_complete,
            "completion must require all input consumed and every minted batch pushed"
        );
    }

    /// An empty input completes immediately: every counter is zero, so both equalities
    /// hold trivially. Worth pinning because the failure mode is a stage that never
    /// declares itself done and hangs on a file with no records.
    #[test]
    fn test_boundary_finding_completes_on_empty_input() {
        let state = create_test_state(0);
        state.read_done.store(true, Ordering::SeqCst);

        assert!(
            boundary_finding_is_complete(&state),
            "an input with no batches must complete rather than hang"
        );
    }

    #[test]
    fn test_validation_detects_counter_mismatch_decompressed() {
        let state = create_test_state(0);
        state.read_done.store(true, Ordering::SeqCst);
        state.group_done.store(true, Ordering::SeqCst);

        // Simulate: read 5 batches, but only decompressed 3
        state.next_read_serial.store(5, Ordering::SeqCst);
        state.batches_decompressed.store(3, Ordering::SeqCst);

        let result = state.validate_completion();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.counter_mismatches.iter().any(|s| s.contains("batches_decompressed")));
    }

    #[test]
    fn test_validation_detects_counter_mismatch_boundary_processed() {
        let state = create_test_state(0);
        state.read_done.store(true, Ordering::SeqCst);
        state.group_done.store(true, Ordering::SeqCst);

        // All decompressed but not all boundary processed
        state.next_read_serial.store(5, Ordering::SeqCst);
        state.batches_decompressed.store(5, Ordering::SeqCst);
        state.batches_boundary_processed.store(3, Ordering::SeqCst);

        let result = state.validate_completion();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.counter_mismatches.iter().any(|s| s.contains("batches_boundary_processed")));
    }

    #[test]
    fn test_validation_detects_counter_mismatch_grouped() {
        let state = create_test_state(0);
        state.read_done.store(true, Ordering::SeqCst);
        state.group_done.store(true, Ordering::SeqCst);

        // Everything processed up to group, but group didn't finish
        state.next_read_serial.store(5, Ordering::SeqCst);
        state.batches_decompressed.store(5, Ordering::SeqCst);
        state.batches_boundary_processed.store(5, Ordering::SeqCst);
        state.batches_boundary_found.store(5, Ordering::SeqCst);
        state.batches_decoded.store(5, Ordering::SeqCst);
        state.batches_grouped.store(3, Ordering::SeqCst);

        let result = state.validate_completion();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.counter_mismatches.iter().any(|s| s.contains("batches_grouped")));
    }

    #[test]
    fn test_validation_error_display() {
        let err = PipelineValidationError {
            non_empty_queues: vec!["q1 (5)".to_string(), "q2 (3)".to_string()],
            counter_mismatches: vec!["batches_x (5) != batches_y (3)".to_string()],
            leaked_heap_bytes: 0,
        };
        let display = err.to_string();
        assert!(display.contains("q1"));
        assert!(display.contains("q2"));
        assert!(display.contains("batches_x"));
    }

    #[test]
    fn test_validation_detects_non_empty_reorder_buffer() {
        let state = create_test_state(0);
        state.read_done.store(true, Ordering::SeqCst);
        state.group_done.store(true, Ordering::SeqCst);

        // Add item to q2_reorder buffer
        {
            let mut q2_reorder = state.q2_reorder.lock();
            let batch = DecompressedBatch { data: vec![] };
            q2_reorder.insert(0, batch);
        }

        let result = state.validate_completion();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.non_empty_queues.iter().any(|s| s.contains("q2_reorder")));
    }

    #[test]
    fn test_validation_detects_non_empty_q3_reorder() {
        let state = create_test_state(0);
        state.read_done.store(true, Ordering::SeqCst);
        state.group_done.store(true, Ordering::SeqCst);

        // Add item to q3_reorder buffer
        {
            let mut q3_reorder = state.q3_reorder.lock();
            let batch: Vec<DecodedRecord> = vec![];
            q3_reorder.insert(0, batch);
        }

        let result = state.validate_completion();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.non_empty_queues.iter().any(|s| s.contains("q3_reorder")));
    }

    #[test]
    fn test_validation_detects_non_empty_output_queue() {
        let state = create_test_state(0);
        state.read_done.store(true, Ordering::SeqCst);
        state.group_done.store(true, Ordering::SeqCst);

        // Add item to output.groups (q4_groups)
        // G = () for create_test_state
        let batch: Vec<()> = vec![()];
        assert!(state.output.groups.push((0, batch)).is_ok());

        let result = state.validate_completion();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.non_empty_queues.iter().any(|s| s.contains("q4_groups")));
    }

    #[test]
    fn test_validation_detects_non_empty_write_reorder() {
        let state = create_test_state(0);
        state.read_done.store(true, Ordering::SeqCst);
        state.group_done.store(true, Ordering::SeqCst);

        // Add item to write_reorder buffer
        {
            let mut write_reorder = state.output.write_reorder.lock();
            let batch =
                CompressedBlockBatch { blocks: vec![], record_count: 0, secondary_data: None };
            write_reorder.insert(0, batch);
        }

        let result = state.validate_completion();
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.non_empty_queues.iter().any(|s| s.contains("write_reorder")));
    }

    /// Create a worker for testing decompress/decode steps
    /// (`compression_level=6`, `thread_id=0`, `num_threads=2`).
    fn create_test_worker() -> WorkerState<()> {
        WorkerState::new(6, 0, 2, SchedulerStrategy::default())
    }

    /// `Decode` finding `q2b_boundaries` empty always counts as Q2b starvation, but
    /// only counts as Q2 starvation while `FindBoundaries` might still deliver more.
    /// Once it is done and the queue has drained, an empty Q2b is the terminal state.
    ///
    /// This is the condition that used to sit in the `else` arm of the `decode_done`
    /// store; the flag is gone, so the behaviour is pinned here instead.
    #[rstest]
    #[case::boundaries_may_still_arrive(false, 1)]
    #[case::boundary_finding_done(true, 0)]
    fn test_decode_records_q2_starvation_only_while_boundaries_may_arrive(
        #[case] boundary_done: bool,
        #[case] expected_q2_empty: u64,
    ) {
        let config = PipelineConfig::new(2, 6).with_stats(true);
        let input: Box<dyn Read + Send> = Box::new(std::io::empty());
        let output: Box<dyn Write + Send> = Box::new(std::io::sink());
        let header = Header::default();
        let library_index = LibraryIndex::from_header(&header);
        let group_key_config = GroupKeyConfig::new(library_index, SamTag::CB.into());
        let state: BamPipelineState<(), ()> =
            BamPipelineState::new(config, input, output, Some(group_key_config));
        let mut worker = create_test_worker();

        state.boundary_done.store(boundary_done, Ordering::SeqCst);

        // Nothing in q2b, nothing held: Decode falls through to the starvation check.
        assert!(!try_step_decode(&state, &mut worker), "no input means no work done");

        let stats = state.stats().expect("stats were enabled");
        assert_eq!(
            stats.q2b_empty.load(Ordering::Relaxed),
            1,
            "an empty Q2b is always recorded as Q2b starvation"
        );
        assert_eq!(
            stats.q2_empty.load(Ordering::Relaxed),
            expected_q2_empty,
            "Q2 starvation must be recorded only while more boundaries may arrive"
        );
    }

    /// Seed a reorder buffer state with a held-item stress scenario: a
    /// non-next_seq serial is waiting to land while the reorder heap is
    /// under partial memory pressure. `next_seq = 0` (so any held serial
    /// other than 0 is "out of order"), and `heap_bytes = 800` — above the
    /// 50% mark of the `memory_limit = 1024` used by these tests but below
    /// the full `is_memory_high` threshold, so the P2 memory gate does NOT
    /// fire. This lets the held-item P1 unconditional-push path be exercised
    /// in isolation: Priority 1 must admit the out-of-order held batch, and
    /// Priority 2 must still allow Q1/Q2b to drain (because memory is not
    /// high).
    fn setup_memory_backpressure(reorder_state: &ReorderBufferState) {
        reorder_state.heap_bytes.store(800, Ordering::SeqCst);
        reorder_state.next_seq.store(0, Ordering::SeqCst);
    }

    /// Held batch pushes unconditionally at Priority 1 when Q2 has physical
    /// capacity, regardless of memory backpressure. This prevents the deadlock
    /// where all workers hold non-next_seq batches and nobody can produce
    /// `next_seq`.
    #[test]
    fn test_decompress_held_pushes_unconditionally_when_q2_has_room() {
        let state = create_test_state(1024);
        setup_memory_backpressure(&state.q2_reorder_state);

        let raw = RawBlockBatch::new();
        assert!(state.q1_raw_blocks.push((0, raw)).is_ok());

        // Worker holds serial 50 — would be blocked by can_proceed (50 != next_seq 0),
        // but Priority 1 now pushes unconditionally.
        let mut worker = create_test_worker();
        worker.held_decompressed = Some((50, DecompressedBatch { data: vec![0u8; 16] }, 16));

        let result = try_step_decompress(&state, &mut worker);
        assert!(result, "should succeed — held batch pushed, then new batch processed");
        assert!(!state.has_error(), "should not error");
        assert!(state.q1_raw_blocks.is_empty(), "Q1 should have been popped");
        // Both the held batch (serial 50) and the new batch were pushed to Q2
        assert_eq!(state.q2_decompressed.len(), 2, "Q2 should have both batches");
        assert!(worker.held_decompressed.is_none(), "held slot should be empty");
    }

    /// When Q2 is physically full, the held batch cannot push at Priority 1.
    /// The function returns false without error and without popping Q1.
    #[test]
    fn test_decompress_held_blocked_by_full_q2() {
        let state = create_test_state(1024);

        let cap = state.q2_decompressed.capacity();
        for i in 0..cap {
            assert!(
                state
                    .q2_decompressed
                    .push((i as u64, DecompressedBatch { data: vec![0u8; 8] }))
                    .is_ok(),
                "failed to fill q2 at serial {i}"
            );
        }
        assert!(state.q2_decompressed.is_full());

        assert!(state.q1_raw_blocks.push((100, RawBlockBatch::new())).is_ok());

        let mut worker = create_test_worker();
        worker.held_decompressed = Some((50, DecompressedBatch { data: vec![0u8; 16] }, 16));

        let result = try_step_decompress(&state, &mut worker);
        assert!(!result, "should return false when Q2 is full");
        assert!(!state.has_error(), "should NOT set error; got: {:?}", state.take_error());
        assert!(worker.held_decompressed.is_some(), "held batch should be preserved");
        assert!(!state.q1_raw_blocks.is_empty(), "Q1 should not have been popped");
    }

    /// Held batch pushes unconditionally at Priority 1 when Q3 has physical
    /// capacity (symmetric to the decompress test).
    #[test]
    fn test_decode_held_pushes_unconditionally_when_q3_has_room() {
        let state = create_test_state(1024);
        setup_memory_backpressure(&state.q3_reorder_state);

        let boundary = BoundaryBatch { buffer: Vec::new(), offsets: vec![0] };
        assert!(state.q2b_push(0, boundary).is_ok());

        let mut worker = create_test_worker();
        worker.held_decoded = Some((50, vec![], 16));

        let result = try_step_decode(&state, &mut worker);
        assert!(result, "should succeed — held batch pushed, then new batch processed");
        assert!(!state.has_error(), "should not error");
        assert!(state.q2b_boundaries.is_empty(), "Q2b should have been popped");
        // Both the held batch and the new batch were pushed to Q3
        assert_eq!(state.q3_decoded.len(), 2, "Q3 should have both batches");
        assert!(worker.held_decoded.is_none(), "held slot should be empty");
    }

    /// When Q3 is physically full, the held batch cannot push at Priority 1.
    /// The function returns false without error.
    #[test]
    fn test_decode_held_blocked_by_full_q3() {
        let state = create_test_state(1024);

        let cap = state.q3_decoded.capacity();
        for i in 0..cap {
            assert!(
                state.q3_decoded.push((i as u64, vec![])).is_ok(),
                "failed to fill q3 at serial {i}"
            );
        }
        assert!(state.q3_decoded.is_full());

        let boundary = BoundaryBatch { buffer: Vec::new(), offsets: vec![0] };
        assert!(state.q2b_push(100, boundary).is_ok());

        let mut worker = create_test_worker();
        worker.held_decoded = Some((50, vec![], 16));

        let result = try_step_decode(&state, &mut worker);
        assert!(!result, "should return false when Q3 is full");
        assert!(!state.has_error(), "should NOT set error; got: {:?}", state.take_error());
        assert!(worker.held_decoded.is_some(), "held batch should be preserved");
        assert!(!state.q2b_boundaries.is_empty(), "Q2b should not have been popped");
    }

    #[rstest]
    #[case::none_falls_back_to_primary(None, "primary")]
    #[case::some_overrides_primary(Some("rejects"), "rejects")]
    fn test_resolve_secondary_header(
        #[case] override_comment: Option<&str>,
        #[case] expected_comment: &str,
    ) {
        let primary = Header::builder().add_comment("primary").build();
        let override_header = override_comment.map(|c| Header::builder().add_comment(c).build());
        let resolved = resolve_secondary_header(override_header.as_ref(), &primary);
        let expected = Header::builder().add_comment(expected_comment).build();
        assert_eq!(resolved, expected);
    }

    // ========================================================================
    // UMI position caching during Decode (issue #334)
    // ========================================================================

    /// Build a raw BAM record carrying an RX UMI tag for cache tests.
    fn raw_record_with_rx(umi: &[u8]) -> fgumi_raw_bam::RawRecord {
        use fgumi_raw_bam::{SamBuilder as RawSamBuilder, flags as raw_flags};

        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGT")
            .qualities(&[30u8; 4])
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(0)
            .mapq(30)
            .cigar_ops(&[4u32 << 4]);
        b.add_string_tag(SamTag::RX, umi);
        b.build()
    }

    #[test]
    fn test_decode_records_populates_umi_cache_when_configured() {
        // Build a BoundaryBatch holding one record (with the 4-byte block_size
        // prefix BAM records carry on disk) and decode it with a GroupKeyConfig
        // that asks for RX caching.
        let raw = raw_record_with_rx(b"CACHEDUMI");
        let mut buffer = Vec::with_capacity(raw.as_ref().len() + 4);
        let block_size = u32::try_from(raw.as_ref().len()).unwrap();
        buffer.extend_from_slice(&block_size.to_le_bytes());
        buffer.extend_from_slice(raw.as_ref());
        let offsets = vec![0usize, buffer.len()];
        let batch = BoundaryBatch { buffer, offsets };

        let header = Header::default();
        let library_index = LibraryIndex::from_header(&header);
        let cfg = GroupKeyConfig::new(library_index, SamTag::CB.into()).with_umi_tag(*SamTag::RX);

        let decoded = decode_records(&batch, Some(&cfg)).expect("decode succeeds");
        assert_eq!(decoded.len(), 1);
        assert_eq!(decoded[0].cached_umi(), Some(b"CACHEDUMI".as_ref()));
    }

    #[test]
    fn test_decode_records_leaves_umi_cache_unset_when_disabled() {
        let raw = raw_record_with_rx(b"CACHEDUMI");
        let mut buffer = Vec::with_capacity(raw.as_ref().len() + 4);
        let block_size = u32::try_from(raw.as_ref().len()).unwrap();
        buffer.extend_from_slice(&block_size.to_le_bytes());
        buffer.extend_from_slice(raw.as_ref());
        let offsets = vec![0usize, buffer.len()];
        let batch = BoundaryBatch { buffer, offsets };

        let header = Header::default();
        let library_index = LibraryIndex::from_header(&header);
        // No `with_umi_tag` → caching is disabled.
        let cfg = GroupKeyConfig::new(library_index, SamTag::CB.into());

        let decoded = decode_records(&batch, Some(&cfg)).expect("decode succeeds");
        assert_eq!(decoded.len(), 1);
        assert!(decoded[0].cached_umi().is_none());
    }

    #[test]
    fn test_decode_records_skips_group_key_when_none() {
        // With `None`, per-record GroupKey computation is skipped and a default
        // placeholder key is attached — for groupers (e.g. SingleRawRecordGrouper)
        // that never read it. The raw record bytes must be byte-identical to the
        // key-computing path so downstream passthrough is unaffected.
        let raw = raw_record_with_rx(b"CACHEDUMI");
        let make_batch = || {
            let mut buffer = Vec::with_capacity(raw.as_ref().len() + 4);
            let block_size = u32::try_from(raw.as_ref().len()).unwrap();
            buffer.extend_from_slice(&block_size.to_le_bytes());
            buffer.extend_from_slice(raw.as_ref());
            let offsets = vec![0usize, buffer.len()];
            BoundaryBatch { buffer, offsets }
        };

        let header = Header::default();
        let library_index = LibraryIndex::from_header(&header);
        let cfg = GroupKeyConfig::new(library_index, SamTag::CB.into());

        let with_key = decode_records(&make_batch(), Some(&cfg)).expect("decode Some");
        let without_key = decode_records(&make_batch(), None).expect("decode None");

        assert_eq!(with_key.len(), 1);
        assert_eq!(without_key.len(), 1);
        // Identical record bytes either way — the transform is pure passthrough.
        assert_eq!(with_key[0].raw_bytes(), without_key[0].raw_bytes());
        // `None` attaches the default placeholder; `Some` computes a real key that
        // differs from it, proving the computation was actually skipped.
        assert_eq!(without_key[0].key, crate::unified_pipeline::base::GroupKey::default());
        assert_ne!(with_key[0].key, without_key[0].key, "Some path computes a real key");
    }

    /// A secondary/supplementary read carries the exact template coordinate of its
    /// primary pair in the `tc` tag (stamped by `fgumi zipper`). dedup must key such
    /// a read into the SAME position group as its primary — rather than an UNKNOWN
    /// position that relies on the read happening to sort adjacent to its primary.
    ///
    /// This asserts the stronger property: the tc-derived key is byte-for-byte equal
    /// to the *real* primary pair's group key (computed independently from the primary
    /// reads), not merely that its fields echo the raw tc values. Both the SECONDARY
    /// and SUPPLEMENTARY flags take the tc branch, so both are exercised. Only fgumi
    /// can do this, because only fgumi persists `tc`.
    #[rstest]
    #[case::secondary(fgumi_raw_bam::flags::SECONDARY)]
    #[case::supplementary(fgumi_raw_bam::flags::SUPPLEMENTARY)]
    fn test_group_key_secondary_supplementary_matches_primary_via_tc(#[case] sec_supp_flag: u16) {
        use crate::unified_pipeline::GroupKey;
        use fgumi_raw_bam::{SamBuilder as RawSamBuilder, flags as raw_flags};

        let header = Header::default();
        let library_index = LibraryIndex::from_header(&header);

        // Build the REAL R1 primary of the pair: forward, ref 0, mate reverse and
        // downstream (its unclipped 5' resolved from the MC tag). Its group key is
        // the template's true position group — the target the sec/supp must join.
        let mut r1 = RawSamBuilder::new();
        r1.read_name(b"tmpl")
            .sequence(b"ACGTACGT")
            .qualities(&[30u8; 8])
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::MATE_REVERSE)
            .ref_id(0)
            .pos(999)
            .mapq(60)
            .cigar_ops(&[8u32 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(1392);
        r1.add_string_tag(SamTag::MC, b"8M"); // mate cigar → mate unclipped 5'
        let r1 = r1.build();
        let (primary_key, _) = compute_group_key_from_raw(r1.as_ref(), &library_index, None, None);

        // Sanity: the primary produced a real two-lane paired key, not a fallback.
        assert_ne!(primary_key.pos2, GroupKey::UNKNOWN_POS, "primary must be a paired key");
        assert_ne!(primary_key.pos1, primary_key.pos2, "the two lanes differ in position");

        // Stamp the primary's template coordinate into the 6-element `tc` array that
        // zipper writes, and hang it on a sec/supp read whose OWN alignment is far
        // away (pos 9000) — exactly the case where per-record keying misplaces it at
        // its own/mate coordinate instead of the template's.
        let tc = [
            primary_key.ref_id1,
            primary_key.pos1,
            i32::from(primary_key.strand1),
            primary_key.ref_id2,
            primary_key.pos2,
            i32::from(primary_key.strand2),
        ];
        let mut ss = RawSamBuilder::new();
        ss.read_name(b"tmpl")
            .sequence(b"ACGTACGT")
            .qualities(&[30u8; 8])
            .flags(raw_flags::PAIRED | sec_supp_flag)
            .ref_id(0)
            .pos(9000)
            .mapq(60)
            .cigar_ops(&[8u32 << 4])
            .mate_ref_id(0)
            .mate_pos(1392);
        ss.add_array_i32(SamTag::TC, &tc);
        let ss = ss.build();
        let (ss_key, _) = compute_group_key_from_raw(ss.as_ref(), &library_index, None, None);

        assert_eq!(
            ss_key, primary_key,
            "sec/supp read keyed via tc must group IDENTICALLY to its primary pair, \
             not at its own (pos=9000) alignment"
        );
    }

    #[test]
    fn test_pipeline_functions_secondary_serialize() {
        let fns = PipelineFunctions::<Vec<u8>, Vec<u8>>::new(Ok, |data, buf| {
            buf.extend_from_slice(&data);
            Ok(1)
        });
        assert!(fns.secondary_serialize_fn.is_none());

        let fns = fns.with_secondary_serialize(|data: &Vec<u8>, buf: &mut Vec<u8>| {
            buf.extend_from_slice(data);
            Ok(1)
        });
        assert!(fns.secondary_serialize_fn.is_some());

        // Verify the secondary serialize function works
        let test_data = vec![1u8, 2, 3, 4];
        let mut buf = Vec::new();
        let count =
            (fns.secondary_serialize_fn.as_ref().expect("secondary_serialize_fn should be set"))(
                &test_data, &mut buf,
            )
            .expect("serialize should succeed");
        assert_eq!(count, 1);
        assert_eq!(buf, vec![1, 2, 3, 4]);
    }
}
