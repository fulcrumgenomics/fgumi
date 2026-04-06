//! Inline record buffer for samtools-style memory layout.
//!
//! This module provides a contiguous buffer that stores BAM records inline
//! with pre-computed sort keys, eliminating per-record heap allocations.
//!
//! Key benefits over `Vec<(Key, Record)>`:
//! - Single large allocation instead of millions of small ones
//! - Better cache locality - sequential memory access
//! - Sort by index array, not by moving actual record data
//! - ~24 bytes overhead per record vs ~110 bytes

use crate::sort::bam_fields;
use crate::sort::keys::{RawCoordinateKey, RawSortKey, SortContext};
use crate::sort::radix::bytes_needed_u64;
use crate::sort::segmented_buf::SegmentedBuf;
use std::cmp::Ordering;
use std::io::{Read, Write};

// ============================================================================
// Packed Sort Keys
// ============================================================================

/// Packed sort key for coordinate ordering.
///
/// Format: `(tid << 34) | ((pos+1) << 1) | reverse`
///
/// This allows single u64 comparison for most records. The +1 on pos
/// ensures pos=0 doesn't collide with unmapped (which uses MAX).
#[repr(transparent)]
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct PackedCoordinateKey(pub u64);

impl PackedCoordinateKey {
    /// Create a packed coordinate key.
    ///
    /// # Arguments
    /// * `tid` - Reference sequence ID (-1 for unmapped)
    /// * `pos` - 0-based alignment position
    /// * `reverse` - True if reverse complemented
    /// * `nref` - Number of reference sequences (for unmapped handling)
    #[inline]
    #[must_use]
    #[allow(clippy::cast_sign_loss)]
    pub fn new(tid: i32, pos: i32, reverse: bool, _nref: u32) -> Self {
        if tid < 0 {
            return Self::unmapped();
        }
        let tid = tid as u32;
        // Pack: tid in high bits, (pos+1) in middle, reverse in LSB
        // Using pos+1 so that pos=0 doesn't become 0 in the key
        #[allow(clippy::cast_lossless)] // Explicit bit packing requires precise control
        let key = (u64::from(tid) << 34)
            | ((i64::from(pos) as u64).wrapping_add(1) << 1)
            | u64::from(reverse);
        Self(key)
    }

    /// Create a key for unmapped records (sorts after all mapped).
    #[inline]
    #[must_use]
    pub fn unmapped() -> Self {
        Self(u64::MAX)
    }
}

// ============================================================================
// Record References (Index for Sorting)
// ============================================================================

/// Reference to a record in the buffer (used for sorting).
///
/// This is a lightweight handle that can be sorted efficiently.
/// The actual record data stays in place in the buffer.
///
/// Note: No read name tie-breaking is used, matching samtools behavior.
/// Equal records maintain their original input order (stable sort).
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct RecordRef {
    /// Packed primary sort key for fast comparison.
    pub sort_key: u64,
    /// Offset into `RecordBuffer` where record header starts.
    pub offset: u64,
    /// Length of raw BAM data (excluding inline header).
    pub len: u32,
    /// Padding for 8-byte alignment.
    padding: u32,
}

impl PartialEq for RecordRef {
    fn eq(&self, other: &Self) -> bool {
        self.sort_key == other.sort_key
    }
}

impl Eq for RecordRef {}

impl PartialOrd for RecordRef {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RecordRef {
    fn cmp(&self, other: &Self) -> Ordering {
        self.sort_key.cmp(&other.sort_key)
    }
}

// ============================================================================
// Inline Header (stored before each record in buffer)
// ============================================================================

/// Inline record header stored in buffer before raw BAM data.
///
/// This header is written once when the record is added and contains
/// pre-computed sort keys for efficient sorting.
#[repr(C)]
#[derive(Copy, Clone, Debug)]
struct InlineHeader {
    /// Pre-computed packed sort key.
    sort_key: u64,
    /// Length of following raw BAM data.
    record_len: u32,
    /// Padding for 8-byte alignment.
    padding: u32,
}

const HEADER_SIZE: usize = std::mem::size_of::<InlineHeader>(); // 16 bytes

// ============================================================================
// RecordBuffer - Main Data Structure
// ============================================================================

/// Contiguous buffer for inline record storage.
///
/// Records are stored sequentially in a single allocation:
/// ```text
/// [Header0][BAM0][Header1][BAM1][Header2][BAM2]...
/// ```
///
/// An index array of `RecordRef` is maintained for sorting.
/// Sorting only reorders the index; records stay in place.
pub struct RecordBuffer {
    /// Segmented byte storage for all records (headers + BAM data).
    data: SegmentedBuf,
    /// Index of record references for sorting.
    refs: Vec<RecordRef>,
    /// Number of reference sequences (for unmapped handling).
    nref: u32,
}

/// Segment size for `RecordBuffer`'s `SegmentedBuf`: 256 MiB.
const RECORD_SEGMENT_SIZE: usize = 256 * 1024 * 1024;

impl RecordBuffer {
    /// Create a new buffer with estimated capacity.
    ///
    /// # Arguments
    /// * `estimated_records` - Expected number of records
    /// * `estimated_bytes` - Expected total bytes of BAM data
    /// * `nref` - Number of reference sequences in header
    #[must_use]
    pub fn with_capacity(estimated_records: usize, estimated_bytes: usize, nref: u32) -> Self {
        Self {
            data: SegmentedBuf::with_capacity(
                estimated_bytes + estimated_records * HEADER_SIZE,
                RECORD_SEGMENT_SIZE,
            ),
            refs: Vec::with_capacity(estimated_records),
            nref,
        }
    }

    /// Push a record for coordinate sorting.
    ///
    /// Extracts the sort key inline from raw BAM bytes (zero-copy).
    ///
    /// # Errors
    ///
    /// Returns an error if the record (plus header) exceeds the segment size
    /// (256 MiB) or if the record length exceeds `u32::MAX`.
    #[inline]
    pub fn push_coordinate(&mut self, record: &[u8]) -> anyhow::Result<()> {
        // BAM fixed-length block is 32 bytes; coordinate fields end at offset 15.
        const MIN_BAM_RECORD_LEN: usize = 16;
        anyhow::ensure!(
            record.len() >= MIN_BAM_RECORD_LEN,
            "BAM record is truncated: need at least {} bytes to extract coordinate fields, got {}",
            MIN_BAM_RECORD_LEN,
            record.len(),
        );

        let total_bytes = HEADER_SIZE + record.len();
        anyhow::ensure!(
            total_bytes <= RECORD_SEGMENT_SIZE,
            "BAM record of {} bytes (+ {} byte header) exceeds segment size of {} bytes; \
             this is likely a malformed BAM file",
            record.len(),
            HEADER_SIZE,
            RECORD_SEGMENT_SIZE,
        );
        let len = u32::try_from(record.len())
            .map_err(|_| anyhow::anyhow!("record length {} exceeds u32::MAX", record.len()))?;

        // Extract sort key from raw BAM bytes
        let sort_key = extract_coordinate_key_inline(record, self.nref);

        // Reserve contiguous space for header + record
        let offset = self.data.reserve_contiguous(total_bytes) as u64;

        // Write inline header (16 bytes)
        let header = InlineHeader { sort_key, record_len: len, padding: 0 };
        self.data.extend_in_place(&header.sort_key.to_le_bytes());
        self.data.extend_in_place(&header.record_len.to_le_bytes());
        self.data.extend_in_place(&header.padding.to_le_bytes());

        // Write raw BAM data
        self.data.extend_in_place(record);

        // Add to index
        self.refs.push(RecordRef { sort_key, offset, len, padding: 0 });
        Ok(())
    }

    /// Sort the index by key (records stay in place).
    ///
    /// Uses radix sort for O(n×k) performance instead of O(n log n) comparison sort.
    /// Falls back to insertion sort for small arrays.
    pub fn sort(&mut self) {
        radix_sort_record_refs(&mut self.refs);
    }

    /// Sort using parallel radix sort (for large arrays).
    ///
    /// Divides data into chunks, sorts each with radix sort, then merges.
    /// For very large arrays this is faster than single-threaded radix sort.
    pub fn par_sort(&mut self) {
        // For parallel sort, we use chunked radix sort
        // Radix sort is already O(n×k), so parallelizing chunks provides
        // linear speedup without the merge overhead of comparison-based parallel sorts
        parallel_radix_sort_record_refs(&mut self.refs);
    }

    /// Sort in parallel and return each sub-array as a separate sorted chunk.
    ///
    /// Instead of merging the parallel sort sub-arrays back into one sorted
    /// buffer (as `par_sort` does), this returns each sub-array as its own
    /// `Vec<(RawCoordinateKey, Vec<u8>)>` so they can be passed as separate
    /// merge sources to the k-way merge, avoiding the intermediate merge step.
    ///
    /// When `threads <= 1`, returns a single chunk.
    pub fn par_sort_into_chunks(
        &mut self,
        threads: usize,
    ) -> Vec<Vec<(RawCoordinateKey, Vec<u8>)>> {
        use rayon::prelude::*;

        let n = self.refs.len();

        if threads <= 1 || n < RADIX_THRESHOLD * 2 || n <= 10_000 {
            // Single-threaded path: sort everything, return one chunk
            radix_sort_record_refs(&mut self.refs);
            let chunk = self
                .refs
                .iter()
                .map(|r| (RawCoordinateKey { sort_key: r.sort_key }, self.get_record(r).to_vec()))
                .collect();
            return vec![chunk];
        }

        let chunk_size = n.div_ceil(threads);

        // Sort each chunk in parallel using radix sort
        self.refs.par_chunks_mut(chunk_size).for_each(|chunk| {
            radix_sort_record_refs(chunk);
        });

        // Collect each sorted sub-array into its own Vec<(RawCoordinateKey, Vec<u8>)>
        self.refs
            .chunks(chunk_size)
            .map(|chunk| {
                chunk
                    .iter()
                    .map(|r| {
                        (RawCoordinateKey { sort_key: r.sort_key }, self.get_record(r).to_vec())
                    })
                    .collect()
            })
            .collect()
    }

    /// Get record bytes by reference.
    #[inline]
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub fn get_record(&self, r: &RecordRef) -> &[u8] {
        self.data.slice(r.offset as usize + HEADER_SIZE, r.len as usize)
    }

    /// Iterate over sorted records.
    pub fn iter_sorted(&self) -> impl Iterator<Item = &[u8]> {
        self.refs.iter().map(|r| self.get_record(r))
    }

    /// Get the record references.
    #[must_use]
    pub fn refs(&self) -> &[RecordRef] {
        &self.refs
    }

    /// Memory usage in bytes (actual data stored, not capacity).
    #[must_use]
    pub fn memory_usage(&self) -> usize {
        self.data.len() + self.refs.len() * std::mem::size_of::<RecordRef>()
    }

    /// Number of records.
    #[must_use]
    pub fn len(&self) -> usize {
        self.refs.len()
    }

    /// Check if buffer is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.refs.is_empty()
    }

    /// Clear the buffer for reuse.
    pub fn clear(&mut self) {
        self.data.clear();
        self.refs.clear();
    }

    /// Get the number of reference sequences.
    #[must_use]
    pub fn nref(&self) -> u32 {
        self.nref
    }
}

// ============================================================================
// Key Extraction Functions
// ============================================================================

/// Extract coordinate sort key directly from raw BAM bytes.
///
/// This is a zero-copy operation that reads fields at fixed offsets.
///
/// For coordinate sorting (following samtools behavior):
/// - Reads with valid tid (>= 0) are sorted by (tid, pos, reverse)
/// - Reads with tid = -1 (no reference) sort at the end
/// - Unmapped reads with a valid tid use that tid for sorting (typically mate's position)
///
/// Note: No read name tie-breaking is used, matching samtools behavior.
/// Equal records maintain their original input order (stable sort).
#[inline]
#[must_use]
pub fn extract_coordinate_key_inline(bam: &[u8], nref: u32) -> u64 {
    let tid = bam_fields::ref_id(bam);
    let pos = bam_fields::pos(bam);
    let reverse = (bam_fields::flags(bam) & bam_fields::flags::REVERSE) != 0;

    // Pack key based on tid (samtools behavior):
    // - tid >= 0: sort by (tid, pos, reverse) even if unmapped flag is set
    // - tid < 0: unmapped with no reference, sort at end
    if tid < 0 {
        PackedCoordinateKey::unmapped().0
    } else {
        PackedCoordinateKey::new(tid, pos, reverse, nref).0
    }
}

// ============================================================================
// Template-Coordinate Support
// ============================================================================

/// Extended key for template-coordinate sorting.
///
/// Template-coordinate requires comparing multiple fields:
/// tid1, tid2, pos1, pos2, neg1, neg2, CB, library, MI, name, `is_upper`
///
/// We pack these into 5 u64 values for efficient comparison.
/// The `name_hash_upper` field packs both `name_hash` and `is_upper`:
/// - Upper 63 bits: name hash (groups same names together)
/// - Lowest bit: `is_upper` (false=0, true=1)
///
/// This ensures reads from the same template stay together (same hash),
/// with `is_upper=false` sorting before `is_upper=true`.
#[repr(C)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, bytemuck::Pod, bytemuck::Zeroable)]
pub struct TemplateKey {
    /// Packed: (tid1 << 48) | (tid2 << 32) | pos1
    /// Comparison order matches samtools: tid1, tid2, pos1
    pub primary: u64,
    /// Packed: (pos2 << 32) | (!neg1 << 1) | !neg2
    /// neg flags inverted so reverse (neg=true) sorts before forward (neg=false)
    pub secondary: u64,
    /// Hash of the CB (cellular barcode) tag value.
    /// Inserted between neg2 and MI to match fgbio's sort order.
    /// 0 when no cell tag is present or `cell_tag` is disabled.
    pub cb_hash: u64,
    /// Packed: (library << 48) | (`mi_value` << 1) | `mi_suffix`
    pub tertiary: u64,
    /// Packed: (`name_hash` << 1) | `is_upper`
    /// This ensures same-name records group together, with `is_upper` as tie-breaker.
    pub name_hash_upper: u64,
}

impl TemplateKey {
    /// Create a new template key from extracted fields.
    #[allow(clippy::too_many_arguments, clippy::cast_sign_loss)]
    #[must_use]
    pub fn new(
        tid1: i32,
        pos1: i32,
        neg1: bool,
        tid2: i32,
        pos2: i32,
        neg2: bool,
        cb_hash: u64,
        library: u32,
        mi: (u64, bool),
        name_hash: u64,
        is_upper: bool,
    ) -> Self {
        // Handle i32::MAX specially (indicates unmapped mate)
        let tid1_packed = if tid1 == i32::MAX { 0xFFFF_u64 } else { (tid1.max(0) as u64) & 0xFFFF };
        let tid2_packed = if tid2 == i32::MAX { 0xFFFF_u64 } else { (tid2.max(0) as u64) & 0xFFFF };
        // Convert signed positions to unsigned preserving sort order (XOR with sign bit).
        // This ensures negative positions sort correctly before positive ones.
        // i32::MAX maps to 0xFFFFFFFF which sorts last (used for unmapped mate sentinel).
        let pos1_packed = u64::from((pos1 as u32) ^ 0x8000_0000) & 0xFFFF_FFFF;
        let pos2_packed = u64::from((pos2 as u32) ^ 0x8000_0000) & 0xFFFF_FFFF;

        // Pack primary: tid1 (bits 63-48), tid2 (bits 47-32), pos1 (bits 31-0)
        // This ensures comparison order matches samtools: tid1, tid2, pos1
        let p1 = (tid1_packed << 48) | (tid2_packed << 32) | pos1_packed;

        // Pack secondary: pos2 (bits 63-32), neg1 (bit 1), neg2 (bit 0)
        // Invert neg flags: samtools sorts reverse (neg=true) BEFORE forward (neg=false)
        // By storing !neg, we get: !true=0 < !false=1, so reverse sorts first
        let p2 = (pos2_packed << 32) | (u64::from(!neg1) << 1) | u64::from(!neg2);

        // Pack tertiary: library (high 16), mi_value (middle), mi_suffix (bit 0)
        // Note: /B suffix should sort after /A, so we use !is_a as the bit
        let p3 = ((u64::from(library) & 0xFFFF) << 48)
            | ((mi.0 & 0xFFFF_FFFF_FFFF) << 1)
            | u64::from(!mi.1);

        // Pack name_hash and is_upper: hash in upper 63 bits, is_upper in bit 0
        // This ensures same-name records group together, with is_upper=false before is_upper=true
        let p4 = (name_hash << 1) | u64::from(is_upper);

        Self { primary: p1, secondary: p2, cb_hash, tertiary: p3, name_hash_upper: p4 }
    }

    /// Create a key for completely unmapped records.
    #[must_use]
    pub fn unmapped(name_hash: u64, cb_hash: u64, is_read2: bool) -> Self {
        Self {
            primary: u64::MAX,
            secondary: u64::MAX,
            cb_hash,
            tertiary: 0,
            name_hash_upper: (name_hash << 1) | u64::from(is_read2),
        }
    }

    /// Create a zeroed key (used as dummy for memory operations).
    #[inline]
    #[must_use]
    pub fn zeroed() -> Self {
        Self { primary: 0, secondary: 0, cb_hash: 0, tertiary: 0, name_hash_upper: 0 }
    }
}

impl Default for TemplateKey {
    fn default() -> Self {
        Self::zeroed()
    }
}

impl TemplateKey {
    /// Serialize to bytes for storage in keyed temp files.
    #[inline]
    #[must_use]
    pub fn to_bytes(&self) -> [u8; 40] {
        let mut buf = [0u8; 40];
        buf[0..8].copy_from_slice(&self.primary.to_le_bytes());
        buf[8..16].copy_from_slice(&self.secondary.to_le_bytes());
        buf[16..24].copy_from_slice(&self.cb_hash.to_le_bytes());
        buf[24..32].copy_from_slice(&self.tertiary.to_le_bytes());
        buf[32..40].copy_from_slice(&self.name_hash_upper.to_le_bytes());
        buf
    }

    /// Deserialize from bytes read from keyed temp files.
    #[inline]
    #[must_use]
    pub fn from_bytes(buf: &[u8; 40]) -> Self {
        Self {
            primary: u64::from_le_bytes([
                buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7],
            ]),
            secondary: u64::from_le_bytes([
                buf[8], buf[9], buf[10], buf[11], buf[12], buf[13], buf[14], buf[15],
            ]),
            cb_hash: u64::from_le_bytes([
                buf[16], buf[17], buf[18], buf[19], buf[20], buf[21], buf[22], buf[23],
            ]),
            tertiary: u64::from_le_bytes([
                buf[24], buf[25], buf[26], buf[27], buf[28], buf[29], buf[30], buf[31],
            ]),
            name_hash_upper: u64::from_le_bytes([
                buf[32], buf[33], buf[34], buf[35], buf[36], buf[37], buf[38], buf[39],
            ]),
        }
    }
}

impl PartialOrd for TemplateKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TemplateKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.primary
            .cmp(&other.primary)
            .then_with(|| self.secondary.cmp(&other.secondary))
            .then_with(|| self.cb_hash.cmp(&other.cb_hash))
            .then_with(|| self.tertiary.cmp(&other.tertiary))
            // name_hash_upper comparison handles both name grouping AND is_upper ordering
            .then_with(|| self.name_hash_upper.cmp(&other.name_hash_upper))
    }
}

impl TemplateKey {
    /// Compare only core fields (tid1, tid2, pos1, pos2, neg1, neg2, CB, library, MI).
    ///
    /// This ignores the `name_hash` tie-breaker, allowing verification to accept
    /// both fgumi and samtools sorted files (which differ only in tie-breaking).
    #[inline]
    #[must_use]
    pub fn core_cmp(&self, other: &Self) -> Ordering {
        self.primary
            .cmp(&other.primary)
            .then_with(|| self.secondary.cmp(&other.secondary))
            .then_with(|| self.cb_hash.cmp(&other.cb_hash))
            .then_with(|| self.tertiary.cmp(&other.tertiary))
    }
}

impl RawSortKey for TemplateKey {
    const SERIALIZED_SIZE: Option<usize> = Some(40);

    /// # Panics
    ///
    /// Always panics. `TemplateKey` extraction requires a [`LibraryLookup`]
    /// context not available through the `RawSortKey` trait interface. All
    /// callers must use `extract_template_key_inline()` instead.
    fn extract(_bam: &[u8], _ctx: &SortContext) -> Self {
        unreachable!(
            "TemplateKey::extract() should not be called directly. \
             Use extract_template_key_inline() with LibraryLookup instead."
        )
    }

    #[inline]
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writer.write_all(&self.to_bytes())
    }

    #[inline]
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let mut buf = [0u8; 40];
        reader.read_exact(&mut buf)?;
        Ok(Self::from_bytes(&buf))
    }
}

/// Inline header stored before each record in the data buffer.
/// This allows us to use a minimal ref structure while still having
/// fast access to the full sort key.
///
/// Layout (48 bytes total):
/// - primary: 8 bytes (u64)
/// - secondary: 8 bytes (u64)
/// - `cb_hash`: 8 bytes (u64)
/// - tertiary: 8 bytes (u64)
/// - `name_hash_upper`: 8 bytes (u64)
/// - `record_len`: 4 bytes (u32)
/// - padding: 4 bytes (u32)
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct TemplateInlineHeader {
    /// Full sort key for comparison.
    pub key: TemplateKey,
    /// Length of following raw BAM data.
    pub record_len: u32,
    /// Padding for 8-byte alignment.
    pub padding: u32,
}

/// Size of `TemplateInlineHeader` in bytes.
pub const TEMPLATE_HEADER_SIZE: usize = 48; // 5 * 8 (key) + 4 + 4
const _: () = assert!(
    std::mem::size_of::<TemplateInlineHeader>() == TEMPLATE_HEADER_SIZE,
    "TEMPLATE_HEADER_SIZE must match size_of::<TemplateInlineHeader>()"
);

impl TemplateInlineHeader {
    /// Serialize header to a byte array.
    #[inline]
    #[must_use]
    pub fn to_bytes(&self) -> [u8; TEMPLATE_HEADER_SIZE] {
        let mut buf = [0u8; TEMPLATE_HEADER_SIZE];
        buf[0..8].copy_from_slice(&self.key.primary.to_le_bytes());
        buf[8..16].copy_from_slice(&self.key.secondary.to_le_bytes());
        buf[16..24].copy_from_slice(&self.key.cb_hash.to_le_bytes());
        buf[24..32].copy_from_slice(&self.key.tertiary.to_le_bytes());
        buf[32..40].copy_from_slice(&self.key.name_hash_upper.to_le_bytes());
        buf[40..44].copy_from_slice(&self.record_len.to_le_bytes());
        buf[44..48].copy_from_slice(&self.padding.to_le_bytes());
        buf
    }

    /// Read header from a byte slice.
    #[inline]
    #[must_use]
    pub fn read_from(data: &[u8]) -> Self {
        let primary = u64::from_le_bytes([
            data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7],
        ]);
        let secondary = u64::from_le_bytes([
            data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15],
        ]);
        let cb_hash = u64::from_le_bytes([
            data[16], data[17], data[18], data[19], data[20], data[21], data[22], data[23],
        ]);
        let tertiary = u64::from_le_bytes([
            data[24], data[25], data[26], data[27], data[28], data[29], data[30], data[31],
        ]);
        let name_hash_upper = u64::from_le_bytes([
            data[32], data[33], data[34], data[35], data[36], data[37], data[38], data[39],
        ]);
        let record_len = u32::from_le_bytes([data[40], data[41], data[42], data[43]]);

        Self {
            key: TemplateKey { primary, secondary, cb_hash, tertiary, name_hash_upper },
            record_len,
            padding: 0,
        }
    }
}

/// Record reference for template-coordinate sorting with cached key.
///
/// Caches the full `TemplateKey` inline for O(1) comparisons during sort.
/// This trades memory (48 bytes vs 16 bytes per ref) for cache locality -
/// all comparison data is in the ref itself, avoiding random access to
/// the multi-GB data buffer during sorting.
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct TemplateRecordRef {
    /// Cached sort key for O(1) comparisons without data buffer access.
    pub key: TemplateKey,
    /// Offset to inline header in data buffer.
    pub offset: u64,
    /// Length of raw BAM data (excluding inline header).
    pub len: u32,
    /// Padding for alignment.
    pub padding: u32,
}

impl PartialEq for TemplateRecordRef {
    fn eq(&self, other: &Self) -> bool {
        self.offset == other.offset
    }
}

impl Eq for TemplateRecordRef {}

impl PartialOrd for TemplateRecordRef {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TemplateRecordRef {
    fn cmp(&self, other: &Self) -> Ordering {
        // Note: This Ord impl is NOT used for sorting - we use sort_by with data access
        self.offset.cmp(&other.offset)
    }
}

/// Template-coordinate record buffer with inline headers.
///
/// Uses inline headers to store full sort keys in the data buffer,
/// allowing minimal refs (16 bytes) while maintaining fast comparison.
///
/// Memory layout:
/// ```text
/// data: [Header0][Record0][Header1][Record1]...
/// refs: [Ref0][Ref1]...
/// ```
///
/// During sorting, we use a custom comparator that:
/// 1. Compares primary keys from refs (fast, O(1))
/// 2. On ties, fetches full keys from inline headers
pub struct TemplateRecordBuffer {
    /// Segmented byte storage: inline headers + record data.
    data: SegmentedBuf,
    /// Minimal index for sorting.
    refs: Vec<TemplateRecordRef>,
}

/// Segment size for `TemplateRecordBuffer`'s `SegmentedBuf`: 256 MiB.
const TEMPLATE_SEGMENT_SIZE: usize = 256 * 1024 * 1024;

impl TemplateRecordBuffer {
    /// Create a new buffer with estimated capacity.
    #[must_use]
    pub fn with_capacity(estimated_records: usize, estimated_bytes: usize) -> Self {
        let header_bytes = estimated_records * TEMPLATE_HEADER_SIZE;
        Self {
            data: SegmentedBuf::with_capacity(
                estimated_bytes + header_bytes,
                TEMPLATE_SEGMENT_SIZE,
            ),
            refs: Vec::with_capacity(estimated_records),
        }
    }

    /// Push a record with a pre-computed template key.
    ///
    /// # Errors
    ///
    /// Returns an error if the record (plus header) exceeds the segment size
    /// (256 MiB) or if the record length exceeds `u32::MAX`.
    #[inline]
    pub fn push(&mut self, record: &[u8], key: TemplateKey) -> anyhow::Result<()> {
        let total_bytes = TEMPLATE_HEADER_SIZE + record.len();
        anyhow::ensure!(
            total_bytes <= TEMPLATE_SEGMENT_SIZE,
            "BAM record of {} bytes (+ {} byte header) exceeds segment size of {} bytes; \
             this is likely a malformed BAM file",
            record.len(),
            TEMPLATE_HEADER_SIZE,
            TEMPLATE_SEGMENT_SIZE,
        );
        let record_len = u32::try_from(record.len())
            .map_err(|_| anyhow::anyhow!("record length {} exceeds u32::MAX", record.len()))?;

        // Reserve contiguous space for header + record
        let offset = self.data.reserve_contiguous(total_bytes) as u64;

        // Write inline header
        let header = TemplateInlineHeader { key, record_len, padding: 0 };
        self.data.extend_in_place(&header.to_bytes());

        // Write raw BAM data
        self.data.extend_in_place(record);

        // Add ref with cached key for O(1) sort comparisons
        self.refs.push(TemplateRecordRef { key, offset, len: record_len, padding: 0 });
        Ok(())
    }

    /// Sort the index by cached key using stable LSD radix sort.
    ///
    /// Uses multi-field radix sort which is stable (preserves relative order
    /// of records with equal keys). This ensures deterministic output that
    /// matches samtools when records have identical sort keys.
    pub fn sort(&mut self) {
        radix_sort_template_refs(&mut self.refs);
    }

    /// Sort using parallel radix sort with stable k-way merge.
    ///
    /// Each chunk is sorted with stable radix sort, then merged with a
    /// heap that uses `chunk_idx` as tie-breaker to preserve input order.
    pub fn par_sort(&mut self) {
        parallel_radix_sort_template_refs(&mut self.refs);
    }

    /// Get record bytes by reference.
    #[inline]
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub fn get_record(&self, r: &TemplateRecordRef) -> &[u8] {
        self.data.slice(r.offset as usize + TEMPLATE_HEADER_SIZE, r.len as usize)
    }

    /// Iterate over sorted records.
    pub fn iter_sorted(&self) -> impl Iterator<Item = &[u8]> {
        self.refs.iter().map(|r| self.get_record(r))
    }

    /// Get the record references.
    #[must_use]
    pub fn refs(&self) -> &[TemplateRecordRef] {
        &self.refs
    }

    /// Memory usage in bytes (actual data stored, not capacity).
    #[must_use]
    pub fn memory_usage(&self) -> usize {
        self.data.len() + self.refs.len() * std::mem::size_of::<TemplateRecordRef>()
    }

    /// Number of records.
    #[must_use]
    pub fn len(&self) -> usize {
        self.refs.len()
    }

    /// Check if buffer is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.refs.is_empty()
    }

    /// Get the key for a record reference (returns cached key from ref).
    #[inline]
    #[must_use]
    pub fn get_key(&self, r: &TemplateRecordRef) -> TemplateKey {
        r.key
    }

    /// Iterate over sorted (key, record) pairs.
    /// Used for writing keyed temp chunks that preserve sort keys.
    pub fn iter_sorted_keyed(&self) -> impl Iterator<Item = (TemplateKey, &[u8])> {
        self.refs.iter().map(|r| (self.get_key(r), self.get_record(r)))
    }

    /// Clear the buffer for reuse.
    pub fn clear(&mut self) {
        self.data.clear();
        self.refs.clear();
    }

    /// Sort in parallel and return each sub-array as a separate sorted chunk.
    ///
    /// Instead of merging the parallel sort sub-arrays back into one sorted
    /// buffer (as `par_sort` does), this returns each sub-array as its own
    /// `Vec<(TemplateKey, Vec<u8>)>` so they can be passed as separate merge
    /// sources to the k-way merge, avoiding the intermediate merge step.
    ///
    /// When `threads <= 1`, returns a single chunk.
    pub fn par_sort_into_chunks(&mut self, threads: usize) -> Vec<Vec<(TemplateKey, Vec<u8>)>> {
        use rayon::prelude::*;

        let n = self.refs.len();

        if threads <= 1 || n < RADIX_THRESHOLD * 2 || n <= 10_000 {
            // Single-threaded path: sort everything, return one chunk
            radix_sort_template_refs(&mut self.refs);
            let chunk = self.refs.iter().map(|r| (r.key, self.get_record(r).to_vec())).collect();
            return vec![chunk];
        }

        let chunk_size = n.div_ceil(threads);

        // Sort each chunk in parallel using radix sort
        self.refs.par_chunks_mut(chunk_size).for_each(|chunk| {
            radix_sort_template_refs(chunk);
        });

        // Collect each sorted sub-array into its own Vec<(TemplateKey, Vec<u8>)>
        self.refs
            .chunks(chunk_size)
            .map(|chunk| chunk.iter().map(|r| (r.key, self.get_record(r).to_vec())).collect())
            .collect()
    }
}

// ============================================================================
// Radix Sort for RecordRef
// ============================================================================

// NOTE: QuerynameRecordBuffer was removed — it had zero callers and proved to
// regress when wired in (reconstruct_record allocates per record on write-out).

/// Threshold below which we use insertion sort instead of radix sort.
const RADIX_THRESHOLD: usize = 256;

/// Radix sort for `RecordRef` arrays using LSD (Least Significant Digit) approach.
///
/// Sorts by the `sort_key` field using 8-bit radix (256 buckets).
/// This is O(n×k) where k is the number of bytes to sort (typically 5-8).
///
/// # Stability
/// Radix sort is inherently stable - records with equal keys maintain their
/// relative input order, matching samtools behavior.
#[allow(clippy::uninit_vec, unsafe_code)]
pub fn radix_sort_record_refs(refs: &mut [RecordRef]) {
    let n = refs.len();
    if n < RADIX_THRESHOLD {
        // Use insertion sort for small arrays
        insertion_sort_refs(refs);
        return;
    }

    // Find max key to determine how many bytes we need to sort
    let max_key = refs.iter().map(|r| r.sort_key).max().unwrap_or(0);
    let bytes_needed =
        if max_key == 0 { 0 } else { ((64 - max_key.leading_zeros()) as usize).div_ceil(8) };

    if bytes_needed == 0 {
        return; // All keys are 0, already sorted
    }

    // Allocate auxiliary buffer
    let mut aux: Vec<RecordRef> = Vec::with_capacity(n);
    unsafe {
        aux.set_len(n);
    }

    let mut src = refs as *mut [RecordRef];
    let mut dst = aux.as_mut_slice() as *mut [RecordRef];

    // LSD radix sort - byte by byte from least significant
    for byte_idx in 0..bytes_needed {
        let src_slice = unsafe { &*src };
        let dst_slice = unsafe { &mut *dst };

        // Count occurrences of each byte value
        let mut counts = [0usize; 256];
        for r in src_slice {
            let byte = ((r.sort_key >> (byte_idx * 8)) & 0xFF) as usize;
            counts[byte] += 1;
        }

        // Convert to cumulative offsets
        let mut total = 0;
        for count in &mut counts {
            let c = *count;
            *count = total;
            total += c;
        }

        // Scatter elements to destination (stable - preserves order within buckets)
        for r in src_slice {
            let byte = ((r.sort_key >> (byte_idx * 8)) & 0xFF) as usize;
            let dest_idx = counts[byte];
            counts[byte] += 1;
            dst_slice[dest_idx] = *r;
        }

        // Swap src and dst
        std::mem::swap(&mut src, &mut dst);
    }

    // If odd number of passes, copy back to original buffer
    if bytes_needed % 2 == 1 {
        let src_slice = unsafe { &*src };
        refs.copy_from_slice(src_slice);
    }
}

/// Parallel radix sort for `RecordRef` arrays.
///
/// Divides the array into chunks, sorts each chunk with radix sort,
/// then performs k-way merge. This provides near-linear speedup.
pub fn parallel_radix_sort_record_refs(refs: &mut [RecordRef]) {
    use rayon::prelude::*;

    let n = refs.len();
    if n < RADIX_THRESHOLD * 2 {
        // Small array - just use single-threaded radix sort
        radix_sort_record_refs(refs);
        return;
    }

    // Get number of threads from rayon
    let n_threads = rayon::current_num_threads();

    // For very large arrays, parallel chunked sort + merge is faster
    if n_threads > 1 && n > 10_000 {
        let chunk_size = n.div_ceil(n_threads);

        // Sort each chunk in parallel using radix sort
        refs.par_chunks_mut(chunk_size).for_each(|chunk| {
            radix_sort_record_refs(chunk);
        });

        // K-way merge the sorted chunks
        // For simplicity, we use a heap-based merge into auxiliary storage
        let chunk_boundaries: Vec<_> = refs
            .chunks(chunk_size)
            .scan(0, |pos, chunk| {
                let start = *pos;
                *pos += chunk.len();
                Some(start..*pos)
            })
            .collect();

        merge_sorted_chunks(refs, &chunk_boundaries);
    } else {
        // Single-threaded radix sort
        radix_sort_record_refs(refs);
    }
}

/// Merge k sorted chunks in place using auxiliary storage.
fn merge_sorted_chunks(refs: &mut [RecordRef], chunk_ranges: &[std::ops::Range<usize>]) {
    use crate::sort::radix::{heap_make, heap_sift_down};

    // Heap entry: (sort_key, chunk_idx, position_in_chunk)
    struct HeapEntry {
        key: u64,
        chunk_idx: usize,
        pos: usize,
    }

    if chunk_ranges.len() <= 1 {
        return;
    }

    let n = refs.len();
    let mut result: Vec<RecordRef> = Vec::with_capacity(n);

    // Initialize heap with first element from each chunk
    let mut heap: Vec<HeapEntry> = Vec::with_capacity(chunk_ranges.len());
    for (chunk_idx, range) in chunk_ranges.iter().enumerate() {
        if !range.is_empty() {
            heap.push(HeapEntry { key: refs[range.start].sort_key, chunk_idx, pos: range.start });
        }
    }

    if heap.is_empty() {
        return;
    }

    // Min-heap: smaller keys should be at the top
    // Use `chunk_idx` as tie-breaker for stability
    let lt = |a: &HeapEntry, b: &HeapEntry| -> bool { (a.key, a.chunk_idx) > (b.key, b.chunk_idx) };

    heap_make(&mut heap, &lt);
    let mut heap_size = heap.len();

    // Merge
    while heap_size > 0 {
        let entry = &heap[0];
        result.push(refs[entry.pos]);

        let chunk_idx = entry.chunk_idx;
        let next_pos = entry.pos + 1;
        let range = &chunk_ranges[chunk_idx];

        if next_pos < range.end {
            // More elements in this chunk
            heap[0] = HeapEntry { key: refs[next_pos].sort_key, chunk_idx, pos: next_pos };
            heap_sift_down(&mut heap, 0, heap_size, &lt);
        } else {
            // Chunk exhausted
            heap_size -= 1;
            if heap_size > 0 {
                heap.swap(0, heap_size);
                heap_sift_down(&mut heap, 0, heap_size, &lt);
            }
        }
    }

    // Copy result back
    refs.copy_from_slice(&result);
}

/// Binary insertion sort for small arrays of `RecordRef`.
fn insertion_sort_refs(refs: &mut [RecordRef]) {
    for i in 1..refs.len() {
        let key = refs[i].sort_key;
        let insert_pos = refs[..i].partition_point(|r| r.sort_key <= key);
        if insert_pos < i {
            refs[insert_pos..=i].rotate_right(1);
        }
    }
}

// ============================================================================
// Radix Sort for TemplateRecordRef (multi-field 4×u64 key)
// ============================================================================

/// MSD (Most Significant Digit) hybrid radix sort for `TemplateRecordRef` arrays.
///
/// Exploits the fact that the `primary` field (tid1/tid2/pos1) is highly discriminating:
/// most records have unique primary values in typical sequencing data.
///
/// Strategy:
/// 1. Radix sort by `primary` field — O(n × k₁) where k₁ = adaptive bytes needed
/// 2. Find runs of equal `primary` values (single O(n) scan)
/// 3. Sub-sort only the runs that need it:
///    - Size 1: already sorted (majority of records)
///    - Size 2–64: insertion sort by remaining fields
///    - Size >64: LSD radix sort remaining fields
///
/// This avoids full O(n) passes for fields 2–5 when most records are already
/// resolved by the primary field alone.
#[allow(clippy::uninit_vec, unsafe_code)]
pub fn radix_sort_template_refs(refs: &mut [TemplateRecordRef]) {
    let n = refs.len();
    if n < RADIX_THRESHOLD {
        insertion_sort_template_refs(refs);
        return;
    }

    // Allocate auxiliary buffer (reused across phases)
    let mut aux: Vec<TemplateRecordRef> = Vec::with_capacity(n);
    unsafe {
        aux.set_len(n);
    }

    // Phase 1: Radix sort by primary field
    let max_primary = refs.iter().map(|r| r.key.primary).max().unwrap_or(0);
    let bytes_needed = bytes_needed_u64(max_primary);
    if bytes_needed > 0 {
        radix_sort_template_field(refs, &mut aux, |r| r.key.primary, bytes_needed);
    }

    // Phase 2: Find equal-primary runs and sub-sort by remaining fields
    sub_sort_runs(refs, &mut aux, |r| r.key.primary, &REMAINING_FIELDS_AFTER_PRIMARY);
}

/// Remaining `TemplateKey` fields after `primary`, in sort precedence order.
const REMAINING_FIELDS_AFTER_PRIMARY: [fn(&TemplateRecordRef) -> u64; 4] =
    [|r| r.key.secondary, |r| r.key.cb_hash, |r| r.key.tertiary, |r| r.key.name_hash_upper];

/// Threshold for sub-sort runs: below this, use insertion sort.
const SUB_SORT_INSERTION_THRESHOLD: usize = 64;

/// Find runs of equal values for `run_field` and sub-sort each run by `remaining_fields`.
fn sub_sort_runs<F>(
    refs: &mut [TemplateRecordRef],
    aux: &mut [TemplateRecordRef],
    run_field: F,
    remaining_fields: &[fn(&TemplateRecordRef) -> u64],
) where
    F: Fn(&TemplateRecordRef) -> u64,
{
    if remaining_fields.is_empty() {
        return;
    }

    let n = refs.len();
    let mut start = 0;
    while start < n {
        let val = run_field(&refs[start]);
        let mut end = start + 1;
        while end < n && run_field(&refs[end]) == val {
            end += 1;
        }

        let run = &mut refs[start..end];
        let run_len = run.len();
        if run_len > 1 {
            if run_len <= SUB_SORT_INSERTION_THRESHOLD {
                insertion_sort_template_refs(run);
            } else {
                let next_field = remaining_fields[0];
                let max_val = run.iter().map(next_field).max().unwrap_or(0);
                let bytes_needed = bytes_needed_u64(max_val);
                if bytes_needed > 0 {
                    let run_aux = &mut aux[start..end];
                    radix_sort_template_field(run, run_aux, next_field, bytes_needed);
                }

                if remaining_fields.len() > 1 {
                    sub_sort_runs(run, &mut aux[start..end], next_field, &remaining_fields[1..]);
                }
            }
        }
        start = end;
    }
}

/// Radix sort a single u64 field of `TemplateRecordRef` using raw pointers.
#[allow(clippy::uninit_vec, unsafe_code)]
fn radix_sort_template_field<F>(
    refs: &mut [TemplateRecordRef],
    aux: &mut [TemplateRecordRef],
    get_field: F,
    bytes_needed: usize,
) where
    F: Fn(&TemplateRecordRef) -> u64,
{
    let n = refs.len();

    // Use raw pointers to avoid borrow checker issues with swapping
    let mut src = refs as *mut [TemplateRecordRef];
    let mut dst = aux as *mut [TemplateRecordRef];

    for byte_idx in 0..bytes_needed {
        let src_slice = unsafe { &*src };
        let dst_slice = unsafe { &mut *dst };

        // Count occurrences of each byte value
        let mut counts = [0usize; 256];
        for r in src_slice {
            let byte = ((get_field(r) >> (byte_idx * 8)) & 0xFF) as usize;
            counts[byte] += 1;
        }

        // Convert to cumulative offsets
        let mut total = 0;
        for count in &mut counts {
            let c = *count;
            *count = total;
            total += c;
        }

        // Scatter elements to destination
        for item in src_slice.iter().take(n) {
            let byte = ((get_field(item) >> (byte_idx * 8)) & 0xFF) as usize;
            let dest_idx = counts[byte];
            counts[byte] += 1;
            dst_slice[dest_idx] = *item;
        }

        // Swap src and dst
        std::mem::swap(&mut src, &mut dst);
    }

    // If odd number of passes, copy back to original buffer
    if bytes_needed % 2 == 1 {
        let src_slice = unsafe { &*src };
        refs.copy_from_slice(src_slice);
    }
}

/// Parallel radix sort for `TemplateRecordRef` arrays.
///
/// Uses stable radix sort per chunk with stable k-way merge (`chunk_idx` tie-breaker).
pub fn parallel_radix_sort_template_refs(refs: &mut [TemplateRecordRef]) {
    use rayon::prelude::*;

    let n = refs.len();
    if n < RADIX_THRESHOLD * 2 {
        radix_sort_template_refs(refs);
        return;
    }

    let n_threads = rayon::current_num_threads();

    if n_threads > 1 && n > 10_000 {
        let chunk_size = n.div_ceil(n_threads);

        // Sort each chunk in parallel
        refs.par_chunks_mut(chunk_size).for_each(|chunk| {
            radix_sort_template_refs(chunk);
        });

        // K-way merge the sorted chunks
        let chunk_boundaries: Vec<_> = refs
            .chunks(chunk_size)
            .scan(0, |pos, chunk| {
                let start = *pos;
                *pos += chunk.len();
                Some(start..*pos)
            })
            .collect();

        merge_sorted_template_chunks(refs, &chunk_boundaries);
    } else {
        radix_sort_template_refs(refs);
    }
}

/// Merge k sorted chunks of `TemplateRecordRef` in place.
///
/// Uses `chunk_idx` as tie-breaker for stable merge (preserves input order for equal keys).
fn merge_sorted_template_chunks(
    refs: &mut [TemplateRecordRef],
    chunk_ranges: &[std::ops::Range<usize>],
) {
    use crate::sort::radix::{heap_make, heap_sift_down};

    struct HeapEntry {
        key: TemplateKey,
        chunk_idx: usize,
        pos: usize,
    }

    if chunk_ranges.len() <= 1 {
        return;
    }

    let n = refs.len();
    let mut result: Vec<TemplateRecordRef> = Vec::with_capacity(n);

    let mut heap: Vec<HeapEntry> = Vec::with_capacity(chunk_ranges.len());
    for (chunk_idx, range) in chunk_ranges.iter().enumerate() {
        if !range.is_empty() {
            heap.push(HeapEntry { key: refs[range.start].key, chunk_idx, pos: range.start });
        }
    }

    if heap.is_empty() {
        return;
    }

    // Min-heap with chunk_idx tie-breaker for stability
    let lt = |a: &HeapEntry, b: &HeapEntry| -> bool {
        match a.key.cmp(&b.key) {
            std::cmp::Ordering::Greater => true,
            std::cmp::Ordering::Less => false,
            std::cmp::Ordering::Equal => a.chunk_idx > b.chunk_idx,
        }
    };

    heap_make(&mut heap, &lt);
    let mut heap_size = heap.len();

    while heap_size > 0 {
        let entry = &heap[0];
        result.push(refs[entry.pos]);

        let chunk_idx = entry.chunk_idx;
        let next_pos = entry.pos + 1;
        let range = &chunk_ranges[chunk_idx];

        if next_pos < range.end {
            heap[0] = HeapEntry { key: refs[next_pos].key, chunk_idx, pos: next_pos };
            heap_sift_down(&mut heap, 0, heap_size, &lt);
        } else {
            heap_size -= 1;
            if heap_size > 0 {
                heap.swap(0, heap_size);
                heap_sift_down(&mut heap, 0, heap_size, &lt);
            }
        }
    }

    refs.copy_from_slice(&result);
}

/// Binary insertion sort for small arrays of `TemplateRecordRef`.
fn insertion_sort_template_refs(refs: &mut [TemplateRecordRef]) {
    for i in 1..refs.len() {
        let key = &refs[i].key;
        let insert_pos = refs[..i].partition_point(|r| r.key <= *key);
        if insert_pos < i {
            refs[insert_pos..=i].rotate_right(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_packed_coord_key_ordering() {
        // Lower tid should come first
        assert!(
            PackedCoordinateKey::new(0, 100, false, 10)
                < PackedCoordinateKey::new(1, 100, false, 10)
        );

        // Lower pos should come first
        assert!(
            PackedCoordinateKey::new(0, 100, false, 10)
                < PackedCoordinateKey::new(0, 200, false, 10)
        );

        // Forward should come before reverse (false < true)
        assert!(
            PackedCoordinateKey::new(0, 100, false, 10)
                < PackedCoordinateKey::new(0, 100, true, 10)
        );

        // Unmapped should come last
        assert!(PackedCoordinateKey::new(9, 1_000_000, true, 10) < PackedCoordinateKey::unmapped());
    }

    #[test]
    fn test_template_key_ordering() {
        let k1 = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 0, false);
        let k2 = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (2, true), 0, false);
        assert!(k1 < k2);

        // /A suffix should come before /B
        let ka = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 0, false);
        let kb = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, false), 0, false);
        assert!(ka < kb);

        // Same name hash: is_upper=false should come before is_upper=true
        let lower = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 12345, false);
        let upper = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 12345, true);
        assert!(lower < upper, "is_upper=false should sort before is_upper=true");

        // Different name hashes should group separately
        let first_hash_lo =
            TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 100, false);
        let first_hash_hi =
            TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 100, true);
        let second_hash =
            TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 200, false);
        // first_hash records should come before second_hash records
        assert!(first_hash_lo < second_hash);
        assert!(first_hash_hi < second_hash);
    }

    #[test]
    fn test_radix_sort_record_refs() {
        // Create test refs with various keys
        let mut refs = vec![
            RecordRef { sort_key: 100, offset: 0, len: 10, padding: 0 },
            RecordRef { sort_key: 50, offset: 100, len: 10, padding: 0 },
            RecordRef { sort_key: 200, offset: 200, len: 10, padding: 0 },
            RecordRef { sort_key: 50, offset: 300, len: 10, padding: 0 }, // Duplicate key
            RecordRef { sort_key: 1, offset: 400, len: 10, padding: 0 },
        ];

        radix_sort_record_refs(&mut refs);

        // Verify sorted order
        assert_eq!(refs[0].sort_key, 1);
        assert_eq!(refs[1].sort_key, 50);
        assert_eq!(refs[2].sort_key, 50);
        assert_eq!(refs[3].sort_key, 100);
        assert_eq!(refs[4].sort_key, 200);

        // Verify stability: duplicate keys should maintain original order
        // offset=100 was before offset=300 in input, should remain so
        assert_eq!(refs[1].offset, 100);
        assert_eq!(refs[2].offset, 300);
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_radix_sort_large() {
        // Test with larger array to trigger radix sort (> RADIX_THRESHOLD)
        let mut refs: Vec<RecordRef> = (0..1000)
            .map(|i| RecordRef {
                sort_key: (999 - i) as u64, // Reverse order
                offset: i as u64 * 100,
                len: 10,
                padding: 0,
            })
            .collect();

        radix_sort_record_refs(&mut refs);

        // Verify sorted
        for (i, r) in refs.iter().enumerate() {
            assert_eq!(r.sort_key, i as u64, "Expected sort_key {i} at index {i}");
        }
    }

    #[test]
    fn test_radix_sort_empty() {
        let mut refs: Vec<RecordRef> = Vec::new();
        radix_sort_record_refs(&mut refs);
        assert!(refs.is_empty());
    }

    #[test]
    fn test_radix_sort_single() {
        let mut refs = vec![RecordRef { sort_key: 42, offset: 0, len: 10, padding: 0 }];
        radix_sort_record_refs(&mut refs);
        assert_eq!(refs.len(), 1);
        assert_eq!(refs[0].sort_key, 42);
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_radix_sort_all_same_keys() {
        // All same keys should maintain original order (stability)
        let mut refs: Vec<RecordRef> = (0..100)
            .map(|i| RecordRef { sort_key: 42, offset: i as u64 * 100, len: 10, padding: 0 })
            .collect();

        radix_sort_record_refs(&mut refs);

        // Verify all keys are 42 and order is preserved
        for (i, r) in refs.iter().enumerate() {
            assert_eq!(r.sort_key, 42);
            assert_eq!(r.offset, i as u64 * 100, "Stability violated at index {i}");
        }
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_radix_sort_all_zero_keys() {
        let mut refs: Vec<RecordRef> = (0..50)
            .map(|i| RecordRef { sort_key: 0, offset: i as u64 * 100, len: 10, padding: 0 })
            .collect();

        radix_sort_record_refs(&mut refs);

        // Should be stable for all-zero keys
        for (i, r) in refs.iter().enumerate() {
            assert_eq!(r.sort_key, 0);
            assert_eq!(r.offset, i as u64 * 100);
        }
    }

    #[test]
    fn test_radix_sort_max_key() {
        // Test with maximum u64 values
        let mut refs = vec![
            RecordRef { sort_key: u64::MAX, offset: 0, len: 10, padding: 0 },
            RecordRef { sort_key: 0, offset: 100, len: 10, padding: 0 },
            RecordRef { sort_key: u64::MAX / 2, offset: 200, len: 10, padding: 0 },
        ];

        radix_sort_record_refs(&mut refs);

        assert_eq!(refs[0].sort_key, 0);
        assert_eq!(refs[1].sort_key, u64::MAX / 2);
        assert_eq!(refs[2].sort_key, u64::MAX);
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_parallel_radix_sort() {
        // Test parallel sort with enough elements to trigger parallelism
        let mut refs: Vec<RecordRef> = (0..50_000)
            .map(|i| RecordRef {
                sort_key: (49_999 - i) as u64, // Reverse order
                offset: i as u64 * 100,
                len: 10,
                padding: 0,
            })
            .collect();

        parallel_radix_sort_record_refs(&mut refs);

        // Verify sorted
        for (i, r) in refs.iter().enumerate() {
            assert_eq!(r.sort_key, i as u64, "Expected sort_key {i} at index {i}");
        }
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_parallel_radix_sort_stability() {
        // Test that parallel sort maintains stability for equal keys
        let mut refs: Vec<RecordRef> = (0..20_000)
            .map(|i| RecordRef {
                sort_key: (i / 100) as u64, // Groups of 100 with same key
                offset: i as u64,           // Use offset to track original order
                len: 10,
                padding: 0,
            })
            .collect();

        parallel_radix_sort_record_refs(&mut refs);

        // Verify sorted and stable within groups
        for i in 1..refs.len() {
            assert!(refs[i - 1].sort_key <= refs[i].sort_key, "Not sorted at index {i}");
            // Within same key group, offsets should be in order
            if refs[i - 1].sort_key == refs[i].sort_key {
                assert!(refs[i - 1].offset < refs[i].offset, "Stability violated at index {i}");
            }
        }
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_radix_sort_template_refs_stability() {
        // Test that template radix sort is stable for equal keys
        // Create refs with identical TemplateKey but different offsets to track order
        let key = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 12345, false);

        let mut refs: Vec<TemplateRecordRef> = (0..500)
            .map(|i| TemplateRecordRef {
                key,              // All same key
                offset: i as u64, // Use offset to track original order
                len: 10,
                padding: 0,
            })
            .collect();

        radix_sort_template_refs(&mut refs);

        // All keys equal, so offsets should maintain original order (stability)
        for i in 1..refs.len() {
            assert!(
                refs[i - 1].offset < refs[i].offset,
                "Template radix sort stability violated at index {}: offset {} should be < {}",
                i,
                refs[i - 1].offset,
                refs[i].offset
            );
        }
    }

    #[test]
    #[allow(clippy::cast_sign_loss)]
    fn test_parallel_radix_sort_template_refs_stability() {
        // Test that parallel template sort maintains stability for equal keys
        // Use groups of records with same key to verify stability within groups
        let mut refs: Vec<TemplateRecordRef> = (0..20_000)
            .map(|i| {
                // Create groups of 100 with same key (different lib_name_hash for each group)
                let group = i / 100;
                let key = TemplateKey::new(
                    0,
                    100,
                    false,
                    0,
                    200,
                    false,
                    0,
                    0,
                    (1, true),
                    group as u64, // Different hash per group
                    false,
                );
                TemplateRecordRef {
                    key,
                    offset: i as u64, // Use offset to track original order
                    len: 10,
                    padding: 0,
                }
            })
            .collect();

        parallel_radix_sort_template_refs(&mut refs);

        // Verify sorted by key and stable within groups
        for i in 1..refs.len() {
            let prev_key = &refs[i - 1].key;
            let curr_key = &refs[i].key;
            assert!(prev_key <= curr_key, "Not sorted at index {i}");
            // Within same key group, offsets should be in order (stability)
            if prev_key == curr_key {
                assert!(
                    refs[i - 1].offset < refs[i].offset,
                    "Parallel template sort stability violated at index {}: offset {} should be < {}",
                    i,
                    refs[i - 1].offset,
                    refs[i].offset
                );
            }
        }
    }

    #[test]
    fn test_template_record_buffer_sort_stability() {
        // Test that TemplateRecordBuffer::sort() is stable
        // This is the method actually used during template-coordinate sorting
        let mut buffer = TemplateRecordBuffer::with_capacity(100, 10000);

        // Add records with identical keys - they should maintain insertion order after sort
        let key = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 12345, false);

        // Create distinct records (different sequence bytes) with same sort key
        for i in 0..100u8 {
            // Minimal valid BAM record: 4 bytes block_size prefix not included in our data
            // ref_id (4) + pos (4) + name_len (1) + mapq (1) + bin (2) + n_cigar (2) + flag (2) + seq_len (4)
            // + mate_ref_id (4) + mate_pos (4) + tlen (4) + name + cigar + seq + qual
            let mut record = vec![
                0, 0, 0, 0, // ref_id = 0
                100, 0, 0, 0, // pos = 100
                2, // name_len = 2 (including null)
                0, // mapq = 0
                0, 0, // bin
                0, 0, // n_cigar_op = 0
                99, 0, // flag = 99 (paired, proper, mate reverse, first)
                1, 0, 0, 0, // seq_len = 1
                0, 0, 0, 0, // mate_ref_id = 0
                200, 0, 0, 0, // mate_pos = 200
                0, 0, 0, 0, // tlen = 0
                b'A', 0, // read name "A\0"
                // no cigar
                i,    // seq (1 byte for 2 bases) - use i to make records distinguishable
                0xFF, // qual
            ];
            // Pad to make it valid
            while record.len() < 40 {
                record.push(0);
            }
            buffer.push(&record, key).expect("push should succeed in tests");
        }

        buffer.sort();

        // Verify records maintain insertion order (tracked by sequence byte value)
        let mut prev_seq_byte = None;
        for rec in buffer.iter_sorted() {
            // Extract the distinguishing byte (seq field at offset 34 after header parsing)
            let seq_byte = rec.get(34).copied().unwrap_or(0);
            if let Some(prev) = prev_seq_byte {
                assert!(
                    prev < seq_byte,
                    "TemplateRecordBuffer::sort() stability violated: {prev} should be < {seq_byte}"
                );
            }
            prev_seq_byte = Some(seq_byte);
        }
    }

    // ========================================================================
    // TemplateKey cb_hash tests
    // ========================================================================

    #[test]
    fn test_template_key_cb_hash_ordering() {
        // Same position but different cb_hash — lower hash sorts first
        let k1 = TemplateKey::new(0, 100, false, 0, 200, false, 10, 0, (1, true), 0, false);
        let k2 = TemplateKey::new(0, 100, false, 0, 200, false, 20, 0, (1, true), 0, false);
        assert!(k1 < k2, "lower cb_hash should sort before higher cb_hash");
    }

    #[test]
    fn test_template_key_cb_hash_zero_sorts_first() {
        // cb_hash=0 (no CB) should sort before non-zero cb_hash
        let k_no_cb = TemplateKey::new(0, 100, false, 0, 200, false, 0, 0, (1, true), 0, false);
        let k_cb = TemplateKey::new(0, 100, false, 0, 200, false, 42, 0, (1, true), 0, false);
        assert!(k_no_cb < k_cb, "cb_hash=0 should sort before non-zero cb_hash");
    }

    #[test]
    fn test_template_key_cb_hash_between_secondary_and_tertiary() {
        // Verify cb_hash is compared AFTER secondary but BEFORE tertiary.
        // Same primary+secondary, different cb_hash and tertiary values
        let k1 = TemplateKey::new(0, 100, false, 0, 200, false, 10, 0, (1, true), 0, false);
        let k2 = TemplateKey::new(0, 100, false, 0, 200, false, 20, 0, (0, true), 0, false);
        // k1 has lower cb_hash, so it should sort first regardless of tertiary
        assert!(k1 < k2, "cb_hash should be compared before tertiary (library/MI)");
    }

    #[test]
    fn test_template_key_serialization_with_cb_hash() {
        let key =
            TemplateKey::new(1, 500, true, 2, 600, false, 0xDEAD_BEEF, 3, (7, true), 999, false);
        let bytes = key.to_bytes();
        assert_eq!(bytes.len(), 40);
        let roundtrip = TemplateKey::from_bytes(&bytes);
        assert_eq!(key, roundtrip, "serialization roundtrip should preserve cb_hash");
        assert_eq!(roundtrip.cb_hash, 0xDEAD_BEEF);
    }

    #[test]
    fn test_template_key_unmapped_with_cb_hash() {
        let key = TemplateKey::unmapped(12345, 0xCAFE, false);
        assert_eq!(key.primary, u64::MAX);
        assert_eq!(key.secondary, u64::MAX);
        assert_eq!(key.cb_hash, 0xCAFE, "unmapped should preserve cb_hash");
        assert_eq!(key.tertiary, 0);
    }

    #[test]
    fn test_template_key_core_cmp_includes_cb_hash() {
        let k1 = TemplateKey::new(0, 100, false, 0, 200, false, 10, 0, (1, true), 0, false);
        let k2 = TemplateKey::new(0, 100, false, 0, 200, false, 20, 0, (1, true), 0, false);
        assert_eq!(k1.core_cmp(&k2), std::cmp::Ordering::Less, "core_cmp should include cb_hash");

        // Same cb_hash, different name_hash — core_cmp should be Equal
        let k3 = TemplateKey::new(0, 100, false, 0, 200, false, 10, 0, (1, true), 100, false);
        let k4 = TemplateKey::new(0, 100, false, 0, 200, false, 10, 0, (1, true), 200, false);
        assert_eq!(k3.core_cmp(&k4), std::cmp::Ordering::Equal, "core_cmp should ignore name_hash");
    }

    #[test]
    fn test_template_key_zeroed_has_zero_cb_hash() {
        let key = TemplateKey::zeroed();
        assert_eq!(key.cb_hash, 0);
    }

    #[test]
    fn test_template_key_default_has_zero_cb_hash() {
        let key = TemplateKey::default();
        assert_eq!(key.cb_hash, 0);
    }

    /// Helper: build a minimal raw BAM record byte array with a distinguishing index.
    fn make_bam_record(index: u16) -> Vec<u8> {
        let mut record = vec![
            0, 0, 0, 0, // ref_id = 0
            100, 0, 0, 0, // pos = 100
            2, // name_len = 2 (including null)
            0, // mapq = 0
            0, 0, // bin
            0, 0, // n_cigar_op = 0
            99, 0, // flag = 99 (paired, proper, mate reverse, first)
            1, 0, 0, 0, // seq_len = 1
            0, 0, 0, 0, // mate_ref_id = 0
            200, 0, 0, 0, // mate_pos = 200
            0, 0, 0, 0, // tlen = 0
            b'A', 0, // read name "A\0"
        ];
        // Encode the index into two bytes so each record is distinguishable
        record.push((index & 0xFF) as u8);
        record.push((index >> 8) as u8);
        record.push(0xFF); // qual
        // Pad to at least 40 bytes
        while record.len() < 40 {
            record.push(0);
        }
        record
    }

    #[test]
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn test_par_sort_into_chunks_single_threaded_fallback() {
        // With 1 thread, par_sort_into_chunks should return exactly 1 chunk.
        let n = 100;
        let mut buffer = TemplateRecordBuffer::with_capacity(n, n * 50);

        for i in 0..n {
            let key = TemplateKey::new(
                0,
                (n - i) as i32,
                false,
                0,
                200,
                false,
                0,
                0,
                (1, true),
                0,
                false,
            );
            buffer.push(&make_bam_record(i as u16), key).expect("push should succeed in tests");
        }

        let chunks = buffer.par_sort_into_chunks(1);
        assert_eq!(chunks.len(), 1, "single-threaded should produce exactly 1 chunk");
        assert_eq!(chunks[0].len(), n, "the single chunk should contain all records");

        // Verify the chunk is sorted
        for i in 1..chunks[0].len() {
            assert!(
                chunks[0][i - 1].0 <= chunks[0][i].0,
                "single chunk should be sorted at index {i}"
            );
        }
    }

    #[test]
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn test_par_sort_into_chunks_parallel_path() {
        // Use a dedicated rayon thread pool with multiple threads so that
        // rayon::current_num_threads() > 1 inside par_sort_into_chunks.
        let n: usize = 10_500; // > 10_000 and > RADIX_THRESHOLD * 2 (512)
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("failed to build rayon thread pool");

        let chunks = pool.install(|| {
            let mut buffer = TemplateRecordBuffer::with_capacity(n, n * 50);

            for i in 0..n {
                // Vary the primary sort field so records are not all equal
                let key = TemplateKey::new(
                    0,
                    (n - i) as i32,
                    false,
                    0,
                    200,
                    false,
                    0,
                    0,
                    (1, true),
                    0,
                    false,
                );
                buffer.push(&make_bam_record(i as u16), key).expect("push should succeed in tests");
            }

            buffer.par_sort_into_chunks(4)
        });

        // With 4 threads and > 10_000 records we should get multiple chunks
        assert!(
            chunks.len() > 1,
            "expected multiple chunks from parallel path, got {}",
            chunks.len()
        );

        // Verify each chunk is individually sorted
        for (ci, chunk) in chunks.iter().enumerate() {
            for i in 1..chunk.len() {
                assert!(chunk[i - 1].0 <= chunk[i].0, "chunk {ci} not sorted at index {i}");
            }
        }

        // Verify total record count across all chunks matches input
        let total: usize = chunks.iter().map(Vec::len).sum();
        assert_eq!(total, n, "total records across chunks should equal input count");
    }

    /// Helper: build a minimal raw BAM record for coordinate sort with a specific position.
    fn make_coordinate_bam_record(tid: i32, pos: i32) -> Vec<u8> {
        let mut record = Vec::with_capacity(40);
        record.extend_from_slice(&tid.to_le_bytes()); // ref_id
        record.extend_from_slice(&pos.to_le_bytes()); // pos
        record.push(2); // name_len = 2 (including null)
        record.push(0); // mapq = 0
        record.extend_from_slice(&0u16.to_le_bytes()); // bin
        record.extend_from_slice(&0u16.to_le_bytes()); // n_cigar_op = 0
        record.extend_from_slice(&0u16.to_le_bytes()); // flags = 0 (forward)
        record.extend_from_slice(&1u32.to_le_bytes()); // seq_len = 1
        record.extend_from_slice(&(-1i32).to_le_bytes()); // mate_ref_id
        record.extend_from_slice(&(-1i32).to_le_bytes()); // mate_pos
        record.extend_from_slice(&0i32.to_le_bytes()); // tlen
        record.push(b'A'); // read name
        record.push(0); // null terminator
        // Pad to at least 34 bytes (bam_fields::flags reads at offset 14-15)
        while record.len() < 40 {
            record.push(0);
        }
        record
    }

    #[test]
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn test_par_sort_into_chunks_coordinate_single_threaded() {
        let nref = 10u32;
        let n = 100;
        let mut buffer = RecordBuffer::with_capacity(n, n * 50, nref);

        for i in 0..n {
            // Reverse order so sorting is non-trivial
            let pos = (n - i) as i32;
            buffer
                .push_coordinate(&make_coordinate_bam_record(0, pos))
                .expect("push_coordinate should succeed in tests");
        }

        let chunks = buffer.par_sort_into_chunks(1);
        assert_eq!(chunks.len(), 1, "single-threaded should produce exactly 1 chunk");
        assert_eq!(chunks[0].len(), n, "the single chunk should contain all records");

        // Verify the chunk is sorted by key
        for i in 1..chunks[0].len() {
            assert!(
                chunks[0][i - 1].0 <= chunks[0][i].0,
                "single chunk should be sorted at index {i}"
            );
        }
    }

    #[test]
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn test_par_sort_into_chunks_coordinate_parallel() {
        let nref = 10u32;
        let n: usize = 10_500; // > 10_000 and > RADIX_THRESHOLD * 2
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("failed to build rayon thread pool");

        let chunks = pool.install(|| {
            let mut buffer = RecordBuffer::with_capacity(n, n * 50, nref);

            for i in 0..n {
                let pos = (n - i) as i32;
                buffer
                    .push_coordinate(&make_coordinate_bam_record(0, pos))
                    .expect("push_coordinate should succeed in tests");
            }

            buffer.par_sort_into_chunks(4)
        });

        // With 4 threads and > 10_000 records we should get multiple chunks
        assert!(
            chunks.len() > 1,
            "expected multiple chunks from parallel path, got {}",
            chunks.len()
        );

        // Verify each chunk is individually sorted
        for (ci, chunk) in chunks.iter().enumerate() {
            for i in 1..chunk.len() {
                assert!(chunk[i - 1].0 <= chunk[i].0, "chunk {ci} not sorted at index {i}");
            }
        }

        // Verify total record count across all chunks matches input
        let total: usize = chunks.iter().map(Vec::len).sum();
        assert_eq!(total, n, "total records across chunks should equal input count");
    }

    mod proptest_msd {
        use super::*;
        use proptest::{prop_assert_eq, proptest};
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        fn make_ref(key: TemplateKey, offset: u64) -> TemplateRecordRef {
            TemplateRecordRef { key, offset, len: 10, padding: 0 }
        }

        fn hash_pair(a: u64, b: u64) -> u64 {
            let mut h = DefaultHasher::new();
            (a, b).hash(&mut h);
            h.finish()
        }

        proptest! {
            /// Oracle test: MSD hybrid sort must produce the same ordering as a
            /// reference sort-by-key on randomized inputs with heavy primary collisions.
            #[test]
            fn msd_sort_matches_reference(
                n_primaries in 1_usize..=4,
                seed in proptest::num::u64::ANY,
            ) {
                let primaries: Vec<u64> = (0..n_primaries)
                    .map(|i| hash_pair(seed, i as u64))
                    .collect();

                // Build 300+ refs (above RADIX_THRESHOLD) with random keys sharing primaries
                let n = 300;
                let mut refs: Vec<TemplateRecordRef> = Vec::with_capacity(n);
                for i in 0..n {
                    let h = hash_pair(seed, (i + n_primaries) as u64);
                    let primary = primaries[i % n_primaries];
                    let key = TemplateKey {
                        primary,
                        secondary: h,
                        cb_hash: h.wrapping_mul(2_654_435_761),
                        tertiary: h.wrapping_mul(40503),
                        name_hash_upper: h.rotate_left(17),
                    };
                    refs.push(make_ref(key, i as u64));
                }

                let mut expected = refs.clone();
                expected.sort_by(|a, b| a.key.cmp(&b.key));

                radix_sort_template_refs(&mut refs);

                for i in 0..n {
                    prop_assert_eq!(
                        refs[i].key, expected[i].key,
                        "Mismatch at index {}: MSD key {:?} != reference {:?}",
                        i, refs[i].key, expected[i].key
                    );
                }
            }

            /// Oracle test with fully random keys (no forced primary collisions).
            #[test]
            fn msd_sort_matches_reference_random_keys(
                seed in proptest::num::u64::ANY,
            ) {
                let n = 500;
                let mut refs: Vec<TemplateRecordRef> = Vec::with_capacity(n);
                for i in 0..n {
                    let h = hash_pair(seed, i as u64);
                    let key = TemplateKey {
                        primary: h,
                        secondary: h.wrapping_mul(6_364_136_223_846_793_005),
                        cb_hash: h.wrapping_mul(2_654_435_761),
                        tertiary: h.wrapping_mul(40503),
                        name_hash_upper: h.rotate_left(17),
                    };
                    refs.push(make_ref(key, i as u64));
                }

                let mut expected = refs.clone();
                expected.sort_by(|a, b| a.key.cmp(&b.key));

                radix_sort_template_refs(&mut refs);

                for i in 0..n {
                    prop_assert_eq!(
                        refs[i].key, expected[i].key,
                        "Mismatch at index {}: MSD key {:?} != reference {:?}",
                        i, refs[i].key, expected[i].key
                    );
                }
            }
        }
    }
}
