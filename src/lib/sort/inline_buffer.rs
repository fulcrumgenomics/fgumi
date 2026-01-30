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

use crate::sort::keys::{RawSortKey, SortContext};
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
pub struct PackedCoordKey(pub u64);

impl PackedCoordKey {
    /// Create a packed coordinate key.
    ///
    /// # Arguments
    /// * `tid` - Reference sequence ID (-1 for unmapped)
    /// * `pos` - 0-based alignment position
    /// * `reverse` - True if reverse complemented
    /// * `nref` - Number of reference sequences (for unmapped handling)
    #[inline]
    #[must_use]
    pub fn new(tid: i32, pos: i32, reverse: bool, nref: u32) -> Self {
        // Map unmapped (tid=-1) to nref for proper sorting (after all mapped)
        let tid = if tid < 0 { nref } else { tid as u32 };
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
    /// Raw byte storage for all records (headers + BAM data).
    data: Vec<u8>,
    /// Index of record references for sorting.
    refs: Vec<RecordRef>,
    /// Number of reference sequences (for unmapped handling).
    nref: u32,
}

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
            data: Vec::with_capacity(estimated_bytes + estimated_records * HEADER_SIZE),
            refs: Vec::with_capacity(estimated_records),
            nref,
        }
    }

    /// Push a record for coordinate sorting.
    ///
    /// Extracts the sort key inline from raw BAM bytes (zero-copy).
    #[inline]
    pub fn push_coordinate(&mut self, record: &[u8]) {
        let offset = self.data.len() as u64;
        let len = record.len() as u32;

        // Extract sort key from raw BAM bytes
        let sort_key = extract_coordinate_key_inline(record, self.nref);

        // Write inline header (16 bytes)
        let header = InlineHeader { sort_key, record_len: len, padding: 0 };

        // Extend data with header bytes
        self.data.extend_from_slice(&header.sort_key.to_le_bytes());
        self.data.extend_from_slice(&header.record_len.to_le_bytes());
        self.data.extend_from_slice(&header.padding.to_le_bytes());

        // Write raw BAM data
        self.data.extend_from_slice(record);

        // Add to index
        self.refs.push(RecordRef { sort_key, offset, len, padding: 0 });
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

    /// Get record bytes by reference.
    #[inline]
    #[must_use]
    pub fn get_record(&self, r: &RecordRef) -> &[u8] {
        let start = r.offset as usize + HEADER_SIZE;
        let end = start + r.len as usize;
        &self.data[start..end]
    }

    /// Iterate over sorted records.
    pub fn iter_sorted(&self) -> impl Iterator<Item = &[u8]> {
        self.refs.iter().map(|r| self.get_record(r))
    }

    /// Get the sorted record references.
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
fn extract_coordinate_key_inline(bam: &[u8], nref: u32) -> u64 {
    // BAM format offsets (all little-endian):
    // 0-3: tid (i32)
    // 4-7: pos (i32)
    // 14-15: flags (u16)

    let tid = i32::from_le_bytes([bam[0], bam[1], bam[2], bam[3]]);
    let pos = i32::from_le_bytes([bam[4], bam[5], bam[6], bam[7]]);
    let flags = u16::from_le_bytes([bam[14], bam[15]]);
    let reverse = (flags & 0x10) != 0;

    // Pack key based on tid (samtools behavior):
    // - tid >= 0: sort by (tid, pos, reverse) even if unmapped flag is set
    // - tid < 0: unmapped with no reference, sort at end
    if tid < 0 {
        PackedCoordKey::unmapped().0
    } else {
        PackedCoordKey::new(tid, pos, reverse, nref).0
    }
}

// ============================================================================
// Template-Coordinate Support
// ============================================================================

/// Extended key for template-coordinate sorting.
///
/// Template-coordinate requires comparing multiple fields:
/// tid1, tid2, pos1, pos2, neg1, neg2, library, MI, name, `is_upper`
///
/// We pack these into 4 u64 values for efficient comparison.
/// The `name_hash_upper` field packs both `name_hash` and `is_upper`:
/// - Upper 63 bits: name hash (groups same names together)
/// - Lowest bit: `is_upper` (false=0, true=1)
///
/// This ensures reads from the same template stay together (same hash),
/// with `is_upper=false` sorting before `is_upper=true`.
#[repr(C)]
#[derive(Copy, Clone, Debug, Eq, PartialEq, bytemuck::Pod, bytemuck::Zeroable)]
pub struct TemplateKey {
    /// Packed: (tid1 << 48) | (pos1 << 16) | (neg1 << 1)
    pub primary: u64,
    /// Packed: (tid2 << 48) | (pos2 << 16) | neg2
    pub secondary: u64,
    /// Packed: (library << 48) | (`mi_value` << 1) | `mi_suffix`
    pub tertiary: u64,
    /// Packed: (`name_hash` << 1) | `is_upper`
    /// This ensures same-name records group together, with `is_upper` as tie-breaker.
    pub name_hash_upper: u64,
}

impl TemplateKey {
    /// Create a new template key from extracted fields.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        tid1: i32,
        pos1: i32,
        neg1: bool,
        tid2: i32,
        pos2: i32,
        neg2: bool,
        library: u32,
        mi: (u64, bool),
        name_hash: u64,
        is_upper: bool,
    ) -> Self {
        // Handle i32::MAX specially (indicates unmapped mate)
        let tid1_packed = if tid1 == i32::MAX { 0xFFFF_u64 } else { (tid1.max(0) as u64) & 0xFFFF };
        let tid2_packed = if tid2 == i32::MAX { 0xFFFF_u64 } else { (tid2.max(0) as u64) & 0xFFFF };
        let pos1_packed =
            if pos1 == i32::MAX { 0xFFFF_FFFF_u64 } else { (pos1.max(0) as u64) & 0xFFFF_FFFF };
        let pos2_packed =
            if pos2 == i32::MAX { 0xFFFF_FFFF_u64 } else { (pos2.max(0) as u64) & 0xFFFF_FFFF };

        // Pack primary: tid1 (high 16), pos1 (middle 32), neg1 (bit 1)
        let p1 = (tid1_packed << 48) | (pos1_packed << 16) | (u64::from(neg1) << 1);

        // Pack secondary: tid2 (high 16), pos2 (middle 32), neg2 (bit 0)
        let p2 = (tid2_packed << 48) | (pos2_packed << 16) | u64::from(neg2);

        // Pack tertiary: library (high 16), mi_value (middle), mi_suffix (bit 0)
        // Note: /B suffix should sort after /A, so we use !is_a as the bit
        let p3 = ((u64::from(library) & 0xFFFF) << 48)
            | ((mi.0 & 0xFFFF_FFFF_FFFF) << 1)
            | u64::from(!mi.1);

        // Pack name_hash and is_upper: hash in upper 63 bits, is_upper in bit 0
        // This ensures same-name records group together, with is_upper=false before is_upper=true
        let p4 = (name_hash << 1) | u64::from(is_upper);

        Self { primary: p1, secondary: p2, tertiary: p3, name_hash_upper: p4 }
    }

    /// Create a key for completely unmapped records.
    #[must_use]
    pub fn unmapped(name_hash: u64, is_read2: bool) -> Self {
        Self {
            primary: u64::MAX,
            secondary: u64::MAX,
            tertiary: 0,
            name_hash_upper: (name_hash << 1) | u64::from(is_read2),
        }
    }

    /// Create a zeroed key (used as dummy for memory operations).
    #[inline]
    #[must_use]
    pub fn zeroed() -> Self {
        Self { primary: 0, secondary: 0, tertiary: 0, name_hash_upper: 0 }
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
    pub fn to_bytes(&self) -> [u8; 32] {
        let mut buf = [0u8; 32];
        buf[0..8].copy_from_slice(&self.primary.to_le_bytes());
        buf[8..16].copy_from_slice(&self.secondary.to_le_bytes());
        buf[16..24].copy_from_slice(&self.tertiary.to_le_bytes());
        buf[24..32].copy_from_slice(&self.name_hash_upper.to_le_bytes());
        buf
    }

    /// Deserialize from bytes read from keyed temp files.
    #[inline]
    #[must_use]
    pub fn from_bytes(buf: &[u8; 32]) -> Self {
        Self {
            primary: u64::from_le_bytes([
                buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7],
            ]),
            secondary: u64::from_le_bytes([
                buf[8], buf[9], buf[10], buf[11], buf[12], buf[13], buf[14], buf[15],
            ]),
            tertiary: u64::from_le_bytes([
                buf[16], buf[17], buf[18], buf[19], buf[20], buf[21], buf[22], buf[23],
            ]),
            name_hash_upper: u64::from_le_bytes([
                buf[24], buf[25], buf[26], buf[27], buf[28], buf[29], buf[30], buf[31],
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
            .then_with(|| self.tertiary.cmp(&other.tertiary))
            // name_hash_upper comparison handles both name grouping AND is_upper ordering
            .then_with(|| self.name_hash_upper.cmp(&other.name_hash_upper))
    }
}

impl RawSortKey for TemplateKey {
    const SERIALIZED_SIZE: Option<usize> = Some(32);

    fn extract(_bam: &[u8], _ctx: &SortContext) -> Self {
        // Template key extraction requires LibraryLookup which is in raw.rs
        // This method should not be called directly - use extract_template_key_inline() in raw.rs
        // For generic pipeline, we pass the key extraction function as a parameter
        unimplemented!(
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
        let mut buf = [0u8; 32];
        reader.read_exact(&mut buf)?;
        Ok(Self::from_bytes(&buf))
    }
}

/// Inline header stored before each record in the data buffer.
/// This allows us to use a minimal ref structure while still having
/// fast access to the full sort key.
///
/// Layout (40 bytes total):
/// - primary: 8 bytes (u64)
/// - secondary: 8 bytes (u64)
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
pub const TEMPLATE_HEADER_SIZE: usize = 40; // 4 * 8 (key) + 4 + 4

impl TemplateInlineHeader {
    /// Write header to a byte buffer.
    #[inline]
    pub fn write_to(&self, buf: &mut Vec<u8>) {
        buf.extend_from_slice(&self.key.primary.to_le_bytes());
        buf.extend_from_slice(&self.key.secondary.to_le_bytes());
        buf.extend_from_slice(&self.key.tertiary.to_le_bytes());
        buf.extend_from_slice(&self.key.name_hash_upper.to_le_bytes());
        buf.extend_from_slice(&self.record_len.to_le_bytes());
        buf.extend_from_slice(&self.padding.to_le_bytes());
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
        let tertiary = u64::from_le_bytes([
            data[16], data[17], data[18], data[19], data[20], data[21], data[22], data[23],
        ]);
        let name_hash_upper = u64::from_le_bytes([
            data[24], data[25], data[26], data[27], data[28], data[29], data[30], data[31],
        ]);
        let record_len = u32::from_le_bytes([data[32], data[33], data[34], data[35]]);

        Self {
            key: TemplateKey { primary, secondary, tertiary, name_hash_upper },
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
    /// Raw byte storage: inline headers + record data.
    data: Vec<u8>,
    /// Minimal index for sorting.
    refs: Vec<TemplateRecordRef>,
}

impl TemplateRecordBuffer {
    /// Create a new buffer with estimated capacity.
    #[must_use]
    pub fn with_capacity(estimated_records: usize, estimated_bytes: usize) -> Self {
        // Account for inline headers in data capacity
        let header_bytes = estimated_records * TEMPLATE_HEADER_SIZE;
        Self {
            data: Vec::with_capacity(estimated_bytes + header_bytes),
            refs: Vec::with_capacity(estimated_records),
        }
    }

    /// Push a record with a pre-computed template key.
    #[inline]
    pub fn push(&mut self, record: &[u8], key: TemplateKey) {
        let offset = self.data.len() as u64;
        let record_len = record.len() as u32;

        // Write inline header using manual byte operations (avoids alignment issues)
        let header = TemplateInlineHeader { key, record_len, padding: 0 };
        header.write_to(&mut self.data);

        // Write raw BAM data
        self.data.extend_from_slice(record);

        // Add ref with cached key for O(1) sort comparisons
        self.refs.push(TemplateRecordRef { key, offset, len: record_len, padding: 0 });
    }

    /// Sort the index by cached key (O(1) comparison, no data buffer access).
    ///
    /// Note: Uses comparison sort rather than radix sort because `TemplateKey`
    /// is 32 bytes (4 × u64), making multi-field radix sort slower than
    /// comparison sort due to the many passes required.
    pub fn sort(&mut self) {
        self.refs.sort_unstable_by(|a, b| a.key.cmp(&b.key));
    }

    /// Sort using parallel sort with cached keys.
    pub fn par_sort(&mut self) {
        use rayon::prelude::*;
        self.refs.par_sort_unstable_by(|a, b| a.key.cmp(&b.key));
    }

    /// Get record bytes by reference.
    #[inline]
    #[must_use]
    pub fn get_record(&self, r: &TemplateRecordRef) -> &[u8] {
        let start = r.offset as usize + TEMPLATE_HEADER_SIZE;
        let end = start + r.len as usize;
        &self.data[start..end]
    }

    /// Iterate over sorted records.
    pub fn iter_sorted(&self) -> impl Iterator<Item = &[u8]> {
        self.refs.iter().map(|r| self.get_record(r))
    }

    /// Get the sorted record references.
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

    /// Get underlying data buffer (for direct access to raw bytes).
    #[must_use]
    pub fn data(&self) -> &[u8] {
        &self.data
    }
}

// ============================================================================
// Radix Sort for RecordRef
// ============================================================================

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
#[allow(clippy::uninit_vec)]
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
        for r in src_slice.iter() {
            let byte = ((r.sort_key >> (byte_idx * 8)) & 0xFF) as usize;
            counts[byte] += 1;
        }

        // Convert to cumulative offsets
        let mut total = 0;
        for count in counts.iter_mut() {
            let c = *count;
            *count = total;
            total += c;
        }

        // Scatter elements to destination (stable - preserves order within buckets)
        for r in src_slice.iter() {
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

    if chunk_ranges.len() <= 1 {
        return;
    }

    let n = refs.len();
    let mut result: Vec<RecordRef> = Vec::with_capacity(n);

    // Heap entry: (sort_key, chunk_idx, position_in_chunk)
    struct HeapEntry {
        key: u64,
        chunk_idx: usize,
        pos: usize,
    }

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
    // Use chunk_idx as tie-breaker for stability
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

/// Radix sort for `TemplateRecordRef` arrays using multi-field LSD approach.
///
/// The `TemplateKey` consists of 4 u64 fields sorted in order:
/// primary → secondary → tertiary → `name_hash_upper`
///
/// For LSD radix sort, we sort from least significant to most significant:
/// 1. Sort by `name_hash_upper`
/// 2. Sort by tertiary (stable, preserves `name_hash_upper` order)
/// 3. Sort by secondary (stable, preserves tertiary + `name_hash_upper` order)
/// 4. Sort by primary (stable, final order)
///
/// This is O(n×k) where k is the total bytes to sort (up to 32 bytes).
///
/// Note: Currently unused because benchmarks showed multi-field radix sort
/// provides minimal benefit over comparison sort for 32-byte keys.
/// Kept for potential future optimization experiments.
#[allow(dead_code)]
#[allow(clippy::uninit_vec)]
pub fn radix_sort_template_refs(refs: &mut [TemplateRecordRef]) {
    let n = refs.len();
    if n < RADIX_THRESHOLD {
        insertion_sort_template_refs(refs);
        return;
    }

    // Allocate auxiliary buffer
    let mut aux: Vec<TemplateRecordRef> = Vec::with_capacity(n);
    unsafe {
        aux.set_len(n);
    }

    // Sort by each field from least significant to most significant
    // This ensures the final order is: primary → secondary → tertiary → name_hash_upper
    let fields = [
        |r: &TemplateRecordRef| r.key.name_hash_upper,
        |r: &TemplateRecordRef| r.key.tertiary,
        |r: &TemplateRecordRef| r.key.secondary,
        |r: &TemplateRecordRef| r.key.primary,
    ];

    for (field_idx, get_field) in fields.iter().enumerate() {
        // Find max value for this field to determine bytes needed
        let max_val = refs.iter().map(get_field).max().unwrap_or(0);
        let bytes_needed =
            if max_val == 0 { 0 } else { ((64 - max_val.leading_zeros()) as usize).div_ceil(8) };

        if bytes_needed == 0 {
            continue; // All values are 0 for this field, skip
        }

        // Radix sort this field
        radix_sort_template_field(refs, &mut aux, get_field, bytes_needed, field_idx);
    }
}

/// Radix sort a single u64 field of `TemplateRecordRef` using raw pointers.
#[allow(dead_code)]
#[allow(clippy::uninit_vec)]
fn radix_sort_template_field<F>(
    refs: &mut [TemplateRecordRef],
    aux: &mut [TemplateRecordRef],
    get_field: F,
    bytes_needed: usize,
    _field_idx: usize,
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
        for r in src_slice.iter() {
            let byte = ((get_field(r) >> (byte_idx * 8)) & 0xFF) as usize;
            counts[byte] += 1;
        }

        // Convert to cumulative offsets
        let mut total = 0;
        for count in counts.iter_mut() {
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
#[allow(dead_code)]
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
#[allow(dead_code)]
fn merge_sorted_template_chunks(
    refs: &mut [TemplateRecordRef],
    chunk_ranges: &[std::ops::Range<usize>],
) {
    use crate::sort::radix::{heap_make, heap_sift_down};

    if chunk_ranges.len() <= 1 {
        return;
    }

    let n = refs.len();
    let mut result: Vec<TemplateRecordRef> = Vec::with_capacity(n);

    struct HeapEntry {
        key: TemplateKey,
        chunk_idx: usize,
        pos: usize,
    }

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
#[allow(dead_code)]
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
        assert!(PackedCoordKey::new(0, 100, false, 10) < PackedCoordKey::new(1, 100, false, 10));

        // Lower pos should come first
        assert!(PackedCoordKey::new(0, 100, false, 10) < PackedCoordKey::new(0, 200, false, 10));

        // Forward should come before reverse (false < true)
        assert!(PackedCoordKey::new(0, 100, false, 10) < PackedCoordKey::new(0, 100, true, 10));

        // Unmapped should come last
        assert!(PackedCoordKey::new(9, 1_000_000, true, 10) < PackedCoordKey::unmapped());
    }

    #[test]
    fn test_template_key_ordering() {
        let k1 = TemplateKey::new(0, 100, false, 0, 200, false, 0, (1, true), 0, false);
        let k2 = TemplateKey::new(0, 100, false, 0, 200, false, 0, (2, true), 0, false);
        assert!(k1 < k2);

        // /A suffix should come before /B
        let ka = TemplateKey::new(0, 100, false, 0, 200, false, 0, (1, true), 0, false);
        let kb = TemplateKey::new(0, 100, false, 0, 200, false, 0, (1, false), 0, false);
        assert!(ka < kb);

        // Same name hash: is_upper=false should come before is_upper=true
        let lower = TemplateKey::new(0, 100, false, 0, 200, false, 0, (1, true), 12345, false);
        let upper = TemplateKey::new(0, 100, false, 0, 200, false, 0, (1, true), 12345, true);
        assert!(lower < upper, "is_upper=false should sort before is_upper=true");

        // Different name hashes should group separately
        let first_hash_lo =
            TemplateKey::new(0, 100, false, 0, 200, false, 0, (1, true), 100, false);
        let first_hash_hi = TemplateKey::new(0, 100, false, 0, 200, false, 0, (1, true), 100, true);
        let second_hash = TemplateKey::new(0, 100, false, 0, 200, false, 0, (1, true), 200, false);
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
            assert_eq!(r.sort_key, i as u64, "Expected sort_key {} at index {}", i, i);
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
    fn test_radix_sort_all_same_keys() {
        // All same keys should maintain original order (stability)
        let mut refs: Vec<RecordRef> = (0..100)
            .map(|i| RecordRef { sort_key: 42, offset: i as u64 * 100, len: 10, padding: 0 })
            .collect();

        radix_sort_record_refs(&mut refs);

        // Verify all keys are 42 and order is preserved
        for (i, r) in refs.iter().enumerate() {
            assert_eq!(r.sort_key, 42);
            assert_eq!(r.offset, i as u64 * 100, "Stability violated at index {}", i);
        }
    }

    #[test]
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
            assert_eq!(r.sort_key, i as u64, "Expected sort_key {} at index {}", i, i);
        }
    }

    #[test]
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
            assert!(refs[i - 1].sort_key <= refs[i].sort_key, "Not sorted at index {}", i);
            // Within same key group, offsets should be in order
            if refs[i - 1].sort_key == refs[i].sort_key {
                assert!(refs[i - 1].offset < refs[i].offset, "Stability violated at index {}", i);
            }
        }
    }
}
