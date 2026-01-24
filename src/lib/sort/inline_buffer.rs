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
use std::hash::{Hash, Hasher};
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
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct RecordRef {
    /// Packed primary sort key for fast comparison.
    pub sort_key: u64,
    /// Hash of read name for secondary sort (tie-breaker).
    pub name_hash: u64,
    /// Offset into `RecordBuffer` where record header starts.
    pub offset: u64,
    /// Length of raw BAM data (excluding inline header).
    pub len: u32,
    /// Padding for 8-byte alignment.
    padding: u32,
}

impl PartialEq for RecordRef {
    fn eq(&self, other: &Self) -> bool {
        self.sort_key == other.sort_key && self.name_hash == other.name_hash
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
        self.sort_key.cmp(&other.sort_key).then_with(|| self.name_hash.cmp(&other.name_hash))
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
    /// Pre-computed name hash.
    name_hash: u64,
    /// Length of following raw BAM data.
    record_len: u32,
    /// Padding for 8-byte alignment.
    padding: u32,
}

const HEADER_SIZE: usize = std::mem::size_of::<InlineHeader>(); // 24 bytes

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
        let (sort_key, name_hash) = extract_coordinate_key_inline(record, self.nref);

        // Write inline header (24 bytes)
        let header = InlineHeader { sort_key, name_hash, record_len: len, padding: 0 };

        // Extend data with header bytes
        self.data.extend_from_slice(&header.sort_key.to_le_bytes());
        self.data.extend_from_slice(&header.name_hash.to_le_bytes());
        self.data.extend_from_slice(&header.record_len.to_le_bytes());
        self.data.extend_from_slice(&header.padding.to_le_bytes());

        // Write raw BAM data
        self.data.extend_from_slice(record);

        // Add to index
        self.refs.push(RecordRef { sort_key, name_hash, offset, len, padding: 0 });
    }

    /// Sort the index by key (records stay in place).
    pub fn sort(&mut self) {
        // Use unstable sort for better performance (no stability needed)
        self.refs.sort_unstable();
    }

    /// Sort using parallel sort (for large arrays).
    pub fn par_sort(&mut self) {
        use rayon::prelude::*;
        self.refs.par_sort_unstable();
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
#[inline]
fn extract_coordinate_key_inline(bam: &[u8], nref: u32) -> (u64, u64) {
    // BAM format offsets (all little-endian):
    // 0-3: tid (i32)
    // 4-7: pos (i32)
    // 8: l_read_name (u8)
    // 14-15: flags (u16)
    // 32+: read name (null-terminated)

    let tid = i32::from_le_bytes([bam[0], bam[1], bam[2], bam[3]]);
    let pos = i32::from_le_bytes([bam[4], bam[5], bam[6], bam[7]]);
    let flags = u16::from_le_bytes([bam[14], bam[15]]);
    let reverse = (flags & 0x10) != 0;

    // Pack key based on tid (samtools behavior):
    // - tid >= 0: sort by (tid, pos, reverse) even if unmapped flag is set
    // - tid < 0: unmapped with no reference, sort at end
    let sort_key = if tid < 0 {
        PackedCoordKey::unmapped().0
    } else {
        PackedCoordKey::new(tid, pos, reverse, nref).0
    };

    // Hash read name for tie-breaking
    let name_len = (bam[8] as usize).saturating_sub(1); // Exclude null terminator
    let name =
        if name_len > 0 && 32 + name_len <= bam.len() { &bam[32..32 + name_len] } else { &[] };
    let name_hash = hash_name(name);

    (sort_key, name_hash)
}

/// Fast name hashing using ahash (same hasher as used elsewhere in fgumi).
#[inline]
fn hash_name(name: &[u8]) -> u64 {
    let mut hasher = ahash::AHasher::default();
    name.hash(&mut hasher);
    hasher.finish()
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
    pub fn sort(&mut self) {
        // Keys are cached in refs - no need to access data buffer during sort
        self.refs.sort_unstable_by(|a, b| a.key.cmp(&b.key));
    }

    /// Sort using parallel sort with cached keys.
    pub fn par_sort(&mut self) {
        use rayon::prelude::*;
        // Keys are cached in refs - no need to access data buffer during sort
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
}
