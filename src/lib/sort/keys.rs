//! Sort key types for BAM record sorting.
//!
//! This module provides lightweight sort keys that can be extracted from BAM records
//! with minimal parsing overhead. Keys are designed for fast comparison and minimal
//! memory footprint.
//!
//! # Key Types
//!
//! - [`CoordinateKey`]: Standard genomic coordinate (tid, pos, strand)
//! - [`QuerynameKey`]: Read name with natural numeric ordering
//! - [`TemplateKey`](super::inline_buffer::TemplateKey): Template-level position for UMI grouping
//!
//! # Generic Sorting Abstraction
//!
//! The [`RawSortKey`] trait provides a unified interface for sort keys that can be:
//! - Extracted directly from raw BAM bytes (zero-copy)
//! - Serialized/deserialized for temp file storage
//! - Compared efficiently (O(1) for fixed-size keys)
//!
//! This design is inspired by:
//! - fgbio's `SamOrder` trait (Scala)
//! - samtools' `bam1_tag` union (C)

use anyhow::Result;
use noodles::sam::Header;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::cmp::Ordering;

use crate::sort::bam_fields;
use std::io::{Read, Write};

/// The PA (Primary Alignment) tag for secondary/supplementary reads.
/// Stores template-coordinate sort key as B:i array (6 int32s).
/// Added by `fgumi zipper` to enable proper template-coordinate sorting.
pub const PA_TAG: Tag = Tag::new(b'p', b'a');

// ============================================================================
// Primary Alignment Info (PA tag)
// ============================================================================

/// Primary alignment info for secondary/supplementary reads.
///
/// Stores the template-coordinate sort key from the primary alignments,
/// enabling secondary/supplementary reads to sort adjacent to their primaries.
///
/// # Binary Format
///
/// Stored as a `B:i` (int32 array) BAM tag with 6 elements:
/// `[tid1, pos1, neg1, tid2, pos2, neg2]`
/// where `neg1`/`neg2` are 0 for forward, 1 for reverse strand.
///
/// This format is faster to parse than a string representation and ensures
/// supplementary reads get the exact same sort key as their primary reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PrimaryAlignmentInfo {
    /// Reference ID of the earlier mate (lower position).
    pub tid1: i32,
    /// Unclipped 5' position of the earlier mate.
    pub pos1: i32,
    /// True if earlier mate is on reverse strand.
    pub neg1: bool,
    /// Reference ID of the later mate.
    pub tid2: i32,
    /// Unclipped 5' position of the later mate.
    pub pos2: i32,
    /// True if later mate is on reverse strand.
    pub neg2: bool,
}

impl PrimaryAlignmentInfo {
    /// Creates a new `PrimaryAlignmentInfo`.
    #[must_use]
    pub const fn new(tid1: i32, pos1: i32, neg1: bool, tid2: i32, pos2: i32, neg2: bool) -> Self {
        Self { tid1, pos1, neg1, tid2, pos2, neg2 }
    }

    /// Serializes to a BAM tag value (B:i array with 6 elements).
    #[must_use]
    pub fn to_tag_value(&self) -> noodles::sam::alignment::record_buf::data::field::Value {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let values: Vec<i32> = vec![
            self.tid1,
            self.pos1,
            i32::from(self.neg1),
            self.tid2,
            self.pos2,
            i32::from(self.neg2),
        ];
        Value::Array(Array::Int32(values))
    }

    /// Deserializes from a BAM tag value (B:i array with 6 int32 elements).
    ///
    /// Optimized with a fast path for Int32 arrays (the expected format) that
    /// avoids heap allocation by directly indexing the array.
    #[must_use]
    pub fn from_tag_value(
        value: &noodles::sam::alignment::record_buf::data::field::Value,
    ) -> Option<Self> {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        match value {
            Value::Array(arr) => {
                // Fast path: Int32 array (expected format from zipper) - no allocation
                if let Array::Int32(v) = arr {
                    if v.len() == 6 {
                        return Some(Self {
                            tid1: v[0],
                            pos1: v[1],
                            neg1: v[2] != 0,
                            tid2: v[3],
                            pos2: v[4],
                            neg2: v[5] != 0,
                        });
                    }
                    return None;
                }

                // Slow path: other array types (rare) - requires allocation
                let values: Vec<i32> = match arr {
                    Array::Int8(v) => v.iter().map(|&x| i32::from(x)).collect(),
                    Array::UInt8(v) => v.iter().map(|&x| i32::from(x)).collect(),
                    Array::Int16(v) => v.iter().map(|&x| i32::from(x)).collect(),
                    Array::UInt16(v) => v.iter().map(|&x| i32::from(x)).collect(),
                    Array::Int32(_) => unreachable!(), // Handled above
                    Array::UInt32(v) => {
                        // Use try_from to avoid wrapping for values > i32::MAX
                        let result: Result<Vec<i32>, _> =
                            v.iter().map(|&x| i32::try_from(x)).collect();
                        match result {
                            Ok(vals) => vals,
                            Err(_) => return None,
                        }
                    }
                    Array::Float(_) => return None,
                };

                if values.len() != 6 {
                    return None;
                }

                Some(Self {
                    tid1: values[0],
                    pos1: values[1],
                    neg1: values[2] != 0,
                    tid2: values[3],
                    pos2: values[4],
                    neg2: values[5] != 0,
                })
            }
            _ => None,
        }
    }

    /// Extracts from a BAM record's PA tag.
    #[must_use]
    pub fn from_record(record: &RecordBuf) -> Option<Self> {
        let value = record.data().get(&PA_TAG)?;
        Self::from_tag_value(value)
    }
}

// ============================================================================
// Generic Sorting Abstraction (Trait-based, inspired by fgbio/samtools)
// ============================================================================

/// Trait for sort keys extracted from raw BAM bytes.
///
/// This trait provides a unified interface for all sort key types, enabling:
/// - Zero-copy key extraction from raw BAM bytes
/// - Efficient serialization for temp file storage
/// - O(1) comparisons during merge phase (for fixed-size keys)
///
/// # Design
///
/// Inspired by fgbio's `SamOrder` trait and samtools' `bam1_tag` approach.
/// Using a trait with monomorphization (`sort_with_keyed<K>`) gives:
/// - Zero-cost abstraction (no runtime dispatch)
/// - Type safety (can't mix key types)
/// - Compile-time optimization (inlining, etc.)
pub trait RawSortKey: Ord + Clone + Send + Sync + Sized {
    /// Fixed byte size when serialized, or `None` for variable-length keys.
    ///
    /// Fixed-size keys enable O(1) reads from temp files during merge.
    const SERIALIZED_SIZE: Option<usize>;

    /// When `true`, the sort key is embedded in the BAM record bytes at a known
    /// offset, so temp files store only raw records without a key prefix.
    /// This saves I/O for variable-length keys (e.g. queryname) by avoiding
    /// writing the name twice — once as key prefix, once inside the BAM record.
    const EMBEDDED_IN_RECORD: bool = false;

    /// Extract a sort key from raw BAM record bytes.
    ///
    /// This is the hot path during sorting - implementations should minimize
    /// parsing overhead by reading only the fields needed for comparison.
    fn extract(bam: &[u8], ctx: &SortContext) -> Self;

    /// Extract a sort key from raw BAM record bytes without a `SortContext`.
    ///
    /// Only valid when `EMBEDDED_IN_RECORD` is `true`. Used during merge to
    /// reconstruct the key from the record itself, avoiding separate key storage.
    ///
    /// Default implementation panics — only override for embedded key types.
    #[must_use]
    fn extract_from_record(_bam: &[u8]) -> Self {
        unimplemented!("extract_from_record only valid when EMBEDDED_IN_RECORD is true")
    }

    /// Serialize the key to a writer for temp file storage.
    ///
    /// Format should be compact and enable fast deserialization.
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the writer fails.
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()>;

    /// Deserialize a key from a reader.
    ///
    /// Must be the inverse of `write_to`.
    ///
    /// # Errors
    ///
    /// Returns an error if reading from the reader fails or data is invalid.
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self>;
}

/// Context needed for sort key extraction (built from BAM header).
///
/// This struct holds header-derived information needed by some sort orders:
/// - `nref`: Number of reference sequences (for coordinate sort unmapped handling)
/// - `lib_lookup`: Library name ordinals (for template-coordinate sort)
#[derive(Clone)]
pub struct SortContext {
    /// Number of reference sequences (unmapped reads map to nref).
    pub nref: u32,
    // Note: LibraryLookup is in raw.rs - we'll use a callback instead
}

impl SortContext {
    /// Create a sort context from a BAM header.
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub fn from_header(header: &Header) -> Self {
        Self { nref: header.reference_sequences().len() as u32 }
    }

    /// Create a context with explicit nref (for testing).
    #[must_use]
    pub fn new(nref: u32) -> Self {
        Self { nref }
    }
}

/// Queryname comparison strategy.
///
/// The SAM spec allows queryname sort with different sub-sort orders
/// specified via the `SS` header tag.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QuerynameComparator {
    /// Lexicographic byte ordering (fast, default).
    ///
    /// Standard byte-by-byte comparison. This is 8-13x faster than natural
    /// ordering and sufficient for all downstream tools that require
    /// queryname-grouped input (e.g., `fgumi zipper`, `fgumi group`).
    #[default]
    Lexicographic,
    /// Natural numeric ordering (samtools-compatible).
    ///
    /// Handles embedded numbers naturally: "read2" < "read10".
    /// Use when output must match `samtools sort -n` ordering.
    Natural,
}

impl std::fmt::Display for QuerynameComparator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Lexicographic => write!(f, "lexicographic"),
            Self::Natural => write!(f, "natural"),
        }
    }
}

/// Sort order enumeration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    /// Coordinate sort: tid → pos → reverse strand
    Coordinate,
    /// Queryname sort: read name ordering with configurable comparator.
    Queryname(QuerynameComparator),
    /// Template-coordinate sort: template position for UMI grouping
    TemplateCoordinate,
}

impl SortOrder {
    /// Returns true if this is any queryname sort order.
    #[must_use]
    pub fn is_queryname(&self) -> bool {
        matches!(self, Self::Queryname(_))
    }

    /// Returns the queryname comparator, if this is a queryname sort order.
    #[must_use]
    pub fn queryname_comparator(&self) -> Option<QuerynameComparator> {
        match self {
            Self::Queryname(cmp) => Some(*cmp),
            _ => None,
        }
    }

    /// Get the SAM header sort order tag value.
    #[must_use]
    pub fn header_so_tag(&self) -> &'static str {
        match self {
            Self::Coordinate => "coordinate",
            Self::Queryname(_) => "queryname",
            Self::TemplateCoordinate => "unsorted",
        }
    }

    /// Get the SAM header group order tag value.
    #[must_use]
    pub fn header_go_tag(&self) -> Option<&'static str> {
        match self {
            Self::Coordinate | Self::Queryname(_) => None,
            Self::TemplateCoordinate => Some("query"),
        }
    }

    /// Get the SAM header sub-sort tag value.
    #[must_use]
    pub fn header_ss_tag(&self) -> Option<&'static str> {
        match self {
            Self::Coordinate => None,
            Self::Queryname(QuerynameComparator::Lexicographic) => Some("lexicographic"),
            Self::Queryname(QuerynameComparator::Natural) => Some("natural"),
            Self::TemplateCoordinate => Some("template-coordinate"),
        }
    }
}

/// Trait for sort keys that can be extracted from BAM records.
///
/// The `Context` associated type allows sort keys to pre-compute data from the header
/// once (e.g., library index mapping) and reuse it for each record extraction.
pub trait SortKey: Ord + Clone + Send + Sync {
    /// Context built once from header (e.g., `LibraryIndex` for template-coordinate).
    /// Use `()` for keys that don't need context.
    type Context: Clone + Send + Sync;

    /// Build context from header (called once at sort start).
    fn build_context(header: &Header) -> Self::Context;

    /// Extract a sort key from a BAM record using pre-built context.
    ///
    /// # Errors
    ///
    /// Returns an error if key extraction fails due to missing or invalid fields.
    fn from_record(record: &RecordBuf, header: &Header, ctx: &Self::Context) -> Result<Self>;
}

// ============================================================================
// Coordinate Sort Key
// ============================================================================

/// Sort key for coordinate ordering.
///
/// Sort order: reference ID → position → reverse strand flag.
/// Unmapped reads (tid = -1) are sorted to the end.
///
/// Note: No read name tie-breaking is used, matching samtools behavior.
/// Equal records maintain their original input order (stable sort).
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct CoordinateKey {
    /// Reference sequence ID (tid), or `i32::MAX` for unmapped.
    pub tid: i32,
    /// 0-based alignment start position.
    pub pos: i64,
    /// True if reverse strand.
    pub reverse: bool,
}

impl CoordinateKey {
    /// Create a coordinate key for an unmapped read.
    #[must_use]
    pub fn unmapped() -> Self {
        Self { tid: i32::MAX, pos: i64::MAX, reverse: false }
    }
}

impl Ord for CoordinateKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.tid
            .cmp(&other.tid)
            .then_with(|| self.pos.cmp(&other.pos))
            .then_with(|| self.reverse.cmp(&other.reverse))
    }
}

impl PartialOrd for CoordinateKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl SortKey for CoordinateKey {
    type Context = ();

    fn build_context(_header: &Header) -> Self::Context {}

    fn from_record(record: &RecordBuf, _header: &Header, _ctx: &Self::Context) -> Result<Self> {
        if record.flags().is_unmapped() {
            return Ok(Self::unmapped());
        }

        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let tid = record.reference_sequence_id().map_or(-1, |id| id as i32);

        #[allow(clippy::cast_possible_wrap)]
        let pos = record.alignment_start().map_or(0, |p| usize::from(p) as i64);

        let reverse = record.flags().is_reverse_complemented();

        Ok(Self { tid, pos, reverse })
    }
}

// ============================================================================
// Raw Coordinate Sort Key (Fixed-size for RawSortKey trait)
// ============================================================================

/// Fixed-size coordinate sort key for raw BAM sorting (8 bytes).
///
/// This key is designed for efficient temp file storage and O(1) comparisons
/// during merge phase. It packs:
/// - `sort_key`: (tid << 34) | ((pos+1) << 1) | reverse
///
/// Note: No read name tie-breaking is used, matching samtools behavior.
/// Equal records maintain their original input order (stable sort).
#[repr(C)]
#[derive(Copy, Clone, Eq, PartialEq, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct RawCoordinateKey {
    /// Packed primary sort key: (tid << 34) | ((pos+1) << 1) | reverse.
    pub sort_key: u64,
}

impl RawCoordinateKey {
    /// Size in bytes when serialized.
    pub const SIZE: usize = 8;

    /// Create a new coordinate key from components.
    ///
    /// # Arguments
    /// * `tid` - Reference sequence ID (-1 for unmapped)
    /// * `pos` - 0-based alignment position
    /// * `reverse` - True if reverse complemented
    /// * `nref` - Number of reference sequences (for unmapped handling)
    #[inline]
    #[must_use]
    #[allow(clippy::cast_sign_loss)]
    pub fn new(tid: i32, pos: i32, reverse: bool, nref: u32) -> Self {
        // Map unmapped (tid=-1) to nref for proper sorting (after all mapped)
        let tid = if tid < 0 { nref } else { tid as u32 };
        // Pack: tid in high bits, (pos+1) in middle, reverse in LSB
        let key = (u64::from(tid) << 34)
            | ((i64::from(pos) as u64).wrapping_add(1) << 1)
            | u64::from(reverse);
        Self { sort_key: key }
    }

    /// Create a key for unmapped records (sorts after all mapped).
    #[inline]
    #[must_use]
    pub fn unmapped() -> Self {
        Self { sort_key: u64::MAX }
    }

    /// Create a zeroed key (for memory operations).
    #[inline]
    #[must_use]
    pub fn zeroed() -> Self {
        Self { sort_key: 0 }
    }
}

impl Default for RawCoordinateKey {
    fn default() -> Self {
        Self::zeroed()
    }
}

impl Ord for RawCoordinateKey {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        self.sort_key.cmp(&other.sort_key)
    }
}

impl PartialOrd for RawCoordinateKey {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl RawSortKey for RawCoordinateKey {
    const SERIALIZED_SIZE: Option<usize> = Some(Self::SIZE);
    const EMBEDDED_IN_RECORD: bool = true;

    #[inline]
    fn extract(bam: &[u8], ctx: &SortContext) -> Self {
        let tid = bam_fields::ref_id(bam);
        let pos = bam_fields::pos(bam);
        let reverse = (bam_fields::flags(bam) & bam_fields::flags::REVERSE) != 0;

        // Create key based on tid (samtools behavior):
        // - tid >= 0: sort by (tid, pos, reverse) even if unmapped flag is set
        // - tid < 0: unmapped with no reference, sort at end
        if tid < 0 { Self::unmapped() } else { Self::new(tid, pos, reverse, ctx.nref) }
    }

    #[inline]
    fn extract_from_record(bam: &[u8]) -> Self {
        let tid = bam_fields::ref_id(bam);
        if tid < 0 {
            return Self::unmapped();
        }
        let pos = bam_fields::pos(bam);
        let reverse = (bam_fields::flags(bam) & bam_fields::flags::REVERSE) != 0;
        // During merge we don't have nref, but for mapped records (tid >= 0)
        // nref is only used for unmapped handling, so any value > tid works.
        // Use tid+1 as a safe nref since tid is non-negative.
        #[allow(clippy::cast_sign_loss)]
        Self::new(tid, pos, reverse, (tid as u32) + 1)
    }

    #[inline]
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writer.write_all(&self.sort_key.to_le_bytes())
    }

    #[inline]
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let mut buf = [0u8; 8];
        reader.read_exact(&mut buf)?;
        Ok(Self { sort_key: u64::from_le_bytes(buf) })
    }
}

// ============================================================================
// Queryname Sort Key
// ============================================================================

/// Transform flags for queryname sort ordering.
///
/// Matches samtools' queryname sort order (`bam_sort.c` lines 245-248):
/// ```c
/// // Sort order is READ1, READ2, (PRIMARY), SUPPLEMENTARY, SECONDARY
/// fa = ((fa&0xc0)<<8)|((fa&0x100)<<3)|((fa&0x800)>>3);
/// ```
///
/// This transforms the relevant flag bits into a value that sorts correctly:
/// - R1 (0x40) and R2 (0x80) bits are shifted to high positions (bits 14-15)
/// - SECONDARY (0x100) is shifted to middle (bit 11)
/// - SUPPLEMENTARY (0x800) is shifted to low position (bit 8)
///
/// The resulting sort order is:
/// 1. NONE (unpaired) - 0x0000
/// 2. R1 PRIMARY      - 0x4000
/// 3. R1 SUPPLEMENTARY - 0x4100
/// 4. R1 SECONDARY    - 0x4800
/// 5. R2 PRIMARY      - 0x8000
/// 6. R2 SUPPLEMENTARY - 0x8100
/// 7. R2 SECONDARY    - 0x8800
#[inline]
#[must_use]
pub const fn queryname_flag_order(flags: u16) -> u16 {
    ((flags & 0xc0) << 8) | ((flags & 0x100) << 3) | ((flags & 0x800) >> 3)
}

/// Sort key for queryname ordering.
///
/// Uses natural string ordering where numeric runs are compared numerically.
/// Example: "read1" < "read2" < "read10" < "read11"
///
/// Names are stored with a null terminator so that `natural_compare_nul` can
/// walk raw pointers without per-byte bounds checks, matching the performance
/// characteristics of samtools' `strnum_cmp`.
#[derive(Clone, Debug)]
pub struct QuerynameKey {
    /// Read name bytes, null-terminated for `natural_compare_nul`.
    name: Vec<u8>,
    /// Read pair flags for ordering R1 before R2.
    flags: u16,
}

impl PartialEq for QuerynameKey {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for QuerynameKey {}

impl QuerynameKey {
    /// Returns the read name bytes (including the null terminator).
    #[must_use]
    pub fn name(&self) -> &[u8] {
        &self.name
    }
}

impl Ord for QuerynameKey {
    #[allow(unsafe_code)]
    fn cmp(&self, other: &Self) -> Ordering {
        // SAFETY: `name` is always null-terminated (see `from_record`).
        unsafe { natural_compare_nul(self.name.as_ptr(), other.name.as_ptr()) }
            .then_with(|| self.flags.cmp(&other.flags))
    }
}

impl PartialOrd for QuerynameKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl SortKey for QuerynameKey {
    type Context = ();

    fn build_context(_header: &Header) -> Self::Context {}

    fn from_record(record: &RecordBuf, _header: &Header, _ctx: &Self::Context) -> Result<Self> {
        let raw = record.name().map_or(&[] as &[u8], |n| <_ as AsRef<[u8]>>::as_ref(n));
        let mut name = Vec::with_capacity(raw.len() + 1);
        name.extend_from_slice(raw);
        name.push(0); // null-terminate for pointer-walking comparison
        let flags = queryname_flag_order(u16::from(record.flags()));
        Ok(Self { name, flags })
    }
}

// Re-export natural comparison functions from fgumi-raw-bam where the unsafe
// implementations live (the main crate enforces `#![deny(unsafe_code)]`).
pub use fgumi_raw_bam::sort::{natural_compare, natural_compare_nul};

// ============================================================================
// Natural Sort Key Normalization
// ============================================================================

/// Encode a read name into a byte string whose lexicographic (`memcmp`) ordering
/// matches `natural_compare` ordering for names that differ in numeric magnitude.
///
/// **Encoding:** Walk the name left to right. Text runs are copied verbatim.
/// Each numeric run is encoded as a length byte (digit count) followed by the
/// significant digits (leading zeros stripped). The length byte sorts first, so
/// `"2"` (length=1) sorts before `"10"` (length=2), matching natural ordering.
///
/// This produces a compact key (one extra byte per numeric field) that can be
/// compared with plain byte comparison instead of the expensive `natural_compare`
/// byte-by-byte scanning with digit detection.
///
/// **Leading-zero caveat:** Because leading zeros are stripped, names that differ
/// only in leading zeros (e.g. `b"01"` vs `b"1"`) normalize to identical bytes.
/// `natural_compare` distinguishes them via a total-length tiebreaker. Callers
/// that need a total order must fall back to `natural_compare` when normalized
/// keys are equal.
///
/// **Not used in production.** The sort pipeline uses [`natural_compare_nul`]
/// (via [`RawQuerynameKey`]) instead. This function is exported for benchmarks
/// that compare sorting strategies.
///
/// # Examples
///
/// ```text
/// "read2"    → [r, e, a, d, 0x01, '2']          (6 bytes)
/// "read10"   → [r, e, a, d, 0x02, '1', '0']     (7 bytes)
/// "SRR123.5" → [S, R, R, 0x03, '1', '2', '3', '.', 0x01, '5']
/// ```
pub fn normalize_natural_key(name: &[u8], out: &mut Vec<u8>) {
    let mut i = 0;
    while i < name.len() {
        if name[i].is_ascii_digit() {
            // Skip leading zeros (matches samtools natural_compare behavior).
            while i < name.len() && name[i] == b'0' {
                i += 1;
            }
            // Count significant digits.
            let sig_start = i;
            while i < name.len() && name[i].is_ascii_digit() {
                i += 1;
            }
            let sig_len = i - sig_start;
            if sig_len == 0 {
                // The run was all zeros — encode as the number 0 (length=1, digit='0').
                out.push(1);
                out.push(b'0');
            } else {
                // Length byte capped at 255; names with 255+ digit numbers are pathological.
                #[allow(clippy::cast_possible_truncation)]
                let len_byte = sig_len.min(255) as u8;
                out.push(len_byte);
                out.extend_from_slice(&name[sig_start..sig_start + sig_len.min(255)]);
            }
        } else {
            out.push(name[i]);
            i += 1;
        }
    }
}

// ============================================================================
// Shared queryname key helpers
// ============================================================================

/// Extract raw queryname bytes and flag-order value from BAM record bytes.
///
/// BAM format offsets: 8 = `l_read_name` (u8), 14-15 = flags (u16), 32+ = name.
#[inline]
fn extract_raw_name_and_flags(bam: &[u8]) -> (&[u8], u16) {
    let name_len = (bam[8] as usize).saturating_sub(1);
    let name =
        if name_len > 0 && 32 + name_len <= bam.len() { &bam[32..32 + name_len] } else { &[] };
    let raw_flags = u16::from_le_bytes([bam[14], bam[15]]);
    (name, queryname_flag_order(raw_flags))
}

/// Serialize a queryname key as `[name_len: u16][name: bytes][flags: u16]`.
#[inline]
fn write_queryname_key<W: Write>(name: &[u8], flags: u16, writer: &mut W) -> std::io::Result<()> {
    let name_len = u16::try_from(name.len()).map_err(|_| {
        std::io::Error::new(std::io::ErrorKind::InvalidInput, "queryname too long for u16")
    })?;
    writer.write_all(&name_len.to_le_bytes())?;
    writer.write_all(name)?;
    writer.write_all(&flags.to_le_bytes())
}

/// Deserialize a queryname key from `[name_len: u16][name: bytes][flags: u16]`.
#[inline]
fn read_queryname_key<R: Read>(reader: &mut R) -> std::io::Result<(Vec<u8>, u16)> {
    let mut len_buf = [0u8; 2];
    reader.read_exact(&mut len_buf)?;
    let name_len = u16::from_le_bytes(len_buf) as usize;

    let mut name = vec![0u8; name_len];
    reader.read_exact(&mut name)?;

    let mut flags_buf = [0u8; 2];
    reader.read_exact(&mut flags_buf)?;
    Ok((name, u16::from_le_bytes(flags_buf)))
}

// ============================================================================
// Raw Queryname Sort Key (Variable-size for RawSortKey trait)
// ============================================================================

/// Variable-size queryname sort key for raw BAM sorting.
///
/// This key stores the full read name for correct natural ordering.
/// Unlike fixed-size keys, serialization size depends on name length.
///
/// Serialization format: `[name_len: u16][name: bytes][flags: u16]`
#[derive(Clone, Debug, Default)]
pub struct RawQuerynameKey {
    /// Read name bytes, null-terminated for `natural_compare_nul`.
    name: Vec<u8>,
    /// Flags for segment ordering (R1 before R2).
    flags: u16,
}

impl PartialEq for RawQuerynameKey {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for RawQuerynameKey {}

impl RawQuerynameKey {
    /// Create a new queryname key.
    ///
    /// The name will be null-terminated for `natural_compare_nul`.
    #[must_use]
    pub fn new(mut name: Vec<u8>, flags: u16) -> Self {
        if name.last() != Some(&0) {
            name.push(0);
        }
        Self { name, flags }
    }

    /// Returns the read name bytes (including the null terminator).
    #[must_use]
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    /// Extract queryname key from raw BAM record bytes.
    /// The name is stored with a null terminator for `natural_compare_nul`.
    #[inline]
    #[must_use]
    fn extract_queryname_key(bam: &[u8]) -> Self {
        let (raw_name, flags) = extract_raw_name_and_flags(bam);
        let mut name = Vec::with_capacity(raw_name.len() + 1);
        name.extend_from_slice(raw_name);
        name.push(0);
        Self { name, flags }
    }
}

impl Ord for RawQuerynameKey {
    #[inline]
    #[allow(unsafe_code)]
    fn cmp(&self, other: &Self) -> Ordering {
        // SAFETY: `name` is always null-terminated (see `extract_queryname_key` and `new`).
        unsafe { natural_compare_nul(self.name.as_ptr(), other.name.as_ptr()) }
            .then_with(|| self.flags.cmp(&other.flags))
    }
}

impl PartialOrd for RawQuerynameKey {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl RawSortKey for RawQuerynameKey {
    const SERIALIZED_SIZE: Option<usize> = None; // Variable-length
    const EMBEDDED_IN_RECORD: bool = true;

    #[inline]
    fn extract(bam: &[u8], _ctx: &SortContext) -> Self {
        Self::extract_queryname_key(bam)
    }

    #[inline]
    fn extract_from_record(bam: &[u8]) -> Self {
        Self::extract_queryname_key(bam)
    }

    #[inline]
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        write_queryname_key(&self.name, self.flags, writer)
    }

    #[inline]
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let (name, flags) = read_queryname_key(reader)?;
        Ok(Self::new(name, flags))
    }
}

// ============================================================================
// Raw Queryname Lexicographic Sort Key
// ============================================================================

/// Variable-size queryname sort key using lexicographic (byte) ordering.
///
/// This is the fast-path queryname sort key. It compares read names as raw bytes
/// (standard lexicographic ordering), which is 8-13x faster than natural ordering.
///
/// Serialization format: `[name_len: u16][name: bytes][flags: u16]`
#[derive(Clone, Eq, PartialEq, Debug, Default)]
pub struct RawQuerynameLexKey {
    /// Read name bytes.
    name: Vec<u8>,
    /// Flags for segment ordering (R1 before R2).
    flags: u16,
}

impl RawQuerynameLexKey {
    /// Create a new lexicographic queryname key.
    #[must_use]
    pub fn new(name: Vec<u8>, flags: u16) -> Self {
        Self { name, flags }
    }

    /// Returns the read name bytes.
    #[must_use]
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    /// Extract queryname key from raw BAM record bytes.
    #[inline]
    #[must_use]
    fn extract_queryname_key(bam: &[u8]) -> Self {
        let (raw_name, flags) = extract_raw_name_and_flags(bam);
        Self { name: raw_name.to_vec(), flags }
    }
}

impl Ord for RawQuerynameLexKey {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        self.name.cmp(&other.name).then_with(|| self.flags.cmp(&other.flags))
    }
}

impl PartialOrd for RawQuerynameLexKey {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl RawSortKey for RawQuerynameLexKey {
    const SERIALIZED_SIZE: Option<usize> = None;
    const EMBEDDED_IN_RECORD: bool = true;

    #[inline]
    fn extract(bam: &[u8], _ctx: &SortContext) -> Self {
        Self::extract_queryname_key(bam)
    }

    #[inline]
    fn extract_from_record(bam: &[u8]) -> Self {
        Self::extract_queryname_key(bam)
    }

    #[inline]
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        write_queryname_key(&self.name, self.flags, writer)
    }

    #[inline]
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let (name, flags) = read_queryname_key(reader)?;
        Ok(Self { name, flags })
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coordinate_key_ordering() {
        let k1 = CoordinateKey { tid: 0, pos: 100, reverse: false };
        let k2 = CoordinateKey { tid: 0, pos: 200, reverse: false };
        let k3 = CoordinateKey { tid: 1, pos: 50, reverse: false };

        assert!(k1 < k2);
        assert!(k2 < k3);
    }

    #[test]
    fn test_coordinate_key_unmapped_last() {
        let mapped = CoordinateKey { tid: 0, pos: 100, reverse: false };
        let unmapped = CoordinateKey::unmapped();

        assert!(mapped < unmapped);
    }

    #[test]
    fn test_primary_alignment_info_roundtrip() {
        let info = PrimaryAlignmentInfo::new(5, 1000, false, 3, 2000, true);
        let value = info.to_tag_value();
        let parsed = PrimaryAlignmentInfo::from_tag_value(&value);

        assert!(parsed.is_some());
        let parsed = parsed.unwrap();
        assert_eq!(parsed.tid1, 5);
        assert_eq!(parsed.pos1, 1000);
        assert!(!parsed.neg1);
        assert_eq!(parsed.tid2, 3);
        assert_eq!(parsed.pos2, 2000);
        assert!(parsed.neg2);
    }

    #[test]
    fn test_primary_alignment_info_from_record() {
        use crate::sam::builder::RecordBuilder;
        use noodles::sam::alignment::record::Flags;

        // Create a supplementary record with pa tag
        let info = PrimaryAlignmentInfo::new(0, 100, false, 0, 200, true);
        let record = RecordBuilder::new()
            .name("test")
            .sequence("ACGT")
            .flags(Flags::SUPPLEMENTARY)
            .tag("pa", info.to_tag_value())
            .build();

        let result = PrimaryAlignmentInfo::from_record(&record);
        assert!(result.is_some());
        let result = result.unwrap();
        assert_eq!(result.tid1, 0);
        assert_eq!(result.pos1, 100);
        assert!(!result.neg1);
        assert_eq!(result.tid2, 0);
        assert_eq!(result.pos2, 200);
        assert!(result.neg2);
    }

    #[test]
    fn test_primary_alignment_info_missing() {
        use crate::sam::builder::RecordBuilder;
        use noodles::sam::alignment::record::Flags;

        let record =
            RecordBuilder::new().name("test").sequence("ACGT").flags(Flags::SUPPLEMENTARY).build();

        let result = PrimaryAlignmentInfo::from_record(&record);
        assert!(result.is_none());
    }

    // ========================================================================
    // from_tag_value tests (fast path coverage)
    // ========================================================================

    #[test]
    fn test_from_tag_value_int32_fast_path() {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        // Create Int32 array (the expected format from zipper)
        let values: Vec<i32> = vec![5, 1000, 0, 3, 2000, 1];
        let value = Value::Array(Array::Int32(values));

        let result = PrimaryAlignmentInfo::from_tag_value(&value);
        assert!(result.is_some());
        let info = result.unwrap();
        assert_eq!(info.tid1, 5);
        assert_eq!(info.pos1, 1000);
        assert!(!info.neg1);
        assert_eq!(info.tid2, 3);
        assert_eq!(info.pos2, 2000);
        assert!(info.neg2);
    }

    #[test]
    fn test_from_tag_value_int32_wrong_length() {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        // Int32 array with wrong number of elements
        let values: Vec<i32> = vec![5, 1000, 0]; // Only 3 elements
        let value = Value::Array(Array::Int32(values));

        let result = PrimaryAlignmentInfo::from_tag_value(&value);
        assert!(result.is_none());
    }

    #[test]
    fn test_from_tag_value_int16_fallback() {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        // Int16 array (rare, but should work via fallback path)
        let values: Vec<i16> = vec![5, 1000, 0, 3, 2000, 1];
        let value = Value::Array(Array::Int16(values));

        let result = PrimaryAlignmentInfo::from_tag_value(&value);
        assert!(result.is_some());
        let info = result.unwrap();
        assert_eq!(info.tid1, 5);
        assert_eq!(info.pos1, 1000);
    }

    #[test]
    fn test_from_tag_value_non_array_returns_none() {
        use noodles::sam::alignment::record_buf::data::field::Value;

        // String value instead of array
        let value = Value::String("not_an_array".into());

        let result = PrimaryAlignmentInfo::from_tag_value(&value);
        assert!(result.is_none());
    }

    #[test]
    fn test_from_tag_value_float_array_returns_none() {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        // Float array is not supported
        let values: Vec<f32> = vec![5.0, 1000.0, 0.0, 3.0, 2000.0, 1.0];
        let value = Value::Array(Array::Float(values));

        let result = PrimaryAlignmentInfo::from_tag_value(&value);
        assert!(result.is_none());
    }

    // ========================================================================
    // PrimaryAlignmentInfo::new edge cases
    // ========================================================================

    #[test]
    fn test_primary_alignment_info_new_stores_values_unchanged() {
        // Verify new() stores values exactly as provided (no normalization)
        let info = PrimaryAlignmentInfo::new(5, 1000, true, 3, 500, false);

        // Values should be stored exactly as provided
        assert_eq!(info.tid1, 5);
        assert_eq!(info.pos1, 1000);
        assert!(info.neg1);
        assert_eq!(info.tid2, 3);
        assert_eq!(info.pos2, 500);
        assert!(!info.neg2);
    }

    #[test]
    fn test_primary_alignment_info_new_with_negative_positions() {
        // Edge case: negative positions (can happen with soft clips before position 0)
        let info = PrimaryAlignmentInfo::new(0, -5, false, 0, -10, true);

        assert_eq!(info.tid1, 0);
        assert_eq!(info.pos1, -5);
        assert!(!info.neg1);
        assert_eq!(info.tid2, 0);
        assert_eq!(info.pos2, -10);
        assert!(info.neg2);
    }

    #[test]
    fn test_primary_alignment_info_new_with_max_values() {
        // Edge case: maximum i32 values
        let info = PrimaryAlignmentInfo::new(i32::MAX, i32::MAX, true, i32::MAX, i32::MAX, true);

        assert_eq!(info.tid1, i32::MAX);
        assert_eq!(info.pos1, i32::MAX);
        assert!(info.neg1);
        assert_eq!(info.tid2, i32::MAX);
        assert_eq!(info.pos2, i32::MAX);
        assert!(info.neg2);
    }

    #[test]
    fn test_primary_alignment_info_new_with_zero_values() {
        // Edge case: all zeros (unmapped or start of reference)
        let info = PrimaryAlignmentInfo::new(0, 0, false, 0, 0, false);

        assert_eq!(info.tid1, 0);
        assert_eq!(info.pos1, 0);
        assert!(!info.neg1);
        assert_eq!(info.tid2, 0);
        assert_eq!(info.pos2, 0);
        assert!(!info.neg2);
    }

    #[test]
    fn test_primary_alignment_info_roundtrip_with_negative_positions() {
        // Verify negative positions survive roundtrip
        let info = PrimaryAlignmentInfo::new(0, -10, false, 0, -5, true);
        let value = info.to_tag_value();
        let parsed = PrimaryAlignmentInfo::from_tag_value(&value).unwrap();

        assert_eq!(parsed.tid1, 0);
        assert_eq!(parsed.pos1, -10);
        assert!(!parsed.neg1);
        assert_eq!(parsed.tid2, 0);
        assert_eq!(parsed.pos2, -5);
        assert!(parsed.neg2);
    }

    #[test]
    fn test_primary_alignment_info_position_order_is_caller_responsibility() {
        // Verify that new() does NOT enforce pos1 < pos2 ordering
        // (ordering is done by the caller in zipper.rs)
        let info = PrimaryAlignmentInfo::new(0, 2000, false, 0, 1000, false);

        // Values are stored as provided, even if "out of order"
        assert_eq!(info.pos1, 2000); // pos1 > pos2 is allowed
        assert_eq!(info.pos2, 1000);
    }

    // ========================================================================
    // queryname_flag_order tests
    // ========================================================================

    /// Test exact transformation values match samtools formula.
    /// Formula: ((flags & 0xc0) << 8) | ((flags & 0x100) << 3) | ((flags & 0x800) >> 3)
    #[test]
    fn test_queryname_flag_order_exact_transformation_values() {
        // Flag bits:
        // 0x40  = READ1 (bit 6)
        // 0x80  = READ2 (bit 7)
        // 0x100 = SECONDARY (bit 8)
        // 0x800 = SUPPLEMENTARY (bit 11)

        // NONE (unpaired primary)
        assert_eq!(queryname_flag_order(0x0000), 0x0000);

        // R1 PRIMARY: 0x40 << 8 = 0x4000
        assert_eq!(queryname_flag_order(0x0040), 0x4000);

        // R2 PRIMARY: 0x80 << 8 = 0x8000
        assert_eq!(queryname_flag_order(0x0080), 0x8000);

        // SECONDARY only: 0x100 << 3 = 0x800
        assert_eq!(queryname_flag_order(0x0100), 0x0800);

        // SUPPLEMENTARY only: 0x800 >> 3 = 0x100
        assert_eq!(queryname_flag_order(0x0800), 0x0100);

        // R1 SUPPLEMENTARY: (0x40 << 8) | (0x800 >> 3) = 0x4000 | 0x100 = 0x4100
        assert_eq!(queryname_flag_order(0x0840), 0x4100);

        // R1 SECONDARY: (0x40 << 8) | (0x100 << 3) = 0x4000 | 0x800 = 0x4800
        assert_eq!(queryname_flag_order(0x0140), 0x4800);

        // R2 SUPPLEMENTARY: (0x80 << 8) | (0x800 >> 3) = 0x8000 | 0x100 = 0x8100
        assert_eq!(queryname_flag_order(0x0880), 0x8100);

        // R2 SECONDARY: (0x80 << 8) | (0x100 << 3) = 0x8000 | 0x800 = 0x8800
        assert_eq!(queryname_flag_order(0x0180), 0x8800);
    }

    /// Test the complete sort order of all 12 categories.
    /// Order: NONE < R1 PRIMARY < R1 SUPP < R1 SEC < R2 PRIMARY < R2 SUPP < R2 SEC
    ///        (and combined R1+R2 categories)
    #[test]
    fn test_queryname_flag_order_complete_sort_order() {
        // All 12 possible combinations sorted in expected order
        let flags_in_order = [
            // NONE (unpaired)
            0x0000u16, // NONE PRIMARY
            0x0800,    // NONE SUPPLEMENTARY
            0x0100,    // NONE SECONDARY
            // R1
            0x0040, // R1 PRIMARY
            0x0840, // R1 SUPPLEMENTARY
            0x0140, // R1 SECONDARY
            // R2
            0x0080, // R2 PRIMARY
            0x0880, // R2 SUPPLEMENTARY
            0x0180, // R2 SECONDARY
            // R1+R2 (unusual but possible in edge cases)
            0x00c0, // R1+R2 PRIMARY
            0x08c0, // R1+R2 SUPPLEMENTARY
            0x01c0, // R1+R2 SECONDARY
        ];

        // Transform all flags
        let transformed: Vec<u16> =
            flags_in_order.iter().map(|&f| queryname_flag_order(f)).collect();

        // Verify sorted order
        for i in 0..transformed.len() - 1 {
            assert!(
                transformed[i] < transformed[i + 1],
                "Expected flags 0x{:04x} (transformed 0x{:04x}) < 0x{:04x} (transformed 0x{:04x})",
                flags_in_order[i],
                transformed[i],
                flags_in_order[i + 1],
                transformed[i + 1]
            );
        }
    }

    /// Test R1 always sorts before R2.
    #[test]
    fn test_queryname_flag_order_r1_before_r2() {
        // R1 variants
        let r1_primary = queryname_flag_order(0x0040);
        let r1_supp = queryname_flag_order(0x0840);
        let r1_sec = queryname_flag_order(0x0140);

        // R2 variants
        let r2_primary = queryname_flag_order(0x0080);
        let r2_supp = queryname_flag_order(0x0880);
        let r2_sec = queryname_flag_order(0x0180);

        // All R1 variants should be less than all R2 variants
        assert!(r1_primary < r2_primary);
        assert!(r1_primary < r2_supp);
        assert!(r1_primary < r2_sec);
        assert!(r1_supp < r2_primary);
        assert!(r1_supp < r2_supp);
        assert!(r1_supp < r2_sec);
        assert!(r1_sec < r2_primary);
        assert!(r1_sec < r2_supp);
        assert!(r1_sec < r2_sec);
    }

    /// Test PRIMARY < SUPPLEMENTARY < SECONDARY within each read type.
    #[test]
    fn test_queryname_flag_order_primary_before_supplementary_before_secondary() {
        // Test for NONE (unpaired)
        let none_pri = queryname_flag_order(0x0000);
        let none_supp = queryname_flag_order(0x0800);
        let none_sec = queryname_flag_order(0x0100);
        assert!(none_pri < none_supp, "NONE: PRIMARY < SUPP");
        assert!(none_supp < none_sec, "NONE: SUPP < SEC");

        // Test for R1
        let r1_pri = queryname_flag_order(0x0040);
        let r1_supp = queryname_flag_order(0x0840);
        let r1_sec = queryname_flag_order(0x0140);
        assert!(r1_pri < r1_supp, "R1: PRIMARY < SUPP");
        assert!(r1_supp < r1_sec, "R1: SUPP < SEC");

        // Test for R2
        let r2_pri = queryname_flag_order(0x0080);
        let r2_supp = queryname_flag_order(0x0880);
        let r2_sec = queryname_flag_order(0x0180);
        assert!(r2_pri < r2_supp, "R2: PRIMARY < SUPP");
        assert!(r2_supp < r2_sec, "R2: SUPP < SEC");
    }

    /// Test NONE (unpaired) sorts before R1 sorts before R2.
    #[test]
    fn test_queryname_flag_order_none_before_r1_before_r2() {
        let none = queryname_flag_order(0x0000);
        let r1 = queryname_flag_order(0x0040);
        let r2 = queryname_flag_order(0x0080);

        assert!(none < r1, "NONE < R1");
        assert!(r1 < r2, "R1 < R2");
    }

    /// Test that irrelevant flags don't affect ordering.
    /// Only 0x40, 0x80, 0x100, 0x800 should matter.
    #[test]
    fn test_queryname_flag_order_irrelevant_flags_do_not_affect_order() {
        // Irrelevant flags: PAIRED(0x1), PROPER_PAIR(0x2), UNMAPPED(0x4), MATE_UNMAPPED(0x8),
        // REVERSE(0x10), MATE_REVERSE(0x20), DUPLICATE(0x400), FAIL_QC(0x200)
        let irrelevant_flags: u16 = 0x01 | 0x02 | 0x04 | 0x08 | 0x10 | 0x20 | 0x200 | 0x400;

        // Base case: R1 PRIMARY
        let r1_base = queryname_flag_order(0x0040);
        let r1_with_irrelevant = queryname_flag_order(0x0040 | irrelevant_flags);
        assert_eq!(r1_base, r1_with_irrelevant, "Irrelevant flags should not change result");

        // R2 SUPPLEMENTARY
        let r2_supp_base = queryname_flag_order(0x0880);
        let r2_supp_with = queryname_flag_order(0x0880 | irrelevant_flags);
        assert_eq!(r2_supp_base, r2_supp_with);

        // NONE SECONDARY
        let none_sec_base = queryname_flag_order(0x0100);
        let none_sec_with = queryname_flag_order(0x0100 | irrelevant_flags);
        assert_eq!(none_sec_base, none_sec_with);
    }

    /// Test with real-world flag combinations commonly seen in BAM files.
    #[test]
    fn test_queryname_flag_order_real_world_flags() {
        // Common real-world flag values (from paired-end sequencing)

        // R1 forward, properly paired: 0x63 = PAIRED|PROPER|MATE_REV|R1 = 0x1|0x2|0x20|0x40
        let r1_fwd = queryname_flag_order(0x0063);
        assert_eq!(r1_fwd, 0x4000, "R1 fwd should extract to R1 PRIMARY");

        // R2 reverse, properly paired: 0x93 = PAIRED|PROPER|REV|R2 = 0x1|0x2|0x10|0x80
        let r2_rev = queryname_flag_order(0x0093);
        assert_eq!(r2_rev, 0x8000, "R2 rev should extract to R2 PRIMARY");

        // R1 supplementary: 0x841 = PAIRED|R1|SUPP = 0x1|0x40|0x800
        let r1_supp = queryname_flag_order(0x0841);
        assert_eq!(r1_supp, 0x4100, "R1 supp should extract to R1 SUPP");

        // R2 secondary: 0x181 = PAIRED|R2|SEC = 0x1|0x80|0x100
        let r2_sec = queryname_flag_order(0x0181);
        assert_eq!(r2_sec, 0x8800, "R2 sec should extract to R2 SEC");

        // Verify ordering
        assert!(r1_fwd < r1_supp);
        assert!(r1_supp < r2_rev);
        assert!(r2_rev < r2_sec);
    }

    /// Test from actual test data showing the bug fix.
    /// samtools order: 113 → 2161 → 177 (R1 primary → R1 supplementary → R2 primary)
    /// fgumi (before fix): 113 → 177 → 2161 (wrong - puts R2 primary before R1 supp)
    #[test]
    fn test_queryname_flag_order_from_test_data() {
        // Flags from actual test data:
        // 113 = 0x71 = PAIRED|PROPER|REV|R1 (R1 PRIMARY on reverse strand)
        // 177 = 0xB1 = PAIRED|PROPER|REV|R2 (R2 PRIMARY on reverse strand)
        // 2161 = 0x871 = PAIRED|PROPER|REV|R1|SUPP (R1 SUPPLEMENTARY)

        let f113 = queryname_flag_order(113); // R1 PRIMARY
        let f177 = queryname_flag_order(177); // R2 PRIMARY
        let f2161 = queryname_flag_order(2161); // R1 SUPPLEMENTARY

        // Correct order: R1 PRIMARY < R1 SUPP < R2 PRIMARY
        assert!(f113 < f2161, "R1 PRIMARY (113) should be < R1 SUPP (2161): {f113} vs {f2161}");
        assert!(f2161 < f177, "R1 SUPP (2161) should be < R2 PRIMARY (177): {f2161} vs {f177}");
    }

    /// Test edge case: both R1 and R2 flags set (unusual but possible).
    #[test]
    fn test_queryname_flag_order_edge_case_both_r1_r2() {
        // Both R1 and R2 set: 0xc0
        let both = queryname_flag_order(0x00c0);

        // Should combine: (0xc0 << 8) = 0xc000
        assert_eq!(both, 0xc000);

        // Should sort after both R1-only and R2-only
        let r1_only = queryname_flag_order(0x0040);
        let r2_only = queryname_flag_order(0x0080);
        assert!(r1_only < both);
        assert!(r2_only < both);
    }

    /// Test that the function is const-evaluable.
    #[test]
    fn test_queryname_flag_order_is_const() {
        // This compiles only if queryname_flag_order is const
        const R1_PRIMARY: u16 = queryname_flag_order(0x0040);
        const R2_PRIMARY: u16 = queryname_flag_order(0x0080);

        assert_eq!(R1_PRIMARY, 0x4000);
        assert_eq!(R2_PRIMARY, 0x8000);
    }

    // ========================================================================
    // Comprehensive natural_compare tests
    // ========================================================================

    use rstest::rstest;

    // --- Basic string comparison (no digits) ---

    #[rstest]
    #[case(b"abc", b"abd", Ordering::Less)]
    #[case(b"abd", b"abc", Ordering::Greater)]
    #[case(b"abc", b"abc", Ordering::Equal)]
    #[case(b"ABC", b"abc", Ordering::Less)] // uppercase < lowercase in ASCII
    #[case(b"z", b"a", Ordering::Greater)]
    #[case(b"aaa", b"aab", Ordering::Less)]
    #[case(b"ZZZ", b"aaa", Ordering::Less)] // 'Z' (90) < 'a' (97)
    fn test_natural_compare_alpha_only(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: Ordering,
    ) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // --- Empty strings and single characters ---

    #[rstest]
    #[case(b"", b"", Ordering::Equal)]
    #[case(b"", b"x", Ordering::Less)]
    #[case(b"x", b"", Ordering::Greater)]
    #[case(b"", b"0", Ordering::Less)]
    #[case(b"0", b"", Ordering::Greater)]
    #[case(b"a", b"a", Ordering::Equal)]
    #[case(b"0", b"0", Ordering::Equal)]
    fn test_natural_compare_empty_and_single(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: Ordering,
    ) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // --- Prefix relationships ---

    #[rstest]
    #[case(b"abc", b"abcd", Ordering::Less)]
    #[case(b"abcd", b"abc", Ordering::Greater)]
    #[case(b"read", b"read1", Ordering::Less)]
    #[case(b"read1", b"read", Ordering::Greater)]
    #[case(b"a1", b"a1b", Ordering::Less)]
    #[case(b"a1b", b"a1", Ordering::Greater)]
    fn test_natural_compare_prefix_relationships(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: Ordering,
    ) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // --- Single-digit numeric comparison ---

    #[rstest]
    #[case(b"0", b"1", Ordering::Less)]
    #[case(b"1", b"2", Ordering::Less)]
    #[case(b"9", b"0", Ordering::Greater)]
    #[case(b"5", b"5", Ordering::Equal)]
    fn test_natural_compare_single_digit(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: Ordering,
    ) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // --- Multi-digit numeric comparison ---

    #[rstest]
    #[case(b"2", b"10", Ordering::Less)] // numeric: 2 < 10
    #[case(b"10", b"2", Ordering::Greater)]
    #[case(b"9", b"10", Ordering::Less)]
    #[case(b"10", b"9", Ordering::Greater)]
    #[case(b"10", b"10", Ordering::Equal)]
    #[case(b"99", b"100", Ordering::Less)]
    #[case(b"100", b"99", Ordering::Greater)]
    #[case(b"123", b"456", Ordering::Less)]
    #[case(b"456", b"123", Ordering::Greater)]
    #[case(b"999", b"1000", Ordering::Less)]
    #[case(b"1000", b"999", Ordering::Greater)]
    fn test_natural_compare_multidigit(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: Ordering,
    ) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // --- Leading zeros ---
    // Leading zeros are stripped during numeric comparison, so "01" and "1" have the
    // same numeric value. However, the implementation uses total string length as the
    // final tiebreaker (alen.cmp(&blen)), so strings with more leading zeros compare
    // as greater when the rest of the string is otherwise equal. This differs from
    // samtools' strnum_cmp which returns Equal for "01" vs "1".
    //
    // Within a numeric run, leading zeros are skipped, so numeric magnitude is compared
    // correctly regardless of leading zeros.

    #[rstest]
    // Same numeric value but different string lengths due to leading zeros:
    // the longer string sorts later due to the alen.cmp(&blen) tiebreaker.
    #[case(b"01", b"1", Ordering::Greater)]
    #[case(b"1", b"01", Ordering::Less)]
    #[case(b"001", b"1", Ordering::Greater)]
    #[case(b"001", b"01", Ordering::Greater)]
    #[case(b"010", b"10", Ordering::Greater)]
    #[case(b"0010", b"10", Ordering::Greater)]
    #[case(b"00", b"0", Ordering::Greater)]
    #[case(b"000", b"0", Ordering::Greater)]
    // Different numeric values with leading zeros: magnitude still wins.
    #[case(b"02", b"1", Ordering::Greater)] // 2 > 1 (and longer)
    #[case(b"02", b"10", Ordering::Less)] // 2 < 10
    #[case(b"009", b"10", Ordering::Less)] // 9 < 10
    #[case(b"0100", b"99", Ordering::Greater)] // 100 > 99
    fn test_natural_compare_leading_zeros(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: Ordering,
    ) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // --- Digits sort before non-digits ---
    // The fgumi implementation has: when one char is a digit and the other is not,
    // the digit sorts first (Less).

    #[rstest]
    #[case(b"0", b"a", Ordering::Less)] // digit < non-digit
    #[case(b"9", b"a", Ordering::Less)]
    #[case(b"a", b"0", Ordering::Greater)]
    #[case(b"a", b"9", Ordering::Greater)]
    #[case(b"1", b"A", Ordering::Less)]
    #[case(b"1", b"Z", Ordering::Less)]
    #[case(b"1", b"z", Ordering::Less)]
    fn test_natural_compare_digits_before_nondigits(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: Ordering,
    ) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // --- Mixed alphanumeric sequences ---

    #[rstest]
    #[case(b"a1b2", b"a1b10", Ordering::Less)]
    #[case(b"a1b10", b"a1b2", Ordering::Greater)]
    #[case(b"x10y20", b"x10y3", Ordering::Greater)]
    #[case(b"x10y3", b"x10y20", Ordering::Less)]
    #[case(b"abc123def", b"abc123deg", Ordering::Less)]
    #[case(b"abc123def", b"abc124def", Ordering::Less)]
    #[case(b"abc124def", b"abc123def", Ordering::Greater)]
    #[case(b"file1part2", b"file1part10", Ordering::Less)]
    #[case(b"file2part1", b"file10part1", Ordering::Less)]
    #[case(b"v1.2.3", b"v1.2.10", Ordering::Less)] // version-like strings
    #[case(b"v1.2.10", b"v1.10.2", Ordering::Less)]
    fn test_natural_compare_mixed_alphanumeric(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: Ordering,
    ) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // --- Samtools sort test data: expected queryname ordering ---
    // From samtools test/sort/name.sort.expected.sam and name2.sort.expected.sam

    #[test]
    fn test_natural_compare_samtools_name_sort_order() {
        // From name.sort.expected.sam: r000, r001, r002, r003, r004, u1, x1..x6
        let names: Vec<&[u8]> = vec![
            b"r000", b"r001", b"r002", b"r003", b"r004", b"u1", b"x1", b"x2", b"x3", b"x4", b"x5",
            b"x6",
        ];
        for i in 0..names.len() - 1 {
            assert_eq!(
                natural_compare(names[i], names[i + 1]),
                Ordering::Less,
                "{} should sort before {}",
                std::str::from_utf8(names[i]).unwrap(),
                std::str::from_utf8(names[i + 1]).unwrap(),
            );
        }
    }

    #[test]
    fn test_natural_compare_samtools_name2_sort_order() {
        // From name2.sort.expected.sam: consecutive pairs should be in natural order.
        // Note: x7 < x8 < x9 < x10 < x11 < x12 in natural order.
        let names: Vec<&[u8]> =
            vec![b"r005", b"r006", b"r007", b"x7", b"x8", b"x9", b"x10", b"x11", b"x12"];
        for i in 0..names.len() - 1 {
            assert_eq!(
                natural_compare(names[i], names[i + 1]),
                Ordering::Less,
                "{} should sort before {}",
                std::str::from_utf8(names[i]).unwrap(),
                std::str::from_utf8(names[i + 1]).unwrap(),
            );
        }
    }

    // --- Illumina read names ---

    #[rstest]
    // Same flowcell, different tiles
    #[case(
        b"A00132:53:HFHJKDSXX:1:1101:12345:67890",
        b"A00132:53:HFHJKDSXX:1:1102:12345:67890",
        Ordering::Less
    )]
    // Same tile, different X coordinates
    #[case(
        b"A00132:53:HFHJKDSXX:1:1101:12345:67890",
        b"A00132:53:HFHJKDSXX:1:1101:12346:67890",
        Ordering::Less
    )]
    // Different lanes
    #[case(
        b"A00132:53:HFHJKDSXX:1:1101:12345:67890",
        b"A00132:53:HFHJKDSXX:2:1101:12345:67890",
        Ordering::Less
    )]
    // Different runs
    #[case(
        b"A00132:53:HFHJKDSXX:1:1101:12345:67890",
        b"A00132:54:HFHJKDSXX:1:1101:12345:67890",
        Ordering::Less
    )]
    // Same everything
    #[case(
        b"A00132:53:HFHJKDSXX:1:1101:12345:67890",
        b"A00132:53:HFHJKDSXX:1:1101:12345:67890",
        Ordering::Equal
    )]
    // Different instruments
    #[case(
        b"A00132:53:HFHJKDSXX:1:1101:12345:67890",
        b"B00132:53:HFHJKDSXX:1:1101:12345:67890",
        Ordering::Less
    )]
    // Numeric tile comparison (1101 vs 2201)
    #[case(b"INST:1:FLOW:1:1101:1:1", b"INST:1:FLOW:1:2201:1:1", Ordering::Less)]
    fn test_natural_compare_illumina_names(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: Ordering,
    ) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // --- SRR read names ---

    #[rstest]
    #[case(b"SRR099966.1", b"SRR099966.2", Ordering::Less)]
    #[case(b"SRR099966.9", b"SRR099966.10", Ordering::Less)]
    #[case(b"SRR099966.99", b"SRR099966.100", Ordering::Less)]
    #[case(b"SRR099966.12345", b"SRR099966.12346", Ordering::Less)]
    #[case(b"SRR099966.12345", b"SRR099966.12345", Ordering::Equal)]
    #[case(b"SRR099966.12345", b"SRR099967.1", Ordering::Less)] // 6 < 7 in accession
    #[case(b"SRR1.1", b"SRR2.1", Ordering::Less)]
    #[case(b"SRR9.1", b"SRR10.1", Ordering::Less)]
    fn test_natural_compare_srr_names(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: Ordering,
    ) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // --- Synthetic read names ---

    #[rstest]
    #[case(b"mol000001_read0001", b"mol000001_read0002", Ordering::Less)]
    #[case(b"mol000001_read0009", b"mol000001_read0010", Ordering::Less)]
    #[case(b"mol000001_read0001", b"mol000002_read0001", Ordering::Less)]
    #[case(b"mol1_read1", b"mol000001_read0001", Ordering::Less)] // same values but shorter string
    #[case(b"mol1_read1", b"mol1_read2", Ordering::Less)]
    #[case(b"mol2_read1", b"mol10_read1", Ordering::Less)]
    fn test_natural_compare_synthetic_names(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: Ordering,
    ) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // --- Very long numeric runs (test overflow robustness) ---
    // The algorithm compares digit-by-digit without parsing to integers,
    // so it should handle arbitrarily long numbers correctly.

    #[test]
    fn test_natural_compare_long_numeric_runs() {
        // Numbers larger than u64::MAX (which is 18446744073709551615, 20 digits)
        let big_a = b"x99999999999999999999999999999999999999"; // 38 digits of 9
        let big_b = b"x100000000000000000000000000000000000000"; // 1 followed by 38 zeros
        assert_eq!(natural_compare(big_a, big_b), Ordering::Less);

        // Equal very long numbers
        let long_num = b"prefix99999999999999999999999999suffix";
        assert_eq!(natural_compare(long_num, long_num), Ordering::Equal);

        // Different very long numbers of same length
        let a = b"x99999999999999999999999999999999999998";
        let b = b"x99999999999999999999999999999999999999";
        assert_eq!(natural_compare(a, b), Ordering::Less);
    }

    #[test]
    fn test_natural_compare_u64_overflow_boundary() {
        // u64::MAX = 18446744073709551615
        let u64_max = b"read18446744073709551615";
        let u64_max_plus_1 = b"read18446744073709551616";
        assert_eq!(natural_compare(u64_max, u64_max_plus_1), Ordering::Less);

        // Two numbers that would overflow u64 if parsed
        let huge_a = b"read99999999999999999999"; // 20 digits, all 9s
        let huge_b = b"read100000000000000000000"; // 21 digits
        assert_eq!(natural_compare(huge_a, huge_b), Ordering::Less);
    }

    // --- Property-based comparator invariant tests ---

    proptest::proptest! {
        #[test]
        fn test_natural_compare_symmetry(a: Vec<u8>, b: Vec<u8>) {
            let ab = natural_compare(&a, &b);
            let ba = natural_compare(&b, &a);
            proptest::prop_assert_eq!(ab, ba.reverse());
        }

        #[test]
        fn test_natural_compare_reflexive(s: Vec<u8>) {
            proptest::prop_assert_eq!(natural_compare(&s, &s), Ordering::Equal);
        }

        #[test]
        fn test_natural_compare_transitivity(
            mut values in proptest::collection::vec(proptest::collection::vec(0..=255u8, 0..20), 3..10),
        ) {
            // Sort values with natural_compare, then verify all pairs are ordered.
            values.sort_by(|a, b| natural_compare(a, b));
            for i in 0..values.len() {
                for j in (i + 1)..values.len() {
                    let cmp = natural_compare(&values[i], &values[j]);
                    proptest::prop_assert!(
                        cmp != Ordering::Greater,
                        "Transitivity violated: sorted[{}] ({:?}) > sorted[{}] ({:?})",
                        i, values[i], j, values[j],
                    );
                }
            }
        }
    }

    // --- Samtools strnum_cmp edge cases ---
    // Ported from samtools bam_sort.c strnum_cmp behavior analysis.

    #[test]
    fn test_natural_compare_leading_zeros_tiebreak_by_length() {
        // Unlike samtools strnum_cmp which returns Equal for "01" vs "1",
        // this implementation uses total string length as tiebreaker.
        // Numerically equal values with different leading zeros compare
        // by string length (longer = greater).
        assert_eq!(natural_compare(b"00", b"0"), Ordering::Greater);
        assert_eq!(natural_compare(b"000", b"00"), Ordering::Greater);
        assert_eq!(natural_compare(b"0", b"00"), Ordering::Less);
        assert_eq!(natural_compare(b"01", b"1"), Ordering::Greater);
        assert_eq!(natural_compare(b"00100", b"100"), Ordering::Greater);
    }

    #[test]
    fn test_natural_compare_strnum_cmp_compat_numeric_vs_alpha_boundary() {
        // When one string has a digit and the other a letter at the same position,
        // in samtools strnum_cmp it compares raw bytes: '0' (48) < 'A' (65) < 'a' (97).
        // In fgumi, digits explicitly sort before non-digits. Same practical result
        // for typical chars, but it's a deliberate design choice.
        assert_eq!(natural_compare(b"a0", b"aa"), Ordering::Less);
        assert_eq!(natural_compare(b"a9", b"aa"), Ordering::Less);
        assert_eq!(natural_compare(b"aa", b"a0"), Ordering::Greater);
    }

    #[test]
    fn test_natural_compare_strnum_cmp_compat_numeric_then_string() {
        // After a numeric run, comparison resumes at the character level.
        // "a1b" vs "a1c": after matching "a" and "1", compare "b" vs "c"
        assert_eq!(natural_compare(b"a1b", b"a1c"), Ordering::Less);
        assert_eq!(natural_compare(b"a10b", b"a10c"), Ordering::Less);
        assert_eq!(natural_compare(b"a10b", b"a2c"), Ordering::Greater); // 10 > 2
    }

    #[test]
    fn test_natural_compare_all_zeros_tiebreak_by_length() {
        // "0" vs "00" vs "000": numerically equal but differ in string length.
        // Longer string sorts later due to alen.cmp(&blen) tiebreaker.
        assert_eq!(natural_compare(b"0", b"00"), Ordering::Less);
        assert_eq!(natural_compare(b"00", b"000"), Ordering::Less);
        assert_eq!(natural_compare(b"0", b"000"), Ordering::Less);

        // Embedded zeros in otherwise equal strings: "a0b" (len=3) vs "a00b" (len=4).
        // The numeric run "0" vs "00" are numerically equal, so comparison continues.
        // Both then hit "b" which matches, then the longer string wins.
        assert_eq!(natural_compare(b"a0b", b"a00b"), Ordering::Less);
        assert_eq!(natural_compare(b"a0b", b"a000b"), Ordering::Less);
    }

    // --- Digit at end of string vs more characters ---

    #[rstest]
    #[case(b"read1", b"read1a", Ordering::Less)] // "read1" is prefix of "read1a"
    #[case(b"read1a", b"read1", Ordering::Greater)]
    #[case(b"read10", b"read10a", Ordering::Less)]
    #[case(b"read10a", b"read10", Ordering::Greater)]
    fn test_natural_compare_numeric_then_suffix(
        #[case] a: &[u8],
        #[case] b: &[u8],
        #[case] expected: Ordering,
    ) {
        assert_eq!(natural_compare(a, b), expected);
    }

    // --- Multiple numeric segments with varying widths ---

    #[test]
    fn test_natural_compare_multiple_numeric_segments() {
        // Simulating version-like or multi-field names
        assert_eq!(natural_compare(b"r1.1.1", b"r1.1.2"), Ordering::Less);
        assert_eq!(natural_compare(b"r1.1.9", b"r1.1.10"), Ordering::Less);
        assert_eq!(natural_compare(b"r1.9.1", b"r1.10.1"), Ordering::Less);
        assert_eq!(natural_compare(b"r9.1.1", b"r10.1.1"), Ordering::Less);
    }

    // --- Consecutive digits without separators ---

    #[test]
    fn test_natural_compare_adjacent_numeric_runs() {
        // When digits are separated by non-digit characters, each run is independent.
        // "a1b2" vs "a1b10": first runs "1"="1", then "2" < "10"
        assert_eq!(natural_compare(b"a1b2", b"a1b10"), Ordering::Less);

        // But a continuous digit run is one number:
        // "a12" vs "a9": "12" > "9"
        assert_eq!(natural_compare(b"a12", b"a9"), Ordering::Greater);
    }

    // --- Non-ASCII bytes ---

    #[test]
    fn test_natural_compare_high_bytes() {
        // Bytes > 127 should still compare correctly by byte value
        assert_eq!(natural_compare(&[0xFF], &[0xFE]), Ordering::Greater);
        assert_eq!(natural_compare(&[0x80], &[0x7F]), Ordering::Greater);
        // Digits still sort before high bytes
        assert_eq!(natural_compare(b"0", &[0x80]), Ordering::Less);
    }

    // --- Stress: sorting a realistic set of read names ---

    #[test]
    fn test_natural_compare_sort_illumina_batch() {
        let names: Vec<&[u8]> = vec![
            b"INST:100:FLOW:2:2201:9999:1",
            b"INST:1:FLOW:1:1101:1:1",
            b"INST:1:FLOW:1:1101:1:2",
            b"INST:1:FLOW:1:1101:2:1",
            b"INST:1:FLOW:1:1102:1:1",
            b"INST:1:FLOW:2:1101:1:1",
            b"INST:2:FLOW:1:1101:1:1",
            b"INST:10:FLOW:1:1101:1:1",
            b"INST:100:FLOW:1:1101:1:1",
        ];

        // Clone and sort
        let mut sorted = names.clone();
        sorted.sort_by(|a, b| natural_compare(a, b));

        // Expected natural order
        let expected: Vec<&[u8]> = vec![
            b"INST:1:FLOW:1:1101:1:1",
            b"INST:1:FLOW:1:1101:1:2",
            b"INST:1:FLOW:1:1101:2:1",
            b"INST:1:FLOW:1:1102:1:1",
            b"INST:1:FLOW:2:1101:1:1",
            b"INST:2:FLOW:1:1101:1:1",
            b"INST:10:FLOW:1:1101:1:1",
            b"INST:100:FLOW:1:1101:1:1",
            b"INST:100:FLOW:2:2201:9999:1",
        ];

        assert_eq!(sorted, expected);
    }

    #[test]
    fn test_natural_compare_sort_mixed_name_formats() {
        let mut names: Vec<&[u8]> = vec![
            b"SRR099966.100",
            b"SRR099966.9",
            b"SRR099966.10",
            b"SRR099966.1",
            b"SRR099966.99",
            b"SRR099966.2",
        ];

        names.sort_by(|a, b| natural_compare(a, b));

        let expected: Vec<&[u8]> = vec![
            b"SRR099966.1",
            b"SRR099966.2",
            b"SRR099966.9",
            b"SRR099966.10",
            b"SRR099966.99",
            b"SRR099966.100",
        ];

        assert_eq!(names, expected);
    }

    // ========================================================================
    // QuerynameComparator tests
    // ========================================================================

    #[test]
    fn test_queryname_comparator_default_is_lexicographic() {
        assert_eq!(QuerynameComparator::default(), QuerynameComparator::Lexicographic);
    }

    #[test]
    fn test_queryname_comparator_display() {
        assert_eq!(QuerynameComparator::Lexicographic.to_string(), "lexicographic");
        assert_eq!(QuerynameComparator::Natural.to_string(), "natural");
    }

    // ========================================================================
    // SortOrder with QuerynameComparator tests
    // ========================================================================

    #[test]
    fn test_sort_order_is_queryname() {
        assert!(!SortOrder::Coordinate.is_queryname());
        assert!(SortOrder::Queryname(QuerynameComparator::Lexicographic).is_queryname());
        assert!(SortOrder::Queryname(QuerynameComparator::Natural).is_queryname());
        assert!(!SortOrder::TemplateCoordinate.is_queryname());
    }

    #[test]
    fn test_sort_order_queryname_comparator() {
        assert_eq!(SortOrder::Coordinate.queryname_comparator(), None);
        assert_eq!(
            SortOrder::Queryname(QuerynameComparator::Lexicographic).queryname_comparator(),
            Some(QuerynameComparator::Lexicographic)
        );
        assert_eq!(
            SortOrder::Queryname(QuerynameComparator::Natural).queryname_comparator(),
            Some(QuerynameComparator::Natural)
        );
        assert_eq!(SortOrder::TemplateCoordinate.queryname_comparator(), None);
    }

    #[test]
    fn test_sort_order_header_so_tag() {
        assert_eq!(SortOrder::Coordinate.header_so_tag(), "coordinate");
        assert_eq!(
            SortOrder::Queryname(QuerynameComparator::Lexicographic).header_so_tag(),
            "queryname"
        );
        assert_eq!(SortOrder::Queryname(QuerynameComparator::Natural).header_so_tag(), "queryname");
        assert_eq!(SortOrder::TemplateCoordinate.header_so_tag(), "unsorted");
    }

    #[test]
    fn test_sort_order_header_ss_tag() {
        assert_eq!(SortOrder::Coordinate.header_ss_tag(), None);
        assert_eq!(
            SortOrder::Queryname(QuerynameComparator::Lexicographic).header_ss_tag(),
            Some("lexicographic")
        );
        assert_eq!(
            SortOrder::Queryname(QuerynameComparator::Natural).header_ss_tag(),
            Some("natural")
        );
        assert_eq!(SortOrder::TemplateCoordinate.header_ss_tag(), Some("template-coordinate"));
    }

    #[test]
    fn test_sort_order_header_go_tag() {
        assert_eq!(SortOrder::Coordinate.header_go_tag(), None);
        assert_eq!(SortOrder::Queryname(QuerynameComparator::Lexicographic).header_go_tag(), None);
        assert_eq!(SortOrder::TemplateCoordinate.header_go_tag(), Some("query"));
    }

    // ========================================================================
    // RawQuerynameLexKey tests
    // ========================================================================

    #[test]
    fn test_lex_key_lexicographic_ordering() {
        // Lexicographic: "read10" < "read2" (because '1' < '2' in ASCII)
        let k1 = RawQuerynameLexKey::new(b"read10".to_vec(), 0);
        let k2 = RawQuerynameLexKey::new(b"read2".to_vec(), 0);
        assert!(k1 < k2, "lexicographic: 'read10' should be < 'read2'");
    }

    #[test]
    fn test_natural_key_natural_ordering() {
        // Natural: "read2" < "read10" (because 2 < 10 numerically)
        let k1 = RawQuerynameKey::new(b"read2".to_vec(), 0);
        let k2 = RawQuerynameKey::new(b"read10".to_vec(), 0);
        assert!(k1 < k2, "natural: 'read2' should be < 'read10'");
    }

    #[test]
    fn test_lex_vs_natural_ordering_difference() {
        // This is the key difference: natural treats "2" < "10", lexicographic treats "10" < "2"
        let lex_10 = RawQuerynameLexKey::new(b"read10".to_vec(), 0);
        let lex_2 = RawQuerynameLexKey::new(b"read2".to_vec(), 0);
        assert!(lex_10 < lex_2, "lexicographic: read10 < read2");

        let nat_2 = RawQuerynameKey::new(b"read2".to_vec(), 0);
        let nat_10 = RawQuerynameKey::new(b"read10".to_vec(), 0);
        assert!(nat_2 < nat_10, "natural: read2 < read10");
    }

    #[test]
    fn test_lex_key_flag_tiebreak() {
        let k1 = RawQuerynameLexKey::new(b"readA".to_vec(), 0);
        let k2 = RawQuerynameLexKey::new(b"readA".to_vec(), 1);
        assert!(k1 < k2);
    }

    #[test]
    fn test_lex_key_empty_names() {
        let k1 = RawQuerynameLexKey::new(Vec::new(), 0);
        let k2 = RawQuerynameLexKey::new(b"a".to_vec(), 0);
        assert!(k1 < k2);
    }

    #[test]
    fn test_lex_key_illumina_names() {
        // Illumina names: lexicographic ordering
        let mut names: Vec<RawQuerynameLexKey> = vec![
            RawQuerynameLexKey::new(b"A00132:53:HFHJKDSXX:2:1100:5000:1000".to_vec(), 0),
            RawQuerynameLexKey::new(b"A00132:53:HFHJKDSXX:1:1100:5000:1000".to_vec(), 0),
            RawQuerynameLexKey::new(b"A00132:54:HFH2JDSXX:1:1100:5000:1000".to_vec(), 0),
        ];
        names.sort();
        assert_eq!(names[0].name(), b"A00132:53:HFHJKDSXX:1:1100:5000:1000");
        assert_eq!(names[1].name(), b"A00132:53:HFHJKDSXX:2:1100:5000:1000");
        assert_eq!(names[2].name(), b"A00132:54:HFH2JDSXX:1:1100:5000:1000");
    }

    #[test]
    fn test_lex_key_srr_names_differ_from_natural() {
        // SRR names: lexicographic puts .100 before .2 (1 < 2 in ASCII)
        let mut lex_names: Vec<RawQuerynameLexKey> = vec![
            RawQuerynameLexKey::new(b"SRR099966.100".to_vec(), 0),
            RawQuerynameLexKey::new(b"SRR099966.2".to_vec(), 0),
            RawQuerynameLexKey::new(b"SRR099966.10".to_vec(), 0),
        ];
        lex_names.sort();
        assert_eq!(lex_names[0].name(), b"SRR099966.10");
        assert_eq!(lex_names[1].name(), b"SRR099966.100");
        assert_eq!(lex_names[2].name(), b"SRR099966.2");

        // Natural puts them in numeric order
        let mut nat_names: Vec<RawQuerynameKey> = vec![
            RawQuerynameKey::new(b"SRR099966.100".to_vec(), 0),
            RawQuerynameKey::new(b"SRR099966.2".to_vec(), 0),
            RawQuerynameKey::new(b"SRR099966.10".to_vec(), 0),
        ];
        nat_names.sort();
        assert_eq!(nat_names[0].name(), b"SRR099966.2\0");
        assert_eq!(nat_names[1].name(), b"SRR099966.10\0");
        assert_eq!(nat_names[2].name(), b"SRR099966.100\0");
    }

    #[test]
    fn test_lex_key_serialization_roundtrip() {
        let key = RawQuerynameLexKey::new(b"test_read".to_vec(), 42);
        let mut buf = Vec::new();
        key.write_to(&mut buf).unwrap();

        let mut cursor = std::io::Cursor::new(&buf);
        let restored = RawQuerynameLexKey::read_from(&mut cursor).unwrap();
        assert_eq!(key, restored);
    }

    #[test]
    fn test_lex_key_variable_length() {
        const { assert!(RawQuerynameLexKey::SERIALIZED_SIZE.is_none()) };
    }

    #[test]
    fn test_queryname_keys_embedded_in_record() {
        const { assert!(RawQuerynameLexKey::EMBEDDED_IN_RECORD) };
        const { assert!(RawQuerynameKey::EMBEDDED_IN_RECORD) };
        const { assert!(RawCoordinateKey::EMBEDDED_IN_RECORD) };
    }
}
