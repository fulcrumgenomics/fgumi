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
//! - [`TemplateCoordinateKey`]: Template-level position for UMI grouping
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
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::cmp::Ordering;

use crate::read_info::LibraryIndex;
use crate::sam::record_utils::{mate_unclipped_end, mate_unclipped_start};
use noodles::sam::alignment::record_buf::data::field::value::Value as DataValue;
use std::io::{Read, Write};

/// The MI (Molecular Identifier) tag.
const MI_TAG: Tag = Tag::new(b'M', b'I');

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

    /// Extract a sort key from raw BAM record bytes.
    ///
    /// This is the hot path during sorting - implementations should minimize
    /// parsing overhead by reading only the fields needed for comparison.
    fn extract(bam: &[u8], ctx: &SortContext) -> Self;

    /// Serialize the key to a writer for temp file storage.
    ///
    /// Format should be compact and enable fast deserialization.
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()>;

    /// Deserialize a key from a reader.
    ///
    /// Must be the inverse of `write_to`.
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
    pub fn from_header(header: &Header) -> Self {
        Self { nref: header.reference_sequences().len() as u32 }
    }

    /// Create a context with explicit nref (for testing).
    #[must_use]
    pub fn new(nref: u32) -> Self {
        Self { nref }
    }
}

/// Sort order enumeration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    /// Coordinate sort: tid → pos → reverse strand
    Coordinate,
    /// Queryname sort: read name with natural ordering
    Queryname,
    /// Template-coordinate sort: template position for UMI grouping
    TemplateCoordinate,
}

impl SortOrder {
    /// Get the SAM header sort order tag value.
    #[must_use]
    pub fn header_so_tag(&self) -> &'static str {
        match self {
            Self::Coordinate => "coordinate",
            Self::Queryname => "queryname",
            Self::TemplateCoordinate => "unsorted",
        }
    }

    /// Get the SAM header group order tag value.
    #[must_use]
    pub fn header_go_tag(&self) -> Option<&'static str> {
        match self {
            Self::Coordinate => None,
            Self::Queryname => None,
            Self::TemplateCoordinate => Some("query"),
        }
    }

    /// Get the SAM header sub-sort tag value.
    #[must_use]
    pub fn header_ss_tag(&self) -> Option<&'static str> {
        match self {
            Self::Coordinate => None,
            Self::Queryname => None,
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

        #[allow(clippy::cast_possible_wrap)]
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

    #[inline]
    fn extract(bam: &[u8], ctx: &SortContext) -> Self {
        // BAM format offsets (all little-endian):
        // 0-3: tid (i32)
        // 4-7: pos (i32)
        // 14-15: flags (u16)

        let tid = i32::from_le_bytes([bam[0], bam[1], bam[2], bam[3]]);
        let pos = i32::from_le_bytes([bam[4], bam[5], bam[6], bam[7]]);
        let flags = u16::from_le_bytes([bam[14], bam[15]]);
        let reverse = (flags & 0x10) != 0;

        // Create key based on tid (samtools behavior):
        // - tid >= 0: sort by (tid, pos, reverse) even if unmapped flag is set
        // - tid < 0: unmapped with no reference, sort at end
        if tid < 0 { Self::unmapped() } else { Self::new(tid, pos, reverse, ctx.nref) }
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
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct QuerynameKey {
    /// Read name bytes.
    pub name: Vec<u8>,
    /// Read pair flags for ordering R1 before R2.
    pub flags: u16,
}

impl Ord for QuerynameKey {
    fn cmp(&self, other: &Self) -> Ordering {
        natural_compare(&self.name, &other.name).then_with(|| self.flags.cmp(&other.flags))
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
        let name =
            record.name().map_or_else(Vec::new, |n| Vec::from(<_ as AsRef<[u8]>>::as_ref(n)));
        let flags = queryname_flag_order(u16::from(record.flags()));
        Ok(Self { name, flags })
    }
}

/// Natural string comparison that handles numeric runs.
///
/// Compares strings such that "read1" < "read2" < "read10".
fn natural_compare(a: &[u8], b: &[u8]) -> Ordering {
    let mut i = 0;
    let mut j = 0;

    while i < a.len() && j < b.len() {
        let a_digit = a[i].is_ascii_digit();
        let b_digit = b[j].is_ascii_digit();

        match (a_digit, b_digit) {
            (true, true) => {
                // Both are digits - compare numerically
                let (a_num, a_end) = parse_number(&a[i..]);
                let (b_num, b_end) = parse_number(&b[j..]);

                match a_num.cmp(&b_num) {
                    Ordering::Equal => {
                        i += a_end;
                        j += b_end;
                    }
                    ord => return ord,
                }
            }
            (true, false) => return Ordering::Less, // Digits before non-digits
            (false, true) => return Ordering::Greater,
            (false, false) => {
                // Both non-digits - compare bytes
                match a[i].cmp(&b[j]) {
                    Ordering::Equal => {
                        i += 1;
                        j += 1;
                    }
                    ord => return ord,
                }
            }
        }
    }

    // Shorter string sorts first if one is prefix of other
    a.len().cmp(&b.len())
}

/// Parse a numeric run from the start of a byte slice.
/// Returns (number, bytes consumed).
fn parse_number(bytes: &[u8]) -> (u64, usize) {
    let mut num: u64 = 0;
    let mut i = 0;

    while i < bytes.len() && bytes[i].is_ascii_digit() {
        num = num.saturating_mul(10).saturating_add(u64::from(bytes[i] - b'0'));
        i += 1;
    }

    (num, i)
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
#[derive(Clone, Eq, PartialEq, Debug, Default)]
pub struct RawQuerynameKey {
    /// Read name bytes.
    pub name: Vec<u8>,
    /// Flags for segment ordering (R1 before R2).
    pub flags: u16,
}

impl RawQuerynameKey {
    /// Create a new queryname key.
    #[must_use]
    pub fn new(name: Vec<u8>, flags: u16) -> Self {
        Self { name, flags }
    }
}

impl Ord for RawQuerynameKey {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        natural_compare(&self.name, &other.name).then_with(|| self.flags.cmp(&other.flags))
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

    #[inline]
    fn extract(bam: &[u8], _ctx: &SortContext) -> Self {
        // BAM format offsets:
        // 8: l_read_name (u8)
        // 14-15: flags (u16)
        // 32+: read name (null-terminated)

        let name_len = (bam[8] as usize).saturating_sub(1); // Exclude null terminator
        let name = if name_len > 0 && 32 + name_len <= bam.len() {
            bam[32..32 + name_len].to_vec()
        } else {
            Vec::new()
        };
        let raw_flags = u16::from_le_bytes([bam[14], bam[15]]);
        let flags = queryname_flag_order(raw_flags);

        Self { name, flags }
    }

    #[inline]
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        // Format: [name_len: u16][name: bytes][flags: u16]
        let name_len = self.name.len() as u16;
        writer.write_all(&name_len.to_le_bytes())?;
        writer.write_all(&self.name)?;
        writer.write_all(&self.flags.to_le_bytes())
    }

    #[inline]
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        // Read name length
        let mut len_buf = [0u8; 2];
        reader.read_exact(&mut len_buf)?;
        let name_len = u16::from_le_bytes(len_buf) as usize;

        // Read name
        let mut name = vec![0u8; name_len];
        reader.read_exact(&mut name)?;

        // Read flags
        let mut flags_buf = [0u8; 2];
        reader.read_exact(&mut flags_buf)?;
        let flags = u16::from_le_bytes(flags_buf);

        Ok(Self { name, flags })
    }
}

// ============================================================================
// Template-Coordinate Sort Key
// ============================================================================

/// Sort key for template-coordinate ordering.
///
/// Matches the sort order used by `samtools sort --template-coordinate` and
/// fgbio's `TemplateCoordinate` order for `fgumi group` input.
///
/// Sort order (matching fgbio):
/// 1. tid1 - reference ID of earlier read
/// 2. tid2 - reference ID of later read
/// 3. pos1 - unclipped 5' position of earlier read
/// 4. pos2 - unclipped 5' position of later read
/// 5. neg1 - strand of earlier read (reverse before forward in samtools)
/// 6. neg2 - strand of later read
/// 7. MI tag (length first, then lexicographic, stripping /A /B suffix)
/// 8. Read name
/// 9. Library (ordinal index from header's @RG LB field)
/// 10. Whether this is the upper read of the pair
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct TemplateCoordinateKey {
    /// Reference ID of earlier read (`i32::MAX` for unmapped).
    pub tid1: i32,
    /// Reference ID of later read (`i32::MAX` for unmapped).
    pub tid2: i32,
    /// Unclipped 5' position of earlier read.
    pub pos1: i64,
    /// Unclipped 5' position of later read.
    pub pos2: i64,
    /// Strand of earlier read (true = reverse).
    pub neg1: bool,
    /// Strand of later read (true = reverse).
    pub neg2: bool,
    /// Molecular identifier (without /A /B suffix).
    pub mid: String,
    /// Read name.
    pub name: Vec<u8>,
    /// Library ordinal index (from `LibraryIndex`, 0 = unknown).
    pub library_idx: u16,
    /// True if this is the "upper" (later) read of the pair.
    pub is_upper: bool,
}

impl TemplateCoordinateKey {
    /// Create a key for an unmapped read.
    #[must_use]
    pub fn unmapped(name: Vec<u8>, mid: String, library_idx: u16) -> Self {
        Self {
            tid1: i32::MAX,
            tid2: i32::MAX,
            pos1: i64::MAX,
            pos2: i64::MAX,
            neg1: false,
            neg2: false,
            mid,
            name,
            library_idx,
            is_upper: false,
        }
    }
}

impl Ord for TemplateCoordinateKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.tid1
            .cmp(&other.tid1)
            .then_with(|| self.tid2.cmp(&other.tid2))
            .then_with(|| self.pos1.cmp(&other.pos1))
            .then_with(|| self.pos2.cmp(&other.pos2))
            // samtools: reverse (neg=true) sorts before forward (neg=false)
            .then_with(|| match (self.neg1, other.neg1) {
                (true, false) => Ordering::Less,
                (false, true) => Ordering::Greater,
                _ => Ordering::Equal,
            })
            .then_with(|| match (self.neg2, other.neg2) {
                (true, false) => Ordering::Less,
                (false, true) => Ordering::Greater,
                _ => Ordering::Equal,
            })
            .then_with(|| compare_mid(&self.mid, &other.mid))
            .then_with(|| self.name.cmp(&other.name))
            .then_with(|| self.library_idx.cmp(&other.library_idx))
            .then_with(|| self.is_upper.cmp(&other.is_upper))
    }
}

impl PartialOrd for TemplateCoordinateKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Compare MI tags using the samtools/fgbio algorithm.
///
/// The comparison:
/// 1. Strips trailing "/X" suffix (e.g., "/A", "/B") if present
/// 2. Compares by length first (shorter sorts before longer)
/// 3. Then compares lexicographically
fn compare_mid(mid1: &str, mid2: &str) -> Ordering {
    // Strip trailing "/X" suffix if present
    let mid1 = strip_mi_suffix(mid1);
    let mid2 = strip_mi_suffix(mid2);

    // Compare by length first, then lexicographically
    mid1.len().cmp(&mid2.len()).then_with(|| mid1.cmp(mid2))
}

/// Strip the trailing /A or /B suffix from an MI tag.
#[inline]
fn strip_mi_suffix(mid: &str) -> &str {
    if mid.len() >= 2 && mid.as_bytes()[mid.len() - 2] == b'/' {
        &mid[..mid.len() - 2]
    } else {
        mid
    }
}

impl SortKey for TemplateCoordinateKey {
    type Context = LibraryIndex;

    fn build_context(header: &Header) -> Self::Context {
        LibraryIndex::from_header(header)
    }

    fn from_record(record: &RecordBuf, _header: &Header, ctx: &Self::Context) -> Result<Self> {
        let name =
            record.name().map_or_else(Vec::new, |n| Vec::from(<_ as AsRef<[u8]>>::as_ref(n)));

        // Extract MI tag
        let mid = extract_mi_tag(record);

        // Extract library index from RG tag (O(1) hash lookup)
        let library_idx = if let Some(DataValue::String(rg_bytes)) = record.data().get(b"RG") {
            let rg_hash = LibraryIndex::hash_rg(rg_bytes);
            ctx.get(rg_hash)
        } else {
            0 // unknown library
        };

        // Handle unmapped reads
        if record.flags().is_unmapped() {
            return Ok(Self::unmapped(name, mid, library_idx));
        }

        let flags = record.flags();

        // For secondary/supplementary reads, check for pa tag to get template sort key
        // This ensures they sort with their primary alignment
        if flags.is_secondary() || flags.is_supplementary() {
            if let Some(pa_info) = PrimaryAlignmentInfo::from_record(record) {
                // Use primary alignment's template sort key
                return Ok(Self {
                    tid1: pa_info.tid1,
                    tid2: pa_info.tid2,
                    pos1: i64::from(pa_info.pos1),
                    pos2: i64::from(pa_info.pos2),
                    neg1: pa_info.neg1,
                    neg2: pa_info.neg2,
                    mid,
                    name,
                    library_idx,
                    is_upper: false,
                });
            }
            // If no pa tag, fall through to use the read's own coordinates
        }

        let is_reverse = flags.is_reverse_complemented();

        // Get reference ID
        #[allow(clippy::cast_possible_wrap)]
        let tid = record.reference_sequence_id().map_or(-1, |id| id as i32);

        // Calculate unclipped 5' position for this read
        let this_pos = get_unclipped_5prime_position(record)?;

        // Get mate information
        #[allow(clippy::cast_possible_wrap)]
        let mate_tid = record.mate_reference_sequence_id().map_or(i32::MAX, |id| id as i32);

        let mate_reverse = flags.is_mate_reverse_complemented();

        // Calculate mate's unclipped 5' position from MC tag
        // This ensures R1 and R2 compute the same template sort key
        #[allow(clippy::cast_possible_wrap)]
        let mate_pos: i64 = if mate_reverse {
            // Reverse strand - 5' is at the end
            mate_unclipped_end(record).map(|p| p as i64).unwrap_or_else(|| {
                // Fallback to raw mate position if MC tag missing
                record.mate_alignment_start().map_or(0, |p| usize::from(p) as i64)
            })
        } else {
            // Forward strand - 5' is at the start
            mate_unclipped_start(record).map(|p| p as i64).unwrap_or_else(|| {
                // Fallback to raw mate position if MC tag missing
                record.mate_alignment_start().map_or(0, |p| usize::from(p) as i64)
            })
        };

        // Determine which read is "earlier" (lower position)
        // For single-end reads or unmapped mates, use this read's position
        let (tid1, tid2, pos1, pos2, neg1, neg2, is_upper) =
            if flags.is_segmented() && !flags.is_mate_unmapped() {
                // Paired read with mapped mate
                if (tid, this_pos) <= (mate_tid, mate_pos) {
                    // This read is earlier
                    (tid, mate_tid, this_pos, mate_pos, is_reverse, mate_reverse, false)
                } else {
                    // Mate is earlier
                    (mate_tid, tid, mate_pos, this_pos, mate_reverse, is_reverse, true)
                }
            } else {
                // Single-end or unmapped mate
                (tid, tid, this_pos, this_pos, is_reverse, is_reverse, false)
            };

        Ok(Self { tid1, tid2, pos1, pos2, neg1, neg2, mid, name, library_idx, is_upper })
    }
}

/// Calculate the unclipped 5' position of a read.
///
/// For forward strand: unclipped start = alignment start - leading soft clips
/// For reverse strand: unclipped end = alignment start + alignment span + trailing soft clips
///
/// This is the position used by template-coordinate sorting to determine read order.
pub fn get_unclipped_5prime_position(record: &RecordBuf) -> Result<i64> {
    if record.flags().is_unmapped() {
        return Ok(0);
    }

    let Some(alignment_start) = record.alignment_start() else {
        return Ok(0);
    };

    #[allow(clippy::cast_possible_wrap)]
    let alignment_start_i64 = usize::from(alignment_start) as i64;

    if record.flags().is_reverse_complemented() {
        // Negative strand: calculate unclipped end (5' position)
        let cigar = record.cigar();

        // Calculate alignment span (M/D/N/=/X operations)
        #[allow(clippy::cast_possible_wrap)]
        let alignment_span: i64 = cigar
            .as_ref()
            .iter()
            .filter_map(|op| match op.kind() {
                Kind::Match
                | Kind::Deletion
                | Kind::Skip
                | Kind::SequenceMatch
                | Kind::SequenceMismatch => Some(op.len() as i64),
                _ => None,
            })
            .sum();

        // Get trailing soft clips (at the end of CIGAR = 5' end for reverse strand)
        #[allow(clippy::cast_possible_wrap)]
        let trailing_soft_clips: i64 = cigar
            .as_ref()
            .iter()
            .rev()
            .take_while(|op| op.kind() == Kind::SoftClip || op.kind() == Kind::HardClip)
            .filter(|op| op.kind() == Kind::SoftClip)
            .map(|op| op.len() as i64)
            .sum();

        Ok(alignment_start_i64 + alignment_span + trailing_soft_clips - 1)
    } else {
        // Positive strand: calculate unclipped start (5' position)
        let cigar = record.cigar();

        // Get leading soft clips
        #[allow(clippy::cast_possible_wrap)]
        let leading_soft_clips: i64 = cigar
            .as_ref()
            .iter()
            .take_while(|op| op.kind() == Kind::SoftClip || op.kind() == Kind::HardClip)
            .filter(|op| op.kind() == Kind::SoftClip)
            .map(|op| op.len() as i64)
            .sum();

        Ok(alignment_start_i64 - leading_soft_clips)
    }
}

/// Extract the MI tag value from a record, returning empty string if not present.
fn extract_mi_tag(record: &RecordBuf) -> String {
    record
        .data()
        .get(&MI_TAG)
        .map(|v| {
            use noodles::sam::alignment::record_buf::data::field::Value;
            match v {
                Value::Int8(v) => v.to_string(),
                Value::UInt8(v) => v.to_string(),
                Value::Int16(v) => v.to_string(),
                Value::UInt16(v) => v.to_string(),
                Value::Int32(v) => v.to_string(),
                Value::UInt32(v) => v.to_string(),
                Value::String(s) => s.to_string(),
                _ => String::new(),
            }
        })
        .unwrap_or_default()
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_natural_compare_basic() {
        assert_eq!(natural_compare(b"a", b"b"), Ordering::Less);
        assert_eq!(natural_compare(b"b", b"a"), Ordering::Greater);
        assert_eq!(natural_compare(b"a", b"a"), Ordering::Equal);
    }

    #[test]
    fn test_natural_compare_numeric() {
        assert_eq!(natural_compare(b"read1", b"read2"), Ordering::Less);
        assert_eq!(natural_compare(b"read2", b"read10"), Ordering::Less);
        assert_eq!(natural_compare(b"read10", b"read11"), Ordering::Less);
        assert_eq!(natural_compare(b"read9", b"read10"), Ordering::Less);
    }

    #[test]
    fn test_natural_compare_mixed() {
        assert_eq!(natural_compare(b"a1b2", b"a1b10"), Ordering::Less);
        assert_eq!(natural_compare(b"x10y20", b"x10y3"), Ordering::Greater);
    }

    #[test]
    fn test_natural_compare_empty() {
        assert_eq!(natural_compare(b"", b""), Ordering::Equal);
        assert_eq!(natural_compare(b"", b"a"), Ordering::Less);
        assert_eq!(natural_compare(b"a", b""), Ordering::Greater);
    }

    #[test]
    fn test_compare_mid_strips_suffix() {
        assert_eq!(compare_mid("123/A", "123/B"), Ordering::Equal);
        assert_eq!(compare_mid("123/A", "123"), Ordering::Equal);
    }

    #[test]
    fn test_compare_mid_length_first() {
        assert_eq!(compare_mid("1", "12"), Ordering::Less);
        assert_eq!(compare_mid("99", "100"), Ordering::Less);
    }

    #[test]
    fn test_compare_mid_lexicographic() {
        assert_eq!(compare_mid("abc", "abd"), Ordering::Less);
        assert_eq!(compare_mid("12", "13"), Ordering::Less);
    }

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
    fn test_template_coord_key_ordering() {
        let k1 = TemplateCoordinateKey {
            tid1: 0,
            tid2: 0,
            pos1: 100,
            pos2: 200,
            neg1: false,
            neg2: true,
            mid: String::new(),
            name: vec![],
            library_idx: 0,
            is_upper: false,
        };
        let k2 = TemplateCoordinateKey {
            tid1: 0,
            tid2: 0,
            pos1: 200,
            pos2: 300,
            neg1: false,
            neg2: true,
            mid: String::new(),
            name: vec![],
            library_idx: 0,
            is_upper: false,
        };

        assert!(k1 < k2);
    }

    #[test]
    fn test_template_coord_key_strand_ordering() {
        // samtools: reverse strand sorts before forward
        let k_rev = TemplateCoordinateKey {
            tid1: 0,
            tid2: 0,
            pos1: 100,
            pos2: 200,
            neg1: true,
            neg2: false,
            mid: String::new(),
            name: vec![],
            library_idx: 0,
            is_upper: false,
        };
        let k_fwd = TemplateCoordinateKey {
            tid1: 0,
            tid2: 0,
            pos1: 100,
            pos2: 200,
            neg1: false,
            neg2: false,
            mid: String::new(),
            name: vec![],
            library_idx: 0,
            is_upper: false,
        };

        assert!(k_rev < k_fwd);
    }

    #[test]
    fn test_template_coord_key_library_ordering() {
        // Library comes after name in sort order
        let k1 = TemplateCoordinateKey {
            tid1: 0,
            tid2: 0,
            pos1: 100,
            pos2: 200,
            neg1: false,
            neg2: false,
            mid: String::new(),
            name: b"read1".to_vec(),
            library_idx: 1,
            is_upper: false,
        };
        let k2 = TemplateCoordinateKey {
            tid1: 0,
            tid2: 0,
            pos1: 100,
            pos2: 200,
            neg1: false,
            neg2: false,
            mid: String::new(),
            name: b"read1".to_vec(),
            library_idx: 2,
            is_upper: false,
        };

        assert!(k1 < k2);
    }

    /// Comprehensive regression test for `TemplateCoordinateKey` comparison order.
    /// Verifies each field is compared in the correct order:
    /// `tid1` → `tid2` → `pos1` → `pos2` → `neg1` → `neg2` → `mid` → `name` → `library_idx` → `is_upper`
    #[test]
    fn test_template_coord_key_comparison_order_regression() {
        // Helper to create a key with specific values
        #[allow(clippy::too_many_arguments)]
        fn key(
            tid1: i32,
            tid2: i32,
            pos1: i64,
            pos2: i64,
            neg1: bool,
            neg2: bool,
            mid: &str,
            name: &[u8],
            library_idx: u16,
            is_upper: bool,
        ) -> TemplateCoordinateKey {
            TemplateCoordinateKey {
                tid1,
                tid2,
                pos1,
                pos2,
                neg1,
                neg2,
                mid: mid.to_string(),
                name: name.to_vec(),
                library_idx,
                is_upper,
            }
        }

        // Baseline key - all zeros/empty/false
        let baseline = key(0, 0, 0, 0, false, false, "", b"", 0, false);

        // 1. tid1 determines ordering when different
        let higher_tid1 = key(1, 0, 0, 0, false, false, "", b"", 0, false);
        assert!(baseline < higher_tid1, "tid1 should determine ordering");

        // 2. tid2 determines ordering when tid1 matches
        let higher_tid2 = key(0, 1, 0, 0, false, false, "", b"", 0, false);
        assert!(baseline < higher_tid2, "tid2 should determine ordering when tid1 matches");

        // 3. pos1 determines ordering when tid1, tid2 match
        let higher_pos1 = key(0, 0, 1, 0, false, false, "", b"", 0, false);
        assert!(baseline < higher_pos1, "pos1 should determine ordering when tids match");

        // 4. pos2 determines ordering when tid1, tid2, pos1 match
        let higher_pos2 = key(0, 0, 0, 1, false, false, "", b"", 0, false);
        assert!(baseline < higher_pos2, "pos2 should determine ordering when tids and pos1 match");

        // 5. neg1 determines ordering when tid1, tid2, pos1, pos2 match
        // samtools ordering: reverse (true) sorts BEFORE forward (false)
        let neg1_true = key(0, 0, 0, 0, true, false, "", b"", 0, false);
        assert!(
            neg1_true < baseline,
            "neg1=true should sort before neg1=false (samtools ordering)"
        );

        // 6. neg2 determines ordering when tid1, tid2, pos1, pos2, neg1 match
        let neg2_true = key(0, 0, 0, 0, false, true, "", b"", 0, false);
        assert!(
            neg2_true < baseline,
            "neg2=true should sort before neg2=false (samtools ordering)"
        );

        // 7. mid determines ordering when positions and strands match
        // MI comparison: length first, then lexicographic
        let longer_mid = key(0, 0, 0, 0, false, false, "aa", b"", 0, false);
        let shorter_mid = key(0, 0, 0, 0, false, false, "a", b"", 0, false);
        assert!(shorter_mid < longer_mid, "shorter MI should sort before longer MI");

        let mid_a = key(0, 0, 0, 0, false, false, "a", b"", 0, false);
        let mid_b = key(0, 0, 0, 0, false, false, "b", b"", 0, false);
        assert!(mid_a < mid_b, "MI comparison should be lexicographic when same length");

        // 8. name determines ordering when mid matches
        let name_a = key(0, 0, 0, 0, false, false, "", b"a", 0, false);
        let name_b = key(0, 0, 0, 0, false, false, "", b"b", 0, false);
        assert!(name_a < name_b, "name should determine ordering when mid matches");

        // 9. library_idx determines ordering when name matches
        let lib_0 = key(0, 0, 0, 0, false, false, "", b"a", 0, false);
        let lib_1 = key(0, 0, 0, 0, false, false, "", b"a", 1, false);
        assert!(lib_0 < lib_1, "library_idx should determine ordering when name matches");

        // 10. is_upper determines ordering when everything else matches
        let is_upper_false = key(0, 0, 0, 0, false, false, "", b"a", 0, false);
        let is_upper_true = key(0, 0, 0, 0, false, false, "", b"a", 0, true);
        assert!(
            is_upper_false < is_upper_true,
            "is_upper should determine ordering when all else matches"
        );

        // Verify that earlier fields take precedence over later fields
        // tid1 difference should override all other fields
        let k_low_tid1 = key(0, 10, 10, 10, true, true, "zzz", b"zzz", 10, true);
        let k_high_tid1 = key(1, 0, 0, 0, false, false, "", b"", 0, false);
        assert!(k_low_tid1 < k_high_tid1, "tid1 difference should override all later fields");

        // pos1 difference should override later fields (but not tid1/tid2)
        let k_low_pos1 = key(0, 0, 0, 10, true, true, "zzz", b"zzz", 10, true);
        let k_high_pos1 = key(0, 0, 1, 0, false, false, "", b"", 0, false);
        assert!(k_low_pos1 < k_high_pos1, "pos1 difference should override later fields");
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
            .pa_tag(info)
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

    #[test]
    fn test_template_coord_key_uses_pa_tag_for_supplementary() {
        use crate::sam::builder::RecordBuilder;
        use noodles::sam::alignment::record::Flags;

        // Create a supplementary record at chr5:5000 with pa tag pointing to chr0:100/200
        let pa_info = PrimaryAlignmentInfo::new(0, 100, false, 0, 200, true);
        let record = RecordBuilder::mapped_read()
            .name("test")
            .sequence("ACGT")
            .reference_sequence_id(5)
            .alignment_start(5000)
            .flags(Flags::SUPPLEMENTARY)
            .pa_tag(pa_info)
            .build();

        let header = Header::default();
        let ctx = TemplateCoordinateKey::build_context(&header);
        let key = TemplateCoordinateKey::from_record(&record, &header, &ctx).unwrap();

        // Key should use pa tag coordinates, not actual coordinates (chr5:5000)
        assert_eq!(key.tid1, 0);
        assert_eq!(key.tid2, 0);
        assert_eq!(key.pos1, 100);
        assert_eq!(key.pos2, 200);
        assert!(!key.neg1);
        assert!(key.neg2);
    }

    #[test]
    fn test_template_coord_key_supplementary_without_pa_tag_uses_own_coords() {
        use crate::sam::builder::RecordBuilder;
        use noodles::sam::alignment::record::Flags;

        // Create a supplementary record without pa tag
        let record = RecordBuilder::mapped_read()
            .name("test")
            .sequence("ACGT")
            .reference_sequence_id(5)
            .alignment_start(5000)
            .flags(Flags::SUPPLEMENTARY)
            .build();

        let header = Header::default();
        let ctx = TemplateCoordinateKey::build_context(&header);
        let key = TemplateCoordinateKey::from_record(&record, &header, &ctx).unwrap();

        // Key should use own coordinates since no pa tag
        assert_eq!(key.tid1, 5);
        // pos1 should be the unclipped position (5000 with 4M cigar = still 5000)
        assert_eq!(key.pos1, 5000);
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
        assert!(f113 < f2161, "R1 PRIMARY (113) should be < R1 SUPP (2161): {} vs {}", f113, f2161);
        assert!(f2161 < f177, "R1 SUPP (2161) should be < R2 PRIMARY (177): {} vs {}", f2161, f177);
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
}
