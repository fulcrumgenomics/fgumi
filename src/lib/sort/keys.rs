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
use std::hash::{Hash, Hasher};
use std::io::{Read, Write};

/// The MI (Molecular Identifier) tag.
const MI_TAG: Tag = Tag::new(b'M', b'I');

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
pub trait SortKey: Ord + Clone + Send + Sync {
    /// Extract a sort key from a BAM record.
    fn from_record(record: &RecordBuf, header: &Header) -> Result<Self>;
}

// ============================================================================
// Coordinate Sort Key
// ============================================================================

/// Sort key for coordinate ordering.
///
/// Sort order: reference ID → position → reverse strand flag.
/// Unmapped reads (tid = -1) are sorted to the end.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct CoordinateKey {
    /// Reference sequence ID (tid), or `i32::MAX` for unmapped.
    pub tid: i32,
    /// 0-based alignment start position.
    pub pos: i64,
    /// True if reverse strand.
    pub reverse: bool,
    /// Read name for tie-breaking (lexicographic).
    pub name: Vec<u8>,
}

impl CoordinateKey {
    /// Create a coordinate key for an unmapped read.
    #[must_use]
    pub fn unmapped(name: Vec<u8>) -> Self {
        Self { tid: i32::MAX, pos: i64::MAX, reverse: false, name }
    }
}

impl Ord for CoordinateKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.tid
            .cmp(&other.tid)
            .then_with(|| self.pos.cmp(&other.pos))
            .then_with(|| self.reverse.cmp(&other.reverse))
            .then_with(|| self.name.cmp(&other.name))
    }
}

impl PartialOrd for CoordinateKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl SortKey for CoordinateKey {
    fn from_record(record: &RecordBuf, _header: &Header) -> Result<Self> {
        let name =
            record.name().map_or_else(Vec::new, |n| Vec::from(<_ as AsRef<[u8]>>::as_ref(n)));

        if record.flags().is_unmapped() {
            return Ok(Self::unmapped(name));
        }

        #[allow(clippy::cast_possible_wrap)]
        let tid = record.reference_sequence_id().map_or(-1, |id| id as i32);

        #[allow(clippy::cast_possible_wrap)]
        let pos = record.alignment_start().map_or(0, |p| usize::from(p) as i64);

        let reverse = record.flags().is_reverse_complemented();

        Ok(Self { tid, pos, reverse, name })
    }
}

// ============================================================================
// Raw Coordinate Sort Key (Fixed-size for RawSortKey trait)
// ============================================================================

/// Fixed-size coordinate sort key for raw BAM sorting (16 bytes).
///
/// This key is designed for efficient temp file storage and O(1) comparisons
/// during merge phase. It packs:
/// - `sort_key`: (tid << 34) | ((pos+1) << 1) | reverse
/// - `name_hash`: Hash of read name for deterministic tie-breaking
///
/// Unlike [`CoordinateKey`], this type doesn't store the actual read name,
/// making it fixed-size (16 bytes) instead of variable-length.
#[repr(C)]
#[derive(Copy, Clone, Eq, PartialEq, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct RawCoordinateKey {
    /// Packed primary sort key: (tid << 34) | ((pos+1) << 1) | reverse.
    pub sort_key: u64,
    /// Hash of read name for deterministic tie-breaking.
    pub name_hash: u64,
}

impl RawCoordinateKey {
    /// Size in bytes when serialized.
    pub const SIZE: usize = 16;

    /// Create a new coordinate key from components.
    ///
    /// # Arguments
    /// * `tid` - Reference sequence ID (-1 for unmapped)
    /// * `pos` - 0-based alignment position
    /// * `reverse` - True if reverse complemented
    /// * `name_hash` - Hash of read name
    /// * `nref` - Number of reference sequences (for unmapped handling)
    #[inline]
    #[must_use]
    pub fn new(tid: i32, pos: i32, reverse: bool, name_hash: u64, nref: u32) -> Self {
        // Map unmapped (tid=-1) to nref for proper sorting (after all mapped)
        let tid = if tid < 0 { nref } else { tid as u32 };
        // Pack: tid in high bits, (pos+1) in middle, reverse in LSB
        let key = (u64::from(tid) << 34)
            | ((i64::from(pos) as u64).wrapping_add(1) << 1)
            | u64::from(reverse);
        Self { sort_key: key, name_hash }
    }

    /// Create a key for unmapped records (sorts after all mapped).
    #[inline]
    #[must_use]
    pub fn unmapped(name_hash: u64) -> Self {
        Self { sort_key: u64::MAX, name_hash }
    }

    /// Create a zeroed key (for memory operations).
    #[inline]
    #[must_use]
    pub fn zeroed() -> Self {
        Self { sort_key: 0, name_hash: 0 }
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
        self.sort_key.cmp(&other.sort_key).then_with(|| self.name_hash.cmp(&other.name_hash))
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
        // 8: l_read_name (u8)
        // 14-15: flags (u16)
        // 32+: read name (null-terminated)

        let tid = i32::from_le_bytes([bam[0], bam[1], bam[2], bam[3]]);
        let pos = i32::from_le_bytes([bam[4], bam[5], bam[6], bam[7]]);
        let flags = u16::from_le_bytes([bam[14], bam[15]]);
        let reverse = (flags & 0x10) != 0;

        // Hash read name for tie-breaking (exclude null terminator)
        let name_len = (bam[8] as usize).saturating_sub(1);
        let name =
            if name_len > 0 && 32 + name_len <= bam.len() { &bam[32..32 + name_len] } else { &[] };
        let name_hash = hash_name(name);

        // Create key based on tid (samtools behavior):
        // - tid >= 0: sort by (tid, pos, reverse) even if unmapped flag is set
        // - tid < 0: unmapped with no reference, sort at end
        if tid < 0 {
            Self::unmapped(name_hash)
        } else {
            Self::new(tid, pos, reverse, name_hash, ctx.nref)
        }
    }

    #[inline]
    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writer.write_all(&self.sort_key.to_le_bytes())?;
        writer.write_all(&self.name_hash.to_le_bytes())
    }

    #[inline]
    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let mut buf = [0u8; 16];
        reader.read_exact(&mut buf)?;
        Ok(Self {
            sort_key: u64::from_le_bytes([
                buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7],
            ]),
            name_hash: u64::from_le_bytes([
                buf[8], buf[9], buf[10], buf[11], buf[12], buf[13], buf[14], buf[15],
            ]),
        })
    }
}

/// Fast name hashing using ahash.
#[inline]
fn hash_name(name: &[u8]) -> u64 {
    let mut hasher = ahash::AHasher::default();
    name.hash(&mut hasher);
    hasher.finish()
}

// ============================================================================
// Queryname Sort Key
// ============================================================================

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
    fn from_record(record: &RecordBuf, _header: &Header) -> Result<Self> {
        let name =
            record.name().map_or_else(Vec::new, |n| Vec::from(<_ as AsRef<[u8]>>::as_ref(n)));
        let flags = u16::from(record.flags());
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
        let flags = u16::from_le_bytes([bam[14], bam[15]]);

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
/// required for `fgumi group` input.
///
/// Sort order:
/// 1. tid1 - reference ID of earlier read
/// 2. tid2 - reference ID of later read
/// 3. pos1 - unclipped 5' position of earlier read
/// 4. pos2 - unclipped 5' position of later read
/// 5. neg1 - strand of earlier read (reverse before forward in samtools)
/// 6. neg2 - strand of later read
/// 7. MI tag (length first, then lexicographic, stripping /A /B suffix)
/// 8. Read name
/// 9. Whether this is the upper read of the pair
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
    /// True if this is the "upper" (later) read of the pair.
    pub is_upper: bool,
}

impl TemplateCoordinateKey {
    /// Create a key for an unmapped read.
    #[must_use]
    pub fn unmapped(name: Vec<u8>, mid: String) -> Self {
        Self {
            tid1: i32::MAX,
            tid2: i32::MAX,
            pos1: i64::MAX,
            pos2: i64::MAX,
            neg1: false,
            neg2: false,
            mid,
            name,
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
    fn from_record(record: &RecordBuf, _header: &Header) -> Result<Self> {
        let name =
            record.name().map_or_else(Vec::new, |n| Vec::from(<_ as AsRef<[u8]>>::as_ref(n)));

        // Extract MI tag
        let mid = extract_mi_tag(record);

        // Handle unmapped reads
        if record.flags().is_unmapped() {
            return Ok(Self::unmapped(name, mid));
        }

        let flags = record.flags();
        let is_reverse = flags.is_reverse_complemented();

        // Get reference ID
        #[allow(clippy::cast_possible_wrap)]
        let tid = record.reference_sequence_id().map_or(-1, |id| id as i32);

        // Calculate unclipped 5' position for this read
        let this_pos = get_unclipped_5prime_position(record)?;

        // Get mate information
        #[allow(clippy::cast_possible_wrap)]
        let mate_tid = record.mate_reference_sequence_id().map_or(i32::MAX, |id| id as i32);

        #[allow(clippy::cast_possible_wrap)]
        let mate_pos = record.mate_alignment_start().map_or(0, |p| usize::from(p) as i64);

        let mate_reverse = flags.is_mate_reverse_complemented();

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

        Ok(Self { tid1, tid2, pos1, pos2, neg1, neg2, mid, name, is_upper })
    }
}

/// Calculate the unclipped 5' position of a read.
///
/// For forward strand: unclipped start = alignment start - leading soft clips
/// For reverse strand: unclipped end = alignment start + alignment span + trailing soft clips
fn get_unclipped_5prime_position(record: &RecordBuf) -> Result<i64> {
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
        let k1 = CoordinateKey { tid: 0, pos: 100, reverse: false, name: vec![] };
        let k2 = CoordinateKey { tid: 0, pos: 200, reverse: false, name: vec![] };
        let k3 = CoordinateKey { tid: 1, pos: 50, reverse: false, name: vec![] };

        assert!(k1 < k2);
        assert!(k2 < k3);
    }

    #[test]
    fn test_coordinate_key_unmapped_last() {
        let mapped = CoordinateKey { tid: 0, pos: 100, reverse: false, name: vec![] };
        let unmapped = CoordinateKey::unmapped(vec![]);

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
            is_upper: false,
        };

        assert!(k_rev < k_fwd);
    }
}
