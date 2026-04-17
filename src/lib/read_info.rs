//! Data structures for tracking read position information.
//!
//! This module provides the `ReadInfo` structure for tracking genomic position information
//! about reads, which is essential for:
//! - Grouping reads by genomic location
//! - Ordering reads for duplicate marking
//! - Organizing reads for UMI-based consensus calling
//!
//! `ReadInfo` captures the unclipped 5' positions of paired-end reads (or single-end reads),
//! along with library and cell barcode information for proper grouping.

use anyhow::Result;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::header::Header;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::sync::Arc;

use bstr::ByteSlice;

use crate::sam::SamTag;
use crate::template::Template;

/// A lookup table mapping read group IDs to library names.
///
/// This is built from the SAM header's @RG lines and used to resolve the library
/// name (LB field) from a record's RG tag. This matches fgbio's behavior where
/// grouping uses the library name, not the read group ID.
///
/// Uses `Arc<str>` for library names to avoid cloning strings for every read.
///
/// # Note: `LibraryLookup` vs `LibraryIndex`
///
/// Both `LibraryLookup` and [`LibraryIndex`] exist for different use cases:
/// - `LibraryLookup`: String-based lookup returning library names. Used by
///   [`ReadInfo::from`] where the actual library name string is needed.
/// - [`LibraryIndex`]: Hash-based lookup returning numeric indices. Used by
///   `compute_group_key_from_raw` in the hot path where only equality comparison
///   matters, avoiding string allocations.
pub type LibraryLookup = Arc<HashMap<String, Arc<str>>>;

/// Shared "unknown" library string to avoid repeated allocations.
static UNKNOWN_LIBRARY: std::sync::LazyLock<Arc<str>> =
    std::sync::LazyLock::new(|| Arc::from("unknown"));

/// Builds a library lookup table from a SAM header.
///
/// Iterates through all @RG lines in the header and creates a mapping from
/// read group ID to library name. If a read group has no LB field, it maps
/// to "unknown" (matching fgbio's behavior).
///
/// # Arguments
///
/// * `header` - The SAM header containing @RG lines
///
/// # Returns
///
/// An `Arc<HashMap>` mapping read group IDs to library names
#[must_use]
pub fn build_library_lookup(header: &Header) -> LibraryLookup {
    let mut lookup = HashMap::new();

    for (id, rg) in header.read_groups() {
        // Get the LB field from the read group's other_fields
        let library: Arc<str> = rg
            .other_fields()
            .get(&rg_tag::LIBRARY)
            .map_or_else(|| Arc::clone(&UNKNOWN_LIBRARY), |s| Arc::from(s.to_string()));
        lookup.insert(id.to_string(), library);
    }

    Arc::new(lookup)
}

// ============================================================================
// LibraryIndex - Fast RG to library index mapping for GroupKey computation
// ============================================================================

/// Fast lookup from RG tag value to library index.
///
/// This provides O(1) library lookup during Decode using string hashing,
/// replacing the O(n) string comparison in the original `LibraryLookup`.
#[derive(Debug, Clone)]
pub struct LibraryIndex {
    /// Map from RG string hash to library index.
    lookup: ahash::AHashMap<u64, u16>,
    /// Library names for each index (for output/debugging).
    names: Vec<Arc<str>>,
    /// Unknown library index (always 0).
    unknown_idx: u16,
}

impl LibraryIndex {
    /// Build a library index from a SAM header.
    ///
    /// Each unique library name gets a sequential index starting from 0.
    /// Index 0 is reserved for "unknown" library.
    ///
    /// # Panics
    ///
    /// Panics if the header contains more than 65,535 distinct libraries.
    #[must_use]
    pub fn from_header(header: &Header) -> Self {
        let mut lookup = ahash::AHashMap::new();
        let mut names = vec![Arc::clone(&UNKNOWN_LIBRARY)]; // Index 0 = unknown
        let mut library_to_idx: ahash::AHashMap<Arc<str>, u16> = ahash::AHashMap::new();
        library_to_idx.insert(Arc::clone(&UNKNOWN_LIBRARY), 0);

        for (id, rg) in header.read_groups() {
            // Get library name from LB field
            let library: Arc<str> = rg
                .other_fields()
                .get(&rg_tag::LIBRARY)
                .map_or_else(|| Arc::clone(&UNKNOWN_LIBRARY), |s| Arc::from(s.to_string()));

            // Get or create library index
            let lib_idx = *library_to_idx.entry(library.clone()).or_insert_with(|| {
                let idx: u16 =
                    names.len().try_into().expect("too many distinct libraries for u16 index");
                names.push(library);
                idx
            });

            // Hash the RG string and map to library index
            let rg_hash = Self::hash_rg(id.as_bytes());
            lookup.insert(rg_hash, lib_idx);
        }

        Self { lookup, names, unknown_idx: 0 }
    }

    /// Get the library index for a read group hash.
    ///
    /// Returns 0 (unknown) if the RG hash is not found.
    #[must_use]
    pub fn get(&self, rg_hash: u64) -> u16 {
        *self.lookup.get(&rg_hash).unwrap_or(&self.unknown_idx)
    }

    /// Get the library name for an index.
    #[must_use]
    pub fn library_name(&self, idx: u16) -> &Arc<str> {
        self.names.get(idx as usize).unwrap_or(&self.names[0])
    }

    /// Hash a byte slice using `AHash`. Returns 0 for `None`.
    ///
    /// This is the single hashing implementation used by all `hash_*` methods.
    #[must_use]
    pub fn hash_bytes(bytes: Option<&[u8]>) -> u64 {
        use ahash::AHasher;
        use std::hash::{Hash, Hasher};
        match bytes {
            Some(b) => {
                let mut hasher = AHasher::default();
                b.hash(&mut hasher);
                hasher.finish()
            }
            None => 0,
        }
    }

    /// Hash an RG tag value for lookup.
    #[must_use]
    pub fn hash_rg(rg_bytes: &[u8]) -> u64 {
        Self::hash_bytes(Some(rg_bytes))
    }

    /// Hash a cell barcode for `GroupKey`.
    #[must_use]
    pub fn hash_cell_barcode(cell_bytes: Option<&[u8]>) -> u64 {
        Self::hash_bytes(cell_bytes)
    }

    /// Hash a read name for `GroupKey`.
    #[must_use]
    pub fn hash_name(name_bytes: Option<&[u8]>) -> u64 {
        Self::hash_bytes(name_bytes)
    }
}

impl Default for LibraryIndex {
    fn default() -> Self {
        Self {
            lookup: ahash::AHashMap::new(),
            names: vec![Arc::clone(&UNKNOWN_LIBRARY)],
            unknown_idx: 0,
        }
    }
}

/// Information about read positions needed for grouping and ordering.
///
/// This structure stores genomic position information for one or two reads (paired-end),
/// along with library and optional cell barcode information. The R1/R2 positions are
/// stored in their original order without normalization, matching fgbio's behavior.
///
/// # Fields
///
/// * `ref_index1` - Reference sequence index for R1 (first read of pair)
/// * `start1` - Unclipped 5' start position for R1
/// * `strand1` - Strand for R1 (0=forward, 1=reverse)
/// * `ref_index2` - Reference sequence index for R2 (or `UNKNOWN_REF` if unpaired)
/// * `start2` - Unclipped 5' start position for R2
/// * `strand2` - Strand for R2
/// * `library` - Library identifier from RG tag
/// * `cell_barcode` - Optional cell barcode for single-cell applications
///
/// # Ordering
///
/// `ReadInfo` implements `Ord` matching fgbio's case class field order:
/// 1. Reference index 1 (R1)
/// 2. Start position 1 (R1)
/// 3. Strand 1 (R1)
/// 4. Reference index 2 (R2)
/// 5. Start position 2 (R2)
/// 6. Strand 2 (R2)
/// 7. Library
/// 8. Cell barcode (optional)
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ReadInfo {
    /// Reference sequence index for first read
    pub ref_index1: i32,
    /// Unclipped 5' position for first read
    pub start1: i32,
    /// Strand for first read (0=forward, 1=reverse)
    pub strand1: u8,
    /// Reference sequence index for second read
    pub ref_index2: i32,
    /// Unclipped 5' position for second read
    pub start2: i32,
    /// Strand for second read (0=forward, 1=reverse)
    pub strand2: u8,
    /// Library identifier (from RG tag). Uses `Arc<str>` to avoid cloning for every read.
    pub library: Arc<str>,
    /// Optional cell barcode (for single-cell data)
    pub cell_barcode: Option<String>,
}

impl ReadInfo {
    // Maximum values for unmapped reads
    const UNKNOWN_REF: i32 = i32::MAX;
    const UNKNOWN_POS: i32 = i32::MAX;
    const UNKNOWN_STRAND: u8 = u8::MAX;

    /// Creates a new `ReadInfo`, automatically ordering by lower coordinate first.
    ///
    /// This constructor automatically normalizes the order so that the read with the
    /// genomically earlier position becomes position 1. This matches fgbio's behavior
    /// where the lower of the two mates' positions comes first.
    ///
    /// # Arguments
    ///
    /// * `ref_index1` - Reference index for first read
    /// * `start1` - Start position for first read
    /// * `strand1` - Strand for first read (0=forward, 1=reverse)
    /// * `ref_index2` - Reference index for second read
    /// * `start2` - Start position for second read
    /// * `strand2` - Strand for second read
    /// * `library` - Library identifier
    /// * `cell_barcode` - Optional cell barcode
    ///
    /// # Returns
    ///
    /// A new `ReadInfo` with reads ordered by genomic position (lower position first)
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ref_index1: i32,
        start1: i32,
        strand1: u8,
        ref_index2: i32,
        start2: i32,
        strand2: u8,
        library: Arc<str>,
        cell_barcode: Option<String>,
    ) -> Self {
        // Normalize order: put lower genomic position first (matching fgbio)
        // Uses tuple comparison for cleaner, potentially faster code
        let r1_earlier = (ref_index1, start1, strand1) <= (ref_index2, start2, strand2);

        if r1_earlier {
            Self { ref_index1, start1, strand1, ref_index2, start2, strand2, library, cell_barcode }
        } else {
            Self {
                ref_index1: ref_index2,
                start1: start2,
                strand1: strand2,
                ref_index2: ref_index1,
                start2: start1,
                strand2: strand1,
                library,
                cell_barcode,
            }
        }
    }

    /// Creates `ReadInfo` for a single-end or fragment read.
    ///
    /// This constructor is used for unpaired reads or when only one read of a pair
    /// is available. The second read fields are set to UNKNOWN values.
    ///
    /// # Arguments
    ///
    /// * `ref_index` - Reference sequence index
    /// * `start` - Unclipped 5' start position
    /// * `strand` - Strand (0=forward, 1=reverse)
    /// * `library` - Library identifier
    /// * `cell_barcode` - Optional cell barcode
    ///
    /// # Returns
    ///
    /// A new `ReadInfo` for a single read
    #[must_use]
    pub fn single(
        ref_index: i32,
        start: i32,
        strand: u8,
        library: Arc<str>,
        cell_barcode: Option<String>,
    ) -> Self {
        Self {
            ref_index1: ref_index,
            start1: start,
            strand1: strand,
            ref_index2: Self::UNKNOWN_REF,
            start2: Self::UNKNOWN_POS,
            strand2: Self::UNKNOWN_STRAND,
            library,
            cell_barcode,
        }
    }

    /// Creates `ReadInfo` for an unmapped read.
    ///
    /// All position and reference fields are set to UNKNOWN values since the
    /// read has no genomic alignment.
    ///
    /// # Arguments
    ///
    /// * `library` - Library identifier
    /// * `cell_barcode` - Optional cell barcode
    ///
    /// # Returns
    ///
    /// A new `ReadInfo` for an unmapped read
    #[must_use]
    pub fn unmapped(library: Arc<str>, cell_barcode: Option<String>) -> Self {
        Self {
            ref_index1: Self::UNKNOWN_REF,
            start1: Self::UNKNOWN_POS,
            strand1: Self::UNKNOWN_STRAND,
            ref_index2: Self::UNKNOWN_REF,
            start2: Self::UNKNOWN_POS,
            strand2: Self::UNKNOWN_STRAND,
            library,
            cell_barcode,
        }
    }

    /// Converts a strand boolean to a byte representation.
    ///
    /// Converts the strand direction to the byte format used by `ReadInfo`:
    /// - `true` (positive/forward strand) -> 0
    /// - `false` (negative/reverse strand) -> 1
    ///
    /// # Arguments
    ///
    /// * `positive` - `true` for forward strand, `false` for reverse
    ///
    /// # Returns
    ///
    /// Byte representation (0 or 1)
    #[must_use]
    pub fn strand_to_byte(positive: bool) -> u8 {
        u8::from(!positive)
    }

    /// Builds a `ReadInfo` from a template.
    ///
    /// Extracts position information from the template's reads, including unclipped
    /// 5' positions, library information from the header's @RG line (via library lookup),
    /// and optional cell barcode.
    ///
    /// # Arguments
    ///
    /// * `template` - The template containing one or two reads
    /// * `cell_tag` - SAM tag to use for extracting cell barcode (e.g., "CB")
    /// * `library_lookup` - Lookup table mapping read group IDs to library names
    ///
    /// # Returns
    ///
    /// A new `ReadInfo` populated from the template
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The template has no records
    /// - The unclipped position could not be computed (e.g., missing alignment data)
    pub fn from(
        template: &Template,
        cell_tag: Tag,
        library_lookup: &LibraryLookup,
    ) -> Result<Self> {
        use crate::sort::bam_fields;
        use fgumi_raw_bam::{RawRecord, RawRecordView, TagValue};

        // Get primary reads (Option<&RawRecord>)
        let r1 = template.r1();
        let r2 = template.r2();
        let record: &RawRecord = r1.or(r2).ok_or_else(|| {
            anyhow::anyhow!(
                "Template '{}' has no records (empty template)",
                String::from_utf8_lossy(&template.name)
            )
        })?;

        // Extract library from RG tag, then look up the library name from the header.
        // This matches fgbio's behavior where grouping uses the library name (LB field),
        // not the read group ID.
        // Uses Arc<str> to avoid cloning strings - Arc::clone is just a reference count increment.
        let library: Arc<str> =
            if let Some(TagValue::String(rg_bytes)) = record.tags().get(SamTag::RG.as_ref()) {
                let rg_id = std::str::from_utf8(rg_bytes).unwrap_or("unknown");
                library_lookup.get(rg_id).cloned().unwrap_or_else(|| Arc::clone(&UNKNOWN_LIBRARY))
            } else {
                Arc::clone(&UNKNOWN_LIBRARY)
            };

        // Extract cell barcode using the raw tag bytes
        let cell_barcode =
            if let Some(TagValue::String(cb_bytes)) = record.tags().get(cell_tag.as_ref()) {
                Some(String::from_utf8_lossy(cb_bytes).into_owned())
            } else {
                None
            };

        // For paired-end, extract both positions using UNCLIPPED positions.
        // Unmapped reads must use UNKNOWN sentinels rather than position 0 so
        // they sort after mapped reads (matching fgbio's case-class ordering).
        let is_unmapped =
            |r: &RawRecord| (RawRecordView::new(r).flags() & bam_fields::flags::UNMAPPED) != 0;
        let raw_info = |r: &RawRecord| -> Result<(i32, i32, u8)> {
            Ok((
                bam_fields::ref_id(r.as_ref()),
                get_unclipped_position_raw(r)?,
                u8::from((RawRecordView::new(r).flags() & bam_fields::flags::REVERSE) != 0),
            ))
        };

        #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
        if let (Some(r1), Some(r2)) = (r1, r2) {
            match (is_unmapped(r1), is_unmapped(r2)) {
                (true, true) => Ok(ReadInfo::unmapped(library, cell_barcode)),
                (true, false) => {
                    let (ref_index, start, strand) = raw_info(r2)?;
                    Ok(ReadInfo::single(ref_index, start, strand, library, cell_barcode))
                }
                (false, true) => {
                    let (ref_index, start, strand) = raw_info(r1)?;
                    Ok(ReadInfo::single(ref_index, start, strand, library, cell_barcode))
                }
                (false, false) => {
                    let (ref_index1, start1, strand1) = raw_info(r1)?;
                    let (ref_index2, start2, strand2) = raw_info(r2)?;
                    Ok(ReadInfo::new(
                        ref_index1,
                        start1,
                        strand1,
                        ref_index2,
                        start2,
                        strand2,
                        library,
                        cell_barcode,
                    ))
                }
            }
        } else {
            // Single-end: use ReadInfo::single() which correctly sets UNKNOWN values for R2
            // This matches fgbio's behavior where missing R2 uses (Int.MaxValue, Int.MaxValue, Byte.MaxValue)
            if is_unmapped(record) {
                return Ok(ReadInfo::unmapped(library, cell_barcode));
            }
            let (ref_index, start, strand) = raw_info(record)?;
            Ok(ReadInfo::single(ref_index, start, strand, library, cell_barcode))
        }
    }
}

impl PartialOrd for ReadInfo {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ReadInfo {
    /// Compares `ReadInfo` in the same order as fgbio's case class field order:
    /// refIndex1, start1, strand1, refIndex2, start2, strand2, library, cellBarcode
    fn cmp(&self, other: &Self) -> Ordering {
        self.ref_index1
            .cmp(&other.ref_index1)
            .then_with(|| self.start1.cmp(&other.start1))
            .then_with(|| self.strand1.cmp(&other.strand1))
            .then_with(|| self.ref_index2.cmp(&other.ref_index2))
            .then_with(|| self.start2.cmp(&other.start2))
            .then_with(|| self.strand2.cmp(&other.strand2))
            .then_with(|| self.library.cmp(&other.library))
            .then_with(|| self.cell_barcode.cmp(&other.cell_barcode))
    }
}

/// Gets the unclipped 5' position of a read.
///
/// For reads with clipped bases (soft or hard), this function calculates the position that would
/// correspond to the 5' end of the read if it were fully aligned (including clipped bases).
///
/// - **Forward strand**: Returns alignment start minus leading clips (soft + hard)
/// - **Reverse strand**: Returns alignment end plus trailing clips (soft + hard)
///
/// This matches htsjdk's `getUnclippedStart` and `getUnclippedEnd` which consider BOTH
/// soft clips and hard clips.
///
/// # Arguments
///
/// * `record` - The SAM/BAM record
///
/// # Returns
///
/// The unclipped 5' position (0 for unmapped reads)
///
/// Returns the unclipped 5' position (1-based) from raw BAM bytes.
///
/// Returns 0 for unmapped reads. Returns an error if the record is mapped but
/// has no CIGAR (returns `i32::MAX` from the raw helper, which we treat as an error).
///
/// # Errors
///
/// Returns an error if the record is mapped but missing required alignment information
fn get_unclipped_position_raw(record: &fgumi_raw_bam::RawRecord) -> Result<i32> {
    let pos = fgumi_raw_bam::unclipped_5prime_from_raw_bam(record.as_ref());
    if pos == i32::MAX {
        let read_name = String::from_utf8_lossy(fgumi_raw_bam::read_name(record.as_ref()));
        Err(anyhow::anyhow!("Mapped read '{read_name}' missing alignment start position"))
    } else {
        Ok(pos)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_info_ordering() {
        let info1 = ReadInfo::new(0, 100, 0, 0, 200, 1, Arc::from("lib1"), None);
        let info2 = ReadInfo::new(0, 100, 0, 0, 200, 1, Arc::from("lib1"), None);
        assert_eq!(info1, info2);
    }

    #[test]
    fn test_read_info_normalizes_order() {
        // R1 at higher position, R2 at lower position
        let info1 = ReadInfo::new(0, 200, 1, 0, 100, 0, Arc::from("lib1"), None);
        // R1 at lower position, R2 at higher position
        let info2 = ReadInfo::new(0, 100, 0, 0, 200, 1, Arc::from("lib1"), None);

        // Normalization puts lower position first (matching fgbio)
        // info1 started with (200, 1) then (100, 0), should be normalized to (100, 0) then (200, 1)
        assert_eq!(info1.ref_index1, 0);
        assert_eq!(info1.start1, 100);
        assert_eq!(info1.strand1, 0);
        assert_eq!(info1.ref_index2, 0);
        assert_eq!(info1.start2, 200);
        assert_eq!(info1.strand2, 1);

        // Both should be equal after normalization (matching fgbio behavior)
        assert_eq!(info1, info2);
    }

    #[test]
    fn test_read_info_unmapped() {
        let info = ReadInfo::unmapped(Arc::from("lib1"), None);
        assert_eq!(info.ref_index1, ReadInfo::UNKNOWN_REF);
        assert_eq!(info.start1, ReadInfo::UNKNOWN_POS);
    }

    #[test]
    fn test_read_info_single() {
        let info = ReadInfo::single(1, 500, 0, Arc::from("lib1"), None);
        assert_eq!(info.ref_index1, 1);
        assert_eq!(info.start1, 500);
        assert_eq!(info.strand1, 0);
        assert_eq!(info.ref_index2, ReadInfo::UNKNOWN_REF);
    }

    #[test]
    fn test_strand_to_byte() {
        assert_eq!(ReadInfo::strand_to_byte(true), 0);
        assert_eq!(ReadInfo::strand_to_byte(false), 1);
    }

    #[test]
    fn test_read_info_comparison() {
        let info1 = ReadInfo::new(0, 100, 0, 0, 200, 1, Arc::from("lib1"), None);
        let info2 = ReadInfo::new(0, 150, 0, 0, 200, 1, Arc::from("lib1"), None);

        assert!(info1 < info2);
    }

    #[test]
    fn test_different_chromosomes() {
        // R1 on chr0, R2 on chr1 - chr0 is earlier, so R1 stays first
        let info1 = ReadInfo::new(0, 100, 0, 1, 200, 1, Arc::from("lib1"), None);
        // R1 on chr1, R2 on chr0 - chr0 is earlier, so R2 becomes first after normalization
        let info2 = ReadInfo::new(1, 100, 0, 0, 200, 1, Arc::from("lib1"), None);

        // Normalization: lower ref_index comes first
        // info1: chr0 < chr1, so keeps original (ref_index1=0)
        assert_eq!(info1.ref_index1, 0);
        assert_eq!(info1.start1, 100);

        // info2: chr1 > chr0, so swaps to put chr0 first (ref_index1=0)
        assert_eq!(info2.ref_index1, 0);
        assert_eq!(info2.start1, 200); // This was R2's position

        // These are different ReadInfo values (different positions within each slot)
        assert_ne!(info1, info2);
    }

    #[test]
    fn test_build_library_lookup() {
        use bstr::BString;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReadGroup;
        use noodles::sam::header::record::value::map::read_group::tag as rg_tag;

        // Build a header with two read groups having different libraries
        let mut header = Header::builder();

        // First read group with library "libA"
        let rg1 = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("libA"))
            .build()
            .expect("building read group with library should succeed");
        header = header.add_read_group(BString::from("RG1"), rg1);

        // Second read group with library "libB"
        let rg2 = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("libB"))
            .build()
            .expect("building read group with library should succeed");
        header = header.add_read_group(BString::from("RG2"), rg2);

        // Third read group with no library (should default to "unknown")
        let rg3 = Map::<ReadGroup>::builder().build().expect("build should succeed");
        header = header.add_read_group(BString::from("RG3"), rg3);

        let header = header.build();
        let lookup = super::build_library_lookup(&header);

        // Verify the lookup - compare Arc<str> values
        assert_eq!(lookup.len(), 3);
        assert_eq!(lookup.get("RG1").map(AsRef::as_ref), Some("libA"));
        assert_eq!(lookup.get("RG2").map(AsRef::as_ref), Some("libB"));
        assert_eq!(lookup.get("RG3").map(AsRef::as_ref), Some("unknown"));
    }

    #[test]
    fn test_single_uses_unknown_values() {
        // Single-end reads should use UNKNOWN values for R2 (matching fgbio)
        let info = ReadInfo::single(1, 500, 0, Arc::from("lib1"), None);
        assert_eq!(info.ref_index1, 1);
        assert_eq!(info.start1, 500);
        assert_eq!(info.strand1, 0);
        // R2 should have UNKNOWN values (i32::MAX, i32::MAX, u8::MAX)
        assert_eq!(info.ref_index2, i32::MAX);
        assert_eq!(info.start2, i32::MAX);
        assert_eq!(info.strand2, u8::MAX);
    }

    #[test]
    fn test_single_sorts_after_paired() {
        // Single-end reads should sort after paired-end reads because they have UNKNOWN (MAX) values
        let paired = ReadInfo::new(0, 100, 0, 0, 200, 1, Arc::from("lib1"), None);
        let single = ReadInfo::single(0, 100, 0, Arc::from("lib1"), None);

        // Single should be greater than paired (sorts after)
        assert!(single > paired);
    }

    #[test]
    fn test_ord_includes_cell_barcode() {
        // Two ReadInfos that differ only in cell_barcode should compare differently
        let lib: Arc<str> = Arc::from("lib1");
        let info1 = ReadInfo::new(0, 100, 0, 0, 200, 1, lib.clone(), Some("cell1".to_string()));
        let info2 = ReadInfo::new(0, 100, 0, 0, 200, 1, lib.clone(), Some("cell2".to_string()));
        let info3 = ReadInfo::new(0, 100, 0, 0, 200, 1, lib, None);

        // cell1 < cell2 (alphabetically)
        assert!(info1 < info2);
        // None < Some(_) in Rust's Option ordering
        assert!(info3 < info1);
        // All three are different
        assert_ne!(info1, info2);
        assert_ne!(info1, info3);
    }

    #[test]
    fn test_normalization_with_strand_tiebreaker() {
        // Same ref and position, but different strands - should use strand as tiebreaker
        // Forward (0) is "earlier" than reverse (1)
        let info1 = ReadInfo::new(0, 100, 1, 0, 100, 0, Arc::from("lib1"), None);
        // After normalization, strand 0 should come first
        assert_eq!(info1.strand1, 0);
        assert_eq!(info1.strand2, 1);
    }

    mod get_unclipped_position_tests {
        use super::*;
        use fgumi_raw_bam::{RawRecord, SamBuilder as RawSamBuilder, flags as raw_flags};

        /// Parse a CIGAR string to raw BAM u32 ops.
        ///
        /// Encoding: `(len << 4) | op_type` where `op_type` is:
        /// 0=M,1=I,2=D,3=N,4=S,5=H,6=P,7==,8=X
        fn parse_cigar_to_raw(cigar: &str) -> Vec<u32> {
            let mut ops = Vec::new();
            let mut num_str = String::new();
            for c in cigar.chars() {
                if c.is_ascii_digit() {
                    num_str.push(c);
                } else {
                    let len: u32 = num_str.parse().expect("Invalid CIGAR: expected number");
                    let op_type: u32 = match c {
                        'M' => 0,
                        'I' => 1,
                        'D' => 2,
                        'N' => 3,
                        'S' => 4,
                        'H' => 5,
                        'P' => 6,
                        '=' => 7,
                        'X' => 8,
                        _ => panic!("Unknown CIGAR operation: {c}"),
                    };
                    ops.push((len << 4) | op_type);
                    num_str.clear();
                }
            }
            ops
        }

        /// Compute the read sequence length from a CIGAR string
        /// (sum of M, I, S, =, X ops — operations that consume the query).
        fn cigar_query_len(cigar: &str) -> usize {
            let mut num_str = String::new();
            let mut total = 0usize;
            for c in cigar.chars() {
                if c.is_ascii_digit() {
                    num_str.push(c);
                } else {
                    let len: usize = num_str.parse().expect("number");
                    if matches!(c, 'M' | 'I' | 'S' | '=' | 'X') {
                        total += len;
                    }
                    num_str.clear();
                }
            }
            total
        }

        /// Helper to build a test record with specific CIGAR and strand
        fn build_record(alignment_start: usize, cigar: &str, reverse: bool) -> RawRecord {
            let raw_cigar = parse_cigar_to_raw(cigar);
            let qlen = cigar_query_len(cigar).max(1);
            let flags: u16 = if reverse { raw_flags::REVERSE } else { 0 };
            let mut b = RawSamBuilder::new();
            b.sequence(&vec![b'A'; qlen])
                .qualities(&vec![30u8; qlen])
                .flags(flags)
                .ref_id(0)
                .pos(i32::try_from(alignment_start).expect("alignment_start fits i32") - 1)
                .cigar_ops(&raw_cigar);
            b.build()
        }

        #[test]
        fn test_forward_strand_soft_clip_only() {
            // Forward strand with 5S at start: alignment_start=100, CIGAR=5S50M
            // unclipped_start = 100 - 5 = 95
            let record = build_record(100, "5S50M", false);
            let pos = get_unclipped_position_raw(&record)
                .expect("get_unclipped_position_raw should succeed");
            assert_eq!(pos, 95);
        }

        #[test]
        fn test_forward_strand_hard_clip_only() {
            // Forward strand with 10H at start: alignment_start=100, CIGAR=10H50M
            // unclipped_start = 100 - 10 = 90
            let record = build_record(100, "10H50M", false);
            let pos = get_unclipped_position_raw(&record)
                .expect("get_unclipped_position_raw should succeed");
            assert_eq!(pos, 90);
        }

        #[test]
        fn test_forward_strand_hard_and_soft_clip() {
            // Forward strand with 10H5S at start: alignment_start=100, CIGAR=10H5S50M
            // unclipped_start = 100 - 10 - 5 = 85
            let record = build_record(100, "10H5S50M", false);
            let pos = get_unclipped_position_raw(&record)
                .expect("get_unclipped_position_raw should succeed");
            assert_eq!(pos, 85);
        }

        #[test]
        fn test_reverse_strand_soft_clip_only() {
            // Reverse strand with 5S at end: alignment_start=100, CIGAR=50M5S
            // alignment_end = 100 + 50 - 1 = 149
            // unclipped_end = 149 + 5 = 154
            let record = build_record(100, "50M5S", true);
            let pos = get_unclipped_position_raw(&record)
                .expect("get_unclipped_position_raw should succeed");
            assert_eq!(pos, 154);
        }

        #[test]
        fn test_reverse_strand_hard_clip_only() {
            // Reverse strand with 10H at end: alignment_start=100, CIGAR=50M10H
            // alignment_end = 100 + 50 - 1 = 149
            // unclipped_end = 149 + 10 = 159
            let record = build_record(100, "50M10H", true);
            let pos = get_unclipped_position_raw(&record)
                .expect("get_unclipped_position_raw should succeed");
            assert_eq!(pos, 159);
        }

        #[test]
        fn test_reverse_strand_hard_and_soft_clip() {
            // Reverse strand with 5S10H at end: alignment_start=100, CIGAR=50M5S10H
            // alignment_end = 100 + 50 - 1 = 149
            // unclipped_end = 149 + 5 + 10 = 164
            let record = build_record(100, "50M5S10H", true);
            let pos = get_unclipped_position_raw(&record)
                .expect("get_unclipped_position_raw should succeed");
            assert_eq!(pos, 164);
        }

        #[test]
        fn test_no_clipping() {
            // No clipping: alignment_start=100, CIGAR=50M
            // Forward: unclipped_start = 100
            let forward_record = build_record(100, "50M", false);
            assert_eq!(
                get_unclipped_position_raw(&forward_record)
                    .expect("get_unclipped_position_raw should succeed"),
                100
            );

            // Reverse: unclipped_end = 100 + 50 - 1 = 149
            let reverse_record = build_record(100, "50M", true);
            assert_eq!(
                get_unclipped_position_raw(&reverse_record)
                    .expect("get_unclipped_position_raw should succeed"),
                149
            );
        }

        #[test]
        fn test_unmapped_read() {
            let mut b = RawSamBuilder::new();
            b.sequence(b"ACGT").qualities(&[30; 4]).flags(raw_flags::UNMAPPED);
            let record = b.build();
            assert_eq!(
                get_unclipped_position_raw(&record)
                    .expect("get_unclipped_position_raw should succeed"),
                0
            );
        }

        #[test]
        fn test_complex_cigar_with_insertions_deletions() {
            // Complex CIGAR: 5H3S10M2I5M3D10M4S6H
            // Forward strand: alignment_start=100
            // Leading clips = 5H + 3S = 8
            // unclipped_start = 100 - 8 = 92
            let forward_record = build_record(100, "5H3S10M2I5M3D10M4S6H", false);
            assert_eq!(
                get_unclipped_position_raw(&forward_record)
                    .expect("get_unclipped_position_raw should succeed"),
                92
            );

            // Reverse strand: same CIGAR
            // alignment_span = 10 + 5 + 3 + 10 = 28 (M + D consume reference)
            // alignment_end = 100 + 28 - 1 = 127
            // Trailing clips = 4S + 6H = 10
            // unclipped_end = 127 + 10 = 137
            let reverse_record = build_record(100, "5H3S10M2I5M3D10M4S6H", true);
            assert_eq!(
                get_unclipped_position_raw(&reverse_record)
                    .expect("get_unclipped_position_raw should succeed"),
                137
            );
        }
    }
}
