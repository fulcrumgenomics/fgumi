//! Shared grouping and decoded-record domain types.
//!
//! Used by both the group command path (`grouper`, `mi_group`,
//! `commands::{group,dedup}` in the main crate) and the typed-step pipeline
//! (`pipeline::steps::group`, `steps::parse::decode`).
//!
//! These are BAM-record domain types — pre-computed grouping keys, the
//! decoded-record representation, the batching-weight and grouper traits.
//! They live here in `fgumi-bam-io`, next to [`crate::MemoryEstimate`] and
//! the [`fgumi_raw_bam`] raw-record helpers they operate on, rather than in
//! the pipeline crate that merely consumes them.

use std::io;
use std::sync::Arc;

use noodles::sam::alignment::record::data::field::Tag;

use crate::library::LibraryIndex;
use fgumi_raw_bam::{RawRecord, RawRecordView};

pub use crate::mem_estimate::MemoryEstimate;

// ============================================================================
// GroupKey - Pre-computed grouping key for fast comparison in Group step
// ============================================================================

/// Pre-computed grouping key for fast comparison in Group step.
///
/// All fields are integers/hashes for O(1) comparison. This is computed during
/// the parallel Decode step so the serial Group step only does integer comparisons.
///
/// For paired-end reads, positions are normalized so the lower position comes first.
/// For single-end reads, the mate fields use `UNKNOWN_*` sentinel values.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GroupKey {
    // Position info (normalized: lower position first)
    /// Reference sequence index for position 1 (lower).
    pub ref_id1: i32,
    /// Unclipped 5' position for position 1.
    pub pos1: i32,
    /// Strand for position 1 (0=forward, 1=reverse).
    pub strand1: u8,
    /// Reference sequence index for position 2 (higher or mate).
    pub ref_id2: i32,
    /// Unclipped 5' position for position 2.
    pub pos2: i32,
    /// Strand for position 2.
    pub strand2: u8,

    // Grouping metadata
    /// Library index (pre-computed from RG tag via header lookup).
    pub library_idx: u16,
    /// Hash of cell barcode (0 if none).
    pub cell_hash: u64,

    // For name-based grouping within position groups
    /// Hash of QNAME for fast name comparison.
    pub name_hash: u64,
}

impl GroupKey {
    /// Sentinel value for unknown reference ID (unpaired reads).
    pub const UNKNOWN_REF: i32 = i32::MAX;
    /// Sentinel value for unknown position (unpaired reads).
    pub const UNKNOWN_POS: i32 = i32::MAX;
    /// Sentinel value for unknown strand (unpaired reads).
    pub const UNKNOWN_STRAND: u8 = u8::MAX;

    /// Create a `GroupKey` for a paired-end read with mate info.
    ///
    /// Positions are automatically normalized so the lower position comes first.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn paired(
        ref_id: i32,
        pos: i32,
        strand: u8,
        mate_ref_id: i32,
        mate_pos: i32,
        mate_strand: u8,
        library_idx: u16,
        cell_hash: u64,
        name_hash: u64,
    ) -> Self {
        // Normalize: put lower position first (matching ReadInfo behavior)
        let (ref_id1, pos1, strand1, ref_id2, pos2, strand2) =
            if (ref_id, pos, strand) <= (mate_ref_id, mate_pos, mate_strand) {
                (ref_id, pos, strand, mate_ref_id, mate_pos, mate_strand)
            } else {
                (mate_ref_id, mate_pos, mate_strand, ref_id, pos, strand)
            };

        Self { ref_id1, pos1, strand1, ref_id2, pos2, strand2, library_idx, cell_hash, name_hash }
    }

    /// Create a `GroupKey` for a single-end/unpaired read.
    #[must_use]
    pub fn single(
        ref_id: i32,
        pos: i32,
        strand: u8,
        library_idx: u16,
        cell_hash: u64,
        name_hash: u64,
    ) -> Self {
        Self {
            ref_id1: ref_id,
            pos1: pos,
            strand1: strand,
            ref_id2: Self::UNKNOWN_REF,
            pos2: Self::UNKNOWN_POS,
            strand2: Self::UNKNOWN_STRAND,
            library_idx,
            cell_hash,
            name_hash,
        }
    }

    /// Returns the position-only key for grouping by genomic position.
    ///
    /// This is used by `RecordPositionGrouper` to determine if records belong to
    /// the same position group (ignoring name).
    #[must_use]
    pub fn position_key(&self) -> (i32, i32, u8, i32, i32, u8, u16, u64) {
        (
            self.ref_id1,
            self.pos1,
            self.strand1,
            self.ref_id2,
            self.pos2,
            self.strand2,
            self.library_idx,
            self.cell_hash,
        )
    }
}

impl PartialOrd for GroupKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for GroupKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.position_key()
            .cmp(&other.position_key())
            .then_with(|| self.name_hash.cmp(&other.name_hash))
    }
}

impl Default for GroupKey {
    fn default() -> Self {
        Self {
            ref_id1: Self::UNKNOWN_REF,
            pos1: Self::UNKNOWN_POS,
            strand1: Self::UNKNOWN_STRAND,
            ref_id2: Self::UNKNOWN_REF,
            pos2: Self::UNKNOWN_POS,
            strand2: Self::UNKNOWN_STRAND,
            library_idx: 0,
            cell_hash: 0,
            name_hash: 0,
        }
    }
}

// ============================================================================
// DecodedRecord - Record with pre-computed grouping key
// ============================================================================

/// A decoded BAM record with its pre-computed grouping key.
///
/// This is the output of the Decode step and input to the Group step.
/// The key is computed during the parallel Decode step so that the
/// serial Group step only needs to do fast integer comparisons.
#[derive(Debug)]
pub struct DecodedRecord {
    /// Pre-computed grouping key.
    pub key: GroupKey,
    /// Raw BAM record bytes.
    pub(crate) data: RawRecord,
    /// Cached record-relative offset of the UMI tag's value bytes (i.e. the
    /// first byte after the 2-byte tag header and the 1-byte type byte). The
    /// slice `data[umi_value_offset..umi_value_offset + umi_value_len]` yields
    /// the UMI bytes without the trailing NUL.
    ///
    /// Set to [`Self::UMI_OFFSET_UNCACHED`] when no UMI position was cached
    /// during decode (UMI tag missing, not Z-typed, or caching disabled).
    pub(crate) umi_value_offset: u32,
    /// Cached UMI value length in bytes, paired with `umi_value_offset`.
    pub(crate) umi_value_len: u16,
}

impl DecodedRecord {
    /// Sentinel value for `umi_value_offset` indicating no cached UMI position.
    /// Chosen as `u32::MAX` so it can never collide with a real BAM record
    /// offset (BAM records are bounded well under 4 GiB).
    pub const UMI_OFFSET_UNCACHED: u32 = u32::MAX;

    /// Create a decoded record from raw bytes, skipping noodles decode.
    ///
    /// Accepts anything that converts `Into<RawRecord>` (e.g. a bare `Vec<u8>` or
    /// an already-constructed `RawRecord`).
    #[must_use]
    pub fn from_raw_bytes(raw: impl Into<RawRecord>, key: GroupKey) -> Self {
        Self {
            key,
            data: raw.into(),
            umi_value_offset: Self::UMI_OFFSET_UNCACHED,
            umi_value_len: 0,
        }
    }

    /// Attach a cached UMI value position to this decoded record.
    ///
    /// `umi_value_offset` is the record-relative offset of the first UMI value
    /// byte (after the tag header), and `umi_value_len` is the value length
    /// (excluding any trailing NUL).
    pub fn set_cached_umi(&mut self, umi_value_offset: u32, umi_value_len: u16) {
        self.umi_value_offset = umi_value_offset;
        self.umi_value_len = umi_value_len;
    }

    /// Returns the cached UMI bytes (without trailing NUL) if a position was
    /// recorded during decode and still falls within the raw record bytes.
    /// Returns `None` if no cache is set or the cached position is out of range.
    #[must_use]
    pub fn cached_umi(&self) -> Option<&[u8]> {
        if self.umi_value_offset == Self::UMI_OFFSET_UNCACHED {
            return None;
        }
        let start = self.umi_value_offset as usize;
        let end = start.checked_add(self.umi_value_len as usize)?;
        self.data.as_ref().get(start..end)
    }

    /// Returns the cached UMI offset ([`Self::UMI_OFFSET_UNCACHED`] when absent)
    /// and length.
    #[must_use]
    pub fn cached_umi_position(&self) -> (u32, u16) {
        (self.umi_value_offset, self.umi_value_len)
    }

    /// Returns the cached UMI `(offset, len)` position, or `None` when no
    /// position was recorded during decode (the [`Self::UMI_OFFSET_UNCACHED`]
    /// sentinel). Keeps the sentinel check at the owning layer so callers can
    /// branch on a plain `Option` instead of comparing against the sentinel.
    #[must_use]
    pub fn cached_umi_position_opt(&self) -> Option<(u32, u16)> {
        if self.umi_value_offset == Self::UMI_OFFSET_UNCACHED {
            None
        } else {
            Some((self.umi_value_offset, self.umi_value_len))
        }
    }

    /// Returns a reference to the raw bytes.
    #[must_use]
    pub fn raw_bytes(&self) -> &[u8] {
        self.data.as_ref()
    }

    /// Mutable access to the underlying [`RawRecord`]. Used by mid-chain
    /// `Parallel` Process steps that need to mutate record bytes (e.g.,
    /// MQ bumping, tag rewriting) without rebuilding the `DecodedRecord`
    /// or invalidating its pre-computed `GroupKey`. Caller must not
    /// change the record's identity (qname, library bytes, cell barcode)
    /// — those feed `key`, which we don't recompute here.
    pub fn raw_bytes_mut(&mut self) -> &mut RawRecord {
        &mut self.data
    }

    /// Takes the [`RawRecord`] out.
    #[must_use]
    pub fn into_raw_bytes(self) -> RawRecord {
        self.data
    }
}

impl MemoryEstimate for DecodedRecord {
    fn estimate_heap_size(&self) -> usize {
        // RawRecord::capacity() returns the inner Vec<u8> capacity.
        self.data.capacity()
    }
}
// Vec<DecodedRecord>, Vec<RecordBuf>, Vec<u8>, RecordBuf, () — all provided
// by the blanket/foreign impls in fgumi_bam_io::mem_estimate.

// ============================================================================
// GroupKeyConfig - Configuration for computing `GroupKey` during Decode
// ============================================================================

/// Configuration for computing `GroupKey` during the Decode step.
///
/// When this is provided to the pipeline, the Decode step will compute
/// full `GroupKey` values for each record. This moves expensive computations
/// (CIGAR parsing, tag extraction) from the serial Group step to the parallel
/// Decode step.
#[derive(Debug, Clone)]
pub struct GroupKeyConfig {
    /// Library index for fast RG → library lookup.
    pub library_index: Arc<LibraryIndex>,
    /// Tag used for cell barcode extraction. None skips cell extraction.
    pub cell_tag: Option<Tag>,
    /// When `true`, the Decode step computes only the read-name hash and
    /// leaves the rest of the [`GroupKey`] at its default. Use for stages
    /// that group by queryname (e.g. `correct`) and read only
    /// [`GroupKey::name_hash`] — it skips the CIGAR 5′-position walk and the
    /// aux-tag (RG/CB/MC) extraction pass entirely. See [`name_hash_key`].
    pub name_hash_only: bool,
    /// UMI tag (raw 2-byte form) whose value position should be cached on each
    /// [`DecodedRecord`] during decode. `None` disables caching — downstream
    /// code must fall back to scanning aux data. This is orthogonal to
    /// `name_hash_only`: the UMI cache scan is gated solely on `umi_tag` being
    /// set (in practice only the Group stage sets it).
    pub umi_tag: Option<[u8; 2]>,
}

impl GroupKeyConfig {
    /// Create a new `GroupKeyConfig` that computes the full position/cell key.
    #[must_use]
    pub fn new(library_index: LibraryIndex, cell_tag: Tag) -> Self {
        Self {
            library_index: Arc::new(library_index),
            cell_tag: Some(cell_tag),
            name_hash_only: false,
            umi_tag: None,
        }
    }

    /// Create a `GroupKeyConfig` without cell barcode extraction.
    #[must_use]
    pub fn new_raw_no_cell(library_index: LibraryIndex) -> Self {
        Self {
            library_index: Arc::new(library_index),
            cell_tag: None,
            name_hash_only: false,
            umi_tag: None,
        }
    }

    /// Create a `GroupKeyConfig` that computes only the read-name hash.
    ///
    /// For queryname-grouping stages (e.g. `correct`) whose grouper reads only
    /// [`GroupKey::name_hash`]; skips the CIGAR position walk and the aux-tag
    /// extraction pass. `library_index` is retained (unused for the key) so
    /// the config shape is uniform across the pipeline.
    #[must_use]
    pub fn name_hash_only(library_index: LibraryIndex) -> Self {
        Self {
            library_index: Arc::new(library_index),
            cell_tag: None,
            name_hash_only: true,
            umi_tag: None,
        }
    }

    /// Enable UMI position caching for the given raw 2-byte UMI tag (e.g. `*b"RX"`).
    ///
    /// The Decode step will record the UMI value's record-relative position on
    /// each [`DecodedRecord`] so downstream UMI lookups can slice it directly
    /// via [`DecodedRecord::cached_umi`] without re-scanning aux data.
    #[must_use]
    pub fn with_umi_tag(mut self, umi_tag: [u8; 2]) -> Self {
        self.umi_tag = Some(umi_tag);
        self
    }
}

impl Default for GroupKeyConfig {
    fn default() -> Self {
        Self {
            library_index: Arc::new(LibraryIndex::default()),
            cell_tag: Some(Tag::from([b'C', b'B'])), // Default cell barcode tag (CB)
            name_hash_only: false,
            umi_tag: None,
        }
    }
}

/// Compute a [`GroupKey`] containing only the read-name hash, leaving all
/// position/strand/library/cell fields at their default.
///
/// Reproduces exactly the `name_hash` that [`compute_group_key_from_raw`]
/// computes (empty name → `hash_name(None)`), so a queryname grouper sees an
/// identical hash. Skips the CIGAR 5′-position walk and the RG/CB/MC aux-tag
/// extraction pass.
///
/// # Panics
///
/// Panics if `raw` is not a validated BAM record payload (same contract as
/// [`compute_group_key_from_raw`]).
#[must_use]
pub fn name_hash_key(raw: &[u8]) -> GroupKey {
    let name = fgumi_raw_bam::read_name(raw);
    let name_hash = if name.is_empty() {
        LibraryIndex::hash_name(None)
    } else {
        LibraryIndex::hash_name(Some(name))
    };
    GroupKey { name_hash, ..GroupKey::default() }
}

/// Compute a `GroupKey` directly from raw BAM bytes, matching `compute_group_key()` exactly.
///
/// Uses 1-based coordinate helpers to produce identical keys to the noodles path.
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
) -> GroupKey {
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
        return GroupKey { name_hash, ..GroupKey::default() };
    }

    // Own position (1-based, matching noodles) — zero-allocation CIGAR iteration
    let reverse = (flg & fgumi_raw_bam::flags::REVERSE) != 0;
    let own_pos = fgumi_raw_bam::unclipped_5prime_from_raw_bam(raw);

    let own_ref_id = fgumi_raw_bam::ref_id(raw);
    let strand = u8::from(reverse);

    // Single-pass aux tag extraction (RG, cell barcode, MC)
    let aux_data = fgumi_raw_bam::aux_data_slice(raw);
    let cell_tag_bytes = cell_tag.map_or([0u8; 2], |t| [t.as_ref()[0], t.as_ref()[1]]);
    let aux_tags = fgumi_raw_bam::extract_aux_string_tags(aux_data, cell_tag_bytes, None);

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
        return GroupKey::single(own_ref_id, own_pos, strand, library_idx, cell_hash, name_hash);
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

    match mate_pos_result {
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
            // No MC tag — fall back to single-end behavior
            GroupKey::single(own_ref_id, own_pos, strand, library_idx, cell_hash, name_hash)
        }
    }
}

/// Groups a stream of in-order [`DecodedRecord`]s into completed groups.
///
/// Implementors maintain partial groups across `add_records` calls and emit
/// completed ones; `finish` flushes any trailing partial group at EOF. Used
/// by the pipeline's Group step and the standalone grouping commands.
pub trait Grouper: Send {
    /// The type of group produced by this grouper.
    type Group: Send;

    /// Add decoded records to the grouper.
    ///
    /// Records are guaranteed to be in order (from template-coordinate sorted BAM).
    /// The grouper maintains partial groups waiting for more records.
    ///
    /// Each `DecodedRecord` contains the record plus a pre-computed `GroupKey`
    /// for fast comparison (position, name hash, library, etc.).
    ///
    /// Returns completed groups (may be empty if more records are needed).
    ///
    /// # Errors
    ///
    /// Returns an I/O error if grouping logic encounters invalid data.
    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>>;

    /// Signal that no more input will arrive (EOF).
    ///
    /// Returns any remaining partial group.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if finalizing the grouper fails.
    fn finish(&mut self) -> io::Result<Option<Self::Group>>;

    /// Returns true if the grouper has a partial group.
    fn has_pending(&self) -> bool;
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::SamBuilder;
    use fgumi_raw_bam::flags;

    /// CIGAR op `(len << 4) | op_code`; op 0 = `M`.
    fn cigar_m(len: u32) -> u32 {
        len << 4
    }

    // ========================================================================
    // GroupKey::paired normalization
    // ========================================================================

    #[test]
    fn paired_normalizes_swapped_own_and_mate_to_equal_keys() {
        // Two reads of the same template: one sees (own=A, mate=B), the other
        // sees (own=B, mate=A). After normalization both must yield the same
        // GroupKey, so they group together.
        let a = GroupKey::paired(0, 100, 0, 0, 200, 1, 7, 99, 42);
        let b = GroupKey::paired(0, 200, 1, 0, 100, 0, 7, 99, 42);
        assert_eq!(a, b, "swapped own/mate positions must normalize to equal keys");

        // The lower (ref_id, pos, strand) tuple is placed in slot 1.
        assert_eq!((a.ref_id1, a.pos1, a.strand1), (0, 100, 0));
        assert_eq!((a.ref_id2, a.pos2, a.strand2), (0, 200, 1));
    }

    #[test]
    fn paired_keys_order_by_normalized_position_then_name() {
        // Position 1 < position 2 orders the keys; equal positions fall back to
        // name_hash ordering.
        let lower = GroupKey::paired(0, 100, 0, 0, 200, 0, 0, 0, 1);
        let higher = GroupKey::paired(0, 150, 0, 0, 200, 0, 0, 0, 1);
        assert!(lower < higher, "lower position-1 must sort first");

        let same_pos_low_name = GroupKey::paired(0, 100, 0, 0, 200, 0, 0, 0, 1);
        let same_pos_high_name = GroupKey::paired(0, 100, 0, 0, 200, 0, 0, 0, 2);
        assert!(same_pos_low_name < same_pos_high_name, "equal positions must order by name_hash",);
    }

    // ========================================================================
    // compute_group_key_from_raw: secondary / supplementary fallback
    // ========================================================================

    #[test]
    fn secondary_and_supplementary_records_yield_name_hash_only_key() {
        let lib = LibraryIndex::default();
        let expected_name_hash = LibraryIndex::hash_name(Some(b"rec1"));

        for flag in [flags::SECONDARY, flags::SUPPLEMENTARY] {
            // Give the record a real mapped position so we can prove the
            // position fields are *not* derived for secondary/supplementary.
            let mut b = SamBuilder::new();
            b.ref_id(0)
                .pos(500)
                .flags(flag)
                .read_name(b"rec1")
                .cigar_ops(&[cigar_m(50)])
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30]);
            let rec = b.build();

            let key = compute_group_key_from_raw(rec.as_ref(), &lib, None);

            // Only the name hash is set; every position field stays at default.
            assert_eq!(key, GroupKey { name_hash: expected_name_hash, ..GroupKey::default() });
            assert_eq!(key.ref_id1, GroupKey::UNKNOWN_REF);
            assert_eq!(key.pos1, GroupKey::UNKNOWN_POS);
            assert_eq!(key.strand1, GroupKey::UNKNOWN_STRAND);
        }
    }

    // ========================================================================
    // compute_group_key_from_raw: paired record missing MC falls back to single
    // ========================================================================

    #[test]
    fn paired_without_mc_tag_falls_back_to_single_semantics() {
        let lib = LibraryIndex::default();

        // Paired + mate mapped, but no MC tag: cannot compute the mate's
        // unclipped 5' position, so the key must use single-end semantics
        // (mate fields left at the UNKNOWN sentinels).
        let mut b = SamBuilder::new();
        b.ref_id(0)
            .pos(1000)
            .flags(flags::PAIRED)
            .mate_ref_id(0)
            .mate_pos(1200)
            .read_name(b"pair_no_mc")
            .cigar_ops(&[cigar_m(50)])
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30]);
        let rec = b.build();

        let key = compute_group_key_from_raw(rec.as_ref(), &lib, None);

        // Mate fields fall back to the single-end sentinels.
        assert_eq!(key.ref_id2, GroupKey::UNKNOWN_REF);
        assert_eq!(key.pos2, GroupKey::UNKNOWN_POS);
        assert_eq!(key.strand2, GroupKey::UNKNOWN_STRAND);

        // The own position is still populated; the key matches the single-end
        // key built from the same fields.
        let expected = GroupKey::single(
            key.ref_id1,
            key.pos1,
            key.strand1,
            key.library_idx,
            key.cell_hash,
            LibraryIndex::hash_name(Some(b"pair_no_mc")),
        );
        assert_eq!(key, expected);
    }

    // ========================================================================
    // name_hash_key parity with compute_group_key_from_raw
    // ========================================================================

    #[test]
    fn name_hash_key_matches_compute_group_key_name_hash() {
        let lib = LibraryIndex::default();

        // A representative spread of records: mapped single-end, secondary,
        // paired-without-MC, and an empty-named record (exercises the
        // hash_name(None) branch shared by both functions).
        let mut single = SamBuilder::new();
        single
            .ref_id(0)
            .pos(10)
            .read_name(b"alpha")
            .cigar_ops(&[cigar_m(8)])
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8]);

        let mut secondary = SamBuilder::new();
        secondary.ref_id(0).pos(20).flags(flags::SECONDARY).read_name(b"beta");

        let mut paired = SamBuilder::new();
        paired
            .ref_id(0)
            .pos(30)
            .flags(flags::PAIRED)
            .mate_ref_id(0)
            .mate_pos(60)
            .read_name(b"gamma");

        let mut empty_name = SamBuilder::new();
        empty_name.ref_id(0).pos(40).read_name(b"");

        for mut builder in [single, secondary, paired, empty_name] {
            let rec = builder.build();
            let raw = rec.as_ref();
            assert_eq!(
                name_hash_key(raw).name_hash,
                compute_group_key_from_raw(raw, &lib, None).name_hash,
                "name_hash parity mismatch",
            );
        }
    }
}
