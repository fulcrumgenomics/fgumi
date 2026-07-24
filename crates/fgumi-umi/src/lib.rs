#![deny(unsafe_code)]

//! UMI (Unique Molecular Identifier) utilities
//!
//! This crate provides functionality for working with UMIs including:
//! - UMI assignment strategies (identity, edit-distance, adjacency, paired)
//! - Tag set management for consensus calling
//! - UMI grouping and family analysis

pub mod assigner;

// Re-export commonly used items
pub use assigner::{
    AdjacencyUmiAssigner, EDIT_INDEX_THRESHOLD, IdentityUmiAssigner, PairedUmiAssigner,
    SimpleErrorUmiAssigner, Strategy, Umi, UmiAssigner,
};

use std::collections::HashSet;

/// Molecule identifier for UMI grouping.
///
/// Represents the assigned molecule ID during UMI-based grouping. The ID is stored
/// as an integer during processing for efficiency, and only converted to a string
/// when writing the final BAM output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum MoleculeId {
    /// Not yet assigned (default state)
    #[default]
    None,
    /// Non-paired strategy: simple integer ID (e.g., "42")
    Single(u64),
    /// Paired strategy, top strand (e.g., "42/A")
    PairedA(u64),
    /// Paired strategy, bottom strand (e.g., "42/B")
    PairedB(u64),
}

impl MoleculeId {
    /// Get the base molecule ID value, if assigned.
    #[inline]
    #[must_use]
    pub fn id(&self) -> Option<u64> {
        match self {
            MoleculeId::None => None,
            MoleculeId::Single(id) | MoleculeId::PairedA(id) | MoleculeId::PairedB(id) => Some(*id),
        }
    }

    /// Check if a molecule ID has been assigned.
    #[inline]
    #[must_use]
    pub fn is_assigned(&self) -> bool {
        !matches!(self, MoleculeId::None)
    }

    /// Convert to string representation with a base offset applied.
    ///
    /// Used when serializing to BAM output, where local IDs (0, 1, 2, ...)
    /// need to be converted to global IDs by adding the base offset.
    #[must_use]
    pub fn to_string_with_offset(&self, base: u64) -> String {
        let mut buf = String::new();
        self.write_with_offset(base, &mut buf);
        buf
    }

    /// Write string representation into a reusable buffer, returning the bytes.
    ///
    /// Avoids per-call String allocation by reusing the caller's buffer.
    ///
    /// # Panics
    ///
    /// Panics consistently in both debug and release builds if `base + id`
    /// would overflow `u64`. With `id` and `base` both bounded by the number
    /// of molecules in the input (which fits comfortably under `u64::MAX` for
    /// every realistic sequencing run), this is unreachable in practice; the
    /// explicit `checked_add` keeps debug and release behavior identical
    /// rather than wrapping silently in release.
    pub fn write_with_offset<'a>(&self, base: u64, buf: &'a mut String) -> &'a [u8] {
        use std::fmt::Write;
        buf.clear();
        match self {
            MoleculeId::None => {}
            MoleculeId::Single(id) => {
                let total = id.checked_add(base).expect("MoleculeId offset overflow");
                write!(buf, "{total}").expect("write to String is infallible");
            }
            MoleculeId::PairedA(id) => {
                let total = id.checked_add(base).expect("MoleculeId offset overflow");
                write!(buf, "{total}/A").expect("write to String is infallible");
            }
            MoleculeId::PairedB(id) => {
                let total = id.checked_add(base).expect("MoleculeId offset overflow");
                write!(buf, "{total}/B").expect("write to String is infallible");
            }
        }
        buf.as_bytes()
    }

    /// Convert to a Vec index for grouping templates by molecule ID.
    ///
    /// Returns `None` for unassigned `MoleculeId`s. For assigned IDs:
    /// - `Single(id)`: returns id (indices: 0, 1, 2, ...)
    /// - `PairedA(id)`: returns id * 2 (indices: 0, 2, 4, ...)
    /// - `PairedB(id)`: returns id * 2 + 1 (indices: 1, 3, 5, ...)
    ///
    /// This preserves the sort order (A before B for same base id) and allows
    /// efficient Vec-based grouping instead of `HashMap`.
    #[inline]
    #[must_use]
    #[expect(clippy::cast_possible_truncation, reason = "molecule IDs never exceed usize::MAX / 2")]
    pub fn to_vec_index(&self) -> Option<usize> {
        match self {
            MoleculeId::None => None,
            MoleculeId::Single(id) => Some(*id as usize),
            MoleculeId::PairedA(id) => Some(*id as usize * 2),
            MoleculeId::PairedB(id) => Some(*id as usize * 2 + 1),
        }
    }

    /// Check if this is a paired molecule ID (`PairedA` or `PairedB`).
    #[inline]
    #[must_use]
    pub fn is_paired(&self) -> bool {
        matches!(self, MoleculeId::PairedA(_) | MoleculeId::PairedB(_))
    }

    /// Get the base ID without suffix for grouping purposes.
    ///
    /// For paired `MoleculeId`s, this returns the base ID (without /A or /B).
    /// Used when grouping templates that share the same molecule but different strands.
    #[inline]
    #[must_use]
    pub fn base_id_string(&self) -> String {
        match self {
            MoleculeId::None => String::new(),
            MoleculeId::Single(id) | MoleculeId::PairedA(id) | MoleculeId::PairedB(id) => {
                id.to_string()
            }
        }
    }

    /// Return a new `MoleculeId` whose numeric `id` is shifted by `offset`,
    /// preserving the variant tag.
    ///
    /// `MoleculeId::None` is returned unchanged. `Single`, `PairedA`, and
    /// `PairedB` each have `offset` added to their inner `id`. Used by the
    /// MI Assign pipeline stage to convert the per-position-group local IDs
    /// produced by `assigner.assign(...)` into globally unique IDs by
    /// folding in a serial-ordered cumulative offset.
    ///
    /// # Panics
    ///
    /// Panics consistently in both debug and release builds if `id + offset`
    /// would overflow `u64`. The cumulative MI counter is bounded by the
    /// number of molecules in the input, which fits comfortably under
    /// `u64::MAX` for every realistic sequencing run — `checked_add` is used
    /// only to keep debug and release behavior identical rather than
    /// wrapping silently in release.
    #[inline]
    #[must_use]
    pub fn with_offset(self, offset: u64) -> Self {
        match self {
            MoleculeId::None => MoleculeId::None,
            MoleculeId::Single(id) => {
                MoleculeId::Single(id.checked_add(offset).expect("MoleculeId offset overflow"))
            }
            MoleculeId::PairedA(id) => {
                MoleculeId::PairedA(id.checked_add(offset).expect("MoleculeId offset overflow"))
            }
            MoleculeId::PairedB(id) => {
                MoleculeId::PairedB(id.checked_add(offset).expect("MoleculeId offset overflow"))
            }
        }
    }
}

impl std::fmt::Display for MoleculeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MoleculeId::None => write!(f, ""),
            MoleculeId::Single(id) => write!(f, "{id}"),
            MoleculeId::PairedA(id) => write!(f, "{id}/A"),
            MoleculeId::PairedB(id) => write!(f, "{id}/B"),
        }
    }
}

/// Named sets of tags for common use cases
///
/// Provides predefined tag sets that can be referenced by name (e.g., "Consensus")
/// instead of listing individual tags.
pub struct TagSets;

impl TagSets {
    /// Consensus per-base tags that should be **reversed** on negative-strand reads.
    ///
    /// Mirrors fgbio `ConsensusTags.PerBase.TagsToReverse`: the per-base
    /// depth/error/quality arrays — `cd`,`ce` (consensus), `ad`,`ae` (top strand),
    /// `bd`,`be` (bottom strand), and `aq`,`bq` (per-strand quals). These are
    /// position-indexed arrays, so they must be reversed (not reverse-complemented)
    /// when the read maps to the negative strand.
    pub const CONSENSUS_REVERSE: &[&str] = &["cd", "ce", "ad", "ae", "bd", "be", "aq", "bq"];

    /// Consensus per-base **base** tags that should be **reverse-complemented** on
    /// negative-strand reads.
    ///
    /// Mirrors fgbio `ConsensusTags.PerBase.TagsToReverseComplement`: the per-base
    /// consensus base arrays `ac` (top strand) and `bc` (bottom strand). Being
    /// sequence bases, they are reverse-complemented (not merely reversed) on the
    /// negative strand. (Note: `aD`/`bD`/`cD` are per-*read* scalar depths — never
    /// reversed or revcomped.)
    pub const CONSENSUS_REVCOMP: &[&str] = &["ac", "bc"];
}

/// Information about which tags to remove, reverse, or reverse complement
///
/// Stores sets of tag names that should be manipulated during the merge operation.
/// Named tag sets (like "Consensus") are expanded into their constituent tags.
#[derive(Debug, Clone)]
pub struct TagInfo {
    /// Tags to remove from mapped reads
    pub remove: HashSet<String>,
    /// Tags to reverse for negative strand reads
    pub reverse: HashSet<String>,
    /// Tags to reverse complement for negative strand reads
    pub revcomp: HashSet<String>,
}

impl TagInfo {
    /// Creates a new `TagInfo` from lists of tag names
    ///
    /// Expands named tag sets (e.g., "Consensus") into their constituent tags.
    /// Individual tag names are added as-is.
    ///
    /// # Arguments
    ///
    /// * `remove` - Tags to remove from mapped reads
    /// * `reverse` - Tags to reverse for negative strand (or "Consensus")
    /// * `revcomp` - Tags to reverse complement for negative strand (or "Consensus")
    ///
    /// # Examples
    ///
    /// ```
    /// use fgumi_umi::TagInfo;
    /// let tag_info = TagInfo::new(
    ///     vec!["AS".to_string()],
    ///     vec!["Consensus".to_string()],
    ///     vec!["Consensus".to_string()],
    /// );
    /// ```
    #[must_use]
    pub fn new(remove: Vec<String>, reverse: Vec<String>, revcomp: Vec<String>) -> Self {
        let mut reverse_set = HashSet::new();
        let mut revcomp_set = HashSet::new();

        for tag in reverse {
            if tag == "Consensus" {
                reverse_set.extend(TagSets::CONSENSUS_REVERSE.iter().map(|&s| s.to_owned()));
            } else {
                reverse_set.insert(tag);
            }
        }

        for tag in revcomp {
            if tag == "Consensus" {
                revcomp_set.extend(TagSets::CONSENSUS_REVCOMP.iter().map(|&s| s.to_owned()));
            } else {
                revcomp_set.insert(tag);
            }
        }

        TagInfo { remove: remove.into_iter().collect(), reverse: reverse_set, revcomp: revcomp_set }
    }

    /// Checks if any tags need reversal or reverse complementation
    ///
    /// # Returns
    ///
    /// `true` if there are any tags to reverse or revcomp, `false` otherwise
    #[must_use]
    pub fn has_revs_or_revcomps(&self) -> bool {
        !self.reverse.is_empty() || !self.revcomp.is_empty()
    }
}

/// Result of validating a UMI string.
///
/// Used by both `group` and `dedup` commands to consistently filter UMIs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UmiValidation {
    /// UMI is valid with the given number of bases (ACGT, case-insensitive)
    Valid(usize),
    /// UMI contains 'N' (uppercase only, matching fgbio behavior)
    ContainsN,
}

/// Validates a UMI string, matching fgbio's behavior.
///
/// This function:
/// - Counts valid DNA bases (A, C, G, T - case insensitive)
/// - Rejects UMIs containing uppercase 'N' (ambiguous base)
/// - Skips dashes (paired UMI separator) and other characters including lowercase 'n'
///
/// This matches fgbio's `GroupReadsByUmi` which:
/// 1. Only rejects uppercase 'N': `.filter(r => !umi.contains('N'))`
/// 2. Counts bases case-insensitively: `umi.toUpperCase.count(c => isUpperACGTN(c))`
///
/// # Arguments
///
/// * `umi` - The UMI string to validate
///
/// # Returns
///
/// * `UmiValidation::Valid(count)` - UMI is valid with `count` DNA bases
/// * `UmiValidation::ContainsN` - UMI contains uppercase 'N'
///
/// # Examples
///
/// ```
/// use fgumi_umi::{validate_umi, UmiValidation};
///
/// assert_eq!(validate_umi(b"ACGT"), UmiValidation::Valid(4));
/// assert_eq!(validate_umi(b"acgt"), UmiValidation::Valid(4));
/// assert_eq!(validate_umi(b"ACGT-TGCA"), UmiValidation::Valid(8));
/// assert_eq!(validate_umi(b"ACNT"), UmiValidation::ContainsN);
/// assert_eq!(validate_umi(b"acnt"), UmiValidation::Valid(3)); // lowercase 'n' skipped
/// ```
#[must_use]
pub fn validate_umi(umi: &[u8]) -> UmiValidation {
    let mut base_count = 0usize;
    for &b in umi {
        match b {
            b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' => base_count += 1,
            b'N' => return UmiValidation::ContainsN,
            _ => {} // Skip dash, lowercase 'n', other chars (matches fgbio)
        }
    }
    UmiValidation::Valid(base_count)
}

/// Extracts the molecular identifier base by truncating at the last `/`.
///
/// MI tags for paired/duplex molecule ids carry a trailing strand (or other)
/// suffix after a `/` — e.g. `1/A` and `1/B` are the two strands of molecule
/// `1` (see the `Paired` assignment strategy in `group`). To correlate reads
/// back to their source molecule, everything from the last `/` onwards is
/// dropped. This mirrors fgbio `ReviewConsensusVariants.toMi`
/// (`mi.substring(0, mi.lastIndexOf('/'))`), which truncates at the last `/`
/// regardless of what the suffix is — not only `/A`/`/B`.
///
/// A leading `/` (index 0) is preserved, matching fgbio's `slash > 0` guard.
///
/// # Examples
///
/// ```
/// use fgumi_umi::extract_mi_base;
///
/// assert_eq!(extract_mi_base("123"), "123");
/// assert_eq!(extract_mi_base("123/A"), "123");
/// assert_eq!(extract_mi_base("123/B"), "123");
/// assert_eq!(extract_mi_base("abc/456/A"), "abc/456"); // truncates at last `/`
/// assert_eq!(extract_mi_base("123/C"), "123");         // any suffix, not just /A,/B
/// assert_eq!(extract_mi_base("123/456"), "123");
/// assert_eq!(extract_mi_base("/A"), "/A");             // leading `/` preserved
/// ```
#[must_use]
pub fn extract_mi_base(mi: &str) -> &str {
    match mi.rfind('/') {
        Some(slash) if slash > 0 => &mi[..slash],
        _ => mi,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_consensus_reverse_returns_correct_tags() {
        // fgbio ConsensusTags.PerBase.TagsToReverse = {cd, ce, ad, ae, bd, be, aq, bq}
        assert_eq!(TagSets::CONSENSUS_REVERSE.len(), 8);
        assert_eq!(TagSets::CONSENSUS_REVERSE, &["cd", "ce", "ad", "ae", "bd", "be", "aq", "bq"]);
    }

    /// Programmatic fgbio baseline for the per-base consensus tag sets, so the
    /// `TagSets::CONSENSUS_*` pins are checked against fgbio's *definition* rather than a
    /// second hand-copied literal.
    ///
    /// This reconstructs fgbio's `ConsensusTags.PerBase` sets from the individually named tag
    /// symbols, in fgbio's exact `Seq(...)` order (fgbio
    /// `src/main/scala/com/fulcrumgenomics/umi/ConsensusTags.scala`, `object PerBase`,
    /// verified against fgbio main). If fgbio adds, removes, or reorders a per-base tag, update
    /// the named symbols below and this test catches any resulting drift in the production
    /// constants. The intentional divergence from fgumi's *superset* `per_base` lists — which
    /// additionally cover the CODEC conversion tags — is asserted separately by
    /// `test_zipper_consensus_tagset_is_subset_of_canonical_per_base` in `src/lib/tag_reversal.rs`.
    #[test]
    fn test_consensus_tag_sets_match_fgbio_perbase_baseline() {
        // fgbio ConsensusTags.PerBase named symbol -> 2-char SAM tag.
        let raw_read_count = "cd"; // RawReadCount    (consensus depth)
        let raw_read_errors = "ce"; // RawReadErrors   (consensus errors)
        let ab_raw_read_count = "ad"; // AbRawReadCount
        let ba_raw_read_count = "bd"; // BaRawReadCount
        let ab_raw_read_errors = "ae"; // AbRawReadErrors
        let ba_raw_read_errors = "be"; // BaRawReadErrors
        let ab_consensus_bases = "ac"; // AbConsensusBases
        let ba_consensus_bases = "bc"; // BaConsensusBases
        let ab_consensus_quals = "aq"; // AbConsensusQuals
        let ba_consensus_quals = "bq"; // BaConsensusQuals

        // fgbio: TagsToReverse = Seq(RawReadCount, RawReadErrors, AbRawReadCount, AbRawReadErrors,
        //                            BaRawReadCount, BaRawReadErrors, AbConsensusQuals, BaConsensusQuals)
        let fgbio_tags_to_reverse = [
            raw_read_count,
            raw_read_errors,
            ab_raw_read_count,
            ab_raw_read_errors,
            ba_raw_read_count,
            ba_raw_read_errors,
            ab_consensus_quals,
            ba_consensus_quals,
        ];
        // fgbio: TagsToReverseComplement = Seq(AbConsensusBases, BaConsensusBases)
        let fgbio_tags_to_reverse_complement = [ab_consensus_bases, ba_consensus_bases];

        assert_eq!(
            TagSets::CONSENSUS_REVERSE,
            &fgbio_tags_to_reverse,
            "CONSENSUS_REVERSE drifted from fgbio ConsensusTags.PerBase.TagsToReverse",
        );
        assert_eq!(
            TagSets::CONSENSUS_REVCOMP,
            &fgbio_tags_to_reverse_complement,
            "CONSENSUS_REVCOMP drifted from fgbio ConsensusTags.PerBase.TagsToReverseComplement",
        );
    }

    // ========================================================================
    // MoleculeId::with_offset
    //
    // The deterministic-MI design's MI Assign stage rewrites local per-batch
    // ids into global ids by calling `with_offset(cumulative_offset)` on
    // every template's MoleculeId. The contract these tests pin down:
    //   * `None` is preserved (an unassigned slot stays unassigned).
    //   * Each `Some` variant adds the offset to its `id`, leaving the
    //     variant tag (Single / PairedA / PairedB) untouched.
    //   * `with_offset(0)` is the identity function.
    //   * Offsets compose by addition: applying `a` then `b` is identical
    //     to applying `a + b` in one step.
    // ========================================================================

    #[test]
    fn molecule_id_with_offset_preserves_none() {
        assert_eq!(MoleculeId::None.with_offset(0), MoleculeId::None);
        assert_eq!(MoleculeId::None.with_offset(42), MoleculeId::None);
        assert_eq!(MoleculeId::None.with_offset(u64::MAX), MoleculeId::None);
    }

    #[test]
    fn molecule_id_with_offset_shifts_single() {
        assert_eq!(MoleculeId::Single(0).with_offset(7), MoleculeId::Single(7));
        assert_eq!(MoleculeId::Single(5).with_offset(100), MoleculeId::Single(105));
    }

    #[test]
    fn molecule_id_with_offset_shifts_paired_a() {
        assert_eq!(MoleculeId::PairedA(0).with_offset(7), MoleculeId::PairedA(7));
        assert_eq!(MoleculeId::PairedA(5).with_offset(100), MoleculeId::PairedA(105));
    }

    #[test]
    fn molecule_id_with_offset_shifts_paired_b() {
        assert_eq!(MoleculeId::PairedB(0).with_offset(7), MoleculeId::PairedB(7));
        assert_eq!(MoleculeId::PairedB(5).with_offset(100), MoleculeId::PairedB(105));
    }

    #[test]
    fn molecule_id_with_offset_zero_is_identity() {
        for mi in [
            MoleculeId::None,
            MoleculeId::Single(0),
            MoleculeId::Single(42),
            MoleculeId::PairedA(7),
            MoleculeId::PairedB(7),
        ] {
            assert_eq!(mi.with_offset(0), mi, "with_offset(0) must be identity for {mi:?}");
        }
    }

    #[test]
    fn molecule_id_with_offset_composes_by_addition() {
        for mi in [MoleculeId::Single(3), MoleculeId::PairedA(11), MoleculeId::PairedB(11)] {
            let a: u64 = 5;
            let b: u64 = 17;
            assert_eq!(
                mi.with_offset(a).with_offset(b),
                mi.with_offset(a + b),
                "composition of with_offset must be additive for {mi:?}",
            );
        }
    }

    /// `with_offset` must panic on `u64` overflow in both debug and release
    /// builds — `checked_add` keeps the failure mode consistent across
    /// profiles instead of wrapping silently in release (which would let
    /// distinct molecules collide on the same MI).
    #[test]
    #[should_panic(expected = "MoleculeId offset overflow")]
    fn molecule_id_with_offset_panics_on_overflow() {
        let _ = MoleculeId::Single(u64::MAX).with_offset(1);
    }

    /// Same overflow guard for `write_with_offset`, the codepath that
    /// actually emits the MI:Z tag bytes during BAM serialization.
    #[test]
    #[should_panic(expected = "MoleculeId offset overflow")]
    fn molecule_id_write_with_offset_panics_on_overflow() {
        let mut buf = String::new();
        let _ = MoleculeId::PairedA(u64::MAX).write_with_offset(1, &mut buf);
    }

    #[test]
    fn test_consensus_revcomp_returns_correct_tags() {
        // fgbio ConsensusTags.PerBase.TagsToReverseComplement = {ac, bc} (per-base base arrays)
        assert_eq!(TagSets::CONSENSUS_REVCOMP.len(), 2);
        assert_eq!(TagSets::CONSENSUS_REVCOMP, &["ac", "bc"]);
    }

    #[test]
    fn test_taginfo_new_with_empty_lists() {
        let tag_info = TagInfo::new(vec![], vec![], vec![]);
        assert!(tag_info.remove.is_empty());
        assert!(tag_info.reverse.is_empty());
        assert!(tag_info.revcomp.is_empty());
    }

    #[test]
    fn test_taginfo_new_with_individual_tags() {
        let tag_info = TagInfo::new(
            vec!["AS".to_string(), "NM".to_string()],
            vec!["BQ".to_string()],
            vec!["E2".to_string(), "U2".to_string()],
        );
        assert_eq!(tag_info.remove.len(), 2);
        assert!(tag_info.remove.contains("AS"));
        assert!(tag_info.remove.contains("NM"));
        assert_eq!(tag_info.reverse.len(), 1);
        assert!(tag_info.reverse.contains("BQ"));
        assert_eq!(tag_info.revcomp.len(), 2);
        assert!(tag_info.revcomp.contains("E2"));
        assert!(tag_info.revcomp.contains("U2"));
    }

    #[test]
    fn test_taginfo_new_with_consensus_reverse() {
        let tag_info = TagInfo::new(vec![], vec!["Consensus".to_string()], vec![]);
        assert_eq!(tag_info.reverse.len(), 8);
        for t in ["cd", "ce", "ad", "ae", "bd", "be", "aq", "bq"] {
            assert!(tag_info.reverse.contains(t), "missing reverse tag {t}");
        }
    }

    #[test]
    fn test_taginfo_new_with_consensus_revcomp() {
        let tag_info = TagInfo::new(vec![], vec![], vec!["Consensus".to_string()]);
        assert_eq!(tag_info.revcomp.len(), 2);
        assert!(tag_info.revcomp.contains("ac"));
        assert!(tag_info.revcomp.contains("bc"));
    }

    #[test]
    fn test_taginfo_new_with_consensus_and_individual_tags() {
        let tag_info = TagInfo::new(
            vec!["AS".to_string()],
            vec!["Consensus".to_string(), "BQ".to_string()],
            vec!["Consensus".to_string(), "E2".to_string()],
        );
        // Remove set
        assert_eq!(tag_info.remove.len(), 1);
        assert!(tag_info.remove.contains("AS"));
        // Reverse set should have the 8 Consensus tags + BQ
        assert_eq!(tag_info.reverse.len(), 9);
        for t in ["cd", "ce", "ad", "ae", "bd", "be", "aq", "bq", "BQ"] {
            assert!(tag_info.reverse.contains(t), "missing reverse tag {t}");
        }
        // Revcomp set should have the 2 Consensus base tags + E2
        assert_eq!(tag_info.revcomp.len(), 3);
        for t in ["ac", "bc", "E2"] {
            assert!(tag_info.revcomp.contains(t), "missing revcomp tag {t}");
        }
    }

    #[test]
    fn test_taginfo_new_with_duplicate_tags() {
        let tag_info = TagInfo::new(
            vec!["AS".to_string(), "AS".to_string()],
            vec!["BQ".to_string(), "BQ".to_string()],
            vec!["E2".to_string(), "E2".to_string()],
        );
        // HashSet should deduplicate
        assert_eq!(tag_info.remove.len(), 1);
        assert_eq!(tag_info.reverse.len(), 1);
        assert_eq!(tag_info.revcomp.len(), 1);
    }

    #[test]
    fn test_taginfo_new_with_multiple_consensus_references() {
        let tag_info = TagInfo::new(
            vec![],
            vec!["Consensus".to_string(), "Consensus".to_string()],
            vec!["Consensus".to_string(), "Consensus".to_string()],
        );
        // Should not duplicate consensus tags
        assert_eq!(tag_info.reverse.len(), 8);
        assert_eq!(tag_info.revcomp.len(), 2);
    }

    #[test]
    fn test_has_revs_or_revcomps_with_both_empty() {
        let tag_info = TagInfo::new(vec!["AS".to_string()], vec![], vec![]);
        assert!(!tag_info.has_revs_or_revcomps());
    }

    #[test]
    fn test_has_revs_or_revcomps_with_reverse_only() {
        let tag_info = TagInfo::new(vec![], vec!["BQ".to_string()], vec![]);
        assert!(tag_info.has_revs_or_revcomps());
    }

    #[test]
    fn test_has_revs_or_revcomps_with_revcomp_only() {
        let tag_info = TagInfo::new(vec![], vec![], vec!["E2".to_string()]);
        assert!(tag_info.has_revs_or_revcomps());
    }

    #[test]
    fn test_has_revs_or_revcomps_with_both() {
        let tag_info = TagInfo::new(vec![], vec!["BQ".to_string()], vec!["E2".to_string()]);
        assert!(tag_info.has_revs_or_revcomps());
    }

    #[test]
    fn test_has_revs_or_revcomps_with_consensus() {
        let tag_info = TagInfo::new(vec![], vec!["Consensus".to_string()], vec![]);
        assert!(tag_info.has_revs_or_revcomps());
    }

    #[test]
    fn test_taginfo_clone() {
        let tag_info =
            TagInfo::new(vec!["AS".to_string()], vec!["BQ".to_string()], vec!["E2".to_string()]);
        let cloned = tag_info.clone();
        assert_eq!(cloned.remove.len(), tag_info.remove.len());
        assert_eq!(cloned.reverse.len(), tag_info.reverse.len());
        assert_eq!(cloned.revcomp.len(), tag_info.revcomp.len());
    }

    #[test]
    fn test_extract_mi_base_simple() {
        assert_eq!(extract_mi_base("123"), "123");
        assert_eq!(extract_mi_base("A"), "A");
    }

    #[test]
    fn test_extract_mi_base_with_duplex_suffix() {
        // Truncates at the last `/`, matching fgbio `ReviewConsensusVariants.toMi`
        // (ReviewConsensusVariantsTest.scala:157) — the suffix is stripped whatever
        // it is, not only `/A` and `/B`.
        assert_eq!(extract_mi_base("123/A"), "123");
        assert_eq!(extract_mi_base("123/B"), "123");
        assert_eq!(extract_mi_base("123/456/A"), "123/456");
        // Any suffix after the last `/` is stripped (previously left intact).
        assert_eq!(extract_mi_base("123/C"), "123");
        assert_eq!(extract_mi_base("123/ReallyLongSuffix"), "123");
        assert_eq!(extract_mi_base("123/456"), "123");
        // A leading `/` (index 0) is preserved (fgbio's `slash > 0` guard).
        assert_eq!(extract_mi_base("/A"), "/A");
    }

    // =========================================================================
    // validate_umi tests - matching fgbio behavior
    // =========================================================================

    #[test]
    fn test_validate_umi_uppercase_acgt() {
        assert_eq!(validate_umi(b"ACGT"), UmiValidation::Valid(4));
        assert_eq!(validate_umi(b"AAAAAAAA"), UmiValidation::Valid(8));
        assert_eq!(validate_umi(b"TTTTTTTT"), UmiValidation::Valid(8));
    }

    #[test]
    fn test_validate_umi_lowercase_acgt() {
        // Lowercase bases should be counted (matches fgbio's toUpperCase behavior)
        assert_eq!(validate_umi(b"acgt"), UmiValidation::Valid(4));
        assert_eq!(validate_umi(b"AcGt"), UmiValidation::Valid(4));
    }

    #[test]
    fn test_validate_umi_uppercase_n_rejected() {
        // Uppercase N should be rejected (matches fgbio's .contains('N'))
        assert_eq!(validate_umi(b"ACNT"), UmiValidation::ContainsN);
        assert_eq!(validate_umi(b"NACGT"), UmiValidation::ContainsN);
        assert_eq!(validate_umi(b"ACGTN"), UmiValidation::ContainsN);
        assert_eq!(validate_umi(b"NNNN"), UmiValidation::ContainsN);
    }

    #[test]
    fn test_validate_umi_lowercase_n_skipped() {
        // Lowercase 'n' should be skipped, not rejected (matches fgbio)
        // fgbio's .contains('N') is case-sensitive
        assert_eq!(validate_umi(b"ACnT"), UmiValidation::Valid(3));
        assert_eq!(validate_umi(b"acnt"), UmiValidation::Valid(3));
        assert_eq!(validate_umi(b"nnnn"), UmiValidation::Valid(0));
    }

    #[test]
    fn test_validate_umi_dash_skipped() {
        // Dash is skipped for paired UMIs
        assert_eq!(validate_umi(b"ACGT-TGCA"), UmiValidation::Valid(8));
        assert_eq!(validate_umi(b"----"), UmiValidation::Valid(0));
    }

    #[test]
    fn test_validate_umi_other_chars_skipped() {
        // Other characters are silently skipped
        assert_eq!(validate_umi(b"ACGT+TGCA"), UmiValidation::Valid(8));
        assert_eq!(validate_umi(b"AC GT"), UmiValidation::Valid(4));
    }

    #[test]
    fn test_validate_umi_empty() {
        assert_eq!(validate_umi(b""), UmiValidation::Valid(0));
    }

    #[test]
    fn test_validate_umi_mixed_case_with_uppercase_n() {
        // Mixed case with uppercase N should still be rejected
        assert_eq!(validate_umi(b"acNt"), UmiValidation::ContainsN);
        assert_eq!(validate_umi(b"AcNt"), UmiValidation::ContainsN);
    }
}
