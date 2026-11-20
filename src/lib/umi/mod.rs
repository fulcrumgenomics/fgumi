//! UMI (Unique Molecular Identifier) utilities
//!
//! This module re-exports functionality from the `fgumi-umi` crate.

// Re-export the assigner module and top-level items from fgumi-umi
pub use fgumi_umi::{
    MoleculeId, TagInfo, TagSets, UmiValidation, assigner, extract_mi_base, validate_umi,
};

// Re-export commonly used items from the assigner module for convenience
pub use fgumi_umi::assigner::{
    AdjacencyUmiAssigner, IdentityUmiAssigner, PairedUmiAssigner, SimpleErrorUmiAssigner, Strategy,
    Umi, UmiAssigner,
};

pub mod parallel_assigner;
pub mod read_name;

use anyhow::{Result, bail};
use fgumi_umi::assigner::assign_by_segment_lengths;

/// Prefix of the error `group` and `dedup` report when UMI assignment fails for a
/// position group (e.g. mixed UMI lengths, or an unparseable UMI for the strategy).
pub(crate) const ASSIGN_UMI_GROUPS_ERROR: &str = "Failed to assign UMI groups";

/// Wraps a UMI-assignment failure as the `io::Error` `group` and `dedup` return from
/// their pipeline steps, keeping the full cause chain (`{e:#}`).
pub(crate) fn assign_umi_groups_error(e: &anyhow::Error) -> std::io::Error {
    std::io::Error::new(
        std::io::ErrorKind::InvalidData,
        format!("{ASSIGN_UMI_GROUPS_ERROR}: {e:#}"),
    )
}

/// Checks `--min-umi-length` and `--truncate` for `group` and `dedup`: `--min-umi-length` must
/// be at least 1, `--truncate` requires `--min-umi-length`, and the paired strategy rejects
/// both.
///
/// # Errors
///
/// Returns an error if `min_umi_length` or `truncate` is set with [`Strategy::Paired`] (checked
/// first, so a paired run with `--truncate` alone gets the paired message), if `min_umi_length`
/// is `Some(0)`, or if `truncate` is set without `min_umi_length`.
pub(crate) fn validate_min_umi_length(
    min_umi_length: Option<usize>,
    truncate: bool,
    strategy: Strategy,
) -> Result<()> {
    if matches!(strategy, Strategy::Paired) && (min_umi_length.is_some() || truncate) {
        bail!(
            "--min-umi-length and --truncate cannot be used with the paired strategy; paired \
             UMIs of different lengths are grouped by the lengths of their two UMIs"
        );
    }
    match min_umi_length {
        Some(0) => bail!("--min-umi-length must be at least 1"),
        None if truncate => bail!("--truncate requires --min-umi-length"),
        _ => Ok(()),
    }
}

/// Logged once per `group` or `dedup` stage when `--min-umi-length` is given without
/// `--truncate`, because earlier releases also truncated UMIs with `--min-umi-length` alone.
pub(crate) const MIN_UMI_LENGTH_WITHOUT_TRUNCATE_NOTICE: &str = "--min-umi-length only discards reads with too few UMI bases; UMIs of different lengths are \
     grouped separately. Add --truncate to group them together.";

/// The notice to log for `--min-umi-length` without `--truncate` (see
/// [`MIN_UMI_LENGTH_WITHOUT_TRUNCATE_NOTICE`]), or `None` when it does not apply: no
/// `--min-umi-length`, `--truncate` given, or `--no-umi` (which ignores UMIs entirely).
pub(crate) fn min_umi_length_notice(
    min_umi_length: Option<usize>,
    truncate: bool,
    no_umi: bool,
) -> Option<&'static str> {
    (min_umi_length.is_some() && !truncate && !no_umi)
        .then_some(MIN_UMI_LENGTH_WITHOUT_TRUNCATE_NOTICE)
}

/// The number of bases to truncate each UMI to before grouping: `--min-umi-length` when
/// `--truncate` is given, otherwise `None` (UMIs are grouped as they are).
pub(crate) fn truncate_length(min_umi_length: Option<usize>, truncate: bool) -> Option<usize> {
    min_umi_length.filter(|_| truncate)
}

/// Whether `b` is kept when truncating UMIs with `--truncate`: `A`, `C`, `G`, `T` or `N` in
/// either case, like fgbio's `basesOnly` (which keeps upper-case `ACGTN` of an upper-cased
/// UMI). `N` cannot reach here from the CLI, since the template filter discards UMIs with one.
const fn is_umi_base(b: u8) -> bool {
    matches!(b, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' | b'N' | b'n')
}

/// Removes every character that is not a UMI base (see [`is_umi_base`]), e.g. the `-` between
/// UMI segments, in place. Does nothing, and scans `umi` only once, when it holds only bases.
fn remove_non_bases(umi: &mut String) {
    if !umi.bytes().all(is_umi_base) {
        umi.retain(|c| u8::try_from(c).is_ok_and(is_umi_base));
    }
}

/// Truncates every UMI to exactly `length` bases for `--truncate`, matching fgbio's
/// `GroupReadsByUmi.umiForGrouping` (as of fulcrumgenomics/fgbio#1198): every character that
/// is not a base (e.g. `-`) is removed, so the position of a dash does not matter, and the
/// first `length` bases are kept.
///
/// # Errors
///
/// Returns an error if any UMI has fewer than `length` bases. The template filter discards
/// such reads (`--min-umi-length`) before this point.
pub(crate) fn truncate_umis(umis: &mut [String], length: usize) -> Result<()> {
    for umi in umis.iter_mut() {
        remove_non_bases(umi);
        if umi.len() < length {
            bail!("UMI found that had shorter length than expected ({} < {length})", umi.len());
        }
        // Bases are ASCII, so every byte offset is a char boundary.
        umi.truncate(length);
    }
    Ok(())
}

/// Assigns molecule ids to one position (and orientation) subgroup's UMIs the way fgbio's
/// `GroupReadsByUmi.assignUmiGroups` does (as of fulcrumgenomics/fgbio#1198).
///
/// `umis` must already be normalized for the strategy (upper-cased, or prefixed for the
/// paired strategy). With `truncate_length` set (`--truncate`), every UMI is first truncated
/// to exactly that many bases (see [`truncate_umis`]), so all UMIs have the same length and
/// are assigned in one call. Otherwise the UMIs are bucketed by segment lengths and each bucket
/// is assigned separately (see [`assign_by_segment_lengths`]).
///
/// # Errors
///
/// Returns an error if any UMI has fewer than `truncate_length` bases, or if the assigner
/// fails.
pub(crate) fn assign_umis(
    assigner: &dyn UmiAssigner,
    mut umis: Vec<String>,
    truncate_length: Option<usize>,
) -> Result<Vec<MoleculeId>> {
    if let Some(length) = truncate_length {
        // Every truncated UMI is `length` bases with no dash: one bucket.
        truncate_umis(&mut umis, length)?;
        return assigner.assign(&umis);
    }
    assign_by_segment_lengths(assigner, umis)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::collections::{BTreeMap, BTreeSet};

    fn strings(umis: &[&str]) -> Vec<String> {
        umis.iter().map(|s| (*s).to_string()).collect()
    }

    #[rstest]
    #[case::empty(&[], 4, &[])]
    #[case::already_that_length(&["ACGTAC", "TTTTTT"], 6, &["ACGTAC", "TTTTTT"])]
    // fgbio#1198: truncate to the option value, not to the shortest UMI.
    #[case::to_option_not_shortest(&["ACGTACGTAA", "ACGTACGT"], 6, &["ACGTAC", "ACGTAC"])]
    #[case::mixed_lengths(&["ACTACT", "ACTAC", "ACTACTA"], 5, &["ACTAC", "ACTAC", "ACTAC"])]
    // Lengths are counted in bases and dashes are removed before truncating.
    #[case::dash_removed(&["ACGTACGT-A", "ACG-TACGTAC"], 8, &["ACGTACGT", "ACGTACGT"])]
    #[case::trailing_dash(&["ACGT-"], 4, &["ACGT"])]
    // Every non-base character is removed, like fgbio's `basesOnly`.
    #[case::other_characters_removed(&["AC+GT.AC", "ACGTACGT"], 4, &["ACGT", "ACGT"])]
    #[case::lowercase_counts(&["acgtac", "ACGTACGT"], 4, &["acgt", "ACGT"])]
    fn test_truncate_umis(#[case] umis: &[&str], #[case] length: usize, #[case] expected: &[&str]) {
        let mut umis = strings(umis);
        truncate_umis(&mut umis, length).expect("truncation should succeed");
        assert_eq!(umis, strings(expected));
    }

    #[rstest]
    #[case::bytes(&["ACGTACGT", "ACGTAC"], 8, "(6 < 8)")]
    // `ACG-TACG` is 8 characters but 7 bases.
    #[case::dash_not_a_base(&["ACGTACGT", "ACG-TACG"], 8, "(7 < 8)")]
    fn test_truncate_umis_too_short(
        #[case] umis: &[&str],
        #[case] length: usize,
        #[case] detail: &str,
    ) {
        let err = truncate_umis(&mut strings(umis), length).expect_err("too short must fail");
        assert!(err.to_string().ends_with(detail), "unexpected error: {err}");
    }

    const PAIRED_ERROR: &str = "--min-umi-length and --truncate cannot be used with the paired \
                                strategy; paired UMIs of different lengths are grouped by the \
                                lengths of their two UMIs";

    #[rstest]
    #[case::unset_identity(None, false, Strategy::Identity, None)]
    #[case::unset_paired(None, false, Strategy::Paired, None)]
    #[case::one(Some(1), false, Strategy::Adjacency, None)]
    #[case::truncate(Some(6), true, Strategy::Edit, None)]
    // fgbio#1185: `--min-umi-length` must be at least 1.
    #[case::zero(Some(0), false, Strategy::Identity, Some("--min-umi-length must be at least 1"))]
    // fgbio#1198: `--truncate` requires `--min-umi-length`.
    #[case::truncate_without_min_length(
        None,
        true,
        Strategy::Edit,
        Some("--truncate requires --min-umi-length")
    )]
    #[case::paired(Some(4), false, Strategy::Paired, Some(PAIRED_ERROR))]
    #[case::paired_truncate(Some(4), true, Strategy::Paired, Some(PAIRED_ERROR))]
    // The paired check comes first, so `--truncate` alone gets the paired message rather than
    // "--truncate requires --min-umi-length".
    #[case::paired_truncate_without_min_length(None, true, Strategy::Paired, Some(PAIRED_ERROR))]
    #[case::paired_zero_truncate(Some(0), true, Strategy::Paired, Some(PAIRED_ERROR))]
    fn test_validate_min_umi_length(
        #[case] min_umi_length: Option<usize>,
        #[case] truncate: bool,
        #[case] strategy: Strategy,
        #[case] expected_error: Option<&str>,
    ) {
        let result = validate_min_umi_length(min_umi_length, truncate, strategy);
        assert_eq!(result.err().map(|e| e.to_string()).as_deref(), expected_error);
    }

    #[rstest]
    #[case::unset(None, false, false, None)]
    #[case::min_length_only(Some(6), false, false, Some(MIN_UMI_LENGTH_WITHOUT_TRUNCATE_NOTICE))]
    #[case::truncate(Some(6), true, false, None)]
    #[case::no_umi(Some(6), false, true, None)]
    fn test_min_umi_length_notice(
        #[case] min_umi_length: Option<usize>,
        #[case] truncate: bool,
        #[case] no_umi: bool,
        #[case] expected: Option<&str>,
    ) {
        assert_eq!(min_umi_length_notice(min_umi_length, truncate, no_umi), expected);
    }

    #[rstest]
    #[case::unset(None, false, None)]
    #[case::min_length_only(Some(6), false, None)]
    #[case::truncate(Some(6), true, Some(6))]
    fn test_truncate_length(
        #[case] min_umi_length: Option<usize>,
        #[case] truncate: bool,
        #[case] expected: Option<usize>,
    ) {
        assert_eq!(truncate_length(min_umi_length, truncate), expected);
    }

    /// The assigners [`assign_umis`] is exercised with, sequential and parallel (two threads).
    #[derive(Debug, Clone, Copy)]
    enum TestAssigner {
        Identity,
        ParallelIdentity,
        Edit(u32),
        ParallelEdit(u32),
        Adjacency(u32),
        ParallelAdjacency(u32),
        Paired(u32),
        ParallelPaired(u32),
    }

    impl TestAssigner {
        fn build(self) -> Box<dyn UmiAssigner> {
            use parallel_assigner::{
                ParallelAdjacencyAssigner, ParallelEditAssigner, ParallelIdentityAssigner,
                ParallelPairedAssigner,
            };
            match self {
                Self::Identity => Strategy::Identity.new_assigner(0),
                Self::ParallelIdentity => Box::new(ParallelIdentityAssigner::new(2)),
                Self::Edit(edits) => Strategy::Edit.new_assigner(edits),
                Self::ParallelEdit(edits) => Box::new(ParallelEditAssigner::new(edits, 2)),
                Self::Adjacency(edits) => Strategy::Adjacency.new_assigner(edits),
                Self::ParallelAdjacency(edits) => {
                    Box::new(ParallelAdjacencyAssigner::new(edits, 2))
                }
                Self::Paired(edits) => Strategy::Paired.new_assigner(edits),
                Self::ParallelPaired(edits) => Box::new(ParallelPairedAssigner::new(edits, 2)),
            }
        }
    }

    /// Runs [`assign_umis`] and returns the partition it induces, as sets of input indices.
    /// The two strands of a paired molecule (`/A` and `/B`) are in one set.
    fn partition(
        assigner: TestAssigner,
        umis: &[&str],
        truncate_length: Option<usize>,
    ) -> BTreeSet<BTreeSet<usize>> {
        let ids = assign_umis(assigner.build().as_ref(), strings(umis), truncate_length)
            .expect("assignment should succeed");
        assert!(ids.iter().all(|id| *id != MoleculeId::None), "every UMI must be assigned");
        let mut by_id: BTreeMap<Option<u64>, BTreeSet<usize>> = BTreeMap::new();
        for (i, id) in ids.iter().enumerate() {
            by_id.entry(id.id()).or_default().insert(i);
        }
        by_id.into_values().collect()
    }

    fn sets(groups: &[Vec<usize>]) -> BTreeSet<BTreeSet<usize>> {
        groups.iter().map(|g| g.iter().copied().collect()).collect()
    }

    /// Port of fgbio#1198's "group single-end reads with UMIs of different lengths
    /// separately": `AAAAAAA` (index 4) is one edit from `AAAAAAAA` but shorter, so it is its
    /// own molecule.
    #[rstest]
    fn test_assign_umis_keeps_umis_of_different_lengths_apart(
        #[values(
            TestAssigner::Edit(1),
            TestAssigner::ParallelEdit(1),
            TestAssigner::Adjacency(1),
            TestAssigner::ParallelAdjacency(1)
        )]
        assigner: TestAssigner,
    ) {
        let umis = ["AAAAAAAA", "AAAAAAAA", "ACACACA", "ACACACC", "AAAAAAA"];
        assert_eq!(partition(assigner, &umis, None), sets(&[vec![0, 1], vec![2, 3], vec![4]]));
    }

    /// Port of fgbio#1198's "not group UMIs with the same number of bases but different segment
    /// lengths": the 8+8 and 7+7 dual UMIs, the 6+3 and 5+4 dual UMIs, and `ACGT` and `ACGT-`
    /// are never compared with each other.
    #[rstest]
    #[case::identity(
        TestAssigner::Identity,
        vec![vec![0, 1], vec![2], vec![3, 4], vec![5], vec![6, 7], vec![8, 9], vec![10], vec![11]]
    )]
    #[case::parallel_identity(
        TestAssigner::ParallelIdentity,
        vec![vec![0, 1], vec![2], vec![3, 4], vec![5], vec![6, 7], vec![8, 9], vec![10], vec![11]]
    )]
    #[case::edit(
        TestAssigner::Edit(1),
        vec![vec![0, 1, 2], vec![3, 4, 5], vec![6, 7], vec![8, 9], vec![10], vec![11]]
    )]
    #[case::parallel_edit(
        TestAssigner::ParallelEdit(1),
        vec![vec![0, 1, 2], vec![3, 4, 5], vec![6, 7], vec![8, 9], vec![10], vec![11]]
    )]
    #[case::adjacency(
        TestAssigner::Adjacency(1),
        vec![vec![0, 1, 2], vec![3, 4, 5], vec![6, 7], vec![8, 9], vec![10], vec![11]]
    )]
    #[case::parallel_adjacency(
        TestAssigner::ParallelAdjacency(1),
        vec![vec![0, 1, 2], vec![3, 4, 5], vec![6, 7], vec![8, 9], vec![10], vec![11]]
    )]
    fn test_assign_umis_keeps_different_segment_lengths_apart(
        #[case] assigner: TestAssigner,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        let umis = [
            "AAAAAAAA-CCCCCCCC",
            "AAAAAAAA-CCCCCCCC",
            "AAAAAAAA-CCCCCCCG",
            "AAAAAAA-CCCCCCC",
            "AAAAAAA-CCCCCCC",
            "AAAAAAA-CCCCCCG",
            "ACTACT-GGG",
            "ACTACT-GGG",
            "ACTAC-TGGG",
            "ACTAC-TGGG",
            "ACGT",
            "ACGT-",
        ];
        assert_eq!(partition(assigner, &umis, None), sets(&expected));
    }

    /// With `--truncate`, UMIs of different lengths are truncated to `--min-umi-length` bases
    /// and grouped together (ports of fgbio#1198's `--truncate` tests).
    #[rstest]
    // "truncate to the specified minimum length with --truncate".
    #[case::mixed_lengths(&["ACTACT", "ACTAC", "ACTACTA"], 5, vec![vec![0, 1, 2]])]
    // "truncate to --min-umi-length even when every UMI is longer".
    #[case::all_longer(&["ACTACTC", "ACTACTC", "ACTACTG", "ACTACTG"], 6, vec![vec![0, 1, 2, 3]])]
    // "truncate UMIs to --min-umi-length with --truncate".
    #[case::to_min_length(&["ACGTA", "ACGTC", "TTTT", "ACGT"], 4, vec![vec![0, 1, 3], vec![2]])]
    // "compare only the first --min-umi-length bases of UMIs containing dashes".
    #[case::dashed(
        &["ACG-TACGTA", "ACG-TACGTC", "ACG-TACGA", "ACGT-ACGTA"],
        8,
        vec![vec![0, 1, 3], vec![2]]
    )]
    fn test_assign_umis_truncates_identical_only(
        #[values(
            TestAssigner::Identity,
            TestAssigner::ParallelIdentity,
            TestAssigner::Edit(0),
            TestAssigner::ParallelEdit(0),
            TestAssigner::Adjacency(0),
            TestAssigner::ParallelAdjacency(0)
        )]
        assigner: TestAssigner,
        #[case] umis: &[&str],
        #[case] length: usize,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        assert_eq!(partition(assigner, umis, Some(length)), sets(&expected));
    }

    /// Without `--truncate`, the identity / `--edits 0` cases above keep every length (and every
    /// dash layout) apart (fgbio#1198's "not truncate UMIs containing dashes without
    /// --truncate" and "group only identical UMIs of the same length").
    #[rstest]
    #[case::dashed(
        &["ACG-TACGTA", "ACG-TACGTC", "ACG-TACGA", "ACGT-ACGTA"],
        vec![vec![0], vec![1], vec![2], vec![3]]
    )]
    #[case::same_length_only(
        &["ACGTA", "ACGTA", "ACGT", "ACGT-", "ACGTAC"],
        vec![vec![0, 1], vec![2], vec![3], vec![4]]
    )]
    fn test_assign_umis_without_truncate_identical_only(
        #[values(
            TestAssigner::Identity,
            TestAssigner::ParallelIdentity,
            TestAssigner::Edit(0),
            TestAssigner::ParallelEdit(0),
            TestAssigner::Adjacency(0),
            TestAssigner::ParallelAdjacency(0)
        )]
        assigner: TestAssigner,
        #[case] umis: &[&str],
        #[case] expected: Vec<Vec<usize>>,
    ) {
        assert_eq!(partition(assigner, umis, None), sets(&expected));
    }

    /// Ports of fgbio#1198's "compare only the bases of UMIs containing dashes when truncating"
    /// and "not group UMIs whose segments have different lengths without --truncate" (edit
    /// and adjacency, `--edits 1`).
    #[rstest]
    // ACGTT vs ACGTA: one mismatch. Comparing characters ("AC-GTT" vs "ACGTA-") would give four.
    #[case::one_mismatch(&["AC-GTT", "ACGTA-T"], Some(5), vec![vec![0, 1]])]
    #[case::one_mismatch_trailing_dash(&["ACGTAC-", "TCGTACG-TT"], Some(6), vec![vec![0, 1]])]
    // Both have six bases, more than `--min-umi-length`, and ACGTTA vs ACGTAC is two
    // mismatches. UMIs are truncated to `--min-umi-length`, not to the shortest UMI, so ACGTT
    // vs ACGTA (one mismatch) are grouped.
    #[case::truncated_to_option(&["AC-GTTA", "ACGTA-C"], Some(5), vec![vec![0, 1]])]
    // One mismatch apart if the dashes are ignored, but the segment lengths differ (2+4 vs 5+1).
    #[case::segment_lengths_differ(&["AC-GTTA", "ACGTA-A"], None, vec![vec![0], vec![1]])]
    fn test_assign_umis_dashed_with_edits(
        #[values(
            TestAssigner::Edit(1),
            TestAssigner::ParallelEdit(1),
            TestAssigner::Adjacency(1),
            TestAssigner::ParallelAdjacency(1)
        )]
        assigner: TestAssigner,
        #[case] umis: &[&str],
        #[case] truncate_length: Option<usize>,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        assert_eq!(partition(assigner, umis, truncate_length), sets(&expected));
    }

    /// The paired strategy groups mixed-length UMI pairs without error. Both orientations of a
    /// pair share their canonical segment lengths, so `A-B` and `B-A` are still the two strands
    /// of one molecule, including when the two halves differ in length.
    #[rstest]
    #[case::two_lengths(
        &["ACTACT-GGGGGG", "GGGGGG-ACTACT", "ACT-GGG", "GGG-ACT"],
        vec![vec![0, 1], vec![2, 3]]
    )]
    #[case::halves_of_different_lengths(&["ACTACT-GGG", "GGG-ACTACT"], vec![vec![0, 1]])]
    fn test_assign_umis_paired_mixed_lengths(
        #[values(TestAssigner::Paired(1), TestAssigner::ParallelPaired(1))] assigner: TestAssigner,
        #[case] umis: &[&str],
        #[case] expected: Vec<Vec<usize>>,
    ) {
        assert_eq!(partition(assigner, umis, None), sets(&expected));
    }

    /// Paired molecule ids across several length buckets are pinned, including the strand, and
    /// are the same for the sequential and parallel paired assigners. The keys carry the
    /// orientation prefixes `group` adds (`aa:`/`bb:` at `--edits 1`). The 3+3 bucket is
    /// numbered first; its two keys are the two strands of one molecule.
    #[rstest]
    #[case::unordered(
        &["aa:ACTACT-bb:GGGGGG", "aa:ACT-bb:GGG", "bb:GGG-aa:ACT"],
        vec![MoleculeId::PairedA(1), MoleculeId::PairedA(0), MoleculeId::PairedB(0)]
    )]
    #[case::halves_of_different_lengths(
        &["aa:ACTACT-bb:GGG", "bb:GGG-aa:ACTACT", "aa:ACT-bb:GGG", "aa:ACT-bb:GGGGGG"],
        vec![
            MoleculeId::PairedA(2),
            MoleculeId::PairedB(2),
            MoleculeId::PairedA(0),
            MoleculeId::PairedA(1),
        ]
    )]
    fn test_assign_umis_paired_ids_across_buckets(
        #[values(TestAssigner::Paired(1), TestAssigner::ParallelPaired(1))] assigner: TestAssigner,
        #[case] umis: &[&str],
        #[case] expected: Vec<MoleculeId>,
    ) {
        let ids = assign_umis(assigner.build().as_ref(), strings(umis), None)
            .expect("assignment should succeed");
        assert_eq!(ids, expected);
    }

    /// Buckets are assigned in ascending order of their segment lengths, from one id counter, so
    /// molecule ids are deterministic and unique across buckets whatever the input order.
    #[rstest]
    fn test_assign_umis_numbers_buckets_in_ascending_length_order(
        #[values(
            TestAssigner::Identity,
            TestAssigner::ParallelIdentity,
            TestAssigner::Edit(1),
            TestAssigner::ParallelEdit(1),
            TestAssigner::Adjacency(1),
            TestAssigner::ParallelAdjacency(1)
        )]
        assigner: TestAssigner,
    ) {
        let ids = assign_umis(
            assigner.build().as_ref(),
            strings(&["AAAAAA", "CCC", "AAAA", "CCC"]),
            None,
        )
        .expect("assignment should succeed");
        let expected: Vec<Option<u64>> = vec![Some(2), Some(0), Some(1), Some(0)];
        assert_eq!(ids.iter().map(MoleculeId::id).collect::<Vec<_>>(), expected);
    }

    /// When every UMI has the same segment lengths, [`assign_umis`] is exactly one `assign`
    /// call, so same-length data gets the same ids it did before UMIs were bucketed.
    #[rstest]
    fn test_assign_umis_same_lengths_match_a_single_assign(
        #[values(
            TestAssigner::Identity,
            TestAssigner::ParallelIdentity,
            TestAssigner::Edit(1),
            TestAssigner::ParallelEdit(1),
            TestAssigner::Adjacency(1),
            TestAssigner::ParallelAdjacency(1)
        )]
        assigner: TestAssigner,
    ) {
        let umis = strings(&["ACGT-TTGA", "ACGT-TTGA", "ACGA-TTGA", "GGCC-AATT"]);
        let direct = assigner.build().assign(&umis).expect("assignment should succeed");
        let bucketed =
            assign_umis(assigner.build().as_ref(), umis, None).expect("assignment should succeed");
        assert_eq!(bucketed, direct);
    }

    #[test]
    fn test_assign_umis_rejects_too_short() {
        let assigner = Strategy::Identity.new_assigner(0);
        let err = assign_umis(assigner.as_ref(), strings(&["ACG-TA", "ACGTAC"]), Some(6))
            .expect_err("a five-base UMI is shorter than 6");
        assert!(err.to_string().ends_with("(5 < 6)"), "unexpected error: {err}");
    }

    #[test]
    fn test_assign_umi_groups_error_keeps_cause_chain() {
        let e = anyhow::anyhow!("Multiple UMI lengths: 3, 4").context("outer");
        let io = assign_umi_groups_error(&e);
        assert_eq!(io.kind(), std::io::ErrorKind::InvalidData);
        assert_eq!(
            io.to_string(),
            "Failed to assign UMI groups: outer: Multiple UMI lengths: 3, 4"
        );
    }
}
