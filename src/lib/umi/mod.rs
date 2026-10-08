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
use std::collections::BTreeMap;

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

/// Checks `--min-umi-length` for `group` and `dedup`: it must be at least 1, and the paired
/// strategy rejects it.
///
/// # Errors
///
/// Returns an error if `min_umi_length` is `Some(0)`, or is set with [`Strategy::Paired`].
pub(crate) fn validate_min_umi_length(
    min_umi_length: Option<usize>,
    strategy: Strategy,
) -> Result<()> {
    match min_umi_length {
        Some(_) if matches!(strategy, Strategy::Paired) => {
            bail!("Paired strategy cannot be used with --min-umi-length")
        }
        Some(0) => bail!("--min-umi-length must be at least 1"),
        _ => Ok(()),
    }
}

/// Whether `b` is kept when comparing UMIs under `--min-umi-length`: `A`, `C`, `G`, `T` or `N`
/// in either case, like fgbio's `basesOnly` (which keeps upper-case `ACGTN` of an upper-cased
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

/// Removes non-bases from every UMI (see [`remove_non_bases`]) and checks that each has at
/// least `min_len` bases.
///
/// # Errors
///
/// Returns an error if any UMI has fewer than `min_len` bases. The template filter normally
/// discards such reads before this point.
fn remove_non_bases_and_check_length(umis: &mut [String], min_len: usize) -> Result<()> {
    for umi in umis.iter_mut() {
        remove_non_bases(umi);
    }
    if let Some(shortest) = umis.iter().map(String::len).min()
        && shortest < min_len
    {
        bail!("UMI found that had shorter length than expected ({shortest} < {min_len})");
    }
    Ok(())
}

/// Truncates bases-only UMIs to the length of the shortest one.
fn truncate_to_shortest(umis: &mut [String]) {
    if let Some(shortest) = umis.iter().map(String::len).min() {
        // Bases are ASCII, so every byte offset is a char boundary.
        for umi in umis.iter_mut() {
            umi.truncate(shortest);
        }
    }
}

/// Truncates a set of UMIs to a common length for `--min-umi-length`.
///
/// With `None` the UMIs are returned unchanged. With `Some(min_len)`, every character that is
/// not a base (e.g. `-`) is removed and every UMI is truncated to the number of bases in the
/// SHORTEST UMI in `umis`, matching fgbio's `GroupReadsByUmi.truncateUmis` (as of
/// fulcrumgenomics/fgbio#1185). `min_len` only bounds that length from below. Removing the
/// non-bases means the bases of a multi-segment UMI (e.g. `ACG-TACGT`) are compared as one
/// sequence, wherever its dash falls.
///
/// # Errors
///
/// Returns an error if any UMI has fewer than `min_len` bases.
pub(crate) fn truncate_umis(
    mut umis: Vec<String>,
    min_umi_length: Option<usize>,
) -> Result<Vec<String>> {
    let Some(min_len) = min_umi_length else { return Ok(umis) };
    remove_non_bases_and_check_length(&mut umis, min_len)?;
    truncate_to_shortest(&mut umis);
    Ok(umis)
}

/// Assigns molecule ids to one position (and orientation) subgroup's UMIs, applying
/// `--min-umi-length` the way fgbio's `GroupReadsByUmi.assignUmiGroups` does.
///
/// `umis` must already be normalized for the strategy (upper-cased, or prefixed for the
/// paired strategy). With `None` they are assigned as given. With `Some(min_len)`:
///
/// - if the assigner groups only identical UMIs (identity, or `--edits 0`), the UMIs are first
///   split by their first `min_len` bases, and each split is truncated to its own shortest UMI
///   and assigned separately, so the result does not depend on the order of the input;
/// - otherwise the whole set is truncated to its shortest UMI (see [`truncate_umis`]).
///
/// Splits are assigned in ascending key order, so molecule ids are deterministic.
///
/// # Errors
///
/// Returns an error if any UMI has fewer than `min_len` bases, or if the assigner fails.
pub(crate) fn assign_umis(
    assigner: &dyn UmiAssigner,
    mut umis: Vec<String>,
    min_umi_length: Option<usize>,
) -> Result<Vec<MoleculeId>> {
    let Some(min_len) = min_umi_length else { return assigner.assign(&umis) };
    if !assigner.groups_identical_umis_only() {
        return assigner.assign(&truncate_umis(umis, min_umi_length)?);
    }

    remove_non_bases_and_check_length(&mut umis, min_len)?;
    let Some(splits) = split_by_prefix(&umis, min_len) else {
        truncate_to_shortest(&mut umis);
        return assigner.assign(&umis);
    };

    let mut assignments = vec![MoleculeId::None; umis.len()];
    for indices in splits {
        let mut split: Vec<String> =
            indices.iter().map(|&i| std::mem::take(&mut umis[i])).collect();
        truncate_to_shortest(&mut split);
        for (i, id) in indices.into_iter().zip(assigner.assign(&split)?) {
            assignments[i] = id;
        }
    }
    Ok(assignments)
}

/// Splits the indices of bases-only `umis`, each at least `len` bases long, by their first
/// `len` bases, in ascending key order. Returns `None`, without allocating, when every UMI
/// shares the same first `len` bases (one split).
fn split_by_prefix(umis: &[String], len: usize) -> Option<Vec<Vec<usize>>> {
    let first = &umis.first()?[..len];
    if umis.iter().all(|u| &u[..len] == first) {
        return None;
    }
    let mut by_prefix: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (i, umi) in umis.iter().enumerate() {
        by_prefix.entry(&umi[..len]).or_default().push(i);
    }
    Some(by_prefix.into_values().collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::collections::BTreeSet;

    fn strings(umis: &[&str]) -> Vec<String> {
        umis.iter().map(|s| (*s).to_string()).collect()
    }

    #[rstest]
    #[case::none(&["ACGTACGT", "ACG"], None, &["ACGTACGT", "ACG"])]
    #[case::none_keeps_dashes(&["ACG-TACGT", "ACG"], None, &["ACG-TACGT", "ACG"])]
    #[case::empty(&[], Some(4), &[])]
    #[case::to_option_when_shortest_equals_it(&["ACGTAC", "ACGTACGT"], Some(6), &["ACGTAC", "ACGTAC"])]
    // fgbio semantics: truncate to the shortest UMI (8), not to the option value (6).
    #[case::to_shortest_not_option(&["ACGTACGTAA", "ACGTACGT"], Some(6), &["ACGTACGT", "ACGTACGT"])]
    #[case::all_longer_than_option(&["ACGTACGT", "ACGTACGA"], Some(6), &["ACGTACGT", "ACGTACGA"])]
    // Lengths are counted in bases and dashes are removed: `ACGTACGT-A` has 9 bases.
    #[case::dash_removed(&["ACGTACGT-A", "ACGTACGTAC"], Some(8), &["ACGTACGTA", "ACGTACGTA"])]
    // Where the dash falls does not matter: both are the bases `ACGTACGTA...`.
    #[case::interleaved_dash(&["ACGTACGTA-T", "ACG-TACGTA"], Some(8), &["ACGTACGTA", "ACGTACGTA"])]
    // Every non-base character is removed, like fgbio's `basesOnly`.
    #[case::other_characters_removed(&["AC+GT.AC", "ACGTACGT"], Some(4), &["ACGTAC", "ACGTAC"])]
    #[case::lowercase_counts(&["acgtac", "ACGTACGT"], Some(4), &["acgtac", "ACGTAC"])]
    fn test_truncate_umis(
        #[case] umis: &[&str],
        #[case] min_len: Option<usize>,
        #[case] expected: &[&str],
    ) {
        let got = truncate_umis(strings(umis), min_len).expect("truncation should succeed");
        assert_eq!(got, strings(expected));
    }

    #[rstest]
    #[case::bytes(&["ACGTACGT", "ACGTAC"], 8, "(6 < 8)")]
    // `ACG-TACG` is 8 characters but 7 bases.
    #[case::dash_not_a_base(&["ACGTACGT", "ACG-TACG"], 8, "(7 < 8)")]
    fn test_truncate_umis_too_short(
        #[case] umis: &[&str],
        #[case] min_len: usize,
        #[case] detail: &str,
    ) {
        let err = truncate_umis(strings(umis), Some(min_len)).expect_err("too short must fail");
        assert!(err.to_string().ends_with(detail), "unexpected error: {err}");
    }

    #[rstest]
    #[case::unset_identity(None, Strategy::Identity, None)]
    #[case::unset_paired(None, Strategy::Paired, None)]
    #[case::one(Some(1), Strategy::Adjacency, None)]
    // fgbio#1185: `--min-umi-length` must be at least 1.
    #[case::zero(Some(0), Strategy::Identity, Some("--min-umi-length must be at least 1"))]
    #[case::paired(
        Some(4),
        Strategy::Paired,
        Some("Paired strategy cannot be used with --min-umi-length")
    )]
    fn test_validate_min_umi_length(
        #[case] min_umi_length: Option<usize>,
        #[case] strategy: Strategy,
        #[case] expected_error: Option<&str>,
    ) {
        let result = validate_min_umi_length(min_umi_length, strategy);
        assert_eq!(result.err().map(|e| e.to_string()).as_deref(), expected_error);
    }

    /// The assigners [`assign_umis`] is exercised with: every non-paired strategy, sequential
    /// and parallel (two threads).
    #[derive(Debug, Clone, Copy)]
    enum TestAssigner {
        Identity,
        ParallelIdentity,
        Edit(u32),
        ParallelEdit(u32),
        Adjacency(u32),
        ParallelAdjacency(u32),
    }

    impl TestAssigner {
        fn build(self) -> Box<dyn UmiAssigner> {
            use parallel_assigner::{
                ParallelAdjacencyAssigner, ParallelEditAssigner, ParallelIdentityAssigner,
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
            }
        }
    }

    /// Runs [`assign_umis`] and returns the partition it induces, as sets of input indices.
    fn partition(
        assigner: TestAssigner,
        umis: &[&str],
        min_len: usize,
    ) -> BTreeSet<BTreeSet<usize>> {
        let ids = assign_umis(assigner.build().as_ref(), strings(umis), Some(min_len))
            .expect("assignment should succeed");
        assert!(ids.iter().all(|id| *id != MoleculeId::None), "every UMI must be assigned");
        let mut by_id: BTreeMap<Option<usize>, BTreeSet<usize>> = BTreeMap::new();
        for (i, id) in ids.iter().enumerate() {
            by_id.entry(id.to_vec_index()).or_default().insert(i);
        }
        by_id.into_values().collect()
    }

    fn sets(groups: &[Vec<usize>]) -> BTreeSet<BTreeSet<usize>> {
        groups.iter().map(|g| g.iter().copied().collect()).collect()
    }

    /// Ports of fgbio#1185's identity / `--edits 0` tests (`GroupReadsByUmiTest.scala`, "... with
    /// the $strategy strategy, edits=$edits and $sortOrder input"): UMIs are split by their first
    /// `--min-umi-length` bases and each split is truncated to its own shortest UMI, so the result
    /// is the same for every input order. Each case is checked in the given and reversed order.
    #[rstest]
    // "truncate UMIs without dashes to the shortest UMI sharing their first bases": TTTT does not
    // share ACGT, so its length does not shorten ACGTA/ACGTC to ACGT.
    #[case::shortest_sharing_prefix(&["ACGTA", "ACGTC", "TTTT"], 4, vec![vec![0], vec![1], vec![2]])]
    // "compare only the bases of UMIs containing dashes when truncating": among the ACGTACGT split
    // the shortest UMI has nine bases, so 0 and 3 are identical wherever their dashes fall.
    #[case::dashed(
        &["ACG-TACGTA", "ACG-TACGTC", "ACG-TACGA", "ACGT-ACGTA"],
        8,
        vec![vec![0, 3], vec![1], vec![2]]
    )]
    // "truncate UMIs without dashes to the minimum UMI length".
    #[case::to_minimum(
        &["ACGTACGTA", "ACGTACGTC", "ACGTACGA", "ACGTACGT"],
        8,
        vec![vec![0, 1, 3], vec![2]]
    )]
    fn test_assign_umis_splits_identical_only_assigners_by_prefix(
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
        #[case] min_len: usize,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        assert_eq!(partition(assigner, umis, min_len), sets(&expected), "given order");

        let reversed: Vec<&str> = umis.iter().rev().copied().collect();
        let last = umis.len() - 1;
        let expected_reversed: Vec<Vec<usize>> =
            expected.iter().map(|g| g.iter().map(|i| last - i).collect()).collect();
        assert_eq!(
            partition(assigner, &reversed, min_len),
            sets(&expected_reversed),
            "reversed order"
        );
    }

    /// Ports of fgbio#1185's "compare only the bases of UMIs containing dashes when truncating
    /// with the $strategy strategy" (edit and adjacency, `--edits 1`). With edits, the whole set
    /// is truncated to its shortest UMI, and dashes are removed first.
    #[rstest]
    // ACGTT vs ACGTA: one mismatch. Comparing characters ("AC-GTT" vs "ACGTA-") would give four.
    #[case::one_mismatch(&["AC-GTT", "ACGTA-T"], 5, vec![vec![0, 1]])]
    #[case::one_mismatch_trailing_dash(&["ACGTAC-", "TCGTACG-TT"], 6, vec![vec![0, 1]])]
    // The shortest UMI has six bases, more than `--min-umi-length`: ACGTTA vs ACGTAC is two
    // mismatches, so they stay apart. Truncating to five bases (ACGTT vs ACGTA) would merge them.
    #[case::shortest_longer_than_option(&["AC-GTTA", "ACGTA-C"], 5, vec![vec![0], vec![1]])]
    fn test_assign_umis_compares_bases_with_edits(
        #[values(
            TestAssigner::Edit(1),
            TestAssigner::ParallelEdit(1),
            TestAssigner::Adjacency(1),
            TestAssigner::ParallelAdjacency(1)
        )]
        assigner: TestAssigner,
        #[case] umis: &[&str],
        #[case] min_len: usize,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        assert_eq!(partition(assigner, umis, min_len), sets(&expected));
    }

    /// With edits, UMIs are not split by prefix: `ACGTA`/`ACGTC`/`TTTT` with `--min-umi-length 4`
    /// truncate together to `ACGT`, `ACGT` and `TTTT`.
    #[rstest]
    fn test_assign_umis_does_not_split_with_edits(
        #[values(TestAssigner::Edit(1), TestAssigner::ParallelEdit(1))] assigner: TestAssigner,
    ) {
        assert_eq!(
            partition(assigner, &["ACGTA", "ACGTC", "TTTT"], 4),
            sets(&[vec![0, 1], vec![2]])
        );
    }

    #[test]
    fn test_assign_umis_without_min_length_assigns_as_given() {
        let assigner = Strategy::Identity.new_assigner(0);
        let ids = assign_umis(assigner.as_ref(), strings(&["ACGTA", "ACGTC"]), None)
            .expect("assignment should succeed");
        assert_ne!(ids[0], ids[1]);
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
