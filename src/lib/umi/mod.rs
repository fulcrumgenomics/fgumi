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

use crate::template::Template;
use anyhow::{Result, bail};
use fgumi_raw_bam::SingleEndInsertEnd;
use log::warn;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, Ordering};

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

/// The single-end insert-end hints (`tc` tags on unpaired primaries) that
/// [`partition_for_umi_assignment`] rejected; see [`Template::single_end_insert_end`].
///
/// A rejected hint is not an error: the read is grouped by its 5' end only, as without a
/// tag. It usually means the upstream tool wrote the tag wrongly (e.g. 0-based), so `group`
/// and `dedup` warn once per run with the first rejected read and log the run's total.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub(crate) struct RejectedInsertEndHints {
    /// Templates whose hint was rejected.
    pub(crate) count: u64,
    /// The first rejected template's read name and the reason it was rejected.
    pub(crate) first: Option<(Vec<u8>, &'static str)>,
}

impl RejectedInsertEndHints {
    /// Records one rejected hint, on the template named `name`, for `reason`.
    fn record(&mut self, name: &[u8], reason: &'static str) {
        self.count += 1;
        if self.first.is_none() {
            self.first = Some((name.to_vec(), reason));
        }
    }

    /// Warns about the first rejected hint, unless `warned` (the run's flag) says a warning
    /// was already logged. Called after each position group, so a run warns at most once.
    pub(crate) fn warn_first_once(&self, warned: &AtomicBool) {
        if let Some((name, reason)) = &self.first
            && !warned.swap(true, Ordering::Relaxed)
        {
            warn!(
                "Ignoring the tc insert-end hint on single-end read {}: {reason}. Reads whose \
                 hint is rejected are grouped by their 5' end only; the run's total is logged \
                 at the end. See the \"Single-End Reads\" section of the UMI grouping guide.",
                String::from_utf8_lossy(name)
            );
        }
    }
}

/// Logs the run's total of rejected single-end insert-end hints, when there were any.
pub(crate) fn log_rejected_insert_end_hints(total: u64) {
    if total > 0 {
        warn!(
            "Rejected the tc insert-end hint on {total} single-end reads; they were grouped by \
             their 5' end only."
        );
    }
}

/// The key [`partition_for_umi_assignment`] splits templates by: the pair orientation (when
/// split by orientation) and the single-end insert end `(tid, 1-based position, is_reverse)`.
type AssignmentSubgroupKey = (Option<(bool, bool)>, Option<(i32, i32, bool)>);

/// Partitions a position group's templates into the subgroups `group` and `dedup` assign
/// molecule ids in separately, returning each subgroup's template indices, in input order.
///
/// Templates are split:
///
/// - by pair orientation ([`Template::pair_strand_orientation`]), when `split_by_orientation`. Each
///   command decides whether to split by orientation (see `group::assign_umi_groups_impl` and
///   `dedup::splits_by_strand_of_origin`);
/// - always, by the single-end insert end ([`Template::single_end_insert_end`]): single-end
///   reads at the same 5' position whose insert ends differ come from different molecules,
///   whatever their UMIs. A template without a valid insert end (every paired template, and
///   every single-end template without a `tc` tag or whose tag is rejected) has insert end
///   `None`, so it is never grouped with a single-end template that has one. Rejected tags
///   are recorded in `rejected`.
///
/// The position key and the template-coordinate sort stay 5'-only; splitting here, inside the
/// position group, is what lets reads whose insert ends differ interleave in the input without
/// splitting molecules.
///
/// Subgroups are returned in ascending `(orientation, insert end)` order (a `None` insert end
/// first), so molecule ids are reproducible. When every template has the same key (the common
/// case, e.g. no orientation split and no insert ends), the result is one subgroup of every
/// index, exactly as without partitioning.
pub(crate) fn partition_for_umi_assignment(
    templates: &[Template],
    split_by_orientation: bool,
    rejected: &mut RejectedInsertEndHints,
) -> Vec<Vec<usize>> {
    let keys: Vec<AssignmentSubgroupKey> = templates
        .iter()
        .map(|template| {
            let insert_end = match template.single_end_insert_end() {
                SingleEndInsertEnd::End(end) => Some(end),
                SingleEndInsertEnd::NotApplicable => None,
                SingleEndInsertEnd::Rejected(reason) => {
                    rejected.record(template.name(), reason);
                    None
                }
            };
            (split_by_orientation.then(|| template.pair_strand_orientation()), insert_end)
        })
        .collect();

    // Common case: one subgroup, without building a map.
    if let Some(first) = keys.first()
        && keys.iter().all(|key| key == first)
    {
        return vec![(0..templates.len()).collect()];
    }

    let mut subgroups: BTreeMap<AssignmentSubgroupKey, Vec<usize>> = BTreeMap::new();
    for (idx, key) in keys.into_iter().enumerate() {
        subgroups.entry(key).or_default().push(idx);
    }
    subgroups.into_values().collect()
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

/// Test helpers shared by the `group` and `dedup` unit tests of single-end insert-end grouping.
#[cfg(test)]
pub(crate) mod test_support {
    use crate::sam::SamTag;
    use crate::template::Template;
    use crate::umi::MoleculeId;
    use fgumi_raw_bam::{SamBuilder, flags, testutil::encode_op};

    /// An unpaired, 50-base (`50M`) single-end template named `name` aligned at 1-based
    /// position 100 on reference 0 with UMI `umi`, on the reverse strand when `reverse`. Its
    /// unclipped 5' position is 100 forward and 149 reverse. When `insert_end` is `Some(end)`,
    /// its primary carries a `tc` tag whose virtual mate (opposite strand) places the
    /// insert's far end at `end`.
    pub(crate) fn single_end_template(
        name: &str,
        umi: &str,
        reverse: bool,
        insert_end: Option<i32>,
    ) -> Template {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .flags(if reverse { flags::REVERSE } else { 0 })
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[encode_op(0, 50)])
            .sequence(&[b'A'; 50])
            .qualities(&[30; 50]);
        b.add_string_tag(SamTag::RX, umi.as_bytes());
        if let Some(end) = insert_end {
            let tc = if reverse { [0, end, 0, 0, 149, 1] } else { [0, 100, 0, 0, end, 1] };
            b.add_array_i32(SamTag::TC, &tc);
        }
        Template::from_records(vec![b.build()]).expect("single-end template")
    }

    /// Relabels molecule ids by first appearance, so `[7, 3, 7]` becomes `[0, 1, 0]`.
    pub(crate) fn molecule_classes(templates: &[Template]) -> Vec<usize> {
        let mut seen: Vec<MoleculeId> = Vec::new();
        templates
            .iter()
            .map(|t| {
                assert!(t.mi.is_assigned(), "every template is assigned a molecule");
                seen.iter().position(|mi| *mi == t.mi).unwrap_or_else(|| {
                    seen.push(t.mi);
                    seen.len() - 1
                })
            })
            .collect()
    }
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

    // ========================================================================
    // partition_for_umi_assignment
    // ========================================================================

    use crate::umi::test_support::single_end_template;

    /// A forward/reverse paired template named `name` at 1-based 100 (R1) and 300 (R2), with
    /// R1 on the reverse strand when `r1_reverse` (R2 then forward).
    fn paired_template(name: &str, r1_reverse: bool) -> Template {
        use fgumi_raw_bam::{SamBuilder, flags, testutil::encode_op};
        let build = |first: bool, reverse: bool, pos: i32, mate_pos: i32, mate_reverse: bool| {
            let mut f =
                flags::PAIRED | if first { flags::FIRST_SEGMENT } else { flags::LAST_SEGMENT };
            if reverse {
                f |= flags::REVERSE;
            }
            if mate_reverse {
                f |= flags::MATE_REVERSE;
            }
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .flags(f)
                .ref_id(0)
                .pos(pos - 1)
                .mapq(60)
                .cigar_ops(&[encode_op(0, 50)])
                .sequence(&[b'A'; 50])
                .qualities(&[30; 50])
                .mate_ref_id(0)
                .mate_pos(mate_pos - 1);
            b.add_string_tag(crate::sam::SamTag::RX, b"AAAA");
            b.add_string_tag(crate::sam::SamTag::MC, b"50M");
            b.build()
        };
        Template::from_records(vec![
            build(true, r1_reverse, 100, 300, !r1_reverse),
            build(false, !r1_reverse, 300, 100, r1_reverse),
        ])
        .expect("paired template")
    }

    /// One input template for [`test_partition_for_umi_assignment`].
    #[derive(Clone, Copy, Debug)]
    enum Input {
        /// A single-end template: reverse strand when `true`, and its `tc` insert end.
        SingleEnd(bool, Option<i32>),
        /// A paired template with R1 on the reverse strand when `true`.
        Paired(bool),
    }

    /// Subgroups come back in ascending `(orientation, insert end)` order: orientation
    /// `(r1_positive, r2_positive)` with `false` (reverse) before `true`, then a `None` insert
    /// end before any `Some`. Each subgroup's indices are in input order.
    #[rstest]
    #[case::one_subgroup_when_all_keys_equal(
        &[Input::SingleEnd(false, None), Input::SingleEnd(false, None)],
        false,
        vec![vec![0, 1]]
    )]
    #[case::insert_ends_ascending_none_first(
        &[
            Input::SingleEnd(false, Some(300)),
            Input::SingleEnd(false, None),
            Input::SingleEnd(false, Some(250)),
            Input::SingleEnd(false, Some(300)),
            Input::SingleEnd(false, None),
        ],
        false,
        vec![vec![1, 4], vec![2], vec![0, 3]]
    )]
    #[case::forward_and_reverse_without_orientation_split(
        &[
            Input::SingleEnd(true, Some(20)),
            Input::SingleEnd(false, Some(300)),
            Input::SingleEnd(true, None),
            Input::SingleEnd(true, Some(50)),
            Input::SingleEnd(false, None),
        ],
        false,
        vec![vec![2, 4], vec![0], vec![3], vec![1]]
    )]
    #[case::forward_and_reverse_with_orientation_split(
        &[
            Input::SingleEnd(true, Some(20)),
            Input::SingleEnd(false, Some(300)),
            Input::SingleEnd(true, None),
            Input::SingleEnd(true, Some(50)),
            Input::SingleEnd(false, None),
        ],
        true,
        vec![vec![2], vec![0], vec![3], vec![4], vec![1]]
    )]
    #[case::paired_and_single_end_with_orientation_split(
        &[
            Input::Paired(false),
            Input::SingleEnd(false, Some(300)),
            Input::Paired(true),
            Input::SingleEnd(false, None),
            Input::Paired(false),
        ],
        true,
        vec![vec![2], vec![0, 4], vec![3], vec![1]]
    )]
    #[case::paired_and_single_end_without_orientation_split(
        &[
            Input::Paired(false),
            Input::SingleEnd(false, Some(300)),
            Input::Paired(true),
            Input::SingleEnd(false, None),
        ],
        false,
        vec![vec![0, 2, 3], vec![1]]
    )]
    #[case::empty(&[], true, vec![])]
    fn test_partition_for_umi_assignment(
        #[case] inputs: &[Input],
        #[case] split_by_orientation: bool,
        #[case] expected: Vec<Vec<usize>>,
    ) {
        let templates: Vec<Template> = inputs
            .iter()
            .enumerate()
            .map(|(i, input)| match *input {
                Input::SingleEnd(reverse, end) => {
                    single_end_template(&format!("q{i}"), "AAAA", reverse, end)
                }
                Input::Paired(r1_reverse) => paired_template(&format!("q{i}"), r1_reverse),
            })
            .collect();
        let mut rejected = RejectedInsertEndHints::default();
        let actual = partition_for_umi_assignment(&templates, split_by_orientation, &mut rejected);
        assert_eq!(actual, expected);
        assert_eq!(rejected, RejectedInsertEndHints::default(), "no hint is rejected");
    }

    /// A rejected hint is counted, its first template named, and the template partitioned as
    /// if it had no insert end.
    #[test]
    fn test_partition_for_umi_assignment_counts_rejected_hints() {
        let templates = vec![
            single_end_template("valid", "AAAA", false, Some(300)),
            // The far end lies inside the read (aligned 100-149), so the hint is rejected.
            single_end_template("inside", "AAAA", false, Some(120)),
            single_end_template("untagged", "AAAA", false, None),
            single_end_template("inside2", "AAAA", false, Some(110)),
        ];
        let mut rejected = RejectedInsertEndHints::default();
        let actual = partition_for_umi_assignment(&templates, false, &mut rejected);
        assert_eq!(actual, vec![vec![1, 2, 3], vec![0]]);
        assert_eq!(rejected.count, 2);
        let (name, reason) = rejected.first.expect("a rejection is recorded");
        assert_eq!(name, b"inside");
        assert!(reason.contains("not at or beyond the read's aligned extent"), "reason: {reason}");
    }
}
