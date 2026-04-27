//! Property: adjacency clustering is deterministic.
//!
//! Two freshly-constructed `AdjacencyUmiAssigner` instances, given the same
//! UMI multiset and the same parameters, must produce byte-for-byte identical
//! `Vec<MoleculeId>` outputs.
//!
//! Regression guard: `HashMap` iteration order non-determinism previously caused
//! equal-count UMIs to be tie-broken differently across runs (fixed by
//! introducing a lexicographic tie-break during sort). This property would
//! have caught that.

use fgumi_lib::umi::MoleculeId;
use fgumi_lib::umi::assigner::{Strategy as UmiStrategy, Umi};
use proptest::prelude::*;

/// Generate a multiset of 8-mer UMIs of varied size.
///
/// Length 8 is representative of realistic UMIs while keeping the proptest
/// input size small enough to run many iterations quickly. The size range
/// (1..64) exercises both small and moderately large inputs.
fn umi_vec_strategy() -> impl Strategy<Value = Vec<Umi>> {
    prop::collection::vec("[ACGT]{8}", 1..64)
}

proptest! {
    /// Two fresh assigners run on the same input produce identical outputs.
    #[test]
    fn adjacency_clustering_is_deterministic(umis in umi_vec_strategy()) {
        // Fresh assigners each run: counters start at 0, so exact `MoleculeId`
        // equality (not just structural equivalence) is the property.
        let assigner_a = UmiStrategy::Adjacency.new_assigner(1);
        let assigner_b = UmiStrategy::Adjacency.new_assigner(1);

        let result_a: Vec<MoleculeId> = assigner_a.assign(&umis);
        let result_b: Vec<MoleculeId> = assigner_b.assign(&umis);

        prop_assert_eq!(result_a, result_b);
    }

    /// Repeated calls on the *same* assigner produce structurally equivalent
    /// groupings (the counter advances, so IDs shift, but the partition over
    /// input indices does not change).
    #[test]
    fn adjacency_clustering_partition_is_stable(umis in umi_vec_strategy()) {
        let assigner = UmiStrategy::Adjacency.new_assigner(1);
        let first = assigner.assign(&umis);
        let second = assigner.assign(&umis);

        // Two indices share a group in `first` iff they share a group in `second`.
        prop_assert_eq!(first.len(), second.len());
        for i in 0..first.len() {
            for j in (i + 1)..first.len() {
                let same_first = first[i] == first[j];
                let same_second = second[i] == second[j];
                prop_assert_eq!(
                    same_first,
                    same_second,
                    "partition changed between repeated assign() calls"
                );
            }
        }
    }
}
