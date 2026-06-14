//! Tests for `SortSpillDecompress` (v4 — per-slot bounded queue design).
//!
//! TODO(P6): rewrite the test suite per spec v4. Pre-v4 tests covered
//! gap-filler admission, `try_claim_one_task` reservation, and
//! `decomp_in_flight` accounting — all removed in v4. The two
//! retained tests below pin construction + registry sharing across
//! Parallel-step worker copies.

use super::*;

#[test]
fn profile_advertises_parallel_no_ordering() {
    let step = SortSpillDecompress::new(16);
    let profile = step.profile();
    assert_eq!(profile.name, "SortSpillDecompress");
    assert_eq!(profile.kind, StepKind::Parallel);
    assert!(!profile.sticky);
    assert_eq!(profile.branch_ordering.len(), 1);
    assert_eq!(profile.branch_ordering[0], BranchOrdering::None);
}

#[test]
fn clone_shares_registry_across_worker_copies() {
    use std::sync::Arc as StdArc;
    let step = SortSpillDecompress::new(16);
    let cloned = step.new_worker_copy();
    // Both clones should point at the same registry Arc.
    assert!(StdArc::ptr_eq(&step.registry, &cloned.registry));
}
