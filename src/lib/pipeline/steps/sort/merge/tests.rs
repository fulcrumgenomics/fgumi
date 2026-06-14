//! Profile-pin tests for `SortMerge`.
//!
//! End-to-end record-byte parity vs `RawExternalSorter::sort()` is covered
//! by the 3-step chain integration tests in
//! `src/lib/pipeline/steps/sort/tests.rs`. These tests pin the
//! step's `StepProfile` and basic affinity / construction behavior.

use fgumi_sort::{QuerynameComparator, SortOrder};

use super::*;

#[test]
fn profile_advertises_serial_byordinal_byte_bounded() {
    let s = SortMerge::new(SortOrder::Coordinate, 64 * 1024);
    let p = s.profile();
    assert_eq!(p.name, "SortMerge");
    assert_eq!(p.kind, StepKind::Serial);
    assert!(!p.sticky);
    assert_eq!(p.output_queues.len(), 1);
    match p.output_queues[0] {
        QueueSpec::ByteBounded { limit_bytes } => assert_eq!(limit_bytes, 64 * 1024),
        other => panic!("expected ByteBounded limit, got {other:?}"),
    }
    assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
}

#[test]
fn affinity_default_is_none_and_overridable() {
    let s = SortMerge::new(SortOrder::Coordinate, 64 * 1024);
    assert_eq!(s.affinity(), Affinity::None);

    let s = SortMerge::new(SortOrder::Coordinate, 64 * 1024).with_affinity(Affinity::Worker(2));
    assert_eq!(s.affinity(), Affinity::Worker(2));
}

#[test]
fn builds_for_all_sort_orders() {
    for order in [
        SortOrder::Coordinate,
        SortOrder::Queryname(QuerynameComparator::Lexicographic),
        SortOrder::Queryname(QuerynameComparator::Natural),
        SortOrder::TemplateCoordinate,
    ] {
        // Construction itself doesn't fail for any order; the typed
        // dispatch happens inside `build_driver` on transition to
        // Merging (which is exercised in the integration tests).
        let _ = SortMerge::new(order, 64 * 1024);
        let _ = SortMerge::with_target_batch_count(order, 64 * 1024, 256);
    }
}

#[test]
fn target_batch_count_floor_is_one() {
    // Defensive: a zero target_batch_count would cause the inner loop to
    // never break by-count, so we clamp to 1 in the constructor.
    let s = SortMerge::with_target_batch_count(SortOrder::Coordinate, 64 * 1024, 0);
    assert_eq!(s.target_batch_count, 1);
}
