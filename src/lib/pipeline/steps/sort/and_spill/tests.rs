//! Profile-pin tests for `SortAndSpill`.
//!
//! The end-to-end behavior (record-byte parity vs `RawExternalSorter::sort()`)
//! is covered by the 3-step chain integration tests in
//! `src/lib/pipeline/steps/sort/tests.rs`. These tests pin the
//! step's `StepProfile` and basic affinity-override behavior — failures
//! here flag a contract change without having to scroll through the
//! integration-test output.

use fgumi_sort::{RawExternalSorter, SortOrder};
use noodles::sam::Header;

use super::*;

#[test]
fn profile_advertises_serial_count_bounded() {
    let header = Header::default();
    let sorter = RawExternalSorter::new(SortOrder::Coordinate).memory_limit(1024 * 1024);
    let s = SortAndSpill::from_sorter(sorter, &header, 32).expect("from_sorter");
    let p = s.profile();
    assert_eq!(p.name, "SortAndSpill");
    assert_eq!(p.kind, StepKind::Serial);
    assert!(!p.sticky);
    assert_eq!(p.output_queues.len(), 1);
    match p.output_queues[0] {
        QueueSpec::CountBounded { capacity } => assert_eq!(capacity, 32),
        other => panic!("expected CountBounded capacity, got {other:?}"),
    }
    assert_eq!(p.branch_ordering, vec![BranchOrdering::None]);
}

#[test]
fn affinity_default_is_none_and_overridable() {
    let header = Header::default();
    let sorter = RawExternalSorter::new(SortOrder::Coordinate).memory_limit(1024 * 1024);
    let s = SortAndSpill::from_sorter(sorter, &header, 32).expect("from_sorter");
    assert_eq!(s.affinity(), Affinity::None);

    let sorter = RawExternalSorter::new(SortOrder::Coordinate).memory_limit(1024 * 1024);
    let s = SortAndSpill::from_sorter(sorter, &header, 32)
        .expect("from_sorter")
        .with_affinity(Affinity::Worker(2));
    assert_eq!(s.affinity(), Affinity::Worker(2));
}

#[test]
fn from_sorter_builds_for_all_sort_orders() {
    use fgumi_sort::QuerynameComparator;

    let header = Header::default();
    for order in [
        SortOrder::Coordinate,
        SortOrder::Queryname(QuerynameComparator::Lexicographic),
        SortOrder::Queryname(QuerynameComparator::Natural),
        SortOrder::TemplateCoordinate,
    ] {
        let sorter = RawExternalSorter::new(order).memory_limit(1024 * 1024);
        let s = SortAndSpill::from_sorter(sorter, &header, 32);
        assert!(s.is_ok(), "from_sorter must succeed for {order:?}");
    }
}
