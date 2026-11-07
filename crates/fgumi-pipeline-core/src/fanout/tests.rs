use std::sync::atomic::Ordering;

use rstest::rstest;

use super::OrderedInbox;
use crate::item::HeapSize;
use crate::runtime::wake::SideWake;

#[derive(Debug)]
struct Item(usize);

impl HeapSize for Item {
    fn heap_size(&self) -> usize {
        self.0
    }
}

#[rstest]
#[case::reverse_then_forward(&[2, 0, 1])]
#[case::in_order(&[0, 1, 2])]
fn takes_in_ordinal_order(#[case] deliveries: &[u64]) {
    let inbox = OrderedInbox::<Item>::new(true);
    assert!(inbox.try_take_next().is_none());
    for &o in deliveries {
        let _ = inbox.deliver(o, Item(1));
        if o != 0 && deliveries.first() == Some(&o) {
            assert!(inbox.try_take_next().is_none(), "0 has not arrived");
        }
    }
    let got: Vec<u64> = std::iter::from_fn(|| inbox.try_take_next().map(|(o, _)| o)).collect();
    assert_eq!(got, vec![0, 1, 2]);
}

#[derive(Clone, Copy, Debug)]
enum Op {
    Reserve(u64),
    Deliver(usize),
    DeliverTake(usize),
    ReReserve(u64, u64),
    UnreservedDeliverTake(usize),
}

#[rstest]
#[case::reserve(Op::Reserve(100), 100)]
#[case::deliver(Op::Deliver(40), 40)]
#[case::take(Op::DeliverTake(40), 0)]
#[case::re_reserve_replaces_estimate(Op::ReReserve(100, 30), 30)]
#[case::unreserved_delivery_then_take(Op::UnreservedDeliverTake(40), 0)]
fn gauge_tracks_reservation_delivery_take(#[case] op: Op, #[case] expected: u64) {
    let inbox = OrderedInbox::<Item>::new(true);
    match op {
        Op::Reserve(e) => inbox.reserve(0, e),
        Op::Deliver(n) => {
            inbox.reserve(0, 100);
            let _ = inbox.deliver(0, Item(n));
        }
        Op::DeliverTake(n) => {
            inbox.reserve(0, 100);
            let _ = inbox.deliver(0, Item(n));
            assert!(inbox.try_take_next().is_some());
        }
        Op::ReReserve(first, second) => {
            inbox.reserve(0, first);
            inbox.reserve(0, second);
        }
        Op::UnreservedDeliverTake(n) => {
            let _ = inbox.deliver(0, Item(n));
            assert_eq!(inbox.in_flight_bytes(), n as u64, "an unreserved delivery is charged");
            assert!(inbox.try_take_next().is_some());
        }
    }
    assert_eq!(inbox.in_flight_bytes(), expected);
}

#[test]
fn deliver_wakes_only_the_awaited_ordinal() {
    let inbox = OrderedInbox::<Item>::new(true);
    let _ = inbox.deliver(0, Item(1));
    assert!(inbox.try_take_next().is_some());
    assert!(!inbox.await_next(), "next = 1 is absent");
    assert_eq!(inbox.deliver(2, Item(1)), SideWake::Unavailable);
    assert_eq!(inbox.wake_attempts.load(Ordering::Relaxed), 0);
    assert_eq!(inbox.deliver(1, Item(1)), SideWake::Unavailable);
    assert_eq!(inbox.wake_attempts.load(Ordering::Relaxed), 1);
}

#[test]
fn unordered_wakes_when_empty_and_waiting() {
    let inbox = OrderedInbox::<Item>::new(false);
    let _ = inbox.deliver(5, Item(1));
    assert_eq!(inbox.wake_attempts.load(Ordering::Relaxed), 0, "owner not waiting");
    assert!(inbox.try_take_any().is_some());
    assert!(!inbox.await_next());
    let _ = inbox.deliver(3, Item(1));
    assert_eq!(inbox.wake_attempts.load(Ordering::Relaxed), 1);
}

#[test]
fn await_next_true_when_present() {
    let inbox = OrderedInbox::<Item>::new(true);
    let _ = inbox.deliver(0, Item(1));
    assert!(inbox.await_next());
    let _ = inbox.deliver(1, Item(1));
    assert_eq!(inbox.wake_attempts.load(Ordering::Relaxed), 0, "not waiting");
}

#[test]
fn unordered_take_any() {
    let inbox = OrderedInbox::<Item>::new(false);
    let _ = inbox.deliver(5, Item(1));
    let _ = inbox.deliver(3, Item(1));
    let got: Vec<u64> = std::iter::from_fn(|| inbox.try_take_any().map(|(o, _)| o)).collect();
    assert_eq!(got, vec![3, 5]);
}

/// A take method that does not match the inbox mode panics, release builds
/// included.
#[rstest]
#[case::any_on_ordered(true, false)]
#[case::next_on_unordered(false, true)]
#[should_panic(expected = "inbox; use")]
fn mismatched_take_method_panics(#[case] ordered: bool, #[case] take_next: bool) {
    let inbox = OrderedInbox::<Item>::new(ordered);
    let _ = if take_next { inbox.try_take_next() } else { inbox.try_take_any() };
}

/// A broken ordinal contract panics, release builds included, instead of
/// dropping a result or leaving a reservation nothing releases.
#[rstest]
#[case::duplicate_delivery(&[Op::Deliver(1), Op::Deliver(2)], "duplicate delivery for ordinal 0")]
#[case::reserve_after_delivery(&[Op::Deliver(1), Op::Reserve(8)], "which is already delivered")]
#[case::delivery_below_next(&[Op::DeliverTake(1), Op::Deliver(1)], "delivery for ordinal 0 below next 1")]
#[case::reserve_below_next(&[Op::DeliverTake(1), Op::Reserve(8)], "reserve for ordinal 0 below next 1")]
fn ordinal_contract_violation_panics(#[case] ops: &[Op], #[case] message: &str) {
    let inbox = OrderedInbox::<Item>::new(true);
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        for &op in ops {
            match op {
                Op::Reserve(e) => inbox.reserve(0, e),
                Op::Deliver(n) => {
                    let _ = inbox.deliver(0, Item(n));
                }
                Op::DeliverTake(n) => {
                    let _ = inbox.deliver(0, Item(n));
                    assert!(inbox.try_take_next().is_some());
                }
                Op::ReReserve(..) | Op::UnreservedDeliverTake(_) => unreachable!("not a case op"),
            }
        }
    }));
    let payload = result.expect_err("the contract violation panics");
    let got = payload
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| payload.downcast_ref::<&str>().map(|s| (*s).to_owned()))
        .unwrap_or_default();
    assert!(got.contains(message), "panic {got:?} does not name {message:?}");
}

#[test]
fn describe_reports_state() {
    use crate::step::SideChannel;
    let inbox = OrderedInbox::<Item>::new(true);
    let _ = inbox.deliver(1, Item(7));
    assert_eq!(inbox.describe(), "inbox(next=0, pending=1)");
}
