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

mod worker {
    use std::io;
    use std::sync::Arc;

    use rstest::rstest;

    use super::super::{FanoutWorker, Job, JobStep, OrderedInbox, SliceBudget};
    use crate::admission::PhaseCap;
    use crate::item::HeapSize;
    use crate::step::{Step, StepOutcome};
    use crate::testing::{StepProbe, assert_admission_contract};

    const BUDGET: SliceBudget = SliceBudget { bytes: 1024, records: 16 };

    /// A job result that reports 8 heap bytes.
    #[derive(Debug, PartialEq)]
    struct Out(u64);

    impl HeapSize for Out {
        fn heap_size(&self) -> usize {
            8
        }
    }

    struct CountJob {
        ordinal: u64,
        slices_left: u32,
        fail: bool,
    }

    impl CountJob {
        fn new(ordinal: u64, slices_left: u32) -> Self {
            Self { ordinal, slices_left, fail: false }
        }
    }

    impl HeapSize for CountJob {
        fn heap_size(&self) -> usize {
            16
        }
    }

    impl Job for CountJob {
        type Output = Out;

        fn ordinal(&self) -> u64 {
            self.ordinal
        }

        fn run_slice(&mut self, _budget: SliceBudget) -> io::Result<JobStep<Out>> {
            if self.fail {
                return Err(io::Error::other("boom"));
            }
            self.slices_left -= 1;
            Ok(if self.slices_left == 0 {
                JobStep::Done(Out(self.ordinal * 10 + 7))
            } else {
                JobStep::Continue
            })
        }
    }

    fn make(cap: Option<Arc<PhaseCap>>) -> FanoutWorker<CountJob> {
        FanoutWorker::new("count", OrderedInbox::new(true), BUDGET)
            .with_step_cap(cap.map(crate::StepCap::from))
    }

    #[test]
    fn admission_contract() {
        assert_admission_contract("fanout-test", make, CountJob::new(0, 1));
    }

    /// One step of a dispatch script.
    #[derive(Clone, Copy, Debug)]
    enum Act {
        /// `try_run` returns this outcome.
        Run(StepOutcome),
        /// Another holder takes the cap's only permit.
        Hold,
        /// That holder releases it.
        Release,
        /// The input is closed (upstream finished).
        Close,
    }

    use Act::{Close, Hold, Release, Run};
    use StepOutcome::{Capped, Finished, NoProgress, Progress};

    #[rstest]
    #[case::single_slice_job(false, &[1], &[Run(Progress)], &[(0, 7)], false, (1, 1))]
    #[case::single_slice_job_capped(true, &[1], &[Run(Progress)], &[(0, 7)], false, (1, 1))]
    #[case::multi_slice_resumes_same_clone(
        false,
        &[3, 1],
        &[Run(Progress), Run(Progress), Run(Progress)],
        &[(0, 7)],
        true,
        (1, 3)
    )]
    #[case::capped_without_permit(true, &[1], &[Hold, Run(Capped)], &[], true, (0, 0))]
    #[case::in_progress_takes_permit_directly(
        true,
        &[2],
        &[Run(Progress), Hold, Run(Capped), Release, Run(Progress)],
        &[(0, 7)],
        false,
        (1, 2)
    )]
    #[case::finished_only_when_idle(
        false,
        &[2],
        &[Close, Run(Progress), Run(Progress), Run(Finished)],
        &[(0, 7)],
        false,
        (1, 2)
    )]
    #[case::open_idle_is_no_progress(false, &[], &[Run(NoProgress)], &[], false, (0, 0))]
    fn fanout_worker_dispatch(
        #[case] capped: bool,
        #[case] jobs: &[u32],
        #[case] script: &[Act],
        #[case] delivered: &[(u64, u64)],
        #[case] input_left: bool,
        #[case] counters: (u64, u64),
    ) {
        let cap = capped.then(|| PhaseCap::new("c", 1));
        let inbox = OrderedInbox::<Out>::new(true);
        let mut w = FanoutWorker::<CountJob>::new("count", Arc::clone(&inbox), BUDGET)
            .with_step_cap(cap.clone().map(crate::StepCap::from));
        let probe = StepProbe::new(&w);
        for (i, &n) in jobs.iter().enumerate() {
            probe.push_input(CountJob::new(i as u64, n));
        }
        let mut held = None;
        for act in script {
            match *act {
                Run(expected) => {
                    assert_eq!(probe.try_run(&mut w).unwrap(), expected, "script {script:?}");
                }
                Hold => held = cap.as_ref().unwrap().try_acquire(),
                Release => held = None,
                Close => probe.close_input(),
            }
        }
        drop(held);
        let mut got = Vec::new();
        while let Some((o, Out(v))) = inbox.try_take_next() {
            got.push((o, v));
        }
        assert_eq!(got, delivered);
        assert_eq!(!probe.input_is_empty(), input_left);
        assert_eq!((probe.counter(0), probe.counter(1)), counters, "jobs, slices");
        assert_eq!(probe.counter(2), 8 * delivered.len() as u64, "output bytes");
        if let Some(cap) = cap {
            assert_eq!(cap.active(), 0, "no permit leaks");
        }
    }

    /// A job error reaches the pipeline with its ordinal, its kind and its
    /// text, which a chain-walking reporter prints once.
    #[rstest]
    fn error_propagates(#[values(true, false)] capped: bool) {
        let cap = capped.then(|| PhaseCap::new("c", 1));
        let mut w = make(cap.clone());
        let probe = StepProbe::new(&w);
        probe.push_input(CountJob { ordinal: 3, slices_left: 1, fail: true });
        let err = probe.try_run(&mut w).unwrap_err();
        assert_eq!(err.to_string(), "job 3: boom");
        assert_eq!(err.kind(), io::ErrorKind::Other);
        let job_error = err
            .get_ref()
            .and_then(|e| e.downcast_ref::<super::super::JobError>())
            .expect("a JobError");
        assert_eq!(job_error.ordinal(), 3);
        assert_eq!(job_error.job_error().to_string(), "boom");
        // A reporter walking the chain prints the job error's text once: the
        // display has it, and the chain continues at the job error's own
        // source (none here).
        let chain: Vec<String> =
            std::iter::successors(Some(&err as &dyn std::error::Error), |e| e.source())
                .map(ToString::to_string)
                .collect();
        assert_eq!(chain, vec!["job 3: boom".to_owned()]);
        if let Some(cap) = cap {
            assert_eq!(cap.active(), 0, "the permit is released after Err");
        }
    }

    #[test]
    fn new_worker_copy_starts_clean() {
        let mut w = make(None);
        let probe = StepProbe::new(&w);
        probe.push_input(CountJob::new(0, 2));
        assert_eq!(probe.try_run(&mut w).unwrap(), Progress);
        let mut copy = w.new_worker_copy();
        // The copy has no in-progress job: an open empty input is NoProgress.
        assert_eq!(probe.try_run(&mut copy).unwrap(), NoProgress);
    }

    #[test]
    fn exposes_inbox_as_side_channel() {
        let inbox = OrderedInbox::<Out>::new(true);
        let w = FanoutWorker::<CountJob>::new("count", Arc::clone(&inbox), BUDGET);
        let channels = w.side_channels();
        assert_eq!(channels.len(), 1);
        assert!(
            std::ptr::addr_eq(Arc::as_ptr(&channels[0]), Arc::as_ptr(&inbox)),
            "the side channel is the worker's inbox"
        );
        inbox.reserve(0, 64);
        assert_eq!(channels[0].in_flight_bytes(), 64, "the channel reads the inbox gauge");
    }
}
