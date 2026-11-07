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

/// A serial owner posts sliced jobs to a [`FanoutWorker`] on the pool and takes
/// the results back in order, through `Pipeline::run` with the deadlock monitor
/// armed.
mod end_to_end {
    use std::io;
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::{Arc, OnceLock};
    use std::time::Duration;

    use parking_lot::Mutex;
    use rstest::rstest;

    use super::super::{FanoutWorker, Job, JobStep, OrderedInbox, SliceBudget};
    use crate::admission::{PhaseCap, StepCap, StepCapRef, admit_input};
    use crate::builder::{Pipeline, PipelineConfig};
    use crate::handles::Unpushed;
    use crate::held::HeldSlot;
    use crate::item::HeapSize;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;
    use crate::signal::{CancelHandle, PipelineError};
    use crate::step::{SideChannel, Step, StepCtx, StepKind, StepOutcome, StepProfile};

    /// Jobs the owner posts.
    const N: u64 = 2000;
    /// Jobs posted but not yet taken, at most.
    const WINDOW: u64 = 16;
    /// The job that fails in the `job_error` case.
    const FAIL_AT: u64 = 777;
    /// Results the sink consumes before the `cancel` case cancels the run.
    const CANCEL_AFTER: usize = 100;
    const WORKER: &str = "FanoutCount";
    const BUDGET: SliceBudget = SliceBudget { bytes: 1024, records: 16 };

    #[derive(Debug)]
    struct Out(u64);

    impl HeapSize for Out {
        fn heap_size(&self) -> usize {
            8
        }
    }

    /// A job of 1-4 slices that yields its own ordinal.
    struct CountJob {
        ordinal: u64,
        slices_left: u64,
        fail: bool,
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
                return Err(io::Error::other("job failed"));
            }
            self.slices_left -= 1;
            Ok(if self.slices_left == 0 {
                JobStep::Done(Out(self.ordinal))
            } else {
                JobStep::Continue
            })
        }
    }

    /// Posts `N` jobs on branch 1 within a window of [`WINDOW`], takes results in
    /// order from the inbox and emits each ordinal on branch 0.
    struct Owner {
        kind: StepKind,
        inbox: Arc<OrderedInbox<Out>>,
        fail_at: Option<u64>,
        bound: bool,
        posted: u64,
        taken: u64,
        held_ordinal: HeldSlot<Unpushed<u32>>,
        held_job: HeldSlot<Unpushed<CountJob>>,
    }

    impl Owner {
        fn new(kind: StepKind, inbox: Arc<OrderedInbox<Out>>, fail_at: Option<u64>) -> Self {
            Self {
                kind,
                inbox,
                fail_at,
                bound: false,
                posted: 0,
                taken: 0,
                held_ordinal: HeldSlot::new(),
                held_job: HeldSlot::new(),
            }
        }
    }

    impl Step for Owner {
        type Input = ();
        type Outputs = (u32, CountJob);

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "FanoutOwner",
                kind: self.kind,
                sticky: false,
                output_queues: vec![
                    QueueSpec::ByteBounded { limit_bytes: 1 << 10 },
                    QueueSpec::ByteBounded { limit_bytes: 1 << 10 },
                ],
                branch_ordering: vec![BranchOrdering::None, BranchOrdering::None],
            }
        }

        fn side_channels(&self) -> Vec<Arc<dyn SideChannel>> {
            vec![Arc::clone(&self.inbox) as Arc<dyn SideChannel>]
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if !self.bound {
                self.inbox.bind_owner(ctx.pool.waker());
                self.bound = true;
            }
            let out = ctx.outputs.view();
            let mut progressed = false;
            // Flush what the last dispatch could not push.
            if let Some(u) = self.held_ordinal.take() {
                match out.a.retry(u) {
                    Ok(()) => progressed = true,
                    Err(u) => self.held_ordinal.put(u),
                }
            }
            if let Some(u) = self.held_job.take() {
                match out.b.retry(u) {
                    Ok(()) => progressed = true,
                    Err(u) => self.held_job.put(u),
                }
            }
            // Take results in order and emit them.
            while !self.held_ordinal.is_held() {
                let Some((ordinal, result)) = self.inbox.try_take_next() else { break };
                assert_eq!(result.0, ordinal, "a result carries its own ordinal");
                self.taken += 1;
                progressed = true;
                if let Err(u) = out.a.push(u32::try_from(ordinal).expect("small ordinal")) {
                    self.held_ordinal.put(u);
                }
            }
            // Post one job per dispatch.
            if !self.held_job.is_held() && self.posted < N && self.posted - self.taken < WINDOW {
                let ordinal = self.posted;
                self.posted += 1;
                self.inbox.reserve(ordinal, 8);
                let job = CountJob {
                    ordinal,
                    slices_left: 1 + ordinal % 4,
                    fail: self.fail_at == Some(ordinal),
                };
                if let Err(u) = out.b.push(job) {
                    self.held_job.put(u);
                }
                progressed = true;
            }
            if progressed {
                return Ok(StepOutcome::Progress);
            }
            if self.held_ordinal.is_held() || self.held_job.is_held() {
                return Ok(StepOutcome::NoProgress);
            }
            if self.taken == N {
                return Ok(StepOutcome::Finished);
            }
            // Nothing to take, and the window is full or the jobs are all posted.
            Ok(if self.inbox.await_next() {
                StepOutcome::Progress
            } else {
                StepOutcome::NoProgress
            })
        }
    }

    /// Records each ordinal it sees, in arrival order; optionally cancels the
    /// run after [`CANCEL_AFTER`] of them.
    #[derive(Clone)]
    struct CountingSink {
        kind: StepKind,
        cap: Option<Arc<PhaseCap>>,
        seen: Arc<Mutex<Vec<u32>>>,
        cancel: Option<Arc<OnceLock<CancelHandle>>>,
    }

    impl Step for CountingSink {
        type Input = u32;
        type Outputs = ();

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "CountingSink",
                kind: self.kind,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }

        fn step_cap(&self) -> Option<StepCapRef<'_>> {
            self.cap.as_deref().map(StepCapRef::Fixed)
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let permit = match admit_input(ctx.input, self.cap.as_deref()) {
                Ok(p) => p,
                Err(outcome) => return Ok(outcome),
            };
            let Some(ordinal) = ctx.input.pop() else {
                return Ok(if ctx.input.is_drained() {
                    StepOutcome::Finished
                } else {
                    StepOutcome::NoProgress
                });
            };
            let len = {
                let mut seen = self.seen.lock();
                seen.push(ordinal);
                seen.len()
            };
            drop(permit);
            if len == CANCEL_AFTER
                && let Some(handle) = self.cancel.as_ref().and_then(|c| c.get())
            {
                handle.cancel();
            }
            Ok(StepOutcome::Progress)
        }

        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum Outcome {
        Complete,
        JobError,
        Cancel,
    }

    /// Fails the test (aborting the process, as a hung `run` cannot be
    /// interrupted) if the run does not return within 30 seconds.
    fn watchdog(name: &'static str) -> Arc<AtomicBool> {
        let done = Arc::new(AtomicBool::new(false));
        let flag = Arc::clone(&done);
        std::thread::spawn(move || {
            for _ in 0..600 {
                std::thread::sleep(Duration::from_millis(50));
                if flag.load(Ordering::SeqCst) {
                    return;
                }
            }
            eprintln!("{name}: WEDGED (deadlock)");
            std::process::abort();
        });
        done
    }

    /// Sums the per-tick deltas of counter `counter_name` of step `step_name` from the
    /// telemetry TSVs written at `stem`.
    fn counter_total(stem: &std::path::Path, step_name: &str, counter_name: &str) -> u64 {
        let read = |suffix: &str| {
            std::fs::read_to_string(crate::runtime::telemetry::ticks_path(stem, suffix)).unwrap()
        };
        let step_idx = read("steps")
            .lines()
            .skip(1)
            .find_map(|l| {
                l.split_once('\t').filter(|(_, n)| *n == step_name).map(|(i, _)| i.to_owned())
            })
            .expect("step row");
        let counter_idx = read("counter_names")
            .lines()
            .skip(1)
            .find_map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[0] == step_idx && f[2] == counter_name).then(|| f[1].to_owned())
            })
            .expect("counter row");
        read("counters")
            .lines()
            .skip(1)
            .filter_map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[2] == step_idx && f[3] == counter_idx).then(|| f[4].parse::<u64>().unwrap())
            })
            .sum()
    }

    #[rstest]
    #[case::uncapped(false, StepKind::Detached, Outcome::Complete)]
    #[case::cap_one(true, StepKind::Detached, Outcome::Complete)]
    #[case::job_error(false, StepKind::Detached, Outcome::JobError)]
    #[case::cancel(false, StepKind::Detached, Outcome::Cancel)]
    #[case::pool_hosted_owner(false, StepKind::Serial, Outcome::Complete)]
    fn fanout_end_to_end(
        #[case] capped: bool,
        #[case] owner_kind: StepKind,
        #[case] outcome: Outcome,
        #[values(1, 2, 4, 8)] threads: usize,
    ) {
        let inbox = OrderedInbox::<Out>::new(true);
        let cap = capped.then(|| PhaseCap::new("fanout-cap", 1));
        let seen = Arc::new(Mutex::new(Vec::<u32>::new()));
        let cancel = Arc::new(OnceLock::<CancelHandle>::new());

        let owner = Owner::new(
            owner_kind,
            Arc::clone(&inbox),
            (outcome == Outcome::JobError).then_some(FAIL_AT),
        );
        let worker = FanoutWorker::<CountJob>::new(WORKER, Arc::clone(&inbox), BUDGET)
            .with_step_cap(cap.clone().map(StepCap::from));
        // A parallel sink shares the cap; its permit covers pop and record, so
        // the single permit keeps its order.
        let sink = CountingSink {
            kind: if capped { StepKind::Parallel } else { StepKind::Serial },
            cap: cap.clone(),
            seen: Arc::clone(&seen),
            cancel: (outcome == Outcome::Cancel).then(|| Arc::clone(&cancel)),
        };

        let builder = Pipeline::builder();
        let multi = builder.chain(owner).into_multi();
        multi.b0.chain(sink).into_sink_marker();
        multi.b1.chain(worker).into_sink_marker();
        let pipeline = builder.build().expect("pipeline build");
        let _ = cancel.set(pipeline.cancel_handle());

        let dag = pipeline.dag_at(threads);
        let worker_line = dag
            .lines()
            .find(|l| l.contains(WORKER) && l.contains("branches="))
            .unwrap_or_else(|| panic!("worker line in {dag}"));
        assert_eq!(worker_line.contains("cap=fanout-cap(1)"), capped, "{worker_line}");
        // One worker is the whole pool: a detached owner names it, and a
        // pool-hosted owner shares its thread, so nothing needs waking.
        let target = match (threads, owner_kind) {
            (1, StepKind::Detached) => "Worker(0)",
            (1, _) => "None",
            _ => "Pool",
        };
        assert!(dag.contains(&format!(".1 → {WORKER} {target}")), "branch 1 wakes the pool: {dag}");

        let dir = std::env::temp_dir().join(format!(
            "fanout-e2e-{}-{threads}-{owner_kind:?}-{outcome:?}-{capped}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let stem = dir.join("run");
        let mut config = PipelineConfig { threads, ..Default::default() }.with_deadlock_timeout(10);
        if outcome == Outcome::Complete {
            config = config.with_telemetry(crate::runtime::telemetry::TelemetryConfig {
                stem: stem.clone(),
                interval: Duration::from_millis(2),
            });
        }

        let done = watchdog("fanout_end_to_end");
        let result = pipeline.run(config);
        done.store(true, Ordering::SeqCst);

        match outcome {
            Outcome::Complete => {
                result.expect("pipeline run");
                let expected: Vec<u32> = (0..u32::try_from(N).unwrap()).collect();
                assert_eq!(*seen.lock(), expected, "ordinals arrive in order, once each");
                assert_eq!(inbox.in_flight_bytes(), 0, "inbox gauge drained");
                assert_eq!(inbox.pending(), 0);
                assert_eq!(counter_total(&stem, WORKER, "jobs"), N, "jobs on the stats row");
            }
            Outcome::JobError => match result {
                Err(PipelineError::Io { step, source }) => {
                    assert_eq!(step, WORKER);
                    assert_eq!(source.to_string(), format!("job {FAIL_AT}: job failed"));
                }
                other => panic!("expected Io from {WORKER}, got {other:?}"),
            },
            Outcome::Cancel => {
                assert!(matches!(result, Err(PipelineError::Cancelled)), "{result:?}");
                assert!(seen.lock().len() >= CANCEL_AFTER);
            }
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
}

/// An owner waiting forever on an ordinal whose job never delivers is a wedge,
/// not starvation: the inbox reservation is the only in-flight charge (every
/// edge is empty), and the deadlock monitor must count it, fail the run, and
/// name the inbox in its snapshot.
mod wedged_owner {
    use std::io;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::time::Duration;

    use rstest::rstest;

    use super::super::OrderedInbox;
    use crate::builder::{Pipeline, PipelineConfig};
    use crate::item::HeapSize;
    use crate::outputs::Single;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;
    use crate::signal::PipelineError;
    use crate::step::{SideChannel, Step, StepCtx, StepKind, StepOutcome, StepProfile};

    /// The reservation the owner charges for its one job.
    const RESERVED: u64 = 64;

    #[derive(Debug)]
    struct Out;

    impl HeapSize for Out {
        fn heap_size(&self) -> usize {
            8
        }
    }

    /// A job descriptor; reports a non-zero size like a real job.
    #[derive(Debug)]
    struct LostJob;

    impl HeapSize for LostJob {
        fn heap_size(&self) -> usize {
            16
        }
    }

    /// Reserves ordinal 0, posts its job, then waits on the inbox forever.
    struct Owner {
        kind: StepKind,
        inbox: Arc<OrderedInbox<Out>>,
        posted: bool,
    }

    impl Step for Owner {
        type Input = ();
        type Outputs = Single<LostJob>;

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "WedgedOwner",
                kind: self.kind,
                sticky: false,
                output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 10 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }

        fn side_channels(&self) -> Vec<Arc<dyn SideChannel>> {
            vec![Arc::clone(&self.inbox) as Arc<dyn SideChannel>]
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if !self.posted {
                self.inbox.bind_owner(ctx.pool.waker());
                self.inbox.reserve(0, RESERVED);
                ctx.outputs.push(LostJob).map_err(|_| io::Error::other("job edge full"))?;
                self.posted = true;
                return Ok(StepOutcome::Progress);
            }
            Ok(if self.inbox.await_next() {
                StepOutcome::Progress
            } else {
                StepOutcome::NoProgress
            })
        }
    }

    /// Pops each job and drops it, so its result is never delivered.
    #[derive(Clone)]
    struct DroppingWorker;

    impl Step for DroppingWorker {
        type Input = LostJob;
        type Outputs = ();

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "DroppingWorker",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(match ctx.input.pop() {
                Some(LostJob) => StepOutcome::Progress,
                None if ctx.input.is_drained() => StepOutcome::Finished,
                None => StepOutcome::NoProgress,
            })
        }

        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    #[rstest]
    #[case::detached_owner(StepKind::Detached)]
    #[case::pool_hosted_owner(StepKind::Serial)]
    fn owner_wedged_on_inbox_is_a_deadlock(#[case] owner_kind: StepKind) {
        crate::test_support::install_log_capture();
        let inbox = OrderedInbox::<Out>::new(true);
        let builder = Pipeline::builder();
        builder
            .chain(Owner { kind: owner_kind, inbox: Arc::clone(&inbox), posted: false })
            .chain(DroppingWorker)
            .into_sink_marker();
        let pipeline = builder.build().expect("pipeline build");

        // A Starving verdict would reset the stall clock forever; bound the run
        // so that regression fails instead of hanging.
        let done = Arc::new(AtomicBool::new(false));
        let flag = Arc::clone(&done);
        std::thread::spawn(move || {
            for _ in 0..600 {
                std::thread::sleep(Duration::from_millis(50));
                if flag.load(Ordering::SeqCst) {
                    return;
                }
            }
            eprintln!("owner_wedged_on_inbox_is_a_deadlock: run never failed (Starving?)");
            std::process::abort();
        });
        // The shortest monitor: warn after 1 s, fail after 6 s.
        let result = pipeline
            .run(PipelineConfig { threads: 2, ..Default::default() }.with_deadlock_timeout(1));
        done.store(true, Ordering::SeqCst);

        assert!(matches!(result, Err(PipelineError::TimedOut { .. })), "{result:?}");
        let snapshot = crate::test_support::captured_logs(
            "fgumi_pipeline_core::builder",
            "side channel: inbox(",
        );
        assert!(!snapshot.is_empty(), "the stall snapshot names the inbox");
        assert!(
            snapshot.iter().any(|r| r.message.contains(&format!("({RESERVED} bytes in flight)"))),
            "the snapshot charges the reservation: {snapshot:?}"
        );
    }
}
