//! Loom model of the whole-cap reservation, driving the REAL `PhaseCap`
//! (`try_acquire`, `acquire_whole`, permit and whole-cap release) and the real
//! cancel path (`PipelineSignal::cancel` waking the cap it was bound to).
//!
//! Under `--cfg loom` the cap's atomics, fences and parking lock/condvar are
//! loom's (see the aliases at the top of `admission.rs`), so loom explores the
//! fence pair in `release`/`acquire_whole` and the untimed wait exactly as
//! production runs them. The run's signal state stays a std atomic: loom does
//! not reorder it, which is sound here because the cancel's wake is delivered
//! under the cap's (loom) lock that the waiter re-checks the signal under.
//!
//! It also drives the cap-release wake of a refused pool worker: the refusal
//! records the worker's wake slot, fences and re-checks; the release fences
//! and claims the slot, and the bound waker unparks it (`loom::thread::park`
//! stands in for the worker's timer park, with no timeout). The waker here
//! unparks unconditionally: these models check the cap's own record, claim and
//! forward protocol against an ideal delivery. `Pipeline::run` binds
//! `WakePlan::deliver_slot`, which unparks only a thread armed for its timer
//! park; that composition — the worker's arm and re-poll pass against a
//! release, and a skipped wake's forward — is model-checked in
//! `tests/loom_wake.rs` (`cap_release_gated_model`,
//! `forwarded_cap_wake_gated_model`).
//!
//! And the owner side of the per-pass record: the pass-start clear
//! (`PhaseCap::clear_own_record`) raced against a release; the forward of a
//! wake a thread was handed and did not spend (`PhaseCap::forward_wake`); and
//! a cap-parked re-poll's skip of a full cap (`PhaseCap::record_skip`) raced
//! against the release that frees it.
//!
//! Negative controls are revert-checks rather than in-suite tests, since every
//! wake is production code with no switch. Each of these makes loom report a
//! deadlock here:
//! - deleting the release-side `notify_all` in `PhaseCap::release`;
//! - deleting the `cap.notify_cancel()` call in
//!   `PipelineSignal::wake_parked_workers`;
//! - deleting the refusal's fence or re-check, or the release's `wake_refused`;
//! - deleting `forward_wake`'s claim (`an_unspent_wake_is_forwarded_to_the_refused_worker`);
//! - deleting `record_skip`'s re-check (`a_full_cap_skip_racing_a_release_loses_no_wake`);
//! - reporting a claimed skip record as not noted
//!   (`a_skip_whose_fresh_record_a_release_claims_forwards_the_wake`).
//!
//! `a_pass_start_clear_racing_a_release_loses_no_wake` checks safety only (no
//! lost wake when a worker clears and then polls the cap); it is not a
//! revert-check of its own, since it passes with or without the clear.
//!
//! And the real `FanoutWorker` (through a `StepProbe`): two clones slicing
//! jobs under a one-permit cap, one permit per slice, with the cap's release
//! wake (`fanout_cap_slices_model`).
//!
//! Exploration is preemption-bounded (3), as in `fgumi-sort`'s `loom_merge_slots`.
//!
//! Run: `RUSTFLAGS="--cfg loom" cargo test -p fgumi-pipeline-core --test loom_admission --release`
#![cfg(loom)]

use fgumi_pipeline_core::signal::PipelineSignal;
use fgumi_pipeline_core::testing::StepProbe;
use fgumi_pipeline_core::{
    FanoutWorker, HeapSize, Job, JobStep, OrderedInbox, PhaseCap, SliceBudget, Step, StepCap,
    StepOutcome,
};
use loom::sync::Arc;
use loom::sync::atomic::{AtomicUsize, Ordering};
use loom::thread;

const MAX: usize = 2;

/// Permits seen held at once, to check a cap of `max` is never exceeded.
struct Inside {
    held: AtomicUsize,
    max: usize,
}

impl Inside {
    fn new(max: usize) -> Self {
        Self { held: AtomicUsize::new(0), max }
    }

    fn enter(&self, n: usize) {
        let now = self.held.fetch_add(n, Ordering::SeqCst) + n;
        assert!(now <= self.max, "{now} permits held in a cap of {}", self.max);
    }

    fn leave(&self, n: usize) {
        self.held.fetch_sub(n, Ordering::SeqCst);
    }
}

/// A whole-cap caller against two pool admitters (each one item) and,
/// optionally, a run cancel: the caller always finishes — it takes the whole
/// cap or, when cancelled, gives up — the cap is never exceeded, and nothing
/// is left held or reserved.
fn model(cancel: bool) {
    // 3 preemptions cover the lost-wake interleavings (holder release vs
    // reservation vs CAS, cancel vs park) while keeping the four-thread model
    // tractable.
    let mut builder = loom::model::Builder::new();
    builder.preemption_bound = Some(3);
    builder.check(move || {
        let cap = PhaseCap::new("loom", MAX);
        let signal = PipelineSignal::new();
        cap.bind_signal_for_test(&signal);
        let inside = Arc::new(Inside::new(MAX));
        let admitters: Vec<_> = (0..2)
            .map(|_| {
                let cap = std::sync::Arc::clone(&cap);
                let inside = Arc::clone(&inside);
                thread::spawn(move || {
                    if let Some(permit) = cap.try_acquire() {
                        inside.enter(1);
                        inside.leave(1);
                        drop(permit);
                    }
                })
            })
            .collect();
        let canceller = cancel.then(|| {
            let signal = std::sync::Arc::clone(&signal);
            thread::spawn(move || signal.cancel())
        });
        if let Some(whole) = cap.acquire_whole() {
            assert_eq!(whole.width(), MAX);
            inside.enter(MAX);
            inside.leave(MAX);
            drop(whole);
        }
        for a in admitters {
            a.join().unwrap();
        }
        if let Some(c) = canceller {
            c.join().unwrap();
        }
        assert_eq!(cap.active(), 0, "every permit released");
        assert_eq!(cap.pending_reservations(), 0, "no reservation left behind");
    });
}

#[test]
fn whole_cap_waiter_finishes_and_never_exceeds_the_cap() {
    model(false);
}

#[test]
fn whole_cap_waiter_gives_up_when_cancelled() {
    model(true);
}

/// A pool holder that keeps its permit until the whole-cap caller has
/// returned: no release can wake the caller, so only the run's cancel can.
/// The caller must give up (`None`) and leave no reservation; without the
/// cancel's wake it parks forever and loom reports the deadlock.
#[test]
fn a_cancel_wakes_a_waiter_no_release_will_wake() {
    let mut builder = loom::model::Builder::new();
    builder.preemption_bound = Some(3);
    builder.check(|| {
        let cap = PhaseCap::new("loom", MAX);
        let signal = PipelineSignal::new();
        cap.bind_signal_for_test(&signal);
        let held = cap.try_acquire().expect("an idle cap admits the holder");
        let waiter = {
            let cap = std::sync::Arc::clone(&cap);
            thread::spawn(move || cap.acquire_whole().is_some())
        };
        let canceller = {
            let signal = std::sync::Arc::clone(&signal);
            thread::spawn(move || signal.cancel())
        };
        assert!(!waiter.join().unwrap(), "cancelled while a holder remains: no permits");
        canceller.join().unwrap();
        drop(held);
        assert_eq!(cap.active(), 0, "every permit released");
        assert_eq!(cap.pending_reservations(), 0, "no reservation left behind");
    });
}

/// A waker that unparks the loom thread registered for wake slot 0.
fn slot0_waker() -> (Arc<loom::sync::Mutex<Option<thread::Thread>>>, fgumi_pipeline_core::CapWaker)
{
    let handle: Arc<loom::sync::Mutex<Option<thread::Thread>>> =
        Arc::new(loom::sync::Mutex::new(None));
    let h = Arc::clone(&handle);
    let waker: fgumi_pipeline_core::CapWaker = std::sync::Arc::new(move |slot: Option<usize>| {
        if slot == Some(0)
            && let Some(t) = h.lock().unwrap().as_ref()
        {
            t.unpark();
        }
    });
    (handle, waker)
}

/// A pool worker in wake slot 0: registers its handle, then takes a permit,
/// parking (untimed) after each refusal until a release wakes it.
fn refused_worker(
    cap: &std::sync::Arc<PhaseCap>,
    handle: &Arc<loom::sync::Mutex<Option<thread::Thread>>>,
) -> thread::JoinHandle<()> {
    let cap = std::sync::Arc::clone(cap);
    let handle = Arc::clone(handle);
    thread::spawn(move || {
        *handle.lock().unwrap() = Some(thread::current());
        loop {
            if let Some(permit) = cap.try_acquire_as(0) {
                drop(permit);
                return;
            }
            thread::park();
        }
    })
}

/// A refused worker is woken by the release of the permit that refused it: it
/// either sees the freed permit on its post-fence re-check, or the release sees
/// its recorded slot and unparks it. No timer exists here, so a lost wake is a
/// loom deadlock.
#[test]
fn a_release_wakes_the_refused_worker() {
    let mut builder = loom::model::Builder::new();
    builder.preemption_bound = Some(3);
    builder.check(|| {
        let cap = PhaseCap::new("loom", 1);
        let (handle, waker) = slot0_waker();
        cap.bind_waker_for_test(waker);
        let held = cap.try_acquire_as(1).expect("an idle cap admits the holder");
        let worker = refused_worker(&cap, &handle);
        drop(held);
        worker.join().unwrap();
        assert_eq!(cap.active(), 0, "every permit released");
    });
}

/// A worker refused while a whole-cap caller reserves the cap is woken by the
/// whole-cap release, not by the pool release that went to the reservation:
/// the release that empties the cap into a reservation wakes no refused worker
/// (it would only be refused again), and the whole holder's release does.
/// The worker is refused (and parked) before the pool holder releases, and
/// the whole holder releases only after that, so exactly one wake is right: a
/// reservation release that also woke the worker would see it refused again
/// and woken a second time.
#[test]
fn a_whole_cap_release_wakes_a_worker_refused_during_the_reservation() {
    let mut builder = loom::model::Builder::new();
    builder.preemption_bound = Some(3);
    builder.check(|| {
        let cap = PhaseCap::new("loom", 1);
        let signal = PipelineSignal::new();
        cap.bind_signal_for_test(&signal);
        let (handle, waker) = slot0_waker();
        let wakes = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let counted: fgumi_pipeline_core::CapWaker = {
            let wakes = std::sync::Arc::clone(&wakes);
            std::sync::Arc::new(move |slot: Option<usize>| {
                if slot == Some(0) {
                    wakes.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
                waker(slot);
            })
        };
        cap.bind_waker_for_test(counted);
        let held = cap.try_acquire_as(1).expect("an idle cap admits the holder");
        let whole = {
            let cap = std::sync::Arc::clone(&cap);
            thread::spawn(move || {
                let whole = cap.acquire_whole().expect("not cancelled");
                while cap.refused() == 0 {
                    thread::yield_now();
                }
                drop(whole);
            })
        };
        while cap.pending_reservations() == 0 {
            thread::yield_now();
        }
        let worker = refused_worker(&cap, &handle);
        while cap.refused() == 0 {
            thread::yield_now();
        }
        drop(held);
        whole.join().unwrap();
        worker.join().unwrap();
        assert_eq!(cap.active(), 0, "every permit released");
        assert_eq!(cap.pending_reservations(), 0, "no reservation left behind");
        assert_eq!(
            wakes.load(std::sync::atomic::Ordering::Relaxed),
            1,
            "only the whole holder's release wakes the refused worker"
        );
    });
}

/// A waker that unparks the loom thread registered for each wake slot.
type Handles = Arc<loom::sync::Mutex<[Option<thread::Thread>; 2]>>;

fn slots_waker() -> (Handles, fgumi_pipeline_core::CapWaker) {
    let handles: Handles = Arc::new(loom::sync::Mutex::new([None, None]));
    let h = Arc::clone(&handles);
    let waker: fgumi_pipeline_core::CapWaker = std::sync::Arc::new(move |slot: Option<usize>| {
        if let Some(t) = slot.and_then(|s| h.lock().unwrap().get(s).cloned().flatten()) {
            t.unpark();
        }
    });
    (handles, waker)
}

/// The pass-start clear raced against a release. Slot 0 left a record in an
/// earlier pass and starts a new one: it clears its own bit
/// (`clear_own_record`), then polls the cap, parking (untimed) after each
/// refusal. Slot 1 is refused and parks the same way. The holder releases the
/// one permit concurrently. Whether the release claims slot 0's stale bit
/// before the clear (a wake for a thread that is running: its token ends its
/// next park at once, and it takes the permit, whose release then wakes slot
/// 1), or the clear comes first (the release claims slot 1), both workers get
/// the permit in turn. A lost wake is a loom deadlock.
#[test]
fn a_pass_start_clear_racing_a_release_loses_no_wake() {
    let mut builder = loom::model::Builder::new();
    builder.preemption_bound = Some(3);
    builder.check(|| {
        let cap = PhaseCap::new("loom", 1);
        let (handles, waker) = slots_waker();
        cap.bind_waker_for_test(waker);
        let held = cap.try_acquire_as(9).expect("an idle cap admits the holder");
        assert!(cap.try_acquire_as(0).is_none(), "slot 0's record from an earlier pass");
        let worker = |slot: usize, new_pass: bool| {
            let cap = std::sync::Arc::clone(&cap);
            let handles = Arc::clone(&handles);
            thread::spawn(move || {
                handles.lock().unwrap()[slot] = Some(thread::current());
                if new_pass {
                    cap.clear_own_record_as(slot);
                }
                loop {
                    if let Some(permit) = cap.try_acquire_as(slot) {
                        drop(permit);
                        return;
                    }
                    thread::park();
                }
            })
        };
        let zero = worker(0, true);
        let one = worker(1, false);
        drop(held);
        zero.join().unwrap();
        one.join().unwrap();
        assert_eq!(cap.active(), 0, "every permit released");
    });
}

/// A wake handed to a thread with no use for it is forwarded. Slot 0 holds a
/// record from an earlier pass, slot 1 is refused and parks (untimed), and the
/// holder releases the one permit. Slot 0 starts its next pass
/// (`clear_own_record`): if its bit is already clear, a release claimed it —
/// it was handed that wake — and since this pass takes no permit, it forwards
/// the wake (`forward_wake`). Whichever of slot 0 and slot 1 the release
/// claims, slot 1 gets the permit. Without the forward, the release that
/// claims slot 0's stale bit leaves slot 1 parked: a loom deadlock.
#[test]
fn an_unspent_wake_is_forwarded_to_the_refused_worker() {
    let mut builder = loom::model::Builder::new();
    builder.preemption_bound = Some(3);
    builder.check(|| {
        let cap = PhaseCap::new("loom", 1);
        let (handles, waker) = slots_waker();
        cap.bind_waker_for_test(waker);
        let held = cap.try_acquire_as(9).expect("an idle cap admits the holder");
        assert!(cap.try_acquire_as(0).is_none(), "slot 0's record from an earlier pass");
        let zero = {
            let cap = std::sync::Arc::clone(&cap);
            thread::spawn(move || {
                if !cap.clear_own_record_as(0) {
                    cap.forward_wake_for_test();
                }
            })
        };
        let one = {
            let cap = std::sync::Arc::clone(&cap);
            let handles = Arc::clone(&handles);
            thread::spawn(move || {
                handles.lock().unwrap()[1] = Some(thread::current());
                loop {
                    if let Some(permit) = cap.try_acquire_as(1) {
                        drop(permit);
                        return;
                    }
                    thread::park();
                }
            })
        };
        drop(held);
        zero.join().unwrap();
        one.join().unwrap();
        assert_eq!(cap.active(), 0, "every permit released");
    });
}

/// A cap-parked re-poll's skip of a full cap, raced against the release that
/// frees it: the skip records the worker, fences and re-checks
/// (`record_skip`), so either the release claims its bit and wakes it, or the
/// re-check sees the free permit and the worker polls the step. The worker
/// parks (untimed) after each skip, so a lost wake is a loom deadlock.
#[test]
fn a_full_cap_skip_racing_a_release_loses_no_wake() {
    let mut builder = loom::model::Builder::new();
    builder.preemption_bound = Some(3);
    builder.check(|| {
        let cap = PhaseCap::new("loom", 1);
        let (handles, waker) = slots_waker();
        cap.bind_waker_for_test(waker);
        let held = cap.try_acquire_as(9).expect("an idle cap admits the holder");
        let worker = {
            let cap = std::sync::Arc::clone(&cap);
            let handles = Arc::clone(&handles);
            thread::spawn(move || {
                handles.lock().unwrap()[0] = Some(thread::current());
                loop {
                    if cap.record_skip_as(0).skip {
                        thread::park();
                    } else if let Some(permit) = cap.try_acquire_as(0) {
                        drop(permit);
                        return;
                    }
                }
            })
        };
        drop(held);
        worker.join().unwrap();
        assert_eq!(cap.active(), 0, "every permit released");
    });
}

/// A skip whose own bit a release claims before its re-check is still a
/// wake handed to the skipping thread. Slot 1 is refused and parks
/// (untimed). Slot 0 skips the full cap (`record_skip`), which records,
/// fences and re-checks; the holder releases concurrently. If the release
/// claims slot 0's fresh bit, slot 0's re-check sees the freed permit and it
/// polls instead of skipping — so the skip reports the record as noted.
/// Slot 0's poll takes nothing; its next pass start (`clear_own_record`)
/// finds the noted bit gone and forwards the wake. Without the note, slot 1
/// stays parked while the permit is free: a loom deadlock.
#[test]
fn a_skip_whose_fresh_record_a_release_claims_forwards_the_wake() {
    let mut builder = loom::model::Builder::new();
    builder.preemption_bound = Some(3);
    builder.check(|| {
        let cap = PhaseCap::new("loom", 1);
        let (handles, waker) = slots_waker();
        cap.bind_waker_for_test(waker);
        let held = cap.try_acquire_as(9).expect("an idle cap admits the holder");
        assert!(cap.try_acquire_as(1).is_none(), "slot 1 is refused");
        let zero = {
            let cap = std::sync::Arc::clone(&cap);
            thread::spawn(move || {
                let noted = cap.record_skip_as(0).noted;
                let still = cap.clear_own_record_as(0);
                if noted && !still {
                    cap.forward_wake_for_test();
                }
            })
        };
        let one = {
            let cap = std::sync::Arc::clone(&cap);
            let handles = Arc::clone(&handles);
            thread::spawn(move || {
                handles.lock().unwrap()[1] = Some(thread::current());
                loop {
                    if let Some(permit) = cap.try_acquire_as(1) {
                        drop(permit);
                        return;
                    }
                    thread::park();
                }
            })
        };
        drop(held);
        zero.join().unwrap();
        one.join().unwrap();
        assert_eq!(cap.active(), 0, "every permit released");
    });
}

/// A job of `slices` slices; each slice runs inside the `Inside` guard and
/// counts itself in `ran`.
struct SliceJob {
    ordinal: u64,
    left: usize,
    inside: Arc<Inside>,
    ran: Arc<AtomicUsize>,
}

impl HeapSize for SliceJob {
    fn heap_size(&self) -> usize {
        64
    }
}

/// A job's result.
struct SliceOut;

impl HeapSize for SliceOut {}

impl Job for SliceJob {
    type Output = SliceOut;

    fn ordinal(&self) -> u64 {
        self.ordinal
    }

    fn run_slice(&mut self, _budget: SliceBudget) -> std::io::Result<JobStep<SliceOut>> {
        self.inside.enter(1);
        self.ran.fetch_add(1, Ordering::SeqCst);
        self.inside.leave(1);
        self.left -= 1;
        Ok(if self.left == 0 { JobStep::Done(SliceOut) } else { JobStep::Continue })
    }
}

/// A waker for threads with no wake slot (`current_slot` is `None` on a loom
/// thread, so a refusal records the cap's anonymous bit): a `None` wake
/// unparks every registered thread, as `WakePlan::deliver_anonymous` does.
fn broadcast_waker() -> (Handles, fgumi_pipeline_core::CapWaker) {
    let handles: Handles = Arc::new(loom::sync::Mutex::new([None, None]));
    let h = Arc::clone(&handles);
    let waker: fgumi_pipeline_core::CapWaker = std::sync::Arc::new(move |slot: Option<usize>| {
        let threads = h.lock().unwrap();
        match slot {
            Some(s) => threads.get(s).cloned().flatten().into_iter().for_each(|t| t.unpark()),
            None => threads.iter().flatten().for_each(thread::Thread::unpark),
        }
    });
    (handles, waker)
}

/// Slices per job.
const SLICES_PER_JOB: usize = 3;

/// Jobs on the edge: one seated on each clone before the race, one left on
/// the edge.
const JOBS: u64 = 3;

/// Two clones of the REAL `FanoutWorker` sharing one input edge (a
/// `StepProbe`) under a one-permit cap with its release wake bound. The edge
/// holds three jobs of three slices each and is drained. Before the race each
/// clone runs one dispatch on the model's main thread, which pops a job and
/// runs its first slice, so each enters the race with a job in progress; the
/// third job stays on the edge. Then each clone loops on `try_run` the way the
/// worker loop dispatches it: `Progress` runs again, `Capped` parks (a loom
/// park, no timer) until a release wakes it, `Finished` returns. A clone with
/// a job in progress takes a permit directly (`try_acquire`); one without
/// admits through `admit_input` and pops; either releases its permit when the
/// slice returns. Checks: never two permits held (the `Inside` guard in every
/// slice), every slice runs, every result is delivered to the inbox, the cap
/// ends empty. A release wake lost to a refused clone is a loom deadlock.
///
/// The jobs are seated before the race because the edge is a
/// `crossbeam_queue` (`UnboundedQueue`), which loom does not instrument: its
/// partial-order reduction treats two pops as independent, so it would never
/// schedule the second clone's pop while the first holds a job. What it does
/// explore is every interleaving of the cap's (loom) protocol: per-slice
/// `try_acquire` with its refusal record, fence and re-check, the release's
/// fence and claim, and the wake.
#[test]
fn fanout_cap_slices_model() {
    let mut builder = loom::model::Builder::new();
    builder.preemption_bound = Some(3);
    builder.check(|| {
        let cap = PhaseCap::new("fanout", 1);
        let (handles, waker) = broadcast_waker();
        cap.bind_waker_for_test(waker);
        let inside = Arc::new(Inside::new(1));
        let ran = Arc::new(AtomicUsize::new(0));
        let inbox = OrderedInbox::<SliceOut>::new(true);
        let budget = SliceBudget { bytes: 1, records: 1 };
        let worker = FanoutWorker::<SliceJob>::new("fanout", std::sync::Arc::clone(&inbox), budget)
            .with_step_cap(StepCap::from(std::sync::Arc::clone(&cap)));
        let probe = std::sync::Arc::new(StepProbe::new(&worker));
        for ordinal in 0..JOBS {
            probe.push_input(SliceJob {
                ordinal,
                left: SLICES_PER_JOB,
                inside: Arc::clone(&inside),
                ran: Arc::clone(&ran),
            });
        }
        probe.close_input();
        // Seat both clones before either thread starts.
        let seated: Vec<_> = (0..2)
            .map(|_| {
                let mut clone = worker.new_worker_copy();
                let first = probe.try_run(&mut clone).expect("a slice never fails here");
                assert_eq!(first, StepOutcome::Progress, "the clone pops a job");
                clone
            })
            .collect();
        let clones: Vec<_> = seated
            .into_iter()
            .enumerate()
            .map(|(i, mut clone)| {
                let probe = std::sync::Arc::clone(&probe);
                let handles = Arc::clone(&handles);
                thread::spawn(move || {
                    handles.lock().unwrap()[i] = Some(thread::current());
                    loop {
                        match probe.try_run(&mut clone).expect("a slice never fails here") {
                            StepOutcome::Progress => {}
                            StepOutcome::Capped => thread::park(),
                            StepOutcome::Finished => return,
                            other => panic!("unexpected outcome on a drained edge: {other:?}"),
                        }
                    }
                })
            })
            .collect();
        for c in clones {
            c.join().unwrap();
        }
        let jobs = usize::try_from(JOBS).expect("small");
        assert_eq!(ran.load(Ordering::SeqCst), jobs * SLICES_PER_JOB, "every slice runs");
        assert_eq!(inbox.pending(), jobs, "every result delivered");
        assert_eq!(cap.active(), 0, "every permit released");
    });
}
