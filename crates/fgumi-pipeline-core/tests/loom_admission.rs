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
//! Negative controls are revert-checks rather than in-suite tests, since both
//! wakes are production code with no switch: deleting the release-side
//! `notify_all` in `PhaseCap::release`, or the `cap.notify_cancel()` call in
//! `PipelineSignal::wake_parked_workers`, makes loom report a deadlock here.
//!
//! Exploration is preemption-bounded (3), as in `fgumi-sort`'s `loom_merge_slots`.
//!
//! Run: `RUSTFLAGS="--cfg loom" cargo test -p fgumi-pipeline-core --test loom_admission --release`
#![cfg(loom)]

use fgumi_pipeline_core::PhaseCap;
use fgumi_pipeline_core::signal::PipelineSignal;
use loom::sync::Arc;
use loom::sync::atomic::{AtomicUsize, Ordering};
use loom::thread;

const MAX: usize = 2;

/// Permits seen held at once, to check the cap is never exceeded.
struct Inside(AtomicUsize);

impl Inside {
    fn enter(&self, n: usize) {
        let now = self.0.fetch_add(n, Ordering::SeqCst) + n;
        assert!(now <= MAX, "{now} permits held in a cap of {MAX}");
    }

    fn leave(&self, n: usize) {
        self.0.fetch_sub(n, Ordering::SeqCst);
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
        let inside = Arc::new(Inside(AtomicUsize::new(0)));
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
