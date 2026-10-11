//! `WorkerCore`: per-thread state carried by the worker loop.

use std::time::Duration;

use crate::topology::StepIdx;

// Timer-idle bounds for pool workers that do not wait on the event-count: a
// pinned worker (Exclusive owner, sticky owner, or the affinity target of a
// Serial step), the lone worker of a one-worker run, and a worker holding an
// item. Under Directed wakes (`runtime::wake`) a producer can `unpark` it
// early; in Legacy mode nothing holds its handle and the timer is the whole
// policy.
//
// The initial (also the post-progress reset floor) is 20µs, NOT 1µs. A surplus
// worker in an oversubscribed pool makes the occasional lucky pop, which resets
// its backoff; a 1µs floor then puts it back on the steep part of the ramp,
// waking ~every microsecond to run a full empty poll pass. 20µs is well below
// any step's per-item service time yet coarse enough to stop that re-poll.
const SLEEP_INITIAL_US: u64 = 20;
const SLEEP_MAX_US: u64 = 50_000; // 50 milliseconds

// Dedicated-driver idle bounds. A driver drives a small step subset off the
// pool; under Directed wakes every producer of its live steps unparks it, so
// the deadline is the self-heal bound for a wake nobody thought of, not the
// wake mechanism. 2 ms bounds the latency of a wake nobody delivers while keeping
// idle timer expiries rare; `timed out` on the detached line of
// `--pipeline-stats` is the counter that says whether a driver is being fed by
// this timer.
const PARK_INITIAL_US: u64 = 10;
const PARK_MAX_US: u64 = 2_000;

/// Which threads a [`TestBackoff`] applies to. Test support (`test-utils`);
/// not part of the supported API.
#[cfg(any(test, feature = "test-utils"))]
#[doc(hidden)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TestBackoffTarget {
    /// Pool worker `w`.
    Worker(usize),
    /// Every pool worker.
    AllWorkers,
    /// Detached driver `d`.
    Driver(usize),
    /// Every detached driver.
    AllDrivers,
}

/// Replace both idle-backoff bounds of the matching threads with `us`
/// microseconds. A liveness test gives exactly the thread whose wake it checks
/// a 10 s timer, so only that wake can finish the run under its watchdog.
/// Test support (`test-utils`); not part of the supported API.
#[cfg(any(test, feature = "test-utils"))]
#[doc(hidden)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TestBackoff {
    /// The threads this entry applies to.
    pub target: TestBackoffTarget,
    /// The pinned backoff, in microseconds.
    pub us: u64,
}

/// The override for `thread` among `entries`, if any entry covers it (last
/// match wins). Test support (`test-utils`); not part of the supported API.
#[cfg(any(test, feature = "test-utils"))]
#[doc(hidden)]
#[must_use]
pub fn backoff_override_for(entries: &[TestBackoff], thread: TestBackoffTarget) -> Option<u64> {
    entries.iter().rev().find(|b| b.target.covers(thread)).map(|b| b.us)
}

#[cfg(any(test, feature = "test-utils"))]
impl TestBackoffTarget {
    /// Whether an entry with this target applies to `thread`.
    #[must_use]
    pub fn covers(self, thread: Self) -> bool {
        match (self, thread) {
            (Self::AllWorkers, Self::Worker(_)) | (Self::AllDrivers, Self::Driver(_)) => true,
            (a, b) => a == b,
        }
    }
}

/// How a thread idles on a no-progress tick, by role: a pool worker idles with
/// [`ParkedSleep`](BackoffPolicy::ParkedSleep) (or blocks on the event-count,
/// using the ramp only as its wait deadline, when it is unpinned and holds
/// nothing); a dedicated driver thread with [`Park`](BackoffPolicy::Park).
/// Both arms park: the only difference is the ramp.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BackoffPolicy {
    /// `thread::park_timeout`, ramp 20µs→50ms. A registered thread can be
    /// unparked early; with no unparker this is a bounded sleep.
    ParkedSleep,
    /// `thread::park_timeout`, ramp 10µs→2ms. For dedicated driver threads.
    Park,
}

impl BackoffPolicy {
    #[inline]
    fn initial_us(self) -> u64 {
        match self {
            Self::ParkedSleep => SLEEP_INITIAL_US,
            Self::Park => PARK_INITIAL_US,
        }
    }

    #[inline]
    fn max_us(self) -> u64 {
        match self {
            Self::ParkedSleep => SLEEP_MAX_US,
            Self::Park => PARK_MAX_US,
        }
    }
}

/// Whether a `run_worker_loop` thread is an N-pool worker or a dedicated driver
/// (the unified "1-thread pool" for a set of off-pool steps). Controls only how
/// the thread's busy/idle time is attributed in `--pipeline-stats`, always
/// excluded from the pool% so the "N + 2" split stays visible:
/// - pool threads sum their whole-pass busy/idle into the N-worker utilisation
///   line, by `thread_id`;
/// - driver threads record each grouped step's own busy on the off-pool
///   "detached" line (by that step's index — so a multi-step `Shared` group
///   shows each member's real time, not the whole thread's under one name), and
///   attribute the thread-level idle/park to the group's `primary_step`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkerRole {
    /// N-worker pool thread; busy/idle keyed by `thread_id`.
    Pool,
    /// Dedicated driver thread. Per-step busy is recorded by each step's own
    /// index in `dispatch_one_step`; thread-level idle/park is keyed to
    /// `primary_step` (the group's representative).
    Driver { primary_step: StepIdx },
}

pub struct WorkerCore {
    /// `0..n_workers` for pool threads; unused (0) for driver threads.
    pub thread_id: usize,
    /// If this worker is the sole eligible dispatcher for a `sticky` step
    /// (either an `Exclusive sticky` step it owns, or a `Serial + sticky`
    /// step whose `Affinity` targets this worker), the step's index.
    /// The driver drives this step in a tight inner loop until it returns
    /// `NoProgress` / `Contention` / `Finished`, then yields to round-
    /// robin. Mirrors the legacy pipeline's sticky read.
    pub sticky_owner: Option<StepIdx>,
    /// Pool worker vs dedicated driver: how the thread's aggregate busy/idle
    /// time is attributed in `--pipeline-stats`, and which ramp it idles on
    /// (its `BackoffPolicy`).
    role: WorkerRole,
    /// `(initial, max)` bounds of the ramp: the role's policy's, or a test
    /// override.
    bounds: (u64, u64),
    /// Backoff duration in microseconds. Doubled on no-progress; reset on progress.
    backoff_us: u64,
}

impl WorkerCore {
    /// A pool worker (`ParkedSleep` policy).
    #[must_use]
    pub fn new(thread_id: usize, sticky_owner: Option<StepIdx>) -> Self {
        Self::with_role(thread_id, sticky_owner, WorkerRole::Pool)
    }

    /// A dedicated driver thread (the unified "1-thread pool"): `Driver` role +
    /// `Park` backoff. Drives a set of `Owned` steps off the pool; its
    /// aggregate busy/idle is attributed on the off-pool detached line.
    /// `thread_id` is unused for drivers; it never owns Exclusive or sticky
    /// steps.
    #[must_use]
    pub fn driver(primary_step: StepIdx) -> Self {
        Self::with_role(0, None, WorkerRole::Driver { primary_step })
    }

    fn with_role(thread_id: usize, sticky_owner: Option<StepIdx>, role: WorkerRole) -> Self {
        let policy = Self::policy_of(role);
        let bounds = (policy.initial_us(), policy.max_us());
        Self { thread_id, sticky_owner, role, bounds, backoff_us: bounds.0 }
    }

    fn policy_of(role: WorkerRole) -> BackoffPolicy {
        match role {
            WorkerRole::Pool => BackoffPolicy::ParkedSleep,
            WorkerRole::Driver { .. } => BackoffPolicy::Park,
        }
    }

    /// Replace both ramp bounds with `us` (a `TestBackoff` entry); `None`
    /// keeps the policy's.
    #[must_use]
    pub fn with_backoff_override(mut self, us: Option<u64>) -> Self {
        if let Some(us) = us {
            self.bounds = (us, us);
            self.backoff_us = us;
        }
        self
    }

    /// This thread's pool/driver role (drives stats attribution in the loop).
    #[must_use]
    pub fn role(&self) -> WorkerRole {
        self.role
    }

    /// This thread's idle policy (by role).
    #[cfg(test)]
    #[must_use]
    pub fn policy(&self) -> BackoffPolicy {
        Self::policy_of(self.role)
    }

    /// Return the ramp to its floor after progress.
    pub fn reset_backoff(&mut self) {
        self.backoff_us = self.bounds.0;
    }

    /// The ramp's cap as a `Duration`: the park of a thread whose wake is a
    /// direct `unpark` that is certain to come (a cap release), for which the
    /// timer is only the self-heal and ramping up from the floor would just
    /// poll.
    #[must_use]
    pub fn max_backoff_deadline(&self) -> Duration {
        Duration::from_micros(self.bounds.1)
    }

    /// Double the backoff, up to the ramp's cap.
    pub fn increase_backoff(&mut self) {
        self.backoff_us = self.backoff_us.saturating_mul(2).min(self.bounds.1);
    }

    /// The current backoff as a `Duration`, used as the event-count `wait`
    /// timeout (the self-heal bound on a missed wakeup) and as the deadline a
    /// timer park is classified against: the worker loop's timer park is
    /// `park_timeout` of this (or of [`Self::max_backoff_deadline`]). The park
    /// token is the thread's std token, which std's own blocking primitives
    /// (`mpsc::recv`, `Mutex`) also use: a token consumed inside a step that
    /// blocks on one of those is not seen there, and that one wake degrades to
    /// the timer.
    #[must_use]
    pub fn backoff_deadline(&self) -> Duration {
        Duration::from_micros(self.backoff_us)
    }

    #[cfg(test)]
    fn current_backoff_us(&self) -> u64 {
        self.backoff_us
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[test]
    fn fresh_backoff_is_initial() {
        let w = WorkerCore::new(0, None);
        assert_eq!(w.policy(), BackoffPolicy::ParkedSleep);
        assert_eq!(w.current_backoff_us(), SLEEP_INITIAL_US);
    }

    // One doubling-then-cap sequence per policy. Asserting the whole sequence
    // (not just the final cap) would fail a buggy `increase_backoff` that jumped
    // straight to the cap or incremented by a constant; the closing reset
    // confirms it returns to *this* policy's initial, not the other's.
    #[rstest]
    #[case::parked_sleep(
        WorkerCore::new(0, None),
        SLEEP_INITIAL_US,
        SLEEP_MAX_US,
        BackoffPolicy::ParkedSleep
    )]
    #[case::park(WorkerCore::driver(StepIdx(0)), PARK_INITIAL_US, PARK_MAX_US, BackoffPolicy::Park)]
    fn backoff_doubles_then_caps(
        #[case] mut w: WorkerCore,
        #[case] initial: u64,
        #[case] max: u64,
        #[case] policy: BackoffPolicy,
    ) {
        assert_eq!(w.policy(), policy);
        assert_eq!(w.current_backoff_us(), initial);
        let mut expected = initial;
        for _ in 0..20 {
            w.increase_backoff();
            expected = (expected * 2).min(max);
            assert_eq!(w.current_backoff_us(), expected);
        }
        assert_eq!(w.current_backoff_us(), max);
        w.reset_backoff();
        assert_eq!(w.current_backoff_us(), initial);
    }

    /// The override pins both bounds for this core only; a core built without it
    /// keeps the normal ramp (no process-global state).
    #[test]
    fn backoff_override_is_per_core() {
        let mut slow = WorkerCore::driver(StepIdx(0)).with_backoff_override(Some(10_000_000));
        let normal = WorkerCore::driver(StepIdx(0));
        assert_eq!(slow.current_backoff_us(), 10_000_000);
        slow.increase_backoff();
        assert_eq!(slow.current_backoff_us(), 10_000_000);
        slow.reset_backoff();
        assert_eq!(slow.current_backoff_us(), 10_000_000);
        assert_eq!(normal.current_backoff_us(), PARK_INITIAL_US);
        assert_eq!(
            WorkerCore::driver(StepIdx(0)).with_backoff_override(None).current_backoff_us(),
            PARK_INITIAL_US
        );
    }

    #[test]
    fn park_cap_is_two_milliseconds() {
        // The driver timer is self-heal only; 2 ms bounds a missed wake at 4×
        // the old cap and cuts timer expiries 4×.
        assert_eq!(PARK_MAX_US, 2_000);
        assert_eq!(PARK_INITIAL_US, 10);
        assert_ne!(SLEEP_MAX_US, PARK_MAX_US);
        assert_ne!(SLEEP_INITIAL_US, PARK_INITIAL_US);
    }

    #[test]
    fn worker_with_sticky_serial_owner_only() {
        let w = WorkerCore::new(2, Some(StepIdx(0)));
        assert_eq!(w.thread_id, 2);
        assert_eq!(w.sticky_owner, Some(StepIdx(0)));
    }

    #[test]
    fn test_backoff_entries_select_by_target() {
        let entries = [
            TestBackoff { target: TestBackoffTarget::Worker(1), us: 7 },
            TestBackoff { target: TestBackoffTarget::AllDrivers, us: 9 },
        ];
        assert_eq!(backoff_override_for(&entries, TestBackoffTarget::Worker(1)), Some(7));
        assert_eq!(backoff_override_for(&entries, TestBackoffTarget::Worker(0)), None);
        assert_eq!(backoff_override_for(&entries, TestBackoffTarget::Driver(3)), Some(9));
        assert_eq!(backoff_override_for(&[], TestBackoffTarget::Worker(0)), None);
    }
}
