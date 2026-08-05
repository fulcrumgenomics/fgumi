//! `WorkerCore`: per-thread state carried by the worker loop.

use std::time::Duration;

use crate::topology::StepIdx;

// Pool worker idle bounds: a pool worker is never unparked, so it sleeps and can
// ramp to a coarse cap without hurting wake latency.
const SLEEP_INITIAL_US: u64 = 1;
const SLEEP_MAX_US: u64 = 50_000; // 50 milliseconds

// Dedicated-driver idle bounds: a driver drives a small step subset off the pool
// and must stay responsive to its peers (the pool filling/draining its edges), so
// it parks with a tight cap. `park_timeout` also lets a producer holding the
// thread handle unpark it early; the cap is the bounded fallback either way.
const PARK_INITIAL_US: u64 = 10;
const PARK_MAX_US: u64 = 500;

/// How a worker idles on a no-progress tick. Selected per thread at construction:
/// the N-worker pool uses [`Sleep`](BackoffPolicy::Sleep); each dedicated driver
/// thread (the unified "1-thread pool") uses [`Park`](BackoffPolicy::Park). This
/// is the *only* behavioral difference between a pool worker and a driver — both
/// run the same [`run_worker_loop`](crate::runtime::run_worker_loop).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BackoffPolicy {
    /// `thread::sleep`, ramp 1µs→50ms. For pool workers (never unparked).
    Sleep,
    /// `thread::park_timeout`, ramp 10µs→500µs. For dedicated driver threads: a
    /// producer *may* unpark early, and the tight cap bounds wake latency when it
    /// does not. With no unparker this behaves as a bounded sleep.
    Park,
}

impl BackoffPolicy {
    #[inline]
    fn initial_us(self) -> u64 {
        match self {
            Self::Sleep => SLEEP_INITIAL_US,
            Self::Park => PARK_INITIAL_US,
        }
    }

    #[inline]
    fn max_us(self) -> u64 {
        match self {
            Self::Sleep => SLEEP_MAX_US,
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
    /// If this worker owns an `Exclusive` step, that step's index.
    /// Used by the runtime to skip `WorkerStepEntry::Skip` placeholders
    /// for Exclusive steps owned by other workers.
    pub exclusive_owner: Option<StepIdx>,
    /// If this worker is the sole eligible dispatcher for a `sticky` step
    /// (either an `Exclusive sticky` step it owns, or a `Serial + sticky`
    /// step whose `Affinity` targets this worker), the step's index.
    /// The driver drives this step in a tight inner loop until it returns
    /// `NoProgress` / `Contention` / `Finished`, then yields to round-
    /// robin. Mirrors the legacy pipeline's sticky read.
    pub sticky_owner: Option<StepIdx>,
    /// Pool worker vs dedicated driver. Determines both the idle backoff policy
    /// (`Pool → Sleep`, `Driver → Park`, via [`Self::policy`]) and how the
    /// thread's aggregate busy/idle time is attributed in `--pipeline-stats`.
    role: WorkerRole,
    /// Backoff duration in microseconds. Doubled on no-progress; reset on progress.
    /// Bounds come from `role`'s policy.
    backoff_us: u64,
}

impl WorkerCore {
    /// A pool worker (`Sleep` backoff). Existing callers are unchanged.
    #[must_use]
    pub fn new(
        thread_id: usize,
        exclusive_owner: Option<StepIdx>,
        sticky_owner: Option<StepIdx>,
    ) -> Self {
        Self {
            thread_id,
            exclusive_owner,
            sticky_owner,
            role: WorkerRole::Pool,
            backoff_us: BackoffPolicy::Sleep.initial_us(),
        }
    }

    /// A dedicated driver thread (the unified "1-thread pool"): `Driver` role +
    /// `Park` backoff (derived from the role). Drives a set of `Owned` steps off
    /// the pool; its aggregate busy/idle is attributed to `primary_step` on the
    /// off-pool detached line. `thread_id` is unused for drivers; it never owns
    /// Exclusive or sticky steps.
    #[must_use]
    pub fn driver(primary_step: StepIdx) -> Self {
        let mut worker = Self::new(0, None, None);
        worker.role = WorkerRole::Driver { primary_step };
        worker.backoff_us = worker.policy().initial_us();
        worker
    }

    /// This thread's pool/driver role (drives stats attribution in the loop).
    #[must_use]
    pub fn role(&self) -> WorkerRole {
        self.role
    }

    /// The idle backoff policy implied by this thread's role: pool workers
    /// `Sleep` (never unparked), driver threads `Park`.
    #[must_use]
    pub fn policy(&self) -> BackoffPolicy {
        match self.role {
            WorkerRole::Pool => BackoffPolicy::Sleep,
            WorkerRole::Driver { .. } => BackoffPolicy::Park,
        }
    }

    pub fn reset_backoff(&mut self) {
        self.backoff_us = self.policy().initial_us();
    }

    pub fn sleep_backoff(&self) {
        let dur = Duration::from_micros(self.backoff_us);
        match self.policy() {
            BackoffPolicy::Sleep => std::thread::sleep(dur),
            BackoffPolicy::Park => std::thread::park_timeout(dur),
        }
    }

    pub fn increase_backoff(&mut self) {
        self.backoff_us = self.backoff_us.saturating_mul(2).min(self.policy().max_us());
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
        let w = WorkerCore::new(0, None, None);
        assert_eq!(w.current_backoff_us(), SLEEP_INITIAL_US);
    }

    // One doubling-then-cap sequence, exercised for both policies: the pool's
    // `Sleep` bounds and a driver's tight `Park` bounds. Asserting the whole
    // sequence (not just the final cap) would fail a buggy `increase_backoff`
    // that jumped straight to the cap or incremented by a constant; the closing
    // reset confirms it returns to *this* policy's initial, not the other's.
    #[rstest]
    #[case::sleep(
        WorkerCore::new(0, None, None),
        SLEEP_INITIAL_US,
        SLEEP_MAX_US,
        BackoffPolicy::Sleep
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
            assert_eq!(
                w.current_backoff_us(),
                expected,
                "backoff must double each step until it saturates at the cap"
            );
        }
        assert_eq!(w.current_backoff_us(), max);
        w.reset_backoff();
        assert_eq!(w.current_backoff_us(), initial);
    }

    #[test]
    fn reset_after_progress() {
        let mut w = WorkerCore::new(0, None, None);
        for _ in 0..5 {
            w.increase_backoff();
        }
        w.reset_backoff();
        assert_eq!(w.current_backoff_us(), SLEEP_INITIAL_US);
    }

    #[test]
    fn sleep_and_park_have_distinct_caps() {
        // Guard against the two policies drifting to the same bound.
        assert_ne!(SLEEP_MAX_US, PARK_MAX_US);
        assert_ne!(SLEEP_INITIAL_US, PARK_INITIAL_US);
        // Policy is derived from role: pool worker → Sleep, driver → Park.
        assert_eq!(WorkerCore::new(0, None, None).policy(), BackoffPolicy::Sleep);
        assert_eq!(WorkerCore::driver(StepIdx(0)).policy(), BackoffPolicy::Park);
    }

    #[test]
    fn worker_with_exclusive_role() {
        let w = WorkerCore::new(2, Some(StepIdx(7)), Some(StepIdx(7)));
        assert_eq!(w.thread_id, 2);
        assert_eq!(w.exclusive_owner, Some(StepIdx(7)));
        assert_eq!(w.sticky_owner, Some(StepIdx(7)));
    }

    #[test]
    fn worker_with_sticky_serial_owner_only() {
        // E.g., the Serial+Affinity::Reader source — sticky_owner set,
        // exclusive_owner None.
        let w = WorkerCore::new(0, None, Some(StepIdx(0)));
        assert_eq!(w.exclusive_owner, None);
        assert_eq!(w.sticky_owner, Some(StepIdx(0)));
    }
}
