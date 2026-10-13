//! `WorkerCore`: per-thread state carried by the worker loop.

use std::time::Duration;

use crate::runtime::placement::DriverIdx;
use crate::topology::StepIdx;

// Idle bounds for an unpinned pool worker (`BackoffPolicy::EventCount`). It
// waits on the pool event-count, and every item it can consume reaches it as a
// directed wake: a pool notify gated on the push or, for an item it holds, an
// `unpark` from the consumer whose pop makes room (a reorder-stash-cap hold
// included). Its deadline therefore does two jobs, neither of them delivering
// work in normal operation: it is the self-heal bound on a wake nobody
// delivers, and the retry interval of a held item (a holder timer-parks on the
// same deadline). 200µs is that retry interval. Each deadline is jittered ±25%
// per worker (`jittered_wait_micros`), so idle waiters do not time out on the
// same microsecond: the ramp runs 200µs → 50 ms, each value ±25%.
const EC_WAIT_INITIAL_US: u64 = 200;
const EC_WAIT_MAX_US: u64 = 50_000; // 50 milliseconds, before the ±25% jitter
const EC_JITTER_DIVISOR: u64 = 4; // ±25%

// Timer-idle bounds for a pinned pool worker (an Exclusive owner, sticky owner,
// or the affinity target of a Serial step) and for the lone worker of a
// one-worker run (`BackoffPolicy::ParkedSleep`). Its producers (and, for a
// holder, the consumer whose pop makes room) `unpark` it through the wake plan
// (`runtime::wake`), so the timer is only the self-heal for a wake nobody
// delivers.
//
// The initial (also the post-progress reset floor) is 20µs, NOT 1µs. A pinned
// worker whose step is fed by an out-of-engine event (a reader's prefetch
// thread) or that also polls pool steps makes the occasional lucky pop, which
// resets its backoff; a 1µs floor then puts it back on the steep part of the
// ramp, waking ~every microsecond to run a full empty poll pass. 20µs is well
// below any step's per-item service time yet coarse enough to stop that
// re-poll.
const SLEEP_INITIAL_US: u64 = 20;
const SLEEP_MAX_US: u64 = 50_000; // 50 milliseconds

// Dedicated-driver idle bounds. A driver drives a small step subset off the
// pool; every producer of its live steps unparks it, so
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

/// How a thread idles on a no-progress tick, chosen per thread by its role in
/// `Pipeline::run`: an unpinned pool worker with
/// [`EventCount`](BackoffPolicy::EventCount); a pinned pool worker or the lone
/// worker with [`ParkedSleep`](BackoffPolicy::ParkedSleep); a dedicated driver
/// thread with [`Park`](BackoffPolicy::Park). `ParkedSleep` and `Park` both
/// `park_timeout` on the ramp value and differ only in the ramp; `EventCount`
/// waits on the pool event-count with a jittered deadline (and, while it holds
/// an item or a cap refused it, timer-parks on that deadline like the others).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BackoffPolicy {
    /// `thread::park_timeout`, ramp 20µs→50ms. A registered thread can be
    /// unparked early; with no unparker this is a bounded sleep.
    ParkedSleep,
    /// `thread::park_timeout`, ramp 10µs→2ms. For dedicated driver threads.
    Park,
    /// The pool event-count wait, ramp 200µs→50ms with each deadline jittered
    /// ±25% per `(worker, idle episode)` (so up to 62.5 ms at the cap). For
    /// unpinned pool workers.
    EventCount,
}

impl BackoffPolicy {
    #[inline]
    fn initial_us(self) -> u64 {
        match self {
            Self::ParkedSleep => SLEEP_INITIAL_US,
            Self::Park => PARK_INITIAL_US,
            Self::EventCount => EC_WAIT_INITIAL_US,
        }
    }

    #[inline]
    fn max_us(self) -> u64 {
        match self {
            Self::ParkedSleep => SLEEP_MAX_US,
            Self::Park => PARK_MAX_US,
            Self::EventCount => EC_WAIT_MAX_US,
        }
    }
}

/// `backoff_us` spread over ±25%, deterministic per `(worker_id, iter)`, so
/// idle event-count waiters do not all time out on the same microsecond. The
/// offset is computed two-sided: in `u64`, `seed % (2r) − r` would saturate
/// every negative offset to 0 and turn ±25% into 0..+25%.
///
/// The per-worker part of the offset is `worker_id` mixed by splitmix64; the
/// episode `iter` is added after mixing. One worker's consecutive episodes
/// therefore step its offset by +1µs and wrap every `2r` episodes (each offset
/// in `0..2r` once per period), while different workers start from unrelated
/// offsets. This departs from the previous engine's seed
/// (`worker_id × golden-ratio + iter`, in `fgumi-sort`'s `worker_pool.rs`):
/// modulo the small `2r` of the 200µs and 400µs levels, that seed gave workers
/// `w` and `w + 8` the same deadline, so 7 of 16 workers still shared one.
///
/// `backoff_us` is an `EventCount` ramp value, at least `EC_WAIT_INITIAL_US`
/// (an override is never jittered), so `2r` is never zero.
fn jittered_wait_micros(backoff_us: u64, worker_id: usize, iter: u64) -> u64 {
    debug_assert!(backoff_us >= EC_WAIT_INITIAL_US, "jitter below the event-count floor");
    let jitter_range = backoff_us / EC_JITTER_DIVISOR;
    let off = splitmix64(u64::try_from(worker_id).unwrap_or(u64::MAX)).wrapping_add(iter)
        % (jitter_range * 2);
    if off < jitter_range {
        backoff_us - (jitter_range - off)
    } else {
        backoff_us + (off - jitter_range)
    }
}

/// splitmix64's finaliser: a well-mixed 64-bit hash of `x`.
fn splitmix64(x: u64) -> u64 {
    let mut z = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
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
    /// `primary_step` (the group's representative). `driver` is its index in
    /// the run's wake plan: the parking bit its timer park arms
    /// (`WakePlan::arm_driver`), without which no directed wake reaches it.
    Driver { primary_step: StepIdx, driver: DriverIdx },
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
    /// time is attributed in `--pipeline-stats`.
    role: WorkerRole,
    /// How the thread idles, named by whoever builds the core.
    policy: BackoffPolicy,
    /// `(initial, max)` bounds of the ramp: the policy's, or a test override.
    bounds: (u64, u64),
    /// A test override pinned the bounds: the deadline is then exactly the
    /// override, un-jittered, whatever the policy.
    overridden: bool,
    /// Event-count idle episodes so far (one per `ec.wait`): the jitter's
    /// per-episode input.
    idle_episodes: u64,
    /// Backoff duration in microseconds. Doubled on no-progress; reset on progress.
    backoff_us: u64,
    /// This run's test hooks (`PipelineConfig::test_hooks`); `None` unless the
    /// test that built the run installed them.
    #[cfg(test)]
    hooks: Option<std::sync::Arc<crate::runtime::driver::test_hooks::TestHooks>>,
}

impl WorkerCore {
    /// A pool worker idling with `policy`: `EventCount` when it is unpinned in
    /// a multi-worker run (the loop then waits on the pool event-count, which
    /// must exist), `ParkedSleep` when it is pinned or the lone worker.
    #[must_use]
    pub fn new(thread_id: usize, sticky_owner: Option<StepIdx>, policy: BackoffPolicy) -> Self {
        Self::with_role(thread_id, sticky_owner, WorkerRole::Pool, policy)
    }

    /// A dedicated driver thread (the unified "1-thread pool"): `Driver` role +
    /// `Park` backoff. Drives a set of `Owned` steps off the pool; its
    /// aggregate busy/idle is attributed on the off-pool detached line.
    /// `thread_id` is unused for drivers; it never owns Exclusive or sticky
    /// steps. `driver` is the thread's index in the run's wake plan.
    #[must_use]
    pub fn driver(primary_step: StepIdx, driver: DriverIdx) -> Self {
        Self::with_role(0, None, WorkerRole::Driver { primary_step, driver }, BackoffPolicy::Park)
    }

    fn with_role(
        thread_id: usize,
        sticky_owner: Option<StepIdx>,
        role: WorkerRole,
        policy: BackoffPolicy,
    ) -> Self {
        let bounds = (policy.initial_us(), policy.max_us());
        Self {
            thread_id,
            sticky_owner,
            role,
            policy,
            bounds,
            overridden: false,
            idle_episodes: 0,
            backoff_us: bounds.0,
            #[cfg(test)]
            hooks: None,
        }
    }

    /// Attach this run's test hooks (`None` keeps the core hook-free).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_test_hooks(
        mut self,
        hooks: Option<std::sync::Arc<crate::runtime::driver::test_hooks::TestHooks>>,
    ) -> Self {
        self.hooks = hooks;
        self
    }

    /// This run's test hooks, if its test installed them.
    #[cfg(test)]
    pub(crate) fn test_hooks(
        &self,
    ) -> Option<&std::sync::Arc<crate::runtime::driver::test_hooks::TestHooks>> {
        self.hooks.as_ref()
    }

    /// Pin both ramp bounds to `us`: the shared body of both override paths.
    fn apply_override(&mut self, us: u64) {
        self.bounds = (us, us);
        self.backoff_us = us;
        self.overridden = true;
    }

    /// Replace both ramp bounds with `us` (a `TestBackoff` entry); `None`
    /// keeps the policy's. An override is exact: it is never jittered.
    #[must_use]
    pub fn with_backoff_override(mut self, us: Option<u64>) -> Self {
        if let Some(us) = us {
            self.apply_override(us);
        }
        self
    }

    /// [`Self::with_backoff_override`]`(Some(us))` on a live core: a test hook
    /// raises a running worker's timer mid-run.
    #[cfg(test)]
    pub(crate) fn set_override_for_test(&mut self, us: u64) {
        self.apply_override(us);
    }

    /// This thread's pool/driver role (drives stats attribution in the loop).
    #[must_use]
    pub fn role(&self) -> WorkerRole {
        self.role
    }

    /// This thread's idle policy.
    #[must_use]
    pub fn policy(&self) -> BackoffPolicy {
        self.policy
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

    /// The idle wait's deadline: the event-count `wait` timeout (the self-heal
    /// bound on a missed wakeup) and the deadline a timer park is classified
    /// against. With a test override: exactly the override, every policy.
    /// Otherwise `EventCount`: the ramp value jittered per
    /// `(thread_id, idle_episodes)` ([`jittered_wait_micros`]); `ParkedSleep`
    /// and `Park`: the ramp value.
    ///
    /// The worker loop's timer park is `park_timeout` of this (or of
    /// [`Self::max_backoff_deadline`]). The park token is the thread's std
    /// token, which std's own blocking primitives (`mpsc::recv`, `Mutex`) also
    /// use: a token consumed inside a step that blocks on one of those is not
    /// seen there, and that one wake degrades to the timer.
    #[must_use]
    pub fn backoff_deadline(&self) -> Duration {
        let us = match self.policy {
            BackoffPolicy::EventCount if !self.overridden => {
                jittered_wait_micros(self.backoff_us, self.thread_id, self.idle_episodes)
            }
            _ => self.backoff_us,
        };
        Duration::from_micros(us)
    }

    /// Count one event-count idle episode (called right before `ec.wait`).
    pub fn note_idle_episode(&mut self) {
        self.idle_episodes = self.idle_episodes.wrapping_add(1);
    }

    /// Event-count idle episodes so far.
    #[cfg(test)]
    #[must_use]
    pub fn idle_episodes(&self) -> u64 {
        self.idle_episodes
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
        let w = WorkerCore::new(0, None, BackoffPolicy::ParkedSleep);
        assert_eq!(w.policy(), BackoffPolicy::ParkedSleep);
        assert_eq!(w.current_backoff_us(), SLEEP_INITIAL_US);
    }

    // One doubling-then-cap sequence per policy. Asserting the whole sequence
    // (not just the final cap) would fail a buggy `increase_backoff` that jumped
    // straight to the cap or incremented by a constant; the closing reset
    // confirms it returns to *this* policy's initial, not the other's.
    #[rstest]
    #[case::parked_sleep(
        WorkerCore::new(0, None, BackoffPolicy::ParkedSleep),
        SLEEP_INITIAL_US,
        SLEEP_MAX_US,
        BackoffPolicy::ParkedSleep
    )]
    #[case::park(
        WorkerCore::driver(StepIdx(0), DriverIdx(0)),
        PARK_INITIAL_US,
        PARK_MAX_US,
        BackoffPolicy::Park
    )]
    #[case::event_count(
        WorkerCore::new(0, None, BackoffPolicy::EventCount),
        EC_WAIT_INITIAL_US,
        EC_WAIT_MAX_US,
        BackoffPolicy::EventCount
    )]
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
        let mut slow =
            WorkerCore::driver(StepIdx(0), DriverIdx(0)).with_backoff_override(Some(10_000_000));
        let normal = WorkerCore::driver(StepIdx(0), DriverIdx(0));
        assert_eq!(slow.current_backoff_us(), 10_000_000);
        slow.increase_backoff();
        assert_eq!(slow.current_backoff_us(), 10_000_000);
        slow.reset_backoff();
        assert_eq!(slow.current_backoff_us(), 10_000_000);
        assert_eq!(normal.current_backoff_us(), PARK_INITIAL_US);
        assert_eq!(
            WorkerCore::driver(StepIdx(0), DriverIdx(0))
                .with_backoff_override(None)
                .current_backoff_us(),
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
        let w = WorkerCore::new(2, Some(StepIdx(0)), BackoffPolicy::ParkedSleep);
        assert_eq!(w.thread_id, 2);
        assert_eq!(w.sticky_owner, Some(StepIdx(0)));
    }

    /// Two-sided and bounded; over one full period of `iter` every offset in
    /// `0..2r` occurs exactly once per worker, so the sum is exactly `2r·b − r`
    /// (mean `b − ½`). Bounds and both halves as `fgumi-sort`'s
    /// `jitter_spreads_waits_both_sides_of_the_backoff`, over more `b` and
    /// workers.
    #[test]
    fn jitter_is_two_sided_and_bounded_over_a_full_period() {
        for &b in &[EC_WAIT_INITIAL_US, 400, 800, 3_200, EC_WAIT_MAX_US] {
            let r = b / EC_JITTER_DIVISOR;
            for worker in 0..16usize {
                let (mut below, mut above, mut sum) = (0u64, 0u64, 0u64);
                let (mut lo, mut hi) = (u64::MAX, 0u64);
                for iter in 0..2 * r {
                    let v = jittered_wait_micros(b, worker, iter);
                    assert!(v >= b - r && v <= b + r, "b={b} w={worker} v={v}");
                    if v < b {
                        below += 1;
                    } else if v > b {
                        above += 1;
                    }
                    lo = lo.min(v);
                    hi = hi.max(v);
                    sum += v;
                }
                assert!(below > 0 && above > 0 && lo < b && b < hi, "both sides at b={b}");
                assert_eq!(sum, 2 * r * b - r, "full period at b={b} w={worker}");
            }
        }
    }

    /// Lockstep is broken at the levels where it matters. The offset is the
    /// worker's hash plus `iter`, so the count of distinct deadlines among `n`
    /// workers is the same for every `iter`. Computed: at b = 200 (2r = 100)
    /// and b = 400 (2r = 200), 8 workers take 8 distinct deadlines and 16 take
    /// 16; 32 take 28 and 29 (the old `worker × golden-ratio` seed gave 9 of 16
    /// and 10 of 32 at both levels: workers `w` and `w + 8` collided).
    #[rstest]
    #[case::floor_8(EC_WAIT_INITIAL_US, 8, 8)]
    #[case::floor_16(EC_WAIT_INITIAL_US, 16, 16)]
    #[case::floor_32(EC_WAIT_INITIAL_US, 32, 28)]
    #[case::second_8(400, 8, 8)]
    #[case::second_16(400, 16, 16)]
    #[case::second_32(400, 32, 29)]
    fn jitter_differs_across_workers_at_the_floor(
        #[case] b: u64,
        #[case] workers: usize,
        #[case] distinct: usize,
    ) {
        for iter in [0, 7, 99] {
            let n = (0..workers)
                .map(|w| jittered_wait_micros(b, w, iter))
                .collect::<std::collections::BTreeSet<u64>>()
                .len();
            assert_eq!(n, distinct, "b={b} workers={workers} iter={iter}");
        }
    }

    /// The policy drives the deadline, per episode.
    #[test]
    fn event_count_deadline_is_jittered_per_episode() {
        let mut w = WorkerCore::new(5, None, BackoffPolicy::EventCount);
        let mut seen = Vec::new();
        for k in 1..=8u64 {
            w.note_idle_episode();
            assert_eq!(w.idle_episodes(), k);
            let d = w.backoff_deadline();
            assert_eq!(d, Duration::from_micros(jittered_wait_micros(EC_WAIT_INITIAL_US, 5, k)));
            seen.push(d);
        }
        assert!(seen.windows(2).any(|p| p[0] != p[1]), "deadlines vary across episodes");
        // No "above nominal" assertion: for one worker, consecutive episodes
        // step the offset by +1µs, so eight episodes of worker 5 land at
        // 169..=176µs, all below nominal. Two-sidedness is pinned over a full
        // period by `jitter_is_two_sided_and_bounded_over_a_full_period`.
    }

    #[rstest]
    #[case::parked_sleep(BackoffPolicy::ParkedSleep)]
    #[case::park(BackoffPolicy::Park)]
    fn only_event_count_is_jittered(#[case] policy: BackoffPolicy) {
        let mut w = WorkerCore::new(3, None, policy);
        for _ in 0..4 {
            w.note_idle_episode();
            assert_eq!(w.backoff_deadline(), Duration::from_micros(w.current_backoff_us()));
            w.increase_backoff();
        }
    }

    /// A test override is exact (un-jittered) for every policy.
    #[rstest]
    #[case::event_count(BackoffPolicy::EventCount)]
    #[case::parked_sleep(BackoffPolicy::ParkedSleep)]
    #[case::park(BackoffPolicy::Park)]
    fn override_deadline_is_exact_for_every_policy(#[case] policy: BackoffPolicy) {
        let mut w = WorkerCore::new(5, None, policy).with_backoff_override(Some(10_000_000));
        for _ in 0..4 {
            w.note_idle_episode();
            assert_eq!(w.backoff_deadline(), Duration::from_secs(10));
        }
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
