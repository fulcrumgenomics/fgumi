//! Shared admission caps for pool steps.
//!
//! A [`PhaseCap`] is a counter shared by every pool step of one sort phase: at
//! most `max` workers are inside those steps' work sections at once. The cap
//! belongs to the *phase*, not the step — `--sort-threads 4` means four workers
//! in total across inflate and spill compression, not four in each — which is
//! what "cede cores to an upstream producer" promises. Nothing here knows what a
//! phase is; the chain builder decides which step instances share which
//! counter, and the step types stay phase-neutral. A step the builder did not
//! cap holds no cap at all and never touches a counter.
//!
//! ## Admission protocol
//!
//! A capped step calls [`admit_input`] in `try_run` after flushing any held
//! output and before popping. An empty input takes no permit: the poll is
//! recorded as empty on the edge (so starvation metrics still see it) and the
//! step reports `Finished` if the input is drained, else `NoProgress`.
//! Otherwise the step takes a permit; a refused clone reports
//! [`StepOutcome::Capped`] — "work may exist but the cap will not let this
//! worker take it" — and is redispatched later, so a full cap never loses the
//! drain hand-off: the refused clone eventually acquires, or finds the input
//! drained and reports `Finished` like any other clone. An uncapped step skips
//! all of this and pops exactly as it would without a cap.
//!
//! A step may scope the permit to its CPU-heavy section only.
//! `SortSpillDecompress` holds it around the slot fill (spill read +
//! decompress), and takes it only when some live slot can be filled; forwarding
//! its input events to the merge is bookkeeping and runs uncapped, so the merge
//! is never starved of announcements behind decompressing clones.
//!
//! ## Off-pool parallel work takes the whole cap
//!
//! Work that runs off the pool but shares the phase's cores — the per-run sort
//! of a run and the in-memory fast-path gather — runs at the phase's full width,
//! so it takes **every** permit of the cap through [`PhaseCap::acquire_whole`]:
//! it registers a pending reservation, which makes [`PhaseCap::try_acquire`]
//! refuse new pool admissions (counted as refusals, reported `Capped`), parks
//! until the last pool holder releases (each releases after one item), takes
//! all `max` permits at once and runs on exactly the phase width. No deadlock:
//! a reserving caller holds no permit while it waits (it takes all of them in
//! one step, or none), and a pool step holds one only inside a single
//! `try_run`, never waiting while it holds it. The wait observes the pipeline's
//! cancel signal (bound by `Pipeline::run`), so a cancelled run does not leave a
//! sealing driver parked.

use std::sync::{Arc, Weak};

// The admission protocol's atomics and its parking lock come from `loom` under
// `--cfg loom`, so `tests/loom_admission.rs` checks the REAL `PhaseCap` (the
// same arrangement as `fgumi-sort`'s `merge_slots.rs`).
#[cfg(loom)]
use loom::sync::atomic::{AtomicU64, AtomicUsize, Ordering, fence};
#[cfg(not(loom))]
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering, fence};

use crate::item::HeapSize;
use crate::padded::Padded;
use crate::signal::PipelineSignal;
use crate::step::{InputHandle, StepOutcome};

/// The lock + condvar a whole-cap caller parks on. `parking_lot` in a normal
/// build, `loom`'s under `--cfg loom`.
struct Parker {
    #[cfg(not(loom))]
    lock: parking_lot::Mutex<()>,
    #[cfg(not(loom))]
    cvar: parking_lot::Condvar,
    #[cfg(loom)]
    lock: loom::sync::Mutex<()>,
    #[cfg(loom)]
    cvar: loom::sync::Condvar,
}

impl Parker {
    fn new() -> Self {
        Self {
            #[cfg(not(loom))]
            lock: parking_lot::Mutex::new(()),
            #[cfg(not(loom))]
            cvar: parking_lot::Condvar::new(),
            #[cfg(loom)]
            lock: loom::sync::Mutex::new(()),
            #[cfg(loom)]
            cvar: loom::sync::Condvar::new(),
        }
    }

    /// Wake every parked caller, under the lock, so a caller between its
    /// check and its wait cannot miss it.
    fn notify_all(&self) {
        #[cfg(not(loom))]
        let _guard = self.lock.lock();
        #[cfg(loom)]
        let _guard = self.lock.lock().expect("phase-cap lock poisoned");
        self.cvar.notify_all();
    }

    /// Park until `ready` returns `Some`, re-checking it under the lock after
    /// every wake.
    fn park_until<T>(&self, mut ready: impl FnMut() -> Option<T>) -> T {
        #[cfg(not(loom))]
        {
            let mut guard = self.lock.lock();
            loop {
                if let Some(v) = ready() {
                    return v;
                }
                self.cvar.wait(&mut guard);
            }
        }
        #[cfg(loom)]
        {
            let mut guard = self.lock.lock().expect("phase-cap lock poisoned");
            loop {
                if let Some(v) = ready() {
                    return v;
                }
                guard = self.cvar.wait(guard).expect("phase-cap lock poisoned");
            }
        }
    }
}

/// A shared admission counter: at most `max` holders of a permit at once.
///
/// Four hot words, each on its own line: `active` (written by every admitted
/// acquire and every release), `peak` (written only when the high-water mark
/// rises), `refused` (written by every refused acquire) and `reserving`
/// (written only by [`Self::acquire_whole`] callers, read by every acquire), so
/// refused workers never invalidate the line the admitted path reads or writes.
/// `max` and `name` are read-only after construction.
///
/// The whole-cap protocol adds to every admission one `Relaxed` load of
/// `reserving` (a line written only when a whole-cap caller arrives or
/// leaves), and to every release one `SeqCst` fence plus that load. Every sort
/// chain attaches a whole-cap caller to both of its caps (the per-run sort,
/// the fast-path gather), so there is no cap the protocol could be skipped
/// for. Measured cost (Apple M2 Max, release build, a `try_acquire` + drop
/// pair, median of 5 runs): 10.7 ns with the protocol vs 10.7 ns with the
/// fence and both `reserving` loads removed, uncontended (10^8 pairs on one
/// thread); with 8 threads on one cap of 8, ~1.8-1.9 µs per pair either way,
/// dominated by `active`'s cache line moving between cores. Within noise in
/// both cases, against a per-item work section of microseconds.
pub struct PhaseCap {
    active: Padded<AtomicUsize>,
    peak: Padded<AtomicUsize>,
    refused: Padded<AtomicU64>,
    /// [`Self::acquire_whole`] callers waiting for the cap to empty. While
    /// non-zero, [`Self::try_acquire`] refuses, so the waiter gets the next
    /// release instead of racing pool steps for it.
    reserving: Padded<AtomicUsize>,
    /// Parks a waiting [`Self::acquire_whole`] caller; a release that empties
    /// the cap while one is waiting, and the run's cancel, notify under it.
    parker: Parker,
    /// The current run's cancel signal, re-bound by every `Pipeline::run`; a
    /// parked [`Self::acquire_whole`] caller gives up when it is done.
    signal: parking_lot::Mutex<Weak<PipelineSignal>>,
    /// This cap, for registering with a run's signal (`Arc::new_cyclic`).
    this: Weak<PhaseCap>,
    max: usize,
    name: &'static str,
}

impl std::fmt::Debug for PhaseCap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PhaseCap")
            .field("name", &self.name)
            .field("max", &self.max)
            .field("active", &self.active())
            .field("peak", &self.peak())
            .field("refused", &self.refused())
            .finish_non_exhaustive()
    }
}

/// RAII permit: one unit of a cap's `active` count, released on drop. Borrows
/// the cap, so taking one costs no reference-count traffic.
#[must_use = "dropping the permit releases the slot immediately"]
#[derive(Debug)]
pub struct CapPermit<'a>(&'a PhaseCap);

impl Drop for CapPermit<'_> {
    fn drop(&mut self) {
        self.0.release(1);
    }
}

/// RAII guard for every permit of a cap (see [`PhaseCap::acquire_whole`]),
/// released together on drop.
#[must_use = "dropping the guard releases the whole cap immediately"]
#[derive(Debug)]
pub struct WholeCap<'a>(&'a PhaseCap);

impl WholeCap<'_> {
    /// How many permits the guard holds: the cap's whole width, which is the
    /// number of threads the off-pool work holding it may run on.
    #[must_use]
    pub fn width(&self) -> usize {
        self.0.max
    }
}

impl Drop for WholeCap<'_> {
    fn drop(&mut self) {
        self.0.release(self.0.max);
    }
}

/// The admission preamble every capped pool step runs before popping (see the
/// module docs). `Ok(Some(permit))` means the step may pop: bind the permit to
/// a named variable (not `_`) and hold it until the capped work is done;
/// `Ok(None)` means the step is uncapped and pops exactly as it would with no
/// cap (its own empty-pop path records the poll). `Err` is the outcome to
/// return from `try_run` without popping: `Finished` or `NoProgress` for an
/// empty input (recorded as an empty poll on the edge), or `Capped` for a
/// refusal.
///
/// # Errors
///
/// Returns the early `StepOutcome` described above.
#[inline]
pub fn admit_input<'a, T: Send + HeapSize + 'static>(
    input: &dyn InputHandle<T>,
    cap: Option<&'a PhaseCap>,
) -> Result<Option<CapPermit<'a>>, StepOutcome> {
    let Some(cap) = cap else {
        return Ok(None);
    };
    if input.is_empty() {
        input.note_empty_poll();
        return Err(if input.is_drained() {
            StepOutcome::Finished
        } else {
            StepOutcome::NoProgress
        });
    }
    cap.try_acquire().map(Some).ok_or(StepOutcome::Capped)
}

impl PhaseCap {
    /// A cap admitting at most `max` holders (`0` is treated as 1: a zero cap
    /// would wedge the phase).
    #[must_use]
    pub fn new(name: &'static str, max: usize) -> Arc<Self> {
        Arc::new_cyclic(|this| Self {
            active: Padded(AtomicUsize::new(0)),
            peak: Padded(AtomicUsize::new(0)),
            refused: Padded(AtomicU64::new(0)),
            reserving: Padded(AtomicUsize::new(0)),
            parker: Parker::new(),
            signal: parking_lot::Mutex::new(Weak::new()),
            this: this.clone(),
            max: max.max(1),
            name,
        })
    }

    /// Bind `signal` as this cap's cancel signal, replacing any earlier run's,
    /// and register the cap with it so the run's cancel or error wakes a
    /// parked [`Self::acquire_whole`] caller. Called by `Pipeline::run` for
    /// every step's cap before any thread starts.
    pub(crate) fn bind_signal(&self, signal: &Arc<PipelineSignal>) {
        *self.signal.lock() = Arc::downgrade(signal);
        signal.register_cap(self.this.clone());
    }

    /// [`Self::bind_signal`] for tests outside the crate (the loom model and
    /// downstream crates' whole-cap cancel tests). Production binding is
    /// `Pipeline::run`'s alone.
    #[cfg(any(loom, feature = "test-utils"))]
    #[doc(hidden)]
    pub fn bind_signal_for_test(&self, signal: &Arc<PipelineSignal>) {
        self.bind_signal(signal);
    }

    fn is_cancelled(&self) -> bool {
        self.signal.lock().upgrade().is_some_and(|s| s.is_done())
    }

    /// Wake a parked [`Self::acquire_whole`] caller so it re-checks the run's
    /// signal. Called by `PipelineSignal` after its terminal transition.
    pub(crate) fn notify_cancel(&self) {
        self.parker.notify_all();
    }

    /// Take one permit, or `None` (counted in [`Self::refused`]) when the cap is
    /// full or an [`Self::acquire_whole`] caller is waiting for it.
    ///
    /// Uncontended cost: two `Relaxed`/`Acquire` loads (`active` and the
    /// whole-cap `reserving` count) and one CAS of `active`, plus one `Relaxed`
    /// load of `peak` (and a `fetch_max` only when the mark rises).
    #[inline]
    #[must_use]
    pub fn try_acquire(&self) -> Option<CapPermit<'_>> {
        let active = &self.active.0;
        let mut cur = active.load(Ordering::Acquire);
        loop {
            // A pending whole-cap reservation wins the next release: refusing
            // here is what lets the waiter's wait end. `Relaxed` is safe: this
            // load only decides fairness, never correctness. Missing a
            // just-registered reservation admits one more pool item, which the
            // waiter then waits out — its own CAS from 0 is what guarantees it
            // never holds the cap beside a pool step, and the release-side
            // fence pair is what guarantees it is woken.
            if cur >= self.max || self.reserving.0.load(Ordering::Relaxed) > 0 {
                self.refused.0.fetch_add(1, Ordering::Relaxed);
                return None;
            }
            match active.compare_exchange_weak(cur, cur + 1, Ordering::AcqRel, Ordering::Acquire) {
                Ok(_) => {
                    self.note_peak(cur + 1);
                    return Some(CapPermit(self));
                }
                Err(observed) => cur = observed,
            }
        }
    }

    #[inline]
    fn note_peak(&self, now: usize) {
        let peak = &self.peak.0;
        if now > peak.load(Ordering::Relaxed) {
            peak.fetch_max(now, Ordering::Relaxed);
        }
    }

    /// Give back `n` permits. A release that empties the cap while a
    /// [`Self::acquire_whole`] caller is waiting wakes it. The `SeqCst` fence
    /// pair (here between the `fetch_sub` and the `reserving` load; in the
    /// waiter between its `reserving` increment and its `active` CAS) means
    /// either the waiter's CAS sees the cap empty or this release sees the
    /// reservation and wakes it (the notify takes the lock the waiter holds
    /// from its CAS until it parks, so it cannot fall in between).
    #[inline]
    fn release(&self, n: usize) {
        let prev = self.active.0.fetch_sub(n, Ordering::AcqRel);
        fence(Ordering::SeqCst);
        if prev == n && self.reserving.0.load(Ordering::Relaxed) > 0 {
            self.parker.notify_all();
        }
    }

    /// Every permit of the cap, for off-pool work that runs at the phase's full
    /// width (the per-run sort, the fast-path gather). Registers a reservation
    /// so pool steps stop being admitted, parks until the last holder releases,
    /// then takes all `max` permits in one step. Returns `None` only when the
    /// run's cancel signal (bound by `Pipeline::run`) is done: the run is
    /// being torn down, and the caller must skip its work.
    ///
    /// Cannot deadlock: the caller holds no permit of this cap while it waits
    /// (it takes all of them at once or none), and every holder — a pool step
    /// inside one `try_run`, or another whole-cap holder doing bounded work —
    /// releases without waiting on anything this caller holds. Its wake comes
    /// from the release that empties the cap (the fence pair in `release`) or
    /// from the run's cancel (`notify_cancel`), each delivered under the lock
    /// the caller checks under.
    ///
    /// The cap must be bound to the current run's signal, which `Pipeline::run`
    /// does for every cap a step reports through [`Step::phase_cap`]. An unbound
    /// cap's wait could never observe a cancel (the wait is untimed), so a step
    /// that takes the whole cap without reporting it would park forever on a
    /// cancelled run. Debug builds assert the binding here, which turns that
    /// contract into a failure in any test that drives the step.
    ///
    /// # Panics
    ///
    /// In debug builds, if no live cancel signal is bound (see above).
    ///
    /// [`Step::phase_cap`]: crate::step::Step::phase_cap
    #[must_use]
    pub fn acquire_whole(&self) -> Option<WholeCap<'_>> {
        debug_assert!(
            self.signal.lock().strong_count() > 0,
            "PhaseCap `{}`: acquire_whole on a cap no live run signal is bound to; a step that \
             takes the whole cap must report it through `Step::phase_cap`",
            self.name
        );
        self.reserving.0.fetch_add(1, Ordering::Relaxed);
        // Pairs with the fence in `release` (see there).
        fence(Ordering::SeqCst);
        let taken = self.parker.park_until(|| {
            if self
                .active
                .0
                .compare_exchange(0, self.max, Ordering::AcqRel, Ordering::Acquire)
                .is_ok()
            {
                Some(true)
            } else if self.is_cancelled() {
                Some(false)
            } else {
                None
            }
        });
        self.reserving.0.fetch_sub(1, Ordering::Relaxed);
        taken.then(|| {
            self.note_peak(self.max);
            WholeCap(self)
        })
    }

    /// [`Self::acquire_whole`] callers currently waiting (tests).
    #[cfg(any(test, loom, feature = "test-utils"))]
    #[must_use]
    pub fn pending_reservations(&self) -> usize {
        self.reserving.0.load(Ordering::Relaxed)
    }

    /// The label passed to [`Self::new`].
    #[must_use]
    pub fn name(&self) -> &'static str {
        self.name
    }

    /// The admission limit.
    #[must_use]
    pub fn max(&self) -> usize {
        self.max
    }

    /// Short label for `Pipeline::dag()` and diagnostics, e.g. `sort-phase1(4)`.
    #[must_use]
    pub fn describe(&self) -> String {
        format!("{}({})", self.name, self.max)
    }

    /// Permits currently held.
    #[must_use]
    pub fn active(&self) -> usize {
        self.active.0.load(Ordering::Acquire)
    }

    /// Highest number of permits held at once so far.
    #[must_use]
    pub fn peak(&self) -> usize {
        self.peak.0.load(Ordering::Relaxed)
    }

    /// Pool polls deferred so far because the cap was full or reserved (a step
    /// reported `Capped` and retried later). Normal when the cap binds.
    #[must_use]
    pub fn refused(&self) -> u64 {
        self.refused.0.load(Ordering::Relaxed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::thread;

    /// N threads hammer one cap of 3; the number of simultaneous holders
    /// (observed from inside the critical section) never exceeds 3, every
    /// permit is released, and the cap's own high-water mark agrees.
    #[test]
    fn phase_cap_admits_at_most_max() {
        let cap = PhaseCap::new("t", 3);
        let inside = Arc::new(AtomicUsize::new(0));
        let observed_peak = Arc::new(AtomicUsize::new(0));
        let handles: Vec<_> = (0..16)
            .map(|_| {
                let cap = Arc::clone(&cap);
                let inside = Arc::clone(&inside);
                let observed_peak = Arc::clone(&observed_peak);
                thread::spawn(move || {
                    for _ in 0..5_000 {
                        if let Some(_permit) = cap.try_acquire() {
                            let now = inside.fetch_add(1, Ordering::AcqRel) + 1;
                            observed_peak.fetch_max(now, Ordering::AcqRel);
                            assert!(now <= 3, "{now} holders inside a cap of 3");
                            inside.fetch_sub(1, Ordering::AcqRel);
                        }
                    }
                })
            })
            .collect();
        for h in handles {
            h.join().expect("worker thread panicked");
        }
        assert_eq!(cap.active(), 0, "every permit must be released");
        assert!(cap.peak() <= 3 && cap.peak() >= observed_peak.load(Ordering::Acquire));
    }

    /// Refusal is deterministic, not timing-dependent: with every permit held,
    /// the next acquire is refused and counted, and a release admits again.
    #[test]
    fn full_cap_refuses_and_counts_until_a_permit_is_released() {
        let cap = PhaseCap::new("r", 2);
        let p1 = cap.try_acquire().expect("1 of 2");
        let _p2 = cap.try_acquire().expect("2 of 2");
        assert!(cap.try_acquire().is_none(), "full cap refuses");
        assert_eq!(cap.refused(), 1);
        assert_eq!(cap.peak(), 2);
        drop(p1);
        assert!(cap.try_acquire().is_some(), "a release admits again");
        assert_eq!(cap.refused(), 1);
    }

    /// Wait until `cap` shows `n` pending whole-cap reservations.
    fn await_reservations(cap: &PhaseCap, n: usize) {
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
        while cap.pending_reservations() != n {
            assert!(std::time::Instant::now() < deadline, "reservation never registered");
            thread::yield_now();
        }
    }

    /// A whole-cap caller waits for every holder, and while it waits a pool
    /// `try_acquire` racing a released permit loses to it (refused, counted),
    /// so the waiter is not starved; it then holds exactly `max` permits.
    #[test]
    fn whole_cap_waiter_wins_released_permits_over_pool_admission() {
        let cap = PhaseCap::new("w", 4);
        let signal = PipelineSignal::new();
        cap.bind_signal(&signal);
        let held1 = cap.try_acquire().expect("1 of 4");
        let held2 = cap.try_acquire().expect("2 of 4");
        let waiter = {
            let cap = Arc::clone(&cap);
            thread::spawn(move || {
                let whole = cap.acquire_whole().expect("not cancelled");
                (whole.width(), cap.active())
            })
        };
        await_reservations(&cap, 1);
        drop(held1);
        assert!(cap.try_acquire().is_none(), "a pending reservation refuses pool admission");
        assert_eq!(cap.refused(), 1, "and the refusal is counted");
        drop(held2);
        assert_eq!(waiter.join().expect("waiter"), (4, 4), "the whole cap, nothing more");
        assert_eq!(cap.active(), 0, "released together");
        assert_eq!(cap.pending_reservations(), 0);
        assert_eq!(cap.peak(), 4);
        assert!(cap.try_acquire().is_some(), "pool admission resumes");
    }

    /// Spawn a whole-cap caller and return a receiver for its result
    /// (`true` = it got the cap), so a test can bound its wait.
    fn whole_cap_caller(cap: &Arc<PhaseCap>) -> std::sync::mpsc::Receiver<bool> {
        let (tx, rx) = std::sync::mpsc::channel();
        let cap = Arc::clone(cap);
        thread::spawn(move || {
            let got = cap.acquire_whole().is_some();
            let _ = tx.send(got);
        });
        rx
    }

    /// A whole-cap caller on a cap with no live run signal fails fast in a
    /// debug build instead of parking where no cancel can reach it: never
    /// bound, or bound to a run whose signal is gone. The cap is idle, so
    /// without the assertion the call would succeed and the test would fail.
    #[rstest::rstest]
    #[case::never_bound(false)]
    #[case::bound_signal_dropped(true)]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "no live run signal is bound")]
    fn whole_cap_requires_a_bound_signal(#[case] bind_then_drop: bool) {
        let cap = PhaseCap::new("unbound", 2);
        if bind_then_drop {
            cap.bind_signal(&PipelineSignal::new());
        }
        let _ = cap.acquire_whole();
    }

    /// A parked whole-cap caller is woken by the run's cancel (not a timer:
    /// the wait is untimed), gives up, and leaves no reservation behind.
    #[test]
    fn whole_cap_wait_observes_cancel() {
        let cap = PhaseCap::new("c", 2);
        let signal = PipelineSignal::new();
        cap.bind_signal(&signal);
        let held = cap.try_acquire().expect("holder");
        let rx = whole_cap_caller(&cap);
        await_reservations(&cap, 1);
        signal.cancel();
        let got = rx.recv_timeout(std::time::Duration::from_secs(10)).expect("woken by the cancel");
        assert!(!got, "cancelled: no permits");
        assert_eq!(cap.pending_reservations(), 0);
        assert_eq!(cap.active(), 1, "only the pool holder's permit");
        drop(held);
    }

    /// A cap reused by a second run follows that run's signal: cancelling the
    /// second run wakes a caller parked during it, even though the first
    /// run's signal is gone.
    #[test]
    fn a_reused_cap_follows_the_latest_runs_signal() {
        let cap = PhaseCap::new("r", 1);
        let first = PipelineSignal::new();
        cap.bind_signal(&first);
        drop(first);
        let second = PipelineSignal::new();
        cap.bind_signal(&second);
        let held = cap.try_acquire().expect("holder");
        let rx = whole_cap_caller(&cap);
        await_reservations(&cap, 1);
        second.cancel();
        let got = rx.recv_timeout(std::time::Duration::from_secs(10)).expect("woken by the cancel");
        assert!(!got);
        drop(held);
    }

    /// A minimal input for `admit_input`: a fixed emptiness/drain state that
    /// counts `note_empty_poll` calls.
    struct FakeInput {
        empty: bool,
        drained: bool,
        noted: AtomicUsize,
    }
    impl InputHandle<u32> for FakeInput {
        fn pop(&self) -> Option<u32> {
            None
        }
        fn is_drained(&self) -> bool {
            self.drained
        }
        fn is_empty(&self) -> bool {
            self.empty
        }
        fn note_empty_poll(&self) {
            self.noted.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// The capped preamble records an empty poll on the edge (so starvation
    /// stats still see it) and takes no permit; an uncapped step is left to
    /// its own pop, which records the poll itself.
    #[test]
    fn admit_input_records_empty_polls_and_takes_no_permit() {
        let cap = PhaseCap::new("e", 1);
        let input = FakeInput { empty: true, drained: false, noted: AtomicUsize::new(0) };
        assert!(matches!(admit_input(&input, Some(&cap)), Err(StepOutcome::NoProgress)));
        assert_eq!(input.noted.load(Ordering::Relaxed), 1, "the idle poll is recorded");
        assert_eq!((cap.active(), cap.refused()), (0, 0), "and takes no permit");
        assert!(matches!(admit_input(&input, None), Ok(None)));
        assert_eq!(
            input.noted.load(Ordering::Relaxed),
            1,
            "uncapped: the step's own pop records it"
        );

        let drained = FakeInput { empty: true, drained: true, noted: AtomicUsize::new(0) };
        assert!(matches!(admit_input(&drained, Some(&cap)), Err(StepOutcome::Finished)));

        let full = FakeInput { empty: false, drained: false, noted: AtomicUsize::new(0) };
        let _held = cap.try_acquire().expect("the only permit");
        assert!(matches!(admit_input(&full, Some(&cap)), Err(StepOutcome::Capped)));
    }

    #[test]
    fn permit_release_on_drop_including_panic() {
        let cap = PhaseCap::new("p", 1);
        let p1 = cap.try_acquire().expect("first permit");
        assert!(cap.try_acquire().is_none(), "cap of 1 is full");
        drop(p1);
        assert!(cap.try_acquire().is_some(), "slot freed on drop");

        // A panic while holding the permit must still release it (RAII).
        let cap2 = Arc::clone(&cap);
        let r = thread::spawn(move || {
            let _held = cap2.try_acquire().expect("permit");
            panic!("boom");
        })
        .join();
        assert!(r.is_err());
        assert_eq!(cap.active(), 0);
        assert!(cap.try_acquire().is_some(), "slot released by the unwinding thread");
        assert_eq!(cap.describe(), "p(1)");
    }

    /// `max == 0` is treated as 1 (a zero cap would wedge the phase).
    #[test]
    fn zero_max_admits_exactly_one() {
        let cap = PhaseCap::new("z", 0);
        let _p = cap.try_acquire().expect("one permit");
        assert!(cap.try_acquire().is_none());
        assert_eq!(cap.max(), 1);
    }
}
