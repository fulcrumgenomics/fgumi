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
use loom::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering, fence};
#[cfg(not(loom))]
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering, fence};

use crate::item::HeapSize;
use crate::padded::Padded;
use crate::runtime::wake_slot::HolderSet;
use crate::signal::PipelineSignal;
use crate::step::{InputHandle, StepOutcome};

/// Wake slots a cap's refused-worker set names one by one: `0..REFUSED_SLOTS`
/// get their own bit (four words, with the anonymous bit). The anonymous bit is
/// not a slot. It stands for every refused thread whose slot does not fit, or
/// that has none, and its wake is a broadcast (`CapWaker(None)`), so it can
/// never be delivered to the thread that happens to own a real slot of that
/// number.
pub(crate) const REFUSED_SLOTS: usize = 4 * 64 - 1;

/// The most distinct phase caps one run may use: a worker notes the caps it
/// is recorded on in one 64-bit thread-local mask, one bit per cap
/// (`WakePlan::start_pass`). `PipelineBuilder::build` refuses a pipeline
/// whose steps report more (`BuildError::TooManyPhaseCaps`), and
/// `WakePlan::build` checks it again. A sort chain uses two.
pub const MAX_PHASE_CAPS: usize = 64;

/// Delivers a cap-release wake: `Some(slot)` to one wake slot, `None` to every
/// thread (a refused thread without a slot in the cap's set). In a pipeline run
/// both go through `WakePlan::deliver_slot`/`deliver_anonymous`.
pub type CapWaker = Arc<dyn Fn(Option<usize>) + Send + Sync>;

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
    /// Wake slots of the threads (pool workers and detached drivers) this cap
    /// refused, or whose cap-parked re-poll skipped it as full, since a
    /// release last woke them: the queues' holder set, one bit per slot plus
    /// the anonymous one (see [`REFUSED_SLOTS`]). Set by the refusal and skip
    /// paths, claimed bit by bit by [`Self::release`] (and by a forwarded
    /// wake), and cleared by its owner when its next dispatch pass starts
    /// ([`Self::clear_own_record`]). Written on refusals, skips, claims and
    /// those clears, read by every release; its words are their own
    /// allocation, away from the hot words above.
    refused_slots: HolderSet,
    /// Whether this run delivers cap-release wakes (a Directed wake plan bound
    /// a waker). Off: refusals are not recorded and nothing changes (an
    /// untracked cap also has no run index, so an admission notes nothing).
    tracking: AtomicBool,
    /// Delivers a cap-release wake to a recorded wake slot. Set per run by
    /// `Pipeline::run`; read only when a release has claimed a slot.
    waker: parking_lot::Mutex<Option<CapWaker>>,
    /// This cap's index among the current run's caps (`WakePlan::build`;
    /// below [`MAX_PHASE_CAPS`]), or `usize::MAX` for a cap no Directed run
    /// tracks: the bit a worker's thread-local record masks use for it.
    /// Written once per run, before any thread starts; a std atomic even under
    /// loom, since it is configuration, not protocol.
    plan_index: std::sync::atomic::AtomicUsize,
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

/// Defines [`SkipOutcome`] at the visibility given: public only for the loom
/// model and tests (through `PhaseCap::record_skip_as`), crate-private
/// otherwise, so it stays out of the production API.
macro_rules! skip_outcome {
    ($vis:vis) => {
        /// What `PhaseCap::record_skip` decided for a full-cap skip.
        #[doc(hidden)]
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        $vis struct SkipOutcome {
            /// The cap is still full: skip the step. The thread's record then
            /// stands if `noted` is set (a tracked cap); an untracked cap
            /// records nothing.
            pub skip: bool,
            /// The thread must note the cap in its record mask (see
            /// `PhaseCap::record_skip`).
            pub noted: bool,
        }
    };
}
#[cfg(any(test, loom, feature = "test-utils"))]
skip_outcome!(pub);
#[cfg(not(any(test, loom, feature = "test-utils")))]
skip_outcome!(pub(crate));

/// Each distinct cap among `caps` once, in order (identity, not equality: two
/// caps with the same name and width are still two). The one dedup of the
/// run's caps: `PipelineBuilder::build` counts them against
/// [`MAX_PHASE_CAPS`], and `WakePlan::build` indexes them.
pub(crate) fn distinct_caps<'a>(caps: impl IntoIterator<Item = &'a PhaseCap>) -> Vec<&'a PhaseCap> {
    let mut distinct: Vec<&PhaseCap> = Vec::new();
    for cap in caps {
        if !distinct.iter().any(|d| std::ptr::eq(*d, cap)) {
            distinct.push(cap);
        }
    }
    distinct
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
            refused_slots: HolderSet::new(REFUSED_SLOTS),
            tracking: AtomicBool::new(false),
            waker: parking_lot::Mutex::new(None),
            plan_index: std::sync::atomic::AtomicUsize::new(usize::MAX),
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

    /// Bind this run's cap-release waker, replacing any earlier run's. `None`
    /// (a Legacy wake plan) turns refusal recording off: a refused worker then
    /// idles exactly as it did before cap-release wakes existed, and the cap
    /// leaves the run's record masks (its run index is cleared). A Directed
    /// run's index was set by `WakePlan::build`, which `Pipeline::run` calls
    /// just before this. Called by `Pipeline::run` before any thread starts.
    pub(crate) fn bind_waker(&self, waker: Option<CapWaker>) {
        self.refused_slots.clear_all();
        if waker.is_none() {
            self.set_plan_index(usize::MAX);
        }
        self.tracking.store(waker.is_some(), Ordering::Relaxed);
        *self.waker.lock() = waker;
    }

    /// Whether this run delivers cap-release wakes (a waker is bound).
    #[cfg(test)]
    pub(crate) fn is_tracking(&self) -> bool {
        self.tracking.load(Ordering::Relaxed)
    }

    /// [`Self::bind_waker`] for tests outside the crate.
    #[cfg(any(loom, feature = "test-utils"))]
    #[doc(hidden)]
    pub fn bind_waker_for_test(&self, waker: Option<CapWaker>) {
        self.bind_waker(waker);
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
    ///
    /// When the run delivers cap-release wakes, a refusal first records the
    /// calling thread's wake slot, fences and re-checks, pairing with the
    /// fence a permit's release issues before it wakes a recorded slot; the
    /// thread is told to idle where that wake reaches it. An admission adds one
    /// `Relaxed` load of the cap's run index and, when the cap has one (a
    /// Directed run tracks it), one thread-local read (a write only when the
    /// calling thread owes this cap a wake);
    /// a call that recorded and whose re-check then admitted it also withdraws
    /// the bit it set (one `fetch_and`).
    ///
    /// A record covers the dispatch pass that made it: the worker loop clears
    /// the thread's own bit at the start of every pass
    /// (`PhaseCap::clear_own_record`), and an admission withdraws only a bit the
    /// same call newly set. So a thread refused for one step and then admitted
    /// on this cap for another, later in the same pass, stays recorded while
    /// it works: a release may spend its wake on that running thread while
    /// another refused thread waits. That wait is bounded by the one item the
    /// running thread is working, since its own release follows and wakes the
    /// next recorded thread; it is not left to a timer.
    ///
    /// A wake that reaches a thread with no use for it — claimed in the window
    /// between its waking and its next pass start, or a full-cap skip's record
    /// for a step with no input — is not lost either: the thread's next pass
    /// start sees its bit claimed, and if that pass takes no permit here while
    /// one is still free, it forwards the wake to the next recorded thread
    /// (`PhaseCap::forward_wake`).
    #[inline]
    #[must_use]
    pub fn try_acquire(&self) -> Option<CapPermit<'_>> {
        self.acquire_or_record(crate::runtime::wake_slot::current_slot, true)
    }

    /// [`Self::try_acquire`] recording `slot` instead of the calling thread's
    /// wake slot, and leaving the thread's idle hint alone — for the loom model,
    /// whose threads have no pipeline slot.
    #[cfg(any(test, loom, feature = "test-utils"))]
    #[doc(hidden)]
    #[must_use]
    pub fn try_acquire_as(&self, slot: usize) -> Option<CapPermit<'_>> {
        self.acquire_or_record(|| Some(slot), false)
    }

    #[inline]
    fn acquire_or_record(
        &self,
        slot: impl FnOnce() -> Option<usize>,
        note_thread: bool,
    ) -> Option<CapPermit<'_>> {
        let active = &self.active.0;
        let mut slot = Some(slot);
        // The slot this call recorded, if it did, and whether that record set
        // the bit (only such a record is withdrawn if the call is admitted).
        let mut recorded: Option<(Option<usize>, bool)> = None;
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
                // First refusal of this call, with release wakes on: record the
                // slot, fence, and re-check once (`HolderSet::latch_slot`, the
                // queues' protocol). Pairs with the fence in `release`
                // (Dekker): either that release's claim sees this bit, or this
                // re-check sees the permit it freed.
                if let Some(slot) = slot.take()
                    && self.tracking.load(Ordering::Relaxed)
                {
                    let slot = slot();
                    let (newly_set, _) = self.refused_slots.latch_slot_reporting(slot, || {
                        cur = active.load(Ordering::Acquire);
                        true
                    });
                    recorded = Some((slot, newly_set));
                    continue;
                }
                if let Some((slot, _)) = recorded
                    && note_thread
                {
                    crate::runtime::wake_slot::note_cap_refused();
                    self.note_recorded(slot);
                }
                self.refused.0.fetch_add(1, Ordering::Relaxed);
                return None;
            }
            match active.compare_exchange_weak(cur, cur + 1, Ordering::AcqRel, Ordering::Acquire) {
                Ok(_) => {
                    // Admitted after this call recorded (its re-check saw a
                    // freed permit): withdraw the record, so a release spends
                    // its one wake on a thread that is really refused — but
                    // only if this call set the bit. A bit already set was set
                    // earlier in this pass, by a refusal of another step on
                    // this cap, and that refusal still stands (a record covers
                    // its pass; the worker loop clears it when the next pass
                    // starts). The anonymous bit is shared and stays. A release
                    // racing the withdrawal costs one spurious unpark at most.
                    if let Some((slot, true)) = recorded {
                        self.refused_slots.clear_slot(slot);
                    }
                    // A wake this thread owes the cap (a release claimed its
                    // record) is paid by this permit: its release wakes the
                    // next recorded thread.
                    if note_thread && let Some(i) = self.plan_bit() {
                        crate::runtime::wake_slot::note_cap_admitted(i);
                    }
                    self.note_peak(cur + 1);
                    return Some(CapPermit(self));
                }
                Err(observed) => cur = observed,
            }
        }
    }

    /// Set this cap's index among the run's caps (`WakePlan::build`).
    pub(crate) fn set_plan_index(&self, i: usize) {
        self.plan_index.store(i, std::sync::atomic::Ordering::Relaxed);
    }

    /// The bit of this cap in a worker's thread-local record masks: `Some`
    /// when a Directed run tracks the cap, its run index (below
    /// [`MAX_PHASE_CAPS`], which `WakePlan::build` enforces); `None` otherwise.
    #[inline]
    pub(crate) fn plan_bit(&self) -> Option<u32> {
        let i = self.plan_index.load(std::sync::atomic::Ordering::Relaxed);
        (i != usize::MAX)
            .then(|| u32::try_from(i).expect("a run index below MAX_PHASE_CAPS fits a u32"))
    }

    /// The calling thread's record on this cap, in slot `slot`, now stands:
    /// note it in the thread's mask, so its next pass start can tell whether a
    /// release claimed it. The anonymous bit is shared and not tracked.
    #[inline]
    fn note_recorded(&self, slot: Option<usize>) {
        let bit = if self.refused_slots.owns_bit(slot) { self.plan_bit() } else { None };
        crate::runtime::wake_slot::note_cap_recorded(bit);
    }

    /// A dispatch pass starts on the thread in wake slot `slot`: clear that
    /// slot's refusal record, so a record covers only the pass that made it
    /// (see [`Self::try_acquire`]). Returns whether the bit was still set: a
    /// thread that recorded here and finds it clear was woken by a release
    /// that claimed it, and owes that wake ([`Self::forward_wake`]). One
    /// `Relaxed` load, plus a `fetch_and` only when the bit is set; nothing on
    /// an untracked cap. The anonymous bit is shared, so it is never cleared
    /// here.
    #[inline]
    pub(crate) fn clear_own_record(&self, slot: Option<usize>) -> bool {
        self.tracking.load(Ordering::Relaxed) && self.refused_slots.clear_slot_reporting(slot)
    }

    /// [`Self::clear_own_record`] for the loom model, whose threads have no
    /// pipeline slot.
    #[cfg(any(test, loom, feature = "test-utils"))]
    #[doc(hidden)]
    pub fn clear_own_record_as(&self, slot: usize) -> bool {
        self.clear_own_record(Some(slot))
    }

    /// Pass on a wake this thread was handed and did not spend: a release
    /// claimed its record (one wake per freed permit), and its pass then took
    /// no permit here. If a permit is still free, wake the next recorded
    /// thread, so the permit is not left idle while a refused thread sleeps.
    /// Fences first (`SeqCst`): the thread saw the release's claim of its bit,
    /// so this fence follows the release's, and the claim below sees every
    /// record whose fence preceded it; a later record sees the free permit on
    /// its own re-check. Terminates: nothing records on a cap with a free
    /// permit (a refusal and a skip both need it full), so forwards only
    /// consume records. Nothing while a whole-cap reservation is pending (the
    /// whole holder's release wakes them).
    pub(crate) fn forward_wake(&self) {
        fence(Ordering::SeqCst);
        if self.tracking.load(Ordering::Relaxed) && self.has_free_permit() {
            self.wake_refused(1);
        }
    }

    /// [`Self::forward_wake`] for the loom model and tests.
    #[cfg(any(test, loom, feature = "test-utils"))]
    #[doc(hidden)]
    pub fn forward_wake_for_test(&self) {
        self.forward_wake();
    }

    /// The thread in wake slot `slot` is about to skip a step on this cap
    /// because the cap is full (a cap-parked worker's re-poll): record it as a
    /// refusal would — record, fence, re-check — so the release that frees a
    /// permit wakes it. `true`: still full, skip the step (the record stands
    /// for this pass). `false`: a permit is free (before or after the record),
    /// so poll the step; a bit this call set is withdrawn. An untracked cap
    /// only reports whether it is full.
    ///
    /// The record does not depend on the step's input: a step empty now may
    /// get input while the cap is still full, and a cap-parked worker is
    /// reached for that only through this record (the `Pool` route does not
    /// claim it for a full cap). Gating it on input would lose that wake to the
    /// timer. The cost of recording anyway is bounded: a release that spends
    /// its wake on a skipping thread with nothing to do costs one unpark and
    /// one pass, after which that thread forwards the wake
    /// ([`Self::forward_wake`]).
    ///
    /// `noted` in the outcome says the thread must note this cap in its
    /// record mask: its record stands (`skip`), or the re-check saw a permit
    /// freed by a release that had already claimed the bit this call set —
    /// the thread was handed that release's wake, and polls instead of
    /// skipping. If that poll takes no permit, the wake is owed, and only a
    /// noted record lets the next pass start see it and forward it. The
    /// anonymous bit is shared, so a claim of it is never reported.
    #[inline]
    pub(crate) fn record_skip(&self, slot: Option<usize>) -> SkipOutcome {
        if self.has_free_permit() {
            return SkipOutcome { skip: false, noted: false };
        }
        if !self.tracking.load(Ordering::Relaxed) {
            return SkipOutcome { skip: true, noted: false };
        }
        let (newly_set, room) =
            self.refused_slots.latch_slot_reporting(slot, || self.has_free_permit());
        if !room {
            return SkipOutcome { skip: true, noted: true };
        }
        // Room: poll the step. Withdraw the bit this call set; if it is gone
        // already, a release claimed it (only an owned bit can tell).
        let claimed = newly_set
            && self.refused_slots.owns_bit(slot)
            && !self.refused_slots.clear_slot_reporting(slot);
        SkipOutcome { skip: false, noted: claimed }
    }

    /// [`Self::record_skip`] for the calling thread's own wake slot, noting
    /// the record in its masks (the cap-parked re-poll's skip). Returns
    /// whether to skip.
    #[inline]
    pub(crate) fn record_skip_current(&self) -> bool {
        let slot = crate::runtime::wake_slot::current_slot();
        let outcome = self.record_skip(slot);
        if outcome.noted {
            self.note_recorded(slot);
        }
        outcome.skip
    }

    /// [`Self::record_skip`] for the loom model, whose threads have no
    /// pipeline slot.
    #[cfg(any(test, loom, feature = "test-utils"))]
    #[doc(hidden)]
    pub fn record_skip_as(&self, slot: usize) -> SkipOutcome {
        self.record_skip(Some(slot))
    }

    /// Claim up to `n` recorded refused slots and wake each through the run's
    /// waker. The caller has fenced (`SeqCst`) after freeing the permits.
    fn wake_refused(&self, n: usize) {
        let refused = &self.refused_slots;
        // Nothing recorded — the common case on a release: one `Relaxed` load
        // per 64 wake slots and no lock.
        if refused.is_empty() {
            return;
        }
        // Cloned out of the lock once per release, before any wake: the waker
        // unparks threads, which must not happen under the lock.
        let Some(waker) = self.waker.lock().clone() else { return };
        refused.take_up_to(n, &mut |bit| waker((!refused.is_anonymous(bit)).then_some(bit)));
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
    ///
    /// The same fence pairs with a refused worker's (see
    /// [`Self::try_acquire`]): with no reservation pending, the release
    /// claims one recorded refused worker per freed permit and wakes it, so a
    /// refused worker is never left waiting on its timer while a permit sits
    /// free: a claimed worker that has no use for the permit forwards the wake
    /// ([`Self::forward_wake`]). With a reservation pending it wakes none —
    /// they would only be refused again — and the whole-cap holder's release
    /// wakes them instead.
    #[inline]
    fn release(&self, n: usize) {
        let prev = self.active.0.fetch_sub(n, Ordering::AcqRel);
        fence(Ordering::SeqCst);
        let reserving = self.reserving.0.load(Ordering::Relaxed);
        if prev == n && reserving > 0 {
            self.parker.notify_all();
        }
        if reserving == 0 && self.tracking.load(Ordering::Relaxed) {
            self.wake_refused(n);
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

    /// Whether wake slot `slot`'s own refusal bit is set (tests).
    #[cfg(test)]
    pub(crate) fn is_recorded(&self, slot: usize) -> bool {
        self.refused_slots.is_set(Some(slot))
    }

    /// [`Self::acquire_whole`] callers currently waiting (tests).
    #[cfg(any(test, loom, feature = "test-utils"))]
    #[must_use]
    pub fn pending_reservations(&self) -> usize {
        self.reserving.0.load(Ordering::Relaxed)
    }

    /// Whether a pool step could be admitted right now: a permit is free and no
    /// whole-cap reservation is pending. Two `Relaxed` loads; a hint for
    /// routing a wake, not an admission (that is [`Self::try_acquire`]'s). The
    /// one spelling of the occupancy rule: the wake plan's `Pool` routing and a
    /// cap-parked worker's re-poll skip ([`Self::record_skip`]) both call it.
    #[inline]
    #[must_use]
    pub(crate) fn has_free_permit(&self) -> bool {
        self.active.0.load(Ordering::Relaxed) < self.max
            && self.reserving.0.load(Ordering::Relaxed) == 0
    }

    /// This cap as a shared handle (the one every step reporting it holds).
    pub(crate) fn arc(&self) -> Option<Arc<Self>> {
        self.this.upgrade()
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

/// Bind a waker to `cap` that records every slot it is asked to wake, and
/// return the record (tests).
#[cfg(test)]
pub(crate) fn bind_recording_waker(cap: &PhaseCap) -> Arc<parking_lot::Mutex<Vec<Option<usize>>>> {
    let woken: Arc<parking_lot::Mutex<Vec<Option<usize>>>> = Arc::default();
    let sink = Arc::clone(&woken);
    cap.bind_waker(Some(Arc::new(move |slot| sink.lock().push(slot))));
    woken
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

    /// A release wakes the refused slot it recorded, and only a slot that has
    /// its own bit is woken by number. A slot past the set (`--threads` above
    /// it), or none, is woken by broadcast: delivering it by number would
    /// unpark whichever thread owns that number and leave the refused one on
    /// its timer.
    #[rstest::rstest]
    #[case::slot_zero(Some(0), Some(0))]
    #[case::last_numbered_slot(Some(REFUSED_SLOTS - 1), Some(REFUSED_SLOTS - 1))]
    #[case::slot_on_the_anonymous_bit(Some(REFUSED_SLOTS), None)]
    #[case::slot_past_the_set(Some(1000), None)]
    #[case::no_slot(None, None)]
    fn a_release_wakes_the_refused_slot(
        #[case] refused: Option<usize>,
        #[case] expected: Option<usize>,
    ) {
        let cap = PhaseCap::new("t", 1);
        let woken = bind_recording_waker(&cap);
        let held = cap.try_acquire_as(0).expect("an idle cap admits");
        assert!(cap.acquire_or_record(|| refused, false).is_none(), "the cap is full");
        drop(held);
        assert_eq!(*woken.lock(), vec![expected]);
    }

    /// A refused thread whose post-fence re-check admits it withdraws its
    /// record: the next release spends its one wake on a thread that is still
    /// refused, not on the one already running. Slot 0 is refused and records
    /// itself, a permit frees between its refusal and its record (the latch
    /// hook), its re-check admits it; slot 1 is then refused for real, and the
    /// release of slot 0's permit must wake slot 1.
    #[test]
    fn a_refusal_admitted_by_its_recheck_leaves_no_record() {
        let cap = PhaseCap::new("t", 1);
        let woken = bind_recording_waker(&cap);
        std::mem::forget(cap.try_acquire_as(9).expect("an idle cap admits"));
        let releaser = Arc::clone(&cap);
        crate::queues::test_hooks::before_latch(move || releaser.release(1));
        let admitted = cap.try_acquire_as(0).expect("the re-check sees the freed permit");
        assert!(woken.lock().is_empty(), "the release ran before the record: nothing to wake");
        assert!(cap.try_acquire_as(1).is_none(), "slot 1 is refused for real");
        drop(admitted);
        assert_eq!(*woken.lock(), vec![Some(1)], "the release wakes slot 1, not slot 0");
    }

    /// A record covers the pass that made it: an admission on a later call
    /// leaves it, and the start of the thread's next pass clears it, so the
    /// next release wakes a thread that is really refused. Slots 1 and 0 are
    /// refused; the release wakes slot 0 and leaves slot 1's bit; slot 1 starts
    /// a new pass and is admitted; slot 2 is refused; the release of slot 1's
    /// permit must wake slot 2, not slot 1.
    #[test]
    fn a_new_pass_clears_a_stale_record() {
        let cap = PhaseCap::new("t", 1);
        let woken = bind_recording_waker(&cap);
        let holder = cap.try_acquire_as(9).expect("an idle cap admits");
        assert!(cap.try_acquire_as(1).is_none() && cap.try_acquire_as(0).is_none());
        drop(holder);
        assert_eq!(*woken.lock(), vec![Some(0)], "one wake per freed permit, lowest slot first");
        cap.clear_own_record_as(1);
        let one = cap.try_acquire_as(1).expect("slot 1 is admitted in its new pass");
        assert!(cap.try_acquire_as(2).is_none(), "slot 2 is really refused");
        drop(one);
        assert_eq!(*woken.lock(), vec![Some(0), Some(2)], "the release wakes slot 2, not slot 1");
    }

    /// Within one pass an admission does not clear a refusal: a thread refused
    /// for one step and admitted on the same cap for another, which then makes
    /// no progress, is still waiting for the first. Slots 1 and 0 are refused
    /// (slot 1 for its step X); the release wakes slot 0; slot 1 is admitted
    /// for its step Y, finds nothing and releases; that release must wake slot
    /// 1, whose next park then returns at once and re-polls X.
    #[test]
    fn a_refusal_earlier_in_the_pass_survives_a_later_admission() {
        let cap = PhaseCap::new("t", 1);
        let woken = bind_recording_waker(&cap);
        let holder = cap.try_acquire_as(9).expect("an idle cap admits");
        assert!(cap.try_acquire_as(1).is_none() && cap.try_acquire_as(0).is_none());
        drop(holder);
        assert_eq!(*woken.lock(), vec![Some(0)]);
        let y = cap.try_acquire_as(1).expect("slot 1 is admitted for step Y");
        drop(y); // Y made no progress
        assert_eq!(*woken.lock(), vec![Some(0), Some(1)], "slot 1 is still recorded for X");
    }

    /// A skip of a full cap records like a refusal: the release that frees a
    /// permit wakes the skipping thread. A skip whose re-check sees a freed
    /// permit polls the step instead and leaves no record.
    #[test]
    fn a_skip_of_a_full_cap_is_recorded_and_a_freed_one_is_polled() {
        let cap = PhaseCap::new("t", 1);
        let woken = bind_recording_waker(&cap);
        let holder = cap.try_acquire_as(9).expect("an idle cap admits");
        assert_eq!(
            cap.record_skip(Some(1)),
            SkipOutcome { skip: true, noted: true },
            "full: skip, recorded"
        );
        drop(holder);
        assert_eq!(*woken.lock(), vec![Some(1)], "the release wakes the skipping thread");

        std::mem::forget(cap.try_acquire_as(9).expect("an idle cap admits"));
        let releaser = Arc::clone(&cap);
        crate::queues::test_hooks::before_latch(move || releaser.release(1));
        assert_eq!(
            cap.record_skip(Some(2)),
            SkipOutcome { skip: false, noted: false },
            "the re-check sees the freed permit: poll"
        );
        assert!(!cap.is_recorded(2), "and the skip's record is withdrawn");
        assert_eq!(*woken.lock(), vec![Some(1)], "the release ran before the record");
    }

    /// A skip whose re-check sees a freed permit withdraws only a record it set
    /// itself: slot 2's earlier skip record (this pass) still stands. Slots 1
    /// and 2 skip the full cap; a permit frees between slot 2's second skip
    /// and its record (the release wakes slot 1, the lower); that skip polls,
    /// and slot 2 stays recorded.
    #[test]
    fn a_skip_that_polls_keeps_an_earlier_record() {
        let cap = PhaseCap::new("t", 1);
        let woken = bind_recording_waker(&cap);
        std::mem::forget(cap.try_acquire_as(9).expect("an idle cap admits"));
        assert!(cap.record_skip(Some(1)).skip && cap.record_skip(Some(2)).skip, "full: both skip");
        let releaser = Arc::clone(&cap);
        crate::queues::test_hooks::before_latch(move || releaser.release(1));
        assert_eq!(
            cap.record_skip(Some(2)),
            SkipOutcome { skip: false, noted: false },
            "the re-check sees the freed permit: poll"
        );
        assert_eq!(*woken.lock(), vec![Some(1)], "the release woke the lower slot");
        assert!(cap.is_recorded(2), "slot 2's earlier record stands");
    }

    /// A skip on the shared anonymous bit never reports a claim: the bit
    /// cannot be withdrawn, so "already gone" says nothing about a release.
    /// A permit frees before the record (the latch hook); the re-check polls.
    #[test]
    fn a_skip_on_the_anonymous_bit_reports_no_claim() {
        let cap = PhaseCap::new("t", 1);
        let _woken = bind_recording_waker(&cap);
        std::mem::forget(cap.try_acquire_as(9).expect("an idle cap admits"));
        let releaser = Arc::clone(&cap);
        crate::queues::test_hooks::before_latch(move || releaser.release(1));
        assert_eq!(
            cap.record_skip(Some(REFUSED_SLOTS + 45)),
            SkipOutcome { skip: false, noted: false }
        );
    }

    /// An untracked cap (a Legacy run) records no skip; the skip only reports
    /// whether the cap is full.
    #[test]
    fn an_untracked_cap_records_no_skip() {
        let cap = PhaseCap::new("t", 1);
        assert_eq!(
            cap.record_skip(Some(1)),
            SkipOutcome { skip: false, noted: false },
            "an idle cap is not skipped"
        );
        let held = cap.try_acquire_as(9).expect("an idle cap admits");
        assert_eq!(
            cap.record_skip(Some(1)),
            SkipOutcome { skip: true, noted: false },
            "a full cap is skipped, unrecorded"
        );
        assert!(cap.refused_slots.is_empty(), "nothing recorded");
        drop(held);
    }

    /// An admission whose re-check withdraws its own record leaves a record the
    /// same thread made earlier on this cap standing (that record set the bit,
    /// this one did not). Slot 3 is refused, then slot 1; a permit is freed
    /// between slot 3's second refusal and its record (the latch hook), and that
    /// release wakes slot 1; slot 3's re-check admits it, and its first record
    /// must still be there for the next release.
    #[test]
    fn a_re_check_admission_keeps_an_earlier_record() {
        let cap = PhaseCap::new("t", 1);
        let woken = bind_recording_waker(&cap);
        std::mem::forget(cap.try_acquire_as(9).expect("an idle cap admits"));
        assert!(cap.try_acquire_as(3).is_none() && cap.try_acquire_as(1).is_none());
        let releaser = Arc::clone(&cap);
        crate::queues::test_hooks::before_latch(move || releaser.release(1));
        let admitted = cap.try_acquire_as(3).expect("the re-check admits slot 3's second call");
        assert_eq!(*woken.lock(), vec![Some(1)], "the freed permit woke slot 1");
        drop(admitted);
        assert_eq!(*woken.lock(), vec![Some(1), Some(3)], "slot 3's first record stood");
    }

    /// A refusal records the thread and marks it for a timer park only when the
    /// run delivers cap-release wakes. A Legacy run binds no waker: its refused
    /// worker must stay unmarked (it idles on the event-count, where nothing
    /// would unpark it from a timer park) and nothing is recorded.
    #[rstest::rstest]
    #[case::legacy_no_waker(false)]
    #[case::directed_with_waker(true)]
    fn only_a_tracked_cap_marks_and_records_a_refused_thread(#[case] tracked: bool) {
        use crate::runtime::wake_slot::{SlotGuard, cap_refused_pending, take_cap_refused};
        let cap = PhaseCap::new("t", 1);
        cap.bind_waker(tracked.then(|| Arc::new(|_: Option<usize>| {}) as CapWaker));
        let held = cap.try_acquire_as(9).expect("an idle cap admits");
        let (marked, recorded) = thread::spawn({
            let cap = Arc::clone(&cap);
            move || {
                let _slot = SlotGuard::enter(2);
                let _ = take_cap_refused();
                assert!(cap.try_acquire().is_none(), "the cap is full");
                (cap_refused_pending(), !cap.refused_slots.is_empty())
            }
        })
        .join()
        .expect("refused thread");
        assert_eq!((marked, recorded), (tracked, tracked));
        drop(held);
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
