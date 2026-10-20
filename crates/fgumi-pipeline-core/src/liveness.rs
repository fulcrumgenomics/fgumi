//! A cheap, always-on liveness signal for the deadlock monitor.
//!
//! # Why this exists separately from [`PipelineStats`](crate::runtime::stats::PipelineStats)
//!
//! The deadlock monitor needs to answer exactly one question: *has anything in
//! the pipeline made progress since I last looked?* It compares one total
//! against the previous total and cares about nothing else.
//!
//! It used to get that total from `PipelineStats`, which coupled the monitor to
//! full instrumentation — and full instrumentation is deliberately not free.
//! `dispatch_one_step` times each dispatch with `Instant::now()` (~20–50 ns on
//! Apple Silicon, ~50–100 ns on `x86_64`) and gates that on `stats.is_some()`
//! precisely to keep the uninstrumented path zero-cost. So the monitor could
//! only be armed by paying a per-dispatch timing cost on every run, which is why
//! it shipped disarmed by default and every wedge hung silently instead.
//!
//! This type breaks that coupling: liveness is a counter, profiling is a
//! separate opt-in. The monitor can therefore run with no `PipelineStats`
//! attached at all — a stats handle now only adds the per-step snapshot to the
//! stall report, and its absence costs the diagnostic, not the detection.
//!
//! That removes the *cost* argument for shipping disarmed; it does not by
//! itself arm anything. `PipelineConfig::default()` still sets
//! `deadlock_timeout_secs: 0`, because an armed monitor additionally requires
//! every output transport to be `ByteBounded` (see
//! [`PipelineError::MonitorBlindTransport`](crate::signal::PipelineError)) and
//! chains such as `Process2` use `CountBounded` today.
//!
//! Both of those are properties of the **scheduled** path. A run that fuses to
//! a single thread returns before the transport check and before the monitor
//! spawns, so it neither gains the monitor nor is rejected for a byte-blind
//! edge; it is bounded by the fused path's own stall budget instead.
//!
//! # Why the counter is sharded
//!
//! A single shared `AtomicU64` bumped by every worker on every productive
//! dispatch is a false-sharing hotspot: each increment is a read-modify-write on
//! one cache line, so N workers ping-pong that line between cores and the cost
//! grows with thread count — exactly the wrong shape, since more threads is when
//! the monitor matters most.
//!
//! Each worker therefore owns a slot padded to its own cache line, so a bump is
//! normally an uncontended increment on a line no other worker touches. (The one
//! exception is a dedicated driver thread, which reuses `worker_slot` 0 and so
//! shares slot 0 with pool worker 0; the bump is an atomic `fetch_add`, so a
//! coincident increment is counted correctly rather than lost — it only costs the
//! two threads a shared cache line, not accuracy.) The monitor sums
//! the slots, which it does once per poll interval (seconds), so the read side's
//! cost is irrelevant.
//!
//! The sum is not a synchronized snapshot — slots are read one at a time under
//! `Relaxed`, so the total may mix values from slightly different instants. That
//! is fine for the only question asked of it: a torn total still *changes* when
//! any worker makes progress, and monotonicity per slot means it can never
//! spuriously appear frozen while work is happening.
//!
//! # External waits
//!
//! The counter is bumped when a dispatch *returns* with progress, so a step
//! blocked inside one dispatch on something outside the pipeline — a
//! `write(2)` to a full pipe whose reader is still loading an index, a spill
//! write stuck on a slow parallel filesystem — looks exactly like a wedge:
//! no progress, and work queued behind it. Such a step is waiting on the
//! outside world, not deadlocked, and no timeout is long enough to tell the
//! two apart from progress alone (issue #1068).
//!
//! Steps therefore bracket every call that may block on an external resource
//! with [`external_wait`]. While one is in flight, the counter's
//! [`LivenessCounter::external_waits`] is non-zero and the monitor treats the
//! stall as a wait rather than a wedge. The bracket finds the running
//! pipeline's counter through a thread-local the worker loop installs
//! (`LivenessCounter::enter_thread`), so a step needs no handle to the
//! runtime; outside a pipeline worker it is a plain call.

use std::cell::RefCell;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use crate::padded::Padded;

/// One slot per worker, each alone on its line (see [`Padded`]).
type PaddedCounter = Padded<AtomicU64>;

/// Per-worker progress counters, summed by the deadlock monitor.
pub struct LivenessCounter {
    slots: Box<[PaddedCounter]>,
    /// `slots.len() - 1`, valid because the slot count is always a power of two.
    /// The bump path masks with this instead of `%`: a runtime modulo is an
    /// integer division (tens of cycles on both `aarch64` and `x86_64`) on a path
    /// that runs once per productive dispatch, which measured as a double-digit
    /// percentage of dispatch cost.
    mask: usize,
    /// Calls currently blocked inside [`external_wait`] on this pipeline's
    /// threads. Shared (rather than inline) so a worker's thread-local can hold
    /// it without borrowing the counter.
    external_waits: Arc<AtomicUsize>,
}

thread_local! {
    /// The external-wait gauge of the pipeline whose worker loop runs on this
    /// thread, installed by [`LivenessCounter::enter_thread`].
    static EXTERNAL_WAITS: RefCell<Option<Arc<AtomicUsize>>> = const { RefCell::new(None) };
}

/// Run `f`, a call that may block on a resource outside the pipeline — a read
/// or write on a pipe, file or socket — and count it as an external wait for
/// as long as it runs, so the deadlock monitor does not mistake a slow reader
/// or filesystem for a wedge (see the [module docs](self#external-waits)).
///
/// Wrap only calls that block on the outside world. A wait on another step of
/// the same pipeline (a condvar, a channel fed by a sibling step) must not be
/// wrapped: that wait *can* be part of a wedge, and wrapping it would hide it.
///
/// Outside a pipeline worker thread this is just `f()`.
pub fn external_wait<R>(f: impl FnOnce() -> R) -> R {
    // Decrement on unwind too, so a panicking call cannot pin the gauge.
    struct Leave(Arc<AtomicUsize>);
    impl Drop for Leave {
        fn drop(&mut self) {
            self.0.fetch_sub(1, Ordering::Relaxed);
        }
    }
    let gauge = EXTERNAL_WAITS.with(|g| g.borrow().clone());
    let Some(gauge) = gauge else {
        return f();
    };
    gauge.fetch_add(1, Ordering::Relaxed);
    let _leave = Leave(gauge);
    f()
}

/// Installs a pipeline's external-wait gauge on the current thread for its
/// lifetime and restores the previous one on drop. See
/// [`LivenessCounter::enter_thread`].
#[must_use = "the gauge is uninstalled when the scope is dropped"]
pub(crate) struct ThreadScope {
    previous: Option<Arc<AtomicUsize>>,
}

impl Drop for ThreadScope {
    fn drop(&mut self) {
        let previous = self.previous.take();
        EXTERNAL_WAITS.with(|g| *g.borrow_mut() = previous);
    }
}

impl LivenessCounter {
    /// One slot per worker. `n_workers` must cover every index passed to
    /// [`Self::bump`]; indices are taken modulo the slot count so an unexpected
    /// worker id degrades to sharing a slot rather than panicking on the hot
    /// path.
    #[must_use]
    pub fn new(n_workers: usize) -> Self {
        // Round up to a power of two so the bump path can mask instead of
        // dividing. The slack is a few unused cache lines — irrelevant next to
        // removing a division from a per-dispatch path.
        let n = n_workers.max(1).next_power_of_two();
        Self {
            slots: (0..n).map(|_| Padded(AtomicU64::new(0))).collect(),
            mask: n - 1,
            external_waits: Arc::new(AtomicUsize::new(0)),
        }
    }

    /// Route [`external_wait`] calls made on the current thread to this
    /// counter until the returned scope is dropped. Called once at the top of
    /// every worker loop; the previous gauge (if any) is restored on drop, so
    /// the single-worker path, which runs on the caller's thread, leaves it as
    /// it found it.
    pub(crate) fn enter_thread(&self) -> ThreadScope {
        let gauge = Some(Arc::clone(&self.external_waits));
        ThreadScope { previous: EXTERNAL_WAITS.with(|g| g.replace(gauge)) }
    }

    /// Number of [`external_wait`] calls currently blocked on this pipeline's
    /// threads. Read by the monitor once per poll.
    #[must_use]
    pub fn external_waits(&self) -> usize {
        self.external_waits.load(Ordering::Relaxed)
    }

    /// Record one unit of progress for `worker`.
    ///
    /// `Relaxed` is sufficient: the value only has to change, and the monitor
    /// re-reads it on a multi-second cadence, so no ordering relative to other
    /// memory is required.
    #[inline]
    pub fn bump(&self, worker: usize) {
        // Mask, not modulo — see `mask`. The slot count is a power of two, so
        // this is exact for in-range ids and wraps harmlessly for anything else.
        let slot = &self.slots[worker & self.mask];
        slot.0.fetch_add(1, Ordering::Relaxed);
    }

    /// Sum every worker's progress. Called once per monitor poll.
    #[must_use]
    pub fn total(&self) -> u64 {
        self.slots.iter().map(|s| s.0.load(Ordering::Relaxed)).sum()
    }

    /// Number of slots. This is the requested worker count rounded up to a
    /// power of two, not the requested count itself.
    #[must_use]
    pub fn len(&self) -> usize {
        self.slots.len()
    }

    /// Always false — `new` clamps to at least one slot. Present because clippy
    /// requires it alongside `len`.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        false
    }
}

impl std::fmt::Debug for LivenessCounter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Hand-written rather than derived: the per-slot atomics are an
        // implementation detail, and the only useful facts are how many slots
        // there are and what they currently sum to.
        f.debug_struct("LivenessCounter")
            .field("slots", &self.slots.len())
            .field("mask", &self.mask)
            .field("total", &self.total())
            .field("external_waits", &self.external_waits())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn total_sums_every_slot() {
        let counter = LivenessCounter::new(4);
        counter.bump(0);
        counter.bump(1);
        counter.bump(1);
        counter.bump(3);
        assert_eq!(counter.total(), 4);
    }

    #[test]
    fn an_out_of_range_worker_wraps_rather_than_panicking() {
        // The hot path must never panic on an unexpected worker id; sharing a
        // slot only costs contention, and the total stays correct.
        let counter = LivenessCounter::new(2);
        counter.bump(7);
        assert_eq!(counter.total(), 1);
    }

    #[test]
    fn a_single_worker_pipeline_still_gets_a_slot() {
        let counter = LivenessCounter::new(0);
        assert_eq!(counter.len(), 1, "the slot count is clamped to at least one");
        counter.bump(0);
        assert_eq!(counter.total(), 1);
    }

    /// Every increment must be observed — the point of the type is that a frozen
    /// total means a frozen pipeline, so a lost bump would be a false wedge
    /// report.
    #[test]
    fn concurrent_bumps_are_all_counted() {
        const WORKERS: usize = 8;
        const PER_WORKER: u64 = 10_000;

        let counter = Arc::new(LivenessCounter::new(WORKERS));
        let handles: Vec<_> = (0..WORKERS)
            .map(|w| {
                let counter = Arc::clone(&counter);
                std::thread::spawn(move || {
                    for _ in 0..PER_WORKER {
                        counter.bump(w);
                    }
                })
            })
            .collect();
        for h in handles {
            h.join().expect("worker thread joins");
        }
        assert_eq!(counter.total(), WORKERS as u64 * PER_WORKER);
    }

    /// Inside an entered thread, an external wait is visible to the monitor for
    /// exactly as long as the call runs.
    #[test]
    fn an_external_wait_is_counted_only_while_it_runs() {
        let counter = LivenessCounter::new(1);
        let _scope = counter.enter_thread();
        let during = external_wait(|| counter.external_waits());
        assert_eq!(during, 1, "the wait is counted while it runs");
        assert_eq!(counter.external_waits(), 0, "and released when it returns");
    }

    /// A call that panics must still release the gauge, or the monitor would
    /// treat every later stall as an external wait.
    #[test]
    fn a_panicking_external_wait_releases_the_gauge() {
        let counter = LivenessCounter::new(1);
        let _scope = counter.enter_thread();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            external_wait(|| panic!("intentional"));
        }));
        assert!(result.is_err());
        assert_eq!(counter.external_waits(), 0);
    }

    /// Off a pipeline thread the bracket is a plain call; once a scope is
    /// dropped the thread no longer reports to that pipeline.
    #[test]
    fn external_wait_outside_an_entered_thread_is_a_plain_call() {
        let counter = LivenessCounter::new(1);
        drop(counter.enter_thread());
        assert_eq!(external_wait(|| 7), 7);
        assert_eq!(counter.external_waits(), 0, "a dropped scope no longer routes waits");
    }

    /// Scopes nest: the inner pipeline's waits go to the inner counter, and
    /// dropping the inner scope restores the outer one.
    #[test]
    fn nested_scopes_restore_the_outer_gauge() {
        let outer = LivenessCounter::new(1);
        let inner = LivenessCounter::new(1);
        let _outer_scope = outer.enter_thread();
        {
            let _inner_scope = inner.enter_thread();
            assert_eq!(external_wait(|| (outer.external_waits(), inner.external_waits())), (0, 1));
        }
        assert_eq!(external_wait(|| (outer.external_waits(), inner.external_waits())), (1, 0));
    }

    /// Each slot must sit on its own cache line — the whole reason for sharding.
    #[test]
    fn slots_are_cache_line_separated() {
        let counter = LivenessCounter::new(4);
        let a = std::ptr::from_ref(&counter.slots[0]) as usize;
        let b = std::ptr::from_ref(&counter.slots[1]) as usize;
        assert_eq!(
            b - a,
            128,
            "adjacent slots must be a full cache line apart or they false-share",
        );
    }
}
