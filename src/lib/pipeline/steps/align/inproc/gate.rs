//! `CohortGate` — the in-process backend's cohort-granularity byte budget, and
//! the [`CohortLease`] that reserves one `-K` cohort against it.
//!
//! The in-process prepare step admits input one whole cohort at a time: it
//! reserves the cohort's bytes with a single non-blocking
//! [`CohortGate::try_acquire_cohort`] (a pool step must never block a worker
//! inside `try_run`), clones the returned lease into every work item of the
//! cohort, and the reservation is released exactly once, when the last clone
//! drops. The gate also drives the pool scheduler's refill signal (see
//! [`CohortGate::with_refill_signal`]).
//!
//! The subprocess backend's blocking writer/reader byte gate is the separate
//! [`InFlightGate`](super::super::InFlightGate); the two share only the
//! per-base byte estimate.

use std::any::Any;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, OnceLock};

use parking_lot::Mutex;

use super::super::IN_FLIGHT_BYTES_PER_BASE;

/// Upper bound (bases) a single template can add on top of a `-K` cohort. A
/// cohort closes only once its accumulated bases have reached `-K` *and* the
/// read count is even, so it can overshoot `-K` by at most the last template
/// admitted. A generous constant (covers long-read templates) since the derived
/// `bound` is only an over-estimate used to size the two-cohort byte budget.
const MAX_TEMPLATE_BASES: u64 = 1_000_000;

/// Number of `-K` cohorts the in-process backend keeps resident at once: one
/// pairing/merging while the next seeds, matching bwa's own `kt_pipeline`
/// double-buffering.
const IN_PROCESS_COHORTS_IN_FLIGHT: u64 = 2;

/// Byte reservation for one in-process `-K` cohort: `(chunk_size +
/// MAX_TEMPLATE_BASES) × IN_FLIGHT_BYTES_PER_BASE` — the cohort's `-K` bases
/// plus a one-template overshoot, priced at the per-base estimate. The
/// [`CohortGate`] budget is [`IN_PROCESS_COHORTS_IN_FLIGHT`] × this.
#[must_use]
pub(crate) fn cohort_bound_for_chunk_size(chunk_size_bases: u64) -> u64 {
    chunk_size_bases.saturating_add(MAX_TEMPLATE_BASES).saturating_mul(IN_FLIGHT_BYTES_PER_BASE)
}

/// The in-process backend's two-cohort byte budget for a given per-cohort
/// `bound` ([`cohort_bound_for_chunk_size`]). Replaces the subprocess
/// `4×`-slack formula (which existed for FASTQ-pipe buffering that the
/// in-process path does not have) with an explicit two-cohort budget of the
/// same magnitude.
#[must_use]
pub(crate) fn in_process_gate_budget(cohort_bound: u64) -> u64 {
    cohort_bound.saturating_mul(IN_PROCESS_COHORTS_IN_FLIGHT)
}

/// Non-blocking cohort-granularity byte gate for the in-process backend.
///
/// Admits a cohort when the gate is empty (an oversized cohort always passes an
/// empty gate, so a single cohort larger than the whole budget cannot deadlock)
/// or when its reservation fits under `budget`. Never blocks: a refused
/// admission returns `None` and the prepare step retries later.
pub(crate) struct CohortGate {
    budget: u64,
    inner: Mutex<CohortGateInner>,
    /// The pool scheduler's refill signal (see [`Self::with_refill_signal`]).
    /// `None` for a gate built with [`Self::new`] (tests).
    refill: Option<RefillSignal>,
}

struct CohortGateInner {
    in_flight: u64,
    /// The prepare step holds an open cohort it is still filling with input.
    filling: bool,
}

/// The gate's "refill wanted" signal, read by
/// [`RefillDrainScheduler`](crate::pipeline::core::runtime::RefillDrainScheduler).
struct RefillSignal {
    raised: Arc<AtomicBool>,
    /// One cohort's reservation, to tell whether another cohort would fit.
    cohort_bound: u64,
}

impl CohortGate {
    /// A gate over `budget` bytes with no refill signal.
    pub(crate) fn new(budget: u64) -> Self {
        Self {
            budget,
            inner: Mutex::new(CohortGateInner { in_flight: 0, filling: false }),
            refill: None,
        }
    }

    /// A gate that also drives a pool-scheduler refill signal. The signal is
    /// raised while the prepare step is filling an open cohort, or while the
    /// gate has room to admit another cohort: in both states the next input must
    /// be read and decoded, so upstream steps should run. It is lowered only
    /// when every admitted cohort has been fully read and the gate is full, when
    /// draining downstream work is all that helps. Without it, drain-first
    /// dispatch keeps every worker on the heavy downstream steps and never reads
    /// ahead, so seed/extend later waits on just-in-time decoding.
    pub(crate) fn with_refill_signal(budget: u64, cohort_bound: u64) -> Self {
        let raised = Arc::new(AtomicBool::new(true));
        Self { refill: Some(RefillSignal { raised, cohort_bound }), ..Self::new(budget) }
    }

    /// The refill signal, if this gate drives one.
    pub(crate) fn refill_signal(&self) -> Option<Arc<AtomicBool>> {
        self.refill.as_ref().map(|r| Arc::clone(&r.raised))
    }

    /// Recompute the refill signal from the gate state `g` (held under the lock,
    /// so transitions publish in order).
    fn update_refill(&self, g: &CohortGateInner) {
        if let Some(r) = &self.refill {
            let has_room =
                g.in_flight == 0 || g.in_flight.saturating_add(r.cohort_bound) <= self.budget;
            r.raised.store(g.filling || has_room, Ordering::Relaxed);
        }
    }

    /// The prepare step closed its open cohort: it no longer needs input for it.
    pub(crate) fn cohort_filled(&self) {
        let mut g = self.inner.lock();
        g.filling = false;
        self.update_refill(&g);
    }

    /// Reserve `bound` bytes for a whole cohort, returning a [`CohortLease`] on
    /// success and `None` (never blocking) on refusal.
    ///
    /// Admits when the gate is empty (`in_flight == 0`) OR when the reservation
    /// fits under budget (`in_flight + bound <= budget`). Non-blocking because
    /// the prepare step is a pool step and must never block a worker inside
    /// `try_run`.
    ///
    /// The returned lease is cloned into every work item of the cohort; the
    /// reservation is released exactly once, when the last clone drops
    /// (including on the error path). Takes `&Arc<Self>` so the lease can hold an
    /// owning handle to the gate for its release-on-drop.
    pub(crate) fn try_acquire_cohort(self: &Arc<Self>, bound: u64) -> Option<CohortLease> {
        let mut g = self.inner.lock();
        if g.in_flight == 0 || g.in_flight.saturating_add(bound) <= self.budget {
            g.in_flight = g.in_flight.saturating_add(bound);
            g.filling = true;
            self.update_refill(&g);
            drop(g);
            Some(CohortLease::new(Arc::clone(self), bound))
        } else {
            None
        }
    }

    /// Release a cohort's `bound` bytes: called exactly once per admitted cohort,
    /// from [`CohortLeaseInner::drop`].
    fn release(&self, bound: u64) {
        let mut g = self.inner.lock();
        g.in_flight = g.in_flight.saturating_sub(bound);
        self.update_refill(&g);
    }

    /// Current in-flight byte reservation. Test-only probe for the cohort-lease
    /// accounting (production code never reads it — admission uses the check
    /// inside the lock in [`Self::try_acquire_cohort`]).
    #[cfg(test)]
    pub(crate) fn in_flight_bytes(&self) -> u64 {
        self.inner.lock().in_flight
    }
}

/// An `Arc`-refcounted reservation of one in-process `-K` cohort's bytes against
/// the [`CohortGate`].
///
/// [`CohortGate::try_acquire_cohort`] mints one when a cohort opens; the prepare
/// step clones it into every `AlignWork` of that cohort, and each later stage
/// moves it forward. The reservation is released exactly once — when the *last*
/// clone drops (the inner `Arc` reaches zero refcount, running
/// [`CohortLeaseInner::drop`]). That covers both the normal path (pair/emit
/// drops each sub-batch's clone once the sub-batch is emitted) and the error
/// path (a step's `Drop` frees any items it still holds).
#[derive(Clone)]
pub(crate) struct CohortLease(Arc<CohortLeaseInner>);

/// The refcounted body of a [`CohortLease`]: it owns a handle to the gate and
/// the reserved byte count, and releases the reservation in its `Drop` — which
/// runs exactly once, when the last `CohortLease` clone (hence the last `Arc`)
/// is dropped.
struct CohortLeaseInner {
    gate: Arc<CohortGate>,
    bound: u64,
    /// The cohort's resident aligner state (the engine's `ResidentCohort`: every
    /// sub-batch's decoded reads and alignment regions, kept resident from
    /// seed/extend through pair/emit). Created by the first sub-batch that needs
    /// it and freed with the lease, i.e. exactly once, after the cohort's last
    /// work item drops. Type-erased so the lease does not depend on the engine;
    /// see [`CohortLease::resident`].
    resident: OnceLock<Box<dyn Any + Send + Sync>>,
}

impl CohortLease {
    /// Wrap a freshly-reserved gate handle. Called only from
    /// [`CohortGate::try_acquire_cohort`], which has already added `bound` to the
    /// gate's in-flight total; the matching subtraction happens in
    /// [`CohortLeaseInner::drop`].
    fn new(gate: Arc<CohortGate>, bound: u64) -> Self {
        Self(Arc::new(CohortLeaseInner { gate, bound, resident: OnceLock::new() }))
    }

    /// This cohort's resident state of type `T`, created by `init` on the first
    /// call from any clone of the lease (concurrent first calls race to `init`,
    /// and exactly one result is kept). Every later call, from any clone,
    /// returns the same value. It lives until the last clone of the lease drops.
    ///
    /// # Errors
    /// Returns an error if the slot already holds a value of a different type,
    /// which would mean two engines shared one lease.
    pub(crate) fn resident<T: Any + Send + Sync>(
        &self,
        init: impl FnOnce() -> T,
    ) -> std::io::Result<&T> {
        self.0.resident.get_or_init(|| Box::new(init())).downcast_ref::<T>().ok_or_else(|| {
            std::io::Error::other(format!(
                "cohort lease already holds resident state of a type other than {}",
                std::any::type_name::<T>()
            ))
        })
    }

    /// A self-contained lease over a throwaway, unbounded gate, for tests that
    /// need to populate a work item's `lease` field without wiring a real gate.
    /// It owns its gate, so dropping it is a self-contained no-op.
    #[cfg(test)]
    pub(crate) fn for_test() -> Self {
        Arc::new(CohortGate::new(u64::MAX))
            .try_acquire_cohort(0)
            .expect("empty gate admits a test lease")
    }
}

impl std::fmt::Debug for CohortLease {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CohortLease")
            .field("bound", &self.0.bound)
            .field("clones", &Arc::strong_count(&self.0))
            .finish()
    }
}

impl Drop for CohortLeaseInner {
    fn drop(&mut self) {
        self.gate.release(self.bound);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The gate's refill signal across a cohort's life, with room for two
    /// cohorts: raised while a cohort is open or another would fit, lowered only
    /// when both admitted cohorts are filled and the gate is full.
    #[test]
    fn refill_signal_follows_cohort_admission_and_release() {
        let gate = Arc::new(CohortGate::with_refill_signal(200, 100));
        let signal = gate.refill_signal().expect("the gate has a signal");
        let raised = || signal.load(Ordering::Relaxed);
        assert!(raised(), "an empty gate wants input");

        let first = gate.try_acquire_cohort(100).expect("first cohort fits");
        assert!(raised(), "filling the first cohort");
        gate.cohort_filled();
        assert!(raised(), "room for a second cohort");

        let second = gate.try_acquire_cohort(100).expect("second cohort fits");
        assert!(raised(), "filling the second cohort");
        gate.cohort_filled();
        assert!(!raised(), "both cohorts filled and the gate full");

        drop(first);
        assert!(raised(), "a released cohort makes room again");
        drop(second);
        assert!(raised(), "an empty gate wants input");
    }

    #[test]
    fn plain_gate_has_no_refill_signal() {
        assert!(CohortGate::new(100).refill_signal().is_none());
    }

    #[test]
    fn try_acquire_cohort_oversized_passes_empty_gate() {
        // A cohort whose `bound` exceeds the whole budget still admits when the
        // gate is empty — it can't be split, so refusing it would deadlock.
        let gate = Arc::new(CohortGate::new(100));
        let lease = gate.try_acquire_cohort(1000).expect("oversized cohort admits an empty gate");
        assert_eq!(gate.in_flight_bytes(), 1000, "the full oversized bound is reserved");
        drop(lease);
        assert_eq!(gate.in_flight_bytes(), 0, "dropping the lease releases the reservation");
    }

    #[test]
    fn try_acquire_cohort_admits_second_only_under_budget() {
        // Budget = two cohorts. The first admits (empty gate); a second admits
        // only while `in_flight + bound <= budget`; a third is refused until one
        // releases.
        let bound = 100;
        let gate = Arc::new(CohortGate::new(2 * bound));
        let first = gate.try_acquire_cohort(bound).expect("first cohort admits");
        assert_eq!(gate.in_flight_bytes(), bound);
        let second = gate.try_acquire_cohort(bound).expect("second cohort fits the budget");
        assert_eq!(gate.in_flight_bytes(), 2 * bound);
        // Gate is now full and non-empty: a third is refused, no blocking.
        assert!(gate.try_acquire_cohort(bound).is_none(), "third cohort exceeds budget → None");

        // Releasing the first frees exactly one cohort's worth, admitting a third.
        drop(first);
        assert_eq!(gate.in_flight_bytes(), bound);
        let third = gate.try_acquire_cohort(bound).expect("third admits after a release");
        assert_eq!(gate.in_flight_bytes(), 2 * bound);
        drop(second);
        drop(third);
        assert_eq!(gate.in_flight_bytes(), 0, "all reservations released");
    }

    #[test]
    fn cohort_lease_releases_only_on_last_clone_drop() {
        // The reservation is released exactly once — when the LAST clone drops —
        // mirroring a cohort whose sub-batches (each holding a clone) drain one
        // by one; only the final sub-batch's drop frees the gate.
        let bound = 100;
        let gate = Arc::new(CohortGate::new(2 * bound));
        let lease = gate.try_acquire_cohort(bound).expect("admit");
        let clone_a = lease.clone();
        let clone_b = lease.clone();
        assert_eq!(gate.in_flight_bytes(), bound, "clones share one reservation, not three");

        drop(lease);
        assert_eq!(gate.in_flight_bytes(), bound, "still reserved: two clones alive");
        drop(clone_a);
        assert_eq!(gate.in_flight_bytes(), bound, "still reserved: one clone alive");
        drop(clone_b);
        assert_eq!(gate.in_flight_bytes(), 0, "last clone drop releases the reservation");
    }

    #[test]
    fn cohort_lease_releases_on_error_path_unwind() {
        // On an error/panic the pipeline drops the in-flight items (and their
        // lease clones) while unwinding; the reservation must still be freed so a
        // stalled Prepare can admit the next cohort. `catch_unwind` stands in for
        // the framework's error teardown.
        let bound = 100;
        let gate = Arc::new(CohortGate::new(2 * bound));
        let lease = gate.try_acquire_cohort(bound).expect("admit");
        let clone = lease.clone();
        drop(lease);
        assert_eq!(gate.in_flight_bytes(), bound, "one clone still alive");

        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _held = clone; // moved in; dropped as the closure unwinds
            panic!("simulated step error");
        }));
        assert!(result.is_err(), "the panic propagates");
        assert_eq!(gate.in_flight_bytes(), 0, "the lease clone was freed during unwind");
    }

    /// Counts drops of the resident state it stands in for.
    struct DropCounted(Arc<std::sync::atomic::AtomicUsize>);

    impl Drop for DropCounted {
        fn drop(&mut self) {
            self.0.fetch_add(1, Ordering::SeqCst);
        }
    }

    /// The cohort's resident state is created once, shared by every clone of
    /// the lease, and freed exactly once: when the last clone drops, together
    /// with the gate reservation, whether that clone drops normally or while a
    /// failing step unwinds.
    #[rstest::rstest]
    #[case::normal_drop(false)]
    #[case::unwind_drop(true)]
    fn cohort_lease_frees_resident_state_once_on_last_clone_drop(#[case] unwind: bool) {
        use std::sync::atomic::AtomicUsize;

        let bound = 100;
        let gate = Arc::new(CohortGate::new(2 * bound));
        let lease = gate.try_acquire_cohort(bound).expect("admit");
        let drops = Arc::new(AtomicUsize::new(0));
        let inits = AtomicUsize::new(0);
        let init = || {
            inits.fetch_add(1, Ordering::SeqCst);
            DropCounted(Arc::clone(&drops))
        };

        let clone = lease.clone();
        let first: *const DropCounted = lease.resident(init).expect("same type");
        let second: *const DropCounted = clone.resident(init).expect("same type");
        assert!(std::ptr::eq(first, second), "every clone sees the one resident state");
        assert_eq!(inits.load(Ordering::SeqCst), 1, "created once");
        assert!(clone.resident(|| 0u8).is_err(), "a different type is refused, not replaced");

        drop(lease);
        assert_eq!(drops.load(Ordering::SeqCst), 0, "one clone still alive");
        if unwind {
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let _held = clone;
                panic!("simulated step error");
            }));
            assert!(result.is_err(), "the panic propagates");
        } else {
            drop(clone);
        }
        assert_eq!(drops.load(Ordering::SeqCst), 1, "freed exactly once, on the last drop");
        assert_eq!(gate.in_flight_bytes(), 0, "released together with the reservation");
    }
}
