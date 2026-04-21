#![deny(unsafe_code)]

//! Sharded per-thread accumulator used to cap metric memory in parallel
//! pipeline stages.
//!
//! Pipeline commands historically collected per-group metric structs into an
//! unbounded `SegQueue`, retaining one entry per position/MI group for the
//! whole run. At real-data scale (hundreds of millions of records, tens of
//! millions of groups) the backing `AHashMap`s and `SegQueue` nodes grew into
//! tens of gigabytes even though downstream reduction only needed a handful of
//! counters (see issue #285).
//!
//! [`PerThreadAccumulator`] replaces that pattern with a fixed number of
//! [`Mutex<A>`] slots — one per worker thread. Each worker claims a slot on
//! first use and merges per-group results into it immediately, so retained
//! memory is `O(threads × distinct keys)` instead of `O(groups)`. After the
//! pipeline completes, callers fold the slots into the final metric output.
//!
//! The slot index is assigned lazily from a process-wide atomic counter
//! ([`SLOT_COUNTER`]) and cached in thread-local storage ([`THREAD_SLOT`]).
//! Indices are taken modulo the instance's slot count, so reused threads
//! across test runs, or long-lived threads that touch multiple accumulators,
//! remain correct — at worst two threads share a slot and merge under the same
//! mutex. Slot distribution is only uniform when `num_slots` is at least the
//! count of distinct global indices ever assigned; with a smaller `num_slots`
//! the modulo mapping deterministically clusters multiple threads onto the
//! same slot and adds lock contention, but does not affect correctness.

use parking_lot::Mutex;
use std::cell::Cell;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

static SLOT_COUNTER: AtomicUsize = AtomicUsize::new(0);
std::thread_local! {
    static THREAD_SLOT: Cell<Option<usize>> = const { Cell::new(None) };
}

/// Per-thread accumulator with a fixed number of mutex-protected slots.
///
/// See [module docs][self] for rationale and the memory-scaling contract.
#[derive(Debug)]
pub struct PerThreadAccumulator<A> {
    slots: Vec<Mutex<A>>,
}

impl<A: Send> PerThreadAccumulator<A> {
    /// Allocates `num_slots` slots seeded by `init` (called once per slot).
    /// `num_slots` is clamped to at least 1; callers typically pass the
    /// pipeline's worker count.
    ///
    /// Use this when `A` has no natural `Default` (e.g. an enum whose variant
    /// is only known at runtime). For types that do implement `Default`, use
    /// [`Self::new`].
    #[must_use]
    pub fn with_init<F>(num_slots: usize, mut init: F) -> Arc<Self>
    where
        F: FnMut() -> A,
    {
        let slots = (0..num_slots.max(1)).map(|_| Mutex::new(init())).collect();
        Arc::new(Self { slots })
    }

    /// Runs `f` against the calling thread's slot.
    ///
    /// The mutex is held only for the duration of `f`; do not perform blocking
    /// I/O inside the closure. Slot assignment is lazy and persisted in TLS,
    /// so repeated calls from the same thread contend on the same slot.
    #[inline]
    pub fn with_slot<F, R>(&self, f: F) -> R
    where
        F: FnOnce(&mut A) -> R,
    {
        let slot = global_thread_index() % self.slots.len();
        let mut guard = self.slots[slot].lock();
        f(&mut guard)
    }

    /// Borrows each slot in order. Callers typically reduce into a final
    /// aggregate after the pipeline has returned.
    ///
    /// Safe to call while other `Arc` holders exist, unlike [`Self::into_slots`].
    #[must_use]
    pub fn slots(&self) -> &[Mutex<A>] {
        &self.slots
    }
}

impl<A: Default + Send> PerThreadAccumulator<A> {
    /// Allocates `num_slots` default-initialized slots. `num_slots` is clamped
    /// to at least 1; callers typically pass the pipeline's worker count.
    #[must_use]
    pub fn new(num_slots: usize) -> Arc<Self> {
        Self::with_init(num_slots, A::default)
    }

    /// Consumes the accumulator and yields each slot's inner value.
    ///
    /// Requires unique `Arc` ownership — i.e. the pipeline closures that held
    /// clones have been dropped. Debug builds assert this invariant. If other
    /// `Arc` holders remain in release builds, falls back to `std::mem::take`
    /// under `A: Default`, which is lossy under concurrent use: any
    /// `with_slot` call racing from a surviving holder after the take will land
    /// in the freshly defaulted slot and be dropped when that holder releases
    /// its `Arc`. Callers must quiesce all writers before invoking.
    pub fn into_slots(self: Arc<Self>) -> Vec<A> {
        debug_assert_eq!(
            Arc::strong_count(&self),
            1,
            "into_slots called with outstanding Arc holders; fallback is lossy under concurrent writes",
        );
        match Arc::try_unwrap(self) {
            Ok(inner) => inner.slots.into_iter().map(Mutex::into_inner).collect(),
            Err(arc) => arc.slots.iter().map(|m| std::mem::take(&mut *m.lock())).collect(),
        }
    }
}

/// Returns a process-wide unique index for the calling thread, lazily
/// assigned from [`SLOT_COUNTER`] and cached in [`THREAD_SLOT`]. Callers map
/// the index into their own slot range with modulo.
#[inline]
fn global_thread_index() -> usize {
    THREAD_SLOT.with(|c| {
        c.get().unwrap_or_else(|| {
            let s = SLOT_COUNTER.fetch_add(1, Ordering::Relaxed);
            c.set(Some(s));
            s
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;

    #[derive(Default, Debug)]
    struct Counter(u64);

    #[test]
    fn single_thread_accumulates_in_one_slot() {
        let acc: Arc<PerThreadAccumulator<Counter>> = PerThreadAccumulator::new(4);
        for _ in 0..100 {
            acc.with_slot(|c| c.0 += 1);
        }
        let total: u64 = acc.slots().iter().map(|s| s.lock().0).sum();
        assert_eq!(total, 100);
    }

    #[test]
    fn parallel_threads_share_total_count() {
        let acc: Arc<PerThreadAccumulator<Counter>> = PerThreadAccumulator::new(8);
        thread::scope(|s| {
            for _ in 0..8 {
                let acc = Arc::clone(&acc);
                s.spawn(move || {
                    for _ in 0..1_000 {
                        acc.with_slot(|c| c.0 += 1);
                    }
                });
            }
        });
        let total: u64 = acc.slots().iter().map(|s| s.lock().0).sum();
        assert_eq!(total, 8_000);
    }

    #[test]
    fn num_slots_zero_clamped_to_one() {
        let acc: Arc<PerThreadAccumulator<Counter>> = PerThreadAccumulator::new(0);
        assert_eq!(acc.slots().len(), 1);
        acc.with_slot(|c| c.0 = 42);
        assert_eq!(acc.slots()[0].lock().0, 42);
    }

    #[test]
    fn into_slots_drains_under_unique_ownership() {
        let acc: Arc<PerThreadAccumulator<Counter>> = PerThreadAccumulator::new(4);
        acc.with_slot(|c| c.0 = 7);
        let slots = acc.into_slots();
        assert_eq!(slots.len(), 4);
        assert_eq!(slots.iter().map(|c| c.0).sum::<u64>(), 7);
    }
}
