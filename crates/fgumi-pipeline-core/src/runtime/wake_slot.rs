//! Pipeline thread identity for wake routing and hold timing.
//!
//! Every pipeline OS thread has a *wake slot*: pool worker `w` is slot `w`,
//! detached driver `d` is slot `n_threads + d` (the numbering the
//! `WorkerStateBoard` already uses). `Pipeline::run` enters the slot on the
//! spawned thread; a queue's reject path reads it to record *which* thread is
//! holding an item. Threads outside a pipeline loop have no slot.

use std::cell::Cell;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use crate::padded::Padded;
use crate::runtime::stats::PipelineStats;
use crate::topology::StepIdx;

thread_local! {
    static SLOT: Cell<Option<usize>> = const { Cell::new(None) };
}

/// Sets this OS thread's wake slot for its lifetime; clears it on drop.
pub(crate) struct SlotGuard(());

impl SlotGuard {
    /// Enter wake slot `slot` on the calling thread.
    pub(crate) fn enter(slot: usize) -> Self {
        SLOT.with(|s| s.set(Some(slot)));
        Self(())
    }
}

impl Drop for SlotGuard {
    fn drop(&mut self) {
        SLOT.with(|s| s.set(None));
    }
}

/// This OS thread's wake slot, if it is a pipeline thread.
#[inline]
pub(crate) fn current_slot() -> Option<usize> {
    SLOT.with(Cell::get)
}

/// Times held items on one byte-bounded edge (stats-on only). It keeps the
/// time of the first rejection of the hold in progress; the next successful
/// push into the edge (or a stash insert that accepts the item) ends the hold,
/// and rejections in between are retries.
///
/// The scope of a hold depends on how the producer runs. A `Parallel` step has
/// one clone per pool worker, each holding its own item and retrying it on its
/// own thread, so the clock keeps one stamp per thread slot
/// ([`Self::per_thread`]). Every other producer holds at most one item per
/// output branch, and a `Serial` step's held item may be retried by whichever
/// worker next wins the step's lock, so the clock keeps one stamp for the step
/// ([`Self::per_step`]): a flush by another worker ends the hold. One padded
/// word per stamp, so producers on different threads never share a line here.
pub struct HoldClock {
    step: StepIdx,
    stats: Arc<PipelineStats>,
    origin: Instant,
    /// Nanoseconds since `origin` (+1, so 0 means "not holding"). Per thread
    /// slot plus a last, shared entry for threads that have none; or a single
    /// entry for the whole step.
    since: Box<[Padded<AtomicU64>]>,
}

impl HoldClock {
    /// A clock that books holds to the `Parallel` producer `step` in `stats`,
    /// with one stamp per pipeline thread (`n_slots`) plus one shared stamp for
    /// threads without a slot.
    #[must_use]
    pub fn per_thread(step: StepIdx, stats: Arc<PipelineStats>, n_slots: usize) -> Self {
        Self::with_stamps(step, stats, n_slots + 1)
    }

    /// A clock that books holds to producer `step` in `stats` with one stamp
    /// for the step, whichever thread rejects or pushes: for every producer
    /// that is not `Parallel`.
    #[must_use]
    pub fn per_step(step: StepIdx, stats: Arc<PipelineStats>) -> Self {
        Self::with_stamps(step, stats, 1)
    }

    fn with_stamps(step: StepIdx, stats: Arc<PipelineStats>, n: usize) -> Self {
        Self {
            step,
            stats,
            origin: Instant::now(),
            since: (0..n).map(|_| Padded(AtomicU64::new(0))).collect(),
        }
    }

    /// The stamp of the calling thread (per-thread clock) or of the step.
    #[inline]
    fn slot(&self) -> &AtomicU64 {
        let last = self.since.len() - 1;
        if last == 0 {
            return &self.since[0].0;
        }
        let s = current_slot().filter(|&s| s < last).unwrap_or(last);
        &self.since[s].0
    }

    #[inline]
    fn now_ns(&self) -> u64 {
        u64::try_from(self.origin.elapsed().as_nanos()).unwrap_or(u64::MAX).saturating_add(1)
    }

    /// A push was rejected: the first rejection of a hold stamps it, later ones
    /// count a retry.
    #[inline]
    pub(crate) fn on_reject(&self) {
        let c = self.slot();
        if c.load(Ordering::Relaxed) == 0 {
            c.store(self.now_ns(), Ordering::Relaxed);
        } else {
            self.stats.record_held_retry(self.step);
        }
    }

    /// A push succeeded (or a stash accepted the item): if a hold was open, it
    /// ends now.
    #[inline]
    pub(crate) fn on_push(&self) {
        let c = self.slot();
        let t = c.load(Ordering::Relaxed);
        if t != 0 {
            c.store(0, Ordering::Relaxed);
            self.stats.record_hold(self.step, self.now_ns().saturating_sub(t));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slot_guard_sets_and_clears_the_thread_slot() {
        assert_eq!(current_slot(), None);
        let g = SlotGuard::enter(3);
        assert_eq!(current_slot(), Some(3));
        let other = std::thread::spawn(current_slot).join().unwrap();
        assert_eq!(other, None, "another thread never sees this thread's slot");
        drop(g);
        assert_eq!(current_slot(), None);
    }
}
