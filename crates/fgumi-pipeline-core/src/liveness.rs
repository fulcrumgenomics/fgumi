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
//! an uncontended increment on a line no other worker touches. The monitor sums
//! the slots, which it does once per poll interval (seconds), so the read side's
//! cost is irrelevant.
//!
//! The sum is not a synchronized snapshot — slots are read one at a time under
//! `Relaxed`, so the total may mix values from slightly different instants. That
//! is fine for the only question asked of it: a torn total still *changes* when
//! any worker makes progress, and monotonicity per slot means it can never
//! spuriously appear frozen while work is happening.

use std::sync::atomic::{AtomicU64, Ordering};

/// Cache-line padding. 128 rather than 64 because Apple Silicon and some `x86_64`
/// prefetchers pair adjacent 64-byte lines, so a 64-byte stride can still share
/// a coherence unit.
const CACHE_LINE_BYTES: usize = 128;

#[repr(align(128))]
struct PaddedCounter(AtomicU64);

const _: () = assert!(std::mem::size_of::<PaddedCounter>() == CACHE_LINE_BYTES);

/// Per-worker progress counters, summed by the deadlock monitor.
pub struct LivenessCounter {
    slots: Box<[PaddedCounter]>,
    /// `slots.len() - 1`, valid because the slot count is always a power of two.
    /// The bump path masks with this instead of `%`: a runtime modulo is an
    /// integer division (tens of cycles on both `aarch64` and `x86_64`) on a path
    /// that runs once per productive dispatch, which measured as a double-digit
    /// percentage of dispatch cost.
    mask: usize,
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
        Self { slots: (0..n).map(|_| PaddedCounter(AtomicU64::new(0))).collect(), mask: n - 1 }
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

    /// Each slot must sit on its own cache line — the whole reason for sharding.
    #[test]
    fn slots_are_cache_line_separated() {
        let counter = LivenessCounter::new(4);
        let a = std::ptr::from_ref(&counter.slots[0]) as usize;
        let b = std::ptr::from_ref(&counter.slots[1]) as usize;
        assert_eq!(
            b - a,
            CACHE_LINE_BYTES,
            "adjacent slots must be a full cache line apart or they false-share",
        );
    }
}
