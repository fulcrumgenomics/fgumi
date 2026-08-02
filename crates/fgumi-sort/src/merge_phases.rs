//! Per-component timing for the final k-way merge.
//!
//! [`SortPhaseTimer`](crate::external::SortPhaseTimer) reports the merge as one
//! bucket, and its own doc says why that is coarse: the span "includes reader
//! decompression + writer compression". That is enough to see the merge dominate
//! a spill-heavy sort, and not enough to say what to do about it. These counters
//! split it.
//!
//! # These are busy-time sums, not a partition of merge wall time
//!
//! The merge is a concurrent pipeline: worker threads read spill blocks from
//! disk, decompress them, and compress output blocks, while the main thread
//! drains the loser tree. All four happen *at the same time* on different
//! threads. So the four counters here overlap each other and routinely sum to
//! more than the merge's wall clock — on an 8-thread sort they can sum to
//! several times it.
//!
//! Read them as "how much thread-busy time went into each component", which is
//! the right question for a phase that is CPU-bound, and compare them against
//! each other rather than against wall clock. The ratio between them is the
//! signal; their sum is not meaningful.
//!
//! Each counter also carries an operation count, so a component's mean cost per
//! unit is available without assuming a block size.
//!
//! # Scope: worker threads only
//!
//! These cover the three things pool workers do during the merge. The main
//! thread's own work -- loser-tree pops, fetching the next record, enqueueing
//! writes -- is already measured by the sampled timers in
//! `merge_chunks_generic`, and is reported alongside these. Duplicating it here
//! would mean an `Instant::now()` pair per record rather than per block, which
//! on a billion-record sort is real overhead for a number we already have.

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

/// Nanosecond/operation counter pair for one merge component.
#[derive(Debug, Default)]
pub(crate) struct ComponentCounter {
    nanos: AtomicU64,
    ops: AtomicU64,
}

impl ComponentCounter {
    /// Run `f`, adding its elapsed time and one operation to this counter.
    ///
    /// `Relaxed` throughout: these are diagnostic totals with no ordering
    /// relationship to the data they describe, and the merge hot path should not
    /// pay for fences to maintain them.
    pub(crate) fn time<R>(&self, f: impl FnOnce() -> R) -> R {
        let start = Instant::now();
        let result = f();
        #[allow(clippy::cast_possible_truncation)]
        let elapsed = start.elapsed().as_nanos() as u64;
        self.nanos.fetch_add(elapsed, Ordering::Relaxed);
        self.ops.fetch_add(1, Ordering::Relaxed);
        result
    }

    #[allow(
        clippy::cast_precision_loss,
        reason = "nanosecond totals stay far below 2^52; a sort would need ~52 days of busy time to lose a nanosecond of precision here"
    )]
    fn snapshot(&self) -> (f64, u64) {
        let nanos = self.nanos.load(Ordering::Relaxed);
        let ops = self.ops.load(Ordering::Relaxed);
        (nanos as f64 / 1e9, ops)
    }
}

/// Busy-time counters for the four components of the k-way merge.
#[derive(Debug, Default)]
pub(crate) struct MergePhaseCounters {
    /// Pulling raw (still-compressed) blocks off disk into the per-file FIFO.
    pub(crate) read: ComponentCounter,
    /// Decompressing spill blocks (BGZF or zstd, whichever the spill used).
    pub(crate) decompress: ComponentCounter,
    /// Compressing output blocks (BGZF at `--compression-level`).
    pub(crate) compress: ComponentCounter,
}

/// A read-only view of the counters, for logging.
#[derive(Debug, Clone, Copy)]
pub struct MergePhaseBreakdown {
    /// Busy seconds spent reading raw spill blocks from disk, and block count.
    pub read: (f64, u64),
    /// Busy seconds spent decompressing spill blocks, and block count.
    pub decompress: (f64, u64),
    /// Busy seconds spent compressing output blocks, and block count.
    pub compress: (f64, u64),
}

impl MergePhaseCounters {
    pub(crate) fn snapshot(&self) -> MergePhaseBreakdown {
        MergePhaseBreakdown {
            read: self.read.snapshot(),
            decompress: self.decompress.snapshot(),
            compress: self.compress.snapshot(),
        }
    }
}

impl MergePhaseBreakdown {
    /// Total busy seconds across all four components.
    ///
    /// Deliberately *not* comparable to merge wall clock — see the module doc.
    /// Provided so each component can be expressed as a share of total merge
    /// effort, which is the comparison that means something.
    #[must_use]
    pub fn total_busy_secs(self) -> f64 {
        self.read.0 + self.decompress.0 + self.compress.0
    }

    /// Whether anything was recorded, so a sort that never merged stays silent.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.read.1 == 0 && self.decompress.1 == 0 && self.compress.1 == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use std::thread;
    use std::time::Duration;

    #[test]
    fn test_counter_accumulates_time_and_ops() {
        let counter = ComponentCounter::default();
        for _ in 0..3 {
            counter.time(|| thread::sleep(Duration::from_millis(5)));
        }
        let (secs, ops) = counter.snapshot();
        assert_eq!(ops, 3);
        assert!(secs >= 0.015, "expected >= 15ms accumulated, got {secs}");
    }

    #[test]
    fn test_counter_returns_the_closure_value() {
        let counter = ComponentCounter::default();
        assert_eq!(counter.time(|| 42), 42);
    }

    /// The counters are written from several worker threads at once, so the
    /// accumulation must not lose updates.
    #[test]
    fn test_counter_accumulates_across_threads() {
        let counter = Arc::new(ComponentCounter::default());
        thread::scope(|scope| {
            for _ in 0..4 {
                let counter = Arc::clone(&counter);
                scope.spawn(move || {
                    for _ in 0..25 {
                        counter.time(|| ());
                    }
                });
            }
        });
        assert_eq!(counter.snapshot().1, 100);
    }

    /// Busy time from concurrent threads exceeds the wall clock they ran in.
    /// This is the property the module doc warns about, pinned so nobody
    /// "fixes" the counters into a wall-clock partition.
    #[test]
    fn test_busy_time_may_exceed_wall_clock() {
        let counters = MergePhaseCounters::default();
        let start = Instant::now();
        thread::scope(|scope| {
            for _ in 0..4 {
                scope.spawn(|| {
                    counters.decompress.time(|| thread::sleep(Duration::from_millis(20)));
                });
            }
        });
        let wall = start.elapsed().as_secs_f64();
        let busy = counters.snapshot().total_busy_secs();
        assert!(busy > wall, "busy {busy} should exceed wall {wall} for concurrent work");
    }

    #[test]
    fn test_empty_breakdown_is_reported_as_empty() {
        assert!(MergePhaseCounters::default().snapshot().is_empty());
    }
}
