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

use crate::merge_trace::{DurationHistogram, HistogramReport};

/// Nanosecond/operation counter pair for one merge component.
///
/// Backed by a [`DurationHistogram`], which already keeps an exact count and
/// sum. Sharing the storage rather than incrementing a second set of atomics
/// keeps the totals in this report and the distributions in
/// [`crate::merge_trace`] describing literally the same observations — two
/// parallel counters over one event are a divergence waiting to happen, and the
/// divergence would be silent.
#[derive(Debug, Default)]
pub(crate) struct ComponentCounter {
    hist: DurationHistogram,
}

impl ComponentCounter {
    /// Run `f`, adding its elapsed time and one operation to this counter.
    ///
    /// `Relaxed` throughout: these are diagnostic totals with no ordering
    /// relationship to the data they describe, and the merge hot path should not
    /// pay for fences to maintain them.
    pub(crate) fn time<R>(&self, f: impl FnOnce() -> R) -> R {
        self.hist.time(f)
    }

    /// The distribution behind the totals, for [`crate::merge_trace`].
    pub(crate) fn histogram(&self) -> HistogramReport {
        self.hist.snapshot()
    }

    fn snapshot(&self) -> (f64, u64) {
        let report = self.hist.snapshot();
        (report.total_secs(), report.count)
    }
}

/// Busy-time counters for the four components of the k-way merge.
#[derive(Debug, Default)]
pub(crate) struct MergePhaseCounters {
    /// Pulling raw (still-compressed) blocks off disk into the per-file FIFO.
    pub(crate) read: ComponentCounter,
    /// Decompressing spill blocks (BGZF or zstd, whichever the spill used).
    pub(crate) decompress: ComponentCounter,
    /// Compressing OUTPUT blocks (BGZF at `--compression-level`), during the
    /// merge.
    pub(crate) output_compress: ComponentCounter,
    /// Compressing SPILL blocks during Phase 1. Counted separately because the
    /// same worker step serves both targets: lumping them together attributes
    /// Phase 1's spill compression to the merge, which overstates the merge's
    /// compression cost by however much was spilled.
    pub(crate) spill_compress: ComponentCounter,
}

/// Latency distributions for the four stages a block passes through, plus the
/// writer, plus how much scanning was wasted getting there.
///
/// [`MergePhaseCounters`] already gives each stage a busy total and a block
/// count, which yields a *mean*. That is not enough to settle an argument: a
/// stage averaging 187 us could be uniform or bimodal with a long tail, and only
/// the tail explains why a worker is unavailable when the consumer needs one.
/// Output compression is 69% of all worker busy on the measured cell and had no
/// distribution at all; the writer had none either, so a jump in consumer
/// backpressure from 0.0s to 29.9s could not be attributed.
#[derive(Debug, Default)]
pub(crate) struct StageLatency {
    /// One batched read of compressed blocks from a spill file.
    pub(crate) read: crate::merge_trace::DurationHistogram,
    /// Decompressing one spill block.
    pub(crate) decompress: crate::merge_trace::DurationHistogram,
    /// Compressing one output block.
    pub(crate) output_compress: crate::merge_trace::DurationHistogram,
    /// Compressing one Phase 1 spill block.
    pub(crate) spill_compress: crate::merge_trace::DurationHistogram,
    // The writer's own histograms live on `PermitPool`: the I/O writer thread
    // already holds that `Arc`, and one pool per writer keeps the output
    // writer's stats separate from a spill writer's.
    /// Files a worker passed over before one gave it work, summed over all
    /// productive scans.
    ///
    /// The scan tally is published only when a scan finds *nothing*, so the
    /// files skipped on the way to a successful claim were invisible. On an
    /// 89-way merge a worker can walk most of the file set before finding work,
    /// and that walk is the cost this counts.
    pub(crate) wasted_visits: std::sync::atomic::AtomicU64,
    /// Scans that ended in a claim. The denominator for `wasted_visits`.
    pub(crate) useful_claims: std::sync::atomic::AtomicU64,
}

impl StageLatency {
    /// Record a productive scan that skipped `skipped` files before claiming.
    pub(crate) fn record_claim(&self, skipped: u64) {
        use std::sync::atomic::Ordering;
        self.wasted_visits.fetch_add(skipped, Ordering::Relaxed);
        self.useful_claims.fetch_add(1, Ordering::Relaxed);
    }

    /// Mean files visited fruitlessly per unit of work claimed.
    #[expect(clippy::cast_precision_loss, reason = "counts stay far below 2^52")]
    pub(crate) fn wasted_visits_per_claim(&self) -> f64 {
        use std::sync::atomic::Ordering;
        let claims = self.useful_claims.load(Ordering::Relaxed);
        if claims == 0 {
            return 0.0;
        }
        self.wasted_visits.load(Ordering::Relaxed) as f64 / claims as f64
    }
}

/// A read-only view of the counters, for logging.
#[derive(Debug, Clone, Copy)]
pub struct MergePhaseBreakdown {
    /// Busy seconds spent reading raw spill blocks from disk, and block count.
    pub read: (f64, u64),
    /// Busy seconds spent decompressing spill blocks, and block count.
    pub decompress: (f64, u64),
    /// Busy seconds spent compressing output blocks, and block count.
    pub output_compress: (f64, u64),
    /// Busy seconds spent compressing Phase 1 spill blocks, and block count.
    /// Not part of the merge -- reported for context, excluded from merge totals.
    pub spill_compress: (f64, u64),
}

impl MergePhaseCounters {
    pub(crate) fn snapshot(&self) -> MergePhaseBreakdown {
        MergePhaseBreakdown {
            read: self.read.snapshot(),
            decompress: self.decompress.snapshot(),
            output_compress: self.output_compress.snapshot(),
            spill_compress: self.spill_compress.snapshot(),
        }
    }
}

impl MergePhaseBreakdown {
    /// Total worker busy seconds attributable to the merge.
    ///
    /// Excludes `spill_compress`, which happens in Phase 1. Deliberately *not*
    /// a partition of merge wall clock — see the module doc — but it is the
    /// numerator of a meaningful utilization figure: divided by
    /// `merge_wall * num_workers` it says what fraction of the available worker
    /// capacity the merge actually used.
    #[must_use]
    pub fn total_busy_secs(self) -> f64 {
        self.read.0 + self.decompress.0 + self.output_compress.0
    }

    /// Fraction of available worker capacity the merge consumed, in `[0, 1]`.
    ///
    /// `busy / (wall * workers)`. Well below 1 means workers idled -- the merge
    /// was bound by something other than worker CPU (the consumer thread, lock
    /// contention, or I/O latency), and adding compression threads would not
    /// help. Near 1 means the pool was saturated and worker CPU is the wall.
    #[must_use]
    #[allow(
        clippy::cast_precision_loss,
        reason = "thread counts are small integers; f64 represents them exactly"
    )]
    pub fn worker_utilization(self, merge_wall_secs: f64, num_workers: usize) -> Option<f64> {
        let capacity = merge_wall_secs * num_workers as f64;
        (capacity > 0.0).then(|| self.total_busy_secs() / capacity)
    }

    /// Whether anything was recorded, so a sort that never merged stays silent.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.read.1 == 0 && self.decompress.1 == 0 && self.output_compress.1 == 0
    }
}

/// What limited the merge, inferred from worker utilization and consumer wait.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MergeVerdict {
    /// The worker pool was saturated; worker CPU is the wall.
    CpuBound,
    /// Workers idled *and* the consumer waited -- neither side is the
    /// bottleneck, which points at storage.
    IoBound,
    /// Neither condition clearly holds.
    Mixed,
}

/// Worker utilization at or above which the pool counts as saturated.
const SATURATED_UTILIZATION: f64 = 0.85;
/// Worker utilization below which the pool counts as idle.
const IDLE_UTILIZATION: f64 = 0.60;
/// Consumer wait fraction above which the consumer counts as starved.
const STARVED_FETCH_FRACTION: f64 = 0.50;

/// Classify what limited the merge.
///
/// The discriminating case is *both* sides being unfavourable: an idle pool
/// alone could just mean a cheap merge, and a waiting consumer alone could mean
/// slow workers. Together they mean neither side is the constraint, which
/// leaves storage.
///
/// Thresholds are calibrated from a handful of measured runs, not a model, so
/// callers should present the result as a hint and show the numbers behind it.
#[must_use]
pub fn classify_merge(utilization: f64, fetch_fraction: f64) -> MergeVerdict {
    if utilization >= SATURATED_UTILIZATION {
        MergeVerdict::CpuBound
    } else if utilization < IDLE_UTILIZATION && fetch_fraction > STARVED_FETCH_FRACTION {
        MergeVerdict::IoBound
    } else {
        MergeVerdict::Mixed
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::sync::Arc;
    use std::thread;
    use std::time::{Duration, Instant};

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

    /// The thresholds are pinned against the runs that calibrated them, so a
    /// future tweak has to confront the evidence rather than just move a number.
    #[rstest]
    // Measured on a 780M-record WGS merge (EBS): workers half idle AND the
    // consumer waiting three quarters of its loop. Both "optimized" arms of that
    // experiment got slower, confirming storage was the constraint.
    #[case::measured_ebs_merge(0.52, 0.74, MergeVerdict::IoBound)]
    // Same code on a laptop with local NVMe: the pool saturates instead.
    #[case::measured_local_nvme(0.93, 0.03, MergeVerdict::CpuBound)]
    // Idle pool but a consumer that is NOT waiting: the merge is simply cheap,
    // which must not be reported as an I/O problem.
    #[case::idle_pool_busy_consumer(0.20, 0.05, MergeVerdict::Mixed)]
    // Waiting consumer but a busy pool: workers are the constraint, not storage.
    #[case::busy_pool_waiting_consumer(0.80, 0.90, MergeVerdict::Mixed)]
    #[case::at_saturation_boundary(0.85, 0.99, MergeVerdict::CpuBound)]
    #[case::just_below_saturation(0.84, 0.10, MergeVerdict::Mixed)]
    fn test_classify_merge(
        #[case] utilization: f64,
        #[case] fetch_fraction: f64,
        #[case] expected: MergeVerdict,
    ) {
        assert_eq!(classify_merge(utilization, fetch_fraction), expected);
    }

    /// A saturated pool is CPU-bound regardless of consumer wait, since the
    /// saturation test is checked first.
    #[test]
    fn test_saturation_takes_precedence_over_consumer_wait() {
        assert_eq!(classify_merge(0.95, 0.95), MergeVerdict::CpuBound);
    }
}
