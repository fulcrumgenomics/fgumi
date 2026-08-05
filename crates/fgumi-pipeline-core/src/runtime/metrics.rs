//! Per-edge instrumentation counters + occupancy histogram.
//!
//! An [`EdgeMetrics`] is the edge-side complement of the step-side
//! [`PipelineStats`](super::stats::PipelineStats): one per instrumented queue
//! edge, shared (`Arc`) between the producer's transport (push/reject counts)
//! and the consumer's input handle (pop/empty counts), and sampled periodically
//! for occupancy. All counters are `Relaxed` atomics — these are statistics, not
//! synchronization (staleness is fine), mirroring `ByteBoundedQueue::current_bytes`.
//!
//! What the histogram alone can classify is [`RawOccupancy`]; the richer
//! `Empty`-vs-`Starved` / `Full`-vs-`Backpressured` split needs the reject/empty
//! *rates* and is done by the renderer, not here.

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

/// Number of occupancy histogram buckets (depth fraction `0.0..=1.0` split into
/// `OCCUPANCY_BUCKETS` equal bands; the top band captures exactly-full).
pub const OCCUPANCY_BUCKETS: usize = 8;

/// Per-edge counters + occupancy histogram. Construct with [`EdgeMetrics::new`]
/// (returns an `Arc` so the producer transport and consumer input handle share
/// one instance). All methods are lock-free `Relaxed` atomic updates.
#[derive(Debug)]
pub struct EdgeMetrics {
    /// Items the producer successfully pushed.
    pushed_items: AtomicU64,
    /// Bytes pushed (only meaningful for byte-bounded edges; `0` otherwise).
    pushed_bytes: AtomicU64,
    /// Items the consumer popped.
    popped_items: AtomicU64,
    /// Bytes popped (byte-bounded edges only).
    popped_bytes: AtomicU64,
    /// `try_push` rejections — backpressure events (producer wanted to push, the
    /// edge was full).
    push_rejections: AtomicU64,
    /// Consumer `pop` on an empty edge — starvation events.
    pop_empties: AtomicU64,
    /// Number of occupancy samples taken (sampler ticks).
    depth_samples: AtomicU64,
    /// Histogram of occupancy fraction at sample time.
    occupancy_buckets: [AtomicU64; OCCUPANCY_BUCKETS],
    /// Σ of `depth_fraction * 1000` over all samples, for the mean.
    occupancy_sum_milli: AtomicU64,
    /// Σ of raw occupied **bytes** over all samples. Kept alongside the fraction
    /// sum so the mean occupancy can be reported in absolute bytes — which,
    /// unlike `fraction × final_limit`, stays correct when the byte limit is
    /// rebalanced at runtime (`queue_memory_total`). For an ordered edge the
    /// sampled bytes include the `ReorderStage` overflow stash.
    occupancy_bytes_sum: AtomicU64,
}

impl EdgeMetrics {
    /// Construct a fresh metrics instance, shared via `Arc`.
    #[must_use]
    pub fn new() -> Arc<Self> {
        Arc::new(Self {
            pushed_items: AtomicU64::new(0),
            pushed_bytes: AtomicU64::new(0),
            popped_items: AtomicU64::new(0),
            popped_bytes: AtomicU64::new(0),
            push_rejections: AtomicU64::new(0),
            pop_empties: AtomicU64::new(0),
            depth_samples: AtomicU64::new(0),
            occupancy_buckets: std::array::from_fn(|_| AtomicU64::new(0)),
            occupancy_sum_milli: AtomicU64::new(0),
            occupancy_bytes_sum: AtomicU64::new(0),
        })
    }

    /// Record a successful producer push of `bytes` (pass `0` for count-bounded
    /// / unbounded edges that don't track bytes).
    pub fn record_push(&self, bytes: u64) {
        self.pushed_items.fetch_add(1, Ordering::Relaxed);
        self.pushed_bytes.fetch_add(bytes, Ordering::Relaxed);
    }

    /// Record a successful consumer pop of `bytes`.
    pub fn record_pop(&self, bytes: u64) {
        self.popped_items.fetch_add(1, Ordering::Relaxed);
        self.popped_bytes.fetch_add(bytes, Ordering::Relaxed);
    }

    /// Record a `try_push` rejection (backpressure).
    pub fn record_reject(&self) {
        self.push_rejections.fetch_add(1, Ordering::Relaxed);
    }

    /// Record a consumer pop that found the edge empty (starvation).
    pub fn record_empty(&self) {
        self.pop_empties.fetch_add(1, Ordering::Relaxed);
    }

    /// Record one occupancy sample from the edge's current `occupied_bytes` and
    /// its byte `limit_bytes`. The histogram bucket and fraction use
    /// `occupied_bytes / limit_bytes` clamped to `0.0..=1.0`; the **raw**
    /// `occupied_bytes` (un-clamped — so it stays truthful when an ordered
    /// edge's reorder stash pushes total buffered bytes past the transport
    /// limit) is summed for the byte-accurate mean the latency estimate uses. A
    /// `limit_bytes` of `0` is treated as fraction `0.0` defensively (the
    /// sampler already skips count/unbounded edges).
    // Casts are bounded statistics: `f ∈ [0,1]` so `f * BUCKETS` ∈ [0,8] and
    // `f * 1000` ∈ [0,1000] — no meaningful truncation/sign loss/precision loss.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss, clippy::cast_precision_loss)]
    pub fn record_depth(&self, occupied_bytes: u64, limit_bytes: u64) {
        let f = if limit_bytes == 0 {
            0.0
        } else {
            (occupied_bytes as f32 / limit_bytes as f32).clamp(0.0, 1.0)
        };
        let bucket = ((f * OCCUPANCY_BUCKETS as f32) as usize).min(OCCUPANCY_BUCKETS - 1);
        self.occupancy_buckets[bucket].fetch_add(1, Ordering::Relaxed);
        self.occupancy_sum_milli.fetch_add((f * 1000.0) as u64, Ordering::Relaxed);
        self.occupancy_bytes_sum.fetch_add(occupied_bytes, Ordering::Relaxed);
        self.depth_samples.fetch_add(1, Ordering::Relaxed);
    }

    /// Take a consistent-enough point-in-time snapshot of the counters.
    // Mean is a display statistic; u64→f32 precision loss past 2^23 samples is
    // irrelevant to a 0..1 occupancy mean.
    #[allow(clippy::cast_precision_loss)]
    #[must_use]
    pub fn snapshot(&self) -> EdgeMetricsSnapshot {
        let buckets: [u64; OCCUPANCY_BUCKETS] =
            std::array::from_fn(|i| self.occupancy_buckets[i].load(Ordering::Relaxed));
        let depth_samples = self.depth_samples.load(Ordering::Relaxed);
        let mean_occupancy = if depth_samples == 0 {
            0.0
        } else {
            (self.occupancy_sum_milli.load(Ordering::Relaxed) as f32 / 1000.0)
                / depth_samples as f32
        };
        let mean_occupancy_bytes = if depth_samples == 0 {
            0.0
        } else {
            self.occupancy_bytes_sum.load(Ordering::Relaxed) as f64 / depth_samples as f64
        };
        EdgeMetricsSnapshot {
            pushed_items: self.pushed_items.load(Ordering::Relaxed),
            pushed_bytes: self.pushed_bytes.load(Ordering::Relaxed),
            popped_items: self.popped_items.load(Ordering::Relaxed),
            popped_bytes: self.popped_bytes.load(Ordering::Relaxed),
            push_rejections: self.push_rejections.load(Ordering::Relaxed),
            pop_empties: self.pop_empties.load(Ordering::Relaxed),
            depth_samples,
            raw_occupancy: RawOccupancy::from_buckets(&buckets, depth_samples),
            mean_occupancy,
            mean_occupancy_bytes,
        }
    }
}

/// A point-in-time read of an [`EdgeMetrics`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeMetricsSnapshot {
    /// Items successfully pushed onto the edge by the producer.
    pub pushed_items: u64,
    /// Bytes successfully pushed onto the edge (0 for count/unbounded edges).
    pub pushed_bytes: u64,
    /// Items successfully popped off the edge by the consumer.
    pub popped_items: u64,
    /// Bytes successfully popped off the edge (0 for count/unbounded edges).
    pub popped_bytes: u64,
    /// Push attempts rejected by backpressure (edge at its limit).
    pub push_rejections: u64,
    /// Pop attempts that found the edge empty (starvation signal).
    pub pop_empties: u64,
    /// Number of occupancy samples the background sampler took on this edge.
    pub depth_samples: u64,
    /// What the occupancy histogram alone says (no rate-based refinement).
    pub raw_occupancy: RawOccupancy,
    /// Mean occupancy fraction `0.0..=1.0` over all samples (`0.0` if none).
    pub mean_occupancy: f32,
    /// Mean occupied **bytes** over all samples (`0.0` if none). Unlike
    /// `mean_occupancy` (a fraction that must be multiplied by a limit to
    /// recover bytes), this is measured directly at sample time, so the
    /// Little's-Law latency estimate stays correct even when `queue_memory_total`
    /// rebalances the byte limit mid-run. For an ordered edge it includes the
    /// `ReorderStage` overflow stash.
    pub mean_occupancy_bytes: f64,
}

/// Occupancy classification derivable from the histogram **alone** (no
/// reject/empty rates). The renderer refines `MostlyEmpty`→`Empty`/`Starved`
/// and `MostlyFull`→`Full`/`Backpressured` using the rates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RawOccupancy {
    /// Insufficient samples to classify.
    Unknown,
    /// ≥80% of samples in the bottom bucket.
    MostlyEmpty,
    /// ≥80% of samples in the top bucket.
    MostlyFull,
    /// Bottom and top buckets each hold >25% of samples (bursty coupling).
    Bimodal,
    /// Spread across the middle — neither side bound.
    Healthy,
}

impl RawOccupancy {
    /// Classify from the bucket histogram and total sample count.
    // Ratios are display statistics; u64→f64 precision loss is irrelevant to the
    // 0.25/0.80 thresholds.
    #[allow(clippy::cast_precision_loss)]
    #[must_use]
    pub fn from_buckets(buckets: &[u64; OCCUPANCY_BUCKETS], samples: u64) -> Self {
        if samples == 0 {
            return Self::Unknown;
        }
        let s = samples as f64;
        let bottom = buckets[0] as f64 / s;
        let top = buckets[OCCUPANCY_BUCKETS - 1] as f64 / s;
        if bottom >= 0.80 {
            Self::MostlyEmpty
        } else if top >= 0.80 {
            Self::MostlyFull
        } else if bottom > 0.25 && top > 0.25 {
            Self::Bimodal
        } else {
            Self::Healthy
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[test]
    fn edge_metrics_counts_are_exact() {
        let m = EdgeMetrics::new();
        for _ in 0..3 {
            m.record_push(100);
        }
        m.record_reject();
        for _ in 0..2 {
            m.record_pop(100);
        }
        m.record_empty();
        let s = m.snapshot();
        assert_eq!(s.pushed_items, 3);
        assert_eq!(s.pushed_bytes, 300);
        assert_eq!(s.push_rejections, 1);
        assert_eq!(s.popped_items, 2);
        assert_eq!(s.popped_bytes, 200);
        assert_eq!(s.pop_empties, 1);
    }

    // Each case records a depth pattern, then asserts the classification (and,
    // where the pattern pins one, the mean occupancy within a tolerance). The
    // `mostly_empty` case reuses the `(mean, tol)` form as `(0.0, 0.02)`, i.e.
    // `mean < 0.02`, since occupancy is non-negative. `fill` is a non-capturing
    // closure coerced to a `fn` pointer so it can ride as an rstest case value.
    #[rstest]
    // Samples are `(occupied_bytes, limit_bytes)`; a 1000-byte limit makes the
    // occupied byte count read directly as the depth fraction ×1000.
    #[case::mostly_full(|m: &EdgeMetrics| for _ in 0..100 { m.record_depth(1000, 1000); }, RawOccupancy::MostlyFull, Some((1.0, 0.02)))]
    #[case::mostly_empty(|m: &EdgeMetrics| for _ in 0..100 { m.record_depth(0, 1000); }, RawOccupancy::MostlyEmpty, Some((0.0, 0.02)))]
    #[case::bimodal(|m: &EdgeMetrics| for _ in 0..50 { m.record_depth(0, 1000); m.record_depth(1000, 1000); }, RawOccupancy::Bimodal, None)]
    // Spread uniformly across the middle buckets → mean ≈ 0.45.
    #[case::healthy(|m: &EdgeMetrics| for i in 0..100u32 { m.record_depth(300 + u64::from(i % 4) * 100, 1000); }, RawOccupancy::Healthy, Some((0.45, 0.05)))]
    #[case::unknown(|_m: &EdgeMetrics| {}, RawOccupancy::Unknown, None)]
    fn occupancy_classification(
        #[case] fill: fn(&EdgeMetrics),
        #[case] expected: RawOccupancy,
        #[case] mean_within: Option<(f32, f32)>,
    ) {
        let m = EdgeMetrics::new();
        fill(&m);
        let s = m.snapshot();
        assert_eq!(s.raw_occupancy, expected);
        if let Some((mean, tol)) = mean_within {
            assert!(
                (s.mean_occupancy - mean).abs() < tol,
                "mean_occupancy {} not within {tol} of {mean}",
                s.mean_occupancy
            );
        }
    }

    #[test]
    fn mean_occupancy_bytes_tracks_absolute_bytes_independent_of_limit() {
        // Sampled occupied bytes are averaged in absolute terms, not as a
        // fraction of the limit — so a later limit change can't distort the mean
        // (the property the derived-latency estimate relies on). Two 600-byte
        // samples average to 600 bytes regardless of the 1000-byte limit used.
        let m = EdgeMetrics::new();
        m.record_depth(600, 1000);
        m.record_depth(600, 1000);
        let s = m.snapshot();
        assert!((s.mean_occupancy_bytes - 600.0).abs() < 1e-9, "{}", s.mean_occupancy_bytes);
        assert!((s.mean_occupancy - 0.6).abs() < 0.01, "{}", s.mean_occupancy);
    }

    #[test]
    fn mean_occupancy_bytes_keeps_raw_bytes_when_stash_exceeds_limit() {
        // An ordered edge's reorder stash can push total buffered bytes past the
        // transport limit; the fraction saturates at 1.0 but the byte mean stays
        // truthful (1500 bytes), so latency is not under-counted.
        let m = EdgeMetrics::new();
        m.record_depth(1500, 1000); // 150% of limit
        let s = m.snapshot();
        assert!((s.mean_occupancy - 1.0).abs() < 1e-6, "fraction clamps to 1.0");
        assert!((s.mean_occupancy_bytes - 1500.0).abs() < 1e-9, "raw bytes preserved");
    }
}
