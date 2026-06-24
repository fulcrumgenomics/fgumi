//! Per-step pipeline statistics: dispatch counts, outcome breakdown, and
//! cumulative `try_run` time. Indexed by `StepIdx` so it works for any
//! chain shape.
//!
//! Stats are opt-in. When `PipelineConfig::stats` is `Some(arc)`, the
//! worker loop times each `dispatch_one_step` call and records the
//! outcome. When `None`, the loop pays no extra cost.
//!
//! ## What's recorded
//!
//! Per step:
//!   - `try_run_total` — total `try_run_erased` dispatches that returned
//!     a result (excluding `Skip` entries).
//!   - `progress_count` — `StepOutcome::Progress`.
//!   - `no_progress_count` — `StepOutcome::NoProgress`.
//!   - `contention_count` — `StepOutcome::Contention` (Serial step mutex
//!     held by another worker; or skipped via `try_lock`). Always 0 under the
//!     fused single-thread driver, which holds no mutex and never contends.
//!   - `finished_count` — `StepOutcome::Finished` (any step on end-of-stream:
//!     source, mid, or sink all record `Finished` once their inputs drain).
//!   - `error_count` — `try_run_erased` returned `Err`.
//!   - `total_run_ns` — cumulative wall time across all dispatches.
//!
//! ## What's *not* recorded yet
//!
//! - Queue-depth samples. The legacy framework periodically samples each
//!   queue's pending count from a monitor thread; porting that needs a
//!   public `InputHandle::depth()` method and a sampler driver. Deferred
//!   to a follow-up commit.
//! - Per-thread step counts. Legacy carries per-`(thread, step)` counters
//!   for bottleneck attribution. Useful for the rebalancer (#17) port,
//!   not for first-cut observability. Deferred until the rebalancer
//!   work lands.

use std::fmt;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use crate::step::StepOutcome;
use crate::topology::StepIdx;

/// Atomic counters for a single step.
#[derive(Debug)]
pub struct StepStats {
    pub try_run_total: AtomicU64,
    pub progress_count: AtomicU64,
    pub no_progress_count: AtomicU64,
    pub contention_count: AtomicU64,
    pub finished_count: AtomicU64,
    pub error_count: AtomicU64,
    pub total_run_ns: AtomicU64,
    /// Wall ns from pipeline start to the start of this step's FIRST
    /// Progress dispatch. `u64::MAX` until set. Lets us see when each
    /// step actually first did useful work.
    pub first_progress_ns: AtomicU64,
    /// Wall ns from pipeline start to the END of this step's LAST
    /// Progress dispatch. Updated on every Progress (monotonic max).
    /// `0` if the step never made progress.
    pub last_progress_ns: AtomicU64,
}

impl Default for StepStats {
    fn default() -> Self {
        Self {
            try_run_total: AtomicU64::new(0),
            progress_count: AtomicU64::new(0),
            no_progress_count: AtomicU64::new(0),
            contention_count: AtomicU64::new(0),
            finished_count: AtomicU64::new(0),
            error_count: AtomicU64::new(0),
            total_run_ns: AtomicU64::new(0),
            first_progress_ns: AtomicU64::new(u64::MAX),
            last_progress_ns: AtomicU64::new(0),
        }
    }
}

impl StepStats {
    fn snapshot(&self) -> StepStatsSnapshot {
        StepStatsSnapshot {
            try_run_total: self.try_run_total.load(Ordering::Relaxed),
            progress_count: self.progress_count.load(Ordering::Relaxed),
            no_progress_count: self.no_progress_count.load(Ordering::Relaxed),
            contention_count: self.contention_count.load(Ordering::Relaxed),
            finished_count: self.finished_count.load(Ordering::Relaxed),
            error_count: self.error_count.load(Ordering::Relaxed),
            total_run_ns: self.total_run_ns.load(Ordering::Relaxed),
            first_progress_ns: self.first_progress_ns.load(Ordering::Relaxed),
            last_progress_ns: self.last_progress_ns.load(Ordering::Relaxed),
        }
    }
}

/// Per-step counter container. Sized to match the pipeline's chain length;
/// callers obtain one via `Pipeline::stats()`.
#[derive(Debug)]
pub struct PipelineStats {
    steps: Box<[StepStats]>,
    step_names: Box<[&'static str]>,
    /// Anchor for first/last-progress timestamps. Set at `PipelineStats`
    /// construction; all `first_progress_ns` / `last_progress_ns` values
    /// are wall-ns elapsed from this `Instant`.
    pipeline_start: Instant,
}

impl PipelineStats {
    /// Construct a stats container sized to `step_names.len()`. Each step's
    /// counters start at zero. Wrap in `Arc` to share across worker threads.
    #[must_use]
    pub fn new(step_names: Vec<&'static str>) -> Self {
        let steps = (0..step_names.len()).map(|_| StepStats::default()).collect::<Vec<_>>();
        Self {
            steps: steps.into_boxed_slice(),
            step_names: step_names.into_boxed_slice(),
            pipeline_start: Instant::now(),
        }
    }

    /// Wall ns elapsed since pipeline start. Used by the driver to stamp
    /// per-step first/last progress timestamps.
    #[must_use]
    pub fn elapsed_ns(&self) -> u64 {
        u64::try_from(self.pipeline_start.elapsed().as_nanos()).unwrap_or(u64::MAX)
    }

    #[must_use]
    pub fn n_steps(&self) -> usize {
        self.steps.len()
    }

    #[must_use]
    pub fn step_name(&self, step: StepIdx) -> &'static str {
        self.step_names[step.0]
    }

    /// Record a successful `try_run_erased` outcome for the given step.
    /// Hot path: relaxed atomics, no allocation, no locking. `start_ns`
    /// is wall-ns at dispatch start (relative to `pipeline_start`);
    /// `elapsed_ns` is the dispatch duration. On `Progress` we stamp
    /// the step's first/last active timestamps.
    #[inline]
    pub fn record(&self, step: StepIdx, outcome: StepOutcome, start_ns: u64, elapsed_ns: u64) {
        let s = &self.steps[step.0];
        s.try_run_total.fetch_add(1, Ordering::Relaxed);
        s.total_run_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
        match outcome {
            StepOutcome::Progress => {
                s.progress_count.fetch_add(1, Ordering::Relaxed);
                // CAS-min first_progress_ns (initially u64::MAX).
                let mut cur = s.first_progress_ns.load(Ordering::Relaxed);
                while start_ns < cur {
                    match s.first_progress_ns.compare_exchange_weak(
                        cur,
                        start_ns,
                        Ordering::Relaxed,
                        Ordering::Relaxed,
                    ) {
                        Ok(_) => break,
                        Err(seen) => cur = seen,
                    }
                }
                // last_progress_ns = max(last, start + elapsed).
                let end_ns = start_ns.saturating_add(elapsed_ns);
                let mut cur = s.last_progress_ns.load(Ordering::Relaxed);
                while end_ns > cur {
                    match s.last_progress_ns.compare_exchange_weak(
                        cur,
                        end_ns,
                        Ordering::Relaxed,
                        Ordering::Relaxed,
                    ) {
                        Ok(_) => break,
                        Err(seen) => cur = seen,
                    }
                }
            }
            StepOutcome::NoProgress => {
                s.no_progress_count.fetch_add(1, Ordering::Relaxed);
            }
            StepOutcome::Contention => {
                s.contention_count.fetch_add(1, Ordering::Relaxed);
            }
            StepOutcome::Finished => {
                s.finished_count.fetch_add(1, Ordering::Relaxed);
            }
        }
    }

    /// Record an error returned by `try_run_erased`. Counted toward
    /// `try_run_total` and `total_run_ns`; outcome buckets are not bumped.
    #[inline]
    pub fn record_error(&self, step: StepIdx, _start_ns: u64, elapsed_ns: u64) {
        let s = &self.steps[step.0];
        s.try_run_total.fetch_add(1, Ordering::Relaxed);
        s.error_count.fetch_add(1, Ordering::Relaxed);
        s.total_run_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
    }

    /// Snapshot all per-step counters into an owned, lock-free struct
    /// suitable for printing or further analysis.
    #[must_use]
    pub fn snapshot(&self) -> StatsSnapshot {
        let steps = self
            .steps
            .iter()
            .zip(self.step_names.iter())
            .map(|(stats, &name)| (name, stats.snapshot()))
            .collect();
        StatsSnapshot { steps }
    }
}

/// Plain (non-atomic) snapshot of `PipelineStats` at a moment in time.
#[derive(Debug, Clone)]
pub struct StatsSnapshot {
    pub steps: Vec<(&'static str, StepStatsSnapshot)>,
}

#[derive(Debug, Clone, Copy)]
pub struct StepStatsSnapshot {
    pub try_run_total: u64,
    pub progress_count: u64,
    pub no_progress_count: u64,
    pub contention_count: u64,
    pub finished_count: u64,
    pub error_count: u64,
    pub total_run_ns: u64,
    /// Wall ns (from pipeline start) of this step's first Progress
    /// dispatch start. `u64::MAX` if no Progress was ever recorded.
    pub first_progress_ns: u64,
    /// Wall ns (from pipeline start) of this step's last Progress
    /// dispatch end. `0` if no Progress was ever recorded.
    pub last_progress_ns: u64,
}

impl StepStatsSnapshot {
    #[must_use]
    pub fn avg_run_ns(&self) -> Option<u64> {
        if self.try_run_total == 0 { None } else { Some(self.total_run_ns / self.try_run_total) }
    }
}

impl fmt::Display for StatsSnapshot {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Pipeline stats ({} step{}):", self.steps.len(), pluralize(self.steps.len()))?;
        writeln!(
            f,
            "  {:<28} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>14} {:>12} {:>12}",
            "step",
            "tries",
            "progress",
            "noprog",
            "content",
            "fin",
            "err",
            "total_ms",
            "first_ms",
            "last_ms",
        )?;
        for (name, s) in &self.steps {
            // ns -> ms with three decimal places. The `as f64` cast can lose
            // precision above 2^52 ns (~52 days of cumulative run time per
            // step), well beyond any real pipeline run.
            #[allow(clippy::cast_precision_loss)]
            let total_ms = (s.total_run_ns as f64) / 1_000_000.0;
            let first_str = if s.first_progress_ns == u64::MAX {
                "-".to_string()
            } else {
                #[allow(clippy::cast_precision_loss)]
                let ms = (s.first_progress_ns as f64) / 1_000_000.0;
                format!("{ms:.1}")
            };
            let last_str = if s.last_progress_ns == 0 {
                "-".to_string()
            } else {
                #[allow(clippy::cast_precision_loss)]
                let ms = (s.last_progress_ns as f64) / 1_000_000.0;
                format!("{ms:.1}")
            };
            writeln!(
                f,
                "  {:<28} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>14.3} {:>12} {:>12}",
                name,
                s.try_run_total,
                s.progress_count,
                s.no_progress_count,
                s.contention_count,
                s.finished_count,
                s.error_count,
                total_ms,
                first_str,
                last_str,
            )?;
        }
        Ok(())
    }
}

fn pluralize(n: usize) -> &'static str {
    if n == 1 { "" } else { "s" }
}

/// Wrapper used by the driver to keep stats access ergonomic. Always pass
/// the `Option<&StatsRecorder>` shape so the hot path is a single
/// `if let Some` branch.
pub struct StatsRecorder<'a> {
    inner: &'a Arc<PipelineStats>,
}

impl<'a> StatsRecorder<'a> {
    #[must_use]
    pub fn new(stats: &'a Arc<PipelineStats>) -> Self {
        Self { inner: stats }
    }

    #[inline]
    pub fn record(&self, step: StepIdx, outcome: StepOutcome, start_ns: u64, elapsed_ns: u64) {
        self.inner.record(step, outcome, start_ns, elapsed_ns);
    }

    #[inline]
    pub fn record_error(&self, step: StepIdx, start_ns: u64, elapsed_ns: u64) {
        self.inner.record_error(step, start_ns, elapsed_ns);
    }

    /// Wall ns elapsed since pipeline start. Used by the driver to stamp
    /// per-step first/last progress timestamps.
    #[inline]
    #[must_use]
    pub fn elapsed_ns(&self) -> u64 {
        self.inner.elapsed_ns()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn record_increments_progress_bucket() {
        let stats = PipelineStats::new(vec!["A", "B"]);
        stats.record(StepIdx(0), StepOutcome::Progress, 0, 100);
        stats.record(StepIdx(0), StepOutcome::Progress, 200, 50);
        stats.record(StepIdx(1), StepOutcome::NoProgress, 0, 25);

        let snap = stats.snapshot();
        assert_eq!(snap.steps.len(), 2);

        let (name_a, a) = &snap.steps[0];
        assert_eq!(*name_a, "A");
        assert_eq!(a.try_run_total, 2);
        assert_eq!(a.progress_count, 2);
        assert_eq!(a.no_progress_count, 0);
        assert_eq!(a.total_run_ns, 150);
        assert_eq!(a.first_progress_ns, 0);
        assert_eq!(a.last_progress_ns, 250);

        let (name_b, b) = &snap.steps[1];
        assert_eq!(*name_b, "B");
        assert_eq!(b.try_run_total, 1);
        assert_eq!(b.no_progress_count, 1);
        assert_eq!(b.total_run_ns, 25);
        assert_eq!(b.first_progress_ns, u64::MAX);
        assert_eq!(b.last_progress_ns, 0);
    }

    #[test]
    fn record_buckets_each_outcome() {
        let stats = PipelineStats::new(vec!["S"]);
        stats.record(StepIdx(0), StepOutcome::Progress, 0, 1);
        stats.record(StepIdx(0), StepOutcome::NoProgress, 10, 2);
        stats.record(StepIdx(0), StepOutcome::Contention, 20, 3);
        stats.record(StepIdx(0), StepOutcome::Finished, 30, 4);
        stats.record_error(StepIdx(0), 40, 5);

        let snap = stats.snapshot();
        let (_name, s) = &snap.steps[0];
        assert_eq!(s.try_run_total, 5);
        assert_eq!(s.progress_count, 1);
        assert_eq!(s.no_progress_count, 1);
        assert_eq!(s.contention_count, 1);
        assert_eq!(s.finished_count, 1);
        assert_eq!(s.error_count, 1);
        assert_eq!(s.total_run_ns, 1 + 2 + 3 + 4 + 5);
    }

    #[test]
    fn avg_run_ns_handles_empty_step() {
        let snap = StepStatsSnapshot {
            try_run_total: 0,
            progress_count: 0,
            no_progress_count: 0,
            contention_count: 0,
            finished_count: 0,
            error_count: 0,
            total_run_ns: 0,
            first_progress_ns: u64::MAX,
            last_progress_ns: 0,
        };
        assert_eq!(snap.avg_run_ns(), None);

        let snap2 = StepStatsSnapshot { try_run_total: 4, total_run_ns: 1000, ..snap };
        assert_eq!(snap2.avg_run_ns(), Some(250));
    }

    #[test]
    fn display_renders_header_and_rows() {
        let stats = PipelineStats::new(vec!["StepA", "StepB"]);
        stats.record(StepIdx(0), StepOutcome::Progress, 0, 10);
        let s = format!("{}", stats.snapshot());
        assert!(s.contains("Pipeline stats (2 steps):"));
        assert!(s.contains("StepA"));
        assert!(s.contains("StepB"));
        assert!(s.contains("step"));
    }
}
