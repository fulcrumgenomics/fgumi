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
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use crate::step::StepOutcome;
use crate::topology::StepIdx;

/// Atomic counters for a single step.
#[derive(Debug)]
pub struct StepStats {
    /// Total `try_run` dispatches (every outcome, including errors).
    pub try_run_total: AtomicU64,
    /// Dispatches that returned `StepOutcome::Progress`.
    pub progress_count: AtomicU64,
    /// Dispatches that returned `StepOutcome::NoProgress`.
    pub no_progress_count: AtomicU64,
    /// Dispatches that returned `StepOutcome::Contention`.
    pub contention_count: AtomicU64,
    /// Dispatches that returned `StepOutcome::Finished`.
    pub finished_count: AtomicU64,
    /// Dispatches that returned an error from `try_run_erased`.
    pub error_count: AtomicU64,
    /// Cumulative wall-ns spent inside `try_run` across all dispatches.
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

/// Upper bound on per-worker utilization slots. Worker `thread_id`s index the
/// `worker_busy_ns` / `worker_idle_ns` arrays; ids at or above this are not
/// tracked (a no-op, never a panic). Sized well above any realistic
/// `--threads`, so a fixed array avoids threading `num_workers` through every
/// `PipelineStats::new` call site.
const MAX_TRACKED_WORKERS: usize = 512;

/// Per-step counter container. Sized to match the pipeline's chain length;
/// callers obtain one via `Pipeline::stats()`.
#[derive(Debug)]
pub struct PipelineStats {
    steps: Box<[StepStats]>,
    step_names: Box<[&'static str]>,
    /// Per-worker wall-ns spent dispatching (the sticky + round-robin work
    /// section of `run_worker_loop`), indexed by `WorkerCore::thread_id`.
    worker_busy_ns: Box<[AtomicU64]>,
    /// Per-worker wall-ns spent in the no-progress backoff sleep (idle/blocked
    /// waiting for upstream work), indexed by `thread_id`. Together with
    /// `worker_busy_ns` this answers "are workers utilised, or blocked".
    worker_idle_ns: Box<[AtomicU64]>,
    /// Per-step wall-ns a `StepKind::Detached` step's dedicated thread spent
    /// inside `try_run` (busy) and parked on its backoff (idle), indexed by
    /// `step_idx`. Detached threads have no `WorkerCore::thread_id`, so they are
    /// tracked here separately and are intentionally EXCLUDED from the pool
    /// utilisation line (which only sums `worker_busy_ns` / `worker_idle_ns`).
    /// Reported on their own line so the legacy "N + 2" split is visible without
    /// diluting the N-worker pool%.
    detached_busy_ns: Box<[AtomicU64]>,
    detached_idle_ns: Box<[AtomicU64]>,
    /// Per-step count of backoff-park events on a `StepKind::Detached` step's
    /// dedicated thread (one per [`backoff_park`](crate::runtime::detached) call),
    /// indexed by `step_idx`. With `detached_idle_ns` this gives the average park
    /// duration and the park-to-progress ratio — the signal for whether the
    /// backoff (vs a precise per-slot condvar) adds latency on the merge's path.
    detached_park_events: Box<[AtomicU64]>,
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
        let worker_busy_ns =
            (0..MAX_TRACKED_WORKERS).map(|_| AtomicU64::new(0)).collect::<Vec<_>>();
        let worker_idle_ns =
            (0..MAX_TRACKED_WORKERS).map(|_| AtomicU64::new(0)).collect::<Vec<_>>();
        let n_steps = steps.len();
        let detached_busy_ns = (0..n_steps).map(|_| AtomicU64::new(0)).collect::<Vec<_>>();
        let detached_idle_ns = (0..n_steps).map(|_| AtomicU64::new(0)).collect::<Vec<_>>();
        let detached_park_events = (0..n_steps).map(|_| AtomicU64::new(0)).collect::<Vec<_>>();
        Self {
            steps: steps.into_boxed_slice(),
            step_names: step_names.into_boxed_slice(),
            worker_busy_ns: worker_busy_ns.into_boxed_slice(),
            worker_idle_ns: worker_idle_ns.into_boxed_slice(),
            detached_busy_ns: detached_busy_ns.into_boxed_slice(),
            detached_idle_ns: detached_idle_ns.into_boxed_slice(),
            detached_park_events: detached_park_events.into_boxed_slice(),
            pipeline_start: Instant::now(),
        }
    }

    /// Accumulate `ns` of `try_run` (busy) time for a `StepKind::Detached`
    /// step's dedicated thread, keyed by `step_idx`. Excluded from pool
    /// utilisation; reported on the Detached line. Out-of-range steps are a
    /// silent no-op.
    #[inline]
    pub fn record_detached_busy(&self, step: StepIdx, ns: u64) {
        if let Some(c) = self.detached_busy_ns.get(step.0) {
            c.fetch_add(ns, Ordering::Relaxed);
        }
    }

    /// Accumulate `ns` of backoff-park (idle) time for a `StepKind::Detached`
    /// step's dedicated thread, keyed by `step_idx`. Excluded from pool
    /// utilisation; reported on the Detached line.
    #[inline]
    pub fn record_detached_idle(&self, step: StepIdx, ns: u64) {
        if let Some(c) = self.detached_idle_ns.get(step.0) {
            c.fetch_add(ns, Ordering::Relaxed);
        }
    }

    /// Increment the backoff-park event count for a `StepKind::Detached` step's
    /// dedicated thread (called once per park). Out-of-range steps are a silent
    /// no-op.
    #[inline]
    pub fn record_detached_park(&self, step: StepIdx) {
        if let Some(c) = self.detached_park_events.get(step.0) {
            c.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// Accumulate `ns` of dispatch (busy) time for `worker` (its `thread_id`).
    /// Ids `>= MAX_TRACKED_WORKERS` are silently dropped.
    #[inline]
    pub fn record_worker_busy(&self, worker: usize, ns: u64) {
        if let Some(c) = self.worker_busy_ns.get(worker) {
            c.fetch_add(ns, Ordering::Relaxed);
        }
    }

    /// Accumulate `ns` of backoff-sleep (idle) time for `worker`.
    #[inline]
    pub fn record_worker_idle(&self, worker: usize, ns: u64) {
        if let Some(c) = self.worker_idle_ns.get(worker) {
            c.fetch_add(ns, Ordering::Relaxed);
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
        // Only workers that recorded any activity (busy or idle) — the array is
        // sized to MAX_TRACKED_WORKERS but typically few slots are live.
        let workers = (0..self.worker_busy_ns.len())
            .map(|w| {
                (
                    w,
                    self.worker_busy_ns[w].load(Ordering::Relaxed),
                    self.worker_idle_ns[w].load(Ordering::Relaxed),
                )
            })
            .filter(|(_, busy, idle)| *busy != 0 || *idle != 0)
            .collect();
        StatsSnapshot { steps, workers, detached: self.detached_snapshot(), edges: Vec::new() }
    }

    /// Collect `(step_name, busy_ns, idle_ns)` for every step whose Detached
    /// thread recorded any activity. Shared by both snapshot builders so the
    /// Detached line renders identically with or without `--pipeline-trace`.
    fn detached_snapshot(&self) -> Vec<(&'static str, u64, u64, u64)> {
        (0..self.detached_busy_ns.len())
            .map(|s| {
                (
                    self.step_names[s],
                    self.detached_busy_ns[s].load(Ordering::Relaxed),
                    self.detached_idle_ns[s].load(Ordering::Relaxed),
                    self.detached_park_events[s].load(Ordering::Relaxed),
                )
            })
            .filter(|(_, busy, idle, _)| *busy != 0 || *idle != 0)
            .collect()
    }

    /// Like [`snapshot`](Self::snapshot) but also derives per-edge throughput /
    /// occupancy / latency from the chain's instrumented `edges` over a
    /// `wall_ns` run. Used by the `--pipeline-trace` end-of-run report.
    #[must_use]
    pub fn snapshot_with_edges(
        &self,
        edges: &[crate::runtime::contexts::RegisteredEdge],
        wall_ns: u64,
    ) -> StatsSnapshot {
        let mut snap = self.snapshot();
        snap.edges = edges
            .iter()
            .map(|e| {
                let ms = e.metrics.snapshot();
                let limit_bytes = e.depth_source.as_ref().map(|s| s.limit_bytes());
                compute_edge_stats(
                    e.producer_name,
                    e.consumer_name,
                    e.producer_step.0,
                    e.consumer_step.map(|s| s.0),
                    &ms,
                    limit_bytes,
                    wall_ns,
                )
            })
            .collect();
        snap
    }
}

/// Plain (non-atomic) snapshot of `PipelineStats` at a moment in time.
#[derive(Debug, Clone)]
pub struct StatsSnapshot {
    pub steps: Vec<(&'static str, StepStatsSnapshot)>,
    /// `(thread_id, busy_ns, idle_ns)` for each worker that did anything.
    pub workers: Vec<(usize, u64, u64)>,
    /// `(step_name, busy_ns, idle_ns, park_events)` for each `StepKind::Detached`
    /// step's dedicated thread that did anything. Tracked separately from
    /// `workers` and EXCLUDED from the pool utilisation line (legacy "N + 2" —
    /// the merge / writer threads are not pool workers); rendered on their own
    /// line. `park_events` is the backoff-park count (latency signal).
    pub detached: Vec<(&'static str, u64, u64, u64)>,
    /// Per-edge throughput / occupancy / latency. Empty unless the snapshot was
    /// built via [`PipelineStats::snapshot_with_edges`] (i.e. instrumentation on).
    pub edges: Vec<EdgeStatsSnapshot>,
}

/// Occupancy classification refined from the histogram's [`RawOccupancy`] plus
/// the reject/empty rates: `MostlyEmpty` + high empty-rate → `Starved`,
/// `MostlyFull` + high reject-rate → `Backpressured`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OccupancyClass {
    /// No occupancy samples (count/unbounded edge, or run too short).
    Unknown,
    /// Mostly empty, low empty-rate — idle, not starved.
    Empty,
    /// Mostly empty AND consumer frequently found it empty — producer-starved.
    Starved,
    /// Spread across the middle — neither side bound.
    Healthy,
    /// Mostly full AND producer frequently rejected — consumer-backpressured.
    Backpressured,
    /// Mostly full, low reject-rate — full but not actively rejecting.
    Full,
    /// Oscillates empty↔full (bursty / batch coupling).
    Bimodal,
}

/// Rate above which a `MostlyEmpty`/`MostlyFull` edge is reclassified as
/// `Starved`/`Backpressured`. 10% of attempts hitting the wall is a clear signal.
const RATE_REFINE_THRESHOLD: f64 = 0.10;

/// Per-edge derived statistics for the end-of-run report.
#[derive(Debug, Clone)]
pub struct EdgeStatsSnapshot {
    pub producer: &'static str,
    pub consumer: Option<&'static str>,
    /// Producer step index. The bottleneck verdict attributes edges to steps by
    /// this identity (not by `producer`/`consumer` name), so fan-out, fan-in, or
    /// duplicate step names are not misattributed to whichever edge matches a
    /// name first.
    pub producer_step: usize,
    /// Consumer step index, `None` for a terminal edge with no consumer.
    pub consumer_step: Option<usize>,
    pub items_per_s: f64,
    /// Consumer throughput in MiB/s (popped bytes / wall time / 2^20).
    pub mibytes_per_s: f64,
    pub class: OccupancyClass,
    pub mean_occupancy: f32,
    /// Fraction of push attempts rejected (backpressure signal).
    pub reject_rate: f64,
    /// Fraction of pop attempts that found the edge empty (starvation signal).
    pub empty_rate: f64,
    /// Derived residence time (Little's Law `L/X`), byte-bounded edges only;
    /// `None` for count/unbounded edges (no occupancy) or zero throughput.
    pub derived_latency_ms: Option<f64>,
}

/// Refine the histogram's raw occupancy class using the reject/empty rates.
#[must_use]
fn refine_occupancy(
    raw: crate::runtime::metrics::RawOccupancy,
    reject_rate: f64,
    empty_rate: f64,
) -> OccupancyClass {
    use crate::runtime::metrics::RawOccupancy;
    match raw {
        RawOccupancy::Unknown => OccupancyClass::Unknown,
        RawOccupancy::Healthy => OccupancyClass::Healthy,
        RawOccupancy::Bimodal => OccupancyClass::Bimodal,
        RawOccupancy::MostlyEmpty if empty_rate >= RATE_REFINE_THRESHOLD => OccupancyClass::Starved,
        RawOccupancy::MostlyEmpty => OccupancyClass::Empty,
        RawOccupancy::MostlyFull if reject_rate >= RATE_REFINE_THRESHOLD => {
            OccupancyClass::Backpressured
        }
        RawOccupancy::MostlyFull => OccupancyClass::Full,
    }
}

/// Little's-Law residence time (`L/X`), byte-bounded edges only. Mean occupancy
/// in ITEMS = `mean_occupancy_bytes / mean_item_bytes`.
///
/// `mean_occupancy_bytes` is sampled directly (absolute bytes at each tick), NOT
/// reconstructed as `fraction × limit`. That matters because `queue_memory_total`
/// can rebalance an edge's byte limit at runtime: the limit at teardown may
/// differ from the limits in force during sampling, so `fraction × final_limit`
/// would yield a materially wrong byte figure — and hence a wrong latency.
/// `limit_bytes` is retained only as the byte-bounded gate (`None` for
/// count/unbounded edges, which have no occupancy and thus no derived latency).
#[must_use]
#[allow(clippy::cast_precision_loss)]
fn derived_latency_ms(
    ms: &crate::runtime::metrics::EdgeMetricsSnapshot,
    limit_bytes: Option<u64>,
    items_per_s: f64,
) -> Option<f64> {
    if limit_bytes.is_none() || ms.pushed_items == 0 || items_per_s <= 0.0 {
        return None;
    }
    let mean_item_bytes = ms.pushed_bytes as f64 / ms.pushed_items as f64;
    if mean_item_bytes <= 0.0 {
        return None;
    }
    let mean_items = ms.mean_occupancy_bytes / mean_item_bytes;
    Some(1000.0 * mean_items / items_per_s)
}

/// Compute one edge's derived stats from its metrics snapshot + (for a
/// byte-bounded edge) its byte budget, over a `wall_ns` run. Pure — unit-tested
/// directly. `limit_bytes` is `None` for count/unbounded edges (no occupancy →
/// no derived latency).
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub(crate) fn compute_edge_stats(
    producer: &'static str,
    consumer: Option<&'static str>,
    producer_step: usize,
    consumer_step: Option<usize>,
    ms: &crate::runtime::metrics::EdgeMetricsSnapshot,
    limit_bytes: Option<u64>,
    wall_ns: u64,
) -> EdgeStatsSnapshot {
    let wall_secs = (wall_ns as f64) / 1e9;
    let items_per_s = if wall_secs > 0.0 { ms.popped_items as f64 / wall_secs } else { 0.0 };
    // Divisor is 1 MiB (2^20), so the field/column is MiB/s, not MB/s.
    let mibytes_per_s =
        if wall_secs > 0.0 { (ms.popped_bytes as f64 / wall_secs) / 1_048_576.0 } else { 0.0 };
    let push_attempts = ms.pushed_items + ms.push_rejections;
    let reject_rate =
        if push_attempts > 0 { ms.push_rejections as f64 / push_attempts as f64 } else { 0.0 };
    let pop_attempts = ms.popped_items + ms.pop_empties;
    let empty_rate =
        if pop_attempts > 0 { ms.pop_empties as f64 / pop_attempts as f64 } else { 0.0 };

    EdgeStatsSnapshot {
        producer,
        consumer,
        producer_step,
        consumer_step,
        items_per_s,
        mibytes_per_s,
        class: refine_occupancy(ms.raw_occupancy, reject_rate, empty_rate),
        mean_occupancy: ms.mean_occupancy,
        reject_rate,
        empty_rate,
        derived_latency_ms: derived_latency_ms(ms, limit_bytes, items_per_s),
    }
}

#[derive(Debug, Clone, Copy)]
pub struct StepStatsSnapshot {
    /// Total `try_run` dispatches (every outcome, including errors).
    pub try_run_total: u64,
    /// Dispatches that returned `StepOutcome::Progress`.
    pub progress_count: u64,
    /// Dispatches that returned `StepOutcome::NoProgress`.
    pub no_progress_count: u64,
    /// Dispatches that returned `StepOutcome::Contention`.
    pub contention_count: u64,
    /// Dispatches that returned `StepOutcome::Finished`.
    pub finished_count: u64,
    /// Dispatches that returned an error from `try_run_erased`.
    pub error_count: u64,
    /// Cumulative wall-ns spent inside `try_run` across all dispatches.
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

/// Severity of a [`Finding`] from the bottleneck verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Severity {
    /// The single rate-limiting step (full input edge + empty output edge).
    Primary,
    /// A contributing issue (spin, starvation) that is not the prime mover.
    Secondary,
    /// Whole-chain observation (e.g. latency-bound, not throughput-bound).
    Info,
}

/// One mechanically-derived diagnosis from [`bottleneck_verdict`].
#[derive(Debug, Clone)]
pub(crate) struct Finding {
    /// Triage classification. Set by [`bottleneck_verdict`] and asserted on by
    /// the crate's tests; the human-readable report renders only `message`, so
    /// the lib build never reads this field.
    #[allow(dead_code)]
    pub severity: Severity,
    pub message: String,
}

/// Rate above which a Serial step's `contention/tries` ratio is flagged as spin.
const SPIN_THRESHOLD: f64 = 0.20;

/// Mechanically locate the chain's bottleneck and contributing issues from a
/// snapshot's step + edge stats. Rules (all from already-collected numbers):
/// - **Primary**: a step whose input edge is `Full`/`Backpressured` AND output
///   edge is `Empty`/`Starved` — work piles into it, it can't fill downstream.
///   Tagged CPU-bound if it dominates `total_run_ns`, else coordination.
/// - **Secondary (spin)**: a step with a high `contention/tries` ratio.
/// - **Secondary (starvation)**: an edge whose consumer frequently finds it
///   empty — the producer (upstream) can't keep up.
/// - **Info**: if no primary and no full/empty edge, the chain is
///   latency/coordination-bound, not throughput-bound.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub(crate) fn bottleneck_verdict(snap: &StatsSnapshot) -> Vec<Finding> {
    use OccupancyClass::{Backpressured, Empty, Full, Starved};
    let mut findings = Vec::new();
    let total_cpu_ns: u64 = snap.steps.iter().map(|(_, s)| s.total_run_ns).sum();
    // Attribute edges to steps by step IDENTITY (index), not by name. `snap.steps`
    // is in `StepIdx` order, so a step's index is its id. A step's input is "full"
    // if ANY input edge (whose consumer is this step) is Full/Backpressured, and
    // its output is "empty" if ANY output edge (whose producer is this step) is
    // Empty/Starved — aggregating over all matching edges rather than the first
    // name match, so fan-out / fan-in / duplicate step names aren't misattributed.
    let input_full = |step: usize| {
        snap.edges
            .iter()
            .filter(|e| e.consumer_step == Some(step))
            .any(|e| matches!(e.class, Full | Backpressured))
    };
    let output_empty = |step: usize| {
        snap.edges
            .iter()
            .filter(|e| e.producer_step == step)
            .any(|e| matches!(e.class, Empty | Starved))
    };

    // Primary bottleneck: full input edge + empty output edge.
    for (step, (name, s)) in snap.steps.iter().enumerate() {
        if input_full(step) && output_empty(step) {
            let cpu_share =
                if total_cpu_ns > 0 { s.total_run_ns as f64 / total_cpu_ns as f64 } else { 0.0 };
            let cause = if cpu_share >= 0.30 {
                format!("CPU-bound ({:.0}% of pipeline try_run time)", cpu_share * 100.0)
            } else {
                "coordination-bound (low CPU share — likely a serialization stall)".to_string()
            };
            findings.push(Finding {
                severity: Severity::Primary,
                message: format!(
                    "BOTTLENECK: step `{name}` (input edge full, output edge empty) — {cause}"
                ),
            });
        }
    }

    // Secondary: Serial-step spin (contention thrash). Skip Detached steps —
    // their dedicated-thread backoff loop records a NoProgress/Contention on
    // every idle poll, which inflates the contention ratio, and "Detach
    // candidate" is meaningless for a step that is already Detached.
    let detached_names: std::collections::HashSet<&str> =
        snap.detached.iter().map(|(name, ..)| *name).collect();
    for (name, s) in &snap.steps {
        if detached_names.contains(name) {
            continue;
        }
        if s.try_run_total > 0 {
            let spin = s.contention_count as f64 / s.try_run_total as f64;
            if spin >= SPIN_THRESHOLD {
                findings.push(Finding {
                    severity: Severity::Secondary,
                    message: format!(
                        "SPIN: step `{name}` contended on {:.0}% of dispatches — affinity / fuse / Detach candidate",
                        spin * 100.0
                    ),
                });
            }
        }
    }

    // Secondary: starvation (consumer of an edge frequently finds it empty).
    for e in &snap.edges {
        if matches!(e.class, Starved) {
            findings.push(Finding {
                severity: Severity::Secondary,
                message: format!(
                    "STARVATION: `{}` is starved by upstream `{}` (empty {:.0}% of pops) — look upstream",
                    e.consumer.unwrap_or("(sink)"),
                    e.producer,
                    e.empty_rate * 100.0
                ),
            });
        }
    }

    // Info: no clear bottleneck and no full/empty edge → latency-bound.
    let any_extreme =
        snap.edges.iter().any(|e| matches!(e.class, Full | Backpressured | Empty | Starved));
    if findings.is_empty() && !snap.edges.is_empty() && !any_extreme {
        findings.push(Finding {
            severity: Severity::Info,
            message: "No throughput bottleneck — edges have headroom; the chain is \
                      latency/coordination-bound (look at residence time, not service time)."
                .to_string(),
        });
    }
    findings
}

impl fmt::Display for StatsSnapshot {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Pipeline stats ({} step{}):", self.steps.len(), pluralize(self.steps.len()))?;
        // Total dispatch time across steps — the denominator for each step's
        // `cpu%` share (the Amdahl ranking: which step to optimize first) and
        // for the headroom line.
        let total_cpu_ns: u64 = self.steps.iter().map(|(_, s)| s.total_run_ns).sum();
        self.write_steps(f, total_cpu_ns)?;
        self.write_utilization(f, total_cpu_ns)?;
        // Per-edge throughput / occupancy / latency + the bottleneck verdict
        // (both only when instrumented).
        self.write_edges(f)?;
        self.write_verdict(f)
    }
}

impl StatsSnapshot {
    /// Render the per-step counter table, including each step's `cpu%` share of
    /// total dispatch time (the per-step cost ranking). `total_cpu_ns` is the
    /// sum of all steps' `total_run_ns`, passed in to avoid recomputing.
    fn write_steps(&self, f: &mut fmt::Formatter<'_>, total_cpu_ns: u64) -> fmt::Result {
        writeln!(
            f,
            "  {:<28} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>14} {:>7} {:>12} {:>12}",
            "step",
            "tries",
            "progress",
            "noprog",
            "content",
            "fin",
            "err",
            "total_ms",
            "cpu%",
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
            // Share of total dispatch time — the per-step cost ranking. The
            // chain-order rows stay readable; this column gives the magnitude.
            #[allow(clippy::cast_precision_loss)]
            let cpu_pct = if total_cpu_ns == 0 {
                0.0
            } else {
                s.total_run_ns as f64 / total_cpu_ns as f64 * 100.0
            };
            writeln!(
                f,
                "  {:<28} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>14.3} {:>6.1}% {:>12} {:>12}",
                name,
                s.try_run_total,
                s.progress_count,
                s.no_progress_count,
                s.contention_count,
                s.finished_count,
                s.error_count,
                total_ms,
                cpu_pct,
                first_str,
                last_str,
            )?;
        }
        Ok(())
    }

    /// Render per-worker utilisation, the pool-level utilisation, and the
    /// squeeze-headroom synthesis. No-op when no worker recorded activity.
    ///
    /// Per-worker: busy = dispatch time, idle = backoff-sleep time. A high
    /// idle% (cores parked while one worker drives a Serial step) is the
    /// signature of pool under-utilisation. The headroom line reads straight off
    /// these numbers: the pool idle% is recoverable via better step overlap;
    /// once the pool saturates (no idle left) the only remaining wins are fewer
    /// cycles/item or more cores, and the hottest step is the Amdahl target.
    fn write_utilization(&self, f: &mut fmt::Formatter<'_>, total_cpu_ns: u64) -> fmt::Result {
        if self.workers.is_empty() {
            // No pool workers recorded activity, but a Detached thread may still
            // have run (e.g. a degenerate chain). Render its line if so.
            self.write_detached(f)?;
            return Ok(());
        }
        #[allow(clippy::cast_precision_loss)]
        let ms = |ns: u64| (ns as f64) / 1_000_000.0;
        let (mut sum_busy, mut sum_idle) = (0u64, 0u64);
        writeln!(f, "  {:<8} {:>14} {:>14} {:>8}", "worker", "busy_ms", "idle_ms", "busy%")?;
        for &(id, busy, idle) in &self.workers {
            sum_busy += busy;
            sum_idle += idle;
            let pct = if busy + idle == 0 {
                0.0
            } else {
                #[allow(clippy::cast_precision_loss)]
                let p = busy as f64 / (busy + idle) as f64 * 100.0;
                p
            };
            writeln!(f, "  {id:<8} {:>14.3} {:>14.3} {pct:>7.1}%", ms(busy), ms(idle))?;
        }
        let pool_pct = if sum_busy + sum_idle == 0 {
            0.0
        } else {
            #[allow(clippy::cast_precision_loss)]
            let p = sum_busy as f64 / (sum_busy + sum_idle) as f64 * 100.0;
            p
        };
        writeln!(
            f,
            "  pool utilisation: {pool_pct:.1}% busy across {} worker{} ({:.3} ms busy / {:.3} ms idle)",
            self.workers.len(),
            pluralize(self.workers.len()),
            ms(sum_busy),
            ms(sum_idle),
        )?;
        if let Some((hot_name, hot)) = self.steps.iter().max_by_key(|(_, s)| s.total_run_ns) {
            #[allow(clippy::cast_precision_loss)]
            let hot_pct = if total_cpu_ns == 0 {
                0.0
            } else {
                hot.total_run_ns as f64 / total_cpu_ns as f64 * 100.0
            };
            writeln!(
                f,
                "  headroom: {:.1}% pool idle (recoverable via overlap); hottest step `{hot_name}` = {hot_pct:.1}% of dispatch time (Amdahl target once pool saturates)",
                100.0 - pool_pct,
            )?;
        }
        self.write_detached(f)
    }

    /// Render the `StepKind::Detached` threads' busy/idle on their own line(s).
    /// These threads are NOT pool workers (legacy "N + 2"), so they are reported
    /// separately and never folded into the pool utilisation %. No-op when no
    /// Detached thread recorded activity (every non-sort chain).
    fn write_detached(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.detached.is_empty() {
            return Ok(());
        }
        #[allow(clippy::cast_precision_loss)]
        let ms = |ns: u64| (ns as f64) / 1_000_000.0;
        for &(name, busy, idle, parks) in &self.detached {
            let pct = if busy + idle == 0 {
                0.0
            } else {
                #[allow(clippy::cast_precision_loss)]
                let p = busy as f64 / (busy + idle) as f64 * 100.0;
                p
            };
            // avg park = idle / parks: a few long parks (idle at end) is benign;
            // many short parks during productive phases is the backoff-latency
            // signal (the deferred per-slot condvar would eliminate them).
            #[allow(clippy::cast_precision_loss)]
            let avg_park_us = if parks == 0 { 0.0 } else { (idle as f64 / parks as f64) / 1000.0 };
            writeln!(
                f,
                "  detached `{name}`: {:.3} ms busy / {:.3} ms idle ({pct:.1}% busy, off pool); {parks} parks (avg {avg_park_us:.1}µs)",
                ms(busy),
                ms(idle),
            )?;
        }
        Ok(())
    }

    /// Render the mechanically-derived bottleneck verdict. No-op when there are
    /// no instrumented edges (nothing to diagnose).
    fn write_verdict(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.edges.is_empty() {
            return Ok(());
        }
        let findings = bottleneck_verdict(self);
        if findings.is_empty() {
            return Ok(());
        }
        writeln!(f, "Bottleneck verdict:")?;
        for finding in findings {
            writeln!(f, "  {}", finding.message)?;
        }
        Ok(())
    }

    /// Render the per-edge table (throughput / occupancy class / rates /
    /// latency). No-op when there are no instrumented edges.
    fn write_edges(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.edges.is_empty() {
            return Ok(());
        }
        writeln!(f, "Pipeline edges ({}):", self.edges.len())?;
        writeln!(
            f,
            "  {:<40} {:>12} {:>10} {:>14} {:>6} {:>6} {:>10}",
            "producer→consumer", "items/s", "MiB/s", "class", "rej%", "empt%", "lat_ms",
        )?;
        for e in &self.edges {
            let edge = format!("{}→{}", e.producer, e.consumer.unwrap_or("(none)"));
            let lat = e.derived_latency_ms.map_or_else(|| "-".to_string(), |l| format!("{l:.2}"));
            writeln!(
                f,
                "  {:<40} {:>12.0} {:>10.1} {:>14} {:>5.1}% {:>5.1}% {:>10}",
                edge,
                e.items_per_s,
                e.mibytes_per_s,
                format!("{:?}", e.class),
                e.reject_rate * 100.0,
                e.empty_rate * 100.0,
                lat,
            )?;
        }
        writeln!(
            f,
            "  (note: throughput is depressed under tracing — confirm wall/RSS with --pipeline-trace off)"
        )
    }
}

fn pluralize(n: usize) -> &'static str {
    if n == 1 { "" } else { "s" }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

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

    #[allow(clippy::too_many_arguments)] // test builder: one arg per snapshot field
    fn edge_ms(
        pushed_items: u64,
        pushed_bytes: u64,
        popped_items: u64,
        pop_empties: u64,
        push_rejections: u64,
        raw: crate::runtime::metrics::RawOccupancy,
        mean_occupancy: f32,
        mean_occupancy_bytes: f64,
    ) -> crate::runtime::metrics::EdgeMetricsSnapshot {
        crate::runtime::metrics::EdgeMetricsSnapshot {
            pushed_items,
            pushed_bytes,
            popped_items,
            popped_bytes: pushed_bytes,
            push_rejections,
            pop_empties,
            depth_samples: 100,
            raw_occupancy: raw,
            mean_occupancy,
            mean_occupancy_bytes,
        }
    }

    #[test]
    fn edge_stats_byte_bounded_latency_is_dimensionally_correct() {
        use crate::runtime::metrics::RawOccupancy;
        // 1000 items/s; mean_item_bytes = 100_000/1000 = 100; sampled mean
        // occupancy 1000 bytes → mean_items 10 → latency = 1000*10/1000 = 10ms.
        // (The byte figure is sampled directly, so the limit only gates that the
        // edge is byte-bounded — its value no longer feeds the latency.)
        let ms = edge_ms(1000, 100_000, 1000, 0, 0, RawOccupancy::Healthy, 0.5, 1000.0);
        let e = compute_edge_stats("p", Some("c"), 0, Some(1), &ms, Some(2000), 1_000_000_000);
        assert!((e.items_per_s - 1000.0).abs() < 1.0, "items/s ≈ 1000, got {}", e.items_per_s);
        let lat = e.derived_latency_ms.expect("byte edge has derived latency");
        assert!((lat - 10.0).abs() < 0.5, "latency ≈ 10ms, got {lat}");
    }

    #[test]
    fn edge_stats_count_edge_has_no_latency() {
        use crate::runtime::metrics::RawOccupancy;
        let ms = edge_ms(500, 0, 500, 0, 0, RawOccupancy::Unknown, 0.0, 0.0);
        let e = compute_edge_stats("p", Some("c"), 0, Some(1), &ms, None, 1_000_000_000);
        assert!(e.derived_latency_ms.is_none(), "count/unbounded edge → no derived latency");
        assert!((e.items_per_s - 500.0).abs() < 1.0);
    }

    // The raw occupancy class is refined by the empty-rate and reject-rate the
    // edge actually observed; each case pins one (raw class, rates, byte bound)
    // combination to its refined `OccupancyClass`.
    #[rstest]
    // MostlyEmpty + high empty-rate (90/100 pops empty) → Starved.
    #[case::starved(
        edge_ms(10, 0, 10, 90, 0, crate::runtime::metrics::RawOccupancy::MostlyEmpty, 0.0, 0.0),
        None,
        OccupancyClass::Starved
    )]
    // MostlyEmpty + low empty-rate → Empty (idle, not starved).
    #[case::empty(
        edge_ms(100, 0, 100, 0, 0, crate::runtime::metrics::RawOccupancy::MostlyEmpty, 0.0, 0.0),
        None,
        OccupancyClass::Empty
    )]
    // MostlyFull + high reject-rate → Backpressured.
    #[case::backpressured(
        edge_ms(10, 0, 10, 0, 90, crate::runtime::metrics::RawOccupancy::MostlyFull, 1.0, 1000.0),
        Some(1000),
        OccupancyClass::Backpressured
    )]
    fn occupancy_class_refined_by_rates(
        #[case] ms: crate::runtime::metrics::EdgeMetricsSnapshot,
        #[case] limit_bytes: Option<u64>,
        #[case] expected: OccupancyClass,
    ) {
        assert_eq!(
            compute_edge_stats("p", None, 0, None, &ms, limit_bytes, 1_000_000_000).class,
            expected
        );
    }

    fn step_stat(
        name: &'static str,
        total_run_ns: u64,
        contention: u64,
        tries: u64,
    ) -> (&'static str, StepStatsSnapshot) {
        (
            name,
            StepStatsSnapshot {
                try_run_total: tries,
                progress_count: tries,
                no_progress_count: 0,
                contention_count: contention,
                finished_count: 1,
                error_count: 0,
                total_run_ns,
                first_progress_ns: 0,
                last_progress_ns: total_run_ns,
            },
        )
    }

    fn edge_stat(
        producer: &'static str,
        consumer: Option<&'static str>,
        producer_step: usize,
        consumer_step: Option<usize>,
        class: OccupancyClass,
        empty_rate: f64,
    ) -> EdgeStatsSnapshot {
        EdgeStatsSnapshot {
            producer,
            consumer,
            producer_step,
            consumer_step,
            items_per_s: 1000.0,
            mibytes_per_s: 1.0,
            class,
            mean_occupancy: 0.5,
            reject_rate: 0.0,
            empty_rate,
            derived_latency_ms: Some(1.0),
        }
    }

    #[test]
    fn verdict_locates_cpu_bound_bottleneck() {
        // `Slow`'s input edge is Full and output edge is Empty, and it dominates
        // CPU → Primary, CPU-bound.
        let snap = StatsSnapshot {
            steps: vec![step_stat("Slow", 1000, 0, 10)],
            workers: vec![],
            detached: vec![],
            // `Slow` is step 0. Its input edge (consumer_step 0) is Full; its
            // output edge (producer_step 0) is Empty. `Up`/`Down` are non-step
            // ids (1/1) so only `Slow` matches by identity.
            edges: vec![
                edge_stat("Up", Some("Slow"), 1, Some(0), OccupancyClass::Full, 0.0),
                edge_stat("Slow", Some("Down"), 0, Some(1), OccupancyClass::Empty, 0.0),
            ],
        };
        let v = bottleneck_verdict(&snap);
        let primary: Vec<_> = v.iter().filter(|f| f.severity == Severity::Primary).collect();
        assert_eq!(primary.len(), 1, "exactly one primary bottleneck");
        assert!(primary[0].message.contains("Slow"));
        assert!(primary[0].message.contains("CPU-bound"));
    }

    #[test]
    fn verdict_matches_edges_by_step_identity_across_fan_in() {
        // Regression for name-based misattribution: `Merge` (step 2) has TWO
        // input edges (fan-in from steps 0 and 1). The first by iteration order
        // is Healthy; the second is Full. Matching the FIRST edge by consumer
        // name would see only the Healthy input and report no bottleneck;
        // matching by step identity aggregates both inputs and correctly flags
        // `Merge` (input full + output empty).
        let snap = StatsSnapshot {
            steps: vec![
                step_stat("A", 100, 0, 10),
                step_stat("B", 100, 0, 10),
                step_stat("Merge", 1000, 0, 10),
            ],
            workers: vec![],
            detached: vec![],
            edges: vec![
                edge_stat("A", Some("Merge"), 0, Some(2), OccupancyClass::Healthy, 0.0),
                edge_stat("B", Some("Merge"), 1, Some(2), OccupancyClass::Full, 0.0),
                edge_stat("Merge", Some("Sink"), 2, Some(3), OccupancyClass::Empty, 0.0),
            ],
        };
        let v = bottleneck_verdict(&snap);
        let primary: Vec<_> = v.iter().filter(|f| f.severity == Severity::Primary).collect();
        assert_eq!(primary.len(), 1, "fan-in bottleneck detected via identity: {v:?}");
        assert!(primary[0].message.contains("Merge"));
    }

    #[test]
    fn verdict_flags_spin_and_starvation() {
        let snap = StatsSnapshot {
            // 5/10 dispatches contended → spin.
            steps: vec![step_stat("Serializer", 100, 5, 10)],
            workers: vec![],
            detached: vec![],
            // A starved edge (steps 1→2, not the `Serializer` step 0): consumer
            // frequently finds it empty.
            edges: vec![edge_stat(
                "Producer",
                Some("Consumer"),
                1,
                Some(2),
                OccupancyClass::Starved,
                0.7,
            )],
        };
        let v = bottleneck_verdict(&snap);
        assert!(
            v.iter().any(|f| f.severity == Severity::Secondary && f.message.contains("SPIN")),
            "spin finding present"
        );
        assert!(
            v.iter().any(|f| f.message.contains("STARVATION") && f.message.contains("Consumer")),
            "starvation finding present"
        );
    }

    #[test]
    fn verdict_skips_spin_for_detached_steps() {
        // A Detached step's backoff loop records Contention on every idle poll,
        // so its contention/tries ratio is high — but "SPIN: … Detach candidate"
        // is nonsensical for an already-Detached step, so it must be skipped.
        let snap = StatsSnapshot {
            // 8/10 "contended" — would trip SPIN if it were a pool step.
            steps: vec![step_stat("SortMerge", 100, 8, 10)],
            workers: vec![],
            detached: vec![("SortMerge", 900_000_000, 100_000_000, 4_200)],
            edges: vec![],
        };
        let v = bottleneck_verdict(&snap);
        assert!(
            !v.iter().any(|f| f.message.contains("SPIN")),
            "no SPIN finding for a Detached step: {v:?}"
        );
    }

    #[test]
    fn verdict_reports_latency_bound_when_no_extremes() {
        let snap = StatsSnapshot {
            steps: vec![step_stat("A", 100, 0, 10)],
            workers: vec![],
            detached: vec![],
            edges: vec![edge_stat("A", Some("B"), 0, Some(1), OccupancyClass::Healthy, 0.0)],
        };
        let v = bottleneck_verdict(&snap);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].severity, Severity::Info);
        assert!(v[0].message.contains("latency"));
    }

    #[test]
    fn snapshot_display_renders_without_edges() {
        // Back-compat: the edge-less snapshot() still renders (empty edges).
        let stats = PipelineStats::new(vec!["A"]);
        let out = format!("{}", stats.snapshot());
        assert!(out.contains("Pipeline stats"));
        assert!(!out.contains("Pipeline edges"), "no edge section when edges empty");
    }

    #[test]
    fn display_renders_cpu_share_and_headroom() {
        // Hot owns 750/1000 = 75% of dispatch time; Cool owns the rest. The pool
        // is 80% busy (800 ms / 1000 ms) → 20% idle headroom.
        let snap = StatsSnapshot {
            steps: vec![step_stat("Hot", 750, 0, 10), step_stat("Cool", 250, 0, 10)],
            workers: vec![(0, 800, 200)],
            detached: vec![],
            edges: vec![],
        };
        let out = format!("{snap}");
        assert!(out.contains("cpu%"), "per-step table has a cpu% column: {out}");
        assert!(out.contains("75.0%"), "Hot step shows its 75% CPU share: {out}");
        assert!(out.contains("pool utilisation: 80.0%"), "pool utilisation line present: {out}");
        assert!(
            out.contains("headroom:") && out.contains("20.0% pool idle") && out.contains("`Hot`"),
            "headroom line names the hottest step and the idle %: {out}"
        );
    }

    #[test]
    fn utilization_section_absent_without_workers() {
        // No worker activity → no pool/headroom lines (fused or stats-off runs).
        let snap = StatsSnapshot {
            steps: vec![step_stat("Solo", 100, 0, 10)],
            workers: vec![],
            detached: vec![],
            edges: vec![],
        };
        let out = format!("{snap}");
        assert!(!out.contains("pool utilisation"), "no pool line without workers: {out}");
        assert!(!out.contains("headroom:"), "no headroom line without workers: {out}");
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

    /// L2.4: a Detached thread's busy/idle is recorded separately, renders on
    /// its own line, and is EXCLUDED from the pool utilisation % (legacy
    /// "N + 2"). The snapshot must render without panicking.
    #[test]
    fn detached_busy_excluded_from_pool_and_on_own_line() {
        let stats = PipelineStats::new(vec!["Source", "SortMerge"]);
        // One pool worker: 800 ms busy / 200 ms idle → pool% = 80%.
        stats.record_worker_busy(0, 800_000_000);
        stats.record_worker_idle(0, 200_000_000);
        // The Detached merge: 900 ms busy / 100 ms idle (off pool).
        stats.record_detached_busy(StepIdx(1), 900_000_000);
        stats.record_detached_idle(StepIdx(1), 100_000_000);

        let snap = stats.snapshot();
        assert_eq!(snap.detached.len(), 1, "one Detached thread recorded");
        assert_eq!(snap.detached[0].0, "SortMerge");

        let out = format!("{snap}");
        // Pool% reflects ONLY the worker (80%), NOT the Detached thread (90%).
        assert!(
            out.contains("pool utilisation: 80.0%"),
            "pool% must exclude the Detached thread: {out}"
        );
        // The Detached thread is reported on its own line.
        assert!(
            out.contains("detached `SortMerge`") && out.contains("90.0% busy, off pool"),
            "Detached line present with its own busy%: {out}"
        );
    }

    /// L2.4: the Detached line renders even when no pool worker recorded
    /// activity (the `write_utilization` early-return path), without panicking.
    #[test]
    fn detached_line_renders_with_no_pool_workers() {
        let stats = PipelineStats::new(vec!["OnlyDetached"]);
        stats.record_detached_busy(StepIdx(0), 5_000_000);
        let out = format!("{}", stats.snapshot());
        assert!(out.contains("detached `OnlyDetached`"), "Detached line present: {out}");
        assert!(!out.contains("pool utilisation"), "no pool line without workers: {out}");
    }
}
