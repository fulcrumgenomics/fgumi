//! Background occupancy sampler for `--pipeline-trace`.
//!
//! When instrumentation is on, `Pipeline::run` spawns one
//! [`run_occupancy_sampler`] thread (same lifecycle slot as the deadlock
//! monitor / queue rebalancer) that periodically reads each byte-bounded edge's
//! depth (`current_bytes / limit_bytes`) and feeds it to the edge's
//! [`EdgeMetrics`](super::metrics::EdgeMetrics) occupancy histogram. Only edges
//! with a `depth_source` (byte-bounded) are sampled — count/unbounded edges still
//! get their push/pop counters, just no occupancy histogram.
//!
//! Reads are cheap but not free. A direct byte-bounded edge costs one `Relaxed`
//! load per tick. An **ordered** edge additionally reads its `ReorderStage`
//! overflow stash, and `ReorderCapHandle::current_buffer_bytes` takes the stage's
//! `state` mutex — the same one every must-accept push and every
//! `try_pop_in_order` holds. So an ordered edge does briefly touch a worker-hot
//! lock once per tick (at the default interval, negligible against per-item
//! traffic, but not zero).
//!
//! Each tick therefore reads every edge's depth **once**, via `read_depths`,
//! and hands the result to both consumers (`record_depths` and the timeline
//! writer). Letting each consumer read for itself would take that mutex twice per
//! ordered edge per tick and could record two different depths for one tick.

use std::fmt::Write as _;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

use super::contexts::{RegisteredEdge, StepCounters};
use super::metrics::{EdgeMetricsSnapshot, RawOccupancy};
use super::telemetry::{EdgeSample, SummarySample, TelemetryWriters, WorkerBins, WorkerSample};
use super::worker_state::{WorkerState, WorkerStateBoard};

/// Default sampling interval. Two milliseconds is cheap (one atomic load per
/// edge) yet fine-grained enough to resolve a chain's phase structure over a
/// multi-second run.
pub const DEFAULT_SAMPLE_INTERVAL: Duration = Duration::from_millis(2);

/// Floor for the effective sampling interval. A zero (or near-zero) interval
/// would turn the sampler's `while !stop { …; sleep(interval) }` loop into a hot
/// spin that pegs a core and grows the `.ticks.*.tsv` files without bound. The
/// floor is applied inside [`run_occupancy_sampler`] so it holds for EVERY
/// caller — the CLI/env parse, and a direct library `TelemetryConfig` whose
/// `interval` bypasses the CLI — making the busy-spin unreachable. 100µs is finer
/// than the 2ms default yet bounded away from zero.
pub const MIN_SAMPLE_INTERVAL: Duration = Duration::from_micros(100);

/// Extra tick telemetry driven by the occupancy sampler loop, additive to the
/// existing per-edge occupancy histogram + `pipeline-trace.tsv` timeline path
/// (that path is untouched whether or not `tele` is passed).
///
/// `source_edge_idxs` / `sink_edge_idxs` are indices into the same `edges`
/// slice passed to [`run_occupancy_sampler`] — the edges whose producer is a
/// chain source and whose consumer is the terminal sink, respectively. They
/// are computed by the caller (`Pipeline::run`); the sampler only consumes
/// them.
pub struct TelemetryArgs<'a> {
    /// Per-OS-thread `(state, step)` board, one slot per pipeline thread.
    pub board: &'a WorkerStateBoard,
    /// Number of worker rows to emit per tick (must match `worker_slots.len()`).
    pub n_workers: usize,
    /// Count of pool (work-stealing) threads, which occupy the low slot range
    /// `0..n_pool`; slots at or above it are detached-driver threads. The sampler
    /// uses this boundary to label each worker row's `role` column
    /// (`pool` / `detached`).
    pub n_pool: usize,
    /// Step names indexed by `StepIdx`, for the `workers` TSV's per-step
    /// fraction columns and the `steps` TSV.
    pub step_names: &'a [&'static str],
    /// Optional RSS probe for the summary row's `rss_bytes` column; `None`
    /// when no probe is wired up.
    pub rss_probe: Option<&'a (dyn Fn() -> Option<u64> + Send + Sync)>,
    /// Total queue-byte budget for the summary row's `queue_bytes_budget`
    /// column, if known.
    pub queue_bytes_budget: Option<u64>,
    /// Open TSV writers for the six `<stem>.ticks.*.tsv` files.
    pub writers: TelemetryWriters,
    /// Indices into `edges` whose producer is a chain source; their
    /// (cumulative) `pushed_items` sum to the summary row's `reads_in`.
    pub source_edge_idxs: Vec<usize>,
    /// Indices into `edges` whose consumer is the terminal sink; their
    /// (cumulative) `popped_items` sum to the summary row's `reads_out`.
    pub sink_edge_idxs: Vec<usize>,
    /// Worker index (row number in the `workers` TSV) -> `WorkerStateBoard`
    /// slot index.
    pub worker_slots: Vec<usize>,
    /// The pipeline's per-step domain counters (`&contexts.step_counters`), one
    /// entry per step index. The sampler reads each slot once per tick and emits
    /// the per-tick delta + cumulative total to the `counters` TSV.
    pub step_counters: &'a [StepCounters],
}

/// A zero-valued edge snapshot, used to seed the tick telemetry's per-edge
/// "previous" store so tick 0's deltas are the totals since run start (i.e.
/// since the edge's counters were created), not since the sampler thread's
/// first tick.
fn zero_edge_snapshot() -> EdgeMetricsSnapshot {
    EdgeMetricsSnapshot {
        pushed_items: 0,
        pushed_bytes: 0,
        popped_items: 0,
        popped_bytes: 0,
        push_rejections: 0,
        pop_empties: 0,
        depth_samples: 0,
        raw_occupancy: RawOccupancy::Unknown,
        mean_occupancy: 0.0,
        mean_occupancy_bytes: 0.0,
    }
}

/// Mutable per-run state for the tick telemetry: one [`WorkerBins`] per
/// worker (per-step occupancy over the emit window), a parallel
/// Running-sample counter per worker (`d_serviced` — see below), the last
/// point-sampled `(state, step)` per worker (reused for the row so the
/// sampler doesn't read the board twice per tick), the previous tick's edge
/// counter snapshots (for computing deltas), the tick index, and the
/// `Instant` of the last emit (for `dt_ms`).
///
/// `d_serviced` is a **sample-derived activity proxy, not an exact per-worker
/// item count**: exact counts would require a step-side counter that doesn't
/// exist yet (a deferred future upgrade). Here it is simply the number of
/// `Running` point-samples taken for that worker during the emit window — at
/// v1's 1:1 sample:emit cadence that is 0 or 1 per tick, but the field stays
/// meaningful if a future cadence samples faster than it emits.
struct TelemetryState<'a> {
    args: TelemetryArgs<'a>,
    bins: Vec<WorkerBins>,
    d_serviced: Vec<u64>,
    last_read: Vec<(WorkerState, Option<crate::topology::StepIdx>)>,
    prev_edges: Vec<EdgeMetricsSnapshot>,
    /// Previous tick's cumulative value per (step, counter), for computing the
    /// per-tick `d_value` deltas. Seeded to 0 (like `prev_edges`) so tick 0's
    /// delta is the total since run start, not since the first sampler tick.
    prev_counters: Vec<Vec<u64>>,
    /// Per-edge row buffer reused across ticks. The `edge` column prefix is
    /// invariant for the run, so it is built once here (in `new`) and left
    /// untouched by `emit`, which only overwrites the numeric fields — sparing
    /// both the per-tick `Vec` allocation and the per-edge `String` format at
    /// the 100 µs floor (~10k ticks/s).
    edge_rows: Vec<EdgeSample>,
    tick: u64,
    last_emit: Instant,
}

impl<'a> TelemetryState<'a> {
    fn new(args: TelemetryArgs<'a>, edges: &[RegisteredEdge]) -> Self {
        let n_steps = args.step_names.len();
        let bins = (0..args.n_workers).map(|_| WorkerBins::new(n_steps)).collect();
        let d_serviced = vec![0u64; args.n_workers];
        let last_read = vec![(WorkerState::Idle, None); args.n_workers];
        let prev_edges = vec![zero_edge_snapshot(); edges.len()];
        let prev_counters = args.step_counters.iter().map(|sc| vec![0u64; sc.len()]).collect();
        // The edge column prefix is invariant for the run, so build each row
        // once with its prefix and zeroed numerics; `emit` overwrites only the
        // numeric fields per tick and never re-formats the prefix.
        let edge_rows = edges
            .iter()
            .map(|e| EdgeSample {
                edge: edge_column_prefix(e),
                depth_bytes: None,
                limit_bytes: None,
                d_pushed: 0,
                d_popped: 0,
                d_push_rej: 0,
                d_pop_empty: 0,
                d_pushed_bytes: None,
                d_popped_bytes: None,
            })
            .collect();
        Self {
            args,
            bins,
            d_serviced,
            last_read,
            prev_edges,
            prev_counters,
            edge_rows,
            tick: 0,
            last_emit: Instant::now(),
        }
    }

    /// Point-sample every worker slot's current `(state, step)` into its
    /// `WorkerBins`, bump the `d_serviced` proxy on `Running`, and stash the
    /// read for reuse by [`Self::emit`]'s worker row (one board read per
    /// worker per tick, shared the same way `read_depths` is shared between
    /// the occupancy histogram and the timeline row).
    fn sample_workers(&mut self) {
        for (w, &slot) in self.args.worker_slots.iter().enumerate() {
            let (state, step) = self.args.board.read(slot);
            if let Some(bins) = self.bins.get_mut(w) {
                bins.record(state, step);
            }
            if state == WorkerState::Running
                && let Some(count) = self.d_serviced.get_mut(w)
            {
                *count += 1;
            }
            if let Some(slot_read) = self.last_read.get_mut(w) {
                *slot_read = (state, step);
            }
        }
    }

    /// Compute this tick's edge deltas / `reads_in` / `reads_out` /
    /// `queue_bytes_used`, write the summary/edge/worker rows, and reset the
    /// per-window accumulators (bins, `d_serviced`, previous edge snapshots).
    ///
    /// `depths` is the SAME per-tick read `run_occupancy_sampler` already
    /// took via `read_depths` for the occupancy histogram / timeline row —
    /// reused here rather than re-read, for the same reason the timeline
    /// writer reuses it (one edge read per tick, not one per consumer).
    fn emit(&mut self, edges: &[RegisteredEdge], depths: &[EdgeDepth], t_ms: f64) {
        self.emit_summary_and_edges(edges, depths, t_ms);
        self.emit_worker_states(t_ms);
        self.emit_counters(t_ms);
        self.tick += 1;
    }

    /// Terminal finalize for the shutdown path (see `run_occupancy_sampler`'s
    /// final-sample guard). Emits the FINAL cumulative summary/edge/counter
    /// values so edge and `StepCounters` activity during the last sleep window
    /// is not dropped from the TSVs, WITHOUT sampling worker state or taking
    /// another occupancy point — a final worker-state row (its per-window bins
    /// reset and never re-sampled after the last tick) and a final occupancy
    /// point would both bias those per-window distributions toward the drained
    /// final state. The caller pairs this with skipping `record_depths`, the
    /// timeline row, and `sample_workers`: only the monotonic totals are
    /// finalized here, not the per-window distributions.
    fn emit_terminal(&mut self, edges: &[RegisteredEdge], depths: &[EdgeDepth], t_ms: f64) {
        self.emit_summary_and_edges(edges, depths, t_ms);
        self.emit_counters(t_ms);
        self.tick += 1;
    }

    /// Write the per-tick summary row plus one edge row per edge (delta metrics
    /// and current depth), advancing the per-window edge snapshots. Shared by
    /// the regular tick (`emit`) and the terminal finalize (`emit_terminal`).
    fn emit_summary_and_edges(
        &mut self,
        edges: &[RegisteredEdge],
        depths: &[EdgeDepth],
        t_ms: f64,
    ) {
        let dt_ms = self.last_emit.elapsed().as_secs_f64() * 1000.0;
        self.last_emit = Instant::now();

        let mut reads_in = 0u64;
        let mut reads_out = 0u64;
        let mut queue_bytes_used = 0u64;
        for (i, e) in edges.iter().enumerate() {
            let cur = e.metrics.snapshot();
            let prev = self.prev_edges[i];
            if let Some(src) = &e.depth_source {
                // `queue_bytes_used` = Σ transport `current_bytes()` only; it does
                // NOT include an ordered edge's reorder-stash bytes, whereas an
                // edge's `depth_bytes` (from `read_depths`) DOES add the stash. So
                // for ordered edges Σ `depth_bytes` ≥ `queue_bytes_used` — the two
                // intentionally won't reconcile exactly.
                queue_bytes_used += src.current_bytes();
            }
            if self.args.source_edge_idxs.contains(&i) {
                reads_in += cur.pushed_items;
            }
            if self.args.sink_edge_idxs.contains(&i) {
                reads_out += cur.popped_items;
            }
            let (depth_bytes, limit_bytes) =
                depths.get(i).copied().flatten().map_or((None, None), |(o, l)| (Some(o), Some(l)));
            // Byte deltas are NA (None) exactly when the edge has no byte
            // accounting — i.e. a count/unbounded edge with no `depth_source`.
            let (d_pushed_bytes, d_popped_bytes) = if e.depth_source.is_some() {
                (
                    Some(cur.pushed_bytes.saturating_sub(prev.pushed_bytes)),
                    Some(cur.popped_bytes.saturating_sub(prev.popped_bytes)),
                )
            } else {
                (None, None)
            };
            // Overwrite only the numeric fields of the reused row; `edge` holds
            // the invariant prefix built once in `new`.
            let row = &mut self.edge_rows[i];
            row.depth_bytes = depth_bytes;
            row.limit_bytes = limit_bytes;
            row.d_pushed = cur.pushed_items.saturating_sub(prev.pushed_items);
            row.d_popped = cur.popped_items.saturating_sub(prev.popped_items);
            row.d_push_rej = cur.push_rejections.saturating_sub(prev.push_rejections);
            row.d_pop_empty = cur.pop_empties.saturating_sub(prev.pop_empties);
            row.d_pushed_bytes = d_pushed_bytes;
            row.d_popped_bytes = d_popped_bytes;
            self.prev_edges[i] = cur;
        }

        let rss_bytes = self.args.rss_probe.and_then(|f| f());
        let summary = SummarySample {
            reads_in,
            reads_out,
            rss_bytes,
            queue_bytes_used,
            queue_bytes_budget: self.args.queue_bytes_budget,
        };
        self.args.writers.write_summary_row(self.tick, t_ms, dt_ms, &summary);
        for row in &self.edge_rows {
            self.args.writers.write_edge_row(self.tick, t_ms, row);
        }
    }

    /// Write one worker-state row per worker from the accumulated per-window
    /// occupancy bins, then reset those bins and the serviced counters for the
    /// next window. Omitted from `emit_terminal` to avoid biasing the
    /// worker-state distribution toward the drained final state.
    fn emit_worker_states(&mut self, t_ms: f64) {
        for w in 0..self.args.n_workers {
            let (state, step) = self.last_read.get(w).copied().unwrap_or((WorkerState::Idle, None));
            let fractions = self.bins[w].fractions();
            let role = if w < self.args.n_pool { "pool" } else { "detached" };
            let sample = WorkerSample {
                worker: w,
                role,
                state,
                step,
                d_serviced: self.d_serviced[w],
                samples: self.bins[w].total(),
                fractions,
            };
            self.args.writers.write_worker_row(self.tick, t_ms, &sample);
            self.bins[w].reset();
            self.d_serviced[w] = 0;
        }
    }

    /// Emit domain counters: one row per (step, counter) whose cumulative
    /// `value` (read from the shared per-step atomic) and per-tick delta
    /// (`d_value`) are written. Ordered step asc then counter asc. Header-only
    /// file when no step declared a counter. Shared by `emit` and
    /// `emit_terminal` so the final cumulative totals are never dropped.
    fn emit_counters(&mut self, t_ms: f64) {
        for (step, counters) in self.args.step_counters.iter().enumerate() {
            for counter in 0..counters.len() {
                let cur = counters.load(counter);
                let prev = self.prev_counters[step][counter];
                let d_value = cur.saturating_sub(prev);
                self.args.writers.write_counter_row(self.tick, t_ms, step, counter, d_value, cur);
                self.prev_counters[step][counter] = cur;
            }
        }
    }

    fn flush(mut self) {
        self.args.writers.flush();
    }
}

/// Poll each byte-bounded edge's occupancy into its histogram until `stop` is
/// set. Edges without a `depth_source` (count/unbounded) are skipped. When
/// `trace_path` is `Some` (the `Timeline` level), also append one TSV row per
/// tick — `t_ms` plus, per edge, its depth fraction and cumulative
/// pushed/popped item counts — so the run's phase structure can be plotted.
///
/// When `tele` is `Some`, ALSO drives the tick telemetry (`.ticks.*.tsv`)
/// each tick, additively — see `TelemetryState`. This never disturbs the
/// `read_depths` / `record_depths` / `TraceWriter` path above: both consumers
/// share the same per-tick `read_depths` result.
pub fn run_occupancy_sampler(
    stop: &AtomicBool,
    edges: &[RegisteredEdge],
    interval: Duration,
    trace_path: Option<PathBuf>,
    tele: Option<TelemetryArgs<'_>>,
) {
    // Clamp to a floor so a zero/near-zero interval cannot busy-spin (see
    // `MIN_SAMPLE_INTERVAL`); this covers direct-library callers too, not just the
    // CLI parse.
    let interval = interval.max(MIN_SAMPLE_INTERVAL);
    let mut trace = trace_path.and_then(|p| TraceWriter::open(&p, edges));
    let mut tele_state = tele.map(|args| TelemetryState::new(args, edges));
    let start = Instant::now();
    let mut sampled_in_loop = false;
    while !stop.load(Ordering::Relaxed) {
        // One read per tick, shared by the histogram and the timeline row.
        let depths = read_depths(edges);
        record_depths(edges, &depths);
        if let Some(t) = trace.as_mut() {
            t.write_row(edges, &depths, start.elapsed());
        }
        if let Some(ts) = tele_state.as_mut() {
            ts.sample_workers();
            let t_ms = start.elapsed().as_secs_f64() * 1000.0;
            ts.emit(edges, &depths, t_ms);
        }
        sampled_in_loop = true;
        // Sleep the interval in bounded slices, re-checking `stop` between them,
        // so a large interval does not delay shutdown by its full duration.
        let mut remaining = interval;
        let slice = Duration::from_millis(50);
        while remaining > Duration::ZERO && !stop.load(Ordering::Relaxed) {
            let nap = remaining.min(slice);
            std::thread::sleep(nap);
            remaining -= nap;
        }
    }
    // Shutdown sampling, in two mutually exclusive shapes:
    //
    // * The loop never sampled (a run so short `stop` was set before the first
    //   iteration): take ONE full sample so the TSVs are not empty.
    // * The loop did sample: DON'T take another full sample — an extra
    //   occupancy point + timeline row + worker-state sample taken AFTER the
    //   pipeline already drained would bias those per-window distributions
    //   toward the empty final state. But edge metrics and `StepCounters` can
    //   still advance during the final sleep, after the last regular tick
    //   emitted, so finalize just those monotonic cumulative totals
    //   (`emit_terminal`) — no occupancy point, timeline row, or worker sample.
    if !sampled_in_loop {
        let depths = read_depths(edges);
        record_depths(edges, &depths);
        if let Some(t) = trace.as_mut() {
            t.write_row(edges, &depths, start.elapsed());
        }
        if let Some(ts) = tele_state.as_mut() {
            ts.sample_workers();
            let t_ms = start.elapsed().as_secs_f64() * 1000.0;
            ts.emit(edges, &depths, t_ms);
        }
    } else if let Some(ts) = tele_state.as_mut() {
        // `read_depths` (current depth for the final edge row) but NOT
        // `record_depths` (histogram) — reading is accurate at drain, recording
        // would bias the occupancy distribution.
        let depths = read_depths(edges);
        let t_ms = start.elapsed().as_secs_f64() * 1000.0;
        ts.emit_terminal(edges, &depths, t_ms);
    }
    if let Some(mut t) = trace {
        t.flush();
    }
    if let Some(ts) = tele_state {
        ts.flush();
    }
}

/// Bytes buffered in an ordered edge's `ReorderStage` overflow stash (0 for a
/// direct/count/unbounded edge). Added to the transport queue's `current_bytes`
/// when sampling depth so an ordered edge reflects total buffered bytes rather
/// than reading empty while items pile in the reorder buffer awaiting an earlier
/// ordinal.
fn reorder_stash_bytes(edge: &RegisteredEdge) -> u64 {
    edge.reorder_depth.as_ref().map_or(0, |r| r.current_buffer_bytes())
}

/// One edge's depth for a single tick: `(occupied_bytes, limit_bytes)`, or `None`
/// for a count/unbounded edge (no `depth_source`) or one whose limit reads 0.
type EdgeDepth = Option<(u64, u64)>;

/// Read every edge's depth once, for one tick.
///
/// Called once per tick and shared by [`record_depths`] and the timeline writer
/// so that (a) an ordered edge's `ReorderStage` mutex is taken once per tick
/// rather than once per consumer, and (b) the histogram sample and the timeline
/// row for a given tick always report the same number.
fn read_depths(edges: &[RegisteredEdge]) -> Vec<EdgeDepth> {
    edges
        .iter()
        .map(|edge| {
            let src = edge.depth_source.as_ref()?;
            let limit = src.limit_bytes();
            if limit == 0 {
                return None;
            }
            Some((src.current_bytes().saturating_add(reorder_stash_bytes(edge)), limit))
        })
        .collect()
}

/// Feed one tick's depths into each edge's occupancy histogram.
fn record_depths(edges: &[RegisteredEdge], depths: &[EdgeDepth]) {
    for (edge, depth) in edges.iter().zip(depths) {
        if let Some((occupied, limit)) = *depth {
            edge.metrics.record_depth(occupied, limit);
        }
    }
}

/// Column-name prefix for one edge's timeline columns. Includes the producer
/// step index and output branch so fan-out edges (one producer, several
/// branches) and repeated step names produce distinct, collision-free headers —
/// a bare `producer__consumer` prefix duplicates columns whenever two edges
/// share both names. `(producer_step, branch)` uniquely identifies an edge.
fn edge_column_prefix(e: &RegisteredEdge) -> String {
    format!(
        "{}__{}#{}b{}",
        e.producer_name,
        e.consumer_name.unwrap_or("sink"),
        e.producer_step.0,
        e.branch.0,
    )
}

/// Per-tick TSV writer for the `Timeline` level. Best-effort: a write error is
/// logged once and further rows are dropped (instrumentation never aborts a run).
struct TraceWriter {
    writer: BufWriter<std::fs::File>,
    failed: bool,
}

impl TraceWriter {
    /// Open `path` and write the header (`t_ms` + three columns per edge).
    /// Returns `None` (with a warning) if the file can't be created.
    fn open(path: &std::path::Path, edges: &[RegisteredEdge]) -> Option<Self> {
        match std::fs::File::create(path) {
            Ok(file) => {
                let mut writer = BufWriter::new(file);
                let mut header = String::from("t_ms");
                for e in edges {
                    let edge = edge_column_prefix(e);
                    let _ = write!(header, "\t{edge}.depth\t{edge}.pushed\t{edge}.popped");
                }
                if writeln!(writer, "{header}").is_err() {
                    log::warn!(
                        "pipeline-trace: failed to write timeline header to {}",
                        path.display()
                    );
                    return None;
                }
                Some(Self { writer, failed: false })
            }
            Err(e) => {
                log::warn!("pipeline-trace: cannot create timeline file {}: {e}", path.display());
                None
            }
        }
    }

    /// `depths` is this tick's depths from [`read_depths`], shared with
    /// [`record_depths`] so the row and the histogram agree and each ordered
    /// edge's reorder mutex is taken once per tick.
    #[allow(clippy::cast_precision_loss)]
    fn write_row(&mut self, edges: &[RegisteredEdge], depths: &[EdgeDepth], elapsed: Duration) {
        if self.failed {
            return;
        }
        let mut row = format!("{}", elapsed.as_millis());
        for (e, depth) in edges.iter().zip(depths) {
            // Count/unbounded edges (no `depth_source`) are unsampled — emit `NA`
            // rather than `0.000`, which would misread as "empty" instead of
            // "not measured". Byte-bounded edges report total buffered depth
            // (transport + reorder stash) as a fraction of the limit.
            let depth = depth.map_or_else(
                || "NA".to_string(),
                |(occupied, limit)| format!("{:.3}", occupied as f32 / limit as f32),
            );
            let ms = e.metrics.snapshot();
            let _ = write!(row, "\t{depth}\t{}\t{}", ms.pushed_items, ms.popped_items);
        }
        if writeln!(self.writer, "{row}").is_err() {
            log::warn!("pipeline-trace: timeline write failed; dropping further rows");
            self.failed = true;
        }
    }

    fn flush(&mut self) {
        if self.failed {
            return;
        }
        // A dropped flush error can silently lose buffered rows after every
        // write appeared to succeed — warn and mark the writer failed, matching
        // `write_row`'s best-effort error handling.
        if self.writer.flush().is_err() {
            log::warn!("pipeline-trace: timeline flush failed; buffered rows may be lost");
            self.failed = true;
        }
    }
}

/// One sampling sweep over all edges. Extracted so tests can drive a single
/// deterministic tick without the sleep loop.
///
/// The sampler loop does not call this — it uses `read_depths` once per tick
/// and shares the result with the timeline writer, so both record the same
/// numbers from a single read (see the module doc).
pub fn sample_once(edges: &[RegisteredEdge]) {
    record_depths(edges, &read_depths(edges));
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use std::sync::atomic::AtomicBool;

    use crate::item::HeapSize;
    use crate::queues::{BoundedQueueHandle, ByteBoundedQueue, ItemQueue};
    use crate::runtime::contexts::StepCounters;
    use crate::runtime::metrics::EdgeMetrics;
    use crate::step::CounterSpec;
    use crate::topology::{BranchIdx, StepIdx};

    #[derive(Debug)]
    struct Heavy(Vec<u8>);
    impl HeapSize for Heavy {
        fn heap_size(&self) -> usize {
            self.0.len()
        }
    }

    fn edge_over(
        metrics: Arc<EdgeMetrics>,
        depth_source: Option<Arc<dyn BoundedQueueHandle>>,
    ) -> RegisteredEdge {
        RegisteredEdge {
            producer_step: StepIdx(0),
            producer_name: "producer",
            consumer_step: Some(StepIdx(1)),
            consumer_name: Some("consumer"),
            branch: BranchIdx(0),
            metrics,
            depth_source,
            reorder_depth: None,
        }
    }

    #[test]
    fn sample_once_records_occupancy_from_depth_source() {
        let m = EdgeMetrics::new();
        let q = Arc::new(ByteBoundedQueue::<Heavy>::new(1000));
        q.try_push(Heavy(vec![0; 500])).unwrap(); // 50% of the 1000-byte budget
        let edge = edge_over(Arc::clone(&m), Some(Arc::clone(&q) as Arc<dyn BoundedQueueHandle>));
        for _ in 0..10 {
            sample_once(std::slice::from_ref(&edge));
        }
        let s = m.snapshot();
        assert_eq!(s.depth_samples, 10);
        assert!((s.mean_occupancy - 0.5).abs() < 0.05, "mean ≈ 0.5, got {}", s.mean_occupancy);
    }

    #[test]
    fn count_edge_without_depth_source_is_skipped() {
        // An edge with no depth_source (count/unbounded) records no occupancy.
        let m = EdgeMetrics::new();
        let edge = edge_over(Arc::clone(&m), None);
        for _ in 0..5 {
            sample_once(std::slice::from_ref(&edge));
        }
        assert_eq!(m.snapshot().depth_samples, 0, "no depth source → no occupancy samples");
    }

    #[test]
    fn run_occupancy_sampler_stops_and_records() {
        let m = EdgeMetrics::new();
        let q = Arc::new(ByteBoundedQueue::<Heavy>::new(1000));
        q.try_push(Heavy(vec![0; 800])).unwrap();
        let edges =
            vec![edge_over(Arc::clone(&m), Some(Arc::clone(&q) as Arc<dyn BoundedQueueHandle>))];
        let stop = Arc::new(AtomicBool::new(false));
        let stop_c = Arc::clone(&stop);
        let handle = std::thread::spawn(move || {
            run_occupancy_sampler(&stop_c, &edges, Duration::from_millis(1), None, None);
        });
        std::thread::sleep(Duration::from_millis(30));
        stop.store(true, Ordering::Relaxed);
        handle.join().unwrap();
        let s = m.snapshot();
        assert!(s.depth_samples > 0, "sampler recorded at least one tick");
        assert!((s.mean_occupancy - 0.8).abs() < 0.1, "mean ≈ 0.8, got {}", s.mean_occupancy);
    }

    #[test]
    fn edge_columns_are_unique_when_step_names_collide() {
        // Regression: two edges that share producer AND consumer names (fan-out,
        // or duplicate step names) must still produce distinct TSV columns — a
        // bare `producer__consumer` prefix would emit duplicate column headers.
        let m0 = EdgeMetrics::new();
        let m1 = EdgeMetrics::new();
        let e0 = RegisteredEdge {
            producer_step: StepIdx(0),
            producer_name: "dup",
            consumer_step: Some(StepIdx(1)),
            consumer_name: Some("sink"),
            branch: BranchIdx(0),
            metrics: m0,
            depth_source: None,
            reorder_depth: None,
        };
        // Same names, different (producer_step, branch): a fan-out sibling.
        let e1 = RegisteredEdge {
            producer_step: StepIdx(0),
            producer_name: "dup",
            consumer_step: Some(StepIdx(2)),
            consumer_name: Some("sink"),
            branch: BranchIdx(1),
            metrics: m1,
            depth_source: None,
            reorder_depth: None,
        };
        let p0 = edge_column_prefix(&e0);
        let p1 = edge_column_prefix(&e1);
        assert_ne!(p0, p1, "colliding names must yield distinct column prefixes");
        assert_eq!(p0, "dup__sink#0b0");
        assert_eq!(p1, "dup__sink#0b1");
    }

    #[test]
    fn ordered_edge_occupancy_includes_reorder_stash() {
        use crate::queues::CountBoundedQueue;
        use crate::reorder::{ReorderCapHandle, ReorderStage, Sequenced};
        // The transport (occupancy depth source) is empty, but the reorder stash
        // holds 400 buffered bytes waiting for an earlier ordinal. Sampled
        // occupancy must reflect the stash (400/1000 = 0.4), not read empty —
        // otherwise a producer-skewed ordered edge looks idle while backed up.
        let m = EdgeMetrics::new();
        let transport = Arc::new(ByteBoundedQueue::<Heavy>::new(1000));

        let reorder_transport: Arc<dyn ItemQueue<Sequenced<Heavy>>> =
            Arc::new(CountBoundedQueue::<Sequenced<Heavy>>::new(8));
        let stage = Arc::new(ReorderStage::new(reorder_transport));
        stage.try_push(1, Heavy(vec![0; 400])).unwrap(); // out-of-order → stashed
        assert!(stage.try_pop_in_order().is_none(), "ordinal 0 absent → nothing pops");
        assert!(stage.current_buffer_bytes() >= 400, "stash holds the buffered bytes");

        let edge = RegisteredEdge {
            producer_step: StepIdx(0),
            producer_name: "p",
            consumer_step: Some(StepIdx(1)),
            consumer_name: Some("c"),
            branch: BranchIdx(0),
            metrics: Arc::clone(&m),
            depth_source: Some(Arc::clone(&transport) as Arc<dyn BoundedQueueHandle>),
            reorder_depth: Some(Arc::clone(&stage) as Arc<dyn ReorderCapHandle>),
        };
        for _ in 0..10 {
            sample_once(std::slice::from_ref(&edge));
        }
        let s = m.snapshot();
        assert_eq!(s.depth_samples, 10);
        assert!(
            (s.mean_occupancy - 0.4).abs() < 0.05,
            "occupancy reflects the reorder stash, got {}",
            s.mean_occupancy
        );
        assert!(
            (s.mean_occupancy_bytes - 400.0).abs() < 1.0,
            "byte mean equals the stashed bytes, got {}",
            s.mean_occupancy_bytes
        );
    }

    #[test]
    fn timeline_tsv_has_header_and_rows() {
        let m = EdgeMetrics::new();
        let q = Arc::new(ByteBoundedQueue::<Heavy>::new(1000));
        q.try_push(Heavy(vec![0; 400])).unwrap();
        let edges =
            vec![edge_over(Arc::clone(&m), Some(Arc::clone(&q) as Arc<dyn BoundedQueueHandle>))];
        let dir = std::env::temp_dir();
        let path = dir.join(format!("fgumi-trace-test-{}.tsv", std::process::id()));
        let stop = Arc::new(AtomicBool::new(false));
        let stop_c = Arc::clone(&stop);
        let path_c = path.clone();
        let handle = std::thread::spawn(move || {
            run_occupancy_sampler(&stop_c, &edges, Duration::from_millis(2), Some(path_c), None);
        });
        std::thread::sleep(Duration::from_millis(30));
        stop.store(true, Ordering::Relaxed);
        handle.join().unwrap();

        let content = std::fs::read_to_string(&path).expect("trace file written");
        let _ = std::fs::remove_file(&path);
        let mut lines = content.lines();
        let header = lines.next().expect("header row");
        assert!(header.starts_with("t_ms"), "header begins with t_ms");
        assert!(header.contains("producer__consumer#0b0.depth"), "per-edge depth column");
        let rows: Vec<&str> = lines.collect();
        assert!(!rows.is_empty(), "at least one data row");
        // First field of a data row is a monotonic t_ms integer.
        let first_t: u128 = rows[0].split('\t').next().unwrap().parse().expect("t_ms is an int");
        let last_t: u128 = rows.last().unwrap().split('\t').next().unwrap().parse().unwrap();
        assert!(last_t >= first_t, "t_ms is monotonic");
    }

    #[test]
    fn sampler_emits_telemetry_files() {
        use crate::runtime::telemetry::TelemetryWriters;
        use crate::runtime::worker_state::{WorkerState, WorkerStateBoard};
        // `const` so the inner slice is `'static` (a const-fn call is not
        // auto-promoted to `'static` in an rvalue context).
        const COUNTER_SPECS: &[&[CounterSpec]] = &[&[CounterSpec::new("records", "records")]];
        let m = EdgeMetrics::new();
        m.record_push(100);
        m.record_pop(40);
        let q = Arc::new(ByteBoundedQueue::<Heavy>::new(1000));
        q.try_push(Heavy(vec![0; 300])).unwrap();
        let edges =
            vec![edge_over(Arc::clone(&m), Some(Arc::clone(&q) as Arc<dyn BoundedQueueHandle>))];
        let dir = std::env::temp_dir().join(format!("fgumi-sampler-tele-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let stem = dir.join("run");
        let board = WorkerStateBoard::new(1);
        board.stamp(0, WorkerState::Running, Some(StepIdx(0)));
        // The producer step declares one counter, pre-bumped to 42 so the sampler
        // reads a non-zero cumulative value on its first tick.
        let step_counters = vec![StepCounters::new(1)];
        step_counters[0].add(0, 42);
        let writers = TelemetryWriters::open(&stem, &["producer"], COUNTER_SPECS).unwrap();
        let tele = crate::runtime::sampler::TelemetryArgs {
            board: &board,
            n_workers: 1,
            n_pool: 1,
            step_names: &["producer"],
            rss_probe: None,
            queue_bytes_budget: Some(1000),
            writers,
            source_edge_idxs: vec![0],
            sink_edge_idxs: vec![0],
            worker_slots: vec![0],
            step_counters: &step_counters,
        };
        let stop = Arc::new(AtomicBool::new(false));
        let stop_c = Arc::clone(&stop);
        let edges_c = edges;
        // `thread::scope` (not `thread::spawn`) because `tele.board` borrows the
        // stack-local `board` — `spawn` would require that borrow to be
        // `'static`. This is a test-driver detail only: it does not change
        // `TelemetryArgs`'s lifetime shape, which stays a plain borrow (the
        // real caller in `builder.rs` constructs it from an owned `Arc` moved
        // into its spawned closure, so the borrow there is local to that
        // closure's body, not held across the `spawn` boundary itself).
        std::thread::scope(|scope| {
            let handle = scope.spawn(move || {
                run_occupancy_sampler(
                    &stop_c,
                    &edges_c,
                    Duration::from_millis(2),
                    None,
                    Some(tele),
                );
            });
            std::thread::sleep(Duration::from_millis(30));
            stop.store(true, Ordering::Relaxed);
            handle.join().unwrap();
        });
        let summary = std::fs::read_to_string(dir.join("run.ticks.summary.tsv")).unwrap();
        assert!(summary.lines().count() >= 2, "header + >=1 data row");
        assert!(std::fs::metadata(dir.join("run.ticks.workers.tsv")).is_ok());
        // Byte-bounded edge → numeric byte deltas at the row tail (never NA here).
        let edges = std::fs::read_to_string(dir.join("run.ticks.edges.tsv")).unwrap();
        assert!(
            edges.lines().next().unwrap().ends_with("d_pushed_bytes\td_popped_bytes"),
            "edges header carries the byte columns"
        );
        // Static counter-name file names the producer's counter.
        let counter_names =
            std::fs::read_to_string(dir.join("run.ticks.counter_names.tsv")).unwrap();
        assert!(
            counter_names.contains("0\t0\trecords\trecords"),
            "counter_names names step 0's counter 0: {counter_names}"
        );
        // Per-tick counters file carries the pre-bumped cumulative value (42).
        let counters = std::fs::read_to_string(dir.join("run.ticks.counters.tsv")).unwrap();
        let last = counters.lines().last().unwrap();
        let fields: Vec<&str> = last.split('\t').collect();
        assert_eq!(fields[2], "0", "step index");
        assert_eq!(fields[3], "0", "counter index");
        assert_eq!(fields[5], "42", "cumulative value column reads the shared counter");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn timeline_tsv_marks_unsampled_edge_na() {
        // A count/unbounded edge (no depth_source) is unsampled: its depth column
        // must read `NA`, not `0.000` (which would misread as an empty byte edge).
        let m = EdgeMetrics::new();
        let edges = vec![edge_over(Arc::clone(&m), None)];
        let dir = std::env::temp_dir();
        let path = dir.join(format!("fgumi-trace-na-{}.tsv", std::process::id()));
        let stop = Arc::new(AtomicBool::new(false));
        let stop_c = Arc::clone(&stop);
        let path_c = path.clone();
        let handle = std::thread::spawn(move || {
            run_occupancy_sampler(&stop_c, &edges, Duration::from_millis(2), Some(path_c), None);
        });
        std::thread::sleep(Duration::from_millis(30));
        stop.store(true, Ordering::Relaxed);
        handle.join().unwrap();

        let content = std::fs::read_to_string(&path).expect("trace file written");
        let _ = std::fs::remove_file(&path);
        let mut lines = content.lines();
        let _header = lines.next().expect("header row");
        let row = lines.next().expect("at least one data row");
        // Columns: t_ms, <edge>.depth, <edge>.pushed, <edge>.popped.
        let depth = row.split('\t').nth(1).expect("depth column");
        assert_eq!(depth, "NA", "unsampled edge's depth column is NA, row: {row}");
    }

    #[test]
    fn terminal_finalize_captures_post_last_tick_counter_updates() {
        // Regression: work can advance `StepCounters` (and edge metrics) during
        // the sampler's final sleep, AFTER the last regular tick has emitted.
        // The shutdown path must finalize those cumulative totals
        // (`emit_terminal`) so the counters TSV's last row reflects the
        // post-update value — not the stale total from the last regular tick.
        use crate::runtime::telemetry::TelemetryWriters;
        use crate::runtime::worker_state::{WorkerState, WorkerStateBoard};
        const COUNTER_SPECS: &[&[CounterSpec]] = &[&[CounterSpec::new("records", "records")]];
        let m = EdgeMetrics::new();
        let q = Arc::new(ByteBoundedQueue::<Heavy>::new(1000));
        q.try_push(Heavy(vec![0; 300])).unwrap();
        let edges =
            vec![edge_over(Arc::clone(&m), Some(Arc::clone(&q) as Arc<dyn BoundedQueueHandle>))];
        let dir =
            std::env::temp_dir().join(format!("fgumi-sampler-terminal-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let stem = dir.join("run");
        let board = WorkerStateBoard::new(1);
        board.stamp(0, WorkerState::Running, Some(StepIdx(0)));
        // Counter starts at 42; the run bumps it by 58 → 100 during the final
        // sleep, after the first (and only) regular tick has emitted 42.
        let step_counters = vec![StepCounters::new(1)];
        step_counters[0].add(0, 42);
        let writers = TelemetryWriters::open(&stem, &["producer"], COUNTER_SPECS).unwrap();
        let tele = crate::runtime::sampler::TelemetryArgs {
            board: &board,
            n_workers: 1,
            n_pool: 1,
            step_names: &["producer"],
            rss_probe: None,
            queue_bytes_budget: Some(1000),
            writers,
            source_edge_idxs: vec![0],
            sink_edge_idxs: vec![0],
            worker_slots: vec![0],
            step_counters: &step_counters,
        };
        let stop = Arc::new(AtomicBool::new(false));
        let stop_c = Arc::clone(&stop);
        let edges_c = edges;
        std::thread::scope(|scope| {
            let handle = scope.spawn(move || {
                // A very large interval: the first tick fires immediately, then
                // the sampler parks in its stop-checked sleep — so exactly one
                // regular tick emits before the terminal finalize runs.
                run_occupancy_sampler(
                    &stop_c,
                    &edges_c,
                    Duration::from_secs(3600),
                    None,
                    Some(tele),
                );
            });
            // Let the single regular tick emit (it fires within microseconds of
            // thread start), THEN advance the shared counter during the final
            // sleep window before signalling stop.
            std::thread::sleep(Duration::from_millis(80));
            step_counters[0].add(0, 58);
            stop.store(true, Ordering::Relaxed);
            handle.join().unwrap();
        });
        // The terminal finalize must have written a final counters row carrying
        // the post-update cumulative total (100), not the last tick's 42.
        let counters = std::fs::read_to_string(dir.join("run.ticks.counters.tsv")).unwrap();
        let last = counters.lines().last().unwrap();
        let fields: Vec<&str> = last.split('\t').collect();
        assert_eq!(fields[2], "0", "step index");
        assert_eq!(fields[3], "0", "counter index");
        assert_eq!(
            fields[5], "100",
            "terminal row carries the post-final-sleep cumulative total: {counters}"
        );
        std::fs::remove_dir_all(&dir).ok();
    }
}
