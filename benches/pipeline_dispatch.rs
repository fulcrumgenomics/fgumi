//! Dispatch-throughput benchmark for the typed-step pipeline runtime.
//!
//! This measures the runtime's *per-dispatch overhead*, not the work steps do:
//! every step here is a near-no-op, so wall time is dominated by the scheduler
//! loop, queue push/pop, and the reorder stage. That is exactly the path any
//! change to the driver's inner loop lands on.
//!
//! It exists because the pipeline had no benchmark at all — no command on this
//! branch drives it yet (the chain builders and command rewiring arrive later),
//! so there was no way to tell whether a hot-path change cost anything.
//!
//! Read the numbers as *relative* only. Absolute throughput here is meaningless
//! as a product metric (real steps decompress BGZF and parse records, which
//! dwarfs dispatch); the point is A/B on the same host.

use std::collections::VecDeque;
use std::hint::black_box;
use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use criterion::{Criterion, criterion_group, criterion_main};
use fgumi_pipeline_core::Unpushed;
use fgumi_pipeline_core::builder::{Pipeline, PipelineConfig};
use fgumi_pipeline_core::held::HeldSlot;
use fgumi_pipeline_core::item::{HeapSize, Ordered};
use fgumi_pipeline_core::outputs::OrderedBytesSingle;
use fgumi_pipeline_core::queues::QueueSpec;
use fgumi_pipeline_core::reorder::BranchOrdering;
use fgumi_pipeline_core::runtime::stats::{PipelineStats, StepStatsSnapshot};
use fgumi_pipeline_core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

/// Per-edge byte budget. Large enough that the queues never bind — this
/// benchmark measures dispatch cost, and backpressure stalls would swamp it
/// with scheduling noise.
const EDGE_LIMIT_BYTES: u64 = 64 * 1024 * 1024;

/// Items pushed through the chain per iteration.
const N_ITEMS: u64 = 50_000;

/// Items per iteration for the per-item-work variants: at 150 µs per item and
/// two Parallel steps, `N_ITEMS` would cost ~15 s per iteration; this keeps a
/// `sample_size(10)` group near 15 s.
const N_ITEMS_WORK: u64 = 5_000;

/// A minimal ordered item. `heap_size` is a small constant rather than 0 so the
/// byte-bounded accounting does real work per item, as it would in production.
#[derive(Clone, Copy)]
struct Item {
    ordinal: u64,
}

impl HeapSize for Item {
    fn heap_size(&self) -> usize {
        64
    }
}

impl Ordered for Item {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

/// Emits `n` items and finishes. `Exclusive` so it owns its cursor.
struct CountingSource {
    remaining: VecDeque<Item>,
    held: HeldSlot<Unpushed<Item>>,
}

impl CountingSource {
    fn new(n: u64) -> Self {
        Self { remaining: (0..n).map(|ordinal| Item { ordinal }).collect(), held: HeldSlot::new() }
    }
}

impl Step for CountingSource {
    type Input = ();
    type Outputs = OrderedBytesSingle<Item>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "BenchSource",
            kind: StepKind::Exclusive,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: EDGE_LIMIT_BYTES }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }
        let Some(item) = self.remaining.pop_front() else {
            return Ok(StepOutcome::Finished);
        };
        if let Err(unpushed) = ctx.outputs.push(item) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

/// A pass-through. As a `Parallel` step it is cloned per worker, so at
/// `threads > 1` it is what puts several workers on the dispatch path at once;
/// as a `Detached` step it runs on its own driver thread. `work_us` busy-loops
/// that long per item (0 = no work).
struct PassThrough {
    held: HeldSlot<Unpushed<Item>>,
    kind: StepKind,
    work_us: u64,
}

impl PassThrough {
    fn new(kind: StepKind, work_us: u64) -> Self {
        Self { held: HeldSlot::new(), kind, work_us }
    }
}

impl Step for PassThrough {
    type Input = Item;
    type Outputs = OrderedBytesSingle<Item>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "BenchPassThrough",
            kind: self.kind,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: EDGE_LIMIT_BYTES }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }
        let Some(item) = ctx.input.pop() else {
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        if self.work_us > 0 {
            let start = Instant::now();
            while start.elapsed() < Duration::from_micros(self.work_us) {
                std::hint::spin_loop();
            }
        }
        // A trivial amount of real work so the compiler cannot fold the step away.
        let out = Item { ordinal: black_box(item).ordinal };
        if let Err(unpushed) = ctx.outputs.push(out) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }

    fn new_worker_copy(&self) -> Self {
        Self::new(self.kind, self.work_us)
    }
}

/// Terminal sink; counts arrivals into a shared atomic so the run is verifiable.
/// `Exclusive` in the plain chain (as the no-op baseline always was); `Serial`
/// in the detached-mid chain, whose source already takes the one exclusive
/// slot a one-worker run has.
struct CountingSink {
    seen: Arc<AtomicU64>,
    kind: StepKind,
}

impl Step for CountingSink {
    type Input = Item;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "BenchSink",
            kind: self.kind,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(item) => {
                black_box(item);
                self.seen.fetch_add(1, Ordering::Relaxed);
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// Which chain shape a group runs.
#[derive(Clone, Copy)]
enum Shape {
    /// Source → Parallel → Sink: the no-op regression check for the hot path.
    Plain,
    /// Source → Parallel → Detached → Parallel → Sink: the sort's two-sided
    /// shape, where the middle step runs on its own driver thread and every
    /// push into / out of it is a directed wake.
    DetachedMid,
}

/// Build the chain; the returned counter is the sink's.
fn build_chain(shape: Shape, work_us: u64, n_items: u64) -> (Pipeline, Arc<AtomicU64>) {
    let seen = Arc::new(AtomicU64::new(0));
    let builder = Pipeline::builder();
    match shape {
        Shape::Plain => builder
            .chain(CountingSource::new(n_items))
            .chain(PassThrough::new(StepKind::Parallel, work_us))
            .chain(CountingSink { seen: Arc::clone(&seen), kind: StepKind::Exclusive })
            .into_sink_marker(),
        Shape::DetachedMid => builder
            .chain(CountingSource::new(n_items))
            .chain(PassThrough::new(StepKind::Parallel, work_us))
            .chain(PassThrough::new(StepKind::Detached, 0))
            .chain(PassThrough::new(StepKind::Parallel, work_us))
            .chain(CountingSink { seen: Arc::clone(&seen), kind: StepKind::Serial })
            .into_sink_marker(),
    }
    (builder.build().expect("bench chain builds"), seen)
}

/// Run a built chain, asserting every item arrived; with `stats`, the stats
/// handle must come from this pipeline.
///
/// The assert is not decoration: a timing harness reports a wall figure for a
/// run that failed or short-circuited just as happily as for a correct one, and
/// a chain that drops items would look like a speedup.
fn run_built(
    pipeline: Pipeline,
    seen: &AtomicU64,
    n_items: u64,
    threads: usize,
    stats: Option<Arc<PipelineStats>>,
) {
    let mut config = PipelineConfig { threads, ..Default::default() };
    if let Some(s) = stats {
        config = config.with_stats(s);
    }
    pipeline.run(config).expect("bench chain runs");
    assert_eq!(
        seen.load(Ordering::Relaxed),
        n_items,
        "every item must reach the sink — a short run is an invalid measurement",
    );
}

/// One instrumented run per variant, printed before the timing groups so the
/// wake counters sit next to the wall numbers in the same log.
fn print_counters() {
    for (shape, label) in [(Shape::Plain, "plain"), (Shape::DetachedMid, "detached_mid")] {
        for threads in [1usize, 4, 8, 16] {
            for work_us in [0u64, 50, 150] {
                let n = if work_us == 0 { N_ITEMS } else { N_ITEMS_WORK };
                let (pipeline, seen) = build_chain(shape, work_us, n);
                let stats = pipeline.stats();
                run_built(pipeline, &seen, n, threads, Some(Arc::clone(&stats)));
                let snap = stats.snapshot();
                let sum = |f: fn(&StepStatsSnapshot) -> u64| {
                    snap.steps.iter().map(|(_, s)| f(s)).sum::<u64>()
                };
                let timed_out: u64 = snap.worker_waits.iter().map(|w| w.3).sum();
                eprintln!(
                    "counters {label} threads={threads} work_us={work_us}: notifies_issued={} \
                     unparks_issued={} reverse_wakes={} direct_fallbacks={} \
                     ec_waits_timed_out={timed_out}",
                    sum(|s| s.notifies_issued),
                    sum(|s| s.unparks_issued),
                    sum(|s| s.reverse_wakes),
                    sum(|s| s.direct_fallbacks),
                );
            }
        }
    }
}

fn bench_dispatch(c: &mut Criterion) {
    print_counters();
    let mut group = c.benchmark_group("pipeline_dispatch");
    group.throughput(criterion::Throughput::Elements(N_ITEMS));
    // The no-op groups are the regression check for the hot path.
    for threads in [1usize, 4, 8] {
        group.bench_function(format!("threads_{threads}"), |b| {
            b.iter(|| {
                let (p, seen) = build_chain(Shape::Plain, 0, N_ITEMS);
                run_built(p, &seen, N_ITEMS, threads, None);
            });
        });
    }
    group.finish();
    let mut group = c.benchmark_group("pipeline_dispatch_detached");
    group.sample_size(10);
    group.throughput(criterion::Throughput::Elements(N_ITEMS_WORK));
    for threads in [1usize, 4, 8, 16] {
        for work_us in [50u64, 150] {
            group.bench_function(format!("detached_mid_threads_{threads}_work_{work_us}us"), |b| {
                b.iter(|| {
                    let (p, seen) = build_chain(Shape::DetachedMid, work_us, N_ITEMS_WORK);
                    run_built(p, &seen, N_ITEMS_WORK, threads, None);
                });
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_dispatch);
criterion_main!(benches);
