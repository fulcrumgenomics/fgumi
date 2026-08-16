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

use criterion::{Criterion, criterion_group, criterion_main};
use fgumi_pipeline_core::Unpushed;
use fgumi_pipeline_core::builder::{Pipeline, PipelineConfig};
use fgumi_pipeline_core::held::HeldSlot;
use fgumi_pipeline_core::item::{HeapSize, Ordered};
use fgumi_pipeline_core::outputs::OrderedBytesSingle;
use fgumi_pipeline_core::queues::QueueSpec;
use fgumi_pipeline_core::reorder::BranchOrdering;
use fgumi_pipeline_core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

/// Per-edge byte budget. Large enough that the queues never bind — this
/// benchmark measures dispatch cost, and backpressure stalls would swamp it
/// with scheduling noise.
const EDGE_LIMIT_BYTES: u64 = 64 * 1024 * 1024;

/// Items pushed through the chain per iteration.
const N_ITEMS: u64 = 50_000;

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

/// A `Parallel` pass-through. Cloned per worker, so at `threads > 1` this is
/// what puts several workers on the dispatch path at once.
struct PassThrough {
    held: HeldSlot<Unpushed<Item>>,
}

impl Step for PassThrough {
    type Input = Item;
    type Outputs = OrderedBytesSingle<Item>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "BenchPassThrough",
            kind: StepKind::Parallel,
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
        // A trivial amount of real work so the compiler cannot fold the step away.
        let out = Item { ordinal: black_box(item).ordinal };
        if let Err(unpushed) = ctx.outputs.push(out) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }

    fn new_worker_copy(&self) -> Self {
        Self { held: HeldSlot::new() }
    }
}

/// Terminal sink; counts arrivals into a shared atomic so the run is verifiable.
struct CountingSink {
    seen: Arc<AtomicU64>,
}

impl Step for CountingSink {
    type Input = Item;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "BenchSink",
            kind: StepKind::Exclusive,
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

/// Run the chain once, asserting every item arrived.
///
/// The assert is not decoration: a timing harness reports a wall figure for a
/// run that failed or short-circuited just as happily as for a correct one, and
/// a chain that drops items would look like a speedup.
fn run_chain(threads: usize) {
    let seen = Arc::new(AtomicU64::new(0));
    let sink_handle = Arc::clone(&seen);

    let builder = Pipeline::builder();
    builder
        .chain(CountingSource::new(N_ITEMS))
        .chain(PassThrough { held: HeldSlot::new() })
        .chain(CountingSink { seen: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("bench chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("bench chain runs");

    assert_eq!(
        seen.load(Ordering::Relaxed),
        N_ITEMS,
        "every item must reach the sink — a short run is an invalid measurement",
    );
}

fn bench_dispatch(c: &mut Criterion) {
    let mut group = c.benchmark_group("pipeline_dispatch");
    group.throughput(criterion::Throughput::Elements(N_ITEMS));
    for threads in [1usize, 4, 8] {
        group.bench_function(format!("threads_{threads}"), |b| {
            b.iter(|| run_chain(threads));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_dispatch);
criterion_main!(benches);
