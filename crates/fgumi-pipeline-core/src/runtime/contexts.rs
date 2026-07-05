//! `ChainContexts`: per-step typed input + outputs handles, constructed
//! at `Pipeline::run` start by walking `ChainGraph` and calling each
//! step's `build_input_handle` / `build_output_set`.
//!
//! Layout:
//!   - `inputs[step_idx]` — `Box<dyn Any>` carrying `BranchInputHandle<S::Input>`
//!     (or a dummy unit handle for sources). The worker loop hands this to
//!     `ErasedStepCtx.input`.
//!   - `outputs[step_idx]` — `Box<dyn Any>` carrying `OutputHandles<S::Outputs>`.
//!     The worker loop hands this to `ErasedStepCtx.outputs` AND uses it
//!     for the typed `mark_outputs_drained` dispatch via the step's
//!     `ErasedStep::mark_outputs_drained` method.
//!
//! Since `OutputsViewAny` no longer carries per-branch drained-flag Arcs
//! (option (c) pivot — drain marking is per-branch via
//! `BranchOutputHandle::mark_drained`, dispatched through the typed view),
//! the framework's drain-propagation path goes through
//! `ErasedStep::mark_outputs_drained(outputs_any.as_ref())`.

use std::any::Any;
use std::sync::Arc;

use crate::erased::ErasedStep;
use crate::handles::{BranchInputHandle, OutputQueueSet};
use crate::item::HeapSize;
use crate::topology::{BranchIdx, ChainGraph, StepIdx};

/// Per-step input + outputs handles. The worker loop hands these to
/// `ErasedStepCtx` for each `try_run_erased` dispatch and to
/// `ErasedStep::mark_outputs_drained` for the drain-propagation path.
pub struct ChainContexts {
    /// `inputs[step_idx]` = boxed `BranchInputHandle<S::Input>` (or a dummy
    /// drained unit handle for sources).
    pub inputs: Vec<Box<dyn Any + Send + Sync>>,
    /// `outputs[step_idx]` = boxed `OutputHandles<S::Outputs>` (the typed
    /// surface step authors push into; also dispatchable via
    /// `ErasedStep::mark_outputs_drained`).
    pub outputs: Vec<Box<dyn Any + Send + Sync>>,
    /// Registry of every byte-bounded queue in the chain. Populated
    /// during `build_chain_contexts` from each `Branch`'s
    /// `bounded_queue_handle`. Empty when no byte-bounded queues
    /// exist (e.g. a chain composed entirely of `CountBounded` /
    /// `Unbounded` steps).
    pub bounded_queues: Vec<RegisteredQueue>,
    /// Registry of every instrumented edge (`--pipeline-trace`), with its
    /// `EdgeMetrics` + producer/consumer + (for byte-bounded edges) a depth
    /// source. Empty when instrumentation is `Off`. Read by the occupancy
    /// sampler and the end-of-run edge report.
    pub edges: Vec<RegisteredEdge>,
}

/// One instrumented edge: its shared [`EdgeMetrics`](crate::runtime::metrics::EdgeMetrics)
/// (push counts from the transport, pop counts from the input handle), its
/// producer + consumer steps, and a depth source for occupancy sampling.
pub struct RegisteredEdge {
    pub producer_step: StepIdx,
    pub producer_name: &'static str,
    /// `None` for a terminal branch with no consumer (e.g. a `--rejects` tail).
    pub consumer_step: Option<StepIdx>,
    pub consumer_name: Option<&'static str>,
    pub branch: BranchIdx,
    /// On an **ordered** edge, the push-side counters (`pushed_items` /
    /// `pushed_bytes` / `push_rejections`) are recorded at the `ReorderStage`
    /// boundary, not on the transport queue: the reorder stage's must-accept path
    /// turns a full-transport reject into an accepted (stashed) push, so recording
    /// on the transport would miscount stashed items as backpressure. Recorded at
    /// the boundary, `pushed_bytes` counts the bare `T` (`heap_size`) and so
    /// matches `popped_bytes` (the ordinal tag adds no heap).
    pub metrics: std::sync::Arc<crate::runtime::metrics::EdgeMetrics>,
    /// `Some` for byte-bounded edges — the transport queue's `current_bytes`
    /// (plus, for an ordered edge, the `reorder_depth` stash bytes) over
    /// `limit_bytes` gives the occupancy fraction the sampler records. `None`
    /// for count/unbounded edges (counters still apply; occupancy histogram is
    /// skipped).
    pub depth_source: Option<std::sync::Arc<dyn crate::queues::BoundedQueueHandle>>,
    /// `Some` for an **ordered** byte-bounded edge — the `ReorderStage`'s
    /// overflow-stash handle. Its buffered bytes are added to the transport's
    /// `current_bytes` when sampling depth, so an ordered edge reflects total
    /// buffered bytes (transport + reorder stash) rather than the transport
    /// queue alone (which under-reports occupancy when items pile in the reorder
    /// buffer under producer skew). `None` for direct/count/unbounded edges.
    pub reorder_depth: Option<std::sync::Arc<dyn crate::reorder::ReorderCapHandle>>,
}

/// One byte-bounded queue's location in the chain plus the handles for
/// updating its budget: the transport-limit setter (`handle`) and, for an
/// ordered branch, the reorder stage's overflow-cap setter (`reorder_cap`).
/// Used by the budget pass (which sets both from one per-edge budget) and
/// the queue-memory rebalancer (which shifts transport budget at runtime).
pub struct RegisteredQueue {
    pub producer_step_name: &'static str,
    pub producer_step: StepIdx,
    pub branch: BranchIdx,
    pub handle: std::sync::Arc<dyn crate::queues::BoundedQueueHandle>,
    /// `Some` for an ordered byte-bounded branch — its reorder overflow stash
    /// is sized from the same per-edge budget as `handle`. `None` for a direct
    /// (unordered) byte-bounded branch (no reorder stage).
    pub reorder_cap: Option<std::sync::Arc<dyn crate::reorder::ReorderCapHandle>>,
}

/// Build the chain's per-step contexts.
///
/// Walks `graph` to discover each step's producer and pulls the typed input
/// handle out of the producer's `OutputQueueSet`. Sources get a dummy
/// pre-drained `BranchInputHandle<()>` (their input is implicitly drained
/// from the start; framework never pops from it).
///
/// # Panics
///
/// Panics if `graph.n_steps() != steps.len()` or if any non-source step
/// has no producer (the builder's all-wired check should prevent this).
#[must_use]
pub fn build_chain_contexts(
    steps: &[Box<dyn ErasedStep>],
    graph: &ChainGraph,
    level: crate::builder::InstrumentationLevel,
) -> ChainContexts {
    build_chain_contexts_inner(steps, graph, false, level)
}

/// Build the chain's per-step contexts with **direct, unbounded** inter-step
/// transports (no byte bounds, no reorder stages) — the wiring the single-
/// thread *fused* driver runs over.
///
/// Identical to [`build_chain_contexts`] except every producer's output set is
/// built via [`ErasedStep::build_fused_output_set`] instead of
/// [`ErasedStep::build_output_set`]. At one worker FIFO push order is already
/// the correct order, so dropping the reorder stage is sound, and the unbounded
/// buffers stay shallow because the consumer is driven in the same pass. Only
/// call this for a linear (single source → single sink, single-input) chain at
/// `--threads 1`; see [`super::run_fused_single_thread`].
#[must_use]
pub fn build_chain_contexts_fused(
    steps: &[Box<dyn ErasedStep>],
    graph: &ChainGraph,
) -> ChainContexts {
    // The fused path uses direct, unbounded transports (no bounded edges), so
    // there is nothing to instrument — force `Off` regardless of the run's level.
    build_chain_contexts_inner(steps, graph, true, crate::builder::InstrumentationLevel::Off)
}

/// Build the producer→consumer edge map (keyed by `(producer_step, branch)`)
/// used to label registered edges with their consumer. Empty when
/// instrumentation is off (no edges are registered, so the map is unused).
fn build_consumer_map(
    steps: &[Box<dyn ErasedStep>],
    graph: &ChainGraph,
    level: crate::builder::InstrumentationLevel,
) -> std::collections::HashMap<(usize, usize), usize> {
    let mut m = std::collections::HashMap::new();
    if !level.is_on() {
        return m;
    }
    for (consumer_idx, step) in steps.iter().enumerate() {
        if step.is_source() {
            continue;
        }
        match step.input_arity() {
            1 => {
                let (p, b) = find_producer(graph, StepIdx(consumer_idx));
                m.insert((p.0, b.0), consumer_idx);
            }
            2 => {
                for (p, b) in find_all_producers(graph, StepIdx(consumer_idx)) {
                    m.insert((p.0, b.0), consumer_idx);
                }
            }
            _ => {}
        }
    }
    m
}

fn build_chain_contexts_inner(
    steps: &[Box<dyn ErasedStep>],
    graph: &ChainGraph,
    direct: bool,
    level: crate::builder::InstrumentationLevel,
) -> ChainContexts {
    assert_eq!(steps.len(), graph.n_steps(), "chain length mismatch");

    let n_steps = steps.len();
    let mut output_sets: Vec<OutputQueueSet> = Vec::with_capacity(n_steps);
    let mut outputs: Vec<Box<dyn Any + Send + Sync>> = Vec::with_capacity(n_steps);
    let mut inputs: Vec<Box<dyn Any + Send + Sync>> = Vec::with_capacity(n_steps);

    // Pass 1: build each step's output set + boxed `OutputHandles<O>`. The
    // fused driver uses direct, unbounded transports (`build_fused_output_set`);
    // the scheduled path honours each step's profiled queue spec + ordering.
    for step in steps {
        let (set, view) =
            if direct { step.build_fused_output_set(level) } else { step.build_output_set(level) };
        let outputs_box = step.wrap_outputs_view(view);
        output_sets.push(set);
        outputs.push(outputs_box);
    }

    // Producer→consumer map for edge registration (inverse of `find_producer`,
    // which Pass 2 uses consumer→producer). Empty when instrumentation is off.
    let consumer_of = build_consumer_map(steps, graph, level);

    // Pass 1.5: collect byte-bounded queue handles into the registry.
    // Must run before Pass 2 because `take_typed_input` (called via
    // `build_input_handle`) replaces the `BranchEntry` with a fresh
    // one whose `bounded_queue_handle` is `None` — by then the
    // handles have been moved out of the chain.
    let mut bounded_queues: Vec<RegisteredQueue> = Vec::new();
    let mut edges: Vec<RegisteredEdge> = Vec::new();
    for (step_idx_usize, set) in output_sets.iter().enumerate() {
        for (branch_idx_usize, entry) in set.branches.iter().enumerate() {
            if let Some(handles) = entry.bounded_queue_handle.as_ref() {
                bounded_queues.push(RegisteredQueue {
                    producer_step_name: steps[step_idx_usize].profile().name,
                    producer_step: StepIdx(step_idx_usize),
                    branch: BranchIdx(branch_idx_usize),
                    handle: std::sync::Arc::clone(&handles.transport),
                    reorder_cap: handles.reorder_cap.clone(),
                });
            }
            // Instrumented edge: register its shared metrics + producer/consumer
            // + (for byte-bounded edges) a depth source for occupancy sampling.
            if let Some(metrics) = entry.metrics.as_ref() {
                let consumer_step =
                    consumer_of.get(&(step_idx_usize, branch_idx_usize)).map(|c| StepIdx(*c));
                edges.push(RegisteredEdge {
                    producer_step: StepIdx(step_idx_usize),
                    producer_name: steps[step_idx_usize].profile().name,
                    consumer_step,
                    consumer_name: consumer_step.map(|c| steps[c.0].profile().name),
                    branch: BranchIdx(branch_idx_usize),
                    metrics: std::sync::Arc::clone(metrics),
                    depth_source: entry
                        .bounded_queue_handle
                        .as_ref()
                        .map(|h| std::sync::Arc::clone(&h.transport)),
                    reorder_depth: entry
                        .bounded_queue_handle
                        .as_ref()
                        .and_then(|h| h.reorder_cap.clone()),
                });
            }
        }
    }

    // Pass 2: input handles. Sources get a dummy drained unit handle;
    // single-input mid-steps and sinks pull one handle from their
    // producer's `OutputQueueSet`; two-input merge steps (`Step2`
    // adapters, identified by `input_arity() == 2`) pull TWO handles
    // — one per consumer-input-slot — and wrap them in a
    // `TwoInputHandles<A, B>`.
    for (consumer_idx, step) in steps.iter().enumerate() {
        let input_box = if step.is_source() {
            // Source step (Input = ()). Build a permanently-drained
            // dummy unit handle — the worker loop never pops from it.
            // Chains with multiple sources (e.g. zipper's mapped +
            // unmapped subchains converging at a Step2 merger) all
            // share this code path; the runtime walks every source
            // independently to Finished.
            dummy_unit_input_handle()
        } else {
            match step.input_arity() {
                1 => {
                    let (producer_idx, branch_idx) = find_producer(graph, StepIdx(consumer_idx));
                    step.build_input_handle(&mut output_sets[producer_idx.0], branch_idx.0)
                }
                2 => {
                    let edges = find_all_producers(graph, StepIdx(consumer_idx));
                    assert_eq!(
                        edges.len(),
                        2,
                        "Step2 consumer {:?} expects 2 input edges, found {}",
                        StepIdx(consumer_idx),
                        edges.len()
                    );
                    let (p0, p0_branch) = edges[0];
                    let (p1, p1_branch) = edges[1];
                    step.build_two_input_handles(
                        &mut output_sets,
                        p0.0,
                        p0_branch.0,
                        p1.0,
                        p1_branch.0,
                    )
                }
                n => panic!("unsupported input_arity {n} for step {:?}", StepIdx(consumer_idx)),
            }
        };
        inputs.push(input_box);
    }

    // `output_sets` is consumed implicitly here — every branch was taken
    // exactly once via `build_input_handle`, leaving placeholder slots.
    drop(output_sets);

    ChainContexts { inputs, outputs, bounded_queues, edges }
}

/// Construct a `BranchInputHandle<()>` that's already drained — for source
/// steps whose input is implicitly empty + drained from t=0. Uses the zero-state
/// `always_drained` handle (no backing queue), so building a source costs no
/// `SegQueue` allocation for a handle whose only job is to report
/// `is_drained() == true` (the worker loop never pops from a source's input).
fn dummy_unit_input_handle() -> Box<dyn Any + Send + Sync> {
    Box::new(BranchInputHandle::<()>::always_drained())
}

/// Find the (producer, branch) that produces the given consumer step.
///
/// # Panics
///
/// Panics if no producer exists (the builder's all-wired check should
/// have prevented this on a built pipeline).
fn find_producer(graph: &ChainGraph, consumer: StepIdx) -> (StepIdx, BranchIdx) {
    for producer_usize in 0..graph.n_steps() {
        let producer = StepIdx(producer_usize);
        let n_branches = graph.branch_count(producer);
        for branch_usize in 0..n_branches {
            let branch = BranchIdx(branch_usize);
            if graph.consumer(producer, branch) == Some(consumer) {
                return (producer, branch);
            }
        }
    }
    panic!(
        "step {consumer:?} has no producer in chain graph; \
         the builder's all-wired check should have caught this"
    );
}

/// Find all `(producer, producer_branch)` edges feeding the given
/// multi-input consumer, sorted by the consumer's input-slot index
/// (`returned[i]` feeds the consumer's input slot `i`). Used by
/// [`build_chain_contexts`] to assemble per-branch input handles
/// for `Step2` and future `StepN` consumers.
///
/// # Panics
///
/// Panics if any of the consumer's input slots is unwired (the
/// builder's all-wired check should have caught this on a built
/// pipeline).
fn find_all_producers(graph: &ChainGraph, consumer: StepIdx) -> Vec<(StepIdx, BranchIdx)> {
    let arity = graph.input_arity(consumer);
    let mut edges: Vec<Option<(StepIdx, BranchIdx)>> = vec![None; arity];
    for producer_usize in 0..graph.n_steps() {
        let producer = StepIdx(producer_usize);
        let n_branches = graph.branch_count(producer);
        for branch_usize in 0..n_branches {
            let branch = BranchIdx(branch_usize);
            if graph.consumer(producer, branch) == Some(consumer) {
                let slot = graph
                    .consumer_input_slot(producer, branch)
                    .expect("consumer_input_slot must be set when consumer is wired");
                assert!(
                    slot < arity,
                    "consumer {consumer:?} edge has input slot {slot} but arity is {arity}"
                );
                assert!(
                    edges[slot].is_none(),
                    "consumer {consumer:?} input slot {slot} wired more than once",
                );
                edges[slot] = Some((producer, branch));
            }
        }
    }
    edges
        .into_iter()
        .enumerate()
        .map(|(slot, e)| {
            e.unwrap_or_else(|| panic!("consumer {consumer:?} input slot {slot} has no producer"))
        })
        .collect()
}

/// Convenience: borrow the typed `BranchInputHandle<T>` for a given step.
///
/// # Panics
///
/// Panics if the step's input handle doesn't downcast to `T` (a framework
/// invariant violation).
#[must_use]
pub fn input_as<T: Send + HeapSize + 'static>(
    contexts: &Arc<ChainContexts>,
    step: StepIdx,
) -> &BranchInputHandle<T> {
    contexts.inputs[step.0]
        .downcast_ref::<BranchInputHandle<T>>()
        .expect("BranchInputHandle<T> downcast failed in input_as")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    use proptest::prelude::*;
    use rstest::rstest;

    use crate::erased::TypedStep;
    use crate::outputs::Single;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;
    use crate::step::{InputHandle, Step, StepCtx, StepKind, StepOutcome, StepProfile};

    #[derive(Clone)]
    struct StubSource;
    impl Step for StubSource {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Source",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::ByOrdinal],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::Finished)
        }
    }

    #[derive(Clone)]
    struct StubSink;
    impl Step for StubSink {
        type Input = u32;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Sink",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    /// A `u32 → u32` pass-through step used to grow a linear chain to an
    /// arbitrary length between the source and sink.
    #[derive(Clone)]
    struct MiddleStep;
    impl Step for MiddleStep {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Middle",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    /// Build an `n`-step linear chain `Source → Middle×(n - 2) → Sink`,
    /// returning the erased step boxes alongside the fully wired graph.
    ///
    /// # Panics
    ///
    /// Panics if `n < 2` (a linear chain needs at least a source and a sink).
    fn linear_chain(n: usize) -> (Vec<Box<dyn ErasedStep>>, ChainGraph) {
        assert!(n >= 2, "linear chain needs at least a source and a sink");

        let mut graph = ChainGraph::new();
        let mut steps: Vec<Box<dyn ErasedStep>> = Vec::with_capacity(n);
        let mut indices: Vec<StepIdx> = Vec::with_capacity(n);

        indices.push(graph.register_step("Source", 1));
        steps.push(Box::new(TypedStep::new(StubSource)));
        for _ in 1..n - 1 {
            indices.push(graph.register_step("Middle", 1));
            steps.push(Box::new(TypedStep::new(MiddleStep)));
        }
        indices.push(graph.register_step("Sink", 0));
        steps.push(Box::new(TypedStep::new(StubSink)));

        for pair in indices.windows(2) {
            graph.wire(pair[0], BranchIdx(0), pair[1]);
        }

        (steps, graph)
    }

    /// Assert the invariants `build_chain_contexts` must uphold for an
    /// `n`-step linear chain: one input/output slot per step, no byte-bounded
    /// queues (the stubs only use `CountBounded` transport), the source's
    /// input is a dummy pre-drained `BranchInputHandle<()>`, and every
    /// downstream step receives a real `BranchInputHandle<u32>` wired from its
    /// producer.
    fn assert_linear_chain_invariants(ctx: &ChainContexts, n: usize) {
        assert_eq!(ctx.inputs.len(), n);
        assert_eq!(ctx.outputs.len(), n);
        assert!(ctx.bounded_queues.is_empty());

        let source = ctx.inputs[0]
            .downcast_ref::<BranchInputHandle<()>>()
            .expect("source must have a dummy unit input handle");
        // Pin the `always_drained()` invariant: a source's input is implicitly
        // drained from t=0, so a regression that swapped it back to a real
        // (never-draining) queue would be caught here, not just the downcast.
        assert!(
            source.is_drained(),
            "source's dummy unit input handle must report drained (always_drained invariant)"
        );
        for i in 1..n {
            assert!(
                ctx.inputs[i].downcast_ref::<BranchInputHandle<u32>>().is_some(),
                "downstream step {i} must have a BranchInputHandle<u32> wired from its producer"
            );
        }
    }

    /// `build_chain_contexts` wires linear chains of varying length: the
    /// three-step case also exercises `find_all_producers` for a middle step.
    #[rstest]
    #[case(2)]
    #[case(3)]
    #[case(4)]
    fn build_chain_contexts_linear(#[case] n: usize) {
        let (steps, graph) = linear_chain(n);
        let ctx = build_chain_contexts(&steps, &graph, crate::builder::InstrumentationLevel::Off);
        assert_linear_chain_invariants(&ctx, n);
    }

    #[test]
    fn registry_covers_edges_with_producer_and_consumer() {
        use crate::builder::InstrumentationLevel as L;
        // Source → Middle → Sink: two producing edges, each with a consumer.
        let (steps, graph) = linear_chain(3);
        let ctx = build_chain_contexts(&steps, &graph, L::Summary);
        assert_eq!(ctx.edges.len(), 2, "every bounded producing edge is registered");
        for e in &ctx.edges {
            assert!(
                e.consumer_step.is_some(),
                "edge from {} has a resolved consumer",
                e.producer_name
            );
            assert!(e.consumer_name.is_some());
            // CountBounded edges expose no byte depth source (occupancy via len only).
            assert!(e.depth_source.is_none(), "count-bounded edge has no byte depth source");
        }
        // Off → no edges registered (zero-overhead path).
        let (steps, graph) = linear_chain(3);
        let off = build_chain_contexts(&steps, &graph, L::Off);
        assert!(off.edges.is_empty(), "level Off registers no edges");
    }

    proptest! {
        /// The typed-handle and bounded-queue invariants hold for linear
        /// chains of any length, not just the hand-picked rstest cases.
        #[test]
        fn build_chain_contexts_linear_invariants(n in 2usize..=8) {
            let (steps, graph) = linear_chain(n);
            let ctx = build_chain_contexts(&steps, &graph, crate::builder::InstrumentationLevel::Off);
            assert_linear_chain_invariants(&ctx, n);
        }
    }
}
