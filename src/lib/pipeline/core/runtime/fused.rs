//! Single-thread *fused* execution mode (issue #330).
//!
//! At `--threads 1` a linear `source → … → sink` chain gains nothing from the
//! scheduled worker pool: there is one worker, so the inter-step bounded
//! queues, round-robin polling, held-slot retries, and reorder bookkeeping are
//! pure overhead (profiling showed ~2/3 of `try_run` calls do no useful work).
//!
//! This module is the structural fix. Fusion is **not** a per-command rewrite
//! — it is an execution mode of the existing pipeline. [`is_linear_chain`]
//! detects a linear chain; [`run_fused_single_thread`] then drives the chain's
//! own type-erased steps inline, in topological order, over **direct,
//! unbounded** buffers (built by [`build_chain_contexts_fused`]). FIFO push
//! order is already the correct order at one worker, so the reorder stage is
//! dropped; each step emits a bounded number of items per `try_run`, so the
//! unbounded buffers stay shallow (bounded-fan-out invariant). The runtime's
//! generic wiring (`build_chain_contexts`) and dispatch (`try_run_erased`) are
//! reused verbatim — no step logic is duplicated.
//!
//! Non-linear chains (two-input `Step2` merges like zipper, multi-output splits
//! like `correct --rejects`) and `--threads ≥ 2` are not eligible and fall back
//! to the scheduled [`run_worker_loop`](super::driver::run_worker_loop).

use std::sync::Arc;
use std::time::Instant;

use super::contexts::build_chain_contexts_fused;
use super::stats::PipelineStats;
use crate::pipeline::core::erased::{ErasedStep, ErasedStepCtx};
use crate::pipeline::core::signal::{PipelineError, PipelineSignal};
use crate::pipeline::core::step::StepOutcome;
use crate::pipeline::core::topology::{BranchIdx, ChainGraph, StepIdx};

/// Returns `true` iff `steps` (in chain-construction order) form a chain the
/// fused driver can run on one worker: at least two steps, exactly one source
/// (at index 0), no two-input (`Step2`) merge, and every output branch wired
/// **forward** to a later step.
///
/// "Forward-wired" allows fan-out — a step may have more than one output branch
/// (e.g. the kept/rejects split of `--rejects`) as long as each branch feeds a
/// later step. Build order equals a topological order (`append_source` then
/// `append_step`/`append_step2`, producers always before consumers), so the
/// driver can walk `steps` by index in a single pass and every producer runs
/// before its consumer. A chain may therefore have **multiple sinks** (one per
/// fan-out leaf), each a zero-branch step at some index; the only structural
/// requirement is that the last step has no unwired-forward branch, which forces
/// at least one terminal sink.
///
/// Excluded: single-step chains (`n < 2` — nothing to fuse; the scheduled
/// single-thread path is already optimal), and any chain with a `Step2` merge
/// (two input streams — a single-worker inline drive assumes one source). Those
/// fall back to the scheduled worker pool.
#[must_use]
pub fn is_fusible_chain(steps: &[Box<dyn ErasedStep>], graph: &ChainGraph) -> bool {
    let n = steps.len();
    if n < 2 {
        return false;
    }
    // Exactly one source, and it must be the first step.
    if !steps[0].is_source() || steps[1..].iter().any(|s| s.is_source()) {
        return false;
    }
    for i in 0..n {
        let idx = StepIdx(i);
        // No two-input (`Step2`) consumers: a single-worker inline drive follows
        // one input stream. Sources are arity 1 in the graph; `Step2` is arity 2.
        if graph.input_arity(idx) != 1 {
            return false;
        }
        // Every output branch (one for a linear step, ≥2 for a fan-out like the
        // `--rejects` split) must be wired to a strictly later step. A sink has
        // zero branches, so its loop body is skipped; the last step necessarily
        // has no forward target, so it must be a sink.
        for b in 0..graph.branch_count(idx) {
            match graph.consumer(idx, BranchIdx(b)) {
                Some(StepIdx(j)) if j > i => {}
                _ => return false,
            }
        }
    }
    true
}

/// Drive a fusible chain to completion on the calling thread, fused.
///
/// Builds the chain's per-step contexts with direct, unbounded inter-step
/// transports, then repeatedly walks the steps in topological order — popping
/// from each step's input, pushing to its output(s) — until **every** step has
/// reported [`StepOutcome::Finished`]. On a step's `Finished` the driver marks
/// all its output branches drained so downstream steps see their inputs closed
/// (the same drain propagation the scheduled driver performs, minus the
/// `Parallel` counter gate — there is exactly one instance of each step at one
/// worker). Waiting for *all* steps (not just the last) is what lets a fan-out
/// chain (e.g. the `--rejects` split) finish both its sink subchains.
///
/// Callers must have confirmed [`is_fusible_chain`] first.
///
/// # Errors
///
/// Returns [`PipelineError::Io`] if any step's `try_run` returned `Err` (the
/// first such error wins, carrying the originating step's name), or
/// [`PipelineError::Cancelled`] if the run was cancelled via the pipeline's
/// [`CancelHandle`](crate::pipeline::core::signal::CancelHandle) — matching
/// [`crate::pipeline::core::builder::Pipeline::run`]'s contract.
pub fn run_fused_single_thread(
    mut steps: Vec<Box<dyn ErasedStep>>,
    graph: &ChainGraph,
    signal: &Arc<PipelineSignal>,
    stats: Option<&Arc<PipelineStats>>,
) -> Result<(), PipelineError> {
    let n = steps.len();
    let contexts = build_chain_contexts_fused(&steps, graph);
    let mut finished = vec![false; n];

    'drive: loop {
        // Bail promptly on an external cancel or a prior-pass error.
        if signal.is_done() {
            break;
        }
        let mut progressed = false;
        for i in 0..n {
            if finished[i] {
                continue;
            }
            let outputs_any = contexts.outputs[i].as_ref();
            let mut ctx =
                ErasedStepCtx { input: contexts.inputs[i].as_ref(), outputs: outputs_any, signal };

            // Time the dispatch only when stats collection is on (mirrors
            // `dispatch_one_step`): `Instant::now()` is non-trivial on the hot
            // path, so gate it on `stats.is_some()`.
            let start = stats.map(|_| Instant::now());
            let result = steps[i].try_run_erased(&mut ctx);
            if let (Some(stats), Some(start)) = (stats, start) {
                let elapsed_ns = u64::try_from(start.elapsed().as_nanos()).unwrap_or(u64::MAX);
                let start_ns = stats.elapsed_ns().saturating_sub(elapsed_ns);
                match &result {
                    Ok(outcome) => stats.record(StepIdx(i), *outcome, start_ns, elapsed_ns),
                    Err(_) => stats.record_error(StepIdx(i), start_ns, elapsed_ns),
                }
            }

            match result {
                Ok(StepOutcome::Progress) => progressed = true,
                Ok(StepOutcome::Finished) => {
                    // Close this step's output branches so downstream drains.
                    // No `StepDrainCounter` gate: one instance per step here.
                    steps[i].mark_outputs_drained(outputs_any);
                    finished[i] = true;
                    progressed = true;
                }
                Ok(StepOutcome::NoProgress | StepOutcome::Contention) => {}
                Err(io_err) => {
                    signal
                        .record_error(PipelineError::Io { step: steps[i].name(), source: io_err });
                    break 'drive;
                }
            }
        }
        // Done when every step has finished. For a fan-out chain (`--rejects`)
        // that means BOTH sink subchains have drained — checking only the last
        // step would break before an earlier-indexed reject sink had flushed.
        if finished.iter().all(|&f| f) {
            break;
        }
        if !progressed {
            // No step advanced and not all are finished. For a well-formed
            // fusible chain this is unreachable: the source either progresses or
            // finishes, and `Finished` cascades drain downstream so some step
            // always advances until every sink closes. Guard against an infinite
            // spin rather than trusting that invariant blindly.
            debug_assert!(
                false,
                "fused single-thread driver stalled: no step progressed and not \
                 all steps have finished"
            );
            // In release builds the `debug_assert!` is compiled out, so record
            // the stall as an error before breaking. Otherwise the loop would
            // exit and map to `Ok`, silently truncating the output instead of
            // surfacing the broken invariant.
            signal.record_error(PipelineError::Io {
                step: "fused-driver",
                source: std::io::Error::other(
                    "fused single-thread driver stalled: no step progressed and \
                     not all steps have finished",
                ),
            });
            break 'drive;
        }
    }

    // Map the recorded outcome to the run result (same shape as `Pipeline::run`;
    // `PipelineError` is not `Clone`, so reconstruct via `PipelineError::reconstruct`).
    match signal.outcome() {
        Some(err) => Err(err.reconstruct()),
        None => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use std::io;
    use std::sync::{Arc, Mutex};

    use super::*;
    use crate::pipeline::core::erased::TypedStep;
    use crate::pipeline::core::outputs::Single;
    use crate::pipeline::core::queues::QueueSpec;
    use crate::pipeline::core::reorder::BranchOrdering;
    use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepProfile};

    // ── Stub steps: source → +100 mid → collecting sink ──────────────────

    /// Source emitting `0, 1, …, count-1` then `Finished`.
    struct CountSource {
        next: u32,
        count: u32,
    }
    impl Step for CountSource {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "CountSource",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::Unbounded],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if self.next >= self.count {
                return Ok(StepOutcome::Finished);
            }
            let _ = ctx.outputs.push(self.next);
            self.next += 1;
            Ok(StepOutcome::Progress)
        }
    }

    /// Mid: pops a `u32`, pushes `+100`; `Finished` once its input is drained.
    #[derive(Clone)]
    struct AddHundred;
    impl Step for AddHundred {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "AddHundred",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::Unbounded],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(v) => {
                    let _ = ctx.outputs.push(v + 100);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// Sink: collects popped values into a shared `Vec`; `Finished` on drain.
    struct CollectSink {
        out: Arc<Mutex<Vec<u32>>>,
    }
    impl Step for CollectSink {
        type Input = u32;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "CollectSink",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(v) => {
                    self.out.lock().unwrap().push(v);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// Build a linear `source → mid → sink` chain (graph + boxed steps).
    fn linear_chain(
        count: u32,
        out: &Arc<Mutex<Vec<u32>>>,
    ) -> (Vec<Box<dyn ErasedStep>>, ChainGraph) {
        let mut graph = ChainGraph::new();
        let s = graph.register_step("CountSource", 1);
        let m = graph.register_step("AddHundred", 1);
        let k = graph.register_step("CollectSink", 0);
        graph.wire(s, BranchIdx(0), m);
        graph.wire(m, BranchIdx(0), k);
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(CountSource { next: 0, count })),
            Box::new(TypedStep::new(AddHundred)),
            Box::new(TypedStep::new(CollectSink { out: Arc::clone(out) })),
        ];
        (steps, graph)
    }

    /// Fan-out split: routes even values to branch 0, odd to branch 1 — the
    /// shape of a `--rejects` kept/rejects split.
    struct EvenOddSplit;
    impl Step for EvenOddSplit {
        type Input = u32;
        type Outputs = (u32, u32);
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "EvenOddSplit",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::Unbounded, QueueSpec::Unbounded],
                branch_ordering: vec![BranchOrdering::None, BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(v) => {
                    let view = ctx.outputs.view();
                    if v % 2 == 0 {
                        let _ = view.a.push(v);
                    } else {
                        let _ = view.b.push(v);
                    }
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// Build a fan-out `source → split → (even sink, odd sink)` chain.
    fn fan_out_chain(
        count: u32,
        even: &Arc<Mutex<Vec<u32>>>,
        odd: &Arc<Mutex<Vec<u32>>>,
    ) -> (Vec<Box<dyn ErasedStep>>, ChainGraph) {
        let mut graph = ChainGraph::new();
        let s = graph.register_step("CountSource", 1);
        let m = graph.register_step("EvenOddSplit", 2);
        let k0 = graph.register_step("CollectSink", 0);
        let k1 = graph.register_step("CollectSink", 0);
        graph.wire(s, BranchIdx(0), m);
        graph.wire(m, BranchIdx(0), k0);
        graph.wire(m, BranchIdx(1), k1);
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(CountSource { next: 0, count })),
            Box::new(TypedStep::new(EvenOddSplit)),
            Box::new(TypedStep::new(CollectSink { out: Arc::clone(even) })),
            Box::new(TypedStep::new(CollectSink { out: Arc::clone(odd) })),
        ];
        (steps, graph)
    }

    #[test]
    fn is_fusible_detects_source_mid_sink() {
        let out = Arc::new(Mutex::new(Vec::new()));
        let (steps, graph) = linear_chain(3, &out);
        assert!(is_fusible_chain(&steps, &graph));
    }

    #[test]
    fn is_fusible_rejects_single_step() {
        // A self-contained source+sink (Input=(), no output branches) is one
        // step — nothing to fuse, so not eligible.
        let mut graph = ChainGraph::new();
        graph.register_step("CountSource", 0);
        let steps: Vec<Box<dyn ErasedStep>> =
            vec![Box::new(TypedStep::new(CountSource { next: 0, count: 0 }))];
        assert!(!is_fusible_chain(&steps, &graph));
    }

    #[test]
    fn is_fusible_rejects_empty() {
        let steps: Vec<Box<dyn ErasedStep>> = vec![];
        let graph = ChainGraph::new();
        assert!(!is_fusible_chain(&steps, &graph));
    }

    #[test]
    fn is_fusible_accepts_fan_out() {
        // A fan-out (the `--rejects` shape: one step with two output branches,
        // each wired forward to its own sink) IS fusible — both branches feed
        // strictly later steps.
        let even = Arc::new(Mutex::new(Vec::new()));
        let odd = Arc::new(Mutex::new(Vec::new()));
        let (steps, graph) = fan_out_chain(0, &even, &odd);
        assert!(is_fusible_chain(&steps, &graph));
    }

    #[test]
    fn drive_runs_chain_to_completion_in_order() {
        let out = Arc::new(Mutex::new(Vec::new()));
        let (steps, graph) = linear_chain(5, &out);
        let signal = PipelineSignal::new();
        run_fused_single_thread(steps, &graph, &signal, None).expect("clean run");
        // Source emits 0..5, mid adds 100, sink collects in FIFO order.
        assert_eq!(*out.lock().unwrap(), vec![100, 101, 102, 103, 104]);
    }

    #[test]
    fn drive_runs_fan_out_to_completion() {
        // Both sink subchains of a fan-out must drain — the driver waits for ALL
        // steps to finish, not just the last-indexed one.
        let even = Arc::new(Mutex::new(Vec::new()));
        let odd = Arc::new(Mutex::new(Vec::new()));
        let (steps, graph) = fan_out_chain(6, &even, &odd);
        let signal = PipelineSignal::new();
        run_fused_single_thread(steps, &graph, &signal, None).expect("clean run");
        // Source emits 0..6; evens route to branch 0, odds to branch 1.
        assert_eq!(*even.lock().unwrap(), vec![0, 2, 4]);
        assert_eq!(*odd.lock().unwrap(), vec![1, 3, 5]);
    }

    #[test]
    fn drive_propagates_step_error() {
        /// Mid that errors on the first item.
        struct Boom;
        impl Step for Boom {
            type Input = u32;
            type Outputs = Single<u32>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "Boom",
                    kind: StepKind::Serial,
                    sticky: false,
                    output_queues: vec![QueueSpec::Unbounded],
                    branch_ordering: vec![BranchOrdering::None],
                }
            }
            fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                if ctx.input.pop().is_some() {
                    return Err(io::Error::other("boom"));
                }
                if ctx.input.is_drained() {
                    Ok(StepOutcome::Finished)
                } else {
                    Ok(StepOutcome::NoProgress)
                }
            }
        }

        let out = Arc::new(Mutex::new(Vec::new()));
        let mut graph = ChainGraph::new();
        let s = graph.register_step("CountSource", 1);
        let m = graph.register_step("Boom", 1);
        let k = graph.register_step("CollectSink", 0);
        graph.wire(s, BranchIdx(0), m);
        graph.wire(m, BranchIdx(0), k);
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(CountSource { next: 0, count: 3 })),
            Box::new(TypedStep::new(Boom)),
            Box::new(TypedStep::new(CollectSink { out })),
        ];
        let signal = PipelineSignal::new();
        let err = run_fused_single_thread(steps, &graph, &signal, None).expect_err("must error");
        assert!(matches!(err, PipelineError::Io { step: "Boom", .. }));
    }
}
