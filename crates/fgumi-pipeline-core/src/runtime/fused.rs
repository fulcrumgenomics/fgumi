//! Single-thread *fused* execution mode (issue #330).
//!
//! At `--threads 1` a forward-wired `source → … → sink` chain gains nothing from
//! the scheduled worker pool: there is one worker, so the inter-step bounded
//! queues, round-robin polling, held-slot retries, and reorder bookkeeping are
//! pure overhead (profiling showed ~2/3 of `try_run` calls do no useful work).
//!
//! This module is the structural fix. Fusion is **not** a per-command rewrite
//! — it is an execution mode of the existing pipeline. [`is_fusible_chain`]
//! detects a fusible chain (forward-wired, fan-out allowed);
//! [`run_fused_single_thread`] then drives the chain's
//! own type-erased steps inline, in topological order, over **direct** buffers
//! (built by [`build_chain_contexts_fused`]). FIFO push order is already the
//! correct order at one worker, so the reorder stage is dropped; the profile's
//! count/byte bound is kept, because the driver runs a producer before its
//! consumer and a step emitting more per `try_run` than its consumer removes
//! would otherwise grow the edge on every pass. The runtime's
//! generic wiring (`build_chain_contexts`) and dispatch (`try_run_erased`) are
//! reused verbatim — no step logic is duplicated.
//!
//! Non-linear chains (two-input `Step2` merges like zipper, multi-output splits
//! like `correct --rejects`) and `--threads ≥ 2` are not eligible and fall back
//! to the scheduled [`run_worker_loop`](super::driver::run_worker_loop).

use std::sync::Arc;
use std::time::{Duration, Instant};

use super::contexts::build_chain_contexts_fused;
use super::stats::PipelineStats;
use crate::builder::InstrumentationLevel;
use crate::erased::{ErasedStep, ErasedStepCtx};
use crate::signal::{PipelineError, PipelineSignal};
use crate::step::StepOutcome;
use crate::topology::{BranchIdx, ChainGraph, StepIdx};

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
        // one input stream. Sources register arity 0 (their input is implicit),
        // single-input steps 1, and `Step2` 2 — so anything above 1 is a merge.
        if graph.input_arity(idx) > 1 {
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

/// Whether `build_chain`'s fused single-thread fast path may be taken for this
/// chain.
///
/// The fast path drives a fusible chain inline over direct buffers, skipping the
/// scheduler — and with it the per-edge instrumentation the scheduled path sets
/// up (edge [`EdgeMetrics`](super::metrics::EdgeMetrics), the occupancy sampler,
/// and the `snapshot_with_edges` bottleneck verdict). An instrumented run
/// (`InstrumentationLevel != Off`) must therefore NOT fuse, or its
/// `--pipeline-stats` / `--pipeline-trace` output would silently omit all edge /
/// occupancy data. When instrumentation is `Off` (the default), fusion is taken
/// whenever the chain is single-thread and fusible (zero-overhead fast path
/// preserved).
#[must_use]
pub fn should_fuse_single_thread(
    n_threads: usize,
    instrumentation: InstrumentationLevel,
    steps: &[Box<dyn ErasedStep>],
    graph: &ChainGraph,
) -> bool {
    n_threads == 1 && !instrumentation.is_on() && is_fusible_chain(steps, graph)
}

/// Drive a fusible chain to completion on the calling thread, fused.
///
/// Builds the chain's per-step contexts with direct inter-step transports (no
/// reorder stage, profile queue bounds retained — see
/// [`ErasedStep::build_fused_output_set`]), then repeatedly walks the steps in
/// topological order — popping
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
/// `queue_memory_total` is the run's `--queue-memory-total`, if set. Because the
/// fused contexts are built here rather than by the caller, this function is the
/// only place that budget can be applied to the fused transports.
///
/// `deadlock_timeout_secs` is the run's stall patience. A pass in which no step
/// progresses is retried with a backoff until this budget is exhausted, and only
/// then reported as a stall; `0` selects a built-in default rather than
/// "unbounded" (see the constant's comment for why).
///
/// # Errors
///
/// Returns [`PipelineError::Io`] if any step's `try_run` returned `Err` (the
/// first such error wins, carrying the originating step's name), or
/// [`PipelineError::Cancelled`] if the run was cancelled via the pipeline's
/// [`CancelHandle`](crate::signal::CancelHandle) — matching
/// [`crate::builder::Pipeline::run`]'s contract.
pub fn run_fused_single_thread(
    mut steps: Vec<Box<dyn ErasedStep>>,
    graph: &ChainGraph,
    signal: &Arc<PipelineSignal>,
    stats: Option<&Arc<PipelineStats>>,
    queue_memory_total: Option<u64>,
    deadlock_timeout_secs: u64,
) -> Result<(), PipelineError> {
    // A single no-progress pass is NOT a stall. `NoProgress` is the transient
    // "input momentarily empty but not drained" outcome — a source waiting on a
    // background reader thread returns it legitimately — so failing on the first
    // idle pass aborts a healthy run and truncates its output. Tolerate idling
    // until this wall-clock budget is exhausted, which is what the scheduled path
    // does via `WorkerCore::sleep_backoff`. That path has no stall limit at all
    // because the deadlock monitor catches wedges for it; the fused path is not
    // monitored, so it needs its own bound rather than hanging forever.
    //
    // Budgeted by TIME, not by a pass count: `thread::sleep` granularity varies by
    // platform (a requested 50µs can round up to ~1ms), so a fixed number of
    // passes would mean wildly different real budgets across hosts.
    //
    // The budget comes from `PipelineConfig::deadlock_timeout_secs`, so a chain
    // with a genuinely slow source (a network-backed reader, say) can raise it
    // instead of having a library-chosen default abort a legitimate run.
    //
    // `deadlock_timeout_secs == 0` — the config default, meaning "monitor
    // disarmed" on the scheduled path — falls back to `DEFAULT_STALL_BUDGET` here
    // rather than meaning "unbounded". That asymmetry is deliberate and is the
    // whole reason this bound exists: the scheduled path has a deadlock monitor to
    // arm, and this one does not, so disabling the bound would leave a fused wedge
    // with nothing at all to detect it. A caller wanting more patience raises the
    // number; there is deliberately no way to remove the bound.
    // Deliberately generous. This bound exists to catch a PERMANENT wedge, not to
    // police a slow source, and it is the value most runs get (`deadlock_timeout_secs`
    // defaults to 0). A source whose background reader blocks on a cold page cache or
    // network-backed input can legitimately idle for tens of seconds, and the idle
    // timer only resets on progress — so a tight default trades a real risk of failing
    // a healthy run for a few seconds off the report of a wedge that has already hung.
    // The costs are asymmetric; err long.
    const DEFAULT_STALL_BUDGET: Duration = Duration::from_secs(60);
    // Long enough to stop pegging the calling thread at 100%, short enough that it
    // adds no meaningful latency to a source that is about to produce.
    const IDLE_BACKOFF: Duration = Duration::from_micros(50);

    let stall_budget = if deadlock_timeout_secs > 0 {
        Duration::from_secs(deadlock_timeout_secs)
    } else {
        DEFAULT_STALL_BUDGET
    };

    let n = steps.len();
    let contexts = build_chain_contexts_fused(&steps, graph);
    // Honour `--queue-memory-total` here as the scheduled path does. The fused
    // transports keep each step's profiled byte bound (see
    // `ErasedStep::build_fused_output_set`), so without this the user's budget
    // would be silently ignored in favour of the per-step defaults. These
    // contexts are local to this call, so the budget cannot be applied by
    // `Pipeline::run` on its behalf.
    if let Some(total) = queue_memory_total {
        crate::builder::apply_initial_queue_budget(&contexts.bounded_queues, total);
    }
    let mut finished = vec![false; n];
    // `Some(t)` while the driver has been idle since `t`; cleared by any pass that
    // makes progress.
    let mut idle_since: Option<Instant> = None;

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
        if progressed {
            idle_since = None;
        } else {
            // Idle pass. Back off and retry the same state until the budget runs
            // out — only a sustained run of idle passes is evidence of a wedge.
            let idle_start = *idle_since.get_or_insert_with(Instant::now);
            if idle_start.elapsed() < stall_budget {
                std::thread::sleep(IDLE_BACKOFF);
                continue 'drive;
            }
            // Budget exhausted. For a well-formed fusible chain this is
            // unreachable: the source either progresses or finishes, and
            // `Finished` cascades drain downstream so some step always advances
            // until every sink closes. Guard against an infinite spin rather than
            // trusting that invariant blindly.
            debug_assert!(
                false,
                "fused single-thread driver stalled: no step progressed for {stall_budget:?} \
                 and not all steps have finished"
            );
            // In release builds the `debug_assert!` is compiled out, so record
            // the stall as an error before breaking. Otherwise the loop would
            // exit and map to `Ok`, silently truncating the output instead of
            // surfacing the broken invariant. Name the still-unfinished steps so
            // the error points at the wedged step(s) rather than just the
            // synthetic "fused-driver" name.
            let stalled_steps: Vec<&'static str> = finished
                .iter()
                .enumerate()
                .filter(|(_, f)| !**f)
                .map(|(i, _)| steps[i].name())
                .collect();
            signal.record_error(PipelineError::Io {
                step: "fused-driver",
                source: std::io::Error::other(format!(
                    "fused single-thread driver stalled: no step progressed for {stall_budget:?} \
                     and not all steps have finished; unfinished step(s): {}",
                    stalled_steps.join(", "),
                )),
            });
            break 'drive;
        }
    }

    // Map the recorded outcome to the run result (same shape as `Pipeline::run`).
    // `to_result` reconstructs the non-`Clone` `PipelineError` and synthesizes
    // `Cancelled` from the state when an external cancel published the terminal
    // state but its `OnceLock` payload is not yet visible to this thread.
    signal.to_result()
}

#[cfg(test)]
mod tests {
    use std::io;
    use std::sync::{Arc, Mutex};

    use rstest::rstest;

    use super::*;
    use crate::erased::TypedStep;
    use crate::outputs::Single;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;
    use crate::step::{Step, StepCtx, StepKind, StepProfile};

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

    /// The `input_arity(idx) > 1` guard: a chain containing a `Step2` merge
    /// (input arity 2) is NOT fusible — the single-worker inline drive follows
    /// one input stream per step. Every other condition holds (one leading
    /// source; all branches wired strictly forward), so the merge arity is the
    /// sole reason fusion is refused.
    #[test]
    fn is_fusible_rejects_a_two_input_merge() {
        let out = Arc::new(Mutex::new(Vec::new()));
        let mut graph = ChainGraph::new();
        let s = graph.register_step("CountSource", 1);
        let m = graph.register_step_with_input_arity("AddHundred", 1, 2);
        let k = graph.register_step("CollectSink", 0);
        graph.wire(s, BranchIdx(0), m);
        graph.wire(m, BranchIdx(0), k);
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(CountSource { next: 0, count: 0 })),
            Box::new(TypedStep::new(AddHundred)),
            Box::new(TypedStep::new(CollectSink { out: Arc::clone(&out) })),
        ];
        assert!(
            !is_fusible_chain(&steps, &graph),
            "a Step2 merge (input_arity 2) must block fusion"
        );
    }

    /// The forward-wiring guard: an output branch wired to an earlier-or-equal
    /// step index is NOT fusible — the inline drive is a single topological pass,
    /// so a back-edge would revisit an already-drained step. Here the mid step's
    /// branch loops back to the source (consumer index 0 ≤ producer index 1);
    /// every other condition holds, so the back-edge is the sole reason.
    #[test]
    fn is_fusible_rejects_a_backward_wired_branch() {
        let out = Arc::new(Mutex::new(Vec::new()));
        let mut graph = ChainGraph::new();
        let s = graph.register_step("CountSource", 1);
        let m = graph.register_step("AddHundred", 1);
        let _k = graph.register_step("CollectSink", 0);
        graph.wire(s, BranchIdx(0), m);
        graph.wire(m, BranchIdx(0), s); // back-edge: consumer index 0 ≤ producer index 1
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(CountSource { next: 0, count: 0 })),
            Box::new(TypedStep::new(AddHundred)),
            Box::new(TypedStep::new(CollectSink { out: Arc::clone(&out) })),
        ];
        assert!(
            !is_fusible_chain(&steps, &graph),
            "a branch wired to an earlier-or-equal step must block fusion"
        );
    }

    // A fusible chain fuses ONLY when it is single-thread AND uninstrumented: the
    // fused fast path skips the scheduled path's edge metrics / occupancy sampler
    // / bottleneck verdict, so any instrumentation level (or ≥2 threads) must fall
    // through to the scheduled path instead of silently dropping that output.
    #[rstest]
    #[case::off_single_thread(1, InstrumentationLevel::Off, true)]
    #[case::summary_single_thread(1, InstrumentationLevel::Summary, false)]
    #[case::timeline_single_thread(1, InstrumentationLevel::Timeline, false)]
    #[case::deep_single_thread(1, InstrumentationLevel::Deep, false)]
    #[case::off_multi_thread(2, InstrumentationLevel::Off, false)]
    fn should_fuse_only_when_uninstrumented_single_thread(
        #[case] n_threads: usize,
        #[case] level: InstrumentationLevel,
        #[case] expected: bool,
    ) {
        let out = Arc::new(Mutex::new(Vec::new()));
        let (steps, graph) = linear_chain(4, &out);
        assert!(is_fusible_chain(&steps, &graph), "linear chain is fusible");
        assert_eq!(should_fuse_single_thread(n_threads, level, &steps, &graph), expected);
    }

    #[test]
    fn drive_runs_chain_to_completion_in_order() {
        let out = Arc::new(Mutex::new(Vec::new()));
        let (steps, graph) = linear_chain(5, &out);
        let signal = PipelineSignal::new();
        run_fused_single_thread(steps, &graph, &signal, None, None, 0).expect("clean run");
        // Source emits 0..5, mid adds 100, sink collects in FIFO order.
        assert_eq!(*out.lock().unwrap(), vec![100, 101, 102, 103, 104]);
    }

    /// Source with a *byte-bounded* output, so a fused chain built from it
    /// registers a real bounded transport. Emits `0, 1, …, count-1`, holding and
    /// retrying whatever the transport rejects — the contract every step already
    /// owes the scheduled path.
    struct ByteBoundedSource {
        next: u32,
        count: u32,
        held: Option<u32>,
    }
    impl Step for ByteBoundedSource {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ByteBoundedSource",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::ByteBounded { limit_bytes: PROFILE_LIMIT_BYTES }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(item) = self.held.take() {
                if let Err(unpushed) = ctx.outputs.push(item) {
                    self.held = Some(unpushed.into_item());
                }
                return Ok(StepOutcome::Progress);
            }
            if self.next >= self.count {
                return Ok(StepOutcome::Finished);
            }
            let item = self.next;
            self.next += 1;
            if let Err(unpushed) = ctx.outputs.push(item) {
                self.held = Some(unpushed.into_item());
            }
            Ok(StepOutcome::Progress)
        }
    }

    /// Deliberately unusual so the assertions below cannot pass by matching a
    /// default or a budget-derived value.
    const PROFILE_LIMIT_BYTES: u64 = 7 * 1024 * 1024;

    fn byte_bounded_chain(
        count: u32,
        out: &Arc<Mutex<Vec<u32>>>,
    ) -> (Vec<Box<dyn ErasedStep>>, ChainGraph) {
        let mut graph = ChainGraph::new();
        let src = graph.register_step("ByteBoundedSource", 1);
        let sink = graph.register_step("CollectSink", 0);
        graph.wire(src, BranchIdx(0), sink);
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(ByteBoundedSource { next: 0, count, held: None })),
            Box::new(TypedStep::new(CollectSink { out: Arc::clone(out) })),
        ];
        (steps, graph)
    }

    /// The fused path must honour `--queue-memory-total`, as the scheduled path
    /// does. Two facts are load-bearing and both were false before fused
    /// transports kept their profile bounds: the fused contexts register a
    /// bounded queue at all, and the budget resizes it off the profile default.
    #[test]
    fn fused_contexts_register_bounded_queues_and_take_the_budget() {
        let out = Arc::new(Mutex::new(Vec::new()));
        let (steps, graph) = byte_bounded_chain(4, &out);
        let contexts = build_chain_contexts_fused(&steps, &graph);

        assert_eq!(
            contexts.bounded_queues.len(),
            1,
            "a byte-bounded fused edge must be registered, else there is nothing \
             for the queue-memory budget to apply to"
        );
        assert_eq!(
            contexts.bounded_queues[0].handle.limit_bytes(),
            PROFILE_LIMIT_BYTES,
            "the fused transport starts at the profile's declared bound"
        );

        // 64 MiB over one queue → 64 MiB per queue, well above the 1 MiB floor,
        // and distinct from the profile default so the assert cannot pass by
        // accident.
        let total = 64 * 1024 * 1024;
        crate::builder::apply_initial_queue_budget(&contexts.bounded_queues, total);
        assert_eq!(
            contexts.bounded_queues[0].handle.limit_bytes(),
            total,
            "the user's budget must override the per-step default on the fused path"
        );
    }

    /// End-to-end: a byte-bounded fused chain run with a budget still delivers
    /// every item in order. Pins that threading the budget through does not wedge
    /// the driver — a bounded fused edge relies on the producer's hold-and-retry.
    #[rstest]
    #[case::no_budget(None)]
    #[case::with_budget(Some(2 * 1024 * 1024))]
    fn fused_byte_bounded_chain_completes(#[case] queue_memory_total: Option<u64>) {
        let out = Arc::new(Mutex::new(Vec::new()));
        let (steps, graph) = byte_bounded_chain(5, &out);
        let signal = PipelineSignal::new();
        run_fused_single_thread(steps, &graph, &signal, None, queue_memory_total, 0)
            .expect("clean run");
        assert_eq!(*out.lock().unwrap(), vec![0, 1, 2, 3, 4]);
    }

    /// Source that reports `NoProgress` for its first few dispatches before it
    /// starts emitting — the shape of a source waiting on a background reader.
    struct IdleThenEmitSource {
        idle_left: u32,
        remaining: u32,
        dispatches: Arc<Mutex<u32>>,
    }
    impl Step for IdleThenEmitSource {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "IdleThenEmitSource",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            *self.dispatches.lock().unwrap() += 1;
            if self.idle_left > 0 {
                self.idle_left -= 1;
                // Legitimate transient: nothing to hand over *yet*.
                return Ok(StepOutcome::NoProgress);
            }
            if self.remaining == 0 {
                return Ok(StepOutcome::Finished);
            }
            let item = self.remaining;
            self.remaining -= 1;
            let _ = ctx.outputs.push(item);
            Ok(StepOutcome::Progress)
        }
    }

    /// A pass in which no step progressed is not a stall — `NoProgress` is the
    /// transient "input momentarily empty but not drained" outcome. The driver
    /// used to fail the run on the FIRST such pass, recording `PipelineError::Io`
    /// and truncating the output (and panicking through the `debug_assert!` in
    /// test builds). It must back off and retry instead, so a source that idles
    /// before producing still completes.
    #[test]
    fn drive_tolerates_transient_no_progress_passes() {
        const IDLE_PASSES: u32 = 5;
        const N_ITEMS: u32 = 3;

        let mut graph = ChainGraph::new();
        let src = graph.register_step("IdleThenEmitSource", 1);
        let sink = graph.register_step("CollectSink", 0);
        graph.wire(src, BranchIdx(0), sink);

        let out = Arc::new(Mutex::new(Vec::new()));
        let dispatches = Arc::new(Mutex::new(0));
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(IdleThenEmitSource {
                idle_left: IDLE_PASSES,
                remaining: N_ITEMS,
                dispatches: Arc::clone(&dispatches),
            })),
            Box::new(TypedStep::new(CollectSink { out: Arc::clone(&out) })),
        ];
        let signal = PipelineSignal::new();
        run_fused_single_thread(steps, &graph, &signal, None, None, 0)
            .expect("transient NoProgress must not fail the run");

        assert_eq!(
            *out.lock().unwrap(),
            vec![3, 2, 1],
            "every item must still be delivered after the idle passes"
        );
        assert!(!signal.is_done(), "no error recorded for a transient idle pass");
        assert!(
            *dispatches.lock().unwrap() > IDLE_PASSES,
            "the driver must have retried past the idle passes, not given up on the first"
        );
    }

    /// Source that never progresses and never finishes — a permanent wedge, the
    /// case the stall bound exists for.
    struct WedgedSource;
    impl Step for WedgedSource {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "WedgedSource",
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

    /// The stall budget is `deadlock_timeout_secs`, and `0` selects the built-in
    /// default rather than "unbounded" — a fused wedge has no deadlock monitor to
    /// fall back on, so the bound must not be disableable. A 1-second budget is
    /// used here to keep the test fast; the assertion is that the run fails inside
    /// its own budget and names the wedged step.
    ///
    /// Gated `#[should_panic]` rather than an unconditional one: the
    /// `debug_assert!` on the stall path fires first in a debug test build, which
    /// is itself the contract (a genuine wedge must be loud in tests). In a
    /// release test build that `debug_assert!` is compiled out, so the driver
    /// records the stall and returns `Err` instead — the in-body assertions cover
    /// that path, and the `should_panic` expectation must not apply there.
    #[test]
    #[cfg_attr(debug_assertions, should_panic(expected = "fused single-thread driver stalled"))]
    fn drive_reports_a_real_stall_within_the_configured_budget() {
        let mut graph = ChainGraph::new();
        let src = graph.register_step("WedgedSource", 1);
        let sink = graph.register_step("CollectSink", 0);
        graph.wire(src, BranchIdx(0), sink);

        let out = Arc::new(Mutex::new(Vec::new()));
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(WedgedSource)),
            Box::new(TypedStep::new(CollectSink { out })),
        ];
        let signal = PipelineSignal::new();
        let started = Instant::now();
        let result = run_fused_single_thread(steps, &graph, &signal, None, None, 1);
        // Only reached in a release test build, where the `debug_assert!` is gone.
        assert!(result.is_err(), "a permanent wedge must not report success");
        assert!(
            started.elapsed() < Duration::from_secs(30),
            "the configured 1s budget must bound the wait, not the default"
        );
    }

    #[test]
    fn drive_runs_fan_out_to_completion() {
        // Both sink subchains of a fan-out must drain — the driver waits for ALL
        // steps to finish, not just the last-indexed one.
        let even = Arc::new(Mutex::new(Vec::new()));
        let odd = Arc::new(Mutex::new(Vec::new()));
        let (steps, graph) = fan_out_chain(6, &even, &odd);
        let signal = PipelineSignal::new();
        run_fused_single_thread(steps, &graph, &signal, None, None, 0).expect("clean run");
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
        let err =
            run_fused_single_thread(steps, &graph, &signal, None, None, 0).expect_err("must error");
        assert!(matches!(err, PipelineError::Io { step: "Boom", .. }));
    }

    /// The `# Errors` contract promises `PipelineError::Cancelled` when the run
    /// is cancelled through the `CancelHandle`, and the driver's only delivery
    /// of it is the top-of-loop `signal.is_done()` break mapped by
    /// `to_result()`. Every other test drives a clean run, a budget, a transient
    /// idle, a stall, or a step error — none a cancel — so a regression that
    /// moved or dropped this check (for example, breaking only once every step
    /// finished) would keep the whole suite green while a cancelled fused run
    /// reported `Ok`.
    #[test]
    fn drive_maps_an_external_cancel_to_cancelled() {
        /// Source that cancels the shared signal on its first dispatch, then
        /// keeps emitting. Models an external `CancelHandle::cancel` landing
        /// mid-run: the typed `StepCtx` deliberately does not expose the signal,
        /// so the step holds its own clone, exactly as an external canceller
        /// (which owns a `CancelHandle` over the same signal) does.
        struct CancelOnFirstDispatch {
            signal: Arc<PipelineSignal>,
            next: u32,
            count: u32,
        }
        impl Step for CancelOnFirstDispatch {
            type Input = ();
            type Outputs = Single<u32>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "CancelOnFirstDispatch",
                    kind: StepKind::Serial,
                    sticky: false,
                    output_queues: vec![QueueSpec::Unbounded],
                    branch_ordering: vec![BranchOrdering::None],
                }
            }
            fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                if self.next == 0 {
                    self.signal.cancel();
                }
                if self.next >= self.count {
                    return Ok(StepOutcome::Finished);
                }
                let _ = ctx.outputs.push(self.next);
                self.next += 1;
                Ok(StepOutcome::Progress)
            }
        }

        const COUNT: u32 = 8;
        let out = Arc::new(Mutex::new(Vec::new()));
        let signal = PipelineSignal::new();
        let mut graph = ChainGraph::new();
        let s = graph.register_step("CancelOnFirstDispatch", 1);
        let k = graph.register_step("CollectSink", 0);
        graph.wire(s, BranchIdx(0), k);
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(CancelOnFirstDispatch {
                signal: Arc::clone(&signal),
                next: 0,
                count: COUNT,
            })),
            Box::new(TypedStep::new(CollectSink { out: Arc::clone(&out) })),
        ];
        let err = run_fused_single_thread(steps, &graph, &signal, None, None, 0)
            .expect_err("a cancelled run must not report success");
        assert!(
            matches!(err, PipelineError::Cancelled),
            "an external cancel must map to PipelineError::Cancelled, got {err:?}"
        );
        assert!(
            out.lock().unwrap().len() < COUNT as usize,
            "the driver must break on the cancel before draining every item"
        );
    }
}
