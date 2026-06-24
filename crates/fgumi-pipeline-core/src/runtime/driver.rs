//! Worker loop body. Each worker thread runs `run_worker_loop` until
//! `signal.is_done()` or all steps are drained.
//!
//! Loop structure (per iteration):
//!   1. Check `signal.is_done()`; bail if true.
//!   2. **Sticky re-entry**: if this worker owns an `Exclusive` step that's
//!      sticky, drive it to a stop (`Progress`→loop; `Finished`→done;
//!      `NoProgress` or `Contention`→exit sticky). Sticky avoids context-
//!      switch overhead for source-style steps that emit in tight bursts.
//!   3. **Round-robin dispatch**: try each step in chain order. On
//!      `Progress`, restart from step 0 (priority). On `Finished`, the step is
//!      complete — `mark_outputs_drained` (counter-gated for `Parallel`) and
//!      remove it from the worklist. `NoProgress`/`Contention` are idle ticks.
//!   4. If no work happened this iteration, exponential-backoff sleep.
//!
//! Completion: every step — source, mid, or sink — terminates by returning
//! `StepOutcome::Finished` from `try_run` once its input edges are drained and
//! it holds no buffered output. The framework then closes its output edges and
//! drops it from the per-worker worklist. For a `Parallel` step the per-step
//! [`StepDrainCounter`] gates `mark_outputs_drained` so only the last clone to
//! finish closes the shared output queue (see `dispatch_one_step`).

use std::any::Any;
use std::sync::Arc;
use std::time::Instant;

use crate::erased::ErasedStepCtx;
use crate::runtime::contexts::ChainContexts;
use crate::runtime::drain::StepDrainCounter;
use crate::runtime::live::LiveSteps;
use crate::runtime::stats::PipelineStats;
use crate::runtime::storage::WorkerStepEntry;
use crate::runtime::worker_core::WorkerCore;
use crate::signal::{PipelineError, PipelineSignal};
use crate::step::StepOutcome;
use crate::topology::StepIdx;

/// Run the worker loop for one worker thread.
///
/// `entries[step_idx] = WorkerStepEntry` — this worker's storage.
/// `contexts` — shared per-step input/output handles.
/// `drain_counters[step_idx]` — the per-step `StepDrainCounter` that gates the
/// output close on `Finished` (init N for Parallel so only the last clone
/// closes the shared output; init 1 for Serial/Exclusive).
/// `signal` — error/cancel broadcast.
pub fn run_worker_loop(
    worker: &mut WorkerCore,
    entries: &mut [WorkerStepEntry],
    contexts: &Arc<ChainContexts>,
    drain_counters: &[Arc<StepDrainCounter>],
    signal: &Arc<PipelineSignal>,
    stats: Option<&Arc<PipelineStats>>,
) {
    // Per-worker worklist of still-dispatchable steps, in chain order. A step
    // is removed when it returns `StepOutcome::Finished`; the worker exits once
    // the list is empty. Build-time `Skip` placeholders (Exclusive steps owned
    // by other workers, Serial steps gated out by affinity) never enter the
    // worklist.
    let mut live = LiveSteps::from_entries(entries);

    // Cache whether this worker's sticky owner (if any) is still dispatchable,
    // so the hot sticky fast-path does not run a linear `live.contains` scan on
    // every outer-loop iteration (the sticky path exists precisely to shave
    // per-iteration overhead). `sticky_owner` is fixed for the worker's
    // lifetime; the step leaves `live` exactly when it returns `Finished`,
    // either via the sticky block below or via round-robin dispatch — both
    // sites flip this flag false. A `None` sticky owner is permanently "not
    // live" so the fast path is skipped entirely.
    let mut sticky_live = worker.sticky_owner.is_some_and(|idx| live.contains(idx));

    loop {
        if signal.is_done() {
            break;
        }
        // Exit when this worker has nothing left to dispatch.
        if live.is_empty() {
            break;
        }

        let mut did_work = false;

        // 1. Sticky re-entry for sticky-owned steps (either an Exclusive
        // sticky step this worker owns, or a Serial+sticky step whose
        // Affinity targets this worker). The Pipeline::run path only
        // sets `sticky_owner` when the step is actually sticky, so no
        // per-iteration profile peek is needed. Re-enter while the
        // step makes Progress; exit on Finished / NoProgress / Contention
        // / Err. Remove from the worklist on Finished (source drain) or
        // observed input drain (mid-step drain). Gated on the step still
        // being live — once removed, the sticky fast-path is disabled.
        if let Some(owned_idx) = worker.sticky_owner.filter(|_| sticky_live) {
            let mut mark_skip = false;
            loop {
                if signal.is_done() {
                    break;
                }
                let entry = &mut entries[owned_idx.0];
                let Some(info) = dispatch_one_step(
                    entry,
                    owned_idx,
                    contexts,
                    &drain_counters[owned_idx.0],
                    signal,
                    stats,
                ) else {
                    break; // Skip
                };
                match info.result {
                    Ok(StepOutcome::Progress) => {
                        did_work = true;
                        // Continue sticky.
                    }
                    Ok(StepOutcome::Finished) => {
                        // Outputs were marked drained under the dispatch guard.
                        did_work = true;
                        mark_skip = true;
                        break;
                    }
                    // Nothing to do this call — yield out of the sticky loop
                    // back to round-robin. The step terminates via `Finished`,
                    // not a drain protocol.
                    Ok(StepOutcome::NoProgress | StepOutcome::Contention) => break,
                    Err(io_err) => {
                        signal.record_error(PipelineError::Io { step: info.name, source: io_err });
                        break;
                    }
                }
            }
            if mark_skip {
                live.remove(owned_idx);
                // The sticky owner finished here; disable the fast path.
                sticky_live = false;
            }
        }

        // 2. Round-robin priority dispatch over all live steps.
        if !signal.is_done() {
            let outcome = round_robin_dispatch(
                entries,
                &mut live,
                worker.sticky_owner,
                contexts,
                drain_counters,
                signal,
                stats,
            );
            did_work |= outcome.did_work;
            if outcome.removed_sticky_owner {
                // The sticky owner finished during round-robin; disable the
                // fast path so subsequent iterations skip the sticky block.
                sticky_live = false;
            }
        }

        // 3. Exponential-backoff sleep on no-progress; reset on progress.
        if did_work {
            worker.reset_backoff();
        } else if signal.is_done() {
            break;
        } else {
            worker.sleep_backoff();
            worker.increase_backoff();
        }
    }
}

/// Result of one round-robin pass: whether any step did useful work (caller
/// resets backoff), and whether the worker's sticky owner finished during the
/// pass (caller clears its `sticky_live` cache so the sticky fast-path is not
/// re-attempted on a removed step).
struct RoundRobinOutcome {
    did_work: bool,
    removed_sticky_owner: bool,
}

/// One pass of the round-robin dispatch over this worker's live steps, in
/// chain order. Finished steps are removed from `live` at end-of-pass (deferred
/// so the in-progress walk over `live.order()` is not mutated underneath it).
/// `sticky_owner` (if any) is reported back via
/// [`RoundRobinOutcome::removed_sticky_owner`] when it finishes here, so the
/// caller can disable the per-iteration sticky fast-path without a linear
/// `live.contains` scan.
fn round_robin_dispatch(
    entries: &mut [WorkerStepEntry],
    live: &mut LiveSteps,
    sticky_owner: Option<StepIdx>,
    contexts: &Arc<ChainContexts>,
    drain_counters: &[Arc<StepDrainCounter>],
    signal: &Arc<PipelineSignal>,
    stats: Option<&Arc<PipelineStats>>,
) -> RoundRobinOutcome {
    let mut did_work = false;
    // Steps that finished this pass, removed from `live` after the walk. A
    // step is visited at most once per pass (the cursor only advances; we
    // `break` on `Progress`/error, never revisit), so deferring removal is
    // safe and avoids reorder-under-iteration.
    let mut finished: Vec<StepIdx> = Vec::new();
    for pos in 0..live.len() {
        if signal.is_done() {
            break;
        }
        let step_idx = live.order()[pos];
        let entry = &mut entries[step_idx.0];
        let mut mark_skip = false;
        let mut restart_priority = false;
        let Some(info) = dispatch_one_step(
            entry,
            step_idx,
            contexts,
            &drain_counters[step_idx.0],
            signal,
            stats,
        ) else {
            continue; // Skip (build-time placeholder; should not appear in `live`)
        };
        match info.result {
            Ok(StepOutcome::Progress) => {
                did_work = true;
                restart_priority = true;
            }
            // Nothing to do this call. The step terminates by returning
            // `Finished` (handled below); `NoProgress`/`Contention` are idle
            // ticks — there is no separate drain protocol.
            Ok(StepOutcome::NoProgress | StepOutcome::Contention) => {}
            Ok(StepOutcome::Finished) => {
                // Any step (source, mid, or sink) may report `Finished` once
                // all its inputs are drained and it holds no buffered output.
                // Outputs were marked drained under the dispatch guard (and, for
                // a Serial step, the shared `finished` latch was set so the
                // other workers stop re-dispatching it — see `dispatch_one_step`).
                did_work = true;
                mark_skip = true;
            }
            Err(io_err) => {
                signal.record_error(PipelineError::Io { step: info.name, source: io_err });
                break;
            }
        }
        if mark_skip {
            finished.push(step_idx);
        }
        if restart_priority {
            break;
        }
    }
    let removed_sticky_owner = sticky_owner.is_some_and(|owner| finished.contains(&owner));
    for step_idx in finished {
        live.remove(step_idx);
    }
    RoundRobinOutcome { did_work, removed_sticky_owner }
}

/// Outcome of dispatching one step, plus the `name` captured *during* the
/// dispatch — under the same `Shared`-mutex guard as the run itself — for
/// error reporting.
struct DispatchInfo {
    result: std::io::Result<StepOutcome>,
    name: &'static str,
}

/// Dispatch one step's `try_run_erased`. Returns:
///   - `Some(DispatchInfo)` — dispatched (the `result` carries the outcome or
///     the step's `Err`); on `Finished`, outputs are already marked drained.
///   - `None` — entry is `Skip` (caller continues to next step).
///
/// For a `Serial` (`Shared`) step the shared `finished` latch on its
/// `DrainGate` is consulted *before* acquiring the step mutex: once any worker
/// has finished the step (returned `Finished`, or completed its cooperative
/// drain), the latch is set and every other worker short-circuits to a synthetic
/// `Finished` here rather than re-`try_lock`-ing and re-running an already-done
/// step. The winning worker sets the latch under the dispatch guard before
/// `mark_outputs_drained`, so a non-idempotent flusher can never be re-entered.
fn dispatch_one_step(
    entry: &mut WorkerStepEntry,
    step_idx: StepIdx,
    contexts: &ChainContexts,
    counter: &StepDrainCounter,
    signal: &Arc<PipelineSignal>,
    stats: Option<&Arc<PipelineStats>>,
) -> Option<DispatchInfo> {
    let outputs_any: &(dyn Any + Send + Sync) = contexts.outputs[step_idx.0].as_ref();
    let mut ctx =
        ErasedStepCtx { input: contexts.inputs[step_idx.0].as_ref(), outputs: outputs_any, signal };

    // Time the dispatch only when stats collection is on. `Instant::now()`
    // is ~20-50ns on Apple Silicon, ~50-100ns on x86_64; gating on
    // `stats.is_some()` keeps the no-stats path zero-cost.
    let start = stats.map(|_| Instant::now());

    // Run the step and capture `name` *while still holding the `Shared` guard*
    // (or with direct `&mut` for Owned/Exclusive). On `Finished`, mark outputs
    // drained here too — under the same guard — so the caller never re-acquires
    // the lock for any post-dispatch inspection.
    //
    // `mark_outputs_drained` is gated behind `counter.observe_drain()` (the
    // per-step `StepDrainCounter`): for a `Parallel` step (counter init N)
    // every clone returns `Finished` independently when the shared input edge
    // is drained, but only the LAST clone to finish (the one that takes the
    // counter to 0) closes the shared output queue — otherwise a clone could
    // `mark_drained` while a sibling is still pushing (`try_push`-after-drained
    // panic). For `Serial`/`Exclusive` (counter init 1) the single finisher
    // wins on its first call, unchanged.
    //
    // INVARIANT: for a `Parallel` step, `counter` init == clone count == the
    // worker count, and a clone leaves its worklist ONLY by returning
    // `Finished`, so the counter reaches 0 exactly when every clone has
    // finished. Any future scheduler change that removes a Parallel clone for
    // another reason (work-stealing, per-worker early exit) — or makes a source
    // `Parallel` — would leave the counter stuck above 0 and never close the
    // shared output, hanging the downstream consumer. Keep the init (builder.rs)
    // and this gate in lockstep.
    let info: Option<DispatchInfo> = match entry {
        WorkerStepEntry::Owned { step } | WorkerStepEntry::Exclusive { step } => {
            let result = step.try_run_erased(&mut ctx);
            // `name()` returns the cached static name — no per-dispatch
            // `StepProfile` (and its two `Vec`s) is built.
            let name = step.name();
            if matches!(result, Ok(StepOutcome::Finished)) && counter.observe_drain() {
                step.mark_outputs_drained(outputs_any);
            }
            Some(DispatchInfo { result, name })
        }
        WorkerStepEntry::Shared { step, drain } => {
            if drain.is_finished() {
                // Another worker already finished this Serial step. Don't
                // re-`try_lock`/re-run it — report a synthetic `Finished` so the
                // caller drops it from this worker's live set. Outputs were
                // already marked drained by the finishing worker.
                Some(DispatchInfo {
                    result: Ok(StepOutcome::Finished),
                    name: "<finished-serial-step>",
                })
            } else {
                match step.try_lock() {
                    None => Some(DispatchInfo {
                        result: Ok(StepOutcome::Contention),
                        name: "<contended-serial-step>",
                    }),
                    Some(mut guard) => {
                        let result = guard.try_run_erased(&mut ctx);
                        let name = guard.name();
                        if matches!(result, Ok(StepOutcome::Finished)) {
                            // Set the shared finished latch under the guard,
                            // before marking outputs drained, so a concurrent
                            // worker that observes the latch never re-runs the
                            // step nor re-marks its outputs.
                            drain.mark_finished();
                            if counter.observe_drain() {
                                guard.mark_outputs_drained(outputs_any);
                            }
                        }
                        Some(DispatchInfo { result, name })
                    }
                }
            }
        }
        WorkerStepEntry::Skip => None,
    };

    if let (Some(stats), Some(start)) = (stats, start) {
        let elapsed_ns = u64::try_from(start.elapsed().as_nanos()).unwrap_or(u64::MAX);
        // Wall ns at dispatch start, relative to pipeline start.
        let start_ns = stats.elapsed_ns().saturating_sub(elapsed_ns);
        // `None` (Skip) attempted no work, so there is nothing to record.
        if let Some(i) = info.as_ref() {
            match &i.result {
                Ok(outcome) => stats.record(step_idx, *outcome, start_ns, elapsed_ns),
                Err(_) => stats.record_error(step_idx, start_ns, elapsed_ns),
            }
        }
    }

    info
}

#[cfg(test)]
mod tests {
    use std::io;
    use std::sync::atomic::{AtomicUsize, Ordering};

    use super::*;
    use crate::erased::{ErasedStep, TypedStep};
    use crate::handles::BranchInputHandle;
    use crate::outputs::Single;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;
    use crate::runtime::contexts::build_chain_contexts;
    use crate::runtime::storage::DrainGate;
    use crate::step::{InputHandle, Step, StepCtx, StepKind, StepOutcome, StepProfile};
    use crate::topology::{BranchIdx, ChainGraph};
    use parking_lot::Mutex;

    #[test]
    fn run_worker_loop_exits_on_signal_done() {
        let signal = PipelineSignal::new();
        let mut entries: Vec<WorkerStepEntry> = vec![];
        let contexts =
            Arc::new(ChainContexts { inputs: vec![], outputs: vec![], bounded_queues: vec![] });
        let drain_counters: Vec<Arc<StepDrainCounter>> = vec![];
        let _ = ChainGraph::new();
        let mut worker = WorkerCore::new(0, None, None);

        signal.cancel();
        run_worker_loop(&mut worker, &mut entries, &contexts, &drain_counters, &signal, None);
        // If we reach this line, the loop exited cleanly.
    }

    // ── Test steps for dispatch-level coverage ──────────────────────────────

    /// `() → u32` source that returns `Finished` immediately.
    #[derive(Clone)]
    struct SrcFinished;
    impl Step for SrcFinished {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Src",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::Finished)
        }
    }

    /// `() → u32` source that returns `NoProgress` on its first `try_run` (so
    /// the sticky fast-path yields back to round-robin without removing it) and
    /// `Finished` on every later call (so it is removed during the round-robin
    /// pass, exercising `RoundRobinOutcome::removed_sticky_owner`).
    #[derive(Clone)]
    struct SrcIdleThenFinish {
        calls: Arc<AtomicUsize>,
    }
    impl Step for SrcIdleThenFinish {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "SrcIdleThenFinish",
                kind: StepKind::Exclusive,
                sticky: true,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let n = self.calls.fetch_add(1, Ordering::Relaxed);
            if n == 0 { Ok(StepOutcome::NoProgress) } else { Ok(StepOutcome::Finished) }
        }
    }

    /// `u32 → u32` step that always returns `Finished`. Used both as a
    /// `Parallel` body (counter-gated output close) and a `Serial` body
    /// (`DrainGate` short-circuit). The `runs` counter records every `try_run`
    /// so the short-circuit test can prove a worker did NOT re-run the step.
    #[derive(Clone)]
    struct FinishStep {
        kind: StepKind,
        runs: Arc<AtomicUsize>,
    }
    impl Step for FinishStep {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Finish",
                kind: self.kind,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            self.runs.fetch_add(1, Ordering::Relaxed);
            Ok(StepOutcome::Finished)
        }
        fn new_worker_copy(&self) -> Self {
            // Clones share the `runs` counter so the test can total runs across
            // every Parallel clone.
            self.clone()
        }
    }

    #[derive(Clone)]
    struct SinkStep;
    impl Step for SinkStep {
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
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            // Drain any inputs, then finish once the upstream edge is drained so
            // the worker loop can terminate (a real sink finishes on drain).
            while ctx.input.pop().is_some() {}
            if ctx.input.is_drained() {
                Ok(StepOutcome::Finished)
            } else {
                Ok(StepOutcome::NoProgress)
            }
        }
    }

    /// A non-sticky `Exclusive` sink that deliberately stays live for one extra
    /// round-robin pass: it ignores its input-drain status and finishes purely
    /// on an internal tick counter — `NoProgress` on the first `try_run`,
    /// `Finished` after. Keeping a second step alive for one more outer
    /// iteration *after* the sticky source is removed is what makes
    /// `sticky_owner_removed_via_round_robin_and_loop_exits` branch-specific: a
    /// plain `SinkStep` finishes in the same round-robin pass as the source
    /// (its input is already drained), emptying `live` so the loop exits via
    /// `live.is_empty()` even if the `removed_sticky_owner` branch had failed to
    /// clear `sticky_live`. Lingering forces the extra iteration on which a
    /// stale `sticky_live` would re-enter the sticky fast-path and re-invoke the
    /// already-removed source.
    struct LingerThenFinish {
        ticks: usize,
    }
    impl Step for LingerThenFinish {
        type Input = u32;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Linger",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            while ctx.input.pop().is_some() {}
            self.ticks += 1;
            // Stay live for exactly one extra round-robin pass before finishing,
            // regardless of input-drain status.
            if self.ticks >= 2 { Ok(StepOutcome::Finished) } else { Ok(StepOutcome::NoProgress) }
        }
    }

    /// Build `Src → Finish → Sink` (Finish having the given kind) and return the
    /// erased steps + the wired graph. The `Finish` step's `try_run` counter is
    /// returned so tests can assert how many times it actually ran.
    fn three_step_chain(
        finish_kind: StepKind,
    ) -> (Vec<Box<dyn ErasedStep>>, ChainGraph, Arc<AtomicUsize>) {
        let runs = Arc::new(AtomicUsize::new(0));
        let mut graph = ChainGraph::new();
        let src = graph.register_step("Src", 1);
        let mid = graph.register_step("Finish", 1);
        let sink = graph.register_step("Sink", 0);
        graph.wire(src, BranchIdx(0), mid);
        graph.wire(mid, BranchIdx(0), sink);
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(SrcFinished)),
            Box::new(TypedStep::new(FinishStep { kind: finish_kind, runs: Arc::clone(&runs) })),
            Box::new(TypedStep::new(SinkStep)),
        ];
        (steps, graph, runs)
    }

    /// A `Parallel` step's shared output queue must be closed exactly once — by
    /// the LAST clone to finish (the one that takes the `StepDrainCounter` to
    /// 0). Earlier finishers must leave the downstream input un-drained so a
    /// sibling could still push.
    #[test]
    fn parallel_last_finisher_closes_shared_output_exactly_once() {
        const N: usize = 4;
        let (steps, graph, _runs) = three_step_chain(StepKind::Parallel);
        let contexts = Arc::new(build_chain_contexts(&steps, &graph));
        let mid = StepIdx(1);
        let counter = StepDrainCounter::new(N);
        let signal = PipelineSignal::new();

        // One `Owned` clone per worker, each sharing `contexts.outputs[mid]`.
        let mut clones: Vec<WorkerStepEntry> =
            (0..N).map(|_| WorkerStepEntry::Owned { step: steps[mid.0].clone_boxed() }).collect();

        let sink_input = contexts.inputs[2].downcast_ref::<BranchInputHandle<u32>>().unwrap();

        for (i, clone) in clones.iter_mut().enumerate() {
            assert!(
                !InputHandle::is_drained(sink_input),
                "downstream input drained before the last clone finished (after {i} of {N})"
            );
            let info = dispatch_one_step(clone, mid, &contexts, &counter, &signal, None).unwrap();
            assert!(matches!(info.result, Ok(StepOutcome::Finished)));
        }
        assert!(
            InputHandle::is_drained(sink_input),
            "downstream input must be drained once the last Parallel clone finished"
        );
    }

    /// Once one worker finishes a `Serial` step (setting the shared `DrainGate`),
    /// a second `dispatch_one_step` short-circuits to a synthetic `Finished`
    /// WITHOUT re-acquiring the mutex or re-running the step.
    #[test]
    fn serial_drain_gate_short_circuits_second_worker() {
        let (steps, graph, runs) = three_step_chain(StepKind::Serial);
        let contexts = Arc::new(build_chain_contexts(&steps, &graph));
        let mid = StepIdx(1);
        let counter = StepDrainCounter::new(1);
        let signal = PipelineSignal::new();

        let shared = Arc::new(Mutex::new(steps.into_iter().nth(1).unwrap()));
        let drain = Arc::new(DrainGate::default());

        // Worker 1 finishes the step: runs once, sets the latch.
        let mut entry1 =
            WorkerStepEntry::Shared { step: Arc::clone(&shared), drain: Arc::clone(&drain) };
        let info1 =
            dispatch_one_step(&mut entry1, mid, &contexts, &counter, &signal, None).unwrap();
        assert!(matches!(info1.result, Ok(StepOutcome::Finished)));
        assert_eq!(runs.load(Ordering::Relaxed), 1, "step ran exactly once on the first worker");
        assert!(drain.is_finished(), "first finisher must set the DrainGate latch");

        // Worker 2 dispatches the same step: short-circuit, no re-run.
        let mut entry2 =
            WorkerStepEntry::Shared { step: Arc::clone(&shared), drain: Arc::clone(&drain) };
        let info2 =
            dispatch_one_step(&mut entry2, mid, &contexts, &counter, &signal, None).unwrap();
        assert!(matches!(info2.result, Ok(StepOutcome::Finished)));
        assert_eq!(
            info2.name, "<finished-serial-step>",
            "second worker takes the latch short-circuit"
        );
        assert_eq!(
            runs.load(Ordering::Relaxed),
            1,
            "the Serial step must NOT be re-run after the DrainGate latch is set"
        );
    }

    /// A sticky-owned source driven through `run_worker_loop` completes and the
    /// loop exits even though the cached `sticky_live` flag (not a per-iteration
    /// `live.contains` scan) gates the fast path. Exercises the S1b-006 cache:
    /// the sticky step is removed once it returns `Finished`, after which the
    /// fast path must be disabled and the loop must terminate.
    #[test]
    fn sticky_owner_completes_and_loop_exits() {
        let mut graph = ChainGraph::new();
        let src = graph.register_step("Src", 1);
        let sink = graph.register_step("Sink", 0);
        graph.wire(src, BranchIdx(0), sink);
        let steps: Vec<Box<dyn ErasedStep>> =
            vec![Box::new(TypedStep::new(SrcFinished)), Box::new(TypedStep::new(SinkStep))];
        let contexts = Arc::new(build_chain_contexts(&steps, &graph));

        // Single worker; the source (idx 0) is its Exclusive sticky owner, the
        // sink (idx 1) is Exclusive owned by the same worker for this 1-worker
        // run.
        let mut entries: Vec<WorkerStepEntry> = vec![
            WorkerStepEntry::Exclusive { step: steps.into_iter().next().unwrap() },
            WorkerStepEntry::Exclusive { step: Box::new(TypedStep::new(SinkStep)) },
        ];
        let drain_counters = vec![StepDrainCounter::new(1), StepDrainCounter::new(1)];
        let signal = PipelineSignal::new();
        let mut worker = WorkerCore::new(0, Some(src), Some(src));

        // The source finishes immediately; the sink then sees its input drained
        // and finishes too. `run_worker_loop` must return (no hang).
        run_worker_loop(&mut worker, &mut entries, &contexts, &drain_counters, &signal, None);
    }

    /// A sticky owner that returns `NoProgress` on its first call (yielding out
    /// of the sticky fast-path back to round-robin) and `Finished` later must be
    /// removed via the round-robin path (`RoundRobinOutcome::removed_sticky_owner`),
    /// after which the next outer iteration skips the sticky re-entry. This pins
    /// the round-robin removal branch (lines around `outcome.removed_sticky_owner`),
    /// not just the sticky fast-path removal exercised by
    /// `sticky_owner_completes_and_loop_exits`.
    ///
    /// A second `LingerThenFinish` step is kept alive for one extra round-robin
    /// pass *after* the source is removed, so the worker loop must run one more
    /// outer iteration. That iteration is where a stale `sticky_live` would
    /// wrongly re-enter the sticky fast-path and re-invoke the
    /// already-removed-from-`live` source — making the `== 2` source-call
    /// assertion below uniquely diagnostic of the `removed_sticky_owner` branch.
    /// (Without the linger, a plain sink would finish in the same pass as the
    /// source, emptying `live` so the loop exits via `live.is_empty()` whether
    /// or not `sticky_live` was cleared — and `== 2` would not be branch-specific.)
    #[test]
    fn sticky_owner_removed_via_round_robin_and_loop_exits() {
        let mut graph = ChainGraph::new();
        let src = graph.register_step("SrcIdleThenFinish", 1);
        let linger = graph.register_step("Linger", 0);
        graph.wire(src, BranchIdx(0), linger);

        let calls = Arc::new(AtomicUsize::new(0));
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(SrcIdleThenFinish { calls: Arc::clone(&calls) })),
            Box::new(TypedStep::new(LingerThenFinish { ticks: 0 })),
        ];
        let contexts = Arc::new(build_chain_contexts(&steps, &graph));

        let mut entries: Vec<WorkerStepEntry> = vec![
            WorkerStepEntry::Exclusive { step: steps.into_iter().next().unwrap() },
            WorkerStepEntry::Exclusive {
                step: Box::new(TypedStep::new(LingerThenFinish { ticks: 0 })),
            },
        ];
        let drain_counters = vec![StepDrainCounter::new(1), StepDrainCounter::new(1)];
        let signal = PipelineSignal::new();
        let mut worker = WorkerCore::new(0, Some(src), Some(src));

        // First sticky call → NoProgress (yield to round-robin); the source then
        // returns Finished during a round-robin pass, which must remove it and
        // disable the sticky fast-path so the loop terminates rather than hangs.
        // The `Linger` step stays alive for one more iteration, forcing the
        // post-removal outer iteration that exercises the cleared fast path.
        run_worker_loop(&mut worker, &mut entries, &contexts, &drain_counters, &signal, None);

        // The source must have been called EXACTLY twice: call 1 = idle in the
        // sticky fast-path (`NoProgress`, which does NOT remove it there — that
        // block only reaps a `Finished`), call 2 = `Finished` during the
        // round-robin pass. The lingering second step guarantees one more outer
        // iteration after that removal, so `== 2` is the branch-specific signal
        // for the `removed_sticky_owner` path: if that branch had failed to clear
        // `sticky_live`, the extra iteration's sticky fast-path would re-invoke
        // the (already-removed-from-`live`) source — the sticky block dispatches
        // `entries[owned_idx]` directly, not gated on `live` membership —
        // producing a third call. `== 2` therefore proves removal happened via
        // round-robin AND that it correctly disabled the fast path.
        assert_eq!(
            calls.load(Ordering::Relaxed),
            2,
            "source must be called exactly twice (sticky idle, then round-robin finish); \
             a different count means the removed_sticky_owner branch did not gate the fast path",
        );
    }
}
