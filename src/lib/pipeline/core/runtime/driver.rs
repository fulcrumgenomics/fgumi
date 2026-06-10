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

use crate::pipeline::core::erased::ErasedStepCtx;
use crate::pipeline::core::runtime::contexts::ChainContexts;
use crate::pipeline::core::runtime::drain::StepDrainCounter;
use crate::pipeline::core::runtime::live::LiveSteps;
use crate::pipeline::core::runtime::stats::PipelineStats;
use crate::pipeline::core::runtime::storage::WorkerStepEntry;
use crate::pipeline::core::runtime::worker_core::WorkerCore;
use crate::pipeline::core::signal::{PipelineError, PipelineSignal};
use crate::pipeline::core::step::StepOutcome;
use crate::pipeline::core::topology::StepIdx;

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
        if let Some(owned_idx) = worker.sticky_owner.filter(|idx| live.contains(*idx)) {
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
            }
        }

        // 2. Round-robin priority dispatch over all live steps.
        if !signal.is_done() {
            did_work |=
                round_robin_dispatch(entries, &mut live, contexts, drain_counters, signal, stats);
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

/// One pass of the round-robin dispatch over this worker's live steps, in
/// chain order. Returns `true` if any step did useful work (caller resets
/// backoff). Finished steps are removed from `live` at end-of-pass (deferred
/// so the in-progress walk over `live.order()` is not mutated underneath it).
fn round_robin_dispatch(
    entries: &mut [WorkerStepEntry],
    live: &mut LiveSteps,
    contexts: &Arc<ChainContexts>,
    drain_counters: &[Arc<StepDrainCounter>],
    signal: &Arc<PipelineSignal>,
    stats: Option<&Arc<PipelineStats>>,
) -> bool {
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
    for step_idx in finished {
        live.remove(step_idx);
    }
    did_work
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
    use super::*;
    use crate::pipeline::core::topology::ChainGraph;

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
}
