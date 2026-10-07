//! Worker loop body. Each worker thread runs `run_worker_loop` until
//! `signal.is_done()` or all steps are drained.
//!
//! Loop structure (per iteration):
//!   1. Check `signal.is_done()`; bail if true.
//!   2. **Sticky re-entry**: if this worker owns an `Exclusive` step that's
//!      sticky, drive it to a stop (`Progress`→loop; `Finished`→done;
//!      `NoProgress` or `Contention`→exit sticky), for at most
//!      `STICKY_BURST_LIMIT` consecutive calls. Sticky avoids context-
//!      switch overhead for source-style steps that emit in tight bursts.
//!   3. **Round-robin dispatch**: try each step in chain order. On
//!      `Progress`, restart from step 0 (priority) — except for this worker's
//!      sticky owner, which already had its burst in step 2. On `Finished`, the
//!      step is complete — `mark_outputs_drained` (counter-gated for `Parallel`)
//!      and remove it from the worklist. `NoProgress`/`Contention` are idle ticks.
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
use crate::liveness::LivenessCounter;
use crate::runtime::contexts::ChainContexts;
use crate::runtime::drain::StepDrainCounter;
use crate::runtime::event_count::{NotifyOutcome, PoolEventCount};
use crate::runtime::live::LiveSteps;
use crate::runtime::scheduler::{Scheduler, WalkDirection};
use crate::runtime::stats::{PipelineStats, WakeCounts};
use crate::runtime::storage::WorkerStepEntry;
use crate::runtime::worker_core::{WorkerCore, WorkerRole};
use crate::runtime::worker_state::WorkerStateBoard;
use crate::signal::{PipelineError, PipelineSignal};
use crate::step::StepOutcome;
use crate::topology::StepIdx;

/// Maximum consecutive sticky re-entries before the worker drops back to a
/// round-robin pass.
///
/// `StepOutcome::Progress` means the step "pushed **or held** an item" (see the
/// [`crate::step`] completion contract), so a sticky step whose output is full
/// reports `Progress` on every call while moving nothing. Re-entering without a
/// bound would spin on that step forever and never dispatch the downstream step
/// that drains the output — at one worker, a hang; at several, a burnt core.
///
/// The bound is high enough that the fast path keeps its point (a source
/// emitting a tight burst skips the outer loop's bookkeeping per item) and low
/// enough that a step stuck on a full output yields promptly. It caps starvation,
/// it is not a tuning knob: forward progress comes from this bound *together
/// with* round-robin declining to restart the walk on the sticky owner's
/// `Progress`.
const STICKY_BURST_LIMIT: usize = 1024;

/// Per-thread diagnostic state carried across loop iterations.
///
/// A driver thread hosting a `Shared` group parks once for the whole group, so
/// its idle must be booked to *one* step. Booking it to the group's first step
/// made a merge's waiting show up as idle of the group's reader step while the
/// merge reported "0 parks". The honest answer for a shared thread is the step
/// that last made progress on it (exact when one step is live; the step the
/// thread is actually driving otherwise), falling back to the primary before
/// any progress. A step that has finished waits for nothing, so its group's
/// later idle goes to a step that is still live. Pool workers do not attribute
/// idle per step, so `note` is a single not-taken branch there.
pub(crate) struct LoopDiag {
    last_progress: Option<StepIdx>,
    track: bool,
}

impl LoopDiag {
    /// Diagnostic state for one thread; tracks progress only on a driver thread.
    pub(crate) fn new(is_driver: bool) -> Self {
        Self { last_progress: None, track: is_driver }
    }

    /// Record `step`'s dispatch outcome: a `Progress`/`Finished` on a driver
    /// thread makes `step` the one this thread's next idle is booked to.
    #[inline]
    fn note(&mut self, step: StepIdx, result: &std::io::Result<StepOutcome>) {
        if self.track && matches!(result, Ok(StepOutcome::Progress | StepOutcome::Finished)) {
            self.last_progress = Some(step);
        }
    }

    /// The step to book this thread's idle to: the step that last made
    /// progress, else the group primary — whichever is still live, since a
    /// finished step no longer waits for anything; else the first live step.
    #[inline]
    fn attributed(&self, primary: StepIdx, live: &LiveSteps) -> StepIdx {
        let live_step = |s: StepIdx| live.contains(s).then_some(s);
        self.last_progress
            .and_then(live_step)
            .or_else(|| live_step(primary))
            .or_else(|| live.order().first().copied())
            .unwrap_or(primary)
    }
}

/// Run the worker loop for one worker thread.
///
/// `entries[step_idx] = WorkerStepEntry` — this worker's storage.
/// `contexts` — shared per-step input/output handles.
/// `drain_counters[step_idx]` — the per-step `StepDrainCounter` that gates the
/// output close on `Finished` (init N for Parallel so only the last clone
/// closes the shared output; init 1 for Serial/Exclusive).
/// `signal` — error/cancel broadcast.
/// `board` — optional per-thread state board for scheduling telemetry; `None`
/// keeps the loop's stamping a no-op (telemetry-off parity).
/// `state_slot` — this thread's slot in `board` (ignored when `board` is `None`).
/// `parker` — the pool's event-count. `Some` only for pool workers when
/// `n_threads > 1`; drives both the idle-branch park (instead of `sleep_backoff`)
/// and the on-`Progress`/`Finished` notify that wakes parked peers. `None`
/// (single-thread paths, and pinned workers — see the idle branch) keeps the
/// existing `sleep_backoff` behaviour.
#[allow(clippy::too_many_arguments)] // per-step shared state plus the liveness shard; a struct would only rename it
#[allow(clippy::too_many_lines)] // the whole-pass sticky + round-robin + backoff discipline is documented inline; splitting it would scatter the invariant across functions
pub fn run_worker_loop(
    worker: &mut WorkerCore,
    entries: &mut [WorkerStepEntry],
    contexts: &Arc<ChainContexts>,
    drain_counters: &[Arc<StepDrainCounter>],
    signal: &Arc<PipelineSignal>,
    stats: Option<&Arc<PipelineStats>>,
    liveness: &LivenessCounter,
    scheduler: &dyn Scheduler,
    board: Option<&WorkerStateBoard>,
    state_slot: usize,
    parker: Option<&PoolEventCount>,
    pinned: bool,
) {
    // A pinned worker (the sole eligible dispatcher of an Exclusive/sticky step,
    // or the affinity target of a Serial step) must NOT deep-park on the shared
    // event-count: `notify_one` cannot target a specific worker, so a push meant
    // to wake *this* worker's pinned step could wake a peer that `Skip`s it,
    // leaving the pinned step's owner asleep to its timeout. Pinned workers keep
    // `sleep_backoff`. Its step is also often driven by out-of-engine events
    // (a reader's prefetch thread, the affinity-pinned grouper) that no push
    // notifies anyway. `parker` is threaded to `None` for the idle branch here,
    // while still being used for the on-Progress notify (a pinned worker's
    // pushes must still wake parked peers).
    let idle_parker = if pinned { None } else { parker };
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

    // A driver thread attributes busy time PER grouped step on the off-pool
    // detached line (so a multi-step `Shared` group shows each step's real
    // busy, not the whole thread's time under one name), recorded inside
    // `dispatch_one_step`; a pool worker records aggregate busy by `thread_id`
    // below. Idle/park stay thread-level (keyed to the group's primary step).
    let is_driver = matches!(worker.role(), WorkerRole::Driver { .. });
    let mut diag = LoopDiag::new(is_driver);

    loop {
        if signal.is_done() {
            break;
        }
        // Exit when this worker has nothing left to dispatch.
        if live.is_empty() {
            break;
        }

        let mut did_work = false;

        // Time the dispatch (busy) section vs the backoff sleep (idle) below,
        // per worker, only when stats are on (`Instant::now()` is otherwise not
        // called — the no-stats loop stays zero-cost).
        let work_start = stats.map(|_| Instant::now());

        // 1. Sticky re-entry for sticky-owned steps (either an Exclusive
        // sticky step this worker owns, or a Serial+sticky step whose
        // Affinity targets this worker). The Pipeline::run path only
        // sets `sticky_owner` when the step is actually sticky, so no
        // per-iteration profile peek is needed. Re-enter while the
        // step makes Progress; exit on Finished / NoProgress / Contention
        // / Err. Remove from the worklist on Finished (source drain) or
        // observed input drain (mid-step drain). Gated on the step still
        // being live — once removed, the sticky fast-path is disabled.
        //
        // Bounded at `STICKY_BURST_LIMIT` calls: `Progress` also covers "held an
        // item", so a step whose output is full reports it indefinitely without
        // moving anything. An unbounded re-entry would then never reach the
        // round-robin pass that dispatches the downstream step draining that
        // output.
        if let Some(owned_idx) = worker.sticky_owner.filter(|_| sticky_live) {
            let mut mark_skip = false;
            for _ in 0..STICKY_BURST_LIMIT {
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
                    liveness,
                    worker.thread_id,
                    is_driver,
                    board,
                    state_slot,
                    parker,
                    &mut diag,
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
                    Ok(StepOutcome::NoProgress | StepOutcome::Contention | StepOutcome::Capped) => {
                        break;
                    }
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
                liveness,
                worker.thread_id,
                scheduler.walk(),
                is_driver,
                board,
                state_slot,
                parker,
                &mut diag,
            );
            did_work |= outcome.did_work;
            if outcome.removed_sticky_owner {
                // The sticky owner finished during round-robin; disable the
                // fast path so subsequent iterations skip the sticky block.
                sticky_live = false;
            }
        }

        // Attribute the dispatch section's wall time to this thread's busy total.
        // Pool workers sum the whole pass into the N-worker utilisation line (by
        // thread_id). Driver threads instead record each grouped step's own busy
        // inside `dispatch_one_step` (on the off-pool detached line, by step) so a
        // multi-step `Shared` group isn't collapsed onto one name — so nothing to
        // record here for a driver.
        if let (Some(stats), Some(ws)) = (stats, work_start)
            && let WorkerRole::Pool = worker.role()
        {
            let ns = u64::try_from(ws.elapsed().as_nanos()).unwrap_or(u64::MAX);
            stats.record_worker_busy(worker.thread_id, ns);
        }

        // 3. Exponential-backoff sleep on no-progress; reset on progress. The
        // sleep is the worker's idle/blocked time — attribute it per worker so
        // pool under-utilisation (cores parked while one worker drives a Serial
        // step) is visible in `--pipeline-stats`.
        if did_work {
            worker.reset_backoff();
        } else if signal.is_done() {
            break;
        } else if let Some(ec) = idle_parker {
            // Event-count park (Layer 2): block instead of spin-poll when the
            // whole pass found no work. The two-phase protocol closes the
            // lost-wakeup race: arm (register + fence) → re-poll the real
            // condition → block only if still empty.
            //
            // Arm and re-poll *before* stamping `Parked` or starting the idle
            // timer: a productive re-poll (work appeared in the arm→wait window)
            // must not be recorded as idle, nor leave the board reading `Parked`
            // while the re-poll's `dispatch_one_step` runs real work.

            // Phase 1: arm. The fence in `prepare_wait` orders the waiter
            // registration before the re-poll below, so a producer that
            // publishes work after this point either sees us as a waiter (and
            // bumps the generation) or we see its item on the re-poll.
            let key = ec.prepare_wait();

            // Phase 2: re-poll the *real* condition — one full dispatch pass.
            // This is one pass per park episode (amortised over the block), not
            // per iteration, so it is not the continuous polling Layer 2 removes.
            let recheck = round_robin_dispatch(
                entries,
                &mut live,
                worker.sticky_owner,
                contexts,
                drain_counters,
                signal,
                stats,
                liveness,
                worker.thread_id,
                scheduler.walk(),
                is_driver,
                board,
                state_slot,
                parker,
                &mut diag,
            );
            if recheck.removed_sticky_owner {
                sticky_live = false;
            }

            if recheck.did_work || signal.is_done() {
                // Work appeared in the arm→wait window (or we are shutting
                // down): do not block, and record no idle — the re-poll was
                // productive. Balance the arm and re-ramp fresh.
                ec.cancel_wait(key);
                worker.reset_backoff();
            } else {
                // No progress remains: only now do we actually park. Stamp
                // `Parked` and start the idle timer immediately before blocking,
                // so the idle attribution covers exactly the wait, not the
                // preceding (possibly productive) re-poll.
                if let Some(b) = board {
                    b.stamp(state_slot, crate::runtime::worker_state::WorkerState::Parked, None);
                }
                let sleep_start = stats.map(|_| Instant::now());

                // Phase 3: block until the generation moves (a peer's notify),
                // the deadline elapses (self-heal), or a spurious wake. Keep the
                // exponential backoff as the wait *timeout* so a missed wakeup
                // cannot hang longer than the cap, and the deadlock monitor
                // still ticks. Re-poll on return regardless of outcome.
                let deadline = worker.backoff_deadline();
                let outcome = ec.wait(key, deadline);
                worker.increase_backoff();

                if let (Some(stats), Some(ss)) = (stats, sleep_start) {
                    let ns = u64::try_from(ss.elapsed().as_nanos()).unwrap_or(u64::MAX);
                    stats.record_ec_wait(worker.thread_id, outcome);
                    match worker.role() {
                        WorkerRole::Pool => stats.record_worker_idle(worker.thread_id, ns),
                        WorkerRole::Driver { primary_step } => {
                            let step = diag.attributed(primary_step, &live);
                            stats.record_detached_idle(step, ns);
                            stats.record_detached_park(step);
                        }
                    }
                }
            }
        } else {
            // No parker (single-thread path or a pinned worker): keep the
            // original exponential-backoff sleep.
            if let Some(b) = board {
                b.stamp(state_slot, crate::runtime::worker_state::WorkerState::Parked, None);
            }
            let deadline = worker.backoff_deadline();
            let sleep_start = stats.map(|_| Instant::now());
            worker.sleep_backoff();
            worker.increase_backoff();
            if let (Some(stats), Some(ss)) = (stats, sleep_start) {
                let elapsed = ss.elapsed();
                let ns = u64::try_from(elapsed.as_nanos()).unwrap_or(u64::MAX);
                // A park that returned before its deadline was unparked (or woke
                // spuriously); one that ran to it was fed by the timer.
                let timed_out = elapsed >= deadline;
                match worker.role() {
                    WorkerRole::Pool => {
                        stats.record_worker_idle(worker.thread_id, ns);
                        stats.record_timer_park(worker.thread_id, timed_out);
                    }
                    WorkerRole::Driver { primary_step } => {
                        let step = diag.attributed(primary_step, &live);
                        stats.record_detached_idle(step, ns);
                        stats.record_detached_park(step);
                        stats.record_detached_park_outcome(step, timed_out);
                    }
                }
            }
        }
    }

    // The loop exits via several `break`s (signal done at the top, empty
    // worklist, or the no-progress branch), and only the no-progress branch
    // stamps `Parked` before sleeping — the others leave this slot at whatever
    // the final `dispatch_one_step` wrote, which is `Running`. The board has no
    // terminal state, so stamp `Parked` once on exit; otherwise the sampler
    // keeps counting an exited worker's slot as a serviced (busy) `Running`.
    if let Some(b) = board {
        b.stamp(state_slot, crate::runtime::worker_state::WorkerState::Parked, None);
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

/// For [`WalkDirection::RefillThenReverse`], how many of the leading live
/// steps (`order`, in chain order) are at or before the refill step and so walk
/// forward; `0` for the other directions.
fn refill_split(walk: WalkDirection, order: &[StepIdx]) -> usize {
    match walk {
        WalkDirection::RefillThenReverse(through) => order.partition_point(|s| s.0 <= through.0),
        WalkDirection::Forward | WalkDirection::Reverse => 0,
    }
}

/// The live-step position the `i`-th visit of a pass lands on, for `n` live
/// steps of which the first `split` walk forward (see [`refill_split`]).
fn walk_position(walk: WalkDirection, i: usize, n: usize, split: usize) -> usize {
    match walk {
        WalkDirection::Forward => i,
        WalkDirection::Reverse => n - 1 - i,
        WalkDirection::RefillThenReverse(_) if i < split => i,
        WalkDirection::RefillThenReverse(_) => n - 1 - (i - split),
    }
}

/// One pass of the round-robin dispatch over this worker's live steps, in
/// chain order. Finished steps are removed from `live` at end-of-pass (deferred
/// so the in-progress walk over `live.order()` is not mutated underneath it).
/// `sticky_owner` (if any) has two roles here: it is reported back via
/// [`RoundRobinOutcome::removed_sticky_owner`] when it finishes in this pass, so
/// the caller can disable the per-iteration sticky fast-path without a linear
/// `live.contains` scan; and its `Progress` does **not** trigger the priority
/// restart, so the walk continues to the steps downstream of it (see the
/// `Progress` arm).
#[allow(clippy::too_many_arguments)] // shared per-step state + the walk policy; a struct would not clarify
fn round_robin_dispatch(
    entries: &mut [WorkerStepEntry],
    live: &mut LiveSteps,
    sticky_owner: Option<StepIdx>,
    contexts: &Arc<ChainContexts>,
    drain_counters: &[Arc<StepDrainCounter>],
    signal: &Arc<PipelineSignal>,
    stats: Option<&Arc<PipelineStats>>,
    liveness: &LivenessCounter,
    worker_slot: usize,
    walk: WalkDirection,
    is_driver: bool,
    board: Option<&WorkerStateBoard>,
    state_slot: usize,
    parker: Option<&PoolEventCount>,
    diag: &mut LoopDiag,
) -> RoundRobinOutcome {
    let mut did_work = false;
    // Steps that finished this pass, removed from `live` after the walk. A
    // step is visited at most once per pass (the cursor only advances; we
    // `break` on `Progress`/error, never revisit), so deferring removal is
    // safe and avoids reorder-under-iteration.
    let mut finished: Vec<StepIdx> = Vec::new();
    let n = live.len();
    let split = refill_split(walk, live.order());
    for i in 0..n {
        if signal.is_done() {
            break;
        }
        // The Scheduler selects the walk DIRECTION over this worker's live
        // steps: `Forward` = chain order (upstream-first, favour production);
        // `Reverse` = downstream-first (favour draining buffered work before
        // producing more). Everything else — skip-on-contention, the sticky
        // source/sink fast-path, Finished handling — is direction-agnostic.
        let pos = walk_position(walk, i, n, split);
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
            liveness,
            worker_slot,
            is_driver,
            board,
            state_slot,
            parker,
            diag,
        ) else {
            continue; // Skip (build-time placeholder; should not appear in `live`)
        };
        match info.result {
            Ok(StepOutcome::Progress) => {
                did_work = true;
                // Restart the walk at the top so the highest-priority step runs
                // again — EXCEPT for this worker's sticky owner. That step just
                // had its dedicated burst in phase 1 of the worker loop, so a
                // priority restart here only re-runs it. Worse, a sticky step
                // that reports `Progress` while sitting on a full output (against
                // the `StepOutcome::Progress` contract) would report it forever:
                // breaking the pass on that outcome means the downstream step
                // that would drain the output is never reached and the pair
                // livelocks at one worker. Walking past the sticky owner is what
                // turns the bounded burst into real forward progress.
                //
                // A non-sticky step gets the restart, and that is safe under
                // every walk (including the forward refill prefix) because of
                // the held-slot convention its steps follow: a *new* hold reports
                // `Progress` once, while a retry that is still rejected reports
                // `NoProgress`/`Contention`. The restart therefore costs one
                // extra visit, and that visit's idle outcome lets the walk reach
                // the consumer that drains the output
                // (`refill_walk_reaches_the_consumer_of_a_holding_step`).
                restart_priority = sticky_owner != Some(step_idx);
            }
            // Nothing to do this call. The step terminates by returning
            // `Finished` (handled below); `NoProgress`/`Contention` are idle
            // ticks — there is no separate drain protocol.
            Ok(StepOutcome::NoProgress | StepOutcome::Contention | StepOutcome::Capped) => {}
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
#[allow(clippy::too_many_arguments)] // per-step shared state plus the liveness shard; a struct would only rename it
fn dispatch_one_step(
    entry: &mut WorkerStepEntry,
    step_idx: StepIdx,
    contexts: &ChainContexts,
    counter: &StepDrainCounter,
    signal: &Arc<PipelineSignal>,
    stats: Option<&Arc<PipelineStats>>,
    // Always-on liveness signal for the deadlock monitor, sharded per worker so
    // the bump is normally an uncontended increment (a dedicated driver reuses
    // `worker_slot` 0 and so shares slot 0 with pool worker 0 — the atomic
    // `fetch_add` counts a coincident bump correctly; see `crate::liveness`).
    // Separate from `stats` on purpose:
    // liveness must be free enough to leave armed, while `stats` pays for
    // per-dispatch timing and stays opt-in. See `crate::liveness`.
    liveness: &LivenessCounter,
    worker_slot: usize,
    // When true (a dedicated driver thread), attribute this dispatch's wall time
    // to the off-pool detached line keyed by `step_idx` — so each grouped step
    // reports its own busy. Pool workers pass `false` and record aggregate busy
    // by `thread_id` in the loop instead.
    is_driver: bool,
    // Optional per-thread state board for scheduling telemetry; `None` keeps
    // stamping a no-op (telemetry-off parity). `state_slot` is this thread's
    // slot in `board` (ignored when `board` is `None`).
    board: Option<&WorkerStateBoard>,
    state_slot: usize,
    // The pool's event-count. On a `Progress` outcome (the step pushed or held
    // an item — including a reorder must-accept stash) wake one parked peer; on
    // `Finished` (an output edge just closed) wake all. `None` for single-thread
    // paths. This is the notify seam: `Progress`/`Finished` are exactly the
    // outcomes that can give a parked consumer new work, and using the outcome
    // (rather than the `BranchOutputHandle::push` site) captures the stash-then-
    // return-`Ok` reorder case without threading the parker through every queue.
    parker: Option<&PoolEventCount>,
    // Per-thread diagnostic state; updated after the outcome is known.
    diag: &mut LoopDiag,
) -> Option<DispatchInfo> {
    let outputs_any: &(dyn Any + Send + Sync) = contexts.outputs[step_idx.0].as_ref();
    let mut ctx = ErasedStepCtx {
        input: contexts.inputs[step_idx.0].as_ref(),
        outputs: outputs_any,
        signal,
        counters: &contexts.step_counters[step_idx.0],
    };

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
    if let Some(b) = board {
        b.stamp(state_slot, crate::runtime::worker_state::WorkerState::Running, Some(step_idx));
    }

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
            } else if step.is_locked() {
                // Test-and-test-and-set: probe with a shared-mode load before
                // attempting the lock. `parking_lot::Mutex::try_lock` is a
                // `compare_exchange` on the lock word that writes (and so
                // invalidates) the holder's cache line *even when it fails*. A
                // surplus worker that `try_lock`s a `Serial + Affinity::None`
                // step held by the productive worker steals that line on every
                // idle pass, and that step is often the throughput ceiling
                // (e.g. `GroupByQueryname` in `correct`). `is_locked()` is a
                // plain relaxed load — no line-ownership transfer — so probing
                // it first turns the common "held by a peer" case into a
                // read-only `Contention` with no write to the holder's line.
                Some(DispatchInfo {
                    result: Ok(StepOutcome::Contention),
                    name: "<contended-serial-step>",
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

    if let Some(i) = info.as_ref() {
        diag.note(step_idx, &i.result);
    }

    // Liveness first, and unconditionally: this is what the deadlock monitor
    // samples, so it must not depend on `stats` being attached. Only productive
    // outcomes count — a wedged pipeline still spins through `NoProgress` and
    // `Contention` dispatches forever, so counting those would make a wedge look
    // alive and defeat the whole detector.
    if let Some(i) = info.as_ref()
        && matches!(i.result, Ok(StepOutcome::Progress | StepOutcome::Finished))
    {
        liveness.bump(worker_slot);
        // Notify seam: a productive dispatch may have given a parked consumer
        // new work. `Progress` = "pushed or held an item" (incl. a reorder
        // must-accept stash that unblocks an ordinal waiter) → wake one parked
        // peer, which cascades (its own downstream push wakes the next).
        // `Finished` closed an output edge (a parked consumer must observe
        // end-of-stream, and a sink's `Finished` closes no edge but its peers
        // may be waiting on the drain latch) → wake all. When nobody is parked,
        // both are a single relaxed load (see `PoolEventCount`).
        if let Some(ec) = parker {
            let outcome = match i.result {
                Ok(StepOutcome::Finished) => {
                    let _ = ec.notify_all();
                    None
                }
                _ => Some(ec.notify_one()),
            };
            if let (Some(stats), Some(o)) = (stats, outcome) {
                let counts = match o {
                    NotifyOutcome::NoWaiters => {
                        WakeCounts { no_waiters: 1, ..WakeCounts::default() }
                    }
                    NotifyOutcome::Suppressed => {
                        WakeCounts { suppressed: 1, ..WakeCounts::default() }
                    }
                    NotifyOutcome::Woken => WakeCounts { notified: 1, ..WakeCounts::default() },
                };
                stats.record_wake(step_idx, counts);
            }
        }
    }

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
            // On a driver thread, this step's try_run wall is its own off-pool
            // busy (excluded from the pool%); each grouped step accrues its own.
            if is_driver {
                stats.record_detached_busy(step_idx, elapsed_ns);
            }
        }
    }

    info
}

#[cfg(test)]
mod tests {
    use std::io;
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

    use super::*;

    /// The visit order for each walk over five live steps `[0, 2, 3, 5, 7]`
    /// (chain order, with gaps where steps finished).
    #[rstest::rstest]
    #[case::forward(WalkDirection::Forward, &[0, 2, 3, 5, 7])]
    #[case::reverse(WalkDirection::Reverse, &[7, 5, 3, 2, 0])]
    #[case::refill_through_live_step(
        WalkDirection::RefillThenReverse(StepIdx(3)),
        &[0, 2, 3, 7, 5]
    )]
    #[case::refill_through_finished_step(
        WalkDirection::RefillThenReverse(StepIdx(4)),
        &[0, 2, 3, 7, 5]
    )]
    #[case::refill_through_first(WalkDirection::RefillThenReverse(StepIdx(0)), &[0, 7, 5, 3, 2])]
    #[case::refill_through_last(WalkDirection::RefillThenReverse(StepIdx(7)), &[0, 2, 3, 5, 7])]
    #[case::refill_before_all(WalkDirection::RefillThenReverse(StepIdx(9)), &[0, 2, 3, 5, 7])]
    fn walk_visits_live_steps_in_order(#[case] walk: WalkDirection, #[case] expected: &[usize]) {
        let order: Vec<StepIdx> = [0, 2, 3, 5, 7].into_iter().map(StepIdx).collect();
        let split = refill_split(walk, &order);
        let visited: Vec<usize> =
            (0..order.len()).map(|i| order[walk_position(walk, i, order.len(), split)].0).collect();
        assert_eq!(visited, expected);
    }
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
        let contexts = Arc::new(ChainContexts {
            inputs: vec![],
            outputs: vec![],
            bounded_queues: vec![],
            edges: vec![],
            step_counters: vec![],
        });
        let drain_counters: Vec<Arc<StepDrainCounter>> = vec![];
        let _ = ChainGraph::new();
        let mut worker = WorkerCore::new(0, None, None);

        signal.cancel();
        run_worker_loop(
            &mut worker,
            &mut entries,
            &contexts,
            &drain_counters,
            &signal,
            None,
            &crate::liveness::LivenessCounter::new(1),
            &crate::runtime::scheduler::ChainOrderScheduler,
            None,
            0,
            None,
            false,
        );
        // If we reach this line, the loop exited cleanly.
    }

    #[test]
    fn run_worker_loop_stamps_parked_on_exit() {
        // A worker whose slot was last stamped `Running` (mid-dispatch) must be
        // re-stamped `Parked` once `run_worker_loop` returns; otherwise the
        // sampler keeps counting the exited worker as a serviced (busy)
        // `Running` slot. The loop here exits immediately via the top-of-loop
        // `signal.is_done()` break, which never reaches the no-progress `Parked`
        // stamp — so only the on-exit stamp can flip Running back to Parked.
        use crate::runtime::worker_state::{WorkerState, WorkerStateBoard};
        let signal = PipelineSignal::new();
        let mut entries: Vec<WorkerStepEntry> = vec![];
        let contexts = Arc::new(ChainContexts {
            inputs: vec![],
            outputs: vec![],
            bounded_queues: vec![],
            edges: vec![],
            step_counters: vec![],
        });
        let drain_counters: Vec<Arc<StepDrainCounter>> = vec![];
        let mut worker = WorkerCore::new(0, None, None);

        let board = WorkerStateBoard::new(1);
        board.stamp(0, WorkerState::Running, Some(StepIdx(0)));

        signal.cancel();
        run_worker_loop(
            &mut worker,
            &mut entries,
            &contexts,
            &drain_counters,
            &signal,
            None,
            &crate::liveness::LivenessCounter::new(1),
            &crate::runtime::scheduler::ChainOrderScheduler,
            Some(&board),
            0,
            None,
            false,
        );
        assert_eq!(board.read(0).0, WorkerState::Parked);
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
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
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
            let info = dispatch_one_step(
                clone,
                mid,
                &contexts,
                &counter,
                &signal,
                None,
                &LivenessCounter::new(1),
                0,
                false,
                None,
                0,
                None,
                &mut LoopDiag::new(false),
            )
            .unwrap();
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
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let mid = StepIdx(1);
        let counter = StepDrainCounter::new(1);
        let signal = PipelineSignal::new();

        let shared = Arc::new(Mutex::new(steps.into_iter().nth(1).unwrap()));
        let drain = Arc::new(DrainGate::default());

        // Worker 1 finishes the step: runs once, sets the latch.
        let mut entry1 =
            WorkerStepEntry::Shared { step: Arc::clone(&shared), drain: Arc::clone(&drain) };
        let info1 = dispatch_one_step(
            &mut entry1,
            mid,
            &contexts,
            &counter,
            &signal,
            None,
            &LivenessCounter::new(1),
            0,
            false,
            None,
            0,
            None,
            &mut LoopDiag::new(false),
        )
        .unwrap();
        assert!(matches!(info1.result, Ok(StepOutcome::Finished)));
        assert_eq!(runs.load(Ordering::Relaxed), 1, "step ran exactly once on the first worker");
        assert!(drain.is_finished(), "first finisher must set the DrainGate latch");

        // Worker 2 dispatches the same step: short-circuit, no re-run.
        let mut entry2 =
            WorkerStepEntry::Shared { step: Arc::clone(&shared), drain: Arc::clone(&drain) };
        let info2 = dispatch_one_step(
            &mut entry2,
            mid,
            &contexts,
            &counter,
            &signal,
            None,
            &LivenessCounter::new(1),
            0,
            false,
            None,
            0,
            None,
            &mut LoopDiag::new(false),
        )
        .unwrap();
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

    /// Layer 0 (test-and-test-and-set): when a `Serial` step's mutex is already
    /// held (by the productive worker), a surplus worker's dispatch reports
    /// `Contention` via the `is_locked()` probe and does NOT run the step. The
    /// probe is a shared-mode load, so it never takes the lock word's write
    /// path (`try_lock`'s `compare_exchange`) that would invalidate the
    /// holder's cache line. We prove the observable contract: the step body is
    /// not entered while the lock is held.
    #[test]
    fn serial_contention_probe_does_not_run_held_step() {
        let (steps, graph, runs) = three_step_chain(StepKind::Serial);
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let mid = StepIdx(1);
        let counter = StepDrainCounter::new(1);
        let signal = PipelineSignal::new();

        let shared = Arc::new(Mutex::new(steps.into_iter().nth(1).unwrap()));
        let drain = Arc::new(DrainGate::default());

        // Simulate the productive worker holding the step's mutex.
        let held = Arc::clone(&shared);
        let guard = held.lock();
        assert!(shared.is_locked(), "precondition: the Serial mutex is held");

        // A surplus worker dispatches the same step while it is held.
        let mut entry =
            WorkerStepEntry::Shared { step: Arc::clone(&shared), drain: Arc::clone(&drain) };
        let info = dispatch_one_step(
            &mut entry,
            mid,
            &contexts,
            &counter,
            &signal,
            None,
            &LivenessCounter::new(1),
            0,
            false,
            None,
            0,
            None,
            &mut LoopDiag::new(false),
        )
        .unwrap();

        assert!(
            matches!(info.result, Ok(StepOutcome::Contention)),
            "a held Serial step must report Contention to the surplus worker"
        );
        assert_eq!(info.name, "<contended-serial-step>", "the contended-step name is reported");
        assert_eq!(
            runs.load(Ordering::Relaxed),
            0,
            "the Serial step body must NOT run while another worker holds its mutex"
        );

        drop(guard);
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
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));

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
        run_worker_loop(
            &mut worker,
            &mut entries,
            &contexts,
            &drain_counters,
            &signal,
            None,
            &crate::liveness::LivenessCounter::new(1),
            &crate::runtime::scheduler::ChainOrderScheduler,
            None,
            0,
            None,
            false,
        );
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
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));

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
        run_worker_loop(
            &mut worker,
            &mut entries,
            &contexts,
            &drain_counters,
            &signal,
            None,
            &crate::liveness::LivenessCounter::new(1),
            &crate::runtime::scheduler::ChainOrderScheduler,
            None,
            0,
            None,
            false,
        );

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

    // ── G1: the load-bearing driver invariant ───────────────────────────────
    //
    // A driver thread (`WorkerCore::driver`) is just `run_worker_loop` over a
    // few `Owned` steps. Its no-deadlock property rests ENTIRELY on the loop
    // trying EVERY live step in a pass before it parks. A naive "drive the first
    // live step until it yields, then park" loop would wedge the coordination
    // driver: it parks on a step whose input isn't ready yet while a *sibling*
    // step on the same driver holds the work that would unblock it (e.g. park on
    // `FindBoundariesAndSort` while `SpillGather` holds the chunk that frees the
    // capacity-1 arena). These two steps pin that the whole-pass discipline holds.

    /// `Owned` step wedged on `NoProgress` until `gate` is flipped by a sibling,
    /// then `Finished`. Placed FIRST in the walk so a park-on-first-NoProgress
    /// loop would never let the sibling run — the gate never flips — and hang.
    #[derive(Clone)]
    struct WedgedUntilGate {
        gate: Arc<AtomicBool>,
    }
    impl Step for WedgedUntilGate {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "WedgedUntilGate",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if self.gate.load(Ordering::Acquire) {
                Ok(StepOutcome::Finished)
            } else {
                Ok(StepOutcome::NoProgress)
            }
        }
    }

    /// `Owned` sibling that flips `gate` and finishes on its first dispatch —
    /// reached only if the loop tries all live steps in a pass rather than
    /// parking on the wedged step's `NoProgress`.
    #[derive(Clone)]
    struct GateOpener {
        gate: Arc<AtomicBool>,
    }
    impl Step for GateOpener {
        type Input = u32;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "GateOpener",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            self.gate.store(true, Ordering::Release);
            Ok(StepOutcome::Finished)
        }
    }

    /// A driver (`WorkerCore::driver`, Park backoff) driving `[Wedged, Opener]`
    /// as two `Owned` steps must drive the sibling that unblocks the wedged step
    /// and terminate. A watchdog aborts the process on a wedge so the failure is
    /// loud rather than a silent hang.
    #[test]
    fn driver_round_robins_all_live_before_parking() {
        let mut graph = ChainGraph::new();
        let wedged = graph.register_step("WedgedUntilGate", 1);
        let opener = graph.register_step("GateOpener", 0);
        graph.wire(wedged, BranchIdx(0), opener);

        let gate = Arc::new(AtomicBool::new(false));
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(WedgedUntilGate { gate: Arc::clone(&gate) })),
            Box::new(TypedStep::new(GateOpener { gate: Arc::clone(&gate) })),
        ];
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));

        // Hand-built driver row: both steps Owned on the one driver thread.
        let mut entries: Vec<WorkerStepEntry> =
            steps.into_iter().map(|step| WorkerStepEntry::Owned { step }).collect();
        let drain_counters = vec![StepDrainCounter::new(1), StepDrainCounter::new(1)];
        let signal = PipelineSignal::new();

        // Watchdog: a wedge parks forever; abort so the test FAILS loudly.
        let done = Arc::new(AtomicBool::new(false));
        {
            let done = Arc::clone(&done);
            std::thread::spawn(move || {
                for _ in 0..200 {
                    std::thread::sleep(std::time::Duration::from_millis(25));
                    if done.load(Ordering::SeqCst) {
                        return;
                    }
                }
                eprintln!("driver_round_robins_all_live_before_parking: WEDGED");
                std::process::abort();
            });
        }

        // Forward walk (ChainOrderScheduler) tries `wedged` (idx 0) first: it
        // yields NoProgress, and the loop MUST proceed to `opener` in the same
        // pass, flip the gate, then finish `wedged` on the next pass.
        let mut worker = WorkerCore::driver(wedged);
        run_worker_loop(
            &mut worker,
            &mut entries,
            &contexts,
            &drain_counters,
            &signal,
            None,
            &crate::liveness::LivenessCounter::new(1),
            &crate::runtime::scheduler::ChainOrderScheduler,
            None,
            0,
            None,
            false,
        );
        done.store(true, Ordering::SeqCst);

        assert!(gate.load(Ordering::Acquire), "the sibling opener must have run");
        assert!(!signal.is_done(), "clean completion, no error");
    }

    // ── Sticky forward-progress ─────────────────────────────────────────────
    //
    // The source below breaks the `StepOutcome::Progress` contract on purpose:
    // it keeps reporting `Progress` while its output is full and it moves
    // nothing. Even so, a sticky pair must not wedge. Two
    // things must hold for the pair below to make progress at ONE worker: the
    // sticky burst is bounded, and round-robin does not restart the walk on the
    // sticky owner's `Progress` (which would break the pass before the sink is
    // ever reached). Removing either one hangs `sticky_holding_source_yields_to_
    // its_draining_consumer`.

    /// Sticky source over a capacity-1 output. Emits `remaining` items; when the
    /// transport rejects a push it *holds* the item and reports `Progress`, even
    /// on a retry that is still rejected (which the contract says should be
    /// `NoProgress`). That is the outcome that makes an unbounded sticky loop
    /// spin forever.
    struct StickyHoldingSource {
        remaining: u32,
        held: Option<u32>,
        calls: Arc<AtomicUsize>,
    }
    impl Step for StickyHoldingSource {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "StickyHoldingSource",
                kind: StepKind::Exclusive,
                sticky: true,
                // Capacity 1 so the second item in a burst always backs up.
                output_queues: vec![QueueSpec::CountBounded { capacity: 1 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            // Flush-first: retry the held item before producing a new one.
            if let Some(item) = self.held.take() {
                return match ctx.outputs.push(item) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(unpushed) => {
                        self.held = Some(unpushed.into_item());
                        Ok(StepOutcome::Progress)
                    }
                };
            }
            if self.remaining == 0 {
                return Ok(StepOutcome::Finished);
            }
            let item = self.remaining;
            self.remaining -= 1;
            match ctx.outputs.push(item) {
                Ok(()) => Ok(StepOutcome::Progress),
                Err(unpushed) => {
                    self.held = Some(unpushed.into_item());
                    Ok(StepOutcome::Progress)
                }
            }
        }
    }

    /// Sink that pops one item per dispatch — the only thing that frees a slot in
    /// the source's capacity-1 output.
    struct DrainingSink {
        received: Arc<Mutex<Vec<u32>>>,
    }
    impl Step for DrainingSink {
        type Input = u32;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "DrainingSink",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(item) => {
                    self.received.lock().push(item);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// One worker owning a sticky source (capacity-1 output) and the sink that
    /// drains it must deliver every item and terminate. A watchdog aborts on a
    /// wedge so the failure is loud rather than a silent hang.
    #[test]
    fn sticky_holding_source_yields_to_its_draining_consumer() {
        const N_ITEMS: u32 = 4;

        let mut graph = ChainGraph::new();
        let src = graph.register_step("StickyHoldingSource", 1);
        let sink = graph.register_step("DrainingSink", 0);
        graph.wire(src, BranchIdx(0), sink);

        let calls = Arc::new(AtomicUsize::new(0));
        let received = Arc::new(Mutex::new(Vec::new()));
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(StickyHoldingSource {
                remaining: N_ITEMS,
                held: None,
                calls: Arc::clone(&calls),
            })),
            Box::new(TypedStep::new(DrainingSink { received: Arc::clone(&received) })),
        ];
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));

        let mut entries: Vec<WorkerStepEntry> =
            steps.into_iter().map(|step| WorkerStepEntry::Exclusive { step }).collect();
        let drain_counters = vec![StepDrainCounter::new(1), StepDrainCounter::new(1)];
        let signal = PipelineSignal::new();

        // Watchdog: an unbounded sticky loop never returns, so abort loudly.
        let done = Arc::new(AtomicBool::new(false));
        {
            let done = Arc::clone(&done);
            std::thread::spawn(move || {
                for _ in 0..400 {
                    std::thread::sleep(std::time::Duration::from_millis(25));
                    if done.load(Ordering::SeqCst) {
                        return;
                    }
                }
                eprintln!(
                    "sticky_holding_source_yields_to_its_draining_consumer: WEDGED — the \
                     sticky fast-path is starving its downstream consumer"
                );
                std::process::abort();
            });
        }

        // Forward walk: the source is idx 0 and is this worker's sticky owner.
        let mut worker = WorkerCore::new(0, Some(src), Some(src));
        run_worker_loop(
            &mut worker,
            &mut entries,
            &contexts,
            &drain_counters,
            &signal,
            None,
            &crate::liveness::LivenessCounter::new(1),
            &crate::runtime::scheduler::ChainOrderScheduler,
            None,
            0,
            None,
            false,
        );
        done.store(true, Ordering::SeqCst);

        assert_eq!(
            *received.lock(),
            (1..=N_ITEMS).rev().collect::<Vec<u32>>(),
            "every item must reach the sink, in emission order"
        );
        assert!(!signal.is_done(), "clean completion, no error");
        // Each item costs at most one full sticky burst plus a round-robin
        // dispatch, so the source cannot have been called an unbounded number of
        // times. Loose on purpose — it pins "bounded", not an exact schedule.
        let n_calls = calls.load(Ordering::Relaxed);
        let ceiling = (usize::try_from(N_ITEMS).unwrap() + 2) * (STICKY_BURST_LIMIT + 2);
        assert!(
            n_calls <= ceiling,
            "source dispatches must stay bounded by the sticky burst limit: \
             {n_calls} calls > {ceiling}"
        );
    }

    /// Non-sticky source over a capacity-1 output that follows the held-slot
    /// convention every production step uses: a *new* hold reports `Progress`
    /// (it claimed an item), while a retry that is still rejected reports
    /// `NoProgress` (it moved nothing).
    struct HoldingSource {
        remaining: u32,
        held: Option<u32>,
    }
    impl Step for HoldingSource {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "HoldingSource",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 1 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(item) = self.held.take() {
                return match ctx.outputs.push(item) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(unpushed) => {
                        self.held = Some(unpushed.into_item());
                        Ok(StepOutcome::NoProgress)
                    }
                };
            }
            if self.remaining == 0 {
                return Ok(StepOutcome::Finished);
            }
            let item = self.remaining;
            self.remaining -= 1;
            if let Err(unpushed) = ctx.outputs.push(item) {
                self.held = Some(unpushed.into_item());
            }
            Ok(StepOutcome::Progress)
        }
    }

    /// Under a permanently raised refill signal, one worker walking a non-sticky
    /// source that holds on a full output must still reach the sink that drains
    /// it, whether the sink falls after the refill prefix (reverse part) or
    /// inside it (forward part). The priority restart after the source's new
    /// hold costs one extra visit, whose still-held retry reports `NoProgress`,
    /// so the walk continues to the sink instead of restarting again.
    #[rstest::rstest]
    #[case::sink_after_refill_prefix(StepIdx(0))]
    #[case::sink_inside_refill_prefix(StepIdx(1))]
    fn refill_walk_reaches_the_consumer_of_a_holding_step(#[case] refill_through: StepIdx) {
        const N_ITEMS: u32 = 4;

        let mut graph = ChainGraph::new();
        let src = graph.register_step("HoldingSource", 1);
        let sink = graph.register_step("DrainingSink", 0);
        graph.wire(src, BranchIdx(0), sink);

        let received = Arc::new(Mutex::new(Vec::new()));
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(HoldingSource { remaining: N_ITEMS, held: None })),
            Box::new(TypedStep::new(DrainingSink { received: Arc::clone(&received) })),
        ];
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));

        let mut entries: Vec<WorkerStepEntry> =
            steps.into_iter().map(|step| WorkerStepEntry::Exclusive { step }).collect();
        let drain_counters = vec![StepDrainCounter::new(1), StepDrainCounter::new(1)];
        let signal = PipelineSignal::new();

        // Never bound, so the read-ahead cap is absent and the refill walk
        // applies for the whole run.
        let scheduler = crate::runtime::scheduler::RefillDrainScheduler::new(
            Arc::new(AtomicBool::new(true)),
            refill_through,
            BranchIdx(0),
            u64::MAX,
        );
        assert_eq!(
            crate::runtime::scheduler::Scheduler::walk(&scheduler),
            WalkDirection::RefillThenReverse(refill_through)
        );

        // Watchdog: a walk that never reaches the sink spins forever.
        let done = Arc::new(AtomicBool::new(false));
        {
            let done = Arc::clone(&done);
            std::thread::spawn(move || {
                for _ in 0..400 {
                    std::thread::sleep(std::time::Duration::from_millis(25));
                    if done.load(Ordering::SeqCst) {
                        return;
                    }
                }
                eprintln!(
                    "refill_walk_reaches_the_consumer_of_a_holding_step: WEDGED — the refill \
                     walk never reaches the holding step's consumer"
                );
                std::process::abort();
            });
        }

        let mut worker = WorkerCore::new(0, None, None);
        run_worker_loop(
            &mut worker,
            &mut entries,
            &contexts,
            &drain_counters,
            &signal,
            None,
            &crate::liveness::LivenessCounter::new(1),
            &scheduler,
            None,
            0,
            None,
            false,
        );
        done.store(true, Ordering::SeqCst);

        assert_eq!(
            *received.lock(),
            (1..=N_ITEMS).rev().collect::<Vec<u32>>(),
            "every item must reach the sink, in emission order"
        );
        assert!(!signal.is_done(), "clean completion, no error");
    }

    /// A step that spends a measurable, nonzero span inside `try_run` so its
    /// dispatch busy-time rounds above 0 ns.
    #[derive(Clone)]
    struct SlowFinish;
    impl Step for SlowFinish {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "SlowFinish",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            std::thread::sleep(std::time::Duration::from_micros(200));
            Ok(StepOutcome::Finished)
        }
    }

    /// A driver dispatch (`is_driver = true`) records the step's busy on the
    /// off-pool detached line keyed by THAT step's own index — so a multi-step
    /// `Shared` group attributes each member's real time, not the whole thread's
    /// under one name. A pool dispatch (`is_driver = false`) records nothing there.
    #[test]
    fn driver_dispatch_records_detached_busy_by_own_step() {
        let mut graph = ChainGraph::new();
        let a = graph.register_step("SlowFinish", 1);
        let sink = graph.register_step("Sink", 0);
        graph.wire(a, BranchIdx(0), sink);
        let steps: Vec<Box<dyn ErasedStep>> =
            vec![Box::new(TypedStep::new(SlowFinish)), Box::new(TypedStep::new(SinkStep))];
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let mid = StepIdx(0);
        let counter = StepDrainCounter::new(1);
        let signal = PipelineSignal::new();

        // Driver dispatch of step 0 → its busy lands on the detached line keyed
        // to step 0 (not some group primary).
        let stats = Arc::new(PipelineStats::new(vec!["SlowFinish", "Sink"]));
        let mut entry = WorkerStepEntry::Owned { step: Box::new(TypedStep::new(SlowFinish)) };
        let _ = dispatch_one_step(
            &mut entry,
            mid,
            &contexts,
            &counter,
            &signal,
            Some(&stats),
            &LivenessCounter::new(1),
            0,
            true,
            None,
            0,
            None,
            &mut LoopDiag::new(false),
        );
        let snap = stats.snapshot();
        assert!(
            snap.detached.iter().any(|&(step, name, busy, ..)| {
                step == mid.0 && name == "SlowFinish" && busy > 0
            }),
            "driver dispatch must record detached busy for the dispatched step itself"
        );

        // Pool dispatch (is_driver=false) records nothing on the detached line.
        // A FRESH counter: the driver dispatch above consumed the first one, so
        // reusing it would make `observe_drain()` return false here and skip
        // `mark_outputs_drained` — the second dispatch would silently stop
        // exercising the same output-close path as the first.
        let counter_pool = StepDrainCounter::new(1);
        let stats_pool = Arc::new(PipelineStats::new(vec!["SlowFinish", "Sink"]));
        let mut entry_pool = WorkerStepEntry::Owned { step: Box::new(TypedStep::new(SlowFinish)) };
        let _ = dispatch_one_step(
            &mut entry_pool,
            mid,
            &contexts,
            &counter_pool,
            &signal,
            Some(&stats_pool),
            &LivenessCounter::new(1),
            0,
            false,
            None,
            0,
            None,
            &mut LoopDiag::new(false),
        );
        assert!(
            stats_pool.snapshot().detached.is_empty(),
            "pool dispatch must not record on the off-pool detached line"
        );
    }

    /// `u32 -> u32` pass-through used as a Detached group member: pops one item,
    /// pushes it on (holds on backpressure), finishes on drained input.
    struct PassDet {
        name: &'static str,
        held: Option<u32>,
    }
    impl Step for PassDet {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: self.name,
                kind: StepKind::Detached,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 64 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn detached_group(&self) -> crate::step::DetachedGroup {
            crate::step::DetachedGroup::Shared("g")
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(v) = self.held.take() {
                if let Err(u) = ctx.outputs.push(v) {
                    self.held = Some(u.into_item());
                    return Ok(StepOutcome::NoProgress);
                }
                return Ok(StepOutcome::Progress);
            }
            match ctx.input.pop() {
                Some(v) => match ctx.outputs.push(v) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(u) => {
                        self.held = Some(u.into_item());
                        Ok(StepOutcome::NoProgress)
                    }
                },
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// `() -> u32` source with the same types and kind as `three_step_chain`'s
    /// step 0: pushes one item and reports `Progress` on its first call,
    /// `Finished` after.
    #[derive(Default)]
    struct OneShotProgress {
        fired: bool,
    }
    impl Step for OneShotProgress {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "OneShot",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if self.fired {
                return Ok(StepOutcome::Finished);
            }
            self.fired = true;
            ctx.outputs.push(1).map_err(|_| io::Error::other("edge has room"))?;
            Ok(StepOutcome::Progress)
        }
    }

    /// Eight heap bytes, so an 8-byte byte-bounded edge holds exactly one.
    #[derive(Debug)]
    struct Item8;
    impl crate::item::HeapSize for Item8 {
        fn heap_size(&self) -> usize {
            8
        }
    }

    /// Emits `n` [`Item8`]s into a 1 MiB byte-bounded edge, so its pushes never
    /// refuse.
    struct Item8Source {
        n: u32,
    }
    impl Step for Item8Source {
        type Input = ();
        type Outputs = Single<Item8>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Src",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 20 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if self.n == 0 {
                return Ok(StepOutcome::Finished);
            }
            self.n -= 1;
            ctx.outputs.push(Item8).map_err(|_| io::Error::other("the 1 MiB edge has room"))?;
            Ok(StepOutcome::Progress)
        }
    }

    /// A Parallel step behind a shared admission cap, shaped like the
    /// Process-family adapters: re-push the held item first, then take a permit,
    /// pop, hold the permit for `work`, and push; a refused push holds the item.
    struct CappedHolder {
        cap: Arc<crate::admission::PhaseCap>,
        edge: QueueSpec,
        work: std::time::Duration,
        held: Option<crate::handles::Unpushed<Item8>>,
    }
    impl Step for CappedHolder {
        type Input = Item8;
        type Outputs = Single<Item8>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "CappedHolder",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![self.edge],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn phase_cap(&self) -> Option<&crate::admission::PhaseCap> {
            Some(&self.cap)
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(u) = self.held.take() {
                return match ctx.outputs.retry(u) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(u) => {
                        self.held = Some(u);
                        Ok(StepOutcome::NoProgress)
                    }
                };
            }
            let _admission = match crate::admission::admit_input(ctx.input, Some(&self.cap)) {
                Ok(a) => a,
                Err(outcome) => return Ok(outcome),
            };
            match ctx.input.pop() {
                Some(item) => {
                    std::thread::sleep(self.work);
                    if let Err(u) = ctx.outputs.push(item) {
                        self.held = Some(u);
                    }
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
        fn new_worker_copy(&self) -> Self {
            Self { cap: Arc::clone(&self.cap), edge: self.edge, work: self.work, held: None }
        }
    }

    /// A Serial sink that sleeps `work` after each pop.
    struct SlowSink {
        work: std::time::Duration,
    }
    impl Step for SlowSink {
        type Input = Item8;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Sink",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(_) => {
                    std::thread::sleep(self.work);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// A finished step is never booked a shared driver's idle: once the step
    /// that last progressed (or the primary) has left the worklist, the idle
    /// goes to a step that is still live.
    #[test]
    fn driver_idle_is_never_attributed_to_a_finished_step() {
        let (a, b, c) = (StepIdx(1), StepIdx(2), StepIdx(3));
        let mut diag = LoopDiag::new(true);
        let all = LiveSteps::from_order(vec![a, b, c]);
        assert_eq!(diag.attributed(a, &all), a, "before any progress: the primary");
        diag.note(b, &Ok(StepOutcome::Progress));
        assert_eq!(diag.attributed(a, &all), b, "the step that last progressed");
        diag.note(b, &Ok(StepOutcome::Finished));
        let without_b = LiveSteps::from_order(vec![a, c]);
        assert_eq!(diag.attributed(a, &without_b), a, "b finished: back to the live primary");
        let only_c = LiveSteps::from_order(vec![c]);
        assert_eq!(diag.attributed(a, &only_c), c, "a and b finished: the live step");
    }

    /// A shared driver's park is booked to the step that last made progress
    /// on that thread (B after an item flowed A→B), and to the group's primary
    /// (A) before any progress. Pinned by driving `[A, B]` as one driver row and
    /// pushing one item from the test thread between two idle windows.
    #[test]
    fn driver_idle_is_attributed_to_the_last_progressing_step() {
        use crate::step::OutputHandles;
        let mut graph = ChainGraph::new();
        let src = graph.register_step("Src", 1);
        let a = graph.register_step("A", 1);
        let b = graph.register_step("B", 1);
        let sink = graph.register_step("Sink", 0);
        graph.wire(src, BranchIdx(0), a);
        graph.wire(a, BranchIdx(0), b);
        graph.wire(b, BranchIdx(0), sink);
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(SrcFinished)),
            Box::new(TypedStep::new(PassDet { name: "A", held: None })),
            Box::new(TypedStep::new(PassDet { name: "B", held: None })),
            Box::new(TypedStep::new(SinkStep)),
        ];
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let stats = Arc::new(PipelineStats::new(vec!["Src", "A", "B", "Sink"]));
        let signal = PipelineSignal::new();
        let drain_counters: Vec<Arc<StepDrainCounter>> =
            (0..4).map(|_| StepDrainCounter::new(1)).collect();

        // Driver row: A and B Owned, everything else Skip.
        let mut row: Vec<WorkerStepEntry> = (0..4).map(|_| WorkerStepEntry::Skip).collect();
        let mut it = steps.into_iter();
        let _src = it.next();
        row[1] = WorkerStepEntry::Owned { step: it.next().unwrap() };
        row[2] = WorkerStepEntry::Owned { step: it.next().unwrap() };

        let driver = {
            let contexts = Arc::clone(&contexts);
            let signal = Arc::clone(&signal);
            let stats = Arc::clone(&stats);
            let drain_counters = drain_counters.clone();
            std::thread::spawn(move || {
                let mut worker = WorkerCore::driver(a);
                run_worker_loop(
                    &mut worker,
                    &mut row,
                    &contexts,
                    &drain_counters,
                    &signal,
                    Some(&stats),
                    &crate::liveness::LivenessCounter::new(1),
                    &crate::runtime::scheduler::DrainFirstScheduler,
                    None,
                    0,
                    None,
                    true,
                );
            })
        };

        let producer = contexts.outputs[0].downcast_ref::<OutputHandles<Single<u32>>>().unwrap();
        let sink_in = contexts.inputs[3].downcast_ref::<BranchInputHandle<u32>>().unwrap();
        let parks = |snap: &crate::runtime::stats::StatsSnapshot, idx: usize| {
            snap.detached.iter().find(|d| d.0 == idx).map_or(0, |d| d.4)
        };
        // Window 1: nothing pushed yet → parks booked to the primary (A). Wait
        // for the first park instead of sleeping a fixed time.
        // Every wait is bounded, so a missing attribution fails the test rather
        // than hanging it.
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
        let wait_until = |what: &str, cond: &dyn Fn() -> bool| {
            while !cond() {
                assert!(std::time::Instant::now() < deadline, "timed out waiting for {what}");
                std::thread::yield_now();
            }
        };
        wait_until("A's first park", &|| parks(&stats.snapshot(), 1) > 0);
        assert_eq!(parks(&stats.snapshot(), 2), 0, "B has not progressed yet");

        // One item flows A → B; window 2: parks now booked to B.
        assert!(producer.push(7).is_ok());
        wait_until("the item at the sink", &|| sink_in.pop().is_some());
        let mid = stats.snapshot();
        wait_until("a park booked to B", &|| parks(&stats.snapshot(), 2) > parks(&mid, 2));
        let after = stats.snapshot();
        assert_eq!(
            parks(&after, 1),
            parks(&mid, 1),
            "A gains no parks once B was the last to progress"
        );

        producer.mark_all_drained();
        driver.join().expect("driver joins");
    }

    /// A step refused by an admission cap returns `Capped`. That is not a held
    /// item: the input item is still queued, and the refusal never reaches a
    /// transport push, so no hold or held retry is recorded (holds are measured
    /// at the transport — see
    /// `queues::tests::hold_clock_times_first_reject_to_next_push`).
    ///
    /// Run as a real two-worker pipeline with `--pipeline-stats`, so every edge
    /// is byte-bounded and carries a hold clock. `cap_refusal` has the two
    /// workers contend for one permit over roomy edges: many refusals, no
    /// holds. `push_refusal` is the control on the same chain: an uncontended
    /// cap and a one-item edge into a slow sink, so the hold clock does book
    /// holds there — the `cap_refusal` zero is not a fixture that cannot count.
    #[rstest::rstest]
    #[case::cap_refusal(1, QueueSpec::ByteBounded { limit_bytes: 1 << 20 }, 200, 0, true)]
    #[case::push_refusal(2, QueueSpec::ByteBounded { limit_bytes: 8 }, 0, 200, false)]
    fn admission_refusal_is_not_a_hold(
        #[case] permits: usize,
        #[case] edge: QueueSpec,
        #[case] capped_work_us: u64,
        #[case] sink_work_us: u64,
        #[case] cap_refuses: bool,
    ) {
        use std::time::Duration;

        use crate::admission::PhaseCap;
        use crate::builder::{Pipeline, PipelineConfig};
        const N: u32 = 400;
        let cap = PhaseCap::new("test-phase", permits);
        let builder = Pipeline::builder();
        builder
            .chain(Item8Source { n: N })
            .chain(CappedHolder {
                cap: Arc::clone(&cap),
                edge,
                work: Duration::from_micros(capped_work_us),
                held: None,
            })
            .chain(SlowSink { work: Duration::from_micros(sink_work_us) })
            .into_sink_marker();
        let pipeline = builder.build().expect("build");
        let stats = pipeline.stats();
        pipeline
            .run(PipelineConfig { threads: 2, ..Default::default() }.with_stats(Arc::clone(&stats)))
            .expect("run");
        let snap = stats.snapshot();
        let holds = |name: &str| {
            let (_, s) = snap.steps.iter().find(|(n, _)| *n == name).expect("step row");
            (s.holds, s.held_retries)
        };
        if cap_refuses {
            assert!(cap.refused() > 0, "two workers on one permit must be refused: {snap:?}");
            for name in ["Src", "CappedHolder", "Sink"] {
                assert_eq!(holds(name), (0, 0), "a cap refusal is not a hold ({name})");
            }
        } else {
            assert_eq!(cap.refused(), 0, "the control's holds come from push refusals only");
            assert!(holds("CappedHolder").0 > 0, "a refused push into the full edge is a hold");
        }
    }

    /// Every `Progress` issues exactly one `notify_one` and every `Finished` one
    /// `notify_all`. Observed through the event-count itself (an armed key
    /// reports `Woken` iff the generation moved), not through the stats the
    /// same branch records.
    #[test]
    fn legacy_mode_notifies_exactly_as_before() {
        use std::time::Duration;

        use crate::runtime::event_count::WaitOutcome;
        let (steps, graph, _runs) = three_step_chain(StepKind::Parallel);
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let stats = Arc::new(PipelineStats::new(vec!["Src", "Finish", "Sink"]));
        let ec = Arc::new(PoolEventCount::new(2));
        let signal = PipelineSignal::new();
        let counter = StepDrainCounter::new(1);

        // Finished → notify_all: an armed key sees the generation move.
        let mut entry = WorkerStepEntry::Owned { step: steps[1].clone_boxed() };
        let key = ec.prepare_wait();
        let info = dispatch_one_step(
            &mut entry,
            StepIdx(1),
            &contexts,
            &counter,
            &signal,
            Some(&stats),
            &LivenessCounter::new(1),
            0,
            false,
            None,
            0,
            Some(&*ec),
            &mut LoopDiag::new(false),
        )
        .unwrap();
        assert!(matches!(info.result, Ok(StepOutcome::Finished)));
        assert_eq!(ec.wait(key, Duration::ZERO), WaitOutcome::Woken, "Finished must notify_all");

        // Progress → exactly one notify_one: one bump, and a second arm sees none.
        let mut progress =
            WorkerStepEntry::Owned { step: Box::new(TypedStep::new(OneShotProgress::default())) };
        let key = ec.prepare_wait();
        let info = dispatch_one_step(
            &mut progress,
            StepIdx(0),
            &contexts,
            &StepDrainCounter::new(1),
            &signal,
            Some(&stats),
            &LivenessCounter::new(1),
            0,
            false,
            None,
            0,
            Some(&*ec),
            &mut LoopDiag::new(false),
        )
        .unwrap();
        assert!(matches!(info.result, Ok(StepOutcome::Progress)));
        assert_eq!(ec.wait(key, Duration::ZERO), WaitOutcome::Woken, "Progress must notify_one");
        let key = ec.prepare_wait();
        assert_eq!(
            ec.wait(key, Duration::ZERO),
            WaitOutcome::TimedOut,
            "exactly one bump per Progress"
        );
        // Rendering check only (the oracle above is the behavioural one).
        assert_eq!(stats.snapshot().steps[0].1.notifies_issued, 1);
    }

    #[test]
    fn dispatch_stamps_running_step_when_board_present() {
        // A dispatch of step idx 3 with a board present must leave the board's
        // slot reading Running(step 3) at the point try_run is entered. We use a
        // step that finishes immediately; after dispatch the slot holds the last
        // stamp (Running, 3) since nothing overwrites it here.
        use crate::runtime::worker_state::{WorkerState, WorkerStateBoard};
        let mut graph = ChainGraph::new();
        let a = graph.register_step("SrcFinished", 1);
        let sink = graph.register_step("Sink", 0);
        graph.wire(a, BranchIdx(0), sink);
        let steps: Vec<Box<dyn ErasedStep>> =
            vec![Box::new(TypedStep::new(SrcFinished)), Box::new(TypedStep::new(SinkStep))];
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let counter = StepDrainCounter::new(1);
        let signal = PipelineSignal::new();
        let board = WorkerStateBoard::new(1);
        let mut entry = WorkerStepEntry::Owned { step: Box::new(TypedStep::new(SrcFinished)) };
        let _ = dispatch_one_step(
            &mut entry,
            StepIdx(0),
            &contexts,
            &counter,
            &signal,
            None,
            &LivenessCounter::new(1),
            0,
            false,
            Some(&board),
            0,
            None,
            &mut LoopDiag::new(false),
        );
        assert_eq!(board.read(0), (WorkerState::Running, Some(StepIdx(0))));
    }

    #[test]
    fn dispatch_with_no_board_is_a_noop_on_state() {
        // The telemetry-off path must not require a board (None) and must behave
        // exactly as before — this dispatch simply must not panic and must run.
        let mut graph = ChainGraph::new();
        let a = graph.register_step("SrcFinished", 1);
        let sink = graph.register_step("Sink", 0);
        graph.wire(a, BranchIdx(0), sink);
        let steps: Vec<Box<dyn ErasedStep>> =
            vec![Box::new(TypedStep::new(SrcFinished)), Box::new(TypedStep::new(SinkStep))];
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let counter = StepDrainCounter::new(1);
        let signal = PipelineSignal::new();
        let mut entry = WorkerStepEntry::Owned { step: Box::new(TypedStep::new(SrcFinished)) };
        let info = dispatch_one_step(
            &mut entry,
            StepIdx(0),
            &contexts,
            &counter,
            &signal,
            None,
            &LivenessCounter::new(1),
            0,
            false,
            None,
            0,
            None,
            &mut LoopDiag::new(false),
        )
        .unwrap();
        assert!(matches!(info.result, Ok(StepOutcome::Finished)));
    }
}
