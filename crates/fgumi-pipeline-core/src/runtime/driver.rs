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
//!   4. If no work happened this iteration, idle: block on the pool
//!      event-count (an unpinned pool worker that holds nothing), or park with
//!      an exponential-backoff `park_timeout` (pinned workers, a worker holding
//!      an item, the lone worker of a one-worker run, and driver threads) that a
//!      directed wake (`runtime::wake`) can end early.
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
use crate::runtime::live::LiveSteps;
use crate::runtime::scheduler::{Scheduler, WalkDirection};
use crate::runtime::stats::{PipelineStats, WakeCounts};
use crate::runtime::storage::WorkerStepEntry;
use crate::runtime::wake::{PoolHandle, PushSnapshot, WakePlan};
use crate::runtime::worker_core::{WorkerCore, WorkerRole};
use crate::runtime::worker_state::WorkerStateBoard;
use crate::signal::{PipelineError, PipelineSignal};
use crate::step::StepOutcome;
use crate::topology::StepIdx;

/// Maximum consecutive sticky re-entries before the worker drops back to a
/// round-robin pass.
///
/// `StepOutcome::Progress` means the step "took, pushed or newly held an item" (see the
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
/// output close on `Finished` (init = the clone count for Parallel —
/// `ParallelHosts::clone_count()` — so only the last clone closes the shared
/// output; init 1 for Serial/Exclusive/Detached).
/// `signal` — error/cancel broadcast.
/// `board` — optional per-thread state board for scheduling telemetry; `None`
/// keeps the loop's stamping a no-op (telemetry-off parity).
/// `state_slot` — this thread's slot in `board` (ignored when `board` is `None`).
/// `wake` — the chain's wake plan: owns the pool event-count (`wake.pool()`,
/// the idle wait of non-pinned workers that hold nothing) and routes every
/// `Progress`/`Finished` to the thread that consumes it.
/// `pinned` — this worker is the sole eligible dispatcher of some step; it idles
/// on its timer (`park_timeout`), never on the event-count.
#[allow(clippy::too_many_arguments)] // per-step shared state plus the liveness shard; a struct would only rename it
#[allow(clippy::too_many_lines)] // the whole-pass sticky + round-robin + backoff discipline is documented inline; splitting it would scatter the invariant across functions
pub(crate) fn run_worker_loop(
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
    wake: &WakePlan,
    pinned: bool,
) {
    // A pinned worker (the sole eligible dispatcher of an Exclusive/sticky step,
    // or the affinity target of a Serial step) must NOT deep-park on the shared
    // event-count: `notify_one` cannot target a specific worker, so a push meant
    // to wake *this* worker's pinned step could wake a peer that `Skip`s it,
    // leaving the pinned step's owner asleep to its timeout. Pinned workers idle
    // on their timer with `park_timeout`; in a Directed chain the wake plan
    // unparks them directly. Its step is also often driven by out-of-engine
    // events (a reader's prefetch thread, the affinity-pinned grouper) that no
    // push notifies anyway. The plan still routes a pinned worker's own
    // `Progress` wakes.
    let idle_parker = if pinned { None } else { wake.pool() };
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

        // The holding and cap-refused idle hints describe this iteration's
        // passes only: a stale one would send a thread that holds nothing, or
        // that no cap recorded, to a timer park.
        let _ = crate::runtime::wake_slot::take_rejected();
        let _ = crate::runtime::wake_slot::take_cap_refused();
        // A cap's record of this thread covers the pass that made it: a new
        // pass starts clean, and re-records on whatever refuses it.
        wake.start_pass();

        // Time the dispatch (busy) section vs the backoff sleep (idle) below,
        // per worker, only when stats are on (`Instant::now()` is otherwise not
        // called — the no-stats loop stays zero-cost). It starts before the
        // resume dispatch, which can run a whole capped item.
        let work_start = stats.map(|_| Instant::now());

        // 0. A phase cap refused this thread on its last pass and a release may
        // have woken it for that permit: poll the refusing step first, before
        // the sticky burst and the round-robin walk can take the thread
        // elsewhere. No input there takes no permit; the pass continues as
        // usual.
        if let Some(step) = crate::runtime::wake_slot::take_resume_step()
            && live.contains(step)
            && let Some(info) = dispatch_one_step(
                &mut entries[step.0],
                step,
                contexts,
                &drain_counters[step.0],
                signal,
                stats,
                liveness,
                worker.thread_id,
                is_driver,
                board,
                state_slot,
                wake,
                &mut diag,
            )
        {
            match info.result {
                Ok(StepOutcome::Progress) => did_work = true,
                Ok(StepOutcome::Finished) => {
                    did_work = true;
                    live.remove(step);
                    if worker.sticky_owner == Some(step) {
                        sticky_live = false;
                    }
                }
                Ok(StepOutcome::NoProgress | StepOutcome::Contention | StepOutcome::Capped) => {}
                Err(io_err) => {
                    signal.record_error(PipelineError::Io { step: info.name, source: io_err });
                }
            }
        }

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
                    wake,
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
                wake,
                &mut diag,
                false,
            );
            did_work |= outcome.did_work;
            if outcome.removed_sticky_owner {
                // The sticky owner finished during round-robin; disable the
                // fast path so subsequent iterations skip the sticky block.
                sticky_live = false;
            }
        }
        // Forward any cap wake this pass was handed and did not spend.
        wake.end_pass();

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
        } else {
            // A thread whose push was rejected since its last idle may be
            // holding an item. Its wake is a direct `unpark` from whoever makes
            // room, which the event-count cannot deliver — so it idles on its
            // timer, like a pinned worker. Always `false` on edges without a
            // holder set (every Legacy chain), so Legacy idle is unchanged.
            let mut holding = crate::runtime::wake_slot::take_rejected();
            // A phase cap refused this thread and recorded it: the release that
            // frees a permit `unpark`s it directly, so it too idles on its timer.
            let mut cap_refused = crate::runtime::wake_slot::take_cap_refused();
            let ec_parker = if holding || cap_refused { None } else { idle_parker };
            if let Some(ec) = ec_parker {
                // Event-count park: block instead of spin-poll when the whole
                // pass found no work. The two-phase protocol closes the
                // lost-wakeup race: arm (register + fence) → re-poll the real
                // condition → block only if still empty.
                //
                // Arm and re-poll *before* stamping `Parked` or starting the
                // idle timer: a productive re-poll (work appeared in the
                // arm→wait window) must not be recorded as idle, nor leave the
                // board reading `Parked` while the re-poll's `dispatch_one_step`
                // runs real work.

                // Phase 1: arm. The fence in `prepare_wait` orders the waiter
                // registration before the re-poll below, so a producer that
                // publishes work after this point either sees us as a waiter
                // (and bumps the generation) or we see its item on the re-poll.
                let key = ec.prepare_wait();
                // The re-poll is a pass of its own, so it starts and ends like
                // every pass: its start clears the thread's records and works out
                // which cap wakes it owes, and its end forwards one it did not
                // spend. In a run no record reaches here (a refusal sends the
                // thread to the timer park instead, and only the cap-parked
                // re-poll skips), so both are cheap no-ops kept for uniformity.
                wake.start_pass();
                // Until this arm is balanced, a `request_worker` from a step in
                // the re-poll below must not count this worker as the parked
                // worker to wake (it would wake only itself).
                crate::runtime::wake_slot::set_ec_armed(true);

                // Phase 2: re-poll the *real* condition — one full dispatch
                // pass. One pass per park episode (amortised over the block),
                // not per iteration.
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
                    wake,
                    &mut diag,
                    false,
                );
                wake.end_pass();
                if recheck.removed_sticky_owner {
                    sticky_live = false;
                }

                if recheck.did_work || signal.is_done() {
                    // Work appeared in the arm→wait window (or we are shutting
                    // down): do not block, and record no idle — the re-poll was
                    // productive. Balance the arm and re-ramp fresh.
                    ec.cancel_wait(key);
                    crate::runtime::wake_slot::set_ec_armed(false);
                    wake.acknowledge_request();
                    worker.reset_backoff();
                    continue;
                }
                let rejected = crate::runtime::wake_slot::take_rejected();
                let refused = crate::runtime::wake_slot::take_cap_refused();
                if rejected || refused {
                    // The re-poll was rejected (holding now) or refused by a
                    // cap that recorded it. Leave the event-count and fall
                    // through to the timer park, where an `unpark` reaches this
                    // thread.
                    holding |= rejected;
                    cap_refused |= refused;
                    ec.cancel_wait(key);
                    crate::runtime::wake_slot::set_ec_armed(false);
                    // Every exit from the event-count acknowledges a pending
                    // request, or one sent while this worker was armed stays
                    // pending past it (timer-parked workers never acknowledge).
                    wake.acknowledge_request();
                } else {
                    // No progress remains: only now do we actually park. Stamp
                    // `Parked` and start the idle timer immediately before
                    // blocking, so the idle attribution covers exactly the wait,
                    // not the preceding (possibly productive) re-poll.
                    if let Some(b) = board {
                        b.stamp(
                            state_slot,
                            crate::runtime::worker_state::WorkerState::Parked,
                            None,
                        );
                    }
                    let sleep_start = stats.map(|_| Instant::now());

                    // Phase 3: block until the generation moves (a peer's
                    // notify), the deadline elapses (self-heal), or a spurious
                    // wake. The exponential backoff is the wait *timeout* so a
                    // missed wakeup cannot hang longer than the cap, and the
                    // deadlock monitor still ticks. Re-poll on return regardless.
                    let deadline = worker.backoff_deadline();
                    let outcome = ec.wait(key, deadline);
                    crate::runtime::wake_slot::set_ec_armed(false);
                    wake.acknowledge_request();
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
                    continue;
                }
            }

            // Timer park: pinned workers, holders, the lone worker, drivers. A
            // Directed pool worker arms and re-polls first, so a `Pool` wake
            // that finds no event-count waiter can reach it.
            //
            // A worker parked only because a cap refused it arms the separate
            // cap-parked set instead: a `Pool` wake claims it only after every
            // direct-parked worker, and only when the consumer's cap admits
            // (uncapped, or a free permit). A wake for a step whose cap is full
            // would only see it refused again; the release that frees a permit
            // on that cap wakes it, through the record its refusal or its
            // re-poll's skip left.
            let cap_only = cap_refused && !holding;
            let armed = matches!(worker.role(), WorkerRole::Pool)
                && if cap_only {
                    wake.arm_cap_parked(worker.thread_id)
                } else {
                    wake.arm_direct(worker.thread_id)
                };
            let disarm = |wake: &WakePlan, w: usize| {
                if cap_only { wake.disarm_cap_parked(w) } else { wake.disarm_direct(w) }
            };
            if armed {
                // A new pass: it re-records on every cap it is refused by or
                // skips as full (whether or not the skipped step has input: a
                // step empty now may get input while its cap is still full, and
                // that record is the only route to this thread for it). So the
                // thread parks recorded on every full cap it would poll; a wake
                // spent on one whose step has nothing for it is forwarded by the
                // thread's next pass.
                wake.start_pass();
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
                    wake,
                    &mut diag,
                    cap_only,
                );
                wake.end_pass();
                if recheck.removed_sticky_owner {
                    sticky_live = false;
                }
                if recheck.did_work || signal.is_done() {
                    disarm(wake, worker.thread_id);
                    worker.reset_backoff();
                    continue;
                }
                // The re-poll left no cap record (the refusing cap had room by
                // then, and its step nothing to take): no release will wake this
                // thread, and the cap-parked set is reachable only for consumers
                // whose cap admits. It is an ordinary idle worker again, so run
                // one more main pass, which, if it also finds nothing, parks on
                // the normal path (the event-count, or direct-parked), where
                // every `Pool` wake reaches it. No idle is recorded for this
                // round.
                if cap_only && !crate::runtime::wake_slot::cap_recorded_pending() {
                    disarm(wake, worker.thread_id);
                    continue;
                }
            }
            if let Some(b) = board {
                b.stamp(state_slot, crate::runtime::worker_state::WorkerState::Parked, None);
            }
            // Parked only because a cap refused it: the cap's release unparks
            // it (one per freed permit), so it waits at the ramp's cap rather
            // than polling the full cap up the ramp from its floor.
            let deadline =
                if cap_only { worker.max_backoff_deadline() } else { worker.backoff_deadline() };
            let sleep_start = stats.map(|_| Instant::now());
            std::thread::park_timeout(deadline);
            worker.increase_backoff();
            if armed {
                disarm(wake, worker.thread_id);
            }
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

    // An exiting thread leaves no cap record behind (a release must not spend
    // a wake on it) and forwards a wake it was handed and never spent.
    wake.start_pass();
    wake.end_pass();

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
    wake: &WakePlan,
    diag: &mut LoopDiag,
    skip_capped: bool,
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
        // The re-poll of a worker that armed as cap-parked covers only what a
        // `Pool` wake may claim it for: consumers that are uncapped or whose
        // cap has a free permit. A step whose cap is full would only refuse it
        // again, so it is skipped — and the skip is recorded on that cap as a
        // refusal is, so the release that frees a permit wakes this thread.
        if skip_capped && wake.skip_full_cap(step_idx) {
            continue;
        }
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
            wake,
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
    // The chain's wake plan. On a `Progress` outcome (the step took, pushed or
    // newly held an item — including a reorder must-accept stash) it wakes the
    // thread that consumes what was produced and the holders this step's pops
    // made room for; on `Finished` (an output edge just closed) it wakes
    // everyone; on an idle outcome of a dispatch that flushed a held item it
    // delivers that push's forward wake (`on_flushed`). This is the wake seam:
    // these are exactly the outcomes that can give a parked thread new work,
    // and using the outcome (rather than the `BranchOutputHandle::push` site)
    // captures the stash-then-return-`Ok` reorder case without threading the
    // plan through every queue.
    wake: &WakePlan,
    // Per-thread diagnostic state; updated after the outcome is known.
    diag: &mut LoopDiag,
) -> Option<DispatchInfo> {
    // Gate snapshot: one Relaxed load per gated branch of a Detached producer,
    // nothing for an ungated step (`is_gated` is a precomputed bool).
    let mut snap = PushSnapshot::default();
    if wake.is_gated(step_idx) {
        wake.snapshot_pushes(step_idx, &mut snap);
    }
    let outputs_any: &(dyn Any + Send + Sync) = contexts.outputs[step_idx.0].as_ref();
    let mut ctx = ErasedStepCtx {
        input: contexts.inputs[step_idx.0].as_ref(),
        outputs: outputs_any,
        signal,
        counters: &contexts.step_counters[step_idx.0],
        pool: PoolHandle::new(wake, stats.map(|s| &**s), step_idx),
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
    // per-step `StepDrainCounter`): for a `Parallel` step (counter init = its
    // clone count)
    // every clone returns `Finished` independently when the shared input edge
    // is drained, but only the LAST clone to finish (the one that takes the
    // counter to 0) closes the shared output queue — otherwise a clone could
    // `mark_drained` while a sibling is still pushing (`try_push`-after-drained
    // panic). For `Serial`/`Exclusive` (counter init 1) the single finisher
    // wins on its first call, unchanged.
    //
    // INVARIANT: for a `Parallel` step, `counter` init == clone count ==
    // `ParallelHosts::clone_count()` — the host workers, or the one driver
    // hosting it — and a clone leaves its worklist ONLY by returning
    // `Finished`, so the counter reaches 0 exactly when every clone has
    // finished. Any future scheduler change that removes a Parallel clone for
    // another reason (work-stealing, per-worker early exit) — or makes a source
    // `Parallel` — would leave the counter stuck above 0 and never close the
    // shared output, hanging the downstream consumer. Storage, driver-group
    // membership and the init (`placement::drain_counter_inits`) all read the
    // one `plan_parallel_hosts` result; keep them and this gate in lockstep.
    if let Some(b) = board {
        b.stamp(state_slot, crate::runtime::worker_state::WorkerState::Running, Some(step_idx));
    }

    let mut info: Option<DispatchInfo> = match entry {
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
    // A held retry inside this dispatch succeeded (`BranchOutputHandle::retry`).
    // Taken on every call, so the next dispatch on this thread starts clear.
    let flushed = crate::runtime::wake_slot::take_flushed();
    upgrade_a_popping_dispatch(&mut info);

    if let Some(i) = info.as_ref() {
        diag.note(step_idx, &i.result);
        // A phase cap refused this thread and recorded it for a release wake:
        // the thread's next pass starts here, so the permit that wake promises
        // is spent on this step rather than left idle while the thread
        // round-robins to other work.
        if matches!(i.result, Ok(StepOutcome::Capped))
            && crate::runtime::wake_slot::cap_refused_pending()
        {
            crate::runtime::wake_slot::set_resume_step(step_idx);
        }
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
        // Wake seam: a productive dispatch may have given a parked thread new
        // work. `Progress` = "took, pushed or newly held an item" (incl. a reorder
        // must-accept stash that unblocks an ordinal waiter) → the plan wakes
        // the consumer's thread (a Legacy plan: one parked pool worker, which
        // cascades) and any holder this step's pops made room for. `Finished`
        // closed an output edge (a parked consumer must observe end-of-stream,
        // and a sink's `Finished` closes no edge but its peers may be waiting on
        // the drain latch) → wake everyone. When nobody is parked, each wake is
        // a fence and a load (see `PoolEventCount`) or an `unpark` of a running
        // thread.
        if let Ok(StepOutcome::Finished) = i.result {
            wake.on_finished(step_idx);
        } else {
            let r = wake.on_progress(step_idx, &snap);
            if let Some(stats) = stats {
                stats.record_wake(step_idx, r);
            }
        }
    } else if flushed && let Some(i) = info.as_ref() {
        deliver_flushed_wake(wake, step_idx, &i.result, &snap, stats);
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

/// A dispatch that took an input item made room on its input edge, so it
/// reports `Progress` whatever the step returned (the `StepOutcome::Progress`
/// contract), and the plan's reverse wake for that room runs on it. The pop
/// flag is taken on every call, so the next dispatch starts clear.
#[inline]
fn upgrade_a_popping_dispatch(info: &mut Option<DispatchInfo>) {
    if crate::runtime::wake_slot::take_popped()
        && let Some(i) = info.as_mut()
        && matches!(
            i.result,
            Ok(StepOutcome::NoProgress | StepOutcome::Contention | StepOutcome::Capped)
        )
    {
        i.result = Ok(StepOutcome::Progress);
    }
}

/// The flushed push's forward wake, for a dispatch whose held retry landed
/// (`BranchOutputHandle::retry`) and which then reported an idle outcome: the
/// flush-first preamble found an empty input, a refused permit or a contended
/// branch. `Progress` and `Finished` already deliver it in
/// [`dispatch_one_step`], and `Err` ends the run (the signal's terminal path
/// unparks everyone). Liveness is not bumped: the consumer's pop of the item
/// will. A Legacy plan returns an empty report and delivers nothing.
#[inline]
fn deliver_flushed_wake(
    wake: &WakePlan,
    step_idx: StepIdx,
    result: &std::io::Result<StepOutcome>,
    snap: &PushSnapshot,
    stats: Option<&Arc<PipelineStats>>,
) {
    if !matches!(
        result,
        Ok(StepOutcome::NoProgress | StepOutcome::Contention | StepOutcome::Capped)
    ) {
        return;
    }
    let r = wake.on_flushed(step_idx, snap);
    if r != WakeCounts::default()
        && let Some(stats) = stats
    {
        stats.record_wake(step_idx, r);
    }
}

#[cfg(test)]
mod tests {
    use std::io;
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

    use super::*;
    use crate::runtime::event_count::PoolEventCount;

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
            holder_only_queues: vec![],
            edges: vec![],
            step_counters: vec![],
        });
        let drain_counters: Vec<Arc<StepDrainCounter>> = vec![];
        let _ = ChainGraph::new();
        let mut worker = WorkerCore::new(0, None);

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
            &WakePlan::legacy(None),
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
            holder_only_queues: vec![],
            edges: vec![],
            step_counters: vec![],
        });
        let drain_counters: Vec<Arc<StepDrainCounter>> = vec![];
        let mut worker = WorkerCore::new(0, None);

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
            &WakePlan::legacy(None),
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

    /// Pops whatever is there and reports `NoProgress` anyway — a step that
    /// breaks the "took input ⇒ `Progress`" contract.
    struct PopsThenIdle;
    impl Step for PopsThenIdle {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "PopsThenIdle",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let _ = ctx.input.pop();
            Ok(StepOutcome::NoProgress)
        }
    }

    /// The framework enforces "a dispatch that took input reports `Progress`":
    /// a step that pops and returns `NoProgress` is recorded as `Progress`,
    /// bumps liveness, and runs the plan's `Progress` wake (here a Legacy
    /// plan's `notify_one`, observed on the event-count), so the wake for the
    /// room its pop made is not lost. A dispatch that popped nothing keeps its
    /// `NoProgress`: no liveness, no wake.
    #[test]
    fn a_dispatch_that_popped_reports_progress_whatever_the_step_said() {
        use std::time::Duration;

        use crate::runtime::event_count::WaitOutcome;
        use crate::step::OutputHandles;
        let (steps, graph, _runs) = three_step_chain(StepKind::Serial);
        let contexts =
            build_chain_contexts(&steps, &graph, crate::builder::InstrumentationLevel::Off, false);
        let producer = contexts.outputs[0].downcast_ref::<OutputHandles<Single<u32>>>().unwrap();
        assert!(producer.push(7).is_ok());
        let stats = Arc::new(PipelineStats::new(vec!["Src", "Finish", "Sink"]));
        let ec = Arc::new(PoolEventCount::new(2));
        let plan = WakePlan::legacy(Some(Arc::clone(&ec)));
        let liveness = LivenessCounter::new(1);
        let mut entry = WorkerStepEntry::Owned { step: Box::new(TypedStep::new(PopsThenIdle)) };
        let _ = crate::runtime::wake_slot::take_popped();
        let mut dispatch = || {
            dispatch_one_step(
                &mut entry,
                StepIdx(1),
                &contexts,
                &StepDrainCounter::new(1),
                &PipelineSignal::new(),
                Some(&stats),
                &liveness,
                0,
                false,
                None,
                0,
                &plan,
                &mut LoopDiag::new(false),
            )
            .expect("dispatched")
            .result
            .expect("ok")
        };
        let key = ec.prepare_wait();
        assert_eq!(dispatch(), StepOutcome::Progress, "it popped the item");
        assert_eq!(ec.wait(key, Duration::ZERO), WaitOutcome::Woken, "its wake ran");
        let key = ec.prepare_wait();
        assert_eq!(dispatch(), StepOutcome::NoProgress, "nothing to pop");
        assert_eq!(ec.wait(key, Duration::ZERO), WaitOutcome::TimedOut, "no wake");
        assert_eq!(liveness.total(), 1, "only the popping dispatch is alive");
        let s = stats.snapshot().steps[1].1;
        assert_eq!((s.progress_count, s.no_progress_count), (1, 1));
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
                &WakePlan::legacy(None),
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
            &WakePlan::legacy(None),
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
            &WakePlan::legacy(None),
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
            &WakePlan::legacy(None),
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
        let mut worker = WorkerCore::new(0, Some(src));

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
            &WakePlan::legacy(None),
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
        let mut worker = WorkerCore::new(0, Some(src));

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
            &WakePlan::legacy(None),
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
            &WakePlan::legacy(None),
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
        let mut worker = WorkerCore::new(0, Some(src));
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
            &WakePlan::legacy(None),
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
            WalkDirection::Reverse,
            vec![crate::runtime::scheduler::RefillSource {
                signal: Arc::new(AtomicBool::new(true)),
                refill_through,
                refill_branch: BranchIdx(0),
                cap_bytes: u64::MAX,
            }],
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

        let mut worker = WorkerCore::new(0, None);
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
            &WakePlan::legacy(None),
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
            &WakePlan::legacy(None),
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
            &WakePlan::legacy(None),
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

    /// Holds once (marks the thread), makes progress in that same pass, then is
    /// idle; on its third call records how many workers are armed on the
    /// event-count — the caller's own arm, if its idle went there.
    struct HeldThenIdle {
        calls: u32,
        plan: Arc<WakePlan>,
        armed_seen: Arc<Mutex<Option<usize>>>,
    }
    impl Step for HeldThenIdle {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "HeldThenIdle",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            self.calls += 1;
            match self.calls {
                1 => {
                    crate::runtime::wake_slot::note_held();
                    Ok(StepOutcome::Progress)
                }
                2 => Ok(StepOutcome::NoProgress),
                _ => {
                    *self.armed_seen.lock() = self.plan.pool().map(PoolEventCount::waiters);
                    Ok(StepOutcome::Finished)
                }
            }
        }
    }

    /// The holding hint describes one iteration's passes: a worker marked as
    /// holding in a pass that also made progress, and idle with nothing held on
    /// a later iteration, idles on the event-count (it is armed there when its
    /// re-poll runs), not on a timer park as if it still held an item.
    #[test]
    fn a_stale_holding_hint_does_not_reach_a_later_idle() {
        let (steps, graph, _runs) = three_step_chain(StepKind::Parallel);
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let (plan, _g) = crate::runtime::wake::tests_support::directed_plan_with_workers(2);
        let armed_seen = Arc::new(Mutex::new(None));
        let mut row: Vec<WorkerStepEntry> = (0..3).map(|_| WorkerStepEntry::Skip).collect();
        row[0] = WorkerStepEntry::Owned {
            step: Box::new(TypedStep::new(HeldThenIdle {
                calls: 0,
                plan: Arc::clone(&plan),
                armed_seen: Arc::clone(&armed_seen),
            })),
        };
        let drain_counters: Vec<Arc<StepDrainCounter>> =
            (0..3).map(|_| StepDrainCounter::new(1)).collect();
        let _ = crate::runtime::wake_slot::take_rejected();
        let mut worker = WorkerCore::new(0, None);
        run_worker_loop(
            &mut worker,
            &mut row,
            &contexts,
            &drain_counters,
            &PipelineSignal::new(),
            None,
            &crate::liveness::LivenessCounter::new(1),
            &crate::runtime::scheduler::ChainOrderScheduler,
            None,
            0,
            &plan,
            false,
        );
        assert_eq!(*armed_seen.lock(), Some(1), "the idle re-poll ran armed on the event-count");
    }

    /// A capped step that notes, at the start of each call, whether its
    /// thread's (slot 1) refusal record on its cap is set. Call 1 is then
    /// refused (recording slot 1) and frees one of the test's permits, whose
    /// release wakes slot 0 (recorded lower) and leaves slot 1's bit standing;
    /// it reports `Capped`. Any later call finishes.
    struct ObservesOwnRecord {
        cap: &'static Arc<crate::admission::PhaseCap>,
        held: Arc<Mutex<Vec<crate::admission::CapPermit<'static>>>>,
        seen: Arc<Mutex<Vec<bool>>>,
    }
    impl Step for ObservesOwnRecord {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ObservesOwnRecord",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn phase_cap(&self) -> Option<&crate::admission::PhaseCap> {
            Some(self.cap)
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let first = {
                let mut seen = self.seen.lock();
                seen.push(self.cap.is_recorded(1));
                seen.len() == 1
            };
            if first {
                assert!(self.cap.try_acquire().is_none(), "the cap is full");
                drop(self.held.lock().pop());
                return Ok(StepOutcome::Capped);
            }
            Ok(StepOutcome::Finished)
        }
    }

    /// A cap's record of a thread covers the pass that made it: the worker loop
    /// clears the thread's own bit when its main pass starts (a record an
    /// earlier pass left) and again when its cap-parked re-poll starts (the
    /// main pass's refusal, which the re-poll re-records if still refused).
    /// The third site, the event-count re-poll, is
    /// `the_event_count_re_poll_starts_a_clean_pass`.
    /// Worker 1 (slot 1) of two; its one live step is capped by a tracked cap
    /// of 2, both permits held.
    #[test]
    fn every_pass_starts_with_the_thread_s_own_records_cleared() {
        use crate::admission::PhaseCap;
        use crate::runtime::wake::{DriverIdx, WakeEdges};
        let cap: &'static Arc<PhaseCap> = Box::leak(Box::new(PhaseCap::new("t", 2)));
        cap.bind_waker(Some(Arc::new(|_: Option<usize>| {})));
        let held: Vec<_> = (0..2).map(|i| cap.try_acquire_as(10 + i).expect("room")).collect();
        // Records left by an earlier pass: slot 1's, and slot 0's (lower, so a
        // release claims it first).
        assert!(cap.try_acquire_as(1).is_none() && cap.try_acquire_as(0).is_none());
        let (steps, graph, _runs) = three_step_chain(StepKind::Parallel);
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let plan = WakePlan::build(
            &graph,
            &[StepKind::Detached, StepKind::Parallel, StepKind::Serial],
            &[None, None, None],
            &[Some(DriverIdx(0)), None, None],
            &[None, Some(Arc::clone(cap)), None],
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        let seen = Arc::new(Mutex::new(Vec::new()));
        let mut row: Vec<WorkerStepEntry> = (0..3).map(|_| WorkerStepEntry::Skip).collect();
        row[1] = WorkerStepEntry::Owned {
            step: Box::new(TypedStep::new(ObservesOwnRecord {
                cap,
                held: Arc::new(Mutex::new(held)),
                seen: Arc::clone(&seen),
            })),
        };
        let drain_counters: Vec<Arc<StepDrainCounter>> =
            (0..3).map(|_| StepDrainCounter::new(1)).collect();
        let _slot = crate::runtime::wake_slot::SlotGuard::enter(1);
        let mut worker = WorkerCore::new(1, None);
        run_worker_loop(
            &mut worker,
            &mut row,
            &contexts,
            &drain_counters,
            &PipelineSignal::new(),
            None,
            &LivenessCounter::new(2),
            &crate::runtime::scheduler::ChainOrderScheduler,
            None,
            1,
            &plan,
            false,
        );
        assert_eq!(
            *seen.lock(),
            vec![false, false],
            "cleared at the main pass's start, then at the re-poll's"
        );
    }

    /// Notes, on each call, whether slot 1's bit on its cap is set. Call 1
    /// records slot 1 the way no refusal of this thread does — by number,
    /// without marking the thread refused — and reports `NoProgress`, so the
    /// thread idles on the event-count; any later call finishes.
    struct RecordsWithoutRefusal {
        cap: &'static Arc<crate::admission::PhaseCap>,
        seen: Arc<Mutex<Vec<bool>>>,
    }
    impl Step for RecordsWithoutRefusal {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "RecordsWithoutRefusal",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let first = {
                let mut seen = self.seen.lock();
                seen.push(self.cap.is_recorded(1));
                seen.len() == 1
            };
            if first {
                assert!(self.cap.try_acquire_as(1).is_none(), "the cap is full");
                return Ok(StepOutcome::NoProgress);
            }
            Ok(StepOutcome::Finished)
        }
    }

    /// The event-count re-poll is a pass of its own and starts clean too (the
    /// third pass-start site). No refusal of the thread can leave a record
    /// there — it would send the thread to the timer park — so the step plants
    /// one by number. Worker 1 of two, unpinned, so it idles on the event-count.
    #[test]
    fn the_event_count_re_poll_starts_a_clean_pass() {
        use crate::admission::PhaseCap;
        use crate::runtime::wake::{DriverIdx, WakeEdges};
        let cap: &'static Arc<PhaseCap> = Box::leak(Box::new(PhaseCap::new("t", 1)));
        cap.bind_waker(Some(Arc::new(|_: Option<usize>| {})));
        let held = cap.try_acquire_as(10).expect("room");
        let (steps, graph, _runs) = three_step_chain(StepKind::Parallel);
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let plan = WakePlan::build(
            &graph,
            &[StepKind::Detached, StepKind::Parallel, StepKind::Serial],
            &[None, None, None],
            &[Some(DriverIdx(0)), None, None],
            &[None, Some(Arc::clone(cap)), None],
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        let seen = Arc::new(Mutex::new(Vec::new()));
        let mut row: Vec<WorkerStepEntry> = (0..3).map(|_| WorkerStepEntry::Skip).collect();
        row[1] = WorkerStepEntry::Owned {
            step: Box::new(TypedStep::new(RecordsWithoutRefusal { cap, seen: Arc::clone(&seen) })),
        };
        let drain_counters: Vec<Arc<StepDrainCounter>> =
            (0..3).map(|_| StepDrainCounter::new(1)).collect();
        let _slot = crate::runtime::wake_slot::SlotGuard::enter(1);
        let mut worker = WorkerCore::new(1, None);
        run_worker_loop(
            &mut worker,
            &mut row,
            &contexts,
            &drain_counters,
            &PipelineSignal::new(),
            None,
            &LivenessCounter::new(2),
            &crate::runtime::scheduler::ChainOrderScheduler,
            None,
            1,
            &plan,
            false,
        );
        assert_eq!(*seen.lock(), vec![false, false], "cleared when the re-poll started");
        drop(held);
    }

    /// A capped step that notes, on each call, whether its worker (1) is armed
    /// as cap-parked — that is, whether the call is the cap-parked re-poll.
    /// Call 1 is refused by its full cap and reports `Capped`; any later call
    /// finishes.
    struct NotesCapParked {
        cap: &'static Arc<crate::admission::PhaseCap>,
        plan: Arc<WakePlan>,
        seen: Arc<Mutex<Vec<bool>>>,
    }
    impl Step for NotesCapParked {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "NotesCapParked",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn phase_cap(&self) -> Option<&crate::admission::PhaseCap> {
            Some(self.cap)
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let first = {
                let mut seen = self.seen.lock();
                seen.push(self.plan.is_cap_parked(1));
                seen.len() == 1
            };
            if first {
                assert!(self.cap.try_acquire().is_none(), "the cap is full");
                return Ok(StepOutcome::Capped);
            }
            Ok(StepOutcome::Finished)
        }
    }

    /// The cap-parked re-poll skips a step whose cap is full: it would only be
    /// refused again (the skip records the worker on the cap instead). So the
    /// step's next call comes from the main pass after the timer park, not
    /// from the re-poll. Worker 1 of two, one live step capped by a full,
    /// tracked cap.
    #[test]
    fn a_cap_parked_re_poll_skips_a_full_cap() {
        use crate::admission::PhaseCap;
        use crate::runtime::wake::{DriverIdx, WakeEdges};
        let cap: &'static Arc<PhaseCap> = Box::leak(Box::new(PhaseCap::new("t", 1)));
        cap.bind_waker(Some(Arc::new(|_: Option<usize>| {})));
        let held = cap.try_acquire_as(10).expect("room");
        let (steps, graph, _runs) = three_step_chain(StepKind::Parallel);
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let plan = WakePlan::build(
            &graph,
            &[StepKind::Detached, StepKind::Parallel, StepKind::Serial],
            &[None, None, None],
            &[Some(DriverIdx(0)), None, None],
            &[None, Some(Arc::clone(cap)), None],
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        let seen = Arc::new(Mutex::new(Vec::new()));
        let mut row: Vec<WorkerStepEntry> = (0..3).map(|_| WorkerStepEntry::Skip).collect();
        row[1] = WorkerStepEntry::Owned {
            step: Box::new(TypedStep::new(NotesCapParked {
                cap,
                plan: Arc::clone(&plan),
                seen: Arc::clone(&seen),
            })),
        };
        let drain_counters: Vec<Arc<StepDrainCounter>> =
            (0..3).map(|_| StepDrainCounter::new(1)).collect();
        let _slot = crate::runtime::wake_slot::SlotGuard::enter(1);
        let mut worker = WorkerCore::new(1, None);
        run_worker_loop(
            &mut worker,
            &mut row,
            &contexts,
            &drain_counters,
            &PipelineSignal::new(),
            None,
            &LivenessCounter::new(2),
            &crate::runtime::scheduler::ChainOrderScheduler,
            None,
            1,
            &plan,
            false,
        );
        assert_eq!(*seen.lock(), vec![false, false], "the re-poll did not call the step");
        assert_eq!(cap.refused(), 1, "only the main pass's call was refused");
        drop(held);
    }

    /// A capped step for the cap-only re-poll that leaves no record. Call 1 is
    /// refused by its full cap with its input non-empty; before returning it
    /// has the test drain the input and free the permit. Later calls go through
    /// `admit_input` again: an empty input takes no permit and records nothing;
    /// once admitted, the step pops its item and finishes.
    struct DrainedBeforeRepoll {
        cap: &'static Arc<crate::admission::PhaseCap>,
        calls: Arc<AtomicUsize>,
        refused: std::sync::mpsc::Sender<()>,
        drained: Arc<Mutex<std::sync::mpsc::Receiver<()>>>,
        /// Set once a later call's `admit_input` has seen the input empty: the
        /// re-poll ran and recorded nothing.
        saw_empty: Arc<AtomicBool>,
    }
    impl Step for DrainedBeforeRepoll {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "DrainedBeforeRepoll",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn phase_cap(&self) -> Option<&crate::admission::PhaseCap> {
            Some(self.cap)
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let call = self.calls.fetch_add(1, Ordering::SeqCst) + 1;
            let admission = crate::admission::admit_input(ctx.input, Some(self.cap));
            if call == 1 {
                assert!(matches!(admission, Err(StepOutcome::Capped)), "refused with input");
                self.refused.send(()).unwrap();
                self.drained.lock().recv().unwrap();
                return Ok(StepOutcome::Capped);
            }
            let _permit = match admission {
                Ok(p) => p,
                Err(outcome) => {
                    if outcome == StepOutcome::NoProgress {
                        self.saw_empty.store(true, Ordering::SeqCst);
                    }
                    return Ok(outcome);
                }
            };
            assert!(ctx.input.pop().is_some(), "admitted with input");
            Ok(StepOutcome::Finished)
        }
    }

    /// A cap-only worker whose re-poll leaves no cap record (the refusing cap
    /// had room by then and its step nothing to take) parks as an ordinary
    /// idle worker, where a `Pool` wake reaches it: a new item that arrives
    /// while the cap is full again is delivered without its timer. Parked
    /// cap-only, it could be reached only by that cap's release, and no record
    /// is left for the release to find. Worker 1 of two has a 10 s timer; the
    /// run must end within 5 s.
    #[test]
    fn a_cap_only_worker_with_no_record_left_parks_where_pool_wakes_reach_it() {
        use std::time::Duration;

        use crate::admission::PhaseCap;
        use crate::runtime::wake::PushSnapshot;
        use crate::step::OutputHandles;
        let cap: &'static Arc<PhaseCap> = Box::leak(Box::new(PhaseCap::new("t", 1)));
        let (contexts, plan) = capped_middle_plan(cap, 2);
        cap.bind_waker(Some({
            let plan = Arc::clone(&plan);
            Arc::new(move |slot: Option<usize>| match slot {
                Some(s) => plan.deliver_slot(s),
                None => plan.deliver_anonymous(),
            })
        }));
        let producer = |c: &ChainContexts| {
            assert!(
                c.outputs[0].downcast_ref::<OutputHandles<Single<u32>>>().unwrap().push(7).is_ok()
            );
        };
        producer(&contexts);
        let held = cap.try_acquire_as(9).expect("room");
        // Slot 0 is recorded too, below slot 1, so the release that frees the
        // permit before the re-poll wakes slot 0 and leaves slot 1's bit for
        // its re-poll to clear.
        assert!(cap.try_acquire_as(0).is_none());
        let calls = Arc::new(AtomicUsize::new(0));
        let saw_empty = Arc::new(AtomicBool::new(false));
        let (refused_tx, refused_rx) = std::sync::mpsc::channel();
        let (drained_tx, drained_rx) = std::sync::mpsc::channel();
        let (done_tx, done_rx) = std::sync::mpsc::channel();
        let worker = {
            let (contexts, plan, calls, saw_empty) = (
                Arc::clone(&contexts),
                Arc::clone(&plan),
                Arc::clone(&calls),
                Arc::clone(&saw_empty),
            );
            std::thread::spawn(move || {
                let step = Box::new(TypedStep::new(DrainedBeforeRepoll {
                    cap,
                    calls,
                    refused: refused_tx,
                    drained: Arc::new(Mutex::new(drained_rx)),
                    saw_empty,
                }));
                let signal = PipelineSignal::new();
                run_middle_step(1, step, &contexts, &plan, &signal, Some(10_000_000));
                let _ = done_tx.send(());
            })
        };
        refused_rx.recv().unwrap();
        let input = contexts.inputs[1].downcast_ref::<BranchInputHandle<u32>>().unwrap();
        assert_eq!(input.pop(), Some(7), "drain the item before the re-poll");
        drop(held); // wakes slot 0, not the worker
        drained_tx.send(()).unwrap();
        // The worker's re-poll has seen the input empty (and so recorded
        // nothing); after it the worker either parks cap-only (and no later
        // release finds it) or, once demoted, idles where `Pool` wakes reach it.
        let deadline = std::time::Instant::now() + Duration::from_secs(5);
        while !saw_empty.load(Ordering::SeqCst) {
            assert!(std::time::Instant::now() < deadline, "the re-poll never ran");
            std::thread::yield_now();
        }
        let full = cap.try_acquire_as(9).expect("the cap is free");
        producer(&contexts);
        let _ = plan.on_progress(StepIdx(0), &PushSnapshot::default());
        drop(full);
        let finished = done_rx.recv_timeout(Duration::from_secs(5)).is_ok();
        assert!(finished, "the new item waited for the worker's 10 s timer");
        worker.join().unwrap();
        assert_eq!(plan.pool().map(PoolEventCount::waiters), Some(0));
    }

    /// A three-step chain's contexts and a Directed plan for it at `n_threads`
    /// pool workers, the middle (`Parallel`) step capped by `cap`.
    fn capped_middle_plan(
        cap: &Arc<crate::admission::PhaseCap>,
        n_threads: usize,
    ) -> (Arc<ChainContexts>, Arc<WakePlan>) {
        use crate::runtime::wake::{DriverIdx, WakeEdges};
        let (steps, graph, _runs) = three_step_chain(StepKind::Parallel);
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let plan = WakePlan::build(
            &graph,
            &[StepKind::Detached, StepKind::Parallel, StepKind::Serial],
            &[None, None, None],
            &[Some(DriverIdx(0)), None, None],
            &[None, Some(Arc::clone(cap)), None],
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(n_threads))),
            n_threads,
        );
        (contexts, plan)
    }

    /// Run worker `slot`'s loop over just the middle step `step`.
    fn run_middle_step(
        slot: usize,
        step: Box<dyn ErasedStep>,
        contexts: &Arc<ChainContexts>,
        plan: &WakePlan,
        signal: &Arc<PipelineSignal>,
        backoff_us: Option<u64>,
    ) {
        let mut row: Vec<WorkerStepEntry> = (0..3).map(|_| WorkerStepEntry::Skip).collect();
        row[1] = WorkerStepEntry::Owned { step };
        let drain_counters: Vec<Arc<StepDrainCounter>> =
            (0..3).map(|_| StepDrainCounter::new(1)).collect();
        let _slot = crate::runtime::wake_slot::SlotGuard::enter(slot);
        plan.register_worker(slot, std::thread::current());
        let mut worker = WorkerCore::new(slot, None).with_backoff_override(backoff_us);
        run_worker_loop(
            &mut worker,
            &mut row,
            contexts,
            &drain_counters,
            signal,
            None,
            &LivenessCounter::new(slot + 1),
            &crate::runtime::scheduler::ChainOrderScheduler,
            None,
            slot,
            plan,
            false,
        );
    }

    /// A capped step whose cap stays full and whose input stays non-empty:
    /// every call is refused. Counts its calls.
    struct AlwaysRefused {
        cap: &'static Arc<crate::admission::PhaseCap>,
        calls: Arc<AtomicUsize>,
    }
    impl Step for AlwaysRefused {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "AlwaysRefused",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn phase_cap(&self) -> Option<&crate::admission::PhaseCap> {
            Some(self.cap)
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            self.calls.fetch_add(1, Ordering::SeqCst);
            match crate::admission::admit_input(ctx.input, Some(self.cap)) {
                Err(outcome) => Ok(outcome),
                Ok(_) => unreachable!("the cap stays full"),
            }
        }
    }

    /// A worker whose record has no bit in its record mask — on the shared
    /// anonymous bit, past the cap's numbered slots — is still recorded, so it
    /// parks cap-only like any refused worker instead of being demoted to
    /// another pass every round. With its cap full and input on its step, it
    /// runs its step, then sits armed in the cap-parked set: its step is called
    /// a handful of times in 200 ms, not hundreds of thousands.
    #[test]
    fn a_worker_recorded_outside_its_mask_parks() {
        use std::time::{Duration, Instant};

        use crate::admission::{PhaseCap, REFUSED_SLOTS};
        use crate::step::OutputHandles;
        let slot = REFUSED_SLOTS + 45;
        let cap: &'static Arc<PhaseCap> = Box::leak(Box::new(PhaseCap::new("t", 1)));
        let (contexts, plan) = capped_middle_plan(cap, slot + 1);
        cap.bind_waker(Some(Arc::new(|_: Option<usize>| {})));
        let producer = contexts.outputs[0].downcast_ref::<OutputHandles<Single<u32>>>().unwrap();
        assert!(producer.push(7).is_ok());
        let held = cap.try_acquire_as(9).expect("room");
        let calls = Arc::new(AtomicUsize::new(0));
        let signal = PipelineSignal::new();
        let worker = {
            let (contexts, plan, calls, signal) =
                (Arc::clone(&contexts), Arc::clone(&plan), Arc::clone(&calls), Arc::clone(&signal));
            std::thread::spawn(move || {
                let step = Box::new(TypedStep::new(AlwaysRefused { cap, calls }));
                run_middle_step(slot, step, &contexts, &plan, &signal, Some(10_000_000));
            })
        };
        let deadline = Instant::now() + Duration::from_secs(5);
        while calls.load(Ordering::SeqCst) == 0 {
            assert!(Instant::now() < deadline, "the worker never ran its step");
            std::thread::yield_now();
        }
        std::thread::sleep(Duration::from_millis(200));
        let n = calls.load(Ordering::SeqCst);
        let parked = plan.is_cap_parked(slot);
        signal.cancel();
        plan.unpark_all();
        worker.join().unwrap();
        drop(held);
        assert!(n < 50, "{n} calls in 200 ms: the worker never parked");
        assert!(parked, "the worker sits armed in the cap-parked set");
    }

    /// How a [`ClaimThenIdle`] step arranges for its worker to be handed a
    /// cap wake it does not spend, and so which pass's end must forward it.
    #[derive(Clone, Copy, Debug)]
    enum ClaimPath {
        /// Refused in a pass that also makes progress: the next main pass
        /// owes the wake.
        MainPass,
        /// Recorded by a skip without a refusal, then idle: the event-count
        /// re-poll owes it. Synthetic: in a run only the cap-parked re-poll
        /// skips, and a refusal sends the thread to the timer park, so no
        /// record reaches an event-count re-poll; the step plants one to pin
        /// that pass's end.
        EventCountRepoll,
        /// Refused, then idle: the cap-parked re-poll owes it.
        CapParkedRepoll,
        /// Recorded by a skip in its last call: the loop's exit owes it.
        Exit,
    }

    /// Call 1 records this thread (slot 1) on its full cap — by a refusal or a
    /// full-cap skip, per `path` — then frees the test's permit, whose release
    /// claims slot 1's bit (the lowest recorded): the thread is handed that
    /// wake. Call 2 takes no permit and finishes. Slot 2, also recorded, is
    /// owed the wake: only the end of the pass that took no permit forwards it.
    struct ClaimThenIdle {
        cap: &'static Arc<crate::admission::PhaseCap>,
        held: Mutex<Option<crate::admission::CapPermit<'static>>>,
        path: ClaimPath,
        calls: u32,
    }
    impl Step for ClaimThenIdle {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ClaimThenIdle",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn phase_cap(&self) -> Option<&crate::admission::PhaseCap> {
            Some(self.cap)
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            self.calls += 1;
            if self.calls > 1 {
                return Ok(StepOutcome::Finished);
            }
            match self.path {
                ClaimPath::MainPass | ClaimPath::CapParkedRepoll => {
                    assert!(self.cap.try_acquire().is_none(), "the cap is full");
                }
                ClaimPath::EventCountRepoll | ClaimPath::Exit => {
                    assert!(self.cap.record_skip_current(), "the cap is full");
                }
            }
            drop(self.held.lock().take());
            Ok(match self.path {
                ClaimPath::MainPass => StepOutcome::Progress,
                ClaimPath::EventCountRepoll => StepOutcome::NoProgress,
                ClaimPath::CapParkedRepoll => StepOutcome::Capped,
                ClaimPath::Exit => StepOutcome::Finished,
            })
        }
    }

    /// Every pass of the worker loop ends by forwarding a cap wake its thread
    /// was handed and did not spend: the main pass, the event-count re-poll,
    /// the cap-parked re-poll, and the loop's exit. Each case hands the wake
    /// to the thread just before a different one of those passes.
    #[rstest::rstest]
    #[case::main_pass(ClaimPath::MainPass)]
    #[case::event_count_repoll(ClaimPath::EventCountRepoll)]
    #[case::cap_parked_repoll(ClaimPath::CapParkedRepoll)]
    #[case::exit(ClaimPath::Exit)]
    fn every_pass_end_forwards_an_unspent_cap_wake(#[case] path: ClaimPath) {
        use crate::admission::PhaseCap;
        let cap: &'static Arc<PhaseCap> = Box::leak(Box::new(PhaseCap::new("t", 1)));
        let (contexts, plan) = capped_middle_plan(cap, 3);
        let woken = crate::admission::bind_recording_waker(cap);
        let held = cap.try_acquire_as(9).expect("room");
        assert!(cap.try_acquire_as(2).is_none(), "slot 2 waits on the cap");
        let step = Box::new(TypedStep::new(ClaimThenIdle {
            cap,
            held: Mutex::new(Some(held)),
            path,
            calls: 0,
        }));
        run_middle_step(1, step, &contexts, &plan, &PipelineSignal::new(), None);
        assert_eq!(*woken.lock(), vec![Some(1), Some(2)], "slot 1's unspent wake reached slot 2");
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
                    &WakePlan::legacy(None),
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
            &WakePlan::legacy(Some(Arc::clone(&ec))),
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
            &WakePlan::legacy(Some(Arc::clone(&ec))),
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

    /// Slot-0 step for `an_event_count_exit_acknowledges_the_request`. Its
    /// first call is idle, so the worker arms the event-count and re-polls. On
    /// the re-poll it requests a pool worker (this worker is the armed waiter,
    /// so the request is sent and becomes pending), then marks the thread as a
    /// refused push (`note_held`) or a recording cap refusal
    /// (`note_cap_refused`) would, so the worker leaves the event-count without
    /// waiting. Its third call finishes.
    struct RequestThenRefuse {
        calls: u32,
        cap_refused: bool,
        requested: Arc<Mutex<Option<crate::runtime::wake::PoolRequest>>>,
    }
    impl Step for RequestThenRefuse {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "RequestThenRefuse",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            self.calls += 1;
            match self.calls {
                1 => Ok(StepOutcome::NoProgress),
                2 => {
                    *self.requested.lock() = Some(ctx.pool.request_worker());
                    if self.cap_refused {
                        crate::runtime::wake_slot::note_cap_refused();
                    } else {
                        crate::runtime::wake_slot::note_held();
                    }
                    Ok(StepOutcome::NoProgress)
                }
                _ => Ok(StepOutcome::Finished),
            }
        }
    }

    /// A worker that leaves the event-count after its re-poll was refused (it
    /// is holding, or a cap recorded it) acknowledges any pending pool request,
    /// as every other exit from the event-count does. Otherwise a request sent
    /// while it was armed stays pending, and every later `request_worker`
    /// reports `Pending` while timer-parked holders never acknowledge it. A
    /// second worker armed on the event-count is the request's real target
    /// (the requester itself is never one), and it is woken.
    #[rstest::rstest]
    #[case::holding(false)]
    #[case::cap_refused(true)]
    fn an_event_count_exit_acknowledges_the_request(#[case] cap_refused: bool) {
        use crate::runtime::event_count::WaitOutcome;
        use crate::runtime::wake::PoolRequest;
        let (steps, graph, _runs) = three_step_chain(StepKind::Parallel);
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let (plan, _g) = crate::runtime::wake::tests_support::directed_plan_with_workers(2);
        let requested = Arc::new(Mutex::new(None));
        let mut row: Vec<WorkerStepEntry> = (0..3).map(|_| WorkerStepEntry::Skip).collect();
        row[0] = WorkerStepEntry::Owned {
            step: Box::new(TypedStep::new(RequestThenRefuse {
                calls: 0,
                cap_refused,
                requested: Arc::clone(&requested),
            })),
        };
        let drain_counters: Vec<Arc<StepDrainCounter>> =
            (0..3).map(|_| StepDrainCounter::new(1)).collect();
        let _ = crate::runtime::wake_slot::take_rejected();
        let _ = crate::runtime::wake_slot::take_cap_refused();
        // The other worker, armed on the event-count for up to 10 s.
        let (armed_tx, armed_rx) = std::sync::mpsc::channel::<()>();
        let other = {
            let plan = Arc::clone(&plan);
            std::thread::spawn(move || {
                let pool = plan.pool().expect("t2 has a pool");
                let key = pool.prepare_wait();
                armed_tx.send(()).unwrap();
                let t = std::time::Instant::now();
                // No acknowledgement here: only the worker loop's refused exit
                // may clear the pending request this test is about.
                let outcome = pool.wait(key, std::time::Duration::from_secs(10));
                (outcome, t.elapsed())
            })
        };
        armed_rx.recv().unwrap();
        let _slot = crate::runtime::wake_slot::SlotGuard::enter(0);
        let mut worker = WorkerCore::new(0, None);
        run_worker_loop(
            &mut worker,
            &mut row,
            &contexts,
            &drain_counters,
            &PipelineSignal::new(),
            None,
            &crate::liveness::LivenessCounter::new(1),
            &crate::runtime::scheduler::ChainOrderScheduler,
            None,
            0,
            &plan,
            false,
        );
        assert_eq!(
            *requested.lock(),
            Some(PoolRequest::Woken),
            "the request was sent while this worker was armed, to the other waiter"
        );
        let (outcome, waited) = other.join().expect("other worker");
        assert_ne!(outcome, WaitOutcome::TimedOut, "the other worker was woken by the request");
        assert!(waited < std::time::Duration::from_secs(5), "...promptly: {waited:?}");
        // A new request with an armed waiter must be sent, not reported pending.
        let (ready_tx, ready_rx) = std::sync::mpsc::channel::<()>();
        let plan_w = Arc::clone(&plan);
        let waiter = std::thread::spawn(move || {
            let pool = plan_w.pool().expect("t2 has a pool");
            let key = pool.prepare_wait();
            ready_tx.send(()).unwrap();
            pool.wait(key, std::time::Duration::from_secs(10))
        });
        ready_rx.recv().unwrap();
        let again = plan.request_pool_worker();
        let _ = plan.pool().expect("t2 has a pool").notify_all(); // release the waiter regardless
        let _ = waiter.join();
        assert_eq!(again, PoolRequest::Woken, "the earlier request was acknowledged");
    }

    /// Slot-0 step for `a_re_polling_worker_is_not_its_own_request_target`:
    /// idle on its first call (the worker arms the event-count and re-polls),
    /// requests a pool worker twice on the re-poll, then finishes.
    struct RequestOnRepoll {
        calls: u32,
        requested: Arc<Mutex<Vec<crate::runtime::wake::PoolRequest>>>,
    }
    impl Step for RequestOnRepoll {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "RequestOnRepoll",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            self.calls += 1;
            match self.calls {
                1 => Ok(StepOutcome::NoProgress),
                2 => {
                    let mut requested = self.requested.lock();
                    requested.push(ctx.pool.request_worker());
                    requested.push(ctx.pool.request_worker());
                    Ok(StepOutcome::NoProgress)
                }
                _ => Ok(StepOutcome::Finished),
            }
        }
    }

    /// The worker loop marks itself armed (`set_ec_armed`) around its
    /// event-count re-poll, so a step that requests a worker during that
    /// re-poll does not count its own worker as the parked waiter to wake.
    /// Here the requester is the only event-count waiter and another worker is
    /// armed for a timer park: the first request must claim that worker, not
    /// notify the requester's own event-count wait. The second request tells
    /// the two apart: with worker 1's bit claimed nobody else is parked
    /// (`AllAwake`); a request sent to the requester's own wait would instead
    /// leave a pending event-count request (`Pending`).
    #[test]
    fn a_re_polling_worker_is_not_its_own_request_target() {
        use crate::runtime::wake::PoolRequest;
        let (steps, graph, _runs) = three_step_chain(StepKind::Parallel);
        let contexts = Arc::new(build_chain_contexts(
            &steps,
            &graph,
            crate::builder::InstrumentationLevel::Off,
            false,
        ));
        let (plan, _g) = crate::runtime::wake::tests_support::directed_plan_with_workers(2);
        let requested = Arc::new(Mutex::new(Vec::new()));
        let mut row: Vec<WorkerStepEntry> = (0..3).map(|_| WorkerStepEntry::Skip).collect();
        row[0] = WorkerStepEntry::Owned {
            step: Box::new(TypedStep::new(RequestOnRepoll {
                calls: 0,
                requested: Arc::clone(&requested),
            })),
        };
        let drain_counters: Vec<Arc<StepDrainCounter>> =
            (0..3).map(|_| StepDrainCounter::new(1)).collect();
        let _ = crate::runtime::wake_slot::take_rejected();
        let _ = crate::runtime::wake_slot::take_cap_refused();
        // Worker 1, armed for a timer park for up to 10 s.
        let (armed_tx, armed_rx) = std::sync::mpsc::channel::<()>();
        let other = {
            let plan = Arc::clone(&plan);
            std::thread::spawn(move || {
                plan.register_worker(1, std::thread::current());
                assert!(plan.arm_direct(1));
                armed_tx.send(()).unwrap();
                let t = std::time::Instant::now();
                std::thread::park_timeout(std::time::Duration::from_secs(10));
                plan.disarm_direct(1);
                t.elapsed()
            })
        };
        armed_rx.recv().unwrap();
        let _slot = crate::runtime::wake_slot::SlotGuard::enter(0);
        let mut worker = WorkerCore::new(0, None);
        run_worker_loop(
            &mut worker,
            &mut row,
            &contexts,
            &drain_counters,
            &PipelineSignal::new(),
            None,
            &crate::liveness::LivenessCounter::new(1),
            &crate::runtime::scheduler::ChainOrderScheduler,
            None,
            0,
            &plan,
            false,
        );
        assert_eq!(
            *requested.lock(),
            vec![PoolRequest::Woken, PoolRequest::AllAwake],
            "the first request claimed worker 1; nothing was left pending"
        );
        let parked_for = other.join().expect("worker 1");
        assert!(
            parked_for < std::time::Duration::from_secs(5),
            "the request reached the timer-parked worker, not the requester: {parked_for:?}"
        );
    }

    /// Slot-0 stand-in for a flush-first step: marks the dispatch flushed (as a
    /// successful `BranchOutputHandle::retry` does; the `handles.rs` table pins
    /// the real setter) and/or cap-refused (as a recording `PhaseCap` refusal
    /// does), then returns `result` (`None` → `Err`).
    struct FlushThen {
        flush: bool,
        cap_refused: bool,
        result: Option<StepOutcome>,
    }
    impl Step for FlushThen {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "FlushThen",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if self.flush {
                crate::runtime::wake_slot::note_flushed();
            }
            if self.cap_refused {
                crate::runtime::wake_slot::note_cap_refused();
            }
            self.result.ok_or_else(|| io::Error::other("the row's error"))
        }
    }

    /// Whether a `dispatch_runs_the_flushed_wake_for_every_idle_outcome` row
    /// checks the capped resume hint, and for which step.
    #[derive(Clone, Copy, Debug)]
    enum ResumeCheck {
        Unchecked,
        Is(Option<StepIdx>),
    }

    /// Clear every thread-local the flushed-wake rows read, so no row depends
    /// on what an earlier test left on a reused thread.
    fn clear_dispatch_flags() {
        let _ = crate::runtime::wake_slot::take_flushed();
        let _ = crate::runtime::wake_slot::take_cap_refused();
        let _ = crate::runtime::wake_slot::take_resume_step();
    }

    /// The dispatcher runs the flushed push's forward wake for every idle
    /// outcome, and only for those: `Progress` already wakes (once, not
    /// twice), `Finished` broadcasts (unrecorded), `Err` ends the run, a clean
    /// idle dispatch wakes nobody, and every call leaves the flag clear (a
    /// stale flag is cleared even by a `Skip`). Step 0's forward target is
    /// `Driver(0)`, registered as this thread, so a wake is one `unpark`.
    #[rstest::rstest]
    #[case::progress(Some(true), false, Some(StepOutcome::Progress), 1, ResumeCheck::Unchecked)]
    #[case::no_progress_flushed(
        Some(true),
        false,
        Some(StepOutcome::NoProgress),
        1,
        ResumeCheck::Unchecked
    )]
    #[case::contention_flushed(
        Some(true),
        false,
        Some(StepOutcome::Contention),
        1,
        ResumeCheck::Unchecked
    )]
    #[case::capped_flushed(Some(true), false, Some(StepOutcome::Capped), 1, ResumeCheck::Is(None))]
    #[case::capped_flushed_with_cap_refused(
        Some(true),
        true,
        Some(StepOutcome::Capped),
        1,
        ResumeCheck::Is(Some(StepIdx(0)))
    )]
    #[case::no_progress_clean(
        Some(false),
        false,
        Some(StepOutcome::NoProgress),
        0,
        ResumeCheck::Unchecked
    )]
    #[case::contention_clean(
        Some(false),
        false,
        Some(StepOutcome::Contention),
        0,
        ResumeCheck::Unchecked
    )]
    #[case::capped_clean(Some(false), false, Some(StepOutcome::Capped), 0, ResumeCheck::Unchecked)]
    #[case::finished_flushed(
        Some(true),
        false,
        Some(StepOutcome::Finished),
        0,
        ResumeCheck::Unchecked
    )]
    #[case::err_flushed(Some(true), false, None, 0, ResumeCheck::Unchecked)]
    #[case::skip_with_a_stale_flag(None, false, None, 0, ResumeCheck::Unchecked)]
    fn dispatch_runs_the_flushed_wake_for_every_idle_outcome(
        // `None`: a `Skip` entry, with the flag pre-set after the clear.
        #[case] flush: Option<bool>,
        #[case] cap_refused: bool,
        #[case] result: Option<StepOutcome>,
        #[case] unparks: u64,
        #[case] resume: ResumeCheck,
    ) {
        use crate::runtime::wake::{DriverIdx, WakeEdges};
        use crate::runtime::wake_slot::{flushed_pending, note_flushed, take_resume_step};
        let (steps, graph, _runs) = three_step_chain(StepKind::Parallel);
        let contexts =
            build_chain_contexts(&steps, &graph, crate::builder::InstrumentationLevel::Off, false);
        let plan = WakePlan::build(
            &graph,
            &[StepKind::Parallel, StepKind::Detached, StepKind::Serial],
            &[None, None, None],
            &[None, Some(DriverIdx(0)), None],
            &[],
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        plan.register_driver(DriverIdx(0), std::thread::current());
        let stats = Arc::new(PipelineStats::new(vec!["Src", "Finish", "Sink"]));
        clear_dispatch_flags();
        let mut entry = if let Some(flush) = flush {
            WorkerStepEntry::Owned {
                step: Box::new(TypedStep::new(FlushThen { flush, cap_refused, result })),
            }
        } else {
            note_flushed(); // a stale flag from outside any dispatch
            WorkerStepEntry::Skip
        };
        let _ = dispatch_one_step(
            &mut entry,
            StepIdx(0),
            &contexts,
            &StepDrainCounter::new(2), // a `Finished` row never closes the shared output
            &PipelineSignal::new(),
            Some(&stats),
            &LivenessCounter::new(1),
            0,
            false,
            None,
            0,
            &plan,
            &mut LoopDiag::new(false),
        );
        assert_eq!(
            stats.snapshot().steps[0].1.unparks_issued,
            unparks,
            "step 0's recorded unparks"
        );
        assert!(!flushed_pending(), "every dispatch returns with the flag clear");
        if let ResumeCheck::Is(want) = resume {
            assert_eq!(take_resume_step(), want, "the resume hint is independent of the flush");
        }
        if unparks == 1 {
            // The counted unpark must have reached the registered driver (this
            // thread): its token makes a long park return at once.
            let t = std::time::Instant::now();
            std::thread::park_timeout(std::time::Duration::from_secs(10));
            assert!(
                t.elapsed() < std::time::Duration::from_secs(5),
                "the flushed wake's unpark never reached the driver"
            );
        } else {
            std::thread::park_timeout(std::time::Duration::ZERO);
        }
        clear_dispatch_flags();
    }

    /// A Legacy plan ignores a flush, exactly as before: no notify (observed
    /// through the event-count itself) and nothing recorded, for every idle
    /// outcome.
    #[rstest::rstest]
    #[case::no_progress(StepOutcome::NoProgress)]
    #[case::contention(StepOutcome::Contention)]
    #[case::capped(StepOutcome::Capped)]
    fn legacy_flushed_dispatch_notifies_nothing(#[case] result: StepOutcome) {
        use std::time::Duration;

        use crate::runtime::event_count::WaitOutcome;
        let (steps, graph, _runs) = three_step_chain(StepKind::Parallel);
        let contexts =
            build_chain_contexts(&steps, &graph, crate::builder::InstrumentationLevel::Off, false);
        let stats = Arc::new(PipelineStats::new(vec!["Src", "Finish", "Sink"]));
        let ec = Arc::new(PoolEventCount::new(2));
        clear_dispatch_flags();
        let mut entry = WorkerStepEntry::Owned {
            step: Box::new(TypedStep::new(FlushThen {
                flush: true,
                cap_refused: false,
                result: Some(result),
            })),
        };
        let key = ec.prepare_wait();
        let _ = dispatch_one_step(
            &mut entry,
            StepIdx(0),
            &contexts,
            &StepDrainCounter::new(2),
            &PipelineSignal::new(),
            Some(&stats),
            &LivenessCounter::new(1),
            0,
            false,
            None,
            0,
            &WakePlan::legacy(Some(Arc::clone(&ec))),
            &mut LoopDiag::new(false),
        );
        assert_eq!(ec.wait(key, Duration::ZERO), WaitOutcome::TimedOut, "no notify");
        assert_eq!(stats.snapshot().steps[0].1.notifies_issued, 0);
        clear_dispatch_flags();
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
            &WakePlan::legacy(None),
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
            &WakePlan::legacy(None),
            &mut LoopDiag::new(false),
        )
        .unwrap();
        assert!(matches!(info.result, Ok(StepOutcome::Finished)));
    }
}
