//! `run_detached_step`: dedicated-thread drive loop for a
//! [`StepKind::Detached`](crate::step::StepKind::Detached) step.
//!
//! A Detached step is excluded from the work-stealing pool
//! ([`build_worker_storage`](crate::runtime::build_worker_storage) gives every
//! pool worker a `Skip` entry for it) and instead runs on its own OS thread,
//! spawned at run start alongside the deadlock-monitor / queue-rebalancer and
//! joined after the workers. This mirrors the legacy sort's "N + 2" threading:
//! N pool workers do the parallel (compression-bound) work while the merge and
//! the output writer each occupy a dedicated thread, so neither steals a pool
//! worker slot.
//!
//! ## Drive model
//!
//! The dedicated thread drives the step through the **same** type-erased
//! `try_run_erased` path the pool dispatch uses
//! ([`dispatch_one_step`](crate::runtime::driver)'s `Owned`/`Exclusive` arm),
//! so the step body — its non-blocking `input.pop()` / `outputs.push()`
//! cooperative dance — behaves byte-for-byte identically to the pool-scheduled
//! version. The only difference is *who* drives it: a private loop on this
//! thread rather than the round-robin scheduler.
//!
//! On each iteration:
//!   - `try_run_erased` → `Progress`: loop immediately (reset backoff).
//!   - `Finished`: mark the step's output edges drained (counter-gated, init 1
//!     for Detached) and return — this is the shutdown signal for downstream
//!     consumers.
//!   - `NoProgress` / `Contention`: the step's input is momentarily empty (and
//!     not drained) or its output is full. Park with bounded exponential
//!     backoff, re-checking `signal.is_done()` so a cancel/error/peer-panic
//!     promptly tears the thread down.
//!   - `Err`: record the error on the signal and return.
//!
//! ## Lock-ordering acyclicity (non-negotiable, the deadlock proof)
//!
//! A Detached thread must NEVER hold one queue's internal lock while parking on
//! another. This loop upholds that by construction:
//!
//!   1. `try_run_erased` is a single non-blocking call. The step body inside it
//!      pops from its input transport (a lock-free `crossbeam ArrayQueue` —
//!      `try_pop`, no lock held across the call) and pushes to its output
//!      transports (`try_push`, likewise). It never *blocks* inside `try_run`;
//!      a full output or empty input is reported back as `NoProgress` /
//!      `Contention`, returning control to this loop.
//!   2. The only blocking this thread ever does is [`backoff_park`]'s
//!      `thread::sleep`, which holds NO queue lock whatsoever.
//!
//! So there is no cycle: this thread parks only *between* `try_run` calls,
//! never while holding a transport lock, so a two-sided Detached step (consuming
//! from the pool AND producing to the pool, both bounded) cannot deadlock — when
//! its output is full it yields (no input lock held) and the pool consumer
//! drains it; when its input is empty it yields (no output lock held) and the
//! pool producer fills it. The `detached_two_sided_no_deadlock` test pins this
//! with a wall-clock watchdog.

use std::any::Any;
use std::sync::Arc;
use std::time::{Duration, Instant};

use crate::erased::{ErasedStep, ErasedStepCtx};
use crate::runtime::contexts::ChainContexts;
use crate::runtime::drain::StepDrainCounter;
use crate::runtime::stats::PipelineStats;
use crate::signal::{PipelineError, PipelineSignal};
use crate::step::{Affinity, OutputsViewAny, StepKind, StepOutcome, StepProfile};
use crate::topology::StepIdx;

/// Sentinel left in the chain's `steps` vec in place of a `Detached` step after
/// its real instance has been extracted for its dedicated thread (see
/// [`extract_detached_steps`]). Keeping a same-position placeholder preserves
/// the `step_idx`-aligned indexing that `build_worker_storage` and
/// `ChainContexts` rely on. `build_worker_storage` only reads `kind()` (matches
/// `Detached` → every worker gets `Skip`) and then drops the box, so the
/// placeholder never has any other method invoked; they panic to catch a
/// framework bug if one ever is.
struct DetachedPlaceholder {
    name: &'static str,
}

impl ErasedStep for DetachedPlaceholder {
    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Detached,
            sticky: false,
            output_queues: Vec::new(),
            branch_ordering: Vec::new(),
        }
    }
    fn name(&self) -> &'static str {
        self.name
    }
    fn kind(&self) -> StepKind {
        StepKind::Detached
    }
    fn sticky(&self) -> bool {
        false
    }
    fn affinity(&self) -> Affinity {
        Affinity::None
    }
    fn try_run_erased(&mut self, _ctx: &mut ErasedStepCtx<'_>) -> std::io::Result<StepOutcome> {
        panic!("DetachedPlaceholder::try_run_erased invoked — the real instance was extracted");
    }
    fn clone_boxed(&self) -> Box<dyn ErasedStep> {
        panic!("DetachedPlaceholder::clone_boxed invoked — placeholder is never cloned");
    }
    fn build_input_handle(
        &self,
        _producer_set: &mut crate::handles::OutputQueueSet,
        _branch_idx: usize,
    ) -> Box<dyn Any + Send + Sync> {
        panic!("DetachedPlaceholder::build_input_handle invoked");
    }
    fn build_output_set(
        &self,
        _level: crate::builder::InstrumentationLevel,
    ) -> (crate::handles::OutputQueueSet, OutputsViewAny) {
        panic!("DetachedPlaceholder::build_output_set invoked");
    }
    fn build_fused_output_set(
        &self,
        _level: crate::builder::InstrumentationLevel,
    ) -> (crate::handles::OutputQueueSet, OutputsViewAny) {
        panic!("DetachedPlaceholder::build_fused_output_set invoked");
    }
    fn wrap_outputs_view(&self, _view: OutputsViewAny) -> Box<dyn Any + Send + Sync> {
        panic!("DetachedPlaceholder::wrap_outputs_view invoked");
    }
    fn mark_outputs_drained(&self, _outputs: &(dyn Any + Send + Sync)) {
        panic!("DetachedPlaceholder::mark_outputs_drained invoked");
    }
    fn is_source(&self) -> bool {
        false
    }
}

/// Remove every [`StepKind::Detached`](crate::step::StepKind::Detached) step's
/// real instance from `steps`, replacing each in place with a
/// [`DetachedPlaceholder`] so the surviving slots keep their `step_idx`
/// positions (which `build_worker_storage` and `ChainContexts` index by).
///
/// Returns `(step_idx, step)` pairs in chain order — the caller spawns one
/// dedicated thread per pair, driving it with [`run_detached_step`]. Called by
/// `Pipeline::run` **before** `build_worker_storage` consumes `steps`, while the
/// (read-only) `ChainContexts` have already been built from `&steps` (so the
/// extracted step's input/output handles live in `contexts[step_idx]`).
#[must_use]
pub fn extract_detached_steps(
    steps: &mut [Box<dyn ErasedStep>],
) -> Vec<(StepIdx, Box<dyn ErasedStep>)> {
    let mut detached = Vec::new();
    for (idx, slot) in steps.iter_mut().enumerate() {
        if slot.kind() == StepKind::Detached {
            let placeholder: Box<dyn ErasedStep> =
                Box::new(DetachedPlaceholder { name: slot.name() });
            let real = std::mem::replace(slot, placeholder);
            detached.push((StepIdx(idx), real));
        }
    }
    detached
}

/// Initial park duration after the first idle (empty input / full output) tick.
const INITIAL_PARK: Duration = Duration::from_micros(10);
/// Maximum park duration the backoff ramps up to. Bounds wake latency so a
/// Detached step resumes promptly once its peer makes room / supplies input,
/// while keeping an otherwise-idle dedicated thread off a hot spin.
const MAX_PARK: Duration = Duration::from_micros(500);

/// Bounded exponential-backoff park used by the Detached drive loop on an idle
/// tick. Returns the next park duration (doubled, capped at [`MAX_PARK`]).
/// Holds no queue lock — see the module-level lock-ordering proof.
fn backoff_park(current: Duration) -> Duration {
    std::thread::sleep(current);
    (current * 2).min(MAX_PARK)
}

/// Drive one [`StepKind::Detached`](crate::step::StepKind::Detached) step to
/// completion on the calling (dedicated) thread.
///
/// `step` is the single shared instance (Detached steps are never
/// `new_worker_copy`'d). `step_idx` indexes its handles in `contexts`.
/// `counter` is the per-step [`StepDrainCounter`] (init 1 for Detached), which
/// gates `mark_outputs_drained` exactly as for `Serial` / `Exclusive`.
///
/// Returns when the step reports `Finished`, records an `Err`, or the run is
/// cancelled (`signal.is_done()`). On `Finished` the step's output edges are
/// marked drained (the downstream consumer's end-of-stream signal); on a cancel
/// before `Finished` they are NOT closed (the run is tearing down and the
/// recorded error/cancel is what propagates).
pub fn run_detached_step(
    mut step: Box<dyn ErasedStep>,
    step_idx: StepIdx,
    contexts: &Arc<ChainContexts>,
    counter: &StepDrainCounter,
    signal: &Arc<PipelineSignal>,
    stats: Option<&Arc<PipelineStats>>,
) {
    let outputs_any = contexts.outputs[step_idx.0].as_ref();
    let mut park = INITIAL_PARK;

    loop {
        if signal.is_done() {
            // Cancelled / errored elsewhere — tear down without closing the
            // output (the recorded outcome is what propagates).
            return;
        }

        let start = stats.map(|_| Instant::now());
        let mut ctx = ErasedStepCtx {
            input: contexts.inputs[step_idx.0].as_ref(),
            outputs: outputs_any,
            signal,
        };
        let result = step.try_run_erased(&mut ctx);

        if let (Some(stats), Some(start)) = (stats, start) {
            let elapsed_ns = u64::try_from(start.elapsed().as_nanos()).unwrap_or(u64::MAX);
            let start_ns = stats.elapsed_ns().saturating_sub(elapsed_ns);
            match &result {
                Ok(outcome) => stats.record(step_idx, *outcome, start_ns, elapsed_ns),
                Err(_) => stats.record_error(step_idx, start_ns, elapsed_ns),
            }
            // Detached busy time: the `try_run` wall, attributed to this thread's
            // own (off-pool) busy account — never folded into the pool%.
            stats.record_detached_busy(step_idx, elapsed_ns);
        }

        match result {
            Ok(StepOutcome::Progress) => {
                park = INITIAL_PARK;
            }
            Ok(StepOutcome::Finished) => {
                // Single finisher (counter init 1) closes the output edges so
                // the downstream pool consumer observes drain and terminates.
                if counter.observe_drain() {
                    step.mark_outputs_drained(outputs_any);
                }
                return;
            }
            Ok(StepOutcome::NoProgress | StepOutcome::Contention) => {
                // Input momentarily empty (not drained) or output full. Park
                // (no lock held) so the pool peer can make progress, then
                // retry. Re-check `is_done` first to avoid a needless sleep on
                // a cancel that landed during `try_run`.
                if signal.is_done() {
                    return;
                }
                let sleep_start = stats.map(|_| Instant::now());
                park = backoff_park(park);
                if let (Some(stats), Some(ss)) = (stats, sleep_start) {
                    stats.record_detached_idle(
                        step_idx,
                        u64::try_from(ss.elapsed().as_nanos()).unwrap_or(u64::MAX),
                    );
                    stats.record_detached_park(step_idx);
                }
            }
            Err(io_err) => {
                signal.record_error(PipelineError::Io { step: step.name(), source: io_err });
                return;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io;
    use std::sync::atomic::{AtomicU32, Ordering};
    use std::time::Duration;

    use super::*;
    use crate::builder::InstrumentationLevel;
    use crate::erased::TypedStep;
    use crate::outputs::Single;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;
    use crate::step::{InputHandle, OutputHandles, Step, StepCtx, StepKind, StepProfile};

    /// `() -> u32` source stub used only so `build_output_set` constructs the
    /// transport that becomes the Detached step's input edge. Never run.
    #[derive(Clone)]
    struct SrcStub {
        capacity: usize,
    }
    impl Step for SrcStub {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Src",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: self.capacity }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    /// `u32 -> u32` pass-through Detached step: pop one item per `try_run`,
    /// push it on (holding it on output-full backpressure), report `Finished`
    /// once input drains and nothing is held.
    #[derive(Clone)]
    struct PassThroughDetached {
        held: Option<u32>,
    }
    impl Step for PassThroughDetached {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "PassThroughDetached",
                kind: StepKind::Detached,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(v) = self.held.take() {
                if ctx.outputs.push(v).is_err() {
                    self.held = Some(v);
                    return Ok(StepOutcome::Contention);
                }
                return Ok(StepOutcome::Progress);
            }
            match ctx.input.pop() {
                Some(v) => match ctx.outputs.push(v) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(unpushed) => {
                        self.held = Some(unpushed.into_item());
                        Ok(StepOutcome::Contention)
                    }
                },
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// Assemble a `Source -> PassThroughDetached -> Sink` shaped context by
    /// hand: build the producer's output set (= the Detached input edge), wire
    /// the Detached step's input from it, build the Detached step's own output
    /// set (= the downstream consumer's input edge), and return the pieces.
    #[allow(clippy::type_complexity)]
    fn build_one_detached(
        src_capacity: usize,
    ) -> (
        Box<dyn ErasedStep>,                       // the detached step
        Arc<ChainContexts>,                        // contexts for step_idx 1
        Arc<Box<dyn std::any::Any + Send + Sync>>, // producer outputs (push side)
        crate::handles::BranchInputHandle<u32>,    // downstream consumer (pop side)
    ) {
        let producer: Box<dyn ErasedStep> =
            Box::new(TypedStep::new(SrcStub { capacity: src_capacity }));
        let (mut producer_set, producer_view) =
            producer.build_output_set(InstrumentationLevel::Off);
        let producer_outputs_any = producer.wrap_outputs_view(producer_view);

        let det: Box<dyn ErasedStep> = Box::new(TypedStep::new(PassThroughDetached { held: None }));
        let det_input = det.build_input_handle(&mut producer_set, 0);
        let (mut det_set, det_view) = det.build_output_set(InstrumentationLevel::Off);
        let det_outputs_any = det.wrap_outputs_view(det_view);
        let det_output_consumer = det_set.take_typed_input::<u32>(0);

        let contexts = Arc::new(ChainContexts {
            inputs: vec![Box::new(()), det_input, Box::new(())],
            outputs: vec![Box::new(()), det_outputs_any, Box::new(())],
            bounded_queues: vec![],
            edges: vec![],
        });
        (det, contexts, Arc::new(producer_outputs_any), det_output_consumer)
    }

    /// Items pushed onto a Detached step's input flow through to its output;
    /// the step finishes once its input is drained and closes the output edge.
    #[test]
    fn detached_step_flows_items_and_finishes() {
        let (det, contexts, producer_outputs_any, consumer) = build_one_detached(64);
        let producer_outputs =
            producer_outputs_any.downcast_ref::<OutputHandles<Single<u32>>>().unwrap();
        producer_outputs.push(10).unwrap();
        producer_outputs.push(20).unwrap();
        producer_outputs.push(30).unwrap();
        producer_outputs.mark_all_drained();

        let counter = StepDrainCounter::new(1);
        let signal = PipelineSignal::new();
        run_detached_step(det, StepIdx(1), &contexts, &counter, &signal, None);

        let mut got = Vec::new();
        while let Some(v) = consumer.pop() {
            got.push(v);
        }
        assert_eq!(got, vec![10, 20, 30]);
        assert!(InputHandle::is_drained(&consumer), "output closed on Finished");
        assert!(!signal.is_done(), "clean completion, no error");
    }

    /// Zero items in (input drained from the start): the Detached step cleanly
    /// drains its output and returns.
    #[test]
    fn detached_step_zero_items_clean_drain() {
        let (det, contexts, producer_outputs_any, consumer) = build_one_detached(64);
        let producer_outputs =
            producer_outputs_any.downcast_ref::<OutputHandles<Single<u32>>>().unwrap();
        producer_outputs.mark_all_drained(); // no items

        let counter = StepDrainCounter::new(1);
        let signal = PipelineSignal::new();
        run_detached_step(det, StepIdx(1), &contexts, &counter, &signal, None);

        assert!(consumer.pop().is_none(), "no items produced");
        assert!(InputHandle::is_drained(&consumer), "output closed on clean empty drain");
    }

    /// A two-sided Detached step — consumer of a bounded pool-fed input AND
    /// producer to a bounded pool-drained output, both tiny — completes without
    /// deadlock. A wall-clock watchdog aborts (fails the test) if it wedges.
    /// This is the sort merge's topology.
    #[test]
    fn detached_two_sided_no_deadlock() {
        const N: u32 = 5_000;

        // cap-2 input AND cap-4 output both force interleaved backpressure.
        let (det, contexts, producer_outputs_any, consumer) = build_one_detached(2);
        let counter = StepDrainCounter::new(1);
        let signal = PipelineSignal::new();

        // Watchdog: a deadlock parks forever; abort so the test FAILS loudly.
        let done = Arc::new(std::sync::atomic::AtomicBool::new(false));
        {
            let done = Arc::clone(&done);
            std::thread::spawn(move || {
                for _ in 0..200 {
                    std::thread::sleep(Duration::from_millis(50));
                    if done.load(Ordering::SeqCst) {
                        return;
                    }
                }
                eprintln!("detached_two_sided_no_deadlock: WEDGED (deadlock)");
                std::process::abort();
            });
        }

        // Producer thread: push N items into the cap-2 input (backpressure),
        // then close it so the Detached step's input drains.
        let pushed = Arc::new(AtomicU32::new(0));
        let producer_handle = {
            let producer_outputs_any = Arc::clone(&producer_outputs_any);
            let pushed = Arc::clone(&pushed);
            std::thread::spawn(move || {
                let outputs =
                    producer_outputs_any.downcast_ref::<OutputHandles<Single<u32>>>().unwrap();
                let mut held: Option<u32> = None;
                let mut next = 0u32;
                loop {
                    if let Some(v) = held.take() {
                        match outputs.push(v) {
                            Ok(()) => {
                                pushed.fetch_add(1, Ordering::Relaxed);
                            }
                            Err(unpushed) => {
                                held = Some(unpushed.into_item());
                                std::thread::yield_now();
                            }
                        }
                        continue;
                    }
                    if next >= N {
                        break;
                    }
                    match outputs.push(next) {
                        Ok(()) => {
                            pushed.fetch_add(1, Ordering::Relaxed);
                            next += 1;
                        }
                        Err(unpushed) => {
                            // The value at `next` is now held for retry; advance
                            // `next` so the fresh-push branch doesn't re-emit it
                            // after `held` flushes (which would double-count).
                            held = Some(unpushed.into_item());
                            next += 1;
                            std::thread::yield_now();
                        }
                    }
                }
                outputs.mark_all_drained();
            })
        };

        // Consumer thread: pop everything the Detached step produces.
        let received = Arc::new(AtomicU32::new(0));
        let consumer_handle = {
            let received = Arc::clone(&received);
            std::thread::spawn(move || {
                loop {
                    if consumer.pop().is_some() {
                        received.fetch_add(1, Ordering::Relaxed);
                    } else if InputHandle::is_drained(&consumer) {
                        break;
                    } else {
                        std::thread::yield_now();
                    }
                }
            })
        };

        // Drive the Detached step on this thread to completion.
        run_detached_step(det, StepIdx(1), &contexts, &counter, &signal, None);

        producer_handle.join().unwrap();
        consumer_handle.join().unwrap();
        done.store(true, Ordering::SeqCst);

        assert_eq!(pushed.load(Ordering::Relaxed), N, "all items pushed");
        assert_eq!(
            received.load(Ordering::Relaxed),
            N,
            "all items flowed through the two-sided Detached step"
        );
    }
}
