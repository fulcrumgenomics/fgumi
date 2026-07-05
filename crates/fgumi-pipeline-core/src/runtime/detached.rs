//! Detached-step extraction and the dedicated **driver thread** that runs them
//! off the work-stealing pool.
//!
//! A [`StepKind::Detached`](crate::step::StepKind::Detached) step is excluded
//! from the pool ([`build_worker_storage`](crate::runtime::build_worker_storage)
//! gives every pool worker a `Skip` entry for it) and instead runs on a dedicated
//! OS thread spawned at run start alongside the deadlock-monitor /
//! queue-rebalancer and joined after the workers. This is the legacy sort's
//! "N + 2" threading: N pool workers do the parallel (compression-bound) work
//! while off-pool driver threads do the serial coordination + I/O, so neither
//! steals a pool worker slot.
//!
//! ## The driver IS a 1-thread pool ([`run_detached_driver`])
//!
//! A driver thread does **not** have a bespoke drive loop. It runs the *same*
//! [`run_worker_loop`](crate::runtime::run_worker_loop) the pool uses — over a
//! purpose-built storage row where its group's steps are `Owned` and every other
//! step is `Skip` ([`build_driver_storage`]) — with a
//! [`WorkerCore::driver`](crate::runtime::WorkerCore) (Park backoff, off-pool
//! stats attribution) and the [`DrainFirstScheduler`]. So a driver is literally a
//! 1-thread (or, for a group, still 1-thread over several `Owned` steps) instance
//! of the pool loop. Several detached steps sharing a
//! [`DetachedGroup::Shared`](crate::step::DetachedGroup) label are driven by ONE
//! thread that round-robins them; [`DetachedGroup::PerStep`] (the default) keeps
//! one thread per step.
//!
//! Because the group's steps and their phases are temporally disjoint (e.g. the
//! sort's phase-1 admit/sort/frame finish and leave the live set before the
//! phase-2 merge runs), one driver thread covers a whole phase's coordination
//! without oversubscription — exactly like main's single main thread.
//!
//! ## Lock-ordering acyclicity (non-negotiable, the deadlock proof)
//!
//! A driver thread must NEVER hold one queue's internal lock while parking on
//! another. `run_worker_loop` upholds that by construction:
//!
//!   1. Each `try_run_erased` is a single non-blocking call. The step body pops
//!      from its input transport (a lock-free `crossbeam ArrayQueue` — `try_pop`,
//!      no lock held across the call) and pushes to its outputs (`try_push`,
//!      likewise); a full output / empty input is reported back as `NoProgress` /
//!      `Contention`. It never *blocks* inside `try_run`.
//!   2. The only blocking a driver does is the loop's `WorkerCore::sleep_backoff`
//!      (`park_timeout` under the Park policy), which holds NO queue lock.
//!   3. The loop tries EVERY live step in a pass before it parks (round-robin,
//!      park only after a full no-progress pass). So when one grouped step is
//!      blocked, a sibling on the same driver still runs — a park-on-first-idle
//!      loop would wedge (see `driver_round_robins_all_live_before_parking`).
//!
//! So there is no cycle: a driver parks only *between* `try_run` calls, never
//! while holding a transport lock, so a two-sided step (consuming from the pool
//! AND producing to it, both bounded) cannot deadlock. The
//! `detached_two_sided_no_deadlock` test pins this with a wall-clock watchdog.

use std::any::Any;
use std::collections::HashMap;
use std::sync::Arc;

use crate::erased::{ErasedStep, ErasedStepCtx};
use crate::runtime::contexts::ChainContexts;
use crate::runtime::drain::StepDrainCounter;
use crate::runtime::driver::run_worker_loop;
use crate::runtime::scheduler::DrainFirstScheduler;
use crate::runtime::stats::PipelineStats;
use crate::runtime::storage::WorkerStepEntry;
use crate::runtime::worker_core::WorkerCore;
use crate::signal::PipelineSignal;
use crate::step::{Affinity, DetachedGroup, OutputsViewAny, StepKind, StepOutcome, StepProfile};
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
    fn detached_group(&self) -> DetachedGroup {
        // Never consulted: `extract_detached_steps` reads the real step's group
        // before swapping in this placeholder.
        DetachedGroup::PerStep
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

/// One dedicated driver thread's worth of extracted detached steps, in chain
/// (`StepIdx`) order. The caller spawns one OS thread per group and drives it
/// with [`run_detached_driver`].
pub struct DetachedDriverGroup {
    /// The steps this one driver thread runs, in chain order. Always non-empty
    /// and all [`StepKind::Detached`]. Private so those invariants — enforced by
    /// [`Self::new`] — cannot be bypassed by an external caller building the
    /// struct directly (which could otherwise trigger `steps[0]` panics in
    /// `primary_step` / `label`, or run non-detached work on an off-pool driver).
    steps: Vec<(StepIdx, Box<dyn ErasedStep>)>,
}

impl DetachedDriverGroup {
    /// Wrap a driver thread's extracted steps, enforcing the invariants every
    /// consumer relies on: the group is **non-empty** (`primary_step` / `label`
    /// index `steps[0]`) and **every step is [`StepKind::Detached`]** (a
    /// non-detached step must not run off-pool on a dedicated driver thread, and
    /// a `Parallel` step would hang its shared output — see
    /// [`build_driver_storage`]).
    ///
    /// # Panics
    ///
    /// Panics if `steps` is empty or contains a non-`Detached` step.
    #[must_use]
    fn new(steps: Vec<(StepIdx, Box<dyn ErasedStep>)>) -> Self {
        assert!(!steps.is_empty(), "detached driver group must be non-empty");
        assert!(
            steps.iter().all(|(_, step)| step.kind() == StepKind::Detached),
            "detached driver group may only contain Detached steps"
        );
        Self { steps }
    }

    /// The group's representative step (first in chain order). Used as the
    /// driver thread's off-pool stats key (`WorkerCore::driver`) and as a stable
    /// label. Non-empty by construction.
    #[must_use]
    pub fn primary_step(&self) -> StepIdx {
        self.steps[0].0
    }

    /// The name of the group's representative step (first in chain order), used
    /// to label a `PerStep` driver thread. Non-empty by construction.
    #[must_use]
    pub fn primary_name(&self) -> &'static str {
        self.steps[0].1.name()
    }

    /// The group label — derived from the steps (every step in the group reports
    /// the same [`DetachedGroup`]), used to name the driver thread. Non-empty by
    /// construction.
    #[must_use]
    pub fn label(&self) -> DetachedGroup {
        self.steps[0].1.detached_group()
    }

    /// Consume the group, yielding its steps for [`build_driver_storage`].
    #[must_use]
    fn into_steps(self) -> Vec<(StepIdx, Box<dyn ErasedStep>)> {
        self.steps
    }
}

/// Remove every [`StepKind::Detached`](crate::step::StepKind::Detached) step's
/// real instance from `steps`, replacing each in place with a
/// [`DetachedPlaceholder`] so the surviving slots keep their `step_idx`
/// positions (which `build_worker_storage` and `ChainContexts` index by).
///
/// Groups the extracted steps by their [`DetachedGroup`]: every
/// [`DetachedGroup::Shared`] label collects onto ONE group (one driver thread);
/// each [`DetachedGroup::PerStep`] step becomes its own singleton group (the
/// legacy one-thread-per-step behavior — the default, so non-sort chains are
/// unchanged). Within a group and across groups, order follows chain order
/// (first appearance). Called by `Pipeline::run` **before** `build_worker_storage`
/// consumes `steps`, while the (read-only) `ChainContexts` have already been
/// built from `&steps`.
#[must_use]
pub fn extract_detached_steps(steps: &mut [Box<dyn ErasedStep>]) -> Vec<DetachedDriverGroup> {
    // Accumulate each driver thread's steps as a raw vec, then wrap through
    // `DetachedDriverGroup::new` so the non-empty / all-Detached invariants are
    // enforced in one place rather than trusting each construction site.
    let mut group_steps: Vec<Vec<(StepIdx, Box<dyn ErasedStep>)>> = Vec::new();
    // Shared(label) -> index into `group_steps`, for O(1) append. PerStep steps
    // never share, so they are not indexed (each starts its own group).
    let mut shared_index: HashMap<&'static str, usize> = HashMap::new();
    for (idx, slot) in steps.iter_mut().enumerate() {
        if slot.kind() != StepKind::Detached {
            continue;
        }
        let group = slot.detached_group();
        let placeholder: Box<dyn ErasedStep> = Box::new(DetachedPlaceholder { name: slot.name() });
        let real = std::mem::replace(slot, placeholder);
        let entry = (StepIdx(idx), real);
        match group {
            DetachedGroup::PerStep => group_steps.push(vec![entry]),
            DetachedGroup::Shared(label) => {
                if let Some(&gi) = shared_index.get(label) {
                    group_steps[gi].push(entry);
                } else {
                    shared_index.insert(label, group_steps.len());
                    group_steps.push(vec![entry]);
                }
            }
        }
    }
    group_steps.into_iter().map(DetachedDriverGroup::new).collect()
}

/// Build a driver thread's storage row: a full-length `Vec<WorkerStepEntry>`
/// (length `n_total_steps`, indexed by global `step_idx` like every other row)
/// where each of `group_steps` is `Owned` and every other slot is `Skip`. The
/// driver thread runs [`run_worker_loop`] over this row exactly as a pool worker
/// runs over its own row.
///
/// # Panics
///
/// - if a group step's kind is `Parallel` — a 1-thread driver's
///   [`StepDrainCounter`] is init 1 (single finisher), which would never close a
///   `Parallel` step's shared output (that needs init N, all clones finishing),
///   hanging the downstream consumer;
/// - if two group steps map to the same `step_idx` (a double registration), or a
///   step index is out of range.
#[must_use]
pub fn build_driver_storage(
    group_steps: Vec<(StepIdx, Box<dyn ErasedStep>)>,
    n_total_steps: usize,
) -> Vec<WorkerStepEntry> {
    let mut row: Vec<WorkerStepEntry> = (0..n_total_steps).map(|_| WorkerStepEntry::Skip).collect();
    for (idx, step) in group_steps {
        assert_ne!(
            step.kind(),
            StepKind::Parallel,
            "driver group step `{}` is Parallel; a 1-thread driver's StepDrainCounter (init 1) \
             would never close its shared output — group only single-runner steps",
            step.name()
        );
        assert!(
            matches!(row[idx.0], WorkerStepEntry::Skip),
            "driver group step index {} registered twice (dual registration)",
            idx.0
        );
        row[idx.0] = WorkerStepEntry::Owned { step };
    }
    row
}

/// Drive one [`DetachedDriverGroup`] to completion on the calling (dedicated)
/// thread — the unified "1-thread pool". Builds the group's `Owned`/`Skip`
/// storage row ([`build_driver_storage`]) and runs the *same*
/// [`run_worker_loop`] the pool uses, with a [`WorkerCore::driver`] (Park
/// backoff, off-pool stats) and the [`DrainFirstScheduler`] (drain/seal
/// downstream before producing more — frees the sort's capacity-1 arena fastest).
///
/// `drain_counters` is the full per-step slice (init 1 for each detached step, so
/// the single finisher closes its output edges — the downstream consumer's
/// end-of-stream signal). On a cancel before `Finished`, outputs are NOT closed
/// (the run is tearing down; the recorded error/cancel is what propagates) —
/// `run_worker_loop`'s top-of-loop `is_done` break upholds this.
pub fn run_detached_driver(
    group: DetachedDriverGroup,
    contexts: &Arc<ChainContexts>,
    drain_counters: &[Arc<StepDrainCounter>],
    signal: &Arc<PipelineSignal>,
    stats: Option<&Arc<PipelineStats>>,
) {
    let primary = group.primary_step();
    let mut row = build_driver_storage(group.into_steps(), contexts.inputs.len());
    let mut worker = WorkerCore::driver(primary);
    run_worker_loop(
        &mut worker,
        &mut row,
        contexts,
        drain_counters,
        signal,
        stats,
        &DrainFirstScheduler,
    );
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

    /// Drive a single detached `step` at `step_idx` through the unified driver —
    /// a `PerStep` group of one — on the calling thread, with a full-length
    /// `drain_counters` slice (init 1 each). Mirrors what `builder.rs` step 4d
    /// does for a one-step group.
    fn drive_single(
        step: Box<dyn ErasedStep>,
        step_idx: StepIdx,
        contexts: &Arc<ChainContexts>,
        signal: &Arc<PipelineSignal>,
    ) {
        let drain_counters: Vec<Arc<StepDrainCounter>> =
            (0..contexts.inputs.len()).map(|_| StepDrainCounter::new(1)).collect();
        let group = DetachedDriverGroup::new(vec![(step_idx, step)]);
        run_detached_driver(group, contexts, &drain_counters, signal, None);
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

        let signal = PipelineSignal::new();
        drive_single(det, StepIdx(1), &contexts, &signal);

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

        let signal = PipelineSignal::new();
        drive_single(det, StepIdx(1), &contexts, &signal);

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

        // Drive the Detached step (as a one-step group) on this thread to
        // completion via the unified driver.
        drive_single(det, StepIdx(1), &contexts, &signal);

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

    /// A `Detached` step declaring a `Shared` group label. Used only to exercise
    /// `extract_detached_steps` grouping — never actually run.
    #[derive(Clone)]
    struct SharedDetached {
        label: &'static str,
    }
    impl Step for SharedDetached {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "SharedDetached",
                kind: StepKind::Detached,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn detached_group(&self) -> crate::step::DetachedGroup {
            crate::step::DetachedGroup::Shared(self.label)
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::Finished)
        }
    }

    /// A `Parallel` step — must never be placed on a 1-thread driver.
    #[derive(Clone)]
    struct ParallelStub;
    impl Step for ParallelStub {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ParallelStub",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    /// `extract_detached_steps` collects every `Shared(label)` onto one group,
    /// keeps each `PerStep` (default) step as its own singleton, preserves chain
    /// order within and across groups, and leaves non-detached steps in place.
    #[test]
    fn extract_groups_shared_together_and_perstep_alone() {
        use crate::step::DetachedGroup;
        // idx0 non-detached; idx1/idx3 Shared("coord"); idx2 PerStep; idx4 Shared("io").
        let mut steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(SrcStub { capacity: 4 })),
            Box::new(TypedStep::new(SharedDetached { label: "coord" })),
            Box::new(TypedStep::new(PassThroughDetached { held: None })),
            Box::new(TypedStep::new(SharedDetached { label: "coord" })),
            Box::new(TypedStep::new(SharedDetached { label: "io" })),
        ];
        let groups = extract_detached_steps(&mut steps);

        assert_eq!(groups.len(), 3, "coord{{1,3}}, perstep{{2}}, io{{4}}");
        // Order is first-appearance in chain order.
        assert_eq!(groups[0].label(), DetachedGroup::Shared("coord"));
        assert_eq!(
            groups[0].steps.iter().map(|(i, _)| i.0).collect::<Vec<_>>(),
            vec![1, 3],
            "shared group keeps both steps in chain order"
        );
        assert_eq!(groups[0].primary_step(), StepIdx(1));
        assert_eq!(groups[1].label(), DetachedGroup::PerStep);
        assert_eq!(groups[1].steps.iter().map(|(i, _)| i.0).collect::<Vec<_>>(), vec![2]);
        assert_eq!(groups[2].label(), DetachedGroup::Shared("io"));

        // Non-detached step survives; detached slots became placeholders that
        // still report `Detached` (so `build_worker_storage` Skips them on pool).
        assert_eq!(steps[0].kind(), StepKind::Exclusive);
        assert_eq!(steps[1].kind(), StepKind::Detached);
        assert_eq!(steps[2].kind(), StepKind::Detached);
    }

    #[test]
    #[should_panic(expected = "must be non-empty")]
    fn driver_group_rejects_empty() {
        // An empty group would panic later in `primary_step`/`label` (`steps[0]`);
        // the constructor rejects it up front.
        let _ = DetachedDriverGroup::new(vec![]);
    }

    #[test]
    #[should_panic(expected = "only contain Detached steps")]
    fn driver_group_rejects_non_detached_step() {
        // `SrcStub` is Exclusive, not Detached — running it off-pool on a driver
        // thread is a bug, so the constructor rejects the group.
        let step: Box<dyn ErasedStep> = Box::new(TypedStep::new(SrcStub { capacity: 1 }));
        let _ = DetachedDriverGroup::new(vec![(StepIdx(0), step)]);
    }

    /// `build_driver_storage` makes the group's steps `Owned` and every other
    /// slot `Skip`, at the correct global indices.
    #[test]
    fn build_driver_storage_owns_group_skips_rest() {
        let det: Box<dyn ErasedStep> = Box::new(TypedStep::new(PassThroughDetached { held: None }));
        let row = build_driver_storage(vec![(StepIdx(2), det)], 5);
        assert_eq!(row.len(), 5);
        assert!(matches!(row[2], WorkerStepEntry::Owned { .. }), "group step is Owned");
        for i in [0usize, 1, 3, 4] {
            assert!(matches!(row[i], WorkerStepEntry::Skip), "non-group slot {i} is Skip");
        }
    }

    /// G3: a `Parallel` step must never be grouped onto a 1-thread driver (its
    /// init-1 counter would never close the shared output).
    #[test]
    #[should_panic(expected = "is Parallel")]
    fn build_driver_storage_rejects_parallel_group_step() {
        let par: Box<dyn ErasedStep> = Box::new(TypedStep::new(ParallelStub));
        let _ = build_driver_storage(vec![(StepIdx(0), par)], 2);
    }

    /// G3: two group steps at the same index is a dual registration — rejected.
    #[test]
    #[should_panic(expected = "registered twice")]
    fn build_driver_storage_rejects_dual_registration() {
        let a: Box<dyn ErasedStep> = Box::new(TypedStep::new(PassThroughDetached { held: None }));
        let b: Box<dyn ErasedStep> = Box::new(TypedStep::new(PassThroughDetached { held: None }));
        let _ = build_driver_storage(vec![(StepIdx(1), a), (StepIdx(1), b)], 3);
    }
}
