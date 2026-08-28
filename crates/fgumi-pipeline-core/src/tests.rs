//! Phase 1 cross-module tests. Validates that the type vocabulary, queue
//! handles, and `TypedStep` adapter compose into a working chain when
//! wired manually.
//!
//! The runtime that wires these automatically (worker pool, worker loop)
//! lives in Phase 2. These tests confirm the Phase 1 pieces support the
//! unified report-`Finished` completion contract end-to-end.

use std::collections::VecDeque;
use std::io;
use std::sync::Arc;

use super::erased::{ErasedStep, ErasedStepCtx};
use super::outputs::Single;
use super::queues::QueueSpec;
use super::reorder::BranchOrdering;
use super::signal::PipelineSignal;
use super::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

// ─────────────────────────────────────────────────────────────────────────────
// F1' — Step2 (multi-input merge) tests.
//
// Validates that:
//   * `TypedStep2<S>` adapter dispatches correctly through `ErasedStep`.
//   * `TwoInputHandles<A, B>` round-trips per-branch handles via the
//     ChainContexts construction path.
//   * Both branches' drain signals propagate so `is_input_drained`
//     reports drained only when BOTH branches are drained.
// ─────────────────────────────────────────────────────────────────────────────

use super::erased::TypedStep2;
use super::step::{InputHandle as _, Step2, StepCtx2};

/// Test `Step2`: sum of `(a, b)` → `u64`. Emits `a + b as u64` once both
/// branches have yielded an item, `NoProgress` otherwise. Mirrors what zipper
/// does at the level of "pop from two queues in lockstep, combine".
///
/// Buffers a popped item whose sibling branch was empty rather than popping both
/// eagerly: `(ctx.a.pop(), ctx.b.pop())` discards any `Some` that lands in the
/// non-pair arm, which is silent data loss. This test's producers happen to be
/// pre-loaded with equal counts so the hazard never fires here — which is
/// exactly why it must not be modelled this way: a step author copying the
/// pattern into real code loses records, and a future change to the test's input
/// counts would turn the bug into a silently passing test. Same shape as
/// `PairSummer` below.
#[derive(Default)]
struct SumPairStep {
    pending_a: Option<u32>,
    pending_b: Option<u32>,
}

impl Step2 for SumPairStep {
    type InputA = u32;
    type InputB = u32;
    type Outputs = Single<u64>;
    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SumPair",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::CountBounded { capacity: 8 }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }
    fn try_run(&mut self, ctx: &mut StepCtx2<'_, Self>) -> io::Result<StepOutcome> {
        if self.pending_a.is_none() {
            self.pending_a = ctx.a.pop();
        }
        if self.pending_b.is_none() {
            self.pending_b = ctx.b.pop();
        }
        match (self.pending_a.take(), self.pending_b.take()) {
            (Some(a), Some(b)) => {
                ctx.outputs
                    .push(u64::from(a) + u64::from(b))
                    .expect("capacity-8 output must accept the pair");
                Ok(StepOutcome::Progress)
            }
            // Hold whichever branch did yield until its sibling catches up.
            (a, b) => {
                self.pending_a = a;
                self.pending_b = b;
                Ok(StepOutcome::NoProgress)
            }
        }
    }
}

#[test]
fn step2_typed_dispatch_pairs_both_branches() {
    // Build two single-output sub-chains feeding into a SumPair Step2.
    // Producer A and producer B each output u32 on branch 0. We manually
    // wire the TypedStep2<SumPairStep>'s TwoInputHandles<u32, u32> via
    // build_two_input_handles, then drive it through ErasedStep.

    let sum_step: Box<dyn ErasedStep> = Box::new(TypedStep2::new(SumPairStep::default()));

    // Producer A's output set + handles.
    let (a_set, a_view) = <Single<u32> as super::outputs::StepOutputs>::build_queues(
        &[QueueSpec::CountBounded { capacity: 4 }],
        &[BranchOrdering::None],
        crate::builder::InstrumentationLevel::Off,
    );
    let a_outputs: super::step::OutputHandles<Single<u32>> =
        super::step::OutputHandles::new(a_view);

    // Producer B's output set + handles.
    let (b_set, b_view) = <Single<u32> as super::outputs::StepOutputs>::build_queues(
        &[QueueSpec::CountBounded { capacity: 4 }],
        &[BranchOrdering::None],
        crate::builder::InstrumentationLevel::Off,
    );
    let b_outputs: super::step::OutputHandles<Single<u32>> =
        super::step::OutputHandles::new(b_view);

    // Build the merge step's TwoInputHandles<u32, u32> by handing the
    // ErasedStep adapter both producer sets.
    let mut producer_sets = vec![a_set, b_set];
    let merge_input_any = sum_step.build_two_input_handles(&mut producer_sets, 0, 0, 1, 0);
    let _ = producer_sets; // both branches taken by build_two_input_handles

    // Sanity: input_arity reports 2.
    assert_eq!(sum_step.input_arity(), 2);

    // Merge step's output set.
    let (mut merge_outset, merge_view) =
        sum_step.build_output_set(crate::builder::InstrumentationLevel::Off);
    let merge_outputs_any = sum_step.wrap_outputs_view(merge_view);
    let merge_consumer_input = merge_outset.take_typed_input::<u64>(0);

    // Push a couple of items on each producer.
    a_outputs.push(10).unwrap();
    a_outputs.push(20).unwrap();
    b_outputs.push(1).unwrap();
    b_outputs.push(2).unwrap();

    let signal = PipelineSignal::new();
    // Drive the merge step twice — once per pair.
    let mut sum_step = sum_step;
    for _ in 0..2 {
        let mut ctx = ErasedStepCtx {
            input: merge_input_any.as_ref(),
            outputs: merge_outputs_any.as_ref(),
            signal: &signal,
        };
        let outcome = sum_step.try_run_erased(&mut ctx).unwrap();
        assert_eq!(outcome, StepOutcome::Progress);
    }
    // No more pairs available — both branches non-drained but empty.
    {
        let mut ctx = ErasedStepCtx {
            input: merge_input_any.as_ref(),
            outputs: merge_outputs_any.as_ref(),
            signal: &signal,
        };
        let outcome = sum_step.try_run_erased(&mut ctx).unwrap();
        assert_eq!(outcome, StepOutcome::NoProgress);
    }

    // Verify the merge step pushed 11, 22 (10+1, 20+2) in order.
    assert_eq!(merge_consumer_input.pop(), Some(11));
    assert_eq!(merge_consumer_input.pop(), Some(22));
    assert_eq!(merge_consumer_input.pop(), None);
}

#[test]
fn step2_build_two_input_handles_rejects_same_producer_and_branch() {
    // When both `Step2` inputs share ONE producer step,
    // `build_two_input_handles` asserts the two branches differ — it does NOT
    // require distinct producer indices. Taking the same branch twice is the
    // illegal case; pin the panic so a refactor that drops the branch check
    // surfaces immediately. (The legal same-producer/distinct-branch case is
    // covered by the test below.)
    let sum_step: Box<dyn ErasedStep> = Box::new(TypedStep2::new(SumPairStep::default()));
    let (a_set, _a_view) = <Single<u32> as super::outputs::StepOutputs>::build_queues(
        &[QueueSpec::CountBounded { capacity: 4 }],
        &[BranchOrdering::None],
        crate::builder::InstrumentationLevel::Off,
    );
    let mut producer_sets = vec![a_set];
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        sum_step.build_two_input_handles(&mut producer_sets, 0, 0, 0, 0)
    }));
    assert!(result.is_err(), "expected panic on same producer AND same branch");
}

/// Like `SumPairStep` but with a **non-commutative** combiner, so a test can tell
/// input A from input B. `SumPairStep` cannot: `a + b` makes `10 + 1` and `1 + 10`
/// both 11, so a swapped branch-to-input wiring produces an identical result and
/// is invisible. `a * 1000 + b` distinguishes them (10001 vs 1010).
#[derive(Default)]
struct AsymmetricPairStep {
    pending_a: Option<u32>,
    pending_b: Option<u32>,
}

impl Step2 for AsymmetricPairStep {
    type InputA = u32;
    type InputB = u32;
    type Outputs = Single<u64>;
    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "AsymmetricPair",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::CountBounded { capacity: 8 }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }
    fn try_run(&mut self, ctx: &mut StepCtx2<'_, Self>) -> io::Result<StepOutcome> {
        if self.pending_a.is_none() {
            self.pending_a = ctx.a.pop();
        }
        if self.pending_b.is_none() {
            self.pending_b = ctx.b.pop();
        }
        match (self.pending_a, self.pending_b) {
            (Some(a), Some(b)) => {
                self.pending_a = None;
                self.pending_b = None;
                ctx.outputs
                    .push(u64::from(a) * 1000 + u64::from(b))
                    .expect("capacity-8 output must accept the pair");
                Ok(StepOutcome::Progress)
            }
            // Hold whichever branch did yield until its sibling catches up.
            _ => Ok(StepOutcome::NoProgress),
        }
    }
}

/// A fan-out producer feeding BOTH of a `Step2`'s inputs from its two distinct
/// branches is legal, and had no coverage — the duplicate-producer test above
/// passes only because both its branches are `0`, so a regression that took the
/// same branch twice would not be caught by it.
#[test]
fn step2_build_two_input_handles_accepts_one_producer_with_distinct_branches() {
    // Asymmetric combiner on purpose: `SumPairStep`'s `a + b` is commutative, so a
    // builder that wired branch 1 to input A and branch 0 to input B would produce
    // the same 11 and this test would still pass while asserting the mapping.
    let sum_step: Box<dyn ErasedStep> = Box::new(TypedStep2::new(AsymmetricPairStep::default()));

    // One producer, two branches: `(A, B)` is the two-branch output shape.
    let (fanout_set, fanout_view) = <(u32, u32) as super::outputs::StepOutputs>::build_queues(
        &[QueueSpec::CountBounded { capacity: 4 }, QueueSpec::CountBounded { capacity: 4 }],
        &[BranchOrdering::None, BranchOrdering::None],
        crate::builder::InstrumentationLevel::Off,
    );
    let fanout_outputs: super::step::OutputHandles<(u32, u32)> =
        super::step::OutputHandles::new(fanout_view);

    // Same producer index (0) for both inputs, distinct branches 0 and 1.
    let mut producer_sets = vec![fanout_set];
    let merge_input_any = sum_step.build_two_input_handles(&mut producer_sets, 0, 0, 0, 1);

    let (mut merge_outset, merge_view) =
        sum_step.build_output_set(crate::builder::InstrumentationLevel::Off);
    let merge_outputs_any = sum_step.wrap_outputs_view(merge_view);
    let merge_consumer_input = merge_outset.take_typed_input::<u64>(0);

    // Branch 0 feeds input A, branch 1 feeds input B.
    let view = fanout_outputs.view();
    view.a.push(10).unwrap();
    view.b.push(1).unwrap();

    let signal = PipelineSignal::new();
    let mut sum_step = sum_step;
    let mut ctx = ErasedStepCtx {
        input: merge_input_any.as_ref(),
        outputs: merge_outputs_any.as_ref(),
        signal: &signal,
    };
    assert_eq!(sum_step.try_run_erased(&mut ctx).unwrap(), StepOutcome::Progress);
    assert_eq!(
        merge_consumer_input.pop(),
        Some(10 * 1000 + 1),
        "branch 0 must feed input A and branch 1 input B — a swapped mapping yields 1010"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// F1' integration test — end-to-end Pipeline::run with TWO sources joined
// at a Step2, demonstrating that the entire chain-builder + runtime path
// (PipelineBuilder.chain ×2 → MultiChain2::from_chains → MultiChain2::join
// → Chain::chain → PipelineBuilder::build → Pipeline::run) works through
// the real worker loop.
// ─────────────────────────────────────────────────────────────────────────────

#[test]
#[allow(clippy::too_many_lines)]
fn step2_end_to_end_pipeline_pairs_two_sources_through_runtime() {
    use super::builder::{MultiChain2, Pipeline, PipelineConfig};

    /// Source A: emits 1..=5 then Finished.
    #[derive(Clone)]
    struct SourceA {
        remaining: u32,
    }
    impl Step for SourceA {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "SourceA",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 8 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if self.remaining == 0 {
                return Ok(StepOutcome::Finished);
            }
            let n = self.remaining;
            self.remaining -= 1;
            ctx.outputs.push(n).expect("capacity-8 output must accept all 5 items");
            Ok(StepOutcome::Progress)
        }
    }

    /// Source B: emits 10..=50 (step 10) then Finished — same count
    /// as Source A so every paired emit consumes one from each.
    #[derive(Clone)]
    struct SourceB {
        remaining: u32,
    }
    impl Step for SourceB {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "SourceB",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 8 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if self.remaining == 0 {
                return Ok(StepOutcome::Finished);
            }
            let n = self.remaining * 10;
            self.remaining -= 1;
            ctx.outputs.push(n).expect("capacity-8 output must accept all 5 items");
            Ok(StepOutcome::Progress)
        }
    }

    /// Step2 merger: sum each pair (a, b) → u64. Buffers a popped
    /// item from one branch if the sibling branch is empty —
    /// otherwise the popped item would be dropped when the sibling
    /// pop returns None, and we'd lose data.
    struct PairSummer {
        pending_a: Option<u32>,
        pending_b: Option<u32>,
    }
    impl Step2 for PairSummer {
        type InputA = u32;
        type InputB = u32;
        type Outputs = Single<u64>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "PairSummer",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 8 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx2<'_, Self>) -> io::Result<StepOutcome> {
            if self.pending_a.is_none() {
                self.pending_a = ctx.a.pop();
            }
            if self.pending_b.is_none() {
                self.pending_b = ctx.b.pop();
            }
            match (self.pending_a, self.pending_b) {
                (Some(a), Some(b)) => {
                    self.pending_a = None;
                    self.pending_b = None;
                    ctx.outputs
                        .push(u64::from(a) + u64::from(b))
                        .expect("capacity-8 output must accept all 5 pairs");
                    Ok(StepOutcome::Progress)
                }
                // Both inputs drained AND nothing buffered: genuinely done.
                _ if ctx.a.is_drained()
                    && ctx.b.is_drained()
                    && self.pending_a.is_none()
                    && self.pending_b.is_none() =>
                {
                    Ok(StepOutcome::Finished)
                }
                // Both inputs drained but one branch item is still buffered —
                // the two branches emitted different counts. Reporting
                // `Finished` here would silently drop that item, exactly the
                // data loss the `PairSummer` doc warns step authors about, so
                // fail loudly instead of modelling the bug.
                _ if ctx.a.is_drained() && ctx.b.is_drained() => {
                    panic!("PairSummer: unpaired item left buffered after both branches drained")
                }
                _ => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// Sink: records received items so the test can assert pairing.
    #[derive(Clone)]
    struct CollectSink {
        received: Arc<parking_lot::Mutex<Vec<u64>>>,
    }
    impl Step for CollectSink {
        type Input = u64;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "CollectSink",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(n) => {
                    self.received.lock().push(n);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    let received: Arc<parking_lot::Mutex<Vec<u64>>> = Arc::new(parking_lot::Mutex::new(Vec::new()));
    let received_for_run = Arc::clone(&received);

    run_with_deadlock_timeout("Step2 PairSummer join at --threads 4", move || {
        let builder = Pipeline::builder();
        let chain_a = builder.chain(SourceA { remaining: 5 });
        let chain_b = builder.chain(SourceB { remaining: 5 });
        MultiChain2::<u32, u32>::from_chains(chain_a, chain_b)
            .join(PairSummer { pending_a: None, pending_b: None })
            .chain(CollectSink { received: received_for_run })
            .into_sink_marker();

        let pipeline = builder.build().expect("pipeline build");
        // 4 threads: 2 Exclusive sources own 2 workers, the Exclusive sink
        // owns a third, leaving 1 free worker for the Serial PairSummer.
        pipeline.run(PipelineConfig { threads: 4, ..Default::default() }).expect("pipeline run");
    });

    // Both sources emit 5 items in step (a: 5,4,3,2,1; b: 50,40,30,20,10), so
    // pair sums are 55, 44, 33, 22, 11 — and that ORDER is guaranteed, not just
    // the multiset. Each branch queue is FIFO, `PairSummer` pops at most one
    // item per branch per tick into `pending_a`/`pending_b`, so A's k-th item
    // always pairs with B's k-th; and the step is `Serial`, so pairs enter the
    // output queue in emit order and the Exclusive sink pops them in that order.
    // Sorting here would discard an ordering guarantee the transport does
    // provide, letting a reordering regression in Serial dispatch pass.
    let collected = received.lock().clone();
    assert_eq!(
        collected,
        vec![55, 44, 33, 22, 11],
        "Serial dispatch must preserve FIFO pairing and output order"
    );
}

/// Sink: records every value it receives so tests can assert no records were
/// dropped, duplicated, or corrupted — not just that the count matched.
#[derive(Clone)]
struct DrainReproSink {
    received: Arc<parking_lot::Mutex<Vec<u32>>>,
}
impl Step for DrainReproSink {
    type Input = u32;
    type Outputs = ();
    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "DrainReproSink",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }
    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(n) => {
                self.received.lock().push(n);
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Unified model (issue #330): a mid-chain step (not just sources) reports
// `StepOutcome::Finished` once its input is drained and it holds no buffered
// output, folding its final flush into `try_run`'s flush-first path. A
// `NoProgress` tick mid-flush (full output queue) is just an idle yield — the
// step keeps flushing across re-dispatches until it reports `Finished`. At >1
// thread the shared `finished` latch must stop the other workers from
// re-running the finished Serial step.
// ─────────────────────────────────────────────────────────────────────────────

/// Source emitting `0..count` then `Finished`, with an output queue wide enough
/// to hold every item so `push` never hits backpressure: a Serial source
/// hammered by N workers must not silently drop on a full queue, which is why
/// the push below is asserted rather than discarded. Capacity 256 covers the
/// test's 64 items.
#[derive(Clone)]
struct WideQueueSource {
    remaining: u32,
}
impl Step for WideQueueSource {
    type Input = ();
    type Outputs = Single<u32>;
    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "WideQueueSource",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::CountBounded { capacity: 256 }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }
    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if self.remaining == 0 {
            return Ok(StepOutcome::Finished);
        }
        self.remaining -= 1;
        ctx.outputs.push(self.remaining).expect("wide source queue must never reject");
        Ok(StepOutcome::Progress)
    }
}

/// Mid (Serial): buffers every input item, then flushes the whole buffer in
/// `try_run`'s flush-first path through a capacity-1 output queue, and reports
/// `Finished` once input is drained and the buffer is empty — the unified
/// completion contract.
#[derive(Clone)]
struct ReportsFinishedBuffer {
    // A `VecDeque` so the flush-first drain pops the FIFO front in O(1). Real
    // steps process millions of records: a `Vec` with `remove(0)` would be O(n)
    // per flushed item (quadratic drain). Step authors copying this worked
    // example must keep the front-drain O(1) — never `Vec::remove(0)`.
    buffered: VecDeque<u32>,
}
impl Step for ReportsFinishedBuffer {
    type Input = u32;
    type Outputs = Single<u32>;
    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ReportsFinishedBuffer",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::CountBounded { capacity: 1 }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }
    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Flush-first: push held items until the queue rejects one.
        let mut pushed = false;
        while let Some(&item) = self.buffered.front() {
            match ctx.outputs.push(item) {
                Ok(()) => {
                    self.buffered.pop_front();
                    pushed = true;
                }
                Err(_unpushed) => break, // queue full — yield, retry next pass
            }
        }
        if pushed {
            return Ok(StepOutcome::Progress);
        }
        // 2. Consume input.
        if let Some(n) = ctx.input.pop() {
            self.buffered.push_back(n);
            return Ok(StepOutcome::Progress);
        }
        // 3. Completion: input drained AND nothing held — never push again.
        if ctx.input.is_drained() && self.buffered.is_empty() {
            return Ok(StepOutcome::Finished);
        }
        // Input drained but buffer still non-empty (queue was full this call):
        // NoProgress is just an idle tick — the round-robin moves on to the
        // sink, which drains the queue, and we flush more on the next pass,
        // eventually emptying the buffer and reporting Finished above. There is
        // no drain protocol to misfire here.
        Ok(StepOutcome::NoProgress)
    }
    // No `new_worker_copy` override: this step is `Serial`, so the framework
    // holds ONE shared instance behind a mutex and never clones it. An override
    // here would be dead code implying `buffered` could be split per worker.
}

/// Run `build_and_run` on its own thread and fail the test if it has not
/// returned within ten seconds.
///
/// Every end-to-end test in this module must go through here. These pipelines'
/// primary failure modes are deadlock and stall, and calling `Pipeline::run`
/// directly on the test thread turns either one into a hung harness with no
/// message — the run never returns, so no assertion is ever reached. `context`
/// names the pipeline in the timeout message.
///
/// A panic inside `build_and_run` is re-raised as itself, so an assertion
/// failure inside the run reports its own message rather than a spurious
/// deadlock.
fn run_with_deadlock_timeout(context: &str, build_and_run: impl FnOnce() + Send + 'static) {
    let (tx, rx) = std::sync::mpsc::channel::<()>();
    let handle = std::thread::spawn(move || {
        build_and_run();
        let _ = tx.send(());
    });
    match rx.recv_timeout(std::time::Duration::from_secs(10)) {
        Ok(()) => handle.join().expect("worker thread panicked"),
        // The sender dropped without sending, so `build_and_run` panicked —
        // typically a failed assertion inside the run. Re-raise that panic
        // instead of blaming a deadlock: reporting "DEADLOCKED or stalled" here
        // would invert the diagnosis for every assertion failure in a pipeline.
        Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => {
            let panic = handle.join().expect_err("a disconnect implies the closure panicked");
            std::panic::resume_unwind(panic);
        }
        // Nothing sent and the sender is still alive: a genuine stall.
        Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
            panic!("{context} DEADLOCKED or stalled")
        }
    }
}

fn run_reports_finished_pipeline(threads: usize, n_items: u32) {
    let received: Arc<parking_lot::Mutex<Vec<u32>>> = Arc::new(parking_lot::Mutex::new(Vec::new()));
    let received_for_run = Arc::clone(&received);

    run_with_deadlock_timeout(
        &format!("reports_finished mid-step pipeline at --threads {threads}"),
        move || {
            use crate::{Pipeline, PipelineConfig};
            let builder = Pipeline::builder();
            builder
                .chain(WideQueueSource { remaining: n_items })
                .chain(ReportsFinishedBuffer { buffered: VecDeque::new() })
                .chain(DrainReproSink { received: received_for_run })
                .into_sink_marker();
            let pipeline = builder.build().expect("pipeline build");
            pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("pipeline run");
        },
    );
    // WideQueueSource emits exactly 0..n_items. Assert the multiset (not just the
    // count) so a dropped, duplicated, or value-corrupted record is caught, not
    // only an off-by-N total.
    let mut got = received.lock().clone();
    if threads == 1 {
        // At one worker the sequence is fully determined, so sorting would throw
        // away a checkable guarantee: `WideQueueSource` emits `n_items-1..0`,
        // every transport is FIFO, and the single `ReportsFinishedBuffer` drains
        // its `VecDeque` front-first — so a reordering regression in the
        // flush-first path is only visible here. Assert identity, not multiset.
        let expected_order: Vec<u32> = (0..n_items).rev().collect();
        assert_eq!(
            got, expected_order,
            "single worker must deliver items in emit order (flush-first is FIFO)"
        );
        return;
    }
    // Above one worker, order is not asserted: the Serial step's re-dispatch
    // interleaving across workers is a valid scheduling detail.
    got.sort_unstable();
    let expected: Vec<u32> = (0..n_items).collect();
    assert_eq!(
        got, expected,
        "sink must receive every item 0..{n_items} exactly once (no drop/dup/corruption) \
         at --threads {threads}"
    );
}

#[test]
fn reports_finished_mid_step_single_thread_no_premature_drain() {
    // Lone worker: returning NoProgress mid-flush must not lose buffered output;
    // the round-robin must reach the sink so the full queue drains (no spin),
    // and the step eventually reports Finished. (The old BUG #5 drain-spin shape.)
    run_reports_finished_pipeline(1, 5);
}

#[test]
fn reports_finished_mid_step_multi_thread_shared_latch_stops_redispatch() {
    // 4 workers: the mid step is Serial + Affinity::None, so it sits in every
    // worker's live set. When one worker finishes it, the shared `finished`
    // latch must stop the others from re-running it (a non-idempotent flusher
    // re-entered would push into a drained queue → panic). A larger item count
    // widens the window for a concurrent re-dispatch.
    run_reports_finished_pipeline(4, 64);
}

// ─────────────────────────────────────────────────────────────────────────────
// `MultiChain2Ordered::{from_chains, join}` smoke test. Mirrors
// `step2_end_to_end_pipeline_pairs_two_sources_through_runtime` but uses
// `OrderedBytesSingle<_>` source outputs — the chain wrapper real BAM/FASTQ
// source subchains produce — so the converge method is exercised on its
// native input shape.
// ─────────────────────────────────────────────────────────────────────────────

#[test]
#[allow(clippy::too_many_lines)]
fn multi_chain2_ordered_pairs_two_byte_bounded_sources() {
    use super::builder::{MultiChain2Ordered, Pipeline, PipelineConfig};
    use super::item::{HeapSize, Ordered};
    use super::outputs::OrderedBytesSingle;
    use super::step::Step2;
    use super::step::StepCtx2;

    /// Bytes of heap payload each test item carries. Non-zero so the
    /// `ByteBounded` transports below actually account for these items:
    /// `ByteBoundedQueue` sums `HeapSize::heap_size()` only (it never counts
    /// `size_of::<T>()`), so a zero-heap item would make every push free, the
    /// byte cap unreachable, and this test blind to a regression that let a
    /// `ByteBounded` queue grow without bound — the memory-bound failure mode it
    /// exists to cover.
    const PAYLOAD_BYTES: usize = 1024;

    /// Byte budget per output edge. Deliberately smaller than
    /// `5 * PAYLOAD_BYTES` so the sources hit the cap partway through and must
    /// hold the rejected item and retry — exercising backpressure rather than
    /// just fitting everything in one go.
    const EDGE_LIMIT_BYTES: usize = 2 * PAYLOAD_BYTES;

    /// Minimal ordered, byte-bounded item type for the test.
    #[derive(Debug, Clone)]
    struct OrderedU32 {
        ordinal: u64,
        value: u32,
        payload: Vec<u8>,
    }

    impl OrderedU32 {
        fn new(ordinal: u64, value: u32) -> Self {
            Self { ordinal, value, payload: vec![0u8; PAYLOAD_BYTES] }
        }
    }

    impl HeapSize for OrderedU32 {
        fn heap_size(&self) -> usize {
            self.payload.capacity()
        }
    }

    impl Ordered for OrderedU32 {
        fn ordinal(&self) -> u64 {
            self.ordinal
        }
    }

    /// Ordered source: emits items with monotonic `ordinal` carrying
    /// values N..=1 then `Finished`.
    ///
    /// `branch_ordering` is `BranchOrdering::None` deliberately, despite the
    /// `OrderedBytesSingle` output shape and the `Ordered` items. The shape is
    /// chosen for its `HeapSize + Ordered` bounds and its byte-bounded transport,
    /// not to engage the reorder stage: `build_branch_ordered_bytes` maps
    /// `(ByteBounded, None)` to a **direct** branch, so these steps exercise
    /// byte-backpressure over a plain FIFO edge, which is exactly what the
    /// ordered-output assertion below reasons from (FIFO branch queues plus one
    /// pop per branch per tick). Declaring `ByOrdinal`/`ByItemOrdinal` here would
    /// change what the test covers, not strengthen it.
    ///
    /// Reorder-stage restoration is covered where it belongs — see
    /// `erased::tests::step2_serial_byitemordinal_output_is_reordered_not_collapsed`,
    /// `handles::handle_tests::ordered_branch_preserves_ordinal_across_retry`, and
    /// the `reorder` module's own tests.
    ///
    /// Holds an item the byte-bounded output
    /// rejected and retries it on a later tick — dropping it would punch a hole
    /// in the ordinal sequence while still reporting `Progress`, so the sink's
    /// count assertion would fail with no indication of the cause.
    #[derive(Clone)]
    struct OrderedSource {
        next_ordinal: u64,
        remaining: u32,
        held: Option<OrderedU32>,
    }

    impl Step for OrderedSource {
        type Input = ();
        type Outputs = OrderedBytesSingle<OrderedU32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "OrderedSource",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::ByteBounded {
                    limit_bytes: EDGE_LIMIT_BYTES as u64,
                }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            // Flush-first: retry the item a previous tick's full queue rejected.
            // A rejected push reports `NoProgress`, not `Contention`: the driver
            // treats the two identically, but `Contention` means "a Serial step's
            // mutex was held by another worker" and feeds `contention_count`, from
            // which the bottleneck verdict derives its SPIN ratio. Reporting it
            // for ordinary output backpressure would invent mutex contention that
            // never happened and can trip a bogus SPIN finding.
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
            let value = self.remaining;
            self.remaining -= 1;
            let ordinal = self.next_ordinal;
            self.next_ordinal += 1;
            match ctx.outputs.push(OrderedU32::new(ordinal, value)) {
                Ok(()) => Ok(StepOutcome::Progress),
                // Hold, never drop: the ordinal was already consumed above, so
                // discarding the item would leave a hole in the sequence.
                Err(unpushed) => {
                    self.held = Some(unpushed.into_item());
                    Ok(StepOutcome::NoProgress)
                }
            }
        }
    }

    /// Step2 merger over two ordered inputs. Sums per-pair, emits
    /// `OrderedU64`. Branch queues are FIFO and this step pops at most one item
    /// per branch per tick, so pair order is deterministic — the test asserts
    /// the exact output sequence, not just the multiset. Backpressure changes
    /// only *when* each item moves, never the order.
    #[derive(Debug, Clone)]
    struct OrderedU64 {
        ordinal: u64,
        value: u64,
        payload: Vec<u8>,
    }
    impl OrderedU64 {
        fn new(ordinal: u64, value: u64) -> Self {
            Self { ordinal, value, payload: vec![0u8; PAYLOAD_BYTES] }
        }
    }
    impl HeapSize for OrderedU64 {
        fn heap_size(&self) -> usize {
            self.payload.capacity()
        }
    }
    impl Ordered for OrderedU64 {
        fn ordinal(&self) -> u64 {
            self.ordinal
        }
    }

    struct OrderedPairSummer {
        pending_a: Option<OrderedU32>,
        pending_b: Option<OrderedU32>,
        /// A summed pair the byte-bounded output rejected, retried on a later
        /// tick. As in `OrderedSource`, dropping it would consume an out-ordinal
        /// and lose a record while still reporting `Progress`.
        held: Option<OrderedU64>,
        next_out_ordinal: u64,
    }
    impl Step2 for OrderedPairSummer {
        type InputA = OrderedU32;
        type InputB = OrderedU32;
        type Outputs = OrderedBytesSingle<OrderedU64>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "OrderedPairSummer",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::ByteBounded {
                    limit_bytes: EDGE_LIMIT_BYTES as u64,
                }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx2<'_, Self>) -> io::Result<StepOutcome> {
            // Flush-first: retry a pair the full output rejected earlier.
            // `NoProgress` rather than `Contention` on a rejected push, for the
            // reason given in `OrderedSource::try_run`.
            if let Some(out) = self.held.take() {
                return match ctx.outputs.push(out) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(unpushed) => {
                        self.held = Some(unpushed.into_item());
                        Ok(StepOutcome::NoProgress)
                    }
                };
            }
            if self.pending_a.is_none() {
                self.pending_a = ctx.a.pop();
            }
            if self.pending_b.is_none() {
                self.pending_b = ctx.b.pop();
            }
            match (self.pending_a.as_ref(), self.pending_b.as_ref()) {
                (Some(_), Some(_)) => {
                    let a = self.pending_a.take().unwrap();
                    let b = self.pending_b.take().unwrap();
                    let out = OrderedU64::new(
                        self.next_out_ordinal,
                        u64::from(a.value) + u64::from(b.value),
                    );
                    self.next_out_ordinal += 1;
                    match ctx.outputs.push(out) {
                        Ok(()) => Ok(StepOutcome::Progress),
                        Err(unpushed) => {
                            self.held = Some(unpushed.into_item());
                            Ok(StepOutcome::NoProgress)
                        }
                    }
                }
                // Same completion guard as `PairSummer`: both inputs drained and
                // neither branch item buffered. (`held` is always `None` here —
                // the flush-first block above returns when it is `Some`.)
                _ if ctx.a.is_drained()
                    && ctx.b.is_drained()
                    && self.pending_a.is_none()
                    && self.pending_b.is_none() =>
                {
                    Ok(StepOutcome::Finished)
                }
                _ if ctx.a.is_drained() && ctx.b.is_drained() => {
                    panic!(
                        "OrderedPairSummer: unpaired item left buffered after both branches drained"
                    )
                }
                _ => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// Records `(ordinal, value)`, not just `value`.
    ///
    /// The out-ordinal is the thing `OrderedPairSummer`'s hold-and-retry path
    /// exists to protect: it assigns `next_out_ordinal` *before* the push, so a
    /// rejected push must retry the SAME ordinal. Recording only `value` left that
    /// invariant unasserted — a regression that reassigned or skipped an ordinal on
    /// the retry path leaves the value sequence intact and the test still passes.
    #[derive(Clone)]
    struct OrderedSink {
        received: Arc<parking_lot::Mutex<Vec<(u64, u64)>>>,
    }
    impl Step for OrderedSink {
        type Input = OrderedU64;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "OrderedSink",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(item) => {
                    self.received.lock().push((item.ordinal, item.value));
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    // Pin the byte budget's teeth deterministically, independent of how the
    // scheduler happens to interleave the run below: an edge must admit at least
    // one item (or the pipeline would wedge) and must reject before all five fit
    // (or the cap would never bind and the test would prove nothing about byte
    // bounding). Both follow from PAYLOAD_BYTES vs EDGE_LIMIT_BYTES, so a future
    // edit to either constant that silently removes backpressure fails here.
    {
        use crate::queues::{ByteBoundedQueue, ItemQueue};
        let q = ByteBoundedQueue::<OrderedU32>::new(EDGE_LIMIT_BYTES as u64);
        let admitted = (0..5u32)
            .take_while(|i| q.try_push(OrderedU32::new(u64::from(*i), *i)).is_ok())
            .count();
        assert!(admitted >= 1, "a byte-bounded edge must always admit one item, else it wedges");
        assert!(
            admitted < 5,
            "the byte cap must bind before all 5 items fit — otherwise this test \
             exercises no byte bounding at all (admitted {admitted})"
        );
    }

    let received: Arc<parking_lot::Mutex<Vec<(u64, u64)>>> =
        Arc::new(parking_lot::Mutex::new(Vec::new()));
    let received_for_run = Arc::clone(&received);

    run_with_deadlock_timeout("MultiChain2Ordered byte-bounded join at --threads 4", move || {
        let builder = Pipeline::builder();
        let chain_a = builder.chain(OrderedSource { next_ordinal: 0, remaining: 5, held: None });
        let chain_b = builder.chain(OrderedSource { next_ordinal: 0, remaining: 5, held: None });
        MultiChain2Ordered::<OrderedU32, OrderedU32>::from_chains(chain_a, chain_b)
            .join(OrderedPairSummer {
                pending_a: None,
                pending_b: None,
                held: None,
                next_out_ordinal: 0,
            })
            .chain(OrderedSink { received: received_for_run })
            .into_sink_marker();

        let pipeline = builder.build().expect("pipeline build");
        pipeline.run(PipelineConfig { threads: 4, ..Default::default() }).expect("pipeline run");
    });

    let collected = received.lock().clone();
    // Both sources emit values 5,4,3,2,1 in step → pair sums 10,8,6,4,2, in that
    // order for the same reason as the `PairSummer` test above: FIFO branch
    // queues, one pop per branch per tick, and a `Serial` merge feeding an
    // `Exclusive` sink. Backpressure (the byte cap is smaller than the total
    // payload) changes only WHEN each item moves, never the order.
    assert_eq!(
        collected,
        vec![(0, 10), (1, 8), (2, 6), (3, 4), (4, 2)],
        "Serial dispatch must preserve FIFO pairing, output order, AND the out-ordinal \
         sequence under backpressure — the retry path must reuse an ordinal, not skip it"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Flattened public API.
//
// `lib.rs` re-exports the step-author surface at the crate root so an
// implementor writes `fgumi_pipeline_core::StepProfile`, not
// `::step::StepProfile`. A type reachable only through its module path is an
// asymmetry a step author trips over, and nothing else catches a dropped
// `pub use` — everything in-crate refers to these by module path anyway, and
// the crate has no in-tree consumer to break.
// ─────────────────────────────────────────────────────────────────────────────

/// Every type a step author needs must resolve at the crate root. This is a
/// compile-time assertion: naming them in type position fails to build if a
/// `pub use` is dropped or a new peer shape is added without one.
///
/// `DetachedGroup`, the ordered tuple shapes and `InstrumentationLevel` are the
/// ones this caught — `Step::detached_group` returns `DetachedGroup` while `Step`
/// itself was flattened, `OrderedBytesTuple2`/`3` sat beside an exported `Single`
/// and `OrderedBytesSingle`, and `InstrumentationLevel` appears in
/// `StepOutputs::build_queues`' signature while the trait itself was flattened.
///
/// `InstrumentationLevel` also turned up a whole class the earlier passes had
/// only sampled: it is the *parameter* of `build_queues`, and sweeping every
/// `pub` type named by a fully-public signature found the **return** side
/// unflattened too — `OutputQueueSet`, the `Tuple*View`s from
/// `OutputHandles::view()`, and `MultiChain2Ordered` from `Chain::into_multi()`.
/// Hence the sweep, not just the one symbol, is what this test now pins.
///
/// Deliberately NOT asserted: `queues::BoundedQueueHandle` and
/// `reorder::ReorderCapHandle`. They are runtime budget plumbing, named only via
/// `runtime::contexts::RegisteredQueue`, which is not flattened either — and
/// `ByteBoundedQueue` exposes inherent `limit_bytes` / `set_limit_bytes` /
/// `current_bytes`, so a direct user never needs the trait in scope. Promoting
/// them would make the root surface *less* coherent, not more.
///
/// Also deliberately NOT asserted: `handles::BranchOutputHandle`. No public
/// signature anywhere in the crate mentions it — a step author reaches its
/// methods through `OutputHandles` / the `Tuple*View`s and never names the type —
/// so there is no signature it has to be nameable in.
///
/// The rule this test encodes is therefore "nameable because a public signature
/// names it", not "`pub`, therefore flattened". `BranchInputHandle` is asserted
/// under exactly that rule: `OutputQueueSet::take_typed_input` returns it, so
/// flattening `OutputQueueSet` here brought it into reach of a root-only import
/// and it had to follow.
#[test]
fn crate_root_reexports_the_step_author_surface() {
    use crate as api;
    use std::marker::PhantomData;

    fn arity_of<O: api::StepOutputs>() -> usize {
        O::arity()
    }

    // Step / Step2 vocabulary.
    let _: fn() -> api::StepKind = || api::StepKind::Serial;
    let _: fn() -> api::StepOutcome = || api::StepOutcome::Progress;
    let _: fn() -> api::Affinity = || api::Affinity::None;
    let _: fn() -> api::DetachedGroup = || api::DetachedGroup::PerStep;
    let _: fn() -> api::QueueSpec = || api::QueueSpec::Unbounded;
    let _: fn() -> api::BranchOrdering = || api::BranchOrdering::None;
    // `StepOutputs::build_queues` names this in its public signature, so anyone
    // implementing that trait by hand has to be able to name it too.
    let _: fn() -> api::InstrumentationLevel = || api::InstrumentationLevel::Off;

    // Return types of public methods, which a step author has to name to store
    // one or to write a helper that takes it. `PhantomData` rather than a value:
    // naming the type in a position that checks its bounds is the whole
    // assertion, and none of these are constructible from outside the runtime.
    // `build_queues` returns `(OutputQueueSet, OutputsViewAny)` — the second half
    // was already flattened, the first was not.
    let _: PhantomData<api::OutputQueueSet> = PhantomData;
    // `OutputHandles::view()`, one per fan-out arity.
    let _: PhantomData<api::Tuple2View<'static, u32, u32>> = PhantomData;
    let _: PhantomData<api::Tuple3View<'static, u32, u32, u32>> = PhantomData;
    let _: PhantomData<api::Tuple4View<'static, u32, u32, u32, u32>> = PhantomData;
    // `Chain::into_multi()`, the ordered sibling of the exported `MultiChain2`.
    let _: PhantomData<api::MultiChain2Ordered<'static, api::Sequenced<u32>, api::Sequenced<u32>>> =
        PhantomData;
    // `OutputQueueSet::take_typed_input()` — flattening `OutputQueueSet` is what
    // put this one within reach of a root-only import.
    let _: PhantomData<api::BranchInputHandle<u32>> = PhantomData;
    // Naming it with its type is the whole assertion — a comparison against a
    // literal would be a constant expression, not a check.
    let _: usize = api::MAX_ARITY;

    // Output shapes: one arity per declared shape, so a new shape added without
    // a re-export shows up here.
    assert_eq!(arity_of::<api::Single<u32>>(), 1);
    assert_eq!(arity_of::<api::OrderedBytesSingle<api::Sequenced<u32>>>(), 1);
    assert_eq!(arity_of::<(u32, u32)>(), 2);
    assert_eq!(
        arity_of::<api::OrderedBytesTuple2<api::Sequenced<u32>, api::Sequenced<u32>>>(),
        2,
        "OrderedBytesTuple2 must be nameable at the crate root"
    );
    assert_eq!(
        arity_of::<
            api::OrderedBytesTuple3<api::Sequenced<u32>, api::Sequenced<u32>, api::Sequenced<u32>>,
        >(),
        3,
        "OrderedBytesTuple3 must be nameable at the crate root"
    );
}

/// The crate docs claim the dependency graph "stays light" and then enumerate
/// it. That list is a promise about the whole graph, so it has to be exhaustive
/// — a reader weighing this crate as a dependency reads the list, not the
/// manifest. Nothing else notices when the two drift: adding a dependency
/// compiles fine, and the prose keeps asserting the old, shorter graph.
///
/// `anyhow` is the one this caught. It backs `FinalizeHook::finalize`'s return
/// type and went undocumented, so the list understated the graph by one crate.
///
/// Dev-dependencies are deliberately out of scope: the claim is about what a
/// consumer links, and `proptest` / `rstest` / `trybuild` are not that.
/// Platform-gated tables are in scope, because a consumer on that platform does
/// link them — `fgumi-sort` already carries `[target.'cfg(unix)'.dependencies]`,
/// so reading only the plain `[dependencies]` table would fail *open* the day
/// this crate grows one.
#[test]
fn crate_docs_enumerate_every_runtime_dependency() {
    // `include_str!` resolves against this file's directory, so both paths are
    // the real files the claim is made in and about — not a copy that can drift.
    let manifest = include_str!("../Cargo.toml");
    let crate_docs = include_str!("lib.rs");

    // `[dependencies]` plus every `[target.'cfg(..)'.dependencies]`.
    // `[dev-dependencies]` and `[build-dependencies]` end in `-dependencies`,
    // so neither form matches them, at top level or under a `target` table.
    let is_runtime_table =
        |header: &str| header == "dependencies" || header.ends_with(".dependencies");

    let mut dependencies: Vec<&str> = Vec::new();
    let mut in_runtime_dependencies = false;
    for line in manifest.lines() {
        if let Some(header) = line.trim().strip_prefix('[').and_then(|h| h.strip_suffix(']')) {
            in_runtime_dependencies = is_runtime_table(header);
            // A crate may instead declare itself in the header, as
            // `[dependencies.log]` or `[target.'cfg(unix)'.dependencies.libc]`,
            // with only its own keys in the body. Recognizing the sub-table by
            // its parent is what keeps both spellings in the check — matching
            // the `dependencies.` prefix alone would take the first and let the
            // target-scoped one fall through undetected. The trim handles the
            // quoted-key spelling, `[dependencies."log"]`, which is also legal.
            if !in_runtime_dependencies
                && let Some((parent, name)) = header.rsplit_once('.')
                && is_runtime_table(parent)
            {
                dependencies.push(name.trim_matches('"'));
            }
            continue;
        }
        // A dependency key sits at column 0. Skipping indented and commented
        // lines is what keeps a multi-line inline table's `features = [..]`
        // continuation, or a commented-out `# tokio = ..`, from being read as a
        // crate name and failing this test under a name that is not a crate.
        if !in_runtime_dependencies || line.starts_with([' ', '\t', '#']) {
            continue;
        }
        if let Some((name, _)) = line.split_once('=') {
            // Cargo allows a dotted dependency key (`log.workspace = true`), whose
            // crate name is the segment before the dot. Take that, so the code-span
            // check below matches the crate name rather than `log.workspace` and
            // never reports a false undocumented dependency.
            let name = name.trim();
            let name = name.split_once('.').map_or(name, |(crate_name, _)| crate_name);
            if !name.is_empty() {
                dependencies.push(name);
            }
        }
    }
    // A renamed section or a reordered manifest would otherwise leave this test
    // asserting over an empty list and passing vacuously.
    assert!(
        !dependencies.is_empty(),
        "parsed no runtime dependencies from Cargo.toml; the section header or layout moved \
         and this test would silently stop checking anything"
    );

    // Odd-index pieces of a backtick split are the code spans. Matching spans
    // rather than raw substrings is what keeps `log` from being "documented" by
    // the word `logging`, and what lets `` `noodles::sam` `` document `noodles`.
    let code_spans: Vec<&str> = crate_docs
        .lines()
        .filter(|line| line.starts_with("//!"))
        .flat_map(|line| line.split('`').skip(1).step_by(2))
        .collect();
    let undocumented: Vec<&str> = dependencies
        .iter()
        .filter(|dependency| {
            !code_spans.iter().any(|span| {
                span == *dependency
                    || span.strip_prefix(*dependency).is_some_and(|rest| rest.starts_with("::"))
            })
        })
        .copied()
        .collect();

    assert!(
        undocumented.is_empty(),
        "crate docs in lib.rs enumerate the dependency graph but omit {undocumented:?}; \
         add them to the list or drop the claim"
    );
}
