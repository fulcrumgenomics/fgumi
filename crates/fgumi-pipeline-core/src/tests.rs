//! Phase 1 cross-module tests. Validates that the type vocabulary, queue
//! handles, and `TypedStep` adapter compose into a working chain when
//! wired manually.
//!
//! The runtime that wires these automatically (worker pool, worker loop)
//! lives in Phase 2. These tests confirm the Phase 1 pieces support the
//! unified report-`Finished` completion contract end-to-end.

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

/// Test `Step2`: sum of `(a, b)` → `u64`. Pops one item from each
/// branch per `try_run` call; emits `Some(a + b as u64)` if both
/// have items, `NoProgress` otherwise. Mirrors what zipper does at
/// the level of "pop from two queues in lockstep, combine".
struct SumPairStep;

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
        match (ctx.a.pop(), ctx.b.pop()) {
            (Some(a), Some(b)) => {
                let _ = ctx.outputs.push(u64::from(a) + u64::from(b));
                Ok(StepOutcome::Progress)
            }
            _ => Ok(StepOutcome::NoProgress),
        }
    }
}

#[test]
fn step2_typed_dispatch_pairs_both_branches() {
    // Build two single-output sub-chains feeding into a SumPair Step2.
    // Producer A and producer B each output u32 on branch 0. We manually
    // wire the TypedStep2<SumPairStep>'s TwoInputHandles<u32, u32> via
    // build_two_input_handles, then drive it through ErasedStep.

    let sum_step: Box<dyn ErasedStep> = Box::new(TypedStep2::new(SumPairStep));

    // Producer A's output set + handles.
    let (a_set, a_view) = <Single<u32> as super::outputs::StepOutputs>::build_queues(
        &[QueueSpec::CountBounded { capacity: 4 }],
        &[BranchOrdering::None],
    );
    let a_outputs: super::step::OutputHandles<Single<u32>> =
        super::step::OutputHandles::new(a_view);

    // Producer B's output set + handles.
    let (b_set, b_view) = <Single<u32> as super::outputs::StepOutputs>::build_queues(
        &[QueueSpec::CountBounded { capacity: 4 }],
        &[BranchOrdering::None],
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
    let (mut merge_outset, merge_view) = sum_step.build_output_set();
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
fn step2_build_two_input_handles_rejects_duplicate_producer() {
    // The TypedStep2 adapter's `build_two_input_handles` asserts that
    // p0_idx != p1_idx (each upstream is its own subchain). Pin this
    // so a future refactor that drops the check surfaces immediately.
    let sum_step: Box<dyn ErasedStep> = Box::new(TypedStep2::new(SumPairStep));
    let (a_set, _a_view) = <Single<u32> as super::outputs::StepOutputs>::build_queues(
        &[QueueSpec::CountBounded { capacity: 4 }],
        &[BranchOrdering::None],
    );
    let mut producer_sets = vec![a_set];
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        sum_step.build_two_input_handles(&mut producer_sets, 0, 0, 0, 0)
    }));
    assert!(result.is_err(), "expected panic on duplicate producer index");
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
            let _ = ctx.outputs.push(n);
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
            let _ = ctx.outputs.push(n);
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
                    let _ = ctx.outputs.push(u64::from(a) + u64::from(b));
                    Ok(StepOutcome::Progress)
                }
                _ if ctx.a.is_drained() && ctx.b.is_drained() => Ok(StepOutcome::Finished),
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

    let builder = Pipeline::builder();
    let chain_a = builder.chain(SourceA { remaining: 5 });
    let chain_b = builder.chain(SourceB { remaining: 5 });
    MultiChain2::<u32, u32>::from_chains(chain_a, chain_b)
        .join(PairSummer { pending_a: None, pending_b: None })
        .chain(CollectSink { received: Arc::clone(&received) })
        .into_sink_marker();

    let pipeline = builder.build().expect("pipeline build");
    // 4 threads: 2 Exclusive sources own 2 workers, the Exclusive sink
    // owns a third, leaving 1 free worker for the Serial PairSummer.
    pipeline.run(PipelineConfig { threads: 4, ..Default::default() }).expect("pipeline run");

    // Both sources emit 5 items in the same order (a: 5,4,3,2,1; b: 50,40,30,20,10).
    // Pair-summer outputs 5+50=55, 4+40=44, 3+30=33, 2+20=22, 1+10=11 — but
    // the pairing is interleave-based (whichever item arrives first on each
    // branch), so order across pairs isn't deterministic without an ordered
    // branch. Assert set equality on counts + sum.
    let collected = received.lock().clone();
    assert_eq!(collected.len(), 5, "expected 5 paired emits, got {collected:?}");
    // Pairs match A's i with B's i (both emit 5..=1 in step), so the multiset of
    // sums is exactly {11, 22, 33, 44, 55}. Assert the multiset, not just the
    // total: a sum check passes for any multiset summing to 165 (e.g. wrong
    // individual pairings), so it can't catch mis-paired values.
    let mut sorted = collected.clone();
    sorted.sort_unstable();
    assert_eq!(sorted, vec![11, 22, 33, 44, 55], "pair sums multiset must match");
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
/// to hold every item so `push` never hits backpressure (a Serial source hammered
/// by N workers must not silently drop on a full queue — unlike the cap-8
/// `DrainReproSource`, which is only safe at --threads 1). Capacity 256 covers
/// the test's 64 items.
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
    buffered: Vec<u32>,
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
        while let Some(&item) = self.buffered.first() {
            match ctx.outputs.push(item) {
                Ok(()) => {
                    self.buffered.remove(0);
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
            self.buffered.push(n);
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
    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

fn run_reports_finished_pipeline(threads: usize, n_items: u32) {
    let received: Arc<parking_lot::Mutex<Vec<u32>>> = Arc::new(parking_lot::Mutex::new(Vec::new()));
    let received_for_run = Arc::clone(&received);

    let (tx, rx) = std::sync::mpsc::channel::<()>();
    let handle = std::thread::spawn(move || {
        use crate::{Pipeline, PipelineConfig};
        let builder = Pipeline::builder();
        builder
            .chain(WideQueueSource { remaining: n_items })
            .chain(ReportsFinishedBuffer { buffered: Vec::new() })
            .chain(DrainReproSink { received: received_for_run })
            .into_sink_marker();
        let pipeline = builder.build().expect("pipeline build");
        pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("pipeline run");
        let _ = tx.send(());
    });

    assert!(
        rx.recv_timeout(std::time::Duration::from_secs(10)).is_ok(),
        "reports_finished mid-step pipeline DEADLOCKED or stalled at --threads {threads}"
    );
    handle.join().expect("worker thread panicked");
    // WideQueueSource emits exactly 0..n_items. Assert the multiset (not just the
    // count) so a dropped, duplicated, or value-corrupted record is caught, not
    // only an off-by-N total. Order isn't asserted: at >1 thread the Serial
    // step's re-dispatch interleaving is a valid scheduling detail.
    let mut got = received.lock().clone();
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
// `OrderedBytesSingle<_>` source outputs (matching the chain wrapper
// real BAM/FASTQ source subchains produce) so the typed-step migration
// validates the new converge method on its native input shape.
// ─────────────────────────────────────────────────────────────────────────────

#[test]
#[allow(clippy::too_many_lines)]
fn multi_chain2_ordered_pairs_two_byte_bounded_sources() {
    use super::builder::{MultiChain2Ordered, Pipeline, PipelineConfig};
    use super::item::{HeapSize, Ordered};
    use super::outputs::OrderedBytesSingle;
    use super::step::Step2;
    use super::step::StepCtx2;

    /// Minimal ordered, byte-bounded item type for the test.
    #[derive(Debug)]
    struct OrderedU32 {
        ordinal: u64,
        value: u32,
    }

    impl HeapSize for OrderedU32 {
        fn heap_size(&self) -> usize {
            0
        }
    }

    impl Ordered for OrderedU32 {
        fn ordinal(&self) -> u64 {
            self.ordinal
        }
    }

    /// Ordered source: emits items with monotonic `ordinal` carrying
    /// values N..=1 then `Finished`.
    #[derive(Clone)]
    struct OrderedSource {
        next_ordinal: u64,
        remaining: u32,
    }

    impl Step for OrderedSource {
        type Input = ();
        type Outputs = OrderedBytesSingle<OrderedU32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "OrderedSource",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 64 * 1024 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if self.remaining == 0 {
                return Ok(StepOutcome::Finished);
            }
            let value = self.remaining;
            self.remaining -= 1;
            let ordinal = self.next_ordinal;
            self.next_ordinal += 1;
            let _ = ctx.outputs.push(OrderedU32 { ordinal, value });
            Ok(StepOutcome::Progress)
        }
    }

    /// Step2 merger over two ordered inputs. Sums per-pair, emits
    /// `OrderedU64`. The cross-branch pairing is interleave-based, so
    /// the test asserts on the multiset of paired sums, not order.
    #[derive(Debug)]
    struct OrderedU64 {
        ordinal: u64,
        value: u64,
    }
    impl HeapSize for OrderedU64 {
        fn heap_size(&self) -> usize {
            0
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
                output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 64 * 1024 }],
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
            match (self.pending_a.as_ref(), self.pending_b.as_ref()) {
                (Some(_), Some(_)) => {
                    let a = self.pending_a.take().unwrap();
                    let b = self.pending_b.take().unwrap();
                    let out = OrderedU64 {
                        ordinal: self.next_out_ordinal,
                        value: u64::from(a.value) + u64::from(b.value),
                    };
                    self.next_out_ordinal += 1;
                    let _ = ctx.outputs.push(out);
                    Ok(StepOutcome::Progress)
                }
                _ if ctx.a.is_drained() && ctx.b.is_drained() => Ok(StepOutcome::Finished),
                _ => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// Records summed values.
    #[derive(Clone)]
    struct OrderedSink {
        received: Arc<parking_lot::Mutex<Vec<u64>>>,
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
                    self.received.lock().push(item.value);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    let received: Arc<parking_lot::Mutex<Vec<u64>>> = Arc::new(parking_lot::Mutex::new(Vec::new()));

    let builder = Pipeline::builder();
    let chain_a = builder.chain(OrderedSource { next_ordinal: 0, remaining: 5 });
    let chain_b = builder.chain(OrderedSource { next_ordinal: 0, remaining: 5 });
    MultiChain2Ordered::<OrderedU32, OrderedU32>::from_chains(chain_a, chain_b)
        .join(OrderedPairSummer { pending_a: None, pending_b: None, next_out_ordinal: 0 })
        .chain(OrderedSink { received: Arc::clone(&received) })
        .into_sink_marker();

    let pipeline = builder.build().expect("pipeline build");
    pipeline.run(PipelineConfig { threads: 4, ..Default::default() }).expect("pipeline run");

    let collected = received.lock().clone();
    assert_eq!(collected.len(), 5, "expected 5 paired emits, got {collected:?}");
    // Both sources emit values 5,4,3,2,1 in step → pair sums are 10,8,6,4,2.
    // Pairing is interleave-based (not order-deterministic), so assert the
    // multiset rather than the sum: a sum check passes for any multiset
    // totaling 30, missing mis-paired individual values.
    let mut sorted = collected.clone();
    sorted.sort_unstable();
    assert_eq!(sorted, vec![2, 4, 6, 8, 10], "pair sums multiset must match");
}
