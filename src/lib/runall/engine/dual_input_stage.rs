//! [`DualInputStage`]: a pipeline stage that consumes from two upstream
//! queues and emits to one downstream queue.
//!
//! This is a generalisation of [`super::stage::Stage`] for pipelines that
//! need fan-in: e.g. zipper-merging an unmapped stream with a mapped
//! stream, dedup with two ordering streams, or duplex consensus that
//! re-joins per-strand paths.
//!
//! The `Stage` trait is enough for almost every pipeline stage; reach for
//! `DualInputStage` only when the stage genuinely needs two upstream
//! inputs that are produced by different upstream chains. A single-input
//! `Stage` whose `Input` is an `Either<A, B>` — fed by a single upstream
//! merge queue — is usually a worse fit because both producers have to
//! agree on a tag-and-wrap protocol; the dual-input model lets each
//! upstream remain typed against its own output.
//!
//! # Parallelism
//!
//! `DualInputStage::parallelism` returns the same [`Parallelism`] enum as
//! `Stage`. The most common case is `Parallelism::Sequential` because
//! coordinated cross-input state (e.g. a name-based zipper that buffers
//! one stream while waiting for the other) cannot trivially shard across
//! workers. `Parallel` is allowed but only when the stage's two inputs
//! contain enough information per item that no cross-item state is
//! needed.
//!
//! # Drain semantics
//!
//! The driver calls [`DualInputStage::input_a_drained`] exactly once when
//! input A's upstream queue is closed and observed empty, and the
//! analogous [`DualInputStage::input_b_drained`] for input B. These hooks
//! let the stage flush any per-input buffer state. After both drain
//! callbacks have fired, the stage's run is complete and its output queue
//! is closed.

use std::marker::PhantomData;
use std::sync::Arc;

use anyhow::Result;

use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::driver::ErasedQueue;
use crate::runall::engine::erased_queue::{ErasingOutputQueue, UnerasingInputQueue};
use crate::runall::engine::sink::InputQueue;
use crate::runall::engine::source::OutputQueue;
use crate::runall::engine::stage::{Parallelism, SequencedItem};

/// A pool-driven pipeline stage with two input queues and one output
/// queue.
///
/// See module-level documentation for when to choose this over [`Stage`].
///
/// [`Stage`]: super::stage::Stage
pub trait DualInputStage: Send {
    /// Item type consumed from input A.
    type InputA: Send + 'static;
    /// Item type consumed from input B.
    type InputB: Send + 'static;
    /// Item type emitted on the single output queue.
    type Output: Send + 'static;

    /// Process one item popped from input A. May emit zero or one output
    /// items via `output(...)`.
    ///
    /// # Errors
    /// Returns an error if processing fails. The driver propagates the
    /// error to the caller and cancels sibling workers.
    fn process_a(
        &mut self,
        input: Self::InputA,
        output: &mut dyn FnMut(Self::Output),
    ) -> anyhow::Result<()>;

    /// Process one item popped from input B. Same emission contract as
    /// [`process_a`](Self::process_a).
    ///
    /// # Errors
    /// Returns an error if processing fails.
    fn process_b(
        &mut self,
        input: Self::InputB,
        output: &mut dyn FnMut(Self::Output),
    ) -> anyhow::Result<()>;

    /// Called exactly once when input A's upstream is observed drained
    /// (closed and empty). The stage may flush any per-A buffer state and
    /// emit additional outputs.
    ///
    /// Default: no-op.
    ///
    /// # Errors
    /// Returns an error if flushing fails.
    fn input_a_drained(
        &mut self,
        _output: &mut dyn FnMut(Self::Output),
    ) -> anyhow::Result<()> {
        Ok(())
    }

    /// Same as [`input_a_drained`](Self::input_a_drained) for input B.
    ///
    /// # Errors
    /// Returns an error if flushing fails.
    fn input_b_drained(
        &mut self,
        _output: &mut dyn FnMut(Self::Output),
    ) -> anyhow::Result<()> {
        Ok(())
    }

    /// Parallelism constraint, mirroring `Stage::parallelism`. Most
    /// `DualInputStage` implementations return `Sequential` because
    /// matching/joining across two streams requires cross-item state.
    fn parallelism(&self) -> Parallelism;

    /// Estimated heap+stack bytes for a single output item. Same
    /// contract as `Stage::output_memory_estimate`. Used for memory-
    /// bounded queue accounting.
    fn output_memory_estimate(&self, output: &Self::Output) -> usize;

    /// Human-readable stage name for logs and diagnostics.
    fn name(&self) -> &'static str;
}

/// Type-erased interface for [`DualInputStage`] so the engine can hold
/// heterogeneous dual-input stages in a `Box<dyn ErasedDualInputStage>`
/// without knowing their concrete `InputA`/`InputB`/`Output` types.
///
/// The input and output queues carry `Box<dyn Any + Send>` items; the
/// [`TypedDualInputStage`] adapter downcasts/boxes internally.
pub trait ErasedDualInputStage: Send {
    /// Execute the stage with type-erased queues.
    ///
    /// # Errors
    /// Returns an error if the stage encounters an unrecoverable failure.
    fn run_erased(
        self: Box<Self>,
        input_a: Arc<ErasedQueue>,
        input_b: Arc<ErasedQueue>,
        output: Arc<ErasedQueue>,
        cancel: CancelToken,
    ) -> Result<()>;

    /// Human-readable stage name for logs and diagnostics.
    fn name(&self) -> &'static str;
}

/// Adapter that wraps a concrete [`DualInputStage`] as an
/// [`ErasedDualInputStage`].
///
/// Mirrors `TypedSpecialStage` for the dual-input case: bridges between
/// type-erased `Box<dyn Any + Send>` queue items and the concrete
/// per-input typed items the stage works with.
pub struct TypedDualInputStage<S: DualInputStage> {
    inner: S,
    _phantom: PhantomData<fn(S::InputA, S::InputB) -> S::Output>,
}

impl<S: DualInputStage> TypedDualInputStage<S> {
    /// Wrap a concrete dual-input stage.
    #[must_use]
    pub fn new(inner: S) -> Self {
        Self { inner, _phantom: PhantomData }
    }
}

impl<S> ErasedDualInputStage for TypedDualInputStage<S>
where
    S: DualInputStage + 'static,
{
    fn run_erased(
        self: Box<Self>,
        input_a: Arc<ErasedQueue>,
        input_b: Arc<ErasedQueue>,
        output: Arc<ErasedQueue>,
        cancel: CancelToken,
    ) -> Result<()> {
        let typed_input_a: Box<dyn InputQueue<S::InputA>> =
            Box::new(UnerasingInputQueue::<S::InputA>::new(input_a));
        let typed_input_b: Box<dyn InputQueue<S::InputB>> =
            Box::new(UnerasingInputQueue::<S::InputB>::new(input_b));
        let typed_output: Box<dyn OutputQueue<S::Output>> =
            Box::new(ErasingOutputQueue::<S::Output>::new(output));
        run_dual_input_stage(self.inner, typed_input_a, typed_input_b, typed_output, cancel)
    }

    fn name(&self) -> &'static str {
        self.inner.name()
    }
}

/// Run a [`DualInputStage`] on the calling thread until both inputs are
/// drained or cancellation fires.
///
/// Worker loop invariants:
///
/// 1. **Round-robin for fairness**: alternate which input is checked
///    first across iterations so neither stream starves the other when
///    both have pending items.
/// 2. **Drain notifications fire exactly once each**, in the order each
///    input is observed drained-and-empty.
/// 3. **Output sequence numbers are stamped monotonically** by the
///    runner — the stage's `process_*` methods don't see ordinals.
/// 4. **Cancellation** is checked at the top of each iteration; the
///    stage can also error to bail.
///
/// This is a standalone helper. Integrating dual-input stages into the
/// main `run_multi_stage_inner` driver is a follow-on change; this
/// function lets stages be tested in isolation today.
///
/// # Errors
/// Returns the first stage error or input/output queue error encountered.
pub fn run_dual_input_stage<S>(
    mut stage: S,
    input_a: Box<dyn InputQueue<S::InputA>>,
    input_b: Box<dyn InputQueue<S::InputB>>,
    output: Box<dyn OutputQueue<S::Output>>,
    cancel: CancelToken,
) -> Result<()>
where
    S: DualInputStage + 'static,
{
    let mut backoff = crate::runall::engine::backoff::Backoff::new();
    let mut iter_count: u64 = 0;
    let mut a_drained_notified = false;
    let mut b_drained_notified = false;
    let mut out_seq: u64 = 0;

    // The output callback closure needs to share &mut output across calls
    // to process_a / process_b without holding it across the receive
    // operations. We wrap each invocation in its own borrow scope.
    macro_rules! call_with_output_capture {
        ($call:expr) => {{
            let mut emitted: Option<S::Output> = None;
            {
                let mut sink = |v: S::Output| {
                    debug_assert!(emitted.is_none(), "DualInputStage::process_* may emit at most one output");
                    emitted = Some(v);
                };
                $call(&mut sink)?;
            }
            if let Some(out_item) = emitted {
                let mem = stage.output_memory_estimate(&out_item);
                let item = SequencedItem::new(out_seq, out_item, mem);
                if output.push_until_cancelled(item, &cancel).is_err() {
                    // cancellation observed; exit cleanly
                    output.close();
                    return Ok(());
                }
                out_seq += 1;
            }
        }};
    }

    loop {
        if cancel.is_cancelled() {
            break;
        }

        let prefer_a = iter_count.is_multiple_of(2);
        iter_count = iter_count.wrapping_add(1);

        // Try preferred input first; fall back to the other.
        let a_pop = if prefer_a { input_a.pop() } else { None };
        if let Some(item) = a_pop {
            call_with_output_capture!(|sink: &mut dyn FnMut(S::Output)| stage.process_a(item.item, sink));
            backoff.reset();
            continue;
        }
        if let Some(item) = input_b.pop() {
            call_with_output_capture!(|sink: &mut dyn FnMut(S::Output)| stage.process_b(item.item, sink));
            backoff.reset();
            continue;
        }
        // If we didn't try A above (because we preferred B), try A now.
        if !prefer_a
            && let Some(item) = input_a.pop()
        {
            call_with_output_capture!(|sink: &mut dyn FnMut(S::Output)| stage.process_a(item.item, sink));
            backoff.reset();
            continue;
        }

        // No items in either input. Check drain status.
        let a_done = input_a.is_drained();
        let b_done = input_b.is_drained();

        if a_done && !a_drained_notified {
            call_with_output_capture!(|sink: &mut dyn FnMut(S::Output)| stage.input_a_drained(sink));
            a_drained_notified = true;
        }
        if b_done && !b_drained_notified {
            call_with_output_capture!(|sink: &mut dyn FnMut(S::Output)| stage.input_b_drained(sink));
            b_drained_notified = true;
        }

        if a_done && b_done {
            break;
        }

        backoff.snooze();
    }

    output.close();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A minimal `DualInputStage` for trait-shape tests: pairs items from
    /// A and B in arrival order, emitting `(a, b)` tuples whenever both
    /// streams have produced at least one item that hasn't been paired
    /// yet. On drain, any unmatched items are silently dropped.
    struct PairStage {
        pending_a: std::collections::VecDeque<u32>,
        pending_b: std::collections::VecDeque<u32>,
    }

    impl PairStage {
        fn new() -> Self {
            Self { pending_a: Default::default(), pending_b: Default::default() }
        }
    }

    impl DualInputStage for PairStage {
        type InputA = u32;
        type InputB = u32;
        type Output = (u32, u32);

        fn process_a(
            &mut self,
            input: Self::InputA,
            output: &mut dyn FnMut(Self::Output),
        ) -> anyhow::Result<()> {
            if let Some(b) = self.pending_b.pop_front() {
                output((input, b));
            } else {
                self.pending_a.push_back(input);
            }
            Ok(())
        }

        fn process_b(
            &mut self,
            input: Self::InputB,
            output: &mut dyn FnMut(Self::Output),
        ) -> anyhow::Result<()> {
            if let Some(a) = self.pending_a.pop_front() {
                output((a, input));
            } else {
                self.pending_b.push_back(input);
            }
            Ok(())
        }

        fn parallelism(&self) -> Parallelism {
            Parallelism::Sequential
        }
        fn output_memory_estimate(&self, _: &Self::Output) -> usize {
            std::mem::size_of::<(u32, u32)>()
        }
        fn name(&self) -> &'static str {
            "PairStage"
        }
    }

    /// Run a sequence of (input, value) operations through the stage and
    /// return the emitted output items. `input` is `'a'` or `'b'` to
    /// route to the corresponding `process_*` method. Avoids holding a
    /// closure-captured `&mut emitted` across calls.
    fn drive(stage: &mut PairStage, ops: &[(char, u32)]) -> Vec<(u32, u32)> {
        let mut emitted: Vec<(u32, u32)> = Vec::new();
        for &(which, v) in ops {
            let mut sink = |out: (u32, u32)| emitted.push(out);
            match which {
                'a' => stage.process_a(v, &mut sink).unwrap(),
                'b' => stage.process_b(v, &mut sink).unwrap(),
                _ => panic!("unknown input"),
            }
        }
        emitted
    }

    #[test]
    fn pair_stage_emits_when_both_streams_have_data() {
        // A1 buffers, B1 pairs with A1, B2 buffers, A2 pairs with B2.
        let mut stage = PairStage::new();
        let emitted = drive(&mut stage, &[('a', 1), ('b', 10), ('b', 20), ('a', 2)]);
        assert_eq!(emitted, vec![(1, 10), (2, 20)]);
    }

    #[test]
    fn run_dual_input_stage_pairs_via_typed_queues() {
        // End-to-end test of `run_dual_input_stage` with real
        // `StageQueue`s as the input/output. Pushes A and B items,
        // closes both inputs, and verifies the runner emits the
        // matched pairs in the expected order.
        use crate::runall::engine::memory::MemoryTracker;
        use crate::runall::engine::queue::StageQueue;

        let tracker = Arc::new(MemoryTracker::new(10_000_000));
        let qa: Arc<StageQueue<u32>> =
            Arc::new(StageQueue::new("test_a", 16, 10_000, tracker.clone()));
        let qb: Arc<StageQueue<u32>> =
            Arc::new(StageQueue::new("test_b", 16, 10_000, tracker.clone()));
        let qout: Arc<StageQueue<(u32, u32)>> =
            Arc::new(StageQueue::new("test_out", 16, 10_000, tracker));

        // Push A items, then B items. Close both.
        qa.push(SequencedItem::new(0, 1, 4)).unwrap();
        qa.push(SequencedItem::new(1, 2, 4)).unwrap();
        qa.push(SequencedItem::new(2, 3, 4)).unwrap();
        qb.push(SequencedItem::new(0, 10, 4)).unwrap();
        qb.push(SequencedItem::new(1, 20, 4)).unwrap();
        qb.push(SequencedItem::new(2, 30, 4)).unwrap();
        qa.close();
        qb.close();

        let stage = PairStage::new();
        let cancel = CancelToken::new();
        let input_a: Box<dyn InputQueue<u32>> = Box::new(qa.clone());
        let input_b: Box<dyn InputQueue<u32>> = Box::new(qb.clone());
        let output_q: Box<dyn OutputQueue<(u32, u32)>> = Box::new(qout.clone());

        run_dual_input_stage(stage, input_a, input_b, output_q, cancel)
            .expect("run_dual_input_stage must succeed");

        // Drain qout and check the pairs. With both A and B fully
        // available before the runner starts, the pair order depends on
        // the runner's per-iteration A-vs-B alternation.
        let mut emitted: Vec<(u32, u32)> = Vec::new();
        while let Some(item) = qout.pop() {
            emitted.push(item.item);
        }
        // Three full pairs expected; matched order can be either
        // (1,10),(2,20),(3,30) (strict round-robin) or some interleaved
        // variant depending on the runner's pop sequence — but exactly
        // three pairs must be emitted using each (a,b) at most once.
        assert_eq!(emitted.len(), 3);
        let (a_set, b_set): (std::collections::BTreeSet<u32>, std::collections::BTreeSet<u32>) =
            emitted.iter().copied().unzip();
        assert_eq!(a_set, [1u32, 2, 3].into_iter().collect());
        assert_eq!(b_set, [10u32, 20, 30].into_iter().collect());
    }

    #[test]
    fn pair_stage_drain_callbacks_are_optional() {
        // The default `input_a_drained` / `input_b_drained` are no-ops;
        // a stage with leftover buffered items doesn't error on drain.
        let mut stage = PairStage::new();
        let _ = drive(&mut stage, &[('a', 1), ('a', 2)]);
        let mut emitted: Vec<(u32, u32)> = Vec::new();
        {
            let mut sink = |out: (u32, u32)| emitted.push(out);
            stage.input_a_drained(&mut sink).unwrap();
            stage.input_b_drained(&mut sink).unwrap();
        }
        // pending_a still holds (1, 2) but that's expected for this
        // minimal stage; no output is emitted.
        assert!(emitted.is_empty());
    }
}
