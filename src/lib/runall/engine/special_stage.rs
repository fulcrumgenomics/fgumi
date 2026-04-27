//! [`SpecialStage`]: self-threaded pipeline stages — the escape hatch for
//! workloads that do not fit the pool-driven [`super::stage::Stage`] model.
//!
//! # When to implement `SpecialStage` instead of `Stage`
//!
//! `Stage` is the default. The pool owns the loop, calls `process()` per
//! item, handles backpressure accounting, records metrics, drives progress
//! bars, and supports both parallel (Clone/factory) and sequential
//! (Mutex-shared) execution. That covers the majority of stages in
//! [`super::stages`].
//!
//! `SpecialStage` exists because some workloads genuinely cannot be
//! expressed as "per-item, emit 0-or-1 outputs":
//!
//! 1. **1:N emission** — a stage that coalesces many inputs into one output
//!    (e.g. [`super::stages::coalesce::Coalesce`],
//!    [`super::stages::position_batch::PositionBatchStage`]) or fans one
//!    input into several outputs. `Stage::process` must call its output
//!    callback *at most once*, which makes 1:N impossible.
//! 2. **Barrier semantics** — the stage must consume *all* input before
//!    emitting anything (e.g. external merge sort in
//!    [`super::stages::sort::SortStage`], where phase 1 accumulates and
//!    spills and phase 2 k-way-merges the spilled chunks).
//! 3. **Owns internal threads or subprocesses** — e.g.
//!    [`super::stages::align_and_merge::AlignAndMerge`] spawns a subprocess
//!    plus three internal threads (stdin writer, stdout reader, merger)
//!    that cannot cohabit with the pool's per-item invocation.
//!
//! If *none* of those apply, implement [`super::stage::Stage`]: it hands
//! you metrics, progress bars, worker scaling, and backpressure accounting
//! for free.
//!
//! # Trade-offs of picking `SpecialStage`
//!
//! The engine gives up integration in exchange for flexibility:
//!
//! - **No automatic metrics** — `PipelineStats` and the per-stage metrics
//!   TSV (`run_multi_stage_with_metrics`) only cover pool stages.
//! - **No automatic progress bar** — indicatif bars are created per pool
//!   stage; special stages do not get one.
//! - **Manual cancellation handling** — `run()` must check
//!   `cancel.is_cancelled()` inside its loop and exit promptly.
//! - **Manual queue close** — `run()` must call `output.close()` on the
//!   happy path; the driver closes defensively too, but relying on that
//!   is an error-handling path, not the contract.
//! - **Manual memory estimates** — every pushed `SequencedItem` must carry
//!   a realistic byte estimate so backpressure still works.
//!
//! # Execution model
//!
//! When the pipeline builder sees a `.special_stage(...)` call, the driver
//! ([`super::driver`]) splits the chain into segments at that boundary:
//! one work-stealing pool per contiguous run of `Stage`s, and one dedicated
//! thread per `SpecialStage`. The special thread is handed type-erased
//! input and output queues (`Arc<ErasedQueue>`), which [`TypedSpecialStage`]
//! un-erases back to `Box<dyn InputQueue<Self::Input>>` /
//! `Box<dyn OutputQueue<Self::Output>>` before calling [`SpecialStage::run`].

use std::marker::PhantomData;
use std::sync::Arc;

use anyhow::Result;

use super::cancel::CancelToken;
use super::driver::ErasedQueue;
use super::erased_queue::{ErasingOutputQueue, UnerasingInputQueue};
use super::sink::InputQueue;
use super::source::OutputQueue;

/// A self-threaded pipeline stage that manages its own loop.
///
/// Use this trait only when the pool-driven [`super::stage::Stage`] model
/// does not fit — 1:N output, barrier semantics, or owned internal
/// threads/subprocesses. See the module-level docs for the full decision
/// tree.
///
/// The [`SpecialStage::run`] method receives ownership of `self` (via
/// `Box<Self>`) and must:
/// 1. Pop items from `input` until `input.is_drained()` (or cancellation).
/// 2. Push results to `output`, respecting backpressure — prefer
///    [`OutputQueue::push_until_cancelled`] over raw `push` so the stage
///    yields when downstream queues are full.
///
/// [`OutputQueue::push_until_cancelled`]: super::source::OutputQueue::push_until_cancelled
/// 3. Call `output.close()` on every exit path (happy and error). The
///    driver closes defensively too, but relying on that is an
///    error-recovery path, not the contract.
/// 4. Check `cancel.is_cancelled()` periodically inside every loop that
///    could otherwise block forever, and exit promptly if set.
///
/// Sequence numbers on input items must be preserved (or re-stamped
/// monotonically) on output items so the downstream sink can reorder.
///
/// Memory estimates on pushed [`SequencedItem`](super::stage::SequencedItem)s
/// must be realistic — typically `output_batch.data.len()` or similar — so
/// the engine's backpressure accounting still functions.
pub trait SpecialStage: Send {
    /// Input item type consumed from the upstream queue.
    type Input: Send + 'static;

    /// Output item type pushed to the downstream queue.
    type Output: Send + 'static;

    /// Execute the stage to completion.
    ///
    /// Blocks until all input is consumed and all output is pushed, or until
    /// `cancel` is set. Must call `output.close()` before returning on the
    /// happy path.
    ///
    /// # Errors
    ///
    /// Returns an error if the stage encounters an unrecoverable failure.
    /// The engine will cancel the pipeline on error.
    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<Self::Input>>,
        output: Box<dyn OutputQueue<Self::Output>>,
        cancel: CancelToken,
    ) -> Result<()>;

    /// Human-readable stage name for logs and diagnostics.
    fn name(&self) -> &'static str;
}

/// Type-erased interface for [`SpecialStage`] so the engine can hold
/// heterogeneous special stages in a `Box<dyn ErasedSpecialStage>` without
/// knowing their concrete `Input`/`Output` types.
///
/// The input and output queues carry `Box<dyn Any + Send>` items; the
/// [`TypedSpecialStage`] adapter downcasts/boxes internally.
pub trait ErasedSpecialStage: Send {
    /// Execute the stage with type-erased queues.
    ///
    /// # Errors
    ///
    /// Returns an error if the stage encounters an unrecoverable failure.
    fn run_erased(
        self: Box<Self>,
        input: Arc<ErasedQueue>,
        output: Arc<ErasedQueue>,
        cancel: CancelToken,
    ) -> Result<()>;

    /// Human-readable stage name for logs and diagnostics.
    fn name(&self) -> &'static str;
}

/// Adapter that wraps a concrete [`SpecialStage`] as an [`ErasedSpecialStage`].
///
/// Bridges between the type-erased `Box<dyn Any + Send>` queue items used
/// internally by the engine and the typed `S::Input`/`S::Output` items that
/// the concrete stage works with. Uses [`UnerasingInputQueue`] and
/// [`ErasingOutputQueue`] for zero-extra-copy conversion.
pub struct TypedSpecialStage<S: SpecialStage> {
    inner: S,
    _phantom: PhantomData<fn(S::Input) -> S::Output>,
}

impl<S: SpecialStage> TypedSpecialStage<S> {
    /// Wrap a concrete special stage.
    #[must_use]
    pub fn new(inner: S) -> Self {
        Self { inner, _phantom: PhantomData }
    }
}

impl<S> ErasedSpecialStage for TypedSpecialStage<S>
where
    S: SpecialStage + 'static,
{
    fn run_erased(
        self: Box<Self>,
        input: Arc<ErasedQueue>,
        output: Arc<ErasedQueue>,
        cancel: CancelToken,
    ) -> Result<()> {
        let typed_input: Box<dyn InputQueue<S::Input>> =
            Box::new(UnerasingInputQueue::<S::Input>::new(input));
        let typed_output: Box<dyn OutputQueue<S::Output>> =
            Box::new(ErasingOutputQueue::<S::Output>::new(output));
        Box::new(self.inner).run(typed_input, typed_output, cancel)
    }

    fn name(&self) -> &'static str {
        self.inner.name()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::queue::StageQueue;
    use crate::runall::engine::stage::SequencedItem;
    use std::sync::Arc;

    /// Passthrough mock: pops input, pushes unchanged to output.
    struct PassthroughSpecial;

    impl SpecialStage for PassthroughSpecial {
        type Input = u64;
        type Output = u64;

        fn run(
            self: Box<Self>,
            input: Box<dyn InputQueue<u64>>,
            output: Box<dyn OutputQueue<u64>>,
            cancel: CancelToken,
        ) -> Result<()> {
            loop {
                if cancel.is_cancelled() {
                    break;
                }
                if let Some(item) = input.pop() {
                    let _ =
                        output.push(SequencedItem::new(item.seq, item.item, item.memory_estimate));
                } else if input.is_drained() {
                    break;
                } else {
                    std::thread::yield_now();
                }
            }
            output.close();
            Ok(())
        }

        fn name(&self) -> &'static str {
            "PassthroughSpecial"
        }
    }

    #[test]
    fn test_special_stage_passthrough() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let input_q: Arc<StageQueue<u64>> =
            Arc::new(StageQueue::new("special_in", 16, 10_000, tracker.clone()));
        let output_q: Arc<StageQueue<u64>> =
            Arc::new(StageQueue::new("special_out", 16, 10_000, tracker));

        for i in 0u64..5 {
            input_q.push(SequencedItem::new(i, i * 10, 8)).unwrap();
        }
        input_q.close();

        let stage = Box::new(PassthroughSpecial);
        stage.run(Box::new(input_q), Box::new(output_q.clone()), CancelToken::new()).unwrap();

        let mut results: Vec<(u64, u64)> = Vec::new();
        while let Some(item) = output_q.pop() {
            results.push((item.seq, item.item));
        }
        results.sort_by_key(|(seq, _)| *seq);
        assert_eq!(results, vec![(0, 0), (1, 10), (2, 20), (3, 30), (4, 40)]);
        assert!(output_q.is_drained());
    }

    #[test]
    fn test_special_stage_respects_cancellation() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let input_q: Arc<StageQueue<u64>> =
            Arc::new(StageQueue::new("special_in", 16, 10_000, tracker.clone()));
        let output_q: Arc<StageQueue<u64>> =
            Arc::new(StageQueue::new("special_out", 16, 10_000, tracker));

        // Don't close input — only cancellation should exit.
        let cancel = CancelToken::new();
        cancel.cancel();

        let stage = Box::new(PassthroughSpecial);
        stage.run(Box::new(input_q), Box::new(output_q), cancel).unwrap();
    }
}
