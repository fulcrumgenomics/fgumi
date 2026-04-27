//! [`Stage`]: the pool-driven pipeline stage trait.
//!
//! # Two flavours of stage
//!
//! The pipeline engine recognises **two** stage traits; picking the right one
//! is the first design decision when adding a new stage:
//!
//! | Trait                            | Execution model                                           | When to pick it                                                                 |
//! |----------------------------------|-----------------------------------------------------------|---------------------------------------------------------------------------------|
//! | [`Stage`] (this file)            | Pool-driven, per-item `process()` call, 0-or-1 output    | Default. Almost every stage is a `Stage`.                                       |
//! | [`super::special_stage::SpecialStage`] | Owns its loop, runs on a dedicated thread, full queue I/O | Only when `Stage` semantics genuinely do not fit — see the decision tree below. |
//!
//! ## Decision tree: which trait should my new stage implement?
//!
//! Pick [`Stage`] unless **all three** of the `SpecialStage` criteria below
//! apply. When in doubt, start with [`Stage`]: it integrates with the pool's
//! scheduler, backpressure accounting, metrics, and progress bars for free.
//!
//! Pick [`super::special_stage::SpecialStage`] when your stage exhibits
//! **any** of:
//!
//! 1. **Emits 0 or many outputs per input** — `Stage::process` must call its
//!    output callback at most once. A stage that collapses multiple inputs
//!    into one output (e.g. position-batching, coalescing) or fans a single
//!    input out to many outputs cannot be expressed as a `Stage`.
//! 2. **Barrier semantics** — the stage must see *all* inputs before
//!    producing *any* output. External merge sort is the canonical example:
//!    phase 1 spills chunks, phase 2 performs a k-way merge across them.
//! 3. **Owns internal threads or subprocesses** — e.g. the aligner stage
//!    spawns a subprocess with stdin/stdout pipes and its own merge thread.
//!    The pool's one-item-at-a-time model cannot cohabit with that structure.
//!
//! Current `SpecialStage` implementors (for reference):
//! [`super::stages::sort::SortStage`] (barrier),
//! [`super::stages::align_and_merge::AlignAndMerge`] (subprocess +
//! threads), [`super::stages::position_batch::PositionBatchStage`] (1:N
//! emission with boundary detection), and
//! [`super::stages::coalesce::Coalesce`] (1:N re-batching with ordinal
//! reorder).
//!
//! Everything else in [`super::stages`] is a `Stage`.

/// Parallelism constraint declared by a [`Stage`] to the work-stealing pool.
///
/// `SpecialStage` implementors do not use this enum — they own their own
/// threading model and bypass the pool entirely.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Parallelism {
    /// Multiple workers may execute this stage concurrently. One instance
    /// of the stage struct is created per worker (via `Clone` or a factory).
    Parallel,
    /// Only one worker at a time may execute this stage. The pool shares a
    /// single `Mutex<Stage>` across workers so access is serialised.
    Sequential,
    /// Accumulates all input before producing output. Phase boundary.
    ///
    /// **Reserved, not currently used by any [`Stage`] implementor.** Barrier
    /// semantics in this codebase are expressed via `SpecialStage`
    /// (e.g. `SortStage`) because they also require owning the accumulate
    /// and merge loops.
    Barrier,
    /// External resource (e.g., subprocess) with dedicated I/O threads.
    /// Scheduler monitors queues but does not assign workers.
    ///
    /// **Reserved, not currently used by any [`Stage`] implementor.** Stages
    /// with dedicated-thread / subprocess semantics implement `SpecialStage`
    /// (e.g. `AlignAndMerge`) instead; the pool refuses to host them.
    Special,
}

/// Wraps an item with a sequence number and memory estimate.
///
/// Sequence numbers are assigned at the source and flow through the pipeline.
/// Intermediate queues are unordered; only the sink reorders by sequence number.
#[derive(Debug)]
pub struct SequencedItem<T> {
    pub seq: u64,
    pub item: T,
    pub memory_estimate: usize,
}

impl<T> SequencedItem<T> {
    pub fn new(seq: u64, item: T, memory_estimate: usize) -> Self {
        Self { seq, item, memory_estimate }
    }
}

/// A pool-driven pipeline stage: processes one `Input` into at most one
/// `Output` per call.
///
/// Stages are constructed with configuration and wired into a pipeline via
/// [`super::builder::PipelineBuilder`]'s `.stage(...)` (or `.stage_factory(...)`
/// for non-`Clone` stages). The work-stealing pool calls [`Stage::process`] from worker
/// threads based on the stage's declared [`Parallelism`]:
///
/// - `Parallel`: each worker holds its own `Clone`/factory-built instance
///   and they run concurrently.
/// - `Sequential`: one shared instance behind a `Mutex`; at most one worker
///   executes at a time.
///
/// For stages that cannot fit this model — 1:N output, barrier semantics,
/// internal threads / subprocesses — implement
/// [`super::special_stage::SpecialStage`] instead. See the module-level
/// docs for the full decision tree.
pub trait Stage: Send {
    /// Input item type. Must be `Send` so items can be moved across threads.
    type Input: Send + 'static;

    /// Output item type. Must be `Send` so items can be moved across threads.
    type Output: Send + 'static;

    /// Process one input item, optionally emitting an output.
    ///
    /// Call `output(value)` to emit an output item. The callback may be called
    /// 0 or 1 times per `process` call. Calling it more than once panics.
    ///
    /// Calling it 0 times means "suppress-on-empty" — the stage produced no
    /// useful output for this input (e.g., a coordination stage waiting for
    /// data from multiple upstream streams).
    ///
    /// # Errors
    ///
    /// May fail — errors propagate to the pipeline driver which cancels
    /// sibling workers.
    fn process(
        &mut self,
        input: Self::Input,
        output: &mut dyn FnMut(Self::Output),
    ) -> anyhow::Result<()>;

    /// Parallelism constraint for this stage.
    fn parallelism(&self) -> Parallelism;

    /// Estimated total bytes — stack AND owned heap — for a single output item.
    ///
    /// Used for queue memory accounting and backpressure. Implementations must
    /// return a realistic estimate; a stack-only value (e.g., `size_of::<T>()`)
    /// under-counts heap-owning types such as `Vec`, `Box`, `String`, or
    /// `HashMap`, causing the pipeline to under-apply backpressure and
    /// potentially exhaust memory.
    ///
    /// Guidelines:
    /// - For plain primitives / `Copy` types, `std::mem::size_of::<Self::Output>()`
    ///   is correct.
    /// - For `Vec<T>` / `String`, add `capacity() * size_of::<T>()` to the
    ///   stack size.
    /// - For composite types, sum the estimates of each field.
    ///
    /// There is deliberately no default: forcing implementors to think about
    /// memory prevents accidental under-counting.
    fn output_memory_estimate(&self, output: &Self::Output) -> usize;

    /// Human-readable stage name for logs and diagnostics.
    fn name(&self) -> &'static str;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequenced_item_construction() {
        let item = SequencedItem::new(42, "hello", 128);
        assert_eq!(item.seq, 42);
        assert_eq!(item.item, "hello");
        assert_eq!(item.memory_estimate, 128);
    }

    #[test]
    fn test_parallelism_equality() {
        assert_eq!(Parallelism::Parallel, Parallelism::Parallel);
        assert_ne!(Parallelism::Parallel, Parallelism::Sequential);
        assert_ne!(Parallelism::Barrier, Parallelism::Special);
    }

    struct DoubleStage;

    impl Stage for DoubleStage {
        type Input = u32;
        type Output = u32;

        fn process(
            &mut self,
            input: u32,
            output: &mut dyn FnMut(Self::Output),
        ) -> anyhow::Result<()> {
            output(input * 2);
            Ok(())
        }

        fn parallelism(&self) -> Parallelism {
            Parallelism::Parallel
        }

        fn output_memory_estimate(&self, _output: &Self::Output) -> usize {
            std::mem::size_of::<u32>()
        }

        fn name(&self) -> &'static str {
            "DoubleStage"
        }
    }

    #[test]
    fn test_stage_process() {
        let mut stage = DoubleStage;
        let mut captured = None;
        stage
            .process(21, &mut |v| {
                captured = Some(v);
            })
            .unwrap();
        assert_eq!(captured.expect("stage should emit"), 42);
        assert_eq!(stage.parallelism(), Parallelism::Parallel);
        assert_eq!(stage.name(), "DoubleStage");
    }

    #[test]
    fn test_explicit_memory_estimate() {
        let stage = DoubleStage;
        let output = 42u32;
        assert_eq!(stage.output_memory_estimate(&output), std::mem::size_of::<u32>());
    }
}
