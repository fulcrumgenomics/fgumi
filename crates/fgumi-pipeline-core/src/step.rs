//! Step trait, profile types, and context.
//!
//! Completion contract (full detail on [`StepOutcome::Finished`]):
//!
//! - Every step — source, mid, or sink — returns `StepOutcome::Finished` from
//!   `try_run` once its input edges are drained (empty + upstream closed) and
//!   it holds no buffered output. The framework then closes the step's output
//!   edges (`mark_outputs_drained`) and drops it from the worklist.
//! - `try_run` returns `Progress` when it pushed or held an item, `NoProgress`
//!   when there's nothing useful to do this call (input empty but not drained,
//!   no held items), and `Contention` when a Serial-step mutex is held by
//!   another worker; the scheduler reroutes.
//!
//! **Last-worker barrier for `Parallel` steps.** A `Parallel` step has N
//! per-worker `Clone`s sharing one output queue and a single drained input
//! edge. When the input drains, every clone returns `Finished`, but only the
//! LAST clone to finish (the one that takes the per-step `StepDrainCounter` to
//! 0) closes the shared output queue. Otherwise a clone could `mark_drained`
//! while a sibling is still pushing — a `try_push`-after-`mark_drained`
//! violation (debug-asserted in `ItemQueue::try_push`). The gate lives in the
//! driver's `dispatch_one_step`.

/// Concurrency profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepKind {
    /// Any number of workers may be inside `try_run` concurrently. Each
    /// worker holds its own `Clone` of the step.
    Parallel,
    /// At most one worker at a time. Framework holds a per-step mutex; any
    /// worker can acquire it; one instance shared.
    Serial,
    /// Exactly one worker (the owner) ever runs this step. Framework
    /// assigns owners at run start in chain declaration order. Other
    /// workers skip the step entirely.
    Exclusive,
    /// Runs on its own **dedicated OS thread**, spawned at run start
    /// alongside the deadlock-monitor / queue-rebalancer — never dispatched
    /// by the work-stealing pool. The dedicated thread drives the step with
    /// a *blocking* loop (park on empty input / full output via the
    /// `BlockingQueue` facet) instead of the pool's cooperative
    /// `try_run`-and-reschedule, so the step consumes no pool worker slot.
    ///
    /// This mirrors the legacy sort's "N + 2" threading: N pool workers do
    /// the parallel (compression-bound) work while the merge and the output
    /// writer each occupy a dedicated thread, keeping the pool saturated.
    ///
    /// **Opt-in.** No existing step declares `Detached`; it is wired only on
    /// the standalone sort chain's merge + writer. A single shared instance
    /// (like `Serial`/`Exclusive`, never `new_worker_copy`'d) runs on the
    /// dedicated thread; every pool worker gets a `Skip` entry for it.
    Detached,
}

/// Optional scheduling hint for `Serial` steps. Restricts which worker(s)
/// are eligible to attempt the step's mutex on each round-robin pass —
/// non-eligible workers `Skip` the step entirely, eliminating the
/// `try_lock` thrash that pure `Serial` exhibits at high thread counts.
///
/// The framework's per-step mutex is still in place (so a step author
/// can rely on single-thread-at-a-time semantics), but with a non-`None`
/// affinity only the hinted worker(s) ever acquire it.
///
/// Mirrors the legacy framework's per-thread `exclusive_step_owned`
/// mapping, where T0 is the reader, T(N-1) is the writer, and interior
/// exclusive steps fan out from both ends.
///
/// Default `None` keeps the existing pure-mutex-shared `Serial` behavior.
/// Ignored for `Parallel` and `Exclusive` kinds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Affinity {
    /// Any worker may attempt this Serial step. (Default; current behavior.)
    None,
    /// Restrict attempts to worker 0. Other workers `Skip` this step in
    /// dispatch — no `try_lock` thrash. Use for I/O sources where keeping
    /// reads on a single thread improves the kernel readahead pattern.
    Reader,
    /// Restrict attempts to worker `N - 1` (the last worker). Other
    /// workers `Skip`. Use for I/O sinks where the writer's `BufWriter`
    /// benefits from thread locality.
    Writer,
    /// Restrict attempts to a specific worker index. Out-of-range values
    /// (`>= n_threads`) trigger an assertion failure (panic) at run start
    /// via an always-on `assert!` in `build_worker_storage`.
    Worker(usize),
}

impl Affinity {
    /// Returns `true` if `worker_id` is eligible to attempt this Serial
    /// step under this affinity hint.
    #[must_use]
    pub fn eligible(self, worker_id: usize, n_workers: usize) -> bool {
        match self {
            Self::None => true,
            Self::Reader => worker_id == 0,
            Self::Writer => worker_id + 1 == n_workers,
            Self::Worker(idx) => worker_id == idx,
        }
    }
}

/// Static description of how a step is scheduled, plus per-output queue
/// configuration.
///
/// `output_queues[i]` selects the transport-layer queue type for branch `i`
/// (count-bounded, byte-bounded, or unbounded). `branch_ordering[i]` selects
/// whether the framework inserts a `ReorderStage` in front of the consumer's
/// input handle (so consumers see items in producer-emitted ordinal order).
///
/// Both vectors must have length equal to `S::Outputs::arity()`.
#[derive(Debug, Clone)]
pub struct StepProfile {
    pub name: &'static str,
    pub kind: StepKind,
    pub sticky: bool,
    pub output_queues: Vec<super::queues::QueueSpec>,
    pub branch_ordering: Vec<super::reorder::BranchOrdering>,
}

/// Outcome of a single `try_run` call.
///
/// **`Finished` contract.** Any step — source, mid, or sink — returns
/// `Finished` once it will never push to its output again: all its input edges
/// are drained (empty + upstream closed) and it holds no buffered output. The
/// framework then closes the step's output edges (`mark_outputs_drained`) and
/// drops it from the worklist. The step must already have flushed every item
/// (via `try_run`'s flush-first path) before returning `Finished` — returning
/// it with work still buffered loses data. For a `Parallel` step the per-step
/// `StepDrainCounter` gates the output close so only the last clone to finish
/// closes the shared queue (a sibling could still be pushing). Sources are the
/// degenerate case: their (unit) input edge is drained from birth, so they
/// return `Finished` on EOF.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepOutcome {
    /// Step did useful work (pushed an item or held one for later).
    Progress,
    /// Step had nothing to do this call (input empty, no held work).
    NoProgress,
    /// Step's Serial-step mutex was contended; scheduler reroutes.
    Contention,
    /// The step has drained all its input and holds no buffered output — it
    /// will never push again. The framework marks its output queues drained
    /// (counter-gated for `Parallel`) and drops it. See the `Finished`
    /// contract in the enum docs above.
    Finished,
}

/// Type-erased holder for the per-step outputs view, downcast to a typed
/// view (`SingleOutputsView`, `Tuple{2,3,4}OutputsView`, `UnitOutputsView`)
/// inside the `TypedStep<S>` adapter. Defined here so the `Step` trait
/// declaration can reference it; concrete views live in `handles.rs`.
///
/// Drain marking is per-branch via `BranchOutputHandle::mark_drained`,
/// exposed by typed `mark_all_drained()` methods on each `*OutputsView`.
/// `TypedStep<S>` (Phase 1 Task 15) downcasts `inner` to the right typed
/// view and calls `mark_all_drained` from the worker loop's drain path.
pub struct OutputsViewAny {
    pub(crate) inner: Box<dyn std::any::Any + Send + Sync>,
}

use std::io;

use super::item::HeapSize;
use super::outputs::StepOutputs;

/// Handle to this step's input queue.
///
/// Implemented by `OrderedQueueInputHandle<T>` in `handles.rs`. Steps see
/// only the trait surface; the concrete handle is constructed by the
/// framework at chain build time.
pub trait InputHandle<T: Send + HeapSize + 'static>: Send + Sync {
    /// Pop the next item, or `None` if the queue is empty right now.
    /// Non-blocking.
    fn pop(&self) -> Option<T>;

    /// Returns `true` if upstream has marked the queue drained AND the
    /// queue is currently empty. Once `is_drained()` returns true, no
    /// further items will arrive. Used by mid-steps and sinks to detect
    /// when to stop pulling.
    fn is_drained(&self) -> bool;
}

/// Handle to this step's output queues, shaped by `S::Outputs`.
///
/// Concrete shape is the per-arity view: `SingleOutputsView<T>` for
/// `Single<T>` outputs, `Tuple2OutputsView<A, B>` for `(A, B)`, etc.
/// Per-arity extension impls in `handles.rs` provide typed `push` methods
/// (e.g., `OutputHandles<Single<T>>::push(&self, item: T) -> Result<(), T>`).
pub struct OutputHandles<O: StepOutputs> {
    pub(crate) inner: OutputsViewAny,
    /// `PhantomData<fn() -> O>` (not `PhantomData<O>`) so `OutputHandles<O>`
    /// is always `Send + Sync` regardless of `O` — required so the runtime
    /// can box it as `Box<dyn Any + Send + Sync>` without forcing a `Sync`
    /// bound on item types.
    pub(crate) _phantom: std::marker::PhantomData<fn() -> O>,
}

impl<O: StepOutputs> OutputHandles<O> {
    /// Wrap a type-erased outputs view. Constructed by `TypedStep<S>`
    /// (Phase 1 Task 15) when assembling the step's `StepCtx`.
    #[allow(dead_code)] // wired up by Phase 1 Task 15 (TypedStep::wrap_outputs_view).
    pub(crate) fn new(inner: OutputsViewAny) -> Self {
        Self { inner, _phantom: std::marker::PhantomData }
    }
}

/// Step trait — every chain link implements this.
///
/// `Send + 'static`:
///   - `Send`: instances move across threads (workers run on dedicated threads).
///   - `'static`: no borrowed references in step state.
///
/// `Clone` is **not** a super-trait. Per-worker copies for `Parallel` steps
/// go through [`Step::new_worker_copy`] instead, which only `Parallel`
/// authors need to implement; `Serial`/`Exclusive` step authors inherit the
/// default panic and pay nothing.
///
/// # Per-worker copy patterns by step kind
///
/// **`Parallel` steps** must override [`Step::new_worker_copy`]. Each worker
/// thread calls it once during `build_worker_storage` to materialize its
/// private instance. The typical impl forwards to a regular `Clone` impl
/// (most parallel steps either `#[derive(Clone)]` or hand-roll a `Clone`
/// that resets per-worker scratch state). For example:
///
/// ```ignore
/// #[derive(Clone)]
/// pub struct ParseBamRecords;
///
/// impl Step for ParseBamRecords {
///     // ...
///     fn new_worker_copy(&self) -> Self { self.clone() }
/// }
/// ```
///
/// Closure-driven `Parallel` steps (`process(fn)`, `serialize(fn)`) wrap
/// their closure in `Arc<dyn Fn>` so the per-worker copy is one atomic
/// increment.
///
/// **`Serial` and `Exclusive` steps** do not override
/// [`Step::new_worker_copy`]. The framework holds a single instance (behind
/// a `Mutex` for `Serial`, pinned to one worker for `Exclusive`) and never
/// asks for additional copies. The default impl panics with a descriptive
/// message — if it ever fires, that's a framework bug.
pub trait Step: Send + Sized + 'static {
    type Input: Send + HeapSize + 'static;
    type Outputs: StepOutputs;

    /// Static description of how this step is scheduled.
    fn profile(&self) -> StepProfile;

    /// Optional scheduling hint for `Serial` kinds. Defaults to
    /// `Affinity::None` (any worker may attempt the step). Override for
    /// I/O sources/sinks that benefit from thread locality — see the
    /// [`Affinity`] enum docs.
    ///
    /// Ignored for `Parallel` and `Exclusive` kinds (those have their
    /// own per-worker dispatch model).
    fn affinity(&self) -> Affinity {
        Affinity::None
    }

    /// Step body. Pop from `ctx.input`, push to `ctx.outputs`. Returns
    /// `Progress` / `NoProgress` / `Contention` / `Finished`. Errors propagate via `Err`.
    ///
    /// # Errors
    ///
    /// Returns the underlying I/O error from any failed read, write, or
    /// codec operation inside the step body. The framework records the
    /// first error via `PipelineSignal::record_error` and broadcasts to
    /// other workers.
    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome>;

    /// Construct a fresh per-worker copy of this step. Only `Parallel`
    /// steps need to override this — the framework calls it during
    /// `build_worker_storage` to materialize one instance per worker.
    /// `Serial`/`Exclusive` steps inherit the default panic; the framework
    /// guarantees it is never invoked on them (one shared instance behind
    /// a `Mutex` for `Serial`, pinned to a single owner worker for
    /// `Exclusive`).
    ///
    /// # Panics
    ///
    /// Default impl panics with the step name + kind. Hitting it indicates
    /// a framework bug (the runtime tried to clone a non-`Parallel` step).
    #[must_use]
    fn new_worker_copy(&self) -> Self {
        let p = self.profile();
        panic!(
            "Step::new_worker_copy invoked on '{}' (kind = {:?}); \
             only Parallel steps need to override this. The framework \
             never clones Serial/Exclusive steps — file a bug.",
            p.name, p.kind
        );
    }
}

/// Context passed to `try_run`.
pub struct StepCtx<'a, S: Step> {
    pub input: &'a dyn InputHandle<S::Input>,
    pub outputs: &'a OutputHandles<S::Outputs>,
}

/// Two-input variant of [`Step`]. Used by merge steps (zipper,
/// AAM-aligner-output + original-record-stream) that need to pop
/// from two upstream queues independently.
///
/// Sits **alongside** [`Step`]; existing single-input steps are
/// unaffected. The framework's [`TypedStep2`] adapter
/// (`crate::erased::TypedStep2`) is the
/// dual of [`TypedStep`], implementing the same [`ErasedStep`]
/// contract so the runtime sees no difference between single- and
/// two-input steps at dispatch time.
///
/// The two input branches can be the same type (`InputA == InputB`
/// — what zipper uses, both `RawRecord`) or different types
/// (heterogeneous — also supported, no special case in the
/// framework).
///
/// See `docs/design/unified-pipeline-typed-step-migration.md` for
/// the broader migration plan this trait participates in.
pub trait Step2: Send + Sized + 'static {
    type InputA: Send + HeapSize + 'static;
    type InputB: Send + HeapSize + 'static;
    type Outputs: StepOutputs;

    fn profile(&self) -> StepProfile;

    /// Same semantics as [`Step::affinity`]. Defaults to
    /// `Affinity::None`.
    fn affinity(&self) -> Affinity {
        Affinity::None
    }

    /// Step body. Pop from `ctx.a` (input branch 0) and `ctx.b`
    /// (input branch 1), push to `ctx.outputs`.
    ///
    /// # Errors
    ///
    /// Same handling as [`Step::try_run`].
    fn try_run(&mut self, ctx: &mut StepCtx2<'_, Self>) -> io::Result<StepOutcome>;

    /// Same semantics as [`Step::new_worker_copy`]. Default panics;
    /// only `Parallel` Step2 impls need to override.
    ///
    /// # Panics
    ///
    /// Default impl panics with the step name + kind. Hitting it
    /// indicates a framework bug (the runtime tried to clone a
    /// non-`Parallel` step).
    #[must_use]
    fn new_worker_copy(&self) -> Self {
        let p = self.profile();
        panic!(
            "Step2::new_worker_copy invoked on '{}' (kind = {:?}); \
             only Parallel steps need to override this. The framework \
             never clones Serial/Exclusive steps — file a bug.",
            p.name, p.kind
        );
    }
}

/// Context passed to [`Step2::try_run`]. Holds independent input handles
/// for each branch.
pub struct StepCtx2<'a, S: Step2> {
    pub a: &'a dyn InputHandle<S::InputA>,
    pub b: &'a dyn InputHandle<S::InputB>,
    pub outputs: &'a OutputHandles<S::Outputs>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;

    #[test]
    fn step_profile_constructs() {
        let p = StepProfile {
            name: "Test",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::CountBounded { capacity: 64 }],
            branch_ordering: vec![BranchOrdering::None],
        };
        assert_eq!(p.name, "Test");
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky);
        assert_eq!(p.output_queues.len(), 1);
        assert!(matches!(p.output_queues[0], QueueSpec::CountBounded { capacity: 64 }));
        assert_eq!(p.branch_ordering, vec![BranchOrdering::None]);
    }

    #[test]
    fn step_kind_is_comparable() {
        assert_eq!(StepKind::Parallel, StepKind::Parallel);
        assert_ne!(StepKind::Parallel, StepKind::Serial);
        assert_ne!(StepKind::Serial, StepKind::Exclusive);
    }

    #[test]
    fn step_outcome_includes_finished() {
        let o = StepOutcome::Finished;
        assert_eq!(o, StepOutcome::Finished);
        assert_ne!(StepOutcome::Finished, StepOutcome::Progress);
        assert_ne!(StepOutcome::Finished, StepOutcome::NoProgress);
        assert_ne!(StepOutcome::Finished, StepOutcome::Contention);
    }
}

#[cfg(test)]
mod step_trait_compile_tests {
    use super::*;
    use crate::outputs::Single;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;

    /// Stub step that does nothing — exercises the trait declaration.
    #[derive(Clone)]
    struct NopStep;

    impl Step for NopStep {
        type Input = u32;
        type Outputs = Single<u32>;

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Nop",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 1 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }

        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    #[test]
    fn nop_step_advertises_profile() {
        let s = NopStep;
        let p = s.profile();
        assert_eq!(p.name, "Nop");
        assert_eq!(p.kind, StepKind::Parallel);
    }

    /// Stub multi-output step — exercises tuple Outputs.
    #[derive(Clone)]
    struct FanOutStep;

    impl Step for FanOutStep {
        type Input = u32;
        type Outputs = (u32, String);

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "FanOut",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![
                    QueueSpec::CountBounded { capacity: 1 },
                    QueueSpec::CountBounded { capacity: 1 },
                ],
                branch_ordering: vec![BranchOrdering::None, BranchOrdering::None],
            }
        }

        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    #[test]
    fn fan_out_step_advertises_profile() {
        let s = FanOutStep;
        assert_eq!(s.profile().name, "FanOut");
    }

    /// Default `new_worker_copy` panics on Serial steps. The framework
    /// guarantees it never invokes this on Serial/Exclusive in practice;
    /// this test pins the safety-net behavior in case a future refactor
    /// accidentally calls it. The panic message must include the step
    /// name so a real-world hit is debuggable.
    #[test]
    #[should_panic(expected = "SerialStub")]
    fn new_worker_copy_default_panics_for_serial() {
        struct SerialStub;
        impl Step for SerialStub {
            type Input = u32;
            type Outputs = Single<u32>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "SerialStub",
                    kind: StepKind::Serial,
                    sticky: false,
                    output_queues: vec![QueueSpec::CountBounded { capacity: 1 }],
                    branch_ordering: vec![BranchOrdering::None],
                }
            }
            fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                Ok(StepOutcome::NoProgress)
            }
        }
        let _ = SerialStub.new_worker_copy();
    }
}
