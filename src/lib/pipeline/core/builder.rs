//! `Pipeline`, `PipelineBuilder`, `Chain<O>`, `MultiChain2/3/4`, `BuildError`.
//!
//! The builder uses `RefCell` for the steps + graph so multiple `Chain`
//! handles (one per branch of a multi-output step) can coexist as `&self`.
//! Move semantics on `Chain` enforce single-consumer per branch at the
//! type-system level; the `ChainGraph` performs the runtime all-wired check.

use std::cell::RefCell;
use std::marker::PhantomData;
use std::sync::Arc;

use super::erased::{ErasedStep, TypedStep, TypedStep2};
use super::item::HeapSize;
use super::outputs::{Single, StepOutputs};
use super::runtime::stats::PipelineStats;
use super::signal::{CancelHandle, PipelineSignal};
use super::step::{Step, Step2};
use super::topology::{BranchIdx, ChainGraph, StepIdx};

/// Errors from `PipelineBuilder::build()`.
#[derive(Debug)]
pub enum BuildError {
    UnwiredOutput { step: &'static str, branch: &'static str },
    Empty,
}

impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnwiredOutput { step, branch } => {
                write!(f, "step {step:?} has unwired output branch {branch:?}")
            }
            Self::Empty => write!(f, "pipeline has no steps"),
        }
    }
}

impl std::error::Error for BuildError {}

#[derive(Debug, Clone)]
pub struct PipelineConfig {
    pub threads: usize,
    /// Optional shared stats collector. Construct via `Pipeline::stats()`,
    /// then pass into the config; the framework writes per-step counters
    /// during the run, and the caller reads them after `run` returns.
    /// `None` (the default) keeps the worker loop on its zero-cost path.
    pub stats: Option<Arc<PipelineStats>>,
    /// Deadlock-detection timeout in seconds. When > 0 and `stats` is set,
    /// `Pipeline::run` spawns a background monitor that polls the stats
    /// snapshot every `timeout / 4` seconds (clamped to ≥1s). If no step
    /// has progressed (or finished) for `timeout` seconds, the monitor
    /// logs the snapshot at warn level so the user has a starting point
    /// for debugging the stall. `0` (default) disables the monitor.
    ///
    /// Mirrors the legacy framework's `--deadlock-timeout` semantics
    /// (default 10s, 0 = disabled). Stats must be enabled for the
    /// monitor to have anything to read; the helpers in
    /// `commands/common.rs` auto-attach a stats handle when this is
    /// non-zero so callers don't have to pair the two flags manually.
    pub deadlock_timeout_secs: u64,
    /// Total byte budget for byte-bounded queues across the chain.
    /// `None` keeps each queue at the per-step limit set by its
    /// `QueueSpec::ByteBounded { limit_bytes }`. `Some(total)`
    /// enables the queue-memory rebalancer: an initial pass evenly
    /// distributes `total` across all byte-bounded queues, and a
    /// background thread periodically reads queue fullness and
    /// shifts budget toward consistently-full queues (producer-bound
    /// bottlenecks) at the expense of consistently-empty ones.
    ///
    /// Floor of 1 MiB per queue is enforced regardless of the
    /// rebalancer's decisions to prevent pathological starvation.
    pub queue_memory_total: Option<u64>,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        Self {
            threads: std::thread::available_parallelism().map_or(1, std::num::NonZero::get),
            stats: None,
            deadlock_timeout_secs: 0,
            queue_memory_total: None,
        }
    }
}

impl PipelineConfig {
    /// Builder-style helper to attach a stats collector.
    #[must_use]
    pub fn with_stats(mut self, stats: Arc<PipelineStats>) -> Self {
        self.stats = Some(stats);
        self
    }

    /// Builder-style helper to set the deadlock-detection timeout.
    /// `0` disables; the monitor only spawns when this is non-zero
    /// AND `stats` is set.
    #[must_use]
    pub fn with_deadlock_timeout(mut self, timeout_secs: u64) -> Self {
        self.deadlock_timeout_secs = timeout_secs;
        self
    }

    /// Builder-style helper to set the total queue-memory budget.
    /// `None` keeps per-step `ByteBounded` limits at their static
    /// values; `Some(total)` enables the rebalancer.
    #[must_use]
    pub fn with_queue_memory_total(mut self, total: Option<u64>) -> Self {
        self.queue_memory_total = total;
        self
    }
}

pub struct PipelineBuilder {
    inner: RefCell<BuilderInner>,
}

struct BuilderInner {
    steps: Vec<Box<dyn ErasedStep>>,
    graph: ChainGraph,
}

impl Default for PipelineBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl PipelineBuilder {
    #[must_use]
    pub fn new() -> Self {
        Self { inner: RefCell::new(BuilderInner { steps: Vec::new(), graph: ChainGraph::new() }) }
    }

    /// Add the first chain link. Requires `step: Step<Input = ()>` (i.e., a source).
    #[must_use = "pipeline branches must be wired to a sink"]
    pub fn chain<S>(&self, step: S) -> Chain<'_, S::Outputs>
    where
        S: Step<Input = ()>,
    {
        let mut inner = self.inner.borrow_mut();
        let producer = inner.graph.register_step(step.profile().name, S::Outputs::arity());
        inner.steps.push(Box::new(TypedStep::new(step)));

        Chain { builder: self, producer, branch: BranchIdx(0), _phantom: PhantomData }
    }

    /// Add a source step (first link) without requiring the typed `Chain`
    /// return value. Used by [`crate::pipeline::chains::ChainBuilder`]
    /// to accumulate steps across `add_source` / `add_<stage>` / `add_sink`
    /// method calls.
    ///
    /// Returns `(StepIdx, BranchIdx(0))` for the newly registered step.
    /// The caller tracks this tail and passes it to
    /// [`Self::append_step`] for all subsequent steps.
    ///
    /// This deliberately bypasses the typed `Chain<'_, S::Outputs>` API.
    /// Type correctness is NOT validated at [`Self::build`] time — `build()`
    /// only checks that every output branch is wired. A mis-typed step
    /// sequence (e.g. a step whose `Input=A` consumes an output of type `B`)
    /// builds successfully and panics at the first dispatch in
    /// `TypedStep::resolve_input` with "input handle downcast failed —
    /// chain topology invariant". That panic is loud and immediate but it
    /// is a runtime check, not a compile-time or build-time one. Callers
    /// are responsible for maintaining type correctness; the `pub(crate)`
    /// scope limits exposure to within `chains::builder`.
    pub(crate) fn append_source<S>(&self, step: S) -> (StepIdx, BranchIdx)
    where
        S: Step<Input = ()>,
    {
        let mut inner = self.inner.borrow_mut();
        let producer = inner.graph.register_step(step.profile().name, S::Outputs::arity());
        inner.steps.push(Box::new(TypedStep::new(step)));
        (producer, BranchIdx(0))
    }

    /// Append a step to the chain by wiring it to the current tail
    /// `(prev_producer, prev_branch)`. Used by
    /// [`crate::pipeline::chains::ChainBuilder`] for the same
    /// reason as [`Self::append_source`] — type-erased incremental chain
    /// assembly across method-boundary calls.
    ///
    /// Returns `(StepIdx, BranchIdx(0))` for the newly registered step.
    pub(crate) fn append_step<S: Step>(
        &self,
        step: S,
        prev: (StepIdx, BranchIdx),
    ) -> (StepIdx, BranchIdx) {
        let mut inner = self.inner.borrow_mut();
        let consumer = inner.graph.register_step(step.profile().name, S::Outputs::arity());
        inner.graph.wire(prev.0, prev.1, consumer);
        inner.steps.push(Box::new(TypedStep::new(step)));
        (consumer, BranchIdx(0))
    }

    /// Append a two-input [`Step2`] step, wiring `prev_a` into input slot 0
    /// and `prev_b` into input slot 1. The type-erased counterpart of
    /// [`MultiChain2Ordered::join`] — used by
    /// [`crate::pipeline::chains::ChainBuilder::add_zipper`] to wire
    /// the unmapped and mapped source chains into [`ZipperMergeStep`] across
    /// method-boundary calls.
    ///
    /// Returns `(StepIdx, BranchIdx(0))` for the newly registered step.
    ///
    /// [`ZipperMergeStep`]: crate::commands::zipper::merge_step::ZipperMergeStep
    pub(crate) fn append_step2<S: Step2>(
        &self,
        step: S,
        prev_a: (StepIdx, BranchIdx),
        prev_b: (StepIdx, BranchIdx),
    ) -> (StepIdx, BranchIdx) {
        let mut inner = self.inner.borrow_mut();
        let consumer =
            inner.graph.register_step_with_input_arity(step.profile().name, S::Outputs::arity(), 2);
        inner.graph.wire_to_slot(prev_a.0, prev_a.1, consumer, 0);
        inner.graph.wire_to_slot(prev_b.0, prev_b.1, consumer, 1);
        inner.steps.push(Box::new(TypedStep2::new(step)));
        (consumer, BranchIdx(0))
    }

    /// Finalize the chain. Returns `Err(UnwiredOutput)` if any output branch
    /// is dangling, `Err(Empty)` if the chain has zero steps.
    ///
    /// # Errors
    ///
    /// See `BuildError`.
    pub fn build(self) -> Result<Pipeline, BuildError> {
        let inner = self.inner.into_inner();
        if inner.steps.is_empty() {
            return Err(BuildError::Empty);
        }
        if let Some((producer, _branch, branch_name)) = inner.graph.first_unwired() {
            return Err(BuildError::UnwiredOutput {
                step: inner.graph.step_name(producer),
                branch: branch_name,
            });
        }
        Ok(Pipeline { steps: inner.steps, graph: inner.graph, signal: PipelineSignal::new() })
    }
}

/// In-progress chain handle. Each `.chain()` call consumes self and returns
/// a fresh `Chain` rooted at the new tail (Rust move semantics enforce
/// single-consumer at the type-system level).
#[must_use = "pipeline branches must be wired to a sink"]
pub struct Chain<'b, O> {
    builder: &'b PipelineBuilder,
    producer: StepIdx,
    branch: BranchIdx,
    _phantom: PhantomData<fn() -> O>,
}

impl<'b, T: Send + HeapSize + 'static> Chain<'b, Single<T>> {
    /// Extend the chain with a step accepting `T`.
    #[must_use = "pipeline branches must be wired to a sink"]
    pub fn chain<S>(self, step: S) -> Chain<'b, S::Outputs>
    where
        S: Step<Input = T>,
    {
        let mut inner = self.builder.inner.borrow_mut();
        let consumer = inner.graph.register_step(step.profile().name, S::Outputs::arity());
        inner.graph.wire(self.producer, self.branch, consumer);
        inner.steps.push(Box::new(TypedStep::new(step)));

        Chain {
            builder: self.builder,
            producer: consumer,
            branch: BranchIdx(0),
            _phantom: PhantomData,
        }
    }
}

impl<'b, T: Send + super::item::HeapSize + super::item::Ordered + 'static>
    Chain<'b, super::outputs::OrderedBytesSingle<T>>
{
    /// Extend the chain with a step accepting `T`. Mirror of
    /// `Chain<Single<T>>::chain` for the heap-aware ordered output shape
    /// used by Phase 3 BAM steps.
    #[must_use = "pipeline branches must be wired to a sink"]
    pub fn chain<S>(self, step: S) -> Chain<'b, S::Outputs>
    where
        S: Step<Input = T>,
    {
        let mut inner = self.builder.inner.borrow_mut();
        let consumer = inner.graph.register_step(step.profile().name, S::Outputs::arity());
        inner.graph.wire(self.producer, self.branch, consumer);
        inner.steps.push(Box::new(TypedStep::new(step)));

        Chain {
            builder: self.builder,
            producer: consumer,
            branch: BranchIdx(0),
            _phantom: PhantomData,
        }
    }
}

impl Chain<'_, ()> {
    /// Convenience to drop a sink-tail chain handle without `let _ = ...`.
    /// Sinks have `Outputs = ()` (zero branches) so the all-wired check
    /// never flags them; this method exists only to suppress the
    /// `#[must_use]` warning at the call site when a chain naturally ends
    /// at a sink.
    pub fn into_sink_marker(self) {
        let _ = (self.builder, self.producer, self.branch);
    }
}

impl<'b, A: Send + HeapSize + 'static, B: Send + HeapSize + 'static> Chain<'b, (A, B)> {
    /// Convert a 2-output chain into per-branch sub-chains.
    #[must_use = "all chain branches must be wired to a sink"]
    pub fn into_multi(self) -> MultiChain2<'b, A, B> {
        MultiChain2 {
            b0: Chain {
                builder: self.builder,
                producer: self.producer,
                branch: BranchIdx(0),
                _phantom: PhantomData,
            },
            b1: Chain {
                builder: self.builder,
                producer: self.producer,
                branch: BranchIdx(1),
                _phantom: PhantomData,
            },
            _phantom: PhantomData,
        }
    }
}

impl<'b, A, B> Chain<'b, super::outputs::OrderedBytesTuple2<A, B>>
where
    A: Send + super::item::HeapSize + super::item::Ordered + 'static,
    B: Send + super::item::HeapSize + super::item::Ordered + 'static,
{
    /// Convert a 2-output ordered + byte-bounded chain into per-branch
    /// sub-chains. Each branch is exposed as
    /// `Chain<OrderedBytesSingle<X>>` so downstream chained steps see
    /// the byte-aware ordered shape (matching what
    /// `Chain<OrderedBytesSingle<T>>::chain` accepts).
    #[must_use = "all chain branches must be wired to a sink"]
    pub fn into_multi(self) -> MultiChain2Ordered<'b, A, B> {
        MultiChain2Ordered {
            b0: Chain {
                builder: self.builder,
                producer: self.producer,
                branch: BranchIdx(0),
                _phantom: PhantomData,
            },
            b1: Chain {
                builder: self.builder,
                producer: self.producer,
                branch: BranchIdx(1),
                _phantom: PhantomData,
            },
            _phantom: PhantomData,
        }
    }
}

impl<'b, A, B, C> Chain<'b, (A, B, C)>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
{
    #[must_use = "all chain branches must be wired to a sink"]
    pub fn into_multi(self) -> MultiChain3<'b, A, B, C> {
        MultiChain3 {
            b0: Chain {
                builder: self.builder,
                producer: self.producer,
                branch: BranchIdx(0),
                _phantom: PhantomData,
            },
            b1: Chain {
                builder: self.builder,
                producer: self.producer,
                branch: BranchIdx(1),
                _phantom: PhantomData,
            },
            b2: Chain {
                builder: self.builder,
                producer: self.producer,
                branch: BranchIdx(2),
                _phantom: PhantomData,
            },
            _phantom: PhantomData,
        }
    }
}

impl<'b, A, B, C, D> Chain<'b, (A, B, C, D)>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
    D: Send + HeapSize + 'static,
{
    #[must_use = "all chain branches must be wired to a sink"]
    pub fn into_multi(self) -> MultiChain4<'b, A, B, C, D> {
        MultiChain4 {
            b0: Chain {
                builder: self.builder,
                producer: self.producer,
                branch: BranchIdx(0),
                _phantom: PhantomData,
            },
            b1: Chain {
                builder: self.builder,
                producer: self.producer,
                branch: BranchIdx(1),
                _phantom: PhantomData,
            },
            b2: Chain {
                builder: self.builder,
                producer: self.producer,
                branch: BranchIdx(2),
                _phantom: PhantomData,
            },
            b3: Chain {
                builder: self.builder,
                producer: self.producer,
                branch: BranchIdx(3),
                _phantom: PhantomData,
            },
            _phantom: PhantomData,
        }
    }
}

#[must_use = "all chain branches must be wired to a sink"]
pub struct MultiChain2<'b, A: Send + HeapSize + 'static, B: Send + HeapSize + 'static> {
    pub b0: Chain<'b, Single<A>>,
    pub b1: Chain<'b, Single<B>>,
    pub(crate) _phantom: PhantomData<&'b PipelineBuilder>,
}

impl<'b, A: Send + HeapSize + 'static, B: Send + HeapSize + 'static> MultiChain2<'b, A, B> {
    /// Construct a `MultiChain2` from two independent source-side
    /// chains that each produce `Single<A>` / `Single<B>`. Used when
    /// two parallel source subchains converge at a `Step2` consumer
    /// (e.g. zipper's mapped + unmapped BAM source subchains, AAM's
    /// aligner-output + original-record-buffer subchains).
    ///
    /// Mirrors [`Chain::into_multi`]'s output-side counterpart: that
    /// method takes one producer with two output branches and splits
    /// it into a `MultiChain2`; this method takes two distinct
    /// producers (each with one output branch) and packages them.
    ///
    /// # Panics
    ///
    /// Panics if `a` and `b` are anchored at different
    /// `PipelineBuilder` instances (i.e. the framework would have to
    /// wire across pipelines, which is meaningless).
    #[must_use = "pipeline branches must be wired to a sink"]
    pub fn from_chains(a: Chain<'b, Single<A>>, b: Chain<'b, Single<B>>) -> Self {
        assert!(
            std::ptr::eq(a.builder, b.builder),
            "MultiChain2::from_chains: both chains must be anchored at the same builder",
        );
        Self { b0: a, b1: b, _phantom: PhantomData }
    }

    /// Join two parallel sub-chains into a single [`Step2`] consumer.
    ///
    /// Wires `b0` into the consumer's input slot 0
    /// ([`StepCtx2::a`]) and `b1` into input slot 1
    /// ([`StepCtx2::b`]), registering the consumer with input arity
    /// 2 in the chain graph. Each branch's typed producer-side
    /// [`OutputQueueSet`] is later (at chain-run time) drained into
    /// the consumer's
    /// [`crate::pipeline::core::handles::TwoInputHandles<A, B>`]
    /// via [`TypedStep2::build_two_input_handles`].
    ///
    /// Returns a single-branch downstream [`Chain`] typed by the
    /// joined step's `S::Outputs`, ready to chain further single-input
    /// steps onto.
    ///
    /// # Type bounds
    ///
    /// `S: Step2<InputA = A, InputB = B>` — the joined step's input
    /// types must match the two upstream branch types exactly.
    ///
    /// # Panics
    ///
    /// `wire_to_slot`'s defensive panic fires if either upstream
    /// branch was already wired (the `Chain` move semantics prevent
    /// this in well-formed code).
    pub fn join<S>(self, step: S) -> Chain<'b, S::Outputs>
    where
        S: Step2<InputA = A, InputB = B>,
    {
        let builder = self.b0.builder;
        let p0 = self.b0.producer;
        let br0 = self.b0.branch;
        let p1 = self.b1.producer;
        let br1 = self.b1.branch;

        let mut inner = builder.inner.borrow_mut();
        let consumer =
            inner.graph.register_step_with_input_arity(step.profile().name, S::Outputs::arity(), 2);
        inner.graph.wire_to_slot(p0, br0, consumer, 0);
        inner.graph.wire_to_slot(p1, br1, consumer, 1);
        inner.steps.push(Box::new(TypedStep2::new(step)));

        Chain { builder, producer: consumer, branch: BranchIdx(0), _phantom: PhantomData }
    }
}

/// Ordered-bytes variant of `MultiChain2`. Each branch is typed as
/// `Chain<OrderedBytesSingle<_>>` so downstream chained steps see the
/// byte-aware ordered representation (matching the queue topology
/// `OrderedBytesTuple2` actually constructed).
#[must_use = "all chain branches must be wired to a sink"]
pub struct MultiChain2Ordered<'b, A, B>
where
    A: Send + super::item::HeapSize + super::item::Ordered + 'static,
    B: Send + super::item::HeapSize + super::item::Ordered + 'static,
{
    pub b0: Chain<'b, super::outputs::OrderedBytesSingle<A>>,
    pub b1: Chain<'b, super::outputs::OrderedBytesSingle<B>>,
    pub(crate) _phantom: PhantomData<&'b PipelineBuilder>,
}

impl<'b, A, B> MultiChain2Ordered<'b, A, B>
where
    A: Send + super::item::HeapSize + super::item::Ordered + 'static,
    B: Send + super::item::HeapSize + super::item::Ordered + 'static,
{
    /// Construct a `MultiChain2Ordered` from two independent
    /// source-side chains that each produce `OrderedBytesSingle<A>` /
    /// `OrderedBytesSingle<B>`. The ordered/byte-bounded counterpart
    /// of [`MultiChain2::from_chains`] — used when two parallel
    /// source subchains coming out of ordered, byte-bounded BAM/FASTQ
    /// step libraries (decompress → boundaries → decode → group) need
    /// to converge at a [`Step2`] consumer.
    ///
    /// Anchored at the same `PipelineBuilder` as both input chains;
    /// panics otherwise.
    ///
    /// # Panics
    ///
    /// Panics if `a` and `b` are anchored at different
    /// `PipelineBuilder` instances.
    #[must_use = "pipeline branches must be wired to a sink"]
    pub fn from_chains(
        a: Chain<'b, super::outputs::OrderedBytesSingle<A>>,
        b: Chain<'b, super::outputs::OrderedBytesSingle<B>>,
    ) -> Self {
        assert!(
            std::ptr::eq(a.builder, b.builder),
            "MultiChain2Ordered::from_chains: both chains must be anchored at the same builder",
        );
        Self { b0: a, b1: b, _phantom: PhantomData }
    }

    /// Join two parallel ordered/byte-bounded sub-chains into a single
    /// [`Step2`] consumer. Ordered counterpart of
    /// [`MultiChain2::join`]: wires `b0` into the consumer's input
    /// slot 0 ([`StepCtx2::a`]) and `b1` into input slot 1
    /// ([`StepCtx2::b`]), registering the consumer with input arity 2.
    ///
    /// Returns a single-branch downstream [`Chain`] typed by the
    /// joined step's `S::Outputs`.
    ///
    /// # Type bounds
    ///
    /// `S: Step2<InputA = A, InputB = B>` — the joined step's input
    /// types must match the two upstream branch element types exactly.
    /// The framework only requires `HeapSize` on `Step2::InputA` /
    /// `Step2::InputB`; the `Ordered` bound carried by
    /// `OrderedBytesSingle` is upstream-side typing and does not flow
    /// into the consumer's input handle (steps see plain
    /// `InputHandle<T>` regardless of upstream queue topology).
    ///
    /// # Panics
    ///
    /// `wire_to_slot`'s defensive panic fires if either upstream
    /// branch was already wired (the `Chain` move semantics prevent
    /// this in well-formed code).
    pub fn join<S>(self, step: S) -> Chain<'b, S::Outputs>
    where
        S: Step2<InputA = A, InputB = B>,
    {
        let builder = self.b0.builder;
        let p0 = self.b0.producer;
        let br0 = self.b0.branch;
        let p1 = self.b1.producer;
        let br1 = self.b1.branch;

        let mut inner = builder.inner.borrow_mut();
        let consumer =
            inner.graph.register_step_with_input_arity(step.profile().name, S::Outputs::arity(), 2);
        inner.graph.wire_to_slot(p0, br0, consumer, 0);
        inner.graph.wire_to_slot(p1, br1, consumer, 1);
        inner.steps.push(Box::new(TypedStep2::new(step)));

        Chain { builder, producer: consumer, branch: BranchIdx(0), _phantom: PhantomData }
    }
}

#[must_use = "all chain branches must be wired to a sink"]
pub struct MultiChain3<'b, A, B, C>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
{
    pub b0: Chain<'b, Single<A>>,
    pub b1: Chain<'b, Single<B>>,
    pub b2: Chain<'b, Single<C>>,
    pub(crate) _phantom: PhantomData<&'b PipelineBuilder>,
}

#[must_use = "all chain branches must be wired to a sink"]
pub struct MultiChain4<'b, A, B, C, D>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
    D: Send + HeapSize + 'static,
{
    pub b0: Chain<'b, Single<A>>,
    pub b1: Chain<'b, Single<B>>,
    pub b2: Chain<'b, Single<C>>,
    pub b3: Chain<'b, Single<D>>,
    pub(crate) _phantom: PhantomData<&'b PipelineBuilder>,
}

/// A built pipeline ready to run.
///
/// `steps` is consumed by `Pipeline::run`; `graph` is read by `run` and
/// `dag()`; `signal` backs both `cancel_handle` and `run`'s outcome
/// plumbing.
pub struct Pipeline {
    pub(crate) steps: Vec<Box<dyn ErasedStep>>,
    pub(crate) graph: ChainGraph,
    pub(crate) signal: Arc<PipelineSignal>,
}

impl Pipeline {
    #[must_use]
    pub fn builder() -> PipelineBuilder {
        PipelineBuilder::new()
    }

    #[must_use]
    pub fn cancel_handle(&self) -> CancelHandle {
        CancelHandle::from_signal(Arc::clone(&self.signal))
    }

    /// Construct a fresh `PipelineStats` collector sized to this pipeline's
    /// chain. Wrap the returned `Arc` and pass it into `PipelineConfig::stats`
    /// (or via `PipelineConfig::with_stats`) before calling `run`. Counters
    /// can be read at any time after the run completes.
    #[must_use]
    pub fn stats(&self) -> Arc<PipelineStats> {
        let names: Vec<&'static str> = self.steps.iter().map(|s| s.profile().name).collect();
        Arc::new(PipelineStats::new(names))
    }

    /// Render the chain shape as a multi-line debug string. Lists each step
    /// in chain order with its profile (kind, sticky, branch count) and
    /// the consumer for each output branch. Used for diagnostics and for
    /// runall-style `--explain` output.
    ///
    /// The output isn't a stable serialization format — it's a developer-
    /// readable summary, intended to be `println!`'d during debugging or
    /// embedded in error messages.
    #[must_use]
    pub fn dag(&self) -> String {
        use std::fmt::Write as _;

        let mut s = String::new();
        let _ = writeln!(
            s,
            "Pipeline DAG ({} step{}):",
            self.steps.len(),
            if self.steps.len() == 1 { "" } else { "s" }
        );
        for (idx, step) in self.steps.iter().enumerate() {
            let profile = step.profile();
            let n_branches = self.graph.branch_count(super::topology::StepIdx(idx));
            let _ = write!(
                s,
                "  [{idx}] {name:<24} {kind:?} sticky={sticky} branches={n_branches}",
                idx = idx,
                name = profile.name,
                kind = profile.kind,
                sticky = profile.sticky,
                n_branches = n_branches,
            );
            if n_branches == 0 {
                let _ = writeln!(s, " (sink)");
            } else {
                let _ = writeln!(s);
                for branch_usize in 0..n_branches {
                    let branch = super::topology::BranchIdx(branch_usize);
                    let consumer = self.graph.consumer(super::topology::StepIdx(idx), branch);
                    let consumer_name = match consumer {
                        Some(c) => self.graph.step_name(c),
                        None => "<unwired>",
                    };
                    let queue_spec = profile
                        .output_queues
                        .get(branch_usize)
                        .copied()
                        .unwrap_or(super::queues::QueueSpec::Unbounded);
                    let ordering = profile
                        .branch_ordering
                        .get(branch_usize)
                        .copied()
                        .unwrap_or(super::reorder::BranchOrdering::None);
                    let _ = writeln!(
                        s,
                        "      .{branch_usize}: {queue_spec:?} {ordering:?} → {consumer_name}",
                    );
                }
            }
        }
        s
    }

    /// Run the pipeline to completion.
    ///
    /// Spawns `config.threads` worker threads, runs each step's `try_run`
    /// until every step has reported `Finished`, joins on completion, and
    /// returns `Ok(())` on clean exit or an `Err(PipelineError)` if any step
    /// returned `Err` or the caller cancelled via the `CancelHandle`.
    ///
    /// # Errors
    ///
    /// Returns `PipelineError::NotEnoughThreads` if the chain has more
    /// `Exclusive` steps than `config.threads`. Returns `PipelineError::Io`
    /// if any step's `try_run` returned `Err`.
    /// Returns `PipelineError::Cancelled` if `cancel_handle().cancel()` was
    /// called during the run.
    ///
    /// # Panics
    ///
    /// Panics if a worker thread panics — the panic is propagated via
    /// `JoinHandle::join`. Worker panics indicate a framework or step bug
    /// (e.g., a contract violation that triggered a `debug_assert!`).
    #[allow(clippy::too_many_lines)]
    pub fn run(self, config: PipelineConfig) -> Result<(), super::signal::PipelineError> {
        use std::thread;

        use super::runtime::{
            StepDrainCounter, WorkerCore, assign_exclusive_owners, assign_sticky_owners,
            build_chain_contexts, build_worker_storage, is_fusible_chain, run_fused_single_thread,
            run_worker_loop,
        };
        use super::signal::PipelineError;
        use super::step::StepKind;
        use super::topology::StepIdx;

        let Self { steps, graph, signal } = self;
        let n_threads = config.threads;
        let stats_arc = config.stats;
        let deadlock_timeout_secs = config.deadlock_timeout_secs;
        assert!(n_threads > 0, "PipelineConfig::threads must be > 0");
        if let Some(stats) = stats_arc.as_ref() {
            assert_eq!(
                stats.n_steps(),
                steps.len(),
                "PipelineConfig::stats was sized for {} steps but pipeline has {}; \
                 obtain the stats handle from `Pipeline::stats()` after `build()`",
                stats.n_steps(),
                steps.len()
            );
        }
        if deadlock_timeout_secs > 0 && stats_arc.is_none() {
            log::warn!(
                "PipelineConfig::deadlock_timeout_secs is {deadlock_timeout_secs} but \
                 PipelineConfig::stats is None; deadlock monitor cannot run without \
                 stats. Disable one or pair them via the helpers in commands/common.rs."
            );
        }

        // 0. Fused single-thread fast path (issue #330). A single-source
        // source→sink chain at one worker is driven inline over direct,
        // unbounded buffers, skipping the scheduler's round-robin poll /
        // contention / reorder overhead (~2/3 of `try_run` calls at t=1 are
        // otherwise wasted). Fan-out is allowed (e.g. the `--rejects` kept/
        // rejects split); only two-input `Step2` merges (zipper, align) and
        // `--threads ≥ 2` fall through to the scheduled worker pool below. The
        // deadlock monitor / queue rebalancer are not spawned: a single-worker
        // inline drive cannot deadlock and has no bounded queues to rebalance.
        if n_threads == 1 && is_fusible_chain(&steps, &graph) {
            log::debug!(
                "Using fused single-thread pipeline ({} steps, direct buffers)",
                steps.len()
            );
            let result = run_fused_single_thread(steps, &graph, &signal, stats_arc.as_ref());
            // Same end-of-run stats snapshot the scheduled path emits, so
            // `--pipeline-stats` shows the fused chain's per-step counters.
            if let Some(stats) = stats_arc.as_ref() {
                let snapshot = stats.snapshot();
                log::info!("Pipeline end-of-run stats:");
                for line in format!("{snapshot}").lines() {
                    log::info!("{line}");
                }
            }
            return result;
        }

        // 1. Assign Exclusive owners; bail if too many.
        let owners = assign_exclusive_owners(&steps, n_threads)?;

        // 1a. Compute per-worker sticky-driven step (Exclusive sticky union
        // with Serial+sticky+Affinity targeting). Indexed by worker id.
        let sticky_owners = assign_sticky_owners(&steps, &owners, n_threads);

        // 2. Build per-step contexts (input + output handles).
        let contexts = Arc::new(build_chain_contexts(&steps, &graph));

        // 2a. If a total queue-memory budget was supplied, evenly
        // redistribute it across all byte-bounded queues now (before
        // workers start) so the initial state matches the user's
        // budget instead of the per-step defaults baked into each
        // step's `QueueSpec::ByteBounded { limit_bytes }`. Floor at
        // 1 MiB per queue to prevent zero-budget queues that would
        // always reject pushes.
        if let Some(total) = config.queue_memory_total {
            apply_initial_queue_budget(&contexts.bounded_queues, total);
        }

        // 3. Per-step drain counter — init N for Parallel, 1 otherwise.
        let drain_counters: Vec<Arc<StepDrainCounter>> = steps
            .iter()
            .map(|step| {
                let initial = match step.profile().kind {
                    StepKind::Parallel => n_threads,
                    StepKind::Serial | StepKind::Exclusive => 1,
                };
                StepDrainCounter::new(initial)
            })
            .collect();

        // 4. Build per-worker step storage (consumes `steps`).
        let mut worker_entries = build_worker_storage(steps, &owners, n_threads);

        let signal_arc = Arc::clone(&signal);

        // 4a. Optional deadlock-detection monitor. Spawns a watcher
        // thread that periodically samples the stats snapshot; if no
        // step's `progress + finished` counter has advanced for
        // `deadlock_timeout_secs` seconds, it logs the snapshot at
        // `warn` level so the user has a starting point. Polls every
        // `max(1, deadlock_timeout_secs / 4)` seconds.
        let (monitor_stop, monitor_handle) = match (deadlock_timeout_secs, stats_arc.as_ref()) {
            (n, Some(stats)) if n > 0 => {
                let stop = Arc::new(std::sync::atomic::AtomicBool::new(false));
                let stop_clone = Arc::clone(&stop);
                let stats_clone = Arc::clone(stats);
                let timeout = std::time::Duration::from_secs(deadlock_timeout_secs);
                let poll_interval =
                    std::time::Duration::from_secs(deadlock_timeout_secs.max(4) / 4);
                let handle = thread::Builder::new()
                    .name("fgumi-deadlock-monitor".to_string())
                    .spawn(move || {
                        run_deadlock_monitor(&stop_clone, &stats_clone, timeout, poll_interval);
                    })
                    .expect("failed to spawn deadlock monitor thread");
                (Some(stop), Some(handle))
            }
            _ => (None, None),
        };

        // 4b. Optional queue-memory rebalancer. Spawns a watcher
        // thread that periodically samples each registered queue's
        // `current_bytes / limit_bytes` ratio and shifts budget
        // toward consistently-full queues at the expense of
        // consistently-empty ones. Total budget is preserved.
        let (rebalancer_stop, rebalancer_handle) =
            if config.queue_memory_total.is_some() && !contexts.bounded_queues.is_empty() {
                let stop = Arc::new(std::sync::atomic::AtomicBool::new(false));
                let stop_clone = Arc::clone(&stop);
                // Capture handles into a contiguous Vec — the monitor
                // doesn't need step indices/names beyond debug logging.
                let handles: Vec<Arc<dyn super::queues::BoundedQueueHandle>> =
                    contexts.bounded_queues.iter().map(|rq| Arc::clone(&rq.handle)).collect();
                let names: Vec<&'static str> =
                    contexts.bounded_queues.iter().map(|rq| rq.producer_step_name).collect();
                let handle = thread::Builder::new()
                    .name("fgumi-queue-rebalancer".to_string())
                    .spawn(move || {
                        run_queue_rebalancer(&stop_clone, &handles, &names);
                    })
                    .expect("failed to spawn queue rebalancer thread");
                (Some(stop), Some(handle))
            } else {
                (None, None)
            };

        if n_threads == 1 {
            // Single-threaded fast path: run the worker loop directly on
            // the caller's thread instead of spawning + joining a fresh
            // OS thread. The framework machinery (`build_worker_storage`,
            // drain counters, stats) is identical to the multi-threaded
            // path — only the spawn/join is skipped.
            //
            // Savings: thread-spawn-and-join (a few ms one-time on Apple
            // Silicon / Linux), and the caller's thread name / TLS is
            // preserved (matters for log correlation in some tools).
            //
            // Per-call cost gap vs the legacy single-threaded path is
            // dominated by the framework's per-record book-keeping
            // (queue byte tracking, ordinal allocation, drain checks),
            // not the spawn — see commit message benchmarks.
            let entries = worker_entries
                .pop()
                .expect("build_worker_storage with n_threads=1 returns one entry vec");
            debug_assert!(worker_entries.is_empty());
            let exclusive_owner = owners
                .iter()
                .enumerate()
                .find_map(|(idx, &own)| if own == Some(0) { Some(StepIdx(idx)) } else { None });
            let sticky_owner = sticky_owners[0];
            let mut worker = WorkerCore::new(0, exclusive_owner, sticky_owner);
            let mut entries_local = entries;
            run_worker_loop(
                &mut worker,
                &mut entries_local,
                &contexts,
                &drain_counters,
                &signal_arc,
                stats_arc.as_ref(),
            );
        } else {
            // 5. Spawn worker threads.
            let mut handles = Vec::with_capacity(n_threads);
            for (worker_id, entries) in worker_entries.into_iter().enumerate() {
                let exclusive_owner = owners.iter().enumerate().find_map(|(idx, &own)| {
                    if own == Some(worker_id) { Some(StepIdx(idx)) } else { None }
                });
                let sticky_owner = sticky_owners[worker_id];

                let contexts_clone = Arc::clone(&contexts);
                let signal_clone = Arc::clone(&signal_arc);
                let drain_counters_clone: Vec<Arc<StepDrainCounter>> =
                    drain_counters.iter().map(Arc::clone).collect();
                let stats_clone = stats_arc.as_ref().map(Arc::clone);

                let handle = thread::Builder::new()
                    .name(format!("fgumi-worker-{worker_id}"))
                    .spawn(move || {
                        let mut worker = WorkerCore::new(worker_id, exclusive_owner, sticky_owner);
                        let mut entries_local = entries;
                        run_worker_loop(
                            &mut worker,
                            &mut entries_local,
                            &contexts_clone,
                            &drain_counters_clone,
                            &signal_clone,
                            stats_clone.as_ref(),
                        );
                    })
                    .expect("failed to spawn worker thread");
                handles.push(handle);
            }

            // 6. Join workers. If a worker panicked, re-raise its original
            // payload via `resume_unwind` so the main thread aborts with the
            // worker's actual panic message and location — not the opaque
            // `Any { .. }` that `join().expect(...)` would print.
            for h in handles {
                if let Err(panic) = h.join() {
                    std::panic::resume_unwind(panic);
                }
            }
        }

        // 6a. Stop and join the deadlock monitor (if spawned). We
        // signal stop *after* workers join so the monitor sees the
        // final stats state and doesn't fire spurious warnings during
        // normal pipeline shutdown (where steps stop progressing
        // because they're done, not stuck).
        if let Some(stop) = monitor_stop {
            stop.store(true, std::sync::atomic::Ordering::Relaxed);
        }
        if let Some(handle) = monitor_handle {
            let _ = handle.join();
        }

        // 6b. Stop and join the queue rebalancer (if spawned).
        if let Some(stop) = rebalancer_stop {
            stop.store(true, std::sync::atomic::Ordering::Relaxed);
        }
        if let Some(handle) = rebalancer_handle {
            let _ = handle.join();
        }

        // 6c. End-of-run stats snapshot — emitted at info-level so
        // bench profiling can see which step accumulated the most
        // wall-clock time inside `try_run`. Only logged when stats
        // were enabled (`with_stats(...)` on the config); cheap
        // either way.
        if let Some(stats) = stats_arc.as_ref() {
            let snapshot = stats.snapshot();
            log::info!("Pipeline end-of-run stats:");
            for line in format!("{snapshot}").lines() {
                log::info!("{line}");
            }
        }

        // 7. Surface error or cancellation. PipelineError isn't Clone
        // (io::Error isn't Clone); reconstruct from the recorded outcome.
        let _ = graph;
        match signal.outcome() {
            Some(PipelineError::Cancelled) => Err(PipelineError::Cancelled),
            Some(PipelineError::Io { step, source }) => Err(PipelineError::Io {
                step,
                source: std::io::Error::new(source.kind(), format!("{source}")),
            }),
            Some(PipelineError::NotEnoughThreads { required, available }) => {
                Err(PipelineError::NotEnoughThreads { required: *required, available: *available })
            }
            None => Ok(()),
        }
    }
}

/// Background deadlock monitor body. Polls `stats` every `poll_interval`
/// and tracks the cumulative `progress + finished` counter across all
/// steps. If no step has advanced for `timeout`, logs the snapshot at
/// `warn` level (once per stall window) and resets the watermark so a
/// long-running stall keeps emitting periodic diagnostics rather than
/// only one. Exits when `stop` is set (typically by the main thread
/// after workers join).
fn run_deadlock_monitor(
    stop: &Arc<std::sync::atomic::AtomicBool>,
    stats: &Arc<PipelineStats>,
    timeout: std::time::Duration,
    poll_interval: std::time::Duration,
) {
    let mut last_advance = std::time::Instant::now();
    let mut last_total = total_progress(stats);
    while !stop.load(std::sync::atomic::Ordering::Relaxed) {
        sleep_until_stop(stop, poll_interval);
        if stop.load(std::sync::atomic::Ordering::Relaxed) {
            break;
        }
        let now_total = total_progress(stats);
        if now_total != last_total {
            last_total = now_total;
            last_advance = std::time::Instant::now();
            continue;
        }
        if last_advance.elapsed() >= timeout {
            // Stalled — emit diagnostics. Reset the watermark so we
            // continue logging on each `timeout` window if the stall
            // persists, instead of going silent after one warning.
            log::warn!(
                "Pipeline stall detected: no step has progressed in {:.0}s. \
                 Snapshot follows.",
                last_advance.elapsed().as_secs_f64()
            );
            let snapshot = stats.snapshot();
            for line in format!("{snapshot}").lines() {
                log::warn!("{line}");
            }
            last_advance = std::time::Instant::now();
        }
    }
}

/// Sum of `progress_count + finished_count` across all steps. The
/// monitor uses this as a single scalar progress watermark; a change
/// means *some* step is making forward progress.
fn total_progress(stats: &PipelineStats) -> u64 {
    let snap = stats.snapshot();
    snap.steps.iter().map(|(_, s)| s.progress_count + s.finished_count).sum()
}

/// Sleep up to `dur`, returning early as soon as `stop` is set. Polls the
/// flag in short slices so a background helper thread (deadlock monitor,
/// queue rebalancer) exits within tens of milliseconds at teardown instead
/// of blocking the main thread's `join()` for a full poll interval after
/// the pipeline has already finished. A plain `thread::sleep(poll_interval)`
/// here adds a fixed dead-time tail (up to `poll_interval`) to every run —
/// negligible on long jobs but a large *relative* regression on short ones
/// (e.g. FASTQ extract), since the worker pool is already idle and waiting.
fn sleep_until_stop(stop: &std::sync::atomic::AtomicBool, dur: std::time::Duration) {
    const SLICE: std::time::Duration = std::time::Duration::from_millis(25);
    let deadline = std::time::Instant::now() + dur;
    loop {
        if stop.load(std::sync::atomic::Ordering::Relaxed) {
            return;
        }
        let now = std::time::Instant::now();
        if now >= deadline {
            return;
        }
        std::thread::sleep(SLICE.min(deadline - now));
    }
}

/// Per-queue floor: never let the rebalancer take a queue below this
/// (`ByteBoundedQueue` panics if `limit_bytes == 0`, and very small
/// limits effectively wedge the producer).
const MIN_PER_QUEUE_BYTES: u64 = 1024 * 1024;

/// Floor for the per-branch reorder overflow stash. Liveness needs no floor
/// (`next_serial` is always exempt — any cap ≥ 0 is deadlock-free), so this
/// is purely a throughput knob: keep enough lookahead headroom that a
/// reorder-heavy edge doesn't thrash one item per round-robin pass.
const MIN_REORDER_OVERFLOW_BYTES: u64 = 4 * 1024 * 1024;

/// Initial-allocation pass for `queue_memory_total`. Distributes
/// `total` evenly across all byte-bounded queues. Floors each queue
/// at `MIN_PER_QUEUE_BYTES` even if the per-queue share would be
/// smaller — in that case the effective total exceeds the user's
/// budget, but starvation is the worse failure mode.
///
/// Each ordered byte-bounded branch's reorder overflow stash is sized from
/// the SAME `per_queue` value (clamped to `[MIN_REORDER_OVERFLOW_BYTES,
/// DEFAULT_REORDER_OVERFLOW_BYTES]`), so the off-budget stash tracks the
/// transport budget instead of a fixed 256 MiB. The clamp ceiling is the
/// prior fixed value, so high thread counts (large `per_queue`) keep today's
/// reorder headroom — no `--threads N` regression — while low thread counts
/// (small `per_queue`) get a streaming-sized stash.
fn apply_initial_queue_budget(
    queues: &[crate::pipeline::core::runtime::contexts::RegisteredQueue],
    total: u64,
) {
    if queues.is_empty() {
        return;
    }
    let per_queue = (total / (queues.len() as u64)).max(MIN_PER_QUEUE_BYTES);
    let reorder_cap = reorder_cap_for(per_queue);
    for rq in queues {
        rq.handle.set_limit_bytes(per_queue);
        if let Some(reorder) = &rq.reorder_cap {
            reorder.set_max_overflow_bytes(reorder_cap);
        }
    }
}

/// The reorder overflow cap for a branch whose transport budget is
/// `per_queue`: track the transport budget, clamped to
/// `[MIN_REORDER_OVERFLOW_BYTES, DEFAULT_REORDER_OVERFLOW_BYTES]`. The
/// ceiling is the prior fixed value, so a large `per_queue` (high thread
/// counts) keeps today's reorder headroom; a small `per_queue` (low thread
/// counts / lean budget) shrinks the off-budget stash to a streaming size.
fn reorder_cap_for(per_queue: u64) -> u64 {
    per_queue.clamp(
        MIN_REORDER_OVERFLOW_BYTES,
        crate::pipeline::core::reorder::DEFAULT_REORDER_OVERFLOW_BYTES,
    )
}

/// Background queue-memory rebalancer body. Polls each queue's
/// `current_bytes / limit_bytes` fullness ratio every 1 second.
/// Identifies the most-full producer (likely bottleneck) and the
/// least-full consumer (over-budget). Shifts a fraction of budget
/// from least to most full, preserving total budget.
///
/// The algorithm is deliberately simple — incremental shifts (10%
/// of the source's limit per tick) converge gradually so transient
/// spikes don't overshoot. Floors each queue at `MIN_PER_QUEUE_BYTES`.
///
/// Exits when `stop` is set (typically after workers join).
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn run_queue_rebalancer(
    stop: &Arc<std::sync::atomic::AtomicBool>,
    handles: &[Arc<dyn super::queues::BoundedQueueHandle>],
    names: &[&'static str],
) {
    if handles.len() < 2 {
        // Nothing to rebalance with one or zero queues.
        return;
    }
    let poll_interval = std::time::Duration::from_secs(1);
    let shift_fraction: f64 = 0.10;

    while !stop.load(std::sync::atomic::Ordering::Relaxed) {
        sleep_until_stop(stop, poll_interval);
        if stop.load(std::sync::atomic::Ordering::Relaxed) {
            break;
        }

        // Snapshot fullness ratios.
        let snapshot: Vec<(usize, u64, u64, f64)> = handles
            .iter()
            .enumerate()
            .map(|(idx, h)| {
                let cur = h.current_bytes();
                let lim = h.limit_bytes();
                let ratio = if lim == 0 { 0.0 } else { (cur as f64) / (lim as f64) };
                (idx, cur, lim, ratio)
            })
            .collect();

        // Find the most-full and least-full queues.
        let max = snapshot
            .iter()
            .max_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal))
            .copied();
        let min = snapshot
            .iter()
            .min_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal))
            .copied();

        let (Some((max_idx, _, max_lim, max_ratio)), Some((min_idx, _, min_lim, min_ratio))) =
            (max, min)
        else {
            continue;
        };
        if max_idx == min_idx {
            continue;
        }
        // Only rebalance when the imbalance is meaningful: the
        // fullest queue is ≥80% full AND the emptiest is ≤20% full.
        // Otherwise the system is in steady state and we shouldn't
        // perturb the limits.
        if max_ratio < 0.80 || min_ratio > 0.20 {
            continue;
        }

        // Shift from min to max.
        let to_shift = ((min_lim as f64) * shift_fraction) as u64;
        if to_shift == 0 {
            continue;
        }
        let new_min = min_lim.saturating_sub(to_shift).max(MIN_PER_QUEUE_BYTES);
        if new_min == min_lim {
            // Floor reached; can't shrink further.
            continue;
        }
        let actual_shift = min_lim - new_min;
        let new_max = max_lim.saturating_add(actual_shift);
        handles[min_idx].set_limit_bytes(new_min);
        handles[max_idx].set_limit_bytes(new_max);
        log::debug!(
            "queue rebalance: shift {} bytes {} ({} -> {}) -> {} ({} -> {})",
            actual_shift,
            names[min_idx],
            min_lim,
            new_min,
            names[max_idx],
            max_lim,
            new_max
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    use crate::pipeline::core::outputs::Single;
    use crate::pipeline::core::queues::QueueSpec;
    use crate::pipeline::core::reorder::BranchOrdering;
    use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

    #[test]
    fn apply_initial_queue_budget_sets_registered_reorder_cap() {
        // End-to-end wiring: a registered ordered byte-bounded branch's reorder
        // cap is re-sized by the budget pass (not left at its construction
        // default). Guards the Pass-1.5 registration + the `set` call together.
        use crate::pipeline::core::queues::{BoundedQueueHandle, ByteBoundedQueue, ItemQueue};
        use crate::pipeline::core::reorder::{
            DEFAULT_REORDER_OVERFLOW_BYTES, ReorderCapHandle, ReorderStage, Sequenced,
        };
        use crate::pipeline::core::runtime::contexts::RegisteredQueue;
        use crate::pipeline::core::topology::{BranchIdx, StepIdx};

        // A reorder stage constructed at the 256 MiB fallback; keep a concrete
        // handle so we can read the cap back after the budget pass.
        let transport: Arc<dyn ItemQueue<Sequenced<u32>>> =
            Arc::new(ByteBoundedQueue::<Sequenced<u32>>::new(1024 * 1024));
        let stage = Arc::new(ReorderStage::<u32>::with_max_overflow_bytes(
            transport,
            DEFAULT_REORDER_OVERFLOW_BYTES,
        ));
        assert_eq!(stage.current_max_overflow_bytes(), DEFAULT_REORDER_OVERFLOW_BYTES);
        let reorder_dyn: Arc<dyn ReorderCapHandle> = stage.clone();

        // A transport-limit handle for the `RegisteredQueue.handle` slot.
        let transport_q = Arc::new(ByteBoundedQueue::<u32>::new(1024 * 1024));
        let transport_handle: Arc<dyn BoundedQueueHandle> = transport_q;

        let registered = vec![RegisteredQueue {
            producer_step_name: "TestStep",
            producer_step: StepIdx(0),
            branch: BranchIdx(0),
            handle: transport_handle,
            reorder_cap: Some(reorder_dyn),
        }];

        // Lean total → per_queue = 8 MiB (1 queue) → reorder clamp = 8 MiB.
        apply_initial_queue_budget(&registered, 8 * 1024 * 1024);
        assert_eq!(
            stage.current_max_overflow_bytes(),
            8 * 1024 * 1024,
            "budget pass must re-size the registered reorder cap to the clamped per_queue"
        );

        // Huge total → per_queue huge → reorder clamped back to the ceiling.
        apply_initial_queue_budget(&registered, 100 * 1024 * 1024 * 1024);
        assert_eq!(
            stage.current_max_overflow_bytes(),
            DEFAULT_REORDER_OVERFLOW_BYTES,
            "high budget clamps the reorder cap to the 256 MiB ceiling (no t>1 regression)"
        );
    }

    #[test]
    fn reorder_cap_tracks_per_queue_clamped_to_floor_and_ceiling() {
        let ceiling = crate::pipeline::core::reorder::DEFAULT_REORDER_OVERFLOW_BYTES;
        // Mid-range per_queue passes through unchanged.
        assert_eq!(reorder_cap_for(32 * 1024 * 1024), 32 * 1024 * 1024);
        // Tiny per_queue (lean / low-thread) is floored — but stays small.
        assert_eq!(reorder_cap_for(1024), MIN_REORDER_OVERFLOW_BYTES);
        assert_eq!(reorder_cap_for(MIN_REORDER_OVERFLOW_BYTES - 1), MIN_REORDER_OVERFLOW_BYTES);
        // Huge per_queue (high thread counts) is capped at today's ceiling →
        // no `--threads N` regression.
        assert_eq!(reorder_cap_for(4 * ceiling), ceiling);
        assert_eq!(reorder_cap_for(ceiling), ceiling);
    }

    // ───── Test stubs ─────

    #[derive(Clone)]
    struct StubSource;
    impl Step for StubSource {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Source",
                kind: StepKind::Exclusive,
                sticky: true,
                output_queues: vec![QueueSpec::CountBounded { capacity: 64 }],
                branch_ordering: vec![BranchOrdering::ByOrdinal],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::Finished)
        }
    }

    #[derive(Clone)]
    struct StubTransform;
    impl Step for StubTransform {
        type Input = u32;
        type Outputs = Single<u64>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Transform",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 64 }],
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

    #[derive(Clone)]
    struct StubFanOut2;
    impl Step for StubFanOut2 {
        type Input = u64;
        type Outputs = (u32, String);
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "FanOut2",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![
                    QueueSpec::CountBounded { capacity: 32 },
                    QueueSpec::CountBounded { capacity: 32 },
                ],
                branch_ordering: vec![BranchOrdering::None, BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    #[derive(Clone)]
    struct StubSinkU32;
    impl Step for StubSinkU32 {
        type Input = u32;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "SinkU32",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    #[derive(Clone)]
    struct StubSinkString;
    impl Step for StubSinkString {
        type Input = String;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "SinkString",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    // ───── Tests ─────

    #[test]
    fn empty_builder_returns_empty_error() {
        let builder = PipelineBuilder::new();
        assert!(matches!(builder.build(), Err(BuildError::Empty)));
    }

    #[test]
    fn unwired_source_returns_unwired_output() {
        let builder = PipelineBuilder::new();
        let _chain = builder.chain(StubSource);
        let result = builder.build();
        assert!(matches!(result, Err(BuildError::UnwiredOutput { step: "Source", branch: "0" })));
    }

    #[test]
    fn source_to_transform_unwired_at_transform_returns_unwired() {
        let builder = PipelineBuilder::new();
        let _chain = builder.chain(StubSource).chain(StubTransform);
        let result = builder.build();
        assert!(matches!(
            result,
            Err(BuildError::UnwiredOutput { step: "Transform", branch: "0" })
        ));
    }

    #[test]
    fn fully_wired_source_to_sink_succeeds() {
        let builder = PipelineBuilder::new();
        builder.chain(StubSource).chain(StubSinkU32).into_sink_marker();
        let result = builder.build();
        assert!(result.is_ok());
        let pipeline = result.unwrap();
        assert_eq!(pipeline.graph.n_steps(), 2);
    }

    #[test]
    fn unwired_fanout_branch_is_detected() {
        let builder = PipelineBuilder::new();
        let after_fanout = builder.chain(StubSource).chain(StubTransform).chain(StubFanOut2);
        let multi = after_fanout.into_multi();
        // Wire branch 0 to a sink, drop branch 1 (unwired).
        multi.b0.chain(StubSinkU32).into_sink_marker();
        drop(multi.b1);

        let result = builder.build();
        assert!(matches!(result, Err(BuildError::UnwiredOutput { step: "FanOut2", branch: "1" })));
    }

    #[test]
    fn fanout_with_both_branches_wired_succeeds() {
        let builder = PipelineBuilder::new();
        let after_fanout = builder.chain(StubSource).chain(StubTransform).chain(StubFanOut2);
        let multi = after_fanout.into_multi();
        multi.b0.chain(StubSinkU32).into_sink_marker();
        multi.b1.chain(StubSinkString).into_sink_marker();

        let result = builder.build();
        assert!(result.is_ok());
        let pipeline = result.unwrap();
        assert_eq!(pipeline.graph.n_steps(), 5);
    }

    #[test]
    fn pipeline_config_default_uses_available_parallelism() {
        let cfg = PipelineConfig::default();
        assert!(cfg.threads >= 1);
    }

    #[test]
    fn dag_renders_chain_shape() {
        let builder = PipelineBuilder::new();
        builder.chain(StubSource).chain(StubSinkU32).into_sink_marker();
        let pipeline = builder.build().unwrap();
        let dag = pipeline.dag();
        // Verify shape rendering — names and (sink) marker.
        assert!(dag.contains("Source"), "DAG missing source name: {dag}");
        assert!(dag.contains("SinkU32"), "DAG missing sink name: {dag}");
        assert!(dag.contains("(sink)"), "DAG missing sink marker: {dag}");
        assert!(dag.contains("→ SinkU32"), "DAG missing source→sink wiring: {dag}");
    }

    // ─────────────────────────────────────────────────────────────────────
    // Pipeline::run smoke tests
    // ─────────────────────────────────────────────────────────────────────

    use std::sync::atomic::{AtomicU32, Ordering as AtomicOrd};

    use crate::pipeline::core::signal::PipelineError;

    /// Source emitting `remaining` items via a shared atomic counter; safe
    /// for both single-worker and multi-worker `Parallel` execution.
    ///
    /// Uses an `Unbounded` output queue so the test never hits backpressure
    /// (which would require the source to use the `HeldSlot<Unpushed<T>>`
    /// retry pattern — exercised in the bigger end-to-end smoke test below).
    #[derive(Clone)]
    struct SharedCountingSource {
        remaining: Arc<AtomicU32>,
    }
    impl Step for SharedCountingSource {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "SharedSource",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::Unbounded],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> std::io::Result<StepOutcome> {
            let n = self.remaining.load(AtomicOrd::Acquire);
            if n == 0 {
                return Ok(StepOutcome::Finished);
            }
            // CAS down to claim this item; only push on successful claim.
            if self
                .remaining
                .compare_exchange(n, n - 1, AtomicOrd::AcqRel, AtomicOrd::Acquire)
                .is_ok()
            {
                ctx.outputs.push(n).map_err(|_| {
                    std::io::Error::other("Unbounded queue rejected push (impossible)")
                })?;
                Ok(StepOutcome::Progress)
            } else {
                Ok(StepOutcome::NoProgress)
            }
        }
        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    /// Sink that pops and counts.
    #[derive(Clone)]
    struct ParallelCountingSink {
        received: Arc<AtomicU32>,
    }
    impl Step for ParallelCountingSink {
        type Input = u32;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ParallelSink",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> std::io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(_) => {
                    self.received.fetch_add(1, AtomicOrd::Relaxed);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    #[test]
    fn pipeline_run_with_threads_1_drains_chain() {
        let remaining = Arc::new(AtomicU32::new(10));
        let received = Arc::new(AtomicU32::new(0));

        let builder = PipelineBuilder::new();
        builder
            .chain(SharedCountingSource { remaining: Arc::clone(&remaining) })
            .chain(ParallelCountingSink { received: Arc::clone(&received) })
            .into_sink_marker();
        let pipeline = builder.build().unwrap();

        let result = pipeline.run(PipelineConfig { threads: 1, ..Default::default() });
        assert!(result.is_ok(), "run failed: {:?}", result.err());
        assert_eq!(received.load(AtomicOrd::Relaxed), 10);
    }

    #[test]
    fn pipeline_run_with_threads_4_drains_chain() {
        let remaining = Arc::new(AtomicU32::new(50));
        let received = Arc::new(AtomicU32::new(0));

        let builder = PipelineBuilder::new();
        builder
            .chain(SharedCountingSource { remaining: Arc::clone(&remaining) })
            .chain(ParallelCountingSink { received: Arc::clone(&received) })
            .into_sink_marker();
        let pipeline = builder.build().unwrap();

        let result = pipeline.run(PipelineConfig { threads: 4, ..Default::default() });
        assert!(result.is_ok(), "run failed: {:?}", result.err());
        assert_eq!(received.load(AtomicOrd::Relaxed), 50);
    }

    #[test]
    fn pipeline_stats_handle_matches_chain_size() {
        let builder = PipelineBuilder::new();
        builder.chain(StubSource).chain(StubSinkU32).into_sink_marker();
        let pipeline = builder.build().unwrap();
        let stats = pipeline.stats();
        assert_eq!(stats.n_steps(), 2);
        assert_eq!(stats.step_name(StepIdx(0)), "Source");
        assert_eq!(stats.step_name(StepIdx(1)), "SinkU32");
    }

    #[test]
    fn pipeline_run_populates_stats_when_enabled() {
        let remaining = Arc::new(AtomicU32::new(20));
        let received = Arc::new(AtomicU32::new(0));

        let builder = PipelineBuilder::new();
        builder
            .chain(SharedCountingSource { remaining: Arc::clone(&remaining) })
            .chain(ParallelCountingSink { received: Arc::clone(&received) })
            .into_sink_marker();
        let pipeline = builder.build().unwrap();
        let stats = pipeline.stats();

        let cfg =
            PipelineConfig { threads: 2, stats: Some(Arc::clone(&stats)), ..Default::default() };
        let result = pipeline.run(cfg);
        assert!(result.is_ok(), "run failed: {:?}", result.err());
        assert_eq!(received.load(AtomicOrd::Relaxed), 20);

        let snap = stats.snapshot();
        assert_eq!(snap.steps.len(), 2);
        assert_eq!(snap.steps[0].0, "SharedSource");
        assert_eq!(snap.steps[1].0, "ParallelSink");

        // Sink saw exactly 20 items (one Progress per pop with item).
        assert_eq!(snap.steps[1].1.progress_count, 20);
        // Source must have made progress at least 20 times to push the items.
        assert!(snap.steps[0].1.progress_count >= 20);
        // Both steps accumulated wall time.
        assert!(snap.steps[0].1.total_run_ns > 0);
        assert!(snap.steps[1].1.total_run_ns > 0);
        // No errors recorded.
        assert_eq!(snap.steps[0].1.error_count, 0);
        assert_eq!(snap.steps[1].1.error_count, 0);
        // The source returned Finished at least once across the workers.
        assert!(snap.steps[0].1.finished_count >= 1);
    }

    #[test]
    fn pipeline_run_process2_fans_out_to_two_sinks() {
        use crate::pipeline::steps::process::{Process2Output, process2};

        let remaining = Arc::new(AtomicU32::new(20));
        let evens_received = Arc::new(AtomicU32::new(0));
        let odds_received = Arc::new(AtomicU32::new(0));

        // Process2 routes even items to branch A, odd items to branch B.
        let split = process2::<u32, u32, u32, _>("EvenOddSplit", 32, 32, |x: u32| {
            if x.is_multiple_of(2) {
                Ok(Process2Output::only_a(x))
            } else {
                Ok(Process2Output::only_b(x))
            }
        });

        let builder = PipelineBuilder::new();
        let after_split =
            builder.chain(SharedCountingSource { remaining: Arc::clone(&remaining) }).chain(split);
        let multi = after_split.into_multi();
        multi
            .b0
            .chain(ParallelCountingSink { received: Arc::clone(&evens_received) })
            .into_sink_marker();
        multi
            .b1
            .chain(ParallelCountingSink { received: Arc::clone(&odds_received) })
            .into_sink_marker();

        let pipeline = builder.build().unwrap();
        let result = pipeline.run(PipelineConfig { threads: 4, ..Default::default() });
        assert!(result.is_ok(), "run failed: {:?}", result.err());

        // Source emits values [20, 19, ..., 1] (20 items total). Of these:
        //   evens: 20, 18, ..., 2  -> 10 items
        //   odds : 19, 17, ..., 1  -> 10 items
        assert_eq!(evens_received.load(AtomicOrd::Relaxed), 10, "evens count");
        assert_eq!(odds_received.load(AtomicOrd::Relaxed), 10, "odds count");
    }

    #[test]
    fn pipeline_run_process2_drops_branches_emit_none() {
        use crate::pipeline::steps::process::{Process2Output, process2};

        let remaining = Arc::new(AtomicU32::new(15));
        let kept = Arc::new(AtomicU32::new(0));
        let dropped = Arc::new(AtomicU32::new(0));

        // Filter: keep multiples of 3 on branch A, route the rest to branch B,
        // and drop entirely if the value is 1 (Process2Output::none).
        let filter = process2::<u32, u32, u32, _>("FilterStep", 16, 16, |x: u32| {
            if x == 1 {
                Ok(Process2Output::none())
            } else if x.is_multiple_of(3) {
                Ok(Process2Output::only_a(x))
            } else {
                Ok(Process2Output::only_b(x))
            }
        });

        let builder = PipelineBuilder::new();
        let after =
            builder.chain(SharedCountingSource { remaining: Arc::clone(&remaining) }).chain(filter);
        let multi = after.into_multi();
        multi.b0.chain(ParallelCountingSink { received: Arc::clone(&kept) }).into_sink_marker();
        multi.b1.chain(ParallelCountingSink { received: Arc::clone(&dropped) }).into_sink_marker();

        let pipeline = builder.build().unwrap();
        let result = pipeline.run(PipelineConfig { threads: 4, ..Default::default() });
        assert!(result.is_ok(), "run failed: {:?}", result.err());

        // Source emits [15, 14, ..., 1]. Of these:
        //   multiples of 3 (kept):  15, 12, 9, 6, 3 -> 5 items
        //   value == 1   (dropped): 1               -> 0 items emitted
        //   the rest     (other):   14,13,11,10,8,7,5,4,2 -> 9 items
        assert_eq!(kept.load(AtomicOrd::Relaxed), 5, "kept count");
        assert_eq!(dropped.load(AtomicOrd::Relaxed), 9, "other-branch count");
    }

    #[test]
    fn pipeline_run_without_stats_succeeds_unchanged() {
        let remaining = Arc::new(AtomicU32::new(5));
        let received = Arc::new(AtomicU32::new(0));

        let builder = PipelineBuilder::new();
        builder
            .chain(SharedCountingSource { remaining: Arc::clone(&remaining) })
            .chain(ParallelCountingSink { received: Arc::clone(&received) })
            .into_sink_marker();
        let pipeline = builder.build().unwrap();

        let result = pipeline.run(PipelineConfig { threads: 2, stats: None, ..Default::default() });
        assert!(result.is_ok(), "run failed: {:?}", result.err());
        assert_eq!(received.load(AtomicOrd::Relaxed), 5);
    }

    #[test]
    fn pipeline_run_responds_to_cancellation() {
        /// Source that emits forever (never returns `Finished`) until cancelled.
        #[derive(Clone)]
        struct InfiniteSource;
        impl Step for InfiniteSource {
            type Input = ();
            type Outputs = Single<u32>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "InfiniteSource",
                    kind: StepKind::Parallel,
                    sticky: false,
                    output_queues: vec![QueueSpec::CountBounded { capacity: 8 }],
                    branch_ordering: vec![BranchOrdering::None],
                }
            }
            fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> std::io::Result<StepOutcome> {
                let _ = ctx.outputs.push(0);
                Ok(StepOutcome::Progress)
            }
            fn new_worker_copy(&self) -> Self {
                self.clone()
            }
        }

        #[derive(Clone)]
        struct DiscardSink;
        impl Step for DiscardSink {
            type Input = u32;
            type Outputs = ();
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "DiscardSink",
                    kind: StepKind::Parallel,
                    sticky: false,
                    output_queues: vec![],
                    branch_ordering: vec![],
                }
            }
            fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> std::io::Result<StepOutcome> {
                match ctx.input.pop() {
                    Some(_) => Ok(StepOutcome::Progress),
                    None => Ok(StepOutcome::NoProgress),
                }
            }
            fn new_worker_copy(&self) -> Self {
                self.clone()
            }
        }

        let builder = PipelineBuilder::new();
        builder.chain(InfiniteSource).chain(DiscardSink).into_sink_marker();
        let pipeline = builder.build().unwrap();
        let cancel = pipeline.cancel_handle();

        std::thread::spawn(move || {
            std::thread::sleep(std::time::Duration::from_millis(100));
            cancel.cancel();
        });

        let result = pipeline.run(PipelineConfig { threads: 4, ..Default::default() });
        assert!(matches!(result, Err(PipelineError::Cancelled)));
    }

    #[test]
    fn pipeline_run_propagates_step_error() {
        /// Source that emits `n` items, then returns `Err`.
        #[derive(Clone)]
        struct FailingSource {
            remaining: Arc<AtomicU32>,
        }
        impl Step for FailingSource {
            type Input = ();
            type Outputs = Single<u32>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "FailingSource",
                    kind: StepKind::Parallel,
                    sticky: false,
                    output_queues: vec![QueueSpec::CountBounded { capacity: 8 }],
                    branch_ordering: vec![BranchOrdering::None],
                }
            }
            fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> std::io::Result<StepOutcome> {
                let n = self.remaining.load(AtomicOrd::Acquire);
                if n == 0 {
                    return Err(std::io::Error::other("source failed"));
                }
                if self
                    .remaining
                    .compare_exchange(n, n - 1, AtomicOrd::AcqRel, AtomicOrd::Acquire)
                    .is_ok()
                {
                    let _ = ctx.outputs.push(n);
                    Ok(StepOutcome::Progress)
                } else {
                    Ok(StepOutcome::NoProgress)
                }
            }
            fn new_worker_copy(&self) -> Self {
                self.clone()
            }
        }

        let received = Arc::new(AtomicU32::new(0));
        let remaining = Arc::new(AtomicU32::new(5));
        let builder = PipelineBuilder::new();
        builder
            .chain(FailingSource { remaining: Arc::clone(&remaining) })
            .chain(ParallelCountingSink { received: Arc::clone(&received) })
            .into_sink_marker();
        let pipeline = builder.build().unwrap();

        let result = pipeline.run(PipelineConfig { threads: 4, ..Default::default() });
        match result {
            Err(PipelineError::Io { step, source }) => {
                assert_eq!(step, "FailingSource");
                assert_eq!(source.kind(), std::io::ErrorKind::Other);
            }
            other => panic!("expected Io error, got {other:?}"),
        }
    }

    // ───── sleep_until_stop ─────

    /// When `stop` is already set, the helper must return effectively
    /// immediately, never sleeping the full duration. This is the teardown
    /// fast-path: the main thread sets `stop` then `join()`s, and the helper
    /// must not block for a poll interval afterward.
    #[test]
    fn sleep_until_stop_returns_immediately_when_already_stopped() {
        let stop = std::sync::atomic::AtomicBool::new(true);
        let start = std::time::Instant::now();
        sleep_until_stop(&stop, std::time::Duration::from_secs(10));
        assert!(
            start.elapsed() < std::time::Duration::from_millis(100),
            "expected near-immediate return, took {:?}",
            start.elapsed()
        );
    }

    /// When `stop` is set partway through the sleep, the helper must wake and
    /// return well before the full duration elapses (within a couple of poll
    /// slices), proving it interrupts a long sleep rather than waiting it out.
    #[test]
    fn sleep_until_stop_wakes_when_stopped_midway() {
        let stop = Arc::new(std::sync::atomic::AtomicBool::new(false));
        let stop_clone = Arc::clone(&stop);
        let setter = std::thread::spawn(move || {
            std::thread::sleep(std::time::Duration::from_millis(50));
            stop_clone.store(true, std::sync::atomic::Ordering::Relaxed);
        });
        let start = std::time::Instant::now();
        sleep_until_stop(&stop, std::time::Duration::from_secs(10));
        let elapsed = start.elapsed();
        setter.join().unwrap();
        // Woke shortly after the 50ms flag flip — far below the 10s budget.
        assert!(elapsed >= std::time::Duration::from_millis(40), "woke too early: {elapsed:?}");
        assert!(elapsed < std::time::Duration::from_millis(500), "woke too late: {elapsed:?}");
    }

    /// When `stop` is never set, the helper sleeps for approximately the full
    /// duration (it does not return early). Generous upper bound keeps it
    /// non-flaky on a loaded CI host.
    #[test]
    fn sleep_until_stop_sleeps_full_duration_when_never_stopped() {
        let stop = std::sync::atomic::AtomicBool::new(false);
        let start = std::time::Instant::now();
        sleep_until_stop(&stop, std::time::Duration::from_millis(100));
        let elapsed = start.elapsed();
        assert!(elapsed >= std::time::Duration::from_millis(95), "returned too early: {elapsed:?}");
        assert!(elapsed < std::time::Duration::from_secs(2), "ran far too long: {elapsed:?}");
    }
}
