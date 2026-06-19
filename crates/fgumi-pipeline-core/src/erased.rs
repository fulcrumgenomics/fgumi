//! Type erasure: `ErasedStep` trait + `TypedStep<S>` adapter.
//!
//! The runtime holds a heterogeneous chain in `Vec<Box<dyn ErasedStep>>`.
//! `TypedStep<S>` is the adapter that bridges between the type-erased
//! dispatch in the worker loop and the concrete `S::try_run` body.
//!
//! Each `ErasedStep` exposes the methods the runtime needs:
//!   - `clone_boxed` — make per-worker copies for `Parallel` steps
//!   - `build_output_set` — construct the producer's queue set + view from
//!     `StepProfile::output_queues` + `branch_ordering`
//!   - `build_input_handle` — pull the consumer's typed input handle out of
//!     the producer's `OutputQueueSet` (mutable: takes ownership of the
//!     branch's slot)
//!   - `wrap_outputs_view` — wrap the type-erased view into a typed
//!     `OutputHandles<S::Outputs>` for the worker to pass into `ctx.outputs`
//!   - `mark_outputs_drained` — close all output branches (called by the
//!     driver when a step returns `StepOutcome::Finished`, counter-gated for
//!     `Parallel` so only the last clone closes the shared output)
//!   - `is_source` — true iff `S::Input == ()` (used by chain-context
//!     construction to pick the source's unit-input path)

use std::any::{Any, TypeId};
use std::io;
use std::marker::PhantomData;
use std::sync::Arc;

use super::handles::{BranchInputHandle, OutputQueueSet};
use super::outputs::StepOutputs;
use super::queues::QueueSpec;
use super::reorder::BranchOrdering;
use super::signal::PipelineSignal;
use super::step::{
    Affinity, OutputHandles, OutputsViewAny, Step, StepCtx, StepKind, StepOutcome, StepProfile,
};

/// Type-erased step interface used by the worker loop.
pub trait ErasedStep: Send + 'static {
    fn profile(&self) -> StepProfile;

    /// The step's static name, returned WITHOUT building a `StepProfile`.
    ///
    /// `dispatch_one_step` reads the name on every dispatch (for stats /
    /// error reporting), so going through `profile()` there would heap-
    /// allocate the profile's two `Vec`s per dispatch (a virtual call, so
    /// the optimizer cannot elide them). Adapters cache the name at
    /// construction and return it here for free.
    fn name(&self) -> &'static str;

    /// Forward `Step::affinity` for worker-eligibility gating. Only
    /// consulted for `Serial` steps; the runtime calls this once during
    /// `build_worker_storage` to decide which workers get a `Shared`
    /// entry vs a `Skip` placeholder.
    fn affinity(&self) -> Affinity;

    /// Dispatch `S::try_run` after downcasting queue handles.
    ///
    /// # Errors
    ///
    /// Forwards any I/O error from the step body.
    fn try_run_erased(&mut self, ctx: &mut ErasedStepCtx<'_>) -> io::Result<StepOutcome>;

    /// Construct a fresh per-worker copy of this step. Used for Parallel
    /// steps. Cheap — implementing types call `S::clone()`, where typical
    /// state is unit-struct or `Arc<dyn Fn>` (one atomic increment).
    fn clone_boxed(&self) -> Box<dyn ErasedStep>;

    /// Take ownership of the consumer's input handle from the producer's
    /// output queue set. Used by the framework when constructing chain
    /// topology; each branch is taken exactly once.
    fn build_input_handle(
        &self,
        producer_set: &mut OutputQueueSet,
        branch_idx: usize,
    ) -> Box<dyn Any + Send + Sync>;

    /// Input arity. Default `1` for single-input
    /// [`crate::step::Step`] impls; multi-input
    /// adapters (`TypedStep2`, future `StepN`) override to return their
    /// arity. Used by [`crate::runtime::contexts`]
    /// to decide which input-construction path to take.
    fn input_arity(&self) -> usize {
        1
    }

    /// Take ownership of TWO consumer input handles from two
    /// upstream output queue sets, paired into a typed
    /// [`crate::handles::TwoInputHandles<A, B>`]
    /// per the consumer's [`crate::step::Step2`]
    /// associated types.
    ///
    /// Default impl panics — single-input steps never have arity 2.
    /// `TypedStep2<S>` overrides to call `take_typed_input` twice
    /// (once per input slot) and wrap the pair.
    ///
    /// # Panics
    ///
    /// Default impl always panics; multi-input adapters override.
    fn build_two_input_handles(
        &self,
        _producer_sets: &mut [OutputQueueSet],
        _p0_idx: usize,
        _p0_branch: usize,
        _p1_idx: usize,
        _p1_branch: usize,
    ) -> Box<dyn Any + Send + Sync> {
        let p = self.profile();
        panic!(
            "build_two_input_handles called on '{}' (kind = {:?}); \
             only Step2 adapters (`TypedStep2`) support arity-2 input \
             construction. This is a framework bug — chain-build code \
             should only dispatch arity-2 input construction to steps \
             that override `input_arity` to return 2.",
            p.name, p.kind
        );
    }

    /// Build this step's output queue set + outputs view from the profile's
    /// per-branch queue specs and ordering directives.
    fn build_output_set(&self) -> (OutputQueueSet, OutputsViewAny);

    /// Like [`Self::build_output_set`] but forces every output branch to a
    /// **direct, unbounded** transport (no byte bound, no reorder stage).
    ///
    /// Used only by the single-thread *fused* driver
    /// ([`crate::runtime::run_fused_single_thread`]). At one
    /// worker, FIFO push order is already the correct order, so the reorder
    /// stage is dead weight, and the unbounded buffer never needs backpressure
    /// because the consumer is driven in the same pass over the chain. Each
    /// sub-step emits a bounded number of items per `try_run`, so the buffer
    /// stays shallow (the fused-chain bounded-fan-out invariant).
    fn build_fused_output_set(&self) -> (OutputQueueSet, OutputsViewAny);

    /// Wrap a typed `OutputsViewAny` into `Box<OutputHandles<S::Outputs>>`
    /// (type-erased as `Any`). The runtime stores this box per step so the
    /// worker loop's `try_run_erased` can downcast it to the typed
    /// `OutputHandles<S::Outputs>` the step expects in its `ctx.outputs`.
    fn wrap_outputs_view(&self, view: OutputsViewAny) -> Box<dyn Any + Send + Sync>;

    /// Mark all of this step's output branches drained. Downcasts the
    /// type-erased outputs handle to `OutputHandles<S::Outputs>` and calls
    /// the typed `mark_all_drained`. Called by the driver when a step returns
    /// `StepOutcome::Finished` (counter-gated for `Parallel`).
    fn mark_outputs_drained(&self, outputs: &(dyn Any + Send + Sync));

    /// Returns `true` iff `S::Input = ()` — i.e., this is a source step.
    /// Used by chain-context construction to pick the source's unit-input
    /// path (a source has no upstream queue to take an input handle from).
    fn is_source(&self) -> bool;
}

/// Context the worker loop hands to `ErasedStep` methods.
pub struct ErasedStepCtx<'a> {
    /// Boxed `BranchInputHandle<S::Input>`. Adapter downcasts.
    pub input: &'a (dyn Any + Send + Sync),
    /// `OutputHandles<S::Outputs>`. Adapter downcasts.
    pub outputs: &'a (dyn Any + Send + Sync),
    /// Shared signal (error/cancel). Workers consult; steps don't directly.
    pub signal: &'a Arc<PipelineSignal>,
}

/// Adapter that wraps a concrete `Step` impl as an `ErasedStep`.
///
/// ## Cached typed handles
///
/// `try_run_erased` needs typed
/// `&BranchInputHandle<S::Input>` / `&OutputHandles<S::Outputs>` views
/// of the type-erased context boxes. The first call resolves them via
/// `Any::downcast_ref` (a `TypeId` compare + transmute); subsequent calls
/// reuse the cached references through the `cached_input` /
/// `cached_outputs` fields. Without the cache, a 4-thread CODEC 8M
/// run pays ≈230 samples (~1.6%) on `downcast_ref` `TypeId` compares
/// across the dispatch hot path; the cache eliminates them entirely.
///
/// The cache is sound because:
///
///   1. The boxes are owned by `ChainContexts` (an `Arc` held alive
///      for the entire `Pipeline::run` call by every worker thread).
///   2. Every `TypedStep<S>` instance — owned (`Parallel`), shared
///      (`Serial`, behind `Mutex`), or pinned (`Exclusive`) — is
///      destroyed before `ChainContexts` goes out of scope (workers
///      exit, then the runtime drops the contexts).
///   3. Every dispatch passes the **same** box reference for a given
///      `step_idx` (see `run_worker_loop`'s
///      `contexts.inputs[step_idx.0].as_ref()`). The cached pointer
///      always refers to that same box.
///
/// `clone_boxed` calls `TypedStep::new(...)` which initializes the
/// cache to `None`, so per-worker clones (`Parallel` steps) start
/// with a fresh cache.
pub struct TypedStep<S: Step> {
    inner: S,
    /// Static step name, cached from `inner.profile().name` at
    /// construction so `ErasedStep::name()` (read per dispatch) never
    /// rebuilds the profile's `Vec`s. See the `ErasedStep::name` doc.
    name: &'static str,
    /// Cached downcast of `ctx.input` after the first dispatch.
    /// `'static` is a lie — the actual lifetime is bounded by
    /// `ChainContexts` — enforced via the `// SAFETY:` comment below.
    cached_input: Option<&'static BranchInputHandle<S::Input>>,
    /// Cached downcast of `ctx.outputs`. Same lifetime story.
    cached_outputs: Option<&'static OutputHandles<S::Outputs>>,
    _phantom: PhantomData<fn() -> S>,
}

impl<S: Step> TypedStep<S> {
    pub fn new(step: S) -> Self {
        let name = step.profile().name;
        Self { inner: step, name, cached_input: None, cached_outputs: None, _phantom: PhantomData }
    }

    /// Resolve the typed input handle from the dispatch context, caching
    /// the result for subsequent dispatches.
    ///
    /// # Panics
    ///
    /// Panics on the first dispatch if the `ctx.input` box doesn't
    /// downcast to `BranchInputHandle<S::Input>` (a chain topology
    /// invariant violation; the builder should have caught this).
    #[allow(unsafe_code)]
    fn resolve_input<'a>(&mut self, ctx: &ErasedStepCtx<'a>) -> &'a BranchInputHandle<S::Input> {
        if let Some(cached) = self.cached_input {
            // SAFETY: Lifetime extension from `'static` (cache slot) back
            // to `'a` (the dispatch context's lifetime). The cached
            // pointer was originally a `&'a BranchInputHandle<S::Input>`
            // pulled out of `ChainContexts.inputs[step_idx]`; that box
            // outlives every `TypedStep<S>` (point 2 in the type-level
            // doc). Dispatches always pass the same box for the same
            // `step_idx` (point 3), so the pointer is still pointing
            // at the live box on every subsequent dispatch.
            return unsafe {
                std::mem::transmute::<&BranchInputHandle<S::Input>, &'a BranchInputHandle<S::Input>>(
                    cached,
                )
            };
        }
        let r: &'a BranchInputHandle<S::Input> = ctx
            .input
            .downcast_ref::<BranchInputHandle<S::Input>>()
            .expect("input handle downcast failed — chain topology invariant");
        // SAFETY: Lifetime extension from `'a` to `'static` for storage.
        // The same point-2/point-3 invariants from `// SAFETY:` above
        // apply: the box outlives `self`, and we only ever read the
        // cache through `resolve_input`, which immediately re-extends
        // back to a bounded `'a` before handing it to user code.
        let cached: &'static BranchInputHandle<S::Input> = unsafe {
            std::mem::transmute::<
                &'a BranchInputHandle<S::Input>,
                &'static BranchInputHandle<S::Input>,
            >(r)
        };
        self.cached_input = Some(cached);
        r
    }

    /// Resolve the typed outputs handle from the dispatch context,
    /// caching the result for subsequent dispatches. See `resolve_input`
    /// for the safety argument.
    #[allow(unsafe_code)]
    fn resolve_outputs<'a>(&mut self, ctx: &ErasedStepCtx<'a>) -> &'a OutputHandles<S::Outputs> {
        if let Some(cached) = self.cached_outputs {
            // SAFETY: see `resolve_input`; same boxes/lifetimes story.
            return unsafe {
                std::mem::transmute::<&OutputHandles<S::Outputs>, &'a OutputHandles<S::Outputs>>(
                    cached,
                )
            };
        }
        let r: &'a OutputHandles<S::Outputs> = ctx
            .outputs
            .downcast_ref::<OutputHandles<S::Outputs>>()
            .expect("outputs handle downcast failed — chain topology invariant");
        // SAFETY: see `resolve_input`; same boxes/lifetimes story.
        let cached: &'static OutputHandles<S::Outputs> = unsafe {
            std::mem::transmute::<&'a OutputHandles<S::Outputs>, &'static OutputHandles<S::Outputs>>(
                r,
            )
        };
        self.cached_outputs = Some(cached);
        r
    }
}

impl<S> ErasedStep for TypedStep<S>
where
    S: Step,
{
    fn profile(&self) -> StepProfile {
        self.inner.profile()
    }

    fn name(&self) -> &'static str {
        self.name
    }

    fn affinity(&self) -> Affinity {
        self.inner.affinity()
    }

    fn try_run_erased(&mut self, ctx: &mut ErasedStepCtx<'_>) -> io::Result<StepOutcome> {
        let input = self.resolve_input(ctx);
        let outputs = self.resolve_outputs(ctx);
        let mut step_ctx = StepCtx { input, outputs };
        self.inner.try_run(&mut step_ctx)
    }

    fn clone_boxed(&self) -> Box<dyn ErasedStep> {
        // `TypedStep::new` initializes `cached_input` / `cached_outputs`
        // to `None` so the new clone resolves them fresh on its first
        // dispatch (the box references in *this* `TypedStep<S>`'s cache
        // are still valid for the new clone too — they refer to the
        // same `ChainContexts` boxes — but resolving fresh is simpler
        // and keeps the cache lifetime story local to each clone).
        Box::new(TypedStep::new(self.inner.new_worker_copy()))
    }

    fn build_input_handle(
        &self,
        producer_set: &mut OutputQueueSet,
        branch_idx: usize,
    ) -> Box<dyn Any + Send + Sync> {
        let handle: BranchInputHandle<S::Input> =
            producer_set.take_typed_input::<S::Input>(branch_idx);
        Box::new(handle)
    }

    fn build_output_set(&self) -> (OutputQueueSet, OutputsViewAny) {
        let profile = self.inner.profile();
        // Producers with `StepKind::Serial` or `StepKind::Exclusive` emit
        // items in arrival order by construction (the framework's mutex
        // serializes pushes; an Exclusive-owned step has a single
        // dispatcher). Inserting a `ReorderStage` on those output edges
        // is pure overhead — every push pays an ordinal allocation +
        // `Sequenced<T>` wrap + reorder mutex hop for items that are
        // already in order. The framework collapses any
        // `BranchOrdering::ByOrdinal` / `ByItemOrdinal` declaration to
        // `None` here, so `build_queues` constructs the direct transport
        // path with no reorder stage. Step authors keep declaring
        // `ByItemOrdinal` (preserves the intent in the profile;
        // documents what the consumer needs); the framework is just
        // smart enough to skip the wrap when the producer can already
        // satisfy that need.
        let effective_orderings: Vec<BranchOrdering> = match profile.kind {
            StepKind::Parallel => profile.branch_ordering.clone(),
            StepKind::Serial | StepKind::Exclusive => {
                profile.branch_ordering.iter().map(|_| BranchOrdering::None).collect()
            }
        };
        <S::Outputs as StepOutputs>::build_queues(&profile.output_queues, &effective_orderings)
    }

    fn build_fused_output_set(&self) -> (OutputQueueSet, OutputsViewAny) {
        // Force a direct, unbounded transport on every branch: at one worker
        // FIFO is already correct (no reorder), and the consumer is driven in
        // the same pass (no backpressure). Arity comes from the profile.
        let n = self.inner.profile().output_queues.len();
        let specs = vec![QueueSpec::Unbounded; n];
        let orderings = vec![BranchOrdering::None; n];
        <S::Outputs as StepOutputs>::build_queues(&specs, &orderings)
    }

    fn wrap_outputs_view(&self, view: OutputsViewAny) -> Box<dyn Any + Send + Sync> {
        let outputs: OutputHandles<S::Outputs> = OutputHandles::new(view);
        Box::new(outputs)
    }

    fn mark_outputs_drained(&self, outputs: &(dyn Any + Send + Sync)) {
        let typed = outputs
            .downcast_ref::<OutputHandles<S::Outputs>>()
            .expect("outputs handle downcast failed in mark_outputs_drained");
        <S::Outputs as StepOutputs>::mark_all_drained(typed);
    }

    fn is_source(&self) -> bool {
        TypeId::of::<S::Input>() == TypeId::of::<()>()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TypedStep2<S> — adapter for `Step2` impls (two-input merge steps).
//
// Same shape as TypedStep<S> for outputs: one `OutputHandles<S::Outputs>`
// cached after the first dispatch. Inputs differ — instead of a single
// `BranchInputHandle<S::Input>`, the adapter caches a single
// `&TwoInputHandles<S::InputA, S::InputB>` and lends per-branch refs
// (`ctx.a`, `ctx.b`) into the wrapped `Step2::try_run`.
//
// Drain detection: a `Step2` consumer is "drained" when **both** input
// branches report drained. The step itself checks this in `try_run`
// (`ctx.a.is_drained() && ctx.b.is_drained()`) to decide when to report
// `Finished`; the typed-erased input box is a `TwoInputHandles<A, B>`.
// ─────────────────────────────────────────────────────────────────────────────

use super::handles::TwoInputHandles;
use super::step::{Step2, StepCtx2};

/// Adapter wrapping a [`Step2`] impl as an [`ErasedStep`].
///
/// Cached typed handles follow the same pattern as [`TypedStep`]:
/// the first dispatch resolves `&TwoInputHandles<S::InputA, S::InputB>`
/// / `&OutputHandles<S::Outputs>` via `Any::downcast_ref`, subsequent
/// dispatches reuse the cached references through unsafe lifetime
/// extension. The safety argument is identical (the boxes are owned
/// by `ChainContexts` which outlives every adapter instance, and
/// every dispatch passes the same box reference for a given
/// `step_idx`).
pub struct TypedStep2<S: Step2> {
    inner: S,
    /// Static step name, cached at construction. See `ErasedStep::name`.
    name: &'static str,
    cached_inputs: Option<&'static TwoInputHandles<S::InputA, S::InputB>>,
    cached_outputs: Option<&'static OutputHandles<S::Outputs>>,
    _phantom: PhantomData<fn() -> S>,
}

impl<S: Step2> TypedStep2<S> {
    pub fn new(step: S) -> Self {
        let name = step.profile().name;
        Self { inner: step, name, cached_inputs: None, cached_outputs: None, _phantom: PhantomData }
    }

    #[allow(unsafe_code)]
    fn resolve_inputs<'a>(
        &mut self,
        ctx: &ErasedStepCtx<'a>,
    ) -> &'a TwoInputHandles<S::InputA, S::InputB> {
        if let Some(cached) = self.cached_inputs {
            // SAFETY: lifetime extension from `'static` (cache slot)
            // back to `'a` (the dispatch context's lifetime). The
            // cached pointer was originally a
            // `&'a TwoInputHandles<S::InputA, S::InputB>` pulled out
            // of `ChainContexts.inputs[step_idx]`; that box outlives
            // every `TypedStep2<S>` instance (same point-2/point-3
            // invariants as `TypedStep::resolve_input`).
            return unsafe {
                std::mem::transmute::<
                    &TwoInputHandles<S::InputA, S::InputB>,
                    &'a TwoInputHandles<S::InputA, S::InputB>,
                >(cached)
            };
        }
        let r: &'a TwoInputHandles<S::InputA, S::InputB> = ctx
            .input
            .downcast_ref::<TwoInputHandles<S::InputA, S::InputB>>()
            .expect("input handle downcast failed — Step2 chain topology invariant");
        // SAFETY: lifetime extension from `'a` to `'static` for
        // storage. The box outlives `self`; we only ever read the
        // cache through `resolve_inputs`, which immediately re-extends
        // back to a bounded `'a` before handing it to user code.
        let cached: &'static TwoInputHandles<S::InputA, S::InputB> = unsafe {
            std::mem::transmute::<
                &'a TwoInputHandles<S::InputA, S::InputB>,
                &'static TwoInputHandles<S::InputA, S::InputB>,
            >(r)
        };
        self.cached_inputs = Some(cached);
        r
    }

    #[allow(unsafe_code)]
    fn resolve_outputs<'a>(&mut self, ctx: &ErasedStepCtx<'a>) -> &'a OutputHandles<S::Outputs> {
        if let Some(cached) = self.cached_outputs {
            return unsafe {
                std::mem::transmute::<&OutputHandles<S::Outputs>, &'a OutputHandles<S::Outputs>>(
                    cached,
                )
            };
        }
        let r: &'a OutputHandles<S::Outputs> = ctx
            .outputs
            .downcast_ref::<OutputHandles<S::Outputs>>()
            .expect("outputs handle downcast failed — Step2 chain topology invariant");
        let cached: &'static OutputHandles<S::Outputs> = unsafe {
            std::mem::transmute::<&'a OutputHandles<S::Outputs>, &'static OutputHandles<S::Outputs>>(
                r,
            )
        };
        self.cached_outputs = Some(cached);
        r
    }
}

impl<S: Step2> ErasedStep for TypedStep2<S> {
    fn profile(&self) -> StepProfile {
        self.inner.profile()
    }

    fn name(&self) -> &'static str {
        self.name
    }

    fn affinity(&self) -> Affinity {
        self.inner.affinity()
    }

    fn try_run_erased(&mut self, ctx: &mut ErasedStepCtx<'_>) -> io::Result<StepOutcome> {
        let inputs = self.resolve_inputs(ctx);
        let outputs = self.resolve_outputs(ctx);
        let mut typed_ctx = StepCtx2::<S> { a: &inputs.a, b: &inputs.b, outputs };
        self.inner.try_run(&mut typed_ctx)
    }

    fn clone_boxed(&self) -> Box<dyn ErasedStep> {
        Box::new(TypedStep2::new(self.inner.new_worker_copy()))
    }

    fn build_input_handle(
        &self,
        _producer_set: &mut OutputQueueSet,
        _branch_idx: usize,
    ) -> Box<dyn Any + Send + Sync> {
        let p = self.profile();
        panic!(
            "build_input_handle called on Step2 adapter '{}' (kind = {:?}); \
             multi-input steps build their inputs via build_two_input_handles. \
             This is a framework bug — chain-build code should dispatch on \
             input_arity().",
            p.name, p.kind
        );
    }

    fn input_arity(&self) -> usize {
        2
    }

    fn build_two_input_handles(
        &self,
        producer_sets: &mut [OutputQueueSet],
        p0_idx: usize,
        p0_branch: usize,
        p1_idx: usize,
        p1_branch: usize,
    ) -> Box<dyn Any + Send + Sync> {
        assert_ne!(
            p0_idx,
            p1_idx,
            "Step2 producers must be distinct steps (each upstream is its \
             own subchain). Got p0_idx == p1_idx == {p0_idx} for step '{}'.",
            self.profile().name
        );
        // Borrow two disjoint elements of `producer_sets` simultaneously.
        // `split_at_mut(lo+1)` puts producer_sets[lo] in the first half;
        // we index into the second half for the hi side.
        let (lo_idx, hi_idx, swap) =
            if p0_idx < p1_idx { (p0_idx, p1_idx, false) } else { (p1_idx, p0_idx, true) };
        let (lo_half, hi_half) = producer_sets.split_at_mut(lo_idx + 1);
        let lo_set: &mut OutputQueueSet = &mut lo_half[lo_idx];
        let hi_set: &mut OutputQueueSet = &mut hi_half[hi_idx - (lo_idx + 1)];
        let (a_set, b_set) = if swap {
            (hi_set, lo_set) // p0 = high, p1 = low
        } else {
            (lo_set, hi_set) // p0 = low,  p1 = high
        };
        let a: BranchInputHandle<S::InputA> = a_set.take_typed_input::<S::InputA>(p0_branch);
        let b: BranchInputHandle<S::InputB> = b_set.take_typed_input::<S::InputB>(p1_branch);
        Box::new(TwoInputHandles::<S::InputA, S::InputB>::new(a, b))
    }

    fn build_output_set(&self) -> (OutputQueueSet, OutputsViewAny) {
        <S::Outputs as StepOutputs>::build_queues(
            &self.inner.profile().output_queues,
            &self.inner.profile().branch_ordering,
        )
    }

    fn build_fused_output_set(&self) -> (OutputQueueSet, OutputsViewAny) {
        // `Step2` never appears in a linear (single-input) fused chain, so this
        // is unreachable in practice; provided for trait completeness with the
        // same direct+unbounded semantics as the single-input adapter.
        let n = self.inner.profile().output_queues.len();
        let specs = vec![QueueSpec::Unbounded; n];
        let orderings = vec![BranchOrdering::None; n];
        <S::Outputs as StepOutputs>::build_queues(&specs, &orderings)
    }

    fn wrap_outputs_view(&self, view: OutputsViewAny) -> Box<dyn Any + Send + Sync> {
        let typed: OutputHandles<S::Outputs> = OutputHandles::<S::Outputs>::new(view);
        Box::new(typed)
    }

    fn mark_outputs_drained(&self, outputs: &(dyn Any + Send + Sync)) {
        let typed = outputs
            .downcast_ref::<OutputHandles<S::Outputs>>()
            .expect("mark_outputs_drained downcast failed — Step2 chain topology invariant");
        <S::Outputs as StepOutputs>::mark_all_drained(typed);
    }

    fn is_source(&self) -> bool {
        // Step2 consumers are never sources by definition (they have
        // two non-unit input branches).
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::handles::BranchInputHandle;
    use crate::outputs::Single;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;
    use crate::step::{
        InputHandle, OutputHandles, Step, StepCtx, StepKind, StepOutcome, StepProfile,
    };

    /// Trivial step: u32 → u32+1 (single output).
    #[derive(Clone)]
    struct AddOne;

    impl Step for AddOne {
        type Input = u32;
        type Outputs = Single<u32>;

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "AddOne",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(n) => {
                    let _ = ctx.outputs.push(n + 1);
                    Ok(StepOutcome::Progress)
                }
                None => Ok(StepOutcome::NoProgress),
            }
        }

        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    /// Build a `(producer_set, producer_outputs)` pair so we can simulate a
    /// chain link by manually constructing the upstream side.
    fn build_addone_outputs() -> (OutputQueueSet, OutputHandles<Single<u32>>) {
        let producer: Box<dyn ErasedStep> = Box::new(TypedStep::new(AddOne));
        let (queue_set, outputs_view) = producer.build_output_set();
        let outputs: OutputHandles<Single<u32>> = OutputHandles::new(outputs_view);
        (queue_set, outputs)
    }

    #[test]
    fn typed_step_round_trips_through_erased_dispatch() {
        // Producer's output set carries the input handle that the consumer
        // (also an AddOne) will pull from. The consumer's own output set is
        // separate.
        let (mut producer_set, producer_outputs) = build_addone_outputs();

        // Push a u32 onto the producer's output (which is the consumer's input).
        producer_outputs.push(41).unwrap();

        // Consumer takes the typed input handle out of the producer's set.
        let mut consumer: Box<dyn ErasedStep> = Box::new(TypedStep::new(AddOne));
        let input_any = consumer.build_input_handle(&mut producer_set, 0);

        // Consumer needs its own output set + view to run.
        let (consumer_set, consumer_view) = consumer.build_output_set();
        let consumer_outputs: OutputHandles<Single<u32>> = OutputHandles::new(consumer_view);

        let signal = PipelineSignal::new();
        let mut ctx = ErasedStepCtx {
            input: input_any.as_ref(),
            outputs: &consumer_outputs as &(dyn Any + Send + Sync),
            signal: &signal,
        };
        let outcome = consumer.try_run_erased(&mut ctx).unwrap();
        assert_eq!(outcome, StepOutcome::Progress);

        // Consumer pushed (41 + 1) = 42 onto its own output. Pull it.
        let mut consumer_set = consumer_set;
        let consumer_input = consumer_set.take_typed_input::<u32>(0);
        assert_eq!(consumer_input.pop(), Some(42));
    }

    #[test]
    fn typed_step_returns_noprogress_on_empty_input() {
        let (mut producer_set, _producer_outputs) = build_addone_outputs();
        let mut consumer: Box<dyn ErasedStep> = Box::new(TypedStep::new(AddOne));
        let input_any = consumer.build_input_handle(&mut producer_set, 0);

        let (_consumer_set, consumer_view) = consumer.build_output_set();
        let consumer_outputs: OutputHandles<Single<u32>> = OutputHandles::new(consumer_view);

        let signal = PipelineSignal::new();
        let mut ctx = ErasedStepCtx {
            input: input_any.as_ref(),
            outputs: &consumer_outputs as &(dyn Any + Send + Sync),
            signal: &signal,
        };
        let outcome = consumer.try_run_erased(&mut ctx).unwrap();
        assert_eq!(outcome, StepOutcome::NoProgress);
    }

    #[test]
    fn clone_boxed_yields_independent_step() {
        let original: Box<dyn ErasedStep> = Box::new(TypedStep::new(AddOne));
        let cloned = original.clone_boxed();
        assert_eq!(original.profile().name, "AddOne");
        assert_eq!(cloned.profile().name, "AddOne");
    }

    /// `ErasedStep::name()` (the per-dispatch, allocation-free accessor)
    /// must agree with `profile().name`, including across `clone_boxed`
    /// (the per-worker clone path that re-caches the name in `new`).
    #[test]
    fn erased_name_matches_profile_name() {
        let step: Box<dyn ErasedStep> = Box::new(TypedStep::new(AddOne));
        assert_eq!(step.name(), step.profile().name);
        assert_eq!(step.name(), "AddOne");
        let cloned = step.clone_boxed();
        assert_eq!(cloned.name(), cloned.profile().name);
        assert_eq!(cloned.name(), "AddOne");
    }

    #[test]
    fn build_output_set_uses_profile_queues_and_ordering() {
        let typed: Box<dyn ErasedStep> = Box::new(TypedStep::new(AddOne));
        let (mut queue_set, _outputs_view) = typed.build_output_set();
        assert_eq!(queue_set.n_branches(), 1);
        // Verify the queue is u32-typed by taking the input handle.
        let input = queue_set.take_typed_input::<u32>(0);
        // Empty initially.
        assert_eq!(input.pop(), None);
    }

    #[test]
    fn build_input_handle_downcasts_to_typed_branch_handle() {
        // Producer pushes a value onto its output.
        let (mut producer_set, producer_outputs) = build_addone_outputs();
        producer_outputs.push(7).unwrap();

        // Consumer grabs the typed input handle.
        let consumer: Box<dyn ErasedStep> = Box::new(TypedStep::new(AddOne));
        let input_any = consumer.build_input_handle(&mut producer_set, 0);
        let input =
            input_any.downcast_ref::<BranchInputHandle<u32>>().expect("input handle downcast");
        assert_eq!(input.pop(), Some(7));
    }

    #[test]
    fn wrap_outputs_view_yields_typed_outputs_handles() {
        let typed: Box<dyn ErasedStep> = Box::new(TypedStep::new(AddOne));
        let (mut queue_set, view) = typed.build_output_set();
        let outputs_any = typed.wrap_outputs_view(view);

        // Downcast back to OutputHandles<Single<u32>> and push.
        let outputs = outputs_any
            .downcast_ref::<OutputHandles<Single<u32>>>()
            .expect("OutputHandles downcast");
        outputs.push(99).unwrap();

        let input = queue_set.take_typed_input::<u32>(0);
        assert_eq!(input.pop(), Some(99));
    }

    #[test]
    fn mark_outputs_drained_propagates_to_consumer_input() {
        // Producer is an AddOne; its outputs feed a consumer's input.
        let producer: Box<dyn ErasedStep> = Box::new(TypedStep::new(AddOne));
        let (mut producer_set, producer_view) = producer.build_output_set();
        let producer_outputs_any = producer.wrap_outputs_view(producer_view);

        let consumer: Box<dyn ErasedStep> = Box::new(TypedStep::new(AddOne));
        let input_any = consumer.build_input_handle(&mut producer_set, 0);
        let consumer_input = input_any
            .downcast_ref::<BranchInputHandle<u32>>()
            .expect("consumer input handle is BranchInputHandle<u32>");

        // Initially not drained.
        assert!(!InputHandle::is_drained(consumer_input));

        // Producer marks its outputs drained via the type-erased path.
        producer.mark_outputs_drained(producer_outputs_any.as_ref());

        // Drained signal reaches consumer.
        assert!(InputHandle::is_drained(consumer_input));
    }

    /// Multi-output step: u32 → (u32, u32) split into low + high bytes.
    #[derive(Clone)]
    struct SplitBytes;

    impl Step for SplitBytes {
        type Input = u32;
        type Outputs = (u32, u32);

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "SplitBytes",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![
                    QueueSpec::CountBounded { capacity: 2 },
                    QueueSpec::CountBounded { capacity: 2 },
                ],
                branch_ordering: vec![BranchOrdering::None, BranchOrdering::None],
            }
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(n) => {
                    let v = ctx.outputs.view();
                    let _ = v.a.push(n & 0xFF);
                    let _ = v.b.push((n >> 8) & 0xFF);
                    Ok(StepOutcome::Progress)
                }
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    #[test]
    fn split_bytes_routes_to_two_branches() {
        // Producer feeding `SplitBytes`: another AddOne (Single<u32> → Single<u32>).
        // We just need a way to feed u32 into SplitBytes, so reuse AddOne's outputs.
        let (mut producer_set, producer_outputs) = build_addone_outputs();
        producer_outputs.push(0xABCD).unwrap();

        // SplitBytes is the consumer; takes the typed input from the producer.
        let mut splitter: Box<dyn ErasedStep> = Box::new(TypedStep::new(SplitBytes));
        let input_any = splitter.build_input_handle(&mut producer_set, 0);

        // Splitter's own output set has two branches.
        let (mut splitter_set, splitter_view) = splitter.build_output_set();
        let splitter_outputs: OutputHandles<(u32, u32)> = OutputHandles::new(splitter_view);

        let signal = PipelineSignal::new();
        let mut ctx = ErasedStepCtx {
            input: input_any.as_ref(),
            outputs: &splitter_outputs as &(dyn Any + Send + Sync),
            signal: &signal,
        };

        let outcome = splitter.try_run_erased(&mut ctx).unwrap();
        assert_eq!(outcome, StepOutcome::Progress);

        let in_a = splitter_set.take_typed_input::<u32>(0);
        let in_b = splitter_set.take_typed_input::<u32>(1);
        assert_eq!(in_a.pop(), Some(0xCD));
        assert_eq!(in_b.pop(), Some(0xAB));
    }

    #[test]
    fn is_source_true_for_unit_input_steps() {
        #[derive(Clone)]
        struct UnitInputStep;
        impl Step for UnitInputStep {
            type Input = ();
            type Outputs = Single<u32>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "UnitInput",
                    kind: StepKind::Exclusive,
                    sticky: false,
                    output_queues: vec![QueueSpec::CountBounded { capacity: 1 }],
                    branch_ordering: vec![BranchOrdering::ByOrdinal],
                }
            }
            fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                Ok(StepOutcome::Finished)
            }
        }

        let source: Box<dyn ErasedStep> = Box::new(TypedStep::new(UnitInputStep));
        let mid: Box<dyn ErasedStep> = Box::new(TypedStep::new(AddOne));
        assert!(source.is_source());
        assert!(!mid.is_source());
    }
}
