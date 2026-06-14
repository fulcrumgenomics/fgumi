//! Closure-driven mid-steps: `Process`, `ProcessOrdered`,
//! `ProcessWithWorkerState`, `MiAssign`.
//!
//! The four are colocated per Phase 0 design line 371:
//! `process.rs # process(fn), process_with_worker_state(init, fn), mi_assign(fn)`.
//!
//! Phase 3 ships single-output Process only. Multi-output (e.g.,
//! `correct`'s `CorrectOutputs`) is Phase 4 work.

use std::io;
use std::marker::PhantomData;
use std::sync::Arc;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::outputs::{OrderedBytesSingle, Single};
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

// ─────────────────────────────────────────────────────────────────────────────
// Process<I, T, F> — Parallel, FIFO output (Single<T>).
// ─────────────────────────────────────────────────────────────────────────────

/// Closure-driven `Parallel` mid-step with `Single<T>` output (FIFO by
/// arrival). Use [`ProcessOrdered`] when consumers need order preserved.
pub struct Process<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + 'static,
    F: Fn(I) -> io::Result<T> + Send + Sync + 'static,
{
    f: Arc<F>,
    held: HeldSlot<Unpushed<T>>,
    name: &'static str,
    capacity: usize,
    _phantom: PhantomData<fn() -> I>,
}

impl<I, T, F> Process<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + 'static,
    F: Fn(I) -> io::Result<T> + Send + Sync + 'static,
{
    pub fn new(name: &'static str, capacity: usize, f: F) -> Self {
        Self { f: Arc::new(f), held: HeldSlot::new(), name, capacity, _phantom: PhantomData }
    }
}

impl<I, T, F> Clone for Process<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + 'static,
    F: Fn(I) -> io::Result<T> + Send + Sync + 'static,
{
    fn clone(&self) -> Self {
        Self {
            f: Arc::clone(&self.f),
            held: HeldSlot::new(),
            name: self.name,
            capacity: self.capacity,
            _phantom: PhantomData,
        }
    }
}

impl<I, T, F> Step for Process<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + 'static,
    F: Fn(I) -> io::Result<T> + Send + Sync + 'static,
{
    type Input = I;
    type Outputs = Single<T>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::CountBounded { capacity: self.capacity }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }
        let Some(item) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        let out = (self.f)(item)?;
        match ctx.outputs.push(out) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

/// Convenience constructor for [`Process`].
pub fn process<I, T, F>(name: &'static str, capacity: usize, f: F) -> Process<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + 'static,
    F: Fn(I) -> io::Result<T> + Send + Sync + 'static,
{
    Process::new(name, capacity, f)
}

// ─────────────────────────────────────────────────────────────────────────────
// ProcessOrdered<I, T, F> — Parallel, ByItemOrdinal output
// (OrderedBytesSingle<T>). The canonical BAM-step shape.
// ─────────────────────────────────────────────────────────────────────────────

/// Closure-driven `Parallel` mid-step with `OrderedBytesSingle<T>` output.
/// Closure must produce `T: HeapSize + Ordered`; the framework's reorder
/// stage uses `T::ordinal()` to preserve global ordering.
///
/// Phase 3's BAM steps that don't fit into a more specific type
/// (`GroupBam`, `MiAssign`) compose via `ProcessOrdered` with a closure
/// that copies the input batch's serial onto the output.
pub struct ProcessOrdered<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(I) -> io::Result<T> + Send + Sync + 'static,
{
    f: Arc<F>,
    held: HeldSlot<Unpushed<T>>,
    name: &'static str,
    limit_bytes: u64,
    _phantom: PhantomData<fn() -> I>,
}

impl<I, T, F> ProcessOrdered<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(I) -> io::Result<T> + Send + Sync + 'static,
{
    pub fn new(name: &'static str, limit_bytes: u64, f: F) -> Self {
        Self { f: Arc::new(f), held: HeldSlot::new(), name, limit_bytes, _phantom: PhantomData }
    }
}

impl<I, T, F> Clone for ProcessOrdered<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(I) -> io::Result<T> + Send + Sync + 'static,
{
    fn clone(&self) -> Self {
        Self {
            f: Arc::clone(&self.f),
            held: HeldSlot::new(),
            name: self.name,
            limit_bytes: self.limit_bytes,
            _phantom: PhantomData,
        }
    }
}

impl<I, T, F> Step for ProcessOrdered<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(I) -> io::Result<T> + Send + Sync + 'static,
{
    type Input = I;
    type Outputs = OrderedBytesSingle<T>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.limit_bytes }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }
        let Some(item) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        let out = (self.f)(item)?;
        match ctx.outputs.push(out) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

/// Convenience constructor for [`ProcessOrdered`].
pub fn process_ordered<I, T, F>(
    name: &'static str,
    limit_bytes: u64,
    f: F,
) -> ProcessOrdered<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(I) -> io::Result<T> + Send + Sync + 'static,
{
    ProcessOrdered::new(name, limit_bytes, f)
}

// ─────────────────────────────────────────────────────────────────────────────
// ProcessWithWorkerState<I, T, F, S, FInit> — Parallel + lazy per-worker
// state. Used by consensus callers (simplex/duplex/codec) where each
// worker reuses one caller across all batches.
// ─────────────────────────────────────────────────────────────────────────────

/// `Parallel` + lazy per-worker state. The framework's `Clone` resets
/// `state` to `None`; first `try_run` initializes via `init()`, subsequent
/// calls reuse the stored `S` via `&mut`. Each worker owns its clone
/// exclusively (Parallel kind), so no synchronization is needed.
///
/// **Implementation note** (deviates from Phase 0 design line 136 which
/// proposed `OnceLock<S>`): `OnceLock` doesn't expose `&mut S`, but
/// consensus callers need mutable access to internal scratch buffers.
/// `Option<S>` + `get_or_insert_with` provides equivalent lazy-init
/// semantics with the required `&mut` access. Race-free because Parallel
/// workers each own their clone.
pub struct ProcessWithWorkerState<I, T, F, S, FInit>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(&mut S, I) -> io::Result<T> + Send + Sync + 'static,
    S: Send + HeapSize + 'static,
    FInit: Fn() -> S + Send + Sync + 'static,
{
    f: Arc<F>,
    init: Arc<FInit>,
    /// Lazy per-worker state. `Clone` resets to `None`.
    state: Option<S>,
    held: HeldSlot<Unpushed<T>>,
    name: &'static str,
    limit_bytes: u64,
    _phantom: PhantomData<fn() -> I>,
}

impl<I, T, F, S, FInit> ProcessWithWorkerState<I, T, F, S, FInit>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(&mut S, I) -> io::Result<T> + Send + Sync + 'static,
    S: Send + HeapSize + 'static,
    FInit: Fn() -> S + Send + Sync + 'static,
{
    pub fn new(name: &'static str, limit_bytes: u64, init: FInit, f: F) -> Self {
        Self {
            f: Arc::new(f),
            init: Arc::new(init),
            state: None,
            held: HeldSlot::new(),
            name,
            limit_bytes,
            _phantom: PhantomData,
        }
    }
}

impl<I, T, F, S, FInit> Clone for ProcessWithWorkerState<I, T, F, S, FInit>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(&mut S, I) -> io::Result<T> + Send + Sync + 'static,
    S: Send + HeapSize + 'static,
    FInit: Fn() -> S + Send + Sync + 'static,
{
    fn clone(&self) -> Self {
        Self {
            f: Arc::clone(&self.f),
            init: Arc::clone(&self.init),
            state: None, // <-- per-worker fresh; lazy on first try_run
            held: HeldSlot::new(),
            name: self.name,
            limit_bytes: self.limit_bytes,
            _phantom: PhantomData,
        }
    }
}

impl<I, T, F, S, FInit> Step for ProcessWithWorkerState<I, T, F, S, FInit>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(&mut S, I) -> io::Result<T> + Send + Sync + 'static,
    S: Send + HeapSize + 'static,
    FInit: Fn() -> S + Send + Sync + 'static,
{
    type Input = I;
    type Outputs = OrderedBytesSingle<T>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.limit_bytes }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }
        let Some(item) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        let state = self.state.get_or_insert_with(|| (self.init)());
        let out = (self.f)(state, item)?;
        match ctx.outputs.push(out) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

/// Convenience constructor for [`ProcessWithWorkerState`].
pub fn process_with_worker_state<I, T, F, S, FInit>(
    name: &'static str,
    limit_bytes: u64,
    init: FInit,
    f: F,
) -> ProcessWithWorkerState<I, T, F, S, FInit>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(&mut S, I) -> io::Result<T> + Send + Sync + 'static,
    S: Send + HeapSize + 'static,
    FInit: Fn() -> S + Send + Sync + 'static,
{
    ProcessWithWorkerState::new(name, limit_bytes, init, f)
}

// ─────────────────────────────────────────────────────────────────────────────
// Process2WithWorkerState<I, A, B, S, Init, F> — Parallel + lazy per-worker
// state + 2-output fan-out. Mirror of `ProcessWithWorkerState` widened to two
// ordered + byte-bounded output branches. Used by `CorrectStep` (kept +
// rejects) where each worker reuses one `LruCache` across all batches without
// `thread_local!` debt at every 2-output worker-state call site.
// ─────────────────────────────────────────────────────────────────────────────

/// `Parallel` + lazy per-worker state + 2-output fan-out. `Clone` resets
/// `state` to `None`; first `try_run` initializes via `init()`; subsequent
/// calls reuse the stored `S` via `&mut`. Each worker owns its clone
/// exclusively (Parallel kind), so no synchronization is needed.
///
/// **Implementation note** (deviates from Phase 0 design line 136 which
/// proposed `OnceLock<S>`): `OnceLock` doesn't expose `&mut S`, but
/// consensus callers need mutable access to internal scratch buffers.
/// `Option<S>` + `get_or_insert_with` provides equivalent lazy-init
/// semantics with the required `&mut` access. Race-free because Parallel
/// workers each own their clone. This applies to this variant for the same
/// reasons as the 1-output sibling [`ProcessWithWorkerState`].
pub struct Process2WithWorkerState<I, A, B, S, Init, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    S: Send + HeapSize + 'static,
    Init: Fn() -> S + Send + Sync + 'static,
    F: Fn(&mut S, I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    f: Arc<F>,
    init: Arc<Init>,
    state: Option<S>,
    held_a: HeldSlot<Unpushed<A>>,
    held_b: HeldSlot<Unpushed<B>>,
    name: &'static str,
    limit_bytes_a: u64,
    limit_bytes_b: u64,
    _phantom: PhantomData<fn() -> I>,
}

impl<I, A, B, S, Init, F> Process2WithWorkerState<I, A, B, S, Init, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    S: Send + HeapSize + 'static,
    Init: Fn() -> S + Send + Sync + 'static,
    F: Fn(&mut S, I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    pub fn new(
        name: &'static str,
        limit_bytes_a: u64,
        limit_bytes_b: u64,
        init: Init,
        f: F,
    ) -> Self {
        Self {
            f: Arc::new(f),
            init: Arc::new(init),
            state: None,
            held_a: HeldSlot::new(),
            held_b: HeldSlot::new(),
            name,
            limit_bytes_a,
            limit_bytes_b,
            _phantom: PhantomData,
        }
    }
}

impl<I, A, B, S, Init, F> Clone for Process2WithWorkerState<I, A, B, S, Init, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    S: Send + HeapSize + 'static,
    Init: Fn() -> S + Send + Sync + 'static,
    F: Fn(&mut S, I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    fn clone(&self) -> Self {
        Self {
            f: Arc::clone(&self.f),
            init: Arc::clone(&self.init),
            state: None, // per-worker fresh; lazy on first try_run
            held_a: HeldSlot::new(),
            held_b: HeldSlot::new(),
            name: self.name,
            limit_bytes_a: self.limit_bytes_a,
            limit_bytes_b: self.limit_bytes_b,
            _phantom: PhantomData,
        }
    }
}

impl<I, A, B, S, Init, F> Step for Process2WithWorkerState<I, A, B, S, Init, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    S: Send + HeapSize + 'static,
    Init: Fn() -> S + Send + Sync + 'static,
    F: Fn(&mut S, I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    type Input = I;
    type Outputs = crate::pipeline::core::outputs::OrderedBytesTuple2<A, B>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![
                QueueSpec::ByteBounded { limit_bytes: self.limit_bytes_a },
                QueueSpec::ByteBounded { limit_bytes: self.limit_bytes_b },
            ],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal, BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let view = ctx.outputs.view();

        if let Some(unpushed) = self.held_a.take() {
            match view.a.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held_a.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }
        if let Some(unpushed) = self.held_b.take() {
            match view.b.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held_b.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(item) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        let state = self.state.get_or_insert_with(|| (self.init)());
        let Process2Output { a, b } = (self.f)(state, item)?;

        if let Some(av) = a {
            if let Err(unpushed) = view.a.push(av) {
                self.held_a.put(unpushed);
            }
        }
        if let Some(bv) = b {
            if let Err(unpushed) = view.b.push(bv) {
                self.held_b.put(unpushed);
            }
        }
        Ok(StepOutcome::Progress)
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

/// Convenience constructor for [`Process2WithWorkerState`].
pub fn process2_with_worker_state<I, A, B, S, Init, F>(
    name: &'static str,
    limit_bytes_a: u64,
    limit_bytes_b: u64,
    init: Init,
    f: F,
) -> Process2WithWorkerState<I, A, B, S, Init, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    S: Send + HeapSize + 'static,
    Init: Fn() -> S + Send + Sync + 'static,
    F: Fn(&mut S, I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    Process2WithWorkerState::new(name, limit_bytes_a, limit_bytes_b, init, f)
}

// ─────────────────────────────────────────────────────────────────────────────
// MiAssign<I, F> — Serial + monotonic MI counter. Used for deterministic
// MoleculeId numbering across input batches.
// ─────────────────────────────────────────────────────────────────────────────

/// `Serial` closure-driven step that maintains a monotonic MI counter.
/// Used by `group` / `dedup` for deterministic `MoleculeId` numbering.
/// The closure receives `(&mut next_mi, input)` and returns the transformed
/// output (with MI tags applied to records).
pub struct MiAssign<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(&mut u64, I) -> io::Result<T> + Send + Sync + 'static,
{
    f: Arc<F>,
    next_mi: u64,
    held: HeldSlot<Unpushed<T>>,
    name: &'static str,
    limit_bytes: u64,
    _phantom: PhantomData<fn() -> I>,
}

impl<I, T, F> MiAssign<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(&mut u64, I) -> io::Result<T> + Send + Sync + 'static,
{
    pub fn new(name: &'static str, limit_bytes: u64, f: F) -> Self {
        Self {
            f: Arc::new(f),
            next_mi: 0,
            held: HeldSlot::new(),
            name,
            limit_bytes,
            _phantom: PhantomData,
        }
    }
}

impl<I, T, F> Step for MiAssign<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(&mut u64, I) -> io::Result<T> + Send + Sync + 'static,
{
    type Input = I;
    type Outputs = OrderedBytesSingle<T>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.limit_bytes }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }
        let Some(item) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        let out = (self.f)(&mut self.next_mi, item)?;
        match ctx.outputs.push(out) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }
}

/// Convenience constructor for [`MiAssign`].
pub fn mi_assign<I, T, F>(name: &'static str, limit_bytes: u64, f: F) -> MiAssign<I, T, F>
where
    I: Send + HeapSize + 'static,
    T: Send + HeapSize + Ordered + 'static,
    F: Fn(&mut u64, I) -> io::Result<T> + Send + Sync + 'static,
{
    MiAssign::new(name, limit_bytes, f)
}

// ─────────────────────────────────────────────────────────────────────────────
// Process2<I, A, B, F> — Parallel, fan-out to two unordered output branches.
// Used by `correct`/`filter`-style commands that emit a primary output plus
// a rejects/secondary output. Both branches are `BranchOrdering::None` (FIFO
// by arrival) so workers can race to push without ordinal coordination.
//
// Backpressure: each branch carries its own `HeldSlot`. If branch A pushes
// fail, the rejected item is parked in `held_a` and the next `try_run`
// drains it before popping new input. Branch B is independent. This means
// a stuck consumer on one side won't drop items — but it also means a
// permanently slow consumer will eventually block the step (held slot full,
// dispatch returns `Contention`).
// ─────────────────────────────────────────────────────────────────────────────

/// Output of a [`Process2`] closure: one optional value per output branch.
/// Either branch may be `None` to indicate "no item produced for this
/// branch this call" (filter / sparse fan-out).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Process2Output<A, B> {
    pub a: Option<A>,
    pub b: Option<B>,
}

impl<A, B> Process2Output<A, B> {
    #[must_use]
    pub fn both(a: A, b: B) -> Self {
        Self { a: Some(a), b: Some(b) }
    }

    #[must_use]
    pub fn only_a(a: A) -> Self {
        Self { a: Some(a), b: None }
    }

    #[must_use]
    pub fn only_b(b: B) -> Self {
        Self { a: None, b: Some(b) }
    }

    #[must_use]
    pub fn none() -> Self {
        Self { a: None, b: None }
    }
}

/// Closure-driven `Parallel` mid-step that fans out to two unordered output
/// branches `(A, B)`. The closure returns a [`Process2Output`] describing
/// what to push to each branch.
///
/// Both output branches use `BranchOrdering::None`. If the consumer of either
/// branch needs ordering, prefer composing single-output ordered steps
/// instead — `Process2` is for the sparse-output / rejects-path shape.
pub struct Process2<I, A, B, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    F: Fn(I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    f: Arc<F>,
    held_a: HeldSlot<Unpushed<A>>,
    held_b: HeldSlot<Unpushed<B>>,
    name: &'static str,
    capacity_a: usize,
    capacity_b: usize,
    _phantom: PhantomData<fn() -> I>,
}

impl<I, A, B, F> Process2<I, A, B, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    F: Fn(I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    pub fn new(name: &'static str, capacity_a: usize, capacity_b: usize, f: F) -> Self {
        Self {
            f: Arc::new(f),
            held_a: HeldSlot::new(),
            held_b: HeldSlot::new(),
            name,
            capacity_a,
            capacity_b,
            _phantom: PhantomData,
        }
    }
}

impl<I, A, B, F> Clone for Process2<I, A, B, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    F: Fn(I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    fn clone(&self) -> Self {
        Self {
            f: Arc::clone(&self.f),
            held_a: HeldSlot::new(),
            held_b: HeldSlot::new(),
            name: self.name,
            capacity_a: self.capacity_a,
            capacity_b: self.capacity_b,
            _phantom: PhantomData,
        }
    }
}

impl<I, A, B, F> Step for Process2<I, A, B, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    F: Fn(I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    type Input = I;
    type Outputs = (A, B);

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![
                QueueSpec::CountBounded { capacity: self.capacity_a },
                QueueSpec::CountBounded { capacity: self.capacity_b },
            ],
            branch_ordering: vec![BranchOrdering::None, BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let view = ctx.outputs.view();

        // Drain held items first. A branch that's still rejecting blocks
        // forward progress on this worker until the consumer drains.
        if let Some(unpushed) = self.held_a.take() {
            match view.a.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held_a.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }
        if let Some(unpushed) = self.held_b.take() {
            match view.b.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held_b.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(item) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        let Process2Output { a, b } = (self.f)(item)?;

        // Push each side. If either is rejected, hold and continue —
        // we still made progress (consumed input + produced what we could).
        if let Some(av) = a {
            if let Err(unpushed) = view.a.push(av) {
                self.held_a.put(unpushed);
            }
        }
        if let Some(bv) = b {
            if let Err(unpushed) = view.b.push(bv) {
                self.held_b.put(unpushed);
            }
        }
        Ok(StepOutcome::Progress)
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

/// Convenience constructor for [`Process2`].
pub fn process2<I, A, B, F>(
    name: &'static str,
    capacity_a: usize,
    capacity_b: usize,
    f: F,
) -> Process2<I, A, B, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    F: Fn(I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    Process2::new(name, capacity_a, capacity_b, f)
}

// ─────────────────────────────────────────────────────────────────────────────
// Process2Ordered<I, A, B, F> — Parallel, fan-out to two ordered + byte-bounded
// branches. Both branches use `BranchOrdering::ByItemOrdinal` so each output
// stream is reorderable into input order downstream. Use this when both
// outputs need to be BAM-sortable (e.g., `filter` emitting kept and rejected
// records that each feed an ordered `BgzfCompress` + `WriteBgzfFile` chain).
//
// Backpressure semantics mirror `Process2`: per-branch `HeldSlot` parks the
// rejected item and the next `try_run` retries before popping new input.
// ─────────────────────────────────────────────────────────────────────────────

/// Closure-driven `Parallel` mid-step that fans out to two ordered output
/// branches `(A, B)`, both `BranchOrdering::ByItemOrdinal`. The closure
/// returns a [`Process2Output`] describing what to push to each branch.
///
/// Both `A` and `B` must implement [`Ordered`] + [`HeapSize`]: each
/// branch's reorder stage uses `item.ordinal()` to release in input
/// order and the byte-bounded queue uses `heap_size()` for admission
/// control.
pub struct Process2Ordered<I, A, B, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    F: Fn(I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    f: Arc<F>,
    held_a: HeldSlot<Unpushed<A>>,
    held_b: HeldSlot<Unpushed<B>>,
    name: &'static str,
    limit_bytes_a: u64,
    limit_bytes_b: u64,
    _phantom: PhantomData<fn() -> I>,
}

impl<I, A, B, F> Process2Ordered<I, A, B, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    F: Fn(I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    pub fn new(name: &'static str, limit_bytes_a: u64, limit_bytes_b: u64, f: F) -> Self {
        Self {
            f: Arc::new(f),
            held_a: HeldSlot::new(),
            held_b: HeldSlot::new(),
            name,
            limit_bytes_a,
            limit_bytes_b,
            _phantom: PhantomData,
        }
    }
}

impl<I, A, B, F> Clone for Process2Ordered<I, A, B, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    F: Fn(I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    fn clone(&self) -> Self {
        Self {
            f: Arc::clone(&self.f),
            held_a: HeldSlot::new(),
            held_b: HeldSlot::new(),
            name: self.name,
            limit_bytes_a: self.limit_bytes_a,
            limit_bytes_b: self.limit_bytes_b,
            _phantom: PhantomData,
        }
    }
}

impl<I, A, B, F> Step for Process2Ordered<I, A, B, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    F: Fn(I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    type Input = I;
    type Outputs = crate::pipeline::core::outputs::OrderedBytesTuple2<A, B>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![
                QueueSpec::ByteBounded { limit_bytes: self.limit_bytes_a },
                QueueSpec::ByteBounded { limit_bytes: self.limit_bytes_b },
            ],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal, BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let view = ctx.outputs.view();

        if let Some(unpushed) = self.held_a.take() {
            match view.a.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held_a.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }
        if let Some(unpushed) = self.held_b.take() {
            match view.b.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held_b.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(item) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        let Process2Output { a, b } = (self.f)(item)?;

        if let Some(av) = a {
            if let Err(unpushed) = view.a.push(av) {
                self.held_a.put(unpushed);
            }
        }
        if let Some(bv) = b {
            if let Err(unpushed) = view.b.push(bv) {
                self.held_b.put(unpushed);
            }
        }
        Ok(StepOutcome::Progress)
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

/// Convenience constructor for [`Process2Ordered`].
pub fn process2_ordered<I, A, B, F>(
    name: &'static str,
    limit_bytes_a: u64,
    limit_bytes_b: u64,
    f: F,
) -> Process2Ordered<I, A, B, F>
where
    I: Send + HeapSize + 'static,
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    F: Fn(I) -> io::Result<Process2Output<A, B>> + Send + Sync + 'static,
{
    Process2Ordered::new(name, limit_bytes_a, limit_bytes_b, f)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, PartialEq, Eq)]
    struct Item {
        ord: u64,
        v: u32,
    }
    impl HeapSize for Item {
        fn heap_size(&self) -> usize {
            std::mem::size_of::<Self>()
        }
    }
    impl Ordered for Item {
        fn ordinal(&self) -> u64 {
            self.ord
        }
    }

    #[test]
    fn process_profile_is_parallel_fifo() {
        let s: Process<u32, u32, _> = process("test", 4, |x: u32| Ok(x + 1));
        let p = s.profile();
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::None]);
    }

    #[test]
    fn process_ordered_profile_is_parallel_byitem() {
        let s: ProcessOrdered<Item, Item, _> =
            process_ordered("test", 1024, |x: Item| Ok(Item { ord: x.ord, v: x.v + 1 }));
        let p = s.profile();
        assert_eq!(p.kind, StepKind::Parallel);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn process_with_worker_state_clone_resets_state() {
        let s = process_with_worker_state::<Item, Item, _, u64, _>(
            "stateful",
            1024,
            || 0u64, // init
            |state, item| {
                *state += 1;
                Ok(Item { ord: item.ord, v: u32::try_from(*state).unwrap_or(0) })
            },
        );
        let cloned = s.clone();
        assert!(s.state.is_none());
        assert!(cloned.state.is_none());
    }

    #[test]
    fn process2_profile_is_parallel_two_unordered_branches() {
        let s: Process2<u32, u32, String, _> =
            process2("test2", 4, 4, |x: u32| Ok(Process2Output::both(x + 1, x.to_string())));
        let p = s.profile();
        assert_eq!(p.name, "test2");
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky);
        assert_eq!(p.output_queues.len(), 2);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::None, BranchOrdering::None]);
    }

    #[test]
    fn process2_output_constructors() {
        let both: Process2Output<u32, &str> = Process2Output::both(1, "a");
        assert_eq!(both.a, Some(1));
        assert_eq!(both.b, Some("a"));

        let only_a: Process2Output<u32, &str> = Process2Output::only_a(2);
        assert_eq!(only_a.a, Some(2));
        assert_eq!(only_a.b, None);

        let only_b: Process2Output<u32, &str> = Process2Output::only_b("b");
        assert_eq!(only_b.a, None);
        assert_eq!(only_b.b, Some("b"));

        let none: Process2Output<u32, &str> = Process2Output::none();
        assert_eq!(none.a, None);
        assert_eq!(none.b, None);
    }

    #[test]
    fn process2_with_worker_state_clone_resets_state() {
        #[derive(Default)]
        struct WS {
            counter: u32,
        }
        impl HeapSize for WS {
            fn heap_size(&self) -> usize {
                0
            }
        }

        let step = process2_with_worker_state::<Item, Item, Item, WS, _, _>(
            "stateful2",
            1024,
            1024,
            WS::default,
            |state, item| {
                state.counter += 1;
                Ok(Process2Output::both(
                    Item { ord: item.ord, v: state.counter },
                    Item { ord: item.ord, v: state.counter },
                ))
            },
        );
        let cloned = step.clone();
        assert!(step.state.is_none());
        assert!(cloned.state.is_none());
    }

    #[test]
    fn process2_with_worker_state_profile_is_parallel_byitem_two_branches() {
        #[derive(Default)]
        struct WS;
        impl HeapSize for WS {
            fn heap_size(&self) -> usize {
                0
            }
        }

        let step = process2_with_worker_state::<Item, Item, Item, WS, _, _>(
            "stateful2-profile",
            512,
            1024,
            WS::default,
            |_state, item| Ok(Process2Output::both(item, Item { ord: 0, v: 0 })),
        );
        let profile = step.profile();
        assert_eq!(profile.kind, StepKind::Parallel);
        assert_eq!(profile.output_queues.len(), 2);
        assert_eq!(
            profile.output_queues,
            vec![
                QueueSpec::ByteBounded { limit_bytes: 512 },
                QueueSpec::ByteBounded { limit_bytes: 1024 },
            ]
        );
        assert_eq!(
            profile.branch_ordering,
            vec![BranchOrdering::ByItemOrdinal, BranchOrdering::ByItemOrdinal]
        );
    }

    #[test]
    fn mi_assign_profile_is_serial() {
        let s: MiAssign<Item, Item, _> = mi_assign("mi", 1024, |next_mi, item: Item| {
            let mi = *next_mi;
            *next_mi += 1;
            Ok(Item { ord: item.ord, v: u32::try_from(mi).unwrap_or(0) })
        });
        let p = s.profile();
        assert_eq!(p.kind, StepKind::Serial);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }
}
