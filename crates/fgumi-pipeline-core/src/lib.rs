//! Core types and traits for the typed-step pipeline framework.
//!
//! A pipeline is a graph of typed [`Step`]s joined by bounded [`queues`] and run
//! by a work-stealing worker pool. This crate carries only the framework
//! primitives — the step traits, the queues, the [`reorder`] stage that restores
//! input order across parallel branches, and the [`runtime`] that schedules and
//! drives them. Nothing here reads or writes sequencing data, so the dependency
//! graph stays light (`ahash` / `crossbeam-queue` / `parking_lot` / `log`, plus
//! `anyhow` for the one [`FinalizeHook::finalize`] return type and one
//! `noodles::sam` type for the shared header handle; `loom` only under
//! `--cfg loom`, for the admission model check) and the crate compiles fast in
//! isolation.
//!
//! The concrete steps that do the I/O and the computation live elsewhere and
//! plug in by implementing [`Step`] (one input) or [`Step2`] (two inputs);
//! [`PipelineBuilder`] wires them into a [`Pipeline`].

#![deny(unsafe_code)]

pub mod admission;
pub mod builder;
pub mod erased;
pub mod finalize;
pub mod handles;
pub mod header;
pub mod held;
pub mod item;
pub mod liveness;
pub mod outputs;
mod padded;
pub mod queues;
pub mod reorder;
pub mod runtime;
pub mod signal;
pub mod step;
#[cfg(any(test, feature = "test-utils"))]
pub mod testing;
pub mod topology;

#[cfg(test)]
mod tests;

pub use admission::{CapPermit, CapWaker, MAX_PHASE_CAPS, PhaseCap, WholeCap, admit_input};
pub use builder::{
    BuildError, Chain, DEFAULT_DEADLOCK_TIMEOUT_SECS, InstrumentationLevel, MultiChain2,
    MultiChain2Ordered, MultiChain3, MultiChain4, Pipeline, PipelineBuilder, PipelineConfig,
};
pub use erased::{ErasedStep, ErasedStepCtx, TypedStep, TypedStep2, TypedStepK};
pub use finalize::FinalizeHook;
pub use handles::{
    BranchInputHandle, HeldRetry, KInputHandles, OutputQueueSet, Tuple2View, Tuple3View,
    Tuple4View, TwoInputHandles, Unpushed,
};
pub use header::{AlreadySetError, HeaderHandle};
pub use held::HeldSlot;
pub use item::{HeapSize, Ordered};
pub use liveness::LivenessCounter;
pub use outputs::{
    MAX_ARITY, OrderedBytesSingle, OrderedBytesTuple2, OrderedBytesTuple3, Single, StepOutputs,
};
pub use queues::{ByteBoundedQueue, CountBoundedQueue, ItemQueue, QueueSpec, UnboundedQueue};
pub use reorder::{BranchOrdering, ReorderStage, Sequenced};
pub use runtime::wake::{PoolHandle, PoolRequest};
// The real wake protocols, for `tests/loom_wake.rs` (`wake_slot` itself stays
// crate-private).
#[cfg(loom)]
#[doc(hidden)]
pub use runtime::{
    wake::{DriverIdx, WakeEdges, WakePlan},
    wake_slot::{HolderSet, ThreadSlots, current_thread, delivery_fence},
};
pub use signal::{CancelHandle, PipelineError, PipelineSignal};
pub use step::{
    Affinity, CounterSpec, DetachedGroup, InputHandle, OutputHandles, OutputsViewAny,
    PoolPlacement, Step, Step2, StepCtx, StepCtx2, StepCtxK, StepK, StepKind, StepOutcome,
    StepProfile,
};
pub use topology::{BranchIdx, ChainGraph, StepIdx};
