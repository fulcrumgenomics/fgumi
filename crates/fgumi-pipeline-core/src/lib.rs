//! Core types and traits for the typed-step pipeline framework.
//!
//! Extracted from the `fgumi` crate's `pipeline::core` module into its own
//! crate so its lightweight dependency graph (no `noodles`-bam, sort, or
//! consensus crates — just `ahash` / `crossbeam-queue` / `parking_lot` / `log`
//! and one `noodles::sam` type) compiles fast in isolation. The `fgumi` crate
//! re-exports this as `pipeline::core`, so existing `crate::pipeline::core::…`
//! paths are unchanged; the `chains` and `steps` layers that drive these
//! primitives still live in `fgumi`.
//!
//! Reference design: `docs/design/unified-pipeline-step-dsl.md`.

#![deny(unsafe_code)]

pub mod builder;
pub mod erased;
pub mod finalize;
pub mod handles;
pub mod header;
pub mod held;
pub mod item;
pub mod outputs;
pub mod queues;
pub mod reorder;
pub mod runtime;
pub mod signal;
pub mod step;
pub mod topology;

#[cfg(test)]
mod tests;

pub use builder::{
    BuildError, Chain, MultiChain2, MultiChain3, MultiChain4, Pipeline, PipelineBuilder,
    PipelineConfig,
};
pub use erased::{ErasedStep, ErasedStepCtx, TypedStep, TypedStep2};
pub use finalize::FinalizeHook;
pub use handles::{HeldRetry, TwoInputHandles, Unpushed};
pub use header::{AlreadySetError, HeaderHandle};
pub use held::HeldSlot;
pub use item::{HeapSize, Ordered};
pub use outputs::{MAX_ARITY, OrderedBytesSingle, Single, StepOutputs};
pub use queues::{ByteBoundedQueue, CountBoundedQueue, ItemQueue, QueueSpec, UnboundedQueue};
pub use reorder::{BranchOrdering, ReorderStage, Sequenced};
pub use signal::{CancelHandle, PipelineError, PipelineSignal};
pub use step::{
    Affinity, InputHandle, OutputHandles, OutputsViewAny, Step, Step2, StepCtx, StepCtx2, StepKind,
    StepOutcome, StepProfile,
};
pub use topology::{BranchIdx, ChainGraph, StepIdx};
