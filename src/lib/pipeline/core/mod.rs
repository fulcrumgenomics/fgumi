//! Core types and traits for the new step-DSL pipeline framework.
//!
//! This module is added by issue #330 PR 1 alongside the existing per-format
//! `bam.rs` / `fastq.rs` machinery. Subsequent PRs incrementally migrate
//! commands from the old framework to this one; PR 6 deletes the old framework.
//!
//! Reference design: `docs/design/unified-pipeline-step-dsl.md`.

pub mod builder;
pub mod erased;
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
