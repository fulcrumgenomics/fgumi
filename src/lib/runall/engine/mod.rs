//! Composable pipeline engine.
//!
//! Provides a linear chain of stages (Source -> Stages -> Sink) connected by
//! bounded memory-tracked queues. One persistent thread pool executes a
//! generalized `BalancedChaseDrain` scheduler.
//!
//! # Stage traits: `Stage` vs `SpecialStage`
//!
//! The engine recognises two stage traits:
//!
//! - [`Stage`] — the default. Pool-driven, one `process()` call per input
//!   with at most one output. Covers ~all stages in [`stages`].
//! - [`SpecialStage`] — the escape hatch. Owns its own loop on a dedicated
//!   thread. Use only for 1:N emission, barrier accumulation, or stages
//!   with internal threads / subprocesses (sort, aligner, coalesce,
//!   position-batch).
//!
//! See the module-level docs on [`stage`] and [`special_stage`] for the
//! decision tree and trade-offs before picking one.
//!
//! # Design reference
//!
//! See `docs/superpowers/specs/2026-04-13-composable-pipeline-design.md`.

pub mod errors;
pub use errors::PipelineError;

pub mod stage;
pub use stage::{Parallelism, SequencedItem, Stage};

pub mod memory;
pub use memory::MemoryTracker;

pub mod queue;
pub use queue::StageQueue;

pub mod backoff;
pub use backoff::Backoff;

pub mod backpressure;
pub use backpressure::{BackpressureState, QueueSummary};

pub mod counting;
pub use counting::{CountingInputQueue, CountingOutputQueue};

pub mod cancel;
pub use cancel::CancelToken;

pub mod sink;
pub mod source;
pub use sink::{InputQueue, Sink};
pub use source::{AnySource, OutputQueue, Source};

pub mod erased_queue;
pub use erased_queue::{ErasingOutputQueue, UnerasingInputQueue};

pub mod worker;

pub mod scheduler;
pub use scheduler::{Scheduler, StageAttemptOutcome, WorkerRole, assign_sequential_ownership};

pub mod deadlock;
pub use deadlock::{QueueProgress, WatchdogConfig, WatchdogOutcome, watch};

pub mod grouping_types;

pub mod reorder;
pub use reorder::ReorderBuffer;

pub mod special_stage;
pub use special_stage::SpecialStage;

pub mod driver;
pub use driver::{
    AUTO_TUNED_QUEUE_CAPACITY_MAX, AUTO_TUNED_QUEUE_CAPACITY_MIN, AUTO_TUNED_SLOTS_PER_THREAD,
    DEFAULT_GLOBAL_MEMORY_LIMIT, DEFAULT_QUEUE_CAPACITY, DEFAULT_QUEUE_MEMORY_LIMIT, ErasedQueue,
    ErasedStage, PipelineConfig, StageEntry, TypedStage, run_multi_stage,
    run_multi_stage_with_metrics, run_multi_stage_with_watchdog, run_single_stage,
};

pub mod builder;
pub use builder::{Pipeline, PipelineBuilder};

pub mod stage_source;
pub use stage_source::{ErasedStageSource, StageSource};

pub mod test_support;

pub mod pool;
pub use pool::Pool;

pub mod stats;
pub use stats::PipelineStats;

pub mod stages;

pub mod fastq_types;
pub use fastq_types::{
    BlockParsed, FastqParsedStream, FastqTemplateV2, OwnedFastqRecord, PerStreamChunk,
    TemplateBatch,
};

pub mod output_types;
pub use output_types::{
    CompressedBatch, CompressedBytes, RawBytes, SerializedBatch, SerializedFastqBatch,
};
