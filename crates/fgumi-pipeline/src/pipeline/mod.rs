//! Pipeline orchestration: builder, scheduler, and progress tracking.

pub mod aligner_io;
pub mod batch_state;
pub mod builder;
pub mod health;
pub mod phase;
pub mod phase_a;
pub mod phase_a_graph;
pub mod phase_a_steps;
pub mod phase_a_worker;
pub mod queue;
pub mod reorder;
pub mod scheduler;
pub mod stats;
pub mod worker_util;
pub use builder::{PipelineBuilder, PipelineConfig, Stage, TwoPhaseConfig};
pub use health::HealthMonitor;
pub use queue::WorkQueue;
pub use reorder::ReorderBuffer;
pub use scheduler::{Scheduler, StageGraph, WorkerState, WriteState};
pub use stats::StageStats;
