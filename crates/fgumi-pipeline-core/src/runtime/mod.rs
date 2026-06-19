//! Runtime: per-worker step storage, chain contexts, drain coordination,
//! worker pool, worker loop body. Built on top of Phase 1's trait surface.

pub mod contexts;
pub mod drain;
pub mod driver;
pub mod fused;
pub mod live;
pub mod pool;
pub mod stats;
pub mod storage;
pub mod worker_core;

pub use contexts::{ChainContexts, build_chain_contexts, build_chain_contexts_fused};
pub use drain::StepDrainCounter;
pub use driver::run_worker_loop;
pub use fused::{is_fusible_chain, run_fused_single_thread};
pub use live::LiveSteps;
pub use pool::{assign_exclusive_owners, assign_sticky_owners};
pub use stats::{PipelineStats, StatsSnapshot, StepStatsSnapshot};
pub use storage::{WorkerStepEntry, build_worker_storage};
pub use worker_core::WorkerCore;
