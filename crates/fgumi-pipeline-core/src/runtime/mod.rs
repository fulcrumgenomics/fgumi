//! Runtime: per-worker step storage, chain contexts, drain coordination,
//! worker pool, worker loop body. Built on top of Phase 1's trait surface.

pub(crate) mod contexts;
pub(crate) mod detached;
pub(crate) mod drain;
pub(crate) mod driver;
pub mod event_count;
pub(crate) mod fused;
pub(crate) mod live;
pub mod metrics;
pub(crate) mod placement;
pub(crate) mod pool;
pub mod sampler;
pub mod scheduler;
pub mod stats;
pub(crate) mod storage;
pub mod telemetry;
pub(crate) mod wake;
pub(crate) mod wake_slot;
pub(crate) mod worker_core;
pub mod worker_state;

pub use contexts::StepCounters;
pub(crate) use contexts::build_chain_contexts;
pub(crate) use detached::extract_detached_steps;
pub(crate) use detached::run_detached_driver;
pub(crate) use drain::StepDrainCounter;
pub(crate) use driver::run_worker_loop;
pub use event_count::{NotifyOutcome, PoolEventCount, WaitKey, WaitOutcome};
pub(crate) use fused::{run_fused_single_thread, should_fuse_single_thread};
pub(crate) use pool::{assign_exclusive_owners, assign_sticky_owners};
pub use scheduler::{
    ChainOrderScheduler, DrainFirstScheduler, RefillDrainScheduler, Scheduler, WalkDirection,
};
pub use stats::{PipelineStats, StatsSnapshot, StepStatsSnapshot};
pub(crate) use storage::build_worker_storage;
pub(crate) use worker_core::WorkerCore;
#[cfg(any(test, feature = "test-utils"))]
#[doc(hidden)]
pub use worker_core::{TestBackoff, TestBackoffTarget};
