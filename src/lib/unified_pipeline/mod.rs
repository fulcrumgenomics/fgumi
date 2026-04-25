//! Unified thread pool pipeline for `--threads N` mode.
//!
//! This module provides a true unified thread pool where all N threads can perform
//! any type of work: reading, processing, or writing. This ensures `--threads N`
//! is a strict thread cap with no separate I/O thread pools.
//!
//! # Architecture
//!
//! ```text
//! ┌─────────────────────────────────────────────────────────────────────────┐
//! │                   UNIFIED THREAD POOL (N threads)                       │
//! │                                                                         │
//! │  All threads can do ANY work: READ → PROCESS → WRITE                    │
//! │                                                                         │
//! │  Priority scheduling based on queue depths:                             │
//! │  - Input queue low → prioritize reading                                 │
//! │  - Output queue high → prioritize writing                               │
//! │  - Otherwise → prioritize processing                                    │
//! │                                                                         │
//! │  File access via parking_lot::Mutex with try_lock()                     │
//! │  Lock-free queues via crossbeam ArrayQueue                              │
//! └─────────────────────────────────────────────────────────────────────────┘
//! ```
//!
//! # Pipeline Steps
//!
//! Both BAM and FASTQ pipelines use a 9-step structure:
//!
//! ```text
//! Read → Decompress → FindBoundaries → Parse → Group → Process → Serialize → Compress → Write
//!        [parallel]   [sequential]     [par]   [seq]   [parallel] [parallel]  [parallel] [seq]
//! ```
//!
//! - **Read**: Read raw bytes from input (sequential, requires file lock)
//! - **Decompress**: Decompress BGZF blocks (parallel)
//! - **`FindBoundaries`**: Find record boundaries in decompressed data (sequential)
//! - **Parse**: Construct record objects from boundaries (parallel - key optimization!)
//! - **Group**: Group records by position/UMI/template (sequential)
//! - **Process**: Apply domain-specific processing (parallel)
//! - **Serialize**: Convert output to bytes (parallel)
//! - **Compress**: Compress to BGZF blocks (parallel)
//! - **Write**: Write to output file (sequential, requires file lock)
//!
//! # FASTQ Parallel Parse Optimization
//!
//! The FASTQ pipeline includes a parallel Parse step that can improve
//! thread scaling, especially from t4 to t8. This was implemented because profiling
//! showed that 99.6% of Group step time was spent parsing FASTQ records.
//!
//! When `FastqPipelineConfig::use_parallel_parse` is true:
//! - `FindBoundaries` step scans for FASTQ record boundaries (fast, O(N) scan)
//! - Parse step constructs `FastqRecord` objects in parallel
//! - Group step receives pre-parsed records (no parsing under lock!)
//!
//! # Module Structure
//!
//! - `base`: Core infrastructure, traits, and shared types
//! - `bam`: BAM pipeline implementation
//! - `fastq`: FASTQ pipeline with multi-stream grouping and parallel Parse
//! - `scheduler`: Thread scheduling strategies
//! - `deadlock`: Deadlock detection and recovery
//! - `queue`: Queue implementations
//! - `rebalancer`: Dynamic memory rebalancing
//!
//! # Parallel Ordered Batch Processing Pattern
//!
//! Both BAM and FASTQ pipelines use a common pattern for steps that require
//! ordered output but can do work in parallel:
//!
//! 1. **Per-worker held state**: Each worker holds its result if output queue is full
//! 2. **Brief reorder lock**: Lock held only for insert/pop, not during work
//! 3. **Work outside lock**: Actual processing happens without holding locks
//! 4. **Priority advancement**: Always try to push held items first
//!
//! This pattern is implemented in:
//! - `bam.rs`: `try_step_find_boundaries()` with `WorkerState.held_boundaries`
//! - `fastq.rs`: `fastq_try_step_find_boundaries()` with `FastqWorkerState.held_boundaries`
//!
//! When modifying either implementation, ensure the pattern stays in sync.
//! The `HasHeldBoundaries` trait in `base.rs` documents the interface.
//!
//! # Adding a New Pipeline Type
//!
//! To add a new pipeline (e.g., for a new input format):
//!
//! 1. **Define your pipeline state struct** implementing:
//!    - [`PipelineLifecycle`] — completion, error, drain mode, validation
//!    - [`MonitorableState`] — if using the shared monitor loop
//!    - [`OutputPipelineState`] — if writing BAM/BGZF output
//!    - [`ProcessPipelineState`] — for the process step
//!    - [`SerializePipelineState`] — for the serialize step
//!    - [`WritePipelineState`] — for the write step
//!
//! 2. **Define your worker state struct** implementing:
//!    - [`WorkerStateCommon`] + [`HasWorkerCore`] — required for all workers
//!    - `HasHeld*` traits — one per non-blocking step your pipeline uses
//!    - [`HasCompressor`] + [`HasRecycledBuffers`] — if writing compressed output
//!
//! 3. **Implement [`StepContext`]** to plug into `generic_worker_loop`, or write
//!    a custom worker loop.
//!
//! 4. **Reuse shared step functions** (`shared_try_step_compress`, etc.) where
//!    possible — they handle non-blocking held-item logic correctly.
//!
//! See `bam.rs` and `fastq.rs` for complete examples of this pattern.

mod bam;
mod base;
pub mod deadlock;
mod fastq;
pub mod queue;
pub mod rebalancer;
pub mod scheduler;
pub mod serial_ordered_array_queue;
pub mod serialize_context;

// Re-export everything from base
pub use base::*;

// Re-export everything from bam
pub use bam::*;

// Re-export everything from fastq
pub use fastq::*;

// Re-export queue types
pub use queue::{OrderedQueue, QueueStats};

// Re-export rebalancer types
pub use rebalancer::{
    DynamicRebalancer, InitialAllocation, RebalancerConfig, initial_allocation_for_command,
    parse_memory_limit,
};

// Re-export scheduler types
pub use scheduler::{
    BackpressureState, ChaseBottleneckScheduler, FixedPriorityScheduler, Scheduler,
    SchedulerStrategy, create_scheduler,
};

// Re-export deadlock detection types
pub use deadlock::{
    DeadlockAction, DeadlockConfig, DeadlockState, QueueSnapshot, check_deadlock_and_restore,
};
