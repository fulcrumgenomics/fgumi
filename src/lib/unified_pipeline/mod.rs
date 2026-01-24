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
//! - `queue`: Memory-bounded queue implementations
//! - `rebalancer`: Dynamic memory rebalancing

mod bam;
mod base;
pub mod deadlock;
mod fastq;
pub mod queue;
pub mod rebalancer;
pub mod scheduler;

// Re-export everything from base
pub use base::*;

// Re-export everything from bam
pub use bam::*;

// Re-export everything from fastq
pub use fastq::*;

// Re-export queue types
pub use queue::{MemoryBoundedQueue, OrderedQueue, QueueStats};

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
