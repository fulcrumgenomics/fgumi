//! The typed-step pipeline framework for `--threads N` mode.
//!
//! A pipeline is a chain of [`core::step::Step`]s, each declaring a
//! [`core::step::StepKind`] (`Serial` / `Parallel` / `Exclusive`); the engine
//! runs all `N` worker threads through a round-robin dispatch loop with
//! `Affinity`-based pinning for I/O sources/sinks. This keeps `--threads N` a
//! strict thread cap with no separate I/O thread pools.
//!
//! # Module structure
//!
//! - [`core`]: typed-step execution engine (worker loop, queues, drain
//!   protocol, reorder buffers).
//! - [`steps`]: the concrete `Step` implementations (decompress, boundaries,
//!   parse, serialize, compress, write, …).
//! - [`chains`]: declarative chain construction — `build_for`, `ChainBuilder`,
//!   per-command step factories that assemble the steps into runnable pipelines.
//!
//! Commands are rewired onto `chains` one at a time; until a given command is
//! rewired it still runs on `unified_pipeline`.

/// The typed-step execution engine, extracted into the `fgumi-pipeline-core`
/// crate so its lightweight dependency graph (no `noodles`-bam / sort /
/// consensus) compiles fast in isolation. Re-exported here so every
/// `crate::pipeline::core::…` path resolves — which is what lets the ported
/// step sources compile unmodified.
pub use fgumi_pipeline_core as core;
pub mod backpressure;
pub mod chains;
pub mod steps;
