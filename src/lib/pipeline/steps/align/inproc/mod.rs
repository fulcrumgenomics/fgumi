//! In-process bwa-mem3 aligner backend (feature `aligner-bwa-mem3`).
//!
//! This module will replace the subprocess `bwa-mem3` invocation with bwa-mem3
//! linked into the fgumi binary through [`bwa_mem3_rs`], run as work items on
//! the pipeline's shared work-stealing pool. The whole module tree is compiled
//! only under `aligner-bwa-mem3`; default builds stay pure Rust with no C
//! toolchain (the module is gated at its declaration in the parent `align`
//! module).
//!
//! This change adds the foundation the pipeline steps build on:
//!
//! - [`cohort`]: the parity-critical pure logic — the `-K` even-parity cohort
//!   cut, the SE/PE template classification, and the per-sub-batch id-base
//!   arithmetic — each verified by a proptest against a literal Rust port of the
//!   corresponding upstream bwa-mem3 rule;
//! - [`engine`]: the `AlignEngine` abstraction over bwa-mem3-rs's three-phase
//!   API (seed/extend, per-cohort pestat, pair/emit), with a test fake;
//! - [`gate`]: the cohort-granularity in-flight gate and its lease;
//! - [`scratch`]: the per-pool-thread aligner scratch pool.
//!
//! The pipeline steps and the backend that wire these into `runall` land in the
//! following change.

// Wired into the pipeline by the in-process backend in the following change;
// until then these items are exercised only by their unit tests.
#![allow(dead_code)]

pub(crate) mod cohort;
pub(crate) mod engine;
pub(crate) mod gate;
pub(crate) mod scratch;
