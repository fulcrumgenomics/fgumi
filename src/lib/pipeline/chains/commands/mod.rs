//! Per-command chain builders. Each module exports a single
//! `build_<command>_chain(spec: ChainSpec) -> Result<BuiltPipeline>`
//! function. [`super::build`] dispatches to these by matching on
//! `spec.stages`.
//!
//! The bodies are extracted from each command's existing inline
//! chain construction during Phase 2 (T2.12–T2.22 = move). Phase 3
//! will inline these into `super::build`'s match arms directly.

pub mod align;
pub mod clip;
pub mod codec;
pub mod correct;
pub mod dedup;
pub mod duplex;
pub mod extract;
pub mod filter;
pub mod group;
pub mod simplex;
pub mod sort;
pub mod zipper;
