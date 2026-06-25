//! Post-pipeline finalize hook contract.
//!
//! [`FinalizeHook`] is the bare trait that lets a chain register
//! heterogeneous post-pipeline cleanup actions — metrics drains, summary
//! logging, gate checks, rejects-file finalization, BAI indexing — in a single
//! `Vec<Box<dyn FinalizeHook>>` that the caller drains after `Pipeline::run`
//! returns.
//!
//! It lives in `fgumi-pipeline-core` (rather than the umbrella `fgumi` crate)
//! so that any CLI crate built directly on the pipeline core — e.g.
//! `fgumi-sort-cli` — can implement the same contract without redefining a
//! parallel trait. Keeping the trait single-sourced here is what lets the sort
//! finalize hooks be defined once and re-exported by the umbrella, instead of
//! duplicated across both crates (X1-005).
//!
//! The richer machinery built on top of this trait — `BuiltPipeline`,
//! `drain_finalize`, and the built-in stats/timing hooks — stays in the
//! umbrella crate, since it references umbrella-only types.

use anyhow::Result;

/// Post-pipeline cleanup. Each stage with metrics, summary logging,
/// `--min-corrected`-style gates, or rejects-file finalization registers one or
/// more hooks during chain build. The caller iterates the registered hooks
/// (typically via `try_for_each`) after `Pipeline::run` returns.
///
/// Trait object so heterogeneous hooks (correct's metrics drain, consensus's
/// summary log, AAM's records-aligned counter, sort's BAI indexer, etc.) can
/// share a `Vec<Box<dyn FinalizeHook>>`.
pub trait FinalizeHook: Send {
    /// Run the post-pipeline action. Called exactly once, after
    /// `Pipeline::run` returns and before the caller exits.
    ///
    /// # Errors
    ///
    /// Returns the underlying I/O / aggregation / gate-check error.
    /// Callers typically chain via `try_for_each`.
    fn finalize(self: Box<Self>) -> Result<()>;
}
