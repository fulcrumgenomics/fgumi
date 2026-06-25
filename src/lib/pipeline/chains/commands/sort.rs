//! Chain builder for `Stage::Sort`.
//!
//! Phase 2 (T2.16) held the full ~200-LOC chain construction here.
//! Phase 3 (T3a.6) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; the sort-specific
//! types and step factory the builder imports (`SortFinalizeHook`,
//! `IndexBamFinalizeHook`, `build_sort_step`, `log_sort_start`) now live in
//! the `fgumi-sort-cli` crate and are re-exported from here. This module
//! itself retains only the `build_sort_chain` delegate and the
//! source/sink path helpers.
//!
//! [`build_sort_chain`] shrinks to a ~10-line delegate.
//!
//! ## Sort pipeline topology
//!
//! Sort is architecturally unique among the `Stage` variants: it is
//! implemented as a single `Exclusive` [`SortBamFile`] step that reads
//! from and writes to files directly. It therefore bypasses the normal
//! source-preamble (`ReadBgzfBlocks → … → DecodeRecords`) and
//! `BgzfCompress → WriteBgzfFile` sink that every other stage uses.
//! `build_sort_chain` consequently skips `add_source` and `add_sink`;
//! `ChainBuilder::add_sort` registers `SortBamFile` via
//! `PipelineBuilder::append_source` (legal since `SortBamFile` has
//! `Input = ()`) and sets `ChainBuilder::override_pipeline_threads = 1`
//! so the framework runs a single driver thread while the sorter's own
//! `SortWorkerPool` / rayon pools provide internal concurrency.
//!
//! [`SortBamFile`]: crate::pipeline::steps::sort::SortBamFile

use anyhow::{Result, bail};

use crate::pipeline::chains::builder::ChainBuilder;
use crate::pipeline::chains::{BuiltPipeline, ChainSpec, SinkSpec, SourceSpec};

// The sort step factory, its captures bundle, the startup-banner logger, AND
// the two finalize hooks (`SortFinalizeHook`, `IndexBamFinalizeHook`) all live
// in the `fgumi-sort-cli` crate. Re-export them so the umbrella's
// `ChainBuilder::add_sort` keeps constructing and registering them from this
// module. The hooks implement the shared `FinalizeHook` trait
// (`fgumi-pipeline-core`, re-exported as `crate::pipeline::chains::FinalizeHook`),
// so the umbrella registers them in its `Vec<Box<dyn FinalizeHook>>` directly —
// no umbrella-local copies (X1-005).
pub(crate) use fgumi_sort_cli::chains::{
    IndexBamFinalizeHook, SortFinalizeHook, SortStepCaptures, build_sort_step, log_sort_start,
};

// ─────────────────────────────────────────────────────────────────────────────
// build_sort_chain — 10-line delegate
// ─────────────────────────────────────────────────────────────────────────────

/// Build a [`BuiltPipeline`] for a sort-only chain.
///
/// `spec.stages == [Stage::Sort]`. Other layouts are caller errors
/// caught by the dispatch match in [`crate::pipeline::chains::build`].
///
/// Delegates to [`ChainBuilder`]. The full step construction lives in
/// [`ChainBuilder`]'s `add_sort` method. Unlike other stage chains, sort
/// bypasses `add_source` / `add_sink` because `SortBamFile` is a
/// self-contained `Exclusive` step that reads/writes files directly.
///
/// # Errors
///
/// Returns validation errors (missing sort options, invalid source/sink
/// types) or any underlying pipeline construction error.
#[allow(clippy::needless_pass_by_value)]
pub fn build_sort_chain(spec: ChainSpec) -> Result<BuiltPipeline> {
    use crate::pipeline::chains::{Stage, builder::StagePosition};
    let mut chain = ChainBuilder::new(&spec)?;
    chain.add_stage(Stage::Sort, StagePosition::Terminal)?;
    chain.build()
}

// ─────────────────────────────────────────────────────────────────────────────
// Source/sink path helpers re-used by add_sort
// ─────────────────────────────────────────────────────────────────────────────

/// Extract the source path from a [`ChainSpec`] for sort (BAM or SAM only).
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_sort`].
pub(crate) fn sort_source_path(spec: &ChainSpec) -> Result<std::path::PathBuf> {
    match &spec.source {
        SourceSpec::Bam(p) | SourceSpec::Sam(p) => Ok(p.clone()),
        other => bail!("sort requires a BAM or SAM source, got {other:?}"),
    }
}

/// Extract the sink path from a [`ChainSpec`] for sort.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_sort`].
pub(crate) fn sort_sink_path(spec: &ChainSpec) -> std::path::PathBuf {
    match &spec.sink {
        SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p.clone(),
    }
}
