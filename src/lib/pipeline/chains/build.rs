//! [`build_for`] — the single chain-construction entry point.
//!
//! `build_for` constructs the chain as a stage-by-stage loop:
//!
//!   1. Validates the spec (progression, options presence, cross-stage constraints).
//!   2. Constructs a [`ChainBuilder`] from the spec.
//!   3. Calls `add_source`.
//!   4. Walks `spec.stages`, calling `add_stage(stage, position)` for each,
//!      where `position` is `Terminal` for the last stage and `Intermediate` for all
//!      earlier ones.
//!   5. Calls `add_sink`.
//!   6. Calls `build()` to produce the [`BuiltPipeline`].
//!
//! Every chain — including the sole-`[Stage::Sort]` standalone sort — runs
//! through this uniform `add_source` → stages → `add_sink` topology on the
//! work-stealing pool. There is no longer a Sort-terminal special case: the
//! former self-contained file→file `SortBamFile` step was removed, and
//! standalone sort now streams (`ParseBamRecords` → arena sort ingest →
//! `SpillGather` → `SpillCompress` → `SpillWrite` → `SortSpillDecompress` →
//! `SortMerge` → sink) like every other command.

use anyhow::{Result, bail};

use crate::pipeline::chains::{
    BuiltPipeline, ChainSpec,
    builder::{ChainBuilder, StagePosition},
    validate::{
        validate_cross_stage_constraints, validate_stage_opts_present, validate_stage_progression,
    },
};

/// Construct a [`BuiltPipeline`] from a [`ChainSpec`]. See the design
/// doc (`docs/design/refactor/2026-05-20-unified-chain-builder-design.md`)
/// for the architectural rationale.
///
/// Does NOT execute the pipeline — the caller picks `PipelineConfig`
/// (threads, scheduler) and calls `built.pipeline.run(built.config)`,
/// then drains `built.finalize` in order. The convenience method
/// [`BuiltPipeline::run`] wraps both steps.
///
/// ## Dispatch
///
/// `build_for` constructs a [`ChainBuilder`] and drives it stage-by-stage:
///
/// ```text
/// chain.add_source()?
/// for (i, stage) in stages.iter().enumerate() {
///     let pos = if i == last { Terminal } else { Intermediate };
///     chain.add_stage(*stage, pos)?;
/// }
/// chain.add_sink()?
/// chain.build()
/// ```
///
/// Each `add_<stage>` method reads its options from `spec.stage_opts` and pushes
/// the appropriate typed-step sequence onto the pipeline builder. The position
/// argument controls whether a serialise-to-bytes step is appended (Terminal)
/// or whether the native typed output is left for the next stage (Intermediate).
///
/// # Errors
///
/// Returns the first validator violation, or any underlying step
/// construction error (header parse, file-open, etc.).
#[allow(clippy::needless_pass_by_value)]
pub fn build_for(spec: ChainSpec) -> Result<BuiltPipeline> {
    validate_stage_progression(&spec)?;
    validate_stage_opts_present(&spec)?;
    validate_cross_stage_constraints(&spec)?;

    if spec.stages.is_empty() {
        bail!("ChainSpec must specify at least one stage");
    }

    let mut chain = ChainBuilder::new(&spec)?;

    // NOTE: `--threads 1` whole-chain fusion is NOT wired here. It is applied
    // generically by the runtime (`Pipeline::run` → `run_fused_single_thread`)
    // to ANY linear source→sink chain at one worker — see
    // `pipeline::core::runtime::fused`. `build_for` therefore constructs the
    // same staged chain regardless of thread count; the runtime fuses it.
    //
    // Every chain — including the sole-`[Stage::Sort]` chain — runs through the
    // normal `add_source` → stages → `add_sink` topology on the work-stealing
    // pool. The former Sort-terminal special case (which skipped source/sink so
    // a self-contained `SortBamFile` could read/write files itself) was removed:
    // standalone sort now streams like every other command.
    chain.add_source()?;

    let last_idx = spec.stages.len() - 1;
    for (i, &stage) in spec.stages.iter().enumerate() {
        let position =
            if i == last_idx { StagePosition::Terminal } else { StagePosition::Intermediate };
        chain.add_stage(stage, position)?;
    }

    chain.add_sink()?;

    chain.build()
}
