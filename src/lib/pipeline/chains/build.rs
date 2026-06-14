//! [`build_for`] — the single chain-construction entry point.
//!
//! T3b.3 converts this from a match-based dispatch (10 single-stage arms +
//! catch-all bail) to a stage-by-stage loop. `build_for` now:
//!
//!   1. Validates the spec (progression, options presence, cross-stage constraints).
//!   2. Constructs a [`ChainBuilder`] from the spec.
//!   3. Calls `add_source` (skipped for the Sort-terminal special case — `SortBamFile`
//!      reads its own input file, so no pipeline source step is needed).
//!   4. Walks `spec.stages`, calling `add_stage(stage, position)` for each,
//!      where `position` is `Terminal` for the last stage and `Intermediate` for all
//!      earlier ones.
//!   5. Calls `add_sink` (skipped for Sort-terminal, same reason as step 3).
//!   6. Calls `build()` to produce the [`BuiltPipeline`].
//!
//! ## Sort-terminal special case
//!
//! When `spec.stages == [Stage::Sort]`, `SortBamFile` is a self-contained
//! `Exclusive` step that reads from and writes to files directly — it bypasses
//! both the normal source preamble and the `BgzfCompress → WriteBgzfFile` sink.
//! `add_sort(Terminal)` registers `SortBamFile` via `PipelineBuilder::append_source`
//! (legal for `Input = ()` steps) and does NOT call `add_source` / `add_sink`.
//! The calling code in `build_for` detects this case (`stages == [Sort]`) and
//! skips `add_source` + `add_sink` accordingly.
//!
//! For all other stage combinations (including Sort-intermediate, i.e. Sort followed
//! by Group, Simplex, Duplex, or Codec), `add_source` and `add_sink` ARE called:
//! the intermediate sort path uses `ParseBamRecords` (not `SortBamFile`) and feeds
//! records through the normal source → sort → downstream → sink topology.

use anyhow::{Result, bail};

use crate::pipeline::chains::{
    BuiltPipeline, ChainSpec, Stage,
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
/// chain.add_source()?                     // skipped for Sort-terminal
/// for (i, stage) in stages.iter().enumerate() {
///     let pos = if i == last { Terminal } else { Intermediate };
///     chain.add_stage(*stage, pos)?;
/// }
/// chain.add_sink()?                       // skipped for Sort-terminal
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

    // Sort-terminal special case: SortBamFile is self-contained (reads/writes
    // files directly). Skip add_source and add_sink.
    let sort_terminal = matches!(spec.stages.as_slice(), [Stage::Sort]);

    let mut chain = ChainBuilder::new(&spec)?;

    // NOTE: `--threads 1` whole-chain fusion is NOT wired here. It is applied
    // generically by the runtime (`Pipeline::run` → `run_fused_single_thread`)
    // to ANY linear source→sink chain at one worker — see
    // `pipeline::core::runtime::fused`. `build_for` therefore constructs the
    // same staged chain regardless of thread count; the runtime fuses it.

    if !sort_terminal {
        chain.add_source()?;
    }

    let last_idx = spec.stages.len() - 1;
    for (i, &stage) in spec.stages.iter().enumerate() {
        let position =
            if i == last_idx { StagePosition::Terminal } else { StagePosition::Intermediate };
        chain.add_stage(stage, position)?;
    }

    if !sort_terminal {
        chain.add_sink()?;
    }

    chain.build()
}
