//! Declarative chain specification — input to [`super::build_for`].

use crate::commands::common::{
    CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
};
use crate::pipeline::chains::{SinkSpec, SourceSpec, Stage, StageOptionsBag};

/// Declarative chain specification. Both standalone commands and
/// `runall` construct one of these and pass it to
/// [`crate::pipeline::chains::build_for`].
///
/// Construction is mechanical: parse CLI → fill in the relevant
/// `stage_opts` slots → set `source`/`sink`/`stages` based on input
/// shape and the requested operation → pass to `build_for`. Validators
/// inside `build_for` enforce ordering, mutual exclusions, and
/// stage-options presence.
pub struct ChainSpec {
    /// Ordered list of stages — must be a valid progression
    /// (validated by
    /// [`crate::pipeline::chains::validate::validate_stage_progression`]).
    pub stages: Vec<Stage>,
    pub source: SourceSpec,
    pub sink: SinkSpec,
    pub stage_opts: StageOptionsBag,
    pub threading: ThreadingOptions,
    pub compression: CompressionOptions,
    pub scheduler: SchedulerOptions,
    pub queue_memory: QueueMemoryOptions,
    /// When true, the BAM/SAM source is opened with a userspace async
    /// prefetch reader (`--async-reader`), overlapping disk I/O with compute.
    pub async_reader: bool,
    /// For `@PG` line injection into the output header.
    pub command_line: String,
}

impl ChainSpec {
    /// Whether this is the standalone sort-terminal chain (`[Stage::Sort]`).
    ///
    /// This layout is special: `build_for` skips `add_source`/`add_sink` and the
    /// terminal `SortBamFile` step opens and reads the input itself. Both the
    /// source/sink skip in `build_for` and the header-open skip in
    /// `ChainBuilder::new` MUST agree on this predicate — if they drift, a
    /// sort-terminal chain could open its source twice (draining a stdin pipe
    /// before the sort engine reads it). Single-source the test here so the two
    /// call sites cannot disagree.
    #[must_use]
    pub fn is_sort_terminal(&self) -> bool {
        matches!(self.stages.as_slice(), [Stage::Sort])
    }
}
