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
    /// Whether this is the standalone sort chain (`[Stage::Sort]` only).
    ///
    /// Standalone sort no longer has a special build path — like every other
    /// chain it runs through `add_source` → sort stage → `add_sink`. This
    /// predicate is the spec classifier for the sole-`[Sort]` layout: `add_sort`
    /// calls it (as `is_standalone_sort`) to honour `--sort::max-memory=auto`
    /// only when sort owns the whole memory budget, and `validate` uses it.
    #[must_use]
    pub fn is_sort_terminal(&self) -> bool {
        matches!(self.stages.as_slice(), [Stage::Sort])
    }
}
