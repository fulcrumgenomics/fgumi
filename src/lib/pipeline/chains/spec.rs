//! Declarative chain specification — input to [`super::build_for`].

use crate::commands::common::{
    BamIoOptions, CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
};
use crate::pipeline::chains::{SinkSpec, SourceSpec, Stage, StageOptionsBag};
use fgumi_bam_io::ReadStreams;

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
    /// Concurrent positional-read policy for a seekable file source. `Fixed(1)`
    /// (the default every command but `sort` uses) is the plain
    /// sequential/async reader; `sort` sets this from its `--read-streams` flag
    /// to raise the device read queue depth (see
    /// [`fgumi_bam_io::scatter_reader`]).
    pub read_streams: ReadStreams,
    /// Whether the BAM source's BGZF decode verifies each block's CRC32. Set
    /// from the command's `--check-crc`/`--no-check-crc` policy (via
    /// [`crate::commands::common::BamIoOptions::effective_check_crc`]) so the
    /// chain reproduces the non-chain path's CRC behavior. Inert for the SAM
    /// source (no BGZF) and for the FASTQ source (which carries its own policy).
    pub verify_crc: bool,
    /// For `@PG` line injection into the output header.
    pub command_line: String,
}

/// The CLI option blocks every BAM-in/BAM-out single-stage command carries.
///
/// Exists so [`ChainSpec::single_stage`] can take them as one argument instead of
/// six positional ones; the fields are borrowed straight off the command struct.
pub struct SingleStageContext<'a> {
    pub io: &'a BamIoOptions,
    pub threading: &'a ThreadingOptions,
    pub compression: &'a CompressionOptions,
    pub scheduler: &'a SchedulerOptions,
    pub queue_memory: &'a QueueMemoryOptions,
    /// For `@PG` line injection into the output header.
    pub command_line: &'a str,
}

impl ChainSpec {
    /// Build the spec for a command that runs exactly one stage, BAM in to BAM out.
    ///
    /// `simplex`, `duplex`, and `codec` are all this shape, differing only in their
    /// [`Stage`] and which [`StageOptionsBag`] slot they fill. Constructing the
    /// full `ChainSpec` here rather than in each `execute` means a new field
    /// is wired once: previously each command hand-built the struct, so a
    /// field added to one and missed in the others left that command silently
    /// running on a default its siblings did not use.
    #[must_use]
    pub fn single_stage(
        stage: Stage,
        stage_opts: StageOptionsBag,
        ctx: &SingleStageContext<'_>,
    ) -> Self {
        Self {
            stages: vec![stage],
            source: SourceSpec::Bam(ctx.io.input.clone()),
            sink: SinkSpec::Bam(ctx.io.output.clone()),
            stage_opts,
            threading: ctx.threading.clone(),
            compression: ctx.compression.clone(),
            scheduler: ctx.scheduler.clone(),
            queue_memory: ctx.queue_memory.clone(),
            async_reader: ctx.io.async_reader,
            // Non-sort single-stage commands do not expose a read-stream knob;
            // keep the plain sequential/async reader. Only sort sets this.
            read_streams: ReadStreams::Fixed(1),
            verify_crc: ctx.io.effective_check_crc(),
            command_line: ctx.command_line.to_string(),
        }
    }

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
