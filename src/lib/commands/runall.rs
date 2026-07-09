//! `runall` — fused multi-stage runall command.
//!
//! Runs a contiguous slice of the UMI pipeline in a single fused, in-memory
//! pipeline with no intermediate files. The slice spans
//! `extract → correct → align → zipper → sort → group → consensus → filter`;
//! the active stages are selected by the required `--start-from` /
//! `--stop-after` flags (see [`RunAllStage`]), and the terminal consensus
//! algorithm is chosen by `--consensus {simplex,duplex,codec}`:
//!
//!   * `--start-from extract` consumes FASTQ input; the other BAM-input
//!     starts (`correct`/`align`/`zipper`/`sort`/`group`/`consensus`/`filter`)
//!     read a BAM.
//!   * `--start-from sort` includes the in-pipeline `Sort` step (raw
//!     unsorted BAM input).
//!   * `--consensus {simplex,duplex,codec}` selects the consensus caller when
//!     the chain reaches the consensus stage.
//!
//! The chain mirrors the equivalent sequence of standalone commands
//! (`fgumi sort` → `fgumi group` → `fgumi {simplex,duplex,codec}`, plus any
//! `extract`/`correct`/`align`/`zipper`/`filter` stages in range), but skips
//! the serialize-to-BAM-bytes / decompress-from-BGZF / parse-records round
//! trip between stages.
//!
//! The diagram below shows the common `--start-from sort` slice; other starts
//! prepend (extract/correct/align/zipper) or append (filter) the corresponding
//! steps:
//!
//! ```text
//! read_bam → BgzfDecompress → FindBamBoundaries → [Sort steps?] →
//!     DecodeRecords → GroupByPosition → ProcessGroups → MiAssign →
//!     TemplatesToMiGroups (rewrite MI tags + group by MI) →
//!     {Simplex|Duplex|Codec}Consensus → BgzfCompress → WriteBgzfFile
//! ```
//!
//! Goals:
//!
//!   * **Parity**: byte-for-byte equal to running standalone
//!     `fgumi sort` (when `--start-from sort`) then `fgumi group`
//!     then `fgumi {simplex,duplex,codec}` on the same input.
//!     Differences in MI numbering would break parity so MI
//!     assignment uses the exact same cumulative counter as group.
//!   * **No intermediate file**: data flows in memory between every
//!     stage. The temp BAMs that a sequential run would produce are
//!     elided.
//!
//! Anti-goals:
//!
//!   * Per-position-group histograms (`--metrics`, `--family-size-histogram`,
//!     `--grouping-metrics`) are not emitted by runall under any
//!     `--start-from`/`--stop-after` combination. The new pipeline
//!     framework does not yet wire histogram state back through a
//!     side channel; runall's group stage construction sets all three
//!     histogram fields to `None`. Users who need these files should
//!     invoke standalone `fgumi group` directly.

use std::path::PathBuf;

use anyhow::{Result, bail};
use clap::Parser;
use log::{info, warn};

use crate::assigner::Strategy;
use crate::commands::command::Command;
use crate::commands::common::{
    BamIoOptions, CompressionOptions, QueueMemoryOptions, ReadGroupOptions, RejectsOptions,
    SchedulerOptions, StatsOptions, ThreadingOptions,
};
use crate::logging::OperationTimer;

/// Consensus mode selector for `runall`. Controls which consensus
/// caller runs after the fused group + MI-grouping stages.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, clap::ValueEnum)]
pub enum RunAllMode {
    /// Single-strand consensus via `VanillaUmiConsensusCaller`. Pair
    /// with `--strategy adjacency` (or `identity`/`edit`).
    Simplex,
    /// CODEC duplex consensus via `CodecConsensusCaller`. Pair with
    /// `--strategy adjacency` (or `identity`/`edit`).
    Codec,
    /// Two-strand duplex consensus via `DuplexConsensusCaller`. Requires
    /// `--strategy paired` so MIs carry `/A`/`/B` suffixes.
    Duplex,
}

impl std::fmt::Display for RunAllMode {
    /// Lower-case the variant name so error messages match the CLI
    /// flag values the user typed (`--consensus simplex`, not
    /// `RunAllMode::Simplex`).
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = match self {
            Self::Simplex => "simplex",
            Self::Codec => "codec",
            Self::Duplex => "duplex",
        };
        f.write_str(name)
    }
}

/// Pipeline stage marker for the `--start-from` / `--stop-after`
/// runall stage flags. Stages are linearly ordered:
/// `Extract < Correct < AlignAndMerge < Zipper < Sort < Group <
/// Consensus < Filter` (see [`RunAllStage::ord`]).
///
/// **Stage semantics — "the input is positioned AT this stage's
/// entry":**
///
/// * `Sort` — raw unsorted BAM, pre-sort. `--start-from sort`
///   includes the sort step. `--stop-after sort` writes the
///   sorted BAM and stops (apples-to-apples replacement for
///   standalone `fgumi sort`).
/// * `Group` — BAM is sorted by template-coordinate (output of
///   `fgumi sort`). `--start-from group` skips the sort step.
///   `--stop-after group` writes the grouped BAM (output of
///   `fgumi group`) and stops.
/// * `Consensus` — BAM is sorted and grouped (MI tags assigned,
///   output of `fgumi group`). `--start-from consensus` skips sort
///   and group. `--stop-after consensus` runs the consensus caller
///   chosen by `--consensus {simplex,duplex,codec}` and writes the
///   consensus BAM (output of `fgumi simplex`/`duplex`/`codec`). The
///   algorithm is selected by `--consensus`, not by a stage variant.
///
/// `AlignAndMerge` is a valid `--start-from` but not a valid
/// `--stop-after` (see [`StopAfter`]): raw aligner output before the
/// zipper-merge has lost every original tag.
///
/// **Case-insensitive parsing**: `sort`, `Sort`, `SORT` all match
/// because the `--start-from` / `--stop-after` flags use
/// `#[arg(value_enum, ignore_case = true)]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, clap::ValueEnum)]
#[clap(rename_all = "kebab-case")]
pub enum RunAllStage {
    /// Extract UMIs from FASTQ files and produce an unmapped BAM.
    /// `--start-from extract` accepts FASTQ files via
    /// `--extract::inputs` and read structures via
    /// `--extract::read-structures`. `--start-from extract --stop-after
    /// extract` runs the extract-only chain and writes the unmapped BAM
    /// (byte-identical to standalone `fgumi extract`). Otherwise Extract
    /// feeds into Correct (when `--correct::umi-files` is set) or Align.
    Extract,
    /// Correct UMIs in an unmapped BAM against a fixed whitelist of
    /// expected sequences. `--start-from correct` accepts an
    /// unmapped BAM via `--input` and requires at least one of
    /// `--correct::umis` / `--correct::umi-files`. `--stop-after
    /// correct` writes the corrected unmapped BAM and exits.
    Correct,
    /// Align an unmapped BAM with UMIs through an external aligner
    /// subprocess (bwa-mem3 / bwa / user-supplied command), then merge
    /// the resulting alignments with the original unmapped reads via
    /// the zipper-merge logic — all streamed end-to-end. `--start-from
    /// align` accepts an unmapped BAM via `--input` and requires
    /// `--ref` (the aligner reference) plus one of
    /// `--aligner::preset` or `--aligner::command`. Earliest stop
    /// point reachable is `--stop-after zipper` (raw aligner output
    /// without merge is not a supported stop point because it loses
    /// every original tag).
    ///
    /// CLI value is `align` (short form for "align-and-merge" — the
    /// stage does both, but the user-facing name is the verb).
    #[clap(name = "align")]
    AlignAndMerge,
    /// Zip an unmapped BAM with an aligned BAM, transferring tags.
    /// `--start-from zipper` accepts a queryname-sorted unmapped BAM
    /// (via `--unmapped`) and a queryname-sorted mapped BAM (via
    /// `--input`); `--stop-after zipper` writes the merged BAM.
    /// Requires `--reference` (the same FASTA + `.dict` the standalone
    /// `fgumi zipper` requires).
    Zipper,
    /// Raw unsorted BAM. Includes the sort step when used as
    /// `--start-from sort`; writes sorted BAM and stops when used
    /// as `--stop-after sort`.
    Sort,
    /// Sorted BAM (template-coordinate). Includes group + MI tag
    /// assignment when used as `--start-from group`; writes
    /// grouped BAM and stops when used as `--stop-after group`.
    Group,
    /// Consensus calling. The algorithm (simplex / duplex / CODEC) is
    /// selected by `--consensus`, not by this stage. `--start-from
    /// consensus` accepts a grouped (MI-tagged) BAM; `--stop-after
    /// consensus` runs the chosen consensus caller and writes the
    /// consensus BAM (byte-identical to `fgumi simplex` / `duplex` /
    /// `codec`). `--stop-after filter` fuses consensus → filter into one
    /// in-memory pipeline (no intermediate consensus BAM).
    Consensus,
    /// Post-consensus filtering. `--start-from filter` accepts a
    /// consensus BAM (output of `fgumi simplex`/`duplex`/`codec`) via
    /// `--input`; `--stop-after filter` runs the filter stage and writes
    /// the filtered BAM (byte-identical to standalone `fgumi filter`).
    /// Tunable via `--filter::*`. Filter is the last stage in the linear
    /// order — no stage follows it. It runs either as the filter self-pair
    /// (`--start-from filter --stop-after filter`, input = a consensus
    /// BAM) or fused after an upstream consensus caller in the same
    /// pipeline (`… → consensus → filter`), where the consensus output is
    /// decoded in-memory and streamed straight into filter.
    Filter,
}

/// The subset of [`RunAllStage`] values that are valid as `--stop-after`.
///
/// `align` (`AlignAndMerge`) is a valid `--start-from` but **not** a valid
/// `--stop-after`: raw aligner output before the zipper-merge has lost
/// every original tag (RX, QX, RG, …), so the earliest post-align stop
/// point is `zipper`. Giving `--stop-after` its own value-enum makes
/// `--help` list only the genuinely-valid stops and rejects
/// `--stop-after align` at the clap parse layer with a clear
/// "invalid value 'align'" error — instead of parsing it and failing
/// later at runtime validation.
///
/// Convert to [`RunAllStage`] with `RunAllStage::from(stop)` /
/// `stop.into()`; the two enums are kept in sync by
/// `stop_after_values_are_runallstage_minus_align`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
#[clap(rename_all = "kebab-case")]
pub enum StopAfter {
    /// `--stop-after extract` — extract-only (FASTQ → unmapped BAM).
    Extract,
    /// `--stop-after correct` — write the corrected unmapped BAM.
    Correct,
    /// `--stop-after zipper` — merged-but-unsorted BAM (earliest
    /// post-align stop).
    Zipper,
    /// `--stop-after sort` — template-coordinate sorted BAM.
    Sort,
    /// `--stop-after group` — grouped (MI-tagged) BAM.
    Group,
    /// `--stop-after consensus` — consensus BAM (algorithm chosen by
    /// `--consensus`); same output as `fgumi simplex`/`duplex`/`codec`.
    Consensus,
    /// `--stop-after filter` — post-consensus filtered BAM (same output
    /// as `fgumi filter`). Tunable via `--filter::*`.
    Filter,
}

impl From<StopAfter> for RunAllStage {
    fn from(stop: StopAfter) -> Self {
        match stop {
            StopAfter::Extract => Self::Extract,
            StopAfter::Correct => Self::Correct,
            StopAfter::Zipper => Self::Zipper,
            StopAfter::Sort => Self::Sort,
            StopAfter::Group => Self::Group,
            StopAfter::Consensus => Self::Consensus,
            StopAfter::Filter => Self::Filter,
        }
    }
}

impl RunAllStage {
    /// Position in the linear stage order. `Extract = 0`,
    /// `Correct = 1`, `AlignAndMerge = 2`, `Zipper = 3`, `Sort = 4`,
    /// `Group = 5`, `Consensus = 6`, `Filter = 7`. Used to validate that
    /// `start_from <= stop_after` and that the chain is buildable.
    #[must_use]
    pub fn ord(self) -> usize {
        match self {
            Self::Extract => 0,
            Self::Correct => 1,
            Self::AlignAndMerge => 2,
            Self::Zipper => 3,
            Self::Sort => 4,
            Self::Group => 5,
            Self::Consensus => 6,
            Self::Filter => 7,
        }
    }

    /// `true` if this stage is the consensus stage. Used by the Group
    /// step-derivation exception (a consensus-start input is already
    /// MI-tagged, so no explicit Group stage is added).
    #[must_use]
    pub fn is_consensus(self) -> bool {
        matches!(self, Self::Consensus)
    }

    /// `true` if this stage is the extract stage.
    #[must_use]
    pub fn is_extract(self) -> bool {
        matches!(self, Self::Extract)
    }

    /// Retained as a no-op shim during the consensus-mode split; the
    /// consensus algorithm now comes from `--consensus`, not the stage.
    /// Validate that `self` (as `--start-from`) and `stop_after` form
    /// a buildable pipeline range.
    ///
    /// # Errors
    ///
    /// Returns `Err` if `self.ord() > stop_after.ord()` — you can't start
    /// later than you stop. With consensus collapsed to a single stage
    /// (the algorithm comes from `--consensus`), the linear order is
    /// strictly ordinal and no consensus-terminal special case is needed.
    pub fn validate_with(self, stop_after: Self) -> Result<()> {
        if self.ord() > stop_after.ord() {
            bail!(
                "--start-from {self} comes after --stop-after {stop_after}; \
                 stop-after must be at or after start-from in the pipeline order \
                 (extract < correct < align < zipper < sort < group < consensus < filter)"
            );
        }
        Ok(())
    }
}

impl std::fmt::Display for RunAllStage {
    /// Lower-case the variant name so error messages match the CLI
    /// flag values the user typed (`--start-from sort`, not
    /// `--start-from Sort`).
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = match self {
            Self::Extract => "extract",
            Self::Correct => "correct",
            Self::AlignAndMerge => "align",
            Self::Zipper => "zipper",
            Self::Sort => "sort",
            Self::Group => "group",
            Self::Consensus => "consensus",
            Self::Filter => "filter",
        };
        f.write_str(name)
    }
}

/// Fused runall command. Stages are selected by `--start-from` and
/// `--stop-after`; the accepted starts are
/// `{Extract, Correct, AlignAndMerge, Zipper, Sort, Group, Consensus,
/// Filter}` (the consensus algorithm is chosen by `--consensus`, not by a
/// per-algorithm start value — `simplex`/`duplex`/`codec` are NOT accepted
/// start values) and the accepted stops are
/// `{Extract, Correct, Zipper, Sort, Group, Consensus, Filter}`
/// (no `--stop-after align` — raw aligner output without zipper-merge
/// would lose every original tag).
#[derive(Debug, Parser)]
#[command(
    name = "runall",
    about = "\x1b[38;5;166m[UTILITIES]\x1b[0m      \x1b[36mFused multi-stage pipeline (correct + align + zipper + sort + group + consensus, no intermediate BAM)\x1b[0m",
    long_about = "\
Fuses extract, UMI correction, alignment + zipper-merge, sort, group, and \
consensus calling (optionally followed by filter) into a single in-memory \
pipeline. Select the slice to run with `--start-from` and `--stop-after` (both \
required); the stages run in the fixed order extract < correct < align < zipper \
< sort < group < consensus < filter. Records flow directly between stages in \
memory, so every intermediate BAM a sequential run would write is elided, and \
the output is byte-for-byte identical to running the equivalent standalone \
commands in turn.

`--start-from` declares the state of the input and selects the upstream work:

* `extract` reads FASTQ via `--extract::inputs` / `--extract::read-structures` \
  (no `--input`).
* `correct` reads an unmapped BAM via `--input` and corrects UMIs against \
  `--correct::umis` / `--correct::umi-files`.
* `align` reads an unmapped BAM via `--input`, runs the aligner subprocess \
  (`--aligner::preset` or `--aligner::command`) against `--ref`, and \
  zipper-merges the result.
* `zipper` reads a queryname-sorted mapped BAM via `--input` plus the matching \
  unmapped BAM via `--unmapped` and the reference (`--ref` + `.dict`).
* `sort` reads a raw unsorted aligned BAM via `--input`.
* `group` reads a template-coordinate-sorted BAM via `--input`.
* `consensus` reads a grouped (MI-tagged) BAM via `--input`.
* `filter` reads a consensus BAM via `--input`.

The consensus algorithm is chosen by `--consensus {simplex,duplex,codec}`, not \
by a stage name; `--stop-after consensus` produces the same BAM as the matching \
`fgumi {simplex,duplex,codec}` command. `--stop-after filter` fuses a \
post-consensus `fgumi filter` pass (tunable via `--filter::*`) into the same \
pipeline, so no intermediate consensus BAM is written. Tune any stage with its \
prefixed flags: `--sort::max-memory`, `--group::strategy`, `--duplex::min-reads`, \
and so on.

`--stop-after align` is not allowed: raw aligner output before the zipper-merge \
has lost every original tag, so the earliest stop after `--start-from align` is \
`zipper`. The validator checks the requested stage range but cannot verify that \
the on-disk input actually matches `--start-from` — make sure it does."
)]
#[allow(clippy::struct_excessive_bools)]
pub struct RunAll {
    /// Input BAM/SAM file. Optional because runall's input contract is
    /// stage-dependent, unlike the standard single-required-BAM commands
    /// that flatten [`BamIoOptions`]: `--start-from extract` reads FASTQ
    /// from `--extract::inputs` (no `--input`); `--start-from zipper`
    /// pairs `--input` (mapped) with `--unmapped`; the BAM-source start
    /// stages (`correct`/`align`/`sort`/`group`/consensus) require
    /// `--input`. Presence is validated per start-stage in
    /// [`RunAll::derive_source_spec`].
    #[arg(short = 'i', long = "input")]
    pub input: Option<PathBuf>,

    /// Output BAM file.
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,

    /// Wrap the input in a userspace async prefetch reader: a background
    /// thread reads ahead so disk I/O overlaps decompression/compute.
    /// Defaults to off.
    #[arg(long = "async-reader", default_value_t = false, hide = true)]
    pub async_reader: bool,

    /// Optional output for rejected reads.
    #[command(flatten)]
    pub rejects_opts: RejectsOptions,

    /// Optional output for statistics.
    #[command(flatten)]
    pub stats_opts: StatsOptions,

    /// Read group and read name prefix options.
    #[command(flatten)]
    pub read_group: ReadGroupOptions,

    /// Threading options.
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Scheduler and pipeline statistics options.
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,

    /// Methylation-aware consensus calling mode (requires --ref).
    #[arg(long = "methylation-mode", value_enum)]
    pub methylation_mode: Option<crate::commands::common::MethylationModeArg>,

    /// Path to the reference FASTA file. Required when `--methylation-mode`
    /// is set (consensus stages), when `--start-from zipper` is used (zipper
    /// requires the FASTA + `.dict` to build the output BAM header from the
    /// reference dictionary), and whenever the derived stage chain includes
    /// `Stage::Align` (the aligner reference; see `validate_align_and_merge`) —
    /// not only `--start-from align`, but also fused runs that reach align from
    /// an upstream start (e.g. `--start-from extract`/`--start-from correct`
    /// with a `--stop-after` past zipper). For those align-bearing upstream
    /// starts the reference is consumed by the align stage and is accepted
    /// without `--methylation-mode` (the `ref_is_aligner_only` exemption in
    /// `execute`).
    #[arg(long = "ref")]
    pub reference: Option<PathBuf>,

    // ───────── extract-side options (used only when --start-from=extract) ─────────
    /// Per-stage extract tuning, exposed as `--extract::inputs`,
    /// `--extract::read-structures`, `--extract::sample`,
    /// `--extract::library`, etc. via the `MultiExtractRunallOptions`
    /// companion struct (generated by `#[multi_options]` on
    /// `ExtractRunallOptions` in `commands::extract`).
    /// Ignored when `--start-from` is not `extract`.
    #[command(flatten)]
    pub extract_opts: crate::commands::extract::MultiExtractRunallOptions,

    // ───────── correct-side options (used only when --start-from <= correct) ─────────
    /// Per-stage UMI-correction tuning, exposed as `--correct::umis`,
    /// `--correct::umi-files`, `--correct::max-mismatches`,
    /// `--correct::min-distance`, `--correct::revcomp`, etc. via the
    /// `MultiCorrectOptions` companion struct (generated by
    /// `#[multi_options]` on `CorrectOptions` in `commands::correct`).
    /// `--correct::min-distance` is required when `--start-from
    /// correct` is selected (no default); at least one of
    /// `--correct::umis` / `--correct::umi-files` is required.
    /// Ignored when `--start-from` is not `correct`.
    #[command(flatten)]
    pub correct_opts: crate::commands::correct::MultiCorrectOptions,

    // ───────── aligner-side options (used whenever the chain includes Stage::Align) ─────────
    /// Per-stage aligner tuning, exposed as `--aligner::preset`,
    /// `--aligner::command`, `--aligner::threads`, `--aligner::chunk-size`
    /// via the `MultiAlignerOptions` companion struct (generated by
    /// `#[multi_options]` on `AlignerOptions` in `crate::aligner`).
    /// Used whenever the derived stage chain includes `Stage::Align` — not
    /// only on `--start-from=align`, but also on fused runs that reach align
    /// from an upstream start (e.g. `--start-from extract` or
    /// `--start-from correct` with a `--stop-after` past zipper). Ignored only
    /// when the chain has no align stage.
    #[command(flatten)]
    pub aligner_opts: crate::aligner::MultiAlignerOptions,

    /// Override path for the aligner binary (preset mode only). When
    /// unset, the preset's binary (`bwa-mem3` / `bwa`) is found via
    /// `which::which()` on `PATH`. Rejected with a clear error if
    /// `--aligner::command` is used (command mode owns its own
    /// binary). Used whenever the chain includes `Stage::Align` (see
    /// `--aligner::*` above); ignored only when the chain has no align stage.
    ///
    /// Note: this flag is at the runall top level (`--aligner-bin`),
    /// NOT inside the `--aligner::*` family. The `--aligner::*` flags
    /// are mode-orthogonal tuning knobs (chunk-size, threads); the
    /// binary identity is mode-shaped — command mode owns it inside
    /// `--aligner::command "..."`, preset mode owns it here. Putting
    /// the override in the `--aligner::*` family would suggest a
    /// non-existent `--aligner::bin` symmetry.
    #[arg(long = "aligner-bin")]
    pub aligner_bin: Option<PathBuf>,

    // ───────── zipper-side options (used only when --start-from=zipper) ─────────
    /// Path to the unmapped BAM (required when `--start-from zipper`).
    /// The standalone `fgumi zipper` reads this via `-u`/`--unmapped`;
    /// runall surfaces it as a top-level flag because the
    /// `BamIoOptions`-shared `--input` already names the mapped BAM.
    /// Ignored when `--start-from` is not `zipper`.
    #[arg(long = "unmapped")]
    pub unmapped: Option<PathBuf>,

    /// Per-stage zipper tuning, exposed as `--zipper::tags-to-remove`,
    /// `--zipper::buffer`, `--zipper::skip-tc-tags`, etc. via the
    /// `MultiZipperOptions` companion struct (generated by
    /// `#[multi_options]` on `ZipperOptions` in `commands::zipper`).
    /// Ignored when `--start-from` is not `zipper`.
    #[command(flatten)]
    pub zipper_opts: crate::commands::zipper::MultiZipperOptions,

    // ───────── group-side options ─────────
    /// Per-stage group tuning, exposed as `--group::strategy`,
    /// `--group::edits`, `--group::min-map-q`, etc. via the
    /// `MultiGroupOptions` companion struct (generated by
    /// `#[multi_options]` on `GroupOptions` in `commands::group`).
    /// `--group::strategy` is required when `--start-from {sort,group}`
    /// is selected (no default).
    #[command(flatten)]
    pub group_opts: crate::commands::group::MultiGroupOptions,

    // ───────── codec-specific options (used only when --stop-after=codec) ─────────
    /// Per-stage codec tuning, exposed as `--codec::min-duplex-length`,
    /// `--codec::single-strand-qual`, `--codec::outer-bases-qual`,
    /// `--codec::outer-bases-length`, `--codec::max-duplex-disagreements`,
    /// `--codec::max-duplex-disagreement-rate`. Ignored when
    /// `--stop-after` is not `codec`.
    #[command(flatten)]
    pub codec_opts: crate::commands::codec::MultiCodecOptions,

    // ───────── simplex-specific options (used only with --consensus simplex) ─────────
    /// Per-stage simplex tuning, exposed as `--simplex::error-rate-pre-umi`,
    /// `--simplex::min-reads`, `--simplex::max-reads`,
    /// `--simplex::min-input-base-quality`, etc. via the
    /// `MultiSimplexOptions` companion struct (generated by
    /// `#[multi_options]` on `SimplexOptions` in `commands::simplex`).
    /// `--simplex::min-reads` is required when `--consensus simplex` is
    /// selected (no default). Ignored unless the chain reaches the
    /// consensus stage with `--consensus simplex`.
    #[command(flatten)]
    pub simplex_opts: crate::commands::simplex::MultiSimplexOptions,

    // ───────── duplex-specific options (used only with --consensus duplex) ─────────
    /// Per-stage duplex tuning, exposed as `--duplex::error-rate-pre-umi`,
    /// `--duplex::min-reads`, `--duplex::max-reads-per-strand`,
    /// `--duplex::consensus-call-overlapping-bases`, etc. via the
    /// `MultiDuplexOptions` companion struct. Ignored unless the chain
    /// reaches the consensus stage with `--consensus duplex`.
    #[command(flatten)]
    pub duplex_opts: crate::commands::duplex::MultiDuplexOptions,

    // ───────── pipeline stage control (#33) ─────────
    /// Pipeline stage to start from: extract, correct, align, zipper, sort, group, consensus, or filter (required, case-insensitive).
    ///
    /// Declares the state of the input and selects the upstream work
    /// runall performs before the terminal `--stop-after` stage. See the
    /// command description for what each start stage reads (FASTQ for
    /// `extract`; `--input` plus, for `zipper`, `--unmapped`; a
    /// template-coordinate-sorted BAM for `group`; a grouped MI-tagged
    /// BAM for `consensus`; and so on).
    ///
    /// * `sort` — input is raw unsorted BAM; the chain runs the
    ///   in-pipeline sort step, then continues into group + consensus
    ///   unless `--stop-after` truncates it earlier (e.g. `--stop-after
    ///   sort` writes the sorted BAM and stops).
    /// * `group` — input is sorted by template-coordinate (output of
    ///   `fgumi sort`); the chain skips sort and runs group, then
    ///   continues into consensus unless `--stop-after` truncates it
    ///   earlier (e.g. `--stop-after group` writes the grouped BAM and
    ///   stops).
    /// * `consensus` — input is sorted and grouped (MI-tagged, output
    ///   of `fgumi group`); the chain runs the consensus caller selected
    ///   by `--consensus {simplex,duplex,codec}` (optionally followed by
    ///   `filter`). The algorithm is chosen by `--consensus`, NOT by the
    ///   start value — `simplex` / `duplex` / `codec` are not valid
    ///   `--start-from` values. See the consensus-start notes below.
    ///
    /// **Consensus-start chain:** `--start-from consensus` (with
    /// `--consensus {simplex,duplex,codec}`) runs the standalone consensus
    /// chain, which groups by the *existing* `MI` tag via `GroupByMi`. It does
    /// NOT add a position-grouping `Stage::Group` (no `GroupByPosition` /
    /// `ProcessGroups` / `MiAssign`) — doing so would double-group an
    /// already-grouped input (see [`derive_stages_for`]'s consensus-self-pair
    /// exception). On already-grouped input the output is record- and
    /// header-equivalent to the standalone consensus command (`fgumi simplex` /
    /// `fgumi duplex` / `fgumi codec`), ignoring `@PG` provenance — runall and the
    /// standalone command record different `@PG` chains, so the output is not
    /// byte-identical.
    ///
    /// **Input-state hazard:** the validator cannot check whether the on-disk
    /// BAM actually matches the declared `start_from`. Because
    /// `--start-from consensus` groups by the existing `MI` tag and adds no
    /// grouping stage, the input must already be grouped/MI-tagged — i.e. the
    /// output of a prior `group`. A merely *sorted-but-not-grouped* BAM (no `MI`
    /// tags), or a *raw-unsorted* BAM, is NOT grouped in-pipeline and will
    /// silently produce wrong/empty output rather than a clear error. (The
    /// consensus algorithm itself is chosen by `--consensus
    /// {simplex,duplex,codec}`, not by the start value.)
    ///
    /// Case-insensitive. Must satisfy `start_from.ord() <=
    /// stop_after.ord()` and the consensus-terminal rules
    /// (see [`RunAllStage::validate_with`]).
    #[arg(long = "start-from", value_enum, ignore_case = true)]
    pub start_from: RunAllStage,

    /// Pipeline stage to stop after: extract, correct, zipper, sort, group, consensus, or filter (required, case-insensitive).
    ///
    /// Determines the terminal stage of the chain. `--stop-after
    /// consensus` writes the consensus BAM produced by the caller chosen
    /// with `--consensus`; `--stop-after filter` fuses a post-consensus
    /// filter pass into the same pipeline.
    ///
    /// `align` is intentionally not a valid stop: raw aligner output
    /// before the zipper-merge has lost every original tag, so the
    /// earliest stop after `--start-from align` is `zipper`. Must satisfy
    /// `start_from.ord() <= stop_after.ord()` (`RunAllStage::validate_with`);
    /// use [`RunAll::stop_after_stage`] to read this as a [`RunAllStage`].
    #[arg(long = "stop-after", value_enum, ignore_case = true)]
    pub stop_after: StopAfter,

    /// Consensus algorithm to run at the consensus stage: `simplex`,
    /// `duplex`, or `codec`. Required whenever the chain reaches the
    /// consensus stage (`--stop-after consensus`, or `--start-from
    /// consensus`); ignored otherwise. `duplex` requires
    /// `--strategy paired`. Selects which standalone caller's output the
    /// run reproduces (`fgumi simplex` / `duplex` / `codec`).
    #[arg(long = "consensus", value_enum, ignore_case = true)]
    pub consensus_mode: Option<RunAllMode>,

    /// Per-stage sort tuning, exposed as `--sort::max-memory`,
    /// `--sort::tmp-dir`, etc. Consumed only when `--start-from sort`
    /// is selected; otherwise ignored. The full standalone option
    /// set of `fgumi sort` is mirrored here via the
    /// `#[multi_options]` macro — see `SortOptions` in
    /// `commands::sort` for the field list.
    #[command(flatten)]
    pub sort_opts: crate::commands::sort::MultiSortOptions,

    /// Per-stage filter tuning, exposed as `--filter::min-reads`,
    /// `--filter::max-read-error-rate`, `--filter::ref`, … Used when the
    /// chain reaches the filter stage (`--stop-after filter` or
    /// `--start-from filter`); otherwise ignored. The full standalone
    /// option set of `fgumi filter` is mirrored here via the
    /// `#[multi_options]` macro — see `FilterOptions` in
    /// `commands::filter` for the field list.
    #[command(flatten)]
    pub filter_opts: crate::commands::filter::MultiFilterOptions,
}

impl RunAll {
    /// Returns `true` if the user passed enough `--correct::*` flags to
    /// make Correct a load-bearing stage in an extract-fed chain.
    /// Concretely: at least one of `--correct::umi-files` or
    /// `--correct::umis` is non-empty. Used by `derive_stages` to
    /// decide whether to splice Correct between Extract and Align.
    #[must_use]
    pub(crate) fn extract_chain_wants_correct(&self) -> bool {
        self.correct_opts.has_umi_source()
    }

    /// Returns the consensus algorithm for this run, from `--consensus`.
    /// Errors (does not panic) if `--consensus` was omitted — it is required
    /// whenever the chain reaches the consensus stage. Only called on chains
    /// that actually reach consensus; `--stop-after sort`/`group`/`filter` are
    /// accepted stops that may not require it.
    pub(crate) fn require_consensus_mode(&self) -> Result<RunAllMode> {
        self.consensus_mode.ok_or_else(|| {
            anyhow::anyhow!(
                "--consensus <simplex|duplex|codec> is required when the chain \
                 reaches the consensus stage"
            )
        })
    }
}

impl Command for RunAll {
    fn execute(&self, command_line: &str) -> Result<()> {
        use crate::pipeline::chains::{ChainSpec, Stage};

        let timer = OperationTimer::new(&format!(
            "runall (--start-from {} --stop-after {})",
            self.start_from,
            self.stop_after_stage()
        ));

        // Validate the input BAM exists when present. `--start-from
        // extract` has no `--input` (FASTQ comes from `--extract::inputs`,
        // validated in `derive_source_spec`), so this is conditional.
        if let Some(input) = &self.input {
            BamIoOptions::new(input, &self.output).validate()?;
        }

        // ── Extract-start guards ─────────────────────────────────────────
        if self.start_from == RunAllStage::Extract {
            // `extract → extract` is the supported extract-only self-pair
            // (FASTQ → unmapped BAM); no guard needed.
            if self.stop_after_stage() == RunAllStage::Zipper {
                bail!(
                    "--start-from extract --stop-after zipper is incompatible \
                     (zipper requires PairedBams, but extract produces Fastqs)"
                );
            }
            if self.stop_after_stage() == RunAllStage::Correct
                && !self.extract_chain_wants_correct()
            {
                bail!(
                    "--start-from extract --stop-after correct requires \
                     --correct::umi-files or --correct::umis"
                );
            }
        }

        // stdin is supported: the input flows into `SourceSpec::Bam`/`Fastqs`
        // and is opened once via `InputSource::open` (BGZF/SAM auto-detected,
        // stdin-aware). runall reads the source exactly once per run, so a
        // streamed `-` works for every start stage.

        self.validate_stages()?;

        info!(
            "Starting runall (--start-from {} --stop-after {})",
            self.start_from,
            self.stop_after_stage()
        );
        if let Some(input) = &self.input {
            info!("Input: {}", input.display());
        }
        info!("Output: {}", self.output.display());

        // Derive the ordered stage list for this (start_from, stop_after) pair.
        let stages = self.derive_stages()?;

        // `--ref` pairs with `--methylation-mode` for the consensus stages, but
        // it is also legitimately consumed *without* methylation when the
        // derived chain includes `Stage::Align` (the aligner reference) or
        // starts at `zipper` (the FASTA + `.dict` for the output BAM header).
        // For any other chain `--ref` is dead, so reject it without
        // `--methylation-mode` rather than silently ignore it. Keyed on the
        // derived chain — not the start stage — so a fused `extract`/`correct →
        // … → align` run is exempted while a non-aligning chain (e.g.
        // `correct → correct`) is not. This must run after `derive_stages()`.
        let ref_is_consumed_without_methylation =
            self.start_from == RunAllStage::Zipper || stages.contains(&Stage::Align);
        if !ref_is_consumed_without_methylation
            && self.reference.is_some()
            && self.methylation_mode.is_none()
        {
            bail!("--ref requires --methylation-mode to be set");
        }

        // Symmetric guard for the *other* flag: in runall, `--methylation-mode`
        // is wired only into the simplex/duplex consensus stages (see the
        // `simplex_opts`/`duplex_opts` bag population below; codec and align
        // reject it with their own messages). On a chain that reaches no
        // consensus stage (e.g. `group → group`, `correct → sort`) it is dead —
        // silently ignored — so reject it rather than mislead. Align chains are
        // exempt here: `validate_align_and_merge` (run just below for any chain
        // containing `Stage::Align`) owns the align-specific EM-seq message,
        // which is more actionable than this generic one.
        let chain_reaches_consensus = stages.iter().any(|s| s.is_consensus());
        let chain_includes_align = stages.contains(&Stage::Align);
        // `--methylation-mode` on a simplex/duplex chain requires `--ref`: the
        // consensus stage consumes the FASTA (via `consensus_reference()`), and
        // the standalone contract (`--methylation-mode` doc, "requires --ref")
        // demands it. Without this guard a non-align `group → simplex` run with
        // `--methylation-mode` but no `--ref` slips through (the dead-consensus
        // guard below does not fire because the chain *does* reach consensus),
        // and `build_stage_options_bag()` then threads `methylation_mode =
        // Some(...)` with `reference = None` into the stage — deferring the
        // failure to a later, less actionable path. Codec is excluded (it
        // rejects `--methylation-mode` outright below), and align chains are
        // exempt (`validate_align_and_merge` owns the more specific message).
        let chain_uses_methylation_capable_consensus =
            stages.iter().any(|s| matches!(s, Stage::Simplex | Stage::Duplex));
        if self.methylation_mode.is_some()
            && chain_uses_methylation_capable_consensus
            && !chain_includes_align
            && self.reference.is_none()
        {
            bail!("--methylation-mode requires --ref to be set");
        }
        if self.methylation_mode.is_some() && !chain_reaches_consensus && !chain_includes_align {
            bail!(
                "--methylation-mode is only consumed by the consensus stage; \
                 it is dead on a runall chain that stops before consensus"
            );
        }

        // D12: log auto-inserted stages for extract-fed chains so users
        // understand the chain expansion rules.
        if self.start_from == RunAllStage::Extract {
            if stages.contains(&Stage::Correct) {
                info!(
                    "Including Stage::Correct between Extract and downstream \
                     (because --correct::umi-files or --correct::umis is set)"
                );
            }
            if stages.contains(&Stage::Align) {
                info!(
                    "Including Stage::Align between Extract and downstream \
                     (FASTQ source requires alignment)"
                );
            }
            if stages.contains(&Stage::Sort) && self.start_from != RunAllStage::Sort {
                info!(
                    "Including Stage::Sort between Align and Group \
                     (template-coordinate sort precedes Group)"
                );
            }
        }

        // Align-stage pre-flight: run `validate_align_and_merge` so that
        // aligner option conflicts (preset vs command, missing --ref,
        // --methylation-mode with AAM) are caught before the pipeline starts.
        // `add_align` in the chain builder independently validates aligner
        // options via `resolve()`, but `--methylation-mode` is a RunAll-level
        // field not stored in the StageOptionsBag, so we must gate it here.
        if stages.contains(&Stage::Align) {
            self.validate_align_and_merge()?;
        }

        // Cross-stage validation: when Group feeds a consensus stage,
        // the group strategy must be compatible with the consensus mode
        // (e.g. duplex requires `--strategy paired`).
        let has_group = stages.contains(&Stage::Group);
        let has_consensus = stages.iter().any(|s| s.is_consensus());
        if has_group && has_consensus {
            let consensus_mode = self.require_consensus_mode()?;
            let group_opts = self.group_opts.clone().validate()?;
            validate_strategy_for_mode(consensus_mode, group_opts.strategy)?;
            // Log group/consensus summary banner (mirrors old fused dispatcher).
            // Per-mode consensus tuning (min-reads, error rates, …) now comes
            // from `--simplex::*` / `--duplex::*` / `--codec::*`, so the
            // min-reads value is logged by the per-mode chain builder, not
            // here.
            info!("Strategy: {:?}, edits: {}", group_opts.strategy, group_opts.edits);
        }

        // Surface a warning if --rejects is set with --start-from correct
        // and a chained --stop-after. The fused chain uses
        // correct_step_kept_only (no UMI-rejects branch) so the rejects
        // file collects only downstream rejects (post-consensus / etc.),
        // not correct's UMI rejects. Users who need UMI rejects should
        // stage a `fgumi correct --rejects` step separately.
        //
        // Surfaced here on the build_for path so the same UX appears
        // for any fused chain that starts at correct.
        if matches!(self.start_from, RunAllStage::Extract | RunAllStage::Correct)
            && self.stop_after_stage() != RunAllStage::Correct
            && self.rejects_opts.rejects.is_some()
        {
            warn!(
                "--rejects with --start-from correct --stop-after {} discards correct's UMI rejects \
                 (the fused chain has no UMI-rejects branch). The rejects file will collect only \
                 downstream rejects. Use a separate `fgumi correct --rejects` step if you need \
                 UMI rejects captured.",
                self.stop_after_stage()
            );
        }

        // T3b.3.1: all chain shapes (including intermediate-sort) route
        // through `build_for`. `ChainBuilder::add_sort` handles both the
        // standalone file-to-file case (Sort is the sole stage) and the
        // streaming in-pipeline case (Sort follows or precedes other
        // stages); downstream stages accept all source variants for
        // their log-path lookup.

        // Construct the ChainSpec from the derived stages + T3b.2 helpers.
        let source = self.derive_source_spec()?;
        let sink = self.derive_sink_spec()?;
        let stage_opts = self.build_stage_options_bag(&stages)?;

        let spec = ChainSpec {
            async_reader: self.async_reader,
            stages,
            source,
            sink,
            stage_opts,
            threading: self.threading.clone(),
            compression: self.compression.clone(),
            scheduler: self.scheduler_opts.clone(),
            queue_memory: self.queue_memory.clone(),
            command_line: command_line.to_string(),
        };

        crate::pipeline::chains::build_for(spec)?.run()?;

        timer.log_completion(0);
        info!("runall completed successfully");
        Ok(())
    }
}

/// Validate that the `--strategy` chosen by the user is compatible
/// with the terminal consensus mode derived from `--stop-after`.
///
/// `Duplex` requires `Strategy::Paired` (MIs need `/A`/`/B` suffixes);
/// `Simplex` and `Codec` require any non-paired strategy. Error
/// messages name `--stop-after` (Display, lowercase) so the user
/// sees the flag they typed.
///
/// Free function so the error-message Display formatting is unit-testable
/// without constructing a full `RunAll`.
fn validate_strategy_for_mode(mode: RunAllMode, strategy: Strategy) -> Result<()> {
    match (mode, strategy) {
        (RunAllMode::Duplex, Strategy::Paired) => Ok(()),
        (RunAllMode::Duplex, _) => {
            bail!("--consensus duplex requires --strategy paired (so MIs carry /A and /B suffixes)")
        }
        (RunAllMode::Simplex | RunAllMode::Codec, Strategy::Paired) => {
            bail!("--consensus {mode} requires a non-paired strategy (identity/edit/adjacency)")
        }
        _ => Ok(()),
    }
}

/// Validate the `--start-from` / `--stop-after` pair against the
/// structural (linear-ordering) rules. Delegates to
/// [`RunAllStage::validate_with`]; consensus is no longer terminal-only
/// (it fuses into `consensus → filter`), so only the ordinal order is
/// enforced.
///
/// Historically this function also enforced "not yet implemented"
/// gates while #33 tasks 4-9 incrementally wired each combination;
/// those gates are all lifted now. The function is kept as the
/// runall-side entry point so future validation rules (e.g.
/// strategy/stage compatibility beyond what `validate_with`
/// expresses) have a clear home.
fn validate_stages_for(start_from: RunAllStage, stop_after: RunAllStage) -> Result<()> {
    start_from.validate_with(stop_after)?;
    // `--stop-after align-and-merge` would mean "stop after the
    // alignment but before zipper-merge". The design dropped that
    // stop point because raw aligner output without the zipper-merge
    // loses every original tag (RX, QX, RG, etc.) and isn't useful.
    // Earliest stop point reachable from `--start-from align-and-merge`
    // is `--stop-after zipper`. `validate_with`'s linear-order check
    // would have let this combo through (ord 0 → ord 0 is structurally
    // valid as a self-pair), so we reject it explicitly here.
    if stop_after == RunAllStage::AlignAndMerge {
        bail!(
            "--stop-after align is not supported; raw aligner output \
             without zipper-merge loses every original tag (RX, QX, RG, ...). \
             Use `--stop-after zipper` for merged-but-unsorted output."
        );
    }
    Ok(())
}

/// Derive the ordered [`Stage`] list for a given `(start_from, stop_after)`
/// pair.
///
/// Free function so the stage-derivation logic is unit-testable without
/// constructing a full [`RunAll`] instance — same pattern as
/// [`validate_stages_for`]. [`RunAll::derive_stages`] is the thin wrapper
/// that reads `self` and calls this.
///
/// # Stage derivation rules
///
/// 1. **Correct**: included when `start_from == Correct`.
///    When `stop_after > Correct`, the chain always chains through AAM
///    (`Correct → Align` is the only supported downstream path; a corrected
///    unmapped BAM must be aligned before it can be sorted or grouped).
///
/// 2. **Align**: included when `start_from <= AlignAndMerge` AND
///    `stop_after >= Zipper`. `Stage::Align` encapsulates both the aligner
///    subprocess AND the zipper-merge; there is no separate
///    `--stop-after align` stop point (the validator rejects it).
///    Included when:
///    - `start_from == AlignAndMerge` (explicit AAM start), OR
///    - `start_from == Correct` and `stop_after > Correct` (correct feeds
///      into AAM as the mandatory next step).
///
/// 3. **Zipper**: included when `start_from == Zipper` (standalone zipper
///    start). Mutually exclusive with `Stage::Align` — both produce a merged
///    BAM from different source shapes.
///
/// 4. **Sort**: included when `start_from <= Sort` OR when `Align` /
///    `Zipper` is in the chain AND `stop_after >= Sort`.
///    **Sort-forcing rule**: `Stage::Align` and `Stage::Zipper` both emit
///    queryname-sorted output. Downstream `Stage::Group` requires
///    template-coordinate order — Sort is therefore mandatory between any
///    Align/Zipper stage and Group, even when `--start-from` is not `sort`.
///
/// 5. **Group**: included when `stop_after >= Group` AND `start_from` is NOT
///    a consensus stage. **Consensus-self-pair exception**: when
///    `start_from.is_consensus()`, the input BAM is already MI-tagged. The
///    standalone consensus chain builders handle MI-tag grouping internally
///    via `GroupByMi` — adding an explicit `Stage::Group` here would
///    double-group an already-grouped input. For all other start stages
///    (correct, align, zipper, sort, group), the Group step is required to
///    produce MI-tagged output before the consensus caller sees it.
///
/// 6. **Simplex / Duplex / Codec**: one of the three is appended when
///    `stop_after` names that consensus stage (terminal, mutually exclusive).
///
/// # Errors
///
/// Returns `Err` for any `(start, stop)` pair that `validate_stages` would
/// reject (should not be reached if `validate_stages` was called first), or
/// if the Stage derivation reaches a branch that is a programming error (e.g.
/// a non-consensus `stop` that makes it through all prior early-returns).
fn derive_stages_for(
    start_from: RunAllStage,
    stop_after: RunAllStage,
    extract_wants_correct: bool,
    consensus_mode: Option<RunAllMode>,
) -> Result<Vec<crate::pipeline::chains::Stage>> {
    use crate::pipeline::chains::Stage;

    // `validate_stages` must have been called before this; the assert
    // catches misuse in tests or future refactors.
    debug_assert!(
        start_from.ord() <= stop_after.ord() && stop_after != RunAllStage::AlignAndMerge,
        "derive_stages_for called with invalid (start={start_from}, stop={stop_after}); \
         validate_stages must run first"
    );

    let mut stages: Vec<Stage> = Vec::with_capacity(7);

    // ── Step 0: Extract ──────────────────────────────────────────────────
    // Included when the chain starts at Extract.
    if start_from == RunAllStage::Extract {
        stages.push(Stage::Extract);
        // extract → extract: extract-only chain (FASTQ → unmapped BAM).
        // The terminal Extract stage serialises and writes the unmapped
        // BAM directly; there are no downstream stages.
        if stop_after == RunAllStage::Extract {
            return Ok(stages);
        }
    }

    // ── Filter self-pair ─────────────────────────────────────────────────
    // `--start-from filter --stop-after filter`: input is a consensus BAM;
    // run only the (terminal) filter stage. No group/consensus steps apply
    // (the standalone filter chain handles its own queryname grouping when
    // `--filter::filter-by-template` is set), so short-circuit here. The
    // validator guarantees `start_from == Filter` implies
    // `stop_after == Filter` (filter is the last stage in the linear order).
    if start_from == RunAllStage::Filter {
        stages.push(Stage::Filter);
        return Ok(stages);
    }

    // ── Step 1: Correct ──────────────────────────────────────────────────
    // Included when the chain starts at Correct. Also included when the
    // chain starts at Extract AND the user supplied --correct::umi-files
    // or --correct::umis (the extract_wants_correct flag).
    // For a Correct self-pair there are no downstream stages; return early.
    if start_from == RunAllStage::Correct
        || (start_from == RunAllStage::Extract && extract_wants_correct)
    {
        stages.push(Stage::Correct);
        if stop_after == RunAllStage::Correct {
            return Ok(stages);
        }
        // `Correct → stop > Correct` always chains through AAM.
        // Any stop beyond Correct (Zipper, Sort, Group, consensus)
        // requires align-and-merge — a corrected unmapped BAM cannot
        // be fed directly to Sort (no query-coordinate BAM).
    }

    // ── Step 2: Align (AAM = align + zipper-merge, fused) ────────────────
    // `Stage::Align` covers both the aligner subprocess and the
    // zipper-merge; there is no separate `--stop-after align` stop point.
    // Included when:
    //   (a) start_from == AlignAndMerge (explicit AAM start), OR
    //   (b) start_from == Correct and stop_after > Correct (correct feeds
    //       into AAM — see Step 1 above).
    //   (c) start_from == Extract and stop_after > Correct (extract feeds
    //       into Align, possibly through Correct first).
    let includes_align = start_from == RunAllStage::AlignAndMerge
        || (start_from == RunAllStage::Correct && stop_after != RunAllStage::Correct)
        || (start_from == RunAllStage::Extract && stop_after.ord() > RunAllStage::Correct.ord());
    if includes_align {
        stages.push(Stage::Align);
        // `--stop-after zipper` stops after the AAM step (which includes
        // the zipper-merge internally). Return early.
        if stop_after == RunAllStage::Zipper {
            return Ok(stages);
        }
    }

    // ── Step 3: Zipper (standalone zipper start) ─────────────────────────
    // Mutually exclusive with `Stage::Align`. Included only when
    // `start_from == Zipper` (the user supplies separate mapped + unmapped
    // BAMs and the zipper-merge runs as the first step).
    let includes_zipper = start_from == RunAllStage::Zipper;
    if includes_zipper {
        stages.push(Stage::Zipper);
        if stop_after == RunAllStage::Zipper {
            return Ok(stages);
        }
    }

    // ── Step 4: Sort ──────────────────────────────────────────────────────
    // Sort-forcing rule: Sort is ALWAYS included when Align or Zipper
    // precede anything past Sort, because both produce queryname-sorted
    // output and downstream Group requires template-coordinate order.
    // Sort is also included when start_from == Sort (explicit sort start).
    let sort_forced =
        (includes_align || includes_zipper) && stop_after.ord() >= RunAllStage::Sort.ord();
    let sort_explicit = start_from == RunAllStage::Sort;
    if sort_forced || sort_explicit {
        stages.push(Stage::Sort);
    }
    if stop_after == RunAllStage::Sort {
        return Ok(stages);
    }

    // ── Step 5: Group ─────────────────────────────────────────────────────
    // Included when stop_after >= Group AND start_from is NOT a consensus
    // stage. The consensus-self-pair exception: when start_from is already
    // a consensus stage (Simplex/Duplex/Codec), the input BAM is already
    // MI-tagged. The standalone consensus chain builder (build_simplex_chain
    // / build_duplex_chain / build_codec_chain) handles MI-tag grouping
    // internally via GroupByMi — it does NOT run GroupByPosition /
    // ProcessGroups / MiAssign. Adding an explicit Group stage here would
    // double-group an already-grouped input.
    //
    // For `--start-from group --stop-after {consensus}` or any upstream
    // start (sort, align, zipper, correct), Group IS required because the
    // input is not yet MI-tagged.
    let includes_group = !start_from.is_consensus()
        && (start_from == RunAllStage::Group || stop_after.ord() >= RunAllStage::Group.ord());
    if includes_group {
        stages.push(Stage::Group);
    }
    if stop_after == RunAllStage::Group {
        return Ok(stages);
    }

    // ── Step 6: Consensus ────────────────────────────────────────────────
    // Reached when stop_after is Consensus or Filter (all earlier stops
    // have returned, and the filter self-pair short-circuited above). The
    // chain stage is chosen by the `--consensus` algorithm. When
    // stop_after == Filter, consensus runs first, then Filter is appended
    // below (consensus → filter chain).
    debug_assert!(
        matches!(stop_after, RunAllStage::Consensus | RunAllStage::Filter),
        "derive_stages_for reached the consensus branch with unexpected \
         stop stage {stop_after}; validate_stages must run first"
    );
    let mode = consensus_mode.ok_or_else(|| {
        anyhow::anyhow!(
            "--consensus <simplex|duplex|codec> is required when --stop-after \
             reaches the consensus stage"
        )
    })?;
    match mode {
        RunAllMode::Simplex => stages.push(Stage::Simplex),
        RunAllMode::Duplex => stages.push(Stage::Duplex),
        RunAllMode::Codec => stages.push(Stage::Codec),
    }

    // ── Step 7: Filter (terminal) ────────────────────────────────────────
    // Appended after the consensus caller when `--stop-after filter` chains
    // consensus → filter. Filter is the last stage in the linear order, so
    // no stage follows it.
    if stop_after == RunAllStage::Filter {
        stages.push(Stage::Filter);
    }

    Ok(stages)
}

impl RunAll {
    /// The `--stop-after` value as a [`RunAllStage`]. The CLI field is the
    /// narrower [`StopAfter`] enum (which excludes `align`); this converts
    /// it back to the full stage enum the pipeline logic operates on.
    fn stop_after_stage(&self) -> RunAllStage {
        self.stop_after.into()
    }

    /// The input BAM path for a BAM-source start stage. Errors if
    /// `--input` is absent — only `--start-from extract` legitimately
    /// omits it (its FASTQ inputs come from `--extract::inputs`).
    fn require_input(&self) -> Result<PathBuf> {
        self.input.clone().ok_or_else(|| {
            anyhow::anyhow!(
                "--input is required with --start-from {} (the input BAM)",
                self.start_from
            )
        })
    }

    /// Reconstruct a [`BamIoOptions`] for delegating to a standalone
    /// command (the consensus fallback paths build a `Simplex`/`Duplex`/
    /// `Codec` and call its executor). Those stages are all BAM-source,
    /// so `--input` must be present.
    fn delegated_io(&self) -> Result<BamIoOptions> {
        Ok(BamIoOptions {
            input: self.require_input()?,
            output: self.output.clone(),
            async_reader: self.async_reader,
        })
    }

    /// [`BamIoOptions`] for a consensus stage's options bag.
    ///
    /// When `--start-from consensus`, the consensus stage **is** the chain
    /// source, so a real `--input` BAM is required (delegate to
    /// [`Self::delegated_io`]). For any pre-consensus start (extract, correct,
    /// align, zipper, sort, group), consensus runs **mid-chain**: the chain
    /// reads its data from `spec.source` (set by `derive_source_spec`) and the
    /// consensus `add_*` builders never touch the per-stage `.io.input`. In
    /// that case requiring `--input` is spurious (it would reject e.g.
    /// `--start-from extract --stop-after filter`, where input comes from
    /// `--extract::inputs`), so we hand back a placeholder IO whose `input` is
    /// never read.
    fn chain_or_delegated_io(&self) -> Result<BamIoOptions> {
        if self.start_from == RunAllStage::Consensus {
            self.delegated_io()
        } else {
            Ok(BamIoOptions {
                input: PathBuf::new(),
                output: self.output.clone(),
                async_reader: self.async_reader,
            })
        }
    }

    /// Thin wrapper that pulls `self.start_from` and `self.stop_after_stage()`
    /// off `self` and delegates to [`validate_stages_for`].
    fn validate_stages(&self) -> Result<()> {
        validate_stages_for(self.start_from, self.stop_after_stage())
    }
    /// Validate the `--start-from align-and-merge` CLI surface
    /// without executing the chain. Runs the same
    /// [`AlignerOptions::resolve`] path the eventual C4 executor will
    /// use, so a misconfigured invocation fails before any input is
    /// read AND produces the same error message users will see at
    /// execute time.
    ///
    /// # Errors
    ///
    /// - `--ref` not set (the aligner needs a reference path).
    /// - `--input` is stdin (already caught upstream, but documented
    ///   here for completeness).
    /// - `--aligner::preset` and `--aligner::command` both unset.
    /// - `--aligner::preset` and `--aligner::command` both set
    ///   (mutually exclusive).
    /// - `--aligner-bin` or `--aligner::threads` paired with
    ///   `--aligner::command` (preset-mode-only knobs).
    ///   `--aligner::chunk-size` is permitted in both modes —
    ///   Step 1 of the AAM chain uses it for BAM-record batching
    ///   regardless of how the aligner is invoked.
    /// - Preset-mode index files missing.
    /// - Command-mode template missing `{ref}`.
    fn validate_align_and_merge(&self) -> Result<()> {
        let reference = self.reference.as_ref().ok_or_else(|| {
            anyhow::anyhow!(
                "a runall chain that includes align requires --ref (the aligner \
                 reference FASTA with its index files alongside)"
            )
        })?;
        // `--methylation-mode` drives the *consensus* stage (which AAM
        // never reaches alone) and would also conflict with the
        // aligner's reference handling. Reject explicitly rather than
        // silently ignore; methylation-aware AAM presets are a
        // follow-up PR per the design doc.
        if self.methylation_mode.is_some() {
            bail!(
                "--methylation-mode is not yet supported for runall chains that include \
                 align. For EM-seq today, use `--aligner::command \"bwameth.py ...\"` (or \
                 `bwa-mem3 --methylation-mode em-seq ...`) in command mode and apply \
                 methylation downstream as a separate step."
            );
        }
        // The downstream zipper-merge step needs a `.dict` file
        // alongside the reference FASTA (used to populate the output
        // BAM header). Validate up front so a missing `.dict` fails
        // before the aligner runs — saves potentially hours of CPU.
        if crate::reference::find_dict_path(reference).is_none() {
            bail!(
                "no sequence-dictionary file found next to --ref {} \
                 (expected `<ref>.dict` or `<ref-without-extension>.dict`). \
                 Generate one with `samtools dict {} -o <ref>.dict` before running \
                 a runall chain that includes align.",
                reference.display(),
                reference.display(),
            );
        }
        let top_threads = self.threading.num_threads();
        // `resolve` consumes the options struct; we clone so the
        // validator can be called multiple times if needed.
        let _resolved = self.aligner_opts.clone().validate()?.resolve(
            reference,
            top_threads,
            self.aligner_bin.as_deref(),
        )?;
        Ok(())
    }

    /// Reference to thread through to a consensus stage. The
    /// standalone consensus commands (`Simplex`, `Duplex`) enforce
    /// `--ref requires --methylation-mode`, so passing the
    /// runall-level `--ref` unconditionally would falsely bail when
    /// the ref came in for upstream stages (e.g. AAM aligner
    /// reference, zipper dict). Return the ref ONLY when methylation
    /// is also requested.
    fn consensus_reference(&self) -> Option<PathBuf> {
        if self.methylation_mode.is_some() { self.reference.clone() } else { None }
    }

    // ── T3b.2 spec-construction helpers ────────────────────────────────────
    //
    // These four methods extract spec-construction logic from the 15-dispatcher
    // tangle so T3b.3 can collapse `RunAll::execute` into a single
    // `ChainSpec` construction + `build_for` call. All four are ADDITIVE —
    // the existing dispatchers remain untouched until T3b.3.
    //
    /// Return the ordered [`Stage`] list for this `--start-from` /
    /// `--stop-after` pair. Thin wrapper that reads `self` and delegates
    /// to the free function [`derive_stages_for`], which is unit-tested
    /// directly without constructing a full [`RunAll`] instance.
    ///
    /// See [`derive_stages_for`] for the complete stage-derivation rules.
    ///
    /// # Errors
    ///
    /// Propagates errors from [`derive_stages_for`].
    pub(crate) fn derive_stages(&self) -> Result<Vec<crate::pipeline::chains::Stage>> {
        derive_stages_for(
            self.start_from,
            self.stop_after_stage(),
            self.extract_chain_wants_correct(),
            self.consensus_mode,
        )
    }

    /// Build the [`SourceSpec`] for the `ChainSpec` derived from
    /// `--start-from`.
    ///
    /// | `--start-from`     | `SourceSpec` variant                                              |
    /// |--------------------|-------------------------------------------------------------------|
    /// | `correct`          | `Bam(self.require_input()?)` — unmapped BAM feeding correct       |
    /// | `align`            | `Bam(self.require_input()?)` — unmapped BAM feeding AAM           |
    /// | `zipper`           | `PairedBams { unmapped, mapped: input, reference }`               |
    /// | `sort` / `group` / consensus | `Bam(self.require_input()?)` — aligned BAM             |
    ///
    /// # Errors
    ///
    /// - `Zipper` start requires `--unmapped` and `--ref`; returns `Err` if
    ///   either is absent.
    pub(crate) fn derive_source_spec(&self) -> Result<crate::pipeline::chains::SourceSpec> {
        use crate::pipeline::chains::SourceSpec;
        match self.start_from {
            RunAllStage::Extract => {
                // When --start-from extract, FASTQ paths come from
                // --extract::inputs, not from -i/--input.
                if self.input.is_some() {
                    bail!(
                        "when --start-from extract, pass --extract::inputs instead of \
                         -i/--input (--input is for BAM-source start-from values)"
                    );
                }
                let opts = self.extract_opts.clone().validate()?;
                anyhow::ensure!(
                    !opts.inputs.is_empty(),
                    "--extract::inputs is required when --start-from extract"
                );
                anyhow::ensure!(
                    !opts.read_structures.is_empty(),
                    "--extract::read-structures is required when --start-from extract"
                );
                anyhow::ensure!(
                    opts.inputs.len() == opts.read_structures.len(),
                    "--extract::inputs and --extract::read-structures must have the same count"
                );
                Ok(SourceSpec::Fastqs { paths: opts.inputs, read_structures: opts.read_structures })
            }

            RunAllStage::Correct
            | RunAllStage::AlignAndMerge
            | RunAllStage::Sort
            | RunAllStage::Group
            | RunAllStage::Consensus
            | RunAllStage::Filter => Ok(SourceSpec::Bam(self.require_input()?)),

            RunAllStage::Zipper => {
                let unmapped = self.unmapped.clone().ok_or_else(|| {
                    anyhow::anyhow!(
                        "--unmapped is required with --start-from zipper (the second \
                         input BAM, queryname-sorted, that carries the source tags)"
                    )
                })?;
                let reference = self.reference.clone().ok_or_else(|| {
                    anyhow::anyhow!(
                        "--ref is required with --start-from zipper (the reference FASTA \
                         with an accompanying `.dict` for the output BAM header)"
                    )
                })?;
                Ok(SourceSpec::PairedBams { unmapped, mapped: self.require_input()?, reference })
            }
        }
    }

    /// Build the [`SinkSpec`] for the `ChainSpec`.
    ///
    /// Runall always writes a plain BAM via `SinkSpec::Bam` — `BamWithIndex`
    /// is never constructed here because runall's sort output is either
    /// intermediate or a final template-coordinate BAM, neither of which
    /// admits a valid BAI. Runall does not expose a `--sort::write-index`
    /// flag.
    ///
    /// # Errors
    ///
    /// Currently infallible; returns `Result` for API symmetry with the
    /// other `derive_*` helpers.
    pub(crate) fn derive_sink_spec(&self) -> Result<crate::pipeline::chains::SinkSpec> {
        Ok(crate::pipeline::chains::SinkSpec::Bam(self.output.clone()))
    }

    /// Build the [`StageOptionsBag`] for the `ChainSpec` derived from
    /// the active stages in `derive_stages()`.
    ///
    /// Each active stage's options slot is populated; inactive slots remain
    /// `None`. Fields annotated `#[arg(skip)]` on the per-stage options
    /// structs are populated here using the same logic the existing
    /// single-stage `execute_*_only` paths apply.
    ///
    /// # Option-population rules by stage
    ///
    /// * **Correct** — validates `MultiCorrectOptions`, then sets
    ///   `rejects_path`:
    ///   - Self-pair (`--stop-after correct`): `rejects_path` = `self.rejects_opts.rejects`.
    ///   - Cross-stage (correct feeds AAM): `rejects_path` = `None` (the fused
    ///     chain uses the kept-only correct step; UMI rejects are discarded —
    ///     `RunAll::execute` emits a warning when `--rejects` is set on a
    ///     chained correct run).
    ///
    /// * **Align** — constructs [`AlignOptions`] from `--aligner::*` +
    ///   `--aligner-bin` + `--ref`. Requires `--ref` to be set.
    ///
    /// * **Zipper** — validates `MultiZipperOptions`. No `#[arg(skip)]` fields.
    ///
    /// * **Sort** — validates `MultiSortOptions`, then forces
    ///   `order = TemplateCoordinate` (runall's sort step always feeds
    ///   downstream group; coordinate order is not applicable). BAI
    ///   indexing is gated separately, at `derive_sink_spec` time, by
    ///   never emitting `SinkSpec::BamWithIndex` — runall does not
    ///   expose a `--sort::write-index` flag and template-coordinate
    ///   sorted BAMs cannot be BAI-indexed.
    ///
    /// * **Group** — validates `MultiGroupOptions`, then computes
    ///   `effective_strategy` / `effective_edits` and clears histogram
    ///   paths (anti-goal documented in the module-level doc).
    ///
    /// * **Simplex** — constructs a `Simplex` struct from runall's
    ///   cross-cutting fields (`io`, `rejects_opts`, `stats_opts`, etc.)
    ///   mirroring what `execute_consensus_only` does.
    ///
    /// * **Duplex** — validates `MultiDuplexOptions`, then populates the
    ///   `#[arg(skip)]` cross-cutting fields from `self`.
    ///
    /// * **Codec** — validates `MultiCodecOptions`, then populates the
    ///   `#[arg(skip)]` cross-cutting fields from `self`.
    ///
    /// # Errors
    ///
    /// Returns `Err` if any per-stage `validate()` call fails (e.g.
    /// missing required option like `--group::strategy` or
    /// `--correct::min-distance`). Also returns `Err` when `Stage::Align`
    /// is present but `--ref` is absent.
    pub(crate) fn build_stage_options_bag(
        &self,
        stages: &[crate::pipeline::chains::Stage],
    ) -> Result<crate::pipeline::chains::StageOptionsBag> {
        use crate::commands::sort::SortOrderArg;
        use crate::pipeline::chains::options_bag::AlignOptions;
        use crate::pipeline::chains::{Stage, StageOptionsBag};

        let mut bag = StageOptionsBag::default();

        for &stage in stages {
            match stage {
                Stage::Correct => {
                    let mut opts = self.correct_opts.clone().validate()?;
                    // Self-pair: wire the rejects path so the chain builder
                    // can open the rejects writer. Cross-stage: leave None
                    // (the fused correct step runs kept-only; UMI rejects
                    // are not captured in fused chains).
                    opts.rejects_path = if self.stop_after_stage() == RunAllStage::Correct {
                        self.rejects_opts.rejects.clone()
                    } else {
                        None
                    };
                    bag.correct = Some(opts);
                }

                Stage::Align => {
                    let reference = self.reference.clone().ok_or_else(|| {
                        anyhow::anyhow!(
                            "--start-from align requires --ref (the aligner reference \
                             FASTA with its index files alongside)"
                        )
                    })?;
                    let aligner = self.aligner_opts.clone().validate()?;
                    bag.aligner = Some(AlignOptions {
                        aligner,
                        reference,
                        aligner_bin: self.aligner_bin.clone(),
                    });
                }

                Stage::Zipper => {
                    let opts = self.zipper_opts.clone().validate()?;
                    // `--zipper::bwa-chunk-size` is the same deprecated, ignored
                    // field as standalone `zipper -K`; the genuinely-used bwa
                    // batch size for the align stage is `--aligner::chunk-size`.
                    crate::commands::zipper::warn_if_bwa_chunk_size_overridden(opts.bwa_chunk_size);
                    bag.zipper = Some(opts);
                }

                Stage::Sort => {
                    let mut sort_opts = self.sort_opts.clone().validate()?;
                    // Runall's sort step always produces template-coordinate
                    // output — the only order compatible with downstream Group.
                    sort_opts.order = SortOrderArg::TemplateCoordinate;
                    // After Phase 4 (T4.4), `--write-index` no longer lives
                    // on `SortOptions` — it routes through `SinkSpec`. Runall
                    // always uses `SinkSpec::Bam` (never `BamWithIndex`)
                    // because BAI is meaningless for template-coordinate
                    // sorted BAMs, so there is nothing to clear here.
                    bag.sort = Some(sort_opts);
                }

                Stage::Group => {
                    use crate::assigner::Strategy;

                    let mut group_opts = self.group_opts.clone().validate()?;
                    // Mirror GroupReadsByUmi::execute: compute
                    // effective_strategy / effective_edits, then populate
                    // the #[arg(skip)] slots before handing GroupOptions
                    // to the ChainSpec.
                    let (effective_strategy, no_umi_edits_override) = if group_opts.no_umi {
                        (Strategy::Identity, true)
                    } else {
                        (group_opts.strategy, false)
                    };
                    let effective_edits = if no_umi_edits_override
                        || matches!(effective_strategy, Strategy::Identity)
                    {
                        0
                    } else {
                        group_opts.edits
                    };
                    group_opts.effective_strategy = effective_strategy;
                    group_opts.effective_edits = effective_edits;
                    // runall never writes per-position metrics (anti-goal
                    // documented in the module-level doc comment).
                    group_opts.family_size_histogram = None;
                    group_opts.grouping_metrics = None;
                    group_opts.metrics_prefix = None;
                    bag.group = Some(group_opts);
                }

                Stage::Simplex => {
                    // Validate the `--simplex::*` flags into SimplexOptions,
                    // then populate the #[arg(skip)] chain-builder slots.
                    // Per-mode consensus tuning is the SOLE source — there is
                    // no shared top-level consensus option to overlay.
                    let mut simplex_opts = self.simplex_opts.clone().validate()?;
                    simplex_opts.io = self.chain_or_delegated_io()?;
                    simplex_opts.rejects_opts.clone_from(&self.rejects_opts);
                    simplex_opts.stats_opts.clone_from(&self.stats_opts);
                    simplex_opts.read_group.clone_from(&self.read_group);
                    simplex_opts.methylation_mode = self.methylation_mode;
                    simplex_opts.reference = self.consensus_reference();
                    bag.simplex = Some(simplex_opts);
                }

                Stage::Duplex => {
                    // Validate the `--duplex::*` flags into DuplexOptions,
                    // then populate the #[arg(skip)] chain-builder slots.
                    let mut duplex_opts = self.duplex_opts.clone().validate()?;
                    duplex_opts.io = self.chain_or_delegated_io()?;
                    duplex_opts.rejects_opts.clone_from(&self.rejects_opts);
                    duplex_opts.stats_opts.clone_from(&self.stats_opts);
                    duplex_opts.read_group.clone_from(&self.read_group);
                    duplex_opts.methylation_mode = self.methylation_mode;
                    duplex_opts.reference = self.consensus_reference();
                    bag.duplex = Some(duplex_opts);
                }

                Stage::Codec => {
                    // Codec does not support methylation calling — fail loud
                    // rather than silently dropping the flag.
                    if self.methylation_mode.is_some() {
                        bail!(
                            "--methylation-mode is not supported with codec consensus; \
                             it is valid only for simplex/duplex"
                        );
                    }
                    // Validate the `--codec::*` flags into CodecOptions,
                    // then run the numeric/semantic checks (min-reads, qual
                    // ceilings, disagreement-rate, …) that standalone `fgumi
                    // codec` runs — otherwise runall would accept degenerate
                    // configs the standalone command rejects (S5c2-001).
                    let mut codec_opts = self.codec_opts.clone().validate()?;
                    codec_opts.validate()?;
                    codec_opts.io = self.chain_or_delegated_io()?;
                    codec_opts.rejects_opts.clone_from(&self.rejects_opts);
                    codec_opts.stats_opts.clone_from(&self.stats_opts);
                    codec_opts.read_group.clone_from(&self.read_group);
                    bag.codec = Some(codec_opts);
                }

                Stage::Extract => {
                    use crate::commands::extract::{QualityEncoding, open_fastq_reader};
                    use fgumi_simd_fastq::SimdFastqReader;

                    let opts = self.extract_opts.clone().validate()?;

                    // Detect quality encoding from the first FASTQ file
                    // (same pre-build computation as standalone extract).
                    //
                    // NOTE (S5c1-004): this opens the first FASTQ and reads its
                    // leading ~1 MiB, then drops the reader; the chain's extract
                    // source re-opens the same file from the top for the real
                    // pass. This REQUIRES a re-openable (seekable) file — it must
                    // not be stdin, or the pre-read would consume it and the
                    // chain would re-open an empty stream. `open_fastq_reader`
                    // uses `File::open`, so `-`/stdin already fails here; do not
                    // wire FASTQ stdin support without first folding detection
                    // into the extract step so the file is opened only once.
                    let first_input = opts.inputs.first().ok_or_else(|| {
                        anyhow::anyhow!("--extract::inputs is required when --start-from extract")
                    })?;
                    // The pre-read below consumes the first FASTQ's leading bytes
                    // and the chain then re-opens the same path from the top, so
                    // the input MUST be a re-openable regular file. A FIFO / named
                    // pipe / socket would `File::open` fine but be *consumed* by
                    // the pre-read, silently truncating the real pass — reject any
                    // non-regular input up front with a clear error.
                    let first_meta = std::fs::metadata(first_input).map_err(|e| {
                        anyhow::anyhow!(
                            "cannot stat --extract::inputs[0] ({}): {e}",
                            first_input.display()
                        )
                    })?;
                    if !first_meta.file_type().is_file() {
                        bail!(
                            "--extract::inputs[0] ({}) must be a regular file: `fgumi runall \
                             --start-from extract` pre-reads it for quality detection and then \
                             re-opens it, which a pipe/FIFO/stdin cannot satisfy",
                            first_input.display()
                        );
                    }
                    let mut sample_quals = Vec::new();
                    let buffer_size = 1024 * 1024;
                    let quality_detection_sample_size = 400;
                    let mut temp_reader = SimdFastqReader::with_capacity(
                        open_fastq_reader(first_input, 1, opts.async_reader)?,
                        buffer_size,
                    );
                    for _i in 0..quality_detection_sample_size {
                        match temp_reader.next() {
                            Some(Ok(rec)) => sample_quals.push(rec.quality),
                            Some(Err(e)) => return Err(e.into()),
                            None => break,
                        }
                    }
                    drop(temp_reader);
                    let encoding = QualityEncoding::detect(&sample_quals)?;

                    let extract_options = opts.to_extract_options(encoding);
                    extract_options.validate()?;
                    bag.extract = Some(extract_options);
                }

                Stage::Filter => {
                    // Validate the `--filter::*` flags into FilterOptions.
                    // The standalone filter chain builder reads `rejects`
                    // and `stats` straight off the bag's FilterOptions, so
                    // `--filter::rejects` / `--filter::stats` flow through
                    // unchanged. No cross-stage rewiring is needed: filter
                    // is always terminal in a runall chain.
                    //
                    // `Vec` fields that the standalone `fgumi filter` defaults
                    // (`--max-read-error-rate` = `[0.025]`,
                    // `--max-base-error-rate` = `[0.1]`) are backfilled from
                    // `FilterOptions::default()` by the `#[multi_options]`
                    // generated `validate()`, so no per-field backfill is
                    // needed here. `--filter::min-reads` has no default (it is
                    // required), so an omitted value stays empty and
                    // `FilterOptions::validate_parameters` surfaces the error.
                    let filter_opts = self.filter_opts.clone().validate()?;

                    bag.filter = Some(filter_opts);
                }

                // Stages that runall does not use. If any of these appear in
                // the stages list it is a programming error (derive_stages
                // never emits them for a runall chain).
                Stage::Clip | Stage::Dedup | Stage::Downsample | Stage::Fastq => {
                    bail!(
                        "internal error: build_stage_options_bag encountered unexpected \
                         stage {stage:?} in a runall chain; this is a bug in derive_stages"
                    );
                }
            }
        }

        Ok(bag)
    }

    // ── end T3b.2 spec-construction helpers ────────────────────────────────
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::ValueEnum;

    #[test]
    fn runall_stage_ord_linear_order() {
        assert_eq!(RunAllStage::Extract.ord(), 0);
        assert_eq!(RunAllStage::Correct.ord(), 1);
        assert_eq!(RunAllStage::AlignAndMerge.ord(), 2);
        assert_eq!(RunAllStage::Zipper.ord(), 3);
        assert_eq!(RunAllStage::Sort.ord(), 4);
        assert_eq!(RunAllStage::Group.ord(), 5);
        assert_eq!(RunAllStage::Consensus.ord(), 6);
        assert_eq!(RunAllStage::Filter.ord(), 7);
    }

    #[test]
    fn runall_stage_is_consensus() {
        assert!(!RunAllStage::Extract.is_consensus());
        assert!(!RunAllStage::Correct.is_consensus());
        assert!(!RunAllStage::AlignAndMerge.is_consensus());
        assert!(!RunAllStage::Zipper.is_consensus());
        assert!(!RunAllStage::Sort.is_consensus());
        assert!(!RunAllStage::Group.is_consensus());
        assert!(RunAllStage::Consensus.is_consensus());
    }

    #[test]
    fn runall_stage_parses_case_insensitive() {
        // clap::ValueEnum::from_str with `ignore_case = true` is the
        // default for derived ValueEnum impls. The CLI value for
        // `AlignAndMerge` is `align` (set via the variant's
        // `#[clap(name = "align")]` override — the variant identifier
        // describes what the stage does, the CLI value is the verb
        // the user types).
        assert_eq!(RunAllStage::from_str("extract", true).unwrap(), RunAllStage::Extract);
        assert_eq!(RunAllStage::from_str("Extract", true).unwrap(), RunAllStage::Extract);
        assert_eq!(RunAllStage::from_str("EXTRACT", true).unwrap(), RunAllStage::Extract);
        assert_eq!(RunAllStage::from_str("correct", true).unwrap(), RunAllStage::Correct);
        assert_eq!(RunAllStage::from_str("Correct", true).unwrap(), RunAllStage::Correct);
        assert_eq!(RunAllStage::from_str("CORRECT", true).unwrap(), RunAllStage::Correct);
        assert_eq!(RunAllStage::from_str("align", true).unwrap(), RunAllStage::AlignAndMerge);
        assert_eq!(RunAllStage::from_str("Align", true).unwrap(), RunAllStage::AlignAndMerge);
        assert_eq!(RunAllStage::from_str("ALIGN", true).unwrap(), RunAllStage::AlignAndMerge);
        // The old kebab-of-identifier form was `align-and-merge`; pin
        // that it no longer parses so a stale doc/example surfaces as
        // a clean error.
        assert!(RunAllStage::from_str("align-and-merge", true).is_err());
        assert_eq!(RunAllStage::from_str("zipper", true).unwrap(), RunAllStage::Zipper);
        assert_eq!(RunAllStage::from_str("Zipper", true).unwrap(), RunAllStage::Zipper);
        assert_eq!(RunAllStage::from_str("ZIPPER", true).unwrap(), RunAllStage::Zipper);
        assert_eq!(RunAllStage::from_str("sort", true).unwrap(), RunAllStage::Sort);
        assert_eq!(RunAllStage::from_str("Sort", true).unwrap(), RunAllStage::Sort);
        assert_eq!(RunAllStage::from_str("SORT", true).unwrap(), RunAllStage::Sort);
        assert_eq!(RunAllStage::from_str("group", true).unwrap(), RunAllStage::Group);
        assert_eq!(RunAllStage::from_str("consensus", true).unwrap(), RunAllStage::Consensus);
        assert_eq!(RunAllStage::from_str("filter", true).unwrap(), RunAllStage::Filter);
        assert_eq!(RunAllStage::from_str("Consensus", true).unwrap(), RunAllStage::Consensus);
        // The old per-algorithm stage values are gone (the algorithm moved
        // to `--consensus`); pin that they no longer parse.
        assert!(RunAllStage::from_str("simplex", true).is_err());
        assert!(RunAllStage::from_str("duplex", true).is_err());
        assert!(RunAllStage::from_str("codec", true).is_err());
    }

    #[test]
    fn runall_stage_parse_rejects_unknown() {
        assert!(RunAllStage::from_str("merge", true).is_err());
        assert!(RunAllStage::from_str("", true).is_err());
        assert!(RunAllStage::from_str("sortgroup", true).is_err());
    }

    #[test]
    fn validate_with_accepts_extract_to_each() {
        // Extract is the earliest stage; every stage at or after it is a
        // structurally valid stop point (the extract self-pair is the
        // supported extract-only run).
        for stop in [
            RunAllStage::Extract,
            RunAllStage::Correct,
            RunAllStage::AlignAndMerge,
            RunAllStage::Zipper,
            RunAllStage::Sort,
            RunAllStage::Group,
            RunAllStage::Consensus,
        ] {
            assert!(
                RunAllStage::Extract.validate_with(stop).is_ok(),
                "extract → {stop} should be valid (structural order)"
            );
        }
    }

    #[test]
    fn validate_with_accepts_correct_to_each() {
        // Correct is the second earliest stage; every later stage is a
        // valid stop point.
        for stop in [
            RunAllStage::Correct,
            RunAllStage::AlignAndMerge,
            RunAllStage::Zipper,
            RunAllStage::Sort,
            RunAllStage::Group,
            RunAllStage::Consensus,
        ] {
            assert!(
                RunAllStage::Correct.validate_with(stop).is_ok(),
                "correct → {stop} should be valid"
            );
        }
    }

    #[test]
    fn validate_with_accepts_align_and_merge_to_each() {
        // align-and-merge is now stage 1 (Correct is stage 0); every
        // stage at or after align is a valid stop point.
        // `validate_with` is the structural-ordering check; the
        // higher-level `validate_stages_for` additionally rejects
        // `--stop-after align-and-merge` (see
        // `validate_stages_for_rejects_align_and_merge_stop`).
        for stop in [
            RunAllStage::AlignAndMerge,
            RunAllStage::Zipper,
            RunAllStage::Sort,
            RunAllStage::Group,
            RunAllStage::Consensus,
        ] {
            assert!(
                RunAllStage::AlignAndMerge.validate_with(stop).is_ok(),
                "align-and-merge → {stop} should be valid"
            );
        }
    }

    /// `--start-from correct --stop-after correct` is the supported
    /// self-pair: writes the corrected unmapped BAM to `--output`.
    /// Pin the acceptance here so EC-C3's executor has a red test to
    /// turn green and any future re-ordering of `validate_stages_for`
    /// surfaces breakage immediately.
    #[test]
    fn validate_stages_for_accepts_correct_stop() {
        validate_stages_for(RunAllStage::Correct, RunAllStage::Correct).unwrap();
    }

    /// `--start-from extract --stop-after extract` is the supported
    /// extract-only self-pair: it runs the FASTQ→unmapped-BAM extract
    /// chain inside runall and writes the unmapped BAM to `--output`
    /// (byte-identical to standalone `fgumi extract`).
    #[test]
    fn validate_stages_for_accepts_extract_self_pair() {
        validate_stages_for(RunAllStage::Extract, RunAllStage::Extract).unwrap();
    }

    /// Consensus is one stage; the algorithm is chosen with `--consensus`.
    /// `--stop-after consensus --consensus simplex` is the new surface for
    /// what used to be `--stop-after simplex`.
    #[test]
    fn runall_stop_after_consensus_with_mode_parses() {
        use clap::Parser;
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "--consensus",
            "simplex",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--threads",
            "1",
        ]);
        assert!(r.is_ok(), "--stop-after consensus --consensus simplex must parse: {r:?}");
    }

    /// `--duplex::min-reads` is a `Vec` whose `default_value` clap cannot
    /// carry to the Multi side, so it arrives empty when omitted. The
    /// `#[multi_options]`-generated `validate()` backfills the standalone
    /// default (`DuplexOptions::default().min_reads == [1]`) so omitting it
    /// yields default behavior, not an empty-Vec panic in the consensus caller.
    #[test]
    fn runall_duplex_min_reads_defaults_to_one_when_omitted() {
        use clap::Parser;
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "--consensus",
            "duplex",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--group::strategy",
            "paired",
            "--threads",
            "1",
        ])
        .expect("parse");
        let stages = r.derive_stages().expect("derive stages");
        let bag = r.build_stage_options_bag(&stages).expect("build bag");
        assert_eq!(
            bag.duplex.expect("duplex bag entry").min_reads,
            vec![1],
            "--duplex::min-reads must default to [1] when omitted"
        );
    }

    /// Codec consensus does not support methylation calling. The
    /// `Stage::Codec` arm of `build_stage_options_bag` rejects a
    /// `--methylation-mode` that survived parsing — fail loud rather than
    /// silently drop the flag. The guard fires before any `--codec::*`
    /// validation, so no codec tuning flags are needed to reach it.
    #[test]
    fn codec_with_methylation_is_rejected() {
        use clap::Parser;
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "--consensus",
            "codec",
            "--methylation-mode",
            "em-seq",
            "--ref",
            "/tmp/fgumi-nonexistent-ref.fa",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--group::strategy",
            "paired",
            "--threads",
            "1",
        ])
        .expect("parse");
        let stages = r.derive_stages().expect("derive stages");
        assert!(
            stages.contains(&crate::pipeline::chains::Stage::Codec),
            "group → consensus with --consensus codec must include the Codec stage: {stages:?}"
        );
        let err = match r.build_stage_options_bag(&stages) {
            Ok(_) => panic!("--consensus codec + --methylation-mode must be rejected"),
            Err(e) => e,
        };
        let msg = err.to_string();
        assert!(
            msg.contains("methylation-mode") && msg.contains("codec"),
            "error must mention methylation-mode and codec, got: {msg}"
        );
    }

    /// `--start-from align` (align-and-merge) does not yet support
    /// methylation-aware aligning. `validate_align_and_merge` rejects a
    /// `--methylation-mode` after confirming `--ref` is present but before
    /// the `.dict` existence check, so a nonexistent reference path still
    /// reaches the methylation guard.
    #[test]
    fn align_start_with_methylation_is_rejected() {
        use clap::Parser;
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "align",
            "--stop-after",
            "consensus",
            "--consensus",
            "duplex",
            "--methylation-mode",
            "em-seq",
            "--ref",
            "/tmp/fgumi-nonexistent-ref.fa",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--group::strategy",
            "paired",
            "--threads",
            "1",
        ])
        .expect("parse");
        let err = r
            .validate_align_and_merge()
            .expect_err("--start-from align + --methylation-mode must be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains("methylation-mode") && msg.contains("align"),
            "error must mention methylation-mode and align, got: {msg}"
        );
    }

    /// On a non-aligner start stage (e.g. `--start-from group`), passing
    /// `--ref` without `--methylation-mode` is rejected by `execute` before
    /// any stage runs: the ref would otherwise be silently ignored. The
    /// guard is checked early — after the (skipped-for-stdin) input
    /// existence check but before the pipeline opens — so `--input -`
    /// reaches the guard without a real input BAM and without consuming
    /// stdin.
    #[test]
    fn ref_without_methylation_on_group_start_is_rejected() {
        use clap::Parser;
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "group",
            "--ref",
            "/tmp/fgumi-nonexistent-ref.fa",
            "--input",
            "-",
            "--output",
            "o.bam",
            "--group::strategy",
            "paired",
            "--threads",
            "1",
        ])
        .expect("parse");
        let err = r
            .execute("test")
            .expect_err("--ref without --methylation-mode on a non-aligner start must be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains("--ref requires --methylation-mode"),
            "error must mention `--ref requires --methylation-mode`, got: {msg}"
        );
    }

    /// Symmetric to `ref_without_methylation_on_group_start_is_rejected`:
    /// `--methylation-mode` on a chain that never reaches consensus (here
    /// `group → group`) is dead — no stage consumes it — so `execute` rejects it
    /// before the pipeline opens rather than silently dropping it. `--ref` is
    /// supplied too (so the existing `--ref requires --methylation-mode` guard
    /// does NOT fire and the methylation guard is the one under test), and
    /// `--input -` reaches the guard without consuming a real BAM. The error
    /// must name `--methylation-mode` and `consensus` (not the align-specific
    /// message, which `validate_align_and_merge` owns for align chains).
    #[test]
    fn methylation_mode_without_consensus_stage_is_rejected() {
        use clap::Parser;
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "group",
            "--methylation-mode",
            "em-seq",
            "--ref",
            "/tmp/fgumi-nonexistent-ref.fa",
            "--input",
            "-",
            "--output",
            "o.bam",
            "--group::strategy",
            "paired",
            "--threads",
            "1",
        ])
        .expect("parse");
        let err = r.execute("test").expect_err(
            "--methylation-mode on a non-consensus group → group chain must be rejected",
        );
        let msg = err.to_string();
        assert!(
            msg.contains("--methylation-mode") && msg.contains("consensus"),
            "error must mention `--methylation-mode` and `consensus`, got: {msg}"
        );
        assert!(
            !msg.contains("align"),
            "non-align chain must not get the align-specific message, got: {msg}"
        );
    }

    /// `--methylation-mode` on a simplex/duplex consensus chain requires `--ref`:
    /// the consensus stage consumes the FASTA. A non-align `group → simplex` run
    /// with `--methylation-mode` but no `--ref` reaches consensus (so the
    /// dead-consensus guard does NOT fire) yet would otherwise thread
    /// `methylation_mode = Some(...)` / `reference = None` into the stage. The
    /// `execute` guard must reject it up front with a clear message. No
    /// `--ref`/`--simplex::*` flags are supplied because the guard fires before
    /// `build_stage_options_bag()`; `--input -` reaches it without opening a BAM.
    #[test]
    fn methylation_mode_without_ref_on_consensus_chain_is_rejected() {
        use clap::Parser;
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "consensus",
            "--consensus",
            "simplex",
            "--methylation-mode",
            "em-seq",
            "--input",
            "-",
            "--output",
            "o.bam",
            "--group::strategy",
            "paired",
            "--threads",
            "1",
        ])
        .expect("parse");
        let err = r.execute("test").expect_err(
            "--methylation-mode on a group → simplex chain without --ref must be rejected",
        );
        let msg = err.to_string();
        assert!(
            msg.contains("--methylation-mode requires --ref"),
            "error must mention `--methylation-mode requires --ref`, got: {msg}"
        );
    }

    /// `--filter::max-read-error-rate` / `--filter::max-base-error-rate` are
    /// `Vec` fields with standalone defaults (`[0.025]` / `[0.1]`) that the
    /// `#[multi_options]`-generated `validate()` backfills from
    /// `FilterOptions::default()` when omitted. `--filter::min-reads` has no
    /// default (it is required), so it must still be supplied; this pins both
    /// behaviors so the runall filter surface matches standalone `fgumi filter`.
    #[test]
    fn runall_filter_vec_defaults_backfilled_when_omitted() {
        use clap::Parser;
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "filter",
            "--stop-after",
            "filter",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--filter::min-reads",
            "1",
            "--threads",
            "1",
        ])
        .expect("parse");
        let stages = r.derive_stages().expect("derive stages");
        let bag = r.build_stage_options_bag(&stages).expect("build bag");
        let filter = bag.filter.expect("filter bag entry");
        assert_eq!(
            filter.max_read_error_rate,
            vec![0.025],
            "--filter::max-read-error-rate must default to [0.025] when omitted"
        );
        assert_eq!(
            filter.max_base_error_rate,
            vec![0.1],
            "--filter::max-base-error-rate must default to [0.1] when omitted"
        );
        assert_eq!(filter.min_reads, vec![1], "--filter::min-reads passes through");
    }

    /// The old per-algorithm `--stop-after simplex|duplex|codec` values are
    /// gone — the algorithm moved to `--consensus`.
    #[test]
    fn runall_old_stop_after_simplex_rejected() {
        use clap::Parser;
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "simplex",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--threads",
            "1",
        ]);
        assert!(r.is_err(), "--stop-after simplex should no longer be a valid value");
    }

    /// `--start-from extract` reads FASTQ from `--extract::inputs`, so
    /// `--input` must NOT be clap-required: runall's input contract is
    /// stage-dependent, unlike the standard single-required-BAM commands.
    #[test]
    fn runall_extract_start_parses_without_input_flag() {
        use clap::Parser;
        let parsed = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "extract",
            "--stop-after",
            "extract",
            "--output",
            "/tmp/o.bam",
            "--extract::inputs",
            "a.fq.gz",
            "b.fq.gz",
            "--extract::read-structures",
            "+T",
            "+T",
            "--extract::sample",
            "s",
            "--extract::library",
            "lib",
            "--threads",
            "1",
        ]);
        assert!(parsed.is_ok(), "extract-start must parse without --input: {parsed:?}");
    }

    /// `--stop-after align` is rejected at the clap parse layer (not just
    /// at runtime validation): `align` is a valid `--start-from` but not a
    /// valid `--stop-after` (raw aligner output before the zipper-merge
    /// loses every tag), so it is absent from the `StopAfter` value-enum.
    #[test]
    fn runall_stop_after_align_rejected_at_parse() {
        use clap::Parser;
        let res = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "align",
            "--stop-after",
            "align",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--ref",
            "r.fa",
            "--threads",
            "1",
        ]);
        let err = res.expect_err("--stop-after align must be rejected at the parse layer");
        let msg = err.to_string().to_lowercase();
        assert!(
            msg.contains("invalid value") && msg.contains("align"),
            "expected clap invalid-value error for align, got: {err}"
        );
    }

    /// `zipper` (the earliest post-align stop) is still an accepted
    /// `--stop-after` value after the `StopAfter` split.
    #[test]
    fn runall_stop_after_zipper_accepted_at_parse() {
        use clap::Parser;
        let parsed = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "align",
            "--stop-after",
            "zipper",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--ref",
            "r.fa",
            "--threads",
            "1",
        ]);
        assert!(parsed.is_ok(), "--stop-after zipper must parse: {parsed:?}");
    }

    /// Pins `StopAfter` to be exactly `RunAllStage` minus `align`, so the
    /// two enums cannot silently drift (e.g. when a new stage is added).
    #[test]
    fn stop_after_values_are_runallstage_minus_align() {
        use clap::ValueEnum;
        for sa in StopAfter::value_variants() {
            let stage: RunAllStage = (*sa).into();
            assert_ne!(stage, RunAllStage::AlignAndMerge, "StopAfter must not include align");
        }
        assert_eq!(
            StopAfter::value_variants().len(),
            RunAllStage::value_variants().len() - 1,
            "StopAfter must be exactly RunAllStage minus align"
        );
    }

    /// BAM start stages still require an input BAM. With `--input`
    /// optional at the clap level, presence is enforced at
    /// `derive_source_spec` with a clear, stage-specific error.
    #[test]
    fn runall_bam_start_without_input_errors_at_source_spec() {
        use clap::Parser;
        let runall = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "sort",
            "--output",
            "/tmp/o.bam",
            "--threads",
            "1",
        ])
        .expect("sort-start must parse without --input (presence checked at runtime)");
        let err = runall.derive_source_spec().unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("--input"), "expected --input-required error, got: {msg}");
    }

    /// `--stop-after align-and-merge` is rejected at the higher-level
    /// `validate_stages_for` (not at `validate_with`, which checks
    /// only structural ordering). The design dropped this stop point
    /// because raw aligner output without zipper-merge is not useful.
    /// Both the AAM self-pair and the new `correct → align` combo
    /// reach this guard.
    #[test]
    fn validate_stages_for_rejects_align_and_merge_stop() {
        for start in [RunAllStage::Extract, RunAllStage::Correct, RunAllStage::AlignAndMerge] {
            let err = validate_stages_for(start, RunAllStage::AlignAndMerge).unwrap_err();
            let msg = format!("{err:#}");
            assert!(
                msg.contains("--stop-after align is not supported"),
                "expected AAM-stop rejection for {start} → align, got: {msg}"
            );
        }
    }

    #[test]
    fn validate_with_rejects_align_and_merge_at_stop() {
        // Only `Correct → align` is structurally valid (and rejected
        // separately by `validate_stages_for`'s
        // `--stop-after align` guard). Every other stage comes after
        // align in the linear order.
        for start in
            [RunAllStage::Zipper, RunAllStage::Sort, RunAllStage::Group, RunAllStage::Consensus]
        {
            let err = start.validate_with(RunAllStage::AlignAndMerge).unwrap_err();
            let msg = format!("{err:#}");
            assert!(
                msg.contains("comes after"),
                "expected order violation for {start} → align-and-merge, got: {msg}"
            );
        }
    }

    #[test]
    fn validate_with_accepts_zipper_to_each() {
        for stop in
            [RunAllStage::Zipper, RunAllStage::Sort, RunAllStage::Group, RunAllStage::Consensus]
        {
            assert!(
                RunAllStage::Zipper.validate_with(stop).is_ok(),
                "zipper → {stop} should be valid"
            );
        }
    }

    #[test]
    fn validate_with_accepts_sort_to_each() {
        for stop in [RunAllStage::Sort, RunAllStage::Group, RunAllStage::Consensus] {
            assert!(RunAllStage::Sort.validate_with(stop).is_ok(), "sort → {stop} should be valid");
        }
    }

    #[test]
    fn validate_with_accepts_group_to_each_consensus() {
        for stop in [RunAllStage::Group, RunAllStage::Consensus] {
            assert!(
                RunAllStage::Group.validate_with(stop).is_ok(),
                "group → {stop} should be valid"
            );
        }
    }

    #[test]
    fn validate_with_accepts_consensus_self_pair() {
        assert!(
            RunAllStage::Consensus.validate_with(RunAllStage::Consensus).is_ok(),
            "consensus → consensus should be valid"
        );
    }

    #[test]
    fn validate_with_rejects_reverse_order() {
        // stop-after earlier than start-from — the order check (the only
        // rule now that consensus is a single stage).
        for (s, e) in [
            (RunAllStage::Sort, RunAllStage::Zipper),
            (RunAllStage::Group, RunAllStage::Zipper),
            (RunAllStage::Group, RunAllStage::Sort),
            (RunAllStage::Consensus, RunAllStage::Sort),
            (RunAllStage::Consensus, RunAllStage::Group),
            (RunAllStage::Consensus, RunAllStage::Zipper),
        ] {
            let err = s.validate_with(e).unwrap_err();
            let msg = format!("{err:#}");
            assert!(msg.contains("comes after"), "expected order-violation message, got {msg}");
            // Lowercase variant names in the error message match the
            // CLI flag values (Display, not Debug).
            assert!(msg.contains(&format!("--start-from {s}")), "got: {msg}");
            assert!(msg.contains(&format!("--stop-after {e}")), "got: {msg}");
        }
    }

    #[test]
    fn validate_with_rejects_consensus_to_group() {
        // Consensus → earlier stage is rejected by the order check.
        let err = RunAllStage::Consensus.validate_with(RunAllStage::Group).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("comes after"), "expected order-violation message, got {msg}");
    }

    #[test]
    fn display_uses_lowercase_to_match_cli() {
        // Error messages, --help text, and any user-visible string
        // must use the lowercase CLI form. Pin to prevent accidental
        // Debug-format leakage.
        assert_eq!(format!("{}", RunAllStage::Extract), "extract");
        assert_eq!(format!("{}", RunAllStage::Sort), "sort");
        assert_eq!(format!("{}", RunAllStage::Group), "group");
        assert_eq!(format!("{}", RunAllStage::Consensus), "consensus");
    }
}

#[cfg(test)]
mod run_all_accessor_tests {
    //! Tests for the derived-accessor MAPPINGS that
    //! `RunAll::is_sort_enabled` and `RunAll::require_consensus_mode`
    //! delegate to (task 3 of #33).
    //!
    //! These tests pin the `RunAllStage → bool` and
    //! `RunAllStage → RunAllMode` mappings — not the accessor methods
    //! themselves, which are trivial one-line wrappers and cannot be
    //! invoked on a `RunAll` instance from a unit test (constructing
    //! one via clap is blocked by an unrelated `-s` conflict). End-to-end
    //! coverage of the accessors comes from the runall parity
    //! integration tests in `tests/integration/test_runall_parity.rs`.
    //!
    //! `RunAllStage::consensus_mode` has standalone coverage in the
    //! The consensus algorithm now comes from `--consensus` (a
    //! `RunAllMode`), not the stage, so there is no stage→mode mapping.
    use super::*;

    #[test]
    fn validate_strategy_for_mode_accepts_compatible_pairs() {
        validate_strategy_for_mode(RunAllMode::Duplex, Strategy::Paired)
            .expect("duplex + paired must be accepted");
        validate_strategy_for_mode(RunAllMode::Simplex, Strategy::Adjacency)
            .expect("simplex + non-paired must be accepted");
        validate_strategy_for_mode(RunAllMode::Codec, Strategy::Adjacency)
            .expect("codec + non-paired must be accepted");
        validate_strategy_for_mode(RunAllMode::Simplex, Strategy::Identity)
            .expect("simplex + identity must be accepted");
    }

    #[test]
    fn validate_strategy_for_mode_duplex_without_paired_message_pins_new_flag_name() {
        let err = validate_strategy_for_mode(RunAllMode::Duplex, Strategy::Adjacency).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("--consensus duplex"), "got: {msg}");
        assert!(msg.contains("--strategy paired"), "got: {msg}");
        assert!(!msg.contains("--mode"), "stale --mode reference: {msg}");
    }

    #[test]
    fn validate_strategy_for_mode_simplex_with_paired_uses_lowercase_display() {
        // Pins `Display for RunAllMode` (lowercase) for the
        // `{mode}` interpolation in the error message — regression
        // guard against accidental `{:?}` (Debug) revert.
        let err = validate_strategy_for_mode(RunAllMode::Simplex, Strategy::Paired).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("--consensus simplex"),
            "expected lowercase '--consensus simplex', got: {msg}"
        );
        assert!(
            !msg.contains("--consensus Simplex"),
            "Debug-format leaked into error message: {msg}"
        );
        assert!(msg.contains("non-paired strategy"), "got: {msg}");
    }

    #[test]
    fn validate_strategy_for_mode_codec_with_paired_uses_lowercase_display() {
        let err = validate_strategy_for_mode(RunAllMode::Codec, Strategy::Paired).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("--consensus codec"), "got: {msg}");
        assert!(!msg.contains("--consensus Codec"), "got: {msg}");
    }

    #[test]
    fn is_sort_enabled_mapping_table() {
        // Pins the `matches!(start_from, Sort)` rule the accessor
        // delegates to. Only `Sort` enables the in-pipeline sort
        // step; all other start stages skip it. If the accessor
        // body is ever rewritten, this table must change in lockstep
        // — the assertion failure highlights the surface.
        for (stage, expected) in [
            (RunAllStage::Sort, true),
            (RunAllStage::Group, false),
            (RunAllStage::Consensus, false),
        ] {
            assert_eq!(
                matches!(stage, RunAllStage::Sort),
                expected,
                "is_sort_enabled mapping for {stage}"
            );
        }
    }
}

#[cfg(test)]
mod derive_stages_tests {
    //! Unit tests for [`derive_stages_for`].
    //!
    //! Every legal `(start_from, stop_after)` combination is tested.
    //! Illegal combinations (`start > stop`, `--stop-after align`) are
    //! exercised by the `validate_stages_for` tests; `derive_stages_for`
    //! is always called after `validate_stages`.
    //!
    //! Test naming convention:
    //!   `derive_stages_{start}_to_{stop}` where both names are lowercase.
    //!   Consensus stops are represented as `simplex`, `duplex`, and `codec`.
    //!   Tests that cover all three consensus variants share a parameterised
    //!   helper.

    use super::*;
    use crate::pipeline::chains::Stage;

    // ── Extract → Correct (with UMI source) ────────────────────────────────

    #[test]
    fn derive_stages_extract_self_pair() {
        // extract → extract: the extract-only chain is a single
        // terminal Extract stage (FASTQ → unmapped BAM), no align.
        let stages = derive_stages_for(RunAllStage::Extract, RunAllStage::Extract, false, None)
            .expect("extract → extract must succeed");
        assert_eq!(stages, vec![Stage::Extract]);
    }

    #[test]
    fn derive_stages_extract_to_correct_with_umi_source() {
        // extract → correct (when extract_wants_correct = true).
        let stages = derive_stages_for(RunAllStage::Extract, RunAllStage::Correct, true, None)
            .expect("extract → correct (with UMI source) must succeed");
        assert_eq!(stages, vec![Stage::Extract, Stage::Correct]);
    }

    // ── Extract → Sort (chains through Align + Sort-forcing) ───────────────

    #[test]
    fn derive_stages_extract_to_sort_without_correct() {
        // extract → align → sort. No Correct stage because
        // extract_wants_correct = false.
        let stages = derive_stages_for(RunAllStage::Extract, RunAllStage::Sort, false, None)
            .expect("extract → sort must succeed");
        assert_eq!(stages, vec![Stage::Extract, Stage::Align, Stage::Sort]);
    }

    #[test]
    fn derive_stages_extract_to_sort_with_correct() {
        // extract → correct → align → sort.
        let stages = derive_stages_for(RunAllStage::Extract, RunAllStage::Sort, true, None)
            .expect("extract → sort (with correct) must succeed");
        assert_eq!(stages, vec![Stage::Extract, Stage::Correct, Stage::Align, Stage::Sort]);
    }

    // ── Extract → Group ─────────────────────────────────────────────────────

    #[test]
    fn derive_stages_extract_to_group() {
        let stages = derive_stages_for(RunAllStage::Extract, RunAllStage::Group, false, None)
            .expect("extract → group must succeed");
        assert_eq!(stages, vec![Stage::Extract, Stage::Align, Stage::Sort, Stage::Group]);
    }

    // ── Extract → consensus ─────────────────────────────────────────────────

    #[test]
    fn derive_stages_extract_to_simplex() {
        let stages = derive_stages_for(
            RunAllStage::Extract,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Simplex),
        )
        .expect("extract → simplex must succeed");
        assert_eq!(
            stages,
            vec![Stage::Extract, Stage::Align, Stage::Sort, Stage::Group, Stage::Simplex]
        );
    }

    #[test]
    fn derive_stages_extract_to_duplex() {
        let stages = derive_stages_for(
            RunAllStage::Extract,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Duplex),
        )
        .expect("extract → duplex must succeed");
        assert_eq!(
            stages,
            vec![Stage::Extract, Stage::Align, Stage::Sort, Stage::Group, Stage::Duplex]
        );
    }

    #[test]
    fn derive_stages_extract_to_codec() {
        let stages = derive_stages_for(
            RunAllStage::Extract,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Codec),
        )
        .expect("extract → codec must succeed");
        assert_eq!(
            stages,
            vec![Stage::Extract, Stage::Align, Stage::Sort, Stage::Group, Stage::Codec]
        );
    }

    #[test]
    fn derive_stages_extract_to_zipper() {
        // NOTE: `execute` rejects `--start-from extract --stop-after zipper`
        // upfront (zipper requires PairedBams), so this combination never
        // reaches `derive_stages_for` in production. This test pins the
        // pure derivation behavior in isolation: AAM includes the
        // zipper-merge internally, so the derived shape is [Extract, Align].
        let stages = derive_stages_for(RunAllStage::Extract, RunAllStage::Zipper, false, None)
            .expect("extract → zipper must succeed");
        assert_eq!(stages, vec![Stage::Extract, Stage::Align]);
    }

    // ── Correct self-pair ───────────────────────────────────────────────────

    #[test]
    fn derive_stages_correct_to_correct() {
        let stages = derive_stages_for(RunAllStage::Correct, RunAllStage::Correct, false, None)
            .expect("correct → correct must succeed");
        assert_eq!(stages, vec![Stage::Correct]);
    }

    // ── Correct → Zipper (chains through AAM) ──────────────────────────────

    #[test]
    fn derive_stages_correct_to_zipper() {
        // correct → align (AAM includes zipper-merge internally).
        let stages = derive_stages_for(RunAllStage::Correct, RunAllStage::Zipper, false, None)
            .expect("correct → zipper must succeed");
        assert_eq!(stages, vec![Stage::Correct, Stage::Align]);
    }

    // ── Correct → Sort (chains through AAM + sort-forcing) ─────────────────

    #[test]
    fn derive_stages_correct_to_sort() {
        // correct → align (AAM) → sort (forced because Align output is
        // queryname-sorted and the stop is at Sort).
        let stages = derive_stages_for(RunAllStage::Correct, RunAllStage::Sort, false, None)
            .expect("correct → sort must succeed");
        assert_eq!(stages, vec![Stage::Correct, Stage::Align, Stage::Sort]);
    }

    // ── Correct → Group ─────────────────────────────────────────────────────

    #[test]
    fn derive_stages_correct_to_group() {
        // correct → align → sort (forced) → group.
        let stages = derive_stages_for(RunAllStage::Correct, RunAllStage::Group, false, None)
            .expect("correct → group must succeed");
        assert_eq!(stages, vec![Stage::Correct, Stage::Align, Stage::Sort, Stage::Group]);
    }

    // ── Correct → consensus (three variants) ────────────────────────────────

    #[test]
    fn derive_stages_correct_to_simplex() {
        let stages = derive_stages_for(
            RunAllStage::Correct,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Simplex),
        )
        .expect("correct → simplex must succeed");
        assert_eq!(
            stages,
            vec![Stage::Correct, Stage::Align, Stage::Sort, Stage::Group, Stage::Simplex]
        );
    }

    #[test]
    fn derive_stages_correct_to_duplex() {
        let stages = derive_stages_for(
            RunAllStage::Correct,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Duplex),
        )
        .expect("correct → duplex must succeed");
        assert_eq!(
            stages,
            vec![Stage::Correct, Stage::Align, Stage::Sort, Stage::Group, Stage::Duplex]
        );
    }

    #[test]
    fn derive_stages_correct_to_codec() {
        let stages = derive_stages_for(
            RunAllStage::Correct,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Codec),
        )
        .expect("correct → codec must succeed");
        assert_eq!(
            stages,
            vec![Stage::Correct, Stage::Align, Stage::Sort, Stage::Group, Stage::Codec]
        );
    }

    // ── AlignAndMerge → Zipper (AAM self-pair at zipper boundary) ───────────

    #[test]
    fn derive_stages_align_to_zipper() {
        // AlignAndMerge includes the zipper-merge: stop at Zipper means
        // just [Align].
        let stages =
            derive_stages_for(RunAllStage::AlignAndMerge, RunAllStage::Zipper, false, None)
                .expect("align → zipper must succeed");
        assert_eq!(stages, vec![Stage::Align]);
    }

    // ── AlignAndMerge → Sort ─────────────────────────────────────────────────

    #[test]
    fn derive_stages_align_to_sort() {
        // Align output is queryname-sorted; Sort is forced.
        let stages = derive_stages_for(RunAllStage::AlignAndMerge, RunAllStage::Sort, false, None)
            .expect("align → sort must succeed");
        assert_eq!(stages, vec![Stage::Align, Stage::Sort]);
    }

    // ── AlignAndMerge → Group ────────────────────────────────────────────────

    #[test]
    fn derive_stages_align_to_group() {
        let stages = derive_stages_for(RunAllStage::AlignAndMerge, RunAllStage::Group, false, None)
            .expect("align → group must succeed");
        assert_eq!(stages, vec![Stage::Align, Stage::Sort, Stage::Group]);
    }

    // ── AlignAndMerge → consensus ────────────────────────────────────────────

    #[test]
    fn derive_stages_align_to_simplex() {
        let stages = derive_stages_for(
            RunAllStage::AlignAndMerge,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Simplex),
        )
        .expect("align → simplex must succeed");
        assert_eq!(stages, vec![Stage::Align, Stage::Sort, Stage::Group, Stage::Simplex]);
    }

    #[test]
    fn derive_stages_align_to_duplex() {
        let stages = derive_stages_for(
            RunAllStage::AlignAndMerge,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Duplex),
        )
        .expect("align → duplex must succeed");
        assert_eq!(stages, vec![Stage::Align, Stage::Sort, Stage::Group, Stage::Duplex]);
    }

    #[test]
    fn derive_stages_align_to_codec() {
        let stages = derive_stages_for(
            RunAllStage::AlignAndMerge,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Codec),
        )
        .expect("align → codec must succeed");
        assert_eq!(stages, vec![Stage::Align, Stage::Sort, Stage::Group, Stage::Codec]);
    }

    // ── Zipper self-pair ─────────────────────────────────────────────────────

    #[test]
    fn derive_stages_zipper_to_zipper() {
        let stages = derive_stages_for(RunAllStage::Zipper, RunAllStage::Zipper, false, None)
            .expect("zipper → zipper must succeed");
        assert_eq!(stages, vec![Stage::Zipper]);
    }

    // ── Zipper → Sort ────────────────────────────────────────────────────────

    #[test]
    fn derive_stages_zipper_to_sort() {
        // Zipper output is queryname-sorted; Sort is forced.
        let stages = derive_stages_for(RunAllStage::Zipper, RunAllStage::Sort, false, None)
            .expect("zipper → sort must succeed");
        assert_eq!(stages, vec![Stage::Zipper, Stage::Sort]);
    }

    // ── Zipper → Group ───────────────────────────────────────────────────────

    #[test]
    fn derive_stages_zipper_to_group() {
        let stages = derive_stages_for(RunAllStage::Zipper, RunAllStage::Group, false, None)
            .expect("zipper → group must succeed");
        assert_eq!(stages, vec![Stage::Zipper, Stage::Sort, Stage::Group]);
    }

    // ── Zipper → consensus ───────────────────────────────────────────────────

    #[test]
    fn derive_stages_zipper_to_simplex() {
        let stages = derive_stages_for(
            RunAllStage::Zipper,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Simplex),
        )
        .expect("zipper → simplex must succeed");
        assert_eq!(stages, vec![Stage::Zipper, Stage::Sort, Stage::Group, Stage::Simplex]);
    }

    #[test]
    fn derive_stages_zipper_to_duplex() {
        let stages = derive_stages_for(
            RunAllStage::Zipper,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Duplex),
        )
        .expect("zipper → duplex must succeed");
        assert_eq!(stages, vec![Stage::Zipper, Stage::Sort, Stage::Group, Stage::Duplex]);
    }

    #[test]
    fn derive_stages_zipper_to_codec() {
        let stages = derive_stages_for(
            RunAllStage::Zipper,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Codec),
        )
        .expect("zipper → codec must succeed");
        assert_eq!(stages, vec![Stage::Zipper, Stage::Sort, Stage::Group, Stage::Codec]);
    }

    // ── Sort self-pair ────────────────────────────────────────────────────────

    #[test]
    fn derive_stages_sort_to_sort() {
        let stages = derive_stages_for(RunAllStage::Sort, RunAllStage::Sort, false, None)
            .expect("sort → sort must succeed");
        assert_eq!(stages, vec![Stage::Sort]);
    }

    // ── Sort → Group ──────────────────────────────────────────────────────────

    #[test]
    fn derive_stages_sort_to_group() {
        let stages = derive_stages_for(RunAllStage::Sort, RunAllStage::Group, false, None)
            .expect("sort → group must succeed");
        assert_eq!(stages, vec![Stage::Sort, Stage::Group]);
    }

    // ── Sort → consensus ──────────────────────────────────────────────────────

    #[test]
    fn derive_stages_sort_to_simplex() {
        let stages = derive_stages_for(
            RunAllStage::Sort,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Simplex),
        )
        .expect("sort → simplex must succeed");
        assert_eq!(stages, vec![Stage::Sort, Stage::Group, Stage::Simplex]);
    }

    #[test]
    fn derive_stages_sort_to_duplex() {
        let stages = derive_stages_for(
            RunAllStage::Sort,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Duplex),
        )
        .expect("sort → duplex must succeed");
        assert_eq!(stages, vec![Stage::Sort, Stage::Group, Stage::Duplex]);
    }

    #[test]
    fn derive_stages_sort_to_codec() {
        let stages = derive_stages_for(
            RunAllStage::Sort,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Codec),
        )
        .expect("sort → codec must succeed");
        assert_eq!(stages, vec![Stage::Sort, Stage::Group, Stage::Codec]);
    }

    // ── Group self-pair ────────────────────────────────────────────────────────

    #[test]
    fn derive_stages_group_to_group() {
        let stages = derive_stages_for(RunAllStage::Group, RunAllStage::Group, false, None)
            .expect("group → group must succeed");
        assert_eq!(stages, vec![Stage::Group]);
    }

    // ── Group → consensus ─────────────────────────────────────────────────────

    #[test]
    fn derive_stages_group_to_simplex() {
        let stages = derive_stages_for(
            RunAllStage::Group,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Simplex),
        )
        .expect("group → simplex must succeed");
        assert_eq!(stages, vec![Stage::Group, Stage::Simplex]);
    }

    #[test]
    fn derive_stages_group_to_duplex() {
        let stages = derive_stages_for(
            RunAllStage::Group,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Duplex),
        )
        .expect("group → duplex must succeed");
        assert_eq!(stages, vec![Stage::Group, Stage::Duplex]);
    }

    #[test]
    fn derive_stages_group_to_codec() {
        let stages = derive_stages_for(
            RunAllStage::Group,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Codec),
        )
        .expect("group → codec must succeed");
        assert_eq!(stages, vec![Stage::Group, Stage::Codec]);
    }

    // ── Consensus self-pairs ───────────────────────────────────────────────────
    //
    // When start_from is a consensus stage, the input BAM is already MI-tagged.
    // The standalone consensus chain builders (build_simplex_chain,
    // build_duplex_chain, build_codec_chain) handle MI-tag grouping internally
    // via GroupByMi. An explicit Group stage must NOT be added — doing so would
    // double-group an already-grouped input.

    #[test]
    fn derive_stages_simplex_self_pair() {
        // Input is already MI-tagged; Group stage must NOT be included.
        // The standalone simplex chain builder handles GroupByMi internally.
        let stages = derive_stages_for(
            RunAllStage::Consensus,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Simplex),
        )
        .expect("simplex → simplex must succeed");
        assert_eq!(stages, vec![Stage::Simplex]);
    }

    #[test]
    fn derive_stages_duplex_self_pair() {
        let stages = derive_stages_for(
            RunAllStage::Consensus,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Duplex),
        )
        .expect("duplex → duplex must succeed");
        assert_eq!(stages, vec![Stage::Duplex]);
    }

    #[test]
    fn derive_stages_codec_self_pair() {
        let stages = derive_stages_for(
            RunAllStage::Consensus,
            RunAllStage::Consensus,
            false,
            Some(RunAllMode::Codec),
        )
        .expect("codec → codec must succeed");
        assert_eq!(stages, vec![Stage::Codec]);
    }

    // ── Filter stage (post-consensus) ──────────────────────────────────────────

    #[test]
    fn derive_stages_group_to_filter_appends_filter() {
        // group → filter: consensus runs first (chosen by --consensus),
        // then Filter is appended as the terminal stage.
        let stages = derive_stages_for(
            RunAllStage::Group,
            RunAllStage::Filter,
            false,
            Some(RunAllMode::Simplex),
        )
        .expect("group → filter must succeed");
        assert_eq!(stages, vec![Stage::Group, Stage::Simplex, Stage::Filter]);
        assert_eq!(
            *stages.last().expect("non-empty"),
            Stage::Filter,
            "filter must be the terminal stage"
        );
    }

    #[test]
    fn derive_stages_consensus_start_to_filter() {
        // consensus → filter: input is a grouped/MI-tagged BAM; consensus
        // runs (no Group, per the consensus-self-pair rule), then Filter.
        let stages = derive_stages_for(
            RunAllStage::Consensus,
            RunAllStage::Filter,
            false,
            Some(RunAllMode::Duplex),
        )
        .expect("consensus → filter must succeed");
        assert_eq!(stages, vec![Stage::Duplex, Stage::Filter]);
    }

    #[test]
    fn derive_stages_filter_self_pair() {
        // filter → filter: input is a consensus BAM; only the (terminal)
        // filter stage runs. No consensus mode is required.
        let stages = derive_stages_for(RunAllStage::Filter, RunAllStage::Filter, false, None)
            .expect("filter → filter must succeed");
        assert_eq!(stages, vec![Stage::Filter]);
    }

    #[test]
    fn runall_stop_after_filter_parses() {
        use clap::Parser;
        // `--stop-after filter` is a valid CLI value (parsing succeeds);
        // the consensus → filter fusion validates at `validate_stages`
        // (see `runall_consensus_to_filter_accepted`).
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "filter",
            "--consensus",
            "simplex",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--threads",
            "1",
        ]);
        assert!(r.is_ok(), "--stop-after filter must parse: {r:?}");
    }

    #[test]
    fn runall_consensus_to_filter_accepted() {
        use clap::Parser;
        // Fusing consensus → filter is supported: the consensus stage runs
        // as an intermediate stage (its DecompressedBlock output is decoded
        // in-memory) feeding the terminal filter stage. `validate_stages`
        // must accept `--stop-after filter` from a pre-consensus start.
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "filter",
            "--consensus",
            "simplex",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--threads",
            "1",
        ])
        .expect("parses");
        assert!(
            r.validate_stages().is_ok(),
            "group → consensus → filter must pass validate_stages"
        );
        // The derived stage list must be the fused [Group, Simplex, Filter].
        let stages = r.derive_stages().expect("derive_stages");
        assert_eq!(
            stages,
            vec![
                crate::pipeline::chains::Stage::Group,
                crate::pipeline::chains::Stage::Simplex,
                crate::pipeline::chains::Stage::Filter,
            ]
        );
    }

    #[test]
    fn chain_or_delegated_io_placeholder_for_mid_chain_consensus() {
        use clap::Parser;
        // group → consensus → filter: consensus runs mid-chain, so its options
        // bag IO must NOT require --input (the chain reads from spec.source).
        // Regression for the fgumi-benchmarks "--input is required" bug on
        // `--start-from extract/sort --stop-after filter`.
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "group",
            "--stop-after",
            "filter",
            "--consensus",
            "codec",
            "--output",
            "o.bam",
            "--group::strategy",
            "adjacency",
            "--filter::min-reads",
            "1",
            "--threads",
            "1",
        ])
        .expect("parse");
        let io = r.chain_or_delegated_io().expect("mid-chain consensus must not require --input");
        assert_eq!(io.input, std::path::PathBuf::new(), "mid-chain consensus IO is a placeholder");
    }

    #[test]
    fn chain_or_delegated_io_requires_input_for_consensus_start() {
        use clap::Parser;
        // --start-from consensus: consensus IS the chain source, so --input is
        // genuinely required.
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "consensus",
            "--stop-after",
            "consensus",
            "--consensus",
            "codec",
            "--output",
            "o.bam",
            "--threads",
            "1",
        ])
        .expect("parse");
        assert!(
            r.chain_or_delegated_io().is_err(),
            "--start-from consensus must still require --input"
        );
    }

    #[test]
    fn runall_validate_stages_still_rejects_reverse_order() {
        use clap::Parser;
        // The linear-order guard must still reject a stop that precedes the
        // start (consensus → group is backwards).
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "consensus",
            "--stop-after",
            "group",
            "--consensus",
            "simplex",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--threads",
            "1",
        ])
        .expect("parses");
        assert!(
            r.validate_stages().is_err(),
            "consensus → group (reverse order) must fail validate_stages"
        );
    }

    #[test]
    fn runall_filter_self_pair_validate_stages_ok() {
        use clap::Parser;
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "filter",
            "--stop-after",
            "filter",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--threads",
            "1",
        ])
        .expect("parses");
        assert!(r.validate_stages().is_ok(), "filter self-pair must pass validate_stages");
    }

    #[test]
    fn runall_start_from_filter_parses() {
        use clap::Parser;
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "filter",
            "--stop-after",
            "filter",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--threads",
            "1",
        ]);
        assert!(r.is_ok(), "--start-from filter --stop-after filter must parse: {r:?}");
    }

    #[test]
    fn runall_filter_prefixed_flag_parses() {
        use clap::Parser;
        // The `--filter::*` Multi-options surface is exposed on runall.
        let r = RunAll::try_parse_from([
            "runall",
            "--start-from",
            "filter",
            "--stop-after",
            "filter",
            "--filter::min-reads",
            "3",
            "--filter::max-read-error-rate",
            "0.05",
            "--input",
            "x.bam",
            "--output",
            "o.bam",
            "--threads",
            "1",
        ])
        .expect("--filter::* flags must parse");
        let opts = r.filter_opts.validate().expect("filter opts validate");
        assert_eq!(opts.min_reads, vec![3]);
        assert_eq!(opts.max_read_error_rate, vec![0.05]);
    }

    // ── Group self-pair must NOT include Sort ──────────────────────────────────

    #[test]
    fn derive_stages_group_self_pair_no_sort() {
        // Group self-pair: input is sorted; no sort step needed.
        let stages =
            derive_stages_for(RunAllStage::Group, RunAllStage::Group, false, None).unwrap();
        assert!(
            !stages.contains(&Stage::Sort),
            "group self-pair must not include Sort; got {stages:?}"
        );
    }

    // ── Sort-forcing rule: guard against Align/Zipper → Group without Sort ─

    #[test]
    fn sort_forcing_rule_align_to_group_includes_sort() {
        // Regression guard: Align → Group MUST include Sort (output is
        // queryname-sorted; Group requires template-coordinate).
        let stages =
            derive_stages_for(RunAllStage::AlignAndMerge, RunAllStage::Group, false, None).unwrap();
        assert!(
            stages.contains(&Stage::Sort),
            "sort-forcing rule: Sort must be present in align → group; got {stages:?}"
        );
    }

    #[test]
    fn sort_forcing_rule_zipper_to_group_includes_sort() {
        // Same assertion for the Zipper → Group path.
        let stages =
            derive_stages_for(RunAllStage::Zipper, RunAllStage::Group, false, None).unwrap();
        assert!(
            stages.contains(&Stage::Sort),
            "sort-forcing rule: Sort must be present in zipper → group; got {stages:?}"
        );
    }

    #[test]
    fn sort_forcing_rule_align_to_consensus_includes_sort() {
        // Sort-forcing applies to the consensus stop as well.
        for mode in [RunAllMode::Simplex, RunAllMode::Duplex, RunAllMode::Codec] {
            let stages = derive_stages_for(
                RunAllStage::AlignAndMerge,
                RunAllStage::Consensus,
                false,
                Some(mode),
            )
            .unwrap_or_else(|e| panic!("align → consensus ({mode}) failed: {e}"));
            assert!(
                stages.contains(&Stage::Sort),
                "sort-forcing rule: Sort must be present in align → consensus; got {stages:?}"
            );
        }
    }
}
