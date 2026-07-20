//! Common CLI options shared across commands.
//!
//! This module provides shared argument structures that can be composed into
//! command structs using `#[command(flatten)]`.

use std::path::PathBuf;
#[cfg(feature = "consensus")]
use std::sync::Arc;

#[cfg(feature = "consensus")]
use crate::logging::OperationTimer;
use crate::validation::validate_file_exists;
use bytesize::ByteSize;
use clap::Args;
use fgumi_bam_io::is_stdin_path;
#[cfg(feature = "consensus")]
use fgumi_consensus::methylation::RefBaseProvider;
#[cfg(feature = "consensus")]
use log::info;
use noodles::sam::Header;

/// CLI argument value for `--methylation-mode`.
///
/// Maps to [`fgumi_consensus::MethylationMode`] variants (excluding `Disabled`,
/// which is represented by the absence of the flag).
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum MethylationModeArg {
    /// EM-Seq (enzymatic methyl-seq): unmethylated C is converted to T.
    /// C in read at ref-C = methylated (protected); T = unmethylated (converted).
    #[value(name = "em-seq")]
    EmSeq,
    /// TAPs/Illumina 5-base: methylated C is converted to T.
    /// C in read at ref-C = unmethylated (not a target); T = methylated (converted).
    #[value(name = "taps")]
    Taps,
}

impl From<MethylationModeArg> for fgumi_consensus::MethylationMode {
    fn from(arg: MethylationModeArg) -> Self {
        match arg {
            MethylationModeArg::EmSeq => Self::EmSeq,
            MethylationModeArg::Taps => Self::Taps,
        }
    }
}

/// Resolves an optional `--methylation-mode` CLI arg to a [`MethylationMode`].
///
/// Returns `Disabled` when `None` (flag not provided).
///
/// [`MethylationMode`]: fgumi_consensus::MethylationMode
pub fn resolve_methylation_mode(
    arg: Option<MethylationModeArg>,
) -> fgumi_consensus::MethylationMode {
    arg.map_or(fgumi_consensus::MethylationMode::Disabled, Into::into)
}

/// Methylation reference pair: reference base provider + contig name mapping.
#[cfg(feature = "consensus")]
pub type MethylationRef = Option<(
    Arc<dyn fgumi_consensus::methylation::RefBaseProvider + Send + Sync>,
    Arc<Vec<String>>,
)>;

/// Loads the reference FASTA and builds contig name mapping for methylation-aware modes.
///
/// Returns `None` if methylation mode is disabled. Errors if enabled but `reference` is `None`.
#[cfg(feature = "consensus")]
pub fn load_methylation_reference(
    methylation_mode: fgumi_consensus::MethylationMode,
    reference: &Option<PathBuf>,
    header: &Header,
) -> anyhow::Result<MethylationRef> {
    if !methylation_mode.is_enabled() {
        return Ok(None);
    }
    let mode_name = match methylation_mode {
        fgumi_consensus::MethylationMode::EmSeq => "EM-Seq",
        fgumi_consensus::MethylationMode::Taps => "TAPs",
        fgumi_consensus::MethylationMode::Disabled => unreachable!(),
    };
    let ref_path = reference
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--ref is required when --methylation-mode is set"))?;
    let ref_timer = OperationTimer::new("Loading reference FASTA");
    let reference = Arc::new(crate::reference::ReferenceReader::new(ref_path)?);
    ref_timer.log_completion(0);

    let ref_names: Vec<String> =
        header.reference_sequences().keys().map(|name| name.to_string()).collect();

    // Fail fast if any BAM header contigs are missing from the reference FASTA
    let missing_contigs: Vec<&String> =
        ref_names.iter().filter(|name| reference.sequence_for(name).is_none()).collect();
    if !missing_contigs.is_empty() {
        anyhow::bail!(
            "Reference FASTA is missing {} contig(s) from the BAM header: {}",
            missing_contigs.len(),
            missing_contigs.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", ")
        );
    }

    info!("{mode_name} mode enabled with {} reference contigs", ref_names.len());
    Ok(Some((reference, Arc::new(ref_names))))
}

/// Add a @PG record to an existing header, using the current fgumi version.
///
/// Wraps [`fgumi_bam_io::header::add_pg_record`] with the binary's version string.
pub fn add_pg_record(header: Header, command_line: &str) -> anyhow::Result<Header> {
    fgumi_bam_io::header::add_pg_record(header, crate::version::VERSION.as_str(), command_line)
}

/// Add a @PG record to a header builder, using the current fgumi version.
///
/// Wraps [`fgumi_bam_io::header::add_pg_to_builder`] with the binary's version string.
pub fn add_pg_to_builder(
    builder: noodles::sam::header::Builder,
    command_line: &str,
) -> anyhow::Result<noodles::sam::header::Builder> {
    fgumi_bam_io::header::add_pg_to_builder(builder, crate::version::VERSION.as_str(), command_line)
}

/// Validate that `header` advertises a sort order the group stage accepts, and
/// emit the accompanying info logging.
///
/// Accepts template-coordinate order always, and queryname order when
/// `allow_unmapped` is set. On rejection, bails with the order-specific
/// remediation hint. Shared by `Group::execute`'s single-threaded path and the
/// chain builder's `add_group` so the two orchestrations of the same stage
/// cannot drift on the accepted orders, the error text, or the info logging —
/// the standalone-vs-runall divergence class. The predicates themselves are the
/// shared `crate::sam::{is_template_coordinate_sorted, is_sorted}`.
pub(crate) fn require_group_input_sort(
    header: &Header,
    allow_unmapped: bool,
) -> anyhow::Result<()> {
    use crate::sam::{is_sorted, is_template_coordinate_sorted};
    use anyhow::bail;
    use log::info;
    use noodles::sam::header::record::value::map::header::sort_order::QUERY_NAME;

    let is_tc_sorted = is_template_coordinate_sorted(header);
    let is_qname_sorted = is_sorted(header, QUERY_NAME);

    if !(is_tc_sorted || allow_unmapped && is_qname_sorted) {
        if allow_unmapped {
            bail!(
                "Input BAM must be template-coordinate sorted or queryname sorted \
                when --allow-unmapped is enabled.\n\n\
                To queryname sort your BAM file, run:\n  \
                samtools sort -n input.bam -o sorted.bam"
            );
        }
        bail!(
            "Input BAM must be template-coordinate sorted (header must advertise \
            SO:unsorted, GO:query, and SS:template-coordinate).\n\n\
            To sort your BAM file, run:\n  \
            fgumi sort -i input.bam -o sorted.bam --order template-coordinate"
        );
    }

    if is_tc_sorted {
        info!("Input is template-coordinate sorted");
    } else {
        info!("Input is queryname sorted (accepted with --allow-unmapped)");
        info!("All unmapped reads will form a single position group per library/cell");
    }
    Ok(())
}

/// EM-Seq methylation-aware consensus calling options.
#[derive(Debug, Clone, Default, Args)]
pub struct EmSeqOptions {
    /// Enable EM-Seq (enzymatic methyl-seq) methylation-aware consensus calling.
    /// Requires --ref. C→T conversions at reference cytosine positions are treated
    /// as bisulfite/enzymatic conversion, and cu/ct per-base count tags
    /// and MM/ML methylation tags are emitted on consensus reads.
    #[arg(long = "em-seq", default_value_t = false, requires = "reference")]
    pub em_seq: bool,

    /// Path to the reference FASTA file (required when --em-seq is enabled)
    #[arg(long = "ref")]
    pub reference: Option<PathBuf>,
}

/// Common input/output options for commands that read a BAM and write a BAM.
#[derive(Debug, Clone, Args)]
pub struct BamIoOptions {
    /// Input BAM file
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Output BAM file
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,

    /// Wrap the input in a userspace async prefetch reader: a background
    /// thread reads ahead so disk I/O overlaps decompression/compute. Helps
    /// on slow or networked storage. Defaults to off.
    #[arg(long = "async-reader", default_value_t = false, hide = true)]
    pub async_reader: bool,
}

impl Default for BamIoOptions {
    fn default() -> Self {
        Self { input: PathBuf::new(), output: PathBuf::new(), async_reader: false }
    }
}

impl BamIoOptions {
    /// Construct a `BamIoOptions` from input and output paths, with
    /// opt-in tuning flags (e.g. `async_reader`) at their default values.
    pub fn new(input: impl Into<PathBuf>, output: impl Into<PathBuf>) -> Self {
        Self { input: input.into(), output: output.into(), async_reader: false }
    }

    /// Build [`fgumi_bam_io::PipelineReaderOpts`] from the `--async-reader`
    /// flag. When set, BAM/SAM readers wrap the input in a `PrefetchReader`
    /// background thread for overlapped I/O.
    pub fn pipeline_reader_opts(&self) -> fgumi_bam_io::PipelineReaderOpts {
        fgumi_bam_io::PipelineReaderOpts { async_reader: self.async_reader }
    }

    /// Validates that the input file exists (skipped for stdin paths).
    ///
    /// # Errors
    ///
    /// Returns an error if the input file does not exist.
    pub fn validate(&self) -> anyhow::Result<()> {
        if !is_stdin_path(&self.input) {
            validate_file_exists(&self.input, "Input BAM")?;
        }
        Ok(())
    }
}

/// Options for writing rejected reads to a separate file.
#[derive(Debug, Clone, Default, Args)]
pub struct RejectsOptions {
    /// Optional output BAM file for rejected reads
    #[arg(short = 'r', long = "rejects")]
    pub rejects: Option<PathBuf>,
}

impl RejectsOptions {
    /// Returns true if rejects tracking is enabled.
    #[must_use]
    pub fn is_enabled(&self) -> bool {
        self.rejects.is_some()
    }
}

/// Options for writing statistics to a file.
#[derive(Debug, Clone, Default, Args)]
pub struct StatsOptions {
    /// Optional output file for statistics
    #[arg(short = 's', long = "stats")]
    pub stats: Option<PathBuf>,
}

impl StatsOptions {
    /// Returns true if stats output is enabled.
    #[must_use]
    pub fn is_enabled(&self) -> bool {
        self.stats.is_some()
    }
}

/// Common options for consensus calling (simplex, duplex, codec).
#[derive(Debug, Clone, Args)]
pub struct ConsensusCallingOptions {
    /// Phred-scaled error rate prior to UMI integration
    #[arg(short = '1', long = "error-rate-pre-umi", default_value = "45")]
    pub error_rate_pre_umi: u8,

    /// Phred-scaled error rate post UMI integration
    #[arg(short = '2', long = "error-rate-post-umi", default_value = "40")]
    pub error_rate_post_umi: u8,

    /// Minimum base quality in raw reads to use for consensus
    #[arg(short = 'm', long = "min-input-base-quality", default_value = "10")]
    pub min_input_base_quality: u8,

    /// Produce per-base tags (cd, ce) in addition to per-read tags
    #[arg(short = 'B', long = "output-per-base-tags", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub output_per_base_tags: bool,

    /// Quality-trim reads before consensus calling (removes low-quality bases from ends)
    #[arg(long = "trim", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub trim: bool,

    /// Minimum consensus base quality (output consensus bases below this are masked to N)
    #[arg(long = "min-consensus-base-quality", default_value = "2")]
    pub min_consensus_base_quality: u8,
}

impl Default for ConsensusCallingOptions {
    fn default() -> Self {
        Self {
            error_rate_pre_umi: 45,
            error_rate_post_umi: 40,
            min_input_base_quality: 10,
            output_per_base_tags: true,
            trim: false,
            min_consensus_base_quality: 2,
        }
    }
}

impl ConsensusCallingOptions {
    /// Maximum valid Phred score (Illumina 1.8+ uses 0-41, but we allow up to 93).
    const MAX_PHRED: u8 = 93;

    /// Validates the consensus calling options.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Any Phred quality value exceeds `MAX_PHRED` (93)
    /// - `min_consensus_base_quality` is less than 2 (`MIN_PHRED`)
    pub fn validate(&self) -> anyhow::Result<()> {
        use anyhow::bail;

        if self.error_rate_pre_umi > Self::MAX_PHRED {
            bail!(
                "error-rate-pre-umi ({}) exceeds maximum Phred score ({})",
                self.error_rate_pre_umi,
                Self::MAX_PHRED
            );
        }
        if self.error_rate_post_umi > Self::MAX_PHRED {
            bail!(
                "error-rate-post-umi ({}) exceeds maximum Phred score ({})",
                self.error_rate_post_umi,
                Self::MAX_PHRED
            );
        }
        if self.min_input_base_quality > Self::MAX_PHRED {
            bail!(
                "min-input-base-quality ({}) exceeds maximum Phred score ({})",
                self.min_input_base_quality,
                Self::MAX_PHRED
            );
        }
        if self.min_consensus_base_quality < 2 {
            bail!(
                "min-consensus-base-quality ({}) must be at least 2 (MIN_PHRED)",
                self.min_consensus_base_quality
            );
        }
        if self.min_consensus_base_quality > Self::MAX_PHRED {
            bail!(
                "min-consensus-base-quality ({}) exceeds maximum Phred score ({})",
                self.min_consensus_base_quality,
                Self::MAX_PHRED
            );
        }

        Ok(())
    }
}

/// Options for read group and read name prefix in consensus output.
#[derive(Debug, Clone, Args)]
pub struct ReadGroupOptions {
    /// Prefix for consensus read names
    #[arg(short = 'p', long = "read-name-prefix")]
    pub read_name_prefix: Option<String>,

    /// Read group ID for consensus reads
    #[arg(short = 'R', long = "read-group-id", default_value = "A")]
    pub read_group_id: String,
}

impl ReadGroupOptions {
    /// Returns the configured read name prefix, or derives one from the header.
    #[must_use]
    pub fn prefix_or_from_header(&self, header: &noodles::sam::Header) -> String {
        self.read_name_prefix
            .clone()
            .unwrap_or_else(|| crate::consensus_caller::make_prefix_from_header(header))
    }
}

impl Default for ReadGroupOptions {
    fn default() -> Self {
        Self { read_name_prefix: None, read_group_id: "A".to_string() }
    }
}

/// Options for overlapping bases consensus calling.
#[derive(Debug, Clone, Args)]
pub struct OverlappingConsensusOptions {
    /// Consensus call overlapping bases in read pairs before UMI consensus calling
    #[arg(long = "consensus-call-overlapping-bases", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub consensus_call_overlapping_bases: bool,
}

impl Default for OverlappingConsensusOptions {
    fn default() -> Self {
        Self { consensus_call_overlapping_bases: true }
    }
}

impl OverlappingConsensusOptions {
    /// Returns true if overlapping consensus calling is enabled.
    #[must_use]
    pub fn is_enabled(&self) -> bool {
        self.consensus_call_overlapping_bases
    }
}

/// Threading mode for parallel processing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThreadingMode {
    /// Single-threaded mode (no parallelism).
    SingleThreaded,
    /// Thread cap mode: never exceed N total threads.
    /// Suitable for HPC/shared systems where predictable resource usage matters.
    Threads(usize),
}

impl ThreadingMode {
    /// Returns true if this mode enables parallel processing.
    #[must_use]
    pub fn is_parallel(&self) -> bool {
        !matches!(self, Self::SingleThreaded)
    }

    /// Returns the target thread count for this mode.
    #[must_use]
    pub fn num_threads(&self) -> usize {
        match self {
            Self::SingleThreaded => 1,
            Self::Threads(n) => *n,
        }
    }
}

/// Common threading options for parallel BAM processing.
///
/// The `--threads N` option caps total threads at N.
/// Use this on HPC/shared systems where you need predictable resource usage.
///
/// # Examples
///
/// ```bash
/// fgumi group --threads 8 ...
/// # Uses up to 8 threads with the round-robin pipeline scheduler
/// ```
#[derive(Debug, Clone, Args)]
pub struct ThreadingOptions {
    /// Number of threads for the multi-threaded pipeline.
    ///
    /// If not specified, uses a single-threaded fast path optimized for
    /// simple streaming. When specified (even with --threads 1), uses the
    /// typed-step work-stealing pipeline.
    #[arg(long = "threads")]
    pub threads: Option<usize>,
}

/// Default `--deadlock-timeout` in seconds (0 would disable detection). Shared
/// so producers that build a [`SchedulerOptions`] without going through clap
/// (e.g. the standalone `fgumi sort` command, which exposes no
/// `--deadlock-timeout` flag) stay aligned with the CLI default rather than
/// hardcoding a bare literal that silently drifts.
pub const DEFAULT_DEADLOCK_TIMEOUT_SECS: u64 = 10;

/// Pipeline debugging/diagnostics options: statistics output, deadlock
/// timeout, and deadlock recovery. Flattened into every pipeline command.
///
/// (Named `SchedulerOptions` for historical reasons — the `--scheduler`
/// strategy flag it once carried was removed in the issue #330 migration,
/// since the typed-step dispatch model has no pluggable strategy. It survives
/// only as a hidden, deprecated no-op alias so legacy invocations degrade to a
/// warning rather than a hard clap error — see the `scheduler` field.)
#[derive(Debug, Clone, Default, Args)]
pub struct SchedulerOptions {
    /// DEPRECATED: no-op. The pluggable scheduler strategy was removed in the
    /// issue #330 typed-step dispatch migration; the typed-step model has no
    /// strategy to select. Kept as a hidden alias that accepts (and ignores)
    /// any value so scripts passing the old `--scheduler <strategy>` warn
    /// instead of hard-failing on an unknown argument (AUDIT-002).
    #[arg(long = "scheduler", hide = true)]
    pub scheduler: Option<String>,

    /// Print detailed pipeline statistics at completion.
    ///
    /// Shows per-step timing, throughput, contention metrics, and
    /// per-thread work distribution.
    #[arg(long = "pipeline-stats", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool, hide = true)]
    pub pipeline_stats: bool,

    /// Per-edge instrumentation level: off | summary | timeline | deep.
    ///
    /// `summary` adds per-edge throughput/occupancy/latency + a bottleneck
    /// verdict to the end-of-run report; `timeline` also writes a per-tick TSV
    /// (see `--pipeline-trace-out`); `deep` adds direct dwell/park-time latency.
    /// Diagnostic only — throughput is slightly depressed under tracing; confirm
    /// final wall/RSS with the flag off. Overridable via `FGUMI_PIPELINE_TRACE`.
    #[arg(long = "pipeline-trace", default_value = "off", value_parser = ["off", "summary", "timeline", "deep"], hide = true)]
    pub pipeline_trace: String,

    /// Path for the `--pipeline-trace timeline` per-tick TSV
    /// (default `pipeline-trace.tsv` in the working directory).
    #[arg(long = "pipeline-trace-out", hide = true)]
    pub pipeline_trace_out: Option<std::path::PathBuf>,

    /// Timeout in seconds for deadlock detection (default: 10, 0 = disabled).
    ///
    /// When no progress is made for this duration, a warning is logged with
    /// diagnostic info (queue depths, memory usage, per-queue timestamps). If
    /// the stall persists with work still in flight, the run fails fast with a
    /// timeout error after 6× this duration (~60s at the default) — long enough
    /// that a single large dispatch (busy-locus group, large sort merge) is not
    /// mistaken for a wedge.
    #[arg(long = "deadlock-timeout", default_value_t = DEFAULT_DEADLOCK_TIMEOUT_SECS, hide = true)]
    pub deadlock_timeout: u64,

    /// Enable automatic deadlock recovery (default: false, detection only).
    ///
    /// Uses progressive doubling: 2x -> 4x -> unbind, with restoration
    /// after 30s of sustained progress.
    #[arg(long = "deadlock-recover", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool, hide = true)]
    pub deadlock_recover: bool,
}

impl SchedulerOptions {
    /// Returns true if pipeline stats should be collected and printed.
    #[must_use]
    pub fn collect_stats(&self) -> bool {
        self.pipeline_stats
    }

    /// Resolve the instrumentation level from `--pipeline-trace`, with the
    /// `FGUMI_PIPELINE_TRACE` env var taking precedence (so standalone commands
    /// can enable it without flattening the flag). Unknown values → `Off`.
    #[must_use]
    pub fn instrumentation_level(&self) -> crate::pipeline::core::builder::InstrumentationLevel {
        let raw =
            std::env::var("FGUMI_PIPELINE_TRACE").unwrap_or_else(|_| self.pipeline_trace.clone());
        parse_instrumentation_level(&raw)
    }

    /// Returns the deadlock detection timeout in seconds (0 = disabled).
    #[must_use]
    pub fn deadlock_timeout_secs(&self) -> u64 {
        self.deadlock_timeout
    }

    /// Returns true if automatic deadlock recovery is enabled.
    #[must_use]
    pub fn deadlock_recover_enabled(&self) -> bool {
        self.deadlock_recover
    }
}

impl ThreadingOptions {
    /// Default batch size for MI group processing.
    ///
    /// This determines how many MI groups are processed together in a single batch
    /// when using parallel processing. Smaller values reduce latency but increase
    /// synchronization overhead; larger values improve throughput but may cause
    /// uneven work distribution.
    pub const DEFAULT_BATCH_SIZE: usize = 100;

    /// Creates threading options with N threads (uses pipeline).
    #[must_use]
    pub fn new(threads: usize) -> Self {
        Self { threads: Some(threads) }
    }

    /// Creates threading options with no threads specified (uses single-threaded fast path).
    #[must_use]
    pub fn none() -> Self {
        Self { threads: None }
    }

    /// Returns the threading mode based on CLI options.
    ///
    /// - `None` -> `SingleThreaded` (fast path, no pipeline)
    /// - `Some(n)` -> `Threads(n)` (uses pipeline, even when n=1)
    #[must_use]
    pub fn mode(&self) -> ThreadingMode {
        match self.threads {
            None => ThreadingMode::SingleThreaded,
            Some(n) => ThreadingMode::Threads(n),
        }
    }

    /// Returns the number of threads.
    #[must_use]
    pub fn num_threads(&self) -> usize {
        self.mode().num_threads()
    }

    /// Returns true if parallel processing should be used.
    #[must_use]
    pub fn is_parallel(&self) -> bool {
        self.mode().is_parallel()
    }

    /// Returns true if running in single-threaded mode.
    #[must_use]
    pub fn is_single_threaded(&self) -> bool {
        matches!(self.mode(), ThreadingMode::SingleThreaded)
    }

    /// Returns the queue length for parallel processing channels.
    ///
    /// The queue length determines how many batches can be buffered between
    /// the reader and worker threads. A value of `2 * num_threads` provides
    /// good overlap between I/O and compute while limiting memory usage.
    #[must_use]
    pub fn queue_len(&self) -> usize {
        self.num_threads() * 2
    }

    /// Returns a log message describing the threading configuration.
    #[must_use]
    pub fn log_message(&self) -> String {
        match self.mode() {
            ThreadingMode::SingleThreaded => "Single-threaded mode".to_string(),
            ThreadingMode::Threads(n) => format!("Using {n} threads"),
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
// Shared memory-budget types
//
// `--max-memory` / `--memory-reserve` are used both by `fgumi sort` (to bound
// the in-memory sort buffer before spilling to disk) and by every pipeline
// command via [`QueueMemoryOptions`] (to bound the inter-stage pipeline queue).
// The types, parsers, and resolution logic live here so the two surfaces stay
// in lockstep.
//////////////////////////////////////////////////////////////////////////////

// The memory-limit types, parsers, and resolution logic now live in the
// `fgumi-cli-common` crate. Re-export them so `crate::commands::common::*`
// paths keep resolving AND so there is ONE `MemoryLimit`/`MemoryReserve`
// type shared with `fgumi-sort-cli`'s `SortOptions` (type unification).
pub use fgumi_cli_common::{
    MIN_MEMORY_PER_THREAD, MemoryLimit, MemoryReserve, parse_bool, parse_memory,
    parse_memory_reserve, resolve_memory_budget, resolve_reserve,
};

/// Options for pipeline queue memory limits.
///
/// Bounds the memory held in the pipeline's inter-stage queues so a command
/// self-throttles instead of OOM-ing. Shared (via `#[command(flatten)]`) by
/// every multi-threaded command, so the `--max-memory` surface matches
/// `fgumi sort`.
///
/// `--max-memory` is the dominant *controllable* consumer; it is a budget, not
/// a hard RSS cap — a single pathological position group is still processed
/// whole, and each worker has transient working-set memory on top of the queue.
#[derive(Debug, Clone, Args)]
pub struct QueueMemoryOptions {
    /// Maximum memory for the pipeline queues.
    ///
    /// When unset (the common case) the default is 768 MiB per thread — EXCEPT
    /// a single-threaded run (`--threads 1` / `--threads` omitted), which uses
    /// a small fixed total (a streaming footprint), since a lone worker doesn't
    /// need per-thread parallel buffering. Pass "auto" to detect host memory
    /// (cgroup-aware) and subtract --memory-reserve, so the budget shrinks to
    /// fit the host and the command self-throttles instead of OOM-ing. Explicit
    /// values like "512MiB", "1GiB", "4GiB" are per-thread when --memory-per-thread
    /// is enabled (default); plain numbers are MiB. Note bare "M"/"G" are decimal
    /// (1000ⁿ); "MiB"/"GiB" are binary (1024ⁿ). Mirrors `fgumi sort`'s --max-memory.
    #[arg(long = "max-memory", value_parser = parse_memory)]
    pub max_memory: Option<MemoryLimit>,

    /// Memory to reserve for other processes when --max-memory=auto.
    ///
    /// "auto" (default) reserves min(10 GiB, 50% of host memory). Explicit
    /// values like "10G", "8GiB" set a fixed reservation. Ignored when
    /// --max-memory is an explicit value.
    #[arg(long = "memory-reserve", default_value = "auto", value_parser = parse_memory_reserve)]
    pub memory_reserve: MemoryReserve,

    /// Scale the memory limit by thread count.
    ///
    /// When enabled (default), --max-memory is per thread, so total memory =
    /// `max_memory` × threads. Disable for a fixed total budget regardless of
    /// thread count (recommended on fixed-RAM hosts).
    #[arg(long = "memory-per-thread", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub memory_per_thread: bool,

    /// DEPRECATED: use --max-memory instead. Kept as a hidden alias; when set it
    /// overrides --max-memory.
    #[arg(long = "queue-memory", hide = true)]
    pub queue_memory: Option<String>,

    /// DEPRECATED: use --memory-per-thread instead. Hidden alias that overrides
    /// --memory-per-thread when set.
    #[arg(long = "queue-memory-per-thread", hide = true, num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub queue_memory_per_thread: Option<bool>,

    /// DEPRECATED: use --max-memory instead. Hidden alias for a fixed total
    /// limit in megabytes.
    #[arg(long = "queue-memory-limit-mb", hide = true)]
    pub queue_memory_limit_mb: Option<u64>,
}

impl Default for QueueMemoryOptions {
    fn default() -> Self {
        Self {
            // `None` = unset → the computed default applies (lean 64 MiB total
            // when single-threaded, else 768 MiB/thread). Matches the CLI's
            // no-flag case.
            max_memory: None,
            memory_reserve: MemoryReserve::Auto,
            memory_per_thread: true,
            queue_memory: None,
            queue_memory_per_thread: None,
            queue_memory_limit_mb: None,
        }
    }
}

impl QueueMemoryOptions {
    /// Default per-thread queue-memory budget when `--max-memory` is unset
    /// (multi-threaded runs): 768 MiB × threads. Byte-identical to the
    /// historical `--queue-memory 768` default.
    const DEFAULT_PER_THREAD_BYTES: usize = 768 * 1024 * 1024;

    /// Default TOTAL queue-memory budget when `--max-memory` is unset AND the
    /// run is single-threaded. A lone worker streams one batch per stage, so a
    /// per-thread parallel budget is wasteful; this small total keeps the
    /// (off-budget) reorder stash and transport queues at a streaming
    /// footprint. Overridden by any explicit `--max-memory` / `--queue-memory`.
    const LEAN_SINGLE_THREAD_TOTAL_BYTES: u64 = 64 * 1024 * 1024;

    /// Resolve the effective `(limit, reserve, per_thread)` triple, honoring the
    /// deprecated `--queue-memory*` aliases over the modern `--max-memory*`
    /// flags. Pure: no warnings, no system queries.
    ///
    /// `limit` is `None` only when neither `--max-memory` nor a deprecated
    /// `--queue-memory*` alias was given, signalling the caller to apply the
    /// thread-count-aware default (lean total at `threads == 1`, else
    /// 768 MiB/thread).
    ///
    /// # Errors
    /// Returns an error if a deprecated `--queue-memory` string fails to parse.
    fn resolve_spec(&self) -> anyhow::Result<(Option<MemoryLimit>, MemoryReserve, bool)> {
        // The deprecated `--queue-memory-per-thread` alias overrides
        // `--memory-per-thread` when set, even on its own (i.e. without
        // `--queue-memory`); otherwise passing it alone would be a silent no-op.
        let per_thread = self.queue_memory_per_thread.unwrap_or(self.memory_per_thread);

        // Legacy fixed-total flag wins if present (matches prior behavior).
        if let Some(legacy_mb) = self.queue_memory_limit_mb {
            let bytes = usize::try_from(legacy_mb)
                .ok()
                .and_then(|mb| mb.checked_mul(1024 * 1024))
                .ok_or_else(|| anyhow::anyhow!("--queue-memory-limit-mb too large: {legacy_mb}"))?;
            return Ok((Some(MemoryLimit::Fixed(bytes)), MemoryReserve::Auto, false));
        }

        // Deprecated --queue-memory overrides --max-memory when set.
        if let Some(queue_memory) = &self.queue_memory {
            let limit = parse_memory(queue_memory).map_err(|e| {
                anyhow::anyhow!("Failed to parse --queue-memory {queue_memory}: {e}")
            })?;
            return Ok((Some(limit), self.memory_reserve, per_thread));
        }

        Ok((self.max_memory, self.memory_reserve, per_thread))
    }

    /// Emit a one-time deprecation warning for any legacy flag in use.
    fn warn_deprecations(&self) {
        if self.queue_memory_limit_mb.is_some() {
            log::warn!("DEPRECATED: --queue-memory-limit-mb; use --max-memory <N> instead.");
        } else if self.queue_memory.is_some() {
            log::warn!(
                "DEPRECATED: --queue-memory / --queue-memory-per-thread; use --max-memory / --memory-per-thread instead."
            );
        } else if self.queue_memory_per_thread.is_some() {
            log::warn!("DEPRECATED: --queue-memory-per-thread; use --memory-per-thread instead.");
        }
    }

    /// Resolve the total queue memory budget in bytes for `num_threads`.
    ///
    /// When neither `--max-memory` nor a deprecated `--queue-memory*` alias is
    /// set, the budget is the lean fixed total at `num_threads == 1` (a lone
    /// worker streams) and 768 MiB/thread otherwise. Under `--max-memory auto`
    /// the budget is detected from (cgroup-aware) host memory minus
    /// `--memory-reserve`, so it shrinks to fit the host. Under an explicit
    /// value it is `max_memory` (× threads when per-thread).
    ///
    /// # Errors
    /// Returns an error if `num_threads` is 0, a deprecated `--queue-memory`
    /// string fails to parse, or the multiplication overflows.
    pub fn calculate_memory_limit(&self, num_threads: usize) -> anyhow::Result<u64> {
        self.warn_deprecations();
        let (limit, reserve, per_thread) = self.resolve_spec()?;
        let limit = if let Some(limit) = limit {
            limit
        } else {
            // Unset: a lone worker streams, so use the lean fixed total at
            // threads == 1; otherwise the historical 768 MiB/thread default.
            if num_threads == 1 {
                return Ok(Self::LEAN_SINGLE_THREAD_TOTAL_BYTES);
            }
            MemoryLimit::Fixed(Self::DEFAULT_PER_THREAD_BYTES)
        };
        let bytes = resolve_memory_budget(limit, reserve, num_threads, per_thread)?;
        Ok(bytes as u64)
    }

    /// Logs the resolved memory configuration.
    ///
    /// # Arguments
    /// * `num_threads` - Number of threads used for the calculation.
    /// * `total_memory` - The resolved total budget from `calculate_memory_limit`.
    pub fn log_memory_config(&self, num_threads: usize, total_memory: u64) {
        let per_thread = self.resolve_spec().map(|(_, _, pt)| pt).unwrap_or(true);
        if per_thread && num_threads > 1 {
            log::info!(
                "Queue memory budget: {} total ({}/thread × {} threads)",
                ByteSize(total_memory),
                ByteSize(total_memory / num_threads as u64),
                num_threads
            );
        } else {
            log::info!("Queue memory budget: {} total", ByteSize(total_memory));
        }
    }
}

// `parse_bool` is re-exported from `fgumi_cli_common` above.

// Output compression options now live in `fgumi-cli-common`; re-export so
// `crate::commands::common::CompressionOptions` keeps resolving and is the
// same type as `fgumi-sort-cli`'s `Sort::compression`.
pub use fgumi_cli_common::CompressionOptions;

// Re-export from the library crate for backward compatibility.
// Now only the module's own tests reference `detect_total_memory` (the memory
// resolvers that used it moved to `fgumi-cli-common`), so gate the import to
// test builds to avoid an unused-import warning in the library build.
#[cfg(test)]
pub(crate) use crate::system::detect_total_memory;
pub use crate::validation::parse_memory_size;

/// Log warnings for `SchedulerOptions` / `QueueMemoryOptions` flags that the
/// typed-step pipeline doesn't honor. Called from each command's multi-
/// threaded dispatch so users see why their flags might appear to be
/// ignored. `--pipeline-stats` is honored separately via
/// `attach_new_pipeline_stats`; this only warns about the others.
pub fn warn_unwired_pipeline_flags(
    scheduler_opts: &SchedulerOptions,
    queue_memory: &QueueMemoryOptions,
) {
    // --scheduler: removed in the issue #330 typed-step dispatch migration
    // (no pluggable strategy to select). Accepted as a hidden no-op alias so
    // legacy scripts warn rather than hard-failing on an unknown argument.
    if scheduler_opts.scheduler.is_some() {
        log::warn!(
            "DEPRECATED: --scheduler is a no-op and is ignored; the pluggable \
             scheduler strategy was removed in the typed-step pipeline migration \
             (#330). Remove it from your invocation."
        );
    }
    // --deadlock-recover: the legacy progressive-doubling recovery
    // addressed a failure mode (a worker pinned on a stuck step under
    // legacy's static scheduler) that doesn't exist in the typed-step
    // dispatch model — every worker round-robins through every step
    // each iteration, so there's nothing to "recover" from.
    if scheduler_opts.deadlock_recover_enabled() {
        log::info!(
            "--deadlock-recover has no effect in the typed-step pipeline: \
             the dispatch model round-robins all workers across all steps, \
             so the failure mode legacy's progressive recovery addressed \
             does not occur"
        );
    }
    // `--queue-memory` is now honored: the total bytes flow into
    // `PipelineConfig::queue_memory_total`, which both seeds the
    // initial per-queue budget AND enables the rebalancer that
    // shifts budget between consistently-full / consistently-empty
    // queues at runtime. No warning needed.
    let _ = queue_memory; // signal intentional use; dead-code lint dampener
}

/// If `--pipeline-stats` is set, construct a `PipelineStats` and attach it
/// to the given `PipelineConfig` so the worker loop populates per-step
/// counters. Returns the optional `Arc<PipelineStats>` so the caller can
/// print it after `pipeline.run` returns; `None` when stats are disabled.
pub fn attach_new_pipeline_stats(
    scheduler_opts: &SchedulerOptions,
    pipeline: &crate::pipeline::core::builder::Pipeline,
    config: &mut crate::pipeline::core::builder::PipelineConfig,
) -> Option<std::sync::Arc<crate::pipeline::core::runtime::stats::PipelineStats>> {
    if !scheduler_opts.collect_stats() {
        return None;
    }
    let stats = pipeline.stats();
    config.stats = Some(std::sync::Arc::clone(&stats));
    Some(stats)
}

/// Print the new-pipeline `PipelineStats` snapshot to the log if any
/// were collected. Pairs with `attach_new_pipeline_stats`.
pub fn log_new_pipeline_stats(
    stats: Option<std::sync::Arc<crate::pipeline::core::runtime::stats::PipelineStats>>,
) {
    if let Some(stats) = stats {
        let snapshot = stats.snapshot();
        log::info!("=== Pipeline statistics ===");
        for line in format!("{snapshot}").lines() {
            log::info!("{line}");
        }
    }
}

/// Run a new-pipeline `Pipeline`, honoring `--pipeline-stats`,
/// `--deadlock-timeout`, and `--queue-memory`. The stats snapshot is
/// logged after the run completes when `--pipeline-stats` is set.
/// When `--deadlock-timeout` is non-zero, stats are auto-attached
/// even if `--pipeline-stats` was off (the monitor reads stats
/// counters to detect stalls). When `--queue-memory` is provided, the
/// total bytes are passed to `PipelineConfig::queue_memory_total` so
/// the rebalancer can spread the budget across byte-bounded queues.
pub fn run_new_pipeline(
    pipeline: crate::pipeline::core::builder::Pipeline,
    num_threads: usize,
    scheduler_opts: &SchedulerOptions,
    queue_memory: &QueueMemoryOptions,
) -> anyhow::Result<()> {
    let mut config = crate::pipeline::core::builder::PipelineConfig {
        threads: num_threads,
        ..Default::default()
    };
    config.deadlock_timeout_secs = scheduler_opts.deadlock_timeout_secs();
    config.queue_memory_total = Some(queue_memory.calculate_memory_limit(num_threads)?);
    let user_wants_stats = scheduler_opts.collect_stats();
    let monitor_needs_stats = config.deadlock_timeout_secs > 0;
    let stats = if user_wants_stats || monitor_needs_stats {
        let s = pipeline.stats();
        config.stats = Some(std::sync::Arc::clone(&s));
        Some(s)
    } else {
        None
    };
    let result = pipeline.run(config).map_err(|e| anyhow::anyhow!("Pipeline::run: {e:?}"));
    if user_wants_stats {
        log_new_pipeline_stats(stats);
    }
    result
}

/// Map a raw `--pipeline-trace` / `FGUMI_PIPELINE_TRACE` string to an
/// [`InstrumentationLevel`](crate::pipeline::core::builder::InstrumentationLevel).
/// Unknown values → `Off`. Pure (no env read), so it is unit-testable in
/// isolation from ambient process state.
fn parse_instrumentation_level(raw: &str) -> crate::pipeline::core::builder::InstrumentationLevel {
    use crate::pipeline::core::builder::InstrumentationLevel as L;
    match raw {
        "summary" => L::Summary,
        "timeline" => L::Timeline,
        "deep" => L::Deep,
        _ => L::Off,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::core::builder::InstrumentationLevel;

    /// Build a SAM header with the given `@HD` line plus one `@SQ`, for the
    /// `require_group_input_sort` matrix below.
    fn header(hd: &str) -> Header {
        format!("{hd}\n@SQ\tSN:chr1\tLN:1000\n").parse().expect("parse SAM header")
    }

    /// B2 (audit): `require_group_input_sort` is the single source of the group
    /// stage's sort-order precondition, shared by the single-threaded command path
    /// and the chain builder's `add_group`. Pin its ACCEPT matrix — template-coordinate
    /// is accepted regardless of `--allow-unmapped`; queryname only with it — so the two
    /// callers cannot drift.
    #[rstest]
    #[case::tc_no_allow("@HD\tVN:1.6\tSO:unsorted\tGO:query\tSS:template-coordinate", false)]
    #[case::tc_allow("@HD\tVN:1.6\tSO:unsorted\tGO:query\tSS:template-coordinate", true)]
    #[case::queryname_allow("@HD\tVN:1.6\tSO:queryname", true)]
    fn require_group_input_sort_accepts(#[case] hd: &str, #[case] allow_unmapped: bool) {
        assert!(require_group_input_sort(&header(hd), allow_unmapped).is_ok());
    }

    /// B2 (audit): the REJECT side of the same precondition — the diagnostic differs by
    /// `--allow-unmapped` (queryname without the flag points at template-coordinate;
    /// coordinate with the flag mentions the queryname escape hatch).
    #[rstest]
    #[case::queryname_no_allow(
        "@HD\tVN:1.6\tSO:queryname",
        false,
        "must be template-coordinate sorted"
    )]
    #[case::coordinate_allow("@HD\tVN:1.6\tSO:coordinate", true, "queryname sorted")]
    fn require_group_input_sort_rejects(
        #[case] hd: &str,
        #[case] allow_unmapped: bool,
        #[case] expected_substr: &str,
    ) {
        let err = require_group_input_sort(&header(hd), allow_unmapped).unwrap_err().to_string();
        assert!(err.contains(expected_substr), "got: {err}");
    }

    // Pure mapping — independent of the ambient `FGUMI_PIPELINE_TRACE`, so every
    // case always runs (no skip-on-set that could silently no-op if another test
    // in the same binary sets the var first).
    #[rstest]
    #[case("", InstrumentationLevel::Off)] // empty → off
    #[case("off", InstrumentationLevel::Off)]
    #[case("bogus", InstrumentationLevel::Off)] // unknown → off
    #[case("summary", InstrumentationLevel::Summary)]
    #[case("timeline", InstrumentationLevel::Timeline)]
    #[case("deep", InstrumentationLevel::Deep)]
    fn parse_instrumentation_level_maps_trace_strings(
        #[case] input: &str,
        #[case] expected: InstrumentationLevel,
    ) {
        assert_eq!(parse_instrumentation_level(input), expected);
    }

    #[test]
    fn test_none_is_single_threaded() {
        let opts = ThreadingOptions::none();
        assert!(opts.is_single_threaded());
        assert!(!opts.is_parallel());
        assert_eq!(opts.mode(), ThreadingMode::SingleThreaded);
        assert_eq!(opts.threads, None);
    }

    #[test]
    fn test_mode_detection() {
        // None -> SingleThreaded (fast path)
        assert_eq!(ThreadingOptions::none().mode(), ThreadingMode::SingleThreaded);
        // Some(1) -> Threads(1) (uses pipeline, even with 1 thread)
        assert_eq!(ThreadingOptions::new(1).mode(), ThreadingMode::Threads(1));
        // Some(8) -> Threads(8)
        assert_eq!(ThreadingOptions::new(8).mode(), ThreadingMode::Threads(8));
    }

    #[test]
    fn test_num_threads() {
        assert_eq!(ThreadingOptions::none().num_threads(), 1);
        assert_eq!(ThreadingOptions::new(1).num_threads(), 1);
        assert_eq!(ThreadingOptions::new(8).num_threads(), 8);
    }

    #[test]
    fn test_queue_len() {
        assert_eq!(ThreadingOptions::new(1).queue_len(), 2);
        assert_eq!(ThreadingOptions::new(8).queue_len(), 16);
    }

    #[test]
    fn test_log_message() {
        let opts = ThreadingOptions::new(8);
        let msg = opts.log_message();
        assert!(msg.contains("8 threads"));

        let opts = ThreadingOptions::none();
        let msg = opts.log_message();
        assert!(msg.contains("Single-threaded"));
    }

    #[test]
    fn test_new_uses_pipeline() {
        // new(1) should use pipeline (Threads mode), not fast path
        let opts = ThreadingOptions::new(1);
        assert!(!opts.is_single_threaded());
        assert!(opts.is_parallel());
        assert_eq!(opts.threads, Some(1));
    }

    // ========== Tests for option struct validation ==========

    #[test]
    fn test_consensus_calling_options_validate_defaults() {
        let opts = ConsensusCallingOptions::default();
        assert!(opts.validate().is_ok());
    }

    #[test]
    fn test_consensus_calling_options_validate_valid() {
        let opts = ConsensusCallingOptions {
            error_rate_pre_umi: 45,
            error_rate_post_umi: 40,
            min_input_base_quality: 10,
            output_per_base_tags: true,
            trim: false,
            min_consensus_base_quality: 13,
        };
        assert!(opts.validate().is_ok());
    }

    #[test]
    fn test_consensus_calling_options_validate_error_rate_pre_umi_too_high() {
        let opts = ConsensusCallingOptions {
            error_rate_pre_umi: 94, // Exceeds MAX_PHRED
            ..ConsensusCallingOptions::default()
        };
        let err = opts.validate().unwrap_err();
        assert!(err.to_string().contains("error-rate-pre-umi"));
    }

    #[test]
    fn test_consensus_calling_options_validate_error_rate_post_umi_too_high() {
        let opts = ConsensusCallingOptions {
            error_rate_post_umi: 94, // Exceeds MAX_PHRED
            ..ConsensusCallingOptions::default()
        };
        let err = opts.validate().unwrap_err();
        assert!(err.to_string().contains("error-rate-post-umi"));
    }

    #[test]
    fn test_consensus_calling_options_validate_min_consensus_too_low() {
        let opts = ConsensusCallingOptions {
            min_consensus_base_quality: 1, // Below MIN_PHRED
            ..ConsensusCallingOptions::default()
        };
        let err = opts.validate().unwrap_err();
        assert!(err.to_string().contains("min-consensus-base-quality"));
    }

    #[test]
    fn test_consensus_calling_options_validate_min_consensus_at_min() {
        let opts = ConsensusCallingOptions {
            min_consensus_base_quality: 2, // Exactly MIN_PHRED
            ..ConsensusCallingOptions::default()
        };
        assert!(opts.validate().is_ok());
    }

    // ========== Tests for SchedulerOptions ==========

    #[test]
    fn test_scheduler_options_default() {
        let opts = SchedulerOptions::default();
        assert!(!opts.collect_stats());
        assert_eq!(opts.deadlock_timeout_secs(), 0);
    }

    #[test]
    fn test_scheduler_options_collect_stats() {
        let opts = SchedulerOptions {
            pipeline_stats: true,
            deadlock_timeout: 10,
            deadlock_recover: false,
            scheduler: None,
            ..Default::default()
        };
        assert!(opts.collect_stats());
    }

    #[test]
    fn test_scheduler_options_deadlock_timeout() {
        let opts = SchedulerOptions {
            pipeline_stats: false,
            deadlock_timeout: 30,
            deadlock_recover: false,
            scheduler: None,
            ..Default::default()
        };
        assert_eq!(opts.deadlock_timeout_secs(), 30);
    }

    #[test]
    fn test_scheduler_options_deadlock_recover() {
        let opts = SchedulerOptions {
            pipeline_stats: false,
            deadlock_timeout: 10,
            deadlock_recover: true,
            scheduler: None,
            ..Default::default()
        };
        assert!(opts.deadlock_recover_enabled());
    }

    #[test]
    fn test_deprecated_scheduler_flag_is_accepted_as_noop() {
        use clap::Parser;

        // A minimal parser that flattens SchedulerOptions, mirroring how every
        // pipeline command embeds it. The deprecated `--scheduler <value>` must
        // parse (not hard-error like an unknown argument) and be captured so the
        // dispatch path can emit the deprecation warning.
        #[derive(Parser)]
        struct TestCli {
            #[command(flatten)]
            scheduler: SchedulerOptions,
        }

        let cli = TestCli::try_parse_from(["test", "--scheduler", "work-stealing"])
            .expect("--scheduler must be accepted as a deprecated no-op, not rejected");
        assert_eq!(cli.scheduler.scheduler.as_deref(), Some("work-stealing"));

        // Absent by default.
        let cli = TestCli::try_parse_from(["test"]).expect("parses with no args");
        assert_eq!(cli.scheduler.scheduler, None);
    }

    // ========== Tests for QueueMemoryOptions ==========

    #[test]
    fn test_queue_memory_default_is_768_mib_per_thread() {
        // The multi-threaded default must stay byte-identical to the historical
        // `--queue-memory 768` (768 MiB per thread).
        let opts = QueueMemoryOptions::default();
        let result = opts
            .calculate_memory_limit(4)
            .expect("calculate_memory_limit should succeed for the default");
        assert_eq!(result, 768 * 1024 * 1024 * 4);
    }

    #[test]
    fn test_queue_memory_max_memory_per_thread_scaling() {
        // 100 MiB × 4 threads = 400 MiB
        let opts = QueueMemoryOptions {
            max_memory: Some(MemoryLimit::Fixed(100 * 1024 * 1024)),
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(4).expect("should succeed");
        assert_eq!(result, 100 * 1024 * 1024 * 4);
    }

    #[test]
    fn test_queue_memory_fixed_total_does_not_scale() {
        let opts = QueueMemoryOptions {
            max_memory: Some(MemoryLimit::Fixed(200 * 1024 * 1024)),
            memory_per_thread: false,
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(8).expect("should succeed");
        assert_eq!(result, 200 * 1024 * 1024);
    }

    #[test]
    fn test_queue_memory_zero_threads_is_error() {
        // Even with an explicit limit, 0 threads must error (not the lean
        // single-thread path, which only triggers for the unset default).
        let opts = QueueMemoryOptions {
            max_memory: Some(MemoryLimit::Fixed(100 * 1024 * 1024)),
            ..QueueMemoryOptions::default()
        };
        assert!(opts.calculate_memory_limit(0).is_err());
    }

    #[test]
    fn test_queue_memory_human_readable_via_parse() {
        // "2GB" parses as decimal gigabytes (matching `fgumi sort` / ByteSize).
        let opts = QueueMemoryOptions {
            max_memory: Some(parse_memory("2GB").expect("parse 2GB")),
            memory_per_thread: false,
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(4).expect("should succeed");
        assert_eq!(result, 2 * 1000 * 1000 * 1000);
    }

    #[test]
    fn test_queue_memory_auto_is_bounded_by_host() {
        // `auto` self-throttles: the budget must never exceed (cgroup-aware)
        // host memory, regardless of thread count.
        let total = detect_total_memory() as u64;
        let opts = QueueMemoryOptions {
            max_memory: Some(MemoryLimit::Auto),
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(4).expect("auto should resolve");
        assert!(result > 0, "auto resolved to a zero budget");
        assert!(result <= total, "auto budget {result} exceeded host total {total}");
    }

    // NOTE: `test_auto_never_oversubscribes_small_host`,
    // `test_auto_uses_floor_when_host_is_ample`, and
    // `test_fixed_budget_independent_of_host` exercised the (now-private)
    // `resolve_memory_budget_with_total` helper, which moved to
    // `fgumi-cli-common` along with the resolver itself. Their coverage lives
    // in that crate's unit tests now.

    #[test]
    fn test_queue_memory_auto_reserve_shrinks_budget() {
        let large_reserve = QueueMemoryOptions {
            max_memory: Some(MemoryLimit::Auto),
            memory_reserve: MemoryReserve::Fixed(512 * 1024 * 1024),
            ..QueueMemoryOptions::default()
        }
        .calculate_memory_limit(4)
        .expect("should succeed");
        let small_reserve = QueueMemoryOptions {
            max_memory: Some(MemoryLimit::Auto),
            memory_reserve: MemoryReserve::Fixed(128 * 1024 * 1024),
            ..QueueMemoryOptions::default()
        }
        .calculate_memory_limit(4)
        .expect("should succeed");
        assert!(large_reserve <= small_reserve);
    }

    #[test]
    fn test_queue_memory_deprecated_queue_memory_alias() {
        // The hidden `--queue-memory` alias overrides `--max-memory`, honoring
        // its companion `--queue-memory-per-thread`.
        let opts = QueueMemoryOptions {
            queue_memory: Some("200".to_string()),
            queue_memory_per_thread: Some(false),
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(8).expect("should succeed");
        assert_eq!(result, 200 * 1024 * 1024); // fixed total, no scaling
    }

    #[test]
    fn test_queue_memory_per_thread_alias_alone_takes_effect() {
        // `--queue-memory-per-thread false` on its own must override
        // `--memory-per-thread` (default true), not be a silent no-op.
        let opts = QueueMemoryOptions {
            max_memory: Some(MemoryLimit::Fixed(200 * 1024 * 1024)),
            queue_memory_per_thread: Some(false),
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(8).expect("should succeed");
        assert_eq!(result, 200 * 1024 * 1024); // fixed total, no per-thread scaling
    }

    #[test]
    fn test_queue_memory_deprecated_alias_defaults_to_per_thread() {
        // When `--queue-memory-per-thread` is unset, the alias falls back to
        // `--memory-per-thread` (default true).
        let opts = QueueMemoryOptions {
            queue_memory: Some("100".to_string()),
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(4).expect("should succeed");
        assert_eq!(result, 100 * 1024 * 1024 * 4);
    }

    #[test]
    fn test_queue_memory_legacy_limit_mb_fixed_total() {
        let opts = QueueMemoryOptions {
            queue_memory_limit_mb: Some(2048),
            // These must be ignored in favor of the legacy fixed total.
            max_memory: Some(MemoryLimit::Fixed(1)),
            memory_per_thread: true,
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(4).expect("should succeed");
        assert_eq!(result, 2048 * 1024 * 1024); // fixed total, no scaling
    }

    #[test]
    fn test_legacy_limit_mb_overflow_is_error() {
        // u64::MAX MiB overflows when multiplied by 1 MiB.
        let opts = QueueMemoryOptions {
            queue_memory_limit_mb: Some(u64::MAX),
            ..QueueMemoryOptions::default()
        };
        assert!(opts.calculate_memory_limit(1).is_err());
    }

    #[test]
    fn test_deprecated_queue_memory_parse_error_is_error() {
        let opts = QueueMemoryOptions {
            queue_memory: Some("not-a-size".to_string()),
            ..QueueMemoryOptions::default()
        };
        assert!(opts.calculate_memory_limit(4).is_err());
    }

    #[test]
    fn test_fixed_per_thread_overflow_is_error() {
        let opts = QueueMemoryOptions {
            max_memory: Some(MemoryLimit::Fixed(usize::MAX)),
            ..QueueMemoryOptions::default()
        };
        assert!(opts.calculate_memory_limit(4).is_err());
    }

    // `test_auto_per_thread_overflow_is_error` moved to `fgumi-cli-common`
    // with the `resolve_memory_budget_with_total` helper it exercised.

    #[test]
    fn test_log_memory_config_exercises_both_branches() {
        // Logging only (no return value); call both paths to keep them covered
        // and to guard against a panic in the formatting/division.
        let per_thread = QueueMemoryOptions::default();
        per_thread.log_memory_config(8, 8 * 768 * 1024 * 1024); // per-thread, threads > 1
        per_thread.log_memory_config(1, 768 * 1024 * 1024); // threads == 1 → total branch

        let fixed_total =
            QueueMemoryOptions { memory_per_thread: false, ..QueueMemoryOptions::default() };
        fixed_total.log_memory_config(8, 1024 * 1024 * 1024); // fixed total branch
    }

    #[test]
    fn unset_max_memory_lean_at_single_thread_default_at_multi() {
        // `--max-memory` unset (`None`): lean 64 MiB TOTAL at threads==1 (a lone
        // worker streams), but the historical 768 MiB/thread at t>1. This is the
        // #330 lean single-thread default, preserved through the #381 host-aware
        // resolution.
        let opts = QueueMemoryOptions::default(); // max_memory: None
        assert_eq!(
            opts.calculate_memory_limit(1).expect("t=1 default"),
            64 * 1024 * 1024,
            "single-threaded default is the lean 64 MiB total"
        );
        assert_eq!(
            opts.calculate_memory_limit(4).expect("t=4 default"),
            768 * 1024 * 1024 * 4,
            "multi-threaded default is unchanged at 768 MiB per thread"
        );
    }

    #[test]
    fn explicit_max_memory_not_clobbered_at_single_thread() {
        // The override-detection guard: an explicit `--max-memory` at threads==1
        // must be honored, NOT silently rewritten to the 64 MiB lean default.
        // (`Option<MemoryLimit>` makes "user typed 768" — `Some` — distinct from
        // the unset default — `None`.)
        let explicit = QueueMemoryOptions {
            max_memory: Some(MemoryLimit::Fixed(768 * 1024 * 1024)),
            ..QueueMemoryOptions::default()
        };
        assert_eq!(
            explicit.calculate_memory_limit(1).expect("explicit t=1"),
            768 * 1024 * 1024,
            "explicit --max-memory at t=1 is honored (per-thread × 1), not the lean default"
        );
        let explicit_total = QueueMemoryOptions {
            max_memory: Some(MemoryLimit::Fixed(512 * 1024 * 1024)),
            memory_per_thread: false,
            ..QueueMemoryOptions::default()
        };
        assert_eq!(
            explicit_total.calculate_memory_limit(1).expect("explicit total t=1"),
            512 * 1024 * 1024,
            "explicit total at t=1 is honored, not the lean default"
        );
    }

    #[test]
    #[cfg(feature = "memory-debug")]
    fn test_sysinfo_returns_reasonable_values() {
        use sysinfo::System;
        let mut system = System::new();
        system.refresh_memory();

        let total = system.total_memory();
        let available = system.available_memory();

        assert!(total > 100_000_000); // > 100MB
        assert!(available > 0);
        assert!(available <= total);
    }

    use clap::Parser;

    /// Test-only wrapper to exercise clap parsing of flattened Args structs.
    #[derive(Debug, Parser)]
    #[command(name = "test")]
    struct TestBoolFlags {
        #[command(flatten)]
        consensus: ConsensusCallingOptions,
        #[command(flatten)]
        overlapping: OverlappingConsensusOptions,
        #[command(flatten)]
        queue_memory: QueueMemoryOptions,
    }

    use rstest::rstest;

    #[rstest]
    // --output-per-base-tags (default true)
    #[case(&["test"], true)]
    #[case(&["test", "--output-per-base-tags"], true)]
    #[case(&["test", "--output-per-base-tags", "true"], true)]
    #[case(&["test", "--output-per-base-tags", "false"], false)]
    #[case(&["test", "--output-per-base-tags=true"], true)]
    #[case(&["test", "--output-per-base-tags=false"], false)]
    fn test_output_per_base_tags_parsing(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.consensus.output_per_base_tags, expected);
    }

    #[rstest]
    // --trim (default false)
    #[case(&["test"], false)]
    #[case(&["test", "--trim"], true)]
    #[case(&["test", "--trim", "true"], true)]
    #[case(&["test", "--trim", "false"], false)]
    #[case(&["test", "--trim=true"], true)]
    #[case(&["test", "--trim=false"], false)]
    fn test_trim_parsing(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.consensus.trim, expected);
    }

    #[rstest]
    // --consensus-call-overlapping-bases (default true)
    #[case(&["test"], true)]
    #[case(&["test", "--consensus-call-overlapping-bases"], true)]
    #[case(&["test", "--consensus-call-overlapping-bases", "true"], true)]
    #[case(&["test", "--consensus-call-overlapping-bases", "false"], false)]
    #[case(&["test", "--consensus-call-overlapping-bases=true"], true)]
    #[case(&["test", "--consensus-call-overlapping-bases=false"], false)]
    fn test_overlapping_bases_parsing(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.overlapping.consensus_call_overlapping_bases, expected);
    }

    #[rstest]
    // --memory-per-thread (default true)
    #[case(&["test"], true)]
    #[case(&["test", "--memory-per-thread"], true)]
    #[case(&["test", "--memory-per-thread", "true"], true)]
    #[case(&["test", "--memory-per-thread", "false"], false)]
    #[case(&["test", "--memory-per-thread=true"], true)]
    #[case(&["test", "--memory-per-thread=false"], false)]
    fn test_memory_per_thread_parsing(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.queue_memory.memory_per_thread, expected);
    }

    #[rstest]
    // --max-memory parses to a MemoryLimit; "auto" and explicit values both work.
    // Unset stays `None` (the thread-count-aware default applies at resolve time).
    #[case(&["test"], None)]
    #[case(&["test", "--max-memory", "auto"], Some(MemoryLimit::Auto))]
    #[case(&["test", "--max-memory", "2GiB"], Some(MemoryLimit::Fixed(2 * 1024 * 1024 * 1024)))]
    #[case(&["test", "--max-memory=512M"], Some(MemoryLimit::Fixed(512 * 1000 * 1000)))]
    fn test_max_memory_parsing(#[case] args: &[&str], #[case] expected: Option<MemoryLimit>) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.queue_memory.max_memory, expected);
    }

    #[test]
    fn test_deprecated_queue_memory_flag_still_parses() {
        // The hidden aliases remain accepted for backward compatibility.
        let cmd = TestBoolFlags::try_parse_from([
            "test",
            "--queue-memory",
            "512",
            "--queue-memory-per-thread",
            "false",
        ])
        .expect("deprecated flags should still parse");
        assert_eq!(cmd.queue_memory.queue_memory.as_deref(), Some("512"));
        assert_eq!(cmd.queue_memory.queue_memory_per_thread, Some(false));
    }

    #[rstest]
    #[case("true", true)]
    #[case("false", false)]
    #[case("yes", true)]
    #[case("no", false)]
    #[case("t", true)]
    #[case("f", false)]
    #[case("y", true)]
    #[case("n", false)]
    #[case("True", true)]
    #[case("TRUE", true)]
    #[case("False", false)]
    #[case("FALSE", false)]
    #[case("Yes", true)]
    #[case("YES", true)]
    #[case("No", false)]
    #[case("NO", false)]
    #[case("T", true)]
    #[case("F", false)]
    #[case("Y", true)]
    #[case("N", false)]
    #[case("tRuE", true)]
    #[case("fAlSe", false)]
    #[case("yEs", true)]
    fn test_parse_bool_valid(#[case] input: &str, #[case] expected: bool) {
        assert_eq!(parse_bool(input).expect("should parse"), expected);
    }

    #[rstest]
    #[case("")]
    #[case("tru")]
    #[case("fals")]
    #[case("truee")]
    #[case("noo")]
    #[case("yess")]
    #[case("maybe")]
    #[case("0")]
    #[case("1")]
    #[case("on")]
    #[case("off")]
    #[case(" true")]
    #[case("true ")]
    fn test_parse_bool_invalid(#[case] input: &str) {
        assert!(parse_bool(input).is_err(), "expected error for input: {input:?}");
    }

    #[rstest]
    #[case(&["test", "--trim", "yes"], true)]
    #[case(&["test", "--trim", "no"], false)]
    #[case(&["test", "--trim", "y"], true)]
    #[case(&["test", "--trim", "n"], false)]
    #[case(&["test", "--trim", "t"], true)]
    #[case(&["test", "--trim", "f"], false)]
    #[case(&["test", "--trim", "YES"], true)]
    #[case(&["test", "--trim", "NO"], false)]
    #[case(&["test", "--trim=yes"], true)]
    #[case(&["test", "--trim=no"], false)]
    fn test_extended_bool_values_in_cli(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.consensus.trim, expected);
    }

    #[rstest]
    #[case(&["test", "--trim", "maybe"])]
    #[case(&["test", "--trim", "0"])]
    #[case(&["test", "--trim", "1"])]
    #[case(&["test", "--trim", "on"])]
    #[case(&["test", "--trim", "off"])]
    fn test_extended_bool_values_in_cli_invalid(#[case] args: &[&str]) {
        assert!(TestBoolFlags::try_parse_from(args).is_err());
    }

    // -------------------------------------------------------------------------
    // cgroup-aware memory detection
    // -------------------------------------------------------------------------

    #[test]
    fn test_detect_total_memory_returns_nonzero() {
        let total = detect_total_memory();
        assert!(total > 0, "detect_total_memory returned 0");
        // Sanity: at least 64 MiB (even the smallest CI runner has more than this).
        assert!(total >= 64 * 1024 * 1024, "detect_total_memory returned < 64 MiB: {total}");
    }

    #[test]
    fn test_detect_total_memory_bounded_by_sysinfo() {
        // cgroup_limits().total_memory is min(cgroup_max, physical_ram), so
        // detect_total_memory() can never exceed what sysinfo reports.
        let total = detect_total_memory();
        let mut system = sysinfo::System::new();
        system.refresh_memory();
        let sysinfo_total = usize::try_from(system.total_memory()).unwrap_or(usize::MAX);
        assert!(
            total <= sysinfo_total,
            "cgroup-limited total {total} exceeded sysinfo total {sysinfo_total}"
        );
    }
}
