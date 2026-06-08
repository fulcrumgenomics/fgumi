//! Common CLI options shared across commands.
//!
//! This module provides shared argument structures that can be composed into
//! command structs using `#[command(flatten)]`.

use std::path::PathBuf;
#[cfg(feature = "simplex")]
use std::sync::Arc;

#[cfg(feature = "simplex")]
use crate::logging::OperationTimer;
use crate::unified_pipeline::{BamPipelineConfig, SchedulerStrategy};
use crate::validation::validate_file_exists;
use bytesize::ByteSize;
use clap::Args;
use fgumi_bam_io::is_stdin_path;
#[cfg(feature = "simplex")]
use fgumi_consensus::methylation::RefBaseProvider;
#[cfg(feature = "simplex")]
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
#[cfg(feature = "simplex")]
pub type MethylationRef = Option<(
    Arc<dyn fgumi_consensus::methylation::RefBaseProvider + Send + Sync>,
    Arc<Vec<String>>,
)>;

/// Loads the reference FASTA and builds contig name mapping for methylation-aware modes.
///
/// Returns `None` if methylation mode is disabled. Errors if enabled but `reference` is `None`.
#[cfg(feature = "simplex")]
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

    /// Enable async userspace prefetch on the input BAM.
    ///
    /// Spawns a dedicated I/O thread that reads raw bytes into a bounded
    /// queue ahead of the decompression step, so processing threads do
    /// not block on disk. Prototype flag; defaults to off.
    #[arg(long = "async-reader", default_value_t = false, hide = true)]
    pub async_reader: bool,
}

impl Default for BamIoOptions {
    fn default() -> Self {
        Self { input: PathBuf::new(), output: PathBuf::new(), async_reader: false }
    }
}

impl BamIoOptions {
    /// Construct a `BamIoOptions` from input and output paths. Leaves
    /// opt-in tuning flags (e.g. `async_reader`) at their default values.
    pub fn new(input: impl Into<PathBuf>, output: impl Into<PathBuf>) -> Self {
        Self { input: input.into(), output: output.into(), async_reader: false }
    }

    /// Build [`fgumi_bam_io::PipelineReaderOpts`] from the async-reader flag.
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
/// # Uses up to 8 threads with work-stealing scheduler
/// ```
#[derive(Debug, Clone, Args)]
pub struct ThreadingOptions {
    /// Number of threads for the multi-threaded pipeline.
    ///
    /// If not specified, uses a single-threaded fast path optimized for
    /// simple streaming. When specified (even with --threads 1), uses the
    /// 7-step parallel pipeline with work-stealing scheduler.
    #[arg(long = "threads")]
    pub threads: Option<usize>,
}

/// Options for output compression.
///
/// Controls BGZF compression level for BAM output files.
#[derive(Debug, Clone, Default, Args)]
pub struct CompressionOptions {
    /// Compression level for output BAM (0-12).
    ///
    /// Level 0 writes uncompressed BGZF (valid BAM, no DEFLATE) — useful for
    /// piped/intermediate outputs where downstream tools will recompress.
    /// Level 1 is fastest DEFLATE with larger files.
    /// Level 12 produces smallest files but is slowest.
    #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u32).range(0..=12))]
    pub compression_level: u32,
}

/// Options for pipeline scheduler configuration.
///
/// Controls which scheduling strategy is used for thread work assignment
/// in the unified pipeline. Also controls pipeline statistics output.
#[derive(Debug, Clone, Default, Args)]
pub struct SchedulerOptions {
    /// Scheduler strategy for thread work assignment.
    ///
    /// - `chase-bottleneck` (default): Threads dynamically follow work through
    ///   the pipeline, moving downstream when output is blocked and upstream
    ///   when input is empty. Shows ~10% improvement at medium thread counts.
    ///
    /// - `fixed-priority`: Assigns fixed thread roles (reader, writer, workers).
    ///   Thread 0 prioritizes reading, Thread N-1 prioritizes writing.
    #[arg(long = "scheduler", value_enum, default_value_t = SchedulerStrategy::default(), hide = true)]
    pub scheduler: SchedulerStrategy,

    /// Print detailed pipeline statistics at completion.
    ///
    /// Shows per-step timing, throughput, contention metrics, and
    /// per-thread work distribution.
    #[arg(long = "pipeline-stats", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool, hide = true)]
    pub pipeline_stats: bool,

    /// Timeout in seconds for deadlock detection (default: 10, 0 = disabled).
    ///
    /// When no progress is made for this duration, a warning is logged with
    /// diagnostic info (queue depths, memory usage, per-queue timestamps).
    #[arg(long = "deadlock-timeout", default_value_t = 10, hide = true)]
    pub deadlock_timeout: u64,

    /// Enable automatic deadlock recovery (default: false, detection only).
    ///
    /// Uses progressive doubling: 2x -> 4x -> unbind, with restoration
    /// after 30s of sustained progress.
    #[arg(long = "deadlock-recover", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool, hide = true)]
    pub deadlock_recover: bool,
}

impl SchedulerOptions {
    /// Returns the scheduler strategy.
    #[must_use]
    pub fn strategy(&self) -> SchedulerStrategy {
        self.scheduler
    }

    /// Returns true if pipeline stats should be collected and printed.
    #[must_use]
    pub fn collect_stats(&self) -> bool {
        self.pipeline_stats
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

/// A memory limit, either auto-detected from the host or a fixed byte count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryLimit {
    /// Detect the (cgroup-aware) host memory and subtract the reserve.
    Auto,
    /// Use a fixed memory limit in bytes.
    Fixed(usize),
}

/// How much memory to reserve for other processes (OS, aligners, etc.) when a
/// memory limit is set to [`MemoryLimit::Auto`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryReserve {
    /// Automatic: `min(10 GiB, 50% of host memory)`.
    Auto,
    /// Reserve a fixed number of bytes.
    Fixed(usize),
}

/// The minimum per-thread memory budget (256 MiB).
pub(crate) const MIN_MEMORY_PER_THREAD: usize = 256 * 1024 * 1024;

/// Default auto-reserve cap: 10 GiB.
pub(crate) const AUTO_RESERVE_CAP: usize = 10 * 1024 * 1024 * 1024;

/// Parse a memory size string into `usize` bytes, suitable for use in clap
/// value parsers.
///
/// Delegates to [`parse_memory_size`] for numeric parsing. Plain numbers are
/// interpreted as MiB (e.g. "768" = 768 MiB). Supports human-readable formats
/// like "2GB", "1GiB", "512MiB". See [`parse_memory_size`] for full details.
fn parse_memory_bytes(s: &str, label: &str) -> Result<usize, String> {
    let bytes = parse_memory_size(s).map_err(|e| e.to_string())?;
    usize::try_from(bytes).map_err(|_| format!("{label} too large: {bytes}"))
}

/// Parse a memory-limit string (e.g. "512M", "1G", "768", "auto").
pub(crate) fn parse_memory(s: &str) -> Result<MemoryLimit, String> {
    let s = s.trim();
    if s.eq_ignore_ascii_case("auto") {
        return Ok(MemoryLimit::Auto);
    }
    Ok(MemoryLimit::Fixed(parse_memory_bytes(s, "Memory size")?))
}

/// Parse a memory-reserve string (e.g. "10G", "auto").
pub(crate) fn parse_memory_reserve(s: &str) -> Result<MemoryReserve, String> {
    let s = s.trim();
    if s.eq_ignore_ascii_case("auto") {
        return Ok(MemoryReserve::Auto);
    }
    Ok(MemoryReserve::Fixed(parse_memory_bytes(s, "Memory reserve")?))
}

/// Resolve a [`MemoryReserve`] to a concrete byte count given total host memory.
pub(crate) fn resolve_reserve(reserve: MemoryReserve, total_memory: usize) -> usize {
    match reserve {
        MemoryReserve::Fixed(bytes) => bytes,
        // min(10 GiB, 50% of host memory)
        MemoryReserve::Auto => AUTO_RESERVE_CAP.min(total_memory / 2),
    }
}

/// Resolve a memory budget to a concrete byte count.
///
/// For [`MemoryLimit::Auto`]: detects total host memory (cgroup-aware via
/// [`detect_total_memory`]), subtracts the reserve, and—when `per_thread` is
/// set—targets each thread's share with a 256 MiB floor. The result is then
/// **capped to the available (post-reserve) memory**, so the floor can never
/// push the total past what the host has. The reserve makes the budget shrink
/// to fit the host, which is what lets pipeline commands self-throttle instead
/// of OOM-ing.
///
/// For [`MemoryLimit::Fixed`]: multiplies by `threads` when `per_thread` is set;
/// the reserve and host size are ignored.
///
/// Calls [`detect_total_memory`] exactly once (it invokes `sysinfo`, which is
/// not free).
pub(crate) fn resolve_memory_budget(
    limit: MemoryLimit,
    reserve: MemoryReserve,
    threads: usize,
    per_thread: bool,
) -> anyhow::Result<usize> {
    // Call once — detect_total_memory() invokes sysinfo, which is not free.
    resolve_memory_budget_with_total(limit, reserve, threads, per_thread, detect_total_memory())
}

/// Pure resolver behind [`resolve_memory_budget`], with `total` (host memory)
/// injected so the `Auto` math is unit-testable on simulated small hosts.
fn resolve_memory_budget_with_total(
    limit: MemoryLimit,
    reserve: MemoryReserve,
    threads: usize,
    per_thread: bool,
    total: usize,
) -> anyhow::Result<usize> {
    if threads == 0 {
        anyhow::bail!("--threads must be at least 1");
    }

    let budget = match limit {
        MemoryLimit::Fixed(bytes) => {
            if per_thread {
                bytes
                    .checked_mul(threads)
                    .ok_or_else(|| anyhow::anyhow!("memory limit × {threads} threads overflowed"))?
            } else {
                bytes
            }
        }
        MemoryLimit::Auto => {
            let margin = resolve_reserve(reserve, total);
            let available = total.saturating_sub(margin);
            // The per-thread floor is a *target*, not a guarantee. On a small
            // host (or high thread count) the floor-based budget can exceed what
            // is actually available; cap it to `available` so `auto` truly
            // self-throttles instead of multiplying the floor past physical
            // memory — the exact OOM this feature exists to prevent (#380).
            let target = if per_thread {
                (available / threads)
                    .max(MIN_MEMORY_PER_THREAD)
                    .checked_mul(threads)
                    .ok_or_else(|| anyhow::anyhow!("auto memory budget overflowed"))?
            } else {
                available.max(MIN_MEMORY_PER_THREAD)
            };
            let budget = target.min(available);
            if budget < target {
                log::warn!(
                    "Auto memory: capping budget to host-available {} ({}/thread target × {} threads \
                     exceeds it after reserve {}); throughput may drop but the run stays within memory",
                    ByteSize(budget as u64),
                    ByteSize(MIN_MEMORY_PER_THREAD as u64),
                    threads,
                    ByteSize(margin as u64),
                );
            }
            log::debug!(
                "Auto memory: {} of {} ({}/thread × {} threads, reserve {})",
                ByteSize(budget as u64),
                ByteSize(total as u64),
                ByteSize((budget / threads) as u64),
                threads,
                ByteSize(margin as u64),
            );
            budget
        }
    };

    if budget > total {
        log::warn!(
            "Memory budget {} exceeds total host memory {}; this may cause OOM (or, for sort, earlier spill-to-disk)",
            ByteSize(budget as u64),
            ByteSize(total as u64),
        );
    }

    Ok(budget)
}

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
    /// Default is "768" (768 MiB) per thread. Pass "auto" to detect host memory
    /// (cgroup-aware) and subtract --memory-reserve, so the budget shrinks to
    /// fit the host and the command self-throttles instead of OOM-ing. Explicit
    /// values like "512M", "1G", "4GiB" are per-thread when --memory-per-thread
    /// is enabled (default); plain numbers are MiB. Mirrors `fgumi sort`'s
    /// --max-memory.
    #[arg(long = "max-memory", default_value = "768", value_parser = parse_memory)]
    pub max_memory: MemoryLimit,

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
    /// --memory-per-thread when --queue-memory is set.
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
            // 768 MiB per thread — byte-identical to the historical
            // --queue-memory 768 default.
            max_memory: MemoryLimit::Fixed(768 * 1024 * 1024),
            memory_reserve: MemoryReserve::Auto,
            memory_per_thread: true,
            queue_memory: None,
            queue_memory_per_thread: None,
            queue_memory_limit_mb: None,
        }
    }
}

impl QueueMemoryOptions {
    /// Resolve the effective `(limit, reserve, per_thread)` triple, honoring the
    /// deprecated `--queue-memory*` aliases over the modern `--max-memory*`
    /// flags. Pure: no warnings, no system queries.
    ///
    /// # Errors
    /// Returns an error if a deprecated `--queue-memory` string fails to parse.
    fn resolve_spec(&self) -> anyhow::Result<(MemoryLimit, MemoryReserve, bool)> {
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
            return Ok((MemoryLimit::Fixed(bytes), MemoryReserve::Auto, false));
        }

        // Deprecated --queue-memory overrides --max-memory when set.
        if let Some(queue_memory) = &self.queue_memory {
            let limit = parse_memory(queue_memory).map_err(|e| {
                anyhow::anyhow!("Failed to parse --queue-memory {queue_memory}: {e}")
            })?;
            return Ok((limit, self.memory_reserve, per_thread));
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
    /// Under `--max-memory auto` the budget is detected from (cgroup-aware) host
    /// memory minus `--memory-reserve`, so it shrinks to fit the host. Under an
    /// explicit value it is `max_memory` (× threads when per-thread).
    ///
    /// # Errors
    /// Returns an error if `num_threads` is 0, a deprecated `--queue-memory`
    /// string fails to parse, or the multiplication overflows.
    pub fn calculate_memory_limit(&self, num_threads: usize) -> anyhow::Result<u64> {
        self.warn_deprecations();
        let (limit, reserve, per_thread) = self.resolve_spec()?;
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

/// Frame raw BAM record bytes into the pipeline's serialize buffer.
///
/// Each element of `records` is treated as the body of a BAM record (no
/// `block_size` prefix). This function writes `<u32 LE block_size><body>`
/// for each record, matching BAM's on-disk record framing. The pipeline's
/// downstream BGZF compression stage operates on the framed bytes verbatim.
///
/// Used by command-level `serialize_fn` and `secondary_serialize_fn`
/// closures that produce raw record bytes (e.g. `commands::filter`,
/// the `--rejects` paths in `commands::correct`/`simplex`/`duplex`/`codec`).
///
/// Generic over `R: AsRef<[u8]>` so callers can pass either
/// `&[fgumi_raw_bam::RawRecord]` or `&[Vec<u8>]` without copying.
///
/// # Errors
///
/// Returns [`std::io::ErrorKind::InvalidData`] if:
/// - a record body exceeds `u32::MAX` bytes (BAM records cannot be larger
///   than ~4 GiB by spec), or
/// - the summed framed size of `records` overflows `usize` (only reachable
///   on pathological inputs; surfaces as a hard error instead of wrapping
///   silently in release builds).
///
/// Both error paths are checked in a single up-front validation pass, so
/// `output` is never partially appended on error: either every record is
/// written or none of `output`'s bytes are touched.
pub(crate) fn serialize_raw_bam_records<R: AsRef<[u8]>>(
    records: &[R],
    output: &mut Vec<u8>,
) -> std::io::Result<u64> {
    // Reserve total framed size up front so the per-record extend_from_slice
    // loop doesn't repeatedly grow `output`. Validates each record's body fits
    // a `u32` `block_size` *and* that the summed framed size fits `usize` in
    // this same pass: by the time we start writing, every record is known
    // good, so the write loop cannot fail and leave `output` partially
    // appended. Either no bytes are written on error or all records are.
    let additional = records.iter().try_fold(0usize, |acc, record| {
        let len = record.as_ref().len();
        u32::try_from(len).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("BAM record too large ({len} bytes) for u32 block_size"),
            )
        })?;
        len.checked_add(4).and_then(|frame| acc.checked_add(frame)).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "serialized BAM batch size overflowed usize",
            )
        })
    })?;
    output.reserve(additional);

    for record in records {
        let body = record.as_ref();
        // Pre-validated in the first pass above; `body.len()` is known to fit
        // a `u32`, so this conversion cannot fail.
        let block_size = u32::try_from(body.len()).expect("body length pre-validated to fit u32");
        output.extend_from_slice(&block_size.to_le_bytes());
        output.extend_from_slice(body);
    }
    Ok(records.len() as u64)
}

/// Parses a boolean value from a string, accepting: true/false, yes/no, y/n, t/f
/// (case-insensitive). Matches sopt/fgbio behavior.
pub(crate) fn parse_bool(s: &str) -> Result<bool, String> {
    match s.to_ascii_lowercase().as_str() {
        "true" | "t" | "yes" | "y" => Ok(true),
        "false" | "f" | "no" | "n" => Ok(false),
        _ => Err(format!("Invalid boolean value '{s}'. Expected: true|false|yes|no|y|n|t|f")),
    }
}

// Re-export from the library crate for backward compatibility.
pub(crate) use crate::system::detect_total_memory;
pub use crate::validation::parse_memory_size;

/// Builds a [`BamPipelineConfig`] from the common CLI option structs.
///
/// This consolidates the pipeline configuration boilerplate that is repeated
/// across all multi-threaded commands: auto-tuning, scheduler strategy,
/// stats collection, deadlock settings, and queue memory limits.
///
/// After calling this, commands can further customize the returned config
/// (e.g. setting `group_key_config` for raw-byte mode).
pub fn build_pipeline_config(
    scheduler_opts: &SchedulerOptions,
    compression: &CompressionOptions,
    queue_memory: &QueueMemoryOptions,
    num_threads: usize,
) -> anyhow::Result<BamPipelineConfig> {
    let mut config = BamPipelineConfig::auto_tuned(num_threads, compression.compression_level);
    config.pipeline.scheduler_strategy = scheduler_opts.strategy();
    if scheduler_opts.collect_stats() {
        config.pipeline = config.pipeline.with_stats(true);
    }
    config.pipeline.deadlock_timeout_secs = scheduler_opts.deadlock_timeout_secs();
    config.pipeline.deadlock_recover_enabled = scheduler_opts.deadlock_recover_enabled();

    let queue_memory_limit_bytes = queue_memory.calculate_memory_limit(num_threads)?;
    config.pipeline.queue_memory_limit = queue_memory_limit_bytes;
    queue_memory.log_memory_config(num_threads, queue_memory_limit_bytes);
    Ok(config)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Enable an at-Trace logger so `log::warn!`/`debug!` macros evaluate their
    /// arguments — without an enabled logger the `log` crate skips argument
    /// evaluation, leaving the formatting expressions inside the memory-budget
    /// warn/debug branches unexecuted under test. nextest runs each test in its
    /// own process, so `try_init` is local and idempotent.
    fn enable_logging() {
        let _ =
            env_logger::builder().is_test(true).filter_level(log::LevelFilter::Trace).try_init();
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
        assert_eq!(opts.strategy(), SchedulerStrategy::BalancedChaseDrain);
        assert!(!opts.collect_stats());
    }

    #[test]
    fn test_scheduler_options_strategy() {
        let opts = SchedulerOptions {
            scheduler: SchedulerStrategy::FixedPriority,
            pipeline_stats: false,
            deadlock_timeout: 10,
            deadlock_recover: false,
        };
        assert_eq!(opts.strategy(), SchedulerStrategy::FixedPriority);
    }

    #[test]
    fn test_scheduler_options_collect_stats() {
        let opts = SchedulerOptions {
            scheduler: SchedulerStrategy::default(),
            pipeline_stats: true,
            deadlock_timeout: 10,
            deadlock_recover: false,
        };
        assert!(opts.collect_stats());
    }

    #[test]
    fn test_scheduler_options_deadlock_timeout() {
        let opts = SchedulerOptions {
            scheduler: SchedulerStrategy::default(),
            pipeline_stats: false,
            deadlock_timeout: 30,
            deadlock_recover: false,
        };
        assert_eq!(opts.deadlock_timeout_secs(), 30);
    }

    #[test]
    fn test_scheduler_options_deadlock_recover() {
        let opts = SchedulerOptions {
            scheduler: SchedulerStrategy::default(),
            pipeline_stats: false,
            deadlock_timeout: 10,
            deadlock_recover: true,
        };
        assert!(opts.deadlock_recover_enabled());
    }

    // ========== Tests for QueueMemoryOptions ==========

    #[test]
    fn test_queue_memory_default_is_768_mib_per_thread() {
        // The default must stay byte-identical to the historical
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
            max_memory: MemoryLimit::Fixed(100 * 1024 * 1024),
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(4).expect("should succeed");
        assert_eq!(result, 100 * 1024 * 1024 * 4);
    }

    #[test]
    fn test_queue_memory_fixed_total_does_not_scale() {
        let opts = QueueMemoryOptions {
            max_memory: MemoryLimit::Fixed(200 * 1024 * 1024),
            memory_per_thread: false,
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(8).expect("should succeed");
        assert_eq!(result, 200 * 1024 * 1024);
    }

    #[test]
    fn test_queue_memory_zero_threads_is_error() {
        assert!(QueueMemoryOptions::default().calculate_memory_limit(0).is_err());
    }

    #[test]
    fn test_queue_memory_human_readable_via_parse() {
        // "2GB" parses as decimal gigabytes (matching `fgumi sort` / ByteSize).
        let opts = QueueMemoryOptions {
            max_memory: parse_memory("2GB").expect("parse 2GB"),
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
        let opts =
            QueueMemoryOptions { max_memory: MemoryLimit::Auto, ..QueueMemoryOptions::default() };
        let result = opts.calculate_memory_limit(4).expect("auto should resolve");
        assert!(result > 0, "auto resolved to a zero budget");
        assert!(result <= total, "auto budget {result} exceeded host total {total}");
    }

    #[test]
    fn test_auto_never_oversubscribes_small_host() {
        enable_logging(); // exercise the cap-warning and auto-debug log branches
        // Simulated 4 GiB host, 16 threads: the 256 MiB/thread floor would want
        // 4 GiB before reserve, which cannot fit after the auto reserve. The
        // budget must be capped to `available`, never `floor × threads`.
        let total = 4 * 1024 * 1024 * 1024; // 4 GiB
        let margin = resolve_reserve(MemoryReserve::Auto, total); // min(10 GiB, 2 GiB) = 2 GiB
        let available = total - margin;
        let budget = resolve_memory_budget_with_total(
            MemoryLimit::Auto,
            MemoryReserve::Auto,
            16,
            true,
            total,
        )
        .expect("should resolve");
        assert!(budget <= available, "budget {budget} oversubscribed available {available}");
        assert!(budget <= total, "budget {budget} oversubscribed host {total}");
    }

    #[test]
    fn test_auto_uses_floor_when_host_is_ample() {
        // Simulated 256 GiB host, 4 threads: plenty of room, so the budget is
        // the per-thread share and stays under available.
        let total = 256 * 1024 * 1024 * 1024;
        let margin = resolve_reserve(MemoryReserve::Auto, total); // 10 GiB cap
        let available = total - margin;
        let budget = resolve_memory_budget_with_total(
            MemoryLimit::Auto,
            MemoryReserve::Auto,
            4,
            true,
            total,
        )
        .expect("should resolve");
        assert!(budget >= MIN_MEMORY_PER_THREAD * 4, "budget {budget} fell below the floor");
        assert!(budget <= available, "budget {budget} exceeded available {available}");
    }

    #[test]
    fn test_fixed_budget_independent_of_host() {
        enable_logging(); // exercise the "budget exceeds host total" warn branch
        // Fixed limits ignore host size entirely (reserve is irrelevant).
        let tiny_host = 512 * 1024 * 1024;
        let budget = resolve_memory_budget_with_total(
            MemoryLimit::Fixed(2 * 1024 * 1024 * 1024),
            MemoryReserve::Auto,
            4,
            false,
            tiny_host,
        )
        .expect("should resolve");
        assert_eq!(budget, 2 * 1024 * 1024 * 1024);
    }

    #[test]
    fn test_queue_memory_auto_reserve_shrinks_budget() {
        let large_reserve = QueueMemoryOptions {
            max_memory: MemoryLimit::Auto,
            memory_reserve: MemoryReserve::Fixed(512 * 1024 * 1024),
            ..QueueMemoryOptions::default()
        }
        .calculate_memory_limit(4)
        .expect("should succeed");
        let small_reserve = QueueMemoryOptions {
            max_memory: MemoryLimit::Auto,
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
            max_memory: MemoryLimit::Fixed(200 * 1024 * 1024),
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
            max_memory: MemoryLimit::Fixed(1),
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
            max_memory: MemoryLimit::Fixed(usize::MAX),
            ..QueueMemoryOptions::default()
        };
        assert!(opts.calculate_memory_limit(4).is_err());
    }

    #[test]
    fn test_auto_per_thread_overflow_is_error() {
        // A pathological thread count makes the per-thread floor × threads
        // overflow; this must surface as an error, not wrap.
        let result = resolve_memory_budget_with_total(
            MemoryLimit::Auto,
            MemoryReserve::Auto,
            usize::MAX,
            true,
            1024 * 1024 * 1024,
        );
        assert!(result.is_err());
    }

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
    #[case(&["test"], MemoryLimit::Fixed(768 * 1024 * 1024))]
    #[case(&["test", "--max-memory", "auto"], MemoryLimit::Auto)]
    #[case(&["test", "--max-memory", "2GiB"], MemoryLimit::Fixed(2 * 1024 * 1024 * 1024))]
    #[case(&["test", "--max-memory=512M"], MemoryLimit::Fixed(512 * 1000 * 1000))]
    fn test_max_memory_parsing(#[case] args: &[&str], #[case] expected: MemoryLimit) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.queue_memory.max_memory, expected);
    }

    #[test]
    fn test_deprecated_queue_memory_flag_still_parses() {
        // The hidden alias remains accepted for backward compatibility.
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

    #[test]
    fn test_serialize_raw_bam_records_empty() {
        let records: Vec<Vec<u8>> = Vec::new();
        let mut output = Vec::new();
        let count = serialize_raw_bam_records(&records, &mut output).unwrap();
        assert_eq!(count, 0);
        assert!(output.is_empty());
    }

    #[test]
    fn test_serialize_raw_bam_records_single_frames_correctly() {
        let records = vec![vec![0xDEu8, 0xAD, 0xBE, 0xEF]];
        let mut output = Vec::new();
        let count = serialize_raw_bam_records(&records, &mut output).unwrap();
        assert_eq!(count, 1);
        // <u32 LE block_size = 4><body = DE AD BE EF>
        assert_eq!(output, vec![0x04, 0x00, 0x00, 0x00, 0xDE, 0xAD, 0xBE, 0xEF]);
    }

    #[test]
    fn test_serialize_raw_bam_records_multiple_frames_concatenated() {
        let records = vec![vec![0x11u8, 0x22], vec![0x33u8, 0x44, 0x55]];
        let mut output = Vec::new();
        let count = serialize_raw_bam_records(&records, &mut output).unwrap();
        assert_eq!(count, 2);
        assert_eq!(
            output,
            vec![0x02, 0x00, 0x00, 0x00, 0x11, 0x22, 0x03, 0x00, 0x00, 0x00, 0x33, 0x44, 0x55],
        );
    }

    #[test]
    fn test_serialize_raw_bam_records_reserves_capacity_upfront() {
        // Build a batch large enough that a naive per-record extend would
        // trigger at least one realloc growth of `output`. We then assert
        // that after one call `output.capacity()` is at least the exact
        // serialized size — proof that the upfront reserve happened.
        let records: Vec<Vec<u8>> = (0..32).map(|i| vec![i as u8; 64]).collect();
        let expected_size: usize = records.iter().map(|r| 4 + r.len()).sum();

        let mut output = Vec::new();
        serialize_raw_bam_records(&records, &mut output).unwrap();
        assert_eq!(output.len(), expected_size);
        assert!(
            output.capacity() >= expected_size,
            "capacity {} should be >= expected serialized size {expected_size}",
            output.capacity(),
        );
    }

    #[test]
    fn test_serialize_raw_bam_records_preserves_existing_output_content() {
        // Reserve must extend, not clear: writing into a non-empty buffer
        // (as the pipeline can do when reusing serialization scratch
        // buffers) must leave the existing prefix untouched.
        let mut output = vec![0xAAu8, 0xBB];
        let records = vec![vec![0x01u8]];
        serialize_raw_bam_records(&records, &mut output).unwrap();
        assert_eq!(output, vec![0xAA, 0xBB, 0x01, 0x00, 0x00, 0x00, 0x01]);
    }
}
