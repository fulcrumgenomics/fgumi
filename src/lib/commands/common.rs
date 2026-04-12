//! Common CLI options shared across commands.
//!
//! This module provides shared argument structures that can be composed into
//! command structs using `#[command(flatten)]`.

use std::path::PathBuf;
use std::sync::Arc;

use crate::bam_io::is_stdin_path;
use crate::logging::OperationTimer;
use crate::unified_pipeline::{BamPipelineConfig, SchedulerStrategy};
use crate::validation::validate_file_exists;
use bytesize::ByteSize;
use clap::Args;
use fgumi_consensus::methylation::RefBaseProvider;
use log::info;
use noodles::sam::Header;

/// EM-Seq reference pair: reference base provider + contig name mapping.
pub type EmSeqRef = Option<(
    Arc<dyn fgumi_consensus::methylation::RefBaseProvider + Send + Sync>,
    Arc<Vec<String>>,
)>;

/// Loads the reference FASTA and builds contig name mapping for EM-Seq mode.
///
/// Returns `None` if `em_seq` is false. Errors if `em_seq` is true but `reference` is `None`.
pub fn load_em_seq_reference(
    em_seq: bool,
    reference: &Option<PathBuf>,
    header: &Header,
) -> anyhow::Result<EmSeqRef> {
    if !em_seq {
        return Ok(None);
    }
    let ref_path = reference
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--ref is required when --em-seq is enabled"))?;
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

    info!("EM-Seq mode enabled with {} reference contigs", ref_names.len());
    Ok(Some((reference, Arc::new(ref_names))))
}

/// Add a @PG record to an existing header, using the current fgumi version.
///
/// Wraps [`crate::header::add_pg_record`] with the binary's version string.
pub fn add_pg_record(header: Header, command_line: &str) -> anyhow::Result<Header> {
    crate::header::add_pg_record(header, crate::version::VERSION.as_str(), command_line)
}

/// Add a @PG record to a header builder, using the current fgumi version.
///
/// Wraps [`crate::header::add_pg_to_builder`] with the binary's version string.
pub fn add_pg_to_builder(
    builder: noodles::sam::header::Builder,
    command_line: &str,
) -> anyhow::Result<noodles::sam::header::Builder> {
    crate::header::add_pg_to_builder(builder, crate::version::VERSION.as_str(), command_line)
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
}

impl BamIoOptions {
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
    /// Compression level for output BAM (1-12).
    ///
    /// Level 1 is fastest with larger files.
    /// Level 12 produces smallest files but is slowest.
    #[arg(long, default_value_t = 1)]
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

/// Options for pipeline queue memory limits.
///
/// Controls memory usage in pipeline queues to prevent out-of-memory conditions.
/// Supports human-readable formats for better UX.
#[derive(Debug, Clone, Args)]
pub struct QueueMemoryOptions {
    /// Pipeline queue memory limit per thread (default) or total.
    ///
    /// Plain numbers are interpreted as MB. Also supports human-readable
    /// formats like "2GB", "1.5GB", or "1024MiB".
    /// By default this value is per-thread, so with --threads 8 the total
    /// memory will be 8x this value. Use --queue-memory-per-thread false
    /// for a fixed total limit.
    #[arg(long = "queue-memory", default_value = "768")]
    pub queue_memory: String,

    /// Interpret --queue-memory as per-thread (true, default) or total (false).
    ///
    /// When true, total memory = queue-memory * threads. For example,
    /// --queue-memory 768 with --threads 16 allocates 12 GB total.
    /// Set to false for a fixed total memory budget regardless of thread count.
    #[arg(long = "queue-memory-per-thread", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub queue_memory_per_thread: bool,

    /// DEPRECATED: Use --queue-memory instead. Memory limit for pipeline queues in megabytes.
    #[arg(long = "queue-memory-limit-mb", hide = true)]
    pub queue_memory_limit_mb: Option<u64>,
}

impl Default for QueueMemoryOptions {
    fn default() -> Self {
        Self {
            queue_memory: "768".to_string(),
            queue_memory_per_thread: true,
            queue_memory_limit_mb: None,
        }
    }
}

impl QueueMemoryOptions {
    /// Maximum reasonable memory per thread (1TB)
    const MAX_MEMORY_PER_THREAD: u64 = 1024 * 1024 * 1024 * 1024;

    /// Minimum reasonable memory (1MB)
    const MIN_MEMORY_TOTAL: u64 = 1024 * 1024;

    /// Calculates the total queue memory limit in bytes with comprehensive validation.
    ///
    /// This includes system memory checks via `sysinfo` (warns if >90% of system memory).
    /// Use [`compute_memory_limit`](Self::compute_memory_limit) for the pure arithmetic
    /// without system queries (e.g. in tests).
    ///
    /// # Errors
    /// Returns an error if the memory size string cannot be parsed, `num_threads` is 0,
    /// or the memory calculation would overflow.
    pub fn calculate_memory_limit(&self, num_threads: usize) -> anyhow::Result<u64> {
        let total_memory = self.compute_memory_limit(num_threads)?;
        Self::validate_against_system_memory(total_memory);
        Ok(total_memory)
    }

    /// Computes the total queue memory limit in bytes (pure arithmetic, no system queries).
    ///
    /// # Errors
    /// Returns an error if:
    /// - The memory size string cannot be parsed
    /// - `num_threads` is 0
    /// - Memory calculation would overflow
    /// - Memory values are unreasonable
    pub fn compute_memory_limit(&self, num_threads: usize) -> anyhow::Result<u64> {
        // Validate thread count
        if num_threads == 0 {
            anyhow::bail!("Number of threads must be greater than 0, got: {num_threads}");
        }

        // Handle migration from old parameter
        let (base_memory_bytes, is_legacy) = if let Some(legacy_mb) = self.queue_memory_limit_mb {
            log::warn!(
                "DEPRECATED: --queue-memory-limit-mb is deprecated. Use --queue-memory instead."
            );
            log::warn!(
                "Migration: --queue-memory-limit-mb {legacy_mb} → --queue-memory {legacy_mb} --queue-memory-per-thread false"
            );
            (
                legacy_mb.checked_mul(1024 * 1024).ok_or_else(|| {
                    anyhow::anyhow!("Legacy memory size overflow: {legacy_mb} MB")
                })?,
                true,
            )
        } else {
            (
                parse_memory_size(&self.queue_memory).map_err(|e| {
                    anyhow::anyhow!("Failed to parse queue memory size: {}: {e}", self.queue_memory)
                })?,
                false,
            )
        };

        // Validate base memory range
        if base_memory_bytes < Self::MIN_MEMORY_TOTAL {
            anyhow::bail!(
                "Memory limit too small: {} (minimum: {})",
                ByteSize(base_memory_bytes),
                ByteSize(Self::MIN_MEMORY_TOTAL)
            );
        }

        // Calculate total memory with overflow checking
        let total_memory = if self.queue_memory_per_thread && !is_legacy {
            // Validate per-thread memory limit
            if base_memory_bytes > Self::MAX_MEMORY_PER_THREAD {
                anyhow::bail!(
                    "Memory per thread too large: {} (maximum: {})",
                    ByteSize(base_memory_bytes),
                    ByteSize(Self::MAX_MEMORY_PER_THREAD)
                );
            }

            // Check for overflow before multiplication
            base_memory_bytes
                .checked_mul(num_threads as u64)
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "Memory calculation overflow: {} × {} threads exceeds maximum addressable memory",
                        ByteSize(base_memory_bytes),
                        num_threads
                    )
                })?
        } else {
            // Fixed total memory (legacy mode or explicit setting)
            base_memory_bytes
        };

        Ok(total_memory)
    }

    /// Warns if the requested memory exceeds reasonable system limits.
    fn validate_against_system_memory(requested_bytes: u64) {
        let mut system = sysinfo::System::new();
        system.refresh_memory();

        // Use cgroup-aware total clamped to physical RAM, matching detect_total_memory().
        // A misconfigured container may report a cgroup limit exceeding physical RAM;
        // clamping keeps warnings consistent with the value used by resolve_memory_limit.
        let physical = system.total_memory();
        let total_memory_bytes =
            system.cgroup_limits().map_or(physical, |c| c.total_memory.min(physical));

        // available_memory reflects free physical pages; no cgroup equivalent.
        let available_memory_bytes = system.available_memory();

        // Calculate 90% limit using integer arithmetic to avoid precision loss
        let memory_limit = total_memory_bytes - (total_memory_bytes / 10); // 90% = total - 10%
        if requested_bytes > memory_limit {
            log::warn!(
                "Requested memory {} exceeds 90% of system memory ({}). System has {} total, {} available. This may cause OOM conditions.",
                ByteSize(requested_bytes),
                ByteSize(memory_limit),
                ByteSize(total_memory_bytes),
                ByteSize(available_memory_bytes)
            );
        }

        // Warn if requesting more than currently available memory
        if requested_bytes > available_memory_bytes {
            log::warn!(
                "Requested memory {} exceeds currently available memory {}. This may cause swapping.",
                ByteSize(requested_bytes),
                ByteSize(available_memory_bytes)
            );
        }
    }

    /// Logs the memory configuration.
    ///
    /// # Arguments
    /// * `num_threads` - Number of threads for the calculation
    /// * `total_memory` - Pre-computed total memory limit in bytes from `calculate_memory_limit`
    pub fn log_memory_config(&self, num_threads: usize, total_memory: u64) {
        if let Some(legacy_mb) = self.queue_memory_limit_mb {
            log::info!(
                "Queue memory limit: {} (LEGACY: {legacy_mb} MB total, per-thread scaling disabled)",
                ByteSize(total_memory)
            );
        } else if self.queue_memory_per_thread && num_threads > 1 {
            log::info!(
                "Queue memory limit: {} total ({} MB/thread × {} threads)",
                ByteSize(total_memory),
                self.queue_memory,
                num_threads
            );
        } else {
            log::info!("Queue memory limit: {} total (fixed)", ByteSize(total_memory));
        }
    }
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
    fn test_queue_memory_options_compute_basic() {
        let opts = QueueMemoryOptions {
            queue_memory: "100".to_string(),
            queue_memory_per_thread: true,
            queue_memory_limit_mb: None,
        };

        // 100MB × 4 threads = 400MB
        let result = opts
            .compute_memory_limit(4)
            .expect("compute_memory_limit should succeed for 100MB x 4 threads");
        assert_eq!(result, 100 * 1024 * 1024 * 4);

        // Fixed memory (no scaling)
        let opts_fixed = QueueMemoryOptions {
            queue_memory: "200".to_string(),
            queue_memory_per_thread: false,
            queue_memory_limit_mb: None,
        };
        let result_fixed = opts_fixed
            .compute_memory_limit(8)
            .expect("compute_memory_limit should succeed for fixed 200MB");
        assert_eq!(result_fixed, 200 * 1024 * 1024); // Should not scale
    }

    #[test]
    fn test_queue_memory_options_validation_errors() {
        let opts = QueueMemoryOptions::default();

        // Zero threads should fail
        assert!(opts.compute_memory_limit(0).is_err());

        // Very large memory per thread should fail
        let large_opts = QueueMemoryOptions {
            queue_memory: "2TB".to_string(),
            queue_memory_per_thread: true,
            queue_memory_limit_mb: None,
        };
        assert!(large_opts.compute_memory_limit(1).is_err());

        // Very small memory should fail
        let tiny_opts = QueueMemoryOptions {
            queue_memory: "1KB".to_string(),
            queue_memory_per_thread: false,
            queue_memory_limit_mb: None,
        };
        assert!(tiny_opts.compute_memory_limit(1).is_err());
    }

    #[test]
    fn test_queue_memory_options_overflow() {
        let opts = QueueMemoryOptions {
            queue_memory: "2TB".to_string(), // 2TB > 1TB limit
            queue_memory_per_thread: true,
            queue_memory_limit_mb: None,
        };

        // Should fail due to per-thread limit
        assert!(opts.compute_memory_limit(1).is_err());

        // Large total (100TB) succeeds at the math level (system check is separate)
        let opts2 = QueueMemoryOptions {
            queue_memory: "100GB".to_string(),
            queue_memory_per_thread: true,
            queue_memory_limit_mb: None,
        };
        assert!(opts2.compute_memory_limit(1000).is_ok());
    }

    #[test]
    fn test_queue_memory_options_legacy_migration() {
        let legacy_opts = QueueMemoryOptions {
            queue_memory: "768".to_string(), // Should be ignored
            queue_memory_per_thread: true,   // Should be ignored
            queue_memory_limit_mb: Some(2048),
        };

        let result = legacy_opts
            .compute_memory_limit(4)
            .expect("compute_memory_limit should succeed for legacy migration");
        assert_eq!(result, 2048 * 1024 * 1024); // Should use legacy value, no scaling
    }

    #[test]
    fn test_queue_memory_options_human_readable() {
        let opts = QueueMemoryOptions {
            queue_memory: "2GB".to_string(),
            queue_memory_per_thread: false,
            queue_memory_limit_mb: None,
        };

        let result = opts
            .compute_memory_limit(4)
            .expect("compute_memory_limit should succeed for 2GB fixed");
        assert_eq!(result, 2 * 1000 * 1000 * 1000); // 2GB in bytes
    }

    #[test]
    fn test_queue_memory_options_small_value() {
        let opts = QueueMemoryOptions {
            queue_memory: "1".to_string(), // 1MB
            queue_memory_per_thread: false,
            queue_memory_limit_mb: None,
        };

        assert!(opts.compute_memory_limit(1).is_ok());
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
    // --queue-memory-per-thread (default true)
    #[case(&["test"], true)]
    #[case(&["test", "--queue-memory-per-thread"], true)]
    #[case(&["test", "--queue-memory-per-thread", "true"], true)]
    #[case(&["test", "--queue-memory-per-thread", "false"], false)]
    #[case(&["test", "--queue-memory-per-thread=true"], true)]
    #[case(&["test", "--queue-memory-per-thread=false"], false)]
    fn test_queue_memory_per_thread_parsing(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.queue_memory.queue_memory_per_thread, expected);
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
