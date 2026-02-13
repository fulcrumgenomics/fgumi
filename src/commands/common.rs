//! Common CLI options shared across commands.
//!
//! This module provides shared argument structures that can be composed into
//! command structs using `#[command(flatten)]`.

use std::path::PathBuf;

use anyhow::Context;
use bytesize::ByteSize;
use clap::Args;
use fgumi_lib::bam_io::is_stdin_path;
use fgumi_lib::unified_pipeline::SchedulerStrategy;
use fgumi_lib::validation::validate_file_exists;
use noodles::sam::Header;

/// Add a @PG record to an existing header, using the current fgumi version.
///
/// Wraps [`fgumi_lib::header::add_pg_record`] with the binary's version string.
pub fn add_pg_record(header: Header, command_line: &str) -> anyhow::Result<Header> {
    fgumi_lib::header::add_pg_record(header, crate::version::VERSION.as_str(), command_line)
}

/// Add a @PG record to a header builder, using the current fgumi version.
///
/// Wraps [`fgumi_lib::header::add_pg_to_builder`] with the binary's version string.
pub fn add_pg_to_builder(
    builder: noodles::sam::header::Builder,
    command_line: &str,
) -> anyhow::Result<noodles::sam::header::Builder> {
    fgumi_lib::header::add_pg_to_builder(builder, crate::version::VERSION.as_str(), command_line)
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
    #[arg(short = 'B', long = "output-per-base-tags", default_value = "true")]
    pub output_per_base_tags: bool,

    /// Quality-trim reads before consensus calling (removes low-quality bases from ends)
    #[arg(long = "trim", default_value = "false")]
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

impl Default for ReadGroupOptions {
    fn default() -> Self {
        Self { read_name_prefix: None, read_group_id: "A".to_string() }
    }
}

/// Options for overlapping bases consensus calling.
#[derive(Debug, Clone, Args)]
pub struct OverlappingConsensusOptions {
    /// Consensus call overlapping bases in read pairs before UMI consensus calling
    #[arg(long = "consensus-call-overlapping-bases", default_value = "true")]
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
    #[arg(long = "pipeline-stats", default_value_t = false, hide = true)]
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
    #[arg(long = "deadlock-recover", default_value_t = false, hide = true)]
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
    #[arg(long = "queue-memory-per-thread", default_value = "true")]
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
        #[cfg(feature = "memory-debug")]
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
                parse_memory_size(&self.queue_memory).with_context(|| {
                    format!("Failed to parse queue memory size: {}", self.queue_memory)
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
    #[cfg(feature = "memory-debug")]
    fn validate_against_system_memory(requested_bytes: u64) {
        use sysinfo::System;
        let mut system = System::new();
        system.refresh_memory();

        let total_memory_bytes = system.total_memory();
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

/// Parses a memory size string into bytes.
///
/// Accepts both plain numbers (interpreted as MB) and human-readable formats like:
/// - "2GB", "2G" -> 2 gigabytes
/// - "1.5GB" -> 1.5 gigabytes
/// - "1024MB", "1024M" -> 1024 megabytes
/// - "512MiB" -> 512 mebibytes
/// - "768" -> 768 megabytes (backward compatibility)
///
/// # Examples
///
/// ```
/// # use anyhow::Result;
/// # use fgumi::commands::common::parse_memory_size;
/// assert_eq!(parse_memory_size("768")?, 768 * 1024 * 1024);
/// assert_eq!(parse_memory_size("2GB")?, 2 * 1000 * 1000 * 1000);
/// assert_eq!(parse_memory_size("1024MiB")?, 1024 * 1024 * 1024);
/// ```
///
/// # Errors
///
/// Returns an error if the string cannot be parsed as a valid size.
pub fn parse_memory_size(size_str: &str) -> anyhow::Result<u64> {
    // Validate input string
    let trimmed = size_str.trim();
    if trimmed.is_empty() {
        anyhow::bail!("Memory size cannot be empty");
    }

    // Handle negative values early
    if trimmed.starts_with('-') {
        anyhow::bail!("Memory size cannot be negative: '{trimmed}'");
    }

    // First try parsing as a plain integer in MB (backward compatibility)
    // Only accept simple integers, not floats or scientific notation
    if let Ok(mb_value) = trimmed.parse::<u64>() {
        // Validate reasonable range for plain numbers
        if mb_value == 0 {
            anyhow::bail!("Memory size cannot be zero");
        }
        if mb_value > 1_000_000 {
            // Sanity guard: >1TB as a plain number likely means the user forgot a unit suffix.
            // Values above this should use human-readable format (e.g. "2TB") which bypasses
            // this check and goes through ByteSize parsing instead.
            anyhow::bail!(
                "Plain number memory size too large: {} MB. Use human-readable format like '{}GB' instead.",
                mb_value,
                mb_value / 1000
            );
        }

        return mb_value
            .checked_mul(1024 * 1024)
            .ok_or_else(|| anyhow::anyhow!("Memory size calculation overflow for {mb_value} MB"));
    }

    // Reject scientific notation (e.g. "1e3") but allow decimals in human-readable sizes (e.g. "1.5GB")
    if trimmed.contains('e') || trimmed.contains('E') {
        anyhow::bail!(
            "Scientific notation not supported: '{trimmed}'. Use integer values or human-readable formats like '2GB'."
        );
    }

    // Reject bare decimal numbers without a unit suffix (e.g. "1.5") since plain numbers are MB
    if trimmed.contains('.') && trimmed.chars().all(|c| c.is_ascii_digit() || c == '.') {
        anyhow::bail!(
            "Plain decimal numbers not supported: '{trimmed}'. Use an integer for MB (e.g. '768') or a human-readable format (e.g. '1.5GB')."
        );
    }

    // Fall back to parsing as a human-readable size (like "2GB", "1024MiB")
    match trimmed.parse::<ByteSize>() {
        Ok(size) => {
            if size.0 == 0 {
                anyhow::bail!("Memory size cannot be zero: '{trimmed}'");
            }
            Ok(size.0)
        }
        Err(_) => {
            anyhow::bail!(
                "Invalid memory size '{trimmed}'. Valid formats:\n\
                 - Plain numbers (interpreted as MB): '768', '4096'\n\
                 - Human-readable (decimal): '2GB', '1024MB'\n\
                 - Human-readable (binary): '1GiB', '512MiB'"
            )
        }
    }
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

    // ========== Tests for memory size parsing ==========

    #[test]
    fn test_parse_memory_size_plain_numbers() {
        assert_eq!(parse_memory_size("768").unwrap(), 768 * 1024 * 1024);
        assert_eq!(parse_memory_size("1").unwrap(), 1024 * 1024);
        assert_eq!(parse_memory_size("4096").unwrap(), 4096 * 1024 * 1024);
    }

    #[test]
    fn test_parse_memory_size_human_readable() {
        assert_eq!(parse_memory_size("2GB").unwrap(), 2 * 1000 * 1000 * 1000);
        assert_eq!(parse_memory_size("2G").unwrap(), 2 * 1000 * 1000 * 1000);
        assert_eq!(parse_memory_size("1024MB").unwrap(), 1024 * 1000 * 1000);
        assert_eq!(parse_memory_size("1024M").unwrap(), 1024 * 1000 * 1000);
        assert_eq!(parse_memory_size("1GiB").unwrap(), 1024 * 1024 * 1024);
        assert_eq!(parse_memory_size("512MiB").unwrap(), 512 * 1024 * 1024);
    }

    #[test]
    fn test_parse_memory_size_invalid() {
        assert!(parse_memory_size("invalid").is_err());
        assert!(parse_memory_size("").is_err());
        assert!(parse_memory_size("GB2").is_err());
    }

    #[test]
    fn test_parse_memory_size_zero() {
        // Zero values should be rejected
        assert!(parse_memory_size("0").is_err());
        assert!(parse_memory_size("0MB").is_err());
        assert!(parse_memory_size("0GB").is_err());
    }

    // ========== Tests for new edge cases and validation ==========

    #[test]
    fn test_parse_memory_size_edge_cases() {
        // Empty and whitespace
        assert!(parse_memory_size("").is_err());
        assert!(parse_memory_size("   ").is_err());

        // Negative values
        assert!(parse_memory_size("-100").is_err());
        assert!(parse_memory_size("-1GB").is_err());

        // Bare floating point without unit suffix (should be rejected)
        assert!(parse_memory_size("1.5").is_err());

        // Floating point with unit suffix (should be accepted via ByteSize)
        assert!(parse_memory_size("1.5GB").is_ok());
        assert!(parse_memory_size("2.5GB").is_ok());

        // Scientific notation (should be rejected for plain numbers)
        assert!(parse_memory_size("1e3").is_err());
        assert!(parse_memory_size("1E6").is_err());

        // Very large plain numbers (should be rejected)
        assert!(parse_memory_size("9999999").is_err());
    }

    #[test]
    fn test_parse_memory_size_whitespace_handling() {
        // Trimmed input should work
        assert_eq!(parse_memory_size("  768  ").unwrap(), 768 * 1024 * 1024);
        assert_eq!(parse_memory_size("\t1GB\n").unwrap(), 1000 * 1000 * 1000);
    }

    #[test]
    fn test_parse_memory_size_overflow() {
        // Very large MB values should be caught
        let very_large = format!("{}", u64::MAX / 1024); // Would overflow when * 1024 * 1024
        assert!(parse_memory_size(&very_large).is_err());
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
        let result = opts.compute_memory_limit(4).unwrap();
        assert_eq!(result, 100 * 1024 * 1024 * 4);

        // Fixed memory (no scaling)
        let opts_fixed = QueueMemoryOptions {
            queue_memory: "200".to_string(),
            queue_memory_per_thread: false,
            queue_memory_limit_mb: None,
        };
        let result_fixed = opts_fixed.compute_memory_limit(8).unwrap();
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

        let result = legacy_opts.compute_memory_limit(4).unwrap();
        assert_eq!(result, 2048 * 1024 * 1024); // Should use legacy value, no scaling
    }

    #[test]
    fn test_queue_memory_options_human_readable() {
        let opts = QueueMemoryOptions {
            queue_memory: "2GB".to_string(),
            queue_memory_per_thread: false,
            queue_memory_limit_mb: None,
        };

        let result = opts.compute_memory_limit(4).unwrap();
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
}
