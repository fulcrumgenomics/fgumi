//! `Duplex` command implementation.
//!
//! This command takes grouped reads (with MI tags ending in /A and /B) and generates
//! duplex consensus reads by performing two-stage consensus calling:
//! 1. Single-strand consensus for /A and /B reads separately
//! 2. Duplex consensus from paired single-strand consensuses

use anyhow::{Context, Result};
use clap::Parser;
use fgoxide::io::DelimFile;
use fgumi_lib::bam_io::{
    create_bam_reader, create_bam_reader_for_pipeline, create_bam_writer,
    create_optional_bam_writer,
};

use super::common::{
    BamIoOptions, CompressionOptions, ConsensusCallingOptions, OverlappingConsensusOptions,
    ReadGroupOptions, RejectsOptions, SchedulerOptions, StatsOptions, ThreadingOptions,
};
use crate::commands::consensus_runner::{
    ConsensusStatsOps, create_unmapped_consensus_header, log_overlapping_stats,
};
use crossbeam_queue::SegQueue;
use fgumi_lib::consensus_caller::{
    ConsensusCaller, ConsensusCallingStats, make_prefix_from_header,
};
use fgumi_lib::duplex_consensus_caller::DuplexConsensusCaller;
use fgumi_lib::grouper::{PositionMiGroup, PositionMiGroupBatch, PositionMiGrouper};
use fgumi_lib::logging::{OperationTimer, log_consensus_summary};
use fgumi_lib::mi_group::{MiGroup, MiGroupBatch, MiGroupIteratorWithTransform, MiGrouper};
use fgumi_lib::overlapping_consensus::{
    AgreementStrategy, CorrectionStats, DisagreementStrategy, OverlappingBasesConsensusCaller,
    apply_overlapping_consensus,
};
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::umi::extract_mi_base;
use fgumi_lib::unified_pipeline::{
    BamPipelineConfig, Grouper, MemoryEstimate, run_bam_pipeline_from_reader,
    serialize_bam_records_into,
};
use fgumi_lib::validation::{optional_string_to_tag, validate_file_exists};
use log::info;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::io;
use std::sync::Arc;

use super::command::Command;

// ============================================================================
// Types for 7-step pipeline processing
// ============================================================================

/// Result from processing a batch of MI groups through duplex consensus calling.
struct DuplexProcessedBatch {
    /// Consensus reads to write to output BAM
    consensus_reads: Vec<RecordBuf>,
    /// Rejected reads (written to rejects file if enabled)
    rejects: Vec<RecordBuf>,
    /// Number of MI groups in this batch
    groups_count: u64,
    /// Consensus calling statistics for this batch
    stats: ConsensusCallingStats,
    /// Overlapping correction stats for this batch (if enabled)
    overlapping_stats: Option<CorrectionStats>,
}

impl MemoryEstimate for DuplexProcessedBatch {
    fn estimate_heap_size(&self) -> usize {
        self.consensus_reads
            .iter()
            .chain(self.rejects.iter())
            .map(MemoryEstimate::estimate_heap_size)
            .sum()
    }
}

/// Metrics collected from each batch during parallel processing.
#[derive(Default)]
struct CollectedDuplexMetrics {
    /// Consensus calling statistics
    stats: ConsensusCallingStats,
    /// Overlapping consensus stats (if enabled)
    overlapping_stats: Option<CorrectionStats>,
    /// Number of MI groups processed
    groups_processed: u64,
    /// Rejected reads for deferred writing
    rejects: Vec<RecordBuf>,
}

/// Family size metrics collected from each batch during parallel processing.
///
/// This struct collects metrics about family sizes for duplex sequencing QC.
/// It is used for lock-free aggregation across threads via SegQueue.
#[derive(Default, Clone)]
struct CollectedDuplexFamilyMetrics {
    /// CS (Coordinate+Strand) family size histogram.
    /// Maps family_size -> count of families with that size.
    cs_family_sizes: std::collections::HashMap<usize, u64>,
    /// SS (Single-Strand) family size histogram.
    /// Maps family_size -> count of families with that size.
    ss_family_sizes: std::collections::HashMap<usize, u64>,
    /// DS (Double-Strand) family size histogram.
    /// Maps family_size -> count of families with that size.
    ds_family_sizes: std::collections::HashMap<usize, u64>,
    /// Duplex family size histogram.
    /// Maps (ab_size, ba_size) -> count of families with those sizes.
    duplex_family_sizes: std::collections::HashMap<(usize, usize), u64>,
    /// Number of read pairs processed.
    read_pairs: u64,
    /// Number of DS families that qualify as duplexes (meet min reads thresholds).
    ds_duplexes: u64,
}

impl CollectedDuplexFamilyMetrics {
    /// Merge another set of family metrics into this one.
    fn merge(&mut self, other: &CollectedDuplexFamilyMetrics) {
        for (size, count) in &other.cs_family_sizes {
            *self.cs_family_sizes.entry(*size).or_insert(0) += count;
        }
        for (size, count) in &other.ss_family_sizes {
            *self.ss_family_sizes.entry(*size).or_insert(0) += count;
        }
        for (size, count) in &other.ds_family_sizes {
            *self.ds_family_sizes.entry(*size).or_insert(0) += count;
        }
        for (sizes, count) in &other.duplex_family_sizes {
            *self.duplex_family_sizes.entry(*sizes).or_insert(0) += count;
        }
        self.read_pairs += other.read_pairs;
        self.ds_duplexes += other.ds_duplexes;
    }

    /// Record a CS (Coordinate+Strand) family.
    fn record_cs_family(&mut self, size: usize) {
        *self.cs_family_sizes.entry(size).or_insert(0) += 1;
    }

    /// Record an SS (Single-Strand) family.
    fn record_ss_family(&mut self, size: usize) {
        *self.ss_family_sizes.entry(size).or_insert(0) += 1;
    }

    /// Record a DS (Double-Strand) family.
    fn record_ds_family(&mut self, size: usize) {
        *self.ds_family_sizes.entry(size).or_insert(0) += 1;
    }

    /// Record a duplex family with AB and BA sizes.
    /// Normalizes so ab_size >= ba_size.
    fn record_duplex_family(&mut self, ab_size: usize, ba_size: usize) {
        let (ab, ba) = if ab_size >= ba_size { (ab_size, ba_size) } else { (ba_size, ab_size) };
        *self.duplex_family_sizes.entry((ab, ba)).or_insert(0) += 1;
    }
}

/// Call duplex consensus reads from grouped reads with /A and /B MI tags
#[derive(Parser, Debug)]
#[command(
    name = "duplex",
    about = "\x1b[38;5;180m[CONSENSUS]\x1b[0m      \x1b[36mCall duplex consensus sequences from UMI-grouped reads\x1b[0m",
    long_about = r#"
Calls duplex consensus sequences from reads generated from the same double-stranded source molecule. Prior
to running this tool, reads must have been grouped with `group` using the `paired` strategy. Doing
so will apply (by default) MI tags to all reads of the form `*/A` and `*/B` where the /A and /B suffixes
with the same identifier denote reads that are derived from opposite strands of the same source duplex molecule.

Reads from the same unique molecule are first partitioned by source strand and assembled into single
strand consensus molecules as described by the simplex command. Subsequently, for molecules that
have at least one observation of each strand, duplex consensus reads are assembled by combining the evidence
from the two single strand consensus reads.

Because of the nature of duplex sequencing, this tool does not support fragment reads - if found in the
input they are ignored. Similarly, read pairs for which consensus reads cannot be generated for one or
other read (R1 or R2) are omitted from the output.

The consensus reads produced are unaligned, due to the difficulty and error-prone nature of inferring the consensus
alignment. Consensus reads should therefore be aligned after, which should not be too expensive as likely there
are far fewer consensus reads than input raw reads.

Consensus reads have a number of additional optional tags set in the resulting BAM file. The tag names follow
a pattern where the first letter (a, b or c) denotes that the tag applies to the first single strand consensus (a),
second single-strand consensus (b) or the final duplex consensus (c). The second letter is intended to capture
the meaning of the tag (e.g. d=depth, m=min depth, e=errors/error-rate) and is upper case for values that are
one per read and lower case for values that are one per base.

The tags break down into those that are single-valued per read:

  consensus depth      [aD,bD,cD] (int)  : the maximum depth of raw reads at any point in the consensus reads
  consensus min depth  [aM,bM,cM] (int)  : the minimum depth of raw reads at any point in the consensus reads
  consensus error rate [aE,bE,cE] (float): the fraction of bases in raw reads disagreeing with the final consensus calls

And those that have a value per base (duplex values are not generated, but can be generated by summing):

  consensus depth  [ad,bd] (short[]): the count of bases contributing to each single-strand consensus read at each position
  consensus errors [ae,be] (short[]): the count of bases from raw reads disagreeing with the final single-strand consensus base
  consensus bases  [ac,bc] (string) : the single-strand consensus bases
  consensus quals  [aq,bq] (string) : the single-strand consensus qualities

The per base depths and errors are both capped at 32,767. In all cases no-calls (Ns) and bases below the
min-input-base-quality are not counted in tag value calculations.

The --min-reads option can take 1-3 values similar to `filter`. For example:

  fgumi duplex ... --min-reads 10,5,3

If fewer than three values are supplied, the last value is repeated (i.e. `5,4` -> `5 4 4` and `1` -> `1 1 1`). The
first value applies to the final consensus read, the second value to one single-strand consensus, and the last
value to the other single-strand consensus. It is required that if values two and three differ,
the more stringent value comes earlier.
"#
)]
pub struct Duplex {
    /// Input and output BAM files
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Optional output BAM file for rejected reads
    #[command(flatten)]
    pub rejects_opts: RejectsOptions,

    /// Optional output file for consensus calling statistics
    #[command(flatten)]
    pub stats_opts: StatsOptions,

    /// Read group options
    #[command(flatten)]
    pub read_group: ReadGroupOptions,

    /// Consensus calling options
    #[command(flatten)]
    pub consensus: ConsensusCallingOptions,

    /// Overlapping consensus options
    #[command(flatten)]
    pub overlapping: OverlappingConsensusOptions,

    /// Threading options for parallel processing
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Minimum reads for consensus calling.
    /// Can specify 1-3 values: \[duplex\] or \[duplex, AB/BA\] or \[duplex, AB, BA\]
    #[arg(short = 'M', long = "min-reads", value_delimiter = ',', default_value = "1")]
    pub min_reads: Vec<usize>,

    /// Maximum reads per strand (downsample if exceeded)
    #[arg(long = "max-reads-per-strand")]
    pub max_reads_per_strand: Option<usize>,

    /// Cellular barcode tag for filtering (e.g., CB for 10X)
    #[arg(short = 'c', long = "cell-tag", default_value = "CB")]
    pub cell_tag: Option<String>,

    /// Scheduler and pipeline stats options
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Output prefix for duplex metrics files. When specified, collects family size
    /// distributions, yield metrics, and UMI statistics during consensus calling.
    /// Generates: PREFIX.family_sizes.txt, PREFIX.duplex_family_sizes.txt,
    /// PREFIX.duplex_yield_metrics.txt, PREFIX.umi_counts.txt, and PREFIX.duplex_qc.pdf
    #[arg(long = "metrics-output")]
    pub metrics_output: Option<std::path::PathBuf>,

    /// Sample description for metric plots (used in PDF titles).
    /// Defaults to the output BAM filename if not specified.
    #[arg(long = "description")]
    pub description: Option<String>,
}

impl Command for Duplex {
    /// Executes the duplex consensus calling pipeline.
    ///
    /// This method processes grouped reads with MI tags ending in /A and /B suffixes to generate
    /// duplex consensus reads. The process involves:
    /// 1. Parsing and validating the cell tag parameter
    /// 2. Creating a duplex consensus caller with specified parameters
    /// 3. Opening input BAM and creating output BAM writer
    /// 4. Processing templates through the template iterator
    /// 5. Calling consensus for each template
    /// 6. Writing consensus reads to output
    /// 7. Optionally writing rejected reads to a separate BAM
    /// 8. Logging statistics
    ///
    /// # Returns
    ///
    /// `Ok(())` on success.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Cell tag is not exactly 2 characters
    /// - Input BAM cannot be opened or read
    /// - Output BAM cannot be created or written
    /// - Template processing encounters errors
    /// - Consensus calling fails
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use fgumi::commands::duplex::Duplex;
    /// # use fgumi::commands::command::Command;
    /// # use fgumi::commands::common::{
    /// #     BamIoOptions, ConsensusCallingOptions, OverlappingConsensusOptions,
    /// #     ReadGroupOptions, RejectsOptions, SchedulerOptions, StatsOptions, ThreadingOptions,
    /// # };
    /// # use std::path::PathBuf;
    /// let duplex = Duplex {
    ///     io: BamIoOptions {
    ///         input: PathBuf::from("grouped.bam"),
    ///         output: PathBuf::from("duplex.bam"),
    ///     },
    ///     rejects_opts: RejectsOptions::default(),
    ///     stats_opts: StatsOptions::default(),
    ///     read_group: ReadGroupOptions {
    ///         read_name_prefix: Some("consensus".to_string()),
    ///         read_group_id: "duplex".to_string(),
    ///     },
    ///     consensus: ConsensusCallingOptions {
    ///         output_per_base_tags: false,
    ///         ..ConsensusCallingOptions::default()
    ///     },
    ///     overlapping: OverlappingConsensusOptions::default(),
    ///     threading: ThreadingOptions::none(),
    ///     compression: CompressionOptions::default(),
    ///     min_reads: vec![1],
    ///     max_reads_per_strand: None,
    ///     cell_tag: None,
    ///     scheduler_opts: SchedulerOptions::default(),
    ///     metrics_output: None,
    ///     description: None,
    /// };
    ///
    /// duplex.execute("test")?;
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    fn execute(&self, command_line: &str) -> Result<()> {
        // Start timing
        let timer = OperationTimer::new("Calling duplex consensus");

        // Validate input file exists
        validate_file_exists(&self.io.input, "Input BAM")?;

        // Get threading configuration (duplex is worker-heavy with consensus calling)
        let reader_threads = self.threading.num_threads();
        let worker_threads = self.threading.num_threads();
        let writer_threads = self.threading.num_threads();

        info!("Duplex");
        info!("  Input: {}", self.io.input.display());
        info!("  Output: {}", self.io.output.display());
        info!("  Min reads: {:?}", self.min_reads);
        info!("  Min base quality: {}", self.consensus.min_input_base_quality);
        info!("  Output per-base tags: {}", self.consensus.output_per_base_tags);
        info!("  Worker threads: {worker_threads}");
        info!("  Reader threads: {reader_threads}");
        info!("  Trim reads: {}", self.consensus.trim);
        info!("  Max reads per strand: {:?}", self.max_reads_per_strand);
        info!("  Cell tag: {:?}", self.cell_tag);
        info!(
            "  Consensus call overlapping bases: {}",
            self.overlapping.consensus_call_overlapping_bases
        );

        // Open input BAM using streaming-capable reader for pipeline use
        let (reader, header) = create_bam_reader_for_pipeline(&self.io.input)?;

        // Add @PG record with PP chaining to input's last program
        let header = fgumi_lib::header::add_pg_record(
            header,
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        // Use library name from header if no prefix is specified (like fgbio)
        let read_name_prefix = self
            .read_group
            .read_name_prefix
            .clone()
            .unwrap_or_else(|| make_prefix_from_header(&header));

        // Parse cell_tag if provided
        let cell_tag = optional_string_to_tag(self.cell_tag.as_deref(), "cell-tag")?;

        // Enable rejects tracking if rejects file is specified
        let track_rejects = self.rejects_opts.is_enabled();

        // Track overlapping consensus settings (callers created per-thread in threaded mode)
        let overlapping_enabled = self.overlapping.consensus_call_overlapping_bases;
        if overlapping_enabled {
            info!("Overlapping consensus calling enabled");
        }

        // Process reads grouped by MI tag (streaming approach)
        info!("Processing reads...");

        // ============================================================
        // --threads N mode: Use 7-step unified pipeline
        // None: Use single-threaded fast path
        // ============================================================
        // IMPORTANT: Check this BEFORE creating any output file handles to avoid
        // file handle conflicts that can corrupt the output.

        // If metrics output is requested but threading is not enabled,
        // automatically enable threading (default to 4 threads) since
        // metrics collection requires position-based grouping
        let effective_threads = if self.metrics_output.is_some() && self.threading.threads.is_none()
        {
            info!("Metrics collection requires threading; using 4 threads");
            Some(4)
        } else {
            self.threading.threads
        };

        if let Some(threads) = effective_threads {
            let result = self.execute_threads_mode(
                threads,
                reader,
                header.clone(),
                read_name_prefix.clone(),
                track_rejects,
                command_line,
            );
            timer.log_completion(0); // Completion logged in execute_threads_mode
            return result;
        }

        // Drop the reader for single-threaded mode - we use MiGroupIterator which needs its own reader
        drop(reader);

        // ============================================================
        // For non-pipeline modes, create output writers here
        // ============================================================

        // Create output BAM writer with multi-threaded BGZF compression
        let mut writer = create_bam_writer(
            &self.io.output,
            &header,
            writer_threads,
            self.compression.compression_level,
        )?;

        // Create optional rejects writer
        let mut rejects_writer = create_optional_bam_writer(
            self.rejects_opts.rejects.as_ref(),
            &header,
            writer_threads,
            self.compression.compression_level,
        )?;

        // Create duplex consensus caller
        // Note: Threading is handled here in duplex.rs, not in DuplexConsensusCaller
        let mut consensus_caller = DuplexConsensusCaller::new(
            read_name_prefix.clone(),
            self.read_group.read_group_id.clone(),
            self.min_reads.clone(),
            self.consensus.min_input_base_quality,
            self.consensus.output_per_base_tags,
            self.consensus.trim,
            self.max_reads_per_strand,
            cell_tag,
            track_rejects,
            self.consensus.error_rate_pre_umi,
            self.consensus.error_rate_post_umi,
        )?;

        // Accumulator for overlapping stats from parallel processing
        let mut merged_overlapping_stats = CorrectionStats::new();

        // Track progress and aggregate stats
        let mut consensus_count = 0usize;
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        // Create the MI group iterator for single-threaded streaming
        let (mut reader, _) = create_bam_reader(&self.io.input, 1)?;
        // Create filtered record iterator
        // Filter out secondary and supplementary reads, keep only mapped or mate-mapped
        let filtered_iter = reader.record_bufs(&header).filter_map(|result| {
            match result {
                Err(e) => Some(Err(e.into())),
                Ok(record) => {
                    let flags = record.flags();
                    // Skip secondary and supplementary reads
                    if flags.is_secondary() || flags.is_supplementary() {
                        return None;
                    }
                    // Only keep mapped reads or reads with mapped mates
                    let is_mapped = !flags.is_unmapped();
                    let has_mapped_mate = flags.is_segmented() && !flags.is_mate_unmapped();
                    if is_mapped || has_mapped_mate { Some(Ok(record)) } else { None }
                }
            }
        });

        // Group by base MI (strip /A, /B suffix) for streaming
        let mi_group_iter = MiGroupIteratorWithTransform::new(filtered_iter, "MI", |mi| {
            extract_mi_base(mi).to_string()
        });
        // Single-threaded processing
        // Create overlapping consensus caller for single-threaded mode
        let mut overlapping_caller = if overlapping_enabled {
            Some(OverlappingBasesConsensusCaller::new(
                AgreementStrategy::Consensus,
                DisagreementStrategy::Consensus,
            ))
        } else {
            None
        };

        for group_result in mi_group_iter {
            let (_base_mi, mut reads) = group_result.context("Failed to read MI group")?;

            // Apply overlapping consensus if enabled (modifies reads in-place)
            // Skip if group doesn't have both strands - no duplex possible anyway
            if let Some(ref mut oc) = overlapping_caller {
                if DuplexConsensusCaller::has_both_strands(&reads) {
                    apply_overlapping_consensus(&mut reads, oc)?;
                }
            }

            // Call consensus for this group
            let consensus_reads = consensus_caller.consensus_reads_from_sam_records(reads)?;

            // Write consensus reads immediately
            let batch_size = consensus_reads.len();
            for consensus_read in &consensus_reads {
                writer.write_alignment_record(&header, consensus_read)?;
                consensus_count += 1;
            }
            progress.log_if_needed(batch_size as u64);
        }

        // Collect stats from single-threaded caller and merge overlapping stats
        let merged_stats = consensus_caller.statistics();
        if let Some(ref oc) = overlapping_caller {
            merged_overlapping_stats.merge(oc.stats());
        }

        // Write rejected reads if tracking is enabled
        if let Some(ref mut rw) = rejects_writer {
            let rejected_reads = consensus_caller.rejected_reads();
            for read in rejected_reads {
                rw.write_alignment_record(&header, read)?;
            }
            info!("Wrote {} rejected reads", rejected_reads.len());
        }

        progress.log_final();

        // Finish the buffered writer (flush remaining records and wait for writer thread)
        writer.into_inner().finish().context("Failed to finish output BAM")?;

        // Log overlapping consensus statistics if enabled
        if overlapping_enabled {
            log_overlapping_stats(&merged_overlapping_stats);
        }

        // Convert stats to metrics and log
        let mut metrics = merged_stats.to_metrics();
        // Use local count since we tracked it during streaming
        metrics.consensus_reads = consensus_count as u64;
        log_consensus_summary(&metrics);

        // Write statistics file if requested
        if let Some(stats_path) = &self.stats_opts.stats {
            DelimFile::default().write_tsv(stats_path, [metrics]).map_err(|e| {
                anyhow::anyhow!("Failed to write statistics to {}: {}", stats_path.display(), e)
            })?;
            info!("Statistics written to: {}", stats_path.display());
        }

        // Log timing completion
        timer.log_completion(consensus_count as u64);

        info!("Done!");
        Ok(())
    }
}

impl Duplex {
    /// Execute using 7-step unified pipeline with --threads.
    ///
    /// This method is called when `--threads N` is specified with N > 1.
    /// It uses the lock-free 7-step unified pipeline for maximum performance.
    fn execute_threads_mode(
        &self,
        num_threads: usize,
        reader: Box<dyn std::io::Read + Send>,
        input_header: Header,
        read_name_prefix: String,
        track_rejects: bool,
        command_line: &str,
    ) -> Result<()> {
        // If metrics output is requested, use the metrics-enabled pipeline path
        if self.metrics_output.is_some() {
            return self.execute_threads_mode_with_metrics(
                num_threads,
                reader,
                input_header,
                read_name_prefix,
                track_rejects,
                command_line,
            );
        }

        info!("Using 7-step unified pipeline with {num_threads} threads");

        // Create output header (for duplex, output is unmapped like simplex)
        let output_header = create_unmapped_consensus_header(
            &input_header,
            &self.read_group.read_group_id,
            "Read group",
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        // Configure pipeline
        let mut pipeline_config =
            BamPipelineConfig::auto_tuned(num_threads, self.compression.compression_level);
        pipeline_config.pipeline.scheduler_strategy = self.scheduler_opts.strategy();
        if self.scheduler_opts.collect_stats() {
            pipeline_config.pipeline = pipeline_config.pipeline.with_stats(true);
        }
        pipeline_config.pipeline.deadlock_timeout_secs =
            self.scheduler_opts.deadlock_timeout_secs();
        pipeline_config.pipeline.deadlock_recover_enabled =
            self.scheduler_opts.deadlock_recover_enabled();

        // Lock-free metrics collection
        let collected_metrics: Arc<SegQueue<CollectedDuplexMetrics>> = Arc::new(SegQueue::new());
        let collected_metrics_for_serialize = Arc::clone(&collected_metrics);

        // Capture configuration for closures
        let min_reads = self.min_reads.clone();
        let min_input_base_quality = self.consensus.min_input_base_quality;
        let output_per_base_tags = self.consensus.output_per_base_tags;
        let trim = self.consensus.trim;
        let max_reads_per_strand = self.max_reads_per_strand;
        let error_rate_pre_umi = self.consensus.error_rate_pre_umi;
        let error_rate_post_umi = self.consensus.error_rate_post_umi;
        let overlapping_enabled = self.overlapping.consensus_call_overlapping_bases;
        let read_group_id = self.read_group.read_group_id.clone();
        let cell_tag = optional_string_to_tag(self.cell_tag.as_deref(), "cell-tag")?;
        let batch_size = 100; // MI groups per batch

        // Record filter for duplex: skip secondary/supplementary, keep mapped or mate-mapped
        let record_filter = |record: &RecordBuf| -> bool {
            let flags = record.flags();
            // Skip secondary and supplementary reads
            if flags.is_secondary() || flags.is_supplementary() {
                return false;
            }
            // Only keep mapped reads or reads with mapped mates
            let is_mapped = !flags.is_unmapped();
            let has_mapped_mate = flags.is_segmented() && !flags.is_mate_unmapped();
            is_mapped || has_mapped_mate
        };

        // MI transform: strip /A and /B suffixes for duplex grouping
        let mi_transform = |mi: &str| extract_mi_base(mi).to_string();

        // Clone input_header before pipeline (needed for rejects writing)
        let rejects_header = input_header.clone();

        // Run the 7-step pipeline with the already-opened reader (supports streaming)
        let groups_processed = run_bam_pipeline_from_reader(
            pipeline_config,
            reader,
            input_header,
            &self.io.output,
            Some(output_header.clone()),
            // ========== grouper_fn: Create MiGrouper with filter and transform ==========
            move |_header: &Header| {
                Box::new(MiGrouper::with_filter_and_transform(
                    "MI",
                    batch_size,
                    record_filter,
                    mi_transform,
                )) as Box<dyn Grouper<Group = MiGroupBatch> + Send>
            },
            // ========== process_fn: Duplex consensus calling ==========
            move |batch: MiGroupBatch| -> io::Result<DuplexProcessedBatch> {
                // Create per-thread duplex consensus caller
                let mut caller = DuplexConsensusCaller::new(
                    read_name_prefix.clone(),
                    read_group_id.clone(),
                    min_reads.clone(),
                    min_input_base_quality,
                    output_per_base_tags,
                    trim,
                    max_reads_per_strand,
                    cell_tag,
                    track_rejects,
                    error_rate_pre_umi,
                    error_rate_post_umi,
                )
                .map_err(|e| {
                    io::Error::other(format!("Failed to create DuplexConsensusCaller: {e}"))
                })?;

                // Create overlapping caller if enabled
                let mut overlapping_caller = if overlapping_enabled {
                    Some(OverlappingBasesConsensusCaller::new(
                        AgreementStrategy::Consensus,
                        DisagreementStrategy::Consensus,
                    ))
                } else {
                    None
                };

                let mut all_consensus = Vec::new();
                let mut all_rejects = Vec::new();
                let mut batch_stats = ConsensusCallingStats::new();
                let mut batch_overlapping = CorrectionStats::new();
                let groups_count = batch.groups.len() as u64;

                for MiGroup { mi, records: mut group_reads } in batch.groups {
                    caller.clear();

                    // Apply overlapping consensus if enabled (modifies in-place)
                    // Skip if group doesn't have both strands - no duplex possible anyway
                    if let Some(ref mut oc) = overlapping_caller {
                        if DuplexConsensusCaller::has_both_strands(&group_reads) {
                            oc.reset_stats();
                            if apply_overlapping_consensus(&mut group_reads, oc).is_err() {
                                continue;
                            }
                            batch_overlapping.merge(oc.stats());
                        }
                    }

                    // Call duplex consensus
                    match caller.consensus_reads_from_sam_records(group_reads) {
                        Ok(consensus_reads) => {
                            all_consensus.extend(consensus_reads);
                            batch_stats.merge(&caller.statistics());
                            if track_rejects {
                                all_rejects.extend(caller.take_rejected_reads());
                            }
                        }
                        Err(e) => {
                            log::warn!("Duplex consensus error for MI {mi}: {e}");
                        }
                    }
                }

                Ok(DuplexProcessedBatch {
                    consensus_reads: all_consensus,
                    rejects: all_rejects,
                    groups_count,
                    stats: batch_stats,
                    overlapping_stats: if overlapping_enabled {
                        Some(batch_overlapping)
                    } else {
                        None
                    },
                })
            },
            // ========== serialize_fn: Serialize + collect metrics ==========
            move |processed: DuplexProcessedBatch,
                  header: &Header,
                  output: &mut Vec<u8>|
                  -> io::Result<u64> {
                // Collect metrics (lock-free)
                collected_metrics_for_serialize.push(CollectedDuplexMetrics {
                    stats: processed.stats,
                    overlapping_stats: processed.overlapping_stats,
                    groups_processed: processed.groups_count,
                    rejects: processed.rejects,
                });

                // Serialize consensus reads
                serialize_bam_records_into(&processed.consensus_reads, header, output)
            },
        )
        .map_err(|e| anyhow::anyhow!("Pipeline error: {e}"))?;

        // ========== Post-pipeline: Aggregate metrics and write rejects ==========
        let mut total_groups = 0u64;
        let mut merged_stats = ConsensusCallingStats::new();
        let mut merged_overlapping_stats = CorrectionStats::new();
        let mut all_rejects: Vec<RecordBuf> = Vec::new();

        while let Some(metrics) = collected_metrics.pop() {
            total_groups += metrics.groups_processed;
            merged_stats.merge(&metrics.stats);
            if let Some(ref ocs) = metrics.overlapping_stats {
                merged_overlapping_stats.merge(ocs);
            }
            if track_rejects {
                all_rejects.extend(metrics.rejects);
            }
        }

        // Write deferred rejects
        if track_rejects && !all_rejects.is_empty() {
            if let Some(rejects_path) = &self.rejects_opts.rejects {
                let writer_threads = self.threading.num_threads();
                let mut rejects_writer = create_optional_bam_writer(
                    Some(rejects_path),
                    &rejects_header,
                    writer_threads,
                    self.compression.compression_level,
                )?;
                if let Some(ref mut rw) = rejects_writer {
                    for record in &all_rejects {
                        rw.write_alignment_record(&rejects_header, record)
                            .context("Failed to write rejected read")?;
                    }
                    rw.finish(&rejects_header).context("Failed to finish rejects file")?;
                    info!("Wrote {} rejected reads", all_rejects.len());
                }
            }
        }

        // Log overlapping consensus statistics if enabled
        if self.overlapping.consensus_call_overlapping_bases {
            log_overlapping_stats(&merged_overlapping_stats);
        }

        // Log statistics
        info!("Duplex consensus calling complete");
        info!("Total MI groups processed: {total_groups}");
        info!("Total groups processed by pipeline: {groups_processed}");

        let metrics = merged_stats.to_metrics();
        let consensus_count = metrics.consensus_reads;
        log_consensus_summary(&metrics);

        // Write statistics file if requested
        if let Some(stats_path) = &self.stats_opts.stats {
            let kv_metrics = metrics.to_kv_metrics();
            DelimFile::default()
                .write_tsv(stats_path, kv_metrics)
                .with_context(|| format!("Failed to write statistics: {}", stats_path.display()))?;
            info!("Wrote statistics to: {}", stats_path.display());
        }

        info!("Wrote {consensus_count} duplex consensus reads");

        Ok(())
    }

    /// Execute using 7-step unified pipeline with metrics collection.
    ///
    /// This method is called when `--metrics-output` is specified.
    /// It uses `PositionMiGrouper` to group by position first, enabling
    /// collection of CS, SS, and DS family size metrics.
    #[allow(clippy::too_many_lines)]
    fn execute_threads_mode_with_metrics(
        &self,
        num_threads: usize,
        reader: Box<dyn std::io::Read + Send>,
        input_header: Header,
        read_name_prefix: String,
        track_rejects: bool,
        command_line: &str,
    ) -> Result<()> {
        info!("Using 7-step unified pipeline with {num_threads} threads (with metrics collection)");

        // Create output header (for duplex, output is unmapped like simplex)
        let output_header = create_unmapped_consensus_header(
            &input_header,
            &self.read_group.read_group_id,
            "Read group",
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        // Configure pipeline
        let mut pipeline_config =
            BamPipelineConfig::auto_tuned(num_threads, self.compression.compression_level);
        pipeline_config.pipeline.scheduler_strategy = self.scheduler_opts.strategy();
        if self.scheduler_opts.collect_stats() {
            pipeline_config.pipeline = pipeline_config.pipeline.with_stats(true);
        }
        pipeline_config.pipeline.deadlock_timeout_secs =
            self.scheduler_opts.deadlock_timeout_secs();
        pipeline_config.pipeline.deadlock_recover_enabled =
            self.scheduler_opts.deadlock_recover_enabled();

        // Lock-free metrics collection
        let collected_metrics: Arc<SegQueue<CollectedDuplexMetrics>> = Arc::new(SegQueue::new());
        let collected_metrics_for_serialize = Arc::clone(&collected_metrics);

        // Lock-free family metrics collection
        let collected_family_metrics: Arc<SegQueue<CollectedDuplexFamilyMetrics>> =
            Arc::new(SegQueue::new());
        let collected_family_metrics_for_serialize = Arc::clone(&collected_family_metrics);

        // Capture configuration for closures
        let min_reads = self.min_reads.clone();
        let min_input_base_quality = self.consensus.min_input_base_quality;
        let output_per_base_tags = self.consensus.output_per_base_tags;
        let trim = self.consensus.trim;
        let max_reads_per_strand = self.max_reads_per_strand;
        let error_rate_pre_umi = self.consensus.error_rate_pre_umi;
        let error_rate_post_umi = self.consensus.error_rate_post_umi;
        let overlapping_enabled = self.overlapping.consensus_call_overlapping_bases;
        let read_group_id = self.read_group.read_group_id.clone();
        let cell_tag = optional_string_to_tag(self.cell_tag.as_deref(), "cell-tag")?;
        let batch_size = 100; // Position groups per batch

        // Clone input_header before pipeline (needed for rejects writing)
        let rejects_header = input_header.clone();

        // MI tag for extracting strand info
        let mi_tag = noodles::sam::alignment::record::data::field::Tag::from([b'M', b'I']);

        // Run the 7-step pipeline with PositionMiGrouper
        let groups_processed = run_bam_pipeline_from_reader(
            pipeline_config,
            reader,
            input_header,
            &self.io.output,
            Some(output_header.clone()),
            // ========== grouper_fn: Create PositionMiGrouper ==========
            move |_header: &Header| {
                Box::new(PositionMiGrouper::new("MI", batch_size))
                    as Box<dyn Grouper<Group = PositionMiGroupBatch> + Send>
            },
            // ========== process_fn: Duplex consensus calling with metrics ==========
            move |batch: PositionMiGroupBatch| -> io::Result<DuplexProcessedBatchWithMetrics> {
                // Create per-thread duplex consensus caller
                let mut caller = DuplexConsensusCaller::new(
                    read_name_prefix.clone(),
                    read_group_id.clone(),
                    min_reads.clone(),
                    min_input_base_quality,
                    output_per_base_tags,
                    trim,
                    max_reads_per_strand,
                    cell_tag,
                    track_rejects,
                    error_rate_pre_umi,
                    error_rate_post_umi,
                )
                .map_err(|e| {
                    io::Error::other(format!("Failed to create DuplexConsensusCaller: {e}"))
                })?;

                // Create overlapping caller if enabled
                let mut overlapping_caller = if overlapping_enabled {
                    Some(OverlappingBasesConsensusCaller::new(
                        AgreementStrategy::Consensus,
                        DisagreementStrategy::Consensus,
                    ))
                } else {
                    None
                };

                let mut all_consensus = Vec::new();
                let mut all_rejects = Vec::new();
                let mut batch_stats = ConsensusCallingStats::new();
                let mut batch_overlapping = CorrectionStats::new();
                let mut family_metrics = CollectedDuplexFamilyMetrics::default();
                let mut groups_count = 0u64;

                for PositionMiGroup { mi_groups, total_templates, .. } in batch.groups {
                    // Record CS family (all templates at this position)
                    family_metrics.record_cs_family(total_templates);

                    // Process each DS family (grouped by base MI)
                    for (_base_mi, group_records) in mi_groups {
                        groups_count += 1;
                        caller.clear();

                        // Count SS families (by strand suffix /A or /B)
                        let mut a_qnames = std::collections::HashSet::new();
                        let mut b_qnames = std::collections::HashSet::new();

                        for record in &group_records {
                            if let Some(val) = record.data().get(&mi_tag) {
                                use noodles::sam::alignment::record_buf::data::field::Value;
                                if let Value::String(s) = val {
                                    let qname_hash = record
                                        .name()
                                        .map(|n| {
                                            use std::hash::{Hash, Hasher};
                                            let mut hasher =
                                                std::collections::hash_map::DefaultHasher::new();
                                            <_ as AsRef<[u8]>>::as_ref(n).hash(&mut hasher);
                                            hasher.finish()
                                        })
                                        .unwrap_or(0);
                                    if s.ends_with(b"/A") {
                                        a_qnames.insert(qname_hash);
                                    } else if s.ends_with(b"/B") {
                                        b_qnames.insert(qname_hash);
                                    }
                                }
                            }
                        }

                        let a_count = a_qnames.len();
                        let b_count = b_qnames.len();
                        let ds_size = a_count + b_count;

                        // Record DS family size
                        if ds_size > 0 {
                            family_metrics.record_ds_family(ds_size);
                        }

                        // Record SS family sizes
                        if a_count > 0 {
                            family_metrics.record_ss_family(a_count);
                        }
                        if b_count > 0 {
                            family_metrics.record_ss_family(b_count);
                        }

                        // Record duplex family sizes (AB x BA)
                        if a_count > 0 && b_count > 0 {
                            family_metrics.record_duplex_family(a_count, b_count);
                            family_metrics.ds_duplexes += 1;
                        }

                        // Track read pairs
                        family_metrics.read_pairs += (a_count + b_count) as u64;

                        // Apply overlapping consensus if enabled
                        let mut group_reads = group_records;
                        if let Some(ref mut oc) = overlapping_caller {
                            if DuplexConsensusCaller::has_both_strands(&group_reads) {
                                oc.reset_stats();
                                if apply_overlapping_consensus(&mut group_reads, oc).is_err() {
                                    continue;
                                }
                                batch_overlapping.merge(oc.stats());
                            }
                        }

                        // Call duplex consensus
                        match caller.consensus_reads_from_sam_records(group_reads) {
                            Ok(consensus_reads) => {
                                all_consensus.extend(consensus_reads);
                                batch_stats.merge(&caller.statistics());
                                if track_rejects {
                                    all_rejects.extend(caller.take_rejected_reads());
                                }
                            }
                            Err(e) => {
                                log::warn!("Duplex consensus error: {e}");
                            }
                        }
                    }
                }

                Ok(DuplexProcessedBatchWithMetrics {
                    consensus_reads: all_consensus,
                    rejects: all_rejects,
                    groups_count,
                    stats: batch_stats,
                    overlapping_stats: if overlapping_enabled {
                        Some(batch_overlapping)
                    } else {
                        None
                    },
                    family_metrics,
                })
            },
            // ========== serialize_fn: Serialize + collect metrics ==========
            move |processed: DuplexProcessedBatchWithMetrics,
                  header: &Header,
                  output: &mut Vec<u8>|
                  -> io::Result<u64> {
                // Collect consensus metrics (lock-free)
                collected_metrics_for_serialize.push(CollectedDuplexMetrics {
                    stats: processed.stats,
                    overlapping_stats: processed.overlapping_stats,
                    groups_processed: processed.groups_count,
                    rejects: processed.rejects,
                });

                // Collect family metrics (lock-free)
                collected_family_metrics_for_serialize.push(processed.family_metrics);

                // Serialize consensus reads
                serialize_bam_records_into(&processed.consensus_reads, header, output)
            },
        )
        .map_err(|e| anyhow::anyhow!("Pipeline error: {e}"))?;

        // ========== Post-pipeline: Aggregate metrics ==========
        let mut total_groups = 0u64;
        let mut merged_stats = ConsensusCallingStats::new();
        let mut merged_overlapping_stats = CorrectionStats::new();
        let mut all_rejects: Vec<RecordBuf> = Vec::new();

        while let Some(metrics) = collected_metrics.pop() {
            total_groups += metrics.groups_processed;
            merged_stats.merge(&metrics.stats);
            if let Some(ref ocs) = metrics.overlapping_stats {
                merged_overlapping_stats.merge(ocs);
            }
            if track_rejects {
                all_rejects.extend(metrics.rejects);
            }
        }

        // Aggregate family metrics
        let mut merged_family_metrics = CollectedDuplexFamilyMetrics::default();
        while let Some(family_metrics) = collected_family_metrics.pop() {
            merged_family_metrics.merge(&family_metrics);
        }

        // Write deferred rejects
        if track_rejects && !all_rejects.is_empty() {
            if let Some(rejects_path) = &self.rejects_opts.rejects {
                let writer_threads = self.threading.num_threads();
                let mut rejects_writer = create_optional_bam_writer(
                    Some(rejects_path),
                    &rejects_header,
                    writer_threads,
                    self.compression.compression_level,
                )?;
                if let Some(ref mut rw) = rejects_writer {
                    for record in &all_rejects {
                        rw.write_alignment_record(&rejects_header, record)
                            .context("Failed to write rejected read")?;
                    }
                    rw.finish(&rejects_header).context("Failed to finish rejects file")?;
                    info!("Wrote {} rejected reads", all_rejects.len());
                }
            }
        }

        // Log overlapping consensus statistics if enabled
        if self.overlapping.consensus_call_overlapping_bases {
            log_overlapping_stats(&merged_overlapping_stats);
        }

        // Log statistics
        info!("Duplex consensus calling complete");
        info!("Total MI groups processed: {total_groups}");
        info!("Total groups processed by pipeline: {groups_processed}");

        let metrics = merged_stats.to_metrics();
        let consensus_count = metrics.consensus_reads;
        log_consensus_summary(&metrics);

        // Write statistics file if requested
        if let Some(stats_path) = &self.stats_opts.stats {
            let kv_metrics = metrics.to_kv_metrics();
            DelimFile::default()
                .write_tsv(stats_path, kv_metrics)
                .with_context(|| format!("Failed to write statistics: {}", stats_path.display()))?;
            info!("Wrote statistics to: {}", stats_path.display());
        }

        // Write family metrics files
        if let Some(ref metrics_prefix) = self.metrics_output {
            self.write_family_metrics(metrics_prefix, &merged_family_metrics)?;
        }

        info!("Wrote {consensus_count} duplex consensus reads");

        Ok(())
    }

    /// Write family metrics files.
    fn write_family_metrics(
        &self,
        prefix: &std::path::Path,
        metrics: &CollectedDuplexFamilyMetrics,
    ) -> Result<()> {
        use std::io::Write;

        // Write family_sizes.txt
        let family_sizes_path = prefix.with_extension("family_sizes.txt");
        let mut file = std::fs::File::create(&family_sizes_path)
            .with_context(|| format!("Failed to create {}", family_sizes_path.display()))?;
        writeln!(file, "family_type\tfamily_size\tcount")?;

        // Sort and write CS family sizes
        let mut cs_sizes: Vec<_> = metrics.cs_family_sizes.iter().collect();
        cs_sizes.sort_by_key(|(k, _)| *k);
        for (size, count) in cs_sizes {
            writeln!(file, "CS\t{size}\t{count}")?;
        }

        // Sort and write SS family sizes
        let mut ss_sizes: Vec<_> = metrics.ss_family_sizes.iter().collect();
        ss_sizes.sort_by_key(|(k, _)| *k);
        for (size, count) in ss_sizes {
            writeln!(file, "SS\t{size}\t{count}")?;
        }

        // Sort and write DS family sizes
        let mut ds_sizes: Vec<_> = metrics.ds_family_sizes.iter().collect();
        ds_sizes.sort_by_key(|(k, _)| *k);
        for (size, count) in ds_sizes {
            writeln!(file, "DS\t{size}\t{count}")?;
        }

        info!("Wrote family sizes to: {}", family_sizes_path.display());

        // Write duplex_family_sizes.txt
        let duplex_sizes_path = prefix.with_extension("duplex_family_sizes.txt");
        let mut file = std::fs::File::create(&duplex_sizes_path)
            .with_context(|| format!("Failed to create {}", duplex_sizes_path.display()))?;
        writeln!(file, "ab_size\tba_size\tcount")?;

        let mut duplex_sizes: Vec<_> = metrics.duplex_family_sizes.iter().collect();
        duplex_sizes.sort_by_key(|((a, b), _)| (*a, *b));
        for ((ab, ba), count) in duplex_sizes {
            writeln!(file, "{ab}\t{ba}\t{count}")?;
        }

        info!("Wrote duplex family sizes to: {}", duplex_sizes_path.display());

        // Log summary
        let total_cs: u64 = metrics.cs_family_sizes.values().sum();
        let total_ss: u64 = metrics.ss_family_sizes.values().sum();
        let total_ds: u64 = metrics.ds_family_sizes.values().sum();
        let total_duplexes: u64 = metrics.duplex_family_sizes.values().sum();

        info!("Family metrics summary:");
        info!("  CS families: {total_cs}");
        info!("  SS families: {total_ss}");
        info!("  DS families: {total_ds}");
        info!("  Duplex families: {total_duplexes}");
        info!("  Read pairs: {}", metrics.read_pairs);

        Ok(())
    }
}

/// Result from processing a batch of position groups with metrics.
struct DuplexProcessedBatchWithMetrics {
    /// Consensus reads to write to output BAM
    consensus_reads: Vec<RecordBuf>,
    /// Rejected reads (written to rejects file if enabled)
    rejects: Vec<RecordBuf>,
    /// Number of MI groups in this batch
    groups_count: u64,
    /// Consensus calling statistics for this batch
    stats: ConsensusCallingStats,
    /// Overlapping correction stats for this batch (if enabled)
    overlapping_stats: Option<CorrectionStats>,
    /// Family metrics collected from this batch
    family_metrics: CollectedDuplexFamilyMetrics,
}

impl MemoryEstimate for DuplexProcessedBatchWithMetrics {
    fn estimate_heap_size(&self) -> usize {
        self.consensus_reads
            .iter()
            .chain(self.rejects.iter())
            .map(MemoryEstimate::estimate_heap_size)
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;
    use fgumi_lib::bam_io::create_bam_writer;
    use fgumi_lib::sam::builder::{RecordBuilder, RecordPairBuilder};
    use noodles::sam;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReferenceSequence;
    use rstest::rstest;
    use std::collections::HashSet;
    use std::num::NonZeroUsize;
    use std::path::PathBuf;
    use tempfile::{NamedTempFile, TempDir};

    /// Helper struct for managing temporary test file paths.
    struct TestPaths {
        #[allow(dead_code)]
        dir: TempDir,
        pub output: PathBuf,
    }

    impl TestPaths {
        fn new() -> Result<Self> {
            let dir = TempDir::new()?;
            Ok(Self { output: dir.path().join("output.bam"), dir })
        }

        fn output_n(&self, n: usize) -> PathBuf {
            self.dir.path().join(format!("output{n}.bam"))
        }
    }

    /// Helper function to create a Duplex command with commonly used test defaults.
    /// Uses sensible test defaults: disabled per-base tags, disabled overlapping consensus.
    fn create_duplex_with_paths(input: PathBuf, output: PathBuf) -> Duplex {
        Duplex {
            io: BamIoOptions { input, output },
            rejects_opts: RejectsOptions::default(),
            stats_opts: StatsOptions::default(),
            read_group: ReadGroupOptions {
                read_name_prefix: Some("consensus".to_string()),
                read_group_id: "duplex".to_string(),
            },
            consensus: ConsensusCallingOptions {
                output_per_base_tags: false,
                min_consensus_base_quality: 0,
                ..ConsensusCallingOptions::default()
            },
            overlapping: OverlappingConsensusOptions::default(),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            min_reads: vec![1],
            max_reads_per_strand: None,
            cell_tag: None,
            scheduler_opts: SchedulerOptions::default(),
            metrics_output: None,
            description: None,
        }
    }

    // ========================================================================
    // Test Helper Functions
    // ========================================================================

    /// Create a test header with reference sequences
    fn create_test_header() -> sam::Header {
        use noodles::sam::header::record::value::map::Program;

        let builder = sam::Header::builder()
            .set_header(Map::default())
            .add_program("fgumi", Map::<Program>::default())
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(248_956_422).unwrap()),
            );

        builder.build()
    }

    /// Build a test read with common parameters
    #[allow(clippy::cast_sign_loss)]
    fn build_test_read(
        name: &str,
        ref_id: usize,
        pos: i32,
        mapq: u8,
        flags: u16,
        bases: &str,
    ) -> sam::alignment::RecordBuf {
        let cigar = format!("{}M", bases.len());
        let qual = vec![40u8; bases.len()];

        let mut builder = RecordBuilder::new()
            .name(name)
            .sequence(bases)
            .qualities(&qual)
            .reference_sequence_id(ref_id)
            .alignment_start(pos as usize)
            .mapping_quality(mapq)
            .cigar(&cigar);

        // Set flags based on raw u16
        let sam_flags = sam::alignment::record::Flags::from(flags);
        if sam_flags.is_segmented() {
            builder = builder.paired(true);
        }
        if sam_flags.is_properly_segmented() {
            builder = builder.properly_paired(true);
        }
        if sam_flags.is_first_segment() {
            builder = builder.first_segment(true);
        } else if sam_flags.is_last_segment() {
            builder = builder.first_segment(false);
        }
        if sam_flags.is_reverse_complemented() {
            builder = builder.reverse_complement(true);
        }
        if sam_flags.is_unmapped() {
            builder = builder.unmapped(true);
        }
        if sam_flags.is_mate_unmapped() {
            builder = builder.mate_unmapped(true);
        }

        builder.build()
    }

    /// Create a pair of reads with MI tag
    #[allow(clippy::cast_sign_loss, clippy::too_many_arguments)]
    fn build_duplex_pair(
        name: &str,
        ref_id: usize,
        pos1: i32,
        pos2: i32,
        mi_tag: &str,
        bases: &str,
        rx_tag: Option<&str>,
        cell_tag: Option<&str>,
    ) -> (sam::alignment::RecordBuf, sam::alignment::RecordBuf) {
        let mut builder = RecordPairBuilder::new()
            .name(name)
            .r1_sequence(bases)
            .r2_sequence(bases)
            .reference_sequence_id(ref_id)
            .r1_start(pos1 as usize)
            .r2_start(pos2 as usize)
            .tag("MI", mi_tag);

        if let Some(rx) = rx_tag {
            builder = builder.tag("RX", rx);
        }
        if let Some(cell) = cell_tag {
            builder = builder.tag("XX", cell);
        }

        builder.build()
    }

    /// Write records to a temporary BAM file
    fn create_test_bam(mut records: Vec<sam::alignment::RecordBuf>) -> Result<NamedTempFile> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        // Sort records by template coordinate as required by duplex consensus calling
        // Sort by: min(R1_pos, R2_pos), max(R1_pos, R2_pos), MI tag, read name
        records.sort_by_key(|r| {
            let pos = r.alignment_start().map_or(usize::MAX, usize::from);
            let mate_pos = r.mate_alignment_start().map_or(usize::MAX, usize::from);
            let min_pos = pos.min(mate_pos);
            let max_pos = pos.max(mate_pos);
            let mi_tag = noodles::sam::alignment::record::data::field::Tag::from([b'M', b'I']);
            let mi =
                if let Some(noodles::sam::alignment::record_buf::data::field::Value::String(s)) =
                    r.data().get(&mi_tag)
                {
                    String::from_utf8_lossy(&s.iter().copied().collect::<Vec<u8>>()).to_string()
                } else {
                    String::new()
                };
            let name = r.name().map(|n| n.to_vec()).unwrap_or_default();
            (min_pos, max_pos, mi, name)
        });

        let mut writer = create_bam_writer(temp_file.path(), &header, 1, 6)?;

        for record in records {
            writer.write_alignment_record(&header, &record)?;
        }

        drop(writer); // Ensure file is flushed
        Ok(temp_file)
    }

    /// Read all records from a BAM file
    fn read_bam_records(path: &std::path::Path) -> Result<Vec<sam::alignment::RecordBuf>> {
        let (mut reader, header) = create_bam_reader(path, 1)?;
        let mut records = Vec::new();

        for result in reader.records() {
            let record = result?;
            let record_buf =
                sam::alignment::RecordBuf::try_from_alignment_record(&header, &record)?;
            records.push(record_buf);
        }

        Ok(records)
    }

    /// Helper function to get a string tag from a record
    fn get_string_tag(record: &sam::alignment::RecordBuf, tag_name: &str) -> Option<String> {
        let tag_bytes = tag_name.as_bytes();
        let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

        record.data().get(&tag).and_then(|v| {
            if let Value::String(s) = v {
                Some(String::from_utf8_lossy(s).to_string())
            } else {
                None
            }
        })
    }

    /// Count unique MI tags (without /A or /B suffix)
    fn count_unique_mi_tags(records: &[sam::alignment::RecordBuf]) -> usize {
        let tags: HashSet<String> =
            records.iter().filter_map(|r| get_string_tag(r, "MI")).collect();
        tags.len()
    }

    // ========================================================================
    // Unit Tests for Duplex Configuration
    // ========================================================================

    #[test]
    fn test_default_parameters() {
        let duplex =
            create_duplex_with_paths(PathBuf::from("test.bam"), PathBuf::from("output.bam"));

        assert_eq!(duplex.min_reads, vec![1]);
        assert_eq!(duplex.consensus.min_input_base_quality, 10);
        assert!(duplex.threading.is_single_threaded());
        assert!(!duplex.consensus.trim);
    }

    #[test]
    fn test_custom_min_reads_single_value() {
        let mut duplex =
            create_duplex_with_paths(PathBuf::from("test.bam"), PathBuf::from("output.bam"));
        duplex.min_reads = vec![5];

        assert_eq!(duplex.min_reads, vec![5]);
    }

    #[test]
    fn test_custom_min_reads_three_values() {
        let mut duplex =
            create_duplex_with_paths(PathBuf::from("test.bam"), PathBuf::from("output.bam"));
        duplex.min_reads = vec![10, 5, 3];

        assert_eq!(duplex.min_reads, vec![10, 5, 3]);
    }

    #[test]
    fn test_output_per_base_tags_enabled() {
        let mut duplex =
            create_duplex_with_paths(PathBuf::from("test.bam"), PathBuf::from("output.bam"));
        duplex.consensus.output_per_base_tags = true;

        assert!(duplex.consensus.output_per_base_tags);
    }

    #[test]
    fn test_multithreaded_configuration() {
        let mut duplex =
            create_duplex_with_paths(PathBuf::from("test.bam"), PathBuf::from("output.bam"));
        duplex.threading = ThreadingOptions::new(8);

        assert_eq!(duplex.threading.threads, Some(8));
    }

    #[test]
    fn test_trim_and_downsample_options() {
        let mut duplex =
            create_duplex_with_paths(PathBuf::from("test.bam"), PathBuf::from("output.bam"));
        duplex.consensus.trim = true;
        duplex.max_reads_per_strand = Some(100);

        assert!(duplex.consensus.trim);
        assert_eq!(duplex.max_reads_per_strand, Some(100));
    }

    // ========================================================================
    // Integration Tests for Duplex Consensus Calling
    // ========================================================================

    #[test]
    fn test_duplex_consensus_basic_ab_ba_pairing() -> Result<()> {
        // Test basic duplex consensus from matching AB and BA strand groups
        // Following Scala test pattern: 1 A pair and 1 B pair
        let mut records = Vec::new();

        // A reads: R1 at 100 (forward), R2 at 200 (reverse) - query name "q1"
        let (r1_a, r2_a) =
            build_duplex_pair("q1", 0, 100, 200, "1/A", "AAAAAAAAAA", Some("AAT-CCG"), None);
        records.push(r1_a);
        records.push(r2_a);

        // B reads: R1 at 200 (reverse), R2 at 100 (forward) - query name "q2", positions swapped
        let flags1 = 0x53; // paired, proper pair, first in pair, reverse
        let flags2 = 0x83; // paired, proper pair, second in pair, forward
        let mut r1_b = build_test_read("q2", 0, 200, 60, flags1, "TTTTTTTTTT"); // Reverse complement of A
        let mut r2_b = build_test_read("q2", 0, 100, 60, flags2, "AAAAAAAAAA");

        // Set mate info
        *r1_b.mate_reference_sequence_id_mut() = Some(0);
        *r1_b.mate_alignment_start_mut() = noodles::core::Position::try_from(100_usize).ok();
        *r2_b.mate_reference_sequence_id_mut() = Some(0);
        *r2_b.mate_alignment_start_mut() = noodles::core::Position::try_from(200_usize).ok();

        let mi = Tag::from([b'M', b'I']);
        r1_b.data_mut().insert(mi, Value::String(b"1/B".into()));
        r2_b.data_mut().insert(mi, Value::String(b"1/B".into()));

        let rx = Tag::from([b'R', b'X']);
        r1_b.data_mut().insert(rx, Value::String(b"CCG-AAT".into()));
        r2_b.data_mut().insert(rx, Value::String(b"CCG-AAT".into()));

        records.push(r1_b);
        records.push(r2_b);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = create_duplex_with_paths(input.path().to_path_buf(), paths.output.clone());

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should produce 2 consensus reads (R1 and R2 of the duplex consensus)
        assert_eq!(output_records.len(), 2, "Should have 2 duplex consensus reads");

        // Check that MI tags no longer have /A or /B suffix
        for record in &output_records {
            let mi = get_string_tag(record, "MI");
            assert!(mi.is_some(), "MI tag should be present");
            let mi = mi.unwrap();
            assert_eq!(mi, "1", "MI tag should be '1' without /A or /B suffix");
        }

        // Check that RX tags are preserved
        for record in &output_records {
            let rx = get_string_tag(record, "RX");
            assert!(rx.is_some(), "RX tag should be preserved");
        }

        Ok(())
    }

    #[test]
    fn test_duplex_no_consensus_when_strands_mismatched_r1() -> Result<()> {
        // Test that no consensus is generated when AB-R1s and BA-R2s are on different strands
        // (not properly complementary)
        let mut records = Vec::new();

        // AB reads: R1 forward (+), R2 forward (+) - WRONG!
        let (r1, mut r2) = build_duplex_pair("ab1", 0, 100, 200, "1/A", "AAAAAAAAAA", None, None);
        // Make R2 forward instead of reverse
        let mut flags2 = r2.flags();
        flags2.set(noodles::sam::alignment::record::Flags::REVERSE_COMPLEMENTED, false);
        *r2.flags_mut() = flags2;
        records.push(r1);
        records.push(r2);

        // BA reads: R1 forward (+), R2 reverse (-) - correct orientation
        let flags1 = 0x53; // paired, first in pair, reverse
        let flags2 = 0x83; // paired, second in pair
        let mut r1 = build_test_read("ba1", 0, 200, 60, flags1, "AAAAAAAAAA");
        let mut r2 = build_test_read("ba1", 0, 100, 60, flags2, "AAAAAAAAAA");

        *r1.mate_reference_sequence_id_mut() = Some(0);
        *r1.mate_alignment_start_mut() = noodles::core::Position::try_from(100_usize).ok();
        *r2.mate_reference_sequence_id_mut() = Some(0);
        *r2.mate_alignment_start_mut() = noodles::core::Position::try_from(200_usize).ok();

        let mi = Tag::from([b'M', b'I']);
        r1.data_mut().insert(mi, Value::String(b"1/B".into()));
        r2.data_mut().insert(mi, Value::String(b"1/B".into()));

        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = create_duplex_with_paths(input.path().to_path_buf(), paths.output.clone());

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should produce 0 consensus reads due to strand mismatch
        assert_eq!(output_records.len(), 0, "Should have 0 consensus reads due to strand mismatch");

        Ok(())
    }

    #[test]
    fn test_duplex_with_min_reads_filtering() -> Result<()> {
        // Test that min-reads filtering works for AB and BA strands independently
        let mut records = Vec::new();

        // AB reads: 5 pairs (should pass min-reads=3)
        for i in 1..=5 {
            let (r1, r2) = build_duplex_pair(
                &format!("ab{i}"),
                0,
                100,
                100,
                "1/A",
                "AAAAAAAAAA",
                Some("AAT-CCG"),
                None,
            );
            records.push(r1);
            records.push(r2);
        }

        // BA reads: 2 pairs (should fail min-reads=3)
        for i in 1..=2 {
            let flags1 = 0x53;
            let flags2 = 0x83;
            let mut r1 = build_test_read(&format!("ba{i}"), 0, 100, 60, flags1, "AAAAAAAAAA");
            let mut r2 = build_test_read(&format!("ba{i}"), 0, 100, 60, flags2, "AAAAAAAAAA");

            *r1.mate_reference_sequence_id_mut() = Some(0);
            *r1.mate_alignment_start_mut() = noodles::core::Position::try_from(100_usize).ok();
            *r2.mate_reference_sequence_id_mut() = Some(0);
            *r2.mate_alignment_start_mut() = noodles::core::Position::try_from(100_usize).ok();

            let mi = Tag::from([b'M', b'I']);
            r1.data_mut().insert(mi, Value::String(b"1/B".into()));
            r2.data_mut().insert(mi, Value::String(b"1/B".into()));

            let rx = Tag::from([b'R', b'X']);
            r1.data_mut().insert(rx, Value::String(b"CCG-AAT".into()));
            r2.data_mut().insert(rx, Value::String(b"CCG-AAT".into()));

            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let mut cmd = create_duplex_with_paths(input.path().to_path_buf(), paths.output.clone());
        cmd.min_reads = vec![3, 3, 3]; // duplex=3, AB=3, BA=3

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // BA strand has only 2 reads (< 3), so no duplex consensus should be generated
        assert_eq!(
            output_records.len(),
            0,
            "Should have 0 consensus reads because BA strand has insufficient reads"
        );

        Ok(())
    }

    #[test]
    fn test_duplex_with_per_base_tags() -> Result<()> {
        // Test that per-base tags are generated when enabled
        let mut records = Vec::new();

        // AB reads
        for i in 1..=3 {
            let (r1, r2) = build_duplex_pair(
                &format!("ab{i}"),
                0,
                100,
                100,
                "1/A",
                "AAAAAAAAAA",
                Some("AAT-CCG"),
                None,
            );
            records.push(r1);
            records.push(r2);
        }

        // BA reads
        for i in 1..=3 {
            let flags1 = 0x53;
            let flags2 = 0x83;
            let mut r1 = build_test_read(&format!("ba{i}"), 0, 100, 60, flags1, "AAAAAAAAAA");
            let mut r2 = build_test_read(&format!("ba{i}"), 0, 100, 60, flags2, "AAAAAAAAAA");

            *r1.mate_reference_sequence_id_mut() = Some(0);
            *r1.mate_alignment_start_mut() = noodles::core::Position::try_from(100_usize).ok();
            *r2.mate_reference_sequence_id_mut() = Some(0);
            *r2.mate_alignment_start_mut() = noodles::core::Position::try_from(100_usize).ok();

            let mi = Tag::from([b'M', b'I']);
            r1.data_mut().insert(mi, Value::String(b"1/B".into()));
            r2.data_mut().insert(mi, Value::String(b"1/B".into()));

            let rx = Tag::from([b'R', b'X']);
            r1.data_mut().insert(rx, Value::String(b"CCG-AAT".into()));
            r2.data_mut().insert(rx, Value::String(b"CCG-AAT".into()));

            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let mut cmd = create_duplex_with_paths(input.path().to_path_buf(), paths.output.clone());
        cmd.consensus.output_per_base_tags = true; // Enable per-base tags

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 2, "Should have 2 consensus reads");

        // Check for per-base tags (ad, bd for depths)
        for record in &output_records {
            // At least one of the per-base depth tags should be present
            let has_ad = record.data().get(&Tag::from([b'a', b'd'])).is_some();
            let has_bd = record.data().get(&Tag::from([b'b', b'd'])).is_some();
            assert!(has_ad || has_bd, "Per-base tags should be present when enabled");
        }

        Ok(())
    }

    #[test]
    fn test_duplex_with_cell_barcode_preservation() -> Result<()> {
        // Test that cell barcodes are preserved in duplex consensus reads
        let mut records = Vec::new();

        // AB reads with cell barcode
        for i in 1..=3 {
            let (r1, r2) = build_duplex_pair(
                &format!("ab{i}"),
                0,
                100,
                100,
                "1/A",
                "AAAAAAAAAA",
                Some("AAT-CCG"),
                Some("CELLBC"),
            );
            records.push(r1);
            records.push(r2);
        }

        // BA reads with same cell barcode
        for i in 1..=3 {
            let flags1 = 0x53;
            let flags2 = 0x83;
            let mut r1 = build_test_read(&format!("ba{i}"), 0, 100, 60, flags1, "AAAAAAAAAA");
            let mut r2 = build_test_read(&format!("ba{i}"), 0, 100, 60, flags2, "AAAAAAAAAA");

            *r1.mate_reference_sequence_id_mut() = Some(0);
            *r1.mate_alignment_start_mut() = noodles::core::Position::try_from(100_usize).ok();
            *r2.mate_reference_sequence_id_mut() = Some(0);
            *r2.mate_alignment_start_mut() = noodles::core::Position::try_from(100_usize).ok();

            let mi = Tag::from([b'M', b'I']);
            r1.data_mut().insert(mi, Value::String(b"1/B".into()));
            r2.data_mut().insert(mi, Value::String(b"1/B".into()));

            let rx = Tag::from([b'R', b'X']);
            r1.data_mut().insert(rx, Value::String(b"CCG-AAT".into()));
            r2.data_mut().insert(rx, Value::String(b"CCG-AAT".into()));

            let xx = Tag::from([b'X', b'X']);
            r1.data_mut().insert(xx, Value::String(b"CELLBC".into()));
            r2.data_mut().insert(xx, Value::String(b"CELLBC".into()));

            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let mut cmd = create_duplex_with_paths(input.path().to_path_buf(), paths.output.clone());
        cmd.cell_tag = Some("XX".to_string()); // Enable cell barcode tracking

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 2, "Should have 2 consensus reads");

        // Check that cell barcode is preserved
        for record in &output_records {
            let cell_bc = get_string_tag(record, "XX");
            assert_eq!(cell_bc, Some("CELLBC".to_string()), "Cell barcode should be preserved");
        }

        Ok(())
    }

    #[test]
    fn test_duplex_multithreading_produces_same_results() -> Result<()> {
        // Test that single-threaded and multi-threaded processing produce the same results
        let mut records = Vec::new();

        // Create multiple duplex groups
        for group in 1..=3 {
            // AB reads
            for i in 1..=3 {
                let (r1, r2) = build_duplex_pair(
                    &format!("g{group}_ab{i}"),
                    0,
                    100 * group,
                    100 * group,
                    &format!("{group}/A"),
                    "AAAAAAAAAA",
                    Some("AAT-CCG"),
                    None,
                );
                records.push(r1);
                records.push(r2);
            }

            // BA reads
            for i in 1..=3 {
                let flags1 = 0x53;
                let flags2 = 0x83;
                let mut r1 = build_test_read(
                    &format!("g{group}_ba{i}"),
                    0,
                    100 * group,
                    60,
                    flags1,
                    "AAAAAAAAAA",
                );
                let mut r2 = build_test_read(
                    &format!("g{group}_ba{i}"),
                    0,
                    100 * group,
                    60,
                    flags2,
                    "AAAAAAAAAA",
                );

                *r1.mate_reference_sequence_id_mut() = Some(0);
                *r1.mate_alignment_start_mut() =
                    noodles::core::Position::try_from((100 * group) as usize).ok();
                *r2.mate_reference_sequence_id_mut() = Some(0);
                *r2.mate_alignment_start_mut() =
                    noodles::core::Position::try_from((100 * group) as usize).ok();

                let mi = Tag::from([b'M', b'I']);
                r1.data_mut().insert(mi, Value::String(format!("{group}/B").as_bytes().into()));
                r2.data_mut().insert(mi, Value::String(format!("{group}/B").as_bytes().into()));

                let rx = Tag::from([b'R', b'X']);
                r1.data_mut().insert(rx, Value::String(b"CCG-AAT".into()));
                r2.data_mut().insert(rx, Value::String(b"CCG-AAT".into()));

                records.push(r1);
                records.push(r2);
            }
        }

        // Test with 1 thread
        let input1 = create_test_bam(records.clone())?;
        let paths = TestPaths::new()?;

        let cmd1 = create_duplex_with_paths(input1.path().to_path_buf(), paths.output.clone());

        cmd1.execute("test")?;
        let output1_records = read_bam_records(&paths.output)?;

        // Test with 4 threads
        let input2 = create_test_bam(records)?;
        let output2_path = paths.output_n(2);

        let mut cmd2 = create_duplex_with_paths(input2.path().to_path_buf(), output2_path.clone());
        cmd2.threading = ThreadingOptions::new(4);

        cmd2.execute("test")?;
        let output2_records = read_bam_records(&output2_path)?;

        // Should produce the same number of consensus reads
        assert_eq!(
            output1_records.len(),
            output2_records.len(),
            "Single-threaded and multi-threaded should produce same number of reads"
        );

        // Should have 3 groups * 2 reads per group = 6 reads
        assert_eq!(output1_records.len(), 6, "Should have 6 consensus reads (3 duplex pairs)");

        // Check that we have 3 unique MI tags in both outputs
        assert_eq!(count_unique_mi_tags(&output1_records), 3, "Should have 3 unique MI tags");
        assert_eq!(count_unique_mi_tags(&output2_records), 3, "Should have 3 unique MI tags");

        Ok(())
    }

    #[test]
    fn test_duplex_only_one_strand_no_consensus() -> Result<()> {
        // Test that no duplex consensus is generated when only one strand (A or B) is present
        let mut records = Vec::new();

        // Only AB reads, no BA reads
        for i in 1..=5 {
            let (r1, r2) = build_duplex_pair(
                &format!("ab{i}"),
                0,
                100,
                100,
                "1/A",
                "AAAAAAAAAA",
                Some("AAT-CCG"),
                None,
            );
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = create_duplex_with_paths(input.path().to_path_buf(), paths.output.clone());

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should produce 0 consensus reads because we need both A and B strands
        assert_eq!(
            output_records.len(),
            0,
            "Should have 0 consensus reads when only one strand is present"
        );

        Ok(())
    }

    /// Parameterized test for all threading modes.
    ///
    /// Tests:
    /// - `None`: Single-threaded fast path, no pipeline
    /// - `Some(1)`: Pipeline with 1 thread
    /// - `Some(2)`: Pipeline with 2 threads
    #[rstest]
    #[case::fast_path(ThreadingOptions::none())]
    #[case::pipeline_1(ThreadingOptions::new(1))]
    #[case::pipeline_2(ThreadingOptions::new(2))]
    fn test_threading_modes(#[case] threading: ThreadingOptions) -> Result<()> {
        // Create test data with A and B strand reads
        let mut records = Vec::new();

        // Create a group with A strand (MI=1/A) and B strand (MI=1/B) for duplex calling
        for i in 0..5 {
            let (r1, r2) = build_duplex_pair(
                &format!("a{i}"),
                0,
                100,
                200,
                "1/A",
                "AAAAAAAAAA",
                Some("AAT-CCG"),
                None,
            );
            records.push(r1);
            records.push(r2);
        }
        for i in 0..5 {
            let (r1, r2) = build_duplex_pair(
                &format!("b{i}"),
                0,
                100,
                200,
                "1/B",
                "AAAAAAAAAA",
                Some("AAT-CCG"),
                None,
            );
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let mut cmd = create_duplex_with_paths(input.path().to_path_buf(), paths.output.clone());
        cmd.threading = threading;
        cmd.execute("test")?;

        // Should have at least some output (duplex consensus)
        assert!(&paths.output.exists());

        Ok(())
    }
}
