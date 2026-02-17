//! # CODEC Consensus Calling Command
//!
//! Implementation of the `codec` command for calling consensus from CODEC sequencing data.
//!
//! CODEC (Bae et al 2023) is a sequencing protocol where each read-pair sequences both strands
//! of the original duplex molecule. R1 comes from one strand, R2 from the opposite strand,
//! allowing even a single read-pair to generate duplex consensus.

use crate::commands::command::Command;
use crate::commands::consensus_runner::{ConsensusStatsOps, create_unmapped_consensus_header};
use anyhow::{Context, Result, bail};
use clap::Parser;
use crossbeam_queue::SegQueue;
use fgoxide::io::DelimFile;
use fgumi_lib::bam_io::{
    create_bam_reader_for_pipeline, create_bam_writer, create_optional_bam_writer,
    create_raw_bam_reader,
};

use super::common::{
    BamIoOptions, CompressionOptions, ConsensusCallingOptions, QueueMemoryOptions,
    ReadGroupOptions, RejectsOptions, SchedulerOptions, StatsOptions, ThreadingOptions,
    build_pipeline_config,
};
use fgumi_lib::consensus::codec_caller::{
    CodecConsensusCaller, CodecConsensusOptions, CodecConsensusStats,
};
use fgumi_lib::consensus_caller::{ConsensusCaller, ConsensusOutput};
use fgumi_lib::logging::{OperationTimer, log_consensus_summary};
use fgumi_lib::mi_group::{RawMiGroup, RawMiGroupBatch, RawMiGroupIterator, RawMiGrouper};
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::read_info::LibraryIndex;
use fgumi_lib::unified_pipeline::{
    GroupKeyConfig, Grouper, MemoryEstimate, run_bam_pipeline_from_reader,
};
use fgumi_lib::vendored::RawRecord;
// RejectionTracker now used via ConsensusStatsOps trait in consensus_runner
use fgumi_lib::validation::{optional_string_to_tag, validate_file_exists};
use log::info;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::io::{self, Write as IoWrite};
use std::sync::Arc;

// ============================================================================
// Types for 7-step pipeline processing
// ============================================================================

/// Result from processing a batch of MI groups through CODEC consensus calling.
struct CodecProcessedBatch {
    /// Consensus reads to write to output BAM
    consensus_output: ConsensusOutput,
    /// Rejected reads as raw BAM bytes (written to rejects file if enabled)
    rejects: Vec<Vec<u8>>,
    /// Number of MI groups in this batch
    groups_count: u64,
    /// CODEC consensus calling statistics for this batch
    stats: CodecConsensusStats,
}

impl MemoryEstimate for CodecProcessedBatch {
    fn estimate_heap_size(&self) -> usize {
        self.consensus_output.data.capacity()
            + self.rejects.iter().map(Vec::capacity).sum::<usize>()
            + self.rejects.capacity() * std::mem::size_of::<Vec<u8>>()
    }
}

/// Metrics collected from each batch during parallel processing.
#[derive(Default)]
struct CollectedCodecMetrics {
    /// CODEC consensus calling statistics
    stats: CodecConsensusStats,
    /// Number of MI groups processed
    groups_processed: u64,
    /// Rejected reads as raw BAM bytes for deferred writing
    rejects: Vec<Vec<u8>>,
}

/// Call CODEC consensus reads from template-coordinate sorted BAM
///
/// CODEC is a sequencing protocol where a single read-pair sequences both strands of the
/// original duplex molecule. R1 sequences one strand, R2 sequences the opposite strand,
/// enabling duplex consensus calling from individual read-pairs.
///
/// The algorithm:
/// 1. Groups reads by molecule ID (MI tag)
/// 2. For each molecule:
///    - Clips read pairs where they extend past mate ends
///    - Filters R1s and R2s for compatible alignments
///    - Checks for sufficient overlap
///    - Calls single-strand consensus from R1s and R2s separately
///    - Combines into duplex consensus
///    - Applies quality masking
/// 3. Outputs unmapped consensus fragments (not paired-end)
///
/// ## Input Requirements
///
/// - BAM must be template-coordinate sorted (or queryname sorted)
/// - Reads must be aligned (mapped)
/// - Required tags:
///   - `RX`: Raw UMI bases
///   - `MI`: Molecule ID (from `group`)
///   - `CB`: Cell barcode (optional, if using --cell-tag)
///
/// ## Grouping Strategy
///
/// Must use `adjacency` or `identity` grouping (NOT `paired`).
/// MI tags should not have /A or /B suffixes.
///
/// ## Output
///
/// - Output BAM contains unmapped consensus fragments (not paired-end)
/// - Each consensus represents a full duplex molecule
/// - Includes rich per-base and per-read tags
///
/// Consensus reads have a number of additional optional tags set in the resulting BAM file.
/// The tag names follow a pattern where the first letter (a, b or c) denotes that the tag
/// applies to the first single strand consensus (a), second single-strand consensus (b) or
/// the final duplex consensus (c). The second letter captures the meaning of the tag
/// (e.g. d=depth, m=min depth, e=errors/error-rate) and is upper case for values that are
/// one per read and lower case for values that are one per base.
///
/// The tags break down into those that are single-valued per read:
///
///   consensus depth      \[aD,bD,cD\] (int)  : the maximum depth of raw reads
///   consensus min depth  \[aM,bM,cM\] (int)  : the minimum depth of raw reads
///   consensus error rate \[aE,bE,cE\] (float): the fraction of bases disagreeing with consensus
///
/// And those that have a value per base:
///
///   consensus depth  \[ad,bd\] (short[]): the count of bases contributing to consensus
///   consensus errors \[ae,be\] (short[]): the count of disagreeing bases
///   consensus bases  \[ac,bc\] (string) : the single-strand consensus bases
///   consensus quals  \[aq,bq\] (string) : the single-strand consensus qualities
///
/// ## Examples
///
/// Basic usage:
/// ```bash
/// fgumi codec -i grouped.bam -o codec_consensus.bam
/// ```
///
/// With quality thresholds and statistics:
/// ```bash
/// fgumi codec \
///   -i grouped.bam \
///   -o codec_consensus.bam \
///   -r rejects.bam \
///   -s stats.txt \
///   -m 15 \
///   -M 3 \
///   -d 10
/// ```
#[derive(Parser, Debug)]
#[command(
    name = "codec",
    about = "\x1b[38;5;180m[CONSENSUS]\x1b[0m      \x1b[36mCall CODEC consensus reads from grouped BAM\x1b[0m",
    long_about = None
)]
pub struct Codec {
    /// Input/output BAM options
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Options for writing rejected reads
    #[command(flatten)]
    pub rejects_opts: RejectsOptions,

    /// Options for writing statistics
    #[command(flatten)]
    pub stats_opts: StatsOptions,

    /// Read group and name prefix options
    #[command(flatten)]
    pub read_group: ReadGroupOptions,

    /// Consensus calling options
    #[command(flatten)]
    pub consensus: ConsensusCallingOptions,

    /// Threading options for parallel processing
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output
    #[command(flatten)]
    pub compression: CompressionOptions,

    // --- CODEC-specific options below ---
    /// Minimum read pairs per strand to form consensus (same as --min-reads)
    #[arg(short = 'M', long = "min-reads", default_value = "1")]
    pub min_reads: usize,

    /// Maximum read pairs per strand (downsample if exceeded)
    #[arg(long = "max-reads")]
    pub max_reads: Option<usize>,

    /// Minimum duplex overlap length in bases
    #[arg(short = 'd', long = "min-duplex-length", default_value = "1")]
    pub min_duplex_length: usize,

    /// Reduce single-strand region quality to this value (0-93).
    /// Note: This uses a different short flag than duplex's -q for min-base-quality.
    #[arg(long = "single-strand-qual")]
    pub single_strand_qual: Option<u8>,

    /// Reduce outer bases quality to this value (0-93)
    #[arg(short = 'Q', long = "outer-bases-qual")]
    pub outer_bases_qual: Option<u8>,

    /// Number of outer bases to reduce quality for
    #[arg(short = 'O', long = "outer-bases-length", default_value = "5")]
    pub outer_bases_length: usize,

    /// Maximum duplex disagreement rate (0.0-1.0)
    #[arg(short = 'x', long = "max-duplex-disagreement-rate", default_value = "1.0")]
    pub max_duplex_disagreement_rate: f64,

    /// Maximum number of duplex disagreements
    #[arg(short = 'X', long = "max-duplex-disagreements")]
    pub max_duplex_disagreements: Option<usize>,

    /// SAM tag containing the cell barcode
    #[arg(short = 'c', long = "cell-tag", default_value = "CB")]
    pub cell_tag: Option<String>,

    /// Output sort order (Unsorted, Queryname, Coordinate, Unknown)
    #[arg(short = 'S', long = "sort-order")]
    pub sort_order: Option<String>,

    /// Scheduler and pipeline statistics options.
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

impl Command for Codec {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate inputs
        self.validate()?;
        validate_file_exists(&self.io.input, "input BAM file")?;

        let timer = OperationTimer::new("Calling CODEC consensus");

        // Get threading configuration (codec is balanced workload)
        let reader_threads = self.threading.num_threads();
        let worker_threads = self.threading.num_threads();
        let writer_threads = self.threading.num_threads();

        info!("Starting CODEC consensus calling");
        info!("Input: {}", self.io.input.display());
        info!("Output: {}", self.io.output.display());
        info!("Min reads: {}", self.min_reads);
        if let Some(max) = self.max_reads {
            info!("Max reads: {max}");
        }
        info!("Error rate pre-UMI: Q{}", self.consensus.error_rate_pre_umi);
        info!("Error rate post-UMI: Q{}", self.consensus.error_rate_post_umi);
        info!("Min duplex length: {}", self.min_duplex_length);
        info!("Worker threads: {worker_threads}");
        info!("Reader threads: {reader_threads}");
        if self.consensus.trim {
            info!("Quality trimming enabled");
        }
        // Note: Unlike simplex/duplex, CODEC does not support overlapping consensus calling
        // (matching fgbio's CallCodecConsensusReads which has no such option).

        // Open input BAM using streaming-capable reader for pipeline use
        let (reader, header) = create_bam_reader_for_pipeline(&self.io.input)?;

        // Create output header for unmapped consensus reads
        let output_header = create_unmapped_consensus_header(
            &header,
            &self.read_group.read_group_id,
            "Read group",
            command_line,
        )?;

        // Use library name from header if no prefix is specified (like fgbio)
        let read_name_prefix = self.read_group.prefix_or_from_header(&header);

        // Parse cell tag
        let cell_tag = optional_string_to_tag(self.cell_tag.as_deref(), "cell-tag")?;

        // Enable rejects tracking if rejects file is specified
        let track_rejects = self.rejects_opts.is_enabled();

        // Process reads using streaming by MI groups
        info!("Processing reads and calling consensus (streaming)...");

        // ============================================================
        // --threads N mode: Use 7-step unified pipeline
        // None: Use single-threaded fast path
        // ============================================================
        // IMPORTANT: Check this BEFORE creating any output file handles to avoid
        // file handle conflicts that can corrupt the output.
        if let Some(threads) = self.threading.threads {
            let result = self.execute_threads_mode(
                threads,
                reader,
                header.clone(),
                output_header,
                read_name_prefix.clone(),
                track_rejects,
            );
            timer.log_completion(0); // Completion logged in execute_threads_mode
            return result;
        }

        // Drop the reader for single-threaded mode - we use MiGroupIterator which needs its own reader
        drop(reader);

        // ============================================================
        // For non-pipeline modes, create output writers here
        // ============================================================

        // Open output BAM writer with multi-threaded BGZF compression
        let mut writer = create_bam_writer(
            &self.io.output,
            &output_header,
            writer_threads,
            self.compression.compression_level,
        )?;

        // Open rejects writer if rejects file is specified
        let mut rejects_writer = create_optional_bam_writer(
            self.rejects_opts.rejects.as_ref(),
            &header,
            writer_threads,
            self.compression.compression_level,
        )?;

        // Create options
        let options = CodecConsensusOptions {
            min_input_base_quality: self.consensus.min_input_base_quality,
            error_rate_pre_umi: self.consensus.error_rate_pre_umi,
            error_rate_post_umi: self.consensus.error_rate_post_umi,
            min_reads_per_strand: self.min_reads,
            max_reads_per_strand: self.max_reads,
            min_duplex_length: self.min_duplex_length,
            single_strand_qual: self.single_strand_qual,
            outer_bases_qual: self.outer_bases_qual,
            outer_bases_length: self.outer_bases_length,
            max_duplex_disagreements: self.max_duplex_disagreements.unwrap_or(usize::MAX),
            max_duplex_disagreement_rate: self.max_duplex_disagreement_rate,
            cell_tag,
            produce_per_base_tags: self.consensus.output_per_base_tags,
            trim: self.consensus.trim,
            min_consensus_base_quality: self.consensus.min_consensus_base_quality,
        };

        // Note: CODEC does not support overlapping consensus (matching fgbio)
        // We keep the infrastructure in place but it's always disabled.

        // Track progress (count records written, not UMI groups)
        let mut record_count: usize = 0;
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        // Create the MI group iterator for single-threaded streaming (raw bytes)
        let (mut raw_reader, _) = create_raw_bam_reader(&self.io.input, 1)?;
        let raw_record_iter = std::iter::from_fn(move || {
            let mut record = RawRecord::new();
            match raw_reader.read_record(&mut record) {
                Ok(0) => None, // EOF
                Ok(_) => Some(Ok(record.into_inner())),
                Err(e) => Some(Err(e.into())),
            }
        });
        let mi_group_iter = RawMiGroupIterator::new(raw_record_iter, "MI")
            .with_cell_tag(cell_tag.map(|ct| *ct.as_ref()));

        let mut caller = CodecConsensusCaller::new_with_rejects_tracking(
            read_name_prefix,
            self.read_group.read_group_id.clone(),
            options,
            track_rejects,
        );

        for result in mi_group_iter {
            let (umi, records) = result.context("Failed to read MI group")?;

            // Call consensus directly — records are already raw bytes!
            let result: anyhow::Result<ConsensusOutput> = caller.consensus_reads(records);
            match result {
                Ok(output) => {
                    let batch_size = output.count;
                    record_count += batch_size;
                    writer
                        .get_mut()
                        .write_all(&output.data)
                        .context("Failed to write consensus read")?;
                    progress.log_if_needed(batch_size as u64);
                }
                Err(e) => {
                    if !e.to_string().contains("duplex disagreement") {
                        return Err(e)
                            .with_context(|| format!("Failed to call consensus for UMI: {umi}"));
                    }
                }
            }

            // Write rejected reads if tracking is enabled (already raw BAM bytes)
            if let Some(ref mut rw) = rejects_writer {
                for raw_record in caller.rejected_reads() {
                    let block_size = raw_record.len() as u32;
                    rw.get_mut()
                        .write_all(&block_size.to_le_bytes())
                        .context("Failed to write rejected read block size")?;
                    rw.get_mut().write_all(raw_record).context("Failed to write rejected read")?;
                }
                caller.clear_rejected_reads();
            }
        }

        // For single-threaded, use the caller's stats
        let merged_stats = caller.statistics().clone();

        progress.log_final();

        // Finish the buffered writer (flush remaining records and wait for writer thread)
        writer.into_inner().finish().context("Failed to finish output BAM")?;

        // Log statistics and write to file
        info!("Consensus calling complete");
        info!("Total records processed: {record_count}");

        let metrics = merged_stats.to_metrics();
        let consensus_count = metrics.consensus_reads;
        log_consensus_summary(&metrics);

        if let Some(stats_path) = &self.stats_opts.stats {
            DelimFile::default()
                .write_tsv(stats_path, [metrics])
                .with_context(|| format!("Failed to write statistics: {}", stats_path.display()))?;
            info!("Wrote statistics to: {}", stats_path.display());
        }

        timer.log_completion(consensus_count);

        // Close rejects writer if it was opened
        if let Some(mut rw) = rejects_writer {
            rw.finish(&header).context("Failed to finish rejects file")?;
            info!("Rejected reads written successfully");
        }

        Ok(())
    }
}
impl Codec {
    /// Validates command-line arguments
    fn validate(&self) -> Result<()> {
        // Validate error rates
        if self.consensus.error_rate_pre_umi == 0 {
            bail!("error-rate-pre-umi must be > 0");
        }
        if self.consensus.error_rate_post_umi == 0 {
            bail!("error-rate-post-umi must be > 0");
        }

        // Validate min/max reads
        if self.min_reads == 0 {
            bail!("min-reads must be >= 1");
        }
        if let Some(max) = self.max_reads {
            if max < self.min_reads {
                bail!("max-reads ({}) must be >= min-reads ({})", max, self.min_reads);
            }
        }

        // Validate duplex length
        if self.min_duplex_length == 0 {
            bail!("min-duplex-length must be >= 1");
        }

        // Validate disagreement rate
        if self.max_duplex_disagreement_rate < 0.0 || self.max_duplex_disagreement_rate > 1.0 {
            bail!("max-duplex-disagreement-rate must be between 0.0 and 1.0");
        }

        // Note: cell_tag validation is done in execute() via optional_string_to_tag()

        Ok(())
    }

    /// Execute using 7-step unified pipeline with --threads.
    ///
    /// This method is called when `--threads N` is specified with N > 1.
    /// It uses the lock-free 7-step unified pipeline for maximum performance.
    fn execute_threads_mode(
        &self,
        num_threads: usize,
        reader: Box<dyn std::io::Read + Send>,
        input_header: Header,
        output_header: Header,
        read_name_prefix: String,
        track_rejects: bool,
    ) -> Result<()> {
        // Configure pipeline
        let mut pipeline_config = build_pipeline_config(
            &self.scheduler_opts,
            &self.compression,
            &self.queue_memory,
            num_threads,
        )?;

        // Lock-free metrics collection
        let collected_metrics: Arc<SegQueue<CollectedCodecMetrics>> = Arc::new(SegQueue::new());
        let collected_metrics_for_serialize = Arc::clone(&collected_metrics);

        // Parse cell tag
        let cell_tag = optional_string_to_tag(self.cell_tag.as_deref(), "cell-tag")?;

        // Create options for CODEC consensus caller
        let options = CodecConsensusOptions {
            min_input_base_quality: self.consensus.min_input_base_quality,
            error_rate_pre_umi: self.consensus.error_rate_pre_umi,
            error_rate_post_umi: self.consensus.error_rate_post_umi,
            min_reads_per_strand: self.min_reads,
            max_reads_per_strand: self.max_reads,
            min_duplex_length: self.min_duplex_length,
            single_strand_qual: self.single_strand_qual,
            outer_bases_qual: self.outer_bases_qual,
            outer_bases_length: self.outer_bases_length,
            max_duplex_disagreements: self.max_duplex_disagreements.unwrap_or(usize::MAX),
            max_duplex_disagreement_rate: self.max_duplex_disagreement_rate,
            cell_tag,
            produce_per_base_tags: self.consensus.output_per_base_tags,
            trim: self.consensus.trim,
            min_consensus_base_quality: self.consensus.min_consensus_base_quality,
        };

        // Capture configuration for closures
        let read_group_id = self.read_group.read_group_id.clone();

        // Use larger batch size for codec (less work per group than simplex)
        let batch_size = 1000;

        // Clone input_header before pipeline (needed for rejects writing)
        let rejects_header = input_header.clone();

        // Enable raw-byte mode: skip noodles decode/encode for CPU savings
        let library_index = LibraryIndex::from_header(&input_header);
        pipeline_config.group_key_config = Some(if let Some(ct) = cell_tag {
            GroupKeyConfig::new_raw(library_index, ct)
        } else {
            GroupKeyConfig::new_raw_no_cell(library_index)
        });

        // Run the 7-step pipeline with the already-opened reader (supports streaming)
        let groups_processed = run_bam_pipeline_from_reader(
            pipeline_config,
            reader,
            input_header,
            &self.io.output,
            Some(output_header.clone()),
            // ========== grouper_fn: Create RawMiGrouper ==========
            move |_header: &Header| {
                Box::new(RawMiGrouper::new("MI", batch_size))
                    as Box<dyn Grouper<Group = RawMiGroupBatch> + Send>
            },
            // ========== process_fn: CODEC consensus calling ==========
            move |batch: RawMiGroupBatch| -> io::Result<CodecProcessedBatch> {
                // Create per-thread CODEC consensus caller
                let mut caller = CodecConsensusCaller::new_with_rejects_tracking(
                    read_name_prefix.clone(),
                    read_group_id.clone(),
                    options.clone(),
                    track_rejects,
                );

                let mut all_output = ConsensusOutput::default();
                let mut all_rejects = Vec::new();
                let mut batch_stats = CodecConsensusStats::default();
                let groups_count = batch.groups.len() as u64;

                for RawMiGroup { mi, records } in batch.groups {
                    caller.clear();

                    // Call CODEC consensus directly — records are already raw bytes!
                    let result: anyhow::Result<ConsensusOutput> = caller.consensus_reads(records);
                    match result {
                        Ok(batch_output) => {
                            all_output.merge(batch_output);
                            batch_stats.merge(caller.statistics());
                            if track_rejects {
                                all_rejects.extend(caller.take_rejected_reads());
                            }
                        }
                        Err(e) => {
                            // Handle duplex disagreement errors by merging stats
                            if e.to_string().contains("duplex disagreement") {
                                batch_stats.merge(caller.statistics());
                            } else {
                                log::warn!("CODEC consensus error for MI {mi}: {e}");
                            }
                        }
                    }
                }

                Ok(CodecProcessedBatch {
                    consensus_output: all_output,
                    rejects: all_rejects,
                    groups_count,
                    stats: batch_stats,
                })
            },
            // ========== serialize_fn: Serialize + collect metrics ==========
            move |processed: CodecProcessedBatch,
                  _header: &Header,
                  output: &mut Vec<u8>|
                  -> io::Result<u64> {
                // Collect metrics (lock-free)
                collected_metrics_for_serialize.push(CollectedCodecMetrics {
                    stats: processed.stats,
                    groups_processed: processed.groups_count,
                    rejects: processed.rejects,
                });

                // Serialize consensus reads
                let count = processed.consensus_output.count as u64;
                output.extend_from_slice(&processed.consensus_output.data);
                Ok(count)
            },
        )
        .map_err(|e| anyhow::anyhow!("Pipeline error: {e}"))?;

        // ========== Post-pipeline: Aggregate metrics and write rejects ==========
        let mut total_groups = 0u64;
        let mut merged_stats = CodecConsensusStats::default();
        let mut all_rejects: Vec<Vec<u8>> = Vec::new();

        while let Some(metrics) = collected_metrics.pop() {
            total_groups += metrics.groups_processed;
            merged_stats.merge(&metrics.stats);
            if track_rejects {
                all_rejects.extend(metrics.rejects);
            }
        }

        // Write deferred rejects as raw BAM bytes
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
                    for raw_record in &all_rejects {
                        let block_size = raw_record.len() as u32;
                        rw.get_mut()
                            .write_all(&block_size.to_le_bytes())
                            .context("Failed to write rejected read block size")?;
                        rw.get_mut()
                            .write_all(raw_record)
                            .context("Failed to write rejected read")?;
                    }
                    rw.finish(&rejects_header).context("Failed to finish rejects file")?;
                    info!("Wrote {} rejected reads", all_rejects.len());
                }
            }
        }

        // Log statistics
        info!("CODEC consensus calling complete");
        info!("Total MI groups processed: {total_groups}");
        info!("Total groups processed by pipeline: {groups_processed}");

        let metrics = merged_stats.to_metrics();
        let consensus_count = metrics.consensus_reads;
        log_consensus_summary(&metrics);

        // Write statistics file if requested
        if let Some(stats_path) = &self.stats_opts.stats {
            DelimFile::default()
                .write_tsv(stats_path, [metrics])
                .with_context(|| format!("Failed to write statistics: {}", stats_path.display()))?;
            info!("Wrote statistics to: {}", stats_path.display());
        }

        info!("Wrote {consensus_count} CODEC consensus reads");

        Ok(())
    }
}

// Helper methods moved to consensus_runner module

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::path::PathBuf;

    /// Helper to create a Codec with specified input/output paths
    fn create_codec_with_paths(input: PathBuf, output: PathBuf) -> Codec {
        Codec {
            io: BamIoOptions { input, output },
            rejects_opts: RejectsOptions::default(),
            stats_opts: StatsOptions::default(),
            read_group: ReadGroupOptions::default(),
            consensus: ConsensusCallingOptions {
                output_per_base_tags: false,
                min_consensus_base_quality: 0,
                ..ConsensusCallingOptions::default()
            },
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            min_reads: 1,
            max_reads: None,
            min_duplex_length: 1,
            single_strand_qual: None,
            outer_bases_qual: None,
            outer_bases_length: 5,
            max_duplex_disagreement_rate: 1.0,
            max_duplex_disagreements: None,
            cell_tag: None,
            sort_order: None,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        }
    }

    fn create_test_codec() -> Codec {
        create_codec_with_paths(PathBuf::from("input.bam"), PathBuf::from("output.bam"))
    }

    #[test]
    fn test_validation() {
        let mut cmd = create_test_codec();

        assert!(cmd.validate().is_ok());

        // Test invalid error rate
        cmd.consensus.error_rate_pre_umi = 0;
        assert!(cmd.validate().is_err());
        cmd.consensus.error_rate_pre_umi = 45;

        // Test invalid min reads
        cmd.min_reads = 0;
        assert!(cmd.validate().is_err());
        cmd.min_reads = 1;

        // Test invalid max < min
        cmd.max_reads = Some(0);
        assert!(cmd.validate().is_err());
        cmd.max_reads = None;

        // Test invalid duplex length
        cmd.min_duplex_length = 0;
        assert!(cmd.validate().is_err());
        cmd.min_duplex_length = 1;

        // Test invalid disagreement rate
        cmd.max_duplex_disagreement_rate = 1.5;
        assert!(cmd.validate().is_err());
        cmd.max_duplex_disagreement_rate = 1.0;

        // Note: cell_tag validation is now done at execute() time via optional_string_to_tag()
        // rather than in validate(), so we don't test it here
        cmd.cell_tag = None;
        assert!(cmd.validate().is_ok());
    }

    #[test]
    fn test_stats_to_metrics() {
        let stats = CodecConsensusStats {
            total_input_reads: 1000,
            consensus_reads_generated: 500,
            reads_filtered: 500,
            ..Default::default()
        };

        let metrics = stats.to_metrics();
        assert_eq!(metrics.total_input_reads, 1000);
        assert_eq!(metrics.consensus_reads, 500);
        assert_eq!(metrics.filtered_reads, 500);
    }

    // Integration tests
    use fgumi_lib::sam::builder::RecordBuilder;
    use noodles::sam::Header;
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use tempfile::TempDir;

    /// Helper to create a CODEC-style FR read pair with proper overlap
    /// R1 forward at start1, R2 reverse at start2, with `read_len` bases each
    fn create_codec_fr_pair_overlapping(
        mi_value: &str,
        start1: usize,
        start2: usize,
        read_len: usize,
        quals: &[u8],
    ) -> (RecordBuf, RecordBuf) {
        // Use a simple reference-matching sequence
        let seq_forward = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let seq = &seq_forward[..read_len];

        // For R2 reverse, we need the reverse complement
        let seq_rc: String = seq
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _ => c,
            })
            .collect();

        let insert_size: i32 = (start2 + read_len) as i32 - start1 as i32;

        // R1: Forward strand, first segment
        let r1 = RecordBuilder::new()
            .name(&format!("read_{mi_value}"))
            .sequence(seq)
            .qualities(quals)
            .flags(
                Flags::SEGMENTED
                    | Flags::PROPERLY_SEGMENTED
                    | Flags::FIRST_SEGMENT
                    | Flags::MATE_REVERSE_COMPLEMENTED,
            )
            .reference_sequence_id(0)
            .alignment_start(start1)
            .cigar(&format!("{read_len}M"))
            .mate_reference_sequence_id(0)
            .mate_alignment_start(start2)
            .template_length(insert_size)
            .tag("MI", mi_value)
            .tag("RG", "A")
            .build();

        // R2: Reverse strand, last segment (stored as revcomp in BAM)
        let r2 = RecordBuilder::new()
            .name(&format!("read_{mi_value}"))
            .sequence(&seq_rc)
            .qualities(quals)
            .flags(
                Flags::SEGMENTED
                    | Flags::PROPERLY_SEGMENTED
                    | Flags::LAST_SEGMENT
                    | Flags::REVERSE_COMPLEMENTED,
            )
            .reference_sequence_id(0)
            .alignment_start(start2)
            .cigar(&format!("{read_len}M"))
            .mate_reference_sequence_id(0)
            .mate_alignment_start(start1)
            .template_length(-insert_size)
            .tag("MI", mi_value)
            .tag("RG", "A")
            .build();

        (r1, r2)
    }

    fn write_codec_bam(path: &std::path::Path, records: Vec<RecordBuf>) -> Result<()> {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::{ReferenceSequence, header::Version};
        use std::num::NonZeroUsize;

        let mut header = Header::builder()
            .set_header(Map::<noodles::sam::header::record::value::map::Header>::new(Version::new(
                1, 6,
            )))
            .build();

        // Add reference sequence
        let rs = Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap());
        header.reference_sequences_mut().insert(bstr::BString::from("chr1"), rs);

        // Add read group
        let rg = Map::<noodles::sam::header::record::value::map::ReadGroup>::default();
        header.read_groups_mut().insert(bstr::BString::from("A"), rg);

        let mut writer = noodles::bam::io::writer::Builder.build_from_path(path)?;
        writer.write_header(&header)?;

        for record in &records {
            writer.write_alignment_record(&header, record)?;
        }

        Ok(())
    }

    fn read_bam_records(path: &std::path::Path) -> Result<Vec<RecordBuf>> {
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(path)?;
        let header = reader.read_header()?;
        let records: Vec<_> = reader.record_bufs(&header).collect::<std::io::Result<Vec<_>>>()?;
        Ok(records)
    }

    #[test]
    fn test_codec_execute_basic() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create FR pairs with same MI (molecule ID)
        // For proper overlap: R1 at pos 100 forward, R2 at pos 105 reverse, both 20bp
        // This gives overlap from pos 105-119 (15 bases overlap)
        let mut records = Vec::new();
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;

        // Execute should complete without errors (coverage is the goal)
        cmd.execute("test")?;

        // Output file should be created
        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_codec_execute_with_rejects() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");
        let rejects_path = dir.path().join("rejects.bam");

        let mut records = Vec::new();
        // Create multiple FR pairs
        for i in 0..3 {
            let (r1, r2) = create_codec_fr_pair_overlapping(
                &format!("UMI{i:03}"),
                100 + i * 50,
                105 + i * 50,
                20,
                &[30; 20],
            );
            records.push(r1);
            records.push(r2);
        }

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.rejects_opts.rejects = Some(rejects_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;

        cmd.execute("test")?;

        // Rejects file should exist
        assert!(rejects_path.exists());

        Ok(())
    }

    #[test]
    fn test_codec_execute_with_stats() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");
        let stats_path = dir.path().join("stats.txt");

        let mut records = Vec::new();
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.stats_opts.stats = Some(stats_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;

        cmd.execute("test")?;

        // Stats file should exist and contain data
        assert!(stats_path.exists());
        let stats_content = std::fs::read_to_string(&stats_path)?;
        assert!(!stats_content.is_empty());

        Ok(())
    }

    #[test]
    fn test_codec_execute_multithreaded() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_single = dir.path().join("output_single.bam");
        let output_multi = dir.path().join("output_multi.bam");
        let rejects_path = dir.path().join("rejects.bam");

        // Create many UMI groups to test parallel processing (25+ to trigger multiple batches with BATCH_SIZE=10)
        let mut records = Vec::new();
        for i in 0..25 {
            let (r1, r2) = create_codec_fr_pair_overlapping(
                &format!("UMI{i:03}"),
                100 + i * 50,
                105 + i * 50,
                20,
                &[30; 20],
            );
            records.push(r1);
            records.push(r2);
        }

        write_codec_bam(&input_path, records)?;

        // Single-threaded
        let mut cmd_single = create_codec_with_paths(input_path.clone(), output_single.clone());
        cmd_single.read_group.read_name_prefix = Some("codec".to_string());
        cmd_single.outer_bases_length = 0;

        cmd_single.execute("test")?;

        // Multi-threaded with rejects
        let mut cmd_multi = create_codec_with_paths(input_path, output_multi.clone());
        cmd_multi.rejects_opts.rejects = Some(rejects_path.clone());
        cmd_multi.read_group.read_name_prefix = Some("codec".to_string());
        cmd_multi.outer_bases_length = 0;
        cmd_multi.threading = ThreadingOptions::new(4);

        cmd_multi.execute("test")?;

        let records_single = read_bam_records(&output_single)?;
        let records_multi = read_bam_records(&output_multi)?;

        // Both should produce the same number of consensus reads
        assert_eq!(
            records_single.len(),
            records_multi.len(),
            "Single and multi-threaded should produce same number of reads"
        );

        Ok(())
    }

    #[test]
    fn test_codec_execute_with_per_base_tags() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut records = Vec::new();
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;
        cmd.consensus.output_per_base_tags = true; // Enable per-base tags

        cmd.execute("test")?;

        let output_records = read_bam_records(&output_path)?;
        assert!(!output_records.is_empty());

        Ok(())
    }

    #[test]
    fn test_codec_execute_with_trim() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut records = Vec::new();
        // Create reads with low quality at ends
        let quals = [5, 5, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 5, 5, 5, 5]; // Low quality at ends
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &quals);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;
        cmd.consensus.trim = true; // Enable trimming

        cmd.execute("test")?;

        // Should complete successfully (may or may not produce output depending on overlap after trim)
        assert!(output_path.exists());

        Ok(())
    }

    // Note: Overlapping consensus tests removed - CODEC does not support overlapping consensus
    // (matching fgbio's CallCodecConsensusReads which has no such option).

    #[test]
    fn test_codec_execute_with_max_reads() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create multiple read pairs for same UMI (to test max_reads downsampling)
        let mut records = Vec::new();
        for i in 0..5 {
            let (r1, r2) =
                create_codec_fr_pair_overlapping("UMI001", 100 + i, 105 + i, 20, &[30; 20]);
            records.push(r1);
            records.push(r2);
        }

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;
        cmd.max_reads = Some(2); // Limit to 2 reads per strand

        cmd.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_codec_execute_validation_errors() {
        // Test various validation errors
        let mut cmd = create_test_codec();

        // Test post-UMI error rate = 0
        cmd.consensus.error_rate_post_umi = 0;
        assert!(cmd.validate().is_err());
        cmd.consensus.error_rate_post_umi = 40;

        // Test max < min reads
        cmd.min_reads = 5;
        cmd.max_reads = Some(2);
        assert!(cmd.validate().is_err());
        cmd.min_reads = 1;
        cmd.max_reads = None;

        // Test invalid disagreement rate
        cmd.max_duplex_disagreement_rate = -0.1;
        assert!(cmd.validate().is_err());
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
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create FR pairs with same MI (molecule ID)
        let mut records = Vec::new();
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.outer_bases_length = 0;
        cmd.threading = threading;
        cmd.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_codec_processed_batch_memory_estimate() {
        let mut data = Vec::with_capacity(1024);
        data.extend_from_slice(&[0u8; 100]);
        let mut reject = Vec::with_capacity(512);
        reject.extend_from_slice(&[0u8; 50]);

        let batch = CodecProcessedBatch {
            consensus_output: ConsensusOutput { data, count: 1 },
            rejects: vec![reject],
            groups_count: 1,
            stats: CodecConsensusStats::default(),
        };

        let estimate = batch.estimate_heap_size();
        assert!(estimate >= 1024 + 512, "estimate {estimate} should be >= 1536 (capacities)");
        assert!(estimate > 1024 + 512, "estimate {estimate} should include Vec overhead");
    }
}
