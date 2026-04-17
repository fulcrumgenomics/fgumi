//! Call molecular consensus reads from UMI-grouped reads.
//!
//! This tool takes reads that have been grouped by UMI (via `group`) and generates
//! consensus reads using a likelihood-based model that accounts for sequencing errors and
//! errors introduced during sample preparation.

use crate::bam_io::{
    create_bam_reader_for_pipeline_with_opts, create_bam_writer, create_optional_bam_writer,
    create_raw_bam_reader_with_opts,
};
use crate::consensus_caller::{
    ConsensusCaller, ConsensusCallingStats, ConsensusOutput, RejectionReason,
};
use crate::logging::{OperationTimer, log_consensus_summary};
use crate::mi_group::{RawMiGroup, RawMiGroupBatch, RawMiGroupIterator, RawMiGrouper};
use crate::overlapping_consensus::{
    AgreementStrategy, CorrectionStats, DisagreementStrategy, OverlappingBasesConsensusCaller,
    apply_overlapping_consensus,
};
use crate::progress::ProgressTracker;
use crate::read_info::LibraryIndex;
use crate::unified_pipeline::{
    GroupKeyConfig, Grouper, MemoryEstimate, run_bam_pipeline_from_reader,
};
use anyhow::{Context, Result, bail};
use clap::Parser;
use fgoxide::io::DelimFile;
use fgumi_raw_bam::RawRecord;
// RejectionTracker now used via ConsensusStatsOps trait in consensus_runner
use crate::sam::SamTag;
use crate::vanilla_consensus_caller::{VanillaUmiConsensusCaller, VanillaUmiConsensusOptions};
use crossbeam_queue::SegQueue;

use log::info;
use noodles::sam::Header;
use noodles::sam::alignment::record::data::field::Tag;
use std::io;
use std::io::Write as IoWrite;
use std::sync::Arc;

use crate::commands::command::Command;
use crate::commands::common::{
    BamIoOptions, CompressionOptions, ConsensusCallingOptions, OverlappingConsensusOptions,
    QueueMemoryOptions, ReadGroupOptions, RejectsOptions, SchedulerOptions, StatsOptions,
    ThreadingOptions, build_pipeline_config,
};
use crate::commands::consensus_runner::{
    ConsensusStatsOps, create_unmapped_consensus_header, log_overlapping_stats,
};

use super::common::{MethylationRef, load_methylation_reference};

// ============================================================================
// Types for 7-step pipeline processing
// ============================================================================

/// Result from processing a batch of MI groups through consensus calling.
///
/// This type is used by the 7-step unified pipeline to pass processed results
/// from the process step to the serialize step.
struct SimplexProcessedBatch {
    /// Pre-serialized consensus reads to write to output BAM
    consensus_output: ConsensusOutput,
    /// Rejected reads as raw BAM bytes (written to rejects file if enabled)
    rejects: Vec<Vec<u8>>,
    /// Number of MI groups in this batch
    groups_count: u64,
    /// Consensus calling statistics for this batch
    stats: ConsensusCallingStats,
    /// Overlapping correction stats for this batch (if enabled)
    overlapping_stats: Option<CorrectionStats>,
}

impl MemoryEstimate for SimplexProcessedBatch {
    fn estimate_heap_size(&self) -> usize {
        self.consensus_output.data.capacity()
            + self.rejects.iter().map(Vec::capacity).sum::<usize>()
            + self.rejects.capacity() * std::mem::size_of::<Vec<u8>>()
    }
}

/// Metrics collected from each batch during parallel processing.
///
/// Used with a lock-free queue (`SegQueue`) to collect metrics from
/// parallel workers for post-pipeline aggregation.
#[derive(Default)]
struct CollectedSimplexMetrics {
    /// Consensus calling statistics
    stats: ConsensusCallingStats,
    /// Overlapping consensus stats (if enabled)
    overlapping_stats: Option<CorrectionStats>,
    /// Number of MI groups processed
    groups_processed: u64,
    /// Rejected reads as raw BAM bytes for deferred writing
    rejects: Vec<Vec<u8>>,
}

/// Calls simplex consensus sequences from reads with the same unique molecular tag.
#[derive(Debug, Parser)]
#[command(
    name = "simplex",
    about = "\x1b[38;5;180m[CONSENSUS]\x1b[0m      \x1b[36mCall simplex consensus sequences from UMI-grouped reads\x1b[0m",
    long_about = r#"
Calls consensus sequences from reads with the same unique molecular tag.

Reads with the same unique molecular tag are examined base-by-base to assess the likelihood of each base in the
source molecule. The likelihood model is as follows:

1. First, the base qualities are adjusted. The base qualities are assumed to represent the probability of a
   sequencing error (i.e. the sequencer observed the wrong base present on the cluster/flowcell/well). The base
   quality scores are converted to probabilities incorporating a probability representing the chance of an error
   from the time the unique molecular tags were integrated to just prior to sequencing. The resulting probability
   is the error rate of all processes from right after integrating the molecular tag through to the end of
   sequencing.
2. Next, a consensus sequence is called for all reads with the same unique molecular tag base-by-base. For a
   given base position in the reads, the likelihoods that an A, C, G, or T is the base for the underlying
   source molecule respectively are computed by multiplying the likelihood of each read observing the base
   position being considered. The probability of error (from 1.) is used when the observed base does not match
   the hypothesized base for the underlying source molecule, while one minus that probability is used otherwise.
   The computed likelihoods are normalized by dividing them by the sum of all four likelihoods to produce a
   posterior probability, namely the probability that the source molecule was an A, C, G, or T from just after
   integrating molecular tag through to sequencing, given the observations. The base with the maximum posterior
   probability as the consensus call, and the posterior probability is used as its raw base quality.
3. Finally, the consensus raw base quality is modified by incorporating the probability of an error prior to
   integrating the unique molecular tags. Therefore, the probability used for the final consensus base
   quality is the posterior probability of the source molecule having the consensus base given the observed
   reads with the same molecular tag, all the way from sample extraction and through sample and library
   preparation, through preparing the library for sequencing (e.g. amplification, target selection), and finally,
   through sequencing.

This tool assumes that reads with the same tag are grouped together (consecutive in the file). Also, this tool
calls each end of a pair independently, and does not jointly call bases that overlap within a pair. Insertion or
deletion errors in the reads are not considered in the consensus model.

The consensus reads produced are unaligned, due to the difficulty and error-prone nature of inferring the consensus
alignment. Consensus reads should therefore be aligned after, which should not be too expensive as likely there
are far fewer consensus reads than input raw reads.

Particular attention should be paid to setting the --min-reads parameter as this can have a dramatic effect on
both results and runtime. For libraries with low duplication rates (e.g. 100-300X exomes libraries) in which it
is desirable to retain singleton reads while making consensus reads from sets of duplicates, --min-reads=1 is
appropriate. For libraries with high duplication rates where it is desirable to only produce consensus reads
supported by 2+ reads to allow error correction, --min-reads=2 or higher is appropriate. After generation,
consensus reads can be further filtered using the filter tool. As such it is always safe to run
with --min-reads=1 and filter later, but filtering at this step can improve performance significantly.

Consensus reads have a number of additional optional tags set in the resulting BAM file. The tags break down into
those that are single-valued per read:

  consensus depth      [cD] (int)  : the maximum depth of raw reads at any point in the consensus read
  consensus min depth  [cM] (int)  : the minimum depth of raw reads at any point in the consensus read
  consensus error rate [cE] (float): the fraction of bases in raw reads disagreeing with the final consensus calls

And those that have a value per base:

  consensus depth  [cd] (short[]): the count of bases contributing to the consensus read at each position
  consensus errors [ce] (short[]): the number of bases from raw reads disagreeing with the final consensus base

The per base depths and errors are both capped at 32,767. In all cases no-calls (Ns) and bases below the
--min-input-base-quality are not counted in tag value calculations.
"#
)]
pub struct Simplex {
    /// Input/output BAM file paths
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Optional output for rejected reads
    #[command(flatten)]
    pub rejects_opts: RejectsOptions,

    /// Optional output for statistics
    #[command(flatten)]
    pub stats_opts: StatsOptions,

    /// Read group and read name prefix options
    #[command(flatten)]
    pub read_group: ReadGroupOptions,

    /// Consensus calling options (error rates, quality thresholds)
    #[command(flatten)]
    pub consensus: ConsensusCallingOptions,

    /// Overlapping bases consensus options
    #[command(flatten)]
    pub overlapping: OverlappingConsensusOptions,

    /// Threading options for parallel processing
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Minimum number of reads to produce a consensus (required, no default)
    /// Matches fgbio's `CallMolecularConsensusReads` which requires this argument
    #[arg(short = 'M', long = "min-reads")]
    pub min_reads: usize,

    /// Maximum reads to use per tag family (downsample if exceeded)
    #[arg(long = "max-reads")]
    pub max_reads: Option<usize>,

    /// Scheduler and pipeline statistics options.
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,

    /// Methylation-aware consensus calling mode.
    /// EM-Seq: C→T at ref-C = unmethylated (enzymatic conversion); TAPs: C→T at ref-C = methylated.
    /// Emits MM/ML methylation tags and cu/ct per-base count tags on consensus reads.
    /// Requires --ref.
    #[arg(long = "methylation-mode", value_enum)]
    pub methylation_mode: Option<crate::commands::common::MethylationModeArg>,

    /// Path to the reference FASTA file (required when --methylation-mode is set).
    #[arg(long = "ref")]
    pub reference: Option<std::path::PathBuf>,
}

impl Command for Simplex {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Start timing
        let timer = OperationTimer::new("Calling simplex consensus");

        // Validate inputs
        self.io.validate()?;

        if let Some(max) = self.max_reads {
            if max < self.min_reads {
                bail!("--max-reads ({}) must be >= --min-reads ({})", max, self.min_reads);
            }
        }

        info!("Starting Simplex");
        info!("Input: {}", self.io.input.display());
        info!("Output: {}", self.io.output.display());
        info!("Min reads: {}", self.min_reads);
        if let Some(max) = self.max_reads {
            info!("Max reads: {max}");
        }
        info!("Error rate pre-UMI: Q{}", self.consensus.error_rate_pre_umi);
        info!("Error rate post-UMI: Q{}", self.consensus.error_rate_post_umi);

        // Open input using streaming-capable reader for pipeline use
        // Get threading configuration
        let writer_threads = self.threading.num_threads();

        let (reader, header) = create_bam_reader_for_pipeline_with_opts(
            &self.io.input,
            self.io.pipeline_reader_opts(),
        )?;

        // Create output header (cleared for unmapped consensus reads)
        let output_header = create_unmapped_consensus_header(
            &header,
            &self.read_group.read_group_id,
            "Read group",
            command_line,
        )?;

        // Use library name from header if no prefix is specified (like fgbio)
        let read_name_prefix = self.read_group.prefix_or_from_header(&header);

        let cell_tag = Tag::from(SamTag::CB);

        // Enable rejects tracking if rejects file is specified
        let track_rejects = self.rejects_opts.is_enabled();

        // Load reference for methylation-aware consensus calling
        if self.reference.is_some() && self.methylation_mode.is_none() {
            bail!("--ref requires --methylation-mode to be set");
        }
        let methylation_mode =
            crate::commands::common::resolve_methylation_mode(self.methylation_mode);
        let methylation_ref: MethylationRef =
            load_methylation_reference(methylation_mode, &self.reference, &header)?;

        // Track overlapping consensus settings (callers created per-thread in threaded mode)
        let overlapping_enabled = self.overlapping.is_enabled();
        if overlapping_enabled {
            info!("Overlapping consensus calling enabled");
        }

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
                output_header.clone(),
                read_name_prefix.clone(),
                track_rejects,
                methylation_ref.clone(),
                methylation_mode,
            );
            timer.log_completion(0); // Completion logged in execute_threads_mode
            return result;
        }

        // Drop the reader for single-threaded mode - we use MiGroupIterator which needs its own reader
        drop(reader);

        // ============================================================
        // For non-pipeline modes, create output writers here
        // ============================================================

        // Open output - use multi-threaded BGZF writer
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

        let options = VanillaUmiConsensusOptions {
            tag: "MI".to_string(),
            error_rate_pre_umi: self.consensus.error_rate_pre_umi,
            error_rate_post_umi: self.consensus.error_rate_post_umi,
            min_input_base_quality: self.consensus.min_input_base_quality,
            min_reads: self.min_reads,
            max_reads: self.max_reads,
            produce_per_base_tags: self.consensus.output_per_base_tags,
            seed: Some(42), // Hard-coded seed for reproducible downsampling
            trim: self.consensus.trim,
            min_consensus_base_quality: self.consensus.min_consensus_base_quality,
            cell_tag: Some(cell_tag),
            methylation_mode,
        };

        // Create a single-threaded caller for stats collection
        let mut caller = VanillaUmiConsensusCaller::new_with_rejects_tracking(
            read_name_prefix.clone(),
            self.read_group.read_group_id.clone(),
            options.clone(),
            track_rejects,
        );

        // Set reference for methylation-aware consensus if enabled
        if let Some((ref reference, ref ref_names)) = methylation_ref {
            caller.set_reference(Arc::clone(reference), Arc::clone(ref_names));
        }

        // Accumulator for overlapping stats from parallel processing
        let mut merged_overlapping_stats = CorrectionStats::new();

        // Track progress (count records written, not UMI groups)
        let mut record_count: usize = 0;
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        // Create the MI group iterator for single-threaded streaming (raw bytes)
        let (mut raw_reader, _) =
            create_raw_bam_reader_with_opts(&self.io.input, 1, self.io.pipeline_reader_opts())?;
        let raw_record_iter = std::iter::from_fn(move || {
            let mut record = RawRecord::new();
            match raw_reader.read_record(&mut record) {
                Ok(0) => None, // EOF
                Ok(_) => Some(Ok(record)),
                Err(e) => Some(Err(e.into())),
            }
        });
        let mi_group_iter =
            RawMiGroupIterator::new(raw_record_iter, "MI").with_cell_tag(Some(*b"CB"));
        // Single-threaded streaming processing
        // Create overlapping consensus caller for single-threaded mode
        let mut overlapping_caller = if overlapping_enabled {
            Some(OverlappingBasesConsensusCaller::new(
                AgreementStrategy::Consensus,
                DisagreementStrategy::Consensus,
            ))
        } else {
            None
        };

        for result in mi_group_iter {
            let (umi, mut records) = result.context("Failed to read MI group")?;

            // Apply overlapping consensus if enabled (modifies raw bytes in-place)
            if let Some(ref mut oc) = overlapping_caller {
                apply_overlapping_consensus(&mut records, oc)?;
            }

            // Call consensus — mi_group now yields Vec<RawRecord> directly, which is
            // what ConsensusCaller::consensus_reads expects.
            let output = caller
                .consensus_reads(records)
                .with_context(|| format!("Failed to call consensus for UMI: {umi}"))?;

            let batch_size = output.count;
            record_count += batch_size;

            // Write pre-serialized consensus reads to output
            writer.get_mut().write_all(&output.data).context("Failed to write consensus read")?;

            // Write rejected reads if tracking is enabled
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

            progress.log_if_needed(batch_size as u64);
        }

        // For single-threaded, use the caller's stats and merge overlapping stats
        let merged_stats = caller.statistics();
        if let Some(ref oc) = overlapping_caller {
            merged_overlapping_stats.merge(oc.stats());
        }

        progress.log_final();

        // Finish the buffered writer (flush remaining records and wait for writer thread)
        writer.into_inner().finish().context("Failed to finish output BAM")?;

        // Log overlapping consensus statistics if enabled
        if overlapping_enabled {
            log_overlapping_stats(&merged_overlapping_stats);
        }

        // Log statistics and write to file
        info!("Consensus calling complete");
        info!("Total records processed: {record_count}");

        let metrics = merged_stats.to_metrics();
        let consensus_count = metrics.consensus_reads;
        log_consensus_summary(&metrics);

        if let Some(stats_path) = &self.stats_opts.stats {
            // Convert to fgbio-compatible vertical key-value-description format
            let kv_metrics = metrics.to_kv_metrics();
            DelimFile::default()
                .write_tsv(stats_path, kv_metrics)
                .with_context(|| format!("Failed to write statistics: {}", stats_path.display()))?;
            info!("Wrote statistics to: {}", stats_path.display());
        }

        timer.log_completion(consensus_count);

        // Close rejects writer if it was opened
        if let Some(rw) = rejects_writer {
            rw.into_inner().finish().context("Failed to finish rejects file")?;
            info!("Rejected reads written successfully");
        }

        Ok(())
    }
}
impl Simplex {
    /// Execute using 7-step unified pipeline with --threads.
    ///
    /// This method is called when `--threads N` is specified with N > 1.
    /// It uses the lock-free 7-step unified pipeline for maximum performance.
    #[expect(clippy::too_many_arguments, reason = "pipeline setup needs all configuration")]
    fn execute_threads_mode(
        &self,
        num_threads: usize,
        reader: Box<dyn std::io::Read + Send>,
        input_header: Header,
        output_header: Header,
        read_name_prefix: String,
        track_rejects: bool,
        methylation_ref: MethylationRef,
        methylation_mode: fgumi_consensus::MethylationMode,
    ) -> Result<()> {
        // Configure pipeline
        let mut pipeline_config = build_pipeline_config(
            &self.scheduler_opts,
            &self.compression,
            &self.queue_memory,
            num_threads,
        )?;

        // Lock-free metrics collection
        let collected_metrics: Arc<SegQueue<CollectedSimplexMetrics>> = Arc::new(SegQueue::new());
        let collected_metrics_for_serialize = Arc::clone(&collected_metrics);

        // Capture configuration for closures
        let tag_str = "MI".to_string();
        let min_reads = self.min_reads;
        let max_reads = self.max_reads;
        let error_rate_pre_umi = self.consensus.error_rate_pre_umi;
        let error_rate_post_umi = self.consensus.error_rate_post_umi;
        let min_input_base_quality = self.consensus.min_input_base_quality;
        let output_per_base_tags = self.consensus.output_per_base_tags;
        let min_consensus_base_quality = self.consensus.min_consensus_base_quality;
        let trim = self.consensus.trim;
        let overlapping_enabled = self.overlapping.is_enabled();
        let read_group_id = self.read_group.read_group_id.clone();
        let cell_tag = Tag::from(SamTag::CB);
        let batch_size = 50; // MI groups per batch (reduced for memory efficiency)

        // Create options for consensus caller
        let options = VanillaUmiConsensusOptions {
            tag: tag_str.clone(),
            error_rate_pre_umi,
            error_rate_post_umi,
            min_input_base_quality,
            min_reads,
            max_reads,
            produce_per_base_tags: output_per_base_tags,
            seed: Some(42),
            trim,
            min_consensus_base_quality,
            cell_tag: Some(cell_tag),
            methylation_mode,
        };

        // Clone input_header before pipeline (needed for rejects writing)
        let rejects_header = input_header.clone();

        // Enable raw-byte mode: skip noodles decode/encode for CPU savings
        let library_index = LibraryIndex::from_header(&input_header);
        pipeline_config.group_key_config = Some(GroupKeyConfig::new_raw(library_index, cell_tag));

        // Run the 7-step pipeline with the already-opened reader (supports streaming)
        let groups_processed = run_bam_pipeline_from_reader(
            pipeline_config,
            reader,
            input_header,
            &self.io.output,
            Some(output_header.clone()),
            // ========== grouper_fn: Create RawMiGrouper ==========
            move |_header: &Header| {
                Box::new(RawMiGrouper::new("MI", batch_size).with_cell_tag(Some(*b"CB")))
                    as Box<dyn Grouper<Group = RawMiGroupBatch> + Send>
            },
            // ========== process_fn: Consensus calling ==========
            move |batch: RawMiGroupBatch| -> io::Result<SimplexProcessedBatch> {
                // Create per-thread consensus caller
                let mut caller = VanillaUmiConsensusCaller::new_with_rejects_tracking(
                    read_name_prefix.clone(),
                    read_group_id.clone(),
                    options.clone(),
                    track_rejects,
                );

                // Set reference for methylation-aware consensus if enabled
                if let Some((ref reference, ref ref_names)) = methylation_ref {
                    caller.set_reference(Arc::clone(reference), Arc::clone(ref_names));
                }

                // Create overlapping caller if enabled
                let mut overlapping_caller = if overlapping_enabled {
                    Some(OverlappingBasesConsensusCaller::new(
                        AgreementStrategy::Consensus,
                        DisagreementStrategy::Consensus,
                    ))
                } else {
                    None
                };

                let mut all_output = ConsensusOutput::default();
                let mut all_rejects = Vec::new();
                let mut batch_stats = ConsensusCallingStats::new();
                let mut batch_overlapping = CorrectionStats::new();
                let groups_count = batch.groups.len() as u64;

                for RawMiGroup { mi, records: mut raw_records } in batch.groups {
                    caller.clear();

                    // Skip if below min_reads threshold
                    if raw_records.len() < min_reads {
                        batch_stats.record_input(raw_records.len());
                        batch_stats.record_rejection(
                            RejectionReason::InsufficientReads,
                            raw_records.len(),
                        );
                        if track_rejects {
                            // TODO(#272): bridge – remove when rejects Vec accepts RawRecord
                            all_rejects.extend(raw_records.into_iter().map(RawRecord::into_inner));
                        }
                        continue;
                    }

                    // Apply overlapping consensus if enabled (modifies raw bytes in-place)
                    if let Some(ref mut oc) = overlapping_caller {
                        oc.reset_stats();
                        if apply_overlapping_consensus(&mut raw_records, oc).is_err() {
                            batch_overlapping.merge(oc.stats());
                            batch_stats.record_input(raw_records.len());
                            batch_stats.record_rejection(RejectionReason::Other, raw_records.len());
                            if track_rejects {
                                // TODO(#272): bridge – remove when rejects Vec accepts RawRecord
                                all_rejects
                                    .extend(raw_records.into_iter().map(RawRecord::into_inner));
                            }
                            continue;
                        }
                        batch_overlapping.merge(oc.stats());
                    }

                    // Call consensus — mi_group now yields Vec<RawRecord> directly.
                    match caller.consensus_reads(raw_records) {
                        Ok(batch_output) => {
                            all_output.merge(batch_output);
                            batch_stats.merge(&caller.statistics());
                            if track_rejects {
                                all_rejects.extend(caller.take_rejected_reads());
                            }
                        }
                        Err(e) => {
                            log::warn!("Consensus error for MI {mi}: {e}");
                        }
                    }
                }

                Ok(SimplexProcessedBatch {
                    consensus_output: all_output,
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
            move |processed: SimplexProcessedBatch,
                  _header: &Header,
                  output: &mut Vec<u8>|
                  -> io::Result<u64> {
                // Collect metrics (lock-free)
                collected_metrics_for_serialize.push(CollectedSimplexMetrics {
                    stats: processed.stats,
                    overlapping_stats: processed.overlapping_stats,
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
        let mut merged_stats = ConsensusCallingStats::new();
        let mut merged_overlapping_stats = CorrectionStats::new();
        let mut all_rejects: Vec<Vec<u8>> = Vec::new();

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

        // Write deferred rejects as raw BAM bytes
        if track_rejects && !all_rejects.is_empty() {
            if let Some(rejects_path) = &self.rejects_opts.rejects {
                let writer_threads = self.threading.num_threads();
                let rejects_writer = create_optional_bam_writer(
                    Some(rejects_path),
                    &rejects_header,
                    writer_threads,
                    self.compression.compression_level,
                )?;
                if let Some(mut rw) = rejects_writer {
                    for raw_record in &all_rejects {
                        let block_size = raw_record.len() as u32;
                        rw.get_mut()
                            .write_all(&block_size.to_le_bytes())
                            .context("Failed to write rejected read block size")?;
                        rw.get_mut()
                            .write_all(raw_record)
                            .context("Failed to write rejected read")?;
                    }
                    rw.into_inner().finish().context("Failed to finish rejects file")?;
                    info!("Wrote {} rejected reads", all_rejects.len());
                }
            }
        }

        // Log overlapping consensus statistics if enabled
        if self.overlapping.is_enabled() {
            log_overlapping_stats(&merged_overlapping_stats);
        }

        // Log statistics and write to file
        info!("Consensus calling complete");
        info!("Total MI groups processed: {total_groups}");
        info!("Total groups processed by pipeline: {groups_processed}");

        let metrics = merged_stats.to_metrics();
        let consensus_count = metrics.consensus_reads;
        log_consensus_summary(&metrics);

        if let Some(stats_path) = &self.stats_opts.stats {
            // Convert to fgbio-compatible vertical key-value-description format
            let kv_metrics = metrics.to_kv_metrics();
            DelimFile::default()
                .write_tsv(stats_path, kv_metrics)
                .with_context(|| format!("Failed to write statistics: {}", stats_path.display()))?;
            info!("Wrote statistics to: {}", stats_path.display());
        }

        info!("Wrote {consensus_count} consensus reads");

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::consensus::ConsensusKvMetric;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use rstest::rstest;
    use std::path::PathBuf;

    /// Creates a test Simplex command instance with default parameters.
    ///
    /// Generates a Simplex command configured with standard test values,
    /// including file paths, UMI tag, error rates, and threading settings.
    ///
    /// # Returns
    ///
    /// A `Simplex` instance configured for testing
    fn create_test_simplex() -> Simplex {
        create_simplex_with_paths(PathBuf::from("test.bam"), PathBuf::from("out.bam"))
    }

    /// Creates a Simplex command with the given input/output paths and default parameters.
    fn create_simplex_with_paths(input: PathBuf, output: PathBuf) -> Simplex {
        Simplex {
            io: BamIoOptions::new(input, output),
            rejects_opts: RejectsOptions::default(),
            stats_opts: StatsOptions::default(),
            read_group: ReadGroupOptions::default(),
            consensus: ConsensusCallingOptions {
                output_per_base_tags: false,
                min_consensus_base_quality: 0,
                ..ConsensusCallingOptions::default()
            },
            overlapping: OverlappingConsensusOptions { consensus_call_overlapping_bases: false },
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            queue_memory: QueueMemoryOptions::default(),
            min_reads: 1,
            max_reads: None,
            scheduler_opts: SchedulerOptions::default(),
            methylation_mode: None,
            reference: None,
        }
    }

    #[test]
    fn test_default_parameters() {
        let cmd = create_test_simplex();
        assert_eq!(cmd.consensus.error_rate_pre_umi, 45);
        assert_eq!(cmd.consensus.error_rate_post_umi, 40);
        assert_eq!(cmd.consensus.min_input_base_quality, 10);
        assert_eq!(cmd.min_reads, 1);
        assert_eq!(cmd.max_reads, None);
        assert!(!cmd.consensus.output_per_base_tags);
        assert!(cmd.threading.is_single_threaded());
        // Note: create_simplex_with_paths disables overlapping by default for tests
        assert!(!cmd.overlapping.consensus_call_overlapping_bases);
    }

    #[test]
    fn test_custom_parameters() {
        let mut cmd = create_test_simplex();
        cmd.min_reads = 3;
        cmd.max_reads = Some(10);
        cmd.consensus.output_per_base_tags = true;
        cmd.threading = ThreadingOptions::new(4);

        assert_eq!(cmd.min_reads, 3);
        assert_eq!(cmd.max_reads, Some(10));
        assert!(cmd.consensus.output_per_base_tags);
        assert_eq!(cmd.threading.threads, Some(4));
    }

    #[test]
    fn test_missing_input_file_fails() {
        let cmd = create_test_simplex();
        let result = cmd.execute("test");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("does not exist"));
    }

    // ========================================================================
    // Integration Tests
    // ========================================================================

    use crate::sam::builder::SamBuilder;
    use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
    use std::collections::HashMap;
    use tempfile::TempDir;

    /// Helper struct for managing temporary test file paths.
    struct TestPaths {
        #[allow(dead_code)]
        dir: TempDir,
        pub input: PathBuf,
        pub output: PathBuf,
        pub rejects: PathBuf,
        pub stats: PathBuf,
    }

    impl TestPaths {
        fn new() -> Result<Self> {
            let dir = TempDir::new()?;
            Ok(Self {
                input: dir.path().join("input.bam"),
                output: dir.path().join("output.bam"),
                rejects: dir.path().join("rejects.bam"),
                stats: dir.path().join("stats.txt"),
                dir,
            })
        }
    }

    /// Helper to read all records from a BAM file
    fn read_bam_records(path: &std::path::Path) -> Result<Vec<RecordBuf>> {
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(path)?;
        let header = reader.read_header()?;
        let mut records = Vec::new();

        for result in reader.records() {
            let record = result?;
            let record_buf = RecordBuf::try_from_alignment_record(&header, &record)?;
            records.push(record_buf);
        }

        Ok(records)
    }

    /// Helper to get a string tag from a record
    fn get_string_tag(record: &RecordBuf, tag_name: &str) -> Option<String> {
        let tag_bytes = tag_name.as_bytes();
        let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

        record.data().get(&tag).and_then(|v| {
            if let noodles::sam::alignment::record_buf::data::field::Value::String(s) = v {
                Some(String::from_utf8_lossy(s).to_string())
            } else {
                None
            }
        })
    }

    /// Helper to get an integer tag from a record
    fn get_int_tag(record: &RecordBuf, tag_name: &str) -> Option<i64> {
        let tag_bytes = tag_name.as_bytes();
        let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

        record.data().get(&tag).and_then(|v| match v {
            noodles::sam::alignment::record_buf::data::field::Value::Int8(i) => Some(*i as i64),
            noodles::sam::alignment::record_buf::data::field::Value::UInt8(i) => Some(*i as i64),
            noodles::sam::alignment::record_buf::data::field::Value::Int16(i) => Some(*i as i64),
            noodles::sam::alignment::record_buf::data::field::Value::UInt16(i) => Some(*i as i64),
            noodles::sam::alignment::record_buf::data::field::Value::Int32(i) => Some(*i as i64),
            noodles::sam::alignment::record_buf::data::field::Value::UInt32(i) => Some(*i as i64),
            _ => None,
        })
    }

    #[test]
    fn test_end_to_end_paired_end_workflow() -> Result<()> {
        // Create test data similar to Scala test: 1000 UMI groups with 2 pairs each
        let mut builder = SamBuilder::new_unmapped();

        for idx in 0..1000 {
            let umi = format!("GATTACA:{idx}");
            let mut attrs = HashMap::new();
            attrs.insert("MI", BufValue::from(umi.clone()));
            attrs.insert("RX", BufValue::from("ACGT-TGCA"));

            // Add 2 pairs per UMI group
            builder.add_pair_with_attrs(
                &format!("READ:{}", 2 * idx),
                None,
                None,
                true,
                true,
                &attrs,
            );
            builder.add_pair_with_attrs(
                &format!("READ:{}", 2 * idx + 1),
                None,
                None,
                true,
                true,
                &attrs,
            );
        }

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.read_group.read_group_id = "ABC".to_string();

        cmd.execute("test")?;

        // Verify output
        let records = read_bam_records(&paths.output)?;

        // Should have 2000 consensus reads (1000 UMI groups × 2 reads per pair)
        assert_eq!(records.len(), 2000, "Should have 2000 consensus reads");

        // Verify first and second of pair counts
        let first_count = records.iter().filter(|r| r.flags().is_first_segment()).count();
        let second_count = records.iter().filter(|r| r.flags().is_last_segment()).count();
        assert_eq!(first_count, 1000, "Should have 1000 first-of-pair reads");
        assert_eq!(second_count, 1000, "Should have 1000 second-of-pair reads");

        // Verify all reads have expected attributes
        for record in &records {
            // Check sequence
            assert_eq!(record.sequence().len(), 100, "Sequence length should be 100");

            // Check read group
            let rg = get_string_tag(record, "RG");
            assert_eq!(rg.as_deref(), Some("ABC"), "Read group should be ABC");

            // Check consensus tags are present
            let cd_tag = get_int_tag(record, "cD");
            assert!(cd_tag.is_some(), "cD tag should be present");
            assert_eq!(
                cd_tag.expect("cD tag should have a value"),
                2,
                "Depth should be 2 (2 reads per UMI)"
            );

            // MI and RX tags should be preserved by the consensus caller
            let mi_tag = get_string_tag(record, "MI");
            assert!(mi_tag.is_some(), "MI tag should be preserved");
        }

        Ok(())
    }

    #[rstest]
    #[case::fast_path(ThreadingOptions::none())]
    #[case::threaded(ThreadingOptions::new(2))]
    fn test_end_to_end_single_end_workflow(#[case] threading: ThreadingOptions) -> Result<()> {
        // Verify CB partitioning: two groups share the same MI but differ in CB.
        // They must produce two separate consensus reads, one per cell barcode.
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs_a = HashMap::new();
        attrs_a.insert("RX", BufValue::from("ACGT"));
        attrs_a.insert("MI", BufValue::from("shared"));
        attrs_a.insert("CB", BufValue::from("AB"));

        // 3 fragments from cell AB
        builder.add_frag_with_attrs("a1", None, true, &attrs_a);
        builder.add_frag_with_attrs("a2", None, true, &attrs_a);
        builder.add_frag_with_attrs("a3", None, true, &attrs_a);

        let mut attrs_b = HashMap::new();
        attrs_b.insert("RX", BufValue::from("ACGT"));
        attrs_b.insert("MI", BufValue::from("shared"));
        attrs_b.insert("CB", BufValue::from("CD"));

        // 2 fragments from cell CD — same MI, different cell
        builder.add_frag_with_attrs("b1", None, true, &attrs_b);
        builder.add_frag_with_attrs("b2", None, true, &attrs_b);

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.read_group.read_group_id = "ABC".to_string();
        cmd.threading = threading;

        cmd.execute("test")?;

        // Verify output
        let records = read_bam_records(&paths.output)?;

        // Should have 2 consensus reads — one per cell barcode
        assert_eq!(records.len(), 2, "Should have 2 consensus reads (one per cell barcode)");

        // All should be unpaired
        assert!(records.iter().all(|r| !r.flags().is_segmented()), "All reads should be unpaired");

        // Both cell barcodes must survive independently
        let observed_cbs: std::collections::HashSet<_> =
            records.iter().filter_map(|r| get_string_tag(r, "CB")).collect();
        assert_eq!(
            observed_cbs,
            ["AB".to_string(), "CD".to_string()].into_iter().collect(),
            "Both cell barcodes should produce independent consensus reads"
        );

        // Verify per-record attributes
        for record in &records {
            let rg = get_string_tag(record, "RG");
            assert_eq!(rg.as_deref(), Some("ABC"), "Read group should be ABC");

            let cd_tag = get_int_tag(record, "cD");
            assert!(cd_tag.is_some(), "cD tag should be present");
            assert!(cd_tag.expect("cD tag should have a value") >= 2, "Depth should be at least 2");

            assert_eq!(record.sequence().len(), 100, "Sequence length should be 100");
        }

        Ok(())
    }

    #[test]
    fn test_min_reads_filtering() -> Result<()> {
        // Test that groups with fewer than min_reads are rejected
        let mut builder = SamBuilder::new_unmapped();

        // Group 1: 3 reads (should pass with min_reads=3)
        let mut attrs1 = HashMap::new();
        attrs1.insert("MI", BufValue::from("group1"));
        builder.add_frag_with_attrs("a1", None, true, &attrs1);
        builder.add_frag_with_attrs("a2", None, true, &attrs1);
        builder.add_frag_with_attrs("a3", None, true, &attrs1);

        // Group 2: 2 reads (should be rejected with min_reads=3)
        let mut attrs2 = HashMap::new();
        attrs2.insert("MI", BufValue::from("group2"));
        builder.add_frag_with_attrs("b1", None, true, &attrs2);
        builder.add_frag_with_attrs("b2", None, true, &attrs2);

        // Group 3: 1 read (should be rejected with min_reads=3)
        let mut attrs3 = HashMap::new();
        attrs3.insert("MI", BufValue::from("group3"));
        builder.add_frag_with_attrs("c1", None, true, &attrs3);

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.min_reads = 3;

        cmd.execute("test")?;

        // Verify output
        let records = read_bam_records(&paths.output)?;

        // Should have only 1 consensus read (from group1)
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        // Check consensus depth tag is present (validates we got the group with 3 reads)
        let cd_tag = get_int_tag(&records[0], "cD");
        assert_eq!(
            cd_tag.expect("cD tag should have a value"),
            3,
            "Depth should be 3 for the passing group"
        );

        Ok(())
    }

    #[test]
    fn test_per_base_tags_generation() -> Result<()> {
        // Test that per-base tags (cd, ce) are generated when requested
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("test_umi"));

        // Add 5 reads to get sufficient depth
        for i in 0..5 {
            builder.add_frag_with_attrs(&format!("read{i}"), None, true, &attrs);
        }

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.consensus.output_per_base_tags = true;

        cmd.execute("test")?;

        // Verify output
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        let record = &records[0];

        // Check per-read tags (cD, cM, cE) are present
        let cd_tag = get_int_tag(record, "cD");
        assert!(cd_tag.is_some(), "cD tag should be present");
        assert_eq!(cd_tag.expect("cD tag should have a value"), 5, "Max depth should be 5");

        let cm_tag = get_int_tag(record, "cM");
        assert!(cm_tag.is_some(), "cM tag should be present");
        assert_eq!(cm_tag.expect("cM tag should have a value"), 5, "Min depth should be 5");

        // Check per-base tags (cd, ce) are present
        let tag_bytes = "cd".as_bytes();
        let cd_array_tag = Tag::from([tag_bytes[0], tag_bytes[1]]);
        assert!(record.data().get(&cd_array_tag).is_some(), "cd per-base tag should be present");

        let tag_bytes = "ce".as_bytes();
        let ce_array_tag = Tag::from([tag_bytes[0], tag_bytes[1]]);
        assert!(record.data().get(&ce_array_tag).is_some(), "ce per-base tag should be present");

        Ok(())
    }

    #[test]
    fn test_multithreading() -> Result<()> {
        // Test that multithreading produces same results as single-threaded
        let mut builder = SamBuilder::new_unmapped();

        // Create 100 UMI groups with 2 reads each
        for idx in 0..100 {
            let mut attrs = HashMap::new();
            attrs.insert("MI", BufValue::from(format!("umi_{idx}")));
            builder.add_frag_with_attrs(&format!("read_{idx}a"), None, true, &attrs);
            builder.add_frag_with_attrs(&format!("read_{idx}b"), None, true, &attrs);
        }

        let paths = TestPaths::new()?;
        let output_multi_path = paths.dir.path().join("output_multi.bam");
        builder.write(&paths.input)?;

        // Run with single thread
        let cmd_single = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd_single.execute("test")?;

        // Run with multiple threads
        let mut cmd_multi =
            create_simplex_with_paths(paths.input.clone(), output_multi_path.clone());
        cmd_multi.threading = ThreadingOptions::new(4);
        cmd_multi.execute("test")?;

        // Verify both outputs have same number of records
        let records_single = read_bam_records(&paths.output)?;
        let records_multi = read_bam_records(&output_multi_path)?;

        assert_eq!(
            records_single.len(),
            records_multi.len(),
            "Single and multi-threaded should produce same number of records"
        );
        assert_eq!(records_single.len(), 100, "Should have 100 consensus reads");

        Ok(())
    }

    #[test]
    fn test_max_reads_downsampling() -> Result<()> {
        // Test that max_reads properly downsamples large UMI groups
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("large_group"));

        // Add 20 reads to this group
        for i in 0..20 {
            builder.add_frag_with_attrs(&format!("read{i}"), None, true, &attrs);
        }

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.max_reads = Some(10);

        cmd.execute("test")?;

        // Verify output
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        // Check depth tag shows downsampling was applied (max depth should be 10)
        let cd_tag = get_int_tag(&records[0], "cD");
        assert!(cd_tag.is_some(), "cD tag should be present");
        assert!(
            cd_tag.expect("cD tag should have a value") <= 10,
            "Max depth should be <= 10 due to downsampling"
        );

        Ok(())
    }

    #[test]
    fn test_max_reads_less_than_min_reads_fails() {
        // Create a dummy BAM file so validation gets to tag check
        let builder = SamBuilder::new_unmapped();
        let paths = TestPaths::new().expect("failed to create test paths");
        builder.write(&paths.input).expect("failed to write test BAM");

        let mut cmd = create_simplex_with_paths(paths.input.clone(), PathBuf::from("out.bam"));
        cmd.min_reads = 5;
        cmd.max_reads = Some(3); // Invalid: less than min_reads

        let result = cmd.execute("test");
        assert!(result.is_err());
        let error_msg = result.unwrap_err().to_string();
        assert!(error_msg.contains("--max-reads"));
        assert!(error_msg.contains("--min-reads"));
    }

    #[test]
    fn test_statistics_file_generation() -> Result<()> {
        let mut builder = SamBuilder::new_unmapped();

        // Group 1: passes (3 reads)
        let mut attrs1 = HashMap::new();
        attrs1.insert("MI", BufValue::from("pass"));
        builder.add_frag_with_attrs("p1", None, true, &attrs1);
        builder.add_frag_with_attrs("p2", None, true, &attrs1);
        builder.add_frag_with_attrs("p3", None, true, &attrs1);

        // Group 2: filtered (1 read, min_reads=2)
        let mut attrs2 = HashMap::new();
        attrs2.insert("MI", BufValue::from("fail"));
        builder.add_frag_with_attrs("f1", None, true, &attrs2);

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.stats_opts.stats = Some(paths.stats.clone());
        cmd.min_reads = 2;

        cmd.execute("test")?;

        // Verify stats file was created and contains data
        assert!(&paths.stats.exists(), "Stats file should exist");

        // Read and verify TSV format (now vertical key-value-description format)
        let kv_metrics: Vec<ConsensusKvMetric> =
            DelimFile::default().read_tsv(&paths.stats).expect("Failed to read metrics file");

        // Should have multiple rows (one per metric)
        assert!(!kv_metrics.is_empty(), "Should have metrics");

        // Verify expected keys are present
        let keys: Vec<&str> = kv_metrics.iter().map(|m| m.key.as_str()).collect();
        assert!(keys.contains(&"raw_reads_considered"), "Should have raw_reads_considered");
        assert!(keys.contains(&"raw_reads_rejected"), "Should have raw_reads_rejected");
        assert!(keys.contains(&"consensus_reads_emitted"), "Should have consensus_reads_emitted");

        // Verify raw_reads_considered has a value
        let raw_reads = kv_metrics
            .iter()
            .find(|m| m.key == "raw_reads_considered")
            .expect("Should have raw_reads_considered");
        let count: u64 = raw_reads.value.parse().expect("Should be a number");
        assert!(count > 0, "Stats should have total reads");

        Ok(())
    }

    #[test]
    fn test_custom_read_name_prefix() -> Result<()> {
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("test"));
        builder.add_frag_with_attrs("read1", None, true, &attrs);
        builder.add_frag_with_attrs("read2", None, true, &attrs);

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.read_group.read_name_prefix = Some("MYCONSENSUS".to_string());

        cmd.execute("test")?;

        // Verify output read name has custom prefix
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        let read_name = records[0].name().map(std::string::ToString::to_string).unwrap_or_default();
        assert!(read_name.starts_with("MYCONSENSUS"), "Read name should start with custom prefix");

        Ok(())
    }

    #[test]
    fn test_rejects_file_generation() -> Result<()> {
        // Test that rejected reads are written to rejects file
        let mut builder = SamBuilder::new_unmapped();

        // Group 1: passes (3 reads)
        let mut attrs1 = HashMap::new();
        attrs1.insert("MI", BufValue::from("pass"));
        builder.add_frag_with_attrs("p1", None, true, &attrs1);
        builder.add_frag_with_attrs("p2", None, true, &attrs1);
        builder.add_frag_with_attrs("p3", None, true, &attrs1);

        // Group 2: filtered (1 read, min_reads=2)
        let mut attrs2 = HashMap::new();
        attrs2.insert("MI", BufValue::from("fail"));
        builder.add_frag_with_attrs("f1", None, true, &attrs2);

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.rejects_opts.rejects = Some(paths.rejects.clone());
        cmd.min_reads = 2;

        cmd.execute("test")?;

        // Verify rejects file exists and contains rejected reads
        assert!(&paths.rejects.exists(), "Rejects file should exist");

        let reject_records = read_bam_records(&paths.rejects)?;
        assert_eq!(reject_records.len(), 1, "Should have 1 rejected read");

        // Verify the rejected read has the correct UMI
        let umi = get_string_tag(&reject_records[0], "MI");
        assert_eq!(umi.as_deref(), Some("fail"), "Rejected read should have 'fail' UMI");

        Ok(())
    }

    #[test]
    fn test_multithreaded_with_rejects() -> Result<()> {
        // Test multi-threaded execution with rejects tracking
        let mut builder = SamBuilder::new_unmapped();

        // Create 50 UMI groups
        for idx in 0..50 {
            let mut attrs = HashMap::new();
            if idx % 3 == 0 {
                // Every 3rd group has only 1 read (will be rejected with min_reads=2)
                attrs.insert("MI", BufValue::from(format!("umi_{idx}")));
                builder.add_frag_with_attrs(&format!("read_{idx}"), None, true, &attrs);
            } else {
                // Other groups have 2 reads (will pass)
                attrs.insert("MI", BufValue::from(format!("umi_{idx}")));
                builder.add_frag_with_attrs(&format!("read_{idx}a"), None, true, &attrs);
                builder.add_frag_with_attrs(&format!("read_{idx}b"), None, true, &attrs);
            }
        }

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.rejects_opts.rejects = Some(paths.rejects.clone());
        cmd.min_reads = 2;
        cmd.threading = ThreadingOptions::new(4);

        cmd.execute("test")?;

        // Verify outputs
        let output_records = read_bam_records(&paths.output)?;
        let reject_records = read_bam_records(&paths.rejects)?;

        // About 33 groups should pass (those with 2 reads)
        assert!(output_records.len() >= 30, "Should have at least 30 consensus reads");

        // About 17 groups should be rejected (those with 1 read)
        assert!(reject_records.len() >= 15, "Should have at least 15 rejected reads");

        Ok(())
    }

    #[test]
    fn test_multithreaded_with_stats() -> Result<()> {
        // Test multi-threaded execution with statistics file generation
        let mut builder = SamBuilder::new_unmapped();

        for idx in 0..100 {
            let mut attrs = HashMap::new();
            attrs.insert("MI", BufValue::from(format!("umi_{idx}")));
            builder.add_frag_with_attrs(&format!("read_{idx}a"), None, true, &attrs);
            builder.add_frag_with_attrs(&format!("read_{idx}b"), None, true, &attrs);
        }

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.stats_opts.stats = Some(paths.stats.clone());
        cmd.threading = ThreadingOptions::new(4);

        cmd.execute("test")?;

        // Verify stats file exists and has expected content
        assert!(&paths.stats.exists(), "Stats file should exist");

        // Read and verify TSV format (now vertical key-value-description format)
        let kv_metrics: Vec<ConsensusKvMetric> =
            DelimFile::default().read_tsv(&paths.stats).expect("Failed to read metrics file");

        // Should have multiple rows (one per metric)
        assert!(!kv_metrics.is_empty(), "Should have metrics");

        // Verify raw_reads_considered has a value
        let raw_reads = kv_metrics
            .iter()
            .find(|m| m.key == "raw_reads_considered")
            .expect("Should have raw_reads_considered");
        let count: u64 = raw_reads.value.parse().expect("Should be a number");
        assert!(count > 0, "Stats should have total reads");

        // Verify consensus_reads_emitted has a value
        let consensus = kv_metrics
            .iter()
            .find(|m| m.key == "consensus_reads_emitted")
            .expect("Should have consensus_reads_emitted");
        let count: u64 = consensus.value.parse().expect("Should be a number");
        assert!(count > 0, "Stats should have consensus count");

        Ok(())
    }

    #[test]
    fn test_overlapping_consensus_calling_paired() -> Result<()> {
        // Test overlapping consensus calling on paired-end reads
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("test_umi"));

        // Add 2 pairs with overlapping bases (both are unmapped in this simple test)
        builder.add_pair_with_attrs("pair1", None, None, true, true, &attrs);
        builder.add_pair_with_attrs("pair2", None, None, true, true, &attrs);

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.overlapping.consensus_call_overlapping_bases = true;

        cmd.execute("test")?;

        // Verify output exists (overlapping consensus was called)
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 2, "Should have 2 consensus reads (R1 and R2)");

        Ok(())
    }

    #[test]
    fn test_overlapping_consensus_unpaired_read() -> Result<()> {
        // Test that unpaired reads in overlapping mode are handled correctly
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("test_umi"));

        // Add 3 unpaired fragments
        builder.add_frag_with_attrs("frag1", None, true, &attrs);
        builder.add_frag_with_attrs("frag2", None, true, &attrs);
        builder.add_frag_with_attrs("frag3", None, true, &attrs);

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.overlapping.consensus_call_overlapping_bases = true;

        cmd.execute("test")?;

        // Verify unpaired reads are processed correctly
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read from unpaired fragments");

        Ok(())
    }

    #[test]
    fn test_trim_enabled() -> Result<()> {
        // Test quality trimming option
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("test_umi"));
        builder.add_frag_with_attrs("read1", None, true, &attrs);
        builder.add_frag_with_attrs("read2", None, true, &attrs);

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.consensus.trim = true;

        cmd.execute("test")?;

        // Verify output exists with trimming enabled
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        Ok(())
    }

    #[test]
    fn test_min_consensus_base_quality() -> Result<()> {
        // Test minimum consensus base quality masking
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("test_umi"));
        builder.add_frag_with_attrs("read1", None, true, &attrs);
        builder.add_frag_with_attrs("read2", None, true, &attrs);

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.consensus.min_consensus_base_quality = 30;

        cmd.execute("test")?;

        // Verify output exists with quality filtering
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        Ok(())
    }

    #[test]
    fn test_reads_without_umi_tag_skipped() -> Result<()> {
        // Test that reads without the UMI tag are skipped
        let mut builder = SamBuilder::new_unmapped();

        let mut attrs_with_umi = HashMap::new();
        attrs_with_umi.insert("MI", BufValue::from("has_umi"));
        builder.add_frag_with_attrs("with_umi1", None, true, &attrs_with_umi);
        builder.add_frag_with_attrs("with_umi2", None, true, &attrs_with_umi);

        // Add reads without UMI tag
        let attrs_no_umi = HashMap::new();
        builder.add_frag_with_attrs("no_umi1", None, true, &attrs_no_umi);
        builder.add_frag_with_attrs("no_umi2", None, true, &attrs_no_umi);

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());

        cmd.execute("test")?;

        // Verify only reads with UMI tag generated consensus
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read from reads with UMI only");

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
        let paths = TestPaths::new()?;

        let mut builder = SamBuilder::with_single_ref("chr1", 100);
        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("1"));
        // Add a few reads with the same MI tag
        builder.add_frag_with_attrs("read1", None, true, &attrs);
        builder.add_frag_with_attrs("read2", None, true, &attrs);
        builder.add_frag_with_attrs("read3", None, true, &attrs);
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.threading = threading;
        cmd.execute("test")?;

        // Should produce consensus output
        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        Ok(())
    }

    #[test]
    fn test_simplex_processed_batch_memory_estimate() {
        let mut data = Vec::with_capacity(1024);
        data.extend_from_slice(&[0u8; 100]);
        let mut reject = Vec::with_capacity(512);
        reject.extend_from_slice(&[0u8; 50]);

        let batch = SimplexProcessedBatch {
            consensus_output: ConsensusOutput { data, count: 1 },
            rejects: vec![reject],
            groups_count: 1,
            stats: ConsensusCallingStats::default(),
            overlapping_stats: None,
        };

        let estimate = batch.estimate_heap_size();
        // Should use capacity (1024 + 512) not len (100 + 50)
        assert!(estimate >= 1024 + 512, "estimate {estimate} should be >= 1536 (capacities)");
        // Should also include Vec<Vec<u8>> overhead
        assert!(estimate > 1024 + 512, "estimate {estimate} should include Vec overhead");
    }

    #[rstest]
    #[case::single_threaded(ThreadingOptions::none())]
    #[case::multi_threaded(ThreadingOptions::new(2))]
    fn test_simplex_em_seq_command(#[case] threading: ThreadingOptions) -> Result<()> {
        use crate::sam::builder::{Strand, create_test_fasta};

        // Create a FASTA with chr1 containing C bases at positions 0..20
        let ref_seq = "C".repeat(200);
        let ref_fasta = create_test_fasta(&[("chr1", &ref_seq)])?;

        // Create mapped reads with MI tag showing C bases (methylated) at ref-C positions
        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        let mut attrs = HashMap::new();
        attrs.insert("MI", BufValue::from("1"));

        // Add 3 fragment reads at position 1 (1-based), all showing C (methylated)
        for i in 0..3 {
            let _ = builder
                .add_frag()
                .name(&format!("r{i}"))
                .start(1)
                .strand(Strand::Plus)
                .bases("CCCCCCCCCC")
                .attr("MI", BufValue::from("1"))
                .build();
        }

        let paths = TestPaths::new()?;
        builder.write(&paths.input)?;

        let mut cmd = create_simplex_with_paths(paths.input.clone(), paths.output.clone());
        cmd.methylation_mode = Some(crate::commands::common::MethylationModeArg::EmSeq);
        cmd.reference = Some(ref_fasta.path().to_path_buf());
        cmd.threading = threading;
        cmd.execute("test")?;

        let records = read_bam_records(&paths.output)?;
        assert_eq!(records.len(), 1, "Should have 1 consensus read");

        // Verify methylation tags are present and correct
        use noodles::sam::alignment::record_buf::data::field::value::{
            Array as BufArray, Value as BufValue,
        };
        let record = &records[0];
        let cu_tag = Tag::from([b'c', b'u']);
        let cu_value = record.data().get(&cu_tag).expect("cu tag should be present with EM-Seq");
        // All 3 reads show C at ref-C positions → unconverted counts should be non-zero
        if let BufValue::Array(BufArray::Int16(cu_vals)) = cu_value {
            assert!(
                cu_vals.iter().any(|&v| v > 0),
                "cu should have non-zero values for methylated reads"
            );
        } else {
            panic!("cu tag should be an i16 array");
        }

        let ct_tag = Tag::from([b'c', b't']);
        let ct_value = record.data().get(&ct_tag).expect("ct tag should be present with EM-Seq");
        // All reads show C (not T) → converted counts should be 0
        if let BufValue::Array(BufArray::Int16(ct_vals)) = ct_value {
            assert!(
                ct_vals.iter().all(|&v| v == 0),
                "ct should be all zeros for fully methylated reads"
            );
        } else {
            panic!("ct tag should be an i16 array");
        }

        Ok(())
    }
}
