//! `Clip` command implementation.
//!
//! Clips reads in a BAM file to remove overlapping portions of read pairs.
//! This is useful for variant calling to avoid double-counting evidence from
//! overlapping portions of paired reads.

use crate::alignment_tags::regenerate_alignment_tags_raw;
use crate::clipper::{ClippingMode, RawRecordClipper};
use crate::grouper::TemplateGrouper;
use crate::logging::OperationTimer;
use crate::metrics::clip::{ClipCounts, ClippingMetricsCollection};
use crate::metrics::writer::write_metrics as write_metrics_tsv;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::reference::ReferenceReader;
use crate::sam::SamTag;
use crate::template::{TemplateBatch, TemplateIterator};
use crate::unified_pipeline::{
    GroupKeyConfig, Grouper, MemoryEstimate, run_bam_pipeline_from_reader,
};
use crate::validation::validate_file_exists;
use anyhow::Result;
use clap::Parser;
use fgumi_bam_io::ProgressTracker;
use fgumi_bam_io::{create_bam_reader_for_pipeline, create_raw_bam_reader, create_raw_bam_writer};
use fgumi_raw_bam::RawRecord;
use log::info;
use noodles::sam::Header;
use std::io;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use super::command::Command;
use super::common::{
    BamIoOptions, CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
    build_pipeline_config, parse_bool,
};

/// Clips reads in a BAM file to remove overlaps
#[derive(Parser, Debug)]
#[command(
    name = "clip",
    about = "\x1b[38;5;173m[POST-CONSENSUS]\x1b[0m \x1b[36mClip overlapping reads in BAM files\x1b[0m",
    long_about = r#"
Clips reads from the same template. Ensures that at least N bases are clipped from any end of the read (i.e.
R1 5' end, R1 3' end, R2 5' end, and R2 3' end). Optionally clips reads from the same template to eliminate overlap
between the reads. This ensures that downstream processes, particularly variant calling, cannot double-count
evidence from the same template when both reads span a variant site in the same template.

Clipping overlapping reads is only performed on FR read pairs, and is implemented by clipping approximately half
the overlapping bases from each read. By default hard clipping is performed; the mode may be changed with
--clipping-mode.

Secondary alignments and supplemental alignments are not clipped, but are passed through into the output.

In order to correctly clip reads by template and update mate information, the input BAM must be either
queryname sorted or query grouped. If your input BAM is not in an appropriate order the sort can be
done in streaming fashion with, for example:

  fgumi sort -i in.bam --order queryname | fgumi clip -i /dev/stdin ...

The output is written in the same order as the input. To produce coordinate-sorted output, pipe the
result through a separate `fgumi sort`.

Any existing NM, UQ and MD tags are repaired, and mate-pair information is updated.

Three clipping modes are supported:
1. `soft` - soft-clip the bases and qualities.
2. `soft-with-mask` - soft-clip and mask the bases and qualities (make bases Ns and qualities the minimum).
3. `hard` - hard-clip the bases and qualities.

The --upgrade-clipping parameter will convert all existing clipping in the input to the given more stringent mode:
from `soft` to either `soft-with-mask` or `hard`, and `soft-with-mask` to `hard`. In all other cases, clipping remains
the same prior to applying any other clipping criteria.
"#
)]
#[allow(clippy::struct_excessive_bools)]
pub struct Clip {
    /// Input/output BAM options
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Reference FASTA file (required for tag regeneration)
    #[arg(short = 'r', long = "reference", alias = "ref", required = true)]
    pub reference: PathBuf,

    /// Clipping mode: soft, soft-with-mask, or hard
    #[arg(short = 'c', long = "clipping-mode", default_value_t = ClippingMode::Hard)]
    pub clipping_mode: ClippingMode,

    /// Clip overlapping read pairs
    #[arg(long = "clip-overlapping-reads", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub clip_overlapping_reads: bool,

    /// Clip reads that extend past their mate's start position
    #[arg(
        long = "clip-bases-past-mate",
        alias = "clip-extending-past-mate",
        default_value = "false",
        num_args = 0..=1,
        default_missing_value = "true",
        action = clap::ArgAction::Set,
        value_parser = parse_bool,
    )]
    pub clip_extending_past_mate: bool,

    /// Minimum bases to clip from 5' end of R1
    #[arg(long = "read-one-five-prime", default_value = "0")]
    pub read_one_five_prime: usize,

    /// Minimum bases to clip from 3' end of R1
    #[arg(long = "read-one-three-prime", default_value = "0")]
    pub read_one_three_prime: usize,

    /// Minimum bases to clip from 5' end of R2
    #[arg(long = "read-two-five-prime", default_value = "0")]
    pub read_two_five_prime: usize,

    /// Minimum bases to clip from 3' end of R2
    #[arg(long = "read-two-three-prime", default_value = "0")]
    pub read_two_three_prime: usize,

    /// Upgrade existing clipping to the specified clipping mode
    #[arg(short = 'H', long = "upgrade-clipping", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub upgrade_clipping: bool,

    /// Automatically clip extended attributes that match read length
    #[arg(short = 'a', long = "auto-clip-attributes", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub auto_clip_attributes: bool,

    /// Output file for clipping metrics (only produced by the single-threaded path; cannot be
    /// combined with --threads)
    #[arg(short = 'm', long = "metrics")]
    pub metrics: Option<PathBuf>,

    /// Threading options for parallel processing
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Scheduler and pipeline stats options
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

// ============================================================================
// Types for 7-step pipeline processing
// ============================================================================

/// Result from processing a batch of templates through clipping.
struct ClipProcessedBatch {
    /// Clipped records to write to output BAM.
    clipped_records: Vec<RawRecord>,
    /// Number of templates processed.
    templates_count: u64,
    /// Number of templates with overlap clipping applied.
    overlap_clipped_count: u64,
    /// Number of templates with mate extension clipping applied.
    extend_clipped_count: u64,
}

impl MemoryEstimate for ClipProcessedBatch {
    fn estimate_heap_size(&self) -> usize {
        self.clipped_records.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
            + self.clipped_records.capacity() * std::mem::size_of::<RawRecord>()
    }
}

/// Metrics collected from clipping processing, aggregated post-pipeline.
#[derive(Default)]
struct CollectedClipMetrics {
    /// Total templates processed.
    total_templates: u64,
    /// Templates with overlap clipping.
    overlap_clipped: u64,
    /// Templates with mate extension clipping.
    extend_clipped: u64,
}

impl Command for Clip {
    #[allow(clippy::too_many_lines)]
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate the input exists (stdin paths are exempt).
        self.io.validate()?;
        validate_file_exists(&self.reference, "Reference FASTA")?;

        info!("Clip");
        info!("  Input: {}", self.io.input.display());
        info!("  Output: {}", self.io.output.display());
        info!("  Clipping mode: {}", self.clipping_mode);
        info!("  Clip overlapping reads: {}", self.clip_overlapping_reads);
        info!("  Clip extending past mate: {}", self.clip_extending_past_mate);
        info!("  {}", self.threading.log_message());

        let timer = OperationTimer::new("Clipping reads");

        let mode = self.clipping_mode;

        // Validate clipping parameters
        if self.upgrade_clipping
            || self.clip_overlapping_reads
            || self.clip_extending_past_mate
            || self.read_one_five_prime > 0
            || self.read_one_three_prime > 0
            || self.read_two_five_prime > 0
            || self.read_two_three_prime > 0
        {
            // At least one clipping option is active
        } else {
            anyhow::bail!("At least one clipping option is required");
        }

        // --metrics is not produced in --threads mode; fail fast rather than
        // silently dropping a user-requested output file.
        if self.threading.threads.is_some()
            && let Some(path) = &self.metrics
        {
            anyhow::bail!(
                "--metrics {} cannot be used with --threads: detailed clipping metrics \
                     are only produced by the single-threaded path",
                path.display()
            );
        }

        // ========================================================================
        // CRITICAL: Check --threads mode BEFORE creating any file handles.
        // The 7-step pipeline manages its own I/O and writer lifecycle, so we
        // must not create a writer here if we're going to use the pipeline.
        // Route: Some(threads) -> pipeline, None -> single-threaded fast path
        // ========================================================================
        let total_records = if let Some(threads) = self.threading.threads {
            // Read header for the 7-step pipeline (supports stdin)
            let (reader, header) = create_bam_reader_for_pipeline(&self.io.input)?;

            // Synthesize @HD VN:1.6 SO:unsorted when the input lacks one (match fgbio).
            // Run before require_query_grouped so both execution modes guard the same
            // normalized header (matching execute_single_threaded's ordering).
            let header = crate::commands::common::ensure_hd_record(header)?;

            // CLIP3-05: fgbio's ClipBam calls Bams.requireQueryGrouped. Clipping is
            // template-based (pair clip, overlap, past-mate, mate-fix), so coordinate-
            // sorted input silently mis-groups mates. Guard the *input* header.
            crate::commands::common::require_query_grouped(
                &header,
                &self.io.input.display().to_string(),
            )?;

            // Load reference (always required for clip)
            let reference = Arc::new(ReferenceReader::new(&self.reference)?);

            // Add @PG record with PP chaining to input's last program
            let header = crate::commands::common::add_pg_record(header, command_line)?;

            self.execute_threads_mode(reader, threads, header, reference)?
        } else {
            self.execute_single_threaded(mode, command_line)?
        };

        timer.log_completion(total_records);
        Ok(())
    }
}

impl Clip {
    /// Execute single-threaded clipping.
    #[allow(clippy::too_many_lines)]
    fn execute_single_threaded(&self, mode: ClippingMode, command_line: &str) -> Result<u64> {
        // Create clipper
        let clipper = if self.auto_clip_attributes {
            RawRecordClipper::with_auto_clip(mode, true)
        } else {
            RawRecordClipper::new(mode)
        };

        // Create metrics collection if metrics output requested
        let mut metrics_collection =
            self.metrics.as_ref().map(|_| ClippingMetricsCollection::new());

        // Load reference (always required and tags are always regenerated to match Scala fgbio)
        let reference_reader = ReferenceReader::new(&self.reference)?;

        // Open input BAM with MT BGZF decompression
        let reader_threads = self.threading.num_threads();
        let (reader, header) = create_raw_bam_reader(&self.io.input, reader_threads)?;

        // Synthesize @HD VN:1.6 SO:unsorted when the input lacks one (match fgbio).
        let header = crate::commands::common::ensure_hd_record(header)?;
        // CLIP3-05: require query-grouped input (fgbio Bams.requireQueryGrouped).
        crate::commands::common::require_query_grouped(
            &header,
            &self.io.input.display().to_string(),
        )?;

        // Add @PG record with PP chaining to input's last program
        let header = crate::commands::common::add_pg_record(header, command_line)?;

        // Create output BAM writer with multi-threaded BGZF compression
        let writer_threads = self.threading.num_threads();
        let mut writer = create_raw_bam_writer(
            &self.io.output,
            &header,
            writer_threads,
            self.compression.compression_level,
        )?;

        // Process templates
        info!("Processing templates...");

        let mut total_records: usize = 0;
        let mut total_clipped_overlap: u64 = 0;
        let mut total_clipped_mate_extension: u64 = 0;
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        // Single-threaded processing: iterate raw templates, clip in place
        let template_iter = TemplateIterator::new(reader);

        for template in template_iter {
            let template = template?;
            let mut records: Vec<RawRecord> = template.into_records();

            // Upgrade existing clipping on *every* read of the template first — including
            // secondary/supplementary alignments — matching fgbio ClipBam, which runs
            // `upgradeAllClipping` over `template.allReads` (ClipBam.scala:123) before clipping the
            // primary pair. Doing it per-template here (rather than per-primary inside
            // clip_pair/clip_fragment) is what lets supplementary reads' clipping be upgraded too.
            if self.upgrade_clipping {
                for record in &mut records {
                    clipper.upgrade_all_clipping_raw(record)?;
                }
            }

            // Clip the primary R1/R2 (or a lone fragment) — never positional records[0]/[1] —
            // so templates with secondary/supplementary reads (len > 2) are clipped too, matching
            // fgbio ClipBam's Template.r1/r2 handling.
            match find_primary_pair_indices(&records)? {
                (Some(i1), Some(i2)) => {
                    let [r1, r2] =
                        records.get_disjoint_mut([i1, i2]).expect("distinct primary indices");

                    let (overlap_clip, extend_clip) =
                        self.clip_pair(&clipper, r1, r2, metrics_collection.as_mut())?;

                    if overlap_clip {
                        total_clipped_overlap += 1;
                    }
                    if extend_clip {
                        total_clipped_mate_extension += 1;
                    }

                    // Fix full mate-pair info (mate coords/strand, mate-unmapped flag, MQ/MC,
                    // TLEN), matching fgbio ClipBam's SamPairUtil.setMateInfo call, then repair
                    // mate info on any supplementary alignments in the template.
                    set_mate_info_raw(r1, r2);
                    fix_supplemental_mate_info(&mut records, i1, i2);
                }
                (Some(i1), None) => {
                    self.clip_fragment(&clipper, &mut records[i1], metrics_collection.as_mut())?;
                }
                _ => {}
            }

            // Regenerate alignment tags (always done to match Scala fgbio behavior)
            for record in &mut records {
                regenerate_alignment_tags_raw(record.as_mut_vec(), &header, &reference_reader)?;
            }

            // Count and write records
            let batch_size = records.len();
            total_records += batch_size;
            for record in &records {
                writer.write_raw_record(record.as_ref())?;
            }
            progress.log_if_needed(batch_size as u64);
        }

        progress.log_final();
        info!("Total records processed: {total_records}");
        info!("Templates with overlap clipping: {total_clipped_overlap}");
        info!("Templates with mate extension clipping: {total_clipped_mate_extension}");

        // Write metrics if requested
        if let Some(metrics_path) = &self.metrics
            && let Some(mut metrics) = metrics_collection
        {
            info!("Writing metrics to {}", metrics_path.display());
            metrics.finalize();
            let all_metrics = metrics.all_metrics().map(Clone::clone);
            write_metrics_tsv(metrics_path, &all_metrics, "clipping")?;
            info!("Metrics written successfully");
        }

        // Flush and finish the writer
        writer.finish()?;

        info!("Done!");
        Ok(total_records as u64)
    }

    /// Clips a fragment (unpaired) read.
    ///
    /// Applies clipping operations to a single fragment read, including:
    /// 1. Applying fixed-position 5' and 3' clipping
    /// 2. Updating metrics if a metrics collector is provided
    ///
    /// Clipping upgrades (`--upgrade-clipping`) are *not* performed here: the caller runs a
    /// template-wide pre-pass that upgrades clipping on every read of the template (including
    /// secondary/supplementary alignments) before this method runs, matching fgbio `ClipBam`
    /// (`ClipBam.scala:123`).
    ///
    /// Fragment reads are treated as R1 for the purposes of fixed-position clipping.
    ///
    /// # Arguments
    ///
    /// * `clipper` - The record clipper instance
    /// * `record` - The fragment record to clip (mutable)
    /// * `metrics` - Optional metrics collector to update
    ///
    /// # Returns
    ///
    /// `Ok(())` on success.
    ///
    /// # Errors
    ///
    /// Returns an error if clipping operations fail.
    fn clip_fragment(
        &self,
        clipper: &RawRecordClipper,
        record: &mut RawRecord,
        metrics: Option<&mut ClippingMetricsCollection>,
    ) -> Result<()> {
        let prior_bases_clipped = clipped_bases_raw(record);

        // Note: clipping upgrades (`--upgrade-clipping`) are applied once per template over all
        // reads by the caller (matching fgbio ClipBam.scala:123), not here.

        // Apply fixed-position clipping
        let num_five_prime = if self.read_one_five_prime > 0 {
            clipper.clip_5_prime_end_of_read_raw(record, self.read_one_five_prime)
        } else {
            0
        };

        let num_three_prime = if self.read_one_three_prime > 0 {
            clipper.clip_3_prime_end_of_read_raw(record, self.read_one_three_prime)
        } else {
            0
        };

        // Update metrics
        if let Some(metrics) = metrics {
            metrics.fragment.update_raw(
                record,
                ClipCounts {
                    prior: prior_bases_clipped,
                    five_prime: num_five_prime,
                    three_prime: num_three_prime,
                    ..ClipCounts::default()
                },
            );
        }

        Ok(())
    }

    /// Clips a pair of reads with comprehensive clipping logic.
    ///
    /// Applies multiple types of clipping to a read pair:
    /// 1. Fixed-position 5' and 3' clipping for each read
    /// 2. Overlap clipping to remove duplicate coverage
    /// 3. Mate-extension clipping to remove reads extending past mate start
    /// 4. Updating metrics for both reads
    ///
    /// Clipping upgrades (`--upgrade-clipping`) are *not* performed here: the caller runs a
    /// template-wide pre-pass that upgrades clipping on every read of the template (including
    /// secondary/supplementary alignments) before this method runs, matching fgbio `ClipBam`
    /// (`ClipBam.scala:123`).
    ///
    /// The method intelligently determines which read is R1 vs R2 based on SAM flags
    /// and applies the appropriate fixed-position clipping thresholds.
    ///
    /// # Arguments
    ///
    /// * `clipper` - The record clipper instance
    /// * `r1` - First read of the pair (mutable)
    /// * `r2` - Second read of the pair (mutable)
    /// * `metrics` - Optional metrics collector to update
    ///
    /// # Returns
    ///
    /// A tuple of `(overlap_clipped, extend_clipped)` booleans indicating whether
    /// overlap or mate-extension clipping was performed.
    ///
    /// # Errors
    ///
    /// Returns an error if clipping operations fail.
    fn clip_pair(
        &self,
        clipper: &RawRecordClipper,
        r1: &mut RawRecord,
        r2: &mut RawRecord,
        metrics: Option<&mut ClippingMetricsCollection>,
    ) -> Result<(bool, bool)> {
        let prior_bases_clipped_r1 = clipped_bases_raw(r1);
        let prior_bases_clipped_r2 = clipped_bases_raw(r2);

        // Note: clipping upgrades (`--upgrade-clipping`) are applied once per template over all
        // reads by the caller (matching fgbio ClipBam.scala:123), not here.

        // Determine read types (raw flags: bit 6 = first segment, bit 7 = last segment)
        let (is_r1_first, is_r2_last) = (r1.is_first_segment(), r2.is_last_segment());

        // Apply fixed-position clipping for R1
        let num_r1_five_prime = if is_r1_first && self.read_one_five_prime > 0 {
            clipper.clip_5_prime_end_of_read_raw(r1, self.read_one_five_prime)
        } else if !is_r1_first && self.read_two_five_prime > 0 {
            clipper.clip_5_prime_end_of_read_raw(r1, self.read_two_five_prime)
        } else {
            0
        };

        let num_r1_three_prime = if is_r1_first && self.read_one_three_prime > 0 {
            clipper.clip_3_prime_end_of_read_raw(r1, self.read_one_three_prime)
        } else if !is_r1_first && self.read_two_three_prime > 0 {
            clipper.clip_3_prime_end_of_read_raw(r1, self.read_two_three_prime)
        } else {
            0
        };

        // Apply fixed-position clipping for R2
        let num_r2_five_prime = if is_r2_last && self.read_two_five_prime > 0 {
            clipper.clip_5_prime_end_of_read_raw(r2, self.read_two_five_prime)
        } else if !is_r2_last && self.read_one_five_prime > 0 {
            clipper.clip_5_prime_end_of_read_raw(r2, self.read_one_five_prime)
        } else {
            0
        };

        let num_r2_three_prime = if is_r2_last && self.read_two_three_prime > 0 {
            clipper.clip_3_prime_end_of_read_raw(r2, self.read_two_three_prime)
        } else if !is_r2_last && self.read_one_three_prime > 0 {
            clipper.clip_3_prime_end_of_read_raw(r2, self.read_one_three_prime)
        } else {
            0
        };

        // Clip overlapping reads
        let (num_overlapping_r1, num_overlapping_r2) = if self.clip_overlapping_reads {
            clipper.clip_overlapping_reads(r1, r2)
        } else {
            (0, 0)
        };

        // Clip reads extending past mate
        let (num_extending_r1, num_extending_r2) = if self.clip_extending_past_mate {
            clipper.clip_extending_past_mate_ends(r1, r2)
        } else {
            (0, 0)
        };

        // Update metrics
        if let Some(metrics) = metrics {
            let r1_counts = ClipCounts {
                prior: prior_bases_clipped_r1,
                five_prime: num_r1_five_prime,
                three_prime: num_r1_three_prime,
                overlapping: num_overlapping_r1,
                extending: num_extending_r1,
            };
            let r2_counts = ClipCounts {
                prior: prior_bases_clipped_r2,
                five_prime: num_r2_five_prime,
                three_prime: num_r2_three_prime,
                overlapping: num_overlapping_r2,
                extending: num_extending_r2,
            };

            // Determine which metric to update based on read flags
            if is_r1_first {
                metrics.read_one.update_raw(r1, r1_counts);
            } else {
                metrics.read_two.update_raw(r1, r1_counts);
            }

            if is_r2_last {
                metrics.read_two.update_raw(r2, r2_counts);
            } else {
                metrics.read_one.update_raw(r2, r2_counts);
            }
        }

        let overlap_clipped = num_overlapping_r1 > 0 || num_overlapping_r2 > 0;
        let extend_clipped = num_extending_r1 > 0 || num_extending_r2 > 0;

        Ok((overlap_clipped, extend_clipped))
    }

    /// Execute using the 7-step unified pipeline with multi-threading.
    ///
    /// Uses `TemplateGrouper` to batch records by template (QNAME) for parallel processing.
    #[allow(clippy::too_many_lines)]
    fn execute_threads_mode(
        &self,
        reader: Box<dyn std::io::Read + Send>,
        num_threads: usize,
        header: Header,
        reference: Arc<ReferenceReader>,
    ) -> Result<u64> {
        // Configure pipeline - clip is writer-heavy workload
        let mut pipeline_config = build_pipeline_config(
            &self.scheduler_opts,
            &self.compression,
            &self.queue_memory,
            num_threads,
        )?;
        // Clip uses raw-byte mode so TemplateGrouper receives RawRecord items.
        pipeline_config.group_key_config =
            Some(GroupKeyConfig::new_raw_no_cell(crate::read_info::LibraryIndex::default()));

        // Per-thread metrics accumulator: bounded memory, no unbounded queue.
        let collected_metrics = PerThreadAccumulator::<CollectedClipMetrics>::new(num_threads);
        let collected_for_serialize = Arc::clone(&collected_metrics);

        // Configuration for closures
        const BATCH_SIZE: usize = 1000;
        let clipping_mode = self.clipping_mode;
        let auto_clip_attributes = self.auto_clip_attributes;
        let upgrade_clipping = self.upgrade_clipping;
        let clip_overlapping_reads = self.clip_overlapping_reads;
        let clip_extending_past_mate = self.clip_extending_past_mate;
        let read_one_five_prime = self.read_one_five_prime;
        let read_one_three_prime = self.read_one_three_prime;
        let read_two_five_prime = self.read_two_five_prime;
        let read_two_three_prime = self.read_two_three_prime;
        let reference_for_process = Arc::clone(&reference);
        let header_for_process = header.clone();

        // Progress tracking
        let progress_counter = Arc::new(AtomicU64::new(0));
        let progress_for_process = Arc::clone(&progress_counter);

        // Grouper: batch records by template (QNAME)
        let grouper_fn = move |_header: &Header| {
            Box::new(TemplateGrouper::new(BATCH_SIZE))
                as Box<dyn Grouper<Group = TemplateBatch> + Send>
        };

        // Process function: clip each template batch
        let process_fn = move |batch: TemplateBatch| -> io::Result<ClipProcessedBatch> {
            // Create per-worker clipper
            let clipper = if auto_clip_attributes {
                RawRecordClipper::with_auto_clip(clipping_mode, true)
            } else {
                RawRecordClipper::new(clipping_mode)
            };

            let mut clipped_records = Vec::new();
            let mut templates_count = 0u64;
            let mut overlap_clipped_count = 0u64;
            let mut extend_clipped_count = 0u64;

            for template in batch {
                templates_count += 1;
                let mut records: Vec<RawRecord> = template.into_records();

                // Upgrade existing clipping on *every* read of the template first — including
                // secondary/supplementary alignments — matching fgbio ClipBam's
                // `template.allReads.foreach(upgradeAllClipping)` (ClipBam.scala:123) before the
                // primary pair is clipped.
                if upgrade_clipping {
                    for record in &mut records {
                        clipper.upgrade_all_clipping_raw(record).map_err(io::Error::other)?;
                    }
                }

                // Clip the primary R1/R2 (or a lone fragment) by flags, never positional
                // records[0]/[1], so templates with secondary/supplementary reads are clipped
                // too (matches fgbio ClipBam's Template.r1/r2 handling).
                match find_primary_pair_indices(&records).map_err(io::Error::other)? {
                    (Some(i1), Some(i2)) => {
                        let [r1, r2] =
                            records.get_disjoint_mut([i1, i2]).expect("distinct primary indices");

                        // Determine read types (raw flags)
                        let is_r1_first = r1.is_first_segment();

                        // Apply fixed-position clipping for R1
                        if is_r1_first && read_one_five_prime > 0 {
                            clipper.clip_5_prime_end_of_read_raw(r1, read_one_five_prime);
                        } else if !is_r1_first && read_two_five_prime > 0 {
                            clipper.clip_5_prime_end_of_read_raw(r1, read_two_five_prime);
                        }
                        if is_r1_first && read_one_three_prime > 0 {
                            clipper.clip_3_prime_end_of_read_raw(r1, read_one_three_prime);
                        } else if !is_r1_first && read_two_three_prime > 0 {
                            clipper.clip_3_prime_end_of_read_raw(r1, read_two_three_prime);
                        }

                        // Apply fixed-position clipping for R2. The primary R2 slot is always a
                        // last-segment read (find_primary_pair_indices only fills it from an
                        // is_last_segment read), so it always uses the read-two thresholds.
                        if read_two_five_prime > 0 {
                            clipper.clip_5_prime_end_of_read_raw(r2, read_two_five_prime);
                        }
                        if read_two_three_prime > 0 {
                            clipper.clip_3_prime_end_of_read_raw(r2, read_two_three_prime);
                        }

                        // Clip overlapping reads
                        if clip_overlapping_reads {
                            let (num_r1, num_r2) = clipper.clip_overlapping_reads(r1, r2);
                            if num_r1 > 0 || num_r2 > 0 {
                                overlap_clipped_count += 1;
                            }
                        }

                        // Clip reads extending past mate
                        if clip_extending_past_mate {
                            let (num_r1, num_r2) = clipper.clip_extending_past_mate_ends(r1, r2);
                            if num_r1 > 0 || num_r2 > 0 {
                                extend_clipped_count += 1;
                            }
                        }

                        // Fix full mate-pair info (mate coords/strand, mate-unmapped flag,
                        // MQ/MC, TLEN), matching fgbio ClipBam's SamPairUtil.setMateInfo call,
                        // then repair mate info on any supplementary alignments.
                        set_mate_info_raw(r1, r2);
                        fix_supplemental_mate_info(&mut records, i1, i2);
                    }
                    (Some(i1), None) => {
                        // Fragment - apply fixed-position clipping (R1 thresholds). Any clipping
                        // upgrade was already applied over all template reads above.
                        let record = &mut records[i1];
                        if read_one_five_prime > 0 {
                            clipper.clip_5_prime_end_of_read_raw(record, read_one_five_prime);
                        }
                        if read_one_three_prime > 0 {
                            clipper.clip_3_prime_end_of_read_raw(record, read_one_three_prime);
                        }
                    }
                    _ => {}
                }

                // Regenerate alignment tags (always done to match fgbio behavior)
                for record in &mut records {
                    regenerate_alignment_tags_raw(
                        record.as_mut_vec(),
                        &header_for_process,
                        &reference_for_process,
                    )
                    .map_err(io::Error::other)?;
                }

                clipped_records.extend(records);
            }

            // Progress logging (count records, not templates)
            let records_count = clipped_records.len() as u64;
            let count = progress_for_process.fetch_add(records_count, Ordering::Relaxed);
            if (count + records_count) / 1_000_000 > count / 1_000_000 {
                info!("Processed {} records", count + records_count);
            }

            Ok(ClipProcessedBatch {
                clipped_records,
                templates_count,
                overlap_clipped_count,
                extend_clipped_count,
            })
        };

        // Serialize function: write raw records directly into the output buffer
        let serialize_fn = move |processed: ClipProcessedBatch,
                                 _header: &Header,
                                 output: &mut Vec<u8>|
              -> io::Result<u64> {
            // Merge per-batch counts into this worker's accumulator slot
            collected_for_serialize.with_slot(|m| {
                m.total_templates += processed.templates_count;
                m.overlap_clipped += processed.overlap_clipped_count;
                m.extend_clipped += processed.extend_clipped_count;
            });

            // Write raw records directly (4-byte block_size + bytes per record)
            let count = processed.clipped_records.len() as u64;
            fgumi_raw_bam::write_raw_records(output, &processed.clipped_records)?;
            Ok(count)
        };

        // Run the 7-step pipeline
        let records_written = run_bam_pipeline_from_reader(
            pipeline_config,
            reader,
            header.clone(),
            &self.io.output,
            None, // Use input header for output
            grouper_fn,
            process_fn,
            serialize_fn,
        )?;

        // ========== Post-pipeline: Aggregate metrics ==========
        let mut total_templates = 0u64;
        let mut total_overlap_clipped = 0u64;
        let mut total_extend_clipped = 0u64;

        for slot in collected_metrics.slots() {
            let m = slot.lock();
            total_templates += m.total_templates;
            total_overlap_clipped += m.overlap_clipped;
            total_extend_clipped += m.extend_clipped;
        }

        info!("Total templates processed: {total_templates}");
        info!("Templates with overlap clipping: {total_overlap_clipped}");
        info!("Templates with mate extension clipping: {total_extend_clipped}");

        // `--metrics` with `--threads` is rejected in `execute()` before we ever
        // reach this path, so no per-pipeline metrics file is written here.

        info!("Done!");
        Ok(records_written)
    }
}

/// Snapshot of the mate-relevant fields of a read, taken before any mutation.
struct MateSnap {
    ref_id: i32,
    pos: i32,
    neg: bool,
    unmapped: bool,
    mapq: u8,
    cigar: String,
    aln_start: Option<usize>,
    aln_end: Option<usize>,
}

impl MateSnap {
    fn of(rec: &RawRecord) -> Self {
        use fgumi_raw_bam::flags as rflags;
        Self {
            ref_id: rec.ref_id(),
            pos: rec.pos(),
            neg: rec.flags() & rflags::REVERSE != 0,
            unmapped: rec.flags() & rflags::UNMAPPED != 0,
            mapq: rec.mapq(),
            cigar: rec.cigar_to_string(),
            aln_start: rec.alignment_start_1based(),
            aln_end: rec.alignment_end_1based(),
        }
    }
}

/// Sets or clears the `MATE_REVERSE` and `MATE_UNMAPPED` flags on a record.
fn set_mate_flags_raw(rec: &mut RawRecord, mate_neg: bool, mate_unmapped: bool) {
    use fgumi_raw_bam::flags as rflags;
    let mut f = rec.flags();
    if mate_neg {
        f |= rflags::MATE_REVERSE;
    } else {
        f &= !rflags::MATE_REVERSE;
    }
    if mate_unmapped {
        f |= rflags::MATE_UNMAPPED;
    } else {
        f &= !rflags::MATE_UNMAPPED;
    }
    rec.set_flags(f);
}

/// Writes the MQ (mate mapping quality) and MC (mate CIGAR) tags for a read.
fn set_mate_mq_mc_raw(rec: &mut RawRecord, mate_mapq: u8, mate_cigar: &str) {
    let mut editor = rec.tags_editor();
    editor.update_int(SamTag::MQ, i32::from(mate_mapq));
    editor.update_string(SamTag::MC, mate_cigar.as_bytes());
}

/// Removes the MQ and MC tags from a read (used when the mate is unmapped).
fn clear_mate_mq_mc_raw(rec: &mut RawRecord) {
    let mut editor = rec.tags_editor();
    editor.remove(SamTag::MQ);
    editor.remove(SamTag::MC);
}

/// Computes the TLEN (inferred insert size) for the first read of a pair, mirroring
/// htsjdk `SamPairUtil.computeInsertSize`. Returns 0 unless both reads are mapped to
/// the same reference; the second read's TLEN is the negation of this value.
fn compute_insert_size_raw(s1: &MateSnap, s2: &MateSnap) -> i32 {
    if s1.unmapped || s2.unmapped || s1.ref_id != s2.ref_id {
        return 0;
    }
    let five_prime = |s: &MateSnap| if s.neg { s.aln_end } else { s.aln_start };
    let (Some(p1), Some(p2)) = (five_prime(s1), five_prime(s2)) else {
        return 0;
    };
    let (p1, p2) = (p1 as i64, p2 as i64);
    let adjustment = if p2 >= p1 { 1 } else { -1 };
    i32::try_from(p2 - p1 + adjustment).unwrap_or(0)
}

/// Sets full mate-pair information on a read pair, mirroring htsjdk
/// `SamPairUtil.setMateInfo(rec1, rec2, setMateCigar=true)`.
///
/// fgbio `ClipBam` calls `SamPairUtil.setMateInfo` on every pair after clipping so that
/// mate reference/position/strand, the mate-unmapped flag, the MQ/MC tags, and TLEN all
/// reflect the post-clip state — including reads that clipping unmapped (see CLIP-01). The
/// three branches (both mapped, both unmapped, one of each) match htsjdk exactly; an
/// unmapped read is relocated to its mapped mate's coordinate.
fn set_mate_info_raw(r1: &mut RawRecord, r2: &mut RawRecord) {
    // Snapshot both reads up front so mutating one never reads stale/updated fields of the other.
    let s1 = MateSnap::of(r1);
    let s2 = MateSnap::of(r2);

    if !s1.unmapped && !s2.unmapped {
        // Both mapped: copy each read's coordinates into the other's mate fields.
        r1.set_mate_ref_id(s2.ref_id);
        r1.set_mate_pos(s2.pos);
        set_mate_flags_raw(r1, s2.neg, false);
        set_mate_mq_mc_raw(r1, s2.mapq, &s2.cigar);

        r2.set_mate_ref_id(s1.ref_id);
        r2.set_mate_pos(s1.pos);
        set_mate_flags_raw(r2, s1.neg, false);
        set_mate_mq_mc_raw(r2, s1.mapq, &s1.cigar);

        let insert_size = compute_insert_size_raw(&s1, &s2);
        r1.set_template_length(insert_size);
        r2.set_template_length(-insert_size);
    } else if s1.unmapped && s2.unmapped {
        // Both unmapped: clear coordinates and mate coordinates, flag each mate unmapped.
        for (rec, mate_neg) in [(&mut *r1, s2.neg), (&mut *r2, s1.neg)] {
            rec.set_ref_id(-1);
            rec.set_pos(-1);
            rec.set_mate_ref_id(-1);
            rec.set_mate_pos(-1);
            set_mate_flags_raw(rec, mate_neg, true);
            clear_mate_mq_mc_raw(rec);
            rec.set_template_length(0);
        }
    } else {
        // Exactly one is unmapped: relocate it to the mapped mate's coordinate.
        let (mapped, unmapped, mapped_snap, unmapped_snap) = if s1.unmapped {
            (&mut *r2, &mut *r1, &s2, &s1)
        } else {
            (&mut *r1, &mut *r2, &s1, &s2)
        };

        // The unmapped read is placed at the mapped read's coordinate.
        unmapped.set_ref_id(mapped_snap.ref_id);
        unmapped.set_pos(mapped_snap.pos);
        unmapped.set_mate_ref_id(mapped_snap.ref_id);
        unmapped.set_mate_pos(mapped_snap.pos);
        set_mate_flags_raw(unmapped, mapped_snap.neg, false);
        set_mate_mq_mc_raw(unmapped, mapped_snap.mapq, &mapped_snap.cigar);
        unmapped.set_template_length(0);

        // The mapped read points its mate fields at the (now co-located) unmapped read.
        // Its mate-reverse flag reflects the unmapped read's *actual* strand (htsjdk
        // `SamPairUtil.java:267`): the clipper's own unmap clears REVERSE, but a read that
        // arrived unmapped-on-input may still carry it, and this runs on every pair.
        mapped.set_mate_ref_id(mapped_snap.ref_id);
        mapped.set_mate_pos(mapped_snap.pos);
        set_mate_flags_raw(mapped, unmapped_snap.neg, true);
        clear_mate_mq_mc_raw(mapped);
        mapped.set_template_length(0);
    }
}

/// Ports htsjdk 5.0.0 `SamPairUtil.setMateInformationOnSupplementalAlignment(supp, matePrimary,
/// setMateCigar=true)`.
///
/// fgbio `ClipBam` calls this on every supplementary alignment after clipping the primary pair,
/// so a supplemental's mate fields point at its mate *primary* read. It sets mate ref/pos/strand,
/// the mate-unmapped flag, TLEN (the negation of the mate primary's already-updated TLEN), the
/// mate CIGAR (MC, only when the mate is mapped) and — as of htsjdk 5.0.0 — the mate mapping
/// quality (MQ, unconditionally). `mate` and `mate_tlen` are snapshotted from the post-clip
/// primary so the caller can fix several supplementals without re-borrowing the primaries.
fn set_supplemental_mate_info_raw(supp: &mut RawRecord, mate: &MateSnap, mate_tlen: i32) {
    supp.set_mate_ref_id(mate.ref_id);
    supp.set_mate_pos(mate.pos);
    set_mate_flags_raw(supp, mate.neg, mate.unmapped);
    supp.set_template_length(-mate_tlen);
    let mut editor = supp.tags_editor();
    if mate.unmapped {
        editor.remove(SamTag::MC);
    } else {
        editor.update_string(SamTag::MC, mate.cigar.as_bytes());
    }
    // htsjdk 5.0.0 sets MQ unconditionally from the mate primary's mapping quality.
    editor.update_int(SamTag::MQ, i32::from(mate.mapq));
}

/// Finds the primary R1 and R2 record indices in a template, following fgbio `Template.r1`/`r2`:
/// the first non-secondary, non-supplementary read that is unpaired or first-of-pair, and the
/// first that is paired and second-of-pair, respectively. Either may be `None`.
///
/// Like fgbio's `Bams.Template` (`Bams.scala:161,167`), a template carrying more than one primary
/// (non-secondary, non-supplementary) R1 — or more than one primary R2 — is malformed and rejected
/// with an error rather than silently keeping the first and passing the extra through unclipped and
/// with unrepaired mate info. The error message mirrors fgbio verbatim.
fn find_primary_pair_indices(records: &[RawRecord]) -> Result<(Option<usize>, Option<usize>)> {
    let mut r1_idx = None;
    let mut r2_idx = None;
    for (i, rec) in records.iter().enumerate() {
        if rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        if !rec.is_paired() || rec.is_first_segment() {
            if r1_idx.is_some() {
                anyhow::bail!(
                    "Multiple non-secondary, non-supplemental R1s for {}",
                    String::from_utf8_lossy(rec.read_name()).trim_end_matches('\0')
                );
            }
            r1_idx = Some(i);
        } else if rec.is_last_segment() {
            if r2_idx.is_some() {
                anyhow::bail!(
                    "Multiple non-secondary, non-supplemental R2s for {}",
                    String::from_utf8_lossy(rec.read_name()).trim_end_matches('\0')
                );
            }
            r2_idx = Some(i);
        }
    }
    Ok((r1_idx, r2_idx))
}

/// Repairs mate information on the template's supplementary alignments after the primary pair
/// (at `r1_idx`/`r2_idx`) has been clipped and had its own mate info set. Mirrors fgbio
/// `ClipBam`: R1 supplementals point at the primary R2, R2 supplementals at the primary R1.
fn fix_supplemental_mate_info(records: &mut [RawRecord], r1_idx: usize, r2_idx: usize) {
    // Snapshot the post-clip primaries so the per-supplemental updates don't re-borrow them.
    let r1_snap = MateSnap::of(&records[r1_idx]);
    let r1_tlen = records[r1_idx].template_length();
    let r2_snap = MateSnap::of(&records[r2_idx]);
    let r2_tlen = records[r2_idx].template_length();

    for rec in records.iter_mut() {
        if !rec.is_supplementary() {
            continue;
        }
        // R1 supplementals (unpaired or first-of-pair) take R2 as their mate; R2 supplementals R1.
        if !rec.is_paired() || rec.is_first_segment() {
            set_supplemental_mate_info_raw(rec, &r2_snap, r2_tlen);
        } else if rec.is_last_segment() {
            set_supplemental_mate_info_raw(rec, &r1_snap, r1_tlen);
        }
    }
}

/// Returns the number of clipped bases (soft + hard) in a raw record's CIGAR.
fn clipped_bases_raw(record: &RawRecord) -> usize {
    record
        .cigar_ops_typed()
        .filter(|op| {
            matches!(
                op.kind(),
                fgumi_raw_bam::CigarKind::SoftClip | fgumi_raw_bam::CigarKind::HardClip
            )
        })
        .map(|op| op.len() as usize)
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::path::PathBuf;

    // R2-CLIP-02: the primary R1/R2 are located by SAM flags, not by position, so templates
    // that carry secondary/supplementary alignments (len > 2) still resolve their primary pair.
    #[test]
    fn test_find_primary_pair_indices_ignores_secondary_and_supplementary() {
        use crate::sam::RecordBuilder;
        use fgumi_raw_bam::encode_record_buf_to_raw;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100_000).expect("ref length must be nonzero"),
        );
        let header =
            noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build();
        let enc = |b: &noodles::sam::alignment::RecordBuf| {
            encode_record_buf_to_raw(b, &header).expect("encode")
        };
        let mapped = |first: bool, secondary: bool, supplementary: bool, start: usize| {
            RecordBuilder::mapped_read()
                .name("q")
                .paired(true)
                .first_segment(first)
                .secondary(secondary)
                .supplementary(supplementary)
                .reference_sequence_id(0)
                .alignment_start(start)
                .cigar("50M")
                .sequence(&"A".repeat(50))
                .build()
        };

        let r1 = mapped(true, false, false, 100);
        let r2 = mapped(false, false, false, 300);
        let supp_r1 = mapped(true, false, true, 700);
        let sec_r2 = mapped(false, true, false, 900);

        // Primaries at positions 0/1 with trailing secondary + supplementary reads.
        let recs = vec![enc(&r1), enc(&r2), enc(&supp_r1), enc(&sec_r2)];
        assert_eq!(find_primary_pair_indices(&recs).unwrap(), (Some(0), Some(1)));

        // Order-independent: a supplementary read first must not be mistaken for a primary.
        let recs2 = vec![enc(&supp_r1), enc(&r2), enc(&r1)];
        assert_eq!(find_primary_pair_indices(&recs2).unwrap(), (Some(2), Some(1)));

        // A lone fragment (unpaired primary) resolves R1 only.
        let frag = RecordBuilder::mapped_read()
            .name("f")
            .paired(false)
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("50M")
            .sequence(&"A".repeat(50))
            .build();
        assert_eq!(find_primary_pair_indices(&[enc(&frag)]).unwrap(), (Some(0), None));
    }

    // A malformed template with two primary (non-secondary, non-supplementary) R1s — or two
    // primary R2s — is rejected loudly, matching fgbio `Bams.Template` (`Bams.scala:161,167`)
    // rather than silently keeping the first and passing the extra through unclipped.
    #[test]
    fn test_find_primary_pair_indices_rejects_duplicate_primaries() {
        use crate::sam::RecordBuilder;
        use fgumi_raw_bam::encode_record_buf_to_raw;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100_000).expect("ref length must be nonzero"),
        );
        let header =
            noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build();
        let enc = |b: &noodles::sam::alignment::RecordBuf| {
            encode_record_buf_to_raw(b, &header).expect("encode")
        };
        let mapped = |first: bool, start: usize| {
            RecordBuilder::mapped_read()
                .name("dup")
                .paired(true)
                .first_segment(first)
                .reference_sequence_id(0)
                .alignment_start(start)
                .cigar("50M")
                .sequence(&"A".repeat(50))
                .build()
        };

        // Two primary R1s (both first-of-pair, neither secondary/supplementary).
        let two_r1 =
            vec![enc(&mapped(true, 100)), enc(&mapped(false, 300)), enc(&mapped(true, 500))];
        let err = find_primary_pair_indices(&two_r1).unwrap_err().to_string();
        assert_eq!(err, "Multiple non-secondary, non-supplemental R1s for dup");

        // Two primary R2s (both last-of-pair).
        let two_r2 =
            vec![enc(&mapped(true, 100)), enc(&mapped(false, 300)), enc(&mapped(false, 500))];
        let err = find_primary_pair_indices(&two_r2).unwrap_err().to_string();
        assert_eq!(err, "Multiple non-secondary, non-supplemental R2s for dup");
    }

    /// Encodes a `RecordBuf` to a `RawRecord` with a shared single-contig header. Used by the
    /// raw mate-info tests below.
    fn encode_raw(rec: &noodles::sam::alignment::RecordBuf) -> RawRecord {
        use fgumi_raw_bam::encode_record_buf_to_raw;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100_000).expect("ref length must be nonzero"),
        );
        let header =
            noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build();
        encode_record_buf_to_raw(rec, &header).expect("encode")
    }

    /// Builds a mapped raw record from the common fields the mate-info tests need.
    fn raw_read(
        first: bool,
        supplementary: bool,
        reverse: bool,
        start: usize,
        mapq: u8,
    ) -> RawRecord {
        use crate::sam::RecordBuilder;
        encode_raw(
            &RecordBuilder::mapped_read()
                .name("q")
                .paired(true)
                .first_segment(first)
                .supplementary(supplementary)
                .reverse_complement(reverse)
                .reference_sequence_id(0)
                .alignment_start(start)
                .mapping_quality(mapq)
                .cigar("50M")
                .sequence(&"A".repeat(50))
                .build(),
        )
    }

    // set_supplemental_mate_info_raw copies a mapped mate's coordinate/strand/MAPQ onto a
    // supplementary read, writes MC from the mate CIGAR, sets MQ, and negates the mate TLEN.
    #[test]
    fn test_set_supplemental_mate_info_raw_mapped_mate() {
        use fgumi_raw_bam::flags as rflags;

        // A reverse-strand primary mate at 1-based 301 (0-based 300), MAPQ 40, CIGAR 50M.
        let mate = raw_read(false, false, true, 301, 40);
        let mut supp = raw_read(true, true, false, 700, 30);

        set_supplemental_mate_info_raw(&mut supp, &MateSnap::of(&mate), 120);

        assert_eq!(supp.mate_ref_id(), 0);
        assert_eq!(supp.mate_pos(), 300);
        assert_eq!(supp.template_length(), -120);
        assert_ne!(supp.flags() & rflags::MATE_REVERSE, 0, "mate is reverse");
        assert_eq!(supp.flags() & rflags::MATE_UNMAPPED, 0, "mate is mapped");
        assert_eq!(supp.tags().find_mc(), Some("50M"));
        assert_eq!(supp.tags().find_int(SamTag::MQ), Some(40));
    }

    // With an unmapped mate, set_supplemental_mate_info_raw flags the mate unmapped and drops MC,
    // but still writes MQ from the mate's mapping quality (htsjdk 5.0.0 sets MQ unconditionally).
    #[test]
    fn test_set_supplemental_mate_info_raw_unmapped_mate() {
        use crate::sam::RecordBuilder;
        use fgumi_raw_bam::flags as rflags;

        // Give the unmapped mate a distinct, non-zero MAPQ so the MQ assertion below proves the
        // value was copied from the mate rather than defaulting to 0.
        let mate = encode_raw(
            &RecordBuilder::mapped_read()
                .name("q")
                .paired(true)
                .first_segment(false)
                .unmapped(true)
                .reference_sequence_id(0)
                .alignment_start(301)
                .mapping_quality(37)
                .cigar("50M")
                .sequence(&"A".repeat(50))
                .build(),
        );
        // Seed an MC tag so we can confirm it is removed when the mate is unmapped.
        let mut supp = raw_read(true, true, false, 700, 30);
        supp.tags_editor().update_string(SamTag::MC, b"10M");
        assert!(supp.tags().contains(SamTag::MC), "MC present before");

        set_supplemental_mate_info_raw(&mut supp, &MateSnap::of(&mate), 0);

        assert_ne!(supp.flags() & rflags::MATE_UNMAPPED, 0, "mate is unmapped");
        assert!(!supp.tags().contains(SamTag::MC), "MC dropped when mate unmapped");
        // MQ is set unconditionally to the mate's mapping quality, even for an unmapped mate.
        assert_eq!(supp.tags().find_int(SamTag::MQ), Some(37), "MQ set from unmapped mate MAPQ");
    }

    // fix_supplemental_mate_info points R1 supplementals at the primary R2 and R2 supplementals
    // at the primary R1, inheriting each primary's coordinate/strand/MAPQ and negated TLEN.
    #[test]
    fn test_fix_supplemental_mate_info() {
        use fgumi_raw_bam::flags as rflags;

        let mut recs = vec![
            raw_read(true, false, false, 101, 60), // 0: primary R1 (forward, MAPQ 60)
            raw_read(false, false, true, 301, 40), // 1: primary R2 (reverse, MAPQ 40)
            raw_read(true, true, false, 701, 30),  // 2: supplementary R1
            raw_read(false, true, false, 901, 20), // 3: supplementary R2
        ];
        recs[0].set_template_length(200);
        recs[1].set_template_length(-200);

        fix_supplemental_mate_info(&mut recs, 0, 1);

        // Supp R1 (idx 2) takes primary R2 (idx 1) as its mate.
        assert_eq!(recs[2].mate_ref_id(), 0);
        assert_eq!(recs[2].mate_pos(), 300);
        assert_ne!(recs[2].flags() & rflags::MATE_REVERSE, 0, "primary R2 is reverse");
        assert_eq!(recs[2].tags().find_int(SamTag::MQ), Some(40));
        assert_eq!(recs[2].template_length(), 200); // -(-200)

        // Supp R2 (idx 3) takes primary R1 (idx 0) as its mate.
        assert_eq!(recs[3].mate_ref_id(), 0);
        assert_eq!(recs[3].mate_pos(), 100);
        assert_eq!(recs[3].flags() & rflags::MATE_REVERSE, 0, "primary R1 is forward");
        assert_eq!(recs[3].tags().find_int(SamTag::MQ), Some(60));
        assert_eq!(recs[3].template_length(), -200);
    }

    // The RecordBuf unsoftclipped_start/end helpers (used by the RecordBuf clip path) subtract or
    // add only *soft* clips, ignore hard clips, and return None for unmapped reads.
    #[test]
    fn test_unsoftclipped_recordbuf_helpers() {
        use crate::sam::RecordBuilder;
        use crate::sam::record_utils::{unsoftclipped_end, unsoftclipped_start};

        // 5H10S30M10S at 1-based 100: start = 100 - 10 (leading soft) = 90; hard clips ignored.
        // end = 100 + 30 (ref span) - 1 + 10 (trailing soft) = 139.
        let mapped = RecordBuilder::mapped_read()
            .name("q")
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("5H10S30M10S")
            .sequence(&"A".repeat(50))
            .build();
        assert_eq!(unsoftclipped_start(&mapped), Some(90));
        assert_eq!(unsoftclipped_end(&mapped), Some(139));

        let unmapped = RecordBuilder::mapped_read()
            .name("q")
            .unmapped(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("50M")
            .sequence(&"A".repeat(50))
            .build();
        assert_eq!(unsoftclipped_start(&unmapped), None);
        assert_eq!(unsoftclipped_end(&unmapped), None);
    }

    #[test]
    fn test_default_clip_parameters() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(clip.clipping_mode, ClippingMode::Hard);
        assert!(!clip.clip_overlapping_reads);
        assert!(!clip.clip_extending_past_mate);
    }

    #[test]
    fn test_clip_with_fixed_positions() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 5,
            read_one_three_prime: 3,
            read_two_five_prime: 7,
            read_two_three_prime: 2,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(clip.read_one_five_prime, 5);
        assert_eq!(clip.read_one_three_prime, 3);
        assert_eq!(clip.read_two_five_prime, 7);
        assert_eq!(clip.read_two_three_prime, 2);
    }

    #[test]
    fn test_clip_with_overlapping_enabled() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(clip.clipping_mode, ClippingMode::Hard);
        assert!(clip.clip_overlapping_reads);
        assert!(clip.clip_extending_past_mate);
    }

    #[test]
    fn test_clip_with_metrics_output() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::SoftWithMask,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: Some(PathBuf::from("metrics.txt")),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(clip.clipping_mode, ClippingMode::SoftWithMask);
        assert!(clip.upgrade_clipping);
        assert_eq!(clip.metrics, Some(PathBuf::from("metrics.txt")));
    }

    #[test]
    fn test_clip_with_tag_regeneration() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: true,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(clip.reference, PathBuf::from("reference.fa"));
        assert!(clip.auto_clip_attributes);
    }

    #[test]
    fn test_clip_all_modes_enabled() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,

            read_one_five_prime: 5,
            read_one_three_prime: 5,
            read_two_five_prime: 5,
            read_two_three_prime: 5,
            upgrade_clipping: true,
            auto_clip_attributes: true,
            metrics: Some(PathBuf::from("metrics.txt")),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // All options enabled
        assert!(clip.clip_overlapping_reads);
        assert!(clip.clip_extending_past_mate);
        assert!(clip.upgrade_clipping);
        assert!(clip.auto_clip_attributes);
        assert!(clip.read_one_five_prime > 0);
    }

    #[test]
    fn test_clipping_mode_enum_values() {
        // Test that clipping_mode enum variants are set properly
        let soft = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(soft.clipping_mode, ClippingMode::Soft);
    }

    #[test]
    fn test_clip_asymmetric_fixed_positions() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 10,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 15,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // R1 5' and R2 3' clipping only
        assert_eq!(clip.read_one_five_prime, 10);
        assert_eq!(clip.read_one_three_prime, 0);
        assert_eq!(clip.read_two_five_prime, 0);
        assert_eq!(clip.read_two_three_prime, 15);
    }

    #[test]
    fn test_clip_with_upgrade_all_clipping() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // upgrade_clipping should upgrade existing soft clips to hard clips
        assert!(clip.upgrade_clipping);
        assert_eq!(clip.clipping_mode, ClippingMode::Hard);
    }

    #[test]
    fn test_clip_extending_past_mate_only() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: true,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Only clip_extending_past_mate is enabled
        assert!(!clip.clip_overlapping_reads);
        assert!(clip.clip_extending_past_mate);
    }

    #[test]
    fn test_clip_overlapping_reads_only() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Only clip_overlapping_reads is enabled
        assert!(clip.clip_overlapping_reads);
        assert!(!clip.clip_extending_past_mate);
    }

    #[test]
    fn test_clip_modes_with_auto_clip_attributes() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: true,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // auto_clip_attributes should work with hard clipping
        assert!(clip.auto_clip_attributes);
        assert_eq!(clip.clipping_mode, ClippingMode::Hard);
    }

    #[test]
    fn test_clip_zero_bases_all_positions() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // All fixed position clipping is zero (no fixed clipping)
        assert_eq!(clip.read_one_five_prime, 0);
        assert_eq!(clip.read_one_three_prime, 0);
        assert_eq!(clip.read_two_five_prime, 0);
        assert_eq!(clip.read_two_three_prime, 0);
    }

    #[test]
    fn test_clip_soft_with_mask_mode() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::SoftWithMask,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 5,
            read_one_three_prime: 5,
            read_two_five_prime: 5,
            read_two_three_prime: 5,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(clip.clipping_mode, ClippingMode::SoftWithMask);
        assert!(clip.clip_overlapping_reads);
        assert_eq!(clip.read_one_five_prime, 5);
    }

    #[test]
    fn test_clip_large_fixed_positions() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 50,
            read_one_three_prime: 50,
            read_two_five_prime: 50,
            read_two_three_prime: 50,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Large fixed clipping values (e.g., for adapter trimming)
        assert_eq!(clip.read_one_five_prime, 50);
        assert_eq!(clip.read_one_three_prime, 50);
        assert_eq!(clip.read_two_five_prime, 50);
        assert_eq!(clip.read_two_three_prime, 50);
    }

    #[test]
    fn test_clip_combination_overlapping_and_fixed() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 10,
            read_one_three_prime: 10,
            read_two_five_prime: 10,
            read_two_three_prime: 10,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Both overlapping read clipping and fixed position clipping
        assert!(clip.clip_overlapping_reads);
        assert!(clip.read_one_five_prime > 0);
        assert!(clip.read_one_three_prime > 0);
    }

    #[test]
    fn test_clip_all_three_modes_comparison() {
        let soft = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        let soft_mask = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::SoftWithMask,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        let hard = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Verify all three modes are distinct
        assert_eq!(soft.clipping_mode, ClippingMode::Soft);
        assert_eq!(soft_mask.clipping_mode, ClippingMode::SoftWithMask);
        assert_eq!(hard.clipping_mode, ClippingMode::Hard);
    }

    #[test]
    fn test_clip_single_read_end_clipping() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 20,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Only R2 3' end clipping
        assert_eq!(clip.read_one_five_prime, 0);
        assert_eq!(clip.read_one_three_prime, 0);
        assert_eq!(clip.read_two_five_prime, 0);
        assert_eq!(clip.read_two_three_prime, 20);
    }

    #[test]
    fn test_clip_with_metrics_and_upgrade() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: Some(PathBuf::from("metrics.txt")),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert!(clip.clip_overlapping_reads);
        assert!(clip.upgrade_clipping);
        assert!(clip.metrics.is_some());
    }

    #[test]
    fn test_clip_both_extending_and_overlapping() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert!(clip.clip_overlapping_reads);
        assert!(clip.clip_extending_past_mate);
    }

    // Integration tests
    use crate::sam::{SamBuilder, Strand};
    use anyhow::Result;
    use tempfile::TempDir;

    fn create_test_reference(dir: &TempDir) -> PathBuf {
        let ref_path = dir.path().join("ref.fa");
        // Create a 200bp reference to accommodate read positions + padding
        let ref_content = ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
        std::fs::write(&ref_path, ref_content).expect("failed to write reference FASTA");
        // Also create index (200bp)
        let fai_content = "chr1\t200\t6\t200\t201\n";
        std::fs::write(dir.path().join("ref.fa.fai"), fai_content)
            .expect("failed to write FASTA index");
        ref_path
    }

    fn read_bam_records(path: &std::path::Path) -> Result<Vec<noodles::sam::alignment::RecordBuf>> {
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(path)?;
        let header = reader.read_header()?;
        let records: Vec<_> = reader.record_bufs(&header).collect::<std::io::Result<Vec<_>>>()?;
        Ok(records)
    }

    /// Renders a record's CIGAR as a string (e.g. `"5S15M"`) for assertions.
    fn cigar_string(record: &noodles::sam::alignment::RecordBuf) -> String {
        use noodles::sam::alignment::record::Cigar as _;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use std::fmt::Write as _;
        record.cigar().iter().filter_map(std::result::Result::ok).fold(
            String::new(),
            |mut acc, op| {
                let kind = match op.kind() {
                    Kind::Match => 'M',
                    Kind::Insertion => 'I',
                    Kind::Deletion => 'D',
                    Kind::Skip => 'N',
                    Kind::SoftClip => 'S',
                    Kind::HardClip => 'H',
                    Kind::Pad => 'P',
                    Kind::SequenceMatch => '=',
                    Kind::SequenceMismatch => 'X',
                };
                let _ = write!(acc, "{}{}", op.len(), kind);
                acc
            },
        )
    }

    /// Rewrites `builder`'s `@HD` line to advertise `SO:unsorted GO:query` (query grouped).
    /// `SamBuilder` writes a template's reads adjacently, so the emitted BAM is genuinely
    /// query grouped; the overlap-clip / clip-past-mate pipeline rejects input whose header
    /// does not advertise queryname sorting or query grouping (`clip` requires a template's
    /// reads to be adjacent).
    fn mark_query_grouped(builder: &mut SamBuilder) {
        builder.header = crate::sam::header_as_query_grouped(&builder.header);
    }

    /// CLIP3-02: fixed-position clipping must ensure *at least* N bases are clipped at
    /// the 5'/3' end *including any existing clipping*, matching fgbio `ClipBam`'s use of
    /// `clip5PrimeEndOfRead`/`clip3PrimeEndOfRead`. A read that already carries >= N
    /// bases of clipping at the requested end must not be clipped further.
    #[rstest]
    #[case::single_threaded(ThreadingOptions::none())]
    #[case::multi_threaded(ThreadingOptions::new(2))]
    fn test_fixed_position_clip_counts_existing_clipping(
        #[case] threading: ThreadingOptions,
    ) -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // R1 forward with 5 bases already soft-clipped at its 5' (start) end.
        // R2 reverse with 4 bases already soft-clipped at its 5' (end) end.
        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        let _ = builder
            .add_pair()
            .name("read1")
            .contig(0)
            .start1(10)
            .cigar1("5S15M")
            .bases1("ACGTACGTACGTACGTACGT") // 20 bases
            .start2(120)
            .cigar2("16M4S")
            .bases2("ACGTACGTACGTACGTACGT") // 20 bases
            .build();
        mark_query_grouped(&mut builder);
        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            // Request fewer 5' bases than already exist on each read: fgbio clips 0 more there
            // (the existing-clipping check). Also request 5 bases at R1's 3' end, which has no
            // existing clipping, so clipping *is* applied there — this makes the command path
            // exercise real clipping rather than a no-op.
            read_one_five_prime: 3,
            read_one_three_prime: 5,
            read_two_five_prime: 2,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading,
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        let records = read_bam_records(&output_path)?;
        assert_eq!(records.len(), 2);
        for record in &records {
            let cigar = cigar_string(record);
            if record.flags().is_first_segment() {
                assert_eq!(
                    cigar, "5S10M5S",
                    "R1: existing 5' clip counted (5S retained, no extra 5' clip), and the \
                     requested 5 bases newly clipped at the 3' end"
                );
            } else {
                assert_eq!(
                    cigar, "16M4S",
                    "R2 already had >= 2 bases 5'-clipped; expected no change"
                );
            }
        }
        Ok(())
    }

    /// CLIP3-01 (command level): overlap clipping must be strand-normalized end to end.
    /// A pair whose first-of-pair read is the reverse strand must yield the same
    /// per-strand output as the mirror pair whose first-of-pair read is the forward
    /// strand — the `clip` command must not clip the wrong (outer) ends. Exercised in
    /// both the single-threaded (`threads == 0`) and multi-threaded pipeline paths.
    #[rstest]
    #[case::single_threaded(0)]
    #[case::multi_threaded(4)]
    fn test_clip_overlap_negative_strand_first_matches_mirror(
        #[case] threads: usize,
    ) -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);

        // Runs `clip --clip-overlapping-reads` on an FR pair with the given strands and
        // returns (forward-read CIGAR, reverse-read CIGAR) from the output.
        let run = |tag: &str,
                   r1_minus: bool,
                   r1_start: usize,
                   r2_start: usize|
         -> Result<(String, String)> {
            let input_path = dir.path().join(format!("in_{tag}.bam"));
            let output_path = dir.path().join(format!("out_{tag}.bam"));
            let mut builder = SamBuilder::with_single_ref("chr1", 200);
            let _ = builder
                .add_pair()
                .name("read1")
                .contig(0)
                .start1(r1_start)
                .start2(r2_start)
                .strand1(if r1_minus { Strand::Minus } else { Strand::Plus })
                .strand2(if r1_minus { Strand::Plus } else { Strand::Minus })
                .build();
            mark_query_grouped(&mut builder);
            builder.write(&input_path)?;

            let threading = if threads == 0 {
                ThreadingOptions::none()
            } else {
                ThreadingOptions::new(threads)
            };
            let clip = Clip {
                io: BamIoOptions {
                    input: input_path,
                    output: output_path.clone(),
                    async_reader: false,
                },
                reference: ref_path.clone(),
                clipping_mode: ClippingMode::Soft,
                clip_overlapping_reads: true,
                clip_extending_past_mate: false,
                read_one_five_prime: 0,
                read_one_three_prime: 0,
                read_two_five_prime: 0,
                read_two_three_prime: 0,
                upgrade_clipping: false,
                auto_clip_attributes: false,
                metrics: None,
                threading,
                compression: CompressionOptions { compression_level: 1 },
                scheduler_opts: SchedulerOptions::default(),
                queue_memory: QueueMemoryOptions::default(),
            };
            clip.execute("test")?;

            let records = read_bam_records(&output_path)?;
            assert_eq!(records.len(), 2);
            let forward = records
                .iter()
                .find(|r| !r.flags().is_reverse_complemented())
                .map(cigar_string)
                .expect("a forward-strand read");
            let reverse = records
                .iter()
                .find(|r| r.flags().is_reverse_complemented())
                .map(cigar_string)
                .expect("a reverse-strand read");
            Ok((forward, reverse))
        };

        // Forward-first mirror: R1 forward on the left, R2 reverse on the right.
        let (fwd_forward, fwd_reverse) = run("fwdfirst", false, 10, 20)?;
        // Negative-strand first: R1 reverse on the right, R2 forward on the left.
        let (rev_forward, rev_reverse) = run("revfirst", true, 20, 10)?;

        // Pin the canonical (forward-first) output so the test fails if the command path
        // becomes a no-op: the 100bp reads overlap over [20, 110), so overlap clipping trims
        // 45bp from each read's 3' end where they meet (forward keeps [10, 65) -> 55M45S;
        // reverse keeps [65, 120) -> 45S55M). Both differ from the unclipped 100M input, so
        // a pipeline that silently skipped clipping would not satisfy these assertions.
        assert_eq!(fwd_forward, "55M45S", "forward-read must be 3'-clipped over the overlap");
        assert_eq!(fwd_reverse, "45S55M", "reverse-read must be 3'-clipped over the overlap");

        assert_eq!(
            fwd_forward, rev_forward,
            "forward-read CIGAR must not depend on R1/R2 strand order"
        );
        assert_eq!(
            fwd_reverse, rev_reverse,
            "reverse-read CIGAR must not depend on R1/R2 strand order"
        );
        Ok(())
    }

    /// CLIP3-03: `--upgrade-clipping` in `soft-with-mask` mode must mask existing
    /// soft-clipped bases to N with minimum quality while leaving the CIGAR intact,
    /// matching fgbio `ClipBam --upgrade-clipping --clipping-mode SoftWithMask`.
    #[test]
    fn test_upgrade_clipping_soft_with_mask_masks_existing_soft_bases() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // A fragment with 5 existing soft-clipped bases at the 5' end.
        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        let _ = builder
            .add_frag()
            .name("maskme")
            .contig(0)
            .start(20)
            .cigar("5S15M")
            .bases("ACGTACGTACGTACGTACGT") // 20 bases
            .strand(Strand::Plus)
            .build();
        mark_query_grouped(&mut builder);
        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::SoftWithMask,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        let records = read_bam_records(&output_path)?;
        assert_eq!(records.len(), 1);
        let record = &records[0];
        assert_eq!(cigar_string(record), "5S15M", "CIGAR unchanged");
        let bases = record.sequence().as_ref().to_vec();
        assert_eq!(&bases[..5], b"NNNNN", "leading soft bases masked to N");
        assert!(bases[5..].iter().all(|&b| b != b'N'), "aligned bases untouched");
        let quals = record.quality_scores().as_ref().to_vec();
        assert!(quals[..5].iter().all(|&q| q == 2), "leading soft quals masked to 2");
        Ok(())
    }

    #[test]
    fn test_clip_execute_basic_overlapping() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create overlapping read pair
        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT") // 20 bases
            .contig(0)
            .start1(10) // R1 at pos 10
            .start2(20) // R2 at pos 20 - overlaps with R1
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        let output_records = read_bam_records(&output_path)?;
        assert_eq!(output_records.len(), 2);

        Ok(())
    }

    #[test]
    fn test_clip_execute_soft_mode() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(20)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_soft_with_mask_mode() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(20)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::SoftWithMask,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_fixed_positions() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(30)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 3,
            read_one_three_prime: 2,
            read_two_five_prime: 2,
            read_two_three_prime: 3,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_extending_past_mate() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(20)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: true,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_upgrade_clipping() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(30)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_metrics() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");
        let metrics_path = dir.path().join("metrics.txt");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(20)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: Some(metrics_path.clone()),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());
        assert!(metrics_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_fragment() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create a fragment (unpaired) read
        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_frag()
            .name("frag1")
            .bases("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start(10)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 3,
            read_one_three_prime: 2,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        let output_records = read_bam_records(&output_path)?;
        assert_eq!(output_records.len(), 1);

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_auto_clip_attributes() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(30)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: true,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_all_clipping_options() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");
        let metrics_path = dir.path().join("metrics.txt");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(20)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,

            read_one_five_prime: 2,
            read_one_three_prime: 2,
            read_two_five_prime: 2,
            read_two_three_prime: 2,
            upgrade_clipping: true,
            auto_clip_attributes: true,
            metrics: Some(metrics_path.clone()),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());
        assert!(metrics_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_no_clipping_option() {
        let dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(30)
            .build();

        builder.write(&input_path).expect("failed to write test BAM");

        let clip = Clip {
            io: BamIoOptions { input: input_path, output: output_path, async_reader: false },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false, // No clipping option
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        let result = clip.execute("test");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("At least one clipping option"));
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
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(20)
            .build();
        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading,
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        let output_records = read_bam_records(&output_path)?;
        assert_eq!(output_records.len(), 2, "Should have 2 records");

        Ok(())
    }

    #[test]
    fn test_clip_processed_batch_memory_estimate() {
        let record =
            fgumi_raw_bam::SamBuilder::new().sequence(b"ACGT").qualities(&[30, 30, 30, 30]).build();
        let mut records = Vec::with_capacity(10);
        records.push(record);

        let batch = ClipProcessedBatch {
            clipped_records: records,
            templates_count: 1,
            overlap_clipped_count: 0,
            extend_clipped_count: 0,
        };

        let estimate = batch.estimate_heap_size();
        // Should include Vec overhead for capacity * size_of::<RawRecord>()
        let vec_overhead = 10 * std::mem::size_of::<RawRecord>();
        assert!(
            estimate >= vec_overhead,
            "estimate {estimate} should include Vec<RawRecord> overhead {vec_overhead}"
        );
    }

    /// Helper to create a `Clip` struct with specified clipping parameters and
    /// all other fields set to sensible defaults.
    fn make_clip(
        read_one_five_prime: usize,
        read_one_three_prime: usize,
        read_two_five_prime: usize,
        read_two_three_prime: usize,
    ) -> Clip {
        Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime,
            read_one_three_prime,
            read_two_five_prime,
            read_two_three_prime,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        }
    }

    #[test]
    fn test_clip_fragment_with_metrics() {
        use crate::metrics::clip::ClippingMetricsCollection;

        let clip = make_clip(3, 2, 0, 0);
        let clipper = RawRecordClipper::new(ClippingMode::Soft);
        let mut metrics = ClippingMetricsCollection::new();

        // Build a mapped fragment record with 20M CIGAR
        // SamBuilder pos is 0-based; alignment_start 100 → pos 99
        let mut record = fgumi_raw_bam::SamBuilder::new()
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4]) // 20M
            .ref_id(0)
            .pos(99) // 0-based (alignment_start 100)
            .mapq(60)
            .build();

        clip.clip_fragment(&clipper, &mut record, Some(&mut metrics))
            .expect("clip_fragment should succeed");

        // Fragment metrics should be updated
        assert_eq!(metrics.fragment.reads, 1);
        assert_eq!(metrics.fragment.bases_clipped_five_prime, 3);
        assert_eq!(metrics.fragment.bases_clipped_three_prime, 2);
        // Remaining aligned bases: 20 - 3 - 2 = 15
        assert_eq!(metrics.fragment.bases, 15);
    }

    #[test]
    fn test_clip_fragment_no_clipping_with_metrics() {
        use crate::metrics::clip::ClippingMetricsCollection;

        let clip = make_clip(0, 0, 0, 0);
        // clip_fragment only uses read_one_*; we just test the method directly.
        let clipper = RawRecordClipper::new(ClippingMode::Soft);
        let mut metrics = ClippingMetricsCollection::new();

        let mut record = fgumi_raw_bam::SamBuilder::new()
            .sequence(b"ACGTACGTACGT")
            .qualities(&[30; 12])
            .cigar_ops(&[12u32 << 4]) // 12M
            .ref_id(0)
            .pos(99) // 0-based (alignment_start 100)
            .mapq(60)
            .build();

        clip.clip_fragment(&clipper, &mut record, Some(&mut metrics))
            .expect("clip_fragment should succeed");

        assert_eq!(metrics.fragment.reads, 1);
        assert_eq!(metrics.fragment.bases, 12);
        assert_eq!(metrics.fragment.bases_clipped_five_prime, 0);
        assert_eq!(metrics.fragment.bases_clipped_three_prime, 0);
    }

    #[test]
    fn test_clip_pair_with_metrics() {
        use crate::metrics::clip::ClippingMetricsCollection;
        use fgumi_raw_bam::flags as raw_flags;

        let clip = make_clip(2, 1, 1, 2);
        let clipper = RawRecordClipper::new(ClippingMode::Soft);
        let mut metrics = ClippingMetricsCollection::new();

        // R1: first segment, forward strand, 20M
        // SamBuilder pos is 0-based; alignment_start 100 → pos 99
        let mut r1 = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"pair1")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4]) // 20M
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(99) // 0-based (alignment_start 100)
            .mapq(60)
            .mate_ref_id(0)
            .mate_pos(199) // 0-based (mate_alignment_start 200)
            .template_length(120)
            .build();

        // R2: last segment, reverse strand, 20M (non-overlapping)
        let mut r2 = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"pair1")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4]) // 20M
            .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
            .ref_id(0)
            .pos(199) // 0-based (alignment_start 200)
            .mapq(60)
            .mate_ref_id(0)
            .mate_pos(99) // 0-based (mate_alignment_start 100)
            .template_length(-120)
            .build();

        let (overlap, extend) = clip
            .clip_pair(&clipper, &mut r1, &mut r2, Some(&mut metrics))
            .expect("clip_pair should succeed");

        assert!(!overlap, "no overlapping clipping expected");
        assert!(!extend, "no extending clipping expected");

        // R1 is first segment -> gets read_one clipping: 5'=2, 3'=1
        assert_eq!(metrics.read_one.reads, 1);
        assert_eq!(metrics.read_one.bases_clipped_five_prime, 2);
        assert_eq!(metrics.read_one.bases_clipped_three_prime, 1);

        // R2 is last segment -> gets read_two clipping: 5'=1, 3'=2
        assert_eq!(metrics.read_two.reads, 1);
        assert_eq!(metrics.read_two.bases_clipped_five_prime, 1);
        assert_eq!(metrics.read_two.bases_clipped_three_prime, 2);
    }

    #[test]
    fn test_clip_pair_with_metrics_swapped_flags() {
        use crate::metrics::clip::ClippingMetricsCollection;
        use fgumi_raw_bam::flags as raw_flags;

        // Test the !is_r1_first and !is_r2_last branches:
        // r1 is NOT first_segment, r2 is NOT last_segment
        let clip = make_clip(2, 1, 1, 2);
        let clipper = RawRecordClipper::new(ClippingMode::Soft);
        let mut metrics = ClippingMetricsCollection::new();

        // r1: last segment (not first)
        let mut r1 = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"pair1")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4]) // 20M
            .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT)
            .ref_id(0)
            .pos(99) // 0-based (alignment_start 100)
            .mapq(60)
            .mate_ref_id(0)
            .mate_pos(199) // 0-based (mate_alignment_start 200)
            .template_length(120)
            .build();

        // r2: first segment (not last)
        let mut r2 = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"pair1")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4]) // 20M
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::REVERSE)
            .ref_id(0)
            .pos(199) // 0-based (alignment_start 200)
            .mapq(60)
            .mate_ref_id(0)
            .mate_pos(99) // 0-based (mate_alignment_start 100)
            .template_length(-120)
            .build();

        let (overlap, extend) = clip
            .clip_pair(&clipper, &mut r1, &mut r2, Some(&mut metrics))
            .expect("clip_pair should succeed");

        assert!(!overlap);
        assert!(!extend);

        // r1 is not first_segment -> goes to read_two metrics
        // r1 gets read_two clipping: 5'=1, 3'=2
        assert_eq!(metrics.read_two.reads, 1);
        assert_eq!(metrics.read_two.bases_clipped_five_prime, 1);
        assert_eq!(metrics.read_two.bases_clipped_three_prime, 2);

        // r2 is not last_segment -> goes to read_one metrics
        // r2 gets read_one clipping: 5'=2, 3'=1
        assert_eq!(metrics.read_one.reads, 1);
        assert_eq!(metrics.read_one.bases_clipped_five_prime, 2);
        assert_eq!(metrics.read_one.bases_clipped_three_prime, 1);
    }

    /// `set_mate_info_raw` mirrors htsjdk `SamPairUtil.setMateInfo` across all three
    /// branches (both mapped, one unmapped, both unmapped) so `fgumi clip` matches fgbio
    /// `ClipBam` when clipping unmaps a read (CLIP-01).
    #[test]
    fn test_set_mate_info_raw_all_branches() {
        use fgumi_raw_bam::flags as raw_flags;

        let has_flag = |rec: &RawRecord, f: u16| rec.flags() & f != 0;
        let has_tag = |rec: &RawRecord, tag: [u8; 2]| {
            fgumi_raw_bam::find_tag_type(fgumi_raw_bam::aux_data_slice(rec.as_ref()), tag).is_some()
        };
        // Forward, mapped read: 20M at 0-based pos 99 (alignment start 100).
        let mapped_fwd = || {
            fgumi_raw_bam::SamBuilder::new()
                .read_name(b"t")
                .sequence(b"ACGTACGTACGTACGTACGT")
                .qualities(&[30; 20])
                .cigar_ops(&[20u32 << 4])
                .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .build()
        };
        // Reverse, mapped read: 20M at 0-based pos 199 (alignment start 200).
        let mapped_rev = || {
            fgumi_raw_bam::SamBuilder::new()
                .read_name(b"t")
                .sequence(b"ACGTACGTACGTACGTACGT")
                .qualities(&[30; 20])
                .cigar_ops(&[20u32 << 4])
                .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
                .ref_id(0)
                .pos(199)
                .mapq(40)
                .build()
        };
        // Unmapped read as produced by makeReadUnmapped: no ref/pos, no strand, no CIGAR.
        let unmapped = |last: bool| {
            let seg = if last { raw_flags::LAST_SEGMENT } else { raw_flags::FIRST_SEGMENT };
            fgumi_raw_bam::SamBuilder::new()
                .read_name(b"t")
                .sequence(b"ACGTACGTACGTACGTACGT")
                .qualities(&[30; 20])
                .cigar_ops(&[])
                .flags(raw_flags::PAIRED | seg | raw_flags::UNMAPPED)
                .ref_id(-1)
                .pos(-1)
                .mapq(0)
                .build()
        };

        // --- Branch 1: both mapped ---
        let (mut r1, mut r2) = (mapped_fwd(), mapped_rev());
        set_mate_info_raw(&mut r1, &mut r2);
        assert_eq!((r1.mate_ref_id(), r1.mate_pos()), (0, 199), "r1 mate coords <- r2");
        assert!(has_flag(&r1, raw_flags::MATE_REVERSE), "r1 mate-reverse (r2 is reverse)");
        assert!(!has_flag(&r1, raw_flags::MATE_UNMAPPED), "r1 mate mapped");
        assert_eq!((r2.mate_ref_id(), r2.mate_pos()), (0, 99), "r2 mate coords <- r1");
        assert!(!has_flag(&r2, raw_flags::MATE_REVERSE), "r2 mate not reverse (r1 forward)");
        // TLEN: fwd 5' = 100, rev 5' = alignmentEnd 219 -> 219 - 100 + 1 = 120.
        assert_eq!(r1.template_length(), 120, "r1 TLEN");
        assert_eq!(r2.template_length(), -120, "r2 TLEN");
        assert!(has_tag(&r1, *SamTag::MQ) && has_tag(&r1, *SamTag::MC), "r1 MQ/MC set");
        assert!(has_tag(&r2, *SamTag::MQ) && has_tag(&r2, *SamTag::MC), "r2 MQ/MC set");

        // --- Branch 2: one mapped (r1), one unmapped (r2) ---
        let (mut r1, mut r2) = (mapped_fwd(), unmapped(true));
        set_mate_info_raw(&mut r1, &mut r2);
        // Unmapped r2 is relocated to r1's coordinate.
        assert_eq!((r2.ref_id(), r2.pos()), (0, 99), "unmapped r2 placed at mate coord");
        assert_eq!((r2.mate_ref_id(), r2.mate_pos()), (0, 99), "r2 mate coords <- r1");
        assert!(!has_flag(&r2, raw_flags::MATE_UNMAPPED), "r2's mate (r1) is mapped");
        assert!(
            has_tag(&r2, *SamTag::MQ) && has_tag(&r2, *SamTag::MC),
            "r2 MQ/MC from mapped mate"
        );
        // Mapped r1 points at the co-located unmapped mate.
        assert_eq!((r1.mate_ref_id(), r1.mate_pos()), (0, 99), "r1 mate coords <- unmapped r2");
        assert!(has_flag(&r1, raw_flags::MATE_UNMAPPED), "r1 mate unmapped");
        assert!(!has_flag(&r1, raw_flags::MATE_REVERSE), "r1 mate-reverse cleared (unmapped)");
        assert!(!has_tag(&r1, *SamTag::MQ) && !has_tag(&r1, *SamTag::MC), "r1 MQ/MC removed");
        assert_eq!(
            (r1.template_length(), r2.template_length()),
            (0, 0),
            "TLEN 0 when a mate unmapped"
        );

        // --- Branch 3: both unmapped ---
        let (mut r1, mut r2) = (unmapped(false), unmapped(true));
        set_mate_info_raw(&mut r1, &mut r2);
        for rec in [&r1, &r2] {
            assert_eq!((rec.ref_id(), rec.pos()), (-1, -1), "unmapped coords cleared");
            assert_eq!(
                (rec.mate_ref_id(), rec.mate_pos()),
                (-1, -1),
                "unmapped mate coords cleared"
            );
            assert!(has_flag(rec, raw_flags::MATE_UNMAPPED), "mate-unmapped set");
            assert!(!has_flag(rec, raw_flags::MATE_REVERSE), "mate-reverse cleared");
            assert!(!has_tag(rec, *SamTag::MQ) && !has_tag(rec, *SamTag::MC), "MQ/MC removed");
            assert_eq!(rec.template_length(), 0, "TLEN 0");
        }
    }

    /// htsjdk's one-mapped branch sets the mapped read's mate-reverse flag from the
    /// unmapped read's *actual current* strand (`SamPairUtil.java:267`), not an assumed
    /// `false`. This only diverges for a read that was already unmapped on input while
    /// still carrying the REVERSE flag — the clipper's own unmap clears it, but
    /// `set_mate_info_raw` runs on every pair, so such a read reaches this branch.
    #[test]
    fn test_set_mate_info_raw_one_mapped_uses_unmapped_reads_strand() {
        use fgumi_raw_bam::flags as raw_flags;
        let has_flag = |rec: &RawRecord, f: u16| rec.flags() & f != 0;

        // Mapped forward r1.
        let mut r1 = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"t")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4])
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .build();
        // Unmapped r2 that still carries the REVERSE flag (already unmapped on input).
        let mut r2 = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"t")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[])
            .flags(
                raw_flags::PAIRED
                    | raw_flags::LAST_SEGMENT
                    | raw_flags::UNMAPPED
                    | raw_flags::REVERSE,
            )
            .ref_id(-1)
            .pos(-1)
            .mapq(0)
            .build();

        set_mate_info_raw(&mut r1, &mut r2);

        assert!(
            has_flag(&r1, raw_flags::MATE_REVERSE),
            "mapped r1 mate-reverse must reflect unmapped r2's actual REVERSE flag"
        );
    }
}
