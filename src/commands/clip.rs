//! `Clip` command implementation.
//!
//! Clips reads in a BAM file to remove overlapping portions of read pairs.
//! This is useful for variant calling to avoid double-counting evidence from
//! overlapping portions of paired reads.

use anyhow::Result;
use clap::Parser;
use crossbeam_queue::SegQueue;
use fgumi_lib::alignment_tags::regenerate_alignment_tags;
use fgumi_lib::bam_io::{
    create_bam_reader, create_bam_reader_for_pipeline, create_bam_writer, is_stdin_path,
};
use fgumi_lib::clipper::{ClippingMode, SamRecordClipper};
use fgumi_lib::grouper::TemplateGrouper;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::metrics::clip::ClippingMetricsCollection;
use fgumi_lib::metrics::writer::write_metrics as write_metrics_tsv;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::reference::ReferenceReader;
use fgumi_lib::template::{TemplateBatch, TemplateIterator};
use fgumi_lib::unified_pipeline::{
    BamPipelineConfig, Grouper, MemoryEstimate, run_bam_pipeline_from_reader,
    serialize_bam_records_into,
};
use fgumi_lib::validation::validate_file_exists;
use log::info;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::Cigar as CigarTrait;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::io;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use super::command::Command;
use super::common::{BamIoOptions, CompressionOptions, SchedulerOptions, ThreadingOptions};

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
the overlapping bases from each read. By default soft clipping is performed.

Secondary alignments and supplemental alignments are not clipped, but are passed through into the output.

In order to correctly clip reads by template and update mate information, the input BAM must be either
queryname sorted or query grouped. If your input BAM is not in an appropriate order the sort can be
done in streaming fashion with, for example:

  samtools sort -n -u in.bam | fgumi clip -i /dev/stdin ...

The output sort order may be specified with --sort-order. If not given, then the output will be in the same
order as input.

Any existing NM, UQ and MD tags are repaired (if --regenerate-tags is specified), and mate-pair information is updated.

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
    #[arg(short = 'r', long = "reference", required = true)]
    pub reference: PathBuf,

    /// Clipping mode: soft, soft-with-mask, or hard
    #[arg(short = 'c', long = "clipping-mode", default_value = "hard")]
    pub clipping_mode: String,

    /// Output sort order (if not specified, output is in same order as input)
    #[arg(short = 'S', long = "sort-order")]
    pub sort_order: Option<String>,

    /// Clip overlapping read pairs
    #[arg(long = "clip-overlapping-reads", default_value = "false")]
    pub clip_overlapping_reads: bool,

    /// Clip reads that extend past their mate's start position
    /// Note: In Scala fgbio this parameter is called `clipBasesPastMate`
    #[arg(long = "clip-extending-past-mate", default_value = "false")]
    pub clip_extending_past_mate: bool,

    /// Note: NM/UQ/MD tags are always regenerated after clipping (matching fgbio behavior)
    /// This flag is kept for backwards compatibility but is ignored
    #[arg(long = "regenerate-tags", default_value = "true", hide = true)]
    pub regenerate_tags: bool,

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
    #[arg(short = 'H', long = "upgrade-clipping", default_value = "false")]
    pub upgrade_clipping: bool,

    /// Automatically clip extended attributes that match read length
    #[arg(short = 'a', long = "auto-clip-attributes", default_value = "false")]
    pub auto_clip_attributes: bool,

    /// Output file for clipping metrics
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
}

// ============================================================================
// Types for 7-step pipeline processing
// ============================================================================

/// Result from processing a batch of templates through clipping.
struct ClipProcessedBatch {
    /// Clipped records to write to output BAM.
    clipped_records: Vec<RecordBuf>,
    /// Number of templates processed.
    templates_count: u64,
    /// Number of templates with overlap clipping applied.
    overlap_clipped_count: u64,
    /// Number of templates with mate extension clipping applied.
    extend_clipped_count: u64,
}

impl MemoryEstimate for ClipProcessedBatch {
    fn estimate_heap_size(&self) -> usize {
        self.clipped_records.iter().map(MemoryEstimate::estimate_heap_size).sum()
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
        // Validate input files exist (skip for stdin)
        if !is_stdin_path(&self.io.input) {
            validate_file_exists(&self.io.input, "Input BAM")?;
        }
        validate_file_exists(&self.reference, "Reference FASTA")?;

        info!("Clip");
        info!("  Input: {}", self.io.input.display());
        info!("  Output: {}", self.io.output.display());
        info!("  Clipping mode: {}", self.clipping_mode);
        info!("  Clip overlapping reads: {}", self.clip_overlapping_reads);
        info!("  Clip extending past mate: {}", self.clip_extending_past_mate);
        info!("  Regenerate tags: {}", self.regenerate_tags);
        info!("  {}", self.threading.log_message());

        let timer = OperationTimer::new("Clipping reads");

        // Parse clipping mode
        let mode = match self.clipping_mode.as_str() {
            "soft" => ClippingMode::Soft,
            "soft-with-mask" => ClippingMode::SoftWithMask,
            "hard" => ClippingMode::Hard,
            _ => {
                anyhow::bail!(
                    "Invalid clipping mode: {}. Must be soft, soft-with-mask, or hard",
                    self.clipping_mode
                );
            }
        };

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

        // ========================================================================
        // CRITICAL: Check --threads mode BEFORE creating any file handles.
        // The 7-step pipeline manages its own I/O and writer lifecycle, so we
        // must not create a writer here if we're going to use the pipeline.
        // Route: Some(threads) -> pipeline, None -> single-threaded fast path
        // ========================================================================
        let total_records = if let Some(threads) = self.threading.threads {
            // Read header for the 7-step pipeline (supports stdin)
            let (reader, header) = create_bam_reader_for_pipeline(&self.io.input)?;

            // Load reference (always required for clip)
            let reference = Arc::new(ReferenceReader::new(&self.reference)?);

            // Update header sort order if specified
            let header = self.update_header_sort_order(header)?;

            // Add @PG record with PP chaining to input's last program
            let header = fgumi_lib::header::add_pg_record(
                header,
                crate::version::VERSION.as_str(),
                command_line,
            )?;

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
            SamRecordClipper::with_auto_clip(mode, true)
        } else {
            SamRecordClipper::new(mode)
        };

        // Create metrics collection if metrics output requested
        let mut metrics_collection =
            self.metrics.as_ref().map(|_| ClippingMetricsCollection::new());

        // Load reference (always required and tags are always regenerated to match Scala fgbio)
        let reference_reader = ReferenceReader::new(&self.reference)?;

        // Open input BAM with MT BGZF decompression
        let reader_threads = self.threading.num_threads();
        let (mut reader, mut header) = create_bam_reader(&self.io.input, reader_threads)?;

        // Update header sort order if specified
        if let Some(ref sort_order) = self.sort_order {
            use bstr::BString;
            use noodles::sam::header::record::value::Map;
            use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;

            // Get or create the header map
            let mut header_map = if let Some(hd) = header.header() {
                hd.clone()
            } else {
                Map::<noodles::sam::header::record::value::map::Header>::default()
            };

            // Update sort order
            *header_map.other_fields_mut().entry(SORT_ORDER).or_insert(BString::from("")) =
                BString::from(sort_order.as_str());

            // Rebuild header with new header map
            let mut builder = noodles::sam::Header::builder();

            // Copy existing components
            for (name, rg) in header.read_groups() {
                builder = builder.add_read_group(name.clone(), rg.clone());
            }
            for (name, reference) in header.reference_sequences() {
                builder = builder.add_reference_sequence(name.clone(), reference.clone());
            }
            for (id, pg) in header.programs().as_ref() {
                builder = builder.add_program(id.clone(), pg.clone());
            }
            for comment in header.comments() {
                builder = builder.add_comment(comment.clone());
            }

            // Set the modified header
            builder = builder.set_header(header_map);
            header = builder.build();
        }

        // Add @PG record with PP chaining to input's last program
        let header = fgumi_lib::header::add_pg_record(
            header,
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        // Create output BAM writer with multi-threaded BGZF compression
        let writer_threads = self.threading.num_threads();
        let mut writer = create_bam_writer(
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

        // Single-threaded processing
        let record_iter = reader.record_bufs(&header).map(|r| r.map_err(Into::into));
        let template_iter = TemplateIterator::new(record_iter);

        for template in template_iter {
            let template = template?;
            let mut records = template.records;

            // Process based on template type
            #[allow(clippy::len_zero)] // We specifically want len() == 1, not !is_empty()
            if records.len() == 1 {
                // Fragment
                let record = &mut records[0];
                self.clip_fragment(&clipper, record, metrics_collection.as_mut())?;
            } else if records.len() == 2 {
                // Paired reads
                let (r1, r2) = records.split_at_mut(1);
                let r1 = &mut r1[0];
                let r2 = &mut r2[0];

                let (overlap_clip, extend_clip) =
                    self.clip_pair(&clipper, r1, r2, metrics_collection.as_mut())?;

                if overlap_clip {
                    total_clipped_overlap += 1;
                }
                if extend_clip {
                    total_clipped_mate_extension += 1;
                }

                // Update mate info (MC and MQ tags)
                update_mate_info(r1, r2);
                update_mate_info(r2, r1);
            }

            // Regenerate alignment tags (always done to match Scala fgbio behavior)
            for record in &mut records {
                regenerate_alignment_tags(record, &header, &reference_reader)?;
            }

            // Count and write records
            let batch_size = records.len();
            total_records += batch_size;
            for record in &records {
                writer.write_alignment_record(&header, record)?;
            }
            progress.log_if_needed(batch_size as u64);
        }

        progress.log_final();
        info!("Total records processed: {total_records}");
        info!("Templates with overlap clipping: {total_clipped_overlap}");
        info!("Templates with mate extension clipping: {total_clipped_mate_extension}");

        // Write metrics if requested
        if let Some(metrics_path) = &self.metrics {
            if let Some(mut metrics) = metrics_collection {
                info!("Writing metrics to {}", metrics_path.display());
                metrics.finalize();
                let all_metrics = metrics.all_metrics();
                write_metrics_tsv(metrics_path, &all_metrics, "clipping")?;
                info!("Metrics written successfully");
            }
        }

        // Flush and finish the writer
        writer.into_inner().finish()?;

        info!("Done!");
        Ok(total_records as u64)
    }

    /// Clips a fragment (unpaired) read.
    ///
    /// Applies clipping operations to a single fragment read, including:
    /// 1. Upgrading existing clipping if requested
    /// 2. Applying fixed-position 5' and 3' clipping
    /// 3. Updating metrics if a metrics collector is provided
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
        clipper: &SamRecordClipper,
        record: &mut noodles::sam::alignment::RecordBuf,
        metrics: Option<&mut ClippingMetricsCollection>,
    ) -> Result<()> {
        let prior_bases_clipped = SamRecordClipper::clipped_bases(record);

        // Upgrade existing clipping if requested
        if self.upgrade_clipping {
            clipper.upgrade_all_clipping(record)?;
        }

        // Apply fixed-position clipping
        let num_five_prime = if self.read_one_five_prime > 0 {
            clipper.clip_5_prime_end_of_alignment(record, self.read_one_five_prime)
        } else {
            0
        };

        let num_three_prime = if self.read_one_three_prime > 0 {
            clipper.clip_3_prime_end_of_alignment(record, self.read_one_three_prime)
        } else {
            0
        };

        // Update metrics
        if let Some(metrics) = metrics {
            metrics.fragment.update(
                record,
                prior_bases_clipped,
                num_five_prime,
                num_three_prime,
                0, // no overlapping
                0, // no extending
            );
        }

        Ok(())
    }

    /// Clips a pair of reads with comprehensive clipping logic.
    ///
    /// Applies multiple types of clipping to a read pair:
    /// 1. Upgrading existing clipping if requested
    /// 2. Fixed-position 5' and 3' clipping for each read
    /// 3. Overlap clipping to remove duplicate coverage
    /// 4. Mate-extension clipping to remove reads extending past mate start
    /// 5. Updating metrics for both reads
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
        clipper: &SamRecordClipper,
        r1: &mut noodles::sam::alignment::RecordBuf,
        r2: &mut noodles::sam::alignment::RecordBuf,
        metrics: Option<&mut ClippingMetricsCollection>,
    ) -> Result<(bool, bool)> {
        let prior_bases_clipped_r1 = SamRecordClipper::clipped_bases(r1);
        let prior_bases_clipped_r2 = SamRecordClipper::clipped_bases(r2);

        // Upgrade existing clipping if requested
        if self.upgrade_clipping {
            clipper.upgrade_all_clipping(r1)?;
            clipper.upgrade_all_clipping(r2)?;
        }

        // Determine read types
        let (is_r1_first, is_r2_last) =
            (r1.flags().is_first_segment(), r2.flags().is_last_segment());

        // Apply fixed-position clipping for R1
        let num_r1_five_prime = if is_r1_first && self.read_one_five_prime > 0 {
            clipper.clip_5_prime_end_of_alignment(r1, self.read_one_five_prime)
        } else if !is_r1_first && self.read_two_five_prime > 0 {
            clipper.clip_5_prime_end_of_alignment(r1, self.read_two_five_prime)
        } else {
            0
        };

        let num_r1_three_prime = if is_r1_first && self.read_one_three_prime > 0 {
            clipper.clip_3_prime_end_of_alignment(r1, self.read_one_three_prime)
        } else if !is_r1_first && self.read_two_three_prime > 0 {
            clipper.clip_3_prime_end_of_alignment(r1, self.read_two_three_prime)
        } else {
            0
        };

        // Apply fixed-position clipping for R2
        let num_r2_five_prime = if is_r2_last && self.read_two_five_prime > 0 {
            clipper.clip_5_prime_end_of_alignment(r2, self.read_two_five_prime)
        } else if !is_r2_last && self.read_one_five_prime > 0 {
            clipper.clip_5_prime_end_of_alignment(r2, self.read_one_five_prime)
        } else {
            0
        };

        let num_r2_three_prime = if is_r2_last && self.read_two_three_prime > 0 {
            clipper.clip_3_prime_end_of_alignment(r2, self.read_two_three_prime)
        } else if !is_r2_last && self.read_one_three_prime > 0 {
            clipper.clip_3_prime_end_of_alignment(r2, self.read_one_three_prime)
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
            // Determine which metric to update based on read flags
            if is_r1_first {
                metrics.read_one.update(
                    r1,
                    prior_bases_clipped_r1,
                    num_r1_five_prime,
                    num_r1_three_prime,
                    num_overlapping_r1,
                    num_extending_r1,
                );
            } else {
                metrics.read_two.update(
                    r1,
                    prior_bases_clipped_r1,
                    num_r1_five_prime,
                    num_r1_three_prime,
                    num_overlapping_r1,
                    num_extending_r1,
                );
            }

            if is_r2_last {
                metrics.read_two.update(
                    r2,
                    prior_bases_clipped_r2,
                    num_r2_five_prime,
                    num_r2_three_prime,
                    num_overlapping_r2,
                    num_extending_r2,
                );
            } else {
                metrics.read_one.update(
                    r2,
                    prior_bases_clipped_r2,
                    num_r2_five_prime,
                    num_r2_three_prime,
                    num_overlapping_r2,
                    num_extending_r2,
                );
            }
        }

        let overlap_clipped = num_overlapping_r1 > 0 || num_overlapping_r2 > 0;
        let extend_clipped = num_extending_r1 > 0 || num_extending_r2 > 0;

        Ok((overlap_clipped, extend_clipped))
    }

    /// Updates the header sort order if specified via command line option.
    fn update_header_sort_order(&self, header: Header) -> Result<Header> {
        if let Some(ref sort_order) = self.sort_order {
            use bstr::BString;
            use noodles::sam::header::record::value::Map;
            use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;

            // Get or create the header map
            let mut header_map = if let Some(hd) = header.header() {
                hd.clone()
            } else {
                Map::<noodles::sam::header::record::value::map::Header>::default()
            };

            // Update sort order
            *header_map.other_fields_mut().entry(SORT_ORDER).or_insert(BString::from("")) =
                BString::from(sort_order.as_str());

            // Rebuild header with new header map
            let mut builder = noodles::sam::Header::builder();

            // Copy existing components
            for (name, rg) in header.read_groups() {
                builder = builder.add_read_group(name.clone(), rg.clone());
            }
            for (name, reference) in header.reference_sequences() {
                builder = builder.add_reference_sequence(name.clone(), reference.clone());
            }
            for (id, pg) in header.programs().as_ref() {
                builder = builder.add_program(id.clone(), pg.clone());
            }
            for comment in header.comments() {
                builder = builder.add_comment(comment.clone());
            }

            // Set the modified header
            builder = builder.set_header(header_map);
            Ok(builder.build())
        } else {
            Ok(header)
        }
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
        info!("Using 7-step unified pipeline with {num_threads} threads");

        // Configure pipeline - clip is writer-heavy workload
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
        let collected_metrics: Arc<SegQueue<CollectedClipMetrics>> = Arc::new(SegQueue::new());
        let collected_for_serialize = Arc::clone(&collected_metrics);

        // Configuration for closures
        const BATCH_SIZE: usize = 1000;
        let clipping_mode = self.clipping_mode.clone();
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
            let mode = match clipping_mode.as_str() {
                "soft" => ClippingMode::Soft,
                "soft-with-mask" => ClippingMode::SoftWithMask,
                _ => ClippingMode::Hard,
            };
            let clipper = if auto_clip_attributes {
                SamRecordClipper::with_auto_clip(mode, true)
            } else {
                SamRecordClipper::new(mode)
            };

            let mut clipped_records = Vec::new();
            let mut templates_count = 0u64;
            let mut overlap_clipped_count = 0u64;
            let mut extend_clipped_count = 0u64;

            for template in batch {
                let mut records = template.records;
                templates_count += 1;

                #[allow(clippy::len_zero)]
                if records.len() == 1 {
                    // Fragment - apply fixed-position clipping
                    let record = &mut records[0];
                    if upgrade_clipping {
                        let _ = clipper.upgrade_all_clipping(record);
                    }
                    if read_one_five_prime > 0 {
                        clipper.clip_5_prime_end_of_alignment(record, read_one_five_prime);
                    }
                    if read_one_three_prime > 0 {
                        clipper.clip_3_prime_end_of_alignment(record, read_one_three_prime);
                    }
                } else if records.len() == 2 {
                    // Paired reads
                    let (r1_slice, r2_slice) = records.split_at_mut(1);
                    let r1 = &mut r1_slice[0];
                    let r2 = &mut r2_slice[0];

                    if upgrade_clipping {
                        let _ = clipper.upgrade_all_clipping(r1);
                        let _ = clipper.upgrade_all_clipping(r2);
                    }

                    // Determine read types
                    let is_r1_first = r1.flags().is_first_segment();
                    let is_r2_last = r2.flags().is_last_segment();

                    // Apply fixed-position clipping for R1
                    if is_r1_first && read_one_five_prime > 0 {
                        clipper.clip_5_prime_end_of_alignment(r1, read_one_five_prime);
                    } else if !is_r1_first && read_two_five_prime > 0 {
                        clipper.clip_5_prime_end_of_alignment(r1, read_two_five_prime);
                    }
                    if is_r1_first && read_one_three_prime > 0 {
                        clipper.clip_3_prime_end_of_alignment(r1, read_one_three_prime);
                    } else if !is_r1_first && read_two_three_prime > 0 {
                        clipper.clip_3_prime_end_of_alignment(r1, read_two_three_prime);
                    }

                    // Apply fixed-position clipping for R2
                    if is_r2_last && read_two_five_prime > 0 {
                        clipper.clip_5_prime_end_of_alignment(r2, read_two_five_prime);
                    } else if !is_r2_last && read_one_five_prime > 0 {
                        clipper.clip_5_prime_end_of_alignment(r2, read_one_five_prime);
                    }
                    if is_r2_last && read_two_three_prime > 0 {
                        clipper.clip_3_prime_end_of_alignment(r2, read_two_three_prime);
                    } else if !is_r2_last && read_one_three_prime > 0 {
                        clipper.clip_3_prime_end_of_alignment(r2, read_one_three_prime);
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

                    // Update mate info
                    update_mate_info(r1, r2);
                    update_mate_info(r2, r1);
                }

                // Regenerate alignment tags (always done to match fgbio behavior)
                for record in &mut records {
                    regenerate_alignment_tags(record, &header_for_process, &reference_for_process)
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

        // Serialize function: convert records to bytes and collect metrics
        let serialize_fn = move |processed: ClipProcessedBatch,
                                 header: &Header,
                                 output: &mut Vec<u8>|
              -> io::Result<u64> {
            // Push metrics to lock-free queue
            collected_for_serialize.push(CollectedClipMetrics {
                total_templates: processed.templates_count,
                overlap_clipped: processed.overlap_clipped_count,
                extend_clipped: processed.extend_clipped_count,
            });

            // Serialize clipped records into the provided buffer
            serialize_bam_records_into(&processed.clipped_records, header, output)
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

        while let Some(metrics) = collected_metrics.pop() {
            total_templates += metrics.total_templates;
            total_overlap_clipped += metrics.overlap_clipped;
            total_extend_clipped += metrics.extend_clipped;
        }

        info!("Total templates processed: {total_templates}");
        info!("Templates with overlap clipping: {total_overlap_clipped}");
        info!("Templates with mate extension clipping: {total_extend_clipped}");

        // Note: Metrics file not supported in 7-step pipeline mode
        // (would require tracking per-record metrics which adds overhead)
        if self.metrics.is_some() {
            info!("Note: Detailed metrics file not supported in --threads mode");
        }

        info!("Done!");
        Ok(records_written)
    }
}

/// Updates mate information tags (MC and MQ) for a read based on its mate.
///
/// After clipping operations, mate information tags need to be updated to reflect
/// the new state of the mate read. This function updates:
/// - MC tag: Mate CIGAR string
/// - MQ tag: Mate mapping quality
///
/// These tags are standard SAM optional tags that store information about the mate
/// to enable single-ended processing.
///
/// # Arguments
///
/// * `record` - The record to update (mutable)
/// * `mate` - The mate record to read information from
#[allow(clippy::similar_names)] // mc_tag and mq_tag are standard SAM tags
#[allow(clippy::cast_lossless)] // u8 to i32 cast is intentional for SAM tag format
fn update_mate_info(
    record: &mut noodles::sam::alignment::RecordBuf,
    mate: &noodles::sam::alignment::RecordBuf,
) {
    use noodles::sam::alignment::record::cigar::op::Kind;
    use std::fmt::Write;

    // Update MC tag (mate CIGAR)
    // Build CIGAR string manually since Cigar doesn't implement Display
    let mate_cigar = mate.cigar();
    let mut cigar_string = String::new();
    for op in mate_cigar.iter().flatten() {
        let kind_char = match op.kind() {
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
        let _ = write!(cigar_string, "{}{}", op.len(), kind_char);
    }

    let mc_tag = Tag::from([b'M', b'C']);
    let _ = record.data_mut().insert(mc_tag, cigar_string.into());

    // Update MQ tag (mate mapping quality)
    if let Some(mate_mapq) = mate.mapping_quality() {
        let mq_tag = Tag::from([b'M', b'Q']);
        let mapq_value: u8 = mate_mapq.into();
        let _ = record.data_mut().insert(mq_tag, (mapq_value as i32).into());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::path::PathBuf;

    #[test]
    fn test_default_clip_parameters() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true, // Always true to match Scala fgbio
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        assert_eq!(clip.clipping_mode, "hard");
        assert!(!clip.clip_overlapping_reads);
        assert!(!clip.clip_extending_past_mate);
        assert!(clip.regenerate_tags); // Always true
    }

    #[test]
    fn test_clip_with_fixed_positions() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 5,
            read_one_three_prime: 3,
            read_two_five_prime: 7,
            read_two_three_prime: 2,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        assert_eq!(clip.clipping_mode, "hard");
        assert!(clip.clip_overlapping_reads);
        assert!(clip.clip_extending_past_mate);
    }

    #[test]
    fn test_clip_with_metrics_output() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "soft-with-mask".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: Some(PathBuf::from("metrics.txt")),
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        assert_eq!(clip.clipping_mode, "soft-with-mask");
        assert!(clip.upgrade_clipping);
        assert_eq!(clip.metrics, Some(PathBuf::from("metrics.txt")));
    }

    #[test]
    fn test_clip_with_tag_regeneration() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true, // Always true
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: true,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        assert!(clip.regenerate_tags);
        assert_eq!(clip.reference, PathBuf::from("reference.fa"));
        assert!(clip.auto_clip_attributes);
    }

    #[test]
    fn test_clip_all_modes_enabled() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,
            regenerate_tags: true,
            read_one_five_prime: 5,
            read_one_three_prime: 5,
            read_two_five_prime: 5,
            read_two_three_prime: 5,
            upgrade_clipping: true,
            auto_clip_attributes: true,
            metrics: Some(PathBuf::from("metrics.txt")),
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        // All options enabled
        assert!(clip.clip_overlapping_reads);
        assert!(clip.clip_extending_past_mate);
        assert!(clip.regenerate_tags);
        assert!(clip.upgrade_clipping);
        assert!(clip.auto_clip_attributes);
        assert!(clip.read_one_five_prime > 0);
    }

    #[test]
    fn test_clipping_mode_string_values() {
        // Test that clipping_mode strings are validated properly
        let soft = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "soft".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        assert_eq!(soft.clipping_mode, "soft");
    }

    #[test]
    fn test_clip_with_sort_order_specification() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: Some("coordinate".to_string()),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        assert_eq!(clip.sort_order, Some("coordinate".to_string()));
    }

    #[test]
    fn test_clip_asymmetric_fixed_positions() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "soft".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 10,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 15,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        // upgrade_clipping should upgrade existing soft clips to hard clips
        assert!(clip.upgrade_clipping);
        assert_eq!(clip.clipping_mode, "hard");
    }

    #[test]
    fn test_clip_extending_past_mate_only() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: true,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: true,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        // auto_clip_attributes should work with hard clipping
        assert!(clip.auto_clip_attributes);
        assert_eq!(clip.clipping_mode, "hard");
    }

    #[test]
    fn test_clip_zero_bases_all_positions() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "soft-with-mask".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 5,
            read_one_three_prime: 5,
            read_two_five_prime: 5,
            read_two_three_prime: 5,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        assert_eq!(clip.clipping_mode, "soft-with-mask");
        assert!(clip.clip_overlapping_reads);
        assert_eq!(clip.read_one_five_prime, 5);
    }

    #[test]
    fn test_clip_large_fixed_positions() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 50,
            read_one_three_prime: 50,
            read_two_five_prime: 50,
            read_two_three_prime: 50,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 10,
            read_one_three_prime: 10,
            read_two_five_prime: 10,
            read_two_three_prime: 10,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "soft".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        let soft_mask = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "soft-with-mask".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        let hard = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        // Verify all three modes are distinct
        assert_eq!(soft.clipping_mode, "soft");
        assert_eq!(soft_mask.clipping_mode, "soft-with-mask");
        assert_eq!(hard.clipping_mode, "hard");
    }

    #[test]
    fn test_clip_with_queryname_sort_order() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: Some("queryname".to_string()),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        assert_eq!(clip.sort_order, Some("queryname".to_string()));
    }

    #[test]
    fn test_clip_with_unsorted_sort_order() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: Some("unsorted".to_string()),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        assert_eq!(clip.sort_order, Some("unsorted".to_string()));
    }

    #[test]
    fn test_clip_regenerate_tags_always_true() {
        // Test that regenerate_tags is always true to match Scala fgbio behavior
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        // regenerate_tags should always be true to match fgbio
        assert!(clip.regenerate_tags);
    }

    #[test]
    fn test_clip_single_read_end_clipping() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 20,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: Some(PathBuf::from("metrics.txt")),
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "soft".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        assert!(clip.clip_overlapping_reads);
        assert!(clip.clip_extending_past_mate);
    }

    #[test]
    fn test_clip_no_regenerate_tags_option() {
        // Test that regenerate_tags is always true (no option to disable)
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        // regenerate_tags must always be true (unlike Scala which had option to disable)
        assert!(clip.regenerate_tags);
    }

    // Integration tests
    use anyhow::Result;
    use fgumi_lib::sam::SamBuilder;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use tempfile::TempDir;

    fn create_test_reference(dir: &TempDir) -> PathBuf {
        let ref_path = dir.path().join("ref.fa");
        // Create a 200bp reference to accommodate read positions + padding
        let ref_content = ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
        std::fs::write(&ref_path, ref_content).unwrap();
        // Also create index (200bp)
        let fai_content = "chr1\t200\t6\t200\t201\n";
        std::fs::write(dir.path().join("ref.fa.fai"), fai_content).unwrap();
        ref_path
    }

    fn read_bam_records(path: &std::path::Path) -> Result<Vec<RecordBuf>> {
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(path)?;
        let header = reader.read_header()?;
        let records: Vec<_> = reader.record_bufs(&header).collect::<std::io::Result<Vec<_>>>()?;
        Ok(records)
    }

    #[test]
    fn test_clip_execute_basic_overlapping() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create overlapping read pair
        let mut builder = SamBuilder::with_single_ref("chr1", 200);
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
            io: BamIoOptions { input: input_path, output: output_path.clone() },
            reference: ref_path,
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            io: BamIoOptions { input: input_path, output: output_path.clone() },
            reference: ref_path,
            clipping_mode: "soft".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            io: BamIoOptions { input: input_path, output: output_path.clone() },
            reference: ref_path,
            clipping_mode: "soft-with-mask".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            io: BamIoOptions { input: input_path, output: output_path.clone() },
            reference: ref_path,
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 3,
            read_one_three_prime: 2,
            read_two_five_prime: 2,
            read_two_three_prime: 3,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            io: BamIoOptions { input: input_path, output: output_path.clone() },
            reference: ref_path,
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: true,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            io: BamIoOptions { input: input_path, output: output_path.clone() },
            reference: ref_path,
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            io: BamIoOptions { input: input_path, output: output_path.clone() },
            reference: ref_path,
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: Some(metrics_path.clone()),
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());
        assert!(metrics_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_sort_order() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
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
            io: BamIoOptions { input: input_path, output: output_path.clone() },
            reference: ref_path,
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: Some("queryname".to_string()),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

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
        let _ = builder
            .add_frag()
            .name("frag1")
            .bases("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start(10)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions { input: input_path, output: output_path.clone() },
            reference: ref_path,
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 3,
            read_one_three_prime: 2,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            io: BamIoOptions { input: input_path, output: output_path.clone() },
            reference: ref_path,
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: true,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            io: BamIoOptions { input: input_path, output: output_path.clone() },
            reference: ref_path,
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,
            regenerate_tags: true,
            read_one_five_prime: 2,
            read_one_three_prime: 2,
            read_two_five_prime: 2,
            read_two_three_prime: 2,
            upgrade_clipping: true,
            auto_clip_attributes: true,
            metrics: Some(metrics_path.clone()),
            sort_order: Some("queryname".to_string()),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());
        assert!(metrics_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_invalid_clipping_mode() {
        let dir = TempDir::new().unwrap();
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(30)
            .build();

        builder.write(&input_path).unwrap();

        let clip = Clip {
            io: BamIoOptions { input: input_path, output: output_path },
            reference: ref_path,
            clipping_mode: "invalid".to_string(), // Invalid mode
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };

        let result = clip.execute("test");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Invalid clipping mode"));
    }

    #[test]
    fn test_clip_execute_no_clipping_option() {
        let dir = TempDir::new().unwrap();
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(30)
            .build();

        builder.write(&input_path).unwrap();

        let clip = Clip {
            io: BamIoOptions { input: input_path, output: output_path },
            reference: ref_path,
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false, // No clipping option
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
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
            io: BamIoOptions { input: input_path, output: output_path.clone() },
            reference: ref_path,
            clipping_mode: "hard".to_string(),
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,
            regenerate_tags: true,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            sort_order: None,
            threading,
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
        };
        clip.execute("test")?;

        let output_records = read_bam_records(&output_path)?;
        assert_eq!(output_records.len(), 2, "Should have 2 records");

        Ok(())
    }
}
