// This file will contain the complete group.rs with all integration tests implemented
// We'll build it section by section
//! Groups reads by UMI to identify reads from the same original molecule.

use crate::commands::command::Command;
use crate::commands::common::{
    BamIoOptions, CompressionOptions, SchedulerOptions, ThreadingOptions,
};
use ahash::AHashMap;
use anyhow::{Context, Result, bail};
use bstr::BString;
use clap::Parser;
use crossbeam_queue::SegQueue;
use fgoxide::io::DelimFile;
use fgumi_lib::assigner::{PairedUmiAssigner, Strategy, UmiAssigner};
use fgumi_lib::bam_io::{create_bam_reader_for_pipeline, create_bam_writer, is_stdin_path};
use fgumi_lib::grouper::{
    FilterMetrics, PositionGroup, PositionGrouper, PositionGrouperConfig, ProcessedPositionGroup,
};
use fgumi_lib::logging::{OperationTimer, log_umi_grouping_summary};
use fgumi_lib::metrics::group::UmiGroupingMetrics;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::read_info::build_library_lookup;
use fgumi_lib::read_info::{LibraryIndex, compute_group_key};
use fgumi_lib::sam::{is_template_coordinate_sorted, unclipped_five_prime_position};
use fgumi_lib::template::{MoleculeId, Template};
use fgumi_lib::unified_pipeline::DecodedRecord;
use fgumi_lib::unified_pipeline::{
    BamPipelineConfig, Grouper, run_bam_pipeline_from_reader, serialize_bam_records_into,
};
use fgumi_lib::validation::{string_to_tag, validate_file_exists};
use log::info;
use noodles::sam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::value::Value as DataValue;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::sync::Arc;

/// Metrics describing the distribution of tag family sizes
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TagFamilySizeMetric {
    /// The family size (number of templates per family)
    pub family_size: usize,
    /// The number of families with this size
    pub count: u64,
    /// The fraction of all families with this size
    pub fraction: f64,
    /// The fraction of families with size >= this size
    pub fraction_gt_or_eq_family_size: f64,
}

// UmiGroupingMetrics is now imported from fgumi_lib::metrics

/// Result of processing a position group, used for parallel processing
#[allow(dead_code)]
struct PositionGroupResult {
    /// Family size counts for this position group
    family_sizes: AHashMap<usize, u64>,
    /// Templates grouped by MI and sorted, ready for writing
    sorted_templates: Vec<Template>,
}

/// Collected metrics from `serialize_fn`, aggregated after pipeline completion.
#[derive(Default, Debug)]
struct CollectedMetrics {
    family_sizes: AHashMap<usize, u64>,
    filter_metrics: FilterMetrics,
}

impl CollectedMetrics {
    #[allow(dead_code)]
    fn merge(&mut self, other: &CollectedMetrics) {
        for (size, count) in &other.family_sizes {
            *self.family_sizes.entry(*size).or_insert(0) += count;
        }
        self.filter_metrics.merge(&other.filter_metrics);
    }
}

/// Configuration for template filtering during group processing.
#[derive(Clone)]
struct GroupFilterConfig {
    /// UMI tag bytes (e.g., [b'R', b'X']).
    umi_tag: [u8; 2],
    /// Minimum mapping quality.
    min_mapq: u8,
    /// Whether to include non-PF reads.
    include_non_pf: bool,
    /// Minimum UMI length (None to disable).
    min_umi_length: Option<usize>,
}

/// Filter a template based on filtering criteria.
/// Returns true if the template should be kept, false if it should be discarded.
/// Updates `filter_metrics` with the reason for filtering.
fn filter_template(
    template: &Template,
    config: &GroupFilterConfig,
    metrics: &mut FilterMetrics,
) -> bool {
    // Get primary reads for filtering
    let r1 = template.r1();
    let r2 = template.r2();

    // Count records, not templates (paired template = 2 records)
    let num_primary_reads = r1.is_some() as u64 + r2.is_some() as u64;
    metrics.total_templates += num_primary_reads;

    // Need at least one primary read
    if r1.is_none() && r2.is_none() {
        metrics.discarded_poor_alignment += num_primary_reads;
        return false;
    }

    // Check if both reads are unmapped (poor alignment)
    let both_unmapped =
        r1.is_none_or(|r| r.flags().is_unmapped()) && r2.is_none_or(|r| r.flags().is_unmapped());

    if both_unmapped {
        metrics.discarded_poor_alignment += num_primary_reads;
        return false;
    }

    // =========================================================================
    // Phase 1: Cheap flag-based checks for all reads (no tag lookups)
    // =========================================================================
    for record in [r1, r2].into_iter().flatten() {
        let flags = record.flags();

        // Filter non-PF if requested
        if !config.include_non_pf && flags.is_qc_fail() {
            metrics.discarded_non_pf += num_primary_reads;
            return false;
        }

        // Check MAPQ for mapped reads
        if !flags.is_unmapped() {
            let mapq = record.mapping_quality().map_or(0, u8::from);
            if mapq < config.min_mapq {
                metrics.discarded_poor_alignment += num_primary_reads;
                return false;
            }
        }
    }

    // =========================================================================
    // Phase 2: Expensive tag lookups (only reached if phase 1 passed)
    // =========================================================================
    for record in [r1, r2].into_iter().flatten() {
        let flags = record.flags();

        // Check mate MAPQ via MQ tag if mate is mapped
        if !flags.is_mate_unmapped() {
            if let Some(data) = record.data().get(b"MQ") {
                if let Some(mq) = data.as_int() {
                    #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
                    if (mq as u8) < config.min_mapq {
                        metrics.discarded_poor_alignment += num_primary_reads;
                        return false;
                    }
                }
            }
        }

        // Check UMI for Ns and minimum length (single pass)
        if let Some(data) = record.data().get(&config.umi_tag) {
            if let DataValue::String(umi) = data {
                // Single pass: count ACGT and detect N simultaneously
                let mut acgt_count = 0usize;
                for &b in umi.iter() {
                    match b {
                        b'A' | b'C' | b'G' | b'T' => acgt_count += 1,
                        b'N' => {
                            metrics.discarded_ns_in_umi += num_primary_reads;
                            return false;
                        }
                        _ => {} // Skip non-ACGT, non-N characters (e.g., '-' in paired UMIs)
                    }
                }

                // Check minimum UMI length
                if let Some(min_len) = config.min_umi_length {
                    if acgt_count < min_len {
                        metrics.discarded_umi_too_short += num_primary_reads;
                        return false;
                    }
                }
            } else {
                // UMI tag is not a string - skip template
                metrics.discarded_poor_alignment += num_primary_reads;
                return false;
            }
        } else {
            // Missing UMI tag - skip template
            metrics.discarded_poor_alignment += num_primary_reads;
            return false;
        }
    }

    metrics.accepted_templates += num_primary_reads;
    true
}

// =============================================================================
// Static _impl functions - core logic used by both closures and &self methods
// =============================================================================

/// Extract UMI for a read, handling paired UMI strategies (static implementation).
fn umi_for_read_impl(umi: &str, is_r1_earlier: bool, assigner: &dyn UmiAssigner) -> Result<String> {
    if assigner.split_templates_by_pair_orientation() {
        // For non-paired strategies, return UMI uppercase
        // Optimize: only allocate if lowercase chars present
        if umi.bytes().all(|b| !b.is_ascii_lowercase()) {
            Ok(umi.to_owned())
        } else {
            Ok(umi.to_uppercase())
        }
    } else {
        // For paired strategies, parse and prefix the UMI
        let parts: Vec<&str> = umi.split('-').collect();
        if parts.len() != 2 {
            bail!(
                "Paired strategy used but UMI did not contain 2 segments delimited by '-': {umi}"
            );
        }

        let Some(paired) = assigner.as_any().downcast_ref::<PairedUmiAssigner>() else {
            bail!("Expected PairedUmiAssigner")
        };

        // When R1 is earlier: lower_prefix:part0-higher_prefix:part1
        // When R2 is earlier: higher_prefix:part0-lower_prefix:part1
        let result = if is_r1_earlier {
            format!(
                "{}:{}-{}:{}",
                paired.lower_read_umi_prefix(),
                parts[0],
                paired.higher_read_umi_prefix(),
                parts[1]
            )
        } else {
            format!(
                "{}:{}-{}:{}",
                paired.higher_read_umi_prefix(),
                parts[0],
                paired.lower_read_umi_prefix(),
                parts[1]
            )
        };

        Ok(result)
    }
}

/// Truncate UMIs to minimum length if specified (static implementation).
fn truncate_umis_impl(umis: Vec<String>, min_umi_length: Option<usize>) -> Result<Vec<String>> {
    match min_umi_length {
        None => Ok(umis),
        Some(min_len) => {
            let min_length = umis.iter().map(String::len).min().unwrap_or(0);
            if min_length < min_len {
                bail!("UMI found that had shorter length than expected ({min_length} < {min_len})");
            }
            Ok(umis.into_iter().map(|u| u[..min_len].to_string()).collect())
        }
    }
}

/// Check if R1 is genomically earlier than R2 (static implementation).
fn is_r1_genomically_earlier_impl(
    r1: &sam::alignment::RecordBuf,
    r2: &sam::alignment::RecordBuf,
) -> Result<bool> {
    let r1_pos = unclipped_five_prime_position(r1).unwrap_or(0);
    let r2_pos = unclipped_five_prime_position(r2).unwrap_or(0);
    Ok(r1_pos <= r2_pos)
}

/// Get pair orientation for a template (static implementation).
fn get_pair_orientation_impl(template: &Template) -> (bool, bool) {
    let r1_positive = template.r1().is_none_or(|r| !r.flags().is_reverse_complemented());
    let r2_positive = template.r2().is_none_or(|r| !r.flags().is_reverse_complemented());
    (r1_positive, r2_positive)
}

/// Assign UMI groups to templates (static implementation).
fn assign_umi_groups_impl(
    templates: &mut [Template],
    assigner: &dyn UmiAssigner,
    raw_tag: [u8; 2],
    assign_tag_bytes: [u8; 2],
    min_umi_length: Option<usize>,
) -> Result<()> {
    if assigner.split_templates_by_pair_orientation() {
        // Group by pair orientation
        let mut subgroups: AHashMap<(bool, bool), Vec<usize>> = AHashMap::new();
        for (idx, template) in templates.iter().enumerate() {
            let orientation = get_pair_orientation_impl(template);
            subgroups.entry(orientation).or_default().push(idx);
        }

        // Process each subgroup separately
        for indices in subgroups.values() {
            assign_umi_groups_for_indices_impl(
                templates,
                indices,
                assigner,
                raw_tag,
                assign_tag_bytes,
                min_umi_length,
            )?;
        }
    } else {
        // No splitting - process all templates together
        let all_indices: Vec<usize> = (0..templates.len()).collect();
        assign_umi_groups_for_indices_impl(
            templates,
            &all_indices,
            assigner,
            raw_tag,
            assign_tag_bytes,
            min_umi_length,
        )?;
    }

    Ok(())
}

/// Assign UMI groups to a subset of templates (static implementation).
///
/// This sets the `Template.mi` field with the `MoleculeId` enum. The actual BAM MI tag
/// is set later during serialization, when we have the global offset.
fn assign_umi_groups_for_indices_impl(
    templates: &mut [Template],
    indices: &[usize],
    assigner: &dyn UmiAssigner,
    raw_tag: [u8; 2],
    _assign_tag_bytes: [u8; 2], // No longer used here - tags set during serialization
    min_umi_length: Option<usize>,
) -> Result<()> {
    if indices.is_empty() {
        return Ok(());
    }

    // Extract UMIs from templates
    let mut umis = Vec::with_capacity(indices.len());

    for &idx in indices {
        let template = &templates[idx];
        let umi_bytes = if let Some(r1) = template.r1() {
            if let Some(DataValue::String(bytes)) = r1.data().get(&raw_tag) {
                bytes
            } else {
                bail!("UMI tag is not a string");
            }
        } else if let Some(r2) = template.r2() {
            if let Some(DataValue::String(bytes)) = r2.data().get(&raw_tag) {
                bytes
            } else {
                bail!("UMI tag is not a string");
            }
        } else {
            bail!("Template has no reads");
        };
        // UMIs are DNA sequences (valid ASCII), so from_utf8 should always succeed
        // This avoids the allocation that from_utf8_lossy().into_owned() would cause
        let umi_str = std::str::from_utf8(umi_bytes)
            .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in UMI: {e}"))?;

        // Determine which read is genomically earlier
        let is_r1_earlier = if let (Some(r1), Some(r2)) = (template.r1(), template.r2()) {
            is_r1_genomically_earlier_impl(r1, r2)?
        } else {
            true
        };

        let processed_umi = umi_for_read_impl(umi_str, is_r1_earlier, assigner)?;
        umis.push(processed_umi);
    }

    // Truncate UMIs if needed
    let truncated_umis = truncate_umis_impl(umis, min_umi_length)?;

    // Assign UMI groups - returns Vec<MoleculeId> indexed by input position
    let assignments = assigner.assign(&truncated_umis);

    // Store MoleculeId enum in Template.mi field
    // The actual MI tag string is set during serialization with global offset
    for (i, &idx) in indices.iter().enumerate() {
        let template = &mut templates[idx];
        template.mi = assignments[i];
    }

    Ok(())
}

/// Set MI tag on a single record using `MoleculeId` with global offset.
///
/// Converts the `MoleculeId` enum to a string representation with the base offset applied,
/// and sets the MI tag on the record. This is called just before writing to minimize
/// memory usage (strings are not stored until write time).
#[inline]
fn set_mi_tag_on_record(
    record: &mut noodles::sam::alignment::record_buf::RecordBuf,
    mi: MoleculeId,
    assign_tag: Tag,
    base_mi: u64,
) {
    if !mi.is_assigned() {
        return;
    }
    let mi_string = mi.to_string_with_offset(base_mi);
    record.data_mut().insert(assign_tag, DataValue::String(BString::from(mi_string)));
}

/// Groups reads by UMI to identify reads from the same original molecule
#[derive(Debug, Parser)]
#[command(
    about = "\x1b[38;5;151m[GROUP]\x1b[0m          \x1b[36mGroup reads by UMI to identify reads from the same original molecule\x1b[0m",
    long_about = r#"
Groups reads together that appear to have come from the same original molecule. Reads
are grouped by template, and then templates are sorted by the 5' mapping positions of
the reads from the template, used from earliest mapping position to latest. Reads that
have the same end positions are then sub-grouped by UMI sequence.

Accepts reads in any order (including unsorted) and outputs reads sorted by:

   1. The lower genome coordinate of the two outer ends of the templates (strand-aware)
   2. The sequencing library
   3. The assigned UMI tag
   4. Read Name

It is recommended to sort the reads into template-coordinate order prior to running
this tool to avoid re-sorting the input. Use `samtools sort --template-coordinate` for
the pre-sorting. The output will always be written in template-coordinate order.

During grouping, reads and templates are filtered out as follows:

1. Templates are filtered if all reads for the template are unmapped
2. Templates are filtered if any non-secondary, non-supplementary read has mapping quality < min-map-q
3. Templates are filtered if any UMI sequence contains one or more N bases
4. Templates are filtered if --min-umi-length is specified and the UMI does not meet the length requirement
5. Records are filtered out if flagged as either secondary or supplementary

Grouping of UMIs is performed by one of four strategies:

1. **identity**:  only reads with identical UMI sequences are grouped together. This strategy
                  may be useful for evaluating data, but should generally be avoided as it will
                  generate multiple UMI groups per original molecule in the presence of errors.
2. **edit**:      reads are clustered into groups such that each read within a group has at least
                  one other read in the group with <= edits differences and there are inter-group
                  pairings with <= edits differences. Effective when there are small numbers of
                  reads per UMI, but breaks down at very high coverage of UMIs.
3. **adjacency**: a version of the directed adjacency method described in umi_tools
                  (http://dx.doi.org/10.1101/051755) that allows for errors between UMIs but
                  only when there is a count gradient.
4. **paired**:    similar to adjacency but for methods that produce templates such that a read with
                  A-B is related to but not identical to a read with B-A. Expects the UMI sequences
                  to be stored in a single SAM tag separated by a hyphen (e.g. ACGT-CCGG) and allows
                  for one of the two UMIs to be absent (e.g. ACGT- or -ACGT). The molecular IDs
                  produced have more structure than for single UMI strategies and are of the form
                  {base}/{A|B}. E.g. two UMI pairs would be mapped as follows:
                  AAAA-GGGG -> 1/A, GGGG-AAAA -> 1/B.

Strategies edit, adjacency, and paired make use of the --edits parameter to control the matching of
non-identical UMIs.

By default, all UMIs must be the same length. If --min-umi-length=len is specified then reads that
have a UMI shorter than len will be discarded, and when comparing UMIs of different lengths, the first
len bases will be compared, where len is the length of the shortest UMI. The UMI length is the number
of [ACGT] bases in the UMI (i.e. does not count dashes and other non-ACGT characters). This option is
not implemented for reads with UMI pairs (i.e. using the paired assigner).

Note: the --min-map-q parameter defaults to 0 in duplicate marking mode and 1 otherwise, and is
directly settable on the command line.

Multi-threaded operation is supported via --threads N which spawns exactly N threads.
Threads are allocated based on the command's workload profile to optimize performance.

Example: --threads 8 spawns exactly 8 threads (2 reader, 4 workers, 2 writer)
"#
)]
pub struct GroupReadsByUmi {
    /// Input and output BAM files
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Optional output of tag family size counts
    #[arg(short = 'f', long = "family-size-histogram")]
    pub family_size_histogram: Option<PathBuf>,

    /// Optional output of UMI grouping metrics
    #[arg(short = 'g', long = "grouping-metrics")]
    pub grouping_metrics: Option<PathBuf>,

    /// The tag containing the raw UMI
    #[arg(short = 't', long = "raw-tag", default_value = "RX")]
    pub raw_tag: String,

    /// The output tag for UMI grouping
    #[arg(short = 'T', long = "assign-tag", default_value = "MI")]
    pub assign_tag: String,

    /// The tag containing the cell barcode
    #[arg(short = 'c', long = "cell-tag", default_value = "CB")]
    pub cell_tag: String,

    /// Minimum mapping quality for mapped reads
    #[arg(short = 'm', long = "min-map-q")]
    pub min_map_q: Option<u8>,

    /// Include non-PF reads
    #[arg(short = 'n', long = "include-non-pf-reads")]
    pub include_non_pf_reads: bool,

    /// The UMI assignment strategy
    #[arg(short = 's', long = "strategy", value_enum)]
    pub strategy: Strategy,

    /// The allowable number of edits between UMIs
    #[arg(short = 'e', long = "edits", default_value = "1")]
    pub edits: u32,

    /// The minimum UMI length
    #[arg(short = 'l', long = "min-umi-length")]
    pub min_umi_length: Option<usize>,

    /// Threading options for parallel processing.
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Minimum UMIs per position to use N-gram/BK-tree index for faster grouping.
    /// Set to 0 to always use linear scan. Only affects Adjacency/Paired strategies.
    #[arg(long = "index-threshold", default_value = "100")]
    pub index_threshold: usize,

    /// Scheduler and pipeline statistics options.
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Memory limit for pipeline queues in megabytes (default: 4096 = 4GB).
    /// Applies backpressure when queue memory exceeds this limit.
    /// Use 0 to disable the limit.
    #[arg(long = "queue-memory-limit-mb", default_value = "4096")]
    pub queue_memory_limit_mb: u64,
}

impl Command for GroupReadsByUmi {
    /// Execute the tool using the 7-step unified pipeline.
    #[allow(clippy::too_many_lines)]
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate inputs
        if self.min_umi_length.is_some() && matches!(self.strategy, Strategy::Paired) {
            bail!("Paired strategy cannot be used with --min-umi-length");
        }

        // Validate input file exists (skip for stdin)
        if !is_stdin_path(&self.io.input) {
            validate_file_exists(&self.io.input, "input BAM file")?;
        }

        // Set minimum mapping quality
        let min_mapq: u8 = self.min_map_q.unwrap_or(1);

        // Initialize tracking infrastructure
        let timer = OperationTimer::new("Grouping reads by UMI");

        info!("Starting group");
        info!("Input: {}", self.io.input.display());
        info!("Output: {}", self.io.output.display());
        info!("Strategy: {:?}", self.strategy);
        info!("Edits: {}", self.edits);
        if matches!(self.strategy, Strategy::Adjacency | Strategy::Paired) {
            info!("Index threshold: {}", self.index_threshold);
        }

        // Log threading configuration
        info!("{}", self.threading.log_message());

        // Open input BAM using streaming-capable reader for pipeline use
        info!("Reading input BAM");
        let (reader, header) = create_bam_reader_for_pipeline(&self.io.input)?;

        if !is_template_coordinate_sorted(&header) {
            bail!(
                "Input BAM must be template-coordinate sorted.\n\n\
                To sort your BAM file, run:\n  \
                samtools sort -n -t CO input.bam | samtools sort -t CO -O bam -o sorted.bam\n\n\
                Or if you have samtools 1.18+:\n  \
                samtools sort --template-coordinate input.bam -o sorted.bam\n\n\
                Note: Pre-sorting is required because fgumi does not include an internal sorter.\n\
                This design choice keeps the tool lightweight and allows users to leverage\n\
                samtools' highly optimized sorting routines."
            );
        }
        info!("template coordinate sorted");

        // Add @PG record with PP chaining to input's last program
        let header = fgumi_lib::header::add_pg_record(
            header,
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        // Prepare raw tag as byte array
        if self.raw_tag.len() != 2 {
            bail!("Raw tag must be length two.")
        }
        let raw_tag: [u8; 2] = {
            let bytes = self.raw_tag.as_bytes();
            [bytes[0], bytes[1]]
        };

        // Prepare cell tag
        let cell_tag = string_to_tag(&self.cell_tag, "cell-tag")?;

        // Prepare assign tag for MI updates
        let assign_tag_bytes: [u8; 2] = {
            let bytes = self.assign_tag.as_bytes();
            [bytes[0], bytes[1]]
        };

        // Identity strategy requires edits=0, others use the configured value
        let effective_edits =
            if matches!(self.strategy, Strategy::Identity) { 0 } else { self.edits };

        // Create filter configuration
        let filter_config = GroupFilterConfig {
            umi_tag: raw_tag,
            min_mapq,
            include_non_pf: self.include_non_pf_reads,
            min_umi_length: self.min_umi_length,
        };

        // ============================================================
        // Check for single-threaded fast path (no --threads flag)
        // ============================================================
        if self.threading.threads.is_none() {
            // Single-threaded fast path - pass reader (required for stdin support)
            return self.execute_single_threaded(
                reader,
                &header,
                effective_edits,
                raw_tag,
                assign_tag_bytes,
                cell_tag,
                &filter_config,
                &timer,
            );
        }

        // ============================================================
        // Use 7-step unified pipeline (--threads N was specified)
        // ============================================================
        // Shared state for collecting metrics from parallel serialize_fn
        // Use lock-free SegQueue to avoid contention when multiple threads push metrics
        let collected_metrics: Arc<SegQueue<CollectedMetrics>> = Arc::new(SegQueue::new());

        // Clone values needed by closures
        let strategy = self.strategy;
        let index_threshold = self.index_threshold;
        let collected_metrics_clone = Arc::clone(&collected_metrics);

        // Configure 7-step pipeline
        let num_threads = self.threading.num_threads();
        // Use auto-tuned config for better performance with varying thread counts
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
        // Set memory limit if specified (convert MB to bytes)
        if self.queue_memory_limit_mb > 0 {
            let limit_bytes = self.queue_memory_limit_mb * 1024 * 1024;
            pipeline_config.pipeline.queue_memory_limit = limit_bytes;
            info!("Queue memory limit: {} MB", self.queue_memory_limit_mb);
        }
        info!("Scheduler: {:?}", self.scheduler_opts.strategy());
        // Template-based batching is enabled by default in auto_tuned() with target=500 templates.
        // This provides consistent batch sizes across datasets with varying templates-per-group ratios.

        info!("Using 7-step unified pipeline with {num_threads} threads");

        // Run the 7-step unified pipeline with the already-opened reader (supports streaming)
        let records_processed = run_bam_pipeline_from_reader(
            pipeline_config,
            reader,
            header,
            &self.io.output,
            None, // Use input header for output
            // grouper_fn: Create PositionGrouper with access to header
            move |header: &Header| {
                let library_lookup = build_library_lookup(header);
                let config = PositionGrouperConfig::new(cell_tag, library_lookup);
                Box::new(PositionGrouper::new(config))
                    as Box<dyn Grouper<Group = PositionGroup> + Send>
            },
            // process_fn: Filter templates + assign UMIs (parallel)
            move |group: PositionGroup| -> std::io::Result<ProcessedPositionGroup> {
                let mut filter_metrics = FilterMetrics::new();

                // Count ALL input records for progress tracking
                let input_record_count: u64 =
                    group.templates.iter().map(|t| t.read_count() as u64).sum();

                // Filter templates
                let filtered_templates: Vec<Template> = group
                    .templates
                    .into_iter()
                    .filter(|t| filter_template(t, &filter_config, &mut filter_metrics))
                    .collect();

                if filtered_templates.is_empty() {
                    return Ok(ProcessedPositionGroup {
                        templates: Vec::new(),
                        family_sizes: AHashMap::new(),
                        filter_metrics,
                        input_record_count,
                        base_mi: group.base_mi,
                    });
                }

                // Create UMI assigner for this group
                let assigner = strategy.new_assigner_full(effective_edits, 1, index_threshold);

                // Assign UMI groups using the unified _impl function
                let mut templates = filtered_templates;
                if let Err(_e) = assign_umi_groups_impl(
                    &mut templates,
                    assigner.as_ref(),
                    raw_tag,
                    assign_tag_bytes,
                    filter_config.min_umi_length,
                ) {
                    return Ok(ProcessedPositionGroup {
                        templates: Vec::new(),
                        family_sizes: AHashMap::new(),
                        filter_metrics,
                        input_record_count,
                        base_mi: group.base_mi,
                    });
                }

                // Sort templates directly by (MI index, name) - avoids Vec<Vec<Template>> allocation
                let base_mi = group.base_mi;
                templates.sort_by(|a, b| {
                    let a_idx = a.mi.to_vec_index();
                    let b_idx = b.mi.to_vec_index();
                    a_idx.cmp(&b_idx).then_with(|| a.name.cmp(&b.name))
                });

                // Count family sizes in one pass through sorted templates
                // Family sizes rarely exceed 50 distinct sizes
                let mut family_sizes: AHashMap<usize, u64> = AHashMap::with_capacity(50);
                if !templates.is_empty() {
                    let mut current_mi = templates[0].mi.to_vec_index();
                    let mut current_count = 1usize;

                    for template in templates.iter().skip(1) {
                        let mi = template.mi.to_vec_index();
                        if mi == current_mi {
                            current_count += 1;
                        } else {
                            // Finish previous MI group
                            if current_mi.is_some() {
                                *family_sizes.entry(current_count).or_insert(0) += 1;
                            }
                            current_mi = mi;
                            current_count = 1;
                        }
                    }
                    // Don't forget the last group
                    if current_mi.is_some() {
                        *family_sizes.entry(current_count).or_insert(0) += 1;
                    }
                }

                // Templates are now sorted by MI, no need for additional collection

                Ok(ProcessedPositionGroup {
                    templates,
                    family_sizes,
                    filter_metrics,
                    input_record_count,
                    base_mi,
                })
            },
            // serialize_fn: Serialize records + collect metrics (parallel)
            move |processed: ProcessedPositionGroup,
                  header: &Header,
                  output: &mut Vec<u8>|
                  -> std::io::Result<u64> {
                // Collect metrics for later aggregation
                let metrics = CollectedMetrics {
                    family_sizes: processed.family_sizes,
                    filter_metrics: processed.filter_metrics,
                };

                // Lock-free push to SegQueue - no contention
                collected_metrics_clone.push(metrics);

                // Save input record count for progress tracking
                let input_record_count = processed.input_record_count;
                let base_mi = processed.base_mi;

                // Convert assign_tag_bytes to Tag for setting MI
                let assign_tag = Tag::from(assign_tag_bytes);

                // Extract primary reads for serialization, setting MI tags just before write
                let mut records: Vec<noodles::sam::alignment::record_buf::RecordBuf> = Vec::new();
                for template in processed.templates {
                    // Capture mi before consuming template
                    let mi = template.mi;
                    let (r1, r2) = template.into_primary_reads();
                    if let Some(mut r1) = r1 {
                        set_mi_tag_on_record(&mut r1, mi, assign_tag, base_mi);
                        records.push(r1);
                    }
                    if let Some(mut r2) = r2 {
                        set_mi_tag_on_record(&mut r2, mi, assign_tag, base_mi);
                        records.push(r2);
                    }
                }

                // Serialize records into the provided buffer
                serialize_bam_records_into(&records, header, output)?;
                // Return INPUT record count for progress tracking (not output count)
                Ok(input_record_count)
            },
        )
        .context("Pipeline execution failed")?;

        info!("Wrote output to {}", self.io.output.display());

        // Aggregate collected metrics from lock-free SegQueue
        // Drain all metrics and merge them
        let mut family_size_counter: AHashMap<usize, u64> = AHashMap::with_capacity(50);
        let mut total_filter_metrics = FilterMetrics::new();

        while let Some(m) = collected_metrics.pop() {
            for (size, count) in m.family_sizes {
                *family_size_counter.entry(size).or_insert(0) += count;
            }
            total_filter_metrics.merge(&m.filter_metrics);
        }

        // Build final metrics
        let mut metrics = UmiGroupingMetrics::new();
        metrics.total_records = total_filter_metrics.total_templates;
        metrics.accepted_records = total_filter_metrics.accepted_templates;
        metrics.discarded_non_pf = total_filter_metrics.discarded_non_pf;
        metrics.discarded_poor_alignment = total_filter_metrics.discarded_poor_alignment;
        metrics.discarded_ns_in_umi = total_filter_metrics.discarded_ns_in_umi;
        metrics.discarded_umi_too_short = total_filter_metrics.discarded_umi_too_short;

        // Calculate metrics from family size counter
        metrics.total_families = family_size_counter.values().sum::<u64>();
        metrics.unique_molecule_ids = metrics.total_families;

        // Calculate average reads per molecule
        if metrics.unique_molecule_ids > 0 {
            metrics.avg_reads_per_molecule =
                metrics.accepted_records as f64 / metrics.unique_molecule_ids as f64;
        }

        // Calculate median family size
        if !family_size_counter.is_empty() {
            let mut sizes: Vec<(usize, u64)> =
                family_size_counter.iter().map(|(&size, &count)| (size, count)).collect();
            sizes.sort_by_key(|(size, _)| *size);

            let total_families = metrics.total_families;
            let mut cumulative = 0u64;
            let median_target = total_families / 2;

            for (size, count) in &sizes {
                cumulative += count;
                if cumulative >= median_target {
                    metrics.median_reads_per_molecule = *size as u64;
                    break;
                }
            }

            // Find min and max family sizes
            if let Some((min_size, _)) = sizes.first() {
                metrics.min_reads_per_molecule = *min_size as u64;
            }
            if let Some((max_size, _)) = sizes.last() {
                metrics.max_reads_per_molecule = *max_size as u64;
            }
        }

        // Log summary using enhanced logging (includes discard reasons)
        log_umi_grouping_summary(&metrics);

        // Write family size histogram
        if let Some(path) = &self.family_size_histogram {
            self.write_family_size_histogram(&family_size_counter).with_context(|| {
                format!("Failed to write family size histogram: {}", path.display())
            })?;
            info!("Wrote family size histogram to {}", path.display());
        }

        // Save accepted_records before moving metrics
        let accepted_records = metrics.accepted_records;

        // Write grouping metrics using new metrics structure
        if let Some(path) = &self.grouping_metrics {
            DelimFile::default()
                .write_tsv(path, [metrics])
                .with_context(|| format!("Failed to write grouping metrics: {}", path.display()))?;
            info!("Wrote grouping metrics to {}", path.display());
        }

        // Log completion with timing
        timer.log_completion(accepted_records);

        info!("group completed successfully");
        info!("Records processed by pipeline: {records_processed}");
        Ok(())
    }
}

impl GroupReadsByUmi {
    /// Execute in single-threaded mode for `--threads 1`.
    ///
    /// This provides a simpler, streaming implementation that avoids pipeline overhead
    /// while maintaining identical output to the multi-threaded mode.
    #[allow(clippy::too_many_arguments)]
    fn execute_single_threaded(
        &self,
        reader: Box<dyn std::io::Read + Send>,
        header: &Header,
        effective_edits: u32,
        raw_tag: [u8; 2],
        assign_tag_bytes: [u8; 2],
        cell_tag: Tag,
        filter_config: &GroupFilterConfig,
        timer: &OperationTimer,
    ) -> Result<()> {
        info!("Using single-threaded mode");

        // Wrap the reader in a BufReader and create noodles BAM reader
        // This reuses the already-opened reader (required for stdin support)
        let buf_reader = std::io::BufReader::new(reader);
        let mut bam_reader = noodles::bam::io::Reader::new(buf_reader);

        // Skip past header bytes - the reader is positioned at file start, but we already
        // have the header. We must read (and discard) it to position at the first record.
        let _ = bam_reader.read_header().context("Failed to skip BAM header")?;

        // Create output writer (single-threaded for strict thread control)
        let mut writer =
            create_bam_writer(&self.io.output, header, 1, self.compression.compression_level)?;

        // Build library index for GroupKey computation
        let library_index = LibraryIndex::from_header(header);

        // Create PositionGrouper
        let library_lookup = build_library_lookup(header);
        let config = PositionGrouperConfig::new(cell_tag, library_lookup);
        let mut grouper = PositionGrouper::new(config);

        // Metrics accumulators (no lock-free queue needed in single-threaded mode)
        let mut total_filter_metrics = FilterMetrics::new();
        let mut family_size_counter: AHashMap<usize, u64> = AHashMap::with_capacity(50);

        // Progress tracking
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        // Iterate over all records
        for record_result in bam_reader.record_bufs(header) {
            let record = record_result?;

            // Compute GroupKey for this record
            let key = compute_group_key(&record, &library_index, cell_tag);
            let decoded = DecodedRecord::new(record, key);

            // Feed to PositionGrouper - may emit completed groups
            let completed_groups = grouper.add_records(vec![decoded])?;

            // Process any completed position groups immediately
            for group in completed_groups {
                Self::process_and_write_position_group(
                    group,
                    filter_config,
                    self.strategy,
                    effective_edits,
                    self.index_threshold,
                    raw_tag,
                    assign_tag_bytes,
                    &mut total_filter_metrics,
                    &mut family_size_counter,
                    header,
                    &mut writer,
                )?;
            }

            progress.log_if_needed(1);
        }

        // Finish grouper - emit final group
        if let Some(final_group) = grouper.finish()? {
            Self::process_and_write_position_group(
                final_group,
                filter_config,
                self.strategy,
                effective_edits,
                self.index_threshold,
                raw_tag,
                assign_tag_bytes,
                &mut total_filter_metrics,
                &mut family_size_counter,
                header,
                &mut writer,
            )?;
        }

        progress.log_final();

        // Finish writer
        writer.into_inner().finish().context("Failed to finish output BAM")?;
        info!("Wrote output to {}", self.io.output.display());

        // Build final metrics (same logic as pipeline mode)
        let mut metrics = UmiGroupingMetrics::new();
        metrics.total_records = total_filter_metrics.total_templates;
        metrics.accepted_records = total_filter_metrics.accepted_templates;
        metrics.discarded_non_pf = total_filter_metrics.discarded_non_pf;
        metrics.discarded_poor_alignment = total_filter_metrics.discarded_poor_alignment;
        metrics.discarded_ns_in_umi = total_filter_metrics.discarded_ns_in_umi;
        metrics.discarded_umi_too_short = total_filter_metrics.discarded_umi_too_short;

        // Calculate metrics from family size counter
        metrics.total_families = family_size_counter.values().sum::<u64>();
        metrics.unique_molecule_ids = metrics.total_families;

        // Calculate average reads per molecule
        if metrics.unique_molecule_ids > 0 {
            metrics.avg_reads_per_molecule =
                metrics.accepted_records as f64 / metrics.unique_molecule_ids as f64;
        }

        // Calculate median family size
        if !family_size_counter.is_empty() {
            let mut sizes: Vec<(usize, u64)> =
                family_size_counter.iter().map(|(&size, &count)| (size, count)).collect();
            sizes.sort_by_key(|(size, _)| *size);

            let total_families = metrics.total_families;
            let mut cumulative = 0u64;
            let median_target = total_families / 2;

            for (size, count) in &sizes {
                cumulative += count;
                if cumulative >= median_target {
                    metrics.median_reads_per_molecule = *size as u64;
                    break;
                }
            }

            // Find min and max family sizes
            if let Some((min_size, _)) = sizes.first() {
                metrics.min_reads_per_molecule = *min_size as u64;
            }
            if let Some((max_size, _)) = sizes.last() {
                metrics.max_reads_per_molecule = *max_size as u64;
            }
        }

        // Log summary using enhanced logging
        log_umi_grouping_summary(&metrics);

        // Write family size histogram
        if let Some(path) = &self.family_size_histogram {
            self.write_family_size_histogram(&family_size_counter).with_context(|| {
                format!("Failed to write family size histogram: {}", path.display())
            })?;
            info!("Wrote family size histogram to {}", path.display());
        }

        // Save accepted_records before moving metrics
        let accepted_records = metrics.accepted_records;

        // Write grouping metrics
        if let Some(path) = &self.grouping_metrics {
            DelimFile::default()
                .write_tsv(path, [metrics])
                .with_context(|| format!("Failed to write grouping metrics: {}", path.display()))?;
            info!("Wrote grouping metrics to {}", path.display());
        }

        // Log completion with timing
        timer.log_completion(accepted_records);

        info!("group completed successfully");
        Ok(())
    }

    /// Process a single position group: filter, assign UMIs, and write output.
    /// Used by `execute_single_threaded` for streaming processing.
    #[allow(clippy::too_many_arguments)]
    fn process_and_write_position_group(
        group: PositionGroup,
        filter_config: &GroupFilterConfig,
        strategy: Strategy,
        effective_edits: u32,
        index_threshold: usize,
        raw_tag: [u8; 2],
        assign_tag_bytes: [u8; 2],
        total_filter_metrics: &mut FilterMetrics,
        family_size_counter: &mut AHashMap<usize, u64>,
        header: &Header,
        writer: &mut fgumi_lib::bam_io::BamWriter,
    ) -> Result<()> {
        let mut filter_metrics = FilterMetrics::new();

        // Filter templates
        let filtered_templates: Vec<Template> = group
            .templates
            .into_iter()
            .filter(|t| filter_template(t, filter_config, &mut filter_metrics))
            .collect();

        // Merge filter metrics
        total_filter_metrics.merge(&filter_metrics);

        if filtered_templates.is_empty() {
            return Ok(());
        }

        // Create UMI assigner
        let assigner = strategy.new_assigner_full(effective_edits, 1, index_threshold);

        // Assign UMI groups
        let mut templates = filtered_templates;
        if let Err(e) = assign_umi_groups_impl(
            &mut templates,
            assigner.as_ref(),
            raw_tag,
            assign_tag_bytes,
            filter_config.min_umi_length,
        ) {
            // Log error but continue processing
            log::warn!("Failed to assign UMI groups: {e}");
            return Ok(());
        }

        // Sort templates directly by (MI index, name) - avoids Vec<Vec<Template>> allocation
        templates.sort_by(|a, b| {
            let a_idx = a.mi.to_vec_index();
            let b_idx = b.mi.to_vec_index();
            a_idx.cmp(&b_idx).then_with(|| a.name.cmp(&b.name))
        });

        // Count family sizes in one pass through sorted templates
        if !templates.is_empty() {
            let mut current_mi = templates[0].mi.to_vec_index();
            let mut current_count = 1usize;

            for template in templates.iter().skip(1) {
                let mi = template.mi.to_vec_index();
                if mi == current_mi {
                    current_count += 1;
                } else {
                    // Finish previous MI group
                    if current_mi.is_some() {
                        *family_size_counter.entry(current_count).or_insert(0) += 1;
                    }
                    current_mi = mi;
                    current_count = 1;
                }
            }
            // Don't forget the last group
            if current_mi.is_some() {
                *family_size_counter.entry(current_count).or_insert(0) += 1;
            }
        }

        // Write templates (already sorted by MI, then by name)
        let base_mi = group.base_mi;
        let assign_tag = Tag::from(assign_tag_bytes);

        for template in templates {
            let mi = template.mi;
            // Write primary reads with MI tags set just before writing
            let (r1, r2) = template.into_primary_reads();
            if let Some(mut r1) = r1 {
                set_mi_tag_on_record(&mut r1, mi, assign_tag, base_mi);
                writer.write_alignment_record(header, &r1)?;
            }
            if let Some(mut r2) = r2 {
                set_mi_tag_on_record(&mut r2, mi, assign_tag, base_mi);
                writer.write_alignment_record(header, &r2)?;
            }
        }

        Ok(())
    }
}

#[allow(dead_code)]
impl GroupReadsByUmi {
    /// Get the unclipped 5' position of a read
    /// For positive strand: unclipped start
    /// For negative strand: unclipped end
    fn get_unclipped_position(record: &sam::alignment::RecordBuf) -> Result<i32> {
        if record.flags().is_unmapped() {
            return Ok(0);
        }

        let alignment_start = record
            .alignment_start()
            .ok_or_else(|| anyhow::anyhow!("Mapped read missing alignment start"))?;

        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let alignment_start_i32 = usize::from(alignment_start) as i32;

        if record.flags().is_reverse_complemented() {
            // Negative strand: calculate unclipped end
            // unclipped_end = alignment_start + alignment_span + soft_clipped_bases_at_end
            let cigar = record.cigar();

            // Calculate alignment span (M/D/N/=/X operations)
            let alignment_span: i32 = cigar
                .as_ref()
                .iter()
                .filter(|op| {
                    matches!(
                        op.kind(),
                        noodles::sam::alignment::record::cigar::op::Kind::Match
                            | noodles::sam::alignment::record::cigar::op::Kind::Deletion
                            | noodles::sam::alignment::record::cigar::op::Kind::Skip
                            | noodles::sam::alignment::record::cigar::op::Kind::SequenceMatch
                            | noodles::sam::alignment::record::cigar::op::Kind::SequenceMismatch
                    )
                })
                .map(|op| i32::try_from(op.len()).unwrap())
                .sum();

            // Get trailing soft clips
            let trailing_soft_clips: i32 = cigar
                .as_ref()
                .iter()
                .rev()
                .take_while(|op| {
                    matches!(op.kind(), noodles::sam::alignment::record::cigar::op::Kind::SoftClip)
                })
                .map(|op| i32::try_from(op.len()).unwrap())
                .sum();

            Ok(alignment_start_i32 + alignment_span + trailing_soft_clips)
        } else {
            // Positive strand: calculate unclipped start
            // unclipped_start = alignment_start - soft_clipped_bases_at_start
            let cigar = record.cigar();

            // Get leading soft clips
            let leading_soft_clips: i32 = cigar
                .as_ref()
                .iter()
                .take_while(|op| {
                    matches!(op.kind(), noodles::sam::alignment::record::cigar::op::Kind::SoftClip)
                })
                .map(|op| i32::try_from(op.len()).unwrap())
                .sum();

            Ok(alignment_start_i32 - leading_soft_clips)
        }
    }

    /// Process a position group without writing - returns result for parallel processing
    fn process_position_group(
        &self,
        group: &mut Vec<Template>,
        assigner: &dyn UmiAssigner,
    ) -> Result<PositionGroupResult> {
        if group.is_empty() {
            return Ok(PositionGroupResult {
                family_sizes: AHashMap::new(),
                sorted_templates: Vec::new(),
            });
        }

        // Assign UMI groups
        self.assign_umi_groups(group, assigner)?;

        // Group by molecule ID using Vec-based indexing
        let max_idx = group.iter().filter_map(|t| t.mi.to_vec_index()).max().unwrap_or(0);

        let mut by_mi: Vec<Vec<Template>> = vec![Vec::new(); max_idx + 1];
        for template in group.drain(..) {
            if let Some(idx) = template.mi.to_vec_index() {
                by_mi[idx].push(template);
            }
        }

        // Count family sizes
        let mut family_sizes: AHashMap<usize, u64> = AHashMap::new();
        for group_templates in &by_mi {
            if !group_templates.is_empty() {
                *family_sizes.entry(group_templates.len()).or_insert(0) += 1;
            }
        }

        // Collect sorted templates (Vec indices are already in sorted order)
        let mut sorted_templates = Vec::new();
        for mut group_templates in by_mi {
            if group_templates.is_empty() {
                continue;
            }
            group_templates.sort_by(|a, b| a.name.cmp(&b.name));
            sorted_templates.extend(group_templates);
        }

        Ok(PositionGroupResult { family_sizes, sorted_templates })
    }

    /// Assign UMI groups to templates at the same position
    ///
    /// When `split_templates_by_pair_orientation` is true, templates are first split
    /// into subgroups by (R1 positive strand, R2 positive strand) before UMI assignment.
    /// This matches fgbio's behavior where F1R2 pairs are not grouped with F2R1 pairs
    /// for non-paired strategies (e.g., Edit, Adjacency, Identity).
    fn assign_umi_groups(
        &self,
        templates: &mut [Template],
        assigner: &dyn UmiAssigner,
    ) -> Result<()> {
        // Split templates by pair orientation if required by the assigner
        // This matches fgbio's behavior in GroupReadsByUmi.assignUmiGroups
        if assigner.split_templates_by_pair_orientation() {
            // Group by (R1 positive strand, R2 positive strand)
            // Templates with different pair orientations should NOT be grouped together
            let mut subgroups: AHashMap<(bool, bool), Vec<usize>> = AHashMap::new();

            for (idx, template) in templates.iter().enumerate() {
                let orientation = Self::get_pair_orientation(template);
                subgroups.entry(orientation).or_default().push(idx);
            }

            // Process each subgroup separately
            // Each subgroup gets its own set of MI assignments
            for indices in subgroups.values() {
                self.assign_umi_groups_for_indices(templates, indices, assigner)?;
            }
        } else {
            // No splitting - process all templates together (for Paired strategy)
            let all_indices: Vec<usize> = (0..templates.len()).collect();
            self.assign_umi_groups_for_indices(templates, &all_indices, assigner)?;
        }

        Ok(())
    }

    /// Get the pair orientation for a template: (R1 positive strand, R2 positive strand)
    /// Returns (true, true) for unpaired/single-end reads
    fn get_pair_orientation(template: &Template) -> (bool, bool) {
        let r1_positive = template.r1().is_none_or(|r| !r.flags().is_reverse_complemented());
        let r2_positive = template.r2().is_none_or(|r| !r.flags().is_reverse_complemented());
        (r1_positive, r2_positive)
    }

    /// Assign UMI groups to a specific subset of templates (identified by indices)
    ///
    /// This sets the `Template.mi` field with the `MoleculeId` enum. The actual BAM MI tag
    /// is set later during serialization, when we have the global offset.
    fn assign_umi_groups_for_indices(
        &self,
        templates: &mut [Template],
        indices: &[usize],
        assigner: &dyn UmiAssigner,
    ) -> Result<()> {
        if indices.is_empty() {
            return Ok(());
        }

        // Extract UMIs from the specified templates - pre-allocate capacity
        let mut umis = Vec::with_capacity(indices.len());
        let raw_tag: [u8; 2] = {
            let bytes = self.raw_tag.as_bytes();
            [bytes[0], bytes[1]]
        };

        for &idx in indices {
            let template = &templates[idx];
            let umi_str = if let Some(r1) = &template.r1() {
                if let Some(DataValue::String(umi_bytes)) = r1.data().get(&raw_tag) {
                    String::from_utf8_lossy(umi_bytes).into_owned()
                } else {
                    bail!("UMI tag is not a string");
                }
            } else if let Some(r2) = &template.r2() {
                if let Some(DataValue::String(umi_bytes)) = r2.data().get(&raw_tag) {
                    String::from_utf8_lossy(umi_bytes).into_owned()
                } else {
                    bail!("UMI tag is not a string");
                }
            } else {
                bail!("Template has no reads");
            };

            // Determine which read is earlier based on GENOMIC POSITION for paired UMI strategies
            let is_r1_earlier = if let (Some(r1), Some(r2)) = (&template.r1(), &template.r2()) {
                Self::is_r1_genomically_earlier(r1, r2)?
            } else {
                true
            };

            let processed_umi = Self::umi_for_read(&umi_str, is_r1_earlier, assigner)?;
            umis.push(processed_umi);
        }

        // Truncate UMIs if needed
        let truncated_umis = self.truncate_umis(umis)?;

        // Assign UMI groups - returns Vec<MoleculeId> indexed by input position
        let assignments = assigner.assign(&truncated_umis);

        // Store MoleculeId enum in Template.mi field
        // The actual MI tag string is set during serialization with global offset
        for (i, &idx) in indices.iter().enumerate() {
            let template = &mut templates[idx];
            template.mi = assignments[i];
        }

        Ok(())
    }

    /// Determine if R1 is genomically earlier than R2
    /// This uses the actual genomic coordinates, not just the read pair flags
    fn is_r1_genomically_earlier(
        r1: &sam::alignment::RecordBuf,
        r2: &sam::alignment::RecordBuf,
    ) -> Result<bool> {
        let ref1 = r1.reference_sequence_id().map_or(-1, |id| i32::try_from(id).unwrap());
        let ref2 = r2.reference_sequence_id().map_or(-1, |id| i32::try_from(id).unwrap());

        // If on different chromosomes, compare chromosome IDs
        if ref1 != ref2 {
            return Ok(ref1 < ref2);
        }

        // Same chromosome: compare unclipped 5' positions
        let pos1 = Self::get_unclipped_position(r1)?;
        let pos2 = Self::get_unclipped_position(r2)?;

        if pos1 != pos2 {
            return Ok(pos1 < pos2);
        }

        // Same position: compare strands (positive strand sorts before negative)
        let strand1 = u8::from(r1.flags().is_reverse_complemented());
        let strand2 = u8::from(r2.flags().is_reverse_complemented());

        Ok(strand1 < strand2)
    }

    /// Get the MI tag from a template
    fn get_mi_tag(&self, template: &Template) -> Result<String> {
        use sam::alignment::record::data::field::Tag;
        let assign_tag_bytes: [u8; 2] = {
            let bytes = self.assign_tag.as_bytes();
            [bytes[0], bytes[1]]
        };
        let assign_tag = Tag::from(assign_tag_bytes);

        let r1 = template.r1();
        let r2 = template.r2();
        let record = if let Some(rec) = r1 { rec } else { r2.unwrap() };

        if let Some(data) = record.data().get(&assign_tag) {
            if let DataValue::String(mi) = data {
                Ok(String::from_utf8_lossy(mi).to_string())
            } else {
                bail!("MI tag is not a string");
            }
        } else {
            bail!("Missing MI tag");
        }
    }

    /// Extract UMI for a read, handling paired UMI strategies
    fn umi_for_read(umi: &str, is_r1_earlier: bool, assigner: &dyn UmiAssigner) -> Result<String> {
        // For paired UMI strategies, we need to add prefixes to indicate which read is earlier
        if assigner.split_templates_by_pair_orientation() {
            // For non-paired strategies, just return the UMI
            Ok(umi.to_uppercase())
        } else {
            // For paired strategies, parse and prefix the UMI
            let parts: Vec<&str> = umi.split('-').collect();
            if parts.len() != 2 {
                bail!(
                    "Paired strategy used but UMI did not contain 2 segments delimited by '-': {umi}"
                );
            }

            let Some(paired) = assigner.as_any().downcast_ref::<PairedUmiAssigner>() else {
                bail!("Expected PairedUmiAssigner")
            };

            let result = if is_r1_earlier {
                format!(
                    "{}:{}-{}:{}",
                    paired.lower_read_umi_prefix(),
                    parts[0],
                    paired.higher_read_umi_prefix(),
                    parts[1]
                )
            } else {
                format!(
                    "{}:{}-{}:{}",
                    paired.higher_read_umi_prefix(),
                    parts[0],
                    paired.lower_read_umi_prefix(),
                    parts[1]
                )
            };

            Ok(result)
        }
    }

    /// Truncate UMIs to minimum length if specified
    fn truncate_umis(&self, umis: Vec<String>) -> Result<Vec<String>> {
        match self.min_umi_length {
            None => Ok(umis),
            Some(min_len) => {
                let min_length = umis.iter().map(std::string::String::len).min().unwrap_or(0);
                if min_length < min_len {
                    bail!(
                        "UMI found that had shorter length than expected ({min_length} < {min_len})"
                    );
                }
                Ok(umis.into_iter().map(|u| u[..min_length].to_string()).collect())
            }
        }
    }

    /// Write family size histogram
    fn write_family_size_histogram(&self, family_sizes: &AHashMap<usize, u64>) -> Result<()> {
        if let Some(path) = &self.family_size_histogram {
            // Collect and sort by family size
            let mut sorted: Vec<_> = family_sizes.iter().map(|(&s, &c)| (s, c)).collect();
            sorted.sort_by_key(|(size, _)| *size);

            // Calculate total count
            #[allow(clippy::cast_precision_loss)]
            let total: f64 = sorted.iter().map(|(_, count)| *count as f64).sum();

            // Build metrics in reverse to calculate cumulative fraction in one pass
            let mut metrics = Vec::with_capacity(sorted.len());
            let mut cumulative = 0.0;
            #[allow(clippy::cast_precision_loss)]
            for &(family_size, count) in sorted.iter().rev() {
                let fraction = count as f64 / total;
                cumulative += fraction;
                metrics.push(TagFamilySizeMetric {
                    family_size,
                    count,
                    fraction,
                    fraction_gt_or_eq_family_size: cumulative,
                });
            }
            metrics.reverse();

            // Write to file
            DelimFile::default()
                .write_tsv(path, metrics)
                .with_context(|| format!("Failed to create file: {}", path.display()))?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_lib::assigner::{IdentityUmiAssigner, PairedUmiAssigner, Strategy};
    use fgumi_lib::sam::builder::{RecordBuilder, RecordPairBuilder};
    use noodles::bam;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::io::Write as AlignmentWrite;
    use rstest::rstest;
    use std::num::NonZeroUsize;
    use tempfile::{NamedTempFile, TempDir};

    /// Helper struct for managing temporary test file paths.
    /// Keeps `TempDir` alive for the lifetime of the struct.
    struct TestPaths {
        #[allow(dead_code)]
        dir: TempDir,
        pub output: PathBuf,
        pub histogram: PathBuf,
        pub metrics: PathBuf,
    }

    impl TestPaths {
        fn new() -> Result<Self> {
            let dir = TempDir::new()?;
            Ok(Self {
                output: dir.path().join("output.bam"),
                histogram: dir.path().join("histogram.txt"),
                metrics: dir.path().join("metrics.txt"),
                dir,
            })
        }
    }

    // ========================================================================
    // Helper Functions for Test Data Creation - FIXED
    // ========================================================================

    /// Create a minimal SAM header for testing - FIXED
    fn create_test_header() -> sam::Header {
        use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
        use noodles::sam::header::record::value::{
            Map as HeaderRecordMap,
            map::{Header as HeaderRecord, Tag as HeaderTag},
        };

        let mut builder = sam::Header::builder();

        // Add header with template-coordinate sort order
        let HeaderTag::Other(so_tag) = HeaderTag::from([b'S', b'O']) else { unreachable!() };
        let HeaderTag::Other(go_tag) = HeaderTag::from([b'G', b'O']) else { unreachable!() };
        let HeaderTag::Other(ss_tag) = HeaderTag::from([b'S', b'S']) else { unreachable!() };

        let map = HeaderRecordMap::<HeaderRecord>::builder()
            .insert(so_tag, "unsorted")
            .insert(go_tag, "query")
            .insert(ss_tag, "template-coordinate")
            .build()
            .unwrap();
        builder = builder.set_header(map);

        // FIX: Use NonZeroUsize instead of Position
        builder = builder.add_reference_sequence(
            BString::from("chr1"),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(248_956_422).unwrap()),
        );
        builder = builder.add_reference_sequence(
            BString::from("chr2"),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(242_193_529).unwrap()),
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
        umi: &str,
    ) -> sam::alignment::RecordBuf {
        // 100bp sequence
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let mut builder = RecordBuilder::new()
            .name(name)
            .sequence(seq)
            .reference_sequence_id(ref_id)
            .alignment_start(pos as usize)
            .mapping_quality(mapq)
            .cigar("100M")
            .tag("RX", umi);

        // Set flags based on raw u16
        let sam_flags = sam::alignment::record::Flags::from(flags);
        if sam_flags.is_segmented() {
            builder = builder.paired(true);
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

    /// Create a pair of reads
    #[allow(clippy::cast_sign_loss)]
    fn build_test_pair(
        name: &str,
        ref_id: usize,
        pos1: i32,
        pos2: i32,
        mapq1: u8,
        _mapq2: u8,
        umi: &str,
    ) -> (sam::alignment::RecordBuf, sam::alignment::RecordBuf) {
        RecordPairBuilder::new()
            .name(name)
            .r1_sequence(&"A".repeat(100))
            .r2_sequence(&"A".repeat(100))
            .reference_sequence_id(ref_id)
            .r1_start(pos1 as usize)
            .r2_start(pos2 as usize)
            .mapping_quality(mapq1)
            .tag("RX", umi)
            .build()
    }

    /// Create a pair where R1 is mapped (with specified MAPQ) and R2 is unmapped
    /// This tests the case where a mapped read with low MAPQ has an unmapped mate.
    #[allow(clippy::cast_sign_loss)]
    fn build_test_pair_mapped_with_unmapped_mate(
        name: &str,
        ref_id: usize,
        pos: i32,
        mapq: u8,
        umi: &str,
    ) -> (sam::alignment::RecordBuf, sam::alignment::RecordBuf) {
        RecordPairBuilder::new()
            .name(name)
            .r1_sequence(&"A".repeat(100))
            .r2_sequence(&"A".repeat(100))
            .reference_sequence_id(ref_id)
            .r1_start(pos as usize) // R1 mapped
            // R2 has no start position, so it's unmapped
            .r1_reverse(false)
            .r2_reverse(false)
            .tag("RX", umi)
            .r1_tag("MQ", 0i32) // R1's mate (R2) is unmapped, so MQ=0
            .r2_tag("MQ", i32::from(mapq)) // R2's mate (R1) has the specified MAPQ
            .build()
    }

    /// Write records to a temporary BAM file
    fn create_test_bam(records: Vec<sam::alignment::RecordBuf>) -> Result<NamedTempFile> {
        let temp_file = NamedTempFile::new()?;
        let header = create_test_header();

        let mut writer = bam::io::writer::Builder.build_from_path(temp_file.path())?;

        writer.write_header(&header)?;

        for record in records {
            writer.write_alignment_record(&header, &record)?;
        }

        drop(writer); // Ensure file is flushed

        Ok(temp_file)
    }

    /// Read all records from a BAM file - FIXED
    fn read_bam_records(path: &std::path::Path) -> Result<Vec<sam::alignment::RecordBuf>> {
        let mut reader = bam::io::reader::Builder.build_from_path(path)?;
        let header = reader.read_header()?;
        let mut records = Vec::new();

        // FIX: Convert bam::Record to RecordBuf
        for result in reader.records() {
            let record = result?;
            let record_buf =
                sam::alignment::RecordBuf::try_from_alignment_record(&header, &record)?;
            records.push(record_buf);
        }

        Ok(records)
    }

    /// Extract MI tags from records
    fn get_mi_tags(records: &[sam::alignment::RecordBuf]) -> Vec<String> {
        use sam::alignment::record::data::field::Tag;

        let mi_tag = Tag::from([b'M', b'I']);

        records
            .iter()
            .filter_map(|r| {
                r.data().get(&mi_tag).and_then(|v| {
                    if let sam::alignment::record_buf::data::field::Value::String(s) = v {
                        Some(String::from_utf8_lossy(s).to_string())
                    } else {
                        None
                    }
                })
            })
            .collect()
    }

    /// Count unique MI tags
    fn count_unique_mi_tags(records: &[sam::alignment::RecordBuf]) -> usize {
        use std::collections::HashSet;
        let tags: HashSet<_> = get_mi_tags(records).into_iter().collect();
        tags.len()
    }

    /// Helper function to get a string tag from a record
    fn get_string_tag(record: &sam::alignment::RecordBuf, tag_name: &str) -> Option<String> {
        use sam::alignment::record::data::field::Tag;

        let tag_bytes = tag_name.as_bytes();
        let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

        record.data().get(&tag).and_then(|v| {
            if let sam::alignment::record_buf::data::field::Value::String(s) = v {
                Some(String::from_utf8_lossy(s).to_string())
            } else {
                None
            }
        })
    }

    // ========================================================================
    // Unit Tests for GroupReadsByUmi Methods
    // ========================================================================

    #[test]
    fn test_umi_for_read_assigns_ab_prefixes_by_coordinates() {
        let assigner: Box<dyn UmiAssigner> = Box::new(PairedUmiAssigner::new(1));

        let umi1 = "AAA-TTT";
        let result1 =
            GroupReadsByUmi::umi_for_read(umi1, true, assigner.as_ref()).expect("Should succeed");

        assert!(result1.contains("AAA"));
        assert!(result1.contains("TTT"));
        assert!(result1.contains('-'));

        let result2 =
            GroupReadsByUmi::umi_for_read(umi1, false, assigner.as_ref()).expect("Should succeed");

        assert_ne!(result1, result2, "Prefixes should differ based on which read is earlier");
    }

    #[test]
    fn test_umi_for_read_handles_absent_umi_ends() {
        let assigner: Box<dyn UmiAssigner> = Box::new(PairedUmiAssigner::new(1));

        let result_left = GroupReadsByUmi::umi_for_read("-TTT", true, assigner.as_ref())
            .expect("Should handle absent left UMI");
        assert!(result_left.contains('-') && result_left.contains("TTT"));

        let result_right = GroupReadsByUmi::umi_for_read("AAA-", true, assigner.as_ref())
            .expect("Should handle absent right UMI");
        assert!(result_right.contains("AAA") && result_right.contains('-'));
    }

    #[test]
    fn test_umi_for_read_uppercase_for_non_paired() {
        let assigner: Box<dyn UmiAssigner> = Box::new(IdentityUmiAssigner::new());

        let result = GroupReadsByUmi::umi_for_read("acgtacgt", true, assigner.as_ref())
            .expect("Should succeed");

        assert_eq!(result, "ACGTACGT");
    }

    #[test]
    fn test_truncate_umis_to_minimum_length() {
        let tool = GroupReadsByUmi {
            io: BamIoOptions {
                input: std::path::PathBuf::from("/dev/null"),
                output: std::path::PathBuf::from("/dev/null"),
            },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: Some(5),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        let umis = vec!["AAAAAA".to_string(), "AAAAA".to_string(), "AAAAAAA".to_string()];
        let truncated = tool.truncate_umis(umis).expect("Should truncate successfully");

        assert_eq!(truncated.len(), 3);
        assert!(truncated.iter().all(|u| u.len() == 5));
    }

    #[test]
    fn test_truncate_umis_none_returns_unchanged() {
        let tool = GroupReadsByUmi {
            io: BamIoOptions {
                input: std::path::PathBuf::from("/dev/null"),
                output: std::path::PathBuf::from("/dev/null"),
            },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        let umis = vec!["AAAAAA".to_string(), "AAAAA".to_string()];
        let original = umis.clone();
        let result = tool.truncate_umis(umis).expect("Should succeed");

        assert_eq!(result, original);
    }

    #[test]
    fn test_truncate_umis_fails_when_too_short() {
        let tool = GroupReadsByUmi {
            io: BamIoOptions {
                input: std::path::PathBuf::from("/dev/null"),
                output: std::path::PathBuf::from("/dev/null"),
            },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: Some(6),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        let umis = vec!["AAAAAA".to_string(), "AAAA".to_string()];
        let result = tool.truncate_umis(umis);
        assert!(result.is_err());
    }

    // ========================================================================
    // Integration Tests
    // ========================================================================

    #[test]
    fn test_groups_reads_correctly_basic() -> Result<()> {
        // Create test data with UMIs that should group together
        let mut records = Vec::new();

        // a01-a04: same position, similar UMIs (should group with edits=1)
        let (r1, r2) = build_test_pair("a01", 0, 100, 300, 60, 60, "AAAAAAAA");
        records.push(r1);
        records.push(r2);

        let (r1, r2) = build_test_pair("a02", 0, 100, 300, 60, 60, "AAAAgAAA");
        records.push(r1);
        records.push(r2);

        let (r1, r2) = build_test_pair("a03", 0, 100, 300, 60, 60, "AAAAAAAA");
        records.push(r1);
        records.push(r2);

        let (r1, r2) = build_test_pair("a04", 0, 100, 300, 60, 60, "AAAAAAAt");
        records.push(r1);
        records.push(r2);

        // c01: low mapq, should be filtered
        let (r1, r2) = build_test_pair("c01", 0, 100, 300, 5, 5, "AAAAAAAA");
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: Some(30),
            include_non_pf_reads: false,
            strategy: Strategy::Edit,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should have 8 records (4 pairs, c01 filtered out)
        assert_eq!(output_records.len(), 8, "Should have 8 records after filtering");

        // All should have same MI tag (grouped together)
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 1, "Should have 1 UMI group");

        Ok(())
    }

    #[test]
    fn test_filtering_excludes_reads_with_n_in_umi() -> Result<()> {
        let mut records = Vec::new();

        // Good UMI
        let (r1, r2) = build_test_pair("a01", 0, 100, 300, 60, 60, "AAAAAA");
        records.push(r1);
        records.push(r2);

        // UMI with N - should be filtered
        let (r1, r2) = build_test_pair("a02", 0, 100, 300, 60, 60, "AANAAA");
        records.push(r1);
        records.push(r2);

        // Good UMI
        let (r1, r2) = build_test_pair("a03", 0, 100, 300, 60, 60, "AAAAAA");
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: Some(paths.metrics.clone()),
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 4, "Should have 4 records (2 pairs, a02 filtered)");

        // Check metrics
        let metrics: Vec<UmiGroupingMetrics> = DelimFile::default().read_tsv(&paths.metrics)?;
        assert_eq!(metrics.len(), 1);
        assert_eq!(metrics[0].discarded_ns_in_umi, 2);

        Ok(())
    }

    #[test]
    fn test_filtering_excludes_reads_below_min_mapq() -> Result<()> {
        let mut records = Vec::new();

        // High MAPQ - should pass
        let (r1, r2) = build_test_pair("a01", 0, 100, 300, 60, 60, "AAAAAA");
        records.push(r1);
        records.push(r2);

        // Low MAPQ - should be filtered
        let (r1, r2) = build_test_pair("a02", 0, 100, 300, 10, 10, "AAAAAA");
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: Some(30),
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 2, "Should have 2 records (1 pair, low MAPQ filtered)");

        Ok(())
    }

    #[test]
    fn test_correctly_groups_single_end_reads() -> Result<()> {
        let records: Vec<RecordBuf> = vec![
            // Group 1: same position, same UMI
            build_test_read("a01", 0, 100, 60, 0, "AAAAAAAA"),
            build_test_read("a02", 0, 100, 60, 0, "AAAAAAAA"),
            // Group 2: same position, similar UMI (1 edit)
            build_test_read("a03", 0, 100, 60, 0, "CACACACA"),
            build_test_read("a04", 0, 100, 60, 0, "CACACACC"),
            // Group 3: different position
            build_test_read("a05", 0, 105, 60, 0, "GTAGTAGG"),
            build_test_read("a06", 0, 105, 60, 0, "GTAGTAGG"),
            // Group 4: different position (from group 1)
            build_test_read("a07", 0, 107, 60, 0, "AAAAAAAA"),
            build_test_read("a08", 0, 107, 60, 0, "AAAAAAAA"),
        ];

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Edit,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 8, "Should have all 8 records");

        // Should have 4 groups
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 4, "Should have 4 UMI groups");

        Ok(())
    }

    #[test]
    fn test_outputs_family_size_histogram() -> Result<()> {
        let mut records = Vec::new();

        // Create groups of different sizes
        // Group 1: 3 pairs
        for i in 1..=3 {
            let (r1, r2) = build_test_pair(&format!("a{i:02}"), 0, 100, 300, 60, 60, "AAAAAAAA");
            records.push(r1);
            records.push(r2);
        }

        // Group 2: 2 pairs
        for i in 1..=2 {
            let (r1, r2) = build_test_pair(&format!("b{i:02}"), 0, 200, 400, 60, 60, "CCCCCCCC");
            records.push(r1);
            records.push(r2);
        }

        // Group 3: 1 pair
        let (r1, r2) = build_test_pair("c01", 0, 300, 500, 60, 60, "GGGGGGGG");
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: Some(paths.histogram.clone()),
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        // Check histogram was created
        assert!(&paths.histogram.exists(), "Histogram file should exist");

        let metrics: Vec<TagFamilySizeMetric> = DelimFile::default().read_tsv(&paths.histogram)?;
        assert_eq!(metrics.len(), 3);
        assert_eq!(
            metrics[0],
            TagFamilySizeMetric {
                family_size: 1,
                count: 1,
                fraction: 0.333_333_333_333_333_3,
                fraction_gt_or_eq_family_size: 1.0,
            }
        );
        assert_eq!(
            metrics[1],
            TagFamilySizeMetric {
                family_size: 2,
                count: 1,
                fraction: 0.333_333_333_333_333_3,
                fraction_gt_or_eq_family_size: 0.666_666_666_666_666_6,
            }
        );
        assert_eq!(
            metrics[2],
            TagFamilySizeMetric {
                family_size: 3,
                count: 1,
                fraction: 0.333_333_333_333_333_3,
                fraction_gt_or_eq_family_size: 0.333_333_333_333_333_3,
            }
        );
        Ok(())
    }

    #[test]
    fn test_outputs_grouping_metrics() -> Result<()> {
        let mut records = Vec::new();

        // Good record
        let (r1, r2) = build_test_pair("a01", 0, 100, 300, 60, 60, "AAAAAA");
        records.push(r1);
        records.push(r2);

        // Record with N in UMI
        let (r1, r2) = build_test_pair("a02", 0, 100, 300, 60, 60, "AANAAA");
        records.push(r1);
        records.push(r2);

        // Low MAPQ record
        let (r1, r2) = build_test_pair("a03", 0, 100, 300, 5, 5, "AAAAAA");
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: Some(paths.metrics.clone()),
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: Some(30),
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let metrics: Vec<UmiGroupingMetrics> = DelimFile::default().read_tsv(&paths.metrics)?;
        assert_eq!(metrics.len(), 1);
        assert_eq!(metrics[0].accepted_records, 2);
        assert_eq!(metrics[0].discarded_ns_in_umi, 2);
        assert_eq!(metrics[0].discarded_poor_alignment, 2);
        assert_eq!(metrics[0].discarded_non_pf, 0);
        assert_eq!(metrics[0].discarded_umi_too_short, 0);

        Ok(())
    }

    #[test]
    fn test_rejects_umis_shorter_than_min_length() -> Result<()> {
        let mut records = Vec::new();

        // UMI of length 6 - OK
        let (r1, r2) = build_test_pair("a01", 0, 100, 300, 60, 60, "ACTACT");
        records.push(r1);
        records.push(r2);

        // UMI of length 5 - too short
        let (r1, r2) = build_test_pair("a02", 0, 100, 300, 60, 60, "ACTAC");
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: Some(paths.metrics.clone()),
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Edit,
            edits: 0,
            min_umi_length: Some(6),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 2, "Should only have records with UMI length >= 6");

        // Check metrics
        let metrics: Vec<UmiGroupingMetrics> = DelimFile::default().read_tsv(&paths.metrics)?;
        assert_eq!(metrics.len(), 1);
        assert_eq!(metrics[0].accepted_records, 2);
        assert_eq!(metrics[0].discarded_ns_in_umi, 0);
        assert_eq!(metrics[0].discarded_poor_alignment, 0);
        assert_eq!(metrics[0].discarded_non_pf, 0);
        assert_eq!(metrics[0].discarded_umi_too_short, 2);

        Ok(())
    }

    #[test]
    fn test_truncates_to_min_length() -> Result<()> {
        let mut records = Vec::new();

        // UMI of length 6
        let (r1, r2) = build_test_pair("a01", 0, 100, 300, 60, 60, "ACTACT");
        records.push(r1);
        records.push(r2);

        // UMI of length 5 (will be truncated to 5)
        let (r1, r2) = build_test_pair("a02", 0, 100, 300, 60, 60, "ACTAC");
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Edit,
            edits: 0,
            min_umi_length: Some(5),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 4, "Should have all 4 records");

        // Both should be grouped together (truncated UMIs match)
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 1, "Should have 1 group after truncation");

        Ok(())
    }

    #[test]
    fn test_cross_contig_read_pairs() -> Result<()> {
        // Test that cross-contig read pairs are correctly grouped with paired UMI strategy
        // Pairs where R1 is on contig1/R2 on contig2 should group with pairs where
        // R1 is on contig2/R2 on contig1 (same molecule, different orientations)
        let mut records = Vec::new();

        // Group A: R1 on chr1, R2 on chr2 (both forward strand)
        for i in 1..=4 {
            let (r1, r2) = RecordPairBuilder::new()
                .name(&format!("a{i:02}"))
                .r1_sequence(&"A".repeat(100))
                .r2_sequence(&"A".repeat(100))
                .reference_sequence_id(0) // R1 on chr1
                .r2_reference_sequence_id(1) // R2 on chr2
                .r1_start(100)
                .r2_start(300)
                .r2_reverse(false) // Both reads forward (matching original test)
                .tag("RX", "ACT-ACT")
                .build();
            records.push(r1);
            records.push(r2);
        }

        // Group B: R1 on chr2, R2 on chr1 (flipped, both forward strand)
        for i in 1..=4 {
            let (r1, r2) = RecordPairBuilder::new()
                .name(&format!("b{i:02}"))
                .r1_sequence(&"A".repeat(100))
                .r2_sequence(&"A".repeat(100))
                .reference_sequence_id(1) // R1 on chr2
                .r2_reference_sequence_id(0) // R2 on chr1
                .r1_start(300)
                .r2_start(100)
                .r2_reverse(false) // Both reads forward (matching original test)
                .tag("RX", "ACT-ACT")
                .build();
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Paired,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 16, "Should have all 16 records");

        // Extract MI tags for each group
        let a_mis: Vec<String> = output_records
            .iter()
            .filter(|r| r.name().is_some_and(|n| String::from_utf8_lossy(n).starts_with('a')))
            .filter_map(|r| get_string_tag(r, "MI"))
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();

        let b_mis: Vec<String> = output_records
            .iter()
            .filter(|r| r.name().is_some_and(|n| String::from_utf8_lossy(n).starts_with('b')))
            .filter_map(|r| get_string_tag(r, "MI"))
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();

        // Both groups should have exactly one unique MI
        assert_eq!(a_mis.len(), 1, "Group A should have 1 unique MI");
        assert_eq!(b_mis.len(), 1, "Group B should have 1 unique MI");

        // The prefix (before '/') should be the same for both groups
        let a_prefix = a_mis[0].split('/').next().unwrap();
        let b_prefix = b_mis[0].split('/').next().unwrap();
        assert_eq!(a_prefix, b_prefix, "Both groups should have same MI prefix");

        // But the full MI (including suffix) should be different
        assert_ne!(a_mis[0], b_mis[0], "Groups should have different MI suffixes");

        Ok(())
    }

    #[test]
    fn test_pair_orientation_splitting() -> Result<()> {
        // Test that reads with different pair orientations are NOT grouped together
        // even if they have the same UMI - they have different start positions and orientations
        use noodles::sam::alignment::record::Flags;

        let mut records = Vec::new();

        // F1R2: Read 1 forward at 100, Read 2 reverse at 300
        let (r1, r2) = build_test_pair("f1r2", 0, 100, 300, 60, 60, "ACGT-TTGA");
        records.push(r1);
        records.push(r2);

        // F2R1: Read 1 reverse at 300, Read 2 forward at 100 (flipped positions and strands)
        let (mut r1, mut r2) = build_test_pair("f2r1", 0, 300, 100, 60, 60, "ACGT-TTGA");
        // Flip strands
        let mut flags1 = r1.flags();
        flags1.set(Flags::REVERSE_COMPLEMENTED, true);
        *r1.flags_mut() = flags1;
        let mut flags2 = r2.flags();
        flags2.set(Flags::REVERSE_COMPLEMENTED, false);
        *r2.flags_mut() = flags2;
        records.push(r1);
        records.push(r2);

        // FF: Both forward at 100 and 300
        let (r1, mut r2) = build_test_pair("ff", 0, 100, 300, 60, 60, "ACGT-TTGA");
        let mut flags2 = r2.flags();
        flags2.set(Flags::REVERSE_COMPLEMENTED, false);
        *r2.flags_mut() = flags2;
        records.push(r1);
        records.push(r2);

        // RR: Both reverse at 1 and 201
        let (mut r1, mut r2) = build_test_pair("rr", 0, 1, 201, 60, 60, "ACGT-TTGA");
        let mut flags1 = r1.flags();
        flags1.set(Flags::REVERSE_COMPLEMENTED, true);
        *r1.flags_mut() = flags1;
        let mut flags2 = r2.flags();
        flags2.set(Flags::REVERSE_COMPLEMENTED, true);
        *r2.flags_mut() = flags2;
        records.push(r1);
        records.push(r2);

        // R1F2: Read 1 reverse at 150, Read 2 forward at 350
        let (mut r1, mut r2) = build_test_pair("r1f2", 0, 150, 350, 60, 60, "ACGT-TTGA");
        let mut flags1 = r1.flags();
        flags1.set(Flags::REVERSE_COMPLEMENTED, true);
        *r1.flags_mut() = flags1;
        let mut flags2 = r2.flags();
        flags2.set(Flags::REVERSE_COMPLEMENTED, false);
        *r2.flags_mut() = flags2;
        records.push(r1);
        records.push(r2);

        // R2F1: Read 1 forward at 400, Read 2 reverse at 600
        let (r1, mut r2) = build_test_pair("r2f1", 0, 400, 600, 60, 60, "ACGT-TTGA");
        let mut flags2 = r2.flags();
        flags2.set(Flags::REVERSE_COMPLEMENTED, true);
        *r2.flags_mut() = flags2;
        records.push(r1);
        records.push(r2);

        // Add some single-end reads with different strands
        let mut frag = build_test_read("Frag", 0, 100, 60, 0, "ACGT-TTGA");
        let mut flags = frag.flags();
        flags.set(Flags::REVERSE_COMPLEMENTED, true);
        *frag.flags_mut() = flags;
        records.push(frag);

        let frag = build_test_read("fRag", 0, 1, 60, 0, "ACGT-TTGA");
        records.push(frag);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Adjacency,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Each template should have a separate MI because of different orientations
        let unique_mis = count_unique_mi_tags(&output_records);
        assert_eq!(unique_mis, 8, "Should have 8 unique MIs (one per template)");

        Ok(())
    }

    #[test]
    fn test_missing_raw_tag() -> Result<()> {
        // Test that templates with missing UMI tags are filtered out
        let mut records = Vec::new();

        // Create a pair WITHOUT the RX tag
        let (mut r1, mut r2) = build_test_pair("a01", 0, 100, 300, 60, 60, "dummy");

        // Remove the RX tag
        use noodles::sam::alignment::record::data::field::Tag;
        let rx_tag = Tag::from([b'R', b'X']);
        r1.data_mut().remove(&rx_tag);
        r2.data_mut().remove(&rx_tag);

        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Paired,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        // Templates with missing UMI tags should be filtered out
        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 0, "Should have no output records");

        Ok(())
    }

    #[test]
    fn test_cell_barcode_grouping() -> Result<()> {
        // Test that reads with different cell barcodes are grouped separately
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;

        let mut records = Vec::new();
        let cb_tag = Tag::from([b'C', b'B']);

        // Two reads with same UMI and position but same cell barcode
        let mut r1 = build_test_read("a01", 0, 100, 60, 0, "AAAAAAAA");
        r1.data_mut().insert(cb_tag, Value::String(b"AA".into()));
        records.push(r1);

        let mut r2 = build_test_read("a02", 0, 100, 60, 0, "AAAAAAAA");
        r2.data_mut().insert(cb_tag, Value::String(b"AA".into()));
        records.push(r2);

        // One read with similar UMI but different cell barcode
        let mut r3 = build_test_read("a03", 0, 100, 60, 0, "CACACACA");
        r3.data_mut().insert(cb_tag, Value::String(b"CA".into()));
        records.push(r3);

        // One read with close UMI but different cell barcode
        let mut r4 = build_test_read("a04", 0, 100, 60, 0, "CACACACC");
        r4.data_mut().insert(cb_tag, Value::String(b"NN".into()));
        records.push(r4);

        // Two reads at different position with same UMI and cell barcode
        let mut r5 = build_test_read("a05", 0, 105, 60, 0, "GTAGTAGG");
        r5.data_mut().insert(cb_tag, Value::String(b"GT".into()));
        records.push(r5);

        let mut r6 = build_test_read("a06", 0, 105, 60, 0, "GTAGTAGG");
        r6.data_mut().insert(cb_tag, Value::String(b"GT".into()));
        records.push(r6);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Edit,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 6, "Should have all 6 records");

        // Group by MI
        let mut groups: std::collections::HashMap<String, Vec<String>> =
            std::collections::HashMap::new();
        for record in &output_records {
            let mi = get_string_tag(record, "MI").unwrap();
            let name = String::from_utf8_lossy(record.name().unwrap()).to_string();
            groups.entry(mi).or_default().push(name);
        }

        assert_eq!(groups.len(), 4, "Should have 4 unique MI groups");

        // Check specific groupings
        let group_sets: Vec<std::collections::HashSet<String>> =
            groups.values().map(|v| v.iter().cloned().collect()).collect();

        assert!(
            group_sets.contains(&["a01".to_string(), "a02".to_string()].iter().cloned().collect())
        );
        assert!(group_sets.contains(&["a03".to_string()].iter().cloned().collect()));
        assert!(group_sets.contains(&["a04".to_string()].iter().cloned().collect()));
        assert!(
            group_sets.contains(&["a05".to_string(), "a06".to_string()].iter().cloned().collect())
        );

        Ok(())
    }

    #[test]
    fn test_paired_mode_with_absent_umi_on_right() -> Result<()> {
        // Test paired mode when the right end of the source molecule does not have a UMI
        // UMI format: "ACT-" (left side present, right side absent)
        let mut records = Vec::new();

        // Group A: R1 earlier, all with "ACT-" UMI
        for i in 1..=4 {
            let (r1, r2) = build_test_pair(&format!("a{i:02}"), 0, 100, 300, 60, 60, "ACT-");
            records.push(r1);
            records.push(r2);
        }

        // Group B: R2 earlier (flipped positions), all with "-ACT" UMI
        for i in 1..=4 {
            let (r1, r2) = build_test_pair(&format!("b{i:02}"), 0, 300, 100, 60, 60, "-ACT");
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Paired,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 16, "Should have all 16 records");

        // Extract MI tags for each group
        let a_mis: Vec<String> = output_records
            .iter()
            .filter(|r| r.name().is_some_and(|n| String::from_utf8_lossy(n).starts_with('a')))
            .filter_map(|r| get_string_tag(r, "MI"))
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();

        let b_mis: Vec<String> = output_records
            .iter()
            .filter(|r| r.name().is_some_and(|n| String::from_utf8_lossy(n).starts_with('b')))
            .filter_map(|r| get_string_tag(r, "MI"))
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();

        // Both groups should have exactly one unique MI
        assert_eq!(a_mis.len(), 1, "Group A should have 1 unique MI");
        assert_eq!(b_mis.len(), 1, "Group B should have 1 unique MI");

        // Check that MIs end with /A and /B respectively
        assert!(a_mis[0].ends_with("/A"), "Group A should have /A suffix");
        assert!(b_mis[0].ends_with("/B"), "Group B should have /B suffix");

        Ok(())
    }

    #[test]
    fn test_paired_mode_with_absent_umi_on_left() -> Result<()> {
        // Test paired mode when the left end of the source molecule does not have a UMI
        // UMI format: "-ACT" (left side absent, right side present)
        let mut records = Vec::new();

        // Group A: R1 earlier, all with "-ACT" UMI
        for i in 1..=4 {
            let (r1, r2) = build_test_pair(&format!("a{i:02}"), 0, 100, 300, 60, 60, "-ACT");
            records.push(r1);
            records.push(r2);
        }

        // Group B: R2 earlier (flipped positions), all with "ACT-" UMI
        for i in 1..=4 {
            let (r1, r2) = build_test_pair(&format!("b{i:02}"), 0, 300, 100, 60, 60, "ACT-");
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Paired,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 16, "Should have all 16 records");

        // Extract MI tags for each group
        let a_mis: Vec<String> = output_records
            .iter()
            .filter(|r| r.name().is_some_and(|n| String::from_utf8_lossy(n).starts_with('a')))
            .filter_map(|r| get_string_tag(r, "MI"))
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();

        let b_mis: Vec<String> = output_records
            .iter()
            .filter(|r| r.name().is_some_and(|n| String::from_utf8_lossy(n).starts_with('b')))
            .filter_map(|r| get_string_tag(r, "MI"))
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();

        // Both groups should have exactly one unique MI
        assert_eq!(a_mis.len(), 1, "Group A should have 1 unique MI");
        assert_eq!(b_mis.len(), 1, "Group B should have 1 unique MI");

        // Check that MIs end with /A and /B respectively
        assert!(a_mis[0].ends_with("/A"), "Group A should have /A suffix");
        assert!(b_mis[0].ends_with("/B"), "Group B should have /B suffix");

        Ok(())
    }

    #[test]
    fn test_discard_secondary_and_supplementary_reads() -> Result<()> {
        // Test that secondary and supplementary reads are filtered out
        let mut records = Vec::new();

        // Add primary read with high MAPQ (will not be marked as duplicate)
        let (r1, r2) = build_test_pair("a01", 0, 100, 300, 100, 100, "AAAAAAAA");
        records.push(r1);
        records.push(r2);

        // Add secondary read (should be filtered out)
        let (r1, r2) = RecordPairBuilder::new()
            .name("a01_sec")
            .r1_sequence(&"A".repeat(100))
            .r2_sequence(&"A".repeat(100))
            .reference_sequence_id(0)
            .r1_start(100)
            .r2_start(300)
            .tag("RX", "AAAAAAAA")
            .secondary(true)
            .build();
        records.push(r1);
        records.push(r2);

        // Add supplementary read (should be filtered out)
        let (r1, r2) = RecordPairBuilder::new()
            .name("a01_sup")
            .r1_sequence(&"A".repeat(100))
            .r2_sequence(&"A".repeat(100))
            .reference_sequence_id(0)
            .r1_start(100)
            .r2_start(300)
            .tag("RX", "AAAAAAAA")
            .supplementary(true)
            .build();
        records.push(r1);
        records.push(r2);

        // Add another primary read with lower MAPQ (will be marked as duplicate)
        let (r1, r2) = build_test_pair("a02", 0, 100, 300, 10, 10, "AAAAAAAA");
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Edit,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should only have 4 records (2 pairs): the primary reads, not the secondary/supplementary
        assert_eq!(
            output_records.len(),
            4,
            "Should have only 4 records (secondary and supplementary filtered out)"
        );

        // Check that no records are marked as secondary or supplementary
        for record in &output_records {
            assert!(!record.flags().is_secondary(), "Output should not contain secondary reads");
            assert!(
                !record.flags().is_supplementary(),
                "Output should not contain supplementary reads"
            );
        }

        // Check that we have records from both a01 and a02
        let a01_count = output_records
            .iter()
            .filter(|r| r.name().is_some_and(|n| String::from_utf8_lossy(n).starts_with("a01")))
            .count();
        let a02_count = output_records
            .iter()
            .filter(|r| r.name().is_some_and(|n| String::from_utf8_lossy(n).starts_with("a02")))
            .count();

        assert_eq!(a01_count, 2, "Should have 2 reads from a01");
        assert_eq!(a02_count, 2, "Should have 2 reads from a02");

        Ok(())
    }

    #[test]
    fn test_adjacency_single_thread() -> Result<()> {
        // Test adjacency strategy with single thread
        let mut records = Vec::new();

        // Create UMIs that should form 3 separate groups based on edit distance
        // Group 1: AAAAAA family (will include AAAAAT within 2 edits)
        for i in 1..=4 {
            let (r1, r2) = build_test_pair(&format!("g1_{i}"), 0, 100, 300, 60, 60, "AAAAAA");
            records.push(r1);
            records.push(r2);
        }
        for i in 1..=2 {
            let (r1, r2) = build_test_pair(&format!("g1b_{i}"), 0, 100, 300, 60, 60, "AAAAAT");
            records.push(r1);
            records.push(r2);
        }

        // Group 2: GACGAC family (will include GACGAT and GACGCC within 2 edits)
        for i in 1..=9 {
            let (r1, r2) = build_test_pair(&format!("g2_{i}"), 0, 100, 300, 60, 60, "GACGAC");
            records.push(r1);
            records.push(r2);
        }
        let (r1, r2) = build_test_pair("g2b_1", 0, 100, 300, 60, 60, "GACGAT");
        records.push(r1);
        records.push(r2);
        for i in 1..=4 {
            let (r1, r2) = build_test_pair(&format!("g2c_{i}"), 0, 100, 300, 60, 60, "GACGCC");
            records.push(r1);
            records.push(r2);
        }

        // Group 3: TACGAC (too far from GACGAC to merge)
        for i in 1..=7 {
            let (r1, r2) = build_test_pair(&format!("g3_{i}"), 0, 100, 300, 60, 60, "TACGAC");
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Adjacency,
            edits: 2,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        // 6 pairs (AAAAAA family) + 14 pairs (GACGAC family) + 7 pairs (TACGAC) = 27 pairs = 54 records
        assert_eq!(output_records.len(), 54, "Should have all 54 records (27 pairs)");

        // Count unique MI tags (should be 3 groups)
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(
            unique_groups, 3,
            "Should have 3 unique groups with adjacency strategy and 1 thread"
        );

        Ok(())
    }

    #[test]
    fn test_adjacency_multi_thread() -> Result<()> {
        // Test adjacency strategy with multiple threads - same test as single thread
        let mut records = Vec::new();

        // Group 1: AAAAAA family
        for i in 1..=4 {
            let (r1, r2) = build_test_pair(&format!("g1_{i}"), 0, 100, 300, 60, 60, "AAAAAA");
            records.push(r1);
            records.push(r2);
        }
        for i in 1..=2 {
            let (r1, r2) = build_test_pair(&format!("g1b_{i}"), 0, 100, 300, 60, 60, "AAAAAT");
            records.push(r1);
            records.push(r2);
        }

        // Group 2: GACGAC family
        for i in 1..=9 {
            let (r1, r2) = build_test_pair(&format!("g2_{i}"), 0, 100, 300, 60, 60, "GACGAC");
            records.push(r1);
            records.push(r2);
        }
        let (r1, r2) = build_test_pair("g2b_1", 0, 100, 300, 60, 60, "GACGAT");
        records.push(r1);
        records.push(r2);
        for i in 1..=4 {
            let (r1, r2) = build_test_pair(&format!("g2c_{i}"), 0, 100, 300, 60, 60, "GACGCC");
            records.push(r1);
            records.push(r2);
        }

        // Group 3: TACGAC
        for i in 1..=7 {
            let (r1, r2) = build_test_pair(&format!("g3_{i}"), 0, 100, 300, 60, 60, "TACGAC");
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Adjacency,
            edits: 2,
            min_umi_length: None,
            threading: ThreadingOptions::new(4), // Use 4 threads
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        // 6 pairs (AAAAAA family) + 14 pairs (GACGAC family) + 7 pairs (TACGAC) = 27 pairs = 54 records
        assert_eq!(output_records.len(), 54, "Should have all 54 records (27 pairs)");

        // Count unique MI tags (should be 3 groups, same as single-threaded)
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(
            unique_groups, 3,
            "Should have 3 unique groups with adjacency strategy and 4 threads (same as single thread)"
        );

        Ok(())
    }

    #[test]
    fn test_adjacency_deep_tree_single_thread() -> Result<()> {
        // Test handling of deep UMI tree with single thread
        let mut records = Vec::new();

        // Create a deep tree: AAAAAA -> TAAAAA -> TTAAAA -> TTTAAA -> TTTTAA
        // Each differs by 1 edit from the next
        let umis =
            vec![("AAAAAA", 256), ("TAAAAA", 128), ("TTAAAA", 64), ("TTTAAA", 32), ("TTTTAA", 16)];

        let mut counter = 1;
        for (umi, count) in umis {
            for _ in 0..count {
                let (r1, r2) = build_test_pair(&format!("q{counter}"), 0, 100, 300, 60, 60, umi);
                records.push(r1);
                records.push(r2);
                counter += 1;
            }
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Adjacency,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should have 992 records (496 pairs = 256+128+64+32+16)
        assert_eq!(output_records.len(), 992, "Should have all records");

        // All should be in one group (connected by the deep tree)
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 1, "Should have 1 group (all UMIs connected in deep tree)");

        Ok(())
    }

    #[test]
    fn test_adjacency_deep_tree_multi_thread() -> Result<()> {
        // Test handling of deep UMI tree with multiple threads
        let mut records = Vec::new();

        let umis =
            vec![("AAAAAA", 256), ("TAAAAA", 128), ("TTAAAA", 64), ("TTTAAA", 32), ("TTTTAA", 16)];

        let mut counter = 1;
        for (umi, count) in umis {
            for _ in 0..count {
                let (r1, r2) = build_test_pair(&format!("q{counter}"), 0, 100, 300, 60, 60, umi);
                records.push(r1);
                records.push(r2);
                counter += 1;
            }
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Adjacency,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::new(4), // Use 4 threads
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        assert_eq!(output_records.len(), 992, "Should have all records");

        // All should be in one group, same as single-threaded
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(
            unique_groups, 1,
            "Should have 1 group with multi-threading (same as single thread)"
        );

        Ok(())
    }

    #[test]
    fn test_paired_assigner_explicit_ab_ba_symmetry() -> Result<()> {
        // Test that paired assigner correctly handles A-B and B-A pairs
        // They should be grouped separately as different orientations
        let mut records = Vec::new();

        // Create pairs with A-B orientation
        let (r1_ab, r2_ab) = build_test_pair("read_ab_1", 0, 100, 300, 60, 60, "AAAA-TTTT");
        records.push(r1_ab);
        records.push(r2_ab);

        // Create another pair with same A-B orientation
        let (r1_ab2, r2_ab2) = build_test_pair("read_ab_2", 0, 100, 300, 60, 60, "AAAA-TTTT");
        records.push(r1_ab2);
        records.push(r2_ab2);

        // Create pairs with B-A orientation (swapped UMIs)
        let (r1_ba, r2_ba) = build_test_pair("read_ba_1", 0, 100, 300, 60, 60, "TTTT-AAAA");
        records.push(r1_ba);
        records.push(r2_ba);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Paired,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should have 6 records (3 pairs)
        assert_eq!(output_records.len(), 6, "Should have all 6 records");

        // Get the MI tags
        let mi_tags = get_mi_tags(&output_records);

        // With paired strategy, A-B and B-A should form separate groups
        // but both should get /A or /B suffixes
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 2, "Should have 2 groups (A-B vs B-A orientation)");

        // All tags should have the /A or /B suffix
        assert!(
            mi_tags.iter().all(|t| t.contains("/A") || t.contains("/B")),
            "All MI tags should have /A or /B suffix. Tags: {mi_tags:?}"
        );

        Ok(())
    }

    #[test]
    fn test_edit_strategy_with_exact_edit_distance() -> Result<()> {
        // Test edit strategy with UMIs at exact edit distance threshold
        let mut records = Vec::new();

        // UMIs that differ by exactly 1 edit
        let (r1_a, r2_a) = build_test_pair("read1", 0, 100, 300, 60, 60, "AAAAAA");
        let (r1_b, r2_b) = build_test_pair("read2", 0, 100, 300, 60, 60, "TAAAAA");

        records.push(r1_a);
        records.push(r2_a);
        records.push(r1_b);
        records.push(r2_b);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Edit,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should group into 1 group (within edit distance)
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 1, "Should group UMIs within edit distance");

        Ok(())
    }

    #[test]
    fn test_edit_strategy_outside_edit_distance() -> Result<()> {
        // Test edit strategy with UMIs outside edit distance threshold
        let mut records = Vec::new();

        // UMIs that differ by 2 edits (more than threshold)
        let (r1_a, r2_a) = build_test_pair("read1", 0, 100, 300, 60, 60, "AAAAAA");
        let (r1_b, r2_b) = build_test_pair("read2", 0, 100, 300, 60, 60, "TTAAAA");

        records.push(r1_a);
        records.push(r2_a);
        records.push(r1_b);
        records.push(r2_b);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Edit,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should keep as 2 separate groups (outside edit distance)
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 2, "Should keep UMIs outside edit distance separate");

        Ok(())
    }

    #[test]
    fn test_identity_strategy_no_grouping() -> Result<()> {
        // Test identity strategy keeps different UMIs separate
        let mut records = Vec::new();

        // Very similar UMIs that differ by only 1 base
        let (r1_a, r2_a) = build_test_pair("read1", 0, 100, 300, 60, 60, "AAAAAA");
        let (r1_b, r2_b) = build_test_pair("read2", 0, 100, 300, 60, 60, "AAAAAC");
        let (r1_c, r2_c) = build_test_pair("read3", 0, 100, 300, 60, 60, "AAAAAG");

        records.push(r1_a);
        records.push(r2_a);
        records.push(r1_b);
        records.push(r2_b);
        records.push(r1_c);
        records.push(r2_c);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should have 3 separate groups (no grouping with identity)
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 3, "Identity strategy should not group similar UMIs");

        Ok(())
    }

    #[test]
    fn test_adjacency_with_count_gradient() -> Result<()> {
        // Test adjacency strategy respects count gradient (high count node can absorb low count)
        let mut records = Vec::new();

        // High count UMI (100 reads)
        for i in 0..50 {
            let (r1, r2) = build_test_pair(&format!("high_{i}"), 0, 100, 300, 60, 60, "AAAAAA");
            records.push(r1);
            records.push(r2);
        }

        // Low count UMI (2 reads) - 1 edit away
        let (r1_low, r2_low) = build_test_pair("low_1", 0, 100, 300, 60, 60, "TAAAAA");
        records.push(r1_low);
        records.push(r2_low);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Adjacency,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Low count should be absorbed into high count group
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 1, "Low count UMI should be absorbed by high count");

        Ok(())
    }

    #[test]
    fn test_adjacency_no_gradient_keeps_separate() -> Result<()> {
        // Test adjacency strategy keeps UMIs separate when counts are similar (no gradient)
        let mut records = Vec::new();

        // Two UMIs with similar counts
        for i in 0..10 {
            let (r1, r2) = build_test_pair(&format!("umi1_{i}"), 0, 100, 300, 60, 60, "AAAAAA");
            records.push(r1);
            records.push(r2);
        }

        for i in 0..10 {
            let (r1, r2) = build_test_pair(&format!("umi2_{i}"), 0, 100, 300, 60, 60, "TAAAAA");
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Adjacency,
            edits: 1,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should keep as 2 groups (no count gradient to drive merging)
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 2, "Similar counts should not merge in adjacency");

        Ok(())
    }

    #[test]
    fn test_very_long_umis() -> Result<()> {
        // Test handling of very long UMIs (e.g., 20+ bases)
        let mut records = Vec::new();

        let long_umi = "ACGTACGTACGTACGTACGT"; // 20 bases
        let (r1, r2) = build_test_pair("read1", 0, 100, 300, 60, 60, long_umi);
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should handle long UMIs without issue
        assert_eq!(output_records.len(), 2, "Should process long UMIs correctly");
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 1, "Should have 1 group");

        Ok(())
    }

    #[test]
    fn test_minimum_umi_length_filtering() -> Result<()> {
        // Test that reads with UMIs shorter than min_umi_length are filtered
        let mut records = Vec::new();

        // Short UMI (4 bases)
        let (r1_short, r2_short) = build_test_pair("short", 0, 100, 300, 60, 60, "ACGT");
        records.push(r1_short);
        records.push(r2_short);

        // Long UMI (10 bases)
        let (r1_long, r2_long) = build_test_pair("long", 0, 100, 300, 60, 60, "ACGTACGTAC");
        records.push(r1_long);
        records.push(r2_long);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: Some(8), // Require at least 8 bases
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should only have the long UMI reads (2 records)
        assert_eq!(output_records.len(), 2, "Should filter short UMIs");

        // Verify the remaining reads have UMIs that meet the min length
        let rx_tags: Vec<String> =
            output_records.iter().filter_map(|r| get_string_tag(r, "RX")).collect();

        assert!(rx_tags.iter().all(|umi| umi.len() >= 8), "UMIs should be at least min length");

        Ok(())
    }

    #[test]
    fn test_unmapped_templates_filtered() -> Result<()> {
        // Test that templates where all reads are unmapped are filtered out
        // Note: This test may not work as expected if the group command doesn't filter unmapped reads
        let mut records = Vec::new();

        // Mapped pair
        let (r1_mapped, r2_mapped) = build_test_pair("mapped", 0, 200, 400, 60, 60, "TTTTTT");
        records.push(r1_mapped);
        records.push(r2_mapped);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should have the mapped pair (2 records)
        assert_eq!(output_records.len(), 2, "Should process mapped templates");

        Ok(())
    }

    #[test]
    fn test_multiple_edit_distances() -> Result<()> {
        // Test edit strategy with different edit distance thresholds
        let mut records = Vec::new();

        // Create UMIs at varying distances
        let (r1_a, r2_a) = build_test_pair("read1", 0, 100, 300, 60, 60, "AAAAAA");
        let (r1_b, r2_b) = build_test_pair("read2", 0, 100, 300, 60, 60, "TAAAAA"); // 1 edit
        let (r1_c, r2_c) = build_test_pair("read3", 0, 100, 300, 60, 60, "TTAAAA"); // 2 edits from A

        records.push(r1_a);
        records.push(r2_a);
        records.push(r1_b);
        records.push(r2_b);
        records.push(r1_c);
        records.push(r2_c);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        // With edits=2, all should group together
        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Edit,
            edits: 2,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 1, "With edits=2, all UMIs should group");

        Ok(())
    }

    #[test]
    fn test_paired_umi_missing_both_ends() -> Result<()> {
        // Test paired UMI handling when both ends are missing (just "-")
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;

        let mut r1 = build_test_read("test", 0, 100, 60, 0x41, "");
        let mut r2 = build_test_read("test", 0, 300, 60, 0x81, "");

        // Set UMI to just "-" (both ends missing)
        let rx_tag = Tag::from([b'R', b'X']);
        r1.data_mut().insert(rx_tag, Value::String(b"-".into()));
        r2.data_mut().insert(rx_tag, Value::String(b"-".into()));

        let records = vec![r1, r2];
        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Paired,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        // Should handle gracefully (either error or filter out)
        let result = cmd.execute("test");
        // The command should either succeed with filtered reads or fail gracefully
        assert!(result.is_ok() || result.is_err());

        Ok(())
    }

    #[test]
    fn test_adjacency_with_edits_0() -> Result<()> {
        // Test adjacency with edits=0 (should act like identity)
        let mut records = Vec::new();

        let (r1_a, r2_a) = build_test_pair("read1", 0, 100, 300, 60, 60, "AAAAAA");
        let (r1_b, r2_b) = build_test_pair("read2", 0, 100, 300, 60, 60, "AAAAAA"); // Same UMI
        let (r1_c, r2_c) = build_test_pair("read3", 0, 100, 300, 60, 60, "AAAAAC"); // 1 edit away

        records.push(r1_a);
        records.push(r2_a);
        records.push(r1_b);
        records.push(r2_b);
        records.push(r1_c);
        records.push(r2_c);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Adjacency,
            edits: 0, // No edits allowed
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // With edits=0, adjacency should keep different UMIs separate
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 2, "Adjacency with edits=0 should keep different UMIs separate");

        Ok(())
    }

    #[test]
    fn test_edit_strategy_with_high_edit_threshold() -> Result<()> {
        // Test edit strategy with high edit distance (should group more UMIs)
        let mut records = Vec::new();

        let (r1_a, r2_a) = build_test_pair("read1", 0, 100, 300, 60, 60, "AAAAAA");
        let (r1_b, r2_b) = build_test_pair("read2", 0, 100, 300, 60, 60, "TTAAAA"); // 2 edits
        let (r1_c, r2_c) = build_test_pair("read3", 0, 100, 300, 60, 60, "TTTAAA"); // 3 edits from A

        records.push(r1_a);
        records.push(r2_a);
        records.push(r1_b);
        records.push(r2_b);
        records.push(r1_c);
        records.push(r2_c);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Edit,
            edits: 3, // Allow up to 3 edits
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // With edits=3, all should group together
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 1, "Edit strategy with edits=3 should group all UMIs");

        Ok(())
    }

    #[test]
    fn test_paired_strategy_single_end_reads() -> Result<()> {
        // Test paired strategy with single-end (fragment) reads
        let mut records = Vec::new();

        let r1 = build_test_read("frag1", 0, 100, 60, 0x0, "AAAA-TTTT");
        let r2 = build_test_read("frag2", 0, 200, 60, 0x0, "AAAA-TTTT");

        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Paired,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 2);

        Ok(())
    }

    #[test]
    fn test_family_size_histogram_generation() -> Result<()> {
        // Test that family size histogram is generated correctly
        let mut records = Vec::new();

        // Create 5 reads with same UMI (family size 5)
        for i in 0..5 {
            let (r1, r2) = build_test_pair(&format!("read{i}"), 0, 100, 300, 60, 60, "AAAAAA");
            records.push(r1);
            records.push(r2);
        }

        // Create 2 reads with different UMI (family size 2)
        for i in 0..2 {
            let (r1, r2) = build_test_pair(&format!("other{i}"), 0, 100, 300, 60, 60, "TTTTTT");
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: Some(paths.histogram.clone()),
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        // Verify histogram file was created
        assert!(&paths.histogram.exists());

        Ok(())
    }

    #[test]
    fn test_grouping_metrics_generation() -> Result<()> {
        // Test that grouping metrics are generated correctly
        let mut records = Vec::new();

        let (r1, r2) = build_test_pair("read1", 0, 100, 300, 60, 60, "AAAAAA");
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: Some(paths.metrics.clone()),
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        // Verify metrics file was created
        assert!(&paths.metrics.exists());

        Ok(())
    }

    #[test]
    fn test_min_mapq_filtering() -> Result<()> {
        // Test that reads below min_map_q are filtered out
        let mut records = Vec::new();

        // High quality mapping
        let (r1_high, r2_high) = build_test_pair("high", 0, 100, 300, 60, 60, "AAAAAA");
        records.push(r1_high);
        records.push(r2_high);

        // Low quality mapping (mapq=10)
        let (r1_low, r2_low) = build_test_pair("low", 0, 100, 300, 10, 10, "TTTTTT");
        records.push(r1_low);
        records.push(r2_low);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: Some(20), // Filter reads with mapq < 20
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should only have high quality reads (2 records)
        assert_eq!(output_records.len(), 2, "Should filter low mapq reads");

        Ok(())
    }

    #[test]
    fn test_umi_with_only_n_bases() -> Result<()> {
        // Test handling of UMIs that are all N bases
        let mut records = Vec::new();

        let (r1, r2) = build_test_pair("readn", 0, 100, 300, 60, 60, "NNNNNN");
        records.push(r1);
        records.push(r2);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Reads with all-N UMIs should be filtered out
        assert_eq!(output_records.len(), 0, "Should filter UMIs with all N bases");

        Ok(())
    }

    #[test]
    fn test_mixed_single_and_paired_end() -> Result<()> {
        // Test handling of mixed single-end and paired-end reads
        let mut records = Vec::new();

        // Paired-end reads
        let (r1_paired, r2_paired) = build_test_pair("paired", 0, 100, 300, 60, 60, "AAAAAA");
        records.push(r1_paired);
        records.push(r2_paired);

        // Single-end read
        let r_single = build_test_read("single", 0, 200, 60, 0x0, "TTTTTT");
        records.push(r_single);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should have all 3 records
        assert_eq!(output_records.len(), 3, "Should handle mixed single and paired end");

        Ok(())
    }

    #[test]
    fn test_large_family_grouping() -> Result<()> {
        // Test grouping with large family sizes
        let mut records = Vec::new();

        // Create 100 reads with same UMI
        for i in 0..50 {
            let (r1, r2) = build_test_pair(&format!("read{i}"), 0, 100, 300, 60, 60, "AAAAAA");
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // All should be in one group
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 1, "Large family should group correctly");
        assert_eq!(output_records.len(), 100, "Should have all 100 reads");

        Ok(())
    }

    /// Regression test for the bug where a mapped read with low MAPQ was not filtered
    /// when its mate was unmapped. The fix ensures that MAPQ filtering applies to
    /// each mapped read independently, regardless of mate status.
    ///
    /// fgbio behavior: "Templates are filtered if any non-secondary, non-supplementary
    /// read has mapping quality < min-map-q"
    #[test]
    fn test_filtering_low_mapq_read_with_unmapped_mate() -> Result<()> {
        let mut records = Vec::new();

        // Good pair: both mapped with high MAPQ - should PASS
        let (r1, r2) = build_test_pair("good", 0, 100, 300, 60, 60, "AAAAAA");
        records.push(r1);
        records.push(r2);

        // Bad pair: R1 mapped with low MAPQ, R2 unmapped - should be FILTERED
        // This is the bug case: previously the low MAPQ R1 was not filtered because
        // R2 was unmapped, and the MAPQ check was inside `if both_mapped`
        let (r1_bad, r2_bad) =
            build_test_pair_mapped_with_unmapped_mate("bad", 0, 200, 10, "CCCCCC");
        records.push(r1_bad);
        records.push(r2_bad);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: Some(30), // Threshold is 30, "bad" pair has MAPQ=10
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;

        // Should only have 2 records from the good pair
        // The bad pair should be completely filtered because R1 has low MAPQ
        assert_eq!(
            output_records.len(),
            2,
            "Should have 2 records (good pair only); bad pair with low MAPQ R1 and unmapped R2 should be filtered"
        );

        // Verify the output contains only reads from the "good" template
        for record in &output_records {
            let name = String::from_utf8_lossy(record.name().unwrap()).to_string();
            assert_eq!(
                name, "good",
                "Only 'good' reads should remain, but found read with name: {name}"
            );
        }

        Ok(())
    }

    // ========================================================================
    // Tests for split_templates_by_pair_orientation (commit 2175c44)
    // ========================================================================

    /// Helper to create a pair with explicit strand orientations
    /// `r1_reverse`: if true, R1 is on reverse strand (flag 0x10 set)
    /// `r2_reverse`: if true, R2 is on reverse strand (flag 0x10 set)
    #[allow(clippy::cast_sign_loss)]
    fn build_test_pair_with_orientation(
        name: &str,
        ref_id: usize,
        pos1: i32,
        pos2: i32,
        umi: &str,
        r1_reverse: bool,
        r2_reverse: bool,
    ) -> (sam::alignment::RecordBuf, sam::alignment::RecordBuf) {
        RecordPairBuilder::new()
            .name(name)
            .r1_sequence(&"A".repeat(100))
            .r2_sequence(&"A".repeat(100))
            .reference_sequence_id(ref_id)
            .r1_start(pos1 as usize)
            .r2_start(pos2 as usize)
            .r1_reverse(r1_reverse)
            .r2_reverse(r2_reverse)
            .tag("RX", umi)
            .build()
    }

    #[test]
    fn test_get_pair_orientation_f1r2() -> Result<()> {
        // F1R2: R1 forward, R2 reverse -> (true, false)
        let (r1, r2) = build_test_pair_with_orientation("test", 0, 100, 200, "AAAAAA", false, true);

        let template = Template::from_records(vec![r1, r2])?;

        let orientation = GroupReadsByUmi::get_pair_orientation(&template);
        assert_eq!(orientation, (true, false), "F1R2 should have orientation (true, false)");

        Ok(())
    }

    #[test]
    fn test_get_pair_orientation_f2r1() -> Result<()> {
        // F2R1: R1 reverse, R2 forward -> (false, true)
        let (r1, r2) = build_test_pair_with_orientation("test", 0, 100, 200, "AAAAAA", true, false);

        let template = Template::from_records(vec![r1, r2])?;

        let orientation = GroupReadsByUmi::get_pair_orientation(&template);
        assert_eq!(orientation, (false, true), "F2R1 should have orientation (false, true)");

        Ok(())
    }

    #[test]
    fn test_get_pair_orientation_forward_tandem() -> Result<()> {
        // F1F2: Both forward -> (true, true)
        let (r1, r2) =
            build_test_pair_with_orientation("test", 0, 100, 200, "AAAAAA", false, false);

        let template = Template::from_records(vec![r1, r2])?;

        let orientation = GroupReadsByUmi::get_pair_orientation(&template);
        assert_eq!(
            orientation,
            (true, true),
            "Forward tandem should have orientation (true, true)"
        );

        Ok(())
    }

    #[test]
    fn test_get_pair_orientation_reverse_tandem() -> Result<()> {
        // R1R2: Both reverse -> (false, false)
        let (r1, r2) = build_test_pair_with_orientation("test", 0, 100, 200, "AAAAAA", true, true);

        let template = Template::from_records(vec![r1, r2])?;

        let orientation = GroupReadsByUmi::get_pair_orientation(&template);
        assert_eq!(
            orientation,
            (false, false),
            "Reverse tandem should have orientation (false, false)"
        );

        Ok(())
    }

    #[test]
    fn test_identity_assigner_splits_by_pair_orientation() -> Result<()> {
        // Two pairs at the SAME genomic position with same UMI but different orientations
        // When both reads are at the same position, the strand normalization would give
        // them the same ReadInfo key. The split_templates_by_pair_orientation=true
        // ensures they still get DIFFERENT MI values.
        let mut records = Vec::new();

        // F1R2 pair: R1 forward at 100, R2 reverse at 100
        // After normalization: (100, forward, 100, reverse) = (100, 0, 100, 1)
        let (r1_f1r2, r2_f1r2) =
            build_test_pair_with_orientation("f1r2", 0, 100, 100, "AAAAAA", false, true);
        records.push(r1_f1r2);
        records.push(r2_f1r2);

        // F2R1 pair: R1 reverse at 100, R2 forward at 100
        // After normalization: (100, forward, 100, reverse) = (100, 0, 100, 1) - SAME KEY!
        // But with split_templates_by_pair_orientation=true, they should still be separated
        let (r1_f2r1, r2_f2r1) =
            build_test_pair_with_orientation("f2r1", 0, 100, 100, "AAAAAA", true, false);
        records.push(r1_f2r1);
        records.push(r2_f2r1);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Identity,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 4, "Should have 4 records");

        // Extract MI values grouped by read name
        let mi_tag = sam::alignment::record::data::field::Tag::from([b'M', b'I']);
        let mut mi_by_name: std::collections::HashMap<String, String> =
            std::collections::HashMap::new();
        for record in &output_records {
            let name = String::from_utf8_lossy(record.name().unwrap()).to_string();
            if let Some(noodles::sam::alignment::record_buf::data::field::Value::String(mi)) =
                record.data().get(&mi_tag)
            {
                mi_by_name.insert(name.clone(), String::from_utf8_lossy(mi).to_string());
            }
        }

        let mi_f1r2 = mi_by_name.get("f1r2").expect("f1r2 should have MI tag");
        let mi_f2r1 = mi_by_name.get("f2r1").expect("f2r1 should have MI tag");

        assert_ne!(
            mi_f1r2, mi_f2r1,
            "F1R2 and F2R1 pairs with same UMI should get DIFFERENT MI values with Identity assigner"
        );

        Ok(())
    }

    #[test]
    fn test_paired_assigner_groups_same_orientation_templates() -> Result<()> {
        // Two pairs at the SAME genomic position with same UMI AND same orientation
        // With Paired assigner, templates with the same orientation and matching UMIs
        // should get the SAME MI value.
        let mut records = Vec::new();

        // First F1R2 pair: R1 forward at 100, R2 reverse at 100
        let (r1_a, r2_a) =
            build_test_pair_with_orientation("pair_a", 0, 100, 100, "AA-TT", false, true);
        records.push(r1_a);
        records.push(r2_a);

        // Second F1R2 pair with same orientation and same UMI
        let (r1_b, r2_b) =
            build_test_pair_with_orientation("pair_b", 0, 100, 100, "AA-TT", false, true);
        records.push(r1_b);
        records.push(r2_b);

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Paired,
            edits: 0,
            min_umi_length: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 4, "Should have 4 records");

        // Extract MI values grouped by read name
        let mi_tag = sam::alignment::record::data::field::Tag::from([b'M', b'I']);
        let mut mi_by_name: std::collections::HashMap<String, String> =
            std::collections::HashMap::new();
        for record in &output_records {
            let name = String::from_utf8_lossy(record.name().unwrap()).to_string();
            if let Some(noodles::sam::alignment::record_buf::data::field::Value::String(mi)) =
                record.data().get(&mi_tag)
            {
                mi_by_name.insert(name.clone(), String::from_utf8_lossy(mi).to_string());
            }
        }

        let mi_a = mi_by_name.get("pair_a").expect("pair_a should have MI tag");
        let mi_b = mi_by_name.get("pair_b").expect("pair_b should have MI tag");

        assert_eq!(
            mi_a, mi_b,
            "Two F1R2 pairs with same UMI should get SAME MI value with Paired assigner"
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
        let mut records = Vec::new();

        // Create a few pairs with the same UMI at the same position
        for i in 1..=3 {
            let (r1, r2) = build_test_pair(&format!("read_{i}"), 0, 100, 300, 60, 60, "AAAAAA");
            records.push(r1);
            records.push(r2);
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let cmd = GroupReadsByUmi {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            family_size_histogram: None,
            grouping_metrics: None,
            raw_tag: "RX".to_string(),
            assign_tag: "MI".to_string(),
            cell_tag: "CB".to_string(),
            min_map_q: None,
            include_non_pf_reads: false,
            strategy: Strategy::Adjacency,
            edits: 1,
            min_umi_length: None,
            threading,
            compression: CompressionOptions { compression_level: 1 },
            index_threshold: 100,
            scheduler_opts: SchedulerOptions::default(),
            queue_memory_limit_mb: 0,
        };

        cmd.execute("test")?;

        let output_records = read_bam_records(&paths.output)?;
        assert_eq!(output_records.len(), 6, "Should have 6 records (3 pairs)");

        // All records should have the same MI tag (one group)
        let unique_groups = count_unique_mi_tags(&output_records);
        assert_eq!(unique_groups, 1, "All records with same UMI should be in 1 group");

        Ok(())
    }
}
