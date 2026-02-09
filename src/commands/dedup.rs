//! UMI-aware duplicate marking command.
//!
//! This command marks or removes PCR duplicates using UMI information.
//! It operates on template-coordinate sorted BAM files and requires
//! the `pa` tag on secondary/supplementary reads (added by `fgumi zipper`).
//!
//! # Algorithm
//!
//! 1. Group templates by genomic position (like `fgumi group`)
//! 2. Within each position group, cluster by UMI using the specified strategy
//! 3. Score each template (sum of base qualities from primary reads)
//! 4. Select the highest-scoring template as the representative
//! 5. Mark all other templates as duplicates
//!
//! # Output Modes
//!
//! - Mark only: Set duplicate flag on non-representative reads (default)
//! - Remove: Exclude duplicate reads from output

use std::io;
use std::path::PathBuf;
use std::sync::Arc;

use ahash::AHashMap;
use anyhow::{Context, Result, bail};
use bstr::BString;
use clap::Parser;
use crossbeam_queue::SegQueue;
use fgoxide::io::DelimFile;
use fgumi_lib::assigner::{PairedUmiAssigner, Strategy, UmiAssigner};
use fgumi_lib::bam_io::{create_bam_reader_for_pipeline, is_stdin_path};
use fgumi_lib::grouper::{FilterMetrics, PositionGroup, PositionGrouper, PositionGrouperConfig};
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::metrics::group::FamilySizeMetrics;
use fgumi_lib::read_info::build_library_lookup;
use fgumi_lib::sam::{is_template_coordinate_sorted, unclipped_five_prime_position};
use fgumi_lib::template::{MoleculeId, Template};
use fgumi_lib::umi::{UmiValidation, validate_umi};
use fgumi_lib::unified_pipeline::{
    BamPipelineConfig, BatchWeight, Grouper, MemoryEstimate, run_bam_pipeline_from_reader,
    serialize_bam_records_into,
};
use fgumi_lib::validation::{string_to_tag, validate_file_exists, validate_tag};
use log::info;
use noodles::sam::Header;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::value::Value as DataValue;
use serde::{Deserialize, Serialize};

use crate::commands::command::Command;
use crate::commands::common::{
    BamIoOptions, CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
};
use fgumi_lib::sort::PA_TAG;

/// Duplicate flag bit in SAM flags (0x400)
const DUPLICATE_FLAG: u16 = 0x400;

//////////////////////////////////////////////////////////////////////////////
// Metrics
//////////////////////////////////////////////////////////////////////////////

/// Metrics collected during deduplication.
#[derive(Debug, Default, Clone)]
pub struct DedupMetrics {
    /// Total templates processed
    pub total_templates: u64,
    /// Templates marked as duplicates
    pub duplicate_templates: u64,
    /// Templates kept (not duplicates)
    pub unique_templates: u64,
    /// Total reads processed
    pub total_reads: u64,
    /// Reads marked as duplicates
    pub duplicate_reads: u64,
    /// Unique reads (not duplicates)
    pub unique_reads: u64,
    /// Secondary reads processed
    pub secondary_reads: u64,
    /// Supplementary reads processed
    pub supplementary_reads: u64,
    /// Secondary/supplementary without pa tag
    pub missing_pa_tag: u64,
    /// Filter metrics from position grouping
    pub filter_metrics: FilterMetrics,
}

impl DedupMetrics {
    /// Merge another `DedupMetrics` into this one.
    pub fn merge(&mut self, other: &DedupMetrics) {
        self.total_templates += other.total_templates;
        self.duplicate_templates += other.duplicate_templates;
        self.unique_templates += other.unique_templates;
        self.total_reads += other.total_reads;
        self.duplicate_reads += other.duplicate_reads;
        self.unique_reads += other.unique_reads;
        self.secondary_reads += other.secondary_reads;
        self.supplementary_reads += other.supplementary_reads;
        self.missing_pa_tag += other.missing_pa_tag;
        self.filter_metrics.merge(&other.filter_metrics);
    }

    /// Calculate duplicate rate.
    #[must_use]
    pub fn duplicate_rate(&self) -> f64 {
        if self.total_templates == 0 {
            0.0
        } else {
            self.duplicate_templates as f64 / self.total_templates as f64
        }
    }
}

/// Serializable version of `DedupMetrics` for file output.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct DedupMetricsOutput {
    total_templates: u64,
    unique_templates: u64,
    duplicate_templates: u64,
    duplicate_rate: f64,
    total_reads: u64,
    unique_reads: u64,
    duplicate_reads: u64,
    secondary_reads: u64,
    supplementary_reads: u64,
    missing_pa_tag: u64,
}

impl From<&DedupMetrics> for DedupMetricsOutput {
    fn from(m: &DedupMetrics) -> Self {
        Self {
            total_templates: m.total_templates,
            unique_templates: m.unique_templates,
            duplicate_templates: m.duplicate_templates,
            duplicate_rate: m.duplicate_rate(),
            total_reads: m.total_reads,
            unique_reads: m.unique_reads,
            duplicate_reads: m.duplicate_reads,
            secondary_reads: m.secondary_reads,
            supplementary_reads: m.supplementary_reads,
            missing_pa_tag: m.missing_pa_tag,
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
// Collected metrics for pipeline aggregation
//////////////////////////////////////////////////////////////////////////////

/// Metrics collected per position group, aggregated after pipeline completion.
#[derive(Default, Debug)]
struct CollectedDedupMetrics {
    /// Dedup-specific metrics
    dedup_metrics: DedupMetrics,
    /// Family size counts
    family_sizes: AHashMap<usize, u64>,
}

//////////////////////////////////////////////////////////////////////////////
// Processed position group for dedup
//////////////////////////////////////////////////////////////////////////////

/// Result of processing a position group for deduplication.
pub struct ProcessedDedupGroup {
    /// Templates with duplicate flags set appropriately.
    pub templates: Vec<Template>,
    /// Family size counts for this group.
    pub family_sizes: AHashMap<usize, u64>,
    /// Dedup metrics for this group.
    pub dedup_metrics: DedupMetrics,
    /// Total input records processed (for progress tracking).
    pub input_record_count: u64,
    /// Base MI offset for global uniqueness.
    pub base_mi: u64,
}

impl BatchWeight for ProcessedDedupGroup {
    fn batch_weight(&self) -> usize {
        self.templates.len()
    }
}

impl MemoryEstimate for ProcessedDedupGroup {
    fn estimate_heap_size(&self) -> usize {
        self.templates.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
            + self.family_sizes.capacity() * std::mem::size_of::<(usize, u64)>()
            + std::mem::size_of::<DedupMetrics>()
    }
}

//////////////////////////////////////////////////////////////////////////////
// Configuration
//////////////////////////////////////////////////////////////////////////////

/// Configuration for template filtering during deduplication.
#[derive(Clone)]
struct DedupFilterConfig {
    /// UMI tag bytes (e.g., [b'R', b'X']).
    umi_tag: [u8; 2],
    /// Minimum mapping quality.
    min_mapq: u8,
    /// Whether to include non-PF reads.
    include_non_pf: bool,
    /// Minimum UMI length.
    min_umi_length: Option<usize>,
}

//////////////////////////////////////////////////////////////////////////////
// Duplicate scoring
//////////////////////////////////////////////////////////////////////////////

/// Scores a template for duplicate selection.
/// Higher score = more likely to be chosen as representative.
///
/// Scoring is based on sum of base qualities from primary reads,
/// matching HTSJDK/Picard's `SUM_OF_BASE_QUALITIES` strategy.
#[inline]
fn score_template(template: &Template) -> i64 {
    let mut score: u32 = 0;

    // Score R1 primary
    if let Some(r1) = template.r1() {
        score += sum_base_qualities(r1);
    }

    // Score R2 primary
    if let Some(r2) = template.r2() {
        score += sum_base_qualities(r2);
    }

    i64::from(score)
}

/// Sums base qualities for a record, capping at 15 per base (like Picard).
///
/// Uses u32 accumulator for faster arithmetic. u32 is sufficient for reads
/// up to ~286 million bases (2^32 / 15), far exceeding any real read length.
#[inline]
fn sum_base_qualities(record: &noodles::sam::alignment::RecordBuf) -> u32 {
    record.quality_scores().as_ref().iter().map(|&q| u32::from(std::cmp::min(q, 15))).sum()
}

//////////////////////////////////////////////////////////////////////////////
// Template filtering (adapted from group command)
//////////////////////////////////////////////////////////////////////////////

/// Filter a template based on filtering criteria.
/// Returns true if the template should be kept.
fn filter_template(
    template: &Template,
    config: &DedupFilterConfig,
    metrics: &mut FilterMetrics,
) -> bool {
    let r1 = template.r1();
    let r2 = template.r2();

    metrics.total_templates += 1;

    // Need at least one primary read
    if r1.is_none() && r2.is_none() {
        metrics.discarded_poor_alignment += 1;
        return false;
    }

    // Check if both reads are unmapped
    let both_unmapped =
        r1.is_none_or(|r| r.flags().is_unmapped()) && r2.is_none_or(|r| r.flags().is_unmapped());

    if both_unmapped {
        metrics.discarded_poor_alignment += 1;
        return false;
    }

    // Phase 1: Flag-based checks
    for record in [r1, r2].into_iter().flatten() {
        let flags = record.flags();

        // Filter non-PF if requested
        if !config.include_non_pf && flags.is_qc_fail() {
            metrics.discarded_non_pf += 1;
            return false;
        }

        // Check MAPQ for mapped reads
        if !flags.is_unmapped() {
            let mapq = record.mapping_quality().map_or(0, u8::from);
            if mapq < config.min_mapq {
                metrics.discarded_poor_alignment += 1;
                return false;
            }
        }
    }

    // Phase 2: Tag-based checks
    for record in [r1, r2].into_iter().flatten() {
        let flags = record.flags();

        // Check mate MAPQ via MQ tag
        if !flags.is_mate_unmapped() {
            if let Some(data) = record.data().get(b"MQ") {
                if let Some(mq) = data.as_int() {
                    #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
                    if (mq as u8) < config.min_mapq {
                        metrics.discarded_poor_alignment += 1;
                        return false;
                    }
                }
            }
        }

        // Check UMI for Ns and minimum length using common validation
        if let Some(DataValue::String(umi)) = record.data().get(&config.umi_tag) {
            match validate_umi(umi) {
                UmiValidation::ContainsN => {
                    metrics.discarded_ns_in_umi += 1;
                    return false;
                }
                UmiValidation::Valid(base_count) => {
                    if let Some(min_len) = config.min_umi_length {
                        if base_count < min_len {
                            metrics.discarded_umi_too_short += 1;
                            return false;
                        }
                    }
                }
            }
        } else {
            metrics.discarded_poor_alignment += 1;
            return false;
        }
    }

    metrics.accepted_templates += 1;
    true
}

//////////////////////////////////////////////////////////////////////////////
// UMI assignment (adapted from group command)
//////////////////////////////////////////////////////////////////////////////

/// Extract UMI for a read, handling paired UMI strategies.
fn umi_for_read(umi: &str, is_r1_earlier: bool, assigner: &dyn UmiAssigner) -> Result<String> {
    if assigner.split_templates_by_pair_orientation() {
        if umi.bytes().all(|b| !b.is_ascii_lowercase()) {
            Ok(umi.to_owned())
        } else {
            Ok(umi.to_uppercase())
        }
    } else {
        let parts: Vec<&str> = umi.split('-').collect();
        if parts.len() != 2 {
            bail!("Paired strategy used but UMI did not contain 2 segments: {umi}");
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

/// Get pair orientation for a template.
fn get_pair_orientation(template: &Template) -> (bool, bool) {
    let r1_positive = template.r1().is_none_or(|r| !r.flags().is_reverse_complemented());
    let r2_positive = template.r2().is_none_or(|r| !r.flags().is_reverse_complemented());
    (r1_positive, r2_positive)
}

/// Check if R1 is genomically earlier than R2.
fn is_r1_genomically_earlier(
    r1: &noodles::sam::alignment::RecordBuf,
    r2: &noodles::sam::alignment::RecordBuf,
) -> Result<bool> {
    let r1_pos = unclipped_five_prime_position(r1).unwrap_or(0);
    let r2_pos = unclipped_five_prime_position(r2).unwrap_or(0);
    Ok(r1_pos <= r2_pos)
}

/// Truncate UMIs to minimum length if specified.
fn truncate_umis(umis: Vec<String>, min_umi_length: Option<usize>) -> Result<Vec<String>> {
    match min_umi_length {
        None => Ok(umis),
        Some(min_len) => {
            let min_length = umis.iter().map(String::len).min().unwrap_or(0);
            if min_length < min_len {
                bail!("UMI found shorter than expected ({min_length} < {min_len})");
            }
            Ok(umis.into_iter().map(|u| u[..min_len].to_string()).collect())
        }
    }
}

/// Assign UMI groups to a subset of templates.
fn assign_umi_groups_for_indices(
    templates: &mut [Template],
    indices: &[usize],
    assigner: &dyn UmiAssigner,
    raw_tag: [u8; 2],
    min_umi_length: Option<usize>,
) -> Result<()> {
    if indices.is_empty() {
        return Ok(());
    }

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

        let umi_str =
            std::str::from_utf8(umi_bytes).map_err(|e| anyhow::anyhow!("Invalid UTF-8: {e}"))?;

        let is_r1_earlier = if let (Some(r1), Some(r2)) = (template.r1(), template.r2()) {
            is_r1_genomically_earlier(r1, r2)?
        } else {
            true
        };

        let processed_umi = umi_for_read(umi_str, is_r1_earlier, assigner)?;
        umis.push(processed_umi);
    }

    let truncated_umis = truncate_umis(umis, min_umi_length)?;
    let assignments = assigner.assign(&truncated_umis);

    for (i, &idx) in indices.iter().enumerate() {
        templates[idx].mi = assignments[i];
    }

    Ok(())
}

/// Assign UMI groups to templates.
fn assign_umi_groups(
    templates: &mut [Template],
    assigner: &dyn UmiAssigner,
    raw_tag: [u8; 2],
    min_umi_length: Option<usize>,
) -> Result<()> {
    if assigner.split_templates_by_pair_orientation() {
        let mut subgroups: AHashMap<(bool, bool), Vec<usize>> = AHashMap::new();
        for (idx, template) in templates.iter().enumerate() {
            let orientation = get_pair_orientation(template);
            subgroups.entry(orientation).or_default().push(idx);
        }

        for indices in subgroups.values() {
            assign_umi_groups_for_indices(templates, indices, assigner, raw_tag, min_umi_length)?;
        }
    } else {
        let all_indices: Vec<usize> = (0..templates.len()).collect();
        assign_umi_groups_for_indices(templates, &all_indices, assigner, raw_tag, min_umi_length)?;
    }

    Ok(())
}

//////////////////////////////////////////////////////////////////////////////
// Duplicate marking
//////////////////////////////////////////////////////////////////////////////

/// Marks duplicates within a UMI family.
/// The highest-scoring template is kept; others are marked as duplicates.
///
/// Optimized with fast paths for small families (size 1-3) to avoid
/// allocation and sorting overhead for the common case.
fn mark_duplicates_in_family(templates: &mut [&mut Template], dedup_metrics: &mut DedupMetrics) {
    match templates.len() {
        0 => return,
        1 => {
            // Single template - not a duplicate
            dedup_metrics.unique_templates += 1;
            return;
        }
        2 => {
            // Fast path for size 2 (very common case) - no allocation needed
            let s0 = score_template(templates[0]);
            let s1 = score_template(templates[1]);
            dedup_metrics.unique_templates += 1;
            dedup_metrics.duplicate_templates += 1;
            if s0 >= s1 {
                mark_template_as_duplicate(templates[1], dedup_metrics);
            } else {
                mark_template_as_duplicate(templates[0], dedup_metrics);
            }
            return;
        }
        3 => {
            // Fast path for size 3 - find max without allocation
            let s0 = score_template(templates[0]);
            let s1 = score_template(templates[1]);
            let s2 = score_template(templates[2]);
            dedup_metrics.unique_templates += 1;
            dedup_metrics.duplicate_templates += 2;
            if s0 >= s1 && s0 >= s2 {
                mark_template_as_duplicate(templates[1], dedup_metrics);
                mark_template_as_duplicate(templates[2], dedup_metrics);
            } else if s1 >= s0 && s1 >= s2 {
                mark_template_as_duplicate(templates[0], dedup_metrics);
                mark_template_as_duplicate(templates[2], dedup_metrics);
            } else {
                mark_template_as_duplicate(templates[0], dedup_metrics);
                mark_template_as_duplicate(templates[1], dedup_metrics);
            }
            return;
        }
        _ => {} // Fall through to general case
    }

    // General case for families of size 4+
    // Score all templates
    let mut scores: Vec<(usize, i64)> =
        templates.iter().enumerate().map(|(i, t)| (i, score_template(t))).collect();

    // Sort by score descending (highest first)
    scores.sort_by(|a, b| b.1.cmp(&a.1));

    // First template (highest score) is representative
    dedup_metrics.unique_templates += 1;

    // Mark all others as duplicates
    for (idx, _score) in scores.iter().skip(1) {
        mark_template_as_duplicate(templates[*idx], dedup_metrics);
        dedup_metrics.duplicate_templates += 1;
    }
}

/// Marks all reads in a template as duplicates.
fn mark_template_as_duplicate(template: &mut Template, dedup_metrics: &mut DedupMetrics) {
    for record in &mut template.records {
        let current_flags = u16::from(record.flags());
        let new_flags = Flags::from(current_flags | DUPLICATE_FLAG);
        *record.flags_mut() = new_flags;
        dedup_metrics.duplicate_reads += 1;
    }
}

//////////////////////////////////////////////////////////////////////////////
// Position group processing
//////////////////////////////////////////////////////////////////////////////

/// Process a position group for deduplication.
fn process_position_group(
    group: PositionGroup,
    filter_config: &DedupFilterConfig,
    assigner: &dyn UmiAssigner,
    raw_tag: [u8; 2],
    min_umi_length: Option<usize>,
) -> io::Result<ProcessedDedupGroup> {
    let mut dedup_metrics = DedupMetrics::default();
    let input_record_count: u64 = group.templates.iter().map(|t| t.read_count() as u64).sum();
    let base_mi = group.base_mi;

    // Filter templates
    let mut filter_metrics = FilterMetrics::new();
    let filtered_templates: Vec<Template> = group
        .templates
        .into_iter()
        .filter(|t| filter_template(t, filter_config, &mut filter_metrics))
        .collect();

    dedup_metrics.filter_metrics = filter_metrics;

    if filtered_templates.is_empty() {
        return Ok(ProcessedDedupGroup {
            templates: Vec::new(),
            family_sizes: AHashMap::new(),
            dedup_metrics,
            input_record_count,
            base_mi,
        });
    }

    // Clear existing duplicate flags before re-marking
    let mut templates: Vec<Template> = filtered_templates
        .into_iter()
        .map(|mut t| {
            for record in &mut t.records {
                let current_flags = u16::from(record.flags());
                let new_flags = Flags::from(current_flags & !DUPLICATE_FLAG);
                *record.flags_mut() = new_flags;
            }
            t
        })
        .collect();
    if let Err(e) = assign_umi_groups(&mut templates, assigner, raw_tag, min_umi_length) {
        return Err(io::Error::new(io::ErrorKind::InvalidData, e));
    }

    // Sort by molecule ID for grouping
    templates.sort_by(|a, b| {
        let a_idx = a.mi.to_vec_index();
        let b_idx = b.mi.to_vec_index();
        a_idx.cmp(&b_idx).then_with(|| a.name.cmp(&b.name))
    });

    // Group by MI and mark duplicates
    let mut family_sizes: AHashMap<usize, u64> = AHashMap::with_capacity(50);

    if !templates.is_empty() {
        // Collect family boundaries
        let mut family_boundaries: Vec<usize> = vec![0];
        let mut current_mi = templates[0].mi.to_vec_index();

        for (i, template) in templates.iter().enumerate().skip(1) {
            let mi = template.mi.to_vec_index();
            if mi != current_mi {
                family_boundaries.push(i);
                current_mi = mi;
            }
        }
        family_boundaries.push(templates.len());

        // Process each family
        for window in family_boundaries.windows(2) {
            let start = window[0];
            let end = window[1];
            let family_size = end - start;

            // Skip unassigned templates
            if templates[start].mi.to_vec_index().is_none() {
                continue;
            }

            *family_sizes.entry(family_size).or_insert(0) += 1;

            // Get mutable references to family templates (single collection, not double)
            let mut family_refs: Vec<&mut Template> = templates[start..end].iter_mut().collect();
            mark_duplicates_in_family(&mut family_refs, &mut dedup_metrics);
        }
    }

    // Count reads and check for missing pa tags
    for template in &templates {
        dedup_metrics.total_templates += 1;
        for record in &template.records {
            dedup_metrics.total_reads += 1;
            let flags = record.flags();
            let is_secondary = flags.is_secondary();
            let is_supplementary = flags.is_supplementary();

            if is_secondary {
                dedup_metrics.secondary_reads += 1;
            }
            if is_supplementary {
                dedup_metrics.supplementary_reads += 1;
            }
            // Single PA_TAG lookup for both secondary and supplementary
            if (is_secondary || is_supplementary) && record.data().get(&PA_TAG).is_none() {
                dedup_metrics.missing_pa_tag += 1;
            }
        }
    }

    dedup_metrics.unique_reads = dedup_metrics.total_reads - dedup_metrics.duplicate_reads;

    Ok(ProcessedDedupGroup { templates, family_sizes, dedup_metrics, input_record_count, base_mi })
}

//////////////////////////////////////////////////////////////////////////////
// Set MI tag on record
//////////////////////////////////////////////////////////////////////////////

/// Set MI tag on a single record.
#[inline]
fn set_mi_tag_on_record(
    record: &mut noodles::sam::alignment::RecordBuf,
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

//////////////////////////////////////////////////////////////////////////////
// Command definition
//////////////////////////////////////////////////////////////////////////////

/// UMI-aware duplicate marking command.
#[derive(Debug, Parser)]
#[command(
    name = "dedup",
    about = "\x1b[38;5;151m[DEDUP]\x1b[0m         \x1b[36mMark or remove PCR duplicates using UMI information\x1b[0m",
    long_about = r#"
Marks or removes PCR duplicates from a BAM file using UMI information.
Requires template-coordinate sorted input with `pa` tags on secondary/supplementary
reads (added by `fgumi zipper`).

Within each UMI family, the template with the highest sum of base qualities
is selected as the representative; all others are marked as duplicates.

# Input Requirements

- Must be processed with `fgumi zipper` (adds `pa` tag for secondary/supplementary reads)
- Must be sorted with `fgumi sort --order template-coordinate`
- UMI tags on reads (default: RX tag)

Note: Using `samtools sort` will NOT work correctly because it doesn't use the
`pa` tag for template-coordinate ordering of secondary/supplementary reads.

# Output Modes

- Mark only (default): Set duplicate flag (0x400) on non-representative reads
- Remove (--remove-duplicates): Exclude duplicate reads from output entirely
"#
)]
pub struct MarkDuplicates {
    /// Input and output BAM files
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Path to write deduplication metrics
    #[arg(short = 'm', long = "metrics")]
    pub metrics: Option<PathBuf>,

    /// Path to write family size histogram
    #[arg(short = 'H', long = "family-size-histogram")]
    pub family_size_histogram: Option<PathBuf>,

    /// Remove duplicates instead of just marking them
    #[arg(short = 'r', long = "remove-duplicates", default_value = "false")]
    pub remove_duplicates: bool,

    /// The tag containing the raw UMI sequence
    #[arg(short = 't', long = "raw-tag", default_value = "RX")]
    pub raw_tag: String,

    /// The output tag for the assigned molecule ID
    #[arg(short = 'T', long = "assign-tag", default_value = "MI")]
    pub assign_tag: String,

    /// The tag containing the cell barcode (for single-cell data)
    #[arg(short = 'c', long = "cell-tag", default_value = "CB")]
    pub cell_tag: String,

    /// Minimum mapping quality for a read to be included
    #[arg(short = 'q', long = "min-map-q")]
    pub min_map_q: Option<u8>,

    /// Include reads flagged as not passing QC
    #[arg(short = 'n', long = "include-non-pf-reads", default_value = "false")]
    pub include_non_pf_reads: bool,

    /// UMI grouping strategy
    #[arg(short = 's', long = "strategy", value_enum, default_value = "adjacency")]
    pub strategy: Strategy,

    /// Maximum edit distance for UMI grouping
    #[arg(short = 'e', long = "edits", default_value = "1")]
    pub edits: u32,

    /// Minimum UMI length (UMIs shorter than this are discarded)
    #[arg(short = 'l', long = "min-umi-length")]
    pub min_umi_length: Option<usize>,

    /// Threading options for parallel processing
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output BAM
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Minimum UMIs per position to use index for faster grouping
    #[arg(long = "index-threshold", default_value = "100")]
    pub index_threshold: usize,

    /// Scheduler and pipeline options
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

impl Command for MarkDuplicates {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate strategy/min-umi-length combination
        if self.min_umi_length.is_some() && matches!(self.strategy, Strategy::Paired) {
            bail!("Paired strategy cannot be used with --min-umi-length");
        }

        // Validate input file exists
        if !is_stdin_path(&self.io.input) {
            validate_file_exists(&self.io.input, "input BAM file")?;
        }

        let min_mapq: u8 = self.min_map_q.unwrap_or(0);

        let timer = OperationTimer::new("Marking duplicates");

        info!("Starting dedup");
        info!("Input: {}", self.io.input.display());
        info!("Output: {}", self.io.output.display());
        info!("Strategy: {:?}", self.strategy);
        info!("Edits: {}", self.edits);
        info!("Remove duplicates: {}", self.remove_duplicates);
        if matches!(self.strategy, Strategy::Adjacency | Strategy::Paired) {
            info!("Index threshold: {}", self.index_threshold);
        }
        info!("{}", self.threading.log_message());

        // Open input BAM
        let (reader, header) = create_bam_reader_for_pipeline(&self.io.input)?;

        if !is_template_coordinate_sorted(&header) {
            bail!(
                "Input BAM must be template-coordinate sorted.\n\n\
                To sort your BAM file, run:\n  \
                samtools sort --template-coordinate input.bam -o sorted.bam"
            );
        }
        info!("Template-coordinate sorted");

        // Add @PG record
        let header = fgumi_lib::header::add_pg_record(
            header,
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        // Prepare tags
        if self.raw_tag.len() != 2 {
            bail!("Raw tag must be exactly 2 characters");
        }
        let raw_tag: [u8; 2] = {
            let bytes = self.raw_tag.as_bytes();
            [bytes[0], bytes[1]]
        };

        let cell_tag = string_to_tag(&self.cell_tag, "cell-tag")?;

        let assign_tag_bytes: [u8; 2] = validate_tag(&self.assign_tag, "assign-tag")?;

        let effective_edits =
            if matches!(self.strategy, Strategy::Identity) { 0 } else { self.edits };

        let filter_config = DedupFilterConfig {
            umi_tag: raw_tag,
            min_mapq,
            include_non_pf: self.include_non_pf_reads,
            min_umi_length: self.min_umi_length,
        };

        // Shared state for collecting metrics
        let collected_metrics: Arc<SegQueue<CollectedDedupMetrics>> = Arc::new(SegQueue::new());

        // Clone values needed by closures
        let strategy = self.strategy;
        let index_threshold = self.index_threshold;
        let min_umi_length = self.min_umi_length;
        let remove_duplicates = self.remove_duplicates;
        let collected_metrics_clone = Arc::clone(&collected_metrics);

        // Configure pipeline
        let num_threads = self.threading.num_threads();
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

        // Calculate and apply queue memory limit
        let queue_memory_limit_bytes = self.queue_memory.calculate_memory_limit(num_threads)?;
        pipeline_config.pipeline.queue_memory_limit = queue_memory_limit_bytes;
        self.queue_memory.log_memory_config(num_threads, queue_memory_limit_bytes);
        info!("Scheduler: {:?}", self.scheduler_opts.strategy());
        info!("Using pipeline with {num_threads} threads");

        // Run the pipeline
        let _records_processed = run_bam_pipeline_from_reader(
            pipeline_config,
            reader,
            header,
            &self.io.output,
            None,
            // Grouper factory - use with_secondary_supplementary to include all reads
            move |header: &Header| {
                let library_lookup = build_library_lookup(header);
                let config =
                    PositionGrouperConfig::with_secondary_supplementary(cell_tag, library_lookup);
                Box::new(PositionGrouper::new(config))
                    as Box<dyn Grouper<Group = PositionGroup> + Send>
            },
            // Process function (parallel)
            move |group: PositionGroup| -> io::Result<ProcessedDedupGroup> {
                let assigner = strategy.new_assigner_full(effective_edits, 1, index_threshold);
                process_position_group(
                    group,
                    &filter_config,
                    assigner.as_ref(),
                    raw_tag,
                    min_umi_length,
                )
            },
            // Serialize function (parallel)
            move |processed: ProcessedDedupGroup,
                  header: &Header,
                  output: &mut Vec<u8>|
                  -> io::Result<u64> {
                // Collect metrics
                let metrics = CollectedDedupMetrics {
                    dedup_metrics: processed.dedup_metrics.clone(),
                    family_sizes: processed.family_sizes.clone(),
                };
                collected_metrics_clone.push(metrics);

                let input_record_count = processed.input_record_count;
                let base_mi = processed.base_mi;
                let assign_tag = Tag::from(assign_tag_bytes);

                // Serialize templates
                for template in processed.templates {
                    let mi = template.mi;
                    for mut record in template.records {
                        // Skip duplicates if remove mode is enabled
                        if remove_duplicates && (u16::from(record.flags()) & DUPLICATE_FLAG) != 0 {
                            continue;
                        }

                        // Set MI tag
                        set_mi_tag_on_record(&mut record, mi, assign_tag, base_mi);

                        // Serialize to output buffer
                        serialize_bam_records_into(&[record], header, output)?;
                    }
                }

                Ok(input_record_count)
            },
        )?;

        // Aggregate metrics
        let mut final_metrics = DedupMetrics::default();
        let mut final_family_sizes: AHashMap<usize, u64> = AHashMap::new();

        while let Some(m) = collected_metrics.pop() {
            final_metrics.merge(&m.dedup_metrics);
            for (size, count) in m.family_sizes {
                *final_family_sizes.entry(size).or_insert(0) += count;
            }
        }

        // Write metrics file
        if let Some(metrics_path) = &self.metrics {
            write_dedup_metrics(&final_metrics, metrics_path)?;
        }

        // Write family size histogram
        if let Some(histogram_path) = &self.family_size_histogram {
            write_family_size_histogram(&final_family_sizes, histogram_path)?;
        }

        // Log summary
        info!(
            "Deduplication complete: {} templates ({} unique, {} duplicates, {:.2}% duplicate rate)",
            final_metrics.total_templates,
            final_metrics.unique_templates,
            final_metrics.duplicate_templates,
            final_metrics.duplicate_rate() * 100.0
        );

        if final_metrics.missing_pa_tag > 0 {
            bail!(
                "{} secondary/supplementary reads are missing the `pa` tag.\n\n\
                The `pa` tag is required for correct UMI-aware deduplication of \
                secondary and supplementary alignments. This tag is added by \
                `fgumi zipper` during the merge of unmapped and mapped BAMs.\n\n\
                To fix this, re-run your pipeline starting from `fgumi zipper`:\n  \
                fgumi zipper --unmapped unmapped.bam --mapped aligned.bam -o merged.bam\n  \
                fgumi sort -i merged.bam -o sorted.bam --order template-coordinate\n  \
                fgumi dedup -i sorted.bam -o deduped.bam",
                final_metrics.missing_pa_tag
            );
        }

        timer.log_completion(final_metrics.total_reads);

        Ok(())
    }
}

//////////////////////////////////////////////////////////////////////////////
// Metrics writing
//////////////////////////////////////////////////////////////////////////////

fn write_dedup_metrics(metrics: &DedupMetrics, path: &PathBuf) -> Result<()> {
    let output: DedupMetricsOutput = metrics.into();
    DelimFile::default()
        .write_tsv(path, [output])
        .with_context(|| format!("Failed to write dedup metrics: {}", path.display()))?;
    Ok(())
}

fn write_family_size_histogram(family_sizes: &AHashMap<usize, u64>, path: &PathBuf) -> Result<()> {
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
        metrics.push(FamilySizeMetrics {
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
        .with_context(|| format!("Failed to write family size histogram: {}", path.display()))?;
    Ok(())
}

//////////////////////////////////////////////////////////////////////////////
// Tests
//////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_lib::sam::builder::RecordBuilder;

    // Helper to create a simple template for testing
    fn create_test_template(name: &str, qualities: &[u8]) -> Template {
        let record = RecordBuilder::new()
            .name(name)
            .sequence("ACGT")
            .qualities(qualities)
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .build();

        Template::from_records(vec![record]).unwrap()
    }

    #[test]
    fn test_score_template() {
        // Each quality capped at 15
        let template = create_test_template("q1", &[20, 20, 20, 20]);
        let score = score_template(&template);
        // 4 bases * 15 (capped) = 60
        assert_eq!(score, 60);
    }

    #[test]
    fn test_score_template_low_quality() {
        let template = create_test_template("q1", &[10, 10, 10, 10]);
        let score = score_template(&template);
        // 4 bases * 10 = 40
        assert_eq!(score, 40);
    }

    /// Creates a paired-end test template with R1 and R2.
    fn create_paired_test_template(
        name: &str,
        r1_qualities: &[u8],
        r2_qualities: &[u8],
    ) -> Template {
        let r1 = RecordBuilder::new()
            .name(name)
            .sequence("ACGT")
            .qualities(r1_qualities)
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .build();

        let r2 = RecordBuilder::new()
            .name(name)
            .sequence("TGCA")
            .qualities(r2_qualities)
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .build();

        Template::from_records(vec![r1, r2]).unwrap()
    }

    #[test]
    fn test_score_template_paired_end() {
        // Both R1 and R2 contribute to score
        let template = create_paired_test_template("q1", &[10, 10, 10, 10], &[10, 10, 10, 10]);
        let score = score_template(&template);
        // R1: 4 * 10 = 40, R2: 4 * 10 = 40, total = 80
        assert_eq!(score, 80);
    }

    #[test]
    fn test_score_template_paired_end_capped() {
        // High quality gets capped at 15 per base
        let template = create_paired_test_template("q1", &[40, 40, 40, 40], &[40, 40, 40, 40]);
        let score = score_template(&template);
        // R1: 4 * 15 = 60, R2: 4 * 15 = 60, total = 120
        assert_eq!(score, 120);
    }

    #[test]
    fn test_score_template_paired_end_asymmetric() {
        // R1 and R2 have different qualities
        let template = create_paired_test_template("q1", &[5, 5, 5, 5], &[15, 15, 15, 15]);
        let score = score_template(&template);
        // R1: 4 * 5 = 20, R2: 4 * 15 = 60, total = 80
        assert_eq!(score, 80);
    }

    #[test]
    fn test_score_template_paired_end_mixed_capping() {
        // Some qualities below cap, some above
        let template = create_paired_test_template("q1", &[10, 20, 5, 30], &[8, 12, 25, 3]);
        let score = score_template(&template);
        // R1: 10 + 15 + 5 + 15 = 45 (20 and 30 capped to 15)
        // R2: 8 + 12 + 15 + 3 = 38 (25 capped to 15)
        // Total: 45 + 38 = 83
        assert_eq!(score, 83);
    }

    #[test]
    fn test_duplicate_rate_calculation() {
        let metrics =
            DedupMetrics { total_templates: 100, duplicate_templates: 25, ..Default::default() };
        assert!((metrics.duplicate_rate() - 0.25).abs() < 0.001);
    }

    #[test]
    fn test_duplicate_rate_zero_templates() {
        let metrics = DedupMetrics::default();
        assert!((metrics.duplicate_rate() - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_metrics_merge() {
        let mut m1 =
            DedupMetrics { total_templates: 10, duplicate_templates: 2, ..Default::default() };

        let m2 = DedupMetrics { total_templates: 20, duplicate_templates: 5, ..Default::default() };

        m1.merge(&m2);
        assert_eq!(m1.total_templates, 30);
        assert_eq!(m1.duplicate_templates, 7);
    }

    // ========================================================================
    // mark_duplicates_in_family tests (fast path coverage)
    // ========================================================================

    /// Helper to check if a template is marked as duplicate
    fn is_duplicate(template: &Template) -> bool {
        template.records.iter().any(|r| (u16::from(r.flags()) & DUPLICATE_FLAG) != 0)
    }

    #[test]
    fn test_mark_duplicates_empty_family() {
        let mut metrics = DedupMetrics::default();
        let mut templates: Vec<&mut Template> = vec![];
        mark_duplicates_in_family(&mut templates, &mut metrics);
        assert_eq!(metrics.unique_templates, 0);
        assert_eq!(metrics.duplicate_templates, 0);
    }

    #[test]
    fn test_mark_duplicates_single_template() {
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[30, 30, 30, 30]);
        let mut templates: Vec<&mut Template> = vec![&mut t1];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.unique_templates, 1);
        assert_eq!(metrics.duplicate_templates, 0);
        assert!(!is_duplicate(&t1));
    }

    #[test]
    fn test_mark_duplicates_size_2_first_higher() {
        // Fast path for size 2: first template has higher score
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[15, 15, 15, 15]); // Score: 60
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 40
        let mut templates: Vec<&mut Template> = vec![&mut t1, &mut t2];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.unique_templates, 1);
        assert_eq!(metrics.duplicate_templates, 1);
        assert!(!is_duplicate(&t1)); // Higher score, not duplicate
        assert!(is_duplicate(&t2)); // Lower score, is duplicate
    }

    #[test]
    fn test_mark_duplicates_size_2_second_higher() {
        // Fast path for size 2: second template has higher score
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[5, 5, 5, 5]); // Score: 20
        let mut t2 = create_test_template("q2", &[15, 15, 15, 15]); // Score: 60
        let mut templates: Vec<&mut Template> = vec![&mut t1, &mut t2];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.unique_templates, 1);
        assert_eq!(metrics.duplicate_templates, 1);
        assert!(is_duplicate(&t1)); // Lower score, is duplicate
        assert!(!is_duplicate(&t2)); // Higher score, not duplicate
    }

    #[test]
    fn test_mark_duplicates_size_3_first_highest() {
        // Fast path for size 3: first template has highest score
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[15, 15, 15, 15]); // Score: 60
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 40
        let mut t3 = create_test_template("q3", &[5, 5, 5, 5]); // Score: 20
        let mut templates: Vec<&mut Template> = vec![&mut t1, &mut t2, &mut t3];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.unique_templates, 1);
        assert_eq!(metrics.duplicate_templates, 2);
        assert!(!is_duplicate(&t1)); // Highest score
        assert!(is_duplicate(&t2));
        assert!(is_duplicate(&t3));
    }

    #[test]
    fn test_mark_duplicates_size_3_second_highest() {
        // Fast path for size 3: second template has highest score
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[5, 5, 5, 5]); // Score: 20
        let mut t2 = create_test_template("q2", &[15, 15, 15, 15]); // Score: 60
        let mut t3 = create_test_template("q3", &[10, 10, 10, 10]); // Score: 40
        let mut templates: Vec<&mut Template> = vec![&mut t1, &mut t2, &mut t3];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.unique_templates, 1);
        assert_eq!(metrics.duplicate_templates, 2);
        assert!(is_duplicate(&t1));
        assert!(!is_duplicate(&t2)); // Highest score
        assert!(is_duplicate(&t3));
    }

    #[test]
    fn test_mark_duplicates_size_3_third_highest() {
        // Fast path for size 3: third template has highest score
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[5, 5, 5, 5]); // Score: 20
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 40
        let mut t3 = create_test_template("q3", &[15, 15, 15, 15]); // Score: 60
        let mut templates: Vec<&mut Template> = vec![&mut t1, &mut t2, &mut t3];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.unique_templates, 1);
        assert_eq!(metrics.duplicate_templates, 2);
        assert!(is_duplicate(&t1));
        assert!(is_duplicate(&t2));
        assert!(!is_duplicate(&t3)); // Highest score
    }

    #[test]
    fn test_mark_duplicates_size_4_general_case() {
        // General case for size 4+: uses Vec and sort
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[5, 5, 5, 5]); // Score: 20
        let mut t2 = create_test_template("q2", &[15, 15, 15, 15]); // Score: 60 (highest)
        let mut t3 = create_test_template("q3", &[10, 10, 10, 10]); // Score: 40
        let mut t4 = create_test_template("q4", &[8, 8, 8, 8]); // Score: 32
        let mut templates: Vec<&mut Template> = vec![&mut t1, &mut t2, &mut t3, &mut t4];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.unique_templates, 1);
        assert_eq!(metrics.duplicate_templates, 3);
        assert!(is_duplicate(&t1));
        assert!(!is_duplicate(&t2)); // Highest score
        assert!(is_duplicate(&t3));
        assert!(is_duplicate(&t4));
    }

    #[test]
    fn test_mark_duplicates_counts_duplicate_reads() {
        // Verify duplicate_reads metric is incremented for each read in duplicate templates
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[15, 15, 15, 15]); // Score: 60 (1 read)
        let mut t2 = create_test_template("q2", &[5, 5, 5, 5]); // Score: 20 (1 read)
        let mut templates: Vec<&mut Template> = vec![&mut t1, &mut t2];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.duplicate_reads, 1); // One read from t2
    }

    // ========================================================================
    // Tie-breaking tests - verify deterministic behavior when scores are equal
    // ========================================================================

    #[test]
    fn test_mark_duplicates_tie_breaking_size_2() {
        // When scores are equal, first template should be kept (deterministic)
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[10, 10, 10, 10]); // Score: 40
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 40 (tie)
        let mut templates: Vec<&mut Template> = vec![&mut t1, &mut t2];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.unique_templates, 1);
        assert_eq!(metrics.duplicate_templates, 1);
        // First template wins tie (s0 >= s1)
        assert!(!is_duplicate(&t1));
        assert!(is_duplicate(&t2));
    }

    #[test]
    fn test_mark_duplicates_tie_breaking_size_3_all_equal() {
        // All three have equal scores - first template should be kept
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[10, 10, 10, 10]); // Score: 40
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 40
        let mut t3 = create_test_template("q3", &[10, 10, 10, 10]); // Score: 40
        let mut templates: Vec<&mut Template> = vec![&mut t1, &mut t2, &mut t3];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.unique_templates, 1);
        assert_eq!(metrics.duplicate_templates, 2);
        // First template wins tie (s0 >= s1 && s0 >= s2)
        assert!(!is_duplicate(&t1));
        assert!(is_duplicate(&t2));
        assert!(is_duplicate(&t3));
    }

    #[test]
    fn test_mark_duplicates_tie_breaking_size_3_second_and_third_tie() {
        // First is lower, second and third tie - second should win
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[5, 5, 5, 5]); // Score: 20
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 40
        let mut t3 = create_test_template("q3", &[10, 10, 10, 10]); // Score: 40 (tie with t2)
        let mut templates: Vec<&mut Template> = vec![&mut t1, &mut t2, &mut t3];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.unique_templates, 1);
        assert_eq!(metrics.duplicate_templates, 2);
        // Second template wins tie (s1 >= s0 && s1 >= s2)
        assert!(is_duplicate(&t1));
        assert!(!is_duplicate(&t2));
        assert!(is_duplicate(&t3));
    }

    #[test]
    fn test_mark_duplicates_tie_breaking_size_4_all_equal() {
        // All four have equal scores - first template should be kept
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[10, 10, 10, 10]); // Score: 40
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 40
        let mut t3 = create_test_template("q3", &[10, 10, 10, 10]); // Score: 40
        let mut t4 = create_test_template("q4", &[10, 10, 10, 10]); // Score: 40
        let mut templates: Vec<&mut Template> = vec![&mut t1, &mut t2, &mut t3, &mut t4];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.unique_templates, 1);
        assert_eq!(metrics.duplicate_templates, 3);
        // First template wins tie (stable sort preserves order)
        assert!(!is_duplicate(&t1));
        assert!(is_duplicate(&t2));
        assert!(is_duplicate(&t3));
        assert!(is_duplicate(&t4));
    }

    #[test]
    fn test_mark_duplicates_tie_breaking_size_4_partial_tie() {
        // Third and fourth tie for highest - third should win
        let mut metrics = DedupMetrics::default();
        let mut t1 = create_test_template("q1", &[5, 5, 5, 5]); // Score: 20
        let mut t2 = create_test_template("q2", &[8, 8, 8, 8]); // Score: 32
        let mut t3 = create_test_template("q3", &[10, 10, 10, 10]); // Score: 40
        let mut t4 = create_test_template("q4", &[10, 10, 10, 10]); // Score: 40 (tie with t3)
        let mut templates: Vec<&mut Template> = vec![&mut t1, &mut t2, &mut t3, &mut t4];
        mark_duplicates_in_family(&mut templates, &mut metrics);

        assert_eq!(metrics.unique_templates, 1);
        assert_eq!(metrics.duplicate_templates, 3);
        // Third template wins tie (stable sort preserves order)
        assert!(is_duplicate(&t1));
        assert!(is_duplicate(&t2));
        assert!(!is_duplicate(&t3));
        assert!(is_duplicate(&t4));
    }

    // ========================================================================
    // filter_template tests
    // ========================================================================

    fn default_filter_config() -> DedupFilterConfig {
        DedupFilterConfig {
            umi_tag: [b'R', b'X'],
            min_mapq: 20,
            include_non_pf: false,
            min_umi_length: None,
        }
    }

    /// Helper to create a mapped template with UMI tag
    fn create_mapped_template_with_umi(name: &str, umi: &str, mapq: u8) -> Template {
        let record = RecordBuilder::new()
            .name(name)
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .mapping_quality(mapq)
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .tag("RX", umi.to_string())
            .build();
        Template::from_records(vec![record]).unwrap()
    }

    #[test]
    fn test_filter_template_accepts_valid_template() {
        let config = default_filter_config();
        let mut metrics = FilterMetrics::new();
        let template = create_mapped_template_with_umi("q1", "ACGTACGT", 30);

        assert!(filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.accepted_templates, 1);
    }

    #[test]
    fn test_filter_template_rejects_low_mapq() {
        let config = default_filter_config(); // min_mapq = 20
        let mut metrics = FilterMetrics::new();
        let template = create_mapped_template_with_umi("q1", "ACGTACGT", 10); // MAPQ 10 < 20

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_poor_alignment, 1);
    }

    #[test]
    fn test_filter_template_rejects_umi_with_n() {
        let config = default_filter_config();
        let mut metrics = FilterMetrics::new();
        let template = create_mapped_template_with_umi("q1", "ACNTACGT", 30); // N in UMI

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_ns_in_umi, 1);
    }

    #[test]
    fn test_filter_template_rejects_short_umi() {
        let mut config = default_filter_config();
        config.min_umi_length = Some(8); // Require 8 bases
        let mut metrics = FilterMetrics::new();
        let template = create_mapped_template_with_umi("q1", "ACGT", 30); // Only 4 bases

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_umi_too_short, 1);
    }

    #[test]
    fn test_filter_template_rejects_missing_umi_tag() {
        let config = default_filter_config();
        let mut metrics = FilterMetrics::new();

        // Create template without RX tag
        let record = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .mapping_quality(30)
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .build();
        let template = Template::from_records(vec![record]).unwrap();

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_poor_alignment, 1);
    }

    #[test]
    fn test_filter_template_rejects_qc_fail() {
        let config = default_filter_config(); // include_non_pf = false
        let mut metrics = FilterMetrics::new();

        let record = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .mapping_quality(30)
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::QC_FAIL)
            .tag("RX", "ACGTACGT".to_string())
            .build();
        let template = Template::from_records(vec![record]).unwrap();

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_non_pf, 1);
    }

    #[test]
    fn test_filter_template_accepts_qc_fail_when_included() {
        let mut config = default_filter_config();
        config.include_non_pf = true; // Include QC-fail reads
        let mut metrics = FilterMetrics::new();

        let record = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .mapping_quality(30)
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::QC_FAIL)
            .tag("RX", "ACGTACGT".to_string())
            .build();
        let template = Template::from_records(vec![record]).unwrap();

        assert!(filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.accepted_templates, 1);
    }

    #[test]
    fn test_filter_template_rejects_unmapped() {
        let config = default_filter_config();
        let mut metrics = FilterMetrics::new();

        let record = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::UNMAPPED)
            .tag("RX", "ACGTACGT".to_string())
            .build();
        let template = Template::from_records(vec![record]).unwrap();

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_poor_alignment, 1);
    }

    #[test]
    fn test_filter_template_accepts_paired_umi_with_dash() {
        // Paired UMIs have format "ACGT-TGCA" with a dash separator
        let config = default_filter_config();
        let mut metrics = FilterMetrics::new();
        let template = create_mapped_template_with_umi("q1", "ACGT-TGCA", 30);

        assert!(filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.accepted_templates, 1);
    }

    // ========================================================================
    // umi_for_read tests
    // ========================================================================

    #[test]
    fn test_umi_for_read_identity_uppercase() {
        let assigner = Strategy::Identity.new_assigner_full(0, 1, 100);
        let result = umi_for_read("ACGTACGT", true, assigner.as_ref()).unwrap();
        assert_eq!(result, "ACGTACGT");
    }

    #[test]
    fn test_umi_for_read_identity_lowercase_gets_uppercased() {
        let assigner = Strategy::Identity.new_assigner_full(0, 1, 100);
        let result = umi_for_read("acgtacgt", true, assigner.as_ref()).unwrap();
        assert_eq!(result, "ACGTACGT");
    }

    #[test]
    fn test_umi_for_read_paired_r1_earlier() {
        let assigner = Strategy::Paired.new_assigner_full(1, 1, 100);
        // With max_mismatches=1, prefix_len=2, so lower_prefix="aa", higher_prefix="bb"
        let result = umi_for_read("ACGT-TGCA", true, assigner.as_ref()).unwrap();
        // is_r1_earlier=true => lower_prefix:parts[0]-higher_prefix:parts[1]
        assert_eq!(result, "aa:ACGT-bb:TGCA");
    }

    #[test]
    fn test_umi_for_read_paired_r2_earlier() {
        let assigner = Strategy::Paired.new_assigner_full(1, 1, 100);
        // With max_mismatches=1, prefix_len=2, so lower_prefix="aa", higher_prefix="bb"
        let result = umi_for_read("ACGT-TGCA", false, assigner.as_ref()).unwrap();
        // is_r1_earlier=false => higher_prefix:parts[0]-lower_prefix:parts[1]
        assert_eq!(result, "bb:ACGT-aa:TGCA");
    }

    #[test]
    fn test_umi_for_read_paired_missing_dash_error() {
        let assigner = Strategy::Paired.new_assigner_full(1, 1, 100);
        let result = umi_for_read("ACGTACGT", true, assigner.as_ref());
        assert!(result.is_err());
        assert!(
            result.unwrap_err().to_string().contains("did not contain 2 segments"),
            "Error message should mention missing segments"
        );
    }

    // ========================================================================
    // get_pair_orientation tests
    // ========================================================================

    #[test]
    fn test_get_pair_orientation_both_forward() {
        // Both R1 and R2 on forward strand (no REVERSE_COMPLEMENTED flag)
        let r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .mapping_quality(30)
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .tag("RX", "ACGT-TGCA".to_string())
            .build();
        let r2 = RecordBuilder::new()
            .name("q1")
            .sequence("TGCA")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("4M")
            .mapping_quality(30)
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .tag("RX", "ACGT-TGCA".to_string())
            .build();
        let template = Template::from_records(vec![r1, r2]).unwrap();
        let (r1_positive, r2_positive) = get_pair_orientation(&template);
        assert!(r1_positive);
        assert!(r2_positive);
    }

    #[test]
    fn test_get_pair_orientation_r1_reverse() {
        // R1 on reverse strand, R2 on forward strand
        let r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .mapping_quality(30)
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::REVERSE_COMPLEMENTED)
            .tag("RX", "ACGT-TGCA".to_string())
            .build();
        let r2 = RecordBuilder::new()
            .name("q1")
            .sequence("TGCA")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("4M")
            .mapping_quality(30)
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .tag("RX", "ACGT-TGCA".to_string())
            .build();
        let template = Template::from_records(vec![r1, r2]).unwrap();
        let (r1_positive, r2_positive) = get_pair_orientation(&template);
        assert!(!r1_positive);
        assert!(r2_positive);
    }

    #[test]
    fn test_get_pair_orientation_both_reverse() {
        // Both R1 and R2 on reverse strand
        let r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .mapping_quality(30)
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::REVERSE_COMPLEMENTED)
            .tag("RX", "ACGT-TGCA".to_string())
            .build();
        let r2 = RecordBuilder::new()
            .name("q1")
            .sequence("TGCA")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("4M")
            .mapping_quality(30)
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT | Flags::REVERSE_COMPLEMENTED)
            .tag("RX", "ACGT-TGCA".to_string())
            .build();
        let template = Template::from_records(vec![r1, r2]).unwrap();
        let (r1_positive, r2_positive) = get_pair_orientation(&template);
        assert!(!r1_positive);
        assert!(!r2_positive);
    }

    // ========================================================================
    // is_r1_genomically_earlier tests
    // ========================================================================

    #[test]
    fn test_is_r1_earlier_true() {
        // R1 at position 100, R2 at position 200
        let r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .build();
        let r2 = RecordBuilder::new()
            .name("q1")
            .sequence("TGCA")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("4M")
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .build();
        assert!(is_r1_genomically_earlier(&r1, &r2).unwrap());
    }

    #[test]
    fn test_is_r1_earlier_false() {
        // R1 at position 200, R2 at position 100
        let r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("4M")
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .build();
        let r2 = RecordBuilder::new()
            .name("q1")
            .sequence("TGCA")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .build();
        assert!(!is_r1_genomically_earlier(&r1, &r2).unwrap());
    }

    #[test]
    fn test_is_r1_earlier_equal_position() {
        // Both at position 100 -> true (<=)
        let r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .build();
        let r2 = RecordBuilder::new()
            .name("q1")
            .sequence("TGCA")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .build();
        assert!(is_r1_genomically_earlier(&r1, &r2).unwrap());
    }

    // ========================================================================
    // truncate_umis tests
    // ========================================================================

    #[test]
    fn test_truncate_umis_no_min_length() {
        let umis = vec!["ACGTACGT".to_string(), "TGCATGCA".to_string()];
        let result = truncate_umis(umis.clone(), None).unwrap();
        assert_eq!(result, umis);
    }

    #[test]
    fn test_truncate_umis_truncates() {
        let umis = vec!["ACGTACGT".to_string(), "TGCATGCA".to_string()];
        let result = truncate_umis(umis, Some(4)).unwrap();
        assert_eq!(result, vec!["ACGT".to_string(), "TGCA".to_string()]);
    }

    #[test]
    fn test_truncate_umis_error_too_short() {
        let umis = vec!["ACG".to_string(), "TGCATGCA".to_string()];
        let result = truncate_umis(umis, Some(4));
        assert!(result.is_err());
        assert!(
            result.unwrap_err().to_string().contains("shorter than expected"),
            "Error message should mention UMI being too short"
        );
    }

    // ========================================================================
    // set_mi_tag_on_record tests
    // ========================================================================

    #[test]
    fn test_set_mi_tag_assigned() {
        let mut record =
            RecordBuilder::new().name("q1").sequence("ACGT").qualities(&[30, 30, 30, 30]).build();
        let mi = MoleculeId::Single(42);
        let assign_tag = Tag::from([b'M', b'I']);
        set_mi_tag_on_record(&mut record, mi, assign_tag, 0);

        let mi_value = record.data().get(&assign_tag);
        assert!(mi_value.is_some());
        if let Some(DataValue::String(val)) = mi_value {
            assert_eq!(AsRef::<[u8]>::as_ref(val), b"42");
        } else {
            panic!("MI tag should be a string value");
        }
    }

    #[test]
    fn test_set_mi_tag_unassigned() {
        let mut record =
            RecordBuilder::new().name("q1").sequence("ACGT").qualities(&[30, 30, 30, 30]).build();
        let mi = MoleculeId::None;
        let assign_tag = Tag::from([b'M', b'I']);
        set_mi_tag_on_record(&mut record, mi, assign_tag, 0);

        let mi_value = record.data().get(&assign_tag);
        assert!(mi_value.is_none(), "Unassigned MoleculeId should not set MI tag");
    }

    #[test]
    fn test_set_mi_tag_with_base_offset() {
        let mut record =
            RecordBuilder::new().name("q1").sequence("ACGT").qualities(&[30, 30, 30, 30]).build();
        let mi = MoleculeId::Single(42);
        let assign_tag = Tag::from([b'M', b'I']);
        set_mi_tag_on_record(&mut record, mi, assign_tag, 100);

        let mi_value = record.data().get(&assign_tag);
        assert!(mi_value.is_some());
        if let Some(DataValue::String(val)) = mi_value {
            // 100 (base) + 42 (id) = 142
            assert_eq!(AsRef::<[u8]>::as_ref(val), b"142");
        } else {
            panic!("MI tag should be a string value");
        }
    }

    // ========================================================================
    // filter_template for paired reads tests
    // ========================================================================

    /// Helper to create a paired mapped template with UMI tag on both reads.
    fn create_paired_mapped_template_with_umi(
        name: &str,
        umi: &str,
        r1_mapq: u8,
        r2_mapq: u8,
    ) -> Template {
        let r1 = RecordBuilder::new()
            .name(name)
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .mapping_quality(r1_mapq)
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .tag("RX", umi.to_string())
            .build();
        let r2 = RecordBuilder::new()
            .name(name)
            .sequence("TGCA")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("4M")
            .mapping_quality(r2_mapq)
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .tag("RX", umi.to_string())
            .build();
        Template::from_records(vec![r1, r2]).unwrap()
    }

    #[test]
    fn test_filter_paired_template_accepts_valid() {
        let config = default_filter_config();
        let mut metrics = FilterMetrics::new();
        let template = create_paired_mapped_template_with_umi("q1", "ACGTACGT", 30, 30);

        assert!(filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.accepted_templates, 1);
    }

    #[test]
    fn test_filter_paired_template_rejects_r2_low_mapq() {
        let config = default_filter_config(); // min_mapq = 20
        let mut metrics = FilterMetrics::new();
        let template = create_paired_mapped_template_with_umi("q1", "ACGTACGT", 30, 10);

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_poor_alignment, 1);
    }

    #[test]
    fn test_filter_paired_template_rejects_r2_umi_with_n() {
        let config = default_filter_config();
        let mut metrics = FilterMetrics::new();
        // R1 has valid UMI, but R2 has N in UMI
        let r1 = RecordBuilder::new()
            .name("q1")
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("4M")
            .mapping_quality(30)
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .tag("RX", "ACGTACGT".to_string())
            .build();
        let r2 = RecordBuilder::new()
            .name("q1")
            .sequence("TGCA")
            .qualities(&[30, 30, 30, 30])
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("4M")
            .mapping_quality(30)
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .tag("RX", "ACNTACGT".to_string())
            .build();
        let template = Template::from_records(vec![r1, r2]).unwrap();

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_ns_in_umi, 1);
    }

    #[test]
    fn test_filter_paired_no_reads_rejected() {
        // Template with no primary reads (empty records list)
        let config = default_filter_config();
        let mut metrics = FilterMetrics::new();
        let template = Template::new(b"empty".to_vec());

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_poor_alignment, 1);
    }

    // ========================================================================
    // DedupMetrics tests
    // ========================================================================

    #[test]
    fn test_metrics_merge_all_fields() {
        let mut m1 = DedupMetrics {
            total_templates: 10,
            duplicate_templates: 2,
            unique_templates: 8,
            total_reads: 20,
            duplicate_reads: 4,
            unique_reads: 16,
            secondary_reads: 3,
            supplementary_reads: 1,
            missing_pa_tag: 1,
            ..Default::default()
        };

        let m2 = DedupMetrics {
            total_templates: 5,
            duplicate_templates: 1,
            unique_templates: 4,
            total_reads: 10,
            duplicate_reads: 2,
            unique_reads: 8,
            secondary_reads: 2,
            supplementary_reads: 3,
            missing_pa_tag: 2,
            ..Default::default()
        };

        m1.merge(&m2);
        assert_eq!(m1.total_templates, 15);
        assert_eq!(m1.duplicate_templates, 3);
        assert_eq!(m1.unique_templates, 12);
        assert_eq!(m1.total_reads, 30);
        assert_eq!(m1.duplicate_reads, 6);
        assert_eq!(m1.unique_reads, 24);
        assert_eq!(m1.secondary_reads, 5);
        assert_eq!(m1.supplementary_reads, 4);
        assert_eq!(m1.missing_pa_tag, 3);
    }

    #[test]
    fn test_duplicate_rate_all_duplicates() {
        let metrics =
            DedupMetrics { total_templates: 50, duplicate_templates: 50, ..Default::default() };
        assert!((metrics.duplicate_rate() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_metrics_default() {
        let metrics = DedupMetrics::default();
        assert_eq!(metrics.total_templates, 0);
        assert_eq!(metrics.duplicate_templates, 0);
        assert_eq!(metrics.unique_templates, 0);
        assert_eq!(metrics.total_reads, 0);
        assert_eq!(metrics.duplicate_reads, 0);
        assert_eq!(metrics.unique_reads, 0);
        assert_eq!(metrics.secondary_reads, 0);
        assert_eq!(metrics.supplementary_reads, 0);
        assert_eq!(metrics.missing_pa_tag, 0);
    }

    #[test]
    fn test_dedup_metrics_output_from() {
        let metrics = DedupMetrics {
            total_templates: 100,
            duplicate_templates: 25,
            unique_templates: 75,
            total_reads: 200,
            duplicate_reads: 50,
            unique_reads: 150,
            secondary_reads: 10,
            supplementary_reads: 5,
            missing_pa_tag: 2,
            ..Default::default()
        };
        let output = DedupMetricsOutput::from(&metrics);

        assert_eq!(output.total_templates, 100);
        assert_eq!(output.duplicate_templates, 25);
        assert_eq!(output.unique_templates, 75);
        assert!((output.duplicate_rate - 0.25).abs() < 0.001);
        assert_eq!(output.total_reads, 200);
        assert_eq!(output.duplicate_reads, 50);
        assert_eq!(output.unique_reads, 150);
        assert_eq!(output.secondary_reads, 10);
        assert_eq!(output.supplementary_reads, 5);
        assert_eq!(output.missing_pa_tag, 2);
    }
}
