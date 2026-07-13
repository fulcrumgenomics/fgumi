//! UMI-aware duplicate marking command.
//!
//! This command marks or removes PCR duplicates using UMI information.
//! It operates on template-coordinate sorted BAM files and requires
//! the `tc` tag on secondary/supplementary reads (added by `fgumi zipper`).
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
use std::sync::atomic::{AtomicU64, Ordering};

use crate::assigner::{PairedUmiAssigner, Strategy, UmiAssigner};
use crate::grouper::{
    FilterMetrics, RawPositionGroup, RecordPositionGrouper, build_templates_from_records,
};
use crate::logging::OperationTimer;
use crate::metrics::group::FamilySizeMetrics;
use crate::read_info::LibraryIndex;
use crate::sam::SamTag;
use crate::sam::is_template_coordinate_sorted;
use crate::template::Template;
use crate::umi::{UmiValidation, validate_umi};
use crate::unified_pipeline::{
    BatchWeight, GroupKeyConfig, Grouper, MemoryEstimate,
    run_bam_pipeline_from_reader_with_mi_assign,
};
use ahash::AHashMap;
use anyhow::{Context, Result, bail};
use clap::Parser;
use fgoxide::io::DelimFile;
use fgumi_bam_io::create_bam_reader_for_pipeline_with_opts;

use log::info;
use noodles::sam::Header;
use noodles::sam::alignment::record::data::field::Tag;
use parking_lot::Mutex;
use serde::{Deserialize, Serialize};

use crate::commands::command::Command;
use crate::commands::common::{
    BamIoOptions, CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
    build_pipeline_config, is_r1_genomically_earlier_raw, parse_bool,
};
use crate::sam::TC_TAG;
use fgumi_raw_bam;
use fgumi_raw_bam::RawRecordView;

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
    /// Secondary/supplementary without `tc` tag
    pub missing_tc_tag: u64,
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
        self.missing_tc_tag += other.missing_tc_tag;
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
    // DXM3-04: format the one fraction column like fgbio's DecimalFormat, matching
    // the float serialization #504 wired onto the fgumi-metrics-crate structs.
    #[serde(with = "fgumi_metrics::float")]
    duplicate_rate: f64,
    total_reads: u64,
    unique_reads: u64,
    duplicate_reads: u64,
    secondary_reads: u64,
    supplementary_reads: u64,
    missing_tc_tag: u64,
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
            missing_tc_tag: m.missing_tc_tag,
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
    /// Number of distinct numeric `MoleculeId`s assigned in this group.
    /// Derived directly from the templates the MI Assign hook will mutate, so
    /// the cumulative MI counter stays in lockstep with what `serialize_fn`
    /// emits even if the family-size histogram semantics change.
    pub distinct_mi_count: u64,
}

impl BatchWeight for ProcessedDedupGroup {
    fn batch_weight(&self) -> usize {
        self.templates.len()
    }
}

impl MemoryEstimate for ProcessedDedupGroup {
    fn estimate_heap_size(&self) -> usize {
        self.templates.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
            + self.templates.capacity() * std::mem::size_of::<Template>()
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
    /// Skip UMI validation (position-only grouping).
    no_umi: bool,
}

//////////////////////////////////////////////////////////////////////////////
// Duplicate scoring
//////////////////////////////////////////////////////////////////////////////

/// Minimum base quality counted by Picard's `SUM_OF_BASE_QUALITIES` scoring.
///
/// HTSJDK `DuplicateScoringStrategy.getSumOfBaseQualities` sums only base
/// qualities `>= 15`; anything below is excluded. This is a threshold, not a cap.
const PICARD_MIN_BASE_QUALITY: u8 = 15;

/// Per-read cap on the summed base-quality score, mirroring HTSJDK
/// `DuplicateScoringStrategy.computeDuplicateScore`, which caps each read at
/// `Short.MAX_VALUE / 2` (`Math.min(getSumOfBaseQualities(rec), Short.MAX_VALUE/2)`)
/// so the two mates' scores can be summed without overflowing a `short`.
const PICARD_MAX_SCORE_PER_READ: i64 = i16::MAX as i64 / 2; // 16383

/// Score penalty for a read failing vendor quality checks (QC-fail), mirroring
/// HTSJDK `computeDuplicateScore`'s `Short.MIN_VALUE / 2` discount so a QC-fail
/// read is never preferred as the duplicate representative. Only reachable with
/// `--include-non-pf` (QC-fail reads are otherwise filtered before scoring).
const PICARD_QC_FAIL_DISCOUNT: i64 = i16::MIN as i64 / 2; // -16384

/// Scores a template for duplicate selection.
/// Higher score = more likely to be chosen as representative.
///
/// Implements Picard/HTSJDK's `SUM_OF_BASE_QUALITIES`
/// `DuplicateScoringStrategy.computeDuplicateScore`, summed over the primary reads
/// (R1 and R2) exactly as Picard `MarkDuplicates` sums the two mates:
/// 1. per read, sum every base quality `>= 15` at its **full value**
///    (`getSumOfBaseQualities`; the 15 is a threshold, not a cap),
/// 2. cap that per-read sum at `Short.MAX_VALUE / 2` (16383), and
/// 3. subtract `Short.MIN_VALUE / 2` (16384) from any read flagged QC-fail.
///
/// No mapping-quality term is added: fgumi `dedup` is a Picard `MarkDuplicates`
/// drop-in, whereas fgbio's `GroupReadsByUmi` additionally adds `mapq`. Omitting
/// `mapq` is intentional (see tracker DEDUP-01).
#[inline]
fn score_template(template: &Template) -> i64 {
    let mut score: i64 = 0;

    for raw in [template.r1(), template.r2()].into_iter().flatten() {
        if raw.len() < 32 {
            continue;
        }
        let qo = fgumi_raw_bam::qual_offset(raw);
        let seq_len = fgumi_raw_bam::l_seq(raw) as usize;
        if qo >= raw.len() {
            continue;
        }
        let max_len = std::cmp::min(seq_len, raw.len() - qo);
        // Slice-iterate for auto-vectorization (compiler can vectorize the
        // threshold filter + sum).
        let read_sum: u32 = raw[qo..qo + max_len]
            .iter()
            .filter(|&&q| q >= PICARD_MIN_BASE_QUALITY)
            .map(|&q| u32::from(q))
            .sum();
        // Cap per read (Short.MAX_VALUE/2), then discount QC-fail reads, matching
        // HTSJDK computeDuplicateScore before the two mates' scores are summed.
        let mut read_score = i64::from(read_sum).min(PICARD_MAX_SCORE_PER_READ);
        if (RawRecordView::new(raw).flags() & fgumi_raw_bam::flags::QC_FAIL) != 0 {
            read_score += PICARD_QC_FAIL_DISCOUNT;
        }
        score += read_score;
    }

    score
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
    metrics.total_templates += 1;

    // Fail closed when *any* record (primary, secondary, or supplementary) is
    // shorter than the minimum BAM record length: a truncated record indicates
    // corrupt input, so the template must be rejected rather than have the
    // malformed record silently dropped (and later panic in RawRecordView::new
    // on the dedup/serialize path).
    if template.records().iter().any(|r| r.len() < fgumi_raw_bam::MIN_BAM_RECORD_LEN) {
        metrics.discarded_poor_alignment += 1;
        return false;
    }

    let raw_r1 = template.r1();
    let raw_r2 = template.r2();

    if raw_r1.is_none() && raw_r2.is_none() {
        metrics.discarded_poor_alignment += 1;
        return false;
    }

    let both_unmapped = raw_r1
        .is_none_or(|r| (RawRecordView::new(r).flags() & fgumi_raw_bam::flags::UNMAPPED) != 0)
        && raw_r2
            .is_none_or(|r| (RawRecordView::new(r).flags() & fgumi_raw_bam::flags::UNMAPPED) != 0);
    if both_unmapped {
        metrics.discarded_poor_alignment += 1;
        return false;
    }

    // Phase 1: Flag-based checks
    for raw in [raw_r1, raw_r2].into_iter().flatten() {
        let flg = RawRecordView::new(raw).flags();

        if !config.include_non_pf && (flg & fgumi_raw_bam::flags::QC_FAIL) != 0 {
            metrics.discarded_non_pf += 1;
            return false;
        }

        if (flg & fgumi_raw_bam::flags::UNMAPPED) == 0 {
            let mapq = fgumi_raw_bam::mapq(raw);
            if mapq < config.min_mapq {
                metrics.discarded_poor_alignment += 1;
                return false;
            }
        }
    }

    // Phase 2: Single-pass tag lookups (MQ + UMI in one aux scan)
    for raw in [raw_r1, raw_r2].into_iter().flatten() {
        let flg = RawRecordView::new(raw).flags();
        let aux = fgumi_raw_bam::aux_data_slice(raw);
        let check_mq = (flg & fgumi_raw_bam::flags::MATE_UNMAPPED) == 0;
        let check_umi = !config.no_umi;

        let mut found_mq: Option<i64> = None;
        let mut found_umi: Option<&[u8]> = None;
        let mut p = 0;
        while p + 3 <= aux.len() {
            let t = [aux[p], aux[p + 1]];
            let val_type = aux[p + 2];

            if check_umi && t == config.umi_tag && val_type == b'Z' {
                let start = p + 3;
                if let Some(end) = aux[start..].iter().position(|&b| b == 0) {
                    found_umi = Some(&aux[start..start + end]);
                    p = start + end + 1;
                } else {
                    break;
                }
                if !check_mq || found_mq.is_some() {
                    break;
                }
                continue;
            }

            if check_mq && t == *SamTag::MQ {
                found_mq = match val_type {
                    b'C' if p + 3 < aux.len() => Some(i64::from(aux[p + 3])),
                    b'c' if p + 3 < aux.len() => Some(i64::from(aux[p + 3] as i8)),
                    b'S' if p + 5 <= aux.len() => {
                        Some(i64::from(u16::from_le_bytes([aux[p + 3], aux[p + 4]])))
                    }
                    b's' if p + 5 <= aux.len() => {
                        Some(i64::from(i16::from_le_bytes([aux[p + 3], aux[p + 4]])))
                    }
                    b'I' if p + 7 <= aux.len() => Some(i64::from(u32::from_le_bytes([
                        aux[p + 3],
                        aux[p + 4],
                        aux[p + 5],
                        aux[p + 6],
                    ]))),
                    b'i' if p + 7 <= aux.len() => Some(i64::from(i32::from_le_bytes([
                        aux[p + 3],
                        aux[p + 4],
                        aux[p + 5],
                        aux[p + 6],
                    ]))),
                    _ => None,
                };
            }

            if let Some(size) = fgumi_raw_bam::tag_value_size(val_type, &aux[p + 3..]) {
                p += 3 + size;
            } else {
                break;
            }
            if (!check_umi || found_umi.is_some()) && (!check_mq || found_mq.is_some()) {
                break;
            }
        }

        if check_mq {
            if let Some(mq) = found_mq {
                // Compare as signed so a negative MQ (e.g., MQ:c:-1) fails the filter
                // rather than wrapping to 255 via `as u8`.
                if mq < i64::from(config.min_mapq) {
                    metrics.discarded_poor_alignment += 1;
                    return false;
                }
            }
        }

        // Skip UMI validation in no-umi mode
        if config.no_umi {
            continue;
        }

        if let Some(umi_bytes) = found_umi {
            match validate_umi(umi_bytes) {
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
    let r1_positive = template
        .r1()
        .is_none_or(|r| (RawRecordView::new(r).flags() & fgumi_raw_bam::flags::REVERSE) == 0);
    let r2_positive = template
        .r2()
        .is_none_or(|r| (RawRecordView::new(r).flags() & fgumi_raw_bam::flags::REVERSE) == 0);
    (r1_positive, r2_positive)
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
    raw_tag: SamTag,
    min_umi_length: Option<usize>,
    no_umi: bool,
) -> Result<()> {
    if indices.is_empty() {
        return Ok(());
    }

    let mut umis = Vec::with_capacity(indices.len());

    for &idx in indices {
        let template = &templates[idx];

        // In no-umi mode, use empty string for all templates
        let processed_umi = if no_umi {
            String::new()
        } else {
            // Prefer the UMI position cached during the parallel Decode step
            // (enabled via GroupKeyConfig::with_umi_tag). The fallback exists
            // for templates built outside the Decode path and for cases where
            // the cache was disabled.
            let umi_bytes = if let Some(cached) = template.cached_umi() {
                cached
            } else if let Some(r1_raw) = template.r1() {
                let aux = fgumi_raw_bam::aux_data_slice(r1_raw);
                fgumi_raw_bam::find_string_tag(aux, raw_tag)
                    .ok_or_else(|| anyhow::anyhow!("Missing UMI tag"))?
            } else if let Some(r2_raw) = template.r2() {
                let aux = fgumi_raw_bam::aux_data_slice(r2_raw);
                fgumi_raw_bam::find_string_tag(aux, raw_tag)
                    .ok_or_else(|| anyhow::anyhow!("Missing UMI tag"))?
            } else {
                bail!("Template has no reads");
            };
            let umi_str = std::str::from_utf8(umi_bytes)
                .map_err(|e| anyhow::anyhow!("Invalid UTF-8: {e}"))?;
            let is_r1_earlier = if let (Some(r1), Some(r2)) = (template.r1(), template.r2()) {
                is_r1_genomically_earlier_raw(r1, r2)
            } else {
                true
            };

            umi_for_read(umi_str, is_r1_earlier, assigner)?
        };

        umis.push(processed_umi);
    }

    // Truncate UMIs if needed (skip in no-umi mode)
    let truncated_umis = if no_umi { umis } else { truncate_umis(umis, min_umi_length)? };
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
    raw_tag: SamTag,
    min_umi_length: Option<usize>,
    no_umi: bool,
) -> Result<()> {
    if assigner.split_templates_by_pair_orientation() {
        let mut subgroups: AHashMap<(bool, bool), Vec<usize>> = AHashMap::new();
        for (idx, template) in templates.iter().enumerate() {
            let orientation = get_pair_orientation(template);
            subgroups.entry(orientation).or_default().push(idx);
        }

        // Iterate orientation subgroups in a deterministic order.
        // `AHashMap`'s per-process random hasher would otherwise walk
        // FR/RF subgroups in different orders across runs and assign
        // different local `MoleculeId`s to the same templates.
        let mut ordered_subgroups: Vec<((bool, bool), Vec<usize>)> =
            subgroups.into_iter().collect();
        ordered_subgroups.sort_by_key(|(orientation, _)| *orientation);
        for (_orientation, indices) in &ordered_subgroups {
            assign_umi_groups_for_indices(
                templates,
                indices,
                assigner,
                raw_tag,
                min_umi_length,
                no_umi,
            )?;
        }
    } else {
        let all_indices: Vec<usize> = (0..templates.len()).collect();
        assign_umi_groups_for_indices(
            templates,
            &all_indices,
            assigner,
            raw_tag,
            min_umi_length,
            no_umi,
        )?;
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
    for raw in template.records_mut().iter_mut() {
        let flg = RawRecordView::new(raw.as_ref()).flags();
        fgumi_raw_bam::set_flags(raw, flg | DUPLICATE_FLAG);
        dedup_metrics.duplicate_reads += 1;
    }
}

//////////////////////////////////////////////////////////////////////////////
// Position group processing
//////////////////////////////////////////////////////////////////////////////

/// Process a position group for deduplication.
fn process_position_group(
    group: RawPositionGroup,
    filter_config: &DedupFilterConfig,
    assigner: &dyn UmiAssigner,
    raw_tag: SamTag,
    min_umi_length: Option<usize>,
    no_umi: bool,
) -> io::Result<ProcessedDedupGroup> {
    let mut dedup_metrics = DedupMetrics::default();
    let input_record_count = group.records.len() as u64;

    // Build templates from raw records (deferred from Group step)
    let all_templates = build_templates_from_records(group.records)?;

    // Filter templates
    let mut filter_metrics = FilterMetrics::new();
    let filtered_templates: Vec<Template> = all_templates
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
            distinct_mi_count: 0,
        });
    }

    // Clear existing duplicate flags before re-marking
    let mut templates: Vec<Template> = filtered_templates
        .into_iter()
        .map(|mut t| {
            for raw in t.records_mut().iter_mut() {
                let flg = RawRecordView::new(raw.as_ref()).flags();
                fgumi_raw_bam::set_flags(raw, flg & !DUPLICATE_FLAG);
            }
            t
        })
        .collect();
    if let Err(e) = assign_umi_groups(&mut templates, assigner, raw_tag, min_umi_length, no_umi) {
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

    // Count reads and check for missing tc tags
    let tc_tag_bytes: [u8; 2] = *TC_TAG.as_ref();
    for template in &templates {
        dedup_metrics.total_templates += 1;
        for raw in template.records() {
            dedup_metrics.total_reads += 1;
            let flg = RawRecordView::new(raw.as_ref()).flags();
            let is_secondary = (flg & fgumi_raw_bam::flags::SECONDARY) != 0;
            let is_supplementary = (flg & fgumi_raw_bam::flags::SUPPLEMENTARY) != 0;

            if is_secondary {
                dedup_metrics.secondary_reads += 1;
            }
            if is_supplementary {
                dedup_metrics.supplementary_reads += 1;
            }
            if is_secondary || is_supplementary {
                let aux = fgumi_raw_bam::aux_data_slice(raw);
                if fgumi_raw_bam::find_tag_type(aux, tc_tag_bytes).is_none() {
                    dedup_metrics.missing_tc_tag += 1;
                }
            }
        }
    }

    dedup_metrics.unique_reads = dedup_metrics.total_reads - dedup_metrics.duplicate_reads;

    // Compute the number of distinct numeric molecule IDs assigned in this group.
    // Mirrors `fgumi group`: assigners hand out numeric IDs 0, 1, 2, ... contiguously,
    // so `max(id) + 1` equals the count of distinct IDs used. Derived from the
    // templates the MI Assign hook will mutate, which decouples the cumulative MI
    // counter from how `family_sizes` is populated.
    let distinct_mi_count: u64 =
        templates.iter().filter_map(|t| t.mi.id()).max().map(|max_id| max_id + 1).unwrap_or(0);

    Ok(ProcessedDedupGroup {
        templates,
        family_sizes,
        dedup_metrics,
        input_record_count,
        distinct_mi_count,
    })
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
Requires template-coordinate sorted input with `tc` tags on secondary/supplementary
reads (added by `fgumi zipper`).

Within each UMI family, the template with the highest sum of base qualities
is selected as the representative; all others are marked as duplicates.

# Input Requirements

- Must be processed with `fgumi zipper` (adds `tc` tag for secondary/supplementary reads)
- Must be sorted with `fgumi sort --order template-coordinate`
- UMI tags on reads (RX tag), unless `--no-umi` is specified

Note: Using `samtools sort` will NOT work correctly because it doesn't use the
`tc` tag for template-coordinate ordering of secondary/supplementary reads.

# Output Modes

- Mark only (default): Set duplicate flag (0x400) on non-representative reads
- Remove (--remove-duplicates): Exclude duplicate reads from output entirely

# Cell Barcodes

If the input data contains cell barcodes (e.g. from single-cell sequencing), reads at the same
genomic position are partitioned by cell barcode before deduplication. This ensures that reads from
different cells are never marked as duplicates of each other, even if they share a UMI sequence and
mapping position. The cell barcode is read from the standard `CB` tag. No
correction or error-handling is performed on cell barcodes; they must be corrected upstream.

# Memory

dedup streams the input and processes one genomic position group at a time, so memory scales with
parallelism, not input size. The two contributors are:

  - Pipeline queue: bounded by --max-memory, which is per thread by default. With the default
    (768 MiB/thread), --threads 16 budgets ~12 GiB of queue alone and --threads 8 budgets ~6 GiB.
  - Per-worker working set: each worker transiently holds the templates of the group it is
    processing (building, scoring, sorting). This is on top of the queue budget, so peak RSS is
    higher than --max-memory; empirically ~1-2 GiB/thread on WGS-scale input.

A practical rule of thumb for WGS-scale input is to budget ~2 GiB per thread of total RSS. On a
fixed-RAM host, prefer one of:

  fgumi dedup ... --threads 16 --max-memory auto              # detect host RAM, self-throttle
  fgumi dedup ... --threads 16 --max-memory 16G --memory-per-thread false   # fixed total budget

--max-memory is a budget for the controllable consumer (the queue), not a hard RSS cap: a single
pathological high-coverage position group is still processed whole.
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
    #[arg(short = 'r', long = "remove-duplicates", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub remove_duplicates: bool,

    /// Minimum mapping quality for a read to be included
    #[arg(short = 'q', long = "min-map-q")]
    pub min_map_q: Option<u8>,

    /// Include reads flagged as not passing QC
    #[arg(short = 'n', long = "include-non-pf-reads", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
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

    /// Skip UMI-based grouping; group by position only. Forces identity strategy
    /// and ignores any existing UMI tags.
    #[arg(long = "no-umi", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub no_umi: bool,

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

        // Validate --no-umi is not used with paired strategy
        if self.no_umi && matches!(self.strategy, Strategy::Paired) {
            bail!("--no-umi cannot be used with --strategy paired");
        }

        // Handle --no-umi mode: force identity strategy
        let (effective_strategy, no_umi_edits_override) = if self.no_umi {
            if !matches!(self.strategy, Strategy::Identity) {
                info!("--no-umi mode: overriding strategy to identity");
            }
            (Strategy::Identity, true)
        } else {
            (self.strategy, false)
        };

        // Validate the input exists (stdin paths are exempt).
        self.io.validate()?;

        let min_mapq: u8 = self.min_map_q.unwrap_or(0);

        // Identity strategy requires edits=0, others use the configured value
        // Also force edits=0 in no-umi mode
        let effective_edits =
            if no_umi_edits_override || matches!(effective_strategy, Strategy::Identity) {
                0
            } else {
                self.edits
            };

        let timer = OperationTimer::new("Marking duplicates");

        info!("Starting dedup");
        info!("Input: {}", self.io.input.display());
        info!("Output: {}", self.io.output.display());
        info!("Strategy: {effective_strategy:?}");
        info!("Edits: {effective_edits}");
        info!("Remove duplicates: {}", self.remove_duplicates);
        if self.no_umi {
            info!("No-UMI mode: deduplicating by position only");
        }
        if matches!(effective_strategy, Strategy::Adjacency | Strategy::Paired) {
            info!("Index threshold: {}", self.index_threshold);
        }
        info!("{}", self.threading.log_message());

        // Open input BAM
        let (reader, header) = create_bam_reader_for_pipeline_with_opts(
            &self.io.input,
            self.io.pipeline_reader_opts(),
        )?;

        if !is_template_coordinate_sorted(&header) {
            bail!(
                "Input BAM must be template-coordinate sorted (header must advertise \
                 SO:unsorted, GO:query, and SS:template-coordinate).\n\n\
                 To prepare your BAM file, run:\n  \
                 fgumi zipper -i mapped.bam -u unmapped.bam -r reference.fa -o merged.bam\n  \
                 fgumi sort -i merged.bam -o sorted.bam --order template-coordinate"
            );
        }
        info!("Template-coordinate sorted");

        // Add @PG record
        let header = crate::commands::common::add_pg_record(header, command_line)?;

        // Tag constants per SAM specification
        let raw_tag = SamTag::RX;
        let cell_tag = Tag::from(SamTag::CB);
        let assign_tag_bytes: [u8; 2] = *SamTag::MI;

        let filter_config = DedupFilterConfig {
            umi_tag: *raw_tag,
            min_mapq,
            include_non_pf: self.include_non_pf_reads,
            min_umi_length: self.min_umi_length,
            no_umi: self.no_umi,
        };

        // Shared state for collecting metrics
        let collected_metrics: Arc<Mutex<CollectedDedupMetrics>> =
            Arc::new(Mutex::new(CollectedDedupMetrics::default()));

        // Clone values needed by closures
        let strategy = effective_strategy;
        let index_threshold = self.index_threshold;
        let min_umi_length = self.min_umi_length;
        let no_umi = self.no_umi;
        let remove_duplicates = self.remove_duplicates;
        let collected_metrics_clone = Arc::clone(&collected_metrics);

        // Configure pipeline
        let num_threads = self.threading.num_threads();
        let mut pipeline_config = build_pipeline_config(
            &self.scheduler_opts,
            &self.compression,
            &self.queue_memory,
            num_threads,
        )?;
        info!("Scheduler: {:?}", self.scheduler_opts.strategy());
        info!("Using pipeline with {num_threads} threads");

        let library_index = LibraryIndex::from_header(&header);
        // Cache the UMI tag's value position on each DecodedRecord so the
        // Process step's UMI-assignment pass can slice the value without
        // re-scanning aux data (issue #334). In `--no-umi` mode the
        // assignment site emits `String::new()` and never reads the cache,
        // so populating it during decode would be pure overhead.
        let group_key_config = GroupKeyConfig::new(library_index, cell_tag);
        pipeline_config.group_key_config = Some(if self.no_umi {
            group_key_config
        } else {
            group_key_config.with_umi_tag(*raw_tag)
        });

        // Cumulative MoleculeId counter, advanced **only** by the
        // serial-ordered MI Assign hook installed below. Mirrors the
        // pattern used by `fgumi group`; see
        // `docs/design/deterministic-mi-numbering.md` for why this lives
        // in the hook rather than in `serialize_fn`. The `Arc` is owned
        // by the hook closure; nothing outside the closure needs to read
        // its final value.
        let next_mi_base_for_hook = Arc::new(AtomicU64::new(0));

        // Run the pipeline
        let _records_processed = run_bam_pipeline_from_reader_with_mi_assign(
            pipeline_config,
            reader,
            header,
            &self.io.output,
            None,
            // Grouper factory - use with_secondary_supplementary to include all reads
            move |_header: &Header| {
                Box::new(RecordPositionGrouper::with_secondary_supplementary())
                    as Box<dyn Grouper<Group = RawPositionGroup> + Send>
            },
            // Process function (parallel) — builds templates from raw records
            move |group: RawPositionGroup| -> io::Result<ProcessedDedupGroup> {
                let assigner = strategy.new_assigner_full(effective_edits, 1, index_threshold);
                process_position_group(
                    group,
                    &filter_config,
                    assigner.as_ref(),
                    raw_tag,
                    min_umi_length,
                    no_umi,
                )
            },
            // Serialize function (parallel, output ordered by serial numbers)
            move |processed: ProcessedDedupGroup,
                  _header: &Header,
                  output: &mut Vec<u8>|
                  -> io::Result<u64> {
                // Collect metrics
                {
                    let mut agg = collected_metrics_clone.lock();
                    agg.dedup_metrics.merge(&processed.dedup_metrics);
                    for (size, count) in &processed.family_sizes {
                        *agg.family_sizes.entry(*size).or_insert(0) += count;
                    }
                }

                let input_record_count = processed.input_record_count;

                // Templates already carry their final global `MoleculeId`s
                // because the MI Assign hook (installed below) ran in
                // serial order before this closure was called, so we pass
                // `base_mi = 0` to `write_with_offset`.
                // Pre-allocate output buffer: ~2 records/template × ~400 bytes/record
                output.reserve(processed.templates.len() * 2 * 400);
                let mut scratch = Vec::with_capacity(512);
                let mut mi_buf = String::with_capacity(16);
                for template in &processed.templates {
                    let mi = template.mi;
                    let has_mi = mi.is_assigned();
                    if has_mi {
                        mi.write_with_offset(0, &mut mi_buf);
                    }
                    for raw in template.records() {
                        // Skip duplicates if remove mode
                        if remove_duplicates
                            && (RawRecordView::new(raw).flags() & DUPLICATE_FLAG) != 0
                        {
                            continue;
                        }
                        if has_mi {
                            scratch.clear();
                            scratch.extend_from_slice(raw);
                            fgumi_raw_bam::update_string_tag(
                                &mut scratch,
                                assign_tag_bytes,
                                mi_buf.as_bytes(),
                            );
                            let block_size = scratch.len() as u32;
                            output.extend_from_slice(&block_size.to_le_bytes());
                            output.extend_from_slice(&scratch);
                        } else {
                            let block_size = raw.len() as u32;
                            output.extend_from_slice(&block_size.to_le_bytes());
                            output.extend_from_slice(raw);
                        }
                    }
                }

                Ok(input_record_count)
            },
            // mi_assign_fn: called in serial order by the MI Assign zone
            // before each item's `serialize_fn`. Folds a cumulative offset
            // into every template's local `MoleculeId`. Advances the counter
            // by `distinct_mi_count` so the offset matches what `serialize_fn`
            // emits even when families contain multiple templates.
            //
            // `fetch_update` + `checked_add` so wraparound is detected and surfaced
            // rather than silently reusing MI integers (which would defeat
            // `MoleculeId::with_offset`'s own overflow check on the per-template add).
            move |_ord, processed: &mut ProcessedDedupGroup| {
                let base = next_mi_base_for_hook
                    .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |current| {
                        current.checked_add(processed.distinct_mi_count)
                    })
                    .map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "MoleculeId offset overflow: cumulative MI counter exceeded u64::MAX",
                        )
                    })?;
                for template in &mut processed.templates {
                    template.mi = template.mi.with_offset(base);
                }
                Ok(())
            },
        )?;

        // Aggregate metrics
        let aggregated = Arc::try_unwrap(collected_metrics)
            .expect("bug: metrics Arc still shared after pipeline join")
            .into_inner();
        let final_metrics = aggregated.dedup_metrics;
        let final_family_sizes = aggregated.family_sizes;

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

        if final_metrics.missing_tc_tag > 0 {
            bail!(
                "{} secondary/supplementary reads are missing the `tc` tag.\n\n\
                The `tc` tag is required for correct UMI-aware deduplication of \
                secondary and supplementary alignments. This tag is added by \
                `fgumi zipper` during the merge of unmapped and mapped BAMs.\n\n\
                To fix this, re-run your pipeline starting from `fgumi zipper`:\n  \
                fgumi zipper -i aligned.bam --unmapped unmapped.bam -r reference.fa -o merged.bam\n  \
                fgumi sort -i merged.bam -o sorted.bam --order template-coordinate\n  \
                fgumi dedup -i sorted.bam -o deduped.bam",
                final_metrics.missing_tc_tag
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
    let metrics = FamilySizeMetrics::from_size_counts(family_sizes.iter().map(|(&s, &c)| (s, c)));
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
    use fgumi_raw_bam::{RawRecord, SamBuilder as RawSamBuilder, flags, testutil::encode_op};
    use rstest::rstest;

    // Helper to create a simple template for testing
    fn create_test_template(name: &str, qualities: &[u8]) -> Template {
        let mut b = RawSamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(b"ACGT")
            .qualities(qualities)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT);
        Template::from_records(vec![b.build()]).expect("test template construction should not fail")
    }

    /// Picard/HTSJDK `SUM_OF_BASE_QUALITIES` (`DuplicateScoringStrategy.getSumOfBaseQualities`):
    /// each base quality strictly below 15 is excluded; 15 and above are summed at full value.
    #[rstest]
    #[case::all_above_threshold(&[20, 20, 20, 20], 80)] // 4 * 20
    #[case::all_below_threshold(&[10, 10, 10, 10], 0)] // every q < 15 excluded
    #[case::threshold_boundary(&[14, 15, 16, 40], 71)] // 14 excluded; 15 + 16 + 40
    fn test_score_template_single(#[case] qualities: &[u8], #[case] expected: i64) {
        let template = create_test_template("q1", qualities);
        assert_eq!(score_template(&template), expected);
    }

    /// Creates a paired-end test template with R1 and R2.
    fn create_paired_test_template(
        name: &str,
        r1_qualities: &[u8],
        r2_qualities: &[u8],
    ) -> Template {
        let r1 = {
            let mut b = RawSamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGT")
                .qualities(r1_qualities)
                .flags(flags::PAIRED | flags::FIRST_SEGMENT);
            b.build()
        };
        let r2 = {
            let mut b = RawSamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"TGCA")
                .qualities(r2_qualities)
                .flags(flags::PAIRED | flags::LAST_SEGMENT);
            b.build()
        };
        Template::from_records(vec![r1, r2]).expect("test template construction should not fail")
    }

    /// Both R1 and R2 contribute to a paired template's score, applying the same
    /// Picard `>= 15` per-base rule across both mates.
    #[rstest]
    #[case::both_below_threshold(&[10, 10, 10, 10], &[10, 10, 10, 10], 0)]
    #[case::both_high_quality(&[40, 40, 40, 40], &[40, 40, 40, 40], 320)] // 160 + 160, no cap
    #[case::asymmetric_r1_excluded(&[5, 5, 5, 5], &[15, 15, 15, 15], 60)] // R1 all < 15; R2 4*15
    #[case::mixed_threshold(&[10, 20, 5, 30], &[8, 12, 25, 3], 75)] // R1 20+30; R2 25
    fn test_score_template_paired(
        #[case] r1_qualities: &[u8],
        #[case] r2_qualities: &[u8],
        #[case] expected: i64,
    ) {
        let template = create_paired_test_template("q1", r1_qualities, r2_qualities);
        assert_eq!(score_template(&template), expected);
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
        template.records.iter().any(|r| (r.flags() & DUPLICATE_FLAG) != 0)
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
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 0 (all q < 15)
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
        let mut t1 = create_test_template("q1", &[5, 5, 5, 5]); // Score: 0 (all q < 15)
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
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 0 (all q < 15)
        let mut t3 = create_test_template("q3", &[5, 5, 5, 5]); // Score: 0 (all q < 15)
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
        let mut t1 = create_test_template("q1", &[5, 5, 5, 5]); // Score: 0 (all q < 15)
        let mut t2 = create_test_template("q2", &[15, 15, 15, 15]); // Score: 60
        let mut t3 = create_test_template("q3", &[10, 10, 10, 10]); // Score: 0 (all q < 15)
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
        let mut t1 = create_test_template("q1", &[5, 5, 5, 5]); // Score: 0 (all q < 15)
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 0 (all q < 15)
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
        let mut t1 = create_test_template("q1", &[5, 5, 5, 5]); // Score: 0 (all q < 15)
        let mut t2 = create_test_template("q2", &[15, 15, 15, 15]); // Score: 60 (highest)
        let mut t3 = create_test_template("q3", &[10, 10, 10, 10]); // Score: 0 (all q < 15)
        let mut t4 = create_test_template("q4", &[8, 8, 8, 8]); // Score: 0 (all q < 15)
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
        let mut t2 = create_test_template("q2", &[5, 5, 5, 5]); // Score: 0 (all q < 15; 1 read)
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
        let mut t1 = create_test_template("q1", &[10, 10, 10, 10]); // Score: 0 (all < 15)
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 0 (tie)
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
        let mut t1 = create_test_template("q1", &[10, 10, 10, 10]); // Score: 0 (all < 15)
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 0
        let mut t3 = create_test_template("q3", &[10, 10, 10, 10]); // Score: 0
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
        let mut t1 = create_test_template("q1", &[10, 10, 10, 10]); // Score: 0 (all < 15)
        let mut t2 = create_test_template("q2", &[20, 20, 20, 20]); // Score: 80
        let mut t3 = create_test_template("q3", &[20, 20, 20, 20]); // Score: 80 (tie with t2)
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
        let mut t1 = create_test_template("q1", &[10, 10, 10, 10]); // Score: 0 (all < 15)
        let mut t2 = create_test_template("q2", &[10, 10, 10, 10]); // Score: 0
        let mut t3 = create_test_template("q3", &[10, 10, 10, 10]); // Score: 0
        let mut t4 = create_test_template("q4", &[10, 10, 10, 10]); // Score: 0
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
        let mut t1 = create_test_template("q1", &[10, 10, 10, 10]); // Score: 0 (all < 15)
        let mut t2 = create_test_template("q2", &[16, 16, 16, 16]); // Score: 64
        let mut t3 = create_test_template("q3", &[20, 20, 20, 20]); // Score: 80
        let mut t4 = create_test_template("q4", &[20, 20, 20, 20]); // Score: 80 (tie with t3)
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
            umi_tag: *SamTag::RX,
            min_mapq: 20,
            include_non_pf: false,
            min_umi_length: None,
            no_umi: false,
        }
    }

    /// Helper to create a mapped template with UMI tag
    fn create_mapped_template_with_umi(name: &str, umi: &str, mapq: u8) -> Template {
        let mut b = RawSamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30])
            .ref_id(0)
            .pos(99)
            .cigar_ops(&[encode_op(0, 4)])
            .mapq(mapq)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT);
        b.add_string_tag(SamTag::RX, umi.as_bytes());
        let record = b.build();
        Template::from_records(vec![record]).expect("test template construction should not fail")
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

    /// Truncated mate must reject the template (fail closed) rather than
    /// silently dropping the malformed record and proceeding on the remaining
    /// valid mate. A truncated record indicates corrupt input.
    #[test]
    fn test_filter_template_rejects_truncated_mate() {
        let config = default_filter_config();
        let mut metrics = FilterMetrics::new();

        // A valid R1 (full record) and a truncated R2 (16 bytes; below MIN_BAM_RECORD_LEN = 32).
        let mut b = RawSamBuilder::new();
        b.read_name(b"q1")
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30])
            .ref_id(0)
            .pos(99)
            .cigar_ops(&[encode_op(0, 4)])
            .mapq(30)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT);
        b.add_string_tag(SamTag::RX, b"ACGTACGT");
        let valid_r1 = b.build();
        let truncated_r2 = RawRecord::from(vec![0u8; 16]);
        let template = Template {
            name: b"q1".to_vec(),
            records: vec![valid_r1, truncated_r2],
            r1: Some((0, 1)),
            r2: Some((1, 2)),
            r1_supplementals: None,
            r2_supplementals: None,
            r1_secondaries: None,
            r2_secondaries: None,
            mi: crate::umi::MoleculeId::None,
            cached_umi_position: None,
        };

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_poor_alignment, 1, "truncated mate must fail closed");
        assert_eq!(metrics.accepted_templates, 0);
    }

    /// A truncated secondary or supplementary alignment must also reject the
    /// whole template. The previous r1()/r2()-only guard let truncated
    /// non-primary records slip through and later panic in `RawRecordView::new`.
    #[test]
    fn test_filter_template_rejects_truncated_secondary() {
        let config = default_filter_config();
        let mut metrics = FilterMetrics::new();

        let mut b1 = RawSamBuilder::new();
        b1.read_name(b"q1")
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30])
            .ref_id(0)
            .pos(99)
            .cigar_ops(&[encode_op(0, 4)])
            .mapq(30)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT);
        b1.add_string_tag(SamTag::RX, b"ACGTACGT");
        let valid_r1 = b1.build();

        let mut b2 = RawSamBuilder::new();
        b2.read_name(b"q1")
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30])
            .ref_id(0)
            .pos(199)
            .cigar_ops(&[encode_op(0, 4)])
            .mapq(30)
            .flags(flags::PAIRED | flags::LAST_SEGMENT);
        b2.add_string_tag(SamTag::RX, b"ACGTACGT");
        let valid_r2 = b2.build();

        // 16 bytes is below MIN_BAM_RECORD_LEN = 32.
        let truncated_secondary = RawRecord::from(vec![0u8; 16]);

        let template = Template {
            name: b"q1".to_vec(),
            records: vec![valid_r1, valid_r2, truncated_secondary],
            r1: Some((0, 1)),
            r2: Some((1, 2)),
            r1_supplementals: None,
            r2_supplementals: None,
            r1_secondaries: Some((2, 3)),
            r2_secondaries: None,
            mi: crate::umi::MoleculeId::None,
            cached_umi_position: None,
        };

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_poor_alignment, 1, "truncated secondary must fail closed");
        assert_eq!(metrics.accepted_templates, 0);
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
        let record = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(30)
                .flags(flags::PAIRED | flags::FIRST_SEGMENT);
            b.build()
        };
        let template = Template::from_records(vec![record])
            .expect("test template construction should not fail");

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_poor_alignment, 1);
    }

    #[test]
    fn test_filter_template_rejects_qc_fail() {
        let config = default_filter_config(); // include_non_pf = false
        let mut metrics = FilterMetrics::new();

        let record = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(30)
                .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::QC_FAIL);
            b.add_string_tag(SamTag::RX, b"ACGTACGT");
            b.build()
        };
        let template = Template::from_records(vec![record])
            .expect("test template construction should not fail");

        assert!(!filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.discarded_non_pf, 1);
    }

    #[test]
    fn test_filter_template_accepts_qc_fail_when_included() {
        let mut config = default_filter_config();
        config.include_non_pf = true; // Include QC-fail reads
        let mut metrics = FilterMetrics::new();

        let record = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(30)
                .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::QC_FAIL);
            b.add_string_tag(SamTag::RX, b"ACGTACGT");
            b.build()
        };
        let template = Template::from_records(vec![record])
            .expect("test template construction should not fail");

        assert!(filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.accepted_templates, 1);
    }

    #[test]
    fn test_filter_template_rejects_unmapped() {
        let config = default_filter_config();
        let mut metrics = FilterMetrics::new();

        let record = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::UNMAPPED);
            b.add_string_tag(SamTag::RX, b"ACGTACGT");
            b.build()
        };
        let template = Template::from_records(vec![record])
            .expect("test template construction should not fail");

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
        let result =
            umi_for_read("ACGTACGT", true, assigner.as_ref()).expect("valid UMI should not fail");
        assert_eq!(result, "ACGTACGT");
    }

    #[test]
    fn test_umi_for_read_identity_lowercase_gets_uppercased() {
        let assigner = Strategy::Identity.new_assigner_full(0, 1, 100);
        let result = umi_for_read("acgtacgt", true, assigner.as_ref())
            .expect("valid lowercase UMI should not fail");
        assert_eq!(result, "ACGTACGT");
    }

    #[test]
    fn test_umi_for_read_paired_r1_earlier() {
        let assigner = Strategy::Paired.new_assigner_full(1, 1, 100);
        // With max_mismatches=1, prefix_len=2, so lower_prefix="aa", higher_prefix="bb"
        let result = umi_for_read("ACGT-TGCA", true, assigner.as_ref())
            .expect("valid paired UMI should not fail");
        // is_r1_earlier=true => lower_prefix:parts[0]-higher_prefix:parts[1]
        assert_eq!(result, "aa:ACGT-bb:TGCA");
    }

    #[test]
    fn test_umi_for_read_paired_r2_earlier() {
        let assigner = Strategy::Paired.new_assigner_full(1, 1, 100);
        // With max_mismatches=1, prefix_len=2, so lower_prefix="aa", higher_prefix="bb"
        let result = umi_for_read("ACGT-TGCA", false, assigner.as_ref())
            .expect("valid paired UMI should not fail");
        // is_r1_earlier=false => higher_prefix:parts[0]-lower_prefix:parts[1]
        assert_eq!(result, "bb:ACGT-aa:TGCA");
    }

    #[test]
    fn test_umi_for_read_paired_missing_dash_error() {
        let assigner = Strategy::Paired.new_assigner_full(1, 1, 100);
        let result = umi_for_read("ACGTACGT", true, assigner.as_ref());
        assert!(result.is_err());
        assert!(
            result
                .expect_err("should fail for missing dash")
                .to_string()
                .contains("did not contain 2 segments"),
            "Error message should mention missing segments"
        );
    }

    // ========================================================================
    // get_pair_orientation tests
    // ========================================================================

    #[test]
    fn test_get_pair_orientation_both_forward() {
        // Both R1 and R2 on forward strand (no REVERSE flag)
        let r1 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(30)
                .flags(flags::PAIRED | flags::FIRST_SEGMENT);
            b.add_string_tag(SamTag::RX, b"ACGT-TGCA");
            b.build()
        };
        let r2 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"TGCA")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(199)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(30)
                .flags(flags::PAIRED | flags::LAST_SEGMENT);
            b.add_string_tag(SamTag::RX, b"ACGT-TGCA");
            b.build()
        };
        let template = Template::from_records(vec![r1, r2])
            .expect("test template construction should not fail");
        let (r1_positive, r2_positive) = get_pair_orientation(&template);
        assert!(r1_positive);
        assert!(r2_positive);
    }

    #[test]
    fn test_get_pair_orientation_r1_reverse() {
        // R1 on reverse strand, R2 on forward strand
        let r1 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(30)
                .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::REVERSE);
            b.add_string_tag(SamTag::RX, b"ACGT-TGCA");
            b.build()
        };
        let r2 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"TGCA")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(199)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(30)
                .flags(flags::PAIRED | flags::LAST_SEGMENT);
            b.add_string_tag(SamTag::RX, b"ACGT-TGCA");
            b.build()
        };
        let template = Template::from_records(vec![r1, r2])
            .expect("test template construction should not fail");
        let (r1_positive, r2_positive) = get_pair_orientation(&template);
        assert!(!r1_positive);
        assert!(r2_positive);
    }

    #[test]
    fn test_get_pair_orientation_both_reverse() {
        // Both R1 and R2 on reverse strand
        let r1 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(30)
                .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::REVERSE);
            b.add_string_tag(SamTag::RX, b"ACGT-TGCA");
            b.build()
        };
        let r2 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"TGCA")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(199)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(30)
                .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE);
            b.add_string_tag(SamTag::RX, b"ACGT-TGCA");
            b.build()
        };
        let template = Template::from_records(vec![r1, r2])
            .expect("test template construction should not fail");
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
        let r1 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT);
            b.build()
        };
        let r2 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"TGCA")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(199)
                .cigar_ops(&[encode_op(0, 4)])
                .flags(flags::PAIRED | flags::LAST_SEGMENT);
            b.build()
        };
        assert!(is_r1_genomically_earlier_raw(&r1, &r2));
    }

    #[test]
    fn test_is_r1_earlier_false() {
        // R1 at position 200, R2 at position 100
        let r1 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(199)
                .cigar_ops(&[encode_op(0, 4)])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT);
            b.build()
        };
        let r2 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"TGCA")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)])
                .flags(flags::PAIRED | flags::LAST_SEGMENT);
            b.build()
        };
        assert!(!is_r1_genomically_earlier_raw(&r1, &r2));
    }

    #[test]
    fn test_is_r1_earlier_equal_position() {
        // Both mates share an unclipped 5' coordinate, so the tie breaks on
        // strand (mirroring fgbio): R1 is "earlier" iff it is on the forward
        // strand. A forward-strand R1 => true; a reverse-strand R1 => false.
        // Assert both halves so the strand-dependent tie-break stays pinned.
        let build_mate = |first_segment: bool, r1_reverse: bool| {
            let mut b = RawSamBuilder::new();
            let strand = if r1_reverse { flags::REVERSE } else { 0 };
            let segment = if first_segment { flags::FIRST_SEGMENT } else { flags::LAST_SEGMENT };
            b.read_name(b"q1")
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)])
                .flags(flags::PAIRED | segment | strand);
            b.build()
        };

        // Forward-strand R1 at the tie => R1 is earlier.
        let r1_fwd = build_mate(true, false);
        let r2_fwd = build_mate(false, false);
        assert!(is_r1_genomically_earlier_raw(&r1_fwd, &r2_fwd));

        // Reverse-strand R1 at the same tie => R1 is NOT earlier. This is the
        // half the old `r1_pos <= r2_pos` tie-break got wrong.
        let r1_rev = build_mate(true, true);
        let r2_rev = build_mate(false, false);
        assert!(!is_r1_genomically_earlier_raw(&r1_rev, &r2_rev));
    }

    // ========================================================================
    // truncate_umis tests
    // ========================================================================

    #[test]
    fn test_truncate_umis_no_min_length() {
        let umis = vec!["ACGTACGT".to_string(), "TGCATGCA".to_string()];
        let result = truncate_umis(umis.clone(), None)
            .expect("truncation with no min length should not fail");
        assert_eq!(result, umis);
    }

    #[test]
    fn test_truncate_umis_truncates() {
        let umis = vec!["ACGTACGT".to_string(), "TGCATGCA".to_string()];
        let result =
            truncate_umis(umis, Some(4)).expect("truncation of 8-base UMIs to 4 should not fail");
        assert_eq!(result, vec!["ACGT".to_string(), "TGCA".to_string()]);
    }

    #[test]
    fn test_truncate_umis_error_too_short() {
        let umis = vec!["ACG".to_string(), "TGCATGCA".to_string()];
        let result = truncate_umis(umis, Some(4));
        assert!(result.is_err());
        assert!(
            result
                .expect_err("should fail for too-short UMI")
                .to_string()
                .contains("shorter than expected"),
            "Error message should mention UMI being too short"
        );
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
        let r1 = {
            let mut b = RawSamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(r1_mapq)
                .flags(flags::PAIRED | flags::FIRST_SEGMENT);
            b.add_string_tag(SamTag::RX, umi.as_bytes());
            b.build()
        };
        let r2 = {
            let mut b = RawSamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"TGCA")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(199)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(r2_mapq)
                .flags(flags::PAIRED | flags::LAST_SEGMENT);
            b.add_string_tag(SamTag::RX, umi.as_bytes());
            b.build()
        };
        Template::from_records(vec![r1, r2]).expect("test template construction should not fail")
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
        let r1 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(30)
                .flags(flags::PAIRED | flags::FIRST_SEGMENT);
            b.add_string_tag(SamTag::RX, b"ACGTACGT");
            b.build()
        };
        let r2 = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"TGCA")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(199)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(30)
                .flags(flags::PAIRED | flags::LAST_SEGMENT);
            b.add_string_tag(SamTag::RX, b"ACNTACGT");
            b.build()
        };
        let template = Template::from_records(vec![r1, r2])
            .expect("test template construction should not fail");

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
            missing_tc_tag: 1,
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
            missing_tc_tag: 2,
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
        assert_eq!(m1.missing_tc_tag, 3);
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
        assert_eq!(metrics.missing_tc_tag, 0);
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
            missing_tc_tag: 2,
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
        assert_eq!(output.missing_tc_tag, 2);
    }

    /// DXM3-04: `duplicate_rate` serializes through the fgbio-style float
    /// formatter (#504's now-public `fgumi_metrics::float`), so a whole-number
    /// rate is written as `1`, not serde's raw `1.0` — consistent with the other
    /// fgumi metric files.
    #[test]
    fn test_dedup_metrics_duplicate_rate_uses_float_formatter() -> Result<()> {
        let metrics = DedupMetrics {
            total_templates: 2,
            duplicate_templates: 2,
            unique_templates: 0,
            ..Default::default()
        };
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("dedup.metrics.txt");
        write_dedup_metrics(&metrics, &path)?;

        let text = std::fs::read_to_string(&path)?;
        let data_row = text.lines().nth(1).expect("metrics data row");
        let columns: Vec<&str> = data_row.split('\t').collect();
        // Field order: total, unique, duplicate templates, then duplicate_rate.
        assert_eq!(
            columns.get(3).copied(),
            Some("1"),
            "duplicate_rate must format as `1`, not `1.0`; row: {data_row:?}"
        );
        Ok(())
    }

    // ========================================================================
    // Bounds check panic tests for raw-byte pipeline
    // ========================================================================

    /// Helper to create a minimal raw BAM record bytes for testing.
    /// The record has mapq=30 so it passes the Phase 1 MAPQ check.
    fn make_raw_bam_record_truncated_aux() -> RawRecord {
        // Create a BAM record where aux_data_offset will be beyond the record length
        let mut rec = vec![0u8; 40]; // Small record

        // Set l_seq to a large value that will make aux_data_offset too big
        // aux_data_offset = 32 + l_read_name + n_cigar_op*4 + l_seq.div_ceil(2) + l_seq
        // With l_seq=1000: offset = 32 + 4 + 0 + 500 + 1000 = 1536 >> 40 (record length)
        rec[8] = 4; // l_read_name = 4 (3 chars + null)
        rec[9] = 30; // mapq = 30 (must be >= min_mapq to reach Phase 2 tag checks)
        rec[12..14].copy_from_slice(&0u16.to_le_bytes()); // n_cigar_op = 0
        rec[16..20].copy_from_slice(&1000u32.to_le_bytes()); // l_seq = 1000 (too large)
        rec[32..36].copy_from_slice(b"tst\0"); // read name

        RawRecord::from(rec)
    }

    #[test]
    fn test_filter_template_raw_truncated_record_no_panic() {
        // After fix: truncated records no longer panic -- they gracefully return false
        // (rejected as poor alignment due to missing UMI tag).

        let raw = make_raw_bam_record_truncated_aux();
        let template =
            Template::from_records(vec![raw]).expect("test template construction should not fail");

        let config = DedupFilterConfig {
            umi_tag: *SamTag::RX,
            min_mapq: 20,
            include_non_pf: false,
            min_umi_length: None,
            no_umi: false,
        };
        let mut metrics = FilterMetrics::new();

        // Should not panic; should reject gracefully due to missing UMI
        let result = filter_template(&template, &config, &mut metrics);
        assert!(!result, "Truncated record should be rejected");
        assert_eq!(metrics.discarded_poor_alignment, 1);
    }

    // Test that the raw-byte mode MQ tag issue actually occurs in filter_template_raw
    #[test]
    fn test_filter_template_raw_mq_tag_type_mismatch_bypasses_filter() {
        // Build a BAM record with exact sizes so aux data is contiguous.
        // Header(32) + name(4) + cigar(4) + seq(2) + qual(4) = 46 bytes before aux
        let mut rec = vec![0u8; 46];

        rec[0..4].copy_from_slice(&0i32.to_le_bytes()); // ref_id = 0 (mapped)
        rec[4..8].copy_from_slice(&100i32.to_le_bytes()); // pos = 100
        rec[8] = 4; // l_read_name = 4
        rec[9] = 30; // mapq = 30 (passes the direct MAPQ check)
        rec[12..14].copy_from_slice(&1u16.to_le_bytes()); // n_cigar_op = 1
        // flags: PAIRED but NOT MATE_UNMAPPED (so mate MQ check triggers)
        rec[14..16].copy_from_slice(&(fgumi_raw_bam::flags::PAIRED).to_le_bytes());
        rec[16..20].copy_from_slice(&4u32.to_le_bytes()); // l_seq = 4
        rec[32..36].copy_from_slice(b"tst\0"); // read name
        rec[36..40].copy_from_slice(&(4u32 << 4).to_le_bytes()); // CIGAR: 4M

        // Aux data starts immediately at byte 46.
        // Add RX tag first (needed to pass UMI check)
        rec.extend_from_slice(b"RXZACGTACGT\x00"); // RX:Z:ACGTACGT

        // Add MQ tag as signed byte with low value (should trigger filtering)
        rec.extend_from_slice(b"MQc"); // tag=MQ, type=c (signed byte)
        rec.push(10); // MAPQ = 10 (< min_mapq of 20)

        let template = Template::from_records(vec![RawRecord::from(rec)])
            .expect("test template construction should not fail");

        let config = DedupFilterConfig {
            umi_tag: *SamTag::RX,
            min_mapq: 20,
            include_non_pf: false,
            min_umi_length: None,
            no_umi: false,
        };
        let mut metrics = FilterMetrics::new();

        // After fix: filter_template_raw now uses find_int_tag which handles all integer types.
        // The template should be REJECTED because mate MAPQ=10 < min_mapq=20.
        let result = filter_template(&template, &config, &mut metrics);

        assert!(!result, "Template with low mate MAPQ should be rejected");
        assert_eq!(metrics.accepted_templates, 0);
        assert_eq!(metrics.discarded_poor_alignment, 1);
    }

    // ========================================================================
    // Raw BAM record builder helper for dedup tests
    // ========================================================================

    /// Build a raw BAM record with specified qualities and UMI for testing.
    fn make_raw_bam_for_dedup(
        name: &[u8],
        flag: u16,
        mapq: u8,
        seq_len: usize,
        qualities: &[u8],
        umi: &[u8],
    ) -> RawRecord {
        use fgumi_raw_bam;

        let l_read_name = (name.len() + 1) as u8;
        let cigar_ops: &[u32] = if (flag & fgumi_raw_bam::flags::UNMAPPED) == 0 {
            &[(seq_len as u32) << 4] // NM cigar
        } else {
            &[]
        };
        let n_cigar_op = cigar_ops.len() as u16;
        let seq_bytes = seq_len.div_ceil(2);
        let total = 32 + l_read_name as usize + cigar_ops.len() * 4 + seq_bytes + seq_len;
        let mut buf = vec![0u8; total];

        buf[0..4].copy_from_slice(&0i32.to_le_bytes()); // tid
        buf[4..8].copy_from_slice(&100i32.to_le_bytes()); // pos
        buf[8] = l_read_name;
        buf[9] = mapq;
        buf[12..14].copy_from_slice(&n_cigar_op.to_le_bytes());
        buf[14..16].copy_from_slice(&flag.to_le_bytes());
        buf[16..20].copy_from_slice(&(seq_len as u32).to_le_bytes());
        buf[20..24].copy_from_slice(&(-1i32).to_le_bytes()); // mate_tid
        buf[24..28].copy_from_slice(&(-1i32).to_le_bytes()); // mate_pos

        let name_start = 32;
        buf[name_start..name_start + name.len()].copy_from_slice(name);
        buf[name_start + name.len()] = 0;

        let cigar_start = name_start + l_read_name as usize;
        for (i, &op) in cigar_ops.iter().enumerate() {
            let off = cigar_start + i * 4;
            buf[off..off + 4].copy_from_slice(&op.to_le_bytes());
        }

        // Write quality scores
        let qo = fgumi_raw_bam::qual_offset(&buf);
        let copy_len = std::cmp::min(qualities.len(), seq_len);
        buf[qo..qo + copy_len].copy_from_slice(&qualities[..copy_len]);

        // Append UMI tag: RX:Z:<umi>\0
        buf.extend_from_slice(b"RXZ");
        buf.extend_from_slice(umi);
        buf.push(0);

        RawRecord::from(buf)
    }

    /// Picard `computeDuplicateScore` caps each read's Q>=15 base-quality sum at
    /// `Short.MAX_VALUE / 2` (16383) before the two mates are summed. A single read
    /// whose sum exceeds the cap must contribute exactly the cap.
    #[test]
    fn test_score_template_picard_per_read_cap() {
        // 500 bases at Q40 -> raw sum 20000, capped to 16383.
        let raw = make_raw_bam_for_dedup(
            b"cap",
            fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::FIRST_SEGMENT,
            60,
            500,
            &[40u8; 500],
            b"ACGT",
        );
        let template =
            Template::from_records(vec![raw]).expect("test template construction should not fail");
        assert_eq!(score_template(&template), PICARD_MAX_SCORE_PER_READ);
    }

    /// Picard `computeDuplicateScore` subtracts `Short.MIN_VALUE / 2` (16384) from a
    /// QC-fail read so it is never preferred as the duplicate representative
    /// (reachable in fgumi only under `--include-non-pf`). An otherwise-identical
    /// QC-fail read must score strictly below the passing read.
    #[test]
    fn test_score_template_qc_fail_discount() {
        let base_flags = fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::FIRST_SEGMENT;
        let passing = make_raw_bam_for_dedup(b"pass", base_flags, 60, 4, &[40u8; 4], b"ACGT");
        let qc_fail = make_raw_bam_for_dedup(
            b"nonpf",
            base_flags | fgumi_raw_bam::flags::QC_FAIL,
            60,
            4,
            &[40u8; 4],
            b"ACGT",
        );
        let pf_t = Template::from_records(vec![passing]).expect("template construction");
        let qc_t = Template::from_records(vec![qc_fail]).expect("template construction");
        // 4 * 40 = 160 for the passing read; QC-fail is discounted by 16384.
        assert_eq!(score_template(&pf_t), 160);
        assert_eq!(score_template(&qc_t), 160 + PICARD_QC_FAIL_DISCOUNT);
        assert!(score_template(&qc_t) < score_template(&pf_t), "QC-fail must score lower");
    }

    // ========================================================================
    // score_template_raw tests
    // ========================================================================

    #[test]
    fn test_score_template_raw_basic() {
        // 4 bases at quality 20 (all >= 15) → 4 * 20 = 80
        let raw = make_raw_bam_for_dedup(
            b"rea",
            fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::FIRST_SEGMENT,
            30,
            4,
            &[20, 20, 20, 20],
            b"ACGT",
        );
        let template =
            Template::from_records(vec![raw]).expect("test template construction should not fail");
        assert_eq!(score_template(&template), 80);
    }

    #[test]
    fn test_score_template_raw_low_quality() {
        // 4 bases at quality 10 (all < 15) → excluded → 0
        let raw = make_raw_bam_for_dedup(
            b"rea",
            fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::FIRST_SEGMENT,
            30,
            4,
            &[10, 10, 10, 10],
            b"ACGT",
        );
        let template =
            Template::from_records(vec![raw]).expect("test template construction should not fail");
        assert_eq!(score_template(&template), 0);
    }

    #[test]
    fn test_score_template_raw_paired() {
        // Paired: R1 (4 bases * 20, all >= 15) + R2 (4 bases * 10, all < 15) = 80 + 0 = 80
        let r1 = make_raw_bam_for_dedup(
            b"rea",
            fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::FIRST_SEGMENT,
            30,
            4,
            &[20, 20, 20, 20],
            b"ACGT",
        );
        let r2 = make_raw_bam_for_dedup(
            b"rea",
            fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::LAST_SEGMENT,
            30,
            4,
            &[10, 10, 10, 10],
            b"ACGT",
        );
        let template = Template::from_records(vec![r1, r2])
            .expect("test template construction should not fail");
        assert_eq!(score_template(&template), 80);
    }

    #[test]
    fn test_score_template_raw_matches_recordbuf() {
        // Verify raw scoring matches RecordBuf scoring
        let qualities = &[20, 15, 10, 5];
        let template_rb = create_test_template("q1", qualities);
        let template_raw = {
            let raw = make_raw_bam_for_dedup(
                b"q01",
                fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::FIRST_SEGMENT,
                30,
                4,
                qualities,
                b"ACGT",
            );
            Template::from_records(vec![raw]).expect("test template construction should not fail")
        };
        assert_eq!(score_template(&template_rb), score_template(&template_raw));
    }

    // ========================================================================
    // filter_template_raw passing/rejecting tests
    // ========================================================================

    /// Which `FilterMetrics` field an `rstest` case expects to advance.
    #[derive(Debug, Clone, Copy)]
    enum FilterExpect {
        Accepted,
        PoorAlignment,
        NonPf,
        NsInUmi,
        UmiTooShort,
    }

    #[rstest]
    #[case::accepts_valid(
        fgumi_raw_bam::flags::PAIRED
            | fgumi_raw_bam::flags::FIRST_SEGMENT
            | fgumi_raw_bam::flags::MATE_UNMAPPED,
        30, b"ACGTACGT", 20, None, true, FilterExpect::Accepted
    )]
    #[case::rejects_low_mapq(
        fgumi_raw_bam::flags::PAIRED
            | fgumi_raw_bam::flags::FIRST_SEGMENT
            | fgumi_raw_bam::flags::MATE_UNMAPPED,
        10, b"ACGT", 20, None, false, FilterExpect::PoorAlignment
    )]
    #[case::rejects_qc_fail(
        fgumi_raw_bam::flags::PAIRED
            | fgumi_raw_bam::flags::FIRST_SEGMENT
            | fgumi_raw_bam::flags::MATE_UNMAPPED
            | fgumi_raw_bam::flags::QC_FAIL,
        30, b"ACGT", 0, None, false, FilterExpect::NonPf
    )]
    #[case::rejects_umi_with_n(
        fgumi_raw_bam::flags::PAIRED
            | fgumi_raw_bam::flags::FIRST_SEGMENT
            | fgumi_raw_bam::flags::MATE_UNMAPPED,
        30, b"ANGT", 0, None, false, FilterExpect::NsInUmi
    )]
    #[case::rejects_short_umi(
        fgumi_raw_bam::flags::PAIRED
            | fgumi_raw_bam::flags::FIRST_SEGMENT
            | fgumi_raw_bam::flags::MATE_UNMAPPED,
        30, b"AC", 0, Some(6), false, FilterExpect::UmiTooShort
    )]
    #[case::rejects_unmapped(
        fgumi_raw_bam::flags::UNMAPPED
            | fgumi_raw_bam::flags::MATE_UNMAPPED,
        0, b"ACGT", 0, None, false, FilterExpect::PoorAlignment
    )]
    fn test_filter_template_raw(
        #[case] flags: u16,
        #[case] mapq: u8,
        #[case] umi: &'static [u8],
        #[case] min_mapq: u8,
        #[case] min_umi_length: Option<usize>,
        #[case] expect_pass: bool,
        #[case] expect: FilterExpect,
    ) {
        let raw = make_raw_bam_for_dedup(b"rea", flags, mapq, 4, &[20, 20, 20, 20], umi);
        let template =
            Template::from_records(vec![raw]).expect("test template construction should not fail");
        let config = DedupFilterConfig {
            umi_tag: *SamTag::RX,
            min_mapq,
            include_non_pf: false,
            min_umi_length,
            no_umi: false,
        };
        let mut metrics = FilterMetrics::new();
        assert_eq!(filter_template(&template, &config, &mut metrics), expect_pass);
        match expect {
            FilterExpect::Accepted => assert_eq!(metrics.accepted_templates, 1),
            FilterExpect::PoorAlignment => assert_eq!(metrics.discarded_poor_alignment, 1),
            FilterExpect::NonPf => assert_eq!(metrics.discarded_non_pf, 1),
            FilterExpect::NsInUmi => assert_eq!(metrics.discarded_ns_in_umi, 1),
            FilterExpect::UmiTooShort => assert_eq!(metrics.discarded_umi_too_short, 1),
        }
    }

    // ========================================================================
    // Duplicate flag manipulation on raw bytes
    // ========================================================================

    #[test]
    fn test_set_duplicate_flag_on_raw_record() {
        use fgumi_raw_bam;
        use fgumi_raw_bam::RawRecordView;

        let raw =
            make_raw_bam_for_dedup(b"rea", fgumi_raw_bam::flags::PAIRED, 30, 4, &[20; 4], b"ACGT");
        let mut rec = raw;
        let orig_flags = RawRecordView::new(&rec).flags();
        assert_eq!(orig_flags & fgumi_raw_bam::flags::DUPLICATE, 0);

        // Set duplicate flag
        fgumi_raw_bam::set_flags(&mut rec, orig_flags | fgumi_raw_bam::flags::DUPLICATE);
        assert_ne!(RawRecordView::new(&rec).flags() & fgumi_raw_bam::flags::DUPLICATE, 0);

        // Clear duplicate flag
        let flags_with_dup = RawRecordView::new(&rec).flags();
        fgumi_raw_bam::set_flags(&mut rec, flags_with_dup & !fgumi_raw_bam::flags::DUPLICATE);
        assert_eq!(RawRecordView::new(&rec).flags() & fgumi_raw_bam::flags::DUPLICATE, 0);
    }

    #[test]
    fn test_processed_dedup_group_memory_estimate() {
        let template = create_test_template("test", &[30, 30, 30, 30]);
        let mut templates = Vec::with_capacity(10);
        templates.push(template);

        let batch = ProcessedDedupGroup {
            templates,
            family_sizes: AHashMap::new(),
            dedup_metrics: DedupMetrics::default(),
            input_record_count: 1,
            distinct_mi_count: 0,
        };

        let estimate = batch.estimate_heap_size();
        // Should include Vec<Template> overhead for capacity * size_of::<Template>()
        let vec_overhead = 10 * std::mem::size_of::<Template>();
        assert!(
            estimate >= vec_overhead,
            "estimate {estimate} should include Vec<Template> overhead {vec_overhead}"
        );
    }

    // ========================================================================
    // no_umi mode tests
    // ========================================================================

    #[test]
    fn test_filter_template_accepts_missing_umi_when_no_umi_mode() {
        let mut config = default_filter_config();
        config.no_umi = true;
        let mut metrics = FilterMetrics::new();

        // Create template without RX tag
        let record = {
            let mut b = RawSamBuilder::new();
            b.read_name(b"q1")
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30])
                .ref_id(0)
                .pos(99)
                .cigar_ops(&[encode_op(0, 4)])
                .mapq(30)
                .flags(flags::PAIRED | flags::FIRST_SEGMENT);
            b.build()
        };
        let template = Template::from_records(vec![record])
            .expect("test template construction should not fail");

        // In no_umi mode, templates without UMI tags should be accepted
        assert!(filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.accepted_templates, 1);
    }

    #[test]
    fn test_filter_template_accepts_umi_with_n_when_no_umi_mode() {
        let mut config = default_filter_config();
        config.no_umi = true;
        let mut metrics = FilterMetrics::new();
        let template = create_mapped_template_with_umi("q1", "ACNTACGT", 30);

        // In no_umi mode, UMIs with N should be accepted (UMI validation is skipped)
        assert!(filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.accepted_templates, 1);
    }

    #[test]
    fn test_filter_template_accepts_short_umi_when_no_umi_mode() {
        let mut config = default_filter_config();
        config.min_umi_length = Some(8);
        config.no_umi = true;
        let mut metrics = FilterMetrics::new();
        let template = create_mapped_template_with_umi("q1", "ACGT", 30);

        // In no_umi mode, short UMIs should be accepted (min_umi_length is not checked)
        assert!(filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.accepted_templates, 1);
    }

    // ========================================================================
    // no_umi mode tests (raw path)
    // ========================================================================

    /// Creates a raw BAM record without any auxiliary tags (no UMI tag).
    fn make_raw_bam_for_dedup_no_tags(
        name: &[u8],
        flag: u16,
        mapq: u8,
        seq_len: usize,
        qualities: &[u8],
    ) -> RawRecord {
        use fgumi_raw_bam;

        let l_read_name = (name.len() + 1) as u8;
        let cigar_ops: &[u32] = if (flag & fgumi_raw_bam::flags::UNMAPPED) == 0 {
            &[(seq_len as u32) << 4]
        } else {
            &[]
        };
        let n_cigar_op = cigar_ops.len() as u16;
        let seq_bytes = seq_len.div_ceil(2);
        let total = 32 + l_read_name as usize + cigar_ops.len() * 4 + seq_bytes + seq_len;
        let mut buf = vec![0u8; total];

        buf[0..4].copy_from_slice(&0i32.to_le_bytes());
        buf[4..8].copy_from_slice(&100i32.to_le_bytes());
        buf[8] = l_read_name;
        buf[9] = mapq;
        buf[12..14].copy_from_slice(&n_cigar_op.to_le_bytes());
        buf[14..16].copy_from_slice(&flag.to_le_bytes());
        buf[16..20].copy_from_slice(&(seq_len as u32).to_le_bytes());
        buf[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        buf[24..28].copy_from_slice(&(-1i32).to_le_bytes());

        let name_start = 32;
        buf[name_start..name_start + name.len()].copy_from_slice(name);
        buf[name_start + name.len()] = 0;

        let cigar_start = name_start + l_read_name as usize;
        for (i, &op) in cigar_ops.iter().enumerate() {
            let off = cigar_start + i * 4;
            buf[off..off + 4].copy_from_slice(&op.to_le_bytes());
        }

        let qo = fgumi_raw_bam::qual_offset(&buf);
        let copy_len = std::cmp::min(qualities.len(), seq_len);
        buf[qo..qo + copy_len].copy_from_slice(&qualities[..copy_len]);

        RawRecord::from(buf)
    }

    #[test]
    fn test_filter_template_raw_accepts_missing_umi_when_no_umi_mode() {
        let raw = make_raw_bam_for_dedup_no_tags(
            b"rea",
            fgumi_raw_bam::flags::PAIRED
                | fgumi_raw_bam::flags::FIRST_SEGMENT
                | fgumi_raw_bam::flags::MATE_UNMAPPED,
            30,
            4,
            &[20, 20, 20, 20],
        );
        let template =
            Template::from_records(vec![raw]).expect("test template construction should not fail");
        let config = DedupFilterConfig {
            umi_tag: *SamTag::RX,
            min_mapq: 20,
            include_non_pf: false,
            min_umi_length: None,
            no_umi: true,
        };
        let mut metrics = FilterMetrics::new();

        // In no_umi mode, templates without UMI tags should be accepted
        assert!(filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.accepted_templates, 1);
    }

    #[test]
    fn test_filter_template_raw_accepts_umi_with_n_when_no_umi_mode() {
        let raw = make_raw_bam_for_dedup(
            b"rea",
            fgumi_raw_bam::flags::PAIRED
                | fgumi_raw_bam::flags::FIRST_SEGMENT
                | fgumi_raw_bam::flags::MATE_UNMAPPED,
            30,
            4,
            &[20, 20, 20, 20],
            b"ANGT",
        );
        let template =
            Template::from_records(vec![raw]).expect("test template construction should not fail");
        let config = DedupFilterConfig {
            umi_tag: *SamTag::RX,
            min_mapq: 0,
            include_non_pf: false,
            min_umi_length: None,
            no_umi: true,
        };
        let mut metrics = FilterMetrics::new();

        // In no_umi mode, UMIs with N should be accepted (UMI validation is skipped)
        assert!(filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.accepted_templates, 1);
    }

    #[test]
    fn test_filter_template_raw_accepts_short_umi_when_no_umi_mode() {
        let raw = make_raw_bam_for_dedup(
            b"rea",
            fgumi_raw_bam::flags::PAIRED
                | fgumi_raw_bam::flags::FIRST_SEGMENT
                | fgumi_raw_bam::flags::MATE_UNMAPPED,
            30,
            4,
            &[20, 20, 20, 20],
            b"AC",
        );
        let template =
            Template::from_records(vec![raw]).expect("test template construction should not fail");
        let config = DedupFilterConfig {
            umi_tag: *SamTag::RX,
            min_mapq: 0,
            include_non_pf: false,
            min_umi_length: Some(6),
            no_umi: true,
        };
        let mut metrics = FilterMetrics::new();

        // In no_umi mode, short UMIs should be accepted (min_umi_length is not checked)
        assert!(filter_template(&template, &config, &mut metrics));
        assert_eq!(metrics.accepted_templates, 1);
    }
}
