//! Consensus read filtering logic.
//!
//! This module provides functionality for filtering consensus reads based on quality, depth,
//! and error rate thresholds. It supports both single-strand and duplex consensus reads.

use ahash::AHashMap;
use anyhow::Result;
use noodles::sam::alignment::record::Sequence as SequenceTrait;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

use crate::tags::{per_base, per_read};
use crate::phred::{MIN_PHRED, NO_CALL_BASE};
use fgumi_metrics::rejection::RejectionReason;

/// Filter thresholds for consensus reads
#[derive(Debug, Clone)]
pub struct FilterThresholds {
    /// Minimum number of raw reads to support a consensus base/read
    pub min_reads: usize,

    /// Maximum raw read error rate (0.0-1.0)
    pub max_read_error_rate: f64,

    /// Maximum base error rate (0.0-1.0)
    pub max_base_error_rate: f64,
}

/// Consensus read type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConsensusType {
    /// Single-strand consensus (from `simplex`)
    SingleStrand,

    /// Duplex consensus (from `duplex`)
    Duplex,
}

/// Filtering configuration for consensus reads
#[derive(Debug, Clone)]
pub struct FilterConfig {
    /// Thresholds for duplex consensus reads (if applicable)
    pub duplex_thresholds: Option<FilterThresholds>,

    /// Thresholds for AB/top-strand consensus
    pub ab_thresholds: Option<FilterThresholds>,

    /// Thresholds for BA/bottom-strand consensus
    pub ba_thresholds: Option<FilterThresholds>,

    /// Thresholds for single-strand consensus
    pub single_strand_thresholds: Option<FilterThresholds>,

    /// Minimum base quality after masking (optional - None means no quality masking)
    pub min_base_quality: Option<u8>,

    /// Minimum mean base quality after masking (optional)
    pub min_mean_base_quality: Option<f64>,

    /// Maximum fraction of no-calls (N bases) allowed (0.0-1.0)
    pub max_no_call_fraction: f64,
}

impl FilterConfig {
    /// Creates a filter configuration for single-strand (simplex) consensus reads.
    ///
    /// This is the simplest configuration - one set of thresholds applied to all reads.
    ///
    /// # Arguments
    /// * `thresholds` - Filter thresholds for single-strand consensus
    /// * `min_base_quality` - Minimum base quality after masking (None means no quality masking)
    /// * `min_mean_base_quality` - Optional minimum mean base quality
    /// * `max_no_call_fraction` - Maximum fraction of N bases allowed
    #[must_use]
    pub fn for_single_strand(
        thresholds: FilterThresholds,
        min_base_quality: Option<u8>,
        min_mean_base_quality: Option<f64>,
        max_no_call_fraction: f64,
    ) -> Self {
        Self {
            duplex_thresholds: Some(thresholds.clone()),
            ab_thresholds: Some(thresholds.clone()),
            ba_thresholds: Some(thresholds.clone()),
            single_strand_thresholds: Some(thresholds),
            min_base_quality,
            min_mean_base_quality,
            max_no_call_fraction,
        }
    }

    /// Creates a filter configuration for symmetric duplex consensus reads.
    ///
    /// Uses the same thresholds for both AB and BA strands.
    ///
    /// # Arguments
    /// * `duplex` - Filter thresholds for final duplex consensus
    /// * `strand` - Filter thresholds for both AB and BA strands (symmetric)
    /// * `min_base_quality` - Minimum base quality after masking (None means no quality masking)
    /// * `min_mean_base_quality` - Optional minimum mean base quality
    /// * `max_no_call_fraction` - Maximum fraction of N bases allowed
    ///
    /// # Panics
    ///
    /// Panics if thresholds violate ordering constraints (`strand.min_reads` <= `duplex.min_reads`, etc.)
    #[must_use]
    pub fn for_duplex(
        duplex: FilterThresholds,
        strand: FilterThresholds,
        min_base_quality: Option<u8>,
        min_mean_base_quality: Option<f64>,
        max_no_call_fraction: f64,
    ) -> Self {
        Self::for_duplex_asymmetric(
            duplex,
            strand.clone(),
            strand,
            min_base_quality,
            min_mean_base_quality,
            max_no_call_fraction,
        )
    }

    /// Creates a filter configuration for asymmetric duplex consensus reads.
    ///
    /// Uses different thresholds for AB (higher depth strand) and BA (lower depth strand).
    ///
    /// # Arguments
    /// * `duplex` - Filter thresholds for final duplex consensus
    /// * `ab` - Filter thresholds for AB strand (typically higher depth)
    /// * `ba` - Filter thresholds for BA strand (typically lower depth)
    /// * `min_base_quality` - Minimum base quality after masking (None means no quality masking)
    /// * `min_mean_base_quality` - Optional minimum mean base quality
    /// * `max_no_call_fraction` - Maximum fraction of N bases allowed
    ///
    /// # Panics
    ///
    /// Panics if thresholds violate ordering constraints:
    /// - `min_reads`: BA <= AB <= duplex
    /// - error rates: AB <= BA (AB more stringent)
    #[must_use]
    pub fn for_duplex_asymmetric(
        duplex: FilterThresholds,
        ab: FilterThresholds,
        ba: FilterThresholds,
        min_base_quality: Option<u8>,
        min_mean_base_quality: Option<f64>,
        max_no_call_fraction: f64,
    ) -> Self {
        // Validate threshold ordering (matching fgbio's validation)
        assert!(
            ab.min_reads <= duplex.min_reads,
            "min-reads values must be specified high to low: AB ({}) > duplex ({})",
            ab.min_reads,
            duplex.min_reads
        );
        assert!(
            ba.min_reads <= ab.min_reads,
            "min-reads values must be specified high to low: BA ({}) > AB ({})",
            ba.min_reads,
            ab.min_reads
        );
        assert!(
            ab.max_read_error_rate <= ba.max_read_error_rate,
            "max-read-error-rate for AB ({}) must be <= BA ({})",
            ab.max_read_error_rate,
            ba.max_read_error_rate
        );
        assert!(
            ab.max_base_error_rate <= ba.max_base_error_rate,
            "max-base-error-rate for AB ({}) must be <= BA ({})",
            ab.max_base_error_rate,
            ba.max_base_error_rate
        );

        Self {
            duplex_thresholds: Some(duplex.clone()),
            ab_thresholds: Some(ab),
            ba_thresholds: Some(ba),
            single_strand_thresholds: Some(duplex),
            min_base_quality,
            min_mean_base_quality,
            max_no_call_fraction,
        }
    }

    /// Returns duplex (CC), AB, and BA thresholds, or `None` if any are missing.
    #[must_use]
    pub fn duplex_thresholds(
        &self,
    ) -> Option<(&FilterThresholds, &FilterThresholds, &FilterThresholds)> {
        match (
            self.duplex_thresholds.as_ref(),
            self.ab_thresholds.as_ref(),
            self.ba_thresholds.as_ref(),
        ) {
            (Some(cc), Some(ab), Some(ba)) => Some((cc, ab, ba)),
            _ => None,
        }
    }

    /// Returns the single-strand thresholds, falling back to duplex thresholds.
    #[must_use]
    pub fn effective_single_strand_thresholds(&self) -> Option<&FilterThresholds> {
        self.single_strand_thresholds.as_ref().or(self.duplex_thresholds.as_ref())
    }

    /// Creates a new filter configuration from parameter vectors
    ///
    /// # Arguments
    /// * `min_reads` - 1-3 values for [duplex, AB, BA] or [single-strand]
    /// * `max_read_error_rate` - 1-3 values for [duplex, AB, BA] or [single-strand]
    /// * `max_base_error_rate` - 1-3 values for [duplex, AB, BA] or [single-strand]
    /// * `min_base_quality` - Minimum base quality after masking (None means no quality masking)
    /// * `min_mean_base_quality` - Optional minimum mean base quality
    /// * `max_no_call_fraction` - Maximum fraction of N bases allowed
    ///
    /// # Panics
    ///
    /// Panics if thresholds violate ordering constraints:
    /// - `min_reads`: BA <= AB <= duplex
    /// - error rates: AB <= BA (AB more stringent)
    #[must_use]
    pub fn new(
        min_reads: &[usize],
        max_read_error_rate: &[f64],
        max_base_error_rate: &[f64],
        min_base_quality: Option<u8>,
        min_mean_base_quality: Option<f64>,
        max_no_call_fraction: f64,
    ) -> Self {
        // Helper to get threshold at index with fallback
        let get_threshold = |values: &[usize], index: usize| -> usize {
            match values.len() {
                2 => {
                    if index == 0 {
                        values[0]
                    } else {
                        values[1]
                    }
                }
                3 => values[index],
                _ => values[0],
            }
        };

        let get_threshold_f64 = |values: &[f64], index: usize| -> f64 {
            match values.len() {
                2 => {
                    if index == 0 {
                        values[0]
                    } else {
                        values[1]
                    }
                }
                3 => values[index],
                _ => values[0],
            }
        };

        // Create thresholds for all levels - matching fgbio which always creates all three
        // when filtering either simplex or duplex reads. Single values are replicated to all levels.
        let duplex_thresholds = Some(FilterThresholds {
            min_reads: get_threshold(min_reads, 0),
            max_read_error_rate: get_threshold_f64(max_read_error_rate, 0),
            max_base_error_rate: get_threshold_f64(max_base_error_rate, 0),
        });

        let ab_thresholds = Some(FilterThresholds {
            min_reads: get_threshold(min_reads, 1),
            max_read_error_rate: get_threshold_f64(max_read_error_rate, 1),
            max_base_error_rate: get_threshold_f64(max_base_error_rate, 1),
        });

        let ba_thresholds = Some(FilterThresholds {
            min_reads: get_threshold(min_reads, 2),
            max_read_error_rate: get_threshold_f64(max_read_error_rate, 2),
            max_base_error_rate: get_threshold_f64(max_base_error_rate, 2),
        });

        // Also create single-strand thresholds using the first value
        let single_strand_thresholds = Some(FilterThresholds {
            min_reads: min_reads[0],
            max_read_error_rate: if max_read_error_rate.is_empty() {
                1.0
            } else {
                max_read_error_rate[0]
            },
            max_base_error_rate: if max_base_error_rate.is_empty() {
                1.0
            } else {
                max_base_error_rate[0]
            },
        });

        // Validate threshold ordering for duplex mode (matching fgbio's validation)
        // For depth thresholds: BA <= AB <= CC (values must be specified high to low)
        // For error rates: AB <= BA (AB must be more stringent than BA)
        if let (Some(cc), Some(ab), Some(ba)) = (&duplex_thresholds, &ab_thresholds, &ba_thresholds)
        {
            // min_reads: CC >= AB >= BA
            assert!(
                ab.min_reads <= cc.min_reads,
                "min-reads values must be specified high to low: AB ({}) > CC ({})",
                ab.min_reads,
                cc.min_reads
            );
            assert!(
                ba.min_reads <= ab.min_reads,
                "min-reads values must be specified high to low: BA ({}) > AB ({})",
                ba.min_reads,
                ab.min_reads
            );

            // max_read_error_rate: AB <= BA (AB more stringent)
            assert!(
                ab.max_read_error_rate <= ba.max_read_error_rate,
                "max-read-error-rate for AB ({}) must be <= BA ({})",
                ab.max_read_error_rate,
                ba.max_read_error_rate
            );

            // max_base_error_rate: AB <= BA (AB more stringent)
            assert!(
                ab.max_base_error_rate <= ba.max_base_error_rate,
                "max-base-error-rate for AB ({}) must be <= BA ({})",
                ab.max_base_error_rate,
                ba.max_base_error_rate
            );
        }

        Self {
            duplex_thresholds,
            ab_thresholds,
            ba_thresholds,
            single_strand_thresholds,
            min_base_quality,
            min_mean_base_quality,
            max_no_call_fraction,
        }
    }
}

/// Result of filtering a consensus read
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterResult {
    /// Read passed all filters
    Pass,

    /// Read failed minimum reads threshold
    InsufficientReads,

    /// Read failed maximum error rate threshold
    ExcessiveErrorRate,

    /// Read failed minimum mean base quality threshold
    LowQuality,

    /// Read failed maximum no-call fraction threshold
    TooManyNoCalls,
}

impl FilterResult {
    /// Converts a `FilterResult` to a `RejectionReason`, if rejected.
    ///
    /// Returns `None` if the result is `Pass`, otherwise returns the corresponding
    /// rejection reason for tracking in metrics and logging.
    #[must_use]
    pub fn to_rejection_reason(&self) -> Option<RejectionReason> {
        match self {
            Self::Pass => None,
            Self::InsufficientReads => Some(RejectionReason::InsufficientSupport),
            Self::ExcessiveErrorRate => Some(RejectionReason::ExcessiveErrorRate),
            Self::LowQuality => Some(RejectionReason::LowMeanQuality),
            Self::TooManyNoCalls => Some(RejectionReason::ExcessiveNBases),
        }
    }
}

/// Helper function to extract integer value from SAM tag, handling different integer types
fn extract_int_tag(data: &noodles::sam::alignment::record_buf::Data, tag_str: &str) -> Option<i32> {
    let tag = per_read::tag(tag_str);
    match data.get(&tag)? {
        Value::Int8(v) => Some(i32::from(*v)),
        Value::UInt8(v) => Some(i32::from(*v)),
        Value::Int16(v) => Some(i32::from(*v)),
        Value::UInt16(v) => Some(i32::from(*v)),
        Value::Int32(v) => Some(*v),
        #[expect(clippy::cast_possible_wrap, reason = "BAM tag values in practice never exceed i32::MAX")]
        Value::UInt32(v) => Some(*v as i32),
        _ => None,
    }
}

/// Helper function to extract float value from SAM tag
fn extract_float_tag(
    data: &noodles::sam::alignment::record_buf::Data,
    tag_str: &str,
) -> Option<f32> {
    let tag = per_read::tag(tag_str);
    match data.get(&tag)? {
        Value::Float(v) => Some(*v),
        _ => None,
    }
}

/// Filters a consensus read based on per-read tags
///
/// For duplex reads, also checks AB and BA strand tags independently.
///
/// # Arguments
/// * `record` - The consensus read to filter
/// * `thresholds` - The filter thresholds to apply
/// * `ab_thresholds` - Optional thresholds for AB strand (duplex only)
/// * `ba_thresholds` - Optional thresholds for BA strand (duplex only)
///
/// # Returns
/// `FilterResult` indicating whether the read passed or why it failed
///
/// # Errors
///
/// Returns an error if the record data cannot be read.
pub fn filter_read(record: &RecordBuf, thresholds: &FilterThresholds) -> Result<FilterResult> {
    // Check minimum reads (cD tag)
    let depth_tag = per_read::tag("cD");
    if let Some(Value::UInt8(depth)) = record.data().get(&depth_tag) {
        if (*depth as usize) < thresholds.min_reads {
            return Ok(FilterResult::InsufficientReads);
        }
    }

    // Check maximum error rate (cE tag)
    // Note: Promote f32 to f64 for comparison (matching fgbio's Float->Double promotion)
    let error_tag = per_read::tag("cE");
    if let Some(Value::Float(error_rate)) = record.data().get(&error_tag) {
        if f64::from(*error_rate) > thresholds.max_read_error_rate {
            return Ok(FilterResult::ExcessiveErrorRate);
        }
    }

    Ok(FilterResult::Pass)
}

/// Filters a duplex consensus read based on per-read tags, checking both AB and BA strands
///
/// # Arguments
/// * `record` - The consensus read to filter
/// * `cc_thresholds` - Thresholds for final consensus
/// * `ab_thresholds` - Thresholds for AB strand
/// * `ba_thresholds` - Thresholds for BA strand
///
/// # Returns
/// `FilterResult` indicating whether the read passed or why it failed
///
/// # Errors
///
/// Returns an error if the record data cannot be read.
pub fn filter_duplex_read(
    record: &RecordBuf,
    cc_thresholds: &FilterThresholds,
    ab_thresholds: &FilterThresholds,
    ba_thresholds: &FilterThresholds,
) -> Result<FilterResult> {
    let data = record.data();

    // First check final consensus thresholds
    let result = filter_read(record, cc_thresholds)?;
    if result != FilterResult::Pass {
        return Ok(result);
    }

    // Extract AB and BA depths and error rates
    let ab_depth = extract_int_tag(data, "aD").or_else(|| extract_int_tag(data, "aM"));
    let ba_depth = extract_int_tag(data, "bD").or_else(|| extract_int_tag(data, "bM"));
    let ab_error = extract_float_tag(data, "aE");
    let ba_error = extract_float_tag(data, "bE");

    // Sort depths and errors to identify which is AB (higher) and BA (lower)
    // Following Scala: val Seq(baMaxDepth, abMaxDepth) = Seq(...).sorted
    let (min_depth, max_depth) = match (ab_depth, ba_depth) {
        (Some(a), Some(b)) => {
            if a < b {
                (a, b)
            } else {
                (b, a)
            }
        }
        (Some(a), None) => (0, a),
        (None, Some(b)) => (0, b),
        (None, None) => return Ok(FilterResult::Pass), // No AB/BA tags, not duplex
    };

    let (min_error, max_error) = match (ab_error, ba_error) {
        (Some(a), Some(b)) => {
            if a < b {
                (a, b)
            } else {
                (b, a)
            }
        }
        (Some(a), None) => (a, a),
        (None, Some(b)) => (b, b),
        (None, None) => (0.0, 0.0),
    };

    // Check AB strand (max depth, min error) against AB thresholds
    // Note: Promote f32 to f64 for comparison (matching fgbio's Float->Double promotion)
    #[expect(clippy::cast_sign_loss, reason = "depth values are non-negative in BAM tags")]
    if (max_depth as usize) < ab_thresholds.min_reads {
        return Ok(FilterResult::InsufficientReads);
    }
    if f64::from(min_error) > ab_thresholds.max_read_error_rate {
        return Ok(FilterResult::ExcessiveErrorRate);
    }

    // Check BA strand (min depth, max error) against BA thresholds
    #[expect(clippy::cast_sign_loss, reason = "depth values are non-negative in BAM tags")]
    if (min_depth as usize) < ba_thresholds.min_reads {
        return Ok(FilterResult::InsufficientReads);
    }
    if f64::from(max_error) > ba_thresholds.max_read_error_rate {
        return Ok(FilterResult::ExcessiveErrorRate);
    }

    Ok(FilterResult::Pass)
}

/// Helper to extract per-base array tag
fn extract_per_base_array(
    data: &noodles::sam::alignment::record_buf::Data,
    tag_str: &str,
    len: usize,
) -> Vec<u16> {
    let tag = per_base::tag(tag_str);
    if let Some(Value::Array(arr)) = data.get(&tag) {
        use noodles::sam::alignment::record_buf::data::field::value::Array;
        match arr {
            #[expect(clippy::cast_sign_loss, reason = "clamped to non-negative via .max(0)")]
            Array::Int16(values) => values.iter().map(|&v| v.max(0) as u16).collect(),
            Array::UInt16(values) => values.clone(),
            Array::UInt8(values) => values.iter().map(|&v| u16::from(v)).collect(),
            _ => vec![0; len],
        }
    } else {
        vec![0; len]
    }
}

/// Helper to extract consensus base string tag (ac/bc) for single-strand agreement checking
fn extract_consensus_bases(
    data: &noodles::sam::alignment::record_buf::Data,
    tag_str: &str,
    len: usize,
) -> Vec<u8> {
    let tag = per_base::tag(tag_str);
    if let Some(value) = data.get(&tag) {
        match value {
            Value::String(s) => {
                let bytes: &[u8] = s.as_ref();
                bytes.to_vec()
            }
            Value::Array(arr) => {
                use noodles::sam::alignment::record_buf::data::field::value::Array;
                match arr {
                    Array::UInt8(bytes) => bytes.clone(),
                    _ => vec![NO_CALL_BASE; len],
                }
            }
            _ => vec![NO_CALL_BASE; len],
        }
    } else {
        vec![NO_CALL_BASE; len]
    }
}

/// Masks bases in a consensus read based on per-base tags and thresholds
///
/// # Arguments
/// * `record` - The consensus read to mask (will be modified in place)
/// * `thresholds` - The filter thresholds to apply
/// * `min_base_quality` - Minimum base quality to keep (None means no quality masking)
///
/// # Returns
/// Number of bases masked
///
/// # Errors
///
/// Returns an error if the record data cannot be read or modified.
pub fn mask_bases(
    record: &mut RecordBuf,
    thresholds: &FilterThresholds,
    min_base_quality: Option<u8>,
) -> Result<usize> {
    // Collect sequence and quality scores directly into mutable vectors (no clone needed)
    let mut new_seq: Vec<u8> = record.sequence().iter().collect();
    let mut new_quals: Vec<u8> = record.quality_scores().iter().collect();
    let mut masked_count = 0;
    let len = new_seq.len();

    // Get per-base depth and error counts using helper
    let data = record.data();
    let depths = extract_per_base_array(data, "cd", len);
    let errors = extract_per_base_array(data, "ce", len);

    // Mask bases that don't meet thresholds
    for i in 0..len {
        let should_mask =
            // Quality too low (only check if min_base_quality is specified)
            min_base_quality.is_some_and(|min_qual| new_quals[i] < min_qual) ||
            // Depth too low
            (depths[i] as usize) < thresholds.min_reads ||
            // Error rate too high
            (depths[i] > 0 && (f64::from(errors[i]) / f64::from(depths[i])) > thresholds.max_base_error_rate);

        if should_mask {
            // Only count as newly masked if not already N
            if new_seq[i] != NO_CALL_BASE {
                masked_count += 1;
            }
            new_seq[i] = NO_CALL_BASE;
            new_quals[i] = MIN_PHRED;
        }
    }

    // Update the record with masked bases
    *record.sequence_mut() = Sequence::from(new_seq);
    *record.quality_scores_mut() = QualityScores::from(new_quals);

    Ok(masked_count)
}

/// Masks bases in a duplex consensus read based on per-base AB/BA tags and thresholds
///
/// This function implements the duplex per-base filtering logic from FilterConsensusReads.scala
/// (DuplexConsensusPerBaseValues.maskBaseAt, lines 422-435).
///
/// For each base position, checks:
/// - Total depth (AB + BA) against `cc_thresholds`
/// - Total error rate against `cc_thresholds`
/// - AB depth (max of two strands) against `ab_thresholds`
/// - AB error rate (min of two strands) against `ab_thresholds`
/// - BA depth (min of two strands) against `ba_thresholds`
/// - BA error rate (max of two strands) against `ba_thresholds`
/// - Optionally, single-strand agreement
///
/// # Arguments
/// * `record` - The consensus read to mask (will be modified in place)
/// * `cc_thresholds` - Filter thresholds for final consensus
/// * `ab_thresholds` - Filter thresholds for AB strand
/// * `ba_thresholds` - Filter thresholds for BA strand
/// * `min_base_quality` - Minimum base quality to keep (None means no quality masking)
/// * `require_ss_agreement` - If true, mask bases where AB and BA disagree
///
/// # Returns
/// Number of bases masked
///
/// # Errors
///
/// Returns an error if the record data cannot be read or modified.
pub fn mask_duplex_bases(
    record: &mut RecordBuf,
    cc_thresholds: &FilterThresholds,
    ab_thresholds: &FilterThresholds,
    ba_thresholds: &FilterThresholds,
    min_base_quality: Option<u8>,
    require_ss_agreement: bool,
) -> Result<usize> {
    // Collect sequence and quality scores directly into mutable vectors (no clone needed)
    let mut new_seq: Vec<u8> = record.sequence().iter().collect();
    let mut new_quals: Vec<u8> = record.quality_scores().iter().collect();
    let len = new_seq.len();
    let mut masked_count = 0;

    let data = record.data();

    // Extract AB per-base tags (ad, ae)
    let ab_depths = extract_per_base_array(data, "ad", len);
    let ab_errors = extract_per_base_array(data, "ae", len);

    // Extract BA per-base tags (bd, be)
    let ba_depths = extract_per_base_array(data, "bd", len);
    let ba_errors = extract_per_base_array(data, "be", len);

    // Extract consensus bases for single-strand agreement checking (ac, bc)
    let ab_bases = if require_ss_agreement {
        extract_consensus_bases(data, "ac", len)
    } else {
        vec![NO_CALL_BASE; len]
    };
    let ba_bases = if require_ss_agreement {
        extract_consensus_bases(data, "bc", len)
    } else {
        vec![NO_CALL_BASE; len]
    };

    // Mask bases that don't meet thresholds
    for i in 0..len {
        if new_seq[i] == NO_CALL_BASE {
            continue;
        }

        // Following Scala's maskBaseAt logic (lines 422-435)
        let ab_depth_i = ab_depths[i];
        let ba_depth_i = ba_depths[i];
        let ab_error_i = ab_errors[i];
        let ba_error_i = ba_errors[i];

        // Sort to get max/min depths and errors
        let max_depth = std::cmp::max(ab_depth_i, ba_depth_i);
        let min_depth = std::cmp::min(ab_depth_i, ba_depth_i);

        // For error rates, we need to calculate them first, then sort
        let ab_error_rate =
            if ab_depth_i > 0 { f64::from(ab_error_i) / f64::from(ab_depth_i) } else { 0.0 };
        let ba_error_rate =
            if ba_depth_i > 0 { f64::from(ba_error_i) / f64::from(ba_depth_i) } else { 0.0 };

        let min_error_rate = ab_error_rate.min(ba_error_rate);
        let max_error_rate = ab_error_rate.max(ba_error_rate);

        // Total depth and error for final consensus check
        let total_depth = ab_depth_i + ba_depth_i;
        let total_error_rate = if total_depth > 0 {
            f64::from(ab_error_i + ba_error_i) / f64::from(total_depth)
        } else {
            0.0
        };

        // Check if base should be masked
        let should_mask = min_base_quality.is_some_and(|min_qual| new_quals[i] < min_qual)
            || (total_depth as usize) < cc_thresholds.min_reads
            || total_error_rate > cc_thresholds.max_base_error_rate
            || (max_depth as usize) < ab_thresholds.min_reads
            || min_error_rate > ab_thresholds.max_base_error_rate
            || (min_depth as usize) < ba_thresholds.min_reads
            || max_error_rate > ba_thresholds.max_base_error_rate;

        // Check single-strand agreement if requested
        let ss_disagree =
            require_ss_agreement && ab_depth_i > 0 && ba_depth_i > 0 && ab_bases[i] != ba_bases[i];

        if should_mask || ss_disagree {
            // Only count as newly masked if not already N
            if new_seq[i] != NO_CALL_BASE {
                masked_count += 1;
            }
            new_seq[i] = NO_CALL_BASE;
            new_quals[i] = MIN_PHRED;
        }
    }

    // Update the record with masked bases
    *record.sequence_mut() = Sequence::from(new_seq);
    *record.quality_scores_mut() = QualityScores::from(new_quals);

    Ok(masked_count)
}

/// Counts the number of N bases in a read
#[must_use]
pub fn count_no_calls(record: &RecordBuf) -> usize {
    record.sequence().iter().filter(|b| *b == NO_CALL_BASE).count()
}

/// Calculates the mean base quality of non-N bases
#[must_use]
pub fn mean_base_quality(record: &RecordBuf) -> f64 {
    // Use iterators directly to avoid allocating vectors
    let (sum, count) = record
        .sequence()
        .iter()
        .zip(record.quality_scores().iter())
        .filter(|&(base, _)| base != NO_CALL_BASE)
        .fold((0u64, 0usize), |(sum, count), (_, qual)| (sum + u64::from(qual), count + 1));

    #[expect(clippy::cast_precision_loss, reason = "precision loss is acceptable for quality averaging")]
    if count == 0 { 0.0 } else { sum as f64 / count as f64 }
}

/// Computes both no-call count and mean base quality in a single pass.
///
/// This is more efficient than calling `count_no_calls()` and `mean_base_quality()`
/// separately when both values are needed, as it only iterates over the sequence once.
///
/// # Returns
/// A tuple of (`no_call_count`, `mean_base_quality`)
#[must_use]
pub fn compute_read_stats(record: &RecordBuf) -> (usize, f64) {
    let (qual_sum, non_n_count, n_count) = record
        .sequence()
        .iter()
        .zip(record.quality_scores().iter())
        .fold((0u64, 0usize, 0usize), |(sum, non_n, n), (base, qual)| {
            if base == NO_CALL_BASE {
                (sum, non_n, n + 1)
            } else {
                (sum + u64::from(qual), non_n + 1, n)
            }
        });

    #[expect(clippy::cast_precision_loss, reason = "precision loss is acceptable for quality averaging")]
    let mean_qual = if non_n_count == 0 { 0.0 } else { qual_sum as f64 / non_n_count as f64 };
    (n_count, mean_qual)
}

/// Checks if all reads in a template pass filtering
///
/// For template-aware filtering, all primary reads (R1/R2) must pass for the template to pass.
/// Secondary and supplementary alignments are allowed to fail independently.
///
/// # Arguments
/// * `template_records` - All records with the same read name
/// * `pass_map` - Map of record index to pass/fail status
///
/// # Returns
/// true if all primary reads pass, false otherwise
#[must_use]
pub fn template_passes(template_records: &[RecordBuf], pass_map: &AHashMap<usize, bool>) -> bool {
    // Find primary reads (non-secondary, non-supplementary)
    let mut has_primary = false;
    let mut all_primary_pass = true;

    for (idx, record) in template_records.iter().enumerate() {
        let flags = record.flags();
        let is_primary = !flags.is_secondary() && !flags.is_supplementary();

        if is_primary {
            has_primary = true;
            if let Some(&passes) = pass_map.get(&idx) {
                if !passes {
                    all_primary_pass = false;
                    break;
                }
            } else {
                // If not in map, assume it failed
                all_primary_pass = false;
                break;
            }
        }
    }

    // If no primary reads found, template fails
    // If all primary reads pass, template passes
    has_primary && all_primary_pass
}

/// Checks whether all primary raw records in a template pass their filters.
///
/// A template passes if it has at least one primary read and all primary reads pass.
#[must_use]
pub fn template_passes_raw(raw_records: &[Vec<u8>], pass_map: &AHashMap<usize, bool>) -> bool {
    let mut has_primary = false;
    let mut all_primary_pass = true;

    for (idx, record) in raw_records.iter().enumerate() {
        let flags = bam_fields::flags(record);
        let is_primary = (flags & bam_fields::flags::SECONDARY) == 0
            && (flags & bam_fields::flags::SUPPLEMENTARY) == 0;

        if is_primary {
            has_primary = true;
            if let Some(&passes) = pass_map.get(&idx) {
                if !passes {
                    all_primary_pass = false;
                    break;
                }
            } else {
                all_primary_pass = false;
                break;
            }
        }
    }

    has_primary && all_primary_pass
}

/// Detects if a record is part of a duplex consensus
///
/// Checks for presence of AB/BA specific tags (aD, bD)
///
/// # Arguments
/// * `record` - The record to check
///
/// # Returns
/// true if duplex tags are present, false otherwise
#[must_use]
pub fn is_duplex_consensus(record: &RecordBuf) -> bool {
    let ad_tag = per_read::tag("aD");
    let bd_tag = per_read::tag("bD");

    record.data().get(&ad_tag).is_some() || record.data().get(&bd_tag).is_some()
}

// ============================================================================
// Raw-byte equivalents (operate on &[u8] / &mut Vec<u8>)
// ============================================================================

use noodles_raw_bam as bam_fields;

/// Detects if a raw BAM record is a duplex consensus.
///
/// Checks for presence of `aD` or `bD` tags in the aux data.
#[must_use]
pub fn is_duplex_consensus_raw(aux_data: &[u8]) -> bool {
    bam_fields::find_tag_type(aux_data, b"aD").is_some()
        || bam_fields::find_tag_type(aux_data, b"bD").is_some()
}

/// Filters a raw consensus read based on per-read tags (cD depth, cE error rate).
///
/// # Errors
///
/// Returns an error if the aux data cannot be parsed.
pub fn filter_read_raw(aux_data: &[u8], thresholds: &FilterThresholds) -> Result<FilterResult> {
    // Check minimum reads (cD tag — UInt8)
    if let Some(depth) = bam_fields::find_uint8_tag(aux_data, b"cD") {
        if (depth as usize) < thresholds.min_reads {
            return Ok(FilterResult::InsufficientReads);
        }
    }

    // Check maximum error rate (cE tag — Float)
    if let Some(error_rate) = bam_fields::find_float_tag(aux_data, b"cE") {
        if f64::from(error_rate) > thresholds.max_read_error_rate {
            return Ok(FilterResult::ExcessiveErrorRate);
        }
    }

    Ok(FilterResult::Pass)
}

/// Filters a raw duplex consensus read, checking CC / AB / BA thresholds.
///
/// # Errors
///
/// Returns an error if the aux data cannot be parsed.
pub fn filter_duplex_read_raw(
    aux_data: &[u8],
    cc_thresholds: &FilterThresholds,
    ab_thresholds: &FilterThresholds,
    ba_thresholds: &FilterThresholds,
) -> Result<FilterResult> {
    // First check final consensus thresholds
    let result = filter_read_raw(aux_data, cc_thresholds)?;
    if result != FilterResult::Pass {
        return Ok(result);
    }

    // Extract AB and BA depths and error rates
    let ab_depth = bam_fields::find_int_tag(aux_data, b"aD")
        .or_else(|| bam_fields::find_int_tag(aux_data, b"aM"));
    let ba_depth = bam_fields::find_int_tag(aux_data, b"bD")
        .or_else(|| bam_fields::find_int_tag(aux_data, b"bM"));
    let ab_error = bam_fields::find_float_tag(aux_data, b"aE");
    let ba_error = bam_fields::find_float_tag(aux_data, b"bE");

    // Sort depths to identify AB (higher) and BA (lower)
    let (min_depth, max_depth) = match (ab_depth, ba_depth) {
        (Some(a), Some(b)) => {
            if a < b {
                (a, b)
            } else {
                (b, a)
            }
        }
        (Some(a), None) => (0, a),
        (None, Some(b)) => (0, b),
        (None, None) => return Ok(FilterResult::Pass),
    };

    let (min_error, max_error) = match (ab_error, ba_error) {
        (Some(a), Some(b)) => {
            if a < b {
                (a, b)
            } else {
                (b, a)
            }
        }
        (Some(a), None) => (a, a),
        (None, Some(b)) => (b, b),
        (None, None) => (0.0, 0.0),
    };

    // Check AB strand (max depth, min error) against AB thresholds
    #[expect(clippy::cast_sign_loss, clippy::cast_possible_truncation, reason = "depth values are non-negative and fit in usize on all supported platforms")]
    if (max_depth as usize) < ab_thresholds.min_reads {
        return Ok(FilterResult::InsufficientReads);
    }
    if f64::from(min_error) > ab_thresholds.max_read_error_rate {
        return Ok(FilterResult::ExcessiveErrorRate);
    }

    // Check BA strand (min depth, max error) against BA thresholds
    #[expect(clippy::cast_sign_loss, clippy::cast_possible_truncation, reason = "depth values are non-negative and fit in usize on all supported platforms")]
    if (min_depth as usize) < ba_thresholds.min_reads {
        return Ok(FilterResult::InsufficientReads);
    }
    if f64::from(max_error) > ba_thresholds.max_read_error_rate {
        return Ok(FilterResult::ExcessiveErrorRate);
    }

    Ok(FilterResult::Pass)
}

/// Computes both no-call count and mean base quality in a single pass over raw BAM bytes.
///
/// Returns (`no_call_count`, `mean_base_quality`).
///
/// # Panics
/// Panics if the record is shorter than `MIN_BAM_HEADER_LEN` (36 bytes).
#[must_use]
pub fn compute_read_stats_raw(bam: &[u8]) -> (usize, f64) {
    assert!(bam.len() >= bam_fields::MIN_BAM_HEADER_LEN, "BAM record too short");
    let seq_off = bam_fields::seq_offset(bam);
    let qual_off = bam_fields::qual_offset(bam);
    let len = bam_fields::l_seq(bam) as usize;

    let mut n_count = 0usize;
    let mut qual_sum = 0u64;
    let mut non_n_count = 0usize;

    for i in 0..len {
        if bam_fields::is_base_n(bam, seq_off, i) {
            n_count += 1;
        } else {
            qual_sum += u64::from(bam_fields::get_qual(bam, qual_off, i));
            non_n_count += 1;
        }
    }

    #[expect(clippy::cast_precision_loss, reason = "precision loss is acceptable for quality averaging")]
    let mean_qual = if non_n_count == 0 { 0.0 } else { qual_sum as f64 / non_n_count as f64 };
    (n_count, mean_qual)
}

/// Reads a tag value as either a Z-type string or a B-type `UInt8` array.
///
/// Returns `Some(Vec<u8>)` if found as either type, `None` otherwise.
fn find_string_or_uint8_array(aux_data: &[u8], tag: [u8; 2]) -> Option<Vec<u8>> {
    if let Some(s) = bam_fields::find_string_tag(aux_data, &tag) {
        Some(s.to_vec())
    } else {
        let arr = bam_fields::find_array_tag(aux_data, &tag)?;
        if matches!(arr.elem_type, b'C' | b'c') {
            #[expect(clippy::cast_possible_truncation, reason = "UInt8 array elements are guaranteed to fit in u8")]
            Some((0..arr.count).map(|i| bam_fields::array_tag_element_u16(&arr, i) as u8).collect())
        } else {
            None
        }
    }
}

/// Masks bases in a raw consensus read based on per-base tags and thresholds.
///
/// Modifies sequence and quality bytes in-place (no Vec allocation for the seq/qual data).
/// Returns the number of newly masked bases.
///
/// # Errors
///
/// Returns an error if the record is too short or the aux data cannot be parsed.
#[allow(clippy::similar_names)]
pub fn mask_bases_raw(
    record: &mut [u8],
    thresholds: &FilterThresholds,
    min_base_quality: Option<u8>,
) -> Result<usize> {
    anyhow::ensure!(record.len() >= bam_fields::MIN_BAM_HEADER_LEN, "BAM record too short");
    let seq_off = bam_fields::seq_offset(record);
    let qual_off = bam_fields::qual_offset(record);
    let len = bam_fields::l_seq(record) as usize;
    let aux_off = bam_fields::aux_data_offset_from_record(record).unwrap_or(record.len());

    // Pre-read per-base arrays into owned Vecs to release the immutable borrow on record
    let cd_vals = bam_fields::find_array_tag(&record[aux_off..], b"cd")
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));
    let ce_vals = bam_fields::find_array_tag(&record[aux_off..], b"ce")
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));

    let mut masked_count = 0;
    for i in 0..len {
        let depth = cd_vals.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let errors = ce_vals.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let qual = bam_fields::get_qual(record, qual_off, i);

        let should_mask = min_base_quality.is_some_and(|min_qual| qual < min_qual)
            || (depth as usize) < thresholds.min_reads
            || (depth > 0
                && (f64::from(errors) / f64::from(depth)) > thresholds.max_base_error_rate);

        if should_mask {
            // Only count as newly masked if not already N
            if !bam_fields::is_base_n(record, seq_off, i) {
                masked_count += 1;
            }
            bam_fields::mask_base(record, seq_off, i);
            bam_fields::set_qual(record, qual_off, i, MIN_PHRED);
        }
    }

    Ok(masked_count)
}

/// Masks bases in a raw duplex consensus read based on per-base AB/BA tags and thresholds.
///
/// Returns the number of newly masked bases.
///
/// # Errors
///
/// Returns an error if the record is too short or the aux data cannot be parsed.
#[allow(clippy::too_many_arguments, clippy::similar_names)]
pub fn mask_duplex_bases_raw(
    record: &mut [u8],
    cc_thresholds: &FilterThresholds,
    ab_thresholds: &FilterThresholds,
    ba_thresholds: &FilterThresholds,
    min_base_quality: Option<u8>,
    require_ss_agreement: bool,
) -> Result<usize> {
    anyhow::ensure!(record.len() >= bam_fields::MIN_BAM_HEADER_LEN, "BAM record too short");
    let seq_off = bam_fields::seq_offset(record);
    let qual_off = bam_fields::qual_offset(record);
    let len = bam_fields::l_seq(record) as usize;
    let aux_off = bam_fields::aux_data_offset_from_record(record).unwrap_or(record.len());

    // Pre-read per-base arrays and strings into owned data to release the immutable borrow
    let ad_vals = bam_fields::find_array_tag(&record[aux_off..], b"ad")
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));
    let ae_vals = bam_fields::find_array_tag(&record[aux_off..], b"ae")
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));
    let bd_vals = bam_fields::find_array_tag(&record[aux_off..], b"bd")
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));
    let be_vals = bam_fields::find_array_tag(&record[aux_off..], b"be")
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));

    // For single-strand agreement checking, get ac/bc tags (copy to owned).
    // These may be Z-type strings or B-type UInt8 arrays.
    let ac_owned: Option<Vec<u8>> = if require_ss_agreement {
        find_string_or_uint8_array(&record[aux_off..], *b"ac")
    } else {
        None
    };
    let bc_owned: Option<Vec<u8>> = if require_ss_agreement {
        find_string_or_uint8_array(&record[aux_off..], *b"bc")
    } else {
        None
    };

    let mut masked_count = 0;
    for i in 0..len {
        // Skip if already N
        if bam_fields::is_base_n(record, seq_off, i) {
            continue;
        }

        let ab_depth = ad_vals.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let ba_depth = bd_vals.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let ab_errors = ae_vals.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let ba_errors = be_vals.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));

        let max_depth = std::cmp::max(ab_depth, ba_depth);
        let min_depth = std::cmp::min(ab_depth, ba_depth);

        let ab_error_rate =
            if ab_depth > 0 { f64::from(ab_errors) / f64::from(ab_depth) } else { 0.0 };
        let ba_error_rate =
            if ba_depth > 0 { f64::from(ba_errors) / f64::from(ba_depth) } else { 0.0 };

        let min_error_rate = ab_error_rate.min(ba_error_rate);
        let max_error_rate = ab_error_rate.max(ba_error_rate);

        let total_depth = u32::from(ab_depth) + u32::from(ba_depth);
        let total_error_rate = if total_depth > 0 {
            f64::from(u32::from(ab_errors) + u32::from(ba_errors)) / f64::from(total_depth)
        } else {
            0.0
        };

        let qual = bam_fields::get_qual(record, qual_off, i);

        let should_mask = min_base_quality.is_some_and(|min_qual| qual < min_qual)
            || (total_depth as usize) < cc_thresholds.min_reads
            || total_error_rate > cc_thresholds.max_base_error_rate
            || (max_depth as usize) < ab_thresholds.min_reads
            || min_error_rate > ab_thresholds.max_base_error_rate
            || (min_depth as usize) < ba_thresholds.min_reads
            || max_error_rate > ba_thresholds.max_base_error_rate;

        // Check single-strand agreement if requested.
        // Use NO_CALL_BASE as default for missing/short tags, matching RecordBuf path behavior.
        let ss_disagree = require_ss_agreement && ab_depth > 0 && ba_depth > 0 && {
            let ac_base =
                ac_owned.as_ref().and_then(|ac| ac.get(i).copied()).unwrap_or(NO_CALL_BASE);
            let bc_base =
                bc_owned.as_ref().and_then(|bc| bc.get(i).copied()).unwrap_or(NO_CALL_BASE);
            ac_base != bc_base
        };

        if should_mask || ss_disagree {
            masked_count += 1;
            bam_fields::mask_base(record, seq_off, i);
            bam_fields::set_qual(record, qual_off, i, MIN_PHRED);
        }
    }

    Ok(masked_count)
}

#[cfg(test)]
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
mod tests {
    use super::*;
    use fgumi_sam::builder::RecordBuilder;
    use noodles::sam::alignment::record::data::field::Tag;

    fn create_test_record() -> RecordBuf {
        // Uses RecordBuilder which auto-generates Q30 qualities
        RecordBuilder::new().sequence("ACGTACGT").build()
    }

    #[test]
    fn test_filter_result_to_rejection_reason() {
        assert_eq!(FilterResult::Pass.to_rejection_reason(), None);
        assert_eq!(
            FilterResult::InsufficientReads.to_rejection_reason(),
            Some(RejectionReason::InsufficientSupport)
        );
        assert_eq!(
            FilterResult::ExcessiveErrorRate.to_rejection_reason(),
            Some(RejectionReason::ExcessiveErrorRate)
        );
        assert_eq!(
            FilterResult::LowQuality.to_rejection_reason(),
            Some(RejectionReason::LowMeanQuality)
        );
        assert_eq!(
            FilterResult::TooManyNoCalls.to_rejection_reason(),
            Some(RejectionReason::ExcessiveNBases)
        );
    }

    #[test]
    fn test_filter_config_single_value() {
        let config = FilterConfig::new(&[3], &[0.1], &[0.2], Some(13), None, 0.1);

        // All threshold types are now always created (matching fgbio behavior)
        assert!(config.single_strand_thresholds.is_some());
        assert!(config.duplex_thresholds.is_some());

        let thresh = config.single_strand_thresholds.unwrap();
        assert_eq!(thresh.min_reads, 3);
        assert!((thresh.max_read_error_rate - 0.1).abs() < f64::EPSILON);
        assert!((thresh.max_base_error_rate - 0.2).abs() < f64::EPSILON);
    }

    #[test]
    fn test_filter_config_three_values() {
        let config =
            FilterConfig::new(&[5, 3, 3], &[0.05, 0.1, 0.1], &[0.1, 0.2, 0.2], Some(13), None, 0.1);

        // All threshold types are now always created (matching fgbio behavior)
        assert!(config.duplex_thresholds.is_some());
        assert!(config.ab_thresholds.is_some());
        assert!(config.ba_thresholds.is_some());
        assert!(config.single_strand_thresholds.is_some());

        let duplex = config.duplex_thresholds.unwrap();
        assert_eq!(duplex.min_reads, 5);
        assert!((duplex.max_read_error_rate - 0.05).abs() < f64::EPSILON);

        let ab = config.ab_thresholds.unwrap();
        assert_eq!(ab.min_reads, 3);
        assert!((ab.max_read_error_rate - 0.1).abs() < f64::EPSILON);
    }

    #[test]
    fn test_filter_config_valid_threshold_ordering() {
        // Valid: CC=10 >= AB=5 >= BA=3 for min_reads
        // Valid: AB=0.05 <= BA=0.1 for error rates (AB more stringent)
        let config = FilterConfig::new(
            &[10, 5, 3],
            &[0.02, 0.05, 0.1],
            &[0.05, 0.1, 0.2],
            Some(13),
            None,
            0.1,
        );

        let cc = config.duplex_thresholds.unwrap();
        let ab = config.ab_thresholds.unwrap();
        let ba = config.ba_thresholds.unwrap();

        assert_eq!(cc.min_reads, 10);
        assert_eq!(ab.min_reads, 5);
        assert_eq!(ba.min_reads, 3);
    }

    // ========== Tests for named constructors ==========

    #[test]
    fn test_filter_config_for_single_strand() {
        let thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };

        let config = FilterConfig::for_single_strand(thresholds.clone(), Some(13), None, 0.1);

        // All threshold types should use the same values
        let ss = config.single_strand_thresholds.unwrap();
        assert_eq!(ss.min_reads, 3);
        assert!((ss.max_read_error_rate - 0.1).abs() < f64::EPSILON);
        assert!((ss.max_base_error_rate - 0.2).abs() < f64::EPSILON);

        // Duplex, AB, BA should also be set (for compatibility)
        assert!(config.duplex_thresholds.is_some());
        assert!(config.ab_thresholds.is_some());
        assert!(config.ba_thresholds.is_some());
    }

    #[test]
    fn test_filter_config_for_duplex_symmetric() {
        let duplex_thresholds =
            FilterThresholds { min_reads: 10, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
        let strand_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };

        let config =
            FilterConfig::for_duplex(duplex_thresholds, strand_thresholds, Some(13), None, 0.1);

        let cc = config.duplex_thresholds.unwrap();
        let ab = config.ab_thresholds.unwrap();
        let ba = config.ba_thresholds.unwrap();

        assert_eq!(cc.min_reads, 10);
        assert_eq!(ab.min_reads, 5);
        assert_eq!(ba.min_reads, 5); // Same as AB (symmetric)
    }

    #[test]
    fn test_filter_config_for_duplex_asymmetric() {
        let duplex = FilterThresholds {
            min_reads: 10,
            max_read_error_rate: 0.02,
            max_base_error_rate: 0.05,
        };
        let ab =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
        let ba =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };

        let config =
            FilterConfig::for_duplex_asymmetric(duplex, ab, ba, Some(13), Some(20.0), 0.15);

        let cc = config.duplex_thresholds.unwrap();
        let ab_t = config.ab_thresholds.unwrap();
        let ba_t = config.ba_thresholds.unwrap();

        assert_eq!(cc.min_reads, 10);
        assert_eq!(ab_t.min_reads, 5);
        assert_eq!(ba_t.min_reads, 3);

        // Verify other config fields
        assert_eq!(config.min_base_quality, Some(13));
        assert_eq!(config.min_mean_base_quality, Some(20.0));
        assert!((config.max_no_call_fraction - 0.15).abs() < f64::EPSILON);
    }

    #[test]
    #[should_panic(expected = "min-reads values must be specified high to low: AB")]
    fn test_filter_config_for_duplex_asymmetric_invalid_ab() {
        let duplex =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
        let ab =
            FilterThresholds { min_reads: 10, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };
        let ba =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };

        // AB (10) > duplex (5) - should panic
        let _ = FilterConfig::for_duplex_asymmetric(duplex, ab, ba, Some(13), None, 0.1);
    }

    #[test]
    #[should_panic(expected = "min-reads values must be specified high to low: AB")]
    fn test_filter_config_invalid_ab_greater_than_cc() {
        // Invalid: AB (6) > CC (5) for min_reads
        let _ =
            FilterConfig::new(&[5, 6, 3], &[0.05, 0.1, 0.1], &[0.1, 0.2, 0.2], Some(13), None, 0.1);
    }

    #[test]
    #[should_panic(expected = "min-reads values must be specified high to low: BA")]
    fn test_filter_config_invalid_ba_greater_than_ab() {
        // Invalid: BA (4) > AB (3) for min_reads
        let _ =
            FilterConfig::new(&[5, 3, 4], &[0.05, 0.1, 0.1], &[0.1, 0.2, 0.2], Some(13), None, 0.1);
    }

    #[test]
    #[should_panic(expected = "max-read-error-rate for AB")]
    fn test_filter_config_invalid_ab_read_error_rate_greater_than_ba() {
        // Invalid: AB read error rate (0.2) > BA read error rate (0.1)
        let _ =
            FilterConfig::new(&[5, 3, 3], &[0.05, 0.2, 0.1], &[0.1, 0.2, 0.2], Some(13), None, 0.1);
    }

    #[test]
    #[should_panic(expected = "max-base-error-rate for AB")]
    fn test_filter_config_invalid_ab_base_error_rate_greater_than_ba() {
        // Invalid: AB base error rate (0.3) > BA base error rate (0.2)
        let _ =
            FilterConfig::new(&[5, 3, 3], &[0.05, 0.1, 0.1], &[0.1, 0.3, 0.2], Some(13), None, 0.1);
    }

    #[test]
    fn test_count_no_calls() {
        let mut record = create_test_record();
        *record.sequence_mut() = Sequence::from(b"ACNNGCNN".to_vec());

        assert_eq!(count_no_calls(&record), 4);
    }

    #[test]
    fn test_mean_base_quality() {
        let mut record = create_test_record();
        *record.quality_scores_mut() = QualityScores::from(vec![10, 20, 30, 40, 10, 20, 30, 40]);

        let mean = mean_base_quality(&record);
        assert!((mean - 25.0).abs() < f64::EPSILON); // (10+20+30+40+10+20+30+40) / 8 = 200/8 = 25
    }

    #[test]
    fn test_mean_base_quality_with_n() {
        let mut record = create_test_record();
        *record.sequence_mut() = Sequence::from(b"ACNNACGT".to_vec());
        *record.quality_scores_mut() = QualityScores::from(vec![10, 20, 0, 0, 30, 40, 10, 20]);

        let mean = mean_base_quality(&record);
        // Only non-N bases: 10+20+30+40+10+20 = 130 / 6 = 21.666...
        assert!((mean - 21.666).abs() < 0.01);
    }

    #[test]
    fn test_compute_read_stats() {
        let mut record = create_test_record();
        *record.sequence_mut() = Sequence::from(b"ACNNACGT".to_vec());
        *record.quality_scores_mut() = QualityScores::from(vec![10, 20, 0, 0, 30, 40, 10, 20]);

        let (n_count, mean_qual) = compute_read_stats(&record);

        // Should match individual functions
        assert_eq!(n_count, count_no_calls(&record));
        assert!((mean_qual - mean_base_quality(&record)).abs() < f64::EPSILON);

        // Verify expected values: 2 Ns, mean of non-N quals = 130/6 = 21.666...
        assert_eq!(n_count, 2);
        assert!((mean_qual - 21.666).abs() < 0.01);
    }

    #[test]
    fn test_compute_read_stats_all_n() {
        let mut record = create_test_record();
        *record.sequence_mut() = Sequence::from(b"NNNNNNNN".to_vec());
        *record.quality_scores_mut() = QualityScores::from(vec![0; 8]);

        let (n_count, mean_qual) = compute_read_stats(&record);

        assert_eq!(n_count, 8);
        assert!((mean_qual - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_compute_read_stats_no_n() {
        let mut record = create_test_record();
        *record.quality_scores_mut() = QualityScores::from(vec![30, 30, 30, 30, 30, 30, 30, 30]);

        let (n_count, mean_qual) = compute_read_stats(&record);

        assert_eq!(n_count, 0);
        assert!((mean_qual - 30.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_filter_result() {
        let thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };

        let mut record = create_test_record();

        // Add cD tag with sufficient depth
        let tag = Tag::from([b'c', b'D']);
        record.data_mut().insert(tag, Value::UInt8(5));

        let result = filter_read(&record, &thresholds).unwrap();
        assert_eq!(result, FilterResult::Pass);
    }

    #[test]
    fn test_filter_insufficient_reads() {
        let thresholds =
            FilterThresholds { min_reads: 10, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };

        let mut record = create_test_record();

        // Add cD tag with insufficient depth
        let tag = Tag::from([b'c', b'D']);
        record.data_mut().insert(tag, Value::UInt8(5));

        let result = filter_read(&record, &thresholds).unwrap();
        assert_eq!(result, FilterResult::InsufficientReads);
    }

    #[test]
    fn test_template_passes_all_pass() {
        use noodles::sam::alignment::record::Flags;

        let mut rec1 = create_test_record();
        *rec1.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT;

        let mut rec2 = create_test_record();
        *rec2.flags_mut() = Flags::SEGMENTED | Flags::LAST_SEGMENT;

        let records = vec![rec1, rec2];
        let mut pass_map = AHashMap::new();
        pass_map.insert(0, true);
        pass_map.insert(1, true);

        assert!(template_passes(&records, &pass_map));
    }

    #[test]
    fn test_template_passes_one_fails() {
        use noodles::sam::alignment::record::Flags;

        let mut rec1 = create_test_record();
        *rec1.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT;

        let mut rec2 = create_test_record();
        *rec2.flags_mut() = Flags::SEGMENTED | Flags::LAST_SEGMENT;

        let records = vec![rec1, rec2];
        let mut pass_map = AHashMap::new();
        pass_map.insert(0, true);
        pass_map.insert(1, false); // R2 fails

        assert!(!template_passes(&records, &pass_map));
    }

    #[test]
    fn test_template_passes_secondary_ignored() {
        use noodles::sam::alignment::record::Flags;

        let mut rec1 = create_test_record();
        *rec1.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT;

        let mut rec2 = create_test_record();
        *rec2.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::SECONDARY;

        let records = vec![rec1, rec2];
        let mut pass_map = AHashMap::new();
        pass_map.insert(0, true); // Primary passes
        pass_map.insert(1, false); // Secondary fails (but should be ignored)

        assert!(template_passes(&records, &pass_map));
    }

    #[test]
    fn test_is_duplex_consensus() {
        let mut record = create_test_record();

        // Not duplex initially
        assert!(!is_duplex_consensus(&record));

        // Add aD tag
        let ad_tag = Tag::from([b'a', b'D']);
        record.data_mut().insert(ad_tag, Value::UInt8(3));

        // Now it's duplex
        assert!(is_duplex_consensus(&record));
    }

    // ========== Tests for filter_duplex_read() ==========

    #[test]
    fn test_filter_duplex_read_pass() {
        let mut record = create_test_record();

        // Add final consensus tags (passing)
        record.data_mut().insert(Tag::from([b'c', b'D']), Value::UInt8(10));
        record.data_mut().insert(Tag::from([b'c', b'E']), Value::Float(0.01));

        // Add AB/BA tags (both passing)
        record.data_mut().insert(Tag::from([b'a', b'D']), Value::UInt8(6));
        record.data_mut().insert(Tag::from([b'b', b'D']), Value::UInt8(4));
        record.data_mut().insert(Tag::from([b'a', b'E']), Value::Float(0.01));
        record.data_mut().insert(Tag::from([b'b', b'E']), Value::Float(0.02));

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.25 };

        let result =
            filter_duplex_read(&record, &cc_thresholds, &ab_thresholds, &ba_thresholds).unwrap();
        assert_eq!(result, FilterResult::Pass);
    }

    #[test]
    fn test_filter_duplex_read_insufficient_ab_reads() {
        let mut record = create_test_record();

        // Add final consensus tags (passing)
        record.data_mut().insert(Tag::from([b'c', b'D']), Value::UInt8(10));
        record.data_mut().insert(Tag::from([b'c', b'E']), Value::Float(0.01));

        // Add AB/BA tags - AB has insufficient reads (max depth = 2 < ab_min_reads = 3)
        record.data_mut().insert(Tag::from([b'a', b'D']), Value::UInt8(2));
        record.data_mut().insert(Tag::from([b'b', b'D']), Value::UInt8(2));
        record.data_mut().insert(Tag::from([b'a', b'E']), Value::Float(0.01));
        record.data_mut().insert(Tag::from([b'b', b'E']), Value::Float(0.02));

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.25 };

        let result =
            filter_duplex_read(&record, &cc_thresholds, &ab_thresholds, &ba_thresholds).unwrap();
        assert_eq!(result, FilterResult::InsufficientReads);
    }

    #[test]
    fn test_filter_duplex_read_insufficient_ba_reads() {
        let mut record = create_test_record();

        // Add final consensus tags (passing)
        record.data_mut().insert(Tag::from([b'c', b'D']), Value::UInt8(10));
        record.data_mut().insert(Tag::from([b'c', b'E']), Value::Float(0.01));

        // Add AB/BA tags - BA has insufficient reads (min depth = 1 < ba_min_reads = 2)
        record.data_mut().insert(Tag::from([b'a', b'D']), Value::UInt8(5));
        record.data_mut().insert(Tag::from([b'b', b'D']), Value::UInt8(1));
        record.data_mut().insert(Tag::from([b'a', b'E']), Value::Float(0.01));
        record.data_mut().insert(Tag::from([b'b', b'E']), Value::Float(0.02));

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.25 };

        let result =
            filter_duplex_read(&record, &cc_thresholds, &ab_thresholds, &ba_thresholds).unwrap();
        assert_eq!(result, FilterResult::InsufficientReads);
    }

    #[test]
    fn test_filter_duplex_read_consensus_fails_first() {
        let mut record = create_test_record();

        // Add final consensus tags (failing - insufficient reads)
        record.data_mut().insert(Tag::from([b'c', b'D']), Value::UInt8(3));
        record.data_mut().insert(Tag::from([b'c', b'E']), Value::Float(0.01));

        // Add AB/BA tags (both would pass if checked)
        record.data_mut().insert(Tag::from([b'a', b'D']), Value::UInt8(6));
        record.data_mut().insert(Tag::from([b'b', b'D']), Value::UInt8(4));
        record.data_mut().insert(Tag::from([b'a', b'E']), Value::Float(0.01));
        record.data_mut().insert(Tag::from([b'b', b'E']), Value::Float(0.02));

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.25 };

        let result =
            filter_duplex_read(&record, &cc_thresholds, &ab_thresholds, &ba_thresholds).unwrap();
        // Should fail on consensus check, not get to AB/BA checks
        assert_eq!(result, FilterResult::InsufficientReads);
    }

    #[test]
    fn test_filter_duplex_read_with_only_one_strand() {
        let mut record = create_test_record();

        // Add final consensus tags (passing)
        record.data_mut().insert(Tag::from([b'c', b'D']), Value::UInt8(10));
        record.data_mut().insert(Tag::from([b'c', b'E']), Value::Float(0.01));

        // Add only AB tags (no BA)
        record.data_mut().insert(Tag::from([b'a', b'D']), Value::UInt8(6));
        record.data_mut().insert(Tag::from([b'a', b'E']), Value::Float(0.01));

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.2 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.25 };

        let result =
            filter_duplex_read(&record, &cc_thresholds, &ab_thresholds, &ba_thresholds).unwrap();
        // Should fail because BA has 0 reads (< ba_min_reads = 2)
        assert_eq!(result, FilterResult::InsufficientReads);
    }

    // ========== Tests for mask_duplex_bases() ==========

    fn create_duplex_record_with_per_base_tags() -> RecordBuf {
        let mut record = create_test_record();

        // Add final consensus per-base tags (all bases have good depth and low error)
        // Total depth = 10, total errors = 1, error rate = 1/10 = 0.1
        record.data_mut().insert(
            Tag::from([b'c', b'd']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt16(
                vec![10, 10, 10, 10, 10, 10, 10, 10],
            )),
        );
        record.data_mut().insert(
            Tag::from([b'c', b'e']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![1, 1, 1, 1, 1, 1, 1, 1], // 1 error out of 10 = 0.1 error rate
            )),
        );

        // Add AB per-base tags (higher depth, lower error rate)
        // AB depth = 6, AB errors = 0, error rate = 0/6 = 0.0
        record.data_mut().insert(
            Tag::from([b'a', b'd']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt16(
                vec![6, 6, 6, 6, 6, 6, 6, 6],
            )),
        );
        record.data_mut().insert(
            Tag::from([b'a', b'e']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![0, 0, 0, 0, 0, 0, 0, 0], // 0 errors out of 6 = 0.0 error rate
            )),
        );

        // Add BA per-base tags (lower depth, higher error rate)
        // BA depth = 4, BA errors = 1, error rate = 1/4 = 0.25
        record.data_mut().insert(
            Tag::from([b'b', b'd']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt16(
                vec![4, 4, 4, 4, 4, 4, 4, 4],
            )),
        );
        record.data_mut().insert(
            Tag::from([b'b', b'e']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![1, 1, 1, 1, 1, 1, 1, 1], // 1 error out of 4 = 0.25 error rate
            )),
        );

        record
    }

    #[test]
    fn test_mask_duplex_bases_all_pass() {
        let mut record = create_duplex_record_with_per_base_tags();

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.15 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.25 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.3 };

        mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(13),
            false,
        )
        .unwrap();

        // All bases should pass - no Ns
        let seq = record.sequence();
        assert_eq!(seq.as_ref(), b"ACGTACGT");
    }

    #[test]
    fn test_mask_duplex_bases_low_ab_depth() {
        let mut record = create_duplex_record_with_per_base_tags();

        // Modify AB depth at position 2 to be too low
        record.data_mut().insert(
            Tag::from([b'a', b'd']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt16(
                vec![6, 6, 2, 6, 6, 6, 6, 6], // Position 2 has low AB depth
            )),
        );

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.15 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.25 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.3 };

        mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(13),
            false,
        )
        .unwrap();

        // Position 2 should be masked
        let seq = record.sequence();
        assert_eq!(seq.as_ref(), b"ACNTACGT");
    }

    #[test]
    fn test_mask_duplex_bases_low_ba_depth() {
        let mut record = create_duplex_record_with_per_base_tags();

        // Modify BA depth at positions 3 and 4 to be too low
        record.data_mut().insert(
            Tag::from([b'b', b'd']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt16(
                vec![4, 4, 4, 1, 1, 4, 4, 4], // Positions 3,4 have low BA depth
            )),
        );

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.15 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.25 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.3 };

        mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(13),
            false,
        )
        .unwrap();

        // Positions 3,4 should be masked
        let seq = record.sequence();
        assert_eq!(seq.as_ref(), b"ACGNNCGT");
    }

    #[test]
    fn test_mask_duplex_bases_high_ab_error() {
        let mut record = create_duplex_record_with_per_base_tags();

        // Modify AB error at position 5 to be too high
        // AB depth = 6, set error = 2, error rate = 2/6 = 0.33 > 0.25 threshold
        record.data_mut().insert(
            Tag::from([b'a', b'e']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![0, 0, 0, 0, 0, 2, 0, 0], // Position 5 has 2 errors: 2/6 = 0.33 > 0.25
            )),
        );

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.15 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.25 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.3 };

        mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(13),
            false,
        )
        .unwrap();

        // Position 5 should be masked
        let seq = record.sequence();
        assert_eq!(seq.as_ref(), b"ACGTANGT");
    }

    #[test]
    fn test_mask_duplex_bases_high_ba_error() {
        let mut record = create_duplex_record_with_per_base_tags();

        // Modify BA error at positions 6 and 7 to be too high
        // BA depth = 4, set error = 2, error rate = 2/4 = 0.5 > 0.3 threshold
        record.data_mut().insert(
            Tag::from([b'b', b'e']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![1, 1, 1, 1, 1, 1, 2, 2], // Positions 6,7 have 2 errors: 2/4 = 0.5 > 0.3
            )),
        );

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.15 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.25 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.3 };

        mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(13),
            false,
        )
        .unwrap();

        // Positions 6,7 should be masked
        let seq = record.sequence();
        assert_eq!(seq.as_ref(), b"ACGTACNN");
    }

    #[test]
    fn test_mask_duplex_bases_low_total_depth() {
        let mut record = create_duplex_record_with_per_base_tags();

        // Modify AB and BA depth at position 1 to make total depth too low
        // AB depth = 2, BA depth = 1, total = 3 < 5 threshold
        record.data_mut().insert(
            Tag::from([b'a', b'd']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt16(
                vec![6, 2, 6, 6, 6, 6, 6, 6],
            )),
        );
        record.data_mut().insert(
            Tag::from([b'b', b'd']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt16(
                vec![4, 1, 4, 4, 4, 4, 4, 4],
            )),
        );

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.15 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.25 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.3 };

        mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(13),
            false,
        )
        .unwrap();

        // Position 1 should be masked
        let seq = record.sequence();
        assert_eq!(seq.as_ref(), b"ANGTACGT");
    }

    #[test]
    fn test_mask_duplex_bases_high_total_error() {
        let mut record = create_duplex_record_with_per_base_tags();

        // Modify AB and BA errors at position 0 to make total error rate too high
        // Total depth = 10 (AB=6, BA=4), total errors = 2 (AB=1, BA=1), rate = 2/10 = 0.2 > 0.15
        record.data_mut().insert(
            Tag::from([b'a', b'e']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![1, 0, 0, 0, 0, 0, 0, 0],
            )),
        );
        record.data_mut().insert(
            Tag::from([b'b', b'e']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![1, 1, 1, 1, 1, 1, 1, 1],
            )),
        );

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.15 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.25 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.3 };

        mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(13),
            false,
        )
        .unwrap();

        // Position 0 should be masked
        let seq = record.sequence();
        assert_eq!(seq.as_ref(), b"NCGTACGT");
    }

    #[test]
    fn test_mask_duplex_bases_multiple_failures() {
        let mut record = create_duplex_record_with_per_base_tags();

        // Create multiple failure conditions
        // Position 0: low AB depth (2 < 3)
        record.data_mut().insert(
            Tag::from([b'a', b'd']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt16(
                vec![2, 6, 6, 6, 6, 6, 6, 6],
            )),
        );
        // Position 1: low BA depth (1 < 2)
        record.data_mut().insert(
            Tag::from([b'b', b'd']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt16(
                vec![4, 1, 4, 4, 4, 4, 4, 4],
            )),
        );
        // Position 2: high AB error (2/6 = 0.33 > 0.25)
        record.data_mut().insert(
            Tag::from([b'a', b'e']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![0, 0, 2, 0, 0, 0, 0, 0],
            )),
        );

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.15 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.25 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.3 };

        mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(13),
            false,
        )
        .unwrap();

        // Positions 0,1,2 should be masked
        let seq = record.sequence();
        assert_eq!(seq.as_ref(), b"NNNTACGT");
    }

    // ========== Tests for single-strand agreement ==========

    #[test]
    fn test_mask_duplex_bases_single_strand_agreement_pass() {
        let mut record = create_duplex_record_with_per_base_tags();

        // Add consensus base tags - all agree
        record.data_mut().insert(
            Tag::from([b'a', b'c']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'],
            )),
        );
        record.data_mut().insert(
            Tag::from([b'b', b'c']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'],
            )),
        );

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.15 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.25 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.3 };

        mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(13),
            true,
        )
        .unwrap();

        // All bases agree - no masking
        let seq = record.sequence();
        assert_eq!(seq.as_ref(), b"ACGTACGT");
    }

    #[test]
    fn test_mask_duplex_bases_single_strand_agreement_fail() {
        let mut record = create_duplex_record_with_per_base_tags();

        // Add consensus base tags - positions 2, 5, 7 disagree
        record.data_mut().insert(
            Tag::from([b'a', b'c']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'],
            )),
        );
        record.data_mut().insert(
            Tag::from([b'b', b'c']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![b'A', b'C', b'T', b'T', b'A', b'G', b'G', b'A'], // Disagree at 2,5,7
            )),
        );

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.15 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.25 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.3 };

        mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(13),
            true,
        )
        .unwrap();

        // Positions 2,5,7 should be masked due to disagreement
        let seq = record.sequence();
        assert_eq!(seq.as_ref(), b"ACNTANGN");
    }

    #[test]
    fn test_mask_duplex_bases_single_strand_agreement_disabled() {
        let mut record = create_duplex_record_with_per_base_tags();

        // Add consensus base tags - positions 2, 5, 7 disagree
        record.data_mut().insert(
            Tag::from([b'a', b'c']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'],
            )),
        );
        record.data_mut().insert(
            Tag::from([b'b', b'c']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![b'A', b'C', b'T', b'T', b'A', b'G', b'G', b'A'], // Disagree at 2,5,7
            )),
        );

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.15 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.25 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.3 };

        // Call with require_ss_agreement = false
        mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(13),
            false,
        )
        .unwrap();

        // No masking due to disagreement when disabled
        let seq = record.sequence();
        assert_eq!(seq.as_ref(), b"ACGTACGT");
    }

    #[test]
    fn test_mask_duplex_bases_combined_failures_with_agreement() {
        let mut record = create_duplex_record_with_per_base_tags();

        // Create depth failure at position 0 (low AB depth: 2 < 3)
        record.data_mut().insert(
            Tag::from([b'a', b'd']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt16(
                vec![2, 6, 6, 6, 6, 6, 6, 6],
            )),
        );

        // Create disagreement at positions 3 and 5
        record.data_mut().insert(
            Tag::from([b'a', b'c']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'],
            )),
        );
        record.data_mut().insert(
            Tag::from([b'b', b'c']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(
                vec![b'A', b'C', b'G', b'A', b'A', b'G', b'G', b'T'], // Disagree at 3,5
            )),
        );

        let cc_thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.15 };
        let ab_thresholds =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.25 };
        let ba_thresholds =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.15, max_base_error_rate: 0.3 };

        mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(13),
            true,
        )
        .unwrap();

        // Positions 0 (depth), 3 (disagreement), 5 (disagreement) should be masked
        let seq = record.sequence();
        assert_eq!(seq.as_ref(), b"NCGNANGT");
    }

    #[test]
    fn test_extract_per_base_array_int16() {
        // Test that Int16 arrays are handled correctly (negative values clamped to 0)
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let mut data = noodles::sam::alignment::record_buf::Data::default();
        let tag = Tag::from([b'c', b'D']);
        data.insert(tag, Value::Array(Array::Int16(vec![5, 10, -1, 15, 0])));

        let result = extract_per_base_array(&data, "cD", 5);
        // Negative values should be clamped to 0
        assert_eq!(result, vec![5, 10, 0, 15, 0]);
    }

    #[test]
    fn test_masked_base_quality_phred_2() {
        // Test that masked bases get quality=2 (Phred MIN_VALUE), not 0
        let mut record = RecordBuilder::new().sequence("ACGT").qualities(&[30, 30, 30, 30]).build();

        // Set up low depth at position 0 and 2 (array tags need direct insertion)
        record.data_mut().insert(
            Tag::from([b'c', b'd']),
            Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt16(
                vec![1, 10, 1, 10],
            )),
        );

        let thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 1.0, max_base_error_rate: 1.0 };

        mask_bases(&mut record, &thresholds, Some(10)).unwrap();

        // Positions 0 and 2 should be masked to N with quality=2
        let seq: Vec<u8> = record.sequence().iter().collect();
        let quals: Vec<u8> = record.quality_scores().iter().collect();
        assert_eq!(seq, b"NCNT");
        assert_eq!(
            quals,
            vec![2, 30, 2, 30],
            "Masked bases should have quality=2 (Phred MIN_VALUE)"
        );
    }

    #[test]
    fn test_error_rate_f64_comparison() {
        // Test that error rates are properly compared using f64
        // This ensures the f32 -> f64 promotion works correctly
        let thresholds = FilterThresholds {
            min_reads: 1,
            max_read_error_rate: 0.1,  // f64
            max_base_error_rate: 0.15, // f64
        };

        // Error rate just under the threshold (as f32)
        let error_rate: f32 = 0.099;
        assert!(
            f64::from(error_rate) <= thresholds.max_read_error_rate,
            "f32 -> f64 comparison should work correctly"
        );

        // Error rate just over the threshold (as f32)
        let error_rate_high: f32 = 0.101;
        assert!(
            f64::from(error_rate_high) > thresholds.max_read_error_rate,
            "f32 -> f64 comparison should catch values over threshold"
        );
    }

    #[test]
    fn test_ac_bc_string_tag_handling() {
        // Test that ac/bc tags can be handled as String type (not just Array)
        // This is important for CODEC consensus reads where these may be stored differently
        let mut record = RecordBuilder::new()
            .sequence("ACGT")
            .qualities(&[30, 30, 30, 30])
            .tag("ac", "ACGT")
            .tag("bc", "ACGT")
            .build();

        // These should be handled gracefully even if not array format
        // The function should not panic
        let cc_thresholds =
            FilterThresholds { min_reads: 1, max_read_error_rate: 1.0, max_base_error_rate: 1.0 };
        let ab_thresholds =
            FilterThresholds { min_reads: 1, max_read_error_rate: 1.0, max_base_error_rate: 1.0 };
        let ba_thresholds =
            FilterThresholds { min_reads: 1, max_read_error_rate: 1.0, max_base_error_rate: 1.0 };

        // Should not panic even with String tags
        let result = mask_duplex_bases(
            &mut record,
            &cc_thresholds,
            &ab_thresholds,
            &ba_thresholds,
            Some(10),
            true,
        );
        assert!(result.is_ok());
    }

    // ========== Tests for find_string_or_uint8_array ==========

    #[test]
    fn test_find_string_or_uint8_array_z_tag() {
        // Build aux data with a Z-type string tag: ac:Z:ACGT
        let mut aux = Vec::new();
        aux.extend_from_slice(b"ac"); // tag
        aux.push(b'Z'); // type
        aux.extend_from_slice(b"ACGT\0"); // value + NUL

        let result = super::find_string_or_uint8_array(&aux, *b"ac");
        assert_eq!(result, Some(b"ACGT".to_vec()));
    }

    #[test]
    fn test_find_string_or_uint8_array_b_uint8_tag() {
        // Build aux data with a B-type UInt8 array tag: ac:B:C,65,67,71,84
        let mut aux = Vec::new();
        aux.extend_from_slice(b"ac"); // tag
        aux.push(b'B'); // type = array
        aux.push(b'C'); // sub-type = UInt8
        aux.extend_from_slice(&4u32.to_le_bytes()); // count = 4
        aux.extend_from_slice(&[65u8, 67, 71, 84]); // A, C, G, T

        let result = super::find_string_or_uint8_array(&aux, *b"ac");
        assert_eq!(result, Some(vec![65u8, 67, 71, 84]));
    }

    #[test]
    fn test_find_string_or_uint8_array_missing_tag() {
        let aux: Vec<u8> = Vec::new();
        let result = super::find_string_or_uint8_array(&aux, *b"ac");
        assert!(result.is_none());
    }

    #[test]
    fn test_find_string_or_uint8_array_wrong_array_type() {
        // Build aux data with a B-type Int16 array — should return None since not UInt8
        let mut aux = Vec::new();
        aux.extend_from_slice(b"ac"); // tag
        aux.push(b'B'); // type = array
        aux.push(b's'); // sub-type = Int16
        aux.extend_from_slice(&2u32.to_le_bytes()); // count = 2
        aux.extend_from_slice(&1i16.to_le_bytes());
        aux.extend_from_slice(&2i16.to_le_bytes());

        let result = super::find_string_or_uint8_array(&aux, *b"ac");
        assert!(result.is_none());
    }
}
