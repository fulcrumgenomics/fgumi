//! Consensus read filtering logic.
//!
//! This module provides functionality for filtering consensus reads based on quality, depth,
//! and error rate thresholds. It supports both single-strand and duplex consensus reads.

use ahash::AHashMap;
#[cfg(feature = "simplex")]
use noodles::sam::alignment::record::cigar::op::Kind;

use crate::error::Result;
use crate::phred::{MIN_PHRED, NO_CALL_BASE};
use fgumi_metrics::rejection::RejectionReason;
use fgumi_raw_bam as bam_fields;
use fgumi_raw_bam::{AsTagBytes, RawRecord, RawRecordView, SamTag};

/// Expands a 1-3 element slice to a 3-element array, filling missing values from the last.
///
/// # Panics
/// Panics if `values` is empty.
fn expand_three_from_last<T: Copy>(values: &[T]) -> [T; 3] {
    match values {
        [a, b, c, ..] => [*a, *b, *c],
        [a, b] => [*a, *b, *b],
        [a] => [*a, *a, *a],
        [] => panic!("at least one value required"),
    }
}

/// Filter thresholds for consensus reads
// PartialEq only (not Eq): contains f64 error-rate fields.
#[derive(Debug, Clone, PartialEq)]
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
// PartialEq only (not Eq): contains f64 max_no_call_fraction and
// min_mean_base_quality fields plus nested FilterThresholds with f64 fields.
#[derive(Debug, Clone, PartialEq)]
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

    // NOTE: default-parity between `fgumi filter` (see `commands::filter::Filter::build_filter_config`)
    // and `fgumi runall` (see `commands::runall::parallel::Runall::build_filter_config`, which uses
    // `for_single_strand`) relies on this constructor's exact empty-slice fallbacks below — keep in sync.
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
        let [cc_reads, ab_reads, ba_reads] = expand_three_from_last(min_reads);
        let [cc_read_err, ab_read_err, ba_read_err] = expand_three_from_last(max_read_error_rate);
        let [cc_base_err, ab_base_err, ba_base_err] = expand_three_from_last(max_base_error_rate);

        // Create thresholds for all levels - matching fgbio which always creates all three
        // when filtering either simplex or duplex reads. Single values are replicated to all levels.
        let duplex_thresholds = Some(FilterThresholds {
            min_reads: cc_reads,
            max_read_error_rate: cc_read_err,
            max_base_error_rate: cc_base_err,
        });

        let ab_thresholds = Some(FilterThresholds {
            min_reads: ab_reads,
            max_read_error_rate: ab_read_err,
            max_base_error_rate: ab_base_err,
        });

        let ba_thresholds = Some(FilterThresholds {
            min_reads: ba_reads,
            max_read_error_rate: ba_read_err,
            max_base_error_rate: ba_base_err,
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
/// Checks whether all primary raw records in a template pass their filters.
///
/// A template passes if it has at least one primary read and all primary reads pass.
#[must_use]
pub fn template_passes(raw_records: &[RawRecord], pass_map: &AHashMap<usize, bool>) -> bool {
    let mut has_primary = false;
    let mut all_primary_pass = true;

    for (idx, record) in raw_records.iter().enumerate() {
        let flags = RawRecordView::new(record).flags();
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

/// Outcome of filtering one template's records.
///
/// Returned by [`filter_template_records`]. Callers iterate `keep_mask` to
/// decide which records to emit (or route to rejects), and use `bases_masked`
/// / `pass_map` for metrics.
#[derive(Debug, Clone)]
pub struct TemplateFilterOutcome {
    /// Per-record keep decision (parallel to the input records slice):
    /// `true` means emit to primary output; `false` means reject.
    pub keep_mask: Vec<bool>,
    /// Per-record pass result from the caller-supplied per-record processor.
    pub pass_map: AHashMap<usize, bool>,
    /// Sum of bases masked across all records in this template.
    pub bases_masked: u64,
    /// Whether the template passed the template-level rule: at least one
    /// primary record present and all primary records passed. Always `true`
    /// when `filter_by_template` is `false` (the rule is skipped entirely).
    pub template_pass: bool,
}

/// Apply per-record filtering to a template's records and decide which to emit.
///
/// The single shared orchestration used by both the standalone `fgumi filter`
/// command and the pipeline `FilterStage`. Iterates `records`, invokes
/// `per_record` on each (which masks in place and returns `(bases_masked,
/// pass)`), then applies the template-level rule:
///
/// - When `filter_by_template` is `true`: if *any* primary record fails, drop
///   *all* records of the template. Non-primary records additionally require
///   their own per-record pass.
/// - When `false`: each record is kept iff it individually passed.
///
/// # Errors
///
/// Propagates any error returned by `per_record`.
pub fn filter_template_records<F, E>(
    records: &mut [RawRecord],
    filter_by_template: bool,
    mut per_record: F,
) -> std::result::Result<TemplateFilterOutcome, E>
where
    F: FnMut(&mut RawRecord) -> std::result::Result<(u64, bool), E>,
{
    let mut pass_map: AHashMap<usize, bool> = AHashMap::with_capacity(records.len());
    let mut bases_masked: u64 = 0;
    for (idx, record) in records.iter_mut().enumerate() {
        let (masked, pass) = per_record(record)?;
        bases_masked += masked;
        pass_map.insert(idx, pass);
    }

    let template_pass = if filter_by_template { template_passes(records, &pass_map) } else { true };

    let keep_mask: Vec<bool> = records
        .iter()
        .enumerate()
        .map(|(idx, rec)| {
            let flags = RawRecordView::new(rec).flags();
            let is_primary = (flags & bam_fields::flags::SECONDARY) == 0
                && (flags & bam_fields::flags::SUPPLEMENTARY) == 0;
            let record_pass = pass_map.get(&idx).copied().unwrap_or(false);
            if filter_by_template {
                if is_primary { template_pass } else { template_pass && record_pass }
            } else {
                record_pass
            }
        })
        .collect();

    Ok(TemplateFilterOutcome { keep_mask, pass_map, bases_masked, template_pass })
}

/// Pre-parsed methylation aux tags from a raw BAM record.
///
/// Avoids repeated linear scans of the aux block when multiple filters
/// need the same tag arrays.
pub struct MethylationTags {
    /// Unconverted counts (combined/simplex).
    pub cu: Option<Vec<u16>>,
    /// Converted counts (combined/simplex).
    pub ct: Option<Vec<u16>>,
    /// AB-strand unconverted counts (duplex).
    pub au: Option<Vec<u16>>,
    /// AB-strand converted counts (duplex).
    pub at: Option<Vec<u16>>,
    /// BA-strand unconverted counts (duplex).
    pub bu: Option<Vec<u16>>,
    /// BA-strand converted counts (duplex).
    pub bt: Option<Vec<u16>>,
}

impl MethylationTags {
    /// Parses all methylation tags from a raw BAM record's aux data.
    #[must_use]
    pub fn from_record(record: &[u8]) -> Self {
        use crate::tags::per_base;

        let aux_off = bam_fields::aux_data_offset_from_record(record).unwrap_or(record.len());
        let aux = &record[aux_off..];
        Self {
            cu: Self::find_tag(aux, per_base::UNCONVERTED_COUNT),
            ct: Self::find_tag(aux, per_base::CONVERTED_COUNT),
            au: Self::find_tag(aux, per_base::AB_UNCONVERTED_COUNT),
            at: Self::find_tag(aux, per_base::AB_CONVERTED_COUNT),
            bu: Self::find_tag(aux, per_base::BA_UNCONVERTED_COUNT),
            bt: Self::find_tag(aux, per_base::BA_CONVERTED_COUNT),
        }
    }

    /// Looks up a 2-character tag in the aux data and returns its values as `Vec<u16>`.
    fn find_tag(aux: &[u8], tag: SamTag) -> Option<Vec<u16>> {
        bam_fields::find_array_tag(aux, tag).map(|r| bam_fields::array_tag_to_vec_u16(&r))
    }
}

/// Detects if a raw BAM record is a duplex consensus.
///
/// Checks for presence of `aD` or `bD` tags in the aux data.
#[must_use]
pub fn is_duplex_consensus(aux_data: &[u8]) -> bool {
    bam_fields::find_tag_type(aux_data, SamTag::AD).is_some()
        || bam_fields::find_tag_type(aux_data, SamTag::BD).is_some()
}

/// Filters a raw consensus read based on per-read tags (cD depth, cE error rate).
///
/// # Errors
///
/// Returns an error if the aux data cannot be parsed.
pub fn filter_read(aux_data: &[u8], thresholds: &FilterThresholds) -> Result<FilterResult> {
    // Check minimum reads (cD tag — any integer type)
    if let Some(depth) = bam_fields::find_int_tag(aux_data, SamTag::CD) {
        let min_reads = i64::try_from(thresholds.min_reads).unwrap_or(i64::MAX);
        if depth < min_reads {
            return Ok(FilterResult::InsufficientReads);
        }
    }

    // Check maximum error rate (cE tag — Float)
    if let Some(error_rate) = bam_fields::find_float_tag(aux_data, SamTag::CE) {
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
pub fn filter_duplex_read(
    aux_data: &[u8],
    cc_thresholds: &FilterThresholds,
    ab_thresholds: &FilterThresholds,
    ba_thresholds: &FilterThresholds,
) -> Result<FilterResult> {
    // First check final consensus thresholds
    let result = filter_read(aux_data, cc_thresholds)?;
    if result != FilterResult::Pass {
        return Ok(result);
    }

    // Extract AB and BA depths and error rates
    let ab_depth = bam_fields::find_int_tag(aux_data, SamTag::AD)
        .or_else(|| bam_fields::find_int_tag(aux_data, SamTag::AM));
    let ba_depth = bam_fields::find_int_tag(aux_data, SamTag::BD)
        .or_else(|| bam_fields::find_int_tag(aux_data, SamTag::BM));
    let ab_error = bam_fields::find_float_tag(aux_data, SamTag::AE);
    let ba_error = bam_fields::find_float_tag(aux_data, SamTag::BE);

    // Pick the "best" and "worst" value per metric, independently. `best`/
    // `worst` are per-metric extremes across the two strands — NOT the
    // biological AB/BA strand values. `ab_thresholds` is the stricter tier
    // (checked against the best); `ba_thresholds` is the lenient tier
    // (checked against the worst). Matches fgbio's `abMaxDepth` / `abError`
    // semantics.
    let (worst_depth, best_depth) = match (ab_depth, ba_depth) {
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

    let (best_error, worst_error) = match (ab_error, ba_error) {
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

    // Stricter AB tier: best-per-metric value must clear the threshold.
    #[expect(
        clippy::cast_sign_loss,
        clippy::cast_possible_truncation,
        reason = "depth values are non-negative and fit in usize on all supported platforms"
    )]
    if (best_depth as usize) < ab_thresholds.min_reads {
        return Ok(FilterResult::InsufficientReads);
    }
    if f64::from(best_error) > ab_thresholds.max_read_error_rate {
        return Ok(FilterResult::ExcessiveErrorRate);
    }

    // Lenient BA tier: worst-per-metric value must still clear the threshold.
    #[expect(
        clippy::cast_sign_loss,
        clippy::cast_possible_truncation,
        reason = "depth values are non-negative and fit in usize on all supported platforms"
    )]
    if (worst_depth as usize) < ba_thresholds.min_reads {
        return Ok(FilterResult::InsufficientReads);
    }
    if f64::from(worst_error) > ba_thresholds.max_read_error_rate {
        return Ok(FilterResult::ExcessiveErrorRate);
    }

    Ok(FilterResult::Pass)
}

/// Computes both no-call count and mean base quality in a single pass over raw BAM bytes.
///
/// Returns (`no_call_count`, `mean_base_quality`).
///
/// # Panics
/// Panics if the record is shorter than `MIN_BAM_RECORD_LEN` (36 bytes).
#[must_use]
pub fn compute_read_stats(bam: &[u8]) -> (usize, f64) {
    assert!(bam.len() >= bam_fields::MIN_BAM_RECORD_LEN, "BAM record too short");
    let seq_off = bam_fields::seq_offset(bam);
    let qual_off = bam_fields::qual_offset(bam);
    let len = RawRecordView::new(bam).l_seq() as usize;

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

    #[expect(
        clippy::cast_precision_loss,
        reason = "precision loss is acceptable for quality averaging"
    )]
    let mean_qual = if non_n_count == 0 { 0.0 } else { qual_sum as f64 / non_n_count as f64 };
    (n_count, mean_qual)
}

/// Counts the number of N bases in a raw BAM record.
///
/// This is a thin wrapper around [`compute_read_stats`] for callers that only need the
/// no-call count.  When both the count and mean quality are needed, prefer calling
/// [`compute_read_stats`] directly to avoid a second pass.
///
/// # Panics
/// Panics if the record is shorter than `MIN_BAM_RECORD_LEN` (36 bytes).
#[must_use]
pub fn count_no_calls(bam: &[u8]) -> usize {
    compute_read_stats(bam).0
}

/// Calculates the mean base quality of non-N bases in a raw BAM record.
///
/// This is a thin wrapper around [`compute_read_stats`] for callers that only need the
/// mean quality.  When both values are needed, prefer calling [`compute_read_stats`]
/// directly to avoid a second pass.
///
/// # Panics
/// Panics if the record is shorter than `MIN_BAM_RECORD_LEN` (36 bytes).
#[must_use]
pub fn mean_base_quality(bam: &[u8]) -> f64 {
    compute_read_stats(bam).1
}

/// Reads a tag value as either a Z-type string or a B-type `UInt8` array.
///
/// Returns `Some(Vec<u8>)` if found as either type, `None` otherwise.
fn find_string_or_uint8_array(aux_data: &[u8], tag: impl AsTagBytes) -> Option<Vec<u8>> {
    if let Some(s) = bam_fields::find_string_tag(aux_data, &tag) {
        Some(s.to_vec())
    } else {
        let arr = bam_fields::find_array_tag(aux_data, &tag)?;
        if matches!(arr.elem_type, b'C' | b'c') {
            #[expect(
                clippy::cast_possible_truncation,
                reason = "UInt8 array elements are guaranteed to fit in u8"
            )]
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
#[expect(
    clippy::similar_names,
    reason = "threshold variable names mirror the consensus tag names they check"
)]
pub fn mask_bases(
    record: &mut [u8],
    thresholds: &FilterThresholds,
    min_base_quality: Option<u8>,
) -> Result<usize> {
    if record.len() < bam_fields::MIN_BAM_RECORD_LEN {
        return Err(crate::error::Error::InvalidRecord("BAM record too short".into()));
    }
    let seq_off = bam_fields::seq_offset(record);
    let qual_off = bam_fields::qual_offset(record);
    let len = RawRecordView::new(record).l_seq() as usize;
    let aux_off = bam_fields::aux_data_offset_from_record(record).unwrap_or(record.len());

    // Pre-read per-base arrays into owned Vecs to release the immutable borrow on record
    let cd_vals = bam_fields::find_array_tag(&record[aux_off..], SamTag::CD_BASES)
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));
    let ce_vals = bam_fields::find_array_tag(&record[aux_off..], SamTag::CE_BASES)
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
#[expect(
    clippy::similar_names,
    reason = "threshold variable names mirror the duplex strand tag names"
)]
pub fn mask_duplex_bases(
    record: &mut [u8],
    cc_thresholds: &FilterThresholds,
    ab_thresholds: &FilterThresholds,
    ba_thresholds: &FilterThresholds,
    min_base_quality: Option<u8>,
    require_ss_agreement: bool,
) -> Result<usize> {
    if record.len() < bam_fields::MIN_BAM_RECORD_LEN {
        return Err(crate::error::Error::InvalidRecord("BAM record too short".into()));
    }
    let seq_off = bam_fields::seq_offset(record);
    let qual_off = bam_fields::qual_offset(record);
    let len = RawRecordView::new(record).l_seq() as usize;
    let aux_off = bam_fields::aux_data_offset_from_record(record).unwrap_or(record.len());

    // Pre-read per-base arrays and strings into owned data to release the immutable borrow
    let ad_vals = bam_fields::find_array_tag(&record[aux_off..], SamTag::AD_BASES)
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));
    let ae_vals = bam_fields::find_array_tag(&record[aux_off..], SamTag::AE_BASES)
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));
    let bd_vals = bam_fields::find_array_tag(&record[aux_off..], SamTag::BD_BASES)
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));
    let be_vals = bam_fields::find_array_tag(&record[aux_off..], SamTag::BE_BASES)
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));

    // For single-strand agreement checking, get ac/bc tags (copy to owned).
    // These may be Z-type strings or B-type UInt8 arrays.
    // AC has no SAM-spec clash; BC does, hence BC_BASES for the per-base tag.
    let ac_owned: Option<Vec<u8>> = if require_ss_agreement {
        find_string_or_uint8_array(&record[aux_off..], SamTag::AC)
    } else {
        None
    };
    let bc_owned: Option<Vec<u8>> = if require_ss_agreement {
        find_string_or_uint8_array(&record[aux_off..], SamTag::BC_BASES)
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

        // Best/worst per metric (see `filter_duplex_read_raw` for the tier
        // semantics): AB tier = stricter, checked against best; BA tier =
        // lenient, checked against worst.
        let best_depth = std::cmp::max(ab_depth, ba_depth);
        let worst_depth = std::cmp::min(ab_depth, ba_depth);

        let ab_error_rate =
            if ab_depth > 0 { f64::from(ab_errors) / f64::from(ab_depth) } else { 0.0 };
        let ba_error_rate =
            if ba_depth > 0 { f64::from(ba_errors) / f64::from(ba_depth) } else { 0.0 };

        let best_error_rate = ab_error_rate.min(ba_error_rate);
        let worst_error_rate = ab_error_rate.max(ba_error_rate);

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
            || (best_depth as usize) < ab_thresholds.min_reads
            || best_error_rate > ab_thresholds.max_base_error_rate
            || (worst_depth as usize) < ba_thresholds.min_reads
            || worst_error_rate > ba_thresholds.max_base_error_rate;

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

// ============================================================================
// Methylation (EM-Seq) filter functions
// ============================================================================

/// Thresholds for methylation depth filtering.
///
/// Mirrors the 1-3 value pattern used by other duplex thresholds:
/// `[duplex, AB, BA]` where missing values are filled from the last provided.
#[derive(Debug, Clone)]
pub struct MethylationDepthThresholds {
    /// Minimum combined depth (cu+ct) for the final consensus.
    pub duplex: usize,
    /// Minimum depth (au+at) for the AB strand.
    pub ab: usize,
    /// Minimum depth (bu+bt) for the BA strand.
    pub ba: usize,
}

impl MethylationDepthThresholds {
    /// Creates thresholds from a 1-3 element slice, filling missing values from the last.
    ///
    /// # Panics
    /// Panics if `values` is empty.
    #[must_use]
    pub fn from_values(values: &[usize]) -> Self {
        assert!(!values.is_empty(), "min-methylation-depth must have at least 1 value");
        let [duplex, ab, ba] = expand_three_from_last(values);
        Self { duplex, ab, ba }
    }
}

/// Masks bases in a raw simplex consensus read where methylation depth is too low.
///
/// At each position, checks if `cu[i] + ct[i] < min_depth` and masks if so.
/// Returns the number of newly masked bases.
///
/// # Errors
/// Returns an error if the record is too short.
pub fn mask_methylation_depth_simplex_raw(record: &mut [u8], min_depth: usize) -> Result<usize> {
    let tags = MethylationTags::from_record(record);
    mask_methylation_depth_simplex_raw_with_tags(record, min_depth, &tags)
}

/// Like [`mask_methylation_depth_simplex_raw`] but accepts pre-parsed methylation tags.
///
/// # Errors
/// Returns an error if the record is too short.
pub fn mask_methylation_depth_simplex_raw_with_tags(
    record: &mut [u8],
    min_depth: usize,
    tags: &MethylationTags,
) -> Result<usize> {
    if record.len() < bam_fields::MIN_BAM_RECORD_LEN {
        return Err(crate::error::Error::InvalidRecord("BAM record too short".into()));
    }
    let seq_off = bam_fields::seq_offset(record);
    let qual_off = bam_fields::qual_offset(record);
    let len = RawRecordView::new(record).l_seq() as usize;

    // If no methylation tags present, nothing to filter
    if tags.cu.is_none() && tags.ct.is_none() {
        return Ok(0);
    }

    let mut masked_count = 0;
    for i in 0..len {
        if bam_fields::is_base_n(record, seq_off, i) {
            continue;
        }
        let cu = tags.cu.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let ct = tags.ct.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let total = (cu as usize) + (ct as usize);
        if total < min_depth {
            masked_count += 1;
            bam_fields::mask_base(record, seq_off, i);
            bam_fields::set_qual(record, qual_off, i, MIN_PHRED);
        }
    }
    Ok(masked_count)
}

/// Masks bases in a raw duplex consensus read where methylation depth is too low.
///
/// Checks combined (cu+ct), AB (au+at), and BA (bu+bt) depths against their
/// respective thresholds.
/// Returns the number of newly masked bases.
///
/// # Errors
/// Returns an error if the record is too short.
pub fn mask_methylation_depth_duplex_raw(
    record: &mut [u8],
    thresholds: &MethylationDepthThresholds,
) -> Result<usize> {
    let tags = MethylationTags::from_record(record);
    mask_methylation_depth_duplex_raw_with_tags(record, thresholds, &tags)
}

/// Like [`mask_methylation_depth_duplex_raw`] but accepts pre-parsed methylation tags.
///
/// # Errors
/// Returns an error if the record is too short.
pub fn mask_methylation_depth_duplex_raw_with_tags(
    record: &mut [u8],
    thresholds: &MethylationDepthThresholds,
    tags: &MethylationTags,
) -> Result<usize> {
    if record.len() < bam_fields::MIN_BAM_RECORD_LEN {
        return Err(crate::error::Error::InvalidRecord("BAM record too short".into()));
    }
    let seq_off = bam_fields::seq_offset(record);
    let qual_off = bam_fields::qual_offset(record);
    let len = RawRecordView::new(record).l_seq() as usize;

    // If no methylation tags present, nothing to filter
    if tags.cu.is_none() && tags.ct.is_none() {
        return Ok(0);
    }

    let mut masked_count = 0;
    for i in 0..len {
        if bam_fields::is_base_n(record, seq_off, i) {
            continue;
        }
        let cu = tags.cu.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let ct = tags.ct.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let au = tags.au.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let at = tags.at.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let bu = tags.bu.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let bt = tags.bt.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));

        let cc_total = (cu as usize) + (ct as usize);
        let ab_total = (au as usize) + (at as usize);
        let ba_total = (bu as usize) + (bt as usize);

        if cc_total < thresholds.duplex || ab_total < thresholds.ab || ba_total < thresholds.ba {
            masked_count += 1;
            bam_fields::mask_base(record, seq_off, i);
            bam_fields::set_qual(record, qual_off, i, MIN_PHRED);
        }
    }
    Ok(masked_count)
}

/// Resolves the reference bases for a raw BAM record's alignment region.
///
/// Returns a vector of `Option<u8>` mapping each query position to its reference base
/// (uppercased), or `None` for insertions/soft-clips. Returns `None` if the record is
/// unmapped or the reference cannot be resolved.
#[cfg(feature = "simplex")]
#[expect(
    clippy::cast_sign_loss,
    reason = "ref_id and pos are non-negative for mapped records (checked above)"
)]
pub fn resolve_ref_bases_for_record(
    record: &[u8],
    reference: &dyn crate::methylation::RefBaseProvider,
    ref_names: &[String],
) -> Option<Vec<Option<u8>>> {
    let flags = RawRecordView::new(record).flags();
    if flags & bam_fields::flags::UNMAPPED != 0 {
        return None;
    }

    let tid = RawRecordView::new(record).ref_id();
    if tid < 0 {
        return None;
    }
    let ref_name = ref_names.get(tid as usize)?;
    let alignment_start = RawRecordView::new(record).pos() as u64; // 0-based

    // Try to get the full sequence slice for O(1) indexed access per base,
    // avoiding a HashMap lookup per position.
    let ref_seq = reference.sequence_for(ref_name);

    let cigar_ops = bam_fields::get_cigar_ops(record);
    let len = RawRecordView::new(record).l_seq() as usize;
    let mut result = Vec::with_capacity(len);
    let mut ref_pos = alignment_start;

    for &op in &cigar_ops {
        let op_len = (op >> 4) as usize;
        let kind = bam_fields::cigar_op_kind(op);
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..op_len {
                    let base = if let Some(seq) = ref_seq {
                        usize::try_from(ref_pos)
                            .ok()
                            .and_then(|p| seq.get(p).copied())
                            .map(|b| b.to_ascii_uppercase())
                    } else {
                        reference.base_at_0based(ref_name, ref_pos).map(|b| b.to_ascii_uppercase())
                    };
                    result.push(base);
                    ref_pos += 1;
                }
            }
            Kind::Insertion | Kind::SoftClip => {
                for _ in 0..op_len {
                    result.push(None);
                }
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos += op_len as u64;
            }
            _ => {}
        }
    }

    result.truncate(len);
    while result.len() < len {
        result.push(None);
    }

    Some(result)
}

/// Masks `CpG` positions in a raw duplex consensus where top and bottom strands disagree
/// on methylation status.
///
/// At a `CpG` dinucleotide (ref positions X and `X+1`):
/// - Position X (`ref=C`): top-strand methylation from `au`/`at` tags
/// - Position `X+1` (`ref=G`): bottom-strand methylation from `bu`/`bt` tags (complement)
///
/// If one strand calls methylated and the other calls unmethylated, both positions
/// are masked.
///
/// Returns the number of newly masked bases.
///
/// # Errors
/// Returns an error if the record is too short.
#[cfg(feature = "simplex")]
pub fn mask_strand_methylation_agreement_raw(
    record: &mut [u8],
    reference: &dyn crate::methylation::RefBaseProvider,
    ref_names: &[String],
) -> Result<usize> {
    let ref_base_map = resolve_ref_bases_for_record(record, reference, ref_names);
    let tags = MethylationTags::from_record(record);
    mask_strand_methylation_agreement_raw_with_ref_bases_and_tags(
        record,
        ref_base_map.as_deref(),
        &tags,
    )
}

/// Like `mask_strand_methylation_agreement_raw` but accepts both pre-resolved reference
/// bases and pre-parsed methylation tags, avoiding all redundant work.
///
/// # Errors
/// Returns an error if the record is too short.
pub fn mask_strand_methylation_agreement_raw_with_ref_bases_and_tags(
    record: &mut [u8],
    ref_base_map: Option<&[Option<u8>]>,
    tags: &MethylationTags,
) -> Result<usize> {
    if record.len() < bam_fields::MIN_BAM_RECORD_LEN {
        return Err(crate::error::Error::InvalidRecord("BAM record too short".into()));
    }

    let Some(ref_base_map) = ref_base_map else {
        return Ok(0);
    };

    let seq_off = bam_fields::seq_offset(record);
    let qual_off = bam_fields::qual_offset(record);
    let len = RawRecordView::new(record).l_seq() as usize;

    // If no per-strand methylation tags, nothing to check
    if tags.au.is_none() && tags.bu.is_none() {
        return Ok(0);
    }

    // Identify CpG dinucleotides and check strand agreement.
    // A CpG is where ref_base[i] = C and ref_base[i+1] = G,
    // and both positions are aligned (Some).
    let mut should_mask = vec![false; len];

    for i in 0..len.saturating_sub(1) {
        let (Some(Some(b'C')), Some(Some(b'G'))) = (ref_base_map.get(i), ref_base_map.get(i + 1))
        else {
            continue;
        };

        // Position i: top-strand C -> check au/at for methylation
        let top_unconverted = tags.au.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let top_converted = tags.at.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));

        // Position i+1: bottom-strand G (complement of C) -> check bu/bt for methylation
        let bot_unconverted = tags.bu.as_ref().map_or(0u16, |v| v.get(i + 1).copied().unwrap_or(0));
        let bot_converted = tags.bt.as_ref().map_or(0u16, |v| v.get(i + 1).copied().unwrap_or(0));

        let top_total = u32::from(top_unconverted) + u32::from(top_converted);
        let bot_total = u32::from(bot_unconverted) + u32::from(bot_converted);

        // Skip if either strand has no evidence
        if top_total == 0 || bot_total == 0 {
            continue;
        }

        // Call methylated if unconverted > converted (majority rule)
        let top_methylated = top_unconverted > top_converted;
        let bot_methylated = bot_unconverted > bot_converted;

        if top_methylated != bot_methylated {
            should_mask[i] = true;
            should_mask[i + 1] = true;
        }
    }

    let mut masked_count = 0;
    for (i, &mask) in should_mask.iter().enumerate() {
        if mask && !bam_fields::is_base_n(record, seq_off, i) {
            masked_count += 1;
            bam_fields::mask_base(record, seq_off, i);
            bam_fields::set_qual(record, qual_off, i, MIN_PHRED);
        }
    }

    Ok(masked_count)
}

/// Checks the bisulfite/enzymatic conversion fraction at non-CpG cytosines.
///
/// Returns `true` if the read passes (conversion fraction >= threshold), or if
/// there are no non-CpG cytosine positions to evaluate.
///
/// At non-CpG ref-C positions, the expected behavior is conversion (C->T).
/// A low conversion rate suggests incomplete enzymatic conversion.
#[cfg(feature = "simplex")]
pub fn check_conversion_fraction_raw(
    record: &[u8],
    min_fraction: f64,
    reference: &dyn crate::methylation::RefBaseProvider,
    ref_names: &[String],
    methylation_mode: crate::MethylationMode,
) -> bool {
    let ref_base_map = resolve_ref_bases_for_record(record, reference, ref_names);
    let meth_tags = MethylationTags::from_record(record);
    check_conversion_fraction_raw_with_ref_bases_and_tags(
        record,
        min_fraction,
        ref_base_map.as_deref(),
        &meth_tags,
        methylation_mode,
    )
}

/// Like `check_conversion_fraction_raw` but accepts both pre-resolved reference bases
/// and pre-parsed methylation tags, avoiding all redundant work.
///
/// For EM-Seq, checks `ct / (cu + ct) >= threshold` at non-CpG ref-C positions.
/// For TAPs, checks `cu / (cu + ct) >= threshold` instead, since
/// non-CpG Cs should NOT be converted in TAPs.
#[must_use]
pub fn check_conversion_fraction_raw_with_ref_bases_and_tags(
    record: &[u8],
    min_fraction: f64,
    ref_base_map: Option<&[Option<u8>]>,
    methylation_tags: &MethylationTags,
    methylation_mode: crate::MethylationMode,
) -> bool {
    // Conversion fraction is meaningless without a methylation mode — pass through
    if methylation_mode == crate::MethylationMode::Disabled {
        return true;
    }

    let Some(ref_base_map) = ref_base_map else {
        return true; // unmapped reads pass
    };

    let len = RawRecordView::new(record).l_seq() as usize;

    // If no methylation tags, pass
    if methylation_tags.cu.is_none() && methylation_tags.ct.is_none() {
        return true;
    }

    let mut total_numerator: u64 = 0;
    let mut total_evidence: u64 = 0;

    for i in 0..len {
        // Only consider non-CpG ref-C positions
        let Some(Some(b'C')) = ref_base_map.get(i) else {
            continue;
        };

        // Check if this is a CpG (next base is G) -- skip CpG sites
        let is_cpg = i + 1 < len && matches!(ref_base_map.get(i + 1), Some(Some(b'G')));
        if is_cpg {
            continue;
        }

        let cu = methylation_tags.cu.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let ct = methylation_tags.ct.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let evidence = u64::from(cu) + u64::from(ct);
        if evidence > 0 {
            // For EM-Seq: count converted (ct) — high conversion = good library quality
            // For TAPs: count unconverted (cu) — high non-conversion at non-CpG = good specificity
            // Safety: Disabled is handled by the early return at the top of this function
            let numerator = match methylation_mode {
                crate::MethylationMode::Taps => u64::from(cu),
                crate::MethylationMode::EmSeq => u64::from(ct),
                crate::MethylationMode::Disabled => return true,
            };
            total_numerator += numerator;
            total_evidence += evidence;
        }
    }

    // If no non-CpG C positions with evidence, pass
    if total_evidence == 0 {
        return true;
    }

    #[expect(
        clippy::cast_precision_loss,
        reason = "precision loss is acceptable for fraction calculation"
    )]
    let fraction = total_numerator as f64 / total_evidence as f64;
    fraction >= min_fraction
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::SamBuilder as RawSamBuilder;

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

        let thresh = config.single_strand_thresholds.expect("failed to get thresh");
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

        let duplex = config.duplex_thresholds.expect("failed to get duplex");
        assert_eq!(duplex.min_reads, 5);
        assert!((duplex.max_read_error_rate - 0.05).abs() < f64::EPSILON);

        let ab = config.ab_thresholds.expect("failed to get ab");
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

        let cc = config.duplex_thresholds.expect("failed to get cc");
        let ab = config.ab_thresholds.expect("failed to get ab");
        let ba = config.ba_thresholds.expect("failed to get ba");

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
        let ss = config.single_strand_thresholds.expect("failed to get ss");
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

        let cc = config.duplex_thresholds.expect("failed to get cc");
        let ab = config.ab_thresholds.expect("failed to get ab");
        let ba = config.ba_thresholds.expect("failed to get ba");

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

        let cc = config.duplex_thresholds.expect("failed to get cc");
        let ab_t = config.ab_thresholds.expect("failed to get ab_t");
        let ba_t = config.ba_thresholds.expect("failed to get ba_t");

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

    // ========== Tests for find_string_or_uint8_array ==========

    #[test]
    fn test_find_string_or_uint8_array_z_tag() {
        // Build aux data with a Z-type string tag: ac:Z:ACGT
        let mut aux = Vec::new();
        aux.extend_from_slice(SamTag::AC.as_ref()); // tag
        aux.push(b'Z'); // type
        aux.extend_from_slice(b"ACGT\0"); // value + NUL

        let result = super::find_string_or_uint8_array(&aux, SamTag::AC);
        assert_eq!(result, Some(b"ACGT".to_vec()));
    }

    #[test]
    fn test_find_string_or_uint8_array_b_uint8_tag() {
        // Build aux data with a B-type UInt8 array tag: ac:B:C,65,67,71,84
        let mut aux = Vec::new();
        aux.extend_from_slice(SamTag::AC.as_ref()); // tag
        aux.push(b'B'); // type = array
        aux.push(b'C'); // sub-type = UInt8
        aux.extend_from_slice(&4u32.to_le_bytes()); // count = 4
        aux.extend_from_slice(&[65u8, 67, 71, 84]); // A, C, G, T

        let result = super::find_string_or_uint8_array(&aux, SamTag::AC);
        assert_eq!(result, Some(vec![65u8, 67, 71, 84]));
    }

    #[test]
    fn test_find_string_or_uint8_array_missing_tag() {
        let aux: Vec<u8> = Vec::new();
        let result = super::find_string_or_uint8_array(&aux, SamTag::AC);
        assert!(result.is_none());
    }

    #[test]
    fn test_find_string_or_uint8_array_wrong_array_type() {
        // Build aux data with a B-type Int16 array — should return None since not UInt8
        let mut aux = Vec::new();
        aux.extend_from_slice(SamTag::AC.as_ref()); // tag
        aux.push(b'B'); // type = array
        aux.push(b's'); // sub-type = Int16
        aux.extend_from_slice(&2u32.to_le_bytes()); // count = 2
        aux.extend_from_slice(&1i16.to_le_bytes());
        aux.extend_from_slice(&2i16.to_le_bytes());

        let result = super::find_string_or_uint8_array(&aux, SamTag::AC);
        assert!(result.is_none());
    }

    // ========================================================================
    // Methylation filter tests
    // ========================================================================

    // -- MethylationDepthThresholds tests --

    #[test]
    fn test_methylation_depth_thresholds_single_value() {
        let t = MethylationDepthThresholds::from_values(&[5]);
        assert_eq!(t.duplex, 5);
        assert_eq!(t.ab, 5);
        assert_eq!(t.ba, 5);
    }

    #[test]
    fn test_methylation_depth_thresholds_two_values() {
        let t = MethylationDepthThresholds::from_values(&[10, 3]);
        assert_eq!(t.duplex, 10);
        assert_eq!(t.ab, 3);
        assert_eq!(t.ba, 3);
    }

    #[test]
    fn test_methylation_depth_thresholds_three_values() {
        let t = MethylationDepthThresholds::from_values(&[10, 5, 2]);
        assert_eq!(t.duplex, 10);
        assert_eq!(t.ab, 5);
        assert_eq!(t.ba, 2);
    }

    // -- mask_methylation_depth_simplex_raw tests --

    #[test]
    fn test_mask_methylation_depth_simplex_all_pass() {
        // cu+ct = 5+3=8 at each position, min_depth=5 -> all pass
        let mut raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0).pos(0).mapq(60).cigar_ops(&[4 << 4]).sequence(b"ACGT").qualities(&[30; 4]);
            b.add_array_i16(SamTag::CU, &[5, 5, 5, 5]).add_array_i16(SamTag::CT, &[3, 3, 3, 3]);
            b.build()
        };
        let masked = mask_methylation_depth_simplex_raw(&mut raw, 5).unwrap();
        assert_eq!(masked, 0);
    }

    #[test]
    fn test_mask_methylation_depth_simplex_some_fail() {
        // Position 0: cu+ct = 5+3=8 -> pass
        // Position 1: cu+ct = 1+1=2 -> fail (< 5)
        // Position 2: cu+ct = 0+0=0 -> fail
        // Position 3: cu+ct = 10+0=10 -> pass
        let mut raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0).pos(0).mapq(60).cigar_ops(&[4 << 4]).sequence(b"ACGT").qualities(&[30; 4]);
            b.add_array_i16(SamTag::CU, &[5, 1, 0, 10]).add_array_i16(SamTag::CT, &[3, 1, 0, 0]);
            b.build()
        };
        let masked = mask_methylation_depth_simplex_raw(&mut raw, 5).unwrap();
        assert_eq!(masked, 2, "Positions 1 and 2 should be masked");
    }

    #[test]
    fn test_mask_methylation_depth_simplex_no_tags_no_masking() {
        // No cu/ct tags
        let mut raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0).pos(0).mapq(60).cigar_ops(&[4 << 4]).sequence(b"ACGT").qualities(&[30; 4]);
            b.build()
        };
        let masked = mask_methylation_depth_simplex_raw(&mut raw, 5).unwrap();
        assert_eq!(masked, 0, "No methylation tags should mean no masking");
    }

    // -- mask_methylation_depth_duplex_raw tests --

    #[test]
    fn test_mask_methylation_depth_duplex_all_pass() {
        let mut raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0).pos(0).mapq(60).cigar_ops(&[4 << 4]).sequence(b"ACGT").qualities(&[30; 4]);
            b.add_array_i16(SamTag::CU, &[10, 10, 10, 10])
                .add_array_i16(SamTag::CT, &[2, 2, 2, 2])
                .add_array_i16(SamTag::AU, &[5, 5, 5, 5])
                .add_array_i16(SamTag::AT, &[1, 1, 1, 1])
                .add_array_i16(SamTag::BU, &[5, 5, 5, 5])
                .add_array_i16(SamTag::BT, &[1, 1, 1, 1]);
            b.build()
        };
        let thresholds = MethylationDepthThresholds { duplex: 5, ab: 3, ba: 3 };
        let masked = mask_methylation_depth_duplex_raw(&mut raw, &thresholds).unwrap();
        assert_eq!(masked, 0);
    }

    #[test]
    fn test_mask_methylation_depth_duplex_ab_fails() {
        // Duplex passes (cu+ct=12), but AB fails at position 1 (au+at=1 < 3)
        let mut raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0).pos(0).mapq(60).cigar_ops(&[4 << 4]).sequence(b"ACGT").qualities(&[30; 4]);
            b.add_array_i16(SamTag::CU, &[10, 10, 10, 10])
                .add_array_i16(SamTag::CT, &[2, 2, 2, 2])
                .add_array_i16(SamTag::AU, &[5, 0, 5, 5])
                .add_array_i16(SamTag::AT, &[1, 1, 1, 1])
                .add_array_i16(SamTag::BU, &[5, 5, 5, 5])
                .add_array_i16(SamTag::BT, &[1, 1, 1, 1]);
            b.build()
        };
        let thresholds = MethylationDepthThresholds { duplex: 5, ab: 3, ba: 3 };
        let masked = mask_methylation_depth_duplex_raw(&mut raw, &thresholds).unwrap();
        assert_eq!(masked, 1, "Position 1 should be masked (AB depth too low)");
    }

    // -- mask_strand_methylation_agreement_raw tests --

    #[cfg(feature = "simplex")]
    #[test]
    fn test_strand_methylation_agreement_concordant() {
        use crate::methylation::tests::TestRef;
        // Reference: ...CG... at positions 5,6
        // Both strands agree: methylated at CpG
        let ref_seq = b"AAAAACGAAAA";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        // au/at at position 5 (C): methylated (au=10, at=1)
        // bu/bt at position 6 (G): methylated (bu=10, bt=1)
        let mut raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[11 << 4])
                .sequence(b"AAAAACGAAAA")
                .qualities(&[30; 11]);
            b.add_array_i16(SamTag::AU, &[0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0])
                .add_array_i16(SamTag::AT, &[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0])
                .add_array_i16(SamTag::BU, &[0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0])
                .add_array_i16(SamTag::BT, &[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]);
            b.build()
        };
        let masked =
            mask_strand_methylation_agreement_raw(&mut raw, &reference, &ref_names).unwrap();
        assert_eq!(masked, 0, "Concordant CpG should not be masked");
    }

    #[cfg(feature = "simplex")]
    #[test]
    fn test_strand_methylation_agreement_discordant() {
        use crate::methylation::tests::TestRef;
        // Reference: ...CG... at positions 5,6
        // Top strand says methylated, bottom says unmethylated
        let ref_seq = b"AAAAACGAAAA";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        // au/at at position 5 (C): methylated (au=10, at=1)
        // bu/bt at position 6 (G): unmethylated (bu=1, bt=10)
        let mut raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[11 << 4])
                .sequence(b"AAAAACGAAAA")
                .qualities(&[30; 11]);
            b.add_array_i16(SamTag::AU, &[0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0])
                .add_array_i16(SamTag::AT, &[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0])
                .add_array_i16(SamTag::BU, &[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0])
                .add_array_i16(SamTag::BT, &[0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0]);
            b.build()
        };
        let masked =
            mask_strand_methylation_agreement_raw(&mut raw, &reference, &ref_names).unwrap();
        assert_eq!(masked, 2, "Both CpG positions should be masked when strands disagree");
    }

    #[cfg(feature = "simplex")]
    #[test]
    fn test_strand_methylation_agreement_non_cpg_ignored() {
        use crate::methylation::tests::TestRef;
        // Reference: no CpG sites, just scattered C and G
        let ref_seq = b"ACAGTGCATGA";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        let mut raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[11 << 4])
                .sequence(b"ACAGTGCATGA")
                .qualities(&[30; 11]);
            b.add_array_i16(SamTag::AU, &[0; 11])
                .add_array_i16(SamTag::AT, &[0; 11])
                .add_array_i16(SamTag::BU, &[0; 11])
                .add_array_i16(SamTag::BT, &[0; 11]);
            b.build()
        };
        let masked =
            mask_strand_methylation_agreement_raw(&mut raw, &reference, &ref_names).unwrap();
        assert_eq!(masked, 0, "Non-CpG sites should not be masked");
    }

    // -- check_conversion_fraction_raw tests --

    #[cfg(feature = "simplex")]
    #[test]
    fn test_conversion_fraction_passes_high_conversion() {
        use crate::methylation::tests::TestRef;
        // Reference: ACATACATA (non-CpG C at positions 1, 5 — each C followed by A, not CpG)
        //            012345678
        let ref_seq = b"ACATACATA";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        // Non-CpG C at positions 1 and 5: high conversion (ct >> cu)
        // cu: unconverted count, ct: converted count
        let raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[9 << 4])
                .sequence(b"ACATACATA")
                .qualities(&[30; 9]);
            b.add_array_i16(SamTag::CU, &[0, 1, 0, 0, 0, 1, 0, 0, 0])
                .add_array_i16(SamTag::CT, &[0, 9, 0, 0, 0, 9, 0, 0, 0]);
            b.build()
        };
        // 18 converted out of 20 = 90% conversion (positions 1 and 5)
        assert!(
            check_conversion_fraction_raw(
                &raw,
                0.8,
                &reference,
                &ref_names,
                crate::MethylationMode::EmSeq
            ),
            "Should pass with high conversion fraction"
        );
    }

    #[cfg(feature = "simplex")]
    #[test]
    fn test_conversion_fraction_fails_low_conversion() {
        use crate::methylation::tests::TestRef;
        let ref_seq = b"ACATACATA";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        // Non-CpG C at positions 1 and 5: low conversion (cu >> ct)
        let raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[9 << 4])
                .sequence(b"ACATACATA")
                .qualities(&[30; 9]);
            b.add_array_i16(SamTag::CU, &[0, 9, 0, 0, 0, 9, 0, 0, 0])
                .add_array_i16(SamTag::CT, &[0, 1, 0, 0, 0, 1, 0, 0, 0]);
            b.build()
        };
        // 2 converted out of 20 = 10% conversion (positions 1 and 5)
        assert!(
            !check_conversion_fraction_raw(
                &raw,
                0.8,
                &reference,
                &ref_names,
                crate::MethylationMode::EmSeq
            ),
            "Should fail with low conversion fraction"
        );
    }

    #[cfg(feature = "simplex")]
    #[test]
    fn test_conversion_fraction_skips_cpg() {
        use crate::methylation::tests::TestRef;
        // Reference with CpG at position 4-5: AAAA CG AAA
        let ref_seq = b"AAAACGAAA";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        // CpG C at position 4: high unconverted (methylated) -- should be EXCLUDED
        // No non-CpG C positions -> should pass vacuously
        let raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[9 << 4])
                .sequence(b"AAAACGAAA")
                .qualities(&[30; 9]);
            b.add_array_i16(SamTag::CU, &[0, 0, 0, 0, 10, 0, 0, 0, 0])
                .add_array_i16(SamTag::CT, &[0, 0, 0, 0, 0, 0, 0, 0, 0]);
            b.build()
        };
        assert!(
            check_conversion_fraction_raw(
                &raw,
                0.9,
                &reference,
                &ref_names,
                crate::MethylationMode::EmSeq
            ),
            "Should pass because CpG sites are excluded from conversion check"
        );
    }

    #[cfg(feature = "simplex")]
    #[test]
    fn test_conversion_fraction_no_methylation_tags_passes() {
        use crate::methylation::tests::TestRef;
        let ref_seq = b"ACATACATA";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        // No cu/ct tags
        let raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[9 << 4])
                .sequence(b"ACATACATA")
                .qualities(&[30; 9]);
            b.build()
        };
        assert!(
            check_conversion_fraction_raw(
                &raw,
                0.9,
                &reference,
                &ref_names,
                crate::MethylationMode::EmSeq
            ),
            "Should pass with no methylation tags"
        );
    }

    #[cfg(feature = "simplex")]
    #[test]
    fn test_conversion_fraction_unmapped_passes() {
        use crate::methylation::tests::TestRef;
        let ref_seq = b"ACATACATA";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        // Unmapped record (no ref_id, no cigar)
        let raw = {
            let mut b = RawSamBuilder::new();
            b.sequence(b"ACATACATA").qualities(&[30; 9]);
            b.build()
        };
        assert!(
            check_conversion_fraction_raw(
                &raw,
                0.9,
                &reference,
                &ref_names,
                crate::MethylationMode::EmSeq
            ),
            "Unmapped reads should pass"
        );
    }

    // -- TAPs conversion fraction tests --

    #[cfg(feature = "simplex")]
    #[test]
    fn test_conversion_fraction_taps_passes_high_non_conversion() {
        use crate::methylation::tests::TestRef;
        let ref_seq = b"ACATACATA";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        // Non-CpG C at positions 1 and 5: high unconverted (cu >> ct) = good TAPs specificity
        let raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[9 << 4])
                .sequence(b"ACATACATA")
                .qualities(&[30; 9]);
            b.add_array_i16(SamTag::CU, &[0, 9, 0, 0, 0, 9, 0, 0, 0])
                .add_array_i16(SamTag::CT, &[0, 1, 0, 0, 0, 1, 0, 0, 0]);
            b.build()
        };
        // TAPs numerator = cu: 18 out of 20 = 90%
        assert!(
            check_conversion_fraction_raw(
                &raw,
                0.8,
                &reference,
                &ref_names,
                crate::MethylationMode::Taps
            ),
            "TAPs should pass with high non-conversion fraction"
        );
    }

    #[cfg(feature = "simplex")]
    #[test]
    fn test_conversion_fraction_taps_fails_low_non_conversion() {
        use crate::methylation::tests::TestRef;
        let ref_seq = b"ACATACATA";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        // Non-CpG C at positions 1 and 5: low unconverted (ct >> cu) = poor TAPs specificity
        let raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[9 << 4])
                .sequence(b"ACATACATA")
                .qualities(&[30; 9]);
            b.add_array_i16(SamTag::CU, &[0, 1, 0, 0, 0, 1, 0, 0, 0])
                .add_array_i16(SamTag::CT, &[0, 9, 0, 0, 0, 9, 0, 0, 0]);
            b.build()
        };
        // TAPs numerator = cu: 2 out of 20 = 10%
        assert!(
            !check_conversion_fraction_raw(
                &raw,
                0.8,
                &reference,
                &ref_names,
                crate::MethylationMode::Taps
            ),
            "TAPs should fail with low non-conversion fraction"
        );
    }

    #[cfg(feature = "simplex")]
    #[test]
    fn test_conversion_fraction_taps_vs_emseq_inverted() {
        use crate::methylation::tests::TestRef;
        let ref_seq = b"ACATACATA";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        // Non-CpG C at positions 1 and 5: cu=9, ct=1 at each
        // EM-Seq numerator = ct: 2/20 = 10%, TAPs numerator = cu: 18/20 = 90%
        let raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[9 << 4])
                .sequence(b"ACATACATA")
                .qualities(&[30; 9]);
            b.add_array_i16(SamTag::CU, &[0, 9, 0, 0, 0, 9, 0, 0, 0])
                .add_array_i16(SamTag::CT, &[0, 1, 0, 0, 0, 1, 0, 0, 0]);
            b.build()
        };
        // EM-Seq should fail (10% < 80%), TAPs should pass (90% >= 80%)
        assert!(
            !check_conversion_fraction_raw(
                &raw,
                0.8,
                &reference,
                &ref_names,
                crate::MethylationMode::EmSeq
            ),
            "EM-Seq should fail with high unconverted fraction"
        );
        assert!(
            check_conversion_fraction_raw(
                &raw,
                0.8,
                &reference,
                &ref_names,
                crate::MethylationMode::Taps
            ),
            "TAPs should pass with same data (inverted numerator)"
        );
    }

    // -- Disabled mode tests --

    #[cfg(feature = "simplex")]
    #[test]
    fn test_conversion_fraction_disabled_mode_passes() {
        use crate::methylation::tests::TestRef;
        let ref_seq = b"ACATACATA";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        // Data that would fail EM-Seq at 80% threshold
        let raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[9 << 4])
                .sequence(b"ACATACATA")
                .qualities(&[30; 9]);
            b.add_array_i16(SamTag::CU, &[0, 9, 0, 0, 0, 9, 0, 0, 0])
                .add_array_i16(SamTag::CT, &[0, 1, 0, 0, 0, 1, 0, 0, 0]);
            b.build()
        };
        // Disabled mode should always pass regardless of data
        assert!(
            check_conversion_fraction_raw(
                &raw,
                0.8,
                &reference,
                &ref_names,
                crate::MethylationMode::Disabled
            ),
            "Disabled mode should always pass conversion fraction check"
        );
    }

    // -- resolve_ref_bases_for_record tests --

    #[cfg(feature = "simplex")]
    #[test]
    fn test_resolve_ref_bases_simple_match() {
        use crate::methylation::tests::TestRef;
        let ref_seq = b"ACGTACGTAC";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        let raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0).pos(0).mapq(60).cigar_ops(&[4 << 4]).sequence(b"ACGT").qualities(&[30; 4]);
            b.build()
        };
        let bases = resolve_ref_bases_for_record(&raw, &reference, &ref_names).unwrap();
        assert_eq!(bases, vec![Some(b'A'), Some(b'C'), Some(b'G'), Some(b'T')]);
    }

    #[cfg(feature = "simplex")]
    #[test]
    fn test_resolve_ref_bases_with_insertion() {
        use crate::methylation::tests::TestRef;
        let ref_seq = b"ACGTACGTAC";
        let reference = TestRef::new(&[("chr1", ref_seq)]);
        let ref_names = vec!["chr1".to_string()];

        let raw = {
            let mut b = RawSamBuilder::new();
            // 2M2I2M -> positions 2,3 are insertions
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[2 << 4, 2 << 4 | 1, 2 << 4])
                .sequence(b"ACNNGT")
                .qualities(&[30; 6]);
            b.build()
        };
        let bases = resolve_ref_bases_for_record(&raw, &reference, &ref_names).unwrap();
        assert_eq!(bases, vec![Some(b'A'), Some(b'C'), None, None, Some(b'G'), Some(b'T')]);
    }
}
