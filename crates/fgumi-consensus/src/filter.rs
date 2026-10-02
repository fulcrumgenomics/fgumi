//! Consensus read filtering logic.
//!
//! This module provides functionality for filtering consensus reads based on quality, depth,
//! and error rate thresholds. It supports both single-strand and duplex consensus reads.

use ahash::AHashMap;
use anyhow::Result;
#[cfg(feature = "simplex")]
use noodles::sam::alignment::record::cigar::op::Kind;

use crate::phred::{MIN_PHRED, NO_CALL_BASE};
use fgumi_metrics::rejection::RejectionReason;
use fgumi_raw_bam as bam_fields;
use fgumi_raw_bam::{AsTagBytes, RawRecord, RawRecordView, SamTag};

pub use crate::modifications::{
    drop_masked_modifications_raw, drop_modifications_at_raw, has_modification_tags,
};

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

    /// Minimum mean base quality over the full read length, computed prior to masking (optional)
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

/// Sums the masked-base counts that fgbio's `maskedBases` statistic includes: only the
/// **primary** reads of a **retained** template contribute; every read of a dropped
/// template and every secondary/supplementary read contributes zero.
///
/// This mirrors `FilterConsensusReads.scala:207-219`, where fgbio accumulates
/// `maskedBases += r1Result.maskedBases + r2Result.maskedBases` **only** inside
/// `if (r1Result.keepRead && r2Result.keepRead)` and **only** over the primary R1/R2 —
/// supplementary/secondary reads are filtered and written out, but their masked bases are
/// never added to the tally ("Masked X of Y bases in retained primary consensus reads").
///
/// `masked_by_record[i]` is the number of bases masked in `raw_records[i]`; the slices are
/// parallel and must be the same length. `template_pass` is the result of
/// [`template_passes`] for this template (for the per-record streaming path, a single-read
/// "template" whose pass value is the record's own filter result).
///
/// # Panics
///
/// Panics if `raw_records` and `masked_by_record` differ in length. A release-mode-silent
/// truncation here would under/over-count the exact "Total bases masked" statistic this
/// helper exists to compute, so the invariant is enforced unconditionally rather than only
/// under `debug_assert`.
#[must_use]
pub fn retained_primary_masked_bases(
    raw_records: &[RawRecord],
    masked_by_record: &[u64],
    template_pass: bool,
) -> u64 {
    assert_eq!(
        raw_records.len(),
        masked_by_record.len(),
        "raw_records and masked_by_record must be parallel"
    );
    if !template_pass {
        return 0;
    }
    raw_records
        .iter()
        .zip(masked_by_record)
        .filter(|(record, _)| {
            let flags = RawRecordView::new(record).flags();
            (flags & bam_fields::flags::SECONDARY) == 0
                && (flags & bam_fields::flags::SUPPLEMENTARY) == 0
        })
        .map(|(_, &masked)| masked)
        .sum()
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
    /// Whether the AB-strand counts (`au`) are present (duplex with an AB strand).
    pub has_ab: bool,
    /// Whether the BA-strand counts (`bu`) are present (duplex with a BA strand).
    pub has_ba: bool,
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
            has_ab: bam_fields::find_tag_type(aux, per_base::AB_UNCONVERTED_COUNT).is_some(),
            has_ba: bam_fields::find_tag_type(aux, per_base::BA_UNCONVERTED_COUNT).is_some(),
        }
    }

    /// Whether `cu` and `ct` (those present) have one entry per base of a SEQ of length
    /// `l_seq`. They do not after SEQ was hard-clipped once they were written.
    #[must_use]
    pub fn fits(&self, l_seq: usize) -> bool {
        [&self.cu, &self.ct].into_iter().flatten().all(|counts| counts.len() == l_seq)
    }

    /// Looks up a 2-character tag in the aux data and returns its values as `Vec<u16>`.
    fn find_tag(aux: &[u8], tag: SamTag) -> Option<Vec<u16>> {
        bam_fields::find_array_tag(aux, tag).map(|r| bam_fields::array_tag_to_vec_u16(&r))
    }
}

/// The NM/UQ scoring that matches a consensus record's SEQ convention.
///
/// [`ConversionScoring::Hidden`](fgumi_sam::alignment_tags::ConversionScoring::Hidden) for a simplex methylation consensus: it carries `cu` (the
/// per-base methylation counts) and is not a duplex consensus ([`is_duplex_consensus`]), and its
/// SEQ keeps the converted bases, as a bisulfite-aware aligner scored them.
/// [`ConversionScoring::Literal`](fgumi_sam::alignment_tags::ConversionScoring::Literal) otherwise, including a duplex consensus, whose SEQ is the
/// molecule's sequence, and any non-methylation record. `filter` and `clip` both use it, so NM
/// means the same whichever of them wrote it.
#[must_use]
pub fn conversion_scoring_for_record(
    record: &[u8],
) -> fgumi_sam::alignment_tags::ConversionScoring {
    use fgumi_sam::alignment_tags::ConversionScoring;
    let aux = bam_fields::aux_data_slice(record);
    if bam_fields::find_tag_type(aux, SamTag::CU).is_some() && !is_duplex_consensus(aux) {
        ConversionScoring::Hidden
    } else {
        ConversionScoring::Literal
    }
}

/// Detects if a raw BAM record is a duplex consensus.
///
/// Matches fgbio `Umis.isFgbioDuplexConsensus` (`Umis.scala:144`), which requires **both**
/// the AB and BA raw-read-count tags (`aD` **and** `bD`). A read carrying only one of the two
/// is treated as a simplex consensus (fgbio routes it to `filterVanillaConsensusRead`); routing
/// it through the duplex filter would spuriously reject it on the BA tier (worst-strand depth 0).
#[must_use]
pub fn is_duplex_consensus(aux_data: &[u8]) -> bool {
    bam_fields::find_tag_type(aux_data, SamTag::AD).is_some()
        && bam_fields::find_tag_type(aux_data, SamTag::BD).is_some()
}

/// Detects if a raw BAM record is a simplex consensus.
///
/// The raw-aux counterpart of [`crate::tags::is_simplex_consensus`]: the `cD`
/// (`RawReadCount`) tag without the duplex pair. A read carrying only one of
/// `aD`/`bD` counts as simplex, matching fgbio — see [`is_duplex_consensus`].
#[must_use]
pub fn is_simplex_consensus(aux_data: &[u8]) -> bool {
    bam_fields::find_tag_type(aux_data, SamTag::CD).is_some() && !is_duplex_consensus(aux_data)
}

/// Detects if a raw BAM record is any kind of consensus read.
///
/// The raw-aux counterpart of [`crate::tags::is_consensus`], which matches
/// fgbio's `Umis.isFgbioStyleConsensus()`. Use this to reject consensus input to
/// tools that require pre-consensus, UMI-grouped reads.
#[must_use]
pub fn is_consensus(aux_data: &[u8]) -> bool {
    is_simplex_consensus(aux_data) || is_duplex_consensus(aux_data)
}

/// The scalar consensus tags the read-level filter needs, collected in a
/// **single** pass over a record's aux block.
///
/// The classification + filter path previously walked the aux block ~9 times
/// per record (`is_duplex_consensus` reads aD+bD; `filter_read` reads cD+cE;
/// `filter_duplex_read` re-reads cD+cE then aD/aM/bD/bM/aE/bE — with aD/bD
/// scanned up to three times). Each `find_*_tag` restarts at offset 0, so the
/// cost is `O(tags · aux_len)` per record. This struct collapses those into one
/// `O(aux_len)` walk that decodes each target tag's value at its first occurrence
/// — mirroring the existing single-pass
/// [`fgumi_raw_bam::extract_aux_string_tags`] pattern.
///
/// Every value is decoded with the **same** primitives the per-tag `find_*`
/// functions use (`extract_int_value` and the `f`-type layout), and the walk
/// mirrors `find_tag_position`'s "resolve on the first occurrence of a tag id,
/// regardless of its value type" semantics via the `seen` bitset:
/// a tag is locked to its first entry even when that entry fails to decode, so a
/// later duplicate cannot override it and a present-but-mistyped `aD`/`bD` still
/// counts as *present* for [`is_duplex`](Self::is_duplex) exactly as
/// `find_tag_type` / [`is_duplex_consensus`] treat it. The classification and
/// threshold decisions are therefore identical to the per-tag path even on
/// malformed aux; `consensus_scalar_tags_match_find` and
/// `consensus_scalar_tags_match_find_malformed` pin that equivalence.
#[derive(Debug, Clone, Copy, Default)]
pub struct ConsensusScalarTags {
    /// `cD` raw-read depth (final consensus), decoded at its first occurrence.
    cd: Option<i64>,
    /// `cE` error rate (final consensus).
    ce: Option<f32>,
    /// `aD` AB raw-read depth (also drives duplex classification via `seen`).
    ad: Option<i64>,
    /// `aM` AB min-depth fallback when `aD` is absent.
    am: Option<i64>,
    /// `bD` BA raw-read depth (also drives duplex classification via `seen`).
    bd: Option<i64>,
    /// `bM` BA min-depth fallback when `bD` is absent.
    bm: Option<i64>,
    /// `aE` AB error rate.
    ae: Option<f32>,
    /// `bE` BA error rate.
    be: Option<f32>,
    /// Which target tags were *seen* (by tag id, at their first occurrence)
    /// during the walk, independent of whether the value decoded. Mirrors
    /// `find_tag_position`'s lock-on-first-occurrence: a set bit both stops a
    /// later duplicate from overriding the slot and records tag *presence*
    /// separately from a successful decode (so `is_duplex` matches
    /// `find_tag_type`'s type-agnostic presence check).
    seen: u8,
}

impl ConsensusScalarTags {
    // `seen` bit per target tag, in the same order as the fields above.
    const SEEN_CD: u8 = 1 << 0;
    const SEEN_CE: u8 = 1 << 1;
    const SEEN_AD: u8 = 1 << 2;
    const SEEN_AM: u8 = 1 << 3;
    const SEEN_BD: u8 = 1 << 4;
    const SEEN_BM: u8 = 1 << 5;
    const SEEN_AE: u8 = 1 << 6;
    const SEEN_BE: u8 = 1 << 7;

    /// Collect the scalar consensus tags from `aux_data` in one pass.
    #[must_use]
    pub fn from_aux(aux_data: &[u8]) -> Self {
        let mut out = Self::default();
        // Decode an `f`-type value at entry offset `p` (value bytes at `p + 3`),
        // matching `find_float_tag`'s layout and its silent-skip on a short slice.
        let float_at = |p: usize| -> Option<f32> {
            aux_data.get(p + 3..p + 7).map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
        };

        let mut p = 0;
        while p + 3 <= aux_data.len() {
            let tag = [aux_data[p], aux_data[p + 1]];
            let val_type = aux_data[p + 2];
            // Resolve each tag on its FIRST occurrence by tag id, exactly like
            // `find_tag_position`. The `seen` bit locks the slot even when the
            // decode below yields `None` (unexpected value type), so a later
            // duplicate cannot override it and presence is recorded independently
            // of a successful decode — keeping the single-pass result identical to
            // the per-tag `find_int_tag`/`find_float_tag`/`find_tag_type` path.
            if tag == *SamTag::CD.as_tag_bytes() && out.seen & Self::SEEN_CD == 0 {
                out.seen |= Self::SEEN_CD;
                out.cd = bam_fields::extract_int_value(aux_data, p, val_type);
            } else if tag == *SamTag::CE.as_tag_bytes() && out.seen & Self::SEEN_CE == 0 {
                out.seen |= Self::SEEN_CE;
                out.ce = if val_type == b'f' { float_at(p) } else { None };
            } else if tag == *SamTag::AD.as_tag_bytes() && out.seen & Self::SEEN_AD == 0 {
                out.seen |= Self::SEEN_AD;
                out.ad = bam_fields::extract_int_value(aux_data, p, val_type);
            } else if tag == *SamTag::AM.as_tag_bytes() && out.seen & Self::SEEN_AM == 0 {
                out.seen |= Self::SEEN_AM;
                out.am = bam_fields::extract_int_value(aux_data, p, val_type);
            } else if tag == *SamTag::BD.as_tag_bytes() && out.seen & Self::SEEN_BD == 0 {
                out.seen |= Self::SEEN_BD;
                out.bd = bam_fields::extract_int_value(aux_data, p, val_type);
            } else if tag == *SamTag::BM.as_tag_bytes() && out.seen & Self::SEEN_BM == 0 {
                out.seen |= Self::SEEN_BM;
                out.bm = bam_fields::extract_int_value(aux_data, p, val_type);
            } else if tag == *SamTag::AE.as_tag_bytes() && out.seen & Self::SEEN_AE == 0 {
                out.seen |= Self::SEEN_AE;
                out.ae = if val_type == b'f' { float_at(p) } else { None };
            } else if tag == *SamTag::BE.as_tag_bytes() && out.seen & Self::SEEN_BE == 0 {
                out.seen |= Self::SEEN_BE;
                out.be = if val_type == b'f' { float_at(p) } else { None };
            }

            match bam_fields::tag_value_size(val_type, &aux_data[p + 3..]) {
                Some(size) => p += 3 + size,
                None => break,
            }
        }
        out
    }

    /// Duplex classification: both `aD` and `bD` **present** (any value type),
    /// matching [`is_duplex_consensus`], which requires both tiers' raw-read-count
    /// tags via a `find_tag_type` presence check. Keyed off the `seen` bitset
    /// rather than a successful integer decode, so a present-but-mistyped `aD`/`bD`
    /// classifies as duplex exactly as the per-tag path does.
    #[must_use]
    pub fn is_duplex(&self) -> bool {
        self.seen & Self::SEEN_AD != 0 && self.seen & Self::SEEN_BD != 0
    }
}

/// Read-level simplex filter over pre-extracted scalar tags — the allocation-
/// free equivalent of [`filter_read`] that reuses a single aux walk.
///
/// # Errors
///
/// Returns an error when the `cD`/`cE` consensus tags are absent (fgbio fails
/// hard rather than silently keeping a non-consensus read).
pub fn filter_read_tags(
    tags: &ConsensusScalarTags,
    thresholds: &FilterThresholds,
) -> Result<FilterResult> {
    anyhow::ensure!(
        tags.cd.is_some() && tags.ce.is_some(),
        "read does not appear to have consensus calling tags (cD/cE) present; \
         FilterConsensusReads requires reads produced by consensus calling"
    );
    if let Some(depth) = tags.cd {
        let min_reads = i64::try_from(thresholds.min_reads).unwrap_or(i64::MAX);
        if depth < min_reads {
            return Ok(FilterResult::InsufficientReads);
        }
    }
    if let Some(error_rate) = tags.ce
        && f64::from(error_rate) > thresholds.max_read_error_rate
    {
        return Ok(FilterResult::ExcessiveErrorRate);
    }
    Ok(FilterResult::Pass)
}

/// Read-level duplex filter over pre-extracted scalar tags — the allocation-
/// free equivalent of [`filter_duplex_read`] that reuses a single aux walk.
///
/// # Errors
///
/// Returns an error when the `cD`/`cE` consensus tags are absent.
pub fn filter_duplex_read_tags(
    tags: &ConsensusScalarTags,
    cc_thresholds: &FilterThresholds,
    ab_thresholds: &FilterThresholds,
    ba_thresholds: &FilterThresholds,
) -> Result<FilterResult> {
    let result = filter_read_tags(tags, cc_thresholds)?;
    if result != FilterResult::Pass {
        return Ok(result);
    }

    let ab_depth = tags.ad.or(tags.am);
    let ba_depth = tags.bd.or(tags.bm);
    let ab_error = tags.ae;
    let ba_error = tags.be;

    Ok(duplex_tier_result(ab_depth, ba_depth, ab_error, ba_error, ab_thresholds, ba_thresholds))
}

/// Filters a raw consensus read based on per-read tags (cD depth, cE error rate).
///
/// # Errors
///
/// Returns an error if the aux data cannot be parsed.
pub fn filter_read(aux_data: &[u8], thresholds: &FilterThresholds) -> Result<FilterResult> {
    // The per-read consensus depth (cD) and error-rate (cE) tags must both be present:
    // FilterConsensusReads only operates on consensus reads. fgbio fails hard when either
    // is absent (FilterConsensusReads.scala:242) rather than silently keeping the read.
    let depth = bam_fields::find_int_tag(aux_data, SamTag::CD);
    let error_rate = bam_fields::find_float_tag(aux_data, SamTag::CE);
    anyhow::ensure!(
        depth.is_some() && error_rate.is_some(),
        "read does not appear to have consensus calling tags (cD/cE) present; \
         FilterConsensusReads requires reads produced by consensus calling"
    );

    // Check minimum reads (cD tag — any integer type)
    if let Some(depth) = depth {
        let min_reads = i64::try_from(thresholds.min_reads).unwrap_or(i64::MAX);
        if depth < min_reads {
            return Ok(FilterResult::InsufficientReads);
        }
    }

    // Check maximum error rate (cE tag — Float)
    if let Some(error_rate) = error_rate
        && f64::from(error_rate) > thresholds.max_read_error_rate
    {
        return Ok(FilterResult::ExcessiveErrorRate);
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

    Ok(duplex_tier_result(ab_depth, ba_depth, ab_error, ba_error, ab_thresholds, ba_thresholds))
}

/// Shared AB/BA tier comparison for the duplex read-level filter.
///
/// Single source of truth for both [`filter_duplex_read`] (per-tag walks) and
/// [`filter_duplex_read_tags`] (single-pass), so the two cannot drift. `best`/
/// `worst` are per-metric extremes across the two strands — NOT biological
/// AB/BA values. `ab_thresholds` is the stricter tier (checked against the
/// best), `ba_thresholds` the lenient tier (checked against the worst). Matches
/// fgbio's `abMaxDepth`/`abError` semantics.
fn duplex_tier_result(
    ab_depth: Option<i64>,
    ba_depth: Option<i64>,
    ab_error: Option<f32>,
    ba_error: Option<f32>,
    ab_thresholds: &FilterThresholds,
    ba_thresholds: &FilterThresholds,
) -> FilterResult {
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
        (None, None) => return FilterResult::Pass,
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
        return FilterResult::InsufficientReads;
    }
    if f64::from(best_error) > ab_thresholds.max_read_error_rate {
        return FilterResult::ExcessiveErrorRate;
    }

    // Lenient BA tier: worst-per-metric value must still clear the threshold.
    #[expect(
        clippy::cast_sign_loss,
        clippy::cast_possible_truncation,
        reason = "depth values are non-negative and fit in usize on all supported platforms"
    )]
    if (worst_depth as usize) < ba_thresholds.min_reads {
        return FilterResult::InsufficientReads;
    }
    if f64::from(worst_error) > ba_thresholds.max_read_error_rate {
        return FilterResult::ExcessiveErrorRate;
    }

    FilterResult::Pass
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

/// Mean base quality over the **full** read length, i.e. the sum of all base qualities
/// divided by the read length, **including** no-call (N) bases.
///
/// This mirrors fgbio's read-level mean-quality filter, which is computed *prior to any
/// masking* over the whole read (`FilterConsensusReads.scala:247`:
/// `rec.quals.foldLeft(0)(_ + _) / rec.length`). It differs from [`compute_read_stats`],
/// whose mean is taken over non-N bases only — the correct denominator for this filter is
/// the full read length, so a read whose low-quality bases are masked (and thus excluded
/// from a non-N mean) is still judged on its original, unmasked quality.
///
/// Returns `0.0` for an empty read.
///
/// # Panics
/// Panics if the record is shorter than `MIN_BAM_RECORD_LEN` (36 bytes).
#[must_use]
pub fn mean_base_quality_full_length(bam: &[u8]) -> f64 {
    assert!(bam.len() >= bam_fields::MIN_BAM_RECORD_LEN, "BAM record too short");
    let qual_off = bam_fields::qual_offset(bam);
    let len = RawRecordView::new(bam).l_seq() as usize;
    if len == 0 {
        return 0.0;
    }
    let mut qual_sum = 0u64;
    for i in 0..len {
        qual_sum += u64::from(bam_fields::get_qual(bam, qual_off, i));
    }
    #[expect(
        clippy::cast_precision_loss,
        reason = "precision loss is acceptable for quality averaging"
    )]
    let mean = qual_sum as f64 / len as f64;
    mean
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
    anyhow::ensure!(record.len() >= bam_fields::MIN_BAM_RECORD_LEN, "BAM record too short");
    let seq_off = bam_fields::seq_offset(record);
    let qual_off = bam_fields::qual_offset(record);
    let len = RawRecordView::new(record).l_seq() as usize;
    let aux_off = bam_fields::aux_data_offset_from_record(record).unwrap_or(record.len());

    // Pre-read per-base arrays into owned Vecs to release the immutable borrow on record
    let cd_vals = bam_fields::find_array_tag(&record[aux_off..], SamTag::CD_BASES)
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));
    let ce_vals = bam_fields::find_array_tag(&record[aux_off..], SamTag::CE_BASES)
        .map(|r| bam_fields::array_tag_to_vec_u16(&r));

    // Per-base depth/error masking only applies when BOTH per-base tags are present, matching
    // fgbio's `pb = depths != null && errors != null` guard (FilterConsensusReads.scala:281,
    // 287-289). When they are absent we must NOT treat depth as 0 and mask everything —
    // only the base-quality mask applies.
    let has_per_base = cd_vals.is_some() && ce_vals.is_some();

    let mut masked_count = 0;
    for i in 0..len {
        let depth = cd_vals.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let errors = ce_vals.as_ref().map_or(0u16, |v| v.get(i).copied().unwrap_or(0));
        let qual = bam_fields::get_qual(record, qual_off, i);

        let should_mask = min_base_quality.is_some_and(|min_qual| qual < min_qual)
            || (has_per_base && (depth as usize) < thresholds.min_reads)
            || (has_per_base
                && depth > 0
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
    anyhow::ensure!(record.len() >= bam_fields::MIN_BAM_RECORD_LEN, "BAM record too short");
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

/// Thresholds for methylation depth filtering (`--min-methylation-depth x[,y[,z]]`).
///
/// Mirrors the 1-3 value pattern of `--min-reads`: missing values are filled from the last
/// provided. As there, the tiers are best/worst, not AB/BA: on a duplex `CpG` the two halves
/// (the reference C and G) are read by different strands, and `best`/`worst` bind to the
/// better- and worse-supported half.
#[derive(Debug, Clone)]
pub struct MethylationDepthThresholds {
    /// Minimum depth at a simplex informative position; minimum summed depth of a duplex `CpG`.
    pub total: usize,
    /// Minimum depth of the better-supported half of a duplex `CpG`, and of a duplex
    /// informative position outside a `CpG` pair.
    pub best: usize,
    /// Minimum depth of the worse-supported half of a duplex `CpG`.
    pub worst: usize,
}

impl MethylationDepthThresholds {
    /// Creates thresholds from a 1-3 element slice, filling missing values from the last.
    ///
    /// # Panics
    /// Panics if `values` is empty.
    #[must_use]
    pub fn from_values(values: &[usize]) -> Self {
        assert!(!values.is_empty(), "min-methylation-depth must have at least 1 value");
        let [total, best, worst] = expand_three_from_last(values);
        Self { total, best, worst }
    }
}

/// The reference base aligned to one query position, with its reference neighbours.
///
/// Built by [`resolve_ref_bases_for_record`]. `prev`/`next` are the reference bases before and
/// after this one whether or not the read covers them, so `CpG` context comes from the reference
/// rather than from neighbouring query positions (which an indel or a read edge separates from
/// their reference neighbours).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RefBase {
    /// Uppercase reference base.
    pub base: u8,
    /// Uppercase reference base before this one, if any.
    pub prev: Option<u8>,
    /// Uppercase reference base after this one, if any.
    pub next: Option<u8>,
}

impl RefBase {
    /// Whether this is the C of a reference `CpG`.
    #[must_use]
    pub fn is_cpg_c(&self) -> bool {
        self.base == b'C' && self.next == Some(b'G')
    }

    /// Whether this is the G of a reference `CpG`.
    #[must_use]
    pub fn is_cpg_g(&self) -> bool {
        self.base == b'G' && self.prev == Some(b'C')
    }
}

/// Where a consensus record's methylation sites come from.
#[derive(Debug, Clone, Copy)]
pub enum MethylationSites<'a> {
    /// A single-strand consensus is reference-anchored (a converted `T` could be a C>T
    /// variant): its sites are the aligned reference bases, from [`resolve_ref_bases_for_record`].
    Reference(&'a [Option<RefBase>]),
    /// A duplex consensus is molecule-based: its sites come from its SEQ, the molecule's
    /// sequence, before any masking (`None`: the record's current SEQ). No reference is needed.
    Sequence(Option<&'a [u8]>),
}

/// What [`mask_methylation_sites_raw`] did to a record.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MethylationSiteMaskOutcome {
    /// Bases newly masked to `N` (single-strand consensus records).
    pub masked: usize,
    /// Positions (stored orientation) of the sites that failed and were masked to `N`
    /// (single-strand consensus records), including any already `N`.
    pub masked_sites: Vec<usize>,
    /// Positions (stored orientation) whose methylation calls must be dropped from
    /// `MM`/`ML`/`am`/`bm` (duplex consensus records); see [`drop_modifications_at_raw`].
    pub dropped_calls: Vec<usize>,
}

/// Per-site methylation masks applied by `filter`.
#[derive(Debug, Clone, Default)]
pub struct MethylationSiteMasks {
    /// `--min-methylation-depth`, if set.
    pub depth: Option<MethylationDepthThresholds>,
    /// `--require-strand-methylation-agreement` (duplex only: an error with
    /// [`MethylationSites::Reference`]).
    pub require_strand_agreement: bool,
}

/// Masks methylation sites in a raw consensus record per `--min-methylation-depth` and
/// `--require-strand-methylation-agreement`.
///
/// Works in the record's stored orientation: the reference bases (from
/// [`resolve_ref_bases_for_record`]) and the `cu`/`ct` arrays must index the same positions as SEQ, i.e. per-base tags must
/// already be reversed for reverse-mapped records. `cu + ct` at a position is the depth of the
/// one strand informative there, so no per-strand array is read for depth or agreement.
///
/// - Single-strand consensus (reference-anchored: a converted `T` could be a C>T variant, so the
///   reference tells where a cytosine was): informative positions are aligned reference bases
///   of the record's class (`C` for top-strand-type records, R1 or fragment forward and R2
///   reverse; `G` for the others). One whose SEQ base is the unconverted or converted base (or
///   `N`) and whose depth is below `total` is masked to `N` (methylation is read from SEQ); a
///   third base (a variant) carries no call and is left alone ([`MethylationSites::Reference`]).
/// - Duplex consensus (molecule-based: SEQ is the molecule's sequence;
///   [`MethylationSites::Sequence`]): a site is a `C` or `G` in SEQ before masking, and a site
///   with counts carries a call; on a record with both strands a
///   site without counts carries none. Methylation is carried by `MM`/`ML`, so a failing call is dropped (returned in
///   [`MethylationSiteMaskOutcome::dropped_calls`]) and SEQ is left alone. A `CpG` is a `C`
///   followed in the read by a `G`; its halves are called by different strands, and both calls
///   are dropped unless the summed depth is at least `total`, the better half at least `best`
///   and the worse half at least `worst`; then, if required, both are dropped when the halves
///   disagree on methylation (majority of `cu` vs `ct`; a tied half calls neither and does not
///   disagree). A half without a call is not tested
///   and its partner is tested alone. Every other call must reach `best`. On a single-strand
///   duplex record a site without counts is the absent strand's half: it has depth 0 in its
///   pair (so the pair fails when `worst` > 0) and is not tested alone.
///   `N` is never a site, so a `C` or `G` next to it is tested alone. Needs no reference.
///
/// Masked bases become `N` with quality [`MIN_PHRED`]; bases already `N` are not counted.
/// Records without `cu`/`ct` are returned unchanged. The counts are not changed here; see
/// [`zero_methylation_counts_at_raw`].
///
/// # Errors
/// Returns an error if the record is too short, or if strand agreement is requested with
/// [`MethylationSites::Reference`] (it compares the two strands of a duplex record).
pub fn mask_methylation_sites_raw(
    record: &mut [u8],
    sites: MethylationSites<'_>,
    masks: &MethylationSiteMasks,
    tags: &MethylationTags,
) -> Result<MethylationSiteMaskOutcome> {
    anyhow::ensure!(record.len() >= bam_fields::MIN_BAM_RECORD_LEN, "BAM record too short");
    anyhow::ensure!(
        !(masks.require_strand_agreement && matches!(sites, MethylationSites::Reference(_))),
        "strand methylation agreement compares the two strands of a duplex record and cannot \
         apply to a single-strand consensus"
    );
    if tags.cu.is_none() && tags.ct.is_none() {
        return Ok(MethylationSiteMaskOutcome::default());
    }
    let pre_mask_seq = match sites {
        MethylationSites::Sequence(seq) => seq,
        MethylationSites::Reference(_) => None,
    };

    let view = RawRecordView::new(record);
    let len = view.l_seq() as usize;
    let current_seq;
    let seq = if let Some(seq) = pre_mask_seq {
        seq
    } else {
        current_seq = view.sequence_vec();
        &current_seq
    };
    let count_at = |counts: &Option<Vec<u16>>, i: usize| {
        counts.as_ref().and_then(|v| v.get(i)).map_or(0usize, |&c| usize::from(c))
    };
    let depth_at = |i: usize| count_at(&tags.cu, i) + count_at(&tags.ct, i);
    let base_at = |i: usize| seq.get(i).map(u8::to_ascii_uppercase);

    if matches!(sites, MethylationSites::Sequence(_)) {
        // A half's call is the majority of its counts; a tie calls neither.
        let call_at = |i: usize| match count_at(&tags.cu, i).cmp(&count_at(&tags.ct, i)) {
            std::cmp::Ordering::Greater => Some(true),
            std::cmp::Ordering::Less => Some(false),
            std::cmp::Ordering::Equal => None,
        };
        // A site with counts carries a call. A site without them carries none on a record with
        // both strands; on a single-strand record it is the absent strand's half, which counts
        // as depth 0 in its `CpG` pair.
        let both_strands = tags.has_ab && tags.has_ba;
        let is_site = |i: usize| matches!(base_at(i), Some(b'C' | b'G'));
        let has_call = |i: usize| is_site(i) && depth_at(i) > 0;
        let in_pair = |i: usize| has_call(i) || (!both_strands && is_site(i));
        let fails_single =
            |i: usize| has_call(i) && masks.depth.as_ref().is_some_and(|t| depth_at(i) < t.best);
        let mut outcome = MethylationSiteMaskOutcome::default();
        let mut i = 0;
        while i < len {
            if base_at(i) == Some(b'C') && base_at(i + 1) == Some(b'G') {
                let j = i + 1;
                if (has_call(i) || has_call(j)) && in_pair(i) && in_pair(j) {
                    let (c_depth, g_depth) = (depth_at(i), depth_at(j));
                    let fails_depth = masks.depth.as_ref().is_some_and(|t| {
                        c_depth + g_depth < t.total
                            || c_depth.max(g_depth) < t.best
                            || c_depth.min(g_depth) < t.worst
                    });
                    let disagrees = masks.require_strand_agreement
                        && matches!((call_at(i), call_at(j)), (Some(c), Some(g)) if c != g);
                    if fails_depth || disagrees {
                        outcome.dropped_calls.extend([i, j]);
                    }
                } else {
                    // At most one half carries a call: it is tested alone.
                    outcome.dropped_calls.extend([i, j].into_iter().filter(|&k| fails_single(k)));
                }
                i = j + 1;
                continue;
            }
            if fails_single(i) {
                outcome.dropped_calls.push(i);
            }
            i += 1;
        }
        return Ok(outcome);
    }

    let MethylationSites::Reference(ref_bases) = sites else { unreachable!("handled above") };
    let (own_class, converted) = if fgumi_sam::alignment_tags::is_top_strand(view.flags()) {
        (b'C', b'T')
    } else {
        (b'G', b'A')
    };
    let should_mask = |i: usize| {
        masks.depth.as_ref().is_some_and(|t| {
            ref_bases.get(i).copied().flatten().is_some_and(|r| r.base == own_class)
                && matches!(base_at(i), Some(b) if b == own_class || b == converted || b == b'N')
                && depth_at(i) < t.total
        })
    };
    let seq_off = bam_fields::seq_offset(record);
    let qual_off = bam_fields::qual_offset(record);
    let mut outcome = MethylationSiteMaskOutcome {
        masked_sites: (0..len).filter(|&i| should_mask(i)).collect(),
        ..MethylationSiteMaskOutcome::default()
    };
    for &i in &outcome.masked_sites {
        if !bam_fields::is_base_n(record, seq_off, i) {
            outcome.masked += 1;
        }
        bam_fields::mask_base(record, seq_off, i);
        bam_fields::set_qual(record, qual_off, i, MIN_PHRED);
    }
    Ok(outcome)
}

/// Zeroes the methylation counts (`cu`/`ct` and the per-strand `au`/`at`/`bu`/`bt`) at
/// `positions` (stored orientation), so a reader of the counts sees no call where `filter`
/// rejected one: the sites [`mask_methylation_sites_raw`] masked or whose calls it dropped.
/// Counts of a type other than `B:s`/`B:S` are left alone.
pub fn zero_methylation_counts_at_raw(record: &mut Vec<u8>, positions: &[usize]) {
    use crate::tags::per_base;
    if positions.is_empty() {
        return;
    }
    let tags = [
        per_base::UNCONVERTED_COUNT,
        per_base::CONVERTED_COUNT,
        per_base::AB_UNCONVERTED_COUNT,
        per_base::AB_CONVERTED_COUNT,
        per_base::BA_UNCONVERTED_COUNT,
        per_base::BA_CONVERTED_COUNT,
    ];
    let arrays: Vec<_> = {
        let aux = bam_fields::aux_data_slice(record);
        tags.into_iter()
            .filter_map(|tag| {
                let array = bam_fields::find_array_tag(aux, tag)?;
                matches!(array.elem_type, b's' | b'S').then(|| {
                    let values: Vec<i16> = array
                        .data
                        .chunks_exact(2)
                        .map(|b| i16::from_le_bytes([b[0], b[1]]))
                        .collect();
                    (tag, array.elem_type, values)
                })
            })
            .collect()
    };
    let mut editor = bam_fields::RawTagsEditor::from_vec(record);
    for (tag, elem_type, mut values) in arrays {
        for &i in positions {
            if let Some(value) = values.get_mut(i) {
                *value = 0;
            }
        }
        if elem_type == b's' {
            editor.update_array_i16(tag, &values);
        } else {
            #[expect(clippy::cast_sign_loss, reason = "reinterprets the stored B:S bytes")]
            let unsigned: Vec<u16> = values.iter().map(|&v| v as u16).collect();
            editor.update_array_u16(tag, &unsigned);
        }
    }
}

/// Resolves the reference bases for a raw BAM record's alignment region.
///
/// Returns one entry per query position: the aligned reference base with its reference
/// neighbours ([`RefBase`]), or `None` for insertions/soft-clips. Returns `None` if
/// the record is unmapped or the reference cannot be resolved.
#[cfg(feature = "simplex")]
#[expect(
    clippy::cast_sign_loss,
    reason = "ref_id and pos are non-negative for mapped records (checked above)"
)]
pub fn resolve_ref_bases_for_record(
    record: &[u8],
    reference: &dyn crate::methylation::RefBaseProvider,
    ref_names: &[String],
) -> Option<Vec<Option<RefBase>>> {
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

    let fetch = |pos: u64| {
        let base = if let Some(seq) = ref_seq {
            usize::try_from(pos).ok().and_then(|p| seq.get(p).copied())
        } else {
            reference.base_at_0based(ref_name, pos)
        };
        base.map(|b| b.to_ascii_uppercase())
    };

    let len = RawRecordView::new(record).l_seq() as usize;
    let mut result = Vec::with_capacity(len);
    let mut ref_pos = alignment_start;

    for op in RawRecordView::new(record).cigar_ops_iter() {
        let op_len = (op >> 4) as usize;
        let kind = bam_fields::cigar_op_kind(op);
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // Slide a prev/current/next window along the run: one fetch per base.
                let mut prev = ref_pos.checked_sub(1).and_then(fetch);
                let mut current = fetch(ref_pos);
                for _ in 0..op_len {
                    let next = fetch(ref_pos + 1);
                    result.push(current.map(|base| RefBase { base, prev, next }));
                    prev = current;
                    current = next;
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

/// Checks the bisulfite/enzymatic conversion fraction at non-CpG cytosines.
///
/// Returns `true` if the read passes (conversion fraction >= threshold), or if
/// there are no non-CpG cytosine positions to evaluate.
///
/// At non-CpG cytosines the expected behavior is conversion. They appear at reference C
/// (next base not G) for OT-derived reads and at reference G (previous base not C) for
/// OB-derived reads; both are evaluated where the read has evidence (`cu + ct > 0`).
/// A low conversion rate suggests incomplete enzymatic conversion.
#[cfg(feature = "simplex")]
pub fn check_conversion_fraction_raw(
    record: &[u8],
    min_fraction: f64,
    reference: &dyn crate::methylation::RefBaseProvider,
    ref_names: &[String],
    methylation_mode: crate::MethylationMode,
) -> bool {
    let meth_tags = MethylationTags::from_record(record);
    let ref_bases;
    let sites = if is_duplex_consensus(bam_fields::aux_data_slice(record)) {
        MethylationSites::Sequence(None)
    } else {
        ref_bases = resolve_ref_bases_for_record(record, reference, ref_names);
        let Some(ref_bases) = ref_bases.as_deref() else {
            return true; // unaligned single-strand reads pass
        };
        MethylationSites::Reference(ref_bases)
    };
    check_conversion_fraction_raw_with_ref_bases_and_tags(
        record,
        min_fraction,
        sites,
        &meth_tags,
        methylation_mode,
    )
}

/// Like `check_conversion_fraction_raw` but accepts the record's methylation sites and
/// pre-parsed methylation tags, avoiding all redundant work.
///
/// For EM-Seq, checks `ct / (cu + ct) >= threshold` at non-CpG cytosines; for TAPs,
/// `cu / (cu + ct)` instead, since non-CpG Cs should NOT be converted in TAPs. On a single-strand
/// consensus ([`MethylationSites::Reference`]) the cytosines are non-CpG reference C and G. On a
/// duplex consensus ([`MethylationSites::Sequence`], SEQ before masking: the molecule's
/// sequence) they are a `C` not followed by `G` and a `G` not preceded by `C` in SEQ; one whose
/// neighbour is `N` or outside the read has unknown context and is not counted.
#[must_use]
pub fn check_conversion_fraction_raw_with_ref_bases_and_tags(
    record: &[u8],
    min_fraction: f64,
    sites: MethylationSites<'_>,
    methylation_tags: &MethylationTags,
    methylation_mode: crate::MethylationMode,
) -> bool {
    // Conversion fraction is meaningless without a methylation mode — pass through
    if methylation_mode == crate::MethylationMode::Disabled {
        return true;
    }
    // If no methylation tags, pass
    if methylation_tags.cu.is_none() && methylation_tags.ct.is_none() {
        return true;
    }

    let len = RawRecordView::new(record).l_seq() as usize;
    let current_seq;
    let is_non_cpg: Box<dyn Fn(usize) -> bool + '_> = match sites {
        MethylationSites::Sequence(seq) => {
            let seq = if let Some(seq) = seq {
                seq
            } else {
                current_seq = RawRecordView::new(record).sequence_vec();
                &current_seq
            };
            let base_at = move |i: usize| seq.get(i).map(u8::to_ascii_uppercase);
            Box::new(move |i: usize| match base_at(i) {
                Some(b'C') => matches!(base_at(i + 1), Some(b'A' | b'C' | b'T')),
                Some(b'G') => i > 0 && matches!(base_at(i - 1), Some(b'A' | b'G' | b'T')),
                _ => false,
            })
        }
        // Non-CpG reference C (next reference base not G) and G (previous reference base not
        // C): OT-derived reads carry evidence at reference C, OB-derived reads at reference G.
        MethylationSites::Reference(ref_bases) => {
            Box::new(move |i: usize| match ref_bases.get(i).copied().flatten() {
                Some(r) if r.base == b'C' => !r.is_cpg_c(),
                Some(r) if r.base == b'G' => !r.is_cpg_g(),
                _ => false,
            })
        }
    };

    let mut total_numerator: u64 = 0;
    let mut total_evidence: u64 = 0;

    for i in (0..len).filter(|&i| is_non_cpg(i)) {
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
    use fgumi_raw_bam::flags;
    use rstest::rstest;

    /// Conversions are hidden from NM/UQ only on a simplex methylation consensus (`cu`, no
    /// duplex `aD`/`bD` pair); a duplex consensus or a non-methylation record is literal.
    #[rstest]
    #[case::simplex_methylation(true, false, false, true)]
    #[case::duplex_methylation(true, true, true, false)]
    #[case::non_methylation(false, false, false, false)]
    #[case::ad_only_is_not_duplex(true, true, false, true)]
    fn test_conversion_scoring_for_record(
        #[case] with_cu: bool,
        #[case] with_ad: bool,
        #[case] with_bd: bool,
        #[case] hidden: bool,
    ) {
        use fgumi_sam::alignment_tags::ConversionScoring;
        let mut b = RawSamBuilder::new();
        b.flags(0)
            .ref_id(0)
            .pos(0)
            .mapq(60)
            .cigar_ops(&[4 << 4])
            .sequence(b"ACGT")
            .qualities(&[30; 4]);
        if with_cu {
            b.add_array_i16(SamTag::CU, &[0; 4]);
        }
        if with_ad {
            b.add_int_tag(SamTag::AD, 3);
        }
        if with_bd {
            b.add_int_tag(SamTag::BD, 2);
        }
        let record = b.build();
        let expected = if hidden { ConversionScoring::Hidden } else { ConversionScoring::Literal };
        assert_eq!(conversion_scoring_for_record(record.as_ref()), expected);
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
        assert_eq!((t.total, t.best, t.worst), (5, 5, 5));
    }

    #[test]
    fn test_methylation_depth_thresholds_two_values() {
        let t = MethylationDepthThresholds::from_values(&[10, 3]);
        assert_eq!((t.total, t.best, t.worst), (10, 3, 3));
    }

    #[test]
    fn test_methylation_depth_thresholds_three_values() {
        let t = MethylationDepthThresholds::from_values(&[10, 5, 2]);
        assert_eq!((t.total, t.best, t.worst), (10, 5, 2));
    }

    // -- mask_methylation_sites_raw tests --
    //
    // Fixtures work in genomic orientation (as after `--reverse-per-base-tags`) and hand the
    // function a reference-base map directly (`-` = unaligned query base). Reference
    // `AACGTCAGCGT`: CpGs at 2-3 and 8-9, a non-CpG C at 5, a non-CpG G at 7.

    const R1_FWD: u16 = flags::PAIRED | flags::FIRST_SEGMENT;
    const R1_REV: u16 = R1_FWD | flags::REVERSE;
    const R2_FWD: u16 = flags::PAIRED | flags::LAST_SEGMENT;
    const R2_REV: u16 = R2_FWD | flags::REVERSE;
    const FRAG_FWD: u16 = 0;
    const FRAG_REV: u16 = flags::REVERSE;

    /// Which per-strand count arrays a duplex fixture carries.
    #[derive(Clone, Copy, Debug)]
    enum Strands {
        Both,
        AbOnly,
        BaOnly,
    }

    /// `cu`/`ct` for a fixture; `ct` defaults to all zeros.
    #[derive(Clone, Copy, Debug)]
    struct Counts {
        cu: &'static [i16],
        ct: Option<&'static [i16]>,
    }

    impl Counts {
        fn depth(cu: &'static [i16]) -> Self {
            Self { cu, ct: None }
        }

        fn cu_ct(cu: &'static [i16], ct: &'static [i16]) -> Self {
            Self { cu, ct: Some(ct) }
        }
    }

    fn depth(values: &[usize]) -> MethylationSiteMasks {
        MethylationSiteMasks {
            depth: Some(MethylationDepthThresholds::from_values(values)),
            require_strand_agreement: false,
        }
    }

    fn agreement(values: Option<&[usize]>) -> MethylationSiteMasks {
        MethylationSiteMasks {
            depth: values.map(MethylationDepthThresholds::from_values),
            require_strand_agreement: true,
        }
    }

    /// Reference context for a read described by `reference`, one character per query
    /// position except deletions: an uppercase letter is the reference base aligned to the
    /// query position, `-` a query base with no reference base (insertion, soft clip), and a
    /// lowercase letter a reference base the read skips (deletion, or flanking reference
    /// outside the read). Reference positions count every letter.
    fn ref_map(reference: &str) -> Vec<Option<RefBase>> {
        let genome: Vec<u8> = reference
            .bytes()
            .filter(u8::is_ascii_alphabetic)
            .map(|b| b.to_ascii_uppercase())
            .collect();
        let at = |pos: u64| usize::try_from(pos).ok().and_then(|p| genome.get(p).copied());
        let mut pos = 0u64;
        let mut map = Vec::new();
        for b in reference.bytes() {
            match b {
                b'-' => map.push(None),
                b if b.is_ascii_lowercase() => pos += 1,
                base => {
                    map.push(Some(RefBase {
                        base,
                        prev: pos.checked_sub(1).and_then(at),
                        next: at(pos + 1),
                    }));
                    pos += 1;
                }
            }
        }
        map
    }

    /// SEQ of a read over `reference` (see [`ref_map`]) that matches the reference, with `N` at
    /// `n_positions` and `A` at unaligned query positions.
    fn ref_seq(reference: &str, n_positions: &[usize]) -> Vec<u8> {
        reference
            .bytes()
            .filter(|b| !b.is_ascii_lowercase())
            .enumerate()
            .map(|(i, b)| {
                if n_positions.contains(&i) {
                    b'N'
                } else if b == b'-' {
                    b'A'
                } else {
                    b
                }
            })
            .collect()
    }

    /// Builds a record with SEQ `seq` carrying `counts` as `cu`/`ct` and, for duplex
    /// fixtures, zeroed per-strand arrays per `strands`.
    fn site_record(
        record_flags: u16,
        strands: Option<Strands>,
        counts: Counts,
        seq: &[u8],
    ) -> Vec<u8> {
        let len = counts.cu.len();
        assert_eq!(seq.len(), len, "fixture SEQ and counts differ in length");
        let zeros = vec![0i16; len];
        let mut b = RawSamBuilder::new();
        b.flags(record_flags)
            .ref_id(0)
            .pos(0)
            .mapq(60)
            .cigar_ops(&[u32::try_from(len).expect("short fixture") << 4])
            .sequence(seq)
            .qualities(&vec![30; len]);
        b.add_array_i16(SamTag::CU, counts.cu)
            .add_array_i16(SamTag::CT, counts.ct.unwrap_or(&zeros));
        if matches!(strands, Some(Strands::Both | Strands::AbOnly)) {
            b.add_array_i16(SamTag::AU, &zeros).add_array_i16(SamTag::AT, &zeros);
        }
        if matches!(strands, Some(Strands::Both | Strands::BaOnly)) {
            b.add_array_i16(SamTag::BU, &zeros).add_array_i16(SamTag::BT, &zeros);
        }
        b.build().as_ref().to_vec()
    }

    /// Runs the masks and returns the affected positions with their count: positions that are
    /// newly `N` for a single-strand record, positions whose calls are dropped for a duplex one.
    fn run_site_masks(
        mut record: Vec<u8>,
        consensus_type: ConsensusType,
        masks: &MethylationSiteMasks,
        ref_bases: Option<&[Option<RefBase>]>,
    ) -> (Vec<usize>, usize) {
        let tags = MethylationTags::from_record(&record);
        let before = RawRecordView::new(&record).sequence_vec();
        let sites = match (consensus_type, ref_bases) {
            (ConsensusType::Duplex, _) => MethylationSites::Sequence(None),
            (ConsensusType::SingleStrand, Some(ref_bases)) => {
                MethylationSites::Reference(ref_bases)
            }
            // An unmapped single-strand record has no sites; the caller skips it.
            (ConsensusType::SingleStrand, None) => return (Vec::new(), 0),
        };
        let outcome =
            mask_methylation_sites_raw(&mut record, sites, masks, &tags).expect("masking succeeds");
        if consensus_type == ConsensusType::Duplex {
            assert_eq!(RawRecordView::new(&record).sequence_vec(), before, "duplex SEQ unchanged");
            let count = outcome.dropped_calls.len();
            return (outcome.dropped_calls, count);
        }
        let seq_off = bam_fields::seq_offset(&record);
        let len = RawRecordView::new(&record).l_seq() as usize;
        let n_positions = (0..len).filter(|&i| bam_fields::is_base_n(&record, seq_off, i));
        (n_positions.collect(), outcome.masked)
    }

    /// Simplex: only informative positions are ever masked — reference C for top-strand-type
    /// records (R1/fragment forward, R2 reverse), reference G for the others — at depth < x.
    /// Depths: C 2→2, 5→3, 8→4; G 3→2, 7→3, 9→4. With x = 3 exactly one of each class fails.
    #[rstest]
    #[case::r1_fwd_ref_c(R1_FWD, &[3], &[2])]
    #[case::r1_rev_ref_g(R1_REV, &[3], &[3])]
    #[case::r2_fwd_ref_g(R2_FWD, &[3], &[3])]
    #[case::r2_rev_ref_c(R2_REV, &[3], &[2])]
    #[case::fragment_fwd_ref_c(FRAG_FWD, &[3], &[2])]
    #[case::fragment_rev_ref_g(FRAG_REV, &[3], &[3])]
    #[case::uses_first_value(R1_FWD, &[3, 1, 1], &[2])]
    #[case::zero_threshold_masks_nothing(R1_FWD, &[0], &[])]
    fn test_mask_methylation_sites_simplex_informative_positions(
        #[case] record_flags: u16,
        #[case] thresholds: &[usize],
        #[case] expected: &[usize],
    ) {
        let counts = Counts::depth(&[0, 0, 2, 2, 0, 3, 0, 3, 4, 4, 0]);
        let record = site_record(record_flags, None, counts, &ref_seq("AACGTCAGCGT", &[]));
        let map = ref_map("AACGTCAGCGT");
        let (masked, count) =
            run_site_masks(record, ConsensusType::SingleStrand, &depth(thresholds), Some(&map));
        assert_eq!(masked, expected);
        assert_eq!(count, expected.len());
    }

    /// Simplex: unaligned query bases (insertions, soft clips) are never informative.
    #[test]
    fn test_mask_methylation_sites_simplex_skips_unaligned_positions() {
        let record = site_record(R1_FWD, None, Counts::depth(&[0; 8]), &ref_seq("C-CG--CA", &[]));
        let map = ref_map("C-CG--CA");
        let (masked, _) =
            run_site_masks(record, ConsensusType::SingleStrand, &depth(&[1]), Some(&map));
        assert_eq!(masked, [0, 2, 6]);
    }

    /// Simplex: an informative position whose SEQ base is neither the unconverted nor the
    /// converted base (a C>A or C>G variant) has no methylation call and is never masked; a
    /// C>T variant cannot be told from a conversion and is masked like one.
    #[test]
    fn test_mask_methylation_sites_simplex_third_base_not_masked() {
        let mut seq = ref_seq("AACGTCAGCGT", &[]);
        seq[2] = b'A';
        seq[5] = b'T';
        let record = site_record(R1_FWD, None, Counts::depth(&[0; 11]), &seq);
        let map = ref_map("AACGTCAGCGT");
        let (masked, count) =
            run_site_masks(record, ConsensusType::SingleStrand, &depth(&[1]), Some(&map));
        assert_eq!(masked, [5, 8]);
        assert_eq!(count, 2);
    }

    /// Duplex rules read SEQ and counts only, so an unaligned duplex record is filtered too.
    #[test]
    fn test_mask_methylation_sites_duplex_needs_no_reference() {
        let counts = Counts::depth(&[0, 0, 5, 4, 0, 5, 0, 5, 6, 4, 0]);
        let record = site_record(R1_FWD, Some(Strands::Both), counts, &ref_seq("AACGTCAGCGT", &[]));
        let (dropped, _) = run_site_masks(record, ConsensusType::Duplex, &depth(&[10, 5, 3]), None);
        assert_eq!(dropped, [2, 3]);
    }

    /// Duplex: each `CpG` of the read (a `C` followed in the read by a `G`, whatever the
    /// reference) is gated as a pair — drop both calls unless C + G ≥ x, the better half ≥ y and
    /// the worse half ≥ z — and every other call (non-`CpG` C or G, or a `CpG` half whose
    /// partner is not in the read) must reach y. A half without counts carries no call, so its
    /// partner is tested alone. The baseline (`healthy_passes`) passes (10,5,3) everywhere.
    /// `reference` only builds SEQ here: the rules read SEQ, so a deletion between a C and a G
    /// joins them into a `CpG` and an insertion splits one.
    #[rstest]
    #[case::healthy_passes(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 6, 4, 0, 5, 0, 5, 6, 4, 0]), depth(&[10, 5, 3]), &[])]
    #[case::cpg_fails_total(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 5, 4, 0, 5, 0, 5, 6, 4, 0]), depth(&[10, 5, 3]), &[2, 3])]
    #[case::cpg_fails_best(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 4, 4, 0, 5, 0, 5, 6, 4, 0]), depth(&[8, 5, 3]), &[2, 3])]
    #[case::cpg_fails_worst(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 6, 2, 0, 5, 0, 5, 6, 4, 0]), depth(&[8, 5, 3]), &[2, 3])]
    #[case::better_half_may_be_g(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 3, 6, 0, 5, 0, 5, 6, 4, 0]), depth(&[8, 5, 3]), &[])]
    #[case::worse_half_may_be_c(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 2, 6, 0, 5, 0, 5, 6, 4, 0]), depth(&[8, 5, 3]), &[2, 3])]
    #[case::pairs_are_independent(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 6, 4, 0, 5, 0, 5, 1, 1, 0]), depth(&[10, 5, 3]), &[8, 9])]
    #[case::adjacent_cpgs_do_not_share_g(R1_FWD, "ACGCGA", Counts::depth(&[0, 6, 4, 1, 1, 0]), depth(&[10, 5, 3]), &[3, 4])]
    #[case::one_value_expands_fail(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 3, 2, 0, 5, 0, 5, 6, 4, 0]), depth(&[3]), &[2, 3])]
    #[case::one_value_expands_pass(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 3, 3, 0, 5, 0, 5, 6, 4, 0]), depth(&[3]), &[])]
    #[case::two_values_expand_fail(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 7, 3, 0, 5, 0, 5, 6, 4, 0]), depth(&[10, 4]), &[2, 3])]
    #[case::two_values_expand_pass(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 6, 4, 0, 5, 0, 5, 6, 4, 0]), depth(&[10, 4]), &[])]
    #[case::non_cpg_c_uses_best(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 6, 4, 0, 4, 0, 5, 6, 4, 0]), depth(&[10, 5, 3]), &[5])]
    #[case::non_cpg_g_uses_best(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 6, 4, 0, 5, 0, 4, 6, 4, 0]), depth(&[10, 5, 3]), &[7])]
    #[case::c_at_read_end_is_single_site_pass(R1_FWD, "AACGTCAGC", Counts::depth(&[0, 0, 6, 4, 0, 5, 0, 5, 5]), depth(&[10, 5, 3]), &[])]
    #[case::c_at_read_end_is_single_site_fail(R1_FWD, "AACGTCAGC", Counts::depth(&[0, 0, 6, 4, 0, 5, 0, 5, 4]), depth(&[10, 5, 3]), &[8])]
    #[case::g_at_read_start_is_single_site(R1_FWD, "GTCAGCGT", Counts::depth(&[4, 0, 5, 0, 5, 6, 4, 0]), depth(&[10, 5, 3]), &[0])]
    #[case::cpg_split_by_insertion(R1_FWD, "AAC-GTCAGCGT", Counts::depth(&[0, 0, 6, 0, 4, 0, 5, 0, 5, 6, 4, 0]), depth(&[10, 5, 3]), &[4])]
    #[case::cpg_split_by_soft_clip(R1_FWD, "AAC--------", Counts::depth(&[0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0]), depth(&[10, 5, 3]), &[2])]
    #[case::cpg_g_deleted(R1_FWD, "AACgTCAGCGT", Counts::depth(&[0, 0, 4, 0, 5, 0, 5, 6, 4, 0]), depth(&[10, 5, 3]), &[2])]
    #[case::deletion_joins_a_cpg(R1_FWD, "AACaGTCAGCGT", Counts::depth(&[0, 0, 6, 4, 0, 5, 0, 5, 6, 4, 0]), depth(&[10, 5, 3]), &[])]
    #[case::cpg_g_outside_read(R1_FWD, "AACGTCAGCg", Counts::depth(&[0, 0, 6, 4, 0, 5, 0, 5, 4]), depth(&[10, 5, 3]), &[8])]
    #[case::zero_depth_half_is_no_call(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 0, 6, 0, 5, 0, 5, 6, 4, 0]), depth(&[10, 5, 3]), &[])]
    #[case::zero_depth_half_partner_fails_alone(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 0, 4, 0, 5, 0, 5, 6, 4, 0]), depth(&[10, 5, 3]), &[3])]
    #[case::r1_rev_same_rules(R1_REV, "AACGTCAGCGT", Counts::depth(&[0, 0, 6, 4, 0, 4, 0, 5, 6, 4, 0]), depth(&[10, 5, 3]), &[5])]
    #[case::r2_fwd_same_rules(R2_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 6, 4, 0, 4, 0, 5, 6, 4, 0]), depth(&[10, 5, 3]), &[5])]
    #[case::r2_rev_same_rules(R2_REV, "AACGTCAGCGT", Counts::depth(&[0, 0, 6, 4, 0, 4, 0, 5, 6, 4, 0]), depth(&[10, 5, 3]), &[5])]
    #[case::hemimethylated_passes_depth(R1_FWD, "AACGTCAGCGT", Counts::cu_ct(&[0, 0, 5, 0, 0, 5, 0, 5, 6, 4, 0], &[0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0]), depth(&[10, 5, 3]), &[])]
    #[case::agreement_masks_discordant(R1_FWD, "AACGTCAGCGT", Counts::cu_ct(&[0, 0, 5, 0, 0, 5, 0, 5, 6, 4, 0], &[0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0]), agreement(Some(&[10, 5, 3])), &[2, 3])]
    #[case::agreement_after_depth_failure(R1_FWD, "AACGTCAGCGT", Counts::cu_ct(&[0, 0, 2, 0, 0, 5, 0, 5, 6, 4, 0], &[0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0]), agreement(Some(&[10, 5, 3])), &[2, 3])]
    #[case::agreement_skips_zero_half(R1_FWD, "AACGTCAGCGT", Counts::depth(&[0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0]), agreement(None), &[])]
    #[case::agreement_only_no_single_site_rule(R1_FWD, "AACGTCAGCGT", Counts::cu_ct(&[0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0], &[0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0]), agreement(None), &[2, 3])]
    #[case::agreement_reverse_record(R1_REV, "AACGTCAGCGT", Counts::cu_ct(&[0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0], &[0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0]), agreement(None), &[2, 3])]
    #[case::agreement_tied_half_vs_methylated(R1_FWD, "AACGTCAGCGT", Counts::cu_ct(&[0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0], &[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]), agreement(None), &[])]
    #[case::agreement_tied_half_vs_unmethylated(R1_FWD, "AACGTCAGCGT", Counts::cu_ct(&[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], &[0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0]), agreement(None), &[])]
    fn test_mask_methylation_sites_duplex(
        #[case] record_flags: u16,
        #[case] reference: &str,
        #[case] counts: Counts,
        #[case] masks: MethylationSiteMasks,
        #[case] expected: &[usize],
    ) {
        let seq = ref_seq(reference, &[]);
        let record = site_record(record_flags, Some(Strands::Both), counts, &seq);
        let map = ref_map(reference);
        let (masked, count) = run_site_masks(record, ConsensusType::Duplex, &masks, Some(&map));
        assert_eq!(masked, expected);
        assert_eq!(count, expected.len());
    }

    /// Duplex: a half without counts carries no methylation call (the caller made none there: a
    /// variant, a strand disagreement, or a base the other strand did not confirm), whichever
    /// base won. It is not tested, so its base survives, and its `CpG` partner is tested alone
    /// (at y). A variant that is not `C`/`G` is no site at all. SEQ starts as `AACGTCAGCGT`:
    /// `CpG` at 2-3 and 8-9, non-`CpG` C at 5, G at 7.
    #[rstest]
    #[case::c_to_t_partner_passes(R1_FWD, Strands::Both, &[(2, b'T')], &[], &[0, 0, 0, 5, 0, 5, 0, 5, 6, 4, 0], depth(&[10, 5, 3]), &[])]
    #[case::c_to_t_partner_fails_alone(R1_FWD, Strands::Both, &[(2, b'T')], &[], &[0, 0, 0, 4, 0, 5, 0, 5, 6, 4, 0], depth(&[10, 5, 3]), &[3])]
    #[case::g_to_a_partner_passes(R1_FWD, Strands::Both, &[(3, b'A')], &[], &[0, 0, 5, 0, 0, 5, 0, 5, 6, 4, 0], depth(&[10, 5, 3]), &[])]
    #[case::both_halves_variant(R1_FWD, Strands::Both, &[(2, b'T'), (3, b'A')], &[], &[0, 0, 0, 0, 0, 5, 0, 5, 6, 4, 0], depth(&[10, 5, 3]), &[])]
    #[case::variant_with_n_partner(R1_FWD, Strands::Both, &[(2, b'T')], &[3], &[0, 0, 0, 0, 0, 5, 0, 5, 6, 4, 0], depth(&[10, 5, 3]), &[])]
    #[case::cleared_reference_base_is_no_call(R1_FWD, Strands::Both, &[], &[], &[0, 0, 0, 6, 0, 5, 0, 5, 6, 4, 0], depth(&[10, 5, 3]), &[])]
    #[case::cleared_reference_base_partner_fails_alone(R1_FWD, Strands::Both, &[], &[], &[0, 0, 0, 4, 0, 5, 0, 5, 6, 4, 0], depth(&[10, 5, 3]), &[3])]
    #[case::quality_won_third_base_skipped(R1_FWD, Strands::Both, &[(2, b'G')], &[], &[0, 0, 0, 5, 0, 5, 0, 5, 6, 4, 0], depth(&[10, 5, 3]), &[])]
    #[case::non_cpg_variant_not_tested(R1_FWD, Strands::Both, &[(5, b'T')], &[], &[0, 0, 6, 4, 0, 0, 0, 5, 6, 4, 0], depth(&[10, 5, 3]), &[])]
    #[case::reverse_record(R1_REV, Strands::Both, &[(2, b'T')], &[], &[0, 0, 0, 4, 0, 5, 0, 5, 6, 4, 0], depth(&[10, 5, 3]), &[3])]
    #[case::agreement_ignores_variant_half(R1_FWD, Strands::Both, &[(2, b'T')], &[], &[0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0], agreement(None), &[])]
    #[case::single_strand_third_base_skipped(R1_FWD, Strands::AbOnly, &[(5, b'G')], &[], &[0, 0, 6, 0, 0, 0, 0, 0, 6, 0, 0], depth(&[3, 3, 0]), &[])]
    fn test_mask_methylation_sites_duplex_no_call_halves(
        #[case] record_flags: u16,
        #[case] strands: Strands,
        #[case] seq_overrides: &[(usize, u8)],
        #[case] n_positions: &[usize],
        #[case] depths: &'static [i16],
        #[case] masks: MethylationSiteMasks,
        #[case] expected: &[usize],
    ) {
        let reference = "AACGTCAGCGT";
        let mut seq = ref_seq(reference, n_positions);
        for &(i, base) in seq_overrides {
            seq[i] = base;
        }
        let record = site_record(record_flags, Some(strands), Counts::depth(depths), &seq);
        let map = ref_map(reference);
        let (dropped, count) = run_site_masks(record, ConsensusType::Duplex, &masks, Some(&map));
        assert_eq!(dropped, expected);
        assert_eq!(count, expected.len());
    }

    /// Duplex records with one strand: a `CpG` half without counts is the absent strand's and
    /// counts as 0 (so z = 0 opts in to keeping them), and a site without counts outside a pair
    /// is not tested.
    #[rstest]
    #[case::ab_only_absent_half_fails(Strands::AbOnly, &[3, 3, 1], &[2, 3, 8, 9])]
    #[case::ab_only_z_zero_passes(Strands::AbOnly, &[3, 3, 0], &[])]
    #[case::ba_only_absent_half_fails(Strands::BaOnly, &[3, 3, 1], &[2, 3, 8, 9])]
    #[case::ba_only_z_zero_passes(Strands::BaOnly, &[3, 3, 0], &[])]
    fn test_mask_methylation_sites_duplex_single_strand(
        #[case] strands: Strands,
        #[case] thresholds: &[usize],
        #[case] expected: &[usize],
    ) {
        // C 2, 5, 8 carry depth; the G at 7 (non-CpG, other class) has none and is untouched.
        let counts = Counts::depth(&[0, 0, 6, 0, 0, 5, 0, 0, 6, 0, 0]);
        let record = site_record(R1_FWD, Some(strands), counts, &ref_seq("AACGTCAGCGT", &[]));
        let map = ref_map("AACGTCAGCGT");
        let (masked, count) =
            run_site_masks(record, ConsensusType::Duplex, &depth(thresholds), Some(&map));
        assert_eq!(masked, expected);
        assert_eq!(count, expected.len());
    }

    /// An `N` is no site, so it forms no `CpG` pair: its neighbour is tested alone.
    #[test]
    fn test_mask_methylation_sites_already_n_not_counted() {
        let counts = Counts::depth(&[0, 0, 1, 1, 0, 5, 0, 5, 6, 4, 0]);
        let record =
            site_record(R1_FWD, Some(Strands::Both), counts, &ref_seq("AACGTCAGCGT", &[2]));
        let map = ref_map("AACGTCAGCGT");
        let (masked, count) =
            run_site_masks(record, ConsensusType::Duplex, &depth(&[10, 5, 3]), Some(&map));
        assert_eq!(masked, [3]);
        assert_eq!(count, 1);
    }

    /// Without `cu`/`ct` nothing is informative, from the reference or from SEQ.
    #[test]
    fn test_mask_methylation_sites_nothing_to_evaluate() {
        let mut record = site_record(
            R1_FWD,
            Some(Strands::Both),
            Counts::depth(&[0; 11]),
            &ref_seq("AACGTCAGCGT", &[]),
        );
        let tags = MethylationTags { cu: None, ct: None, has_ab: false, has_ba: false };
        let map = ref_map("AACGTCAGCGT");
        for (sites, masks) in [
            (MethylationSites::Reference(&map), depth(&[3])),
            (MethylationSites::Sequence(None), agreement(Some(&[3]))),
        ] {
            let outcome =
                mask_methylation_sites_raw(&mut record, sites, &masks, &tags).expect("masking");
            assert_eq!(outcome, MethylationSiteMaskOutcome::default(), "{sites:?}");
        }
    }

    /// Strand agreement compares the two strands of a duplex record: asking for it with a
    /// single-strand record's reference sites is an error, not a silent no-op.
    #[test]
    fn test_mask_methylation_sites_agreement_needs_duplex_sites() {
        let mut record = site_record(
            R1_FWD,
            Some(Strands::Both),
            Counts::depth(&[3; 11]),
            &ref_seq("AACGTCAGCGT", &[]),
        );
        let tags = MethylationTags::from_record(&record);
        let map = ref_map("AACGTCAGCGT");
        let result = mask_methylation_sites_raw(
            &mut record,
            MethylationSites::Reference(&map),
            &agreement(None),
            &tags,
        );
        assert!(result.is_err(), "agreement with reference sites must be rejected");
    }

    /// Zeroing the counts at rejected sites touches every count tag present, only there.
    #[test]
    fn test_zero_methylation_counts_at_raw() {
        let mut b = RawSamBuilder::new();
        b.flags(R1_FWD).sequence(b"CGCG").qualities(&[30; 4]);
        b.add_array_i16(SamTag::CU, &[1, 2, 3, 4])
            .add_array_i16(SamTag::CT, &[5, 6, 7, 8])
            .add_array_i16(SamTag::AU, &[1, 0, 3, 0])
            .add_array_i16(SamTag::BT, &[0, 6, 0, 8]);
        let mut record = b.build().as_ref().to_vec();

        zero_methylation_counts_at_raw(&mut record, &[1, 2, 9]);

        let aux = bam_fields::aux_data_slice(&record);
        let values = |tag| {
            let array = bam_fields::find_array_tag(aux, tag).expect("tag kept");
            assert_eq!(array.elem_type, b's', "type kept");
            array.data.chunks_exact(2).map(|b| i16::from_le_bytes([b[0], b[1]])).collect::<Vec<_>>()
        };
        assert_eq!(values(SamTag::CU), vec![1, 0, 0, 4]);
        assert_eq!(values(SamTag::CT), vec![5, 0, 0, 8]);
        assert_eq!(values(SamTag::AU), vec![1, 0, 0, 0]);
        assert_eq!(values(SamTag::BT), vec![0, 0, 0, 8]);
        assert!(bam_fields::find_array_tag(aux, SamTag::AT).is_none(), "absent tags not added");
    }

    /// The single-pass `ConsensusScalarTags` + `_tags` filters must be
    /// byte-identical to the per-tag `find_*` path (`is_duplex_consensus`,
    /// `filter_read`, `filter_duplex_read`) across every tag combination — this
    /// is the correctness contract that lets the hot path drop ~9 aux walks to 1.
    #[test]
    fn consensus_scalar_tags_match_find() {
        // (cd, ce, ad, am, bd, bm, ae, be) as Option, covering: simplex-only,
        // duplex full, duplex with aM/bM fallback, missing pairs, absent-cD, and
        // out-of-order tag layouts.
        type C = (
            Option<i32>,
            Option<f32>,
            Option<i32>,
            Option<i32>,
            Option<i32>,
            Option<i32>,
            Option<f32>,
            Option<f32>,
        );
        let cases: &[C] = &[
            (Some(10), Some(0.01), None, None, None, None, None, None), // simplex
            (Some(10), Some(0.01), Some(6), None, Some(4), None, Some(0.02), Some(0.03)), // duplex full
            (Some(2), Some(0.5), Some(1), None, Some(1), None, Some(0.4), Some(0.4)), // low depth
            (Some(10), Some(0.01), None, Some(5), None, Some(3), Some(0.0), Some(0.0)), // aM/bM fallback
            (Some(10), Some(0.01), Some(6), None, None, None, Some(0.02), None), // only aD (simplex-classed)
            (None, None, None, None, None, None, None, None), // non-consensus (cD/cE absent)
            (Some(100), Some(0.0), Some(50), Some(9), Some(40), Some(7), Some(0.0), Some(0.001)), // both aD+aM present
        ];

        let cc =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.1 };
        let ab =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
        let ba =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.2, max_base_error_rate: 0.1 };

        for (i, &(cd, ce, ad, am, bd, bm, ae, be)) in cases.iter().enumerate() {
            let mut b = RawSamBuilder::new();
            b.ref_id(0).pos(0).mapq(60).cigar_ops(&[4 << 4]).sequence(b"ACGT").qualities(&[30; 4]);
            // Interleave with an unrelated tag so the walk must step over a
            // variable-width entry between targets.
            b.add_int_tag(SamTag::NM, 3);
            if let Some(v) = cd {
                b.add_int_tag(SamTag::CD, v);
            }
            if let Some(v) = ce {
                b.add_float_tag(SamTag::CE, v);
            }
            if let Some(v) = ad {
                b.add_int_tag(SamTag::AD, v);
            }
            if let Some(v) = am {
                b.add_int_tag(SamTag::AM, v);
            }
            if let Some(v) = bd {
                b.add_int_tag(SamTag::BD, v);
            }
            if let Some(v) = bm {
                b.add_int_tag(SamTag::BM, v);
            }
            if let Some(v) = ae {
                b.add_float_tag(SamTag::AE, v);
            }
            if let Some(v) = be {
                b.add_float_tag(SamTag::BE, v);
            }
            let rec = b.build();
            let aux = fgumi_raw_bam::aux_data_slice(rec.as_ref());

            let tags = ConsensusScalarTags::from_aux(aux);

            // Classification parity.
            assert_eq!(tags.is_duplex(), is_duplex_consensus(aux), "case {i}: is_duplex");

            // Simplex filter parity (compare Ok/Err shape + value).
            let want_s = filter_read(aux, &cc);
            let got_s = filter_read_tags(&tags, &cc);
            assert_eq!(want_s.is_err(), got_s.is_err(), "case {i}: filter_read err-shape");
            if let (Ok(w), Ok(g)) = (&want_s, &got_s) {
                assert_eq!(w, g, "case {i}: filter_read result");
            }

            // Duplex filter parity.
            let want_d = filter_duplex_read(aux, &cc, &ab, &ba);
            let got_d = filter_duplex_read_tags(&tags, &cc, &ab, &ba);
            assert_eq!(want_d.is_err(), got_d.is_err(), "case {i}: filter_duplex err-shape");
            if let (Ok(w), Ok(g)) = (&want_d, &got_d) {
                assert_eq!(w, g, "case {i}: filter_duplex result");
            }
        }
    }

    /// A tag value paired with the BAM type it is written as, so a case table can
    /// build deliberately mistyped or duplicated aux entries.
    #[derive(Clone, Copy)]
    enum ScalarTagVal {
        Int(i32),
        Float(f32),
    }

    /// The single-pass path must stay identical to the per-tag `find_*` path even
    /// on **malformed** aux: a present-but-mistyped `aD`/`bD` (which
    /// `is_duplex_consensus` still treats as present), and a duplicated tag whose
    /// first occurrence fails to decode (which `find_tag_position` locks onto).
    /// These are the exact inputs `consensus_scalar_tags_match_find` cannot build
    /// — it only emits well-formed, single, correctly-typed tags — so they pin the
    /// `seen`-bitset first-match/presence semantics against the reference path.
    #[rstest]
    // Duplicate `cD`: first `f`-typed (fails int decode), then a valid int. The
    // per-tag path locks onto the first (undecodable) entry, so both `Err`.
    #[case::dup_cd_bad_first(&[
        (SamTag::CD, ScalarTagVal::Float(2.5)),
        (SamTag::CD, ScalarTagVal::Int(10)),
        (SamTag::CE, ScalarTagVal::Float(0.01)),
    ])]
    // Duplicate `cE`: first int-typed (not `f`), then a valid float.
    #[case::dup_ce_bad_first(&[
        (SamTag::CD, ScalarTagVal::Int(10)),
        (SamTag::CE, ScalarTagVal::Int(99)),
        (SamTag::CE, ScalarTagVal::Float(0.01)),
    ])]
    // `aD`/`bD` present but float-typed: presence -> duplex, but int decode fails.
    #[case::float_typed_ad_bd(&[
        (SamTag::CD, ScalarTagVal::Int(10)),
        (SamTag::CE, ScalarTagVal::Float(0.01)),
        (SamTag::AD, ScalarTagVal::Float(6.0)),
        (SamTag::BD, ScalarTagVal::Float(4.0)),
        (SamTag::AE, ScalarTagVal::Float(0.02)),
        (SamTag::BE, ScalarTagVal::Float(0.03)),
    ])]
    // Only `aD` present and float-typed (simplex either way).
    #[case::float_typed_ad_only(&[
        (SamTag::CD, ScalarTagVal::Int(10)),
        (SamTag::CE, ScalarTagVal::Float(0.01)),
        (SamTag::AD, ScalarTagVal::Float(6.0)),
    ])]
    // Duplicate `aD`: first float-typed (fails int decode), then a valid int; the
    // slot must stay locked to the first occurrence like `find_int_tag(AD)`.
    #[case::dup_ad_bad_first(&[
        (SamTag::CD, ScalarTagVal::Int(10)),
        (SamTag::CE, ScalarTagVal::Float(0.01)),
        (SamTag::AD, ScalarTagVal::Float(9.9)),
        (SamTag::AD, ScalarTagVal::Int(6)),
        (SamTag::BD, ScalarTagVal::Int(4)),
        (SamTag::AE, ScalarTagVal::Float(0.02)),
        (SamTag::BE, ScalarTagVal::Float(0.03)),
    ])]
    fn consensus_scalar_tags_match_find_malformed(#[case] tags: &[(SamTag, ScalarTagVal)]) {
        let cc =
            FilterThresholds { min_reads: 3, max_read_error_rate: 0.1, max_base_error_rate: 0.1 };
        let ab =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
        let ba =
            FilterThresholds { min_reads: 2, max_read_error_rate: 0.2, max_base_error_rate: 0.1 };

        let mut b = RawSamBuilder::new();
        b.ref_id(0).pos(0).mapq(60).cigar_ops(&[4 << 4]).sequence(b"ACGT").qualities(&[30; 4]);
        for &(tag, val) in tags {
            match val {
                ScalarTagVal::Int(v) => b.add_int_tag(tag, v),
                ScalarTagVal::Float(v) => b.add_float_tag(tag, v),
            };
        }
        let rec = b.build();
        let aux = fgumi_raw_bam::aux_data_slice(rec.as_ref());

        let scalar = ConsensusScalarTags::from_aux(aux);

        // Classification parity against the presence-based reference predicate.
        assert_eq!(scalar.is_duplex(), is_duplex_consensus(aux), "is_duplex parity");

        // Simplex read-level filter parity (Err-shape + Ok value).
        let want_s = filter_read(aux, &cc);
        let got_s = filter_read_tags(&scalar, &cc);
        assert_eq!(want_s.is_err(), got_s.is_err(), "filter_read err-shape");
        if let (Ok(w), Ok(g)) = (&want_s, &got_s) {
            assert_eq!(w, g, "filter_read result");
        }

        // Duplex read-level filter parity.
        let want_d = filter_duplex_read(aux, &cc, &ab, &ba);
        let got_d = filter_duplex_read_tags(&scalar, &cc, &ab, &ba);
        assert_eq!(want_d.is_err(), got_d.is_err(), "filter_duplex err-shape");
        if let (Ok(w), Ok(g)) = (&want_d, &got_d) {
            assert_eq!(w, g, "filter_duplex result");
        }
    }

    // -- drop_masked_modifications_raw tests --

    /// Builds a mapped record with genomic SEQ `seq` and the given string/ML tags.
    fn modification_record(record_flags: u16, seq: &[u8], mm: &str, ml: &[u8]) -> Vec<u8> {
        let mut b = RawSamBuilder::new();
        b.flags(record_flags)
            .ref_id(0)
            .pos(0)
            .mapq(60)
            .cigar_ops(&[u32::try_from(seq.len()).expect("short fixture") << 4])
            .sequence(seq)
            .qualities(&vec![30; seq.len()]);
        b.add_string_tag(SamTag::MM, mm.as_bytes())
            .add_array_u8(SamTag::ML, ml)
            .add_int_tag(SamTag::MN, i32::try_from(seq.len()).expect("short fixture"))
            .add_string_tag(SamTag::AM_BASES, mm.as_bytes());
        b.build().as_ref().to_vec()
    }

    fn string_tag(record: &[u8], tag: SamTag) -> Option<String> {
        bam_fields::find_string_tag(bam_fields::aux_data_slice(record), tag)
            .map(|v| String::from_utf8(v.to_vec()).expect("UTF-8 tag"))
    }

    /// MM/ML (and the per-strand `am`) index SEQ in original read orientation, so for a
    /// reverse-mapped record genomic position `i` is MM position `len - 1 - i`. Genomic
    /// `AAGTAAGT` reads `ACTTACTT` as sequenced, with C at 1 and 5 (genomic 6 and 2).
    #[rstest]
    #[case::forward_record(R1_FWD, b"ACTTACTT", 1, "C+m?,0;", &[2])]
    #[case::reverse_record(R1_REV, b"AAGTAAGT", 2, "C+m?,0;", &[1])]
    #[case::untracked_base(R1_REV, b"AAGTAAGT", 0, "C+m?,0,0;", &[1, 2])]
    fn test_drop_masked_modifications_raw(
        #[case] record_flags: u16,
        #[case] seq: &[u8],
        #[case] masked_position: usize,
        #[case] expected_mm: &str,
        #[case] expected_ml: &[u8],
    ) {
        let mut record = modification_record(record_flags, seq, "C+m?,0,0;", &[1, 2]);
        let before = RawRecordView::new(&record).sequence_vec();
        let seq_off = bam_fields::seq_offset(&record);
        bam_fields::mask_base(&mut record, seq_off, masked_position);

        drop_masked_modifications_raw(&mut record, &before);

        let aux = bam_fields::aux_data_slice(&record);
        assert_eq!(string_tag(&record, SamTag::MM).as_deref(), Some(expected_mm), "MM");
        assert_eq!(
            bam_fields::find_array_tag(aux, SamTag::ML)
                .map(|a| bam_fields::array_tag_to_vec_u16(&a)),
            Some(expected_ml.iter().map(|&v| u16::from(v)).collect()),
            "ML"
        );
        assert_eq!(string_tag(&record, SamTag::AM_BASES).as_deref(), Some(expected_mm), "am");
        assert_eq!(bam_fields::find_int_tag(aux, SamTag::MN), Some(8), "MN");
    }

    /// Dropping calls leaves SEQ alone and mirrors genomic positions for reverse-mapped records:
    /// genomic `AAGTAAGT` reads `ACTTACTT` as sequenced, with C at 1 and 5 (genomic 6 and 2).
    #[rstest]
    #[case::forward_first(R1_FWD, b"ACTTACTT", &[1], "C+m?,1;", &[2])]
    #[case::forward_second(R1_FWD, b"ACTTACTT", &[5], "C+m?,0;", &[1])]
    #[case::reverse_genomic_2_is_read_5(R1_REV, b"AAGTAAGT", &[2], "C+m?,0;", &[1])]
    #[case::reverse_genomic_6_is_read_1(R1_REV, b"AAGTAAGT", &[6], "C+m?,1;", &[2])]
    #[case::untracked_position(R1_FWD, b"ACTTACTT", &[0], "C+m?,0,0;", &[1, 2])]
    fn test_drop_modifications_at_raw(
        #[case] record_flags: u16,
        #[case] seq: &[u8],
        #[case] positions: &[usize],
        #[case] expected_mm: &str,
        #[case] expected_ml: &[u8],
    ) {
        let mut record = modification_record(record_flags, seq, "C+m?,0,0;", &[1, 2]);
        drop_modifications_at_raw(&mut record, positions);

        let aux = bam_fields::aux_data_slice(&record);
        assert_eq!(RawRecordView::new(&record).sequence_vec(), seq, "SEQ unchanged");
        assert_eq!(string_tag(&record, SamTag::MM).as_deref(), Some(expected_mm), "MM");
        assert_eq!(
            bam_fields::find_array_tag(aux, SamTag::ML)
                .map(|a| bam_fields::array_tag_to_vec_u16(&a)),
            Some(expected_ml.iter().map(|&v| u16::from(v)).collect()),
            "ML"
        );
        assert_eq!(string_tag(&record, SamTag::AM_BASES).as_deref(), Some(expected_mm), "am");
        assert_eq!(bam_fields::find_int_tag(aux, SamTag::MN), Some(8), "MN");
    }

    /// Tags that no longer fit SEQ cannot be edited and are removed rather than left wrong.
    #[test]
    fn test_drop_masked_modifications_raw_removes_unfixable_tags() {
        let mut record = modification_record(R1_FWD, b"ACTTACTT", "C+m?,0,7;", &[1, 2]);
        let before = RawRecordView::new(&record).sequence_vec();
        let seq_off = bam_fields::seq_offset(&record);
        bam_fields::mask_base(&mut record, seq_off, 1);

        assert!(drop_masked_modifications_raw(&mut record, &before), "removal reported");

        let aux = bam_fields::aux_data_slice(&record);
        for tag in [SamTag::MM, SamTag::ML, SamTag::MN, SamTag::AM_BASES] {
            assert!(bam_fields::find_tag_type(aux, tag).is_none(), "{tag:?} removed");
        }
    }

    /// An `MN` that no longer matches SEQ (SEQ hard-clipped after the tags were written) means no
    /// modification tag fits SEQ: all are removed, and the removal is reported.
    #[rstest]
    #[case::masking(true)]
    #[case::dropping_calls(false)]
    fn test_modification_tags_removed_when_mn_is_stale(#[case] masking: bool) {
        let mut record = modification_record(R1_FWD, b"ACTTACTT", "C+m?,0,0;", &[1, 2]);
        {
            let mut editor = bam_fields::RawTagsEditor::from_vec(&mut record);
            editor.update_int(SamTag::MN, 10);
        }
        let removed = if masking {
            let before = RawRecordView::new(&record).sequence_vec();
            let seq_off = bam_fields::seq_offset(&record);
            bam_fields::mask_base(&mut record, seq_off, 1);
            drop_masked_modifications_raw(&mut record, &before)
        } else {
            drop_modifications_at_raw(&mut record, &[1])
        };
        assert!(removed, "removal reported");
        let aux = bam_fields::aux_data_slice(&record);
        for tag in [SamTag::MM, SamTag::ML, SamTag::MN, SamTag::AM_BASES] {
            assert!(bam_fields::find_tag_type(aux, tag).is_none(), "{tag:?} removed");
        }
    }

    /// Builds a forward R1 record over `ACTTACTT` with `MM:C+m?,0,0;`, MN, and `ml_tag` applied
    /// to the builder (or no ML at all), then masks the C at 1.
    fn masked_record_with_ml(ml_tag: impl FnOnce(&mut RawSamBuilder)) -> Vec<u8> {
        let mut b = RawSamBuilder::new();
        b.flags(R1_FWD)
            .ref_id(0)
            .pos(0)
            .mapq(60)
            .cigar_ops(&[8 << 4])
            .sequence(b"ACTTACTT")
            .qualities(&[30; 8]);
        b.add_string_tag(SamTag::MM, b"C+m?,0,0;").add_int_tag(SamTag::MN, 8);
        ml_tag(&mut b);
        let mut record = b.build().as_ref().to_vec();
        let before = RawRecordView::new(&record).sequence_vec();
        let seq_off = bam_fields::seq_offset(&record);
        bam_fields::mask_base(&mut record, seq_off, 1);
        drop_masked_modifications_raw(&mut record, &before);
        record
    }

    /// MM without ML is rewritten on its own; no (empty) ML is invented.
    #[test]
    fn test_drop_masked_modifications_raw_mm_without_ml() {
        let record = masked_record_with_ml(|_| {});
        let aux = bam_fields::aux_data_slice(&record);
        assert_eq!(string_tag(&record, SamTag::MM).as_deref(), Some("C+m?,0;"), "MM");
        assert!(bam_fields::find_tag_type(aux, SamTag::ML).is_none(), "no ML added");
        assert_eq!(bam_fields::find_int_tag(aux, SamTag::MN), Some(8), "MN");
    }

    /// An ML that is not `B:C` cannot be edited in step with MM, so MM/ML/MN are removed rather
    /// than ML being overwritten with a different type.
    #[test]
    fn test_drop_masked_modifications_raw_non_u8_ml_is_removed() {
        let record = masked_record_with_ml(|b| {
            b.add_array_u16(SamTag::ML, &[1, 2]);
        });
        let aux = bam_fields::aux_data_slice(&record);
        for tag in [SamTag::MM, SamTag::ML, SamTag::MN] {
            assert!(bam_fields::find_tag_type(aux, tag).is_none(), "{tag:?} removed");
        }
    }

    // -- check_conversion_fraction_raw tests --

    /// OB-derived reads (R1 reverse, R2 forward) carry their evidence at reference G, so the
    /// conversion fraction also counts non-`CpG` reference G (previous base not C) and skips
    /// `CpG` G. Reference `TGAACGAG`: non-`CpG` G at 1 and 7, `CpG` G at 5.
    #[rstest]
    #[case::non_cpg_g_unconverted_fails(&[0, 5, 0, 0, 0, 0, 0, 5], &[0; 8], false)]
    #[case::non_cpg_g_converted_passes(&[0; 8], &[0, 5, 0, 0, 0, 0, 0, 5], true)]
    #[case::cpg_g_excluded(&[0, 0, 0, 0, 0, 20, 0, 0], &[0, 5, 0, 0, 0, 0, 0, 5], true)]
    fn test_conversion_fraction_counts_reference_g(
        #[case] cu: &[i16],
        #[case] ct: &[i16],
        #[case] expected_pass: bool,
    ) {
        let raw = {
            let mut b = RawSamBuilder::new();
            b.flags(R1_REV)
                .ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[8 << 4])
                .sequence(b"TAAACAAA")
                .qualities(&[30; 8]);
            b.add_array_i16(SamTag::CU, cu).add_array_i16(SamTag::CT, ct);
            b.build()
        };
        let tags = MethylationTags::from_record(&raw);
        let map = ref_map("TGAACGAG");
        let pass = check_conversion_fraction_raw_with_ref_bases_and_tags(
            &raw,
            0.9,
            MethylationSites::Reference(&map),
            &tags,
            crate::MethylationMode::EmSeq,
        );
        assert_eq!(pass, expected_pass);
    }

    /// `CpG` context comes from the reference, not the read: a C at the last aligned base whose
    /// reference successor is G, or a G at the first aligned base whose reference predecessor is
    /// C, is a `CpG` half and is excluded; a C and G next to each other in the read but separated
    /// by a deletion are both non-`CpG` and counted. Each read has one unconverted site (20
    /// reads, which fails 0.9 unless excluded) and one converted site (50 reads).
    #[rstest]
    #[case::cpg_c_at_last_base_excluded(R1_FWD, "TCTTTCg", &[0, 0, 0, 0, 0, 20], 1, true)]
    #[case::non_cpg_c_at_last_base_counted(R1_FWD, "TCTTTCa", &[0, 0, 0, 0, 0, 20], 1, false)]
    #[case::cpg_g_at_first_base_excluded(R1_REV, "cGAAGA", &[20, 0, 0, 0, 0], 3, true)]
    #[case::deletion_split_counts_both(R1_FWD, "TCaGTT", &[0, 20, 0, 0, 0], 2, false)]
    fn test_conversion_fraction_uses_reference_context(
        #[case] record_flags: u16,
        #[case] reference: &str,
        #[case] cu: &[i16],
        #[case] converted_site: usize,
        #[case] expected_pass: bool,
    ) {
        let mut ct = vec![0i16; cu.len()];
        ct[converted_site] = 50;
        let seq = ref_seq(reference, &[]);
        let raw = {
            let mut b = RawSamBuilder::new();
            b.flags(record_flags)
                .ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[u32::try_from(seq.len()).expect("short fixture") << 4])
                .sequence(&seq)
                .qualities(&vec![30; seq.len()]);
            b.add_array_i16(SamTag::CU, cu).add_array_i16(SamTag::CT, &ct);
            b.build()
        };
        let tags = MethylationTags::from_record(&raw);
        let map = ref_map(reference);
        let pass = check_conversion_fraction_raw_with_ref_bases_and_tags(
            &raw,
            0.9,
            MethylationSites::Reference(&map),
            &tags,
            crate::MethylationMode::EmSeq,
        );
        assert_eq!(pass, expected_pass);
    }

    /// Duplex: SEQ is the molecule's sequence, so non-`CpG` cytosines come from SEQ (a `C` not
    /// followed by `G`, a `G` not preceded by `C`) and need no reference. A cytosine next to `N`
    /// (or the read edge) has unknown context and is not counted. SEQ `CAGTCGC?`: converted
    /// non-`CpG` C at 0 and G at 2, a methylated `CpG` at 4-5 (excluded), and an unconverted C
    /// at 6 that fails 0.9 when its context is known.
    #[rstest]
    #[case::unknown_context_not_counted(b"CAGTCGCN", true)]
    #[case::known_context_counted(b"CAGTCGCA", false)]
    fn test_conversion_fraction_duplex_uses_seq_context(
        #[case] seq: &[u8],
        #[case] expected_pass: bool,
    ) {
        let raw = {
            let mut b = RawSamBuilder::new();
            b.flags(R1_FWD).ref_id(-1).pos(-1).sequence(seq).qualities(&[30; 8]);
            b.add_array_i16(SamTag::CU, &[0, 0, 0, 0, 10, 10, 10, 0])
                .add_array_i16(SamTag::CT, &[10, 0, 10, 0, 0, 0, 0, 0]);
            b.build()
        };
        let tags = MethylationTags::from_record(&raw);
        let pass = check_conversion_fraction_raw_with_ref_bases_and_tags(
            &raw,
            0.9,
            MethylationSites::Sequence(None),
            &tags,
            crate::MethylationMode::EmSeq,
        );
        assert_eq!(pass, expected_pass);
    }

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
        let bases: Vec<Option<u8>> = bases.iter().map(|r| r.map(|r| r.base)).collect();
        assert_eq!(bases, vec![Some(b'A'), Some(b'C'), Some(b'G'), Some(b'T')]);
    }

    /// Each entry carries its reference position and reference neighbours, even across a
    /// deletion and past the ends of the alignment: `ACGTACGTAC` aligned at 1 with `2M1D2M`.
    #[cfg(feature = "simplex")]
    #[test]
    fn test_resolve_ref_bases_positions_and_neighbours() {
        use crate::methylation::tests::TestRef;
        let reference = TestRef::new(&[("chr1", b"ACGTACGTAC")]);
        let ref_names = vec!["chr1".to_string()];
        let raw = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(1)
                .mapq(60)
                .cigar_ops(&[2 << 4, 1 << 4 | 2, 2 << 4])
                .sequence(b"CGAC")
                .qualities(&[30; 4]);
            b.build()
        };
        let bases = resolve_ref_bases_for_record(&raw, &reference, &ref_names).unwrap();
        // Reference positions 1, 2, 4 and 5 (a 1-base deletion skips 3).
        let entry = |base, prev, next| Some(RefBase { base, prev, next });
        assert_eq!(
            bases,
            vec![
                entry(b'C', Some(b'A'), Some(b'G')),
                entry(b'G', Some(b'C'), Some(b'T')),
                entry(b'A', Some(b'T'), Some(b'C')),
                entry(b'C', Some(b'A'), Some(b'G')),
            ]
        );
        assert!(bases[0].unwrap().is_cpg_c() && bases[1].unwrap().is_cpg_g());
        assert!(bases[3].unwrap().is_cpg_c(), "context comes from the reference past the read");
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
        let bases: Vec<Option<u8>> = bases.iter().map(|r| r.map(|r| r.base)).collect();
        assert_eq!(bases, vec![Some(b'A'), Some(b'C'), None, None, Some(b'G'), Some(b'T')]);
    }

    // -- FILT-02: hard-fail when per-read consensus tags (cD/cE) are absent --

    /// Both `cD` and `cE` consensus tags are required: missing either one (or both) must
    /// error (fgbio `FilterConsensusReads.scala:242`); only both-present yields a normal
    /// filter result. `expected` is `None` when the record must be rejected.
    #[rstest]
    #[case::both_absent(false, false, None)]
    #[case::only_cd_present(true, false, None)]
    #[case::only_ce_present(false, true, None)]
    #[case::both_present(true, true, Some(FilterResult::Pass))]
    fn test_filter_read_requires_consensus_tags(
        #[case] with_cd: bool,
        #[case] with_ce: bool,
        #[case] expected: Option<FilterResult>,
    ) {
        let thresholds =
            FilterThresholds { min_reads: 1, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
        let mut b = RawSamBuilder::new();
        b.ref_id(0).pos(0).mapq(0).cigar_ops(&[4 << 4]).sequence(b"ACGT").qualities(&[30; 4]);
        if with_cd {
            b.add_int_tag(SamTag::CD, 20);
        }
        if with_ce {
            b.add_float_tag(SamTag::CE, 0.0);
        }
        let rec = b.build();
        let aux = fgumi_raw_bam::aux_data_slice(rec.as_ref());

        match expected {
            None => assert!(
                filter_read(aux, &thresholds).is_err(),
                "cD={with_cd} cE={with_ce}: absent consensus tag must be rejected, not passed"
            ),
            Some(result) => assert_eq!(filter_read(aux, &thresholds).unwrap(), result),
        }
    }

    // -- FILT-01: mean base quality is over the FULL read (pre-masking), including N bases --

    #[test]
    fn test_mean_base_quality_full_length_includes_all_bases() {
        // 5 bases at q40 and 5 no-call (N) bases at q2.
        let rec = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(0)
                .cigar_ops(&[10 << 4])
                .sequence(b"AAAAANNNNN")
                .qualities(&[40, 40, 40, 40, 40, 2, 2, 2, 2, 2]);
            b.build()
        };
        // fgbio: sum(all quals) / length = (5*40 + 5*2)/10 = 21.0 (FilterConsensusReads.scala:247).
        assert!((mean_base_quality_full_length(rec.as_ref()) - 21.0).abs() < 1e-9);
        // The non-N mean (what the read-level check must NOT use) is 40.0 — proving they differ.
        assert!((compute_read_stats(rec.as_ref()).1 - 40.0).abs() < 1e-9);
    }

    // -- FILT-04: when per-base cd/ce tags are absent, mask only by base quality --

    /// Per-base depth/error masking engages only when BOTH the `cd` and `ce` per-base arrays
    /// are present (`pb = depths != null && errors != null`, fgbio
    /// FilterConsensusReads.scala:281,287-289). With neither — or only one — present, `pb` is
    /// false: fgumi must NOT treat an absent array as depth 0 and mask everything; only the
    /// base-quality mask applies, so exactly the 5 q10 bases (below `--min-base-quality` 30)
    /// are masked. The one-present-one-absent cases guard the `&&` against a regression to
    /// `||` / a single `is_some()`: the deliberately low `cd` depths (1 < `min_reads`=5) would
    /// mask all 10 bases if `pb` were wrongly true.
    #[rstest]
    #[case::neither_tag(false, false)]
    #[case::only_cd_present(true, false)]
    #[case::only_ce_present(false, true)]
    fn test_mask_bases_depth_mask_requires_both_per_base_tags(
        #[case] with_cd: bool,
        #[case] with_ce: bool,
    ) {
        let mut rec = {
            let mut b = RawSamBuilder::new();
            b.ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[10 << 4])
                .sequence(b"AAAAAAAAAA")
                .qualities(&[40, 40, 40, 40, 40, 10, 10, 10, 10, 10]);
            if with_cd {
                b.add_array_i16(SamTag::CD_BASES, &[1; 10]);
            }
            if with_ce {
                b.add_array_i16(SamTag::CE_BASES, &[100; 10]);
            }
            b.build().into_inner()
        };
        let thresholds =
            FilterThresholds { min_reads: 5, max_read_error_rate: 0.05, max_base_error_rate: 0.1 };
        let masked = mask_bases(&mut rec, &thresholds, Some(30)).unwrap();
        assert_eq!(
            masked, 5,
            "with_cd={with_cd} with_ce={with_ce}: pb=false, only the 5 low-quality bases masked \
             (depth/error mask skipped unless BOTH per-base tags present)"
        );
    }

    // -- CONS-02: duplex detection requires BOTH aD and bD (fgbio Umis.isFgbioDuplexConsensus) --

    // fgbio Umis.isFgbioDuplexConsensus = contains(aD) && contains(bD): only both -> duplex.
    #[rstest]
    #[case::neither(false, false, false)]
    #[case::ad_only(true, false, false)]
    #[case::bd_only(false, true, false)]
    #[case::both(true, true, true)]
    fn test_is_duplex_consensus_requires_both_ad_and_bd(
        #[case] with_ad: bool,
        #[case] with_bd: bool,
        #[case] expected: bool,
    ) {
        let mut b = RawSamBuilder::new();
        b.ref_id(0).pos(0).mapq(60).cigar_ops(&[4 << 4]).sequence(b"ACGT").qualities(&[30; 4]);
        if with_ad {
            b.add_int_tag(SamTag::AD, 10);
        }
        if with_bd {
            b.add_int_tag(SamTag::BD, 8);
        }
        let rec = b.build();
        assert_eq!(
            is_duplex_consensus(fgumi_raw_bam::aux_data_slice(rec.as_ref())),
            expected,
            "aD={with_ad} bD={with_bd}: duplex iff both present"
        );
    }

    /// The raw-aux consensus predicates over the `{cD, aD, bD}` presence lattice.
    ///
    /// The asymmetry is deliberate and fgbio-compatible: `aD`/`bD` must *both* be
    /// present for duplex, and a read carrying only one of them counts as simplex
    /// (fgbio routes it to the vanilla filter). That is exactly the property a
    /// future "simplification" would flatten, so every combination is pinned.
    #[rstest]
    #[case::nothing(&[], false, false, false)]
    #[case::cd_only(&[SamTag::CD], true, false, true)]
    #[case::ad_only(&[SamTag::AD], false, false, false)]
    #[case::bd_only(&[SamTag::BD], false, false, false)]
    #[case::ad_and_bd(&[SamTag::AD, SamTag::BD], false, true, true)]
    #[case::cd_and_ad_only(&[SamTag::CD, SamTag::AD], true, false, true)]
    #[case::cd_and_bd_only(&[SamTag::CD, SamTag::BD], true, false, true)]
    #[case::all_three(&[SamTag::CD, SamTag::AD, SamTag::BD], false, true, true)]
    fn raw_aux_consensus_predicates(
        #[case] tags: &[SamTag],
        #[case] simplex: bool,
        #[case] duplex: bool,
        #[case] any: bool,
    ) {
        let mut b = fgumi_raw_bam::SamBuilder::new();
        b.read_name(b"r").sequence(b"A").qualities(&[30]);
        for tag in tags {
            b.add_int_tag(*tag, 1);
        }
        let record = b.build();
        let aux = fgumi_raw_bam::aux_data_slice(record.as_ref());

        assert_eq!(is_simplex_consensus(aux), simplex, "is_simplex_consensus");
        assert_eq!(is_duplex_consensus(aux), duplex, "is_duplex_consensus");
        assert_eq!(is_consensus(aux), any, "is_consensus");
    }

    // -- FILT-05: the masked-bases statistic counts retained primary reads only --

    /// Builds a minimal one-base raw record carrying only the given SAM flag bits.
    fn record_with_flags(flag: u16) -> RawRecord {
        let mut b = RawSamBuilder::new();
        b.ref_id(0).pos(0).mapq(0).flags(flag).cigar_ops(&[1 << 4]).sequence(b"A").qualities(&[30]);
        b.build()
    }

    /// fgbio accumulates `maskedBases` only for the **primary** reads of a **retained**
    /// template (`FilterConsensusReads.scala:207-219`): the increment lives inside
    /// `if (r1.keepRead && r2.keepRead)` and sums only R1/R2, never secondary/supplementary.
    /// [`retained_primary_masked_bases`] must reproduce that: a dropped template contributes
    /// 0, secondary/supplementary reads contribute 0 even when the template is kept, and a
    /// kept template contributes exactly the sum over its primary reads.
    #[rstest]
    #[case::kept_primary(vec![(bam_fields::flags::PAIRED | bam_fields::flags::FIRST_SEGMENT, 7)], true, 7)]
    #[case::dropped_template(vec![(bam_fields::flags::PAIRED | bam_fields::flags::FIRST_SEGMENT, 7)], false, 0)]
    #[case::secondary_excluded(vec![(bam_fields::flags::SECONDARY, 7)], true, 0)]
    #[case::supplementary_excluded(vec![(bam_fields::flags::SUPPLEMENTARY, 7)], true, 0)]
    #[case::primary_pair_sums(vec![
        (bam_fields::flags::PAIRED | bam_fields::flags::FIRST_SEGMENT, 3),
        (bam_fields::flags::PAIRED | bam_fields::flags::LAST_SEGMENT, 2),
    ], true, 5)]
    #[case::supplementary_excluded_from_kept_template(vec![
        (bam_fields::flags::PAIRED | bam_fields::flags::FIRST_SEGMENT, 3),
        (bam_fields::flags::SUPPLEMENTARY, 5),
        (bam_fields::flags::PAIRED | bam_fields::flags::LAST_SEGMENT, 2),
    ], true, 5)]
    #[case::mixed_template_dropped(vec![
        (bam_fields::flags::PAIRED | bam_fields::flags::FIRST_SEGMENT, 3),
        (bam_fields::flags::SUPPLEMENTARY, 5),
        (bam_fields::flags::PAIRED | bam_fields::flags::LAST_SEGMENT, 2),
    ], false, 0)]
    // Boundary: an empty family contributes nothing even when "retained".
    #[case::empty_family(vec![], true, 0)]
    // An unpaired (fragment) primary read is not secondary/supplementary, so a
    // retained single read still contributes its masked bases.
    #[case::unpaired_single_primary(vec![(0, 4)], true, 4)]
    fn test_retained_primary_masked_bases(
        #[case] flags_and_masked: Vec<(u16, u64)>,
        #[case] template_pass: bool,
        #[case] expected: u64,
    ) {
        let records: Vec<RawRecord> =
            flags_and_masked.iter().map(|&(flag, _)| record_with_flags(flag)).collect();
        let masked_by_record: Vec<u64> =
            flags_and_masked.iter().map(|&(_, masked)| masked).collect();

        assert_eq!(
            retained_primary_masked_bases(&records, &masked_by_record, template_pass),
            expected,
            "template_pass={template_pass}, records={flags_and_masked:?}"
        );
    }

    /// The parallel-slice invariant is enforced with `assert_eq!` (not `debug_assert_eq!`)
    /// so a length mismatch cannot silently truncate the masked-base tally in release
    /// builds. `template_pass` is `true` here to prove the panic guards the counting path
    /// itself, not merely the early `!template_pass` return.
    #[rstest]
    #[should_panic(expected = "raw_records and masked_by_record must be parallel")]
    fn test_retained_primary_masked_bases_length_mismatch_panics() {
        let records =
            vec![record_with_flags(bam_fields::flags::PAIRED | bam_fields::flags::FIRST_SEGMENT)];
        let masked_by_record: Vec<u64> = vec![3, 4];
        let _ = retained_primary_masked_bases(&records, &masked_by_record, true);
    }
}
