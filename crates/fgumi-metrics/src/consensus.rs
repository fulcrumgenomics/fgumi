//! Metrics for consensus calling commands.
//!
//! This module provides metrics for tracking consensus calling operations,
//! shared by the `simplex`, `duplex`, and `codec` commands.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::rejection::RejectionReason;

use crate::{Metric, format_float, frac_u64};

/// A key-value-description metric row matching fgbio's `ConsensusKvMetric` format.
///
/// This struct represents a single row in the vertical metrics output format,
/// where each metric is output as a separate row with key, value, and description.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConsensusKvMetric {
    /// The metric key/name
    pub key: String,
    /// The metric value (as a string for TSV output)
    pub value: String,
    /// Human-readable description of the metric
    pub description: String,
}

impl ConsensusKvMetric {
    /// Creates a new key-value-description metric.
    #[must_use]
    pub fn new(key: impl Into<String>, value: String, description: impl Into<String>) -> Self {
        Self { key: key.into(), value, description: description.into() }
    }
}

/// Consensus calling metrics with rejection tracking.
///
/// These metrics track the consensus calling process, including how many reads
/// were accepted, filtered, and the reasons for rejection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConsensusMetrics {
    /// Total input reads processed
    pub total_input_reads: u64,

    /// Number of consensus reads generated
    pub consensus_reads: u64,

    /// Number of input reads filtered out
    pub filtered_reads: u64,

    /// Total number of UMI groups processed
    pub total_umi_groups: u64,

    /// UMI groups that generated consensus
    pub umi_groups_with_consensus: u64,

    /// UMI groups that failed to generate consensus
    pub umi_groups_failed: u64,

    /// Average input reads per consensus read
    #[serde(with = "crate::float")]
    pub avg_input_reads_per_consensus: f64,

    /// Average raw read depth per consensus read
    #[serde(with = "crate::float")]
    pub avg_raw_read_depth: f64,

    /// Minimum raw read depth
    pub min_raw_read_depth: u64,

    /// Maximum raw read depth
    pub max_raw_read_depth: u64,

    /// Reads rejected due to insufficient support
    pub rejected_insufficient_support: u64,

    /// Reads rejected due to minority alignment
    pub rejected_minority_alignment: u64,

    /// Reads rejected due to insufficient strand support
    pub rejected_insufficient_strand_support: u64,

    /// Reads rejected due to low base quality
    pub rejected_low_base_quality: u64,

    /// Reads rejected due to excessive N bases
    pub rejected_excessive_n_bases: u64,

    /// Reads rejected due to no valid alignment
    pub rejected_no_valid_alignment: u64,

    /// Reads rejected due to low mapping quality
    pub rejected_low_mapping_quality: u64,

    /// Reads rejected due to N bases in UMI
    pub rejected_n_bases_in_umi: u64,

    /// Reads rejected due to missing UMI tag
    pub rejected_missing_umi: u64,

    /// Reads rejected due to not passing filter
    pub rejected_not_passing_filter: u64,

    /// Reads rejected due to low mean quality
    pub rejected_low_mean_quality: u64,

    /// Reads rejected due to insufficient min depth
    pub rejected_insufficient_min_depth: u64,

    /// Reads rejected due to excessive error rate
    pub rejected_excessive_error_rate: u64,

    /// Reads rejected due to UMI too short
    pub rejected_umi_too_short: u64,

    /// Reads rejected due to same strand only
    pub rejected_same_strand_only: u64,

    /// Reads rejected because they were unpaired/fragment reads (duplex/codec)
    pub rejected_non_paired_reads: u64,

    /// Reads rejected due to duplicate UMI
    pub rejected_duplicate_umi: u64,

    /// Reads rejected due to orphan consensus (only R1 or R2 had consensus)
    pub rejected_orphan_consensus: u64,

    /// Reads rejected due to zero bases after trimming
    pub rejected_zero_bases_post_trimming: u64,

    /// Reads rejected because the template was not a single primary FR pair (codec)
    pub rejected_not_primary_fr_pair: u64,

    /// Reads rejected because the R1/R2 overlap was too short for CODEC (fgbio
    /// `r1_r2_overlap_too_short`; codec)
    pub rejected_r1_r2_overlap_too_short: u64,

    /// Reads rejected due to an indel error between the top/bottom strands (fgbio
    /// `indel_error_between_strands`; codec)
    pub rejected_indel_error_between_strands: u64,

    /// Reads rejected due to too many errors between the top/bottom strands (fgbio
    /// `high_duplex_disagreement`; codec)
    pub rejected_high_duplex_disagreement: u64,

    /// Reads rejected because overlap clipping failed (fgbio `clip_overlap_failed`;
    /// codec)
    pub rejected_clip_overlap_failed: u64,
}

/// The consensus caller whose statistics are being rendered.
///
/// Used to seed the always-emitted KV rejection rows to match fgbio's per-caller
/// `initializeRejectCounts` (`UmiConsensusCaller.scala:223-224`): the vanilla
/// caller seeds `usedByVanilla` reasons, duplex additionally seeds `usedByDuplex`,
/// and codec additionally seeds `usedByCodec`. Seeded rows are emitted even when
/// their count is zero, exactly as fgbio does.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConsensusCallerKind {
    /// Vanilla / simplex (`CallMolecularConsensusReads`).
    Vanilla,
    /// Duplex (`CallDuplexConsensusReads`).
    Duplex,
    /// Codec (`CallCodecConsensusReads`).
    Codec,
}

/// Reasons fgbio always emits for every caller (`usedByVanilla`), in fgbio order.
const VANILLA_SEEDED: [RejectionReason; 4] = [
    RejectionReason::InsufficientSupport,
    RejectionReason::MinorityAlignment,
    RejectionReason::OrphanConsensus,
    RejectionReason::ZeroBasesPostTrimming,
];

/// Additional reasons the duplex/codec callers always emit (`usedByDuplex`), in fgbio order.
const DUPLEX_EXTRA_SEEDED: [RejectionReason; 3] = [
    RejectionReason::NonPairedReads,
    RejectionReason::SameStrandOnly,
    RejectionReason::DuplicateUmi,
];

/// Additional reasons the codec caller always emits (`usedByCodec`), in fgbio's
/// `RejectionReason.values` order (`UmiConsensusCaller.scala:84-88`) — the four
/// codec-specific reasons followed by `not_primary_fr_pair`.
const CODEC_EXTRA_SEEDED: [RejectionReason; 5] = [
    RejectionReason::R1R2OverlapTooShort,
    RejectionReason::IndelErrorBetweenStrands,
    RejectionReason::HighDuplexDisagreement,
    RejectionReason::ClipOverlapFailed,
    RejectionReason::NotPrimaryFrPair,
];

impl ConsensusCallerKind {
    /// The rejection reasons this caller seeds (always emits), in fgbio's order.
    fn seeded_reasons(self) -> Vec<RejectionReason> {
        let mut reasons = VANILLA_SEEDED.to_vec();
        if matches!(self, Self::Duplex | Self::Codec) {
            reasons.extend(DUPLEX_EXTRA_SEEDED);
        }
        if matches!(self, Self::Codec) {
            reasons.extend(CODEC_EXTRA_SEEDED);
        }
        reasons
    }
}

/// All rejection reasons fgumi tracks, in a stable order (seeded reasons first).
///
/// fgumi's taxonomy is finer-grained than fgbio's; the reasons not in any caller's
/// seeded set (e.g. `LowBaseQuality`) are fgumi-specific and are emitted only when
/// non-zero.
const ALL_REJECTIONS: [RejectionReason; 24] = [
    RejectionReason::InsufficientSupport,
    RejectionReason::MinorityAlignment,
    RejectionReason::OrphanConsensus,
    RejectionReason::ZeroBasesPostTrimming,
    RejectionReason::NonPairedReads,
    RejectionReason::SameStrandOnly,
    RejectionReason::DuplicateUmi,
    RejectionReason::R1R2OverlapTooShort,
    RejectionReason::IndelErrorBetweenStrands,
    RejectionReason::HighDuplexDisagreement,
    RejectionReason::ClipOverlapFailed,
    RejectionReason::NotPrimaryFrPair,
    RejectionReason::InsufficientStrandSupport,
    RejectionReason::LowBaseQuality,
    RejectionReason::ExcessiveNBases,
    RejectionReason::NoValidAlignment,
    RejectionReason::LowMappingQuality,
    RejectionReason::NBasesInUmi,
    RejectionReason::MissingUmi,
    RejectionReason::NotPassingFilter,
    RejectionReason::LowMeanQuality,
    RejectionReason::InsufficientMinDepth,
    RejectionReason::ExcessiveErrorRate,
    RejectionReason::UmiTooShort,
];

/// Returns an iterator over all rejection reasons fgumi tracks.
fn all_rejections() -> impl Iterator<Item = RejectionReason> {
    ALL_REJECTIONS.into_iter()
}

impl ConsensusMetrics {
    /// Creates a new consensus metrics struct with all counts initialized to zero.
    #[must_use]
    pub fn new() -> Self {
        Self {
            total_input_reads: 0,
            consensus_reads: 0,
            filtered_reads: 0,
            total_umi_groups: 0,
            umi_groups_with_consensus: 0,
            umi_groups_failed: 0,
            avg_input_reads_per_consensus: 0.0,
            avg_raw_read_depth: 0.0,
            min_raw_read_depth: 0,
            max_raw_read_depth: 0,
            rejected_insufficient_support: 0,
            rejected_minority_alignment: 0,
            rejected_insufficient_strand_support: 0,
            rejected_low_base_quality: 0,
            rejected_excessive_n_bases: 0,
            rejected_no_valid_alignment: 0,
            rejected_low_mapping_quality: 0,
            rejected_n_bases_in_umi: 0,
            rejected_missing_umi: 0,
            rejected_not_passing_filter: 0,
            rejected_low_mean_quality: 0,
            rejected_insufficient_min_depth: 0,
            rejected_excessive_error_rate: 0,
            rejected_umi_too_short: 0,
            rejected_same_strand_only: 0,
            rejected_non_paired_reads: 0,
            rejected_duplicate_umi: 0,
            rejected_orphan_consensus: 0,
            rejected_zero_bases_post_trimming: 0,
            rejected_not_primary_fr_pair: 0,
            rejected_r1_r2_overlap_too_short: 0,
            rejected_indel_error_between_strands: 0,
            rejected_high_duplex_disagreement: 0,
            rejected_clip_overlap_failed: 0,
        }
    }

    /// Returns the rejection count for a specific reason.
    #[must_use]
    fn rejection_count(&self, reason: RejectionReason) -> u64 {
        match reason {
            RejectionReason::InsufficientSupport => self.rejected_insufficient_support,
            RejectionReason::MinorityAlignment => self.rejected_minority_alignment,
            RejectionReason::InsufficientStrandSupport => self.rejected_insufficient_strand_support,
            RejectionReason::LowBaseQuality => self.rejected_low_base_quality,
            RejectionReason::ExcessiveNBases => self.rejected_excessive_n_bases,
            RejectionReason::NoValidAlignment => self.rejected_no_valid_alignment,
            RejectionReason::LowMappingQuality => self.rejected_low_mapping_quality,
            RejectionReason::NBasesInUmi => self.rejected_n_bases_in_umi,
            RejectionReason::MissingUmi => self.rejected_missing_umi,
            RejectionReason::NotPassingFilter => self.rejected_not_passing_filter,
            RejectionReason::LowMeanQuality => self.rejected_low_mean_quality,
            RejectionReason::InsufficientMinDepth => self.rejected_insufficient_min_depth,
            RejectionReason::ExcessiveErrorRate => self.rejected_excessive_error_rate,
            RejectionReason::UmiTooShort => self.rejected_umi_too_short,
            RejectionReason::SameStrandOnly => self.rejected_same_strand_only,
            RejectionReason::NonPairedReads => self.rejected_non_paired_reads,
            RejectionReason::DuplicateUmi => self.rejected_duplicate_umi,
            RejectionReason::OrphanConsensus => self.rejected_orphan_consensus,
            RejectionReason::ZeroBasesPostTrimming => self.rejected_zero_bases_post_trimming,
            RejectionReason::NotPrimaryFrPair => self.rejected_not_primary_fr_pair,
            RejectionReason::R1R2OverlapTooShort => self.rejected_r1_r2_overlap_too_short,
            RejectionReason::IndelErrorBetweenStrands => self.rejected_indel_error_between_strands,
            RejectionReason::HighDuplexDisagreement => self.rejected_high_duplex_disagreement,
            RejectionReason::ClipOverlapFailed => self.rejected_clip_overlap_failed,
        }
    }

    /// Adds a rejection count for a specific reason.
    ///
    /// This is used to populate rejection statistics from caller stats.
    pub fn add_rejection(&mut self, reason: RejectionReason, count: u64) {
        match reason {
            RejectionReason::InsufficientSupport => self.rejected_insufficient_support += count,
            RejectionReason::MinorityAlignment => self.rejected_minority_alignment += count,
            RejectionReason::InsufficientStrandSupport => {
                self.rejected_insufficient_strand_support += count;
            }
            RejectionReason::LowBaseQuality => self.rejected_low_base_quality += count,
            RejectionReason::ExcessiveNBases => self.rejected_excessive_n_bases += count,
            RejectionReason::NoValidAlignment => self.rejected_no_valid_alignment += count,
            RejectionReason::LowMappingQuality => self.rejected_low_mapping_quality += count,
            RejectionReason::NBasesInUmi => self.rejected_n_bases_in_umi += count,
            RejectionReason::MissingUmi => self.rejected_missing_umi += count,
            RejectionReason::NotPassingFilter => self.rejected_not_passing_filter += count,
            RejectionReason::LowMeanQuality => self.rejected_low_mean_quality += count,
            RejectionReason::InsufficientMinDepth => self.rejected_insufficient_min_depth += count,
            RejectionReason::ExcessiveErrorRate => self.rejected_excessive_error_rate += count,
            RejectionReason::UmiTooShort => self.rejected_umi_too_short += count,
            RejectionReason::SameStrandOnly => self.rejected_same_strand_only += count,
            RejectionReason::NonPairedReads => self.rejected_non_paired_reads += count,
            RejectionReason::DuplicateUmi => self.rejected_duplicate_umi += count,
            RejectionReason::OrphanConsensus => self.rejected_orphan_consensus += count,
            RejectionReason::ZeroBasesPostTrimming => {
                self.rejected_zero_bases_post_trimming += count;
            }
            RejectionReason::NotPrimaryFrPair => self.rejected_not_primary_fr_pair += count,
            RejectionReason::R1R2OverlapTooShort => self.rejected_r1_r2_overlap_too_short += count,
            RejectionReason::IndelErrorBetweenStrands => {
                self.rejected_indel_error_between_strands += count;
            }
            RejectionReason::HighDuplexDisagreement => {
                self.rejected_high_duplex_disagreement += count;
            }
            RejectionReason::ClipOverlapFailed => self.rejected_clip_overlap_failed += count,
        }
    }

    /// Returns the total number of rejections across all reasons.
    #[must_use]
    pub fn total_rejections(&self) -> u64 {
        all_rejections().map(|r| self.rejection_count(r)).sum()
    }

    /// Returns a map of rejection reasons to their counts.
    ///
    /// Only includes rejection reasons with non-zero counts.
    #[must_use]
    pub fn rejection_summary(&self) -> HashMap<RejectionReason, u64> {
        all_rejections()
            .filter_map(|reason| {
                let count = self.rejection_count(reason);
                if count > 0 { Some((reason, count)) } else { None }
            })
            .collect()
    }

    /// Converts metrics to fgbio-compatible key-value-description format.
    ///
    /// Returns a vector of `ConsensusKvMetric` rows matching fgbio's vertical
    /// metrics output format (`UmiConsensusCaller.statistics`). The `kind`
    /// selects which rejection rows are seeded (always emitted, even at zero),
    /// matching fgbio's per-caller `initializeRejectCounts`: vanilla emits four
    /// rows, duplex additionally emits `non_paired_reads`, `single_strand_only`,
    /// and `potential_umi_collision`, and codec additionally emits the four
    /// codec-specific reasons (`r1_r2_overlap_too_short`,
    /// `indel_error_between_strands`, `high_duplex_disagreement`,
    /// `clip_overlap_failed`) plus `not_primary_fr_pair`.
    /// fgumi-specific finer-grained reasons are emitted only when non-zero.
    #[must_use]
    pub fn to_kv_metrics(&self, kind: ConsensusCallerKind) -> Vec<ConsensusKvMetric> {
        let mut metrics = Vec::new();

        let raw_reads_used = self.total_input_reads.saturating_sub(self.filtered_reads);
        let frac_used = frac_u64(raw_reads_used, self.total_input_reads);

        // Core pipeline metrics
        metrics.push(ConsensusKvMetric::new(
            "raw_reads_considered",
            self.total_input_reads.to_string(),
            "Total raw reads considered from input file",
        ));
        metrics.push(ConsensusKvMetric::new(
            "raw_reads_rejected",
            self.filtered_reads.to_string(),
            "Total number of raw reads rejected before consensus calling",
        ));
        metrics.push(ConsensusKvMetric::new(
            "raw_reads_used",
            raw_reads_used.to_string(),
            "Total count of raw reads used in consensus reads",
        ));
        metrics.push(ConsensusKvMetric::new(
            "frac_raw_reads_used",
            format_float(frac_used),
            "Fraction of raw reads used in consensus reads",
        ));

        // Seeded rejection reasons for this caller: always emitted, in fgbio order.
        let seeded = kind.seeded_reasons();
        for reason in &seeded {
            metrics.push(ConsensusKvMetric::new(
                reason.tsv_key(),
                self.rejection_count(*reason).to_string(),
                reason.kv_description(),
            ));
        }

        // fgumi-specific reasons not in the seeded set: emit only when non-zero.
        for reason in all_rejections() {
            if !seeded.contains(&reason) {
                let count = self.rejection_count(reason);
                if count > 0 {
                    metrics.push(ConsensusKvMetric::new(
                        reason.tsv_key(),
                        count.to_string(),
                        reason.kv_description(),
                    ));
                }
            }
        }

        // Final output metric
        metrics.push(ConsensusKvMetric::new(
            "consensus_reads_emitted",
            self.consensus_reads.to_string(),
            "Total number of consensus reads (R1+R2=2) emitted.",
        ));

        metrics
    }
}

impl Default for ConsensusMetrics {
    fn default() -> Self {
        Self::new()
    }
}

impl Metric for ConsensusMetrics {
    fn metric_name() -> &'static str {
        "consensus"
    }
}

impl crate::ProcessingMetrics for ConsensusMetrics {
    fn total_input(&self) -> u64 {
        self.total_input_reads
    }

    fn total_output(&self) -> u64 {
        self.consensus_reads
    }

    fn total_filtered(&self) -> u64 {
        self.filtered_reads
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_consensus_metrics_new() {
        let metrics = ConsensusMetrics::new();
        assert_eq!(metrics.total_input_reads, 0);
        assert_eq!(metrics.consensus_reads, 0);
        assert_eq!(metrics.filtered_reads, 0);
    }

    #[test]
    fn test_consensus_metrics_rejection_summary() {
        let mut metrics = ConsensusMetrics::new();
        metrics.rejected_insufficient_support = 10;
        metrics.rejected_low_base_quality = 5;
        metrics.rejected_excessive_n_bases = 3;

        let summary = metrics.rejection_summary();
        assert_eq!(summary.len(), 3);
        assert_eq!(summary.get(&RejectionReason::InsufficientSupport), Some(&10));
        assert_eq!(summary.get(&RejectionReason::LowBaseQuality), Some(&5));
        assert_eq!(summary.get(&RejectionReason::ExcessiveNBases), Some(&3));
    }

    #[test]
    fn test_consensus_metrics_total_rejections() {
        let mut metrics = ConsensusMetrics::new();
        metrics.rejected_insufficient_support = 1;
        metrics.rejected_minority_alignment = 2;
        metrics.rejected_low_base_quality = 3;
        metrics.rejected_excessive_n_bases = 4;

        assert_eq!(metrics.total_rejections(), 10);
    }

    #[test]
    fn test_consensus_metrics_default() {
        let metrics = ConsensusMetrics::default();
        assert_eq!(metrics.total_input_reads, 0);
        assert_eq!(metrics.consensus_reads, 0);
    }

    #[test]
    fn test_consensus_metrics_rejection_summary_all_types() {
        let mut metrics = ConsensusMetrics::new();
        metrics.rejected_insufficient_support = 1;
        metrics.rejected_minority_alignment = 2;
        metrics.rejected_insufficient_strand_support = 3;
        metrics.rejected_low_base_quality = 4;
        metrics.rejected_excessive_n_bases = 5;
        metrics.rejected_no_valid_alignment = 6;
        metrics.rejected_low_mapping_quality = 7;
        metrics.rejected_n_bases_in_umi = 8;
        metrics.rejected_not_passing_filter = 9;
        metrics.rejected_low_mean_quality = 10;
        metrics.rejected_insufficient_min_depth = 11;
        metrics.rejected_excessive_error_rate = 12;
        metrics.rejected_umi_too_short = 13;
        metrics.rejected_same_strand_only = 14;
        metrics.rejected_duplicate_umi = 15;

        let summary = metrics.rejection_summary();
        assert_eq!(summary.len(), 15);
        assert_eq!(summary.get(&RejectionReason::InsufficientSupport), Some(&1));
        assert_eq!(summary.get(&RejectionReason::MinorityAlignment), Some(&2));
        assert_eq!(summary.get(&RejectionReason::InsufficientStrandSupport), Some(&3));
        assert_eq!(summary.get(&RejectionReason::LowBaseQuality), Some(&4));
        assert_eq!(summary.get(&RejectionReason::ExcessiveNBases), Some(&5));
        assert_eq!(summary.get(&RejectionReason::NoValidAlignment), Some(&6));
        assert_eq!(summary.get(&RejectionReason::LowMappingQuality), Some(&7));
        assert_eq!(summary.get(&RejectionReason::NBasesInUmi), Some(&8));
        assert_eq!(summary.get(&RejectionReason::NotPassingFilter), Some(&9));
        assert_eq!(summary.get(&RejectionReason::LowMeanQuality), Some(&10));
        assert_eq!(summary.get(&RejectionReason::InsufficientMinDepth), Some(&11));
        assert_eq!(summary.get(&RejectionReason::ExcessiveErrorRate), Some(&12));
        assert_eq!(summary.get(&RejectionReason::UmiTooShort), Some(&13));
        assert_eq!(summary.get(&RejectionReason::SameStrandOnly), Some(&14));
        assert_eq!(summary.get(&RejectionReason::DuplicateUmi), Some(&15));
    }

    #[test]
    fn test_consensus_metrics_rejection_summary_empty() {
        let metrics = ConsensusMetrics::new();
        let summary = metrics.rejection_summary();
        assert!(summary.is_empty());
    }

    #[test]
    fn test_consensus_metrics_total_rejections_all() {
        let mut metrics = ConsensusMetrics::new();
        // Set all rejection counts to 1
        metrics.rejected_insufficient_support = 1;
        metrics.rejected_minority_alignment = 1;
        metrics.rejected_insufficient_strand_support = 1;
        metrics.rejected_low_base_quality = 1;
        metrics.rejected_excessive_n_bases = 1;
        metrics.rejected_no_valid_alignment = 1;
        metrics.rejected_low_mapping_quality = 1;
        metrics.rejected_n_bases_in_umi = 1;
        metrics.rejected_not_passing_filter = 1;
        metrics.rejected_low_mean_quality = 1;
        metrics.rejected_insufficient_min_depth = 1;
        metrics.rejected_excessive_error_rate = 1;
        metrics.rejected_umi_too_short = 1;
        metrics.rejected_same_strand_only = 1;
        metrics.rejected_duplicate_umi = 1;

        assert_eq!(metrics.total_rejections(), 15);
    }

    #[test]
    fn test_add_rejection_round_trips_every_reason() {
        use strum::IntoEnumIterator;

        // `add_rejection(reason, n)` must route to a counter that `rejection_count(reason)` reads
        // back as `n`, and must not touch any other reason's counter. Iterating every variant
        // (including the codec-specific reasons) keeps this exhaustive as new reasons are added.
        for reason in RejectionReason::iter() {
            let mut metrics = ConsensusMetrics::new();
            metrics.add_rejection(reason, 7);
            assert_eq!(
                metrics.rejection_count(reason),
                7,
                "add_rejection/rejection_count disagree for {reason:?}"
            );
            for other in RejectionReason::iter() {
                if other != reason {
                    assert_eq!(
                        metrics.rejection_count(other),
                        0,
                        "adding {reason:?} leaked into {other:?}'s counter"
                    );
                }
            }
        }
    }

    #[test]
    fn test_all_rejections_is_exhaustive() {
        use std::collections::HashSet;
        use strum::IntoEnumIterator;

        // `ALL_REJECTIONS` is hand-maintained and drives `all_rejections()`, which
        // `total_rejections()`, `rejection_summary()`, and the fgumi-specific seeding pass in
        // `to_kv_metrics()` all iterate. Unlike the exhaustive `match` in
        // `add_rejection`/`rejection_count`, nothing forces this list to stay in sync when a new
        // `RejectionReason` variant is added — a forgotten entry would silently under-report that
        // reason forever. Assert it covers exactly the variant set.
        let listed: HashSet<_> = ALL_REJECTIONS.iter().copied().collect();
        let all: HashSet<_> = RejectionReason::iter().collect();
        assert_eq!(
            listed, all,
            "ALL_REJECTIONS is out of sync with RejectionReason variants; \
             total_rejections()/rejection_summary()/to_kv_metrics() would silently miss the diff"
        );
        assert_eq!(
            ALL_REJECTIONS.len(),
            listed.len(),
            "ALL_REJECTIONS contains a duplicate rejection reason"
        );
    }

    #[test]
    fn test_metric_trait_impl() {
        assert_eq!(ConsensusMetrics::metric_name(), "consensus");
    }

    #[test]
    fn test_to_kv_metrics_basic() {
        let mut metrics = ConsensusMetrics::new();
        metrics.total_input_reads = 1000;
        metrics.filtered_reads = 100;
        metrics.consensus_reads = 450;
        metrics.rejected_insufficient_support = 50;
        metrics.rejected_minority_alignment = 50;

        let kv_metrics = metrics.to_kv_metrics(ConsensusCallerKind::Vanilla);

        // Should have core metrics
        let raw_reads = kv_metrics.iter().find(|m| m.key == "raw_reads_considered");
        assert!(raw_reads.is_some());
        assert_eq!(raw_reads.expect("raw_reads_considered metric should be present").value, "1000");
        assert_eq!(
            raw_reads.expect("raw_reads_considered metric should be present").description,
            "Total raw reads considered from input file"
        );

        let rejected = kv_metrics.iter().find(|m| m.key == "raw_reads_rejected");
        assert!(rejected.is_some());
        assert_eq!(rejected.expect("raw_reads_rejected metric should be present").value, "100");

        let used = kv_metrics.iter().find(|m| m.key == "raw_reads_used");
        assert!(used.is_some());
        assert_eq!(used.expect("raw_reads_used metric should be present").value, "900");

        let frac = kv_metrics.iter().find(|m| m.key == "frac_raw_reads_used");
        assert!(frac.is_some());
        assert_eq!(frac.expect("frac_raw_reads_used metric should be present").value, "0.900000");

        let consensus = kv_metrics.iter().find(|m| m.key == "consensus_reads_emitted");
        assert!(consensus.is_some());
        assert_eq!(
            consensus.expect("consensus_reads_emitted metric should be present").value,
            "450"
        );

        // Should have rejection reasons
        let insuff =
            kv_metrics.iter().find(|m| m.key == "raw_reads_rejected_for_insufficient_support");
        assert!(insuff.is_some());
        assert_eq!(
            insuff.expect("insufficient_support rejection metric should be present").value,
            "50"
        );

        let minority =
            kv_metrics.iter().find(|m| m.key == "raw_reads_rejected_for_minority_alignment");
        assert!(minority.is_some());
        assert_eq!(
            minority.expect("minority_alignment rejection metric should be present").value,
            "50"
        );
    }

    #[test]
    fn test_to_kv_metrics_zero_reads() {
        let metrics = ConsensusMetrics::new();
        let kv_metrics = metrics.to_kv_metrics(ConsensusCallerKind::Vanilla);

        // frac_raw_reads_used should be 0.0 when total_input_reads is 0
        let frac = kv_metrics.iter().find(|m| m.key == "frac_raw_reads_used");
        assert!(frac.is_some());
        assert_eq!(frac.expect("frac_raw_reads_used metric should be present").value, "0.000000");
    }

    #[test]
    fn test_to_kv_metrics_only_nonzero_optional_rejections() {
        let mut metrics = ConsensusMetrics::new();
        metrics.total_input_reads = 100;
        metrics.rejected_low_base_quality = 5; // This is an optional rejection

        let kv_metrics = metrics.to_kv_metrics(ConsensusCallerKind::Vanilla);

        // Should have low_base_quality since it's non-zero
        let lbq = kv_metrics.iter().find(|m| m.key == "raw_reads_rejected_for_low_base_quality");
        assert!(lbq.is_some());
        assert_eq!(lbq.expect("low_base_quality rejection metric should be present").value, "5");

        // Should NOT have excessive_n_bases since it's zero
        let en = kv_metrics.iter().find(|m| m.key == "raw_reads_rejected_for_excessive_n_bases");
        assert!(en.is_none());
    }

    /// R2-MET-05: the duplex caller seeds (always emits) `non_paired_reads`,
    /// `single_strand_only`, and `potential_umi_collision` even at zero, matching
    /// fgbio's `initializeRejectCounts(usedByVanilla || usedByDuplex)`. The vanilla
    /// caller emits none of them at zero.
    ///
    /// Assert the *complete ordered* key sequence for both callers (not just presence), so
    /// row-order drift is caught. The two expected sequences also encode the semantics
    /// checked before: the duplex-only rows appear (in fgbio order) for Duplex and are absent
    /// for Vanilla.
    #[test]
    fn test_to_kv_metrics_duplex_seeds_duplex_rows() {
        // Vanilla: four core metrics, the four seeded vanilla rejections in fgbio order, then
        // the final emitted count. No duplex-only rows.
        const EXPECTED_VANILLA: &[&str] = &[
            "raw_reads_considered",
            "raw_reads_rejected",
            "raw_reads_used",
            "frac_raw_reads_used",
            "raw_reads_rejected_for_insufficient_support",
            "raw_reads_rejected_for_minority_alignment",
            "raw_reads_rejected_for_orphan_consensus",
            "raw_reads_rejected_for_zero_bases_post_trimming",
            "consensus_reads_emitted",
        ];
        // Duplex: the same sequence plus the three duplex-only rows seeded (in fgbio order)
        // immediately after the vanilla rejections and before the final emitted count.
        const EXPECTED_DUPLEX: &[&str] = &[
            "raw_reads_considered",
            "raw_reads_rejected",
            "raw_reads_used",
            "frac_raw_reads_used",
            "raw_reads_rejected_for_insufficient_support",
            "raw_reads_rejected_for_minority_alignment",
            "raw_reads_rejected_for_orphan_consensus",
            "raw_reads_rejected_for_zero_bases_post_trimming",
            "raw_reads_rejected_for_non_paired_reads",
            "raw_reads_rejected_for_single_strand_only",
            "raw_reads_rejected_for_potential_umi_collision",
            "consensus_reads_emitted",
        ];

        let metrics = ConsensusMetrics::new();

        let vanilla = metrics.to_kv_metrics(ConsensusCallerKind::Vanilla);
        let vanilla_keys: Vec<&str> = vanilla.iter().map(|m| m.key.as_str()).collect();
        assert_eq!(vanilla_keys.as_slice(), EXPECTED_VANILLA, "vanilla KV row order drifted");

        let duplex = metrics.to_kv_metrics(ConsensusCallerKind::Duplex);
        let duplex_keys: Vec<&str> = duplex.iter().map(|m| m.key.as_str()).collect();
        assert_eq!(duplex_keys.as_slice(), EXPECTED_DUPLEX, "duplex KV row order drifted");
    }

    /// R2-MET-05: the codec caller additionally seeds `not_primary_fr_pair` and
    /// the four codec-specific `usedByCodec` reasons on top of the duplex rows,
    /// matching fgbio's `CodecConsensusCaller.initializeRejectCounts`.
    #[test]
    fn test_to_kv_metrics_codec_seeds_codec_row() {
        let metrics = ConsensusMetrics::new();
        let codec = metrics.to_kv_metrics(ConsensusCallerKind::Codec);
        // All codec-only rows (`usedByCodec` && !`usedByDuplex`) are seeded, at zero.
        let codec_only_keys = [
            "raw_reads_rejected_for_r1_r2_overlap_too_short",
            "raw_reads_rejected_for_indel_error_between_strands",
            "raw_reads_rejected_for_high_duplex_disagreement",
            "raw_reads_rejected_for_clip_overlap_failed",
            "raw_reads_rejected_for_not_primary_fr_pair",
        ];
        for key in codec_only_keys {
            let rows: Vec<_> = codec.iter().filter(|m| m.key == key).collect();
            assert_eq!(rows.len(), 1, "codec KV metrics must emit {key} exactly once");
            assert_eq!(rows[0].value, "0", "{key} must be seeded at zero");
        }
        // The codec-only rows must appear in fgbio's `RejectionReason.values` order
        // (`UmiConsensusCaller.scala:84-88`: the four codec reasons then
        // `not_primary_fr_pair`), not merely be present.
        let emitted_codec_order: Vec<&str> =
            codec.iter().map(|m| m.key.as_str()).filter(|k| codec_only_keys.contains(k)).collect();
        assert_eq!(
            emitted_codec_order, codec_only_keys,
            "codec-only seeded rows must be emitted in fgbio order"
        );
        // ...and it still seeds the duplex rows.
        assert!(codec.iter().any(|m| m.key == "raw_reads_rejected_for_non_paired_reads"));
        // Duplex must NOT seed the codec-only rows.
        let duplex = metrics.to_kv_metrics(ConsensusCallerKind::Duplex);
        for key in codec_only_keys {
            assert!(
                !duplex.iter().any(|m| m.key == key),
                "duplex KV metrics must not seed the codec-only row {key} at zero"
            );
        }
    }

    /// R2-MET-05: seeded duplex rows carry fgbio's exact keys and descriptions.
    #[test]
    fn test_duplex_seeded_row_descriptions_match_fgbio() {
        let metrics = ConsensusMetrics::new();
        let kv = metrics.to_kv_metrics(ConsensusCallerKind::Duplex);
        let find = |key: &str| kv.iter().find(|m| m.key == key).map(|m| m.description.clone());
        assert_eq!(
            find("raw_reads_rejected_for_non_paired_reads").as_deref(),
            Some("Unpaired/fragment reads not supported by Duplex caller")
        );
        assert_eq!(
            find("raw_reads_rejected_for_single_strand_only").as_deref(),
            Some("Only Generating One Strand of Duplex Consensus")
        );
        assert_eq!(
            find("raw_reads_rejected_for_potential_umi_collision").as_deref(),
            Some("Potential collision between independent duplex molecules")
        );
    }

    #[test]
    fn test_consensus_kv_metric_new() {
        let kv = ConsensusKvMetric::new("test_key", 42u64.to_string(), "A test metric");
        assert_eq!(kv.key, "test_key");
        assert_eq!(kv.value, "42");
        assert_eq!(kv.description, "A test metric");
    }

    #[test]
    fn test_nan_infinity_serialization() -> anyhow::Result<()> {
        use crate::writer::{read_metrics, write_metrics};
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new()?;
        let path = temp_file.path();

        // Both averaged-depth fields serialize through `crate::float`, so non-finite values
        // (0/0 → NaN, n/0 → Infinity) must round-trip through fgumi's own writer/reader the
        // same way `correct.rs::test_nan_infinity_serialization` pins them for the sibling
        // metric type. Plain ryu `inf`/`-inf` tokens would fail fgbio's `Double.parseDouble`.
        let mut metrics = ConsensusMetrics::new();
        metrics.avg_input_reads_per_consensus = f64::NAN;
        metrics.avg_raw_read_depth = f64::INFINITY;

        write_metrics(path, &[metrics], "consensus")?;
        let read_back: Vec<ConsensusMetrics> = read_metrics(path, "consensus")?;

        assert_eq!(read_back.len(), 1);
        assert!(read_back[0].avg_input_reads_per_consensus.is_nan());
        assert!(read_back[0].avg_raw_read_depth.is_infinite());
        assert!(read_back[0].avg_raw_read_depth > 0.0); // Positive infinity

        Ok(())
    }
}
