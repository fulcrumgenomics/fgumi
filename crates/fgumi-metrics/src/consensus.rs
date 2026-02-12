//! Metrics for consensus calling commands.
//!
//! This module provides metrics for tracking consensus calling operations,
//! shared by the `simplex`, `duplex`, and `codec` commands.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::rejection::RejectionReason;

use crate::{Metric, format_float};

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
    pub avg_input_reads_per_consensus: f64,

    /// Average raw read depth per consensus read
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

    /// Reads rejected due to duplicate UMI
    pub rejected_duplicate_umi: u64,

    /// Reads rejected due to orphan consensus (only R1 or R2 had consensus)
    pub rejected_orphan_consensus: u64,

    /// Reads rejected due to zero bases after trimming
    pub rejected_zero_bases_post_trimming: u64,
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
            rejected_duplicate_umi: 0,
            rejected_orphan_consensus: 0,
            rejected_zero_bases_post_trimming: 0,
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
            RejectionReason::DuplicateUmi => self.rejected_duplicate_umi += count,
            RejectionReason::OrphanConsensus => self.rejected_orphan_consensus += count,
            RejectionReason::ZeroBasesPostTrimming => {
                self.rejected_zero_bases_post_trimming += count;
            }
        }
    }

    /// Returns the total number of rejections across all reasons.
    #[must_use]
    pub fn total_rejections(&self) -> u64 {
        self.rejected_insufficient_support
            + self.rejected_minority_alignment
            + self.rejected_insufficient_strand_support
            + self.rejected_low_base_quality
            + self.rejected_excessive_n_bases
            + self.rejected_no_valid_alignment
            + self.rejected_low_mapping_quality
            + self.rejected_n_bases_in_umi
            + self.rejected_missing_umi
            + self.rejected_not_passing_filter
            + self.rejected_low_mean_quality
            + self.rejected_insufficient_min_depth
            + self.rejected_excessive_error_rate
            + self.rejected_umi_too_short
            + self.rejected_same_strand_only
            + self.rejected_duplicate_umi
            + self.rejected_orphan_consensus
            + self.rejected_zero_bases_post_trimming
    }

    /// Returns a map of rejection reasons to their counts.
    ///
    /// Only includes rejection reasons with non-zero counts.
    #[must_use]
    pub fn rejection_summary(&self) -> HashMap<RejectionReason, u64> {
        let mut summary = HashMap::new();

        if self.rejected_insufficient_support > 0 {
            summary
                .insert(RejectionReason::InsufficientSupport, self.rejected_insufficient_support);
        }
        if self.rejected_minority_alignment > 0 {
            summary.insert(RejectionReason::MinorityAlignment, self.rejected_minority_alignment);
        }
        if self.rejected_insufficient_strand_support > 0 {
            summary.insert(
                RejectionReason::InsufficientStrandSupport,
                self.rejected_insufficient_strand_support,
            );
        }
        if self.rejected_low_base_quality > 0 {
            summary.insert(RejectionReason::LowBaseQuality, self.rejected_low_base_quality);
        }
        if self.rejected_excessive_n_bases > 0 {
            summary.insert(RejectionReason::ExcessiveNBases, self.rejected_excessive_n_bases);
        }
        if self.rejected_no_valid_alignment > 0 {
            summary.insert(RejectionReason::NoValidAlignment, self.rejected_no_valid_alignment);
        }
        if self.rejected_low_mapping_quality > 0 {
            summary.insert(RejectionReason::LowMappingQuality, self.rejected_low_mapping_quality);
        }
        if self.rejected_n_bases_in_umi > 0 {
            summary.insert(RejectionReason::NBasesInUmi, self.rejected_n_bases_in_umi);
        }
        if self.rejected_missing_umi > 0 {
            summary.insert(RejectionReason::MissingUmi, self.rejected_missing_umi);
        }
        if self.rejected_not_passing_filter > 0 {
            summary.insert(RejectionReason::NotPassingFilter, self.rejected_not_passing_filter);
        }
        if self.rejected_low_mean_quality > 0 {
            summary.insert(RejectionReason::LowMeanQuality, self.rejected_low_mean_quality);
        }
        if self.rejected_insufficient_min_depth > 0 {
            summary.insert(
                RejectionReason::InsufficientMinDepth,
                self.rejected_insufficient_min_depth,
            );
        }
        if self.rejected_excessive_error_rate > 0 {
            summary.insert(RejectionReason::ExcessiveErrorRate, self.rejected_excessive_error_rate);
        }
        if self.rejected_umi_too_short > 0 {
            summary.insert(RejectionReason::UmiTooShort, self.rejected_umi_too_short);
        }
        if self.rejected_same_strand_only > 0 {
            summary.insert(RejectionReason::SameStrandOnly, self.rejected_same_strand_only);
        }
        if self.rejected_duplicate_umi > 0 {
            summary.insert(RejectionReason::DuplicateUmi, self.rejected_duplicate_umi);
        }
        if self.rejected_orphan_consensus > 0 {
            summary.insert(RejectionReason::OrphanConsensus, self.rejected_orphan_consensus);
        }
        if self.rejected_zero_bases_post_trimming > 0 {
            summary.insert(
                RejectionReason::ZeroBasesPostTrimming,
                self.rejected_zero_bases_post_trimming,
            );
        }

        summary
    }

    /// Converts metrics to fgbio-compatible key-value-description format.
    ///
    /// Returns a vector of `ConsensusKvMetric` rows matching fgbio's vertical
    /// metrics output format used by `CallMolecularConsensusReads`.
    #[must_use]
    pub fn to_kv_metrics(&self) -> Vec<ConsensusKvMetric> {
        let mut metrics = Vec::new();

        let raw_reads_used = self.total_input_reads.saturating_sub(self.filtered_reads);
        let frac_used = if self.total_input_reads == 0 {
            0.0
        } else {
            #[expect(clippy::cast_precision_loss, reason = "read counts never exceed 2^53")]
            let numerator = raw_reads_used as f64;
            #[expect(clippy::cast_precision_loss, reason = "read counts never exceed 2^53")]
            let denominator = self.total_input_reads as f64;
            numerator / denominator
        };

        // Core metrics matching fgbio's output
        let core = [
            ("raw_reads_considered", self.total_input_reads.to_string(), "Total raw reads considered from input file"),
            ("raw_reads_rejected", self.filtered_reads.to_string(), "Total number of raw reads rejected before consensus calling"),
            ("raw_reads_used", raw_reads_used.to_string(), "Total count of raw reads used in consensus reads"),
            ("frac_raw_reads_used", format_float(frac_used), "Fraction of raw reads used in consensus reads"),
            ("raw_reads_rejected_for_insufficient_support", self.rejected_insufficient_support.to_string(), "Insufficient reads to generate a consensus"),
            ("raw_reads_rejected_for_minority_alignment", self.rejected_minority_alignment.to_string(), "Read has a different, and minority, set of indels"),
            ("raw_reads_rejected_for_orphan_consensus", self.rejected_orphan_consensus.to_string(), "Only one of R1 or R2 consensus generated"),
            ("raw_reads_rejected_for_zero_bases_post_trimming", self.rejected_zero_bases_post_trimming.to_string(), "Read or mate had zero bases post trimming"),
        ];
        for (key, value, description) in core {
            metrics.push(ConsensusKvMetric::new(key, value, description));
        }

        // Optional rejection reasons (only included when non-zero)
        self.push_optional_rejections(&mut metrics);

        // Final output metric
        metrics.push(ConsensusKvMetric::new(
            "consensus_reads_emitted",
            self.consensus_reads.to_string(),
            "Total number of consensus reads (R1+R2=2) emitted.",
        ));

        metrics
    }

    /// Appends fgumi-specific rejection metrics with non-zero counts.
    fn push_optional_rejections(&self, metrics: &mut Vec<ConsensusKvMetric>) {
        let optional = [
            ("raw_reads_rejected_for_insufficient_strand_support", self.rejected_insufficient_strand_support, "Insufficient strand support for consensus"),
            ("raw_reads_rejected_for_low_base_quality", self.rejected_low_base_quality, "Low base quality"),
            ("raw_reads_rejected_for_excessive_n_bases", self.rejected_excessive_n_bases, "Excessive N bases in read"),
            ("raw_reads_rejected_for_no_valid_alignment", self.rejected_no_valid_alignment, "No valid alignment found"),
            ("raw_reads_rejected_for_low_mapping_quality", self.rejected_low_mapping_quality, "Low mapping quality"),
            ("raw_reads_rejected_for_n_bases_in_umi", self.rejected_n_bases_in_umi, "N bases in UMI sequence"),
            ("raw_reads_rejected_for_missing_umi", self.rejected_missing_umi, "Read lacks required UMI tag"),
            ("raw_reads_rejected_for_not_passing_filter", self.rejected_not_passing_filter, "Read did not pass vendor filter"),
            ("raw_reads_rejected_for_low_mean_quality", self.rejected_low_mean_quality, "Low mean base quality"),
            ("raw_reads_rejected_for_insufficient_min_depth", self.rejected_insufficient_min_depth, "Insufficient minimum read depth"),
            ("raw_reads_rejected_for_excessive_error_rate", self.rejected_excessive_error_rate, "Excessive error rate"),
            ("raw_reads_rejected_for_umi_too_short", self.rejected_umi_too_short, "UMI sequence too short"),
            ("raw_reads_rejected_for_single_strand_only", self.rejected_same_strand_only, "Only generating one strand of duplex consensus"),
            ("raw_reads_rejected_for_duplicate_umi", self.rejected_duplicate_umi, "Duplicate UMI detected"),
        ];
        for (key, count, description) in optional {
            if count > 0 {
                metrics.push(ConsensusKvMetric::new(key, count.to_string(), description));
            }
        }
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

        let kv_metrics = metrics.to_kv_metrics();

        // Should have core metrics
        let raw_reads = kv_metrics.iter().find(|m| m.key == "raw_reads_considered");
        assert!(raw_reads.is_some());
        assert_eq!(raw_reads.unwrap().value, "1000");
        assert_eq!(raw_reads.unwrap().description, "Total raw reads considered from input file");

        let rejected = kv_metrics.iter().find(|m| m.key == "raw_reads_rejected");
        assert!(rejected.is_some());
        assert_eq!(rejected.unwrap().value, "100");

        let used = kv_metrics.iter().find(|m| m.key == "raw_reads_used");
        assert!(used.is_some());
        assert_eq!(used.unwrap().value, "900");

        let frac = kv_metrics.iter().find(|m| m.key == "frac_raw_reads_used");
        assert!(frac.is_some());
        assert_eq!(frac.unwrap().value, "0.900000");

        let consensus = kv_metrics.iter().find(|m| m.key == "consensus_reads_emitted");
        assert!(consensus.is_some());
        assert_eq!(consensus.unwrap().value, "450");

        // Should have rejection reasons
        let insuff =
            kv_metrics.iter().find(|m| m.key == "raw_reads_rejected_for_insufficient_support");
        assert!(insuff.is_some());
        assert_eq!(insuff.unwrap().value, "50");

        let minority =
            kv_metrics.iter().find(|m| m.key == "raw_reads_rejected_for_minority_alignment");
        assert!(minority.is_some());
        assert_eq!(minority.unwrap().value, "50");
    }

    #[test]
    fn test_to_kv_metrics_zero_reads() {
        let metrics = ConsensusMetrics::new();
        let kv_metrics = metrics.to_kv_metrics();

        // frac_raw_reads_used should be 0.0 when total_input_reads is 0
        let frac = kv_metrics.iter().find(|m| m.key == "frac_raw_reads_used");
        assert!(frac.is_some());
        assert_eq!(frac.unwrap().value, "0.000000");
    }

    #[test]
    fn test_to_kv_metrics_only_nonzero_optional_rejections() {
        let mut metrics = ConsensusMetrics::new();
        metrics.total_input_reads = 100;
        metrics.rejected_low_base_quality = 5; // This is an optional rejection

        let kv_metrics = metrics.to_kv_metrics();

        // Should have low_base_quality since it's non-zero
        let lbq = kv_metrics.iter().find(|m| m.key == "raw_reads_rejected_for_low_base_quality");
        assert!(lbq.is_some());
        assert_eq!(lbq.unwrap().value, "5");

        // Should NOT have excessive_n_bases since it's zero
        let en = kv_metrics.iter().find(|m| m.key == "raw_reads_rejected_for_excessive_n_bases");
        assert!(en.is_none());
    }

    #[test]
    fn test_consensus_kv_metric_new() {
        let kv = ConsensusKvMetric::new("test_key", 42u64.to_string(), "A test metric");
        assert_eq!(kv.key, "test_key");
        assert_eq!(kv.value, "42");
        assert_eq!(kv.description, "A test metric");
    }
}
