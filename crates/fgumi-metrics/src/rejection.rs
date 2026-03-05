//! Rejection reason tracking for reads and templates.
//!
//! This module provides rejection reason types for tracking why reads or templates
//! are rejected during processing, enabling detailed metrics and debugging.

use serde::{Deserialize, Serialize};
use std::fmt;

/// Reasons why a read or template was rejected during processing.
///
/// Each variant represents a specific reason for rejection, allowing for
/// detailed tracking and reporting of why data was filtered out.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum RejectionReason {
    /// Insufficient reads to generate a consensus
    InsufficientSupport,
    /// Read has a different, and minority, set of indels
    MinorityAlignment,
    /// Too few reads agreed on the strand orientation
    InsufficientStrandSupport,
    /// Base quality scores were below threshold
    LowBaseQuality,
    /// Read group had too many N bases
    ExcessiveNBases,
    /// Template had no valid alignments
    NoValidAlignment,
    /// Reads failed mapping quality threshold
    LowMappingQuality,
    /// UMI contained N bases
    NBasesInUmi,
    /// Read lacks required UMI tag
    MissingUmi,
    /// Reads were marked as not passing filter (PF flag)
    NotPassingFilter,
    /// Consensus read had too low mean quality
    LowMeanQuality,
    /// Consensus read had insufficient minimum depth
    InsufficientMinDepth,
    /// Consensus read had excessive error rate
    ExcessiveErrorRate,
    /// UMI was too short
    UmiTooShort,
    /// Template had reads on same strand only (no proper pair)
    SameStrandOnly,
    /// Duplicate UMI at same genomic position
    DuplicateUmi,
    /// Only one of R1 or R2 consensus was generated (orphan)
    OrphanConsensus,
    /// Read had zero bases after quality trimming
    ZeroBasesPostTrimming,
}

impl RejectionReason {
    /// Returns a human-readable description.
    #[must_use]
    pub fn description(&self) -> &'static str {
        match self {
            Self::InsufficientSupport => "Insufficient reads to generate a consensus",
            Self::MinorityAlignment => "Read has a different, and minority, set of indels",
            Self::InsufficientStrandSupport => "Too few reads agreed on the strand orientation",
            Self::LowBaseQuality => "Base quality scores were below threshold",
            Self::ExcessiveNBases => "Read group had too many N bases",
            Self::NoValidAlignment => "Template had no valid alignments",
            Self::LowMappingQuality => "Reads failed mapping quality threshold",
            Self::NBasesInUmi => "UMI contained N bases",
            Self::MissingUmi => "Read lacks required UMI tag",
            Self::NotPassingFilter => "Reads were marked as not passing filter",
            Self::LowMeanQuality => "Consensus read had too low mean quality",
            Self::InsufficientMinDepth => "Consensus read had insufficient minimum depth",
            Self::ExcessiveErrorRate => "Consensus read had excessive error rate",
            Self::UmiTooShort => "UMI was too short",
            Self::SameStrandOnly => "Template had reads on same strand only",
            Self::DuplicateUmi => "Duplicate UMI at same genomic position",
            Self::OrphanConsensus => "Only one of R1 or R2 consensus generated",
            Self::ZeroBasesPostTrimming => "Read or mate had zero bases post trimming",
        }
    }

    /// Returns the TSV metric key for this rejection reason.
    #[must_use]
    pub fn tsv_key(&self) -> &'static str {
        match self {
            Self::InsufficientSupport => "raw_reads_rejected_for_insufficient_support",
            Self::MinorityAlignment => "raw_reads_rejected_for_minority_alignment",
            Self::InsufficientStrandSupport => "raw_reads_rejected_for_insufficient_strand_support",
            Self::LowBaseQuality => "raw_reads_rejected_for_low_base_quality",
            Self::ExcessiveNBases => "raw_reads_rejected_for_excessive_n_bases",
            Self::NoValidAlignment => "raw_reads_rejected_for_no_valid_alignment",
            Self::LowMappingQuality => "raw_reads_rejected_for_low_mapping_quality",
            Self::NBasesInUmi => "raw_reads_rejected_for_n_bases_in_umi",
            Self::MissingUmi => "raw_reads_rejected_for_missing_umi",
            Self::NotPassingFilter => "raw_reads_rejected_for_not_passing_filter",
            Self::LowMeanQuality => "raw_reads_rejected_for_low_mean_quality",
            Self::InsufficientMinDepth => "raw_reads_rejected_for_insufficient_min_depth",
            Self::ExcessiveErrorRate => "raw_reads_rejected_for_excessive_error_rate",
            Self::UmiTooShort => "raw_reads_rejected_for_umi_too_short",
            Self::SameStrandOnly => "raw_reads_rejected_for_single_strand_only",
            Self::DuplicateUmi => "raw_reads_rejected_for_duplicate_umi",
            Self::OrphanConsensus => "raw_reads_rejected_for_orphan_consensus",
            Self::ZeroBasesPostTrimming => "raw_reads_rejected_for_zero_bases_post_trimming",
        }
    }

    /// Returns a short description for key-value metrics output.
    #[must_use]
    pub fn kv_description(&self) -> &'static str {
        match self {
            Self::InsufficientSupport => "Insufficient reads to generate a consensus",
            Self::MinorityAlignment => "Read has a different, and minority, set of indels",
            Self::InsufficientStrandSupport => "Insufficient strand support for consensus",
            Self::LowBaseQuality => "Low base quality",
            Self::ExcessiveNBases => "Excessive N bases in read",
            Self::NoValidAlignment => "No valid alignment found",
            Self::LowMappingQuality => "Low mapping quality",
            Self::NBasesInUmi => "N bases in UMI sequence",
            Self::MissingUmi => "Read lacks required UMI tag",
            Self::NotPassingFilter => "Read did not pass vendor filter",
            Self::LowMeanQuality => "Low mean base quality",
            Self::InsufficientMinDepth => "Insufficient minimum read depth",
            Self::ExcessiveErrorRate => "Excessive error rate",
            Self::UmiTooShort => "UMI sequence too short",
            Self::SameStrandOnly => "Only generating one strand of duplex consensus",
            Self::DuplicateUmi => "Duplicate UMI detected",
            Self::OrphanConsensus => "Only one of R1 or R2 consensus generated",
            Self::ZeroBasesPostTrimming => "Read or mate had zero bases post trimming",
        }
    }
}

impl fmt::Display for RejectionReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.description())
    }
}

/// Formats a count with thousands separators.
///
/// # Examples
///
/// ```
/// use fgumi_metrics::rejection::format_count;
///
/// assert_eq!(format_count(1234567), "1,234,567");
/// assert_eq!(format_count(123), "123");
/// ```
#[must_use]
pub fn format_count(n: u64) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let len = bytes.len();
    let num_commas = if len > 3 { (len - 1) / 3 } else { 0 };
    let mut result = String::with_capacity(len + num_commas);
    for (i, &byte) in bytes.iter().enumerate() {
        if i > 0 && (len - i).is_multiple_of(3) {
            result.push(',');
        }
        result.push(byte as char);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rejection_reason_description() {
        assert!(RejectionReason::LowBaseQuality.description().contains("quality"));
        assert!(RejectionReason::InsufficientSupport.description().contains("Insufficient"));
        assert_eq!(
            RejectionReason::MinorityAlignment.to_string(),
            "Read has a different, and minority, set of indels"
        );
    }

    #[test]
    fn test_tsv_key_prefix() {
        let all_reasons = [
            RejectionReason::InsufficientSupport,
            RejectionReason::MinorityAlignment,
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
            RejectionReason::SameStrandOnly,
            RejectionReason::DuplicateUmi,
            RejectionReason::OrphanConsensus,
            RejectionReason::ZeroBasesPostTrimming,
        ];
        for reason in &all_reasons {
            assert!(
                reason.tsv_key().starts_with("raw_reads_rejected_for_"),
                "tsv_key for {:?} does not have expected prefix: {}",
                reason,
                reason.tsv_key()
            );
        }
    }

    #[test]
    fn test_kv_description_non_empty() {
        let all_reasons = [
            RejectionReason::InsufficientSupport,
            RejectionReason::MinorityAlignment,
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
            RejectionReason::SameStrandOnly,
            RejectionReason::DuplicateUmi,
            RejectionReason::OrphanConsensus,
            RejectionReason::ZeroBasesPostTrimming,
        ];
        for reason in &all_reasons {
            assert!(
                !reason.kv_description().is_empty(),
                "kv_description for {:?} is empty",
                reason
            );
        }
    }

    #[test]
    fn test_format_count() {
        assert_eq!(format_count(0), "0");
        assert_eq!(format_count(123), "123");
        assert_eq!(format_count(1234), "1,234");
        assert_eq!(format_count(1_234_567), "1,234,567");
        assert_eq!(format_count(1_000_000_000), "1,000,000,000");
    }
}
