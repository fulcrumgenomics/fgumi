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
// `EnumIter` (test-only) auto-generates `RejectionReason::iter()` so coverage assertions
// enumerate every variant without a hand-maintained list that can silently drift.
#[cfg_attr(test, derive(strum::EnumIter))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum RejectionReason {
    /// Insufficient reads to generate a consensus
    InsufficientSupport,
    /// Reads has a different, and minority, set of indels
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
    /// Unpaired/fragment reads supplied to the duplex/codec caller
    NonPairedReads,
    /// Potential collision between independent duplex molecules sharing an MI
    DuplicateUmi,
    /// Only one of R1 or R2 consensus was generated (orphan)
    OrphanConsensus,
    /// Read had zero bases after quality trimming
    ZeroBasesPostTrimming,
    /// Template did not have a single primary FR pair of reads (codec)
    NotPrimaryFrPair,
    /// Overlap between R1s and R2s too short for CODEC calling (fgbio
    /// `r1_r2_overlap_too_short`; used by codec)
    R1R2OverlapTooShort,
    /// Indel error between the top/bottom strands (fgbio
    /// `indel_error_between_strands`; used by codec)
    IndelErrorBetweenStrands,
    /// Too many errors between the top/bottom strands (fgbio
    /// `high_duplex_disagreement`; used by codec)
    HighDuplexDisagreement,
    /// Overlap clipping failed (fgbio `clip_overlap_failed`; used by codec)
    ClipOverlapFailed,
}

impl RejectionReason {
    /// Returns a human-readable description.
    #[must_use]
    pub fn description(&self) -> &'static str {
        match self {
            Self::InsufficientSupport => "Insufficient reads to generate a consensus",
            // fgbio spells this "Reads has …" (sic) — matched verbatim for metric parity.
            Self::MinorityAlignment => "Reads has a different, and minority, set of indels",
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
            Self::NonPairedReads => "Unpaired/fragment reads not supported by Duplex caller",
            Self::DuplicateUmi => "Potential collision between independent duplex molecules",
            Self::OrphanConsensus => "Only one of R1 or R2 consensus generated",
            Self::ZeroBasesPostTrimming => "Read or mate had zero bases post trimming",
            Self::NotPrimaryFrPair => "Template did not have a single primary FR pair of reads",
            Self::R1R2OverlapTooShort => "Overlap between R1s and R2s too short for CODEC calling",
            Self::IndelErrorBetweenStrands => "Indel error between top/bottom strands",
            // fgbio's exact string (note the "top/bottoms" typo), matched for parity.
            Self::HighDuplexDisagreement => "Too many errors between top/bottoms strands",
            Self::ClipOverlapFailed => "See https://github.com/fulcrumgenomics/fgbio/issues/1090",
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
            Self::NonPairedReads => "raw_reads_rejected_for_non_paired_reads",
            Self::DuplicateUmi => "raw_reads_rejected_for_potential_umi_collision",
            Self::OrphanConsensus => "raw_reads_rejected_for_orphan_consensus",
            Self::ZeroBasesPostTrimming => "raw_reads_rejected_for_zero_bases_post_trimming",
            Self::NotPrimaryFrPair => "raw_reads_rejected_for_not_primary_fr_pair",
            Self::R1R2OverlapTooShort => "raw_reads_rejected_for_r1_r2_overlap_too_short",
            Self::IndelErrorBetweenStrands => "raw_reads_rejected_for_indel_error_between_strands",
            Self::HighDuplexDisagreement => "raw_reads_rejected_for_high_duplex_disagreement",
            Self::ClipOverlapFailed => "raw_reads_rejected_for_clip_overlap_failed",
        }
    }

    /// Returns a short description for key-value metrics output.
    #[must_use]
    pub fn kv_description(&self) -> &'static str {
        match self {
            Self::InsufficientSupport => "Insufficient reads to generate a consensus",
            Self::MinorityAlignment => "Reads has a different, and minority, set of indels",
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
            // fgbio's exact title-cased string for `single_strand_only`, matched for parity.
            Self::SameStrandOnly => "Only Generating One Strand of Duplex Consensus",
            Self::NonPairedReads => "Unpaired/fragment reads not supported by Duplex caller",
            Self::DuplicateUmi => "Potential collision between independent duplex molecules",
            Self::OrphanConsensus => "Only one of R1 or R2 consensus generated",
            Self::ZeroBasesPostTrimming => "Read or mate had zero bases post trimming",
            Self::NotPrimaryFrPair => "Template did not have a single primary FR pair of reads",
            Self::R1R2OverlapTooShort => "Overlap between R1s and R2s too short for CODEC calling",
            Self::IndelErrorBetweenStrands => "Indel error between top/bottom strands",
            Self::HighDuplexDisagreement => "Too many errors between top/bottoms strands",
            Self::ClipOverlapFailed => "See https://github.com/fulcrumgenomics/fgbio/issues/1090",
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
        // Matches fgbio's exact (grammatically-odd) wording for metric parity.
        assert_eq!(
            RejectionReason::MinorityAlignment.to_string(),
            "Reads has a different, and minority, set of indels"
        );
    }

    use strum::IntoEnumIterator;

    /// The four codec-only reasons carry fgbio-verbatim strings (`UmiConsensusCaller.scala:84-87`),
    /// including the deliberate "top/bottoms" typo and the issue-URL description. Pin the exact
    /// `description()`, `kv_description()`, and `tsv_key()` so a well-meaning edit can't silently
    /// break fgbio parity.
    #[rstest::rstest]
    #[case::r1_r2_overlap_too_short(
        RejectionReason::R1R2OverlapTooShort,
        "Overlap between R1s and R2s too short for CODEC calling",
        "raw_reads_rejected_for_r1_r2_overlap_too_short"
    )]
    #[case::indel_error_between_strands(
        RejectionReason::IndelErrorBetweenStrands,
        "Indel error between top/bottom strands",
        "raw_reads_rejected_for_indel_error_between_strands"
    )]
    #[case::high_duplex_disagreement(
        RejectionReason::HighDuplexDisagreement,
        "Too many errors between top/bottoms strands",
        "raw_reads_rejected_for_high_duplex_disagreement"
    )]
    #[case::clip_overlap_failed(
        RejectionReason::ClipOverlapFailed,
        "See https://github.com/fulcrumgenomics/fgbio/issues/1090",
        "raw_reads_rejected_for_clip_overlap_failed"
    )]
    fn test_codec_reason_strings_match_fgbio(
        #[case] reason: RejectionReason,
        #[case] expected_description: &str,
        #[case] expected_tsv_key: &str,
    ) {
        // fgbio uses one description string per reason for both the human-readable and KV forms.
        assert_eq!(reason.description(), expected_description);
        assert_eq!(reason.kv_description(), expected_description);
        assert_eq!(reason.tsv_key(), expected_tsv_key);
    }

    #[test]
    fn test_tsv_key_prefix() {
        // `RejectionReason::iter()` (via `EnumIter`) covers every variant automatically, so a
        // newly added reason is exercised here without updating a hand-maintained list.
        for reason in RejectionReason::iter() {
            assert!(
                reason.tsv_key().starts_with("raw_reads_rejected_for_"),
                "tsv_key for {:?} does not have expected prefix: {}",
                reason,
                reason.tsv_key()
            );
        }
    }

    #[test]
    fn test_duplex_row_keys_match_fgbio() {
        // R2-MET-05 / R2-UCC-01: the duplex/codec rejection rows carry fgbio's exact keys.
        assert_eq!(
            RejectionReason::NonPairedReads.tsv_key(),
            "raw_reads_rejected_for_non_paired_reads"
        );
        assert_eq!(
            RejectionReason::SameStrandOnly.tsv_key(),
            "raw_reads_rejected_for_single_strand_only"
        );
        assert_eq!(
            RejectionReason::DuplicateUmi.tsv_key(),
            "raw_reads_rejected_for_potential_umi_collision"
        );
    }

    #[test]
    fn test_kv_description_non_empty() {
        for reason in RejectionReason::iter() {
            assert!(!reason.kv_description().is_empty(), "kv_description for {reason:?} is empty");
        }
    }

    /// Pin the exact fgbio-parity `tsv_key` + `kv_description` for the two reasons whose strings
    /// are deliberately matched to fgbio (see the `for parity` comments above). A non-empty check
    /// would let these silently drift; exact equality is what keeps them aligned with fgbio.
    #[test]
    fn test_kv_description_matches_fgbio_verbatim() {
        assert_eq!(
            RejectionReason::MinorityAlignment.tsv_key(),
            "raw_reads_rejected_for_minority_alignment"
        );
        assert_eq!(
            RejectionReason::MinorityAlignment.kv_description(),
            "Reads has a different, and minority, set of indels"
        );

        assert_eq!(
            RejectionReason::SameStrandOnly.tsv_key(),
            "raw_reads_rejected_for_single_strand_only"
        );
        assert_eq!(
            RejectionReason::SameStrandOnly.kv_description(),
            "Only Generating One Strand of Duplex Consensus"
        );
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
