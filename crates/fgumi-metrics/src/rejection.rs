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
}

impl fmt::Display for RejectionReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.description())
    }
}

/// Formats a count with thousands separators.
///
/// # Panics
///
/// Cannot panic: input is always valid UTF-8 since it comes from `u64::to_string()`.
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

    bytes
        .rchunks(3)
        .rev()
        .map(|chunk| std::str::from_utf8(chunk).unwrap())
        .collect::<Vec<_>>()
        .join(",")
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
    fn test_format_count() {
        assert_eq!(format_count(0), "0");
        assert_eq!(format_count(123), "123");
        assert_eq!(format_count(1234), "1,234");
        assert_eq!(format_count(1_234_567), "1,234,567");
        assert_eq!(format_count(1_000_000_000), "1,000,000,000");
    }
}
