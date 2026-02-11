//! # Consensus Calling Infrastructure
//!
//! This module provides the core trait and supporting structures for consensus calling algorithms
//! that generate high-quality consensus reads from groups of raw reads sharing the same UMI.
//!
//! ## Overview
//!
//! Consensus calling is the process of combining multiple raw sequencing reads from the same
//! source molecule (identified by a shared UMI) into a single, higher-quality consensus read.
//! By comparing bases across all reads at each position, consensus calling can:
//!
//! - Correct random sequencing errors
//! - Assign quality scores that reflect the confidence in each consensus base
//! - Filter out reads that don't meet quality thresholds
//! - Track detailed statistics about the consensus process
//!
//! ## Consensus Calling Process
//!
//! The typical workflow for consensus calling:
//!
//! 1. **Grouping**: Reads are grouped by UMI (molecule ID), typically after UMI assignment
//! 2. **Filtering**: Reads within each group are filtered based on quality criteria:
//!    - Remove unmapped reads (if requiring mapped reads)
//!    - Remove secondary/supplementary alignments
//!    - Remove reads with too many N bases
//!    - Remove reads failing QC flags
//!    - Filter minority alignments (reads mapping to different locations than consensus)
//! 3. **Consensus Calling**: For each position in the read:
//!    - Collect all observed bases and their quality scores
//!    - Use a likelihood-based model to determine the most likely true base
//!    - Calculate a consensus quality score reflecting confidence
//! 4. **Output**: Generate consensus read(s) with:
//!    - Consensus sequence and quality scores
//!    - Metadata tags describing depth, error rates, and per-base statistics
//!
//! ## Likelihood-Based Quality Model
//!
//! Consensus calling uses Bayesian inference to call bases and assign quality scores:
//!
//! **For each candidate base** (A, C, G, T), we calculate the likelihood of observing all
//! the data given that the true base is that candidate:
//!
//! ```text
//! L(B|data) = P(obs_1|B) * P(obs_2|B) * ... * P(obs_n|B)
//! ```
//!
//! Where for each observation:
//! - If observed base matches candidate: `P(obs|B) = 1 - error_rate`
//! - If observed base differs: `P(obs|B) = error_rate / 3`
//!
//! **Quality score calculation**: After normalizing likelihoods to get posterior probabilities,
//! the consensus quality is derived from the Phred-scaled error probability:
//!
//! ```text
//! Q = -10 * log10(P_error)
//! ```
//!
//! Where `P_error = 1 - P(most_likely_base)`.
//!
//! **Error rate priors**: The model incorporates two error rates:
//! - **Pre-UMI error rate**: The base error rate in the original molecule (before amplification)
//! - **Post-UMI error rate**: Additional errors introduced during sequencing
//!
//! These are combined to model the total error rate for each observation.
//!
//! ## Rejection Reasons
//!
//! Reads can be filtered out for various reasons during consensus calling. The `RejectionReason`
//! enum documents why each read was excluded:
//!
//! - **`InsufficientReads`**: Group has fewer than minimum required reads for consensus
//! - **`QualityTooLow`**: Read has average quality below threshold
//! - **`FragmentRead`**: Unpaired read when paired reads are required
//! - **Unmapped**: Read is unmapped when mapped reads are required
//! - **`SecondaryOrSupplementary`**: Read is not a primary alignment
//! - **`TooManyNs`**: Read contains excessive N bases
//! - **`MinorityAlignment`**: Read maps to different location than consensus alignment
//! - **`FailedQC`**: Read has failed QC flag set
//! - **`MissingUmi`**: Read lacks required UMI tag
//! - **`QualityTrimmed`**: Read was trimmed below acceptable length
//!
//! These rejection reasons enable detailed QC reporting and help diagnose library preparation
//! or sequencing issues.
//!
//! ## Implementing Custom Consensus Callers
//!
//! To implement a custom consensus caller, implement the `ConsensusCaller` trait:
//!
//! ```rust,ignore
//! use fgumi_lib::consensus::caller::{ConsensusCaller, ConsensusCallingStats, ConsensusOutput};
//! use anyhow::Result;
//!
//! pub struct MyConsensusCaller {
//!     stats: ConsensusCallingStats,
//!     // ... other fields
//! }
//!
//! impl ConsensusCaller for MyConsensusCaller {
//!     fn consensus_reads(&mut self, records: Vec<Vec<u8>>) -> Result<ConsensusOutput> {
//!         // Each Vec<u8> is a raw BAM record (without block_size prefix)
//!         self.stats.record_input(records.len());
//!
//!         // ... filter reads, call consensus, write into output ...
//!
//!         self.stats.record_consensus();
//!         Ok(output)
//!     }
//!
//!     fn total_reads(&self) -> usize { self.stats.total_reads }
//!     fn total_filtered(&self) -> usize { self.stats.filtered_reads }
//!     fn consensus_reads_constructed(&self) -> usize { self.stats.consensus_reads }
//!     fn statistics(&self) -> ConsensusCallingStats { self.stats.clone() }
//!     fn log_statistics(&self) { /* ... */ }
//! }
//! ```
//!
//! ## Usage Patterns
//!
//! ### Basic Consensus Calling
//!
//! ```rust,ignore
//! use fgumi_lib::consensus::vanilla_consensus_caller::VanillaUmiConsensusCaller;
//! use fgumi_lib::consensus::caller::ConsensusCaller;
//!
//! // Create consensus caller
//! let mut caller = VanillaUmiConsensusCaller::new(
//!     "consensus".to_string(),  // read name prefix
//!     "RG1".to_string(),         // read group ID
//!     options,                   // configuration options
//! );
//!
//! // Process a group of raw-byte BAM records with the same UMI
//! let records: Vec<Vec<u8>> = vec![raw_record1, raw_record2, raw_record3];
//! let consensus = caller.consensus_reads(records)?;
//!
//! // Check statistics
//! println!("Consensus reads: {}", caller.consensus_reads_constructed());
//! println!("Filtered reads: {}", caller.total_filtered());
//! caller.log_statistics();
//! ```
//!
//! ### Tracking Rejection Reasons
//!
//! ```rust,ignore
//! let stats = caller.statistics();
//! println!("Total reads processed: {}", stats.total_reads);
//! println!("Consensus reads generated: {}", stats.consensus_reads);
//!
//! // Break down rejections by reason
//! for (reason, count) in &stats.rejection_reasons {
//!     println!("  {}: {}", reason.description(), count);
//! }
//! ```
//!
//! ## See Also
//!
//! - `vanilla_consensus_caller`: Standard single-strand consensus implementation
//! - `duplex_caller`: Two-stage consensus for duplex sequencing
//! - `base_builder`: Core likelihood-based consensus base calling logic

use anyhow::Result;
use noodles::sam::Header;
use std::collections::HashMap;

/// Pre-serialized consensus output.
///
/// Contains concatenated BAM records (each prefixed with 4-byte `block_size`)
/// ready for direct BGZF compression, bypassing noodles encoding.
#[derive(Default)]
pub struct ConsensusOutput {
    /// Concatenated BAM records with `block_size` prefixes.
    pub data: Vec<u8>,
    /// Number of consensus records in the buffer.
    pub count: usize,
}

impl ConsensusOutput {
    /// Creates a new empty `ConsensusOutput`.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Appends another `ConsensusOutput`'s data to this one (copies).
    pub fn extend(&mut self, other: &ConsensusOutput) {
        self.data.extend_from_slice(&other.data);
        self.count += other.count;
    }

    /// Merges another `ConsensusOutput` by moving its data (avoids copy).
    pub fn merge(&mut self, mut other: ConsensusOutput) {
        self.data.append(&mut other.data);
        self.count += other.count;
    }
}

/// The main trait for consensus callers that generate consensus reads from groups of raw reads.
///
/// All consensus callers operate on raw BAM byte records (`Vec<u8>`) to avoid the overhead
/// of noodles `RecordBuf` parsing. Each `Vec<u8>` is a complete BAM record without the
/// 4-byte `block_size` prefix.
pub trait ConsensusCaller: Send + Sync {
    /// Takes a group of raw-byte BAM records with the same UMI and generates consensus reads.
    ///
    /// # Arguments
    /// * `records` - Raw BAM byte records from the same source molecule (same MI tag)
    ///
    /// # Returns
    /// Pre-serialized consensus output (may be empty if the group doesn't meet minimum requirements)
    fn consensus_reads(&mut self, records: Vec<Vec<u8>>) -> Result<ConsensusOutput>;

    /// Returns the total number of input reads examined by the consensus caller
    fn total_reads(&self) -> usize;

    /// Returns the total number of reads filtered/rejected for any reason
    fn total_filtered(&self) -> usize;

    /// Returns the number of consensus reads constructed by this caller
    fn consensus_reads_constructed(&self) -> usize;

    /// Returns statistics about the consensus calling process
    fn statistics(&self) -> ConsensusCallingStats;

    /// Logs statistics about consensus calling to the logger
    fn log_statistics(&self);
}

/// Statistics tracked during consensus calling
#[derive(Debug, Clone)]
pub struct ConsensusCallingStats {
    /// Total number of input reads processed
    pub total_reads: usize,
    /// Number of consensus reads generated
    pub consensus_reads: usize,
    /// Number of reads filtered/rejected
    pub filtered_reads: usize,
    /// Breakdown of rejection reasons
    pub rejection_reasons: HashMap<RejectionReason, usize>,
}

impl ConsensusCallingStats {
    /// Creates a new empty statistics object
    #[must_use]
    pub fn new() -> Self {
        Self {
            total_reads: 0,
            consensus_reads: 0,
            filtered_reads: 0,
            rejection_reasons: HashMap::new(),
        }
    }

    /// Records that reads were rejected for a given reason
    pub fn record_rejection(&mut self, reason: RejectionReason, count: usize) {
        self.filtered_reads += count;
        *self.rejection_reasons.entry(reason).or_insert(0) += count;
    }

    /// Records that a consensus read was created
    pub fn record_consensus(&mut self) {
        self.consensus_reads += 1;
    }

    /// Records that input reads were processed
    pub fn record_input(&mut self, count: usize) {
        self.total_reads += count;
    }

    /// Merges statistics from another instance into this one
    pub fn merge(&mut self, other: &ConsensusCallingStats) {
        self.total_reads += other.total_reads;
        self.consensus_reads += other.consensus_reads;
        self.filtered_reads += other.filtered_reads;
        for (reason, count) in &other.rejection_reasons {
            *self.rejection_reasons.entry(*reason).or_insert(0) += count;
        }
    }
}

impl Default for ConsensusCallingStats {
    fn default() -> Self {
        Self::new()
    }
}

/// Logs standard consensus calling statistics.
///
/// This is a helper function to reduce duplication across consensus caller implementations.
pub fn log_consensus_statistics(caller_name: &str, stats: &ConsensusCallingStats) {
    log::info!("{caller_name} Consensus Calling Statistics:");
    log::info!("  Total input reads: {}", stats.total_reads);
    log::info!("  Consensus reads generated: {}", stats.consensus_reads);
    log::info!("  Reads filtered: {}", stats.filtered_reads);

    if !stats.rejection_reasons.is_empty() {
        log::info!("  Rejection reasons:");
        for (reason, count) in &stats.rejection_reasons {
            log::info!("    {:?}: {}", reason, count);
        }
    }
}

/// Common options shared across consensus callers.
///
/// These fields are duplicated across `VanillaUmiConsensusOptions`, `CodecConsensusOptions`,
/// and `DuplexConsensusOptions`. This struct provides a common base to reduce duplication.
#[derive(Debug, Clone)]
pub struct ConsensusOptionsBase {
    /// Pre-UMI error rate (Phred scale)
    pub error_rate_pre_umi: crate::phred::PhredScore,

    /// Post-UMI error rate (Phred scale)
    pub error_rate_post_umi: crate::phred::PhredScore,

    /// Minimum base quality to include in consensus
    pub min_input_base_quality: crate::phred::PhredScore,

    /// Whether to produce per-base tags (cd, ce, etc.)
    pub produce_per_base_tags: bool,

    /// Whether to quality-trim reads before consensus calling
    pub trim: bool,

    /// Minimum consensus base quality (output bases below this are masked to N)
    pub min_consensus_base_quality: crate::phred::PhredScore,
}

impl Default for ConsensusOptionsBase {
    fn default() -> Self {
        Self {
            error_rate_pre_umi: 45,
            error_rate_post_umi: 40,
            min_input_base_quality: 10,
            produce_per_base_tags: true,
            trim: false,
            min_consensus_base_quality: 40,
        }
    }
}

/// Calculates error rate from depths and errors arrays.
///
/// This is the standard formula used across all consensus callers:
/// `total_errors / total_depth`
///
/// Returns 0.0 if total depth is zero.
#[must_use]
pub fn calculate_error_rate(depths: &[u16], errors: &[u16]) -> f32 {
    let total_depth: u32 = depths.iter().map(|&d| u32::from(d)).sum();
    if total_depth == 0 {
        return 0.0;
    }
    let total_errors: u32 = errors.iter().map(|&e| u32::from(e)).sum();
    total_errors as f32 / total_depth as f32
}

/// Tracker for rejected reads during consensus calling.
///
/// This struct encapsulates the common pattern of tracking rejected reads
/// that's duplicated across `VanillaUmiConsensusCaller` and `CodecConsensusCaller`.
#[derive(Debug, Default)]
pub struct RejectionTracker {
    /// Whether to track rejected reads
    pub track_rejects: bool,

    /// Rejected reads as raw bytes (only populated if `track_rejects` is true)
    rejected_records: Vec<Vec<u8>>,
}

impl RejectionTracker {
    /// Creates a new rejection tracker.
    #[must_use]
    pub fn new(track_rejects: bool) -> Self {
        Self { track_rejects, rejected_records: Vec::new() }
    }

    /// Adds rejected raw-byte records to the tracker (if tracking is enabled).
    pub fn add_rejected(&mut self, records: impl IntoIterator<Item = Vec<u8>>) {
        if self.track_rejects {
            self.rejected_records.extend(records);
        }
    }

    /// Returns a reference to the rejected records.
    #[must_use]
    pub fn rejected_records(&self) -> &[Vec<u8>] {
        &self.rejected_records
    }

    /// Takes ownership of the rejected records, leaving an empty vector.
    pub fn take_rejected_records(&mut self) -> Vec<Vec<u8>> {
        std::mem::take(&mut self.rejected_records)
    }

    /// Clears the rejected records without returning them.
    pub fn clear(&mut self) {
        self.rejected_records.clear();
    }

    /// Returns whether rejection tracking is enabled.
    #[must_use]
    pub fn is_tracking(&self) -> bool {
        self.track_rejects
    }
}

/// Reasons why reads might be rejected and not used in consensus calling
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RejectionReason {
    /// Fragment read (unpaired) when paired reads are required
    FragmentRead,
    /// Fewer than minimum required reads
    InsufficientReads,
    /// Read failed quality filters
    QualityTooLow,
    /// Read is unmapped when mapped reads are required
    Unmapped,
    /// Read is mapped when unmapped reads are required
    Mapped,
    /// Read has too many Ns in the sequence
    TooManyNs,
    /// Reads do not share the same alignment (minority alignment)
    MinorityAlignment,
    /// Read is a secondary or supplementary alignment
    SecondaryOrSupplementary,
    /// Read failed QC
    FailedQC,
    /// Read has no UMI tag
    MissingUmi,
    /// Read was quality-trimmed below acceptable length
    QualityTrimmed,
    /// Read has zero length after masking low quality bases (all bases became N)
    ZeroLengthAfterTrimming,
    /// Insufficient overlap between reads (e.g., R1/R2 overlap for duplex)
    InsufficientOverlap,
    /// Only one of R1/R2 consensus was generated (orphan)
    OrphanConsensus,
    /// Overlap boundary lands in an indel between strands
    IndelErrorBetweenStrands,
    /// Potential collision
    PotentialCollision,
    /// Other unspecified reason
    Other,
}

impl RejectionReason {
    /// Returns a single character code for the rejection reason (for tagging rejected reads)
    #[must_use]
    pub fn code(&self) -> char {
        match self {
            Self::FragmentRead => 'F',
            Self::InsufficientReads => 'I',
            Self::QualityTooLow => 'Q',
            Self::Unmapped => 'U',
            Self::Mapped => 'M',
            Self::TooManyNs => 'N',
            Self::MinorityAlignment => 'A',
            Self::SecondaryOrSupplementary => 'S',
            Self::FailedQC => 'f',
            Self::MissingUmi => 'u',
            Self::QualityTrimmed => 'T',
            Self::ZeroLengthAfterTrimming => 'Z',
            Self::InsufficientOverlap => 'V',
            Self::OrphanConsensus => 'P',
            Self::IndelErrorBetweenStrands => 'D',
            Self::PotentialCollision => 'C',
            Self::Other => 'O',
        }
    }

    /// Returns a human-readable description of the rejection reason
    #[must_use]
    pub fn description(&self) -> &'static str {
        match self {
            Self::FragmentRead => "Fragment read when paired required",
            Self::InsufficientReads => "Insufficient reads for consensus",
            Self::QualityTooLow => "Quality too low",
            Self::Unmapped => "Read is unmapped",
            Self::Mapped => "Read is mapped",
            Self::TooManyNs => "Too many Ns in sequence",
            Self::MinorityAlignment => "Minority alignment",
            Self::SecondaryOrSupplementary => "Secondary or supplementary alignment",
            Self::FailedQC => "Failed QC",
            Self::QualityTrimmed => "Quality-trimmed below acceptable length",
            Self::ZeroLengthAfterTrimming => "Zero length after masking low quality bases",
            Self::MissingUmi => "Missing UMI tag",
            Self::InsufficientOverlap => "Insufficient overlap between reads",
            Self::OrphanConsensus => "Only one of R1 or R2 consensus generated",
            Self::IndelErrorBetweenStrands => "Overlap boundary lands in indel between strands",
            Self::PotentialCollision => {
                "Potential collisions (reads with the same MI but different strands)"
            }
            Self::Other => "Other reason",
        }
    }

    /// Converts this caller-specific rejection reason to the centralized rejection tracking reason.
    ///
    /// This enables unified rejection tracking across all modules while maintaining backward
    /// compatibility with existing code.
    #[must_use]
    pub fn to_centralized(&self) -> crate::rejection::RejectionReason {
        use crate::rejection::RejectionReason as CentralReason;
        match self {
            Self::InsufficientReads => CentralReason::InsufficientSupport,
            Self::TooManyNs => CentralReason::ExcessiveNBases,
            Self::MinorityAlignment => CentralReason::MinorityAlignment,
            Self::QualityTooLow => CentralReason::LowBaseQuality,
            Self::FailedQC => CentralReason::NotPassingFilter,
            Self::MissingUmi => CentralReason::NBasesInUmi, // Closest match
            // For reasons that don't have a direct centralized equivalent,
            // we use the closest semantic match
            Self::FragmentRead => CentralReason::SameStrandOnly,
            Self::Unmapped | Self::Mapped => CentralReason::NoValidAlignment,
            Self::SecondaryOrSupplementary => CentralReason::NoValidAlignment,
            Self::QualityTrimmed => CentralReason::LowBaseQuality,
            Self::ZeroLengthAfterTrimming => CentralReason::ZeroBasesPostTrimming,
            Self::InsufficientOverlap => CentralReason::InsufficientSupport, // Close match
            Self::OrphanConsensus => CentralReason::OrphanConsensus,
            Self::IndelErrorBetweenStrands => CentralReason::NoValidAlignment, // Indel boundary error
            Self::PotentialCollision => CentralReason::DuplicateUmi,
            Self::Other => CentralReason::InsufficientSupport,
        }
    }
}

/// Creates a read name prefix from the SAM header read groups.
///
/// This matches fgbio's `makePrefixFromSamHeader` behavior:
/// - Extracts library names from all read groups (or read group ID if no library)
/// - Sorts and deduplicates the library names
/// - Joins them with `|` if total length <= 200 characters
/// - Otherwise uses a hash of the joined string
#[must_use]
pub fn make_prefix_from_header(header: &Header) -> String {
    use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    let mut ids: Vec<String> = header
        .read_groups()
        .iter()
        .map(|(rg_id, rg_map)| {
            // Try to get library name, fall back to read group ID
            rg_map
                .other_fields()
                .get(&rg_tag::LIBRARY)
                .map_or_else(|| rg_id.to_string(), std::string::ToString::to_string)
        })
        .collect();

    ids.sort();
    ids.dedup();

    // Calculate total length including separators
    let total_len: usize = ids.iter().map(|s| s.len() + 1).sum();

    if total_len <= 200 {
        ids.join("|")
    } else {
        // Use hash if too long (similar to fgbio's Murmur3 hash)
        let joined = ids.join("|");
        let mut hasher = DefaultHasher::new();
        joined.hash(&mut hasher);
        format!("{:x}", hasher.finish())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rejection_reason_codes() {
        assert_eq!(RejectionReason::FragmentRead.code(), 'F');
        assert_eq!(RejectionReason::InsufficientReads.code(), 'I');
        assert_eq!(RejectionReason::QualityTooLow.code(), 'Q');
    }

    #[test]
    fn test_rejection_reason_to_centralized() {
        use crate::rejection::RejectionReason as CentralReason;

        assert_eq!(
            RejectionReason::InsufficientReads.to_centralized(),
            CentralReason::InsufficientSupport
        );
        assert_eq!(RejectionReason::TooManyNs.to_centralized(), CentralReason::ExcessiveNBases);
        assert_eq!(
            RejectionReason::MinorityAlignment.to_centralized(),
            CentralReason::MinorityAlignment
        );
        assert_eq!(RejectionReason::QualityTooLow.to_centralized(), CentralReason::LowBaseQuality);
        assert_eq!(RejectionReason::FailedQC.to_centralized(), CentralReason::NotPassingFilter);
    }

    #[test]
    fn test_stats_tracking() {
        let mut stats = ConsensusCallingStats::new();
        assert_eq!(stats.total_reads, 0);
        assert_eq!(stats.consensus_reads, 0);
        assert_eq!(stats.filtered_reads, 0);

        stats.record_input(10);
        assert_eq!(stats.total_reads, 10);

        stats.record_rejection(RejectionReason::InsufficientReads, 3);
        assert_eq!(stats.filtered_reads, 3);
        assert_eq!(stats.rejection_reasons[&RejectionReason::InsufficientReads], 3);

        stats.record_consensus();
        assert_eq!(stats.consensus_reads, 1);
    }

    #[test]
    fn test_all_rejection_reason_codes() {
        assert_eq!(RejectionReason::FragmentRead.code(), 'F');
        assert_eq!(RejectionReason::InsufficientReads.code(), 'I');
        assert_eq!(RejectionReason::QualityTooLow.code(), 'Q');
        assert_eq!(RejectionReason::Unmapped.code(), 'U');
        assert_eq!(RejectionReason::Mapped.code(), 'M');
        assert_eq!(RejectionReason::TooManyNs.code(), 'N');
        assert_eq!(RejectionReason::MinorityAlignment.code(), 'A');
        assert_eq!(RejectionReason::SecondaryOrSupplementary.code(), 'S');
        assert_eq!(RejectionReason::FailedQC.code(), 'f');
        assert_eq!(RejectionReason::MissingUmi.code(), 'u');
        assert_eq!(RejectionReason::QualityTrimmed.code(), 'T');
        assert_eq!(RejectionReason::Other.code(), 'O');
    }

    #[test]
    fn test_all_rejection_reason_descriptions() {
        assert_eq!(
            RejectionReason::FragmentRead.description(),
            "Fragment read when paired required"
        );
        assert_eq!(
            RejectionReason::InsufficientReads.description(),
            "Insufficient reads for consensus"
        );
        assert_eq!(RejectionReason::QualityTooLow.description(), "Quality too low");
        assert_eq!(RejectionReason::Unmapped.description(), "Read is unmapped");
        assert_eq!(RejectionReason::Mapped.description(), "Read is mapped");
        assert_eq!(RejectionReason::TooManyNs.description(), "Too many Ns in sequence");
        assert_eq!(RejectionReason::MinorityAlignment.description(), "Minority alignment");
        assert_eq!(
            RejectionReason::SecondaryOrSupplementary.description(),
            "Secondary or supplementary alignment"
        );
        assert_eq!(RejectionReason::FailedQC.description(), "Failed QC");
        assert_eq!(
            RejectionReason::QualityTrimmed.description(),
            "Quality-trimmed below acceptable length"
        );
        assert_eq!(RejectionReason::MissingUmi.description(), "Missing UMI tag");
        assert_eq!(RejectionReason::Other.description(), "Other reason");
    }

    #[test]
    fn test_stats_default() {
        let stats = ConsensusCallingStats::default();
        assert_eq!(stats.total_reads, 0);
        assert_eq!(stats.consensus_reads, 0);
        assert_eq!(stats.filtered_reads, 0);
        assert!(stats.rejection_reasons.is_empty());
    }

    #[test]
    fn test_stats_merge() {
        let mut stats1 = ConsensusCallingStats::new();
        stats1.record_input(10);
        stats1.record_consensus();
        stats1.record_rejection(RejectionReason::QualityTooLow, 3);
        stats1.record_rejection(RejectionReason::InsufficientReads, 2);

        let mut stats2 = ConsensusCallingStats::new();
        stats2.record_input(20);
        stats2.record_consensus();
        stats2.record_consensus();
        stats2.record_rejection(RejectionReason::QualityTooLow, 5);
        stats2.record_rejection(RejectionReason::TooManyNs, 4);

        stats1.merge(&stats2);

        assert_eq!(stats1.total_reads, 30);
        assert_eq!(stats1.consensus_reads, 3);
        assert_eq!(stats1.filtered_reads, 14);
        assert_eq!(stats1.rejection_reasons[&RejectionReason::QualityTooLow], 8);
        assert_eq!(stats1.rejection_reasons[&RejectionReason::InsufficientReads], 2);
        assert_eq!(stats1.rejection_reasons[&RejectionReason::TooManyNs], 4);
    }

    #[test]
    fn test_stats_merge_empty() {
        let mut stats1 = ConsensusCallingStats::new();
        stats1.record_input(10);
        stats1.record_consensus();

        let stats2 = ConsensusCallingStats::new();
        stats1.merge(&stats2);

        assert_eq!(stats1.total_reads, 10);
        assert_eq!(stats1.consensus_reads, 1);
        assert_eq!(stats1.filtered_reads, 0);
    }

    #[test]
    fn test_multiple_rejections_same_reason() {
        let mut stats = ConsensusCallingStats::new();
        stats.record_rejection(RejectionReason::QualityTooLow, 5);
        stats.record_rejection(RejectionReason::QualityTooLow, 3);
        stats.record_rejection(RejectionReason::QualityTooLow, 2);

        assert_eq!(stats.filtered_reads, 10);
        assert_eq!(stats.rejection_reasons[&RejectionReason::QualityTooLow], 10);
    }

    #[test]
    fn test_multiple_consensus_records() {
        let mut stats = ConsensusCallingStats::new();
        for _ in 0..5 {
            stats.record_consensus();
        }
        assert_eq!(stats.consensus_reads, 5);
    }

    #[test]
    fn test_multiple_input_records() {
        let mut stats = ConsensusCallingStats::new();
        stats.record_input(10);
        stats.record_input(20);
        stats.record_input(30);
        assert_eq!(stats.total_reads, 60);
    }

    #[test]
    fn test_rejection_reason_equality() {
        assert_eq!(RejectionReason::QualityTooLow, RejectionReason::QualityTooLow);
        assert_ne!(RejectionReason::QualityTooLow, RejectionReason::InsufficientReads);
    }

    #[test]
    fn test_stats_new_vs_default() {
        let stats_new = ConsensusCallingStats::new();
        let stats_default = ConsensusCallingStats::default();

        assert_eq!(stats_new.total_reads, stats_default.total_reads);
        assert_eq!(stats_new.consensus_reads, stats_default.consensus_reads);
        assert_eq!(stats_new.filtered_reads, stats_default.filtered_reads);
    }

    #[test]
    fn test_all_rejection_reasons_unique_codes() {
        let reasons = [
            RejectionReason::FragmentRead,
            RejectionReason::InsufficientReads,
            RejectionReason::QualityTooLow,
            RejectionReason::Unmapped,
            RejectionReason::Mapped,
            RejectionReason::TooManyNs,
            RejectionReason::MinorityAlignment,
            RejectionReason::SecondaryOrSupplementary,
            RejectionReason::FailedQC,
            RejectionReason::MissingUmi,
            RejectionReason::QualityTrimmed,
            RejectionReason::ZeroLengthAfterTrimming,
            RejectionReason::InsufficientOverlap,
            RejectionReason::OrphanConsensus,
            RejectionReason::IndelErrorBetweenStrands,
            RejectionReason::PotentialCollision,
            RejectionReason::Other,
        ];

        let codes: Vec<char> = reasons.iter().map(super::RejectionReason::code).collect();
        let mut unique_codes = codes.clone();
        unique_codes.sort_unstable();
        unique_codes.dedup();

        assert_eq!(codes.len(), unique_codes.len(), "All rejection reason codes should be unique");
    }

    #[test]
    fn test_stats_clone() {
        let mut stats1 = ConsensusCallingStats::new();
        stats1.record_input(10);
        stats1.record_consensus();
        stats1.record_rejection(RejectionReason::QualityTooLow, 3);

        let stats2 = stats1.clone();

        assert_eq!(stats1.total_reads, stats2.total_reads);
        assert_eq!(stats1.consensus_reads, stats2.consensus_reads);
        assert_eq!(stats1.filtered_reads, stats2.filtered_reads);
        assert_eq!(stats1.rejection_reasons.len(), stats2.rejection_reasons.len());
    }

    #[test]
    fn test_rejection_reason_copy() {
        let reason1 = RejectionReason::QualityTooLow;
        let reason2 = reason1;
        assert_eq!(reason1, reason2);
    }

    /// Test `make_prefix_from_header` with a single read group with library
    #[test]
    fn test_make_prefix_from_header_single_rg_with_library() {
        use bstr::BString;
        use noodles::sam::header::Builder as HeaderBuilder;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReadGroup;
        use noodles::sam::header::record::value::map::read_group::tag as rg_tag;

        let rg = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("MyLibrary"))
            .build()
            .unwrap();
        let header = HeaderBuilder::default().add_read_group(BString::from("RG1"), rg).build();

        let prefix = make_prefix_from_header(&header);
        assert_eq!(prefix, "MyLibrary");
    }

    /// Test `make_prefix_from_header` with a single read group without library
    #[test]
    fn test_make_prefix_from_header_single_rg_no_library() {
        use bstr::BString;
        use noodles::sam::header::Builder as HeaderBuilder;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReadGroup;

        let rg = Map::<ReadGroup>::builder().build().unwrap();
        let header = HeaderBuilder::default().add_read_group(BString::from("RG1"), rg).build();

        let prefix = make_prefix_from_header(&header);
        assert_eq!(prefix, "RG1");
    }

    /// Test `make_prefix_from_header` with multiple read groups with libraries
    #[test]
    fn test_make_prefix_from_header_multiple_rg_with_libraries() {
        use bstr::BString;
        use noodles::sam::header::Builder as HeaderBuilder;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReadGroup;
        use noodles::sam::header::record::value::map::read_group::tag as rg_tag;

        let rg1 = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibB"))
            .build()
            .unwrap();
        let rg2 = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .unwrap();
        let header = HeaderBuilder::default()
            .add_read_group(BString::from("RG1"), rg1)
            .add_read_group(BString::from("RG2"), rg2)
            .build();

        let prefix = make_prefix_from_header(&header);
        // Libraries should be sorted and joined with |
        assert_eq!(prefix, "LibA|LibB");
    }

    /// Test `make_prefix_from_header` with duplicate libraries (should deduplicate)
    #[test]
    fn test_make_prefix_from_header_deduplicate_libraries() {
        use bstr::BString;
        use noodles::sam::header::Builder as HeaderBuilder;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReadGroup;
        use noodles::sam::header::record::value::map::read_group::tag as rg_tag;

        let rg1 = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("SameLib"))
            .build()
            .unwrap();
        let rg2 = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("SameLib"))
            .build()
            .unwrap();
        let header = HeaderBuilder::default()
            .add_read_group(BString::from("RG1"), rg1)
            .add_read_group(BString::from("RG2"), rg2)
            .build();

        let prefix = make_prefix_from_header(&header);
        // Duplicates should be removed
        assert_eq!(prefix, "SameLib");
    }

    /// Test `make_prefix_from_header` with empty header (no read groups)
    #[test]
    fn test_make_prefix_from_header_empty() {
        let header = Header::default();
        let prefix = make_prefix_from_header(&header);
        // Empty string when no read groups
        assert_eq!(prefix, "");
    }
}
