#![deny(unsafe_code)]

//! Consensus calling and filtering for UMI-based molecular consensus reads
//!
//! This crate provides comprehensive functionality for generating consensus sequences
//! from reads grouped by Unique Molecular Identifiers (UMIs). It includes:
//!
//! - **Base-level consensus calling**: Building consensus bases from multiple reads
//! - **Duplex consensus**: Calling duplex consensus from paired single-strand consensuses
//! - **Simple UMI consensus**: Fast consensus for non-overlapping reads
//! - **Vanilla consensus**: Standard consensus without special features
//! - **Overlapping consensus**: Handling overlapping read pairs
//! - **Consensus filtering**: Quality-based filtering and masking of consensus reads
//! - **Consensus tags**: SAM tags for tracking consensus metrics

pub mod base_builder;
pub mod caller;
pub mod filter;
pub mod overlapping;
pub mod phred;
pub mod sequence;
pub mod simple_umi;
pub mod tags;

#[cfg(feature = "simplex")]
pub mod vanilla_caller;

#[cfg(feature = "duplex")]
pub mod duplex_caller;

#[cfg(feature = "codec")]
pub mod codec_caller;

#[cfg(feature = "simplex")]
pub mod methylation;

mod vendored;

/// Methylation chemistry mode for consensus calling.
///
/// Controls how C→T conversions at reference cytosine positions are interpreted
/// during consensus calling and MM/ML tag generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MethylationMode {
    /// No methylation-aware processing. C→T changes are treated as errors.
    #[default]
    Disabled,
    /// EM-Seq (enzymatic methyl-seq): unmethylated C is converted to T.
    /// C in read = methylated (protected), T in read = unmethylated (converted).
    EmSeq,
    /// TAPs/Illumina 5-base: methylated C is converted to T.
    /// C in read = unmethylated (not a target), T in read = methylated (converted).
    Taps,
}

impl MethylationMode {
    /// Returns true if any methylation mode is enabled.
    #[must_use]
    pub fn is_enabled(&self) -> bool {
        !matches!(self, Self::Disabled)
    }
}

// Re-export commonly used items
pub use base_builder::ConsensusBaseBuilder;
pub use caller::{ConsensusCaller, calculate_error_rate, log_consensus_statistics};
pub use filter::{
    ConsensusType, FilterConfig, FilterResult, FilterThresholds, compute_read_stats,
    count_no_calls, filter_duplex_read, filter_read, is_duplex_consensus, mask_bases,
    mask_duplex_bases, mean_base_quality, template_passes,
};
pub use overlapping::{
    AgreementStrategy, CorrectionStats, DisagreementStrategy, OverlappingBasesConsensusCaller,
    apply_overlapping_consensus,
};
pub use sequence::ConsensusSequence;

#[cfg(feature = "simplex")]
pub use vanilla_caller::{
    VanillaConsensusRead, VanillaUmiConsensusCaller, VanillaUmiConsensusOptions,
};

#[cfg(feature = "simplex")]
pub(crate) use vanilla_caller::{
    IndexedSourceRead, ReadType, SourceRead, select_most_common_alignment_group,
};

#[cfg(feature = "duplex")]
pub use duplex_caller::{DuplexConsensusCaller, DuplexConsensusRead};

#[cfg(feature = "codec")]
pub use codec_caller::{CodecConsensusCaller, CodecConsensusOptions, CodecConsensusStats};
