#![deny(unsafe_code)]
#![deny(missing_docs)]

//! Consensus calling and filtering for UMI-tagged sequencing reads.
//!
//! This crate generates consensus sequences from reads grouped by Unique Molecular
//! Identifiers (UMIs). It is consumed by the fgumi CLI and may also be used
//! directly as a library.
//!
//! # Overview
//!
//! The crate is organized around a small number of concepts:
//!
//! - **Base-level consensus** ([`ConsensusBaseBuilder`], [`caller::ConsensusCaller`]):
//!   aggregates evidence at each position to choose a consensus base and quality.
//! - **Consensus callers**: higher-level wrappers that consume groups of reads and
//!   produce consensus reads:
//!   - [`VanillaUmiConsensusCaller`] and [`methylation`] (simplex)
//!   - [`DuplexConsensusCaller`] (duplex)
//!   - [`CodecConsensusCaller`] (CODEC)
//! - **Overlap resolution** ([`overlapping`]): harmonise overlapping read pairs
//!   before consensus.
//! - **Filtering** ([`filter`]): apply quality and agreement thresholds to
//!   consensus output.
//! - **Tags** ([`tags`]): SAM/BAM tag identifiers used throughout.
//!
//! # Example
//!
//! ```
//! use fgumi_consensus::MethylationMode;
//!
//! // Methylation mode defaults to disabled; callers opt into EM-seq or TAPS.
//! assert!(!MethylationMode::default().is_enabled());
//! assert!(MethylationMode::EmSeq.is_enabled());
//! ```
//!
//! # Errors
//!
//! All fallible public functions return [`enum@Error`] via [`Result`]. See
//! [`enum@Error`] for the variants consumers may pattern-match on.

pub mod base_builder;
pub mod caller;
pub mod error;
pub mod filter;
pub mod overlapping;
pub mod phred;
pub mod sequence;
pub mod simple_umi;
pub mod tags;

pub mod vanilla_caller;

pub mod duplex_caller;

pub mod codec_caller;

pub mod methylation;

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
pub use error::{Error, Result};
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

pub use vanilla_caller::{
    VanillaConsensusRead, VanillaUmiConsensusCaller, VanillaUmiConsensusOptions,
};

pub(crate) use vanilla_caller::{
    IndexedSourceRead, ReadType, SourceRead, select_most_common_alignment_group,
};

pub use duplex_caller::{DuplexConsensusCaller, DuplexConsensusRead};

pub use codec_caller::{CodecConsensusCaller, CodecConsensusOptions, CodecConsensusStats};
