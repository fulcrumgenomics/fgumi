//! Consensus calling and filtering for UMI-based molecular consensus reads.
//!
//! This module re-exports from the `fgumi-consensus` crate for backward compatibility.

pub use fgumi_consensus::{base_builder, caller, filter, overlapping, sequence, simple_umi, tags};

#[cfg(feature = "codec")]
pub use fgumi_consensus::codec_caller;
#[cfg(feature = "duplex")]
pub use fgumi_consensus::duplex_caller;
#[cfg(feature = "simplex")]
pub use fgumi_consensus::vanilla_caller;

// Re-export commonly used items
pub use fgumi_consensus::{
    AgreementStrategy, ConsensusBaseBuilder, ConsensusCaller, ConsensusOptionsBase,
    ConsensusSequence, ConsensusType, CorrectionStats, DisagreementStrategy, FilterConfig,
    FilterResult, FilterThresholds, OverlappingBasesConsensusCaller, RejectionTracker,
    apply_overlapping_consensus, calculate_error_rate, compute_read_stats, count_no_calls,
    filter_duplex_read, filter_read, is_duplex_consensus, log_consensus_statistics, mask_bases,
    mask_duplex_bases, mean_base_quality, template_passes,
};

#[cfg(feature = "simplex")]
pub use fgumi_consensus::{
    VanillaConsensusRead, VanillaUmiConsensusCaller, VanillaUmiConsensusOptions,
};

#[cfg(feature = "duplex")]
pub use fgumi_consensus::{DuplexConsensusCaller, DuplexConsensusRead};

#[cfg(feature = "codec")]
pub use fgumi_consensus::{CodecConsensusCaller, CodecConsensusOptions, CodecConsensusStats};
