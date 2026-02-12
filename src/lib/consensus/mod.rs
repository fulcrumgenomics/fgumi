//! Consensus calling and filtering for UMI-based molecular consensus reads.
//!
//! This module re-exports from the `fgumi-consensus` crate for backward compatibility.

pub use fgumi_consensus::{base_builder, caller, filter, overlapping, sequence, simple_umi, tags};

pub use fgumi_consensus::{codec_caller, duplex_caller, vanilla_caller};

// Re-export commonly used items
pub use fgumi_consensus::{
    AgreementStrategy, CodecConsensusCaller, CodecConsensusOptions, CodecConsensusStats,
    ConsensusBaseBuilder, ConsensusCaller, ConsensusOptionsBase, ConsensusSequence, ConsensusType,
    CorrectionStats, DisagreementStrategy, DuplexConsensusCaller, DuplexConsensusRead,
    FilterConfig, FilterResult, FilterThresholds, OverlappingBasesConsensusCaller,
    RejectionTracker, VanillaConsensusRead, VanillaUmiConsensusCaller, VanillaUmiConsensusOptions,
    apply_overlapping_consensus, calculate_error_rate, compute_read_stats, count_no_calls,
    filter_duplex_read, filter_read, is_duplex_consensus, log_consensus_statistics, mask_bases,
    mask_duplex_bases, mean_base_quality, template_passes,
};
