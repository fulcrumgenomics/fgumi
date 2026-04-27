//! Consensus calling and filtering for UMI-based molecular consensus reads.
//!
//! This module re-exports from the `fgumi-consensus` crate for backward compatibility.

pub use fgumi_consensus::{
    MethylationMode, base_builder, caller, filter, methylation, overlapping, sequence, simple_umi,
    tags,
};

pub use fgumi_consensus::codec_caller;
pub use fgumi_consensus::duplex_caller;
pub use fgumi_consensus::vanilla_caller;

// Re-export commonly used items
pub use fgumi_consensus::{
    AgreementStrategy, ConsensusBaseBuilder, ConsensusCaller, ConsensusSequence, ConsensusType,
    CorrectionStats, DisagreementStrategy, FilterConfig, FilterResult, FilterThresholds,
    OverlappingBasesConsensusCaller, apply_overlapping_consensus, calculate_error_rate,
    compute_read_stats, count_no_calls, filter_duplex_read, filter_read, is_duplex_consensus,
    log_consensus_statistics, mask_bases, mask_duplex_bases, mean_base_quality, template_passes,
};

pub use fgumi_consensus::{
    VanillaConsensusRead, VanillaUmiConsensusCaller, VanillaUmiConsensusOptions,
};

pub use fgumi_consensus::{DuplexConsensusCaller, DuplexConsensusRead};

pub use fgumi_consensus::{CodecConsensusCaller, CodecConsensusOptions, CodecConsensusStats};
