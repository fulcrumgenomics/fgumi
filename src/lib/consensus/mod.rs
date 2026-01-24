//! Consensus calling and filtering for UMI-based molecular consensus reads
//!
//! This module provides comprehensive functionality for generating consensus sequences
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
pub mod codec_caller;
pub mod duplex_caller;
pub mod filter;
pub mod overlapping;
pub mod sequence;
pub mod simple_umi;
pub mod tags;
pub mod vanilla_caller;

// Re-export commonly used items
pub use base_builder::ConsensusBaseBuilder;
pub use caller::{
    ConsensusCaller, ConsensusOptionsBase, RejectionTracker, calculate_error_rate,
    log_consensus_statistics,
};
pub use codec_caller::{CodecConsensusCaller, CodecConsensusOptions, CodecConsensusStats};
pub use duplex_caller::{DuplexConsensusCaller, DuplexConsensusRead};
pub use filter::{
    ConsensusType, FilterConfig, FilterResult, FilterThresholds, count_no_calls,
    filter_duplex_read, filter_read, is_duplex_consensus, mask_bases, mask_duplex_bases,
    mean_base_quality, template_passes,
};
pub use overlapping::{
    AgreementStrategy, CorrectionStats, DisagreementStrategy, OverlappingBasesConsensusCaller,
    apply_overlapping_consensus,
};
pub use sequence::ConsensusSequence;
pub use vanilla_caller::{
    VanillaConsensusRead, VanillaUmiConsensusCaller, VanillaUmiConsensusOptions,
};

// Re-export crate-internal items for use between consensus modules
pub(crate) use vanilla_caller::{
    IndexedSourceRead, ReadType, SourceRead, select_most_common_alignment_group,
};
