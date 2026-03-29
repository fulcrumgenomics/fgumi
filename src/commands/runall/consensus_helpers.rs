//! Consensus caller construction helpers.

use anyhow::Result;
use fgumi_lib::consensus_caller::ConsensusCaller;
use fgumi_lib::vanilla_consensus_caller::{VanillaUmiConsensusCaller, VanillaUmiConsensusOptions};

use super::Runall;
use super::options::ConsensusMode;

impl Runall {
    pub(super) fn consensus_options(&self) -> VanillaUmiConsensusOptions {
        VanillaUmiConsensusOptions {
            tag: self.consensus_opts.consensus_mi_tag.clone(),
            error_rate_pre_umi: self.consensus_opts.consensus_error_rate_pre_umi,
            error_rate_post_umi: self.consensus_opts.consensus_error_rate_post_umi,
            min_input_base_quality: self.filter_opts.filter_min_base_quality,
            min_reads: self.filter_opts.filter_min_reads,
            max_reads: None,
            produce_per_base_tags: true,
            seed: Some(42),
            trim: false,
            min_consensus_base_quality: 40,
            cell_tag: None,
        }
    }

    /// Build the appropriate consensus caller based on the selected consensus mode.
    ///
    /// Returns a boxed trait object that handles consensus calling regardless of mode.
    pub(super) fn build_consensus_caller(
        &self,
        read_name_prefix: String,
    ) -> Result<Box<dyn ConsensusCaller>> {
        match self.consensus_mode {
            ConsensusMode::Simplex => {
                let options = self.consensus_options();
                let caller = VanillaUmiConsensusCaller::new_with_rejects_tracking(
                    read_name_prefix,
                    self.consensus_opts.consensus_read_group_id.clone(),
                    options,
                    false,
                );
                Ok(Box::new(caller))
            }
            #[cfg(feature = "duplex")]
            ConsensusMode::Duplex => {
                use fgumi_lib::duplex_consensus_caller::DuplexConsensusCaller;
                let caller = DuplexConsensusCaller::new(
                    read_name_prefix,
                    self.consensus_opts.consensus_read_group_id.clone(),
                    self.duplex_opts.duplex_min_reads.clone().unwrap_or_else(|| vec![1]),
                    self.filter_opts.filter_min_base_quality,
                    true,  // produce_per_base_tags
                    false, // trim
                    self.duplex_opts.duplex_max_reads_per_strand,
                    None,  // cell_tag
                    false, // track_rejects
                    self.consensus_opts.consensus_error_rate_pre_umi,
                    self.consensus_opts.consensus_error_rate_post_umi,
                )?;
                Ok(Box::new(caller))
            }
            #[cfg(feature = "codec")]
            ConsensusMode::Codec => {
                use fgumi_lib::consensus::codec_caller::{
                    CodecConsensusCaller, CodecConsensusOptions,
                };
                let options = CodecConsensusOptions {
                    min_input_base_quality: self.filter_opts.filter_min_base_quality,
                    error_rate_pre_umi: self.consensus_opts.consensus_error_rate_pre_umi,
                    error_rate_post_umi: self.consensus_opts.consensus_error_rate_post_umi,
                    min_reads_per_strand: self.codec_opts.codec_min_reads,
                    max_reads_per_strand: self.codec_opts.codec_max_reads,
                    min_duplex_length: self.codec_opts.codec_min_duplex_length,
                    produce_per_base_tags: true,
                    trim: false,
                    min_consensus_base_quality: 0,
                    ..CodecConsensusOptions::default()
                };
                let caller = CodecConsensusCaller::new_with_rejects_tracking(
                    read_name_prefix,
                    self.consensus_opts.consensus_read_group_id.clone(),
                    options,
                    false,
                );
                Ok(Box::new(caller))
            }
        }
    }
}
