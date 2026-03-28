//! Parallel pipeline configuration and factory builders.

use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result, anyhow};
use fgumi_consensus::filter::{FilterConfig, FilterThresholds};
use fgumi_lib::consensus_caller::ConsensusCaller;
use fgumi_lib::grouper::GroupFilterConfig;
use fgumi_lib::vanilla_consensus_caller::VanillaUmiConsensusCaller;
use fgumi_pipeline::pipeline::builder::{PipelineConfig, Stage as PipelineStage};
use fgumi_pipeline::stages::consensus::CallerFactory;

use super::Runall;
use super::options::{ConsensusMode, StopAfter};

impl Runall {
    pub(super) fn should_use_parallel_pipeline(&self) -> bool {
        // Parallel path requires enough threads to overcome queue overhead.
        // With larger queue capacity (4096 slots) and exponential backoff in
        // push, the push-spin contention at low thread counts is eliminated,
        // so 4 threads is sufficient for the parallel path to outperform serial.
        let enough_threads = self.threads >= 4;
        let stop_ok =
            matches!(self.stop_after, None | Some(StopAfter::Consensus | StopAfter::Filter));
        enough_threads && stop_ok
    }

    /// Build a [`FilterConfig`] from the CLI filter options for the current consensus mode.
    pub(super) fn build_filter_config(&self) -> FilterConfig {
        let thresholds = FilterThresholds {
            min_reads: self.filter_opts.filter_min_reads,
            max_read_error_rate: 1.0, // no read-level error rate filter by default
            max_base_error_rate: 1.0, // no base-level error rate filter by default
        };
        FilterConfig::for_single_strand(
            thresholds,
            Some(self.filter_opts.filter_min_base_quality),
            None, // no min_mean_base_quality
            1.0,  // max_no_call_fraction = 1.0 (no filter)
        )
    }

    /// Build a [`GroupFilterConfig`] from the CLI group options.
    ///
    /// This controls template filtering before UMI assignment (MAPQ, unmapped,
    /// non-PF, UMI validation), matching the standalone `fgumi group` command.
    pub(super) fn build_group_filter_config(&self) -> GroupFilterConfig {
        let umi_tag = super::tag_bytes(&self.group_opts.group_umi_tag);
        GroupFilterConfig {
            umi_tag,
            min_mapq: self.group_opts.group_min_mapq,
            include_non_pf: false,
            min_umi_length: None,
            no_umi: false,
            allow_unmapped: self.group_opts.group_allow_unmapped,
        }
    }

    /// Build a [`CallerFactory`] closure that creates per-thread consensus callers.
    ///
    /// The factory captures the CLI options and creates a fresh caller for each worker
    /// thread in the parallel pipeline. The `read_name_prefix` should be derived from
    /// the input BAM header via [`fgumi_lib::consensus_caller::make_prefix_from_header`]
    /// to match the serial path.
    pub(super) fn build_caller_factory(&self, read_name_prefix: String) -> CallerFactory {
        let consensus = self.consensus_mode;
        let read_group_id = self.consensus_opts.consensus_read_group_id.clone();
        let simplex_options = self.consensus_options();
        #[cfg(feature = "duplex")]
        let error_rate_pre_umi = self.consensus_opts.consensus_error_rate_pre_umi;
        #[cfg(feature = "duplex")]
        let error_rate_post_umi = self.consensus_opts.consensus_error_rate_post_umi;
        #[cfg(feature = "duplex")]
        let min_input_base_quality = self.filter_opts.filter_min_base_quality;
        #[cfg(feature = "duplex")]
        let duplex_min_reads = self.duplex_opts.duplex_min_reads.clone().unwrap_or_else(|| vec![1]);
        #[cfg(feature = "duplex")]
        let duplex_max_reads_per_strand = self.duplex_opts.duplex_max_reads_per_strand;
        #[cfg(feature = "codec")]
        let codec_min_reads = self.codec_opts.codec_min_reads;
        #[cfg(feature = "codec")]
        let codec_max_reads = self.codec_opts.codec_max_reads;
        #[cfg(feature = "codec")]
        let codec_min_duplex_length = self.codec_opts.codec_min_duplex_length;
        #[cfg(feature = "codec")]
        let codec_min_input_base_quality = self.filter_opts.filter_min_base_quality;
        #[cfg(feature = "codec")]
        let codec_error_rate_pre_umi = self.consensus_opts.consensus_error_rate_pre_umi;
        #[cfg(feature = "codec")]
        let codec_error_rate_post_umi = self.consensus_opts.consensus_error_rate_post_umi;

        Arc::new(move || -> Box<dyn ConsensusCaller> {
            match consensus {
                ConsensusMode::Simplex => {
                    Box::new(VanillaUmiConsensusCaller::new_with_rejects_tracking(
                        read_name_prefix.clone(),
                        read_group_id.clone(),
                        simplex_options.clone(),
                        false,
                    ))
                }
                #[cfg(feature = "duplex")]
                ConsensusMode::Duplex => {
                    use fgumi_lib::duplex_consensus_caller::DuplexConsensusCaller;
                    Box::new(
                        DuplexConsensusCaller::new(
                            read_name_prefix.clone(),
                            read_group_id.clone(),
                            duplex_min_reads.clone(),
                            min_input_base_quality,
                            true,
                            false,
                            duplex_max_reads_per_strand,
                            None,
                            false,
                            error_rate_pre_umi,
                            error_rate_post_umi,
                        )
                        .expect("Failed to create duplex consensus caller in factory"),
                    )
                }
                #[cfg(feature = "codec")]
                ConsensusMode::Codec => {
                    use fgumi_lib::consensus::codec_caller::{
                        CodecConsensusCaller, CodecConsensusOptions,
                    };
                    let options = CodecConsensusOptions {
                        min_input_base_quality: codec_min_input_base_quality,
                        error_rate_pre_umi: codec_error_rate_pre_umi,
                        error_rate_post_umi: codec_error_rate_post_umi,
                        min_reads_per_strand: codec_min_reads,
                        max_reads_per_strand: codec_max_reads,
                        min_duplex_length: codec_min_duplex_length,
                        produce_per_base_tags: true,
                        trim: false,
                        min_consensus_base_quality: 0,
                        ..CodecConsensusOptions::default()
                    };
                    Box::new(CodecConsensusCaller::new_with_rejects_tracking(
                        read_name_prefix.clone(),
                        read_group_id.clone(),
                        options,
                        false,
                    ))
                }
            }
        })
    }

    /// Build a [`PipelineConfig`] for the parallel pipeline builder.
    ///
    /// Reads the input BAM header to derive the read name prefix (matching the serial path),
    /// then assembles the full pipeline configuration.
    ///
    /// # Arguments
    ///
    /// * `input` - Input BAM file path (unsorted or template-coordinate sorted).
    /// * `output` - Output BAM file path.
    /// * `stage` - Which stage to start from (`Sort` or `Group`).
    pub(super) fn build_pipeline_config(
        &self,
        input: &Path,
        output: &Path,
        stage: PipelineStage,
    ) -> Result<PipelineConfig> {
        use fgumi_lib::bam_io::create_raw_bam_reader;

        let (_reader, header) = create_raw_bam_reader(input, 1)
            .with_context(|| format!("opening input BAM for header: {}", input.display()))?;

        self.build_pipeline_config_from_header(&header, input, output, stage)
    }

    /// Build a [`PipelineConfig`] from a pre-read header.
    ///
    /// Like [`build_pipeline_config`](Self::build_pipeline_config) but takes the header
    /// directly, avoiding the need to open a BAM file. The `input` path is stored in the
    /// config but may be unused if the caller feeds sorted chunks directly via
    /// [`PipelineBuilder::run_from_sorted_chunks`].
    pub(super) fn build_pipeline_config_from_header(
        &self,
        header: &noodles::sam::Header,
        input: &Path,
        output: &Path,
        stage: PipelineStage,
    ) -> Result<PipelineConfig> {
        let reference = self
            .reference
            .as_ref()
            .ok_or_else(|| anyhow!("--reference is required for parallel pipeline"))?;

        let read_name_prefix = self
            .consensus_opts
            .consensus_read_name_prefix
            .clone()
            .unwrap_or_else(|| fgumi_lib::consensus_caller::make_prefix_from_header(header));

        let umi_tag_bytes = super::tag_bytes(&self.group_opts.group_umi_tag);
        let assign_tag_bytes = super::tag_bytes(&self.consensus_opts.consensus_mi_tag);

        let assigner = self.group_opts.group_strategy.new_assigner(self.group_opts.group_max_edits);

        Ok(PipelineConfig {
            input: input.to_path_buf(),
            output: output.to_path_buf(),
            reference: reference.clone(),
            threads: self.threads,
            start_from: stage,
            assigner: Arc::from(assigner),
            umi_tag: umi_tag_bytes,
            assign_tag: assign_tag_bytes,
            group_filter_config: self.build_group_filter_config(),
            caller_factory: self.build_caller_factory(read_name_prefix),
            filter_config: self.build_filter_config(),
            min_base_quality: self.filter_opts.filter_min_base_quality,
            compression_level: self.compression_level,
            sort_memory_limit: self.sort_opts.sort_memory_limit,
            sort_temp_dir: self.sort_opts.sort_temp_dir.clone(),
        })
    }
}
