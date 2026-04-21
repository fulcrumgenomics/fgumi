//! Consensus caller options and factory dispatch for the runall pipeline.
//!
//! Hosts two concerns decoupled from the `Runall` CLI struct:
//!
//! 1. `Runall::consensus_options` — builds [`VanillaUmiConsensusOptions`] from
//!    CLI flags so runall and standalone `fgumi simplex` produce identical
//!    output for the same inputs.
//! 2. [`ConsensusFactoryOptions`] and `ConsensusMode::build_factory` — per-mode
//!    construction of a per-worker [`CallerFactory`]. The `ConsensusMode` enum
//!    lives next to the CLI but dispatches without touching the CLI type, so
//!    adding a new consensus mode is a single-file change to this module and
//!    the `ConsensusMode` enum.

use std::sync::Arc;

use anyhow::Result;

use crate::runall::engine::stages::consensus::CallerFactory;
use crate::vanilla_consensus_caller::VanillaUmiConsensusOptions;

use super::Runall;
use super::options::ConsensusMode;

impl Runall {
    /// Resolve the methylation mode from the CLI flag.
    pub(super) fn resolved_methylation_mode(&self) -> crate::consensus::MethylationMode {
        crate::commands::common::resolve_methylation_mode(self.methylation_mode)
    }

    /// Build [`VanillaUmiConsensusOptions`] from CLI flags, mirroring the
    /// defaults the standalone `fgumi simplex` command uses so runall and
    /// standalone produce identical output for the same inputs.
    pub(crate) fn consensus_options(&self) -> VanillaUmiConsensusOptions {
        VanillaUmiConsensusOptions {
            tag: crate::sam::SamTag::MI.to_string(),
            error_rate_pre_umi: self.consensus_opts.consensus_error_rate_pre_umi,
            error_rate_post_umi: self.consensus_opts.consensus_error_rate_post_umi,
            min_input_base_quality: self.consensus_opts.consensus_min_input_base_quality,
            // Single-strand consensus uses the first (most-stringent) value;
            // duplex/AB/BA splits are handled downstream in the filter stage.
            // Validate ensures a value exists for real runs; fallback 0 is a
            // placeholder for `--explain` which never executes the stage.
            min_reads: self
                .filter_opts
                .filter_min_reads
                .as_deref()
                .and_then(<[usize]>::first)
                .copied()
                .unwrap_or(0),
            max_reads: None,
            produce_per_base_tags: self.consensus_opts.consensus_output_per_base_tags,
            seed: Some(42),
            trim: self.consensus_opts.consensus_trim,
            min_consensus_base_quality: self.consensus_opts.consensus_min_consensus_base_quality,
            cell_tag: None,
            methylation_mode: self.resolved_methylation_mode(),
        }
    }

    /// Resolve the header-derived read-name prefix (falling back to a CLI
    /// override) and the read-group ID used by every consensus caller.
    pub(crate) fn build_consensus_factory_options(
        &self,
        header: &noodles::sam::Header,
    ) -> ConsensusFactoryOptions {
        use crate::consensus_caller::make_prefix_from_header;

        let read_name_prefix = self
            .consensus_opts
            .consensus_read_name_prefix
            .clone()
            .unwrap_or_else(|| make_prefix_from_header(header));
        let read_group_id = self.consensus_opts.consensus_read_group_id.clone();

        ConsensusFactoryOptions {
            read_name_prefix,
            read_group_id,
            simplex: SimplexFactoryOptions { options: self.consensus_options() },
            duplex: DuplexFactoryOptions {
                min_reads: self.duplex_opts.duplex_min_reads.clone().unwrap_or_else(|| vec![1]),
                min_input_base_quality: self.consensus_opts.consensus_min_input_base_quality,
                output_per_base_tags: self.consensus_opts.consensus_output_per_base_tags,
                trim: self.consensus_opts.consensus_trim,
                max_reads_per_strand: self.duplex_opts.duplex_max_reads_per_strand,
                error_rate_pre_umi: self.consensus_opts.consensus_error_rate_pre_umi,
                error_rate_post_umi: self.consensus_opts.consensus_error_rate_post_umi,
            },
            codec: CodecFactoryOptions {
                min_input_base_quality: self.consensus_opts.consensus_min_input_base_quality,
                error_rate_pre_umi: self.consensus_opts.consensus_error_rate_pre_umi,
                error_rate_post_umi: self.consensus_opts.consensus_error_rate_post_umi,
                min_reads_per_strand: self.codec_opts.codec_min_reads,
                max_reads_per_strand: self.codec_opts.codec_max_reads,
                min_duplex_length: self.codec_opts.codec_min_duplex_length,
                output_per_base_tags: self.consensus_opts.consensus_output_per_base_tags,
                trim: self.consensus_opts.consensus_trim,
                min_consensus_base_quality: self
                    .consensus_opts
                    .consensus_min_consensus_base_quality,
            },
        }
    }

    /// Build a per-worker [`CallerFactory`] for the pipeline.
    ///
    /// Thin wrapper: aggregates CLI args into [`ConsensusFactoryOptions`] and
    /// dispatches via `ConsensusMode::build_factory` so the CLI type stays
    /// decoupled from per-mode construction.
    pub(crate) fn build_consensus_factory(
        &self,
        header: &noodles::sam::Header,
    ) -> Result<CallerFactory> {
        let opts = self.build_consensus_factory_options(header);
        self.consensus_mode.build_factory(&opts)
    }
}

/// Per-mode inputs for constructing a [`CallerFactory`], decoupled from the
/// runall CLI struct. Aggregating into one struct keeps
/// [`ConsensusMode::build_factory`]'s signature stable as new per-mode knobs
/// are added.
#[derive(Debug, Clone)]
pub(crate) struct ConsensusFactoryOptions {
    /// Prefix applied to consensus read names.
    pub read_name_prefix: String,
    /// Read-group ID stamped on every consensus record.
    pub read_group_id: String,
    /// Simplex-mode inputs.
    pub simplex: SimplexFactoryOptions,
    /// Duplex-mode inputs.
    pub duplex: DuplexFactoryOptions,
    /// CODEC-mode inputs.
    pub codec: CodecFactoryOptions,
}

/// Simplex-caller construction inputs.
#[derive(Debug, Clone)]
pub(crate) struct SimplexFactoryOptions {
    pub options: VanillaUmiConsensusOptions,
}

/// Duplex-caller construction inputs.
#[derive(Debug, Clone)]
pub(crate) struct DuplexFactoryOptions {
    pub min_reads: Vec<usize>,
    pub min_input_base_quality: u8,
    pub output_per_base_tags: bool,
    pub trim: bool,
    pub max_reads_per_strand: Option<usize>,
    pub error_rate_pre_umi: u8,
    pub error_rate_post_umi: u8,
}

/// CODEC-caller construction inputs.
#[derive(Debug, Clone)]
pub(crate) struct CodecFactoryOptions {
    pub min_input_base_quality: u8,
    pub error_rate_pre_umi: u8,
    pub error_rate_post_umi: u8,
    pub min_reads_per_strand: usize,
    pub max_reads_per_strand: Option<usize>,
    pub min_duplex_length: usize,
    pub output_per_base_tags: bool,
    pub trim: bool,
    pub min_consensus_base_quality: u8,
}

impl ConsensusMode {
    /// Build a per-worker [`CallerFactory`] for this consensus mode.
    ///
    /// Dispatches on the enum variant and constructs a closure that produces a
    /// fresh caller per pool worker. Adding a new consensus mode only requires
    /// adding a variant and a branch here.
    pub(crate) fn build_factory(self, opts: &ConsensusFactoryOptions) -> Result<CallerFactory> {
        match self {
            Self::Simplex => {
                use crate::vanilla_consensus_caller::VanillaUmiConsensusCaller;
                let read_name_prefix = opts.read_name_prefix.clone();
                let read_group_id = opts.read_group_id.clone();
                let options = opts.simplex.options.clone();
                let factory: CallerFactory = Arc::new(move || {
                    let rnp = read_name_prefix.clone();
                    let rgid = read_group_id.clone();
                    let opts = options.clone();
                    Box::new(VanillaUmiConsensusCaller::new(rnp, rgid, opts))
                });
                Ok(factory)
            }
            Self::Duplex => {
                use crate::duplex_consensus_caller::DuplexConsensusCaller;
                let read_name_prefix = opts.read_name_prefix.clone();
                let read_group_id = opts.read_group_id.clone();
                let duplex = opts.duplex.clone();
                // Probe-build once to surface any constructor error before
                // the factory is handed to worker threads. `min_reads` is the
                // only input `new` validates, and it's captured by value —
                // the per-worker rebuild below cannot see different inputs.
                validate_duplex_min_reads(&duplex.min_reads)?;
                let factory: CallerFactory = Arc::new(move || {
                    Box::new(
                        DuplexConsensusCaller::new(
                            read_name_prefix.clone(),
                            read_group_id.clone(),
                            duplex.min_reads.clone(),
                            duplex.min_input_base_quality,
                            duplex.output_per_base_tags,
                            duplex.trim,
                            duplex.max_reads_per_strand,
                            None,
                            false,
                            duplex.error_rate_pre_umi,
                            duplex.error_rate_post_umi,
                        )
                        .expect("validate_duplex_min_reads already rejected invalid input"),
                    )
                });
                Ok(factory)
            }
            Self::Codec => {
                use crate::consensus::codec_caller::{CodecConsensusCaller, CodecConsensusOptions};
                let read_name_prefix = opts.read_name_prefix.clone();
                let read_group_id = opts.read_group_id.clone();
                let codec = &opts.codec;
                let codec_opts = CodecConsensusOptions {
                    min_input_base_quality: codec.min_input_base_quality,
                    error_rate_pre_umi: codec.error_rate_pre_umi,
                    error_rate_post_umi: codec.error_rate_post_umi,
                    min_reads_per_strand: codec.min_reads_per_strand,
                    max_reads_per_strand: codec.max_reads_per_strand,
                    min_duplex_length: codec.min_duplex_length,
                    produce_per_base_tags: codec.output_per_base_tags,
                    trim: codec.trim,
                    min_consensus_base_quality: codec.min_consensus_base_quality,
                    ..CodecConsensusOptions::default()
                };
                let factory: CallerFactory = Arc::new(move || {
                    Box::new(CodecConsensusCaller::new(
                        read_name_prefix.clone(),
                        read_group_id.clone(),
                        codec_opts.clone(),
                    ))
                });
                Ok(factory)
            }
        }
    }
}

/// Validate a duplex `--min-reads` vector the way [`DuplexConsensusCaller::new`]
/// does, but without constructing a caller instance. Used to hoist validation
/// out of the per-worker factory closure so the closure's caller-construction
/// can panic-assert correctness instead of returning a `Result`.
///
/// [`DuplexConsensusCaller::new`]: crate::duplex_consensus_caller::DuplexConsensusCaller::new
pub(crate) fn validate_duplex_min_reads(min_reads: &[usize]) -> Result<()> {
    if min_reads.is_empty() {
        anyhow::bail!("min_reads parameter must have at least 1 value");
    }
    if min_reads.len() > 3 {
        anyhow::bail!(
            "min_reads parameter must have 1-3 values (total, [XY, [YX]]), got {} values",
            min_reads.len()
        );
    }
    let last = *min_reads.last().expect("checked non-empty above");
    let total = min_reads[0];
    let xy = min_reads.get(1).copied().unwrap_or(last);
    let yx = min_reads.get(2).copied().unwrap_or(last);
    if xy > total {
        anyhow::bail!("min-reads values must be specified high to low (total >= XY)");
    }
    if yx > xy {
        anyhow::bail!("min-reads values must be specified high to low (XY >= YX)");
    }
    Ok(())
}
