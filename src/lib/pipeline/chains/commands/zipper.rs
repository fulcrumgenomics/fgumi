//! Chain builder for `Stage::Zipper`.
//!
//! Body extracted from `commands::zipper::Zipper::execute_new_pipeline`
//! in Phase 2 (T2.22). Zipper is a two-input command, so the spec carries
//! `SourceSpec::PairedBams { unmapped, mapped, reference }`; this function
//! destructures it and bails with a clear message on any other variant.
//!
//! In Phase 3a T3a.12 the bulk of the chain-construction logic moved into
//! [`ChainBuilder::add_zipper`] and the step-construction details were
//! extracted into [`build_zipper_merge_step`]. `build_zipper_chain` is now
//! the canonical 10-line entry point.

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::info;

use crate::commands::zipper::{ZipperOptions, merge_step};
use crate::logging::OperationTimer;
use crate::pipeline::chains::builder::StagePosition;
use crate::pipeline::chains::{BuiltPipeline, ChainSpec, FinalizeHook, Stage};
use crate::pipeline::steps::tuning::BamPipelineTuning;
use crate::reference::ReferenceReader;

// ─────────────────────────────────────────────────────────────────────────────
// ZipperFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline finalize hook for zipper. Logs the optional
/// `--exclude-missing-reads` summary, logs "zipper completed
/// successfully", and calls `timer.log_completion`.
pub(crate) struct ZipperFinalizeHook {
    pub(crate) missing_count: Arc<AtomicU64>,
    pub(crate) records_emitted: Arc<AtomicU64>,
    pub(crate) exclude_missing_reads: bool,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for ZipperFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let ZipperFinalizeHook { missing_count, records_emitted, exclude_missing_reads, timer } =
            *self;

        let missing = missing_count.load(Ordering::Relaxed);
        if exclude_missing_reads && missing > 0 {
            info!("Excluded {missing} templates that were not present in the aligned BAM.");
        }

        info!("zipper completed successfully");
        timer.log_completion(records_emitted.load(Ordering::Relaxed));
        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// build_zipper_merge_step factory
// ─────────────────────────────────────────────────────────────────────────────

/// Captures for [`build_zipper_merge_step`] construction.
pub(crate) struct ZipperMergeCaptures {
    pub(crate) zipper_opts: ZipperOptions,
    pub(crate) output_header: Arc<noodles::sam::Header>,
    pub(crate) reference_path: std::path::PathBuf,
    pub(crate) tuning: BamPipelineTuning,
    pub(crate) missing_count: Arc<AtomicU64>,
    pub(crate) records_emitted: Arc<AtomicU64>,
}

/// Build a [`ZipperMergeStep`] from the supplied captures.
///
/// Logs tag-manipulation summary lines, loads the reference FASTA if
/// `--restore-unconverted-bases` is set, constructs the
/// [`ZipperMergeConfig`], and returns the configured step.
///
/// [`ZipperMergeStep`]: merge_step::ZipperMergeStep
/// [`ZipperMergeConfig`]: merge_step::ZipperMergeConfig
pub(crate) fn build_zipper_merge_step(
    caps: ZipperMergeCaptures,
) -> Result<merge_step::ZipperMergeStep> {
    use crate::umi::TagInfo;

    let ZipperMergeCaptures {
        zipper_opts,
        output_header,
        reference_path,
        tuning,
        missing_count,
        records_emitted,
    } = caps;

    let tag_info = TagInfo::new(
        zipper_opts.tags_to_remove.clone(),
        zipper_opts.tags_to_reverse.clone(),
        zipper_opts.tags_to_revcomp.clone(),
    );
    if !tag_info.remove.is_empty() {
        info!("Tags for removal: {:?}", tag_info.remove);
    }
    if !tag_info.reverse.is_empty() {
        info!("Tags being reversed: {:?}", tag_info.reverse);
    }
    if !tag_info.revcomp.is_empty() {
        info!("Tags being reverse complemented: {:?}", tag_info.revcomp);
    }

    let reference = if zipper_opts.restore_unconverted_bases {
        info!("Loading reference FASTA for unconverted base restoration");
        Some(Arc::new(ReferenceReader::new(&reference_path)?))
    } else {
        None
    };

    let cfg = merge_step::ZipperMergeConfig {
        tag_info: Arc::new(tag_info),
        skip_tc_tags: zipper_opts.skip_tc_tags,
        exclude_missing_reads: zipper_opts.exclude_missing_reads,
        reference,
        output_header,
        missing_count,
        records_emitted,
        target_batch_count: tuning.template_batch_size,
        output_byte_limit: tuning.per_step_byte_limit,
    };
    Ok(merge_step::ZipperMergeStep::new(cfg))
}

// ─────────────────────────────────────────────────────────────────────────────
// build_zipper_chain — 10-line entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Build a [`BuiltPipeline`] for a zipper-only chain.
///
/// `spec.stages == [Stage::Zipper]`. Other layouts are caller errors
/// caught by the dispatch match in [`crate::pipeline::chains::build`].
///
/// `spec.source` must be `SourceSpec::PairedBams { unmapped, mapped, reference }`.
/// The mapped source supports both BAM and SAM input; the unmapped source
/// is BAM-only (as `fgumi extract` always produces BAM). The function
/// bails with a clear message if the unmapped source turns out to be SAM.
///
/// # Errors
///
/// Returns input-validation errors or any underlying pipeline
/// construction error.
#[allow(clippy::needless_pass_by_value)]
pub fn build_zipper_chain(spec: ChainSpec) -> Result<BuiltPipeline> {
    use crate::pipeline::chains::builder::ChainBuilder;

    let mut chain = ChainBuilder::new(&spec)?;
    chain.add_source()?;
    chain.add_stage(Stage::Zipper, StagePosition::Terminal)?;
    chain.add_sink()?;
    chain.build()
}
