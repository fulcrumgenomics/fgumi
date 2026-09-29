//! Chain builder for `Stage::Zipper`.
//!
//! Body extracted from `commands::zipper::Zipper::execute_new_pipeline`
//! in Phase 2 (T2.22). Zipper is a two-input command, so the spec carries
//! `SourceSpec::PairedBams { unmapped, mapped, reference }`; this function
//! destructures it and bails with a clear message on any other variant.
//!
//! In Phase 3a T3a.12 the bulk of the chain-construction logic moved into
//! `ChainBuilder::add_zipper` and the step-construction details were
//! extracted into `build_zipper_merge_config`.

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::info;

use crate::commands::zipper::{ZipperMergeRules, ZipperOptions, merge_step};
use crate::logging::OperationTimer;
use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::steps::tuning::BamPipelineTuning;

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
// build_zipper_merge_config factory
// ─────────────────────────────────────────────────────────────────────────────

/// Captures for [`build_zipper_merge_config`] construction.
pub(crate) struct ZipperMergeCaptures {
    pub(crate) zipper_opts: ZipperOptions,
    pub(crate) output_header: Arc<noodles::sam::Header>,
    pub(crate) reference_path: std::path::PathBuf,
    pub(crate) tuning: BamPipelineTuning,
    pub(crate) missing_count: Arc<AtomicU64>,
    pub(crate) records_emitted: Arc<AtomicU64>,
}

/// Build the [`ZipperMergeConfig`] from the supplied captures.
///
/// Logs tag-manipulation summary lines, loads the reference FASTA if
/// `--restore-unconverted-bases` is set, and constructs the
/// [`ZipperMergeConfig`] shared by the `ZipperZipStep` (pairing) and
/// `ZipperMerge` (per-template merge) steps that `add_zipper` wires up.
///
/// [`ZipperMergeConfig`]: merge_step::ZipperMergeConfig
pub(crate) fn build_zipper_merge_config(
    caps: ZipperMergeCaptures,
) -> Result<merge_step::ZipperMergeConfig> {
    let ZipperMergeCaptures {
        zipper_opts,
        output_header,
        reference_path,
        tuning,
        missing_count,
        records_emitted,
    } = caps;

    let ZipperMergeRules { tag_info, skip_tc_tags, reference } =
        zipper_opts.merge_rules(&reference_path)?;

    let cfg = merge_step::ZipperMergeConfig {
        tag_info,
        skip_tc_tags,
        exclude_missing_reads: zipper_opts.exclude_missing_reads,
        reference,
        output_header,
        missing_count,
        records_emitted,
        target_batch_count: tuning.template_batch_size,
        output_byte_limit: tuning.per_step_byte_limit,
    };
    Ok(cfg)
}
