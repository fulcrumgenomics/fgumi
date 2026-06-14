//! Chain builder for `Stage::Correct`.
//!
//! Phase 2 (T2.21) held the full ~400-LOC chain construction here.
//! Phase 3 (T3a.11) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the correct-specific types and step factories that the builder
//! imports: [`CorrectFinalizeHook`] used by `ChainBuilder::add_correct`.
//!
//! [`build_correct_chain`] shrinks to a ~10-line delegate.
//!
//! ## Four chain shapes
//!
//! Correct dispatches on `(input_source, track_rejects)`:
//!
//! - `(Bam,  false)`: `[read preamble] → GroupByQueryname → correct_step_kept_only → SerializeBamRecords`
//! - `(Bam,  true)`:  `[read preamble] → GroupByQueryname → correct_step_with_rejects` (branch 0 kept → `SerializeBamRecords`; branch 1 rejects bytes → `BgzfCompress`+Write)
//! - `(Sam,  false)`: `[sam preamble]  → GroupByQueryname → correct_step_kept_only  → SerializeBamRecords`
//! - `(Sam,  true)`:  `[sam preamble]  → GroupByQueryname → correct_step_with_rejects` (same as Bam/true)
//!
//! The source preamble differences (BAM vs SAM) are handled by `add_source`
//! in [`ChainBuilder`]; `add_correct` always receives a
//! `BamTemplateBatch`-typed tail (after `GroupByQueryname`).

use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use ahash::AHashMap;
use anyhow::Result;
use log::{info, warn};

use crate::commands::correct::{
    CollectedCorrectMetrics, CorrectUmis, EncodedUmiSet, merge_umi_counts,
};
use crate::logging::OperationTimer;
use crate::metrics::correct::UmiCorrectionMetrics;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::builder::ChainBuilder;
use crate::pipeline::chains::{BuiltPipeline, ChainSpec, FinalizeHook};

/// Post-pipeline finalize hook for correct. Reduces per-thread metrics,
/// writes the optional metrics TSV, logs summary + warn banners, enforces
/// the `--min-corrected` ratio gate, and calls `timer.log_completion`.
///
/// Mirrors the body of `CorrectUmis::finalize_correct_run`, extracted here so
/// `BuiltPipeline::run` can call it after `Pipeline::run` returns.
///
/// `pub(crate)` so [`ChainBuilder`] can construct and register it in
/// `add_correct`.
pub(crate) struct CorrectFinalizeHook {
    pub(crate) metrics: Arc<PerThreadAccumulator<CollectedCorrectMetrics>>,
    pub(crate) records_emitted: Arc<AtomicU64>,
    pub(crate) encoded_umi_set: Arc<EncodedUmiSet>,
    pub(crate) unmatched_umi: String,
    /// Path for the optional UMI-metrics TSV output.
    pub(crate) metrics_path: Option<PathBuf>,
    /// `--min-corrected` gate value (if set).
    pub(crate) min_corrected: Option<f64>,
    /// Used by `finalize_metrics` to compute and write the metrics TSV.
    pub(crate) correct: CorrectUmis,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for CorrectFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let CorrectFinalizeHook {
            metrics,
            records_emitted,
            encoded_umi_set,
            unmatched_umi,
            metrics_path: _metrics_path,
            min_corrected,
            correct,
            timer,
        } = *self;

        // Drain all PerThreadAccumulator slots into per-counter totals +
        // a single AHashMap that finalize_metrics consumes.
        let mut total_templates = 0u64;
        let mut total_missing = 0u64;
        let mut total_wrong_length = 0u64;
        let mut total_mismatched = 0u64;
        let mut merged_umi_matches: AHashMap<String, UmiCorrectionMetrics> = AHashMap::new();
        for slot in metrics.slots() {
            let mut m = slot.lock();
            total_templates += m.templates_processed;
            total_missing += m.missing_umis;
            total_wrong_length += m.wrong_length;
            total_mismatched += m.mismatched;
            for (umi, counts) in m.umi_matches.drain() {
                merge_umi_counts(&mut merged_umi_matches, umi, &counts);
            }
        }

        // Ensure ALL UMI rows present (matches fgbio behavior).
        for umi in encoded_umi_set.strings.iter().chain(std::iter::once(&unmatched_umi.clone())) {
            merged_umi_matches
                .entry(umi.clone())
                .or_insert_with(|| UmiCorrectionMetrics::new(umi.clone()));
        }
        correct.finalize_metrics(&mut merged_umi_matches, &unmatched_umi)?;

        // Log summary.
        let records_written = records_emitted.load(Ordering::Relaxed);
        let rejected = total_missing + total_wrong_length + total_mismatched;
        let total_records = records_written + rejected;
        info!("Read {total_records}; kept {records_written} and rejected {rejected}");
        info!("Total templates processed: {total_templates}");

        // Warn banner on missing/wrong-length UMIs.
        if total_missing > 0 || total_wrong_length > 0 {
            warn!("###################################################################");
            if total_missing > 0 {
                warn!("# {total_missing} were missing UMI attributes in the BAM file!");
            }
            if total_wrong_length > 0 {
                warn!(
                    "# {total_wrong_length} had unexpected UMIs of differing lengths in the BAM file!"
                );
            }
            warn!("###################################################################");
        }

        // Check minimum correction ratio.
        if let Some(min) = min_corrected {
            check_min_corrected(records_written, total_records, min)?;
        }

        timer.log_completion(records_written);

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// build_correct_chain — thin delegate
// ─────────────────────────────────────────────────────────────────────────────

/// Enforce the `--min-corrected` floor on the ratio of kept records.
///
/// Errors when `records_written / total_records` falls below `min`. When no
/// records were processed (`total_records == 0`) the ratio is undefined, so a
/// positive `min` errors (it cannot be met) while a zero `min` passes.
fn check_min_corrected(records_written: u64, total_records: u64, min: f64) -> Result<()> {
    // No reads processed: the kept ratio is undefined (0 / 0 = NaN), and
    // `NaN < min` is always false — which would silently bypass the gate. A
    // positive minimum cannot be satisfied with zero reads, so fail explicitly;
    // a zero minimum imposes no requirement and passes.
    if total_records == 0 {
        if min > 0.0 {
            anyhow::bail!(
                "No reads were processed, so the minimum ratio of reads kept (user specified minimum was {min:.2}) \
                could not be met. This could indicate empty input or a mismatch between library \
                preparation and the provided UMI file."
            );
        }
        return Ok(());
    }
    #[allow(clippy::cast_precision_loss)]
    let ratio_kept = records_written as f64 / total_records as f64;
    if ratio_kept < min {
        anyhow::bail!(
            "Final ratio of reads kept / total was {ratio_kept:.2} (user specified minimum was {min:.2}). \
            This could indicate a mismatch between library preparation and the provided UMI file."
        );
    }
    Ok(())
}

/// Build a [`BuiltPipeline`] for a correct-only chain.
///
/// Delegates to [`ChainBuilder`]. The full step sequence lives in
/// [`ChainBuilder`]'s `add_correct` method in
/// [`crate::pipeline::chains::builder`].
///
/// # Errors
///
/// Returns input-validation errors or any underlying pipeline construction error.
#[allow(clippy::needless_pass_by_value)]
pub fn build_correct_chain(spec: ChainSpec) -> Result<BuiltPipeline> {
    use crate::pipeline::chains::{Stage, builder::StagePosition};
    let mut chain = ChainBuilder::new(&spec)?;
    chain.add_source()?;
    chain.add_stage(Stage::Correct, StagePosition::Terminal)?;
    chain.add_sink()?;
    chain.build()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_records_with_positive_min_is_an_error() {
        // Regression: with no reads processed, `records_written / total_records`
        // is 0/0 = NaN and `NaN < min` is always false, which silently bypassed
        // the --min-corrected gate. A positive minimum cannot be satisfied with
        // zero reads, so this must fail explicitly.
        assert!(
            check_min_corrected(0, 0, 0.9).is_err(),
            "zero records with a positive minimum must error, not silently pass"
        );
    }

    #[test]
    fn zero_records_with_zero_min_is_ok() {
        // A zero minimum imposes no requirement, so an empty run passes.
        assert!(check_min_corrected(0, 0, 0.0).is_ok());
    }

    #[test]
    fn ratio_below_min_is_an_error() {
        // 5 / 100 = 0.05 < 0.90.
        assert!(check_min_corrected(5, 100, 0.9).is_err());
    }

    #[test]
    fn ratio_at_or_above_min_is_ok() {
        // 95 / 100 = 0.95 >= 0.90, and an exact match passes.
        assert!(check_min_corrected(95, 100, 0.9).is_ok());
        assert!(check_min_corrected(90, 100, 0.9).is_ok());
    }
}
