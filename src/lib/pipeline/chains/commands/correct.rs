//! Chain builder for `Stage::Correct`.
//!
//! Phase 2 (T2.21) held the full ~400-LOC chain construction here.
//! Phase 3 (T3a.11) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the correct-specific types and step factories that the builder
//! imports: `CorrectFinalizeHook` used by `ChainBuilder::add_correct`.
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
//! in `ChainBuilder`; `add_correct` always receives a
//! `BamTemplateBatch`-typed tail (after `GroupByQueryname`).

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use ahash::AHashMap;
use anyhow::Result;
use log::{error, info};

use crate::commands::correct::{
    CollectedCorrectMetrics, CorrectOptions, EncodedUmiSet, merge_umi_counts,
};
use crate::logging::OperationTimer;
use crate::metrics::correct::UmiCorrectionMetrics;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::FinalizeHook;

/// Post-pipeline finalize hook for correct. Reduces per-thread metrics,
/// writes the optional metrics TSV, logs the summary and the
/// missing/wrong-length-UMI error banner, enforces the `--min-corrected`
/// ratio gate, and calls `timer.log_completion`.
///
/// Mirrors the body of `CorrectUmis::finalize_correct_run`, extracted here so
/// `BuiltPipeline::run` can call it after `Pipeline::run` returns.
///
/// `pub(crate)` so `ChainBuilder` can construct and register it in
/// `add_correct`.
pub(crate) struct CorrectFinalizeHook {
    pub(crate) metrics: Arc<PerThreadAccumulator<CollectedCorrectMetrics>>,
    pub(crate) records_emitted: Arc<AtomicU64>,
    pub(crate) encoded_umi_set: Arc<EncodedUmiSet>,
    pub(crate) unmatched_umi: String,
    /// `--min-corrected` gate value (if set).
    pub(crate) min_corrected: Option<f64>,
    /// Used by `finalize_metrics` to compute and write the metrics TSV.
    pub(crate) correct: CorrectOptions,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for CorrectFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let CorrectFinalizeHook {
            metrics,
            records_emitted,
            encoded_umi_set,
            unmatched_umi,
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

        // fgbio logs this summary at error level (CorrectUmis.scala:275-280).
        if total_missing > 0 || total_wrong_length > 0 {
            error!("###################################################################");
            error!("# {total_missing} were missing UMI attributes in the BAM file!");
            error!(
                "# {total_wrong_length} had unexpected UMIs of differing lengths in the BAM file!"
            );
            error!("###################################################################");
        }

        // Check minimum correction ratio.
        if let Some(min) = min_corrected {
            check_min_corrected(records_written, total_records, min)?;
        }

        timer.log_completion(records_written);

        Ok(())
    }
}

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

    // Shares the crate-wide capturing logger (see
    // `crate::commands::common::test_log_capture`) so this test does not
    // install a competing process-global logger under plain `cargo t`.
    use crate::commands::common::test_log_capture::{capture_logs, captured_with_level};
    use crate::commands::correct::Target;
    use rstest::rstest;

    /// Content needles identifying the missing/wrong-length-UMI banner lines,
    /// deliberately excluding the `###`-only separator lines: those would
    /// also match an unrelated banner, so anchoring on the two message lines
    /// keeps a false match from masking a real regression.
    const MISSING_UMI_NEEDLE: &str = "were missing UMI attributes";
    const WRONG_LENGTH_NEEDLE: &str = "had unexpected UMIs of differing lengths";

    /// fgbio logs the missing/wrong-length-UMI summary banner at error level
    /// (`CorrectUmis.scala:275-280`), and the legacy single-threaded path
    /// (`src/lib/commands/correct.rs`) matches that. Pin the chain
    /// (`--threads`) path to the same level and presence:
    ///
    /// - `missing_and_wrong_length_present`: regression coverage for the
    ///   banner silently being downgraded to `warn!` here while the legacy
    ///   path stayed at `error!`, which let `correct --threads` users'
    ///   banners slip past error-level log filters that the legacy path's
    ///   banners would not.
    /// - `missing_only` / `wrong_length_only`: independent single-trigger
    ///   cases. Each on its own would flip to no-banner under a `||`->`&&`
    ///   guard regression (`2 > 0 && 0 > 0` is false), which the combined
    ///   `(2, 1)` case cannot catch since both counters are nonzero there.
    /// - `none_missing_or_wrong_length`: without this case, a mutant that
    ///   drops the `if total_missing > 0 || total_wrong_length > 0` guard and
    ///   emits the banner unconditionally would still pass the positive case
    ///   alone.
    ///
    /// Whenever a banner is expected, both message needles must be present
    /// (the guard emits both lines together), not merely one -- requiring only
    /// one would let a dropped message line slip through.
    #[rstest]
    #[case::missing_and_wrong_length_present(2, 1, true)]
    #[case::missing_only(2, 0, true)]
    #[case::wrong_length_only(0, 1, true)]
    #[case::none_missing_or_wrong_length(0, 0, false)]
    fn missing_or_wrong_length_umi_banner_level_and_presence(
        #[case] missing_umis: u64,
        #[case] wrong_length: u64,
        #[case] expect_banner: bool,
    ) {
        let _session = capture_logs();

        let metrics = PerThreadAccumulator::new(1);
        metrics.with_slot(|m: &mut CollectedCorrectMetrics| {
            m.templates_processed = 3;
            m.missing_umis = missing_umis;
            m.wrong_length = wrong_length;
        });

        let hook = CorrectFinalizeHook {
            metrics,
            records_emitted: Arc::new(AtomicU64::new(5)),
            encoded_umi_set: Arc::new(EncodedUmiSet::new(&["AAA".to_string(), "CCC".to_string()])),
            unmatched_umi: "NNN".to_string(),
            min_corrected: None,
            correct: CorrectOptions {
                metrics: None,
                target: Target::Umi,
                max_mismatches: 1,
                min_distance_diff: 1,
                umis: vec!["AAA".to_string(), "CCC".to_string()],
                umi_files: vec![],
                dont_store_original_umis: false,
                cache_size: 10,
                min_corrected: None,
                revcomp: false,
                rejects_path: None,
            },
            timer: OperationTimer::new("Correct"),
        };

        Box::new(hook).finalize().expect("finalize must succeed");

        let logs = captured_with_level();
        // The banner emits both message lines together whenever it fires, so
        // locate each independently and require both when a banner is expected
        // (requiring only one would let a dropped line pass).
        let missing_line = logs.iter().find(|(_, msg)| msg.contains(MISSING_UMI_NEEDLE));
        let wrong_length_line = logs.iter().find(|(_, msg)| msg.contains(WRONG_LENGTH_NEEDLE));

        if expect_banner {
            // `assert!` + `unwrap` rather than `unwrap_or_else(|| panic!(..))`: the
            // never-taken closure of the latter is an uncovered region on the
            // (always-passing) happy path, whereas an `assert!`'s panic branch
            // folds onto its covered line. Same rigor -- both lines must be present.
            assert!(
                missing_line.is_some(),
                "expected the missing-UMI banner line ({MISSING_UMI_NEEDLE:?}); got: {logs:?}"
            );
            assert!(
                wrong_length_line.is_some(),
                "expected the wrong-length-UMI banner line ({WRONG_LENGTH_NEEDLE:?}); got: {logs:?}"
            );
            for (level, msg) in [missing_line.unwrap(), wrong_length_line.unwrap()] {
                assert_eq!(
                    *level,
                    log::Level::Error,
                    "banner line {msg:?} logged at {level:?}; fgbio parity requires error level \
                    (CorrectUmis.scala:275-280)"
                );
            }
        } else {
            assert!(
                missing_line.is_none() && wrong_length_line.is_none(),
                "expected no missing/wrong-length-UMI banner when neither counter is nonzero \
                (the guard must suppress it); got: {logs:?}"
            );
        }
    }
}
