//! Chain builder for `Stage::Correct`.
//!
//! Phase 2 (T2.21) held the full ~400-LOC chain construction here.
//! Phase 3 (T3a.11) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the correct-specific types and step factories that the builder
//! imports: `CorrectFinalizeHook` (always-run: summary/banner) and
//! `CorrectMetricsFinalizeHook` (success-only: `--metrics` TSV write, THEN the
//! `--min-corrected` ratio gate), both used by `ChainBuilder::add_correct`.
//!
//! ## fgbio parity: metrics vs. the `--min-corrected` gate
//!
//! fgbio's `CorrectUmis` writes the metrics file once processing completes,
//! *before* it checks the `--min-corrected` ratio — so a ratio-below-floor
//! failure is not treated as a processing failure and still leaves the
//! metrics (and other output) files in place. Only a genuine processing error
//! (I/O, malformed input) withholds them. `CorrectMetricsFinalizeHook`
//! reproduces that ordering directly: it writes the `--metrics` TSV first,
//! then (only when `--min-corrected` was given) runs the ratio gate — both
//! within the same success-only hook, so the write always happens before the
//! gate can error. A genuine pipeline processing failure skips the whole
//! `finalize_on_success` list, withholding both.
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
    CollectedCorrectMetrics, CorrectOptions, EncodedUmiSet, check_min_corrected, merge_umi_counts,
};
use crate::logging::OperationTimer;
use crate::metrics::correct::UmiCorrectionMetrics;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::FinalizeHook;

/// Aggregated per-thread totals shared by [`CorrectFinalizeHook`]'s summary
/// and [`CorrectMetricsFinalizeHook`]'s `--min-corrected` gate. Extracted into
/// one struct + reducer ([`sum_correct_counts`]) so the two call sites read
/// the exact same counters and cannot silently diverge if the set of
/// rejection reasons on [`CollectedCorrectMetrics`] ever changes.
#[derive(Default)]
pub(crate) struct CorrectCounts {
    pub(crate) templates: u64,
    pub(crate) missing: u64,
    pub(crate) wrong_length: u64,
    pub(crate) mismatched: u64,
}

impl CorrectCounts {
    /// Total rejected records: every rejection-reason counter summed.
    pub(crate) fn rejected(&self) -> u64 {
        self.missing + self.wrong_length + self.mismatched
    }
}

/// Sum the plain per-thread counters (everything except `umi_matches`) across
/// every `PerThreadAccumulator` slot. Reads only `Copy` fields, so this is
/// non-destructive and safe to call from more than one finalize hook sharing
/// the same accumulator.
pub(crate) fn sum_correct_counts(
    metrics: &PerThreadAccumulator<CollectedCorrectMetrics>,
) -> CorrectCounts {
    let mut counts = CorrectCounts::default();
    for slot in metrics.slots() {
        let m = slot.lock();
        counts.templates += m.templates_processed;
        counts.missing += m.missing_umis;
        counts.wrong_length += m.wrong_length;
        counts.mismatched += m.mismatched;
    }
    counts
}

/// Post-pipeline finalize hook for correct. Reduces per-thread counters (via
/// [`sum_correct_counts`]), logs the summary and the missing/wrong-length-UMI
/// error banner, and calls `timer.log_completion`.
///
/// Registered on the **always-run** `finalize` list, so the summary/banner are
/// reported even after a pipeline processing failure. The `--metrics` TSV
/// write and the `--min-corrected` ratio gate live on the success-only
/// [`CorrectMetricsFinalizeHook`] instead — see the module-level fgbio-parity
/// note for why the gate is success-gated rather than living here. Shares the
/// `metrics` accumulator with that hook via `Arc`; only the per-counter
/// totals (not `umi_matches`) are read here, non-destructively, so both hooks
/// can independently read the same accumulator.
///
/// `pub(crate)` so `ChainBuilder` can construct and register it in
/// `add_correct`.
pub(crate) struct CorrectFinalizeHook {
    pub(crate) metrics: Arc<PerThreadAccumulator<CollectedCorrectMetrics>>,
    pub(crate) records_emitted: Arc<AtomicU64>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for CorrectFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let CorrectFinalizeHook { metrics, records_emitted, timer } = *self;

        let counts = sum_correct_counts(&metrics);

        // Log summary.
        let records_written = records_emitted.load(Ordering::Relaxed);
        let rejected = counts.rejected();
        let total_records = records_written + rejected;
        info!("Read {total_records}; kept {records_written} and rejected {rejected}");
        info!("Total templates processed: {}", counts.templates);

        // fgbio logs this summary at error level (CorrectUmis.scala:275-280).
        if counts.missing > 0 || counts.wrong_length > 0 {
            error!("###################################################################");
            error!("# {} were missing UMI attributes in the BAM file!", counts.missing);
            error!(
                "# {} had unexpected UMIs of differing lengths in the BAM file!",
                counts.wrong_length
            );
            error!("###################################################################");
        }

        timer.log_completion(records_written);

        Ok(())
    }
}

/// Success-only finalize hook that writes the `--metrics` TSV and then, only
/// when `--min-corrected` was given, enforces the ratio gate. Registered on
/// `finalize_on_success` (not the always-run `finalize`), so a genuine
/// pipeline processing failure withholds both — but, per fgbio parity, a
/// `--min-corrected` ratio failure does NOT withhold the metrics write: the
/// TSV is written first, inside this same hook, before the gate can error
/// (see the module-level fgbio-parity note). Owns the sole read of the shared
/// `metrics` accumulator's `umi_matches` map, so it drains it (moves out);
/// the ratio gate instead calls [`sum_correct_counts`], the same reducer
/// [`CorrectFinalizeHook`] uses for its summary.
pub(crate) struct CorrectMetricsFinalizeHook {
    pub(crate) metrics: Arc<PerThreadAccumulator<CollectedCorrectMetrics>>,
    pub(crate) records_emitted: Arc<AtomicU64>,
    pub(crate) encoded_umi_set: Arc<EncodedUmiSet>,
    pub(crate) unmatched_umi: String,
    /// Used by `finalize_metrics` to compute and write the metrics TSV, and
    /// to read the `--min-corrected` gate value (if any).
    pub(crate) correct: CorrectOptions,
}

impl FinalizeHook for CorrectMetricsFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let CorrectMetricsFinalizeHook {
            metrics,
            records_emitted,
            encoded_umi_set,
            unmatched_umi,
            correct,
        } = *self;

        // Drain every PerThreadAccumulator slot's `umi_matches` into a single
        // AHashMap that `finalize_metrics` consumes. This hook is the sole
        // consumer of `umi_matches` (`CorrectFinalizeHook` only reads the
        // plain counters via `sum_correct_counts`), so draining (move, not
        // clone) is safe here.
        let mut merged_umi_matches: AHashMap<String, UmiCorrectionMetrics> = AHashMap::new();
        for slot in metrics.slots() {
            let mut m = slot.lock();
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
        // Write the metrics TSV BEFORE the ratio gate below, so a
        // --min-corrected failure still leaves it on disk (fgbio parity).
        correct.finalize_metrics(&mut merged_umi_matches, &unmatched_umi)?;

        // Enforce the --min-corrected ratio gate, if requested. Uses the same
        // `sum_correct_counts` reducer as CorrectFinalizeHook's summary, so
        // the two can't diverge on which counters count as "rejected".
        if let Some(min) = correct.min_corrected {
            let counts = sum_correct_counts(&metrics);
            let records_written = records_emitted.load(Ordering::Relaxed);
            let total_records = records_written + counts.rejected();
            check_min_corrected(records_written, total_records, min)?;
        }

        Ok(())
    }
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
