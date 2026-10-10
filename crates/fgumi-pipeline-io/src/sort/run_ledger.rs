//! Per-sort run facts that size each arena run so that its records **and**
//! its sort refs fit `--max-memory`.
//!
//! The split front's boundary step publishes the ref width when it provisions
//! the sort's key context, and its run sort publishes each run's exact byte and
//! record totals as it seals the run. `ReadBlocks` reads a snapshot when it
//! acquires the next run's arena and seals that run at
//! [`run_cap_for`]`(memory_limit, snapshot)` arena bytes.
//!
//! All stores happen-before the reads that size the next run: run k's arena is
//! acquired only after run k−1's chunk released it, which is after the run
//! sort sealed run k−1. So every run's cap is a deterministic function of the
//! previous run's totals — the same input splits at the same records at every
//! thread count.

use std::sync::Arc;
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};

use fgumi_sort::TemplateKeyVariant;

/// The mean record size the first run is charged at, before any run has been
/// sealed: `TEMPLATE_HEADER_SIZE + EST_BAM_BYTES_PER_TEMPLATE_RECORD`, the
/// estimator v0.7.0 used. It over-charges inputs with longer records (run 0
/// comes out a few percent smaller than steady state) and under-charges inputs
/// with shorter ones; later runs use the previous run's exact mean.
pub const RUN0_MEAN_RECORD_BYTES: u64 = 258;

/// Per-sort facts the run sort publishes at each seal and `ReadBlocks` reads
/// at the next run's first admit.
#[derive(Debug, Default)]
pub struct RunLedger {
    /// Bytes per ref; 0 = unknown (template `--key-types auto` before the
    /// first record provisions the lanes).
    ref_width: AtomicU32,
    /// Runs sealed so far.
    runs_sealed: AtomicU32,
    /// Arena bytes of the run just sealed.
    last_run_bytes: AtomicU64,
    /// Records of the run just sealed.
    last_run_records: AtomicU64,
}

/// What the next run is sized from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunPrior {
    /// Bytes per ref, once known.
    pub ref_width: Option<usize>,
    /// The previous run's mean record bytes, once a run with records sealed.
    pub mean_record_bytes: Option<u64>,
}

impl RunLedger {
    /// A ledger for one sort; `ref_width` is known up front for every order
    /// but template order under `--key-types auto`.
    #[must_use]
    pub fn new(ref_width: Option<usize>) -> Arc<Self> {
        let ledger = Self::default();
        if let Some(width) = ref_width {
            ledger.publish_width(width);
        }
        Arc::new(ledger)
    }

    /// The ref width, once the sort's key lanes are provisioned.
    ///
    /// # Panics
    ///
    /// Panics if `width` exceeds `u32::MAX` (refs are at most 56 bytes).
    pub fn publish_width(&self, width: usize) {
        let width = u32::try_from(width).expect("a ref width fits u32");
        self.ref_width.store(width, Ordering::Release);
    }

    /// A run sealed holding `bytes` arena bytes over `records` records.
    pub fn publish_seal(&self, bytes: u64, records: u64) {
        self.last_run_bytes.store(bytes, Ordering::Release);
        self.last_run_records.store(records, Ordering::Release);
        self.runs_sealed.fetch_add(1, Ordering::AcqRel);
    }

    /// Runs sealed so far.
    #[must_use]
    pub fn runs_sealed(&self) -> u32 {
        self.runs_sealed.load(Ordering::Acquire)
    }

    /// The facts the next run is sized from.
    #[must_use]
    pub fn snapshot(&self) -> RunPrior {
        let width = self.ref_width.load(Ordering::Acquire);
        let records = self.last_run_records.load(Ordering::Acquire);
        let bytes = self.last_run_bytes.load(Ordering::Acquire);
        RunPrior {
            ref_width: (width > 0).then_some(width as usize),
            mean_record_bytes: (records > 0).then(|| bytes / records),
        }
    }
}

/// Arena bytes a run may hold so that its records and its refs fit
/// `memory_limit`: `memory_limit × m / (m + w)`, with `m` the prior's mean
/// record bytes ([`RUN0_MEAN_RECORD_BYTES`] before any run sealed) and `w` the
/// ref width (the widest template ref while template lanes are unknown). A
/// zero width charges nothing.
#[must_use]
pub fn run_cap_for(memory_limit: usize, prior: RunPrior) -> usize {
    let m = u128::from(prior.mean_record_bytes.unwrap_or(RUN0_MEAN_RECORD_BYTES).max(1));
    let w = prior.ref_width.unwrap_or_else(|| TemplateKeyVariant::WIDEST.ref_width()) as u128;
    let cap = memory_limit as u128 * m / (m + w);
    // The cap never exceeds `memory_limit`, so it fits.
    usize::try_from(cap).unwrap_or(memory_limit)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    /// Expected caps are written as literals, computed by hand.
    #[rstest]
    // ⌊1e6 × 449 / 489⌋
    #[case::steady_state(Some(40), Some(449), 918_200)]
    #[case::no_refs(Some(0), Some(449), 1_000_000)]
    // ⌊1e6 × 258 / 298⌋
    #[case::run0_configured_width(Some(40), None, 865_771)]
    // ⌊1e6 × 258 / 314⌋ — WIDEST = 56
    #[case::template_auto_run0(None, None, 821_656)]
    // ⌊1e6 × 449 / 505⌋
    #[case::template_auto_later_run(None, Some(449), 889_108)]
    // ⌊1e6 × 300 / 324⌋
    #[case::coordinate(Some(24), Some(300), 925_925)]
    fn run_cap_for_charges_the_refs(
        #[case] ref_width: Option<usize>,
        #[case] mean_record_bytes: Option<u64>,
        #[case] expected: usize,
    ) {
        let prior = RunPrior { ref_width, mean_record_bytes };
        assert_eq!(run_cap_for(1_000_000, prior), expected);
    }

    #[test]
    fn widest_template_ref_is_56_bytes() {
        assert_eq!(TemplateKeyVariant::WIDEST.ref_width(), 56);
    }

    /// The ledger reports the configured width, then a published one; no mean
    /// until a run with records sealed, then the last run's exact mean.
    #[test]
    fn ledger_reports_the_last_run() {
        let ledger = RunLedger::new(None);
        assert_eq!(ledger.snapshot(), RunPrior { ref_width: None, mean_record_bytes: None });
        ledger.publish_width(48);
        assert_eq!(ledger.snapshot().ref_width, Some(48));
        ledger.publish_seal(0, 0);
        assert_eq!(ledger.snapshot().mean_record_bytes, None, "an empty run has no mean");
        ledger.publish_seal(1_000, 3);
        assert_eq!(ledger.snapshot().mean_record_bytes, Some(333));
        ledger.publish_seal(4_500, 10);
        assert_eq!(
            ledger.snapshot().mean_record_bytes,
            Some(450),
            "the last run, not a running mean"
        );
        assert_eq!(ledger.runs_sealed(), 3);
        let fixed = RunLedger::new(Some(24));
        assert_eq!(fixed.snapshot().ref_width, Some(24));
    }
}
