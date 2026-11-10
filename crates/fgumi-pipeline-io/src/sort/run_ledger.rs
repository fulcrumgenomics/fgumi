//! Per-sort run facts that size each arena run so that its records **and**
//! its sort refs fit `--max-memory`.
//!
//! The split front's boundary step publishes the ref width and the spill key
//! kind when it provisions the sort's key context, and its run sort publishes
//! each run's exact byte and record totals as it seals the run (adding them to
//! the sort's spilled totals). `ReadBlocks` reads a snapshot when it acquires
//! the next run's arena and seals that run at
//! [`run_cap_for`]`(memory_limit, snapshot)` arena bytes. A run index is not
//! charged here: closed runs' indexes live in their spill files' footers, and
//! the merge's budget pays for them.
//!
//! All stores happen-before the reads that size the next run: run k's arena is
//! acquired only after run k−1's chunk released it, which is after the run
//! sort sealed run k−1. So every run's cap is a deterministic function of the
//! previous run's totals — the same input splits at the same records at every
//! thread count.

use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};
use std::sync::{Arc, OnceLock};

use fgumi_sort::{SpillKeyKind, TemplateKeyVariant};

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
    /// Bytes of every run sealed so far (the sum of `last_run_bytes`).
    spilled_bytes: AtomicU64,
    /// Records of every run sealed so far (the sum of `last_run_records`).
    spilled_records: AtomicU64,
    /// The largest record (`u32` length prefix included) of every run sealed
    /// so far; 0 = none known.
    spilled_max_record: AtomicU64,
    /// The key kind runs spill with, once published.
    spill_key_kind: OnceLock<SpillKeyKind>,
}

/// What the next run is sized from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunPrior {
    /// Bytes per ref, once known.
    pub ref_width: Option<usize>,
    /// The previous run's mean record bytes, once a run with records sealed.
    pub mean_record_bytes: Option<u64>,
    /// The previous run's record bytes (0 before any run sealed).
    pub last_run_bytes: u64,
    /// The previous run's records (0 before any run sealed).
    pub last_run_records: u64,
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

    /// A run sealed holding `bytes` arena bytes over `records` records. The
    /// totals are added to the spilled totals before `runs_sealed` is bumped,
    /// so a reader that observes the seal count also observes them.
    pub fn publish_seal(&self, bytes: u64, records: u64) {
        self.last_run_bytes.store(bytes, Ordering::Release);
        self.last_run_records.store(records, Ordering::Release);
        self.spilled_bytes.fetch_add(bytes, Ordering::AcqRel);
        self.spilled_records.fetch_add(records, Ordering::AcqRel);
        self.runs_sealed.fetch_add(1, Ordering::AcqRel);
    }

    /// The largest record (its `u32` length prefix included) of a run about
    /// to be sealed, published before [`Self::publish_seal`] (so before the
    /// run's chunk can release its arena). Keeps the maximum over runs.
    pub fn publish_run_max_record(&self, len: u64) {
        self.spilled_max_record.fetch_max(len, Ordering::AcqRel);
    }

    /// The largest record of every run sealed so far (length prefix
    /// included); `None` when none was published.
    #[must_use]
    pub fn spilled_max_record(&self) -> Option<u64> {
        Some(self.spilled_max_record.load(Ordering::Acquire)).filter(|&m| m > 0)
    }

    /// The kind every run of the sort spills with, published once when the
    /// boundary step provisions the key context (before any run can seal).
    /// Later publishes are ignored.
    pub fn publish_spill_key_kind(&self, kind: SpillKeyKind) {
        let _ = self.spill_key_kind.set(kind);
    }

    /// The spill key kind; `None` until the boundary step published it.
    #[must_use]
    pub fn spill_key_kind(&self) -> Option<SpillKeyKind> {
        self.spill_key_kind.get().copied()
    }

    /// `(bytes, records)` of every run sealed so far.
    #[must_use]
    pub fn spilled_totals(&self) -> (u64, u64) {
        (self.spilled_bytes.load(Ordering::Acquire), self.spilled_records.load(Ordering::Acquire))
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
            last_run_bytes: bytes,
            last_run_records: records,
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
        let prior =
            RunPrior { ref_width, mean_record_bytes, last_run_bytes: 0, last_run_records: 0 };
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
        assert_eq!(
            ledger.snapshot(),
            RunPrior {
                ref_width: None,
                mean_record_bytes: None,
                last_run_bytes: 0,
                last_run_records: 0
            }
        );
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

    /// Every seal adds to the spilled totals; the last run is still the
    /// last seal's alone; the largest record is the maximum over runs.
    #[test]
    fn ledger_accumulates_seals() {
        let ledger = RunLedger::new(Some(24));
        assert_eq!(ledger.spilled_totals(), (0, 0));
        assert_eq!(ledger.spilled_max_record(), None);
        ledger.publish_run_max_record(400);
        ledger.publish_run_max_record(300);
        assert_eq!(ledger.spilled_max_record(), Some(400));
        ledger.publish_seal(900_000, 3_000);
        ledger.publish_seal(800_000, 2_000);
        assert_eq!(ledger.spilled_totals(), (1_700_000, 5_000));
        assert_eq!(ledger.runs_sealed(), 2);
        assert_eq!(ledger.snapshot().mean_record_bytes, Some(400));
        assert_eq!(
            (ledger.snapshot().last_run_bytes, ledger.snapshot().last_run_records),
            (800_000, 2_000)
        );
    }

    /// The spill key kind is unknown until published, and the first publish
    /// wins.
    #[test]
    fn spill_key_kind_is_published_once() {
        let ledger = RunLedger::new(None);
        assert_eq!(ledger.spill_key_kind(), None);
        ledger.publish_spill_key_kind(SpillKeyKind::TemplateK40);
        assert_eq!(ledger.spill_key_kind(), Some(SpillKeyKind::TemplateK40));
        ledger.publish_spill_key_kind(SpillKeyKind::Coordinate);
        assert_eq!(ledger.spill_key_kind(), Some(SpillKeyKind::TemplateK40));
    }

    /// Every kind round-trips through the cell.
    #[rstest]
    #[case::coordinate(SpillKeyKind::Coordinate)]
    #[case::qn_lex(SpillKeyKind::QuerynameLex)]
    #[case::qn_natural(SpillKeyKind::QuerynameNatural)]
    #[case::k24(SpillKeyKind::TemplateK24)]
    #[case::cb32(SpillKeyKind::TemplateCb32)]
    #[case::tert32(SpillKeyKind::TemplateTert32)]
    #[case::k40(SpillKeyKind::TemplateK40)]
    fn spill_key_kind_round_trips(#[case] kind: SpillKeyKind) {
        let ledger = RunLedger::new(None);
        ledger.publish_spill_key_kind(kind);
        assert_eq!(ledger.spill_key_kind(), Some(kind));
    }
}
