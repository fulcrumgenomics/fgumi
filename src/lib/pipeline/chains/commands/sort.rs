//! Chain builder support for `Stage::Sort`.
//!
//! The stage-by-stage construction lives in
//! [`crate::pipeline::chains::builder::ChainBuilder`]'s `add_sort` method.
//! This module provides the standalone-sort summary finalize hook
//! (`SortSummaryFinalizeHook`), registered by `add_sort` for a sole-
//! `[Stage::Sort]` chain. A `SinkSpec::BamWithIndex` request needs no
//! finalize hook here: `ChainBuilder::add_sink` attaches an inline BAI
//! indexer directly to the `WriteBgzfFile` sink
//! (`WriteBgzfFile::with_bai_index`), which builds the `.bai` from the
//! `BamIndexManifest`s each `BgzfCompress`-produced block carries, as the sink
//! drains — no post-pipeline re-read of the finished BAM.
//!
//! ## Sort pipeline topology
//!
//! A sole-`[Stage::Sort]` chain runs through the same streaming source → arena
//! sort ingest → `SpillGather` → `SpillBlockCompress` → `SpillWrite` →
//! `SortSpillDecompress` → `SortMerge` → sink pipeline as a
//! fused sort stage, via the normal `add_source` / `add_sink` flow that
//! [`crate::pipeline::chains::build::build_for`] drives for every stage; there
//! is no longer a self-contained file→file sort step or a sort-only chain
//! builder.

use std::path::PathBuf;
use std::sync::Arc;

use crate::logging::OperationTimer;
use anyhow::Result;
use log::info;
use parking_lot::Mutex;

use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::core::runtime::stats::{PipelineStats, StatsSnapshot};

/// Post-pipeline summary for standalone `fgumi sort`.
///
/// Reads the `SortMerge` stats slot (records processed/written + spill-chunk
/// count) and logs the `=== Summary ===` block, then the timer's
/// records-per-second completion line. Registered by
/// `ChainBuilder::add_sort` only for a sole-`[Stage::Sort]` chain; the fused
/// `runall` path leaves the slot unset and gets no summary block.
pub(crate) struct SortSummaryFinalizeHook {
    pub(crate) stats_slot: Arc<Mutex<Option<fgumi_sort::SortStats>>>,
    pub(crate) output_path: PathBuf,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for SortSummaryFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let SortSummaryFinalizeHook { stats_slot, output_path, timer } = *self;
        let stats = stats_slot.lock().take().unwrap_or_default();
        info!("=== Summary ===");
        info!("Records processed: {}", stats.total_records);
        info!("Records written: {}", stats.output_records);
        if stats.runs_written > 0 {
            info!("Spill runs: {}", stats.runs_written);
        }
        info!("Output: {}", output_path.display());
        timer.log_completion(stats.total_records);
        Ok(())
    }
}

// ─────────────────────────── Sort phase timing ───────────────────────────
//
// Pre-cutover, the owned `RawExternalSorter` engine emitted a
// `=== Sort Phase Timing ===` per-phase wall-time breakdown at `info!` on every
// sort. Standalone `fgumi sort` now runs through the declarative chain, whose
// sort phases are modeled as discrete pipeline steps
// (`ReadBlocks`/`InflateToArena` → `FindBoundariesAndSort`/`SortBuffer` →
// `SpillGather`/`SpillBlockCompress`/`SpillWrite` → `SortSpillDecompress` →
// `SortMerge` → `BgzfCompress`/`WriteBgzfFile`). The pipeline already records
// each step's cumulative `try_run` busy time
// (`StepStatsSnapshot::total_run_ns`) whenever a `PipelineStats` collector is
// attached, so this hook re-derives the same per-phase breakdown from that
// snapshot rather than re-running the retired engine. That restores the
// `=== Sort Phase Timing ===` grep string and the project's documented
// sort-profiling signal (see `CLAUDE.md`, "Benchmarking Notes"). It is
// registered only when stats collection is on (`--pipeline-stats` where the
// command exposes it, else `FGUMI_PIPELINE_STATS=1` — standalone `fgumi sort`
// uses the env var), so the default sort path keeps its zero-overhead behavior
// — which the external `fg-labs/mako` consumer relies on.

/// The six sort phases, in pipeline order, that `=== Sort Phase Timing ===`
/// reports. The declaration order is also the accumulator index order used by
/// [`summarize_sort_phases`] (`phase as usize`) and the report order in
/// [`SortPhase::ORDER`] — keep the three in sync.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SortPhase {
    ReadDecompress,
    InMemorySort,
    SpillWrite,
    Consolidation,
    KWayMerge,
    WriteOutput,
}

impl SortPhase {
    /// Phases in report order. The array index equals each variant's
    /// `as usize`, i.e. its slot in the per-phase nanosecond accumulator.
    const ORDER: [SortPhase; 6] = [
        SortPhase::ReadDecompress,
        SortPhase::InMemorySort,
        SortPhase::SpillWrite,
        SortPhase::Consolidation,
        SortPhase::KWayMerge,
        SortPhase::WriteOutput,
    ];

    /// Human-readable phase label used in the log block. Mirrors the phase
    /// names the retired owned-engine `SortPhaseTimer` printed.
    fn label(self) -> &'static str {
        match self {
            SortPhase::ReadDecompress => "read + decompress",
            SortPhase::InMemorySort => "in-memory sort",
            SortPhase::SpillWrite => "spill write",
            SortPhase::Consolidation => "consolidation",
            SortPhase::KWayMerge => "k-way merge",
            SortPhase::WriteOutput => "write output",
        }
    }

    /// Classify a pipeline step name into its sort phase, or `None` for a step
    /// that is not part of the sort engine (record/format adapters such as a SAM
    /// source's parse steps, whose time is therefore not counted here).
    ///
    /// The arms are the `StepProfile::name` literals of the sort steps in
    /// `crates/fgumi-pipeline-io/src/sort/` plus the shared BGZF sink, covering
    /// the steps a sole-`[Stage::Sort]` chain runs (the only chain the timing
    /// hook is registered on — see `ChainBuilder::build`). A sort step added or
    /// renamed there without a matching arm here falls to `None` and is silently
    /// omitted from the roll-up (and its `total`) — so when adding a sort step,
    /// add its name here too.
    ///
    /// `sort_phase_classifier_covers_known_sort_steps` pins the arm→phase mapping
    /// for this known set (a maintainer checklist), but note what it does NOT do:
    /// because it asserts against the same string literals, it cannot by itself
    /// detect an *upstream* rename or a brand-new step (those are caught at the
    /// io layer by each step's own `assert_eq!(profile.name, …)` test, and here
    /// only when the new name is added to both this match and that case table).
    fn from_step_name(name: &str) -> Option<SortPhase> {
        match name {
            "ReadBlocks" | "InflateToArena" => Some(SortPhase::ReadDecompress),
            "FindBoundariesAndSort" | "SortBuffer" => Some(SortPhase::InMemorySort),
            // `CompressSpill` is the fused compress-and-write spill variant
            // (`SortBuffer → CompressSpill → …`); `add_sort` currently wires the
            // `SpillGather`/`SpillBlockCompress`/`SpillWrite` split instead, but
            // classify it so the roll-up stays correct if that wiring changes.
            "SpillGather" | "SpillBlockCompress" | "SpillWrite" | "CompressSpill" => {
                Some(SortPhase::SpillWrite)
            }
            "SortSpillDecompress" => Some(SortPhase::Consolidation),
            "SortMerge" => Some(SortPhase::KWayMerge),
            "BgzfCompress" | "WriteBgzfFile" => Some(SortPhase::WriteOutput),
            _ => None,
        }
    }
}

/// Sum each sort step's cumulative `try_run` busy time into its phase bucket,
/// indexed by each variant's `as usize` (see [`SortPhase::ORDER`]).
///
/// Reads only `snapshot.steps` (`StepStatsSnapshot::total_run_ns`), which the
/// driver populates for **every** dispatched step — the detached writer
/// included, since `PipelineStats::record` runs on the detached thread's
/// dispatches too (`runtime/driver.rs`). The separate `snapshot.detached` view
/// re-counts that same time under the "N + 2" pool-vs-detached reporting split,
/// so it is deliberately NOT added here (doing so would double-count the
/// writer).
fn summarize_sort_phases(snapshot: &StatsSnapshot) -> [u64; 6] {
    let mut phase_ns = [0u64; 6];
    for (name, stats) in &snapshot.steps {
        if let Some(phase) = SortPhase::from_step_name(name) {
            phase_ns[phase as usize] = phase_ns[phase as usize].saturating_add(stats.total_run_ns);
        }
    }
    phase_ns
}

/// Render the `=== Sort Phase Timing ===` block from per-phase nanoseconds, or
/// `None` when no sort phase did any work (empty or pre-work-failure stats), so
/// a run that recorded nothing prints no misleading block.
// Nanosecond counts → seconds/percentages for a human-readable log line only;
// `f64` precision loss above 2^52 ns (~52 days per phase) is irrelevant to a
// displayed figure. Mirrors `runtime::stats`'s own display-math allows.
#[allow(clippy::cast_precision_loss)]
fn format_sort_phase_timing(phase_ns: &[u64; 6]) -> Option<Vec<String>> {
    let total_ns: u64 = phase_ns.iter().copied().sum();
    if total_ns == 0 {
        return None;
    }
    let total_secs = total_ns as f64 / 1e9;
    let mut lines = Vec::with_capacity(SortPhase::ORDER.len() + 3);
    lines.push("=== Sort Phase Timing ===".to_owned());
    for phase in SortPhase::ORDER {
        let ns = phase_ns[phase as usize];
        let secs = ns as f64 / 1e9;
        let pct = 100.0 * ns as f64 / total_ns as f64;
        lines.push(format!("  {:<18} {secs:>9.3}s  {pct:>5.1}%", phase.label()));
    }
    lines.push(format!("  {:<18} {total_secs:>9.3}s  100.0%", "total"));
    lines.push(
        "  (cumulative per-step busy time; sort steps run concurrently, so each \
         percentage is a share of total sort work, not of wall-clock time)"
            .to_owned(),
    );
    Some(lines)
}

/// Log the `=== Sort Phase Timing ===` block for a stats snapshot, when any
/// sort phase did work. Split out from [`SortPhaseTimingFinalizeHook`] so the
/// roll-up + rendering can be unit-tested against a hand-built snapshot.
fn log_sort_phase_timing(snapshot: &StatsSnapshot) {
    if let Some(lines) = format_sort_phase_timing(&summarize_sort_phases(snapshot)) {
        for line in lines {
            info!("{line}");
        }
    }
}

/// Post-pipeline `=== Sort Phase Timing ===` diagnostic for a chain that
/// contains a sort stage.
///
/// Registered by `ChainBuilder::build` only when a `PipelineStats` collector is
/// attached AND stats output was requested (`--pipeline-stats` where exposed,
/// else `FGUMI_PIPELINE_STATS=1`), so the default sort path pays nothing. Emits the
/// per-phase breakdown re-derived from the end-of-run stats snapshot.
pub(crate) struct SortPhaseTimingFinalizeHook {
    pub(crate) stats: Arc<PipelineStats>,
}

impl FinalizeHook for SortPhaseTimingFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        log_sort_phase_timing(&self.stats.snapshot());
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    // Shares the crate-wide capturing logger (see
    // `crate::commands::common::test_log_capture`) so this test and the
    // memory-budget capture tests in `commands::common` do not each install a
    // competing process-global logger — the second install panics under plain
    // `cargo t`, which runs every test in one process.
    use crate::commands::common::test_log_capture::{capture_logs, captured};

    /// `SortSummaryFinalizeHook::finalize` must log the "Spill runs:" wording
    /// the owned `execute_sort` engine uses (`src/lib/commands/sort.rs`), not
    /// the chain's former "Temporary runs:" wording -- `test_streaming_output`
    /// asserts the former through the sort command once the chain owns the
    /// summary.
    #[test]
    fn sort_summary_uses_spill_runs_wording() {
        let _session = capture_logs();

        let hook = SortSummaryFinalizeHook {
            stats_slot: Arc::new(Mutex::new(Some(fgumi_sort::SortStats {
                runs_written: 3,
                ..Default::default()
            }))),
            output_path: PathBuf::from("out.bam"),
            timer: OperationTimer::new("Sort"),
        };
        Box::new(hook).finalize().expect("finalize must succeed");

        let logs = captured();
        assert!(
            logs.iter().any(|line| line.contains("Spill runs: 3")),
            "expected a 'Spill runs: 3' log line; got: {logs:?}"
        );
        assert!(
            !logs.iter().any(|line| line.contains("Temporary runs:")),
            "must not emit the old 'Temporary runs:' wording; got: {logs:?}"
        );
    }

    // ── Sort phase timing (`=== Sort Phase Timing ===`) ──────────────────────

    use crate::pipeline::core::runtime::stats::StepStatsSnapshot;

    /// A `StepStatsSnapshot` carrying only the `total_run_ns` the phase roll-up
    /// reads; every other counter is zero.
    fn step_snap(total_run_ns: u64) -> StepStatsSnapshot {
        StepStatsSnapshot {
            try_run_total: 0,
            progress_count: 0,
            no_progress_count: 0,
            contention_count: 0,
            finished_count: 0,
            error_count: 0,
            total_run_ns,
            first_progress_ns: u64::MAX,
            last_progress_ns: 0,
        }
    }

    /// A snapshot whose only populated field is `steps` (the field the roll-up
    /// reads); `workers`/`detached`/`edges` empty.
    fn snapshot_with_steps(steps: Vec<(&'static str, u64)>) -> StatsSnapshot {
        StatsSnapshot {
            steps: steps.into_iter().map(|(n, ns)| (n, step_snap(ns))).collect(),
            workers: Vec::new(),
            detached: Vec::new(),
            edges: Vec::new(),
        }
    }

    #[test]
    fn summarize_sort_phases_buckets_every_step_into_its_phase() {
        let snap = snapshot_with_steps(vec![
            ("ReadBlocks", 100),
            ("InflateToArena", 200),        // read + decompress = 300
            ("FindBoundariesAndSort", 400), // in-memory sort = 400
            ("SpillGather", 10),
            ("SpillBlockCompress", 20),
            ("SpillWrite", 30),          // spill write = 60
            ("SortSpillDecompress", 50), // consolidation = 50
            ("SortMerge", 500),          // k-way merge = 500
            ("BgzfCompress", 5),
            ("WriteBgzfFile", 15), // write output = 20
        ]);
        let phase_ns = summarize_sort_phases(&snap);
        assert_eq!(phase_ns[SortPhase::ReadDecompress as usize], 300);
        assert_eq!(phase_ns[SortPhase::InMemorySort as usize], 400);
        assert_eq!(phase_ns[SortPhase::SpillWrite as usize], 60);
        assert_eq!(phase_ns[SortPhase::Consolidation as usize], 50);
        assert_eq!(phase_ns[SortPhase::KWayMerge as usize], 500);
        assert_eq!(phase_ns[SortPhase::WriteOutput as usize], 20);
    }

    /// The `SortBuffer` ingest step (the SAM/RecordBatch arena front) buckets
    /// into in-memory sort, exercising the branch the BAM path's
    /// `FindBoundariesAndSort` does not.
    #[test]
    fn summarize_maps_sort_buffer_ingest_to_in_memory_sort() {
        let phase_ns = summarize_sort_phases(&snapshot_with_steps(vec![("SortBuffer", 42)]));
        assert_eq!(phase_ns[SortPhase::InMemorySort as usize], 42);
        assert_eq!(phase_ns.iter().sum::<u64>(), 42);
    }

    #[test]
    fn summarize_ignores_non_sort_steps() {
        // Adapter / other-stage step names (e.g. a fused sort→group chain) must
        // not land in any sort phase bucket.
        let phase_ns = summarize_sort_phases(&snapshot_with_steps(vec![
            ("SortMerge", 100),
            ("DecodeFromRecords", 999),
            ("GroupBam", 999),
            ("TemplatesToRecordBatch", 999),
        ]));
        assert_eq!(phase_ns[SortPhase::KWayMerge as usize], 100);
        assert_eq!(
            phase_ns.iter().sum::<u64>(),
            100,
            "only the SortMerge step is part of the sort engine here"
        );
    }

    #[test]
    fn summarize_reads_steps_total_run_ns_not_detached_busy() {
        // The detached writer's time is recorded in BOTH steps[].total_run_ns
        // and detached[].busy_ns; the roll-up must read steps only, or it would
        // double-count the writer.
        let snap = StatsSnapshot {
            steps: vec![("WriteBgzfFile", step_snap(100))],
            workers: Vec::new(),
            detached: vec![(0, "WriteBgzfFile", 999_999, 0, 0)],
            edges: Vec::new(),
        };
        let phase_ns = summarize_sort_phases(&snap);
        assert_eq!(
            phase_ns[SortPhase::WriteOutput as usize],
            100,
            "must count steps[].total_run_ns, not detached[].busy_ns"
        );
    }

    /// Pins the arm→phase mapping for every sort step name the classifier
    /// recognizes, as a maintainer checklist of the covered set. It asserts the
    /// classifier's own literals, so it does NOT by itself catch an upstream
    /// rename or a newly-added step (add a case here when adding a `from_step_name`
    /// arm); it does catch an accidental change to which phase an existing name
    /// maps to.
    #[rstest]
    #[case::read_blocks("ReadBlocks", SortPhase::ReadDecompress)]
    #[case::inflate("InflateToArena", SortPhase::ReadDecompress)]
    #[case::find_and_sort("FindBoundariesAndSort", SortPhase::InMemorySort)]
    #[case::sort_buffer("SortBuffer", SortPhase::InMemorySort)]
    #[case::spill_gather("SpillGather", SortPhase::SpillWrite)]
    #[case::spill_block_compress("SpillBlockCompress", SortPhase::SpillWrite)]
    #[case::spill_write("SpillWrite", SortPhase::SpillWrite)]
    #[case::compress_spill("CompressSpill", SortPhase::SpillWrite)]
    #[case::spill_decompress("SortSpillDecompress", SortPhase::Consolidation)]
    #[case::merge("SortMerge", SortPhase::KWayMerge)]
    #[case::bgzf_compress("BgzfCompress", SortPhase::WriteOutput)]
    #[case::write_file("WriteBgzfFile", SortPhase::WriteOutput)]
    fn sort_phase_classifier_covers_known_sort_steps(
        #[case] step_name: &str,
        #[case] expected: SortPhase,
    ) {
        assert_eq!(SortPhase::from_step_name(step_name), Some(expected));
    }

    #[test]
    fn format_sort_phase_timing_emits_header_all_phases_and_percentages() {
        let mut phase_ns = [0u64; 6];
        phase_ns[SortPhase::ReadDecompress as usize] = 300_000_000; // 0.3s → 30%
        phase_ns[SortPhase::KWayMerge as usize] = 700_000_000; // 0.7s → 70%
        let lines = format_sort_phase_timing(&phase_ns).expect("nonzero total must render");

        assert_eq!(lines[0], "=== Sort Phase Timing ===");
        for phase in SortPhase::ORDER {
            assert!(
                lines.iter().any(|l| l.contains(phase.label())),
                "missing a line for phase {phase:?}; got: {lines:?}"
            );
        }
        assert!(
            lines.iter().any(|l| l.contains("read + decompress") && l.contains("30.0%")),
            "read+decompress should be 30.0%; got: {lines:?}"
        );
        assert!(
            lines.iter().any(|l| l.contains("k-way merge") && l.contains("70.0%")),
            "k-way merge should be 70.0%; got: {lines:?}"
        );
        assert!(
            lines.iter().any(|l| l.contains("total") && l.contains("100.0%")),
            "a total line summing to 100% is expected; got: {lines:?}"
        );
    }

    #[test]
    fn format_sort_phase_timing_is_none_when_no_phase_did_work() {
        assert!(
            format_sort_phase_timing(&[0u64; 6]).is_none(),
            "an all-zero snapshot must render no block"
        );
    }

    #[test]
    fn log_sort_phase_timing_emits_the_block_for_a_populated_snapshot() {
        let _session = capture_logs();
        log_sort_phase_timing(&snapshot_with_steps(vec![
            ("SortMerge", 500_000_000),
            ("WriteBgzfFile", 500_000_000),
        ]));
        let logs = captured();
        assert!(
            logs.iter().any(|l| l.contains("=== Sort Phase Timing ===")),
            "expected the phase-timing header; got: {logs:?}"
        );
        assert!(logs.iter().any(|l| l.contains("k-way merge")), "got: {logs:?}");
        assert!(logs.iter().any(|l| l.contains("write output")), "got: {logs:?}");
    }

    #[test]
    fn log_sort_phase_timing_is_silent_for_an_empty_snapshot() {
        let _session = capture_logs();
        log_sort_phase_timing(&snapshot_with_steps(vec![]));
        assert!(
            !captured().iter().any(|l| l.contains("Sort Phase Timing")),
            "no block should be logged when no sort work was recorded"
        );
    }
}
