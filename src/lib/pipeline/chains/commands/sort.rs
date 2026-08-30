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

#[cfg(test)]
mod tests {
    use super::*;

    /// A `log` sink that records every emitted message so a test can assert a
    /// specific `info!` line actually fired.
    ///
    /// This mirrors the `CaptureLogger` pattern in
    /// `crate::commands::common::tests` (that one is private to its own
    /// module, so it cannot be reused here). Installing a process-global
    /// logger is safe under `cargo nextest run` (the repo's standard test
    /// runner), which isolates each test in its own process; a plain `cargo
    /// test` run sharing one process across multiple such loggers would
    /// panic on the second install, same as the existing pattern.
    struct CaptureLogger;

    static CAPTURED_LOGS: std::sync::Mutex<Vec<String>> = std::sync::Mutex::new(Vec::new());
    static CAPTURE_LOGGER: CaptureLogger = CaptureLogger;

    impl log::Log for CaptureLogger {
        fn enabled(&self, _metadata: &log::Metadata) -> bool {
            true
        }

        fn log(&self, record: &log::Record) {
            CAPTURED_LOGS
                .lock()
                .expect("captured-log lock poisoned")
                .push(record.args().to_string());
        }

        fn flush(&self) {}
    }

    /// Install [`CaptureLogger`] as the process-global logger and clear any
    /// records left by an earlier test, so the caller asserts only on what
    /// its own operation emits.
    fn capture_logs() {
        use std::sync::Once;
        static INSTALL: Once = Once::new();
        INSTALL.call_once(|| {
            log::set_logger(&CAPTURE_LOGGER).expect("no logger installed yet in this test process");
            log::set_max_level(log::LevelFilter::Trace);
        });
        CAPTURED_LOGS.lock().expect("captured-log lock poisoned").clear();
    }

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

        let logs = CAPTURED_LOGS.lock().expect("captured-log lock poisoned");
        assert!(
            logs.iter().any(|line| line.contains("Spill runs: 3")),
            "expected a 'Spill runs: 3' log line; got: {logs:?}"
        );
        assert!(
            !logs.iter().any(|line| line.contains("Temporary runs:")),
            "must not emit the old 'Temporary runs:' wording; got: {logs:?}"
        );
    }
}
