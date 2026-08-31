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
}
