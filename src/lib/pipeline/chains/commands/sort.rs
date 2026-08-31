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
            info!("Temporary runs: {}", stats.runs_written);
        }
        info!("Output: {}", output_path.display());
        timer.log_completion(stats.total_records);
        Ok(())
    }
}
