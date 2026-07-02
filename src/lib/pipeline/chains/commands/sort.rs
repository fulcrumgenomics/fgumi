//! Chain builder support for `Stage::Sort`.
//!
//! The stage-by-stage construction lives in
//! [`crate::pipeline::chains::builder::ChainBuilder`]'s `add_sort` method.
//! This module provides the standalone-sort summary finalize hook
//! ([`SortSummaryFinalizeHook`]) and re-exports the `IndexBamFinalizeHook`
//! BAI indexer from the framework-light `fgumi-sort-cli` crate so `add_sort`
//! can register it for a `SinkSpec::BamWithIndex` request.
//!
//! ## Sort pipeline topology
//!
//! A sole-`[Stage::Sort]` chain runs through the same streaming source →
//! `SortAndSpill` → `SortSpillDecompress` → `SortMerge` → sink pipeline as a
//! fused sort stage, via the normal `add_source` / `add_sink` flow that
//! [`crate::pipeline::chains::build::build_for`] drives for every stage; there
//! is no longer a self-contained file→file sort step or a sort-only chain
//! builder.

use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use anyhow::Result;
use fgumi_cli_common::OperationTimer;
use log::info;

use crate::pipeline::chains::FinalizeHook;

/// Post-pipeline summary for standalone `fgumi sort`.
///
/// Reads the `SortMerge` stats slot (records processed/written + spill-chunk
/// count) and logs the `=== Summary ===` block, then the timer's
/// records-per-second completion line. Registered by
/// [`ChainBuilder::add_sort`] only for a sole-`[Stage::Sort]` chain; the fused
/// `runall` path leaves the slot unset and gets no summary block.
pub(crate) struct SortSummaryFinalizeHook {
    pub(crate) stats_slot: Arc<Mutex<Option<fgumi_sort::SortStats>>>,
    pub(crate) output_path: PathBuf,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for SortSummaryFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let SortSummaryFinalizeHook { stats_slot, output_path, timer } = *self;
        let stats = stats_slot.lock().expect("sort stats slot poisoned").take().unwrap_or_default();
        info!("=== Summary ===");
        info!("Records processed: {}", stats.total_records);
        info!("Records written: {}", stats.output_records);
        if stats.chunks_written > 0 {
            info!("Temporary chunks: {}", stats.chunks_written);
        }
        info!("Output: {}", output_path.display());
        timer.log_completion(stats.total_records);
        Ok(())
    }
}

// `IndexBamFinalizeHook` lives in the `fgumi-sort-cli` crate. Re-export it so
// the umbrella's `ChainBuilder::add_sort` registers it from this module. It
// implements the shared `FinalizeHook` trait (`fgumi-pipeline-core`,
// re-exported as `crate::pipeline::chains::FinalizeHook`), so the umbrella
// registers it in its `Vec<Box<dyn FinalizeHook>>` directly — no umbrella-local
// copy (X1-005).
pub(crate) use fgumi_sort_cli::chains::IndexBamFinalizeHook;
