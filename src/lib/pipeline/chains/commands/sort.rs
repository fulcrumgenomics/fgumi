//! Chain builder support for `Stage::Sort`.
//!
//! The stage-by-stage construction lives in
//! [`crate::pipeline::chains::builder::ChainBuilder`]'s `add_sort` method.
//! This module provides the standalone-sort summary finalize hook
//! ([`SortSummaryFinalizeHook`]) and the [`IndexBamFinalizeHook`] BAI indexer
//! (defined locally — see its doc), which `add_sort` registers for a
//! `SinkSpec::BamWithIndex` request.
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

// ─────────────────────────────────────────────────────────────────────────────
// IndexBamFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline action that builds a BAI index for a finished
/// coordinate-sorted BAM, writing the `.bai` sidecar next to the BAM.
///
/// Reads the BAM via `noodles::bam::fs::index` (walking BGZF block offsets) and
/// writes the sidecar via [`fgumi_bam_io::write_bai_sidecar`], which appends
/// `.bai` to the full output path (samtools convention), so `foo.bam` →
/// `foo.bam.bai` and `foo.sorted` → `foo.sorted.bai`.
///
/// **Invariant:** the BAM at `output_path` must be writer-closed before
/// `finalize` runs — guaranteed by `Pipeline::run` returning (every step's
/// writer has dropped). The hook does not `fsync`: the same process re-reading
/// via the OS page cache sees all flushed bytes on every supported local
/// filesystem.
///
/// `pub(crate)` so [`crate::pipeline::chains::builder::ChainBuilder::add_sort`]
/// can construct and register it. Defined here (rather than re-exported from a
/// sort-cli crate as upstream did) because `main-runall` keeps the sort command
/// local to `crate::commands::sort`.
pub(crate) struct IndexBamFinalizeHook {
    /// Path of the finished coordinate-sorted BAM to index.
    pub(crate) output_path: PathBuf,
}

impl FinalizeHook for IndexBamFinalizeHook {
    /// Build and write the BAI index alongside the BAM.
    ///
    /// # Errors
    ///
    /// Returns an error if indexing or writing the BAI sidecar fails.
    fn finalize(self: Box<Self>) -> Result<()> {
        use std::time::Instant;

        use crate::logging::format_duration;

        let IndexBamFinalizeHook { output_path } = *self;

        info!("Indexing BAM: {}", output_path.display());
        let start = Instant::now();

        // Index + write the `<output>.bai` sidecar in one call; the path is
        // derived by appending `.bai` to the full BAM path (samtools
        // convention), correct for any path — not just `*.bam`.
        let index_path = fgumi_bam_io::write_bai_sidecar(&output_path)?;
        info!("Wrote BAM index: {} ({})", index_path.display(), format_duration(start.elapsed()));

        Ok(())
    }
}
