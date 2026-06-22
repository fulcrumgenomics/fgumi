//! Sort typed-steps for the unified pipeline.

pub mod and_spill;
pub mod merge;
pub mod protocol;
pub mod spill_decompress;

pub use and_spill::SortAndSpill;
pub use merge::SortMerge;
pub use spill_decompress::SortSpillDecompress;

use std::io;
use std::path::PathBuf;
use std::sync::Arc;

use fgumi_sort::RawExternalSorter;
use parking_lot::Mutex;

use fgumi_pipeline_core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

// ─────────────────────────────────────────────────────────────────────────────
// SortBamFile — Exclusive single-step wrapping legacy sort end-to-end.
// ─────────────────────────────────────────────────────────────────────────────

/// `Exclusive` step that drives a complete legacy
/// `RawExternalSorter::sort(input, output)` call to completion in a
/// single `try_run`.
pub struct SortBamFile {
    sorter: Option<RawExternalSorter>,
    input: PathBuf,
    output: PathBuf,
    stats_out: Arc<Mutex<Option<fgumi_sort::SortStats>>>,
}

impl SortBamFile {
    /// Build a `SortBamFile` step.
    #[must_use]
    pub fn new(
        sorter: RawExternalSorter,
        input: PathBuf,
        output: PathBuf,
        stats_out: Arc<Mutex<Option<fgumi_sort::SortStats>>>,
    ) -> Self {
        Self { sorter: Some(sorter), input, output, stats_out }
    }
}

impl Step for SortBamFile {
    type Input = ();
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SortBamFile",
            kind: StepKind::Exclusive,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let Some(sorter) = self.sorter.take() else {
            return Ok(StepOutcome::Finished);
        };
        let stats = sorter
            .sort(&self.input, &self.output)
            .map_err(|e| io::Error::other(format!("SortBamFile: sort failed: {e:#}")))?;
        *self.stats_out.lock() = Some(stats);
        Ok(StepOutcome::Finished)
    }
}

#[cfg(test)]
pub mod tests;
