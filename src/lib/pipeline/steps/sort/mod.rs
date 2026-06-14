//! Sort typed-steps for the unified pipeline.
//!
//! Two flavors:
//!
//! 1. **[`SortBamFile`]** — `Exclusive` single-step that wraps the
//!    `RawExternalSorter::sort(input, output)` end-to-end. Used by
//!    standalone `fgumi sort`. The chain is literally `[SortBamFile]`
//!    — the framework provides one thread to run `sort()` on.
//!    Concurrency profile is one driver thread + the sort engine's
//!    own `SortWorkerPool` and rayon pools.
//!
//! 2. **Three-step chain** (`SortAndSpill` → `SortSpillDecompress` →
//!    `SortMerge`) — middle-of-chain composition for `runall --sort`
//!    fusion. The three steps split the legacy monolithic `Sort` body
//!    so the spill decompress phase runs on the framework's
//!    work-stealing pool (eliminating the +N `SortWorkerPool::Phase2`
//!    OS-thread oversubscription during merge). Records stream in via
//!    `SortAndSpill`, sorted records stream out the back of `SortMerge`,
//!    and downstream consumers (Decode → Group → ...) run concurrently
//!    with the merge tail.
//!
//! See [`and_spill`], [`spill_decompress`], and [`merge`] for the
//! per-step state machines, and `docs/design/sort-step-split.md` for
//! the locked design.

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

use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

// ─────────────────────────────────────────────────────────────────────────────
// SortBamFile — Exclusive single-step wrapping legacy sort end-to-end.
// ─────────────────────────────────────────────────────────────────────────────

/// `Exclusive` step that drives a complete legacy
/// `RawExternalSorter::sort(input, output)` call to completion in a
/// single `try_run`. The framework gives this step a dedicated owner
/// worker, on which the entire sort runs. Concurrency profile is
/// identical to invoking the legacy sort directly: one driver thread
/// + the sort engine's `SortWorkerPool` and rayon pools.
///
/// The chain at standalone use is just `[SortBamFile]` — there are no
/// other typed-steps. The pipeline framework adds *zero* concurrency
/// vs legacy. The wrapping exists purely so the operation is
/// composable.
///
/// Supports `--write-index` transparently via the `write_index` flag
/// on the underlying `RawExternalSorter`.
pub struct SortBamFile {
    /// `Some` until the first `try_run` consumes it; `None` after.
    /// Guards against a second `try_run` doing the sort again — once the step
    /// reports `Finished` it is dropped from the worklist and not re-dispatched,
    /// but the option makes the no-op path explicit.
    sorter: Option<RawExternalSorter>,
    input: PathBuf,
    output: PathBuf,
    /// Out-parameter slot for the `SortStats` produced by
    /// `RawExternalSorter::sort`. Filled in by `try_run` after the
    /// sort completes; the caller (e.g. standalone `Sort::execute`)
    /// holds an `Arc` clone and reads the totals once the pipeline
    /// returns so we can log the Records-processed / Records-written
    /// / Temporary-chunks Summary block.
    stats_out: Arc<Mutex<Option<fgumi_sort::SortStats>>>,
}

impl SortBamFile {
    /// Build a `SortBamFile` step. The `sorter` should be fully
    /// configured (`memory_limit`, `threads`, `temp_compression`,
    /// `output_compression`, `temp_dirs`, `cell_tag`, `write_index`,
    /// `pg_info`, `initial_capacity`, `async_reader`) — the same way
    /// the standalone `Sort::execute_sort` path configures its sorter.
    ///
    /// `stats_out` is the slot the step writes its `SortStats` into
    /// after `sorter.sort()` returns. Pass `Arc::clone(&slot)` here and
    /// take the inner value from the same slot after `Pipeline::run`
    /// returns to recover the per-run counters.
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
mod tests;
