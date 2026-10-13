//! Sort typed-steps for the declarative chain.

/// `DetachedGroup::Shared` label for the sort's **coordination** driver thread:
/// the serial phase-1 coordination steps (`ReadBlocks` admit, the split front's
/// boundary scan and run sort — `Find*Boundaries` / `Sort*Run`, with key
/// extraction on the pool between them — and `SpillGather` range planning — the pool frames and
/// compresses the ranges) plus the phase-2 `SortMerge`. One dedicated
/// thread runs all of them off the pool — the true N+2 model (mirrors main's
/// single main thread). Phase 1 and phase 2 are temporally disjoint, so the
/// coordination steps finish and leave the driver's live set before the merge
/// runs, giving the merge a dedicated thread in phase 2.
pub const SORT_COORD_GROUP: &str = "sort-coord";

/// `DetachedGroup::Shared` label for the sort's **I/O writer** driver thread:
/// `SpillWrite` (phase 1) and `WriteBgzfFile` (phase 2). Isolated from the
/// coordination driver so a write flush never stalls coordination (main's reason
/// for the second dedicated thread).
pub const SORT_IO_GROUP: &str = "sort-io";

pub mod arena_ingest;
pub mod arena_split;
pub mod consolidate_runs;
pub mod consolidation;
pub mod merge;
pub mod protocol;
pub(crate) mod read_ahead_budget;
pub mod run_ledger;
pub(crate) mod run_stack;
pub mod sort_buffer;
pub mod spill_block_compress;
pub mod spill_decompress;
#[cfg(test)]
pub(crate) mod spill_fixtures;
pub mod spill_gather;
pub mod spill_read_planner;
pub mod spill_write;
pub mod supply_ledger;
pub mod write_consolidated_run;

pub use arena_ingest::{
    ArenaBlock, ArenaSortStrategy, CoordinateStrategy, InflateToArena, InflatedBlock,
    QuerynameStrategy, ReadBlocks, TemplateStrategy,
};
pub use arena_split::{
    ArenaKeyExtractor, CoordinateCtx, CoordinateKeys, ExtractCoordinateKeys, ExtractKeys,
    ExtractQuerynameKeys, ExtractTemplateKeys, FindCoordinateBoundaries, FindQuerynameBoundaries,
    FindRecordBoundaries, FindTemplateBoundaries, KEY_BATCH_RECORDS, KeyBatch, KeyOverlapCensus,
    QnChunkOrder, QuerynameKeys, RefBatch, RunEnd, SortCoordinateRun, SortQuerynameRun, SortRun,
    SortTemplateRun, TemplateKeyBatch, TemplateKeys, TemplateRefBatch,
};
pub use consolidate_runs::{ConsolidateRuns, NoConsolidation};
pub use consolidation::ConsolidationRegistry;
pub use merge::{BlockOutput, MergeBatchBuilder, MergeOutput, RecordBatchOutput, SortMerge};
pub use run_ledger::{RUN0_MEAN_RECORD_BYTES, RunLedger, RunPrior, run_cap_for};
pub use run_stack::consolidation_can_trigger;
pub use sort_buffer::{SortBuffer, SortBufferConfig};
pub use spill_block_compress::{SpillBlockCompress, SpillBufPool};
pub use spill_decompress::SortSpillDecompress;
pub use spill_gather::SpillGather;
pub use spill_read_planner::SpillReadPlanner;
pub use spill_write::{SpillRunStats, SpillWrite};
pub use supply_ledger::{SpillSupply, SupplyDiagnostics, SupplyLedger};
pub use write_consolidated_run::WriteConsolidatedRun;

#[cfg(all(test, feature = "stress-tests"))]
mod slot_stash_soak;
/// The sort's test oracle (the stable reference sort and its tie-heavy
/// fixture), shared with other crates' tests.
#[cfg(any(test, feature = "test-utils"))]
pub mod testing;
#[cfg(test)]
pub mod tests;
