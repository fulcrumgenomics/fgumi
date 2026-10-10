//! Sort typed-steps for the declarative chain.

/// `DetachedGroup::Shared` label for the sort's **coordination** driver thread:
/// the serial phase-1 coordination steps (`ReadBlocks` admit, the split front's
/// boundary scan and run sort — `FindTemplateBoundaries` / `SortTemplateRun`, with
/// key extraction on the pool between them — `FindBoundariesAndSort` for the
/// orders not yet split, `SpillGather` range planning — the pool frames and
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
pub mod merge;
pub mod protocol;
pub(crate) mod read_ahead_budget;
pub(crate) mod run_stack;
pub mod sort_buffer;
pub mod spill_block_compress;
pub mod spill_decompress;
pub mod spill_gather;
pub mod spill_read_planner;
pub mod spill_write;
pub mod supply_ledger;

pub use arena_ingest::{
    ArenaBlock, ArenaSortStrategy, CoordinateStrategy, FindBoundariesAndSort, InflateToArena,
    InflatedBlock, QuerynameStrategy, ReadBlocks, TemplateStrategy,
};
pub use arena_split::{
    ArenaKeyExtractor, ExtractKeys, ExtractQuerynameKeys, ExtractTemplateKeys,
    FindQuerynameBoundaries, FindRecordBoundaries, FindTemplateBoundaries, KEY_BATCH_RECORDS,
    KeyBatch, KeyOverlapCensus, QnChunkOrder, QuerynameKeys, RefBatch, RunEnd, SortQuerynameRun,
    SortRun, SortTemplateRun, TemplateKeyBatch, TemplateKeys, TemplateRefBatch,
};
pub use merge::{BlockOutput, MergeBatchBuilder, MergeOutput, RecordBatchOutput, SortMerge};
pub use sort_buffer::{SortBuffer, SortBufferConfig};
pub use spill_block_compress::{SpillBlockCompress, SpillBufPool};
pub use spill_decompress::SortSpillDecompress;
pub use spill_gather::SpillGather;
pub use spill_read_planner::SpillReadPlanner;
pub use spill_write::{SpillRunStats, SpillWrite};
pub use supply_ledger::{SpillSupply, SupplyDiagnostics, SupplyLedger};

#[cfg(all(test, feature = "stress-tests"))]
mod slot_stash_soak;
#[cfg(test)]
pub mod tests;
