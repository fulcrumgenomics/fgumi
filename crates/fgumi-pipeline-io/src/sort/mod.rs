//! Sort typed-steps for the unified pipeline.

/// `DetachedGroup::Shared` label for the sort's **coordination** driver thread:
/// the serial phase-1 coordination steps (`ReadBlocks` admit, `FindBoundariesAndSort`
/// sort/seal, `SpillGather` framing) plus the phase-2 `SortMerge`. One dedicated
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
pub mod compress_spill;
pub mod merge;
pub mod protocol;
pub mod sort_buffer;
pub mod spill_block_compress;
pub mod spill_decompress;
pub mod spill_gather;
pub mod spill_write;

pub use arena_ingest::{
    ArenaBlock, ArenaSortStrategy, CoordinateStrategy, FindBoundariesAndSort, InflateToArena,
    InflatedBlock, QuerynameStrategy, ReadBlocks, TemplateStrategy,
};
pub use compress_spill::CompressSpill;
pub use merge::{BlockOutput, MergeBatchBuilder, MergeOutput, RecordBatchOutput, SortMerge};
pub use sort_buffer::SortBuffer;
pub use spill_block_compress::SpillBlockCompress;
pub use spill_decompress::{SortDecompressTuning, SortSpillDecompress};
pub use spill_gather::SpillGather;
pub use spill_write::SpillWrite;

#[cfg(test)]
pub mod tests;
