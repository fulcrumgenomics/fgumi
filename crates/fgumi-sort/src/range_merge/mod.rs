//! Range-parallel final merge of a sort's spill runs.
//!
//! The final merge is a k-way merge whose leaves are the spill slots (in
//! `file_id` order); range mode runs only when the final run was spilled too,
//! so no residual chunk is held in memory.
//! Range mode partitions the merge's record order (key, then leaf index, then
//! position in the leaf) into disjoint, increasing ranges so each range can be
//! merged on its own and the range outputs concatenated in range order: every
//! record lands in the one range holding it, and within a range ties break by
//! leaf index exactly as in the serial merge. A range bound may fall inside a
//! group of equal keys, so no key, however common, makes one range unbounded.

pub(crate) mod plan;
pub(crate) mod task;
#[cfg(any(test, feature = "test-utils"))]
pub mod test_support;

pub use plan::{MAX_RANGES, RangeBound, RangePlan, RangeSizing};
pub use task::{RangeMerge, RangeOut, RangeStats, RangeTask, RecordSink};
