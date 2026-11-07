//! The work `SortMerge` posts for pool threads and the results they return.
//!
//! `SortMerge` runs on one detached thread; work it can hand off travels on its
//! second output branch as a [`MergeJob`] to the `MergeWork` step (a
//! [`fgumi_pipeline_core::FanoutWorker`] on the pool), which delivers each
//! [`MergeResult`] to the `OrderedInbox` the merge takes them back from, in
//! ordinal order.

use std::io;
use std::marker::PhantomData;
use std::sync::Arc;

use fgumi_pipeline_core::outputs::Single;
use fgumi_pipeline_core::{Chain, FanoutWorker, HeapSize, Job, JobStep, OrderedInbox, SliceBudget};

use crate::sort::merge::{FastBlock, MergeOutput, gather_fast_block};
use crate::sort::protocol::MemoryChunkErased;

/// The slice budget `MergeWork` runs jobs under. A gather is one slice and
/// ignores it.
pub(crate) const MERGE_SLICE: SliceBudget =
    SliceBudget { bytes: 4 * 1024 * 1024, records: 16 * 1024 };

/// The `MergeWork` pool step: runs the jobs on `SortMerge`'s branch 1 and
/// delivers their results to `inbox`, the one the merge was built `with_inbox`.
/// Crate-private: production code gets the matched pair from
/// [`SortMerge::with_pooled_gather`](crate::sort::SortMerge::with_pooled_gather).
/// Every `SortMerge`'s branch 1 is wired to one; without a fan-out it never
/// sees a job and finishes when the merge does.
#[must_use]
pub(crate) fn merge_work<O: MergeOutput>(
    inbox: Arc<OrderedInbox<MergeResult<O>>>,
) -> FanoutWorker<MergeJob<O>> {
    FanoutWorker::new("MergeWork", inbox, MERGE_SLICE)
}

/// Wire a `SortMerge` that never fans out (built without
/// [`SortMerge::with_pooled_gather`](crate::sort::SortMerge::with_pooled_gather),
/// so it has no inbox and posts no job) to an idle `MergeWork` on its job
/// branch, and return its merged-output branch.
///
/// Not for a merge built with `with_pooled_gather`: wire the `MergeWork` that
/// call returns instead. This function cannot see the merge inside `merge`,
/// so it does not check; such a merge's gather results would go to its own
/// inbox, which nothing here delivers to.
pub fn serial_merge_output<O: MergeOutput>(
    merge: Chain<'_, (O::Item, MergeJob<O>)>,
) -> Chain<'_, Single<O::Item>> {
    let fan = merge.into_multi();
    fan.b1.chain(merge_work::<O>(OrderedInbox::new(true))).into_sink_marker();
    fan.b0
}

/// How a gather sizes its output builder: the serial fast path's initial
/// buffer reservation and record capacity, so a pooled gather allocates as the
/// serial one does.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GatherSizing {
    /// Initial byte reservation of the block's builder.
    pub initial_bytes: usize,
    /// Record capacity of the block's builder (the batch's record cap).
    pub target_batch_count: usize,
}

/// A unit of merge work for the pool.
pub enum MergeJob<O: MergeOutput> {
    /// Single-chunk fast path: gather planned block `block` of `chunk` into one
    /// `O::Item` with `batch_serial == ordinal`.
    Gather {
        /// The sorted in-memory chunk, shared by every gather job of the run.
        chunk: Arc<MemoryChunkErased>,
        /// The record range this job frames.
        block: FastBlock,
        /// The block's position in the output (its index in the plan).
        ordinal: u64,
        /// The builder's sizing.
        sizing: GatherSizing,
        /// Carries the output type `O` (zero-sized).
        _output: PhantomData<fn() -> O>,
    },
}

/// The result of a [`MergeJob`].
pub enum MergeResult<O: MergeOutput> {
    /// A gathered output block.
    Block(O::Item),
}

impl<O: MergeOutput> HeapSize for MergeJob<O> {
    /// The descriptor itself: the chunk is shared and already resident, so a
    /// queued job holds only its own bytes (never zero, so the job edge's byte
    /// bound counts every queued job).
    fn heap_size(&self) -> usize {
        std::mem::size_of::<Self>()
    }
}

impl<O: MergeOutput> HeapSize for MergeResult<O> {
    fn heap_size(&self) -> usize {
        match self {
            Self::Block(item) => item.heap_size().max(1),
        }
    }
}

impl<O: MergeOutput> Job for MergeJob<O> {
    type Output = MergeResult<O>;

    fn ordinal(&self) -> u64 {
        match self {
            Self::Gather { ordinal, .. } => *ordinal,
        }
    }

    fn run_slice(&mut self, _budget: SliceBudget) -> io::Result<JobStep<Self::Output>> {
        match self {
            Self::Gather { chunk, block, ordinal, sizing, .. } => {
                let item = gather_fast_block::<O>(
                    chunk,
                    *block,
                    *ordinal,
                    sizing.initial_bytes,
                    sizing.target_batch_count,
                )?;
                Ok(JobStep::Done(MergeResult::Block(item)))
            }
        }
    }
}
