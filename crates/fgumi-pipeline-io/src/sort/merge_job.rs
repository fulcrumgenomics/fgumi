//! The work `SortMerge` posts for pool threads and the results they return.
//!
//! `SortMerge` runs on one detached thread; work it can hand off travels on its
//! second output branch as a [`MergeTicket`] to the `MergeWork` step (a
//! [`fgumi_pipeline_core::FanoutWorker`] on the pool), which delivers each
//! [`MergeResult`] to the `OrderedInbox` the merge takes them back from, in
//! ordinal order. A gather ticket carries its [`MergeJob`]; a range-mode
//! ticket claims the oldest unstarted job from the merge's
//! claim table when a worker runs it, and runs that job to completion.

use std::io;
use std::marker::PhantomData;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use fgumi_pipeline_core::outputs::Single;
use fgumi_pipeline_core::{
    Chain, Dispatched, FanoutWorker, HeapSize, Job, JobStep, OrderedInbox, SliceBudget,
    external_wait,
};
use fgumi_sort::{RangeOut, RangeStats, RangeTask, RecordSink};

use crate::sort::merge::claim::ClaimTable;
use crate::sort::merge::{FastBlock, MergeBatchBuilder, MergeOutput, gather_fast_block};
use crate::sort::protocol::MemoryChunkErased;
use crate::sort::range_cut::{BatchPlan, Segment};

/// The slice budget a merge job runs under. A gather or an assembly is one
/// slice and ignores it; a range merge takes `records` records per slice (the
/// merge's own thread runs one slice per dispatch; a `MergeWork` worker runs
/// every slice of a claimed job in one dispatch).
pub(crate) const MERGE_SLICE: SliceBudget =
    SliceBudget { bytes: 4 * 1024 * 1024, records: 16 * 1024 };

/// The `MergeWork` pool step: runs the tickets on `SortMerge`'s branch 1 and
/// delivers their results to `inbox`, the one the merge was built `with_inbox`.
/// Crate-private: production code gets the matched pair from
/// [`SortMerge::with_pooled_gather`](crate::sort::SortMerge::with_pooled_gather).
/// Every `SortMerge`'s branch 1 is wired to one; without a fan-out it never
/// sees a ticket and finishes when the merge does.
#[must_use]
pub(crate) fn merge_work<O: MergeOutput>(
    inbox: Arc<OrderedInbox<MergeResult<O>>>,
) -> FanoutWorker<MergeTicket<O>> {
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
    merge: Chain<'_, (O::Item, MergeTicket<O>)>,
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

/// A unit of merge work.
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
    /// Range-parallel merge: merge key range `range` of the spill runs (range
    /// mode spills the final run, so nothing is merged from memory). The first
    /// slice reads the range's spans from disk; later slices merge up to the
    /// budget's record count each.
    Range {
        /// The range's merge.
        task: Box<dyn RangeTask>,
        /// The range's position in range order.
        range: usize,
        /// The job's inbox ordinal.
        job: u64,
        /// Whether `task.read()` has run.
        read: bool,
        /// Diagnostic timing (`--sort-stats` only).
        timing: Option<JobTiming>,
        /// Recycled output buffers to build the range's output in, taken
        /// when the job starts.
        buffers: Option<Arc<RangeBufPool>>,
    },
    /// Range-parallel merge: copy the records of the output batches closed
    /// while one merged range was cut out of the merged ranges, one `O::Item`
    /// per batch with `batch_serial == batch.ordinal`, in batch order.
    Assemble {
        /// The batches, in output order: their ordinals and record spans.
        batches: Vec<BatchPlan>,
        /// Records over every batch.
        records: u64,
        /// Framed bytes over every batch.
        framed_bytes: u64,
        /// The job's inbox ordinal (not a batch's output ordinal).
        job: u64,
        /// Record capacity of each block's builder (the batch's record cap).
        target_batch_count: usize,
        /// The builder holding the first batch's leading records, when that
        /// batch continues a fused range's open batch (`BatchPlan::lead`).
        lead: Option<O::Builder>,
        /// Diagnostic timing (`--sort-stats` only).
        timing: Option<JobTiming>,
        /// Where a merged range's output buffers go once this job drops the
        /// range's last segment.
        recycle: Option<Arc<RangeBufPool>>,
        /// Carries the output type `O` (zero-sized).
        _output: PhantomData<fn() -> O>,
    },
}

/// Which range-mode job a [`JobTimes`] bucket or a claim counts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum JobKind {
    /// A range merge.
    Range = 0,
    /// A range's batch assembly.
    Assemble = 1,
}

/// Who ran a range-mode job.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Runner {
    /// A `MergeWork` worker.
    Worker = 0,
    /// The merge's own thread.
    Merge = 1,
}

/// Diagnostic job latencies of a range-mode merge, summed per kind and
/// runner: jobs, nanoseconds from post to claim (queue wait), from claim to
/// done (run, including any time between slices), and the largest queue
/// wait.
#[derive(Default)]
pub(crate) struct JobTimes {
    cells: [[JobTimesCell; 2]; 2],
}

#[derive(Default)]
struct JobTimesCell {
    jobs: AtomicU64,
    queue_ns: AtomicU64,
    run_ns: AtomicU64,
    max_queue_ns: AtomicU64,
}

impl JobTimes {
    fn record(&self, kind: JobKind, runner: Runner, queue_ns: u64, run_ns: u64) {
        let c = &self.cells[kind as usize][runner as usize];
        c.jobs.fetch_add(1, Ordering::Relaxed);
        c.queue_ns.fetch_add(queue_ns, Ordering::Relaxed);
        c.run_ns.fetch_add(run_ns, Ordering::Relaxed);
        c.max_queue_ns.fetch_max(queue_ns, Ordering::Relaxed);
    }

    /// One `--sort-stats` line per kind and runner that ran any job.
    #[allow(clippy::cast_precision_loss, reason = "diagnostic means")]
    pub(crate) fn log(&self) {
        for (kind, kname) in [(JobKind::Range, "range"), (JobKind::Assemble, "assembly")] {
            for (runner, rname) in [(Runner::Worker, "workers"), (Runner::Merge, "merge thread")] {
                let c = &self.cells[kind as usize][runner as usize];
                let jobs = c.jobs.load(Ordering::Relaxed);
                if jobs == 0 {
                    continue;
                }
                let ms = |ns: u64| ns as f64 / 1e6;
                let q = c.queue_ns.load(Ordering::Relaxed);
                let r = c.run_ns.load(Ordering::Relaxed);
                log::info!(
                    "Range jobs ({kname}, {rname}): {jobs} jobs, queue wait (post to claim) {:.1} \
                     s total ({:.3} ms mean, {:.1} ms max), run {:.1} s total ({:.3} ms mean)",
                    ms(q) / 1e3,
                    ms(q) / jobs as f64,
                    ms(c.max_queue_ns.load(Ordering::Relaxed)),
                    ms(r) / 1e3,
                    ms(r) / jobs as f64,
                );
            }
        }
    }
}

/// A job's diagnostic timestamps, recorded into `sink` when it finishes.
pub struct JobTiming {
    posted: Instant,
    claimed: Option<Instant>,
    runner: Runner,
    sink: Arc<JobTimes>,
}

impl JobTiming {
    /// A job posted now.
    pub(crate) fn posted(sink: &Arc<JobTimes>) -> Self {
        Self {
            posted: Instant::now(),
            claimed: None,
            runner: Runner::Worker,
            sink: Arc::clone(sink),
        }
    }

    fn finish(&self, kind: JobKind) {
        let now = Instant::now();
        let claimed = self.claimed.unwrap_or(now);
        let ns = |d: std::time::Duration| u64::try_from(d.as_nanos()).unwrap_or(u64::MAX);
        self.sink.record(kind, self.runner, ns(claimed - self.posted), ns(now - claimed));
    }
}

/// The result of a [`MergeJob`].
pub enum MergeResult<O: MergeOutput> {
    /// A gathered output block.
    Block(O::Item),
    /// A range's assembled output blocks, in batch order.
    Blocks(Vec<O::Item>),
    /// A merged key range.
    Range {
        /// The range's position in range order.
        range: usize,
        /// The range's merged records.
        out: Arc<RangeOut>,
    },
    /// A key range the merge's own thread merged straight into output
    /// batches (a fused range): the blocks it closed, in batch order.
    Fused {
        /// The range's position in range order.
        range: usize,
        /// The range's counters.
        stats: RangeStats,
        /// The batches the range closed.
        blocks: Vec<O::Item>,
        /// Records in `blocks`.
        records: u64,
    },
}

impl<O: MergeOutput> MergeJob<O> {
    /// The job's range-mode kind (`None` for a gather).
    pub(crate) fn kind(&self) -> Option<JobKind> {
        match self {
            Self::Gather { .. } => None,
            Self::Range { .. } => Some(JobKind::Range),
            Self::Assemble { .. } => Some(JobKind::Assemble),
        }
    }

    /// Run one slice of this range job merging into `sink` instead of into
    /// the range's output buffer (the merge thread's fused path): the first
    /// slice reads the range's spans, later slices merge up to the budget's
    /// record count each. Returns the range's position and counters once the
    /// range is done.
    ///
    /// # Errors
    ///
    /// The range's read or merge error, a sink error, or a job that is not a
    /// range merge.
    pub(crate) fn run_fused_slice(
        &mut self,
        budget: SliceBudget,
        sink: &mut dyn RecordSink,
    ) -> io::Result<Option<(usize, RangeStats)>> {
        let Self::Range { task, range, read, timing, buffers, .. } = self else {
            return Err(io::Error::other("only a range merge can run fused"));
        };
        if !*read {
            // A fused range builds no output buffer: leave the pool's
            // buffers to the ranges that do.
            *buffers = None;
            external_wait(|| task.read())?;
            *read = true;
            return Ok(None);
        }
        let Some(stats) = task.merge_into(sink, budget.records.max(1))? else {
            return Ok(None);
        };
        if let Some(t) = timing.as_ref() {
            t.finish(JobKind::Range);
        }
        Ok(Some((*range, stats)))
    }

    /// `runner` claimed this job now: its queue wait ends.
    pub(crate) fn claimed(&mut self, runner: Runner) {
        if let Self::Range { timing: Some(t), .. } | Self::Assemble { timing: Some(t), .. } = self {
            t.claimed = Some(Instant::now());
            t.runner = runner;
        }
    }
}

impl<O: MergeOutput> HeapSize for MergeJob<O> {
    /// The descriptor itself, never zero, so an edge or table that counts it
    /// sees every queued job. A gather's chunk and an assembly's ranges are
    /// shared and already charged elsewhere; an assembly adds its batch and
    /// segment lists. A queued range job has not read anything yet.
    fn heap_size(&self) -> usize {
        match self {
            Self::Gather { .. } | Self::Range { .. } => std::mem::size_of::<Self>(),
            Self::Assemble { batches, .. } => {
                std::mem::size_of::<Self>()
                    + batches
                        .iter()
                        .map(|b| {
                            std::mem::size_of::<BatchPlan>()
                                + std::mem::size_of::<Segment>() * b.segments.len()
                        })
                        .sum::<usize>()
            }
        }
    }
}

impl<O: MergeOutput> HeapSize for MergeResult<O> {
    fn heap_size(&self) -> usize {
        match self {
            Self::Block(item) => item.heap_size().max(1),
            Self::Blocks(items) | Self::Fused { blocks: items, .. } => {
                items.iter().map(HeapSize::heap_size).sum::<usize>().max(1)
            }
            Self::Range { out, .. } => {
                (out.framed.capacity() + std::mem::size_of::<u32>() * out.lens.capacity()).max(1)
            }
        }
    }
}

impl<O: MergeOutput> Job for MergeJob<O> {
    type Output = MergeResult<O>;

    fn ordinal(&self) -> u64 {
        match self {
            Self::Gather { ordinal, .. } => *ordinal,
            Self::Range { job, .. } | Self::Assemble { job, .. } => *job,
        }
    }

    fn run_slice(&mut self, budget: SliceBudget) -> io::Result<JobStep<Self::Output>> {
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
            Self::Range { task, range, read, timing, buffers, .. } => {
                if !*read {
                    if let Some((framed, lens)) = buffers.take().and_then(|pool| pool.take()) {
                        task.reuse_output_buffers(framed, lens);
                    }
                    external_wait(|| task.read())?;
                    *read = true;
                    return Ok(JobStep::Continue);
                }
                // `merge(0)` makes no progress, so a zero budget still merges one.
                Ok(match task.merge(budget.records.max(1))? {
                    Some(out) => {
                        if let Some(t) = timing.as_ref() {
                            t.finish(JobKind::Range);
                        }
                        JobStep::Done(MergeResult::Range { range: *range, out: Arc::new(out) })
                    }
                    None => JobStep::Continue,
                })
            }
            Self::Assemble { batches, target_batch_count, lead, timing, recycle, .. } => {
                let items = batches
                    .iter()
                    .map(|batch| {
                        let lead = if batch.lead {
                            Some(lead.take().ok_or_else(|| {
                                io::Error::other(format!(
                                    "batch {} continues a lead builder the job does not hold",
                                    batch.ordinal
                                ))
                            })?)
                        } else {
                            None
                        };
                        assemble_batch::<O>(batch, *target_batch_count, lead)
                    })
                    .collect::<io::Result<Vec<_>>>()?;
                if let Some(pool) = recycle.as_ref() {
                    // A range's segments end in assembly jobs only; whichever
                    // job drops the last one recycles the range's buffers.
                    for segment in std::mem::take(batches).into_iter().flat_map(|b| b.segments) {
                        if let Some(out) = Arc::into_inner(segment.range) {
                            pool.put(out);
                        }
                    }
                }
                if let Some(t) = timing.as_ref() {
                    t.finish(JobKind::Assemble);
                }
                Ok(JobStep::Done(MergeResult::Blocks(items)))
            }
        }
    }
}

/// Merged ranges' output buffers (`RangeOut::framed` / `lens`) kept for
/// the next range jobs to build their output in. A range's output is tens of
/// MiB, above the allocator's large-object size, so a fresh buffer per range
/// is returned to the OS and page-faulted in again every time.
///
/// A buffer is created only when the pool is empty, so the buffers in
/// existence never outnumber the ranges that held output at once (bounded
/// by the window); the pool keeps at most `max_bufs`, and drops any buffer
/// whose capacity is past `max_framed_bytes` (range mode passes the merged
/// output `range_fit` charges per live range, so a buffer from a range larger
/// than planned does not keep its excess resident).
pub struct RangeBufPool {
    bufs: parking_lot::Mutex<Vec<(Vec<u8>, Vec<u32>)>>,
    max_bufs: usize,
    max_framed_bytes: usize,
    /// Range jobs that took a pooled buffer.
    reused: AtomicU64,
    /// Buffers returned to the pool.
    recycled: AtomicU64,
}

impl RangeBufPool {
    /// An empty pool keeping up to `max_bufs` buffers of at most
    /// `max_framed_bytes` framed capacity each.
    pub(crate) fn new(max_bufs: usize, max_framed_bytes: usize) -> Self {
        Self {
            bufs: parking_lot::Mutex::new(Vec::new()),
            max_bufs,
            max_framed_bytes,
            reused: AtomicU64::new(0),
            recycled: AtomicU64::new(0),
        }
    }

    /// A pooled pair of buffers, if any.
    pub(crate) fn take(&self) -> Option<(Vec<u8>, Vec<u32>)> {
        let bufs = self.bufs.lock().pop();
        if bufs.is_some() {
            self.reused.fetch_add(1, Ordering::Relaxed);
        }
        bufs
    }

    /// Keep `out`'s buffers for reuse when there is room.
    pub(crate) fn put(&self, out: RangeOut) {
        if out.framed.capacity() > self.max_framed_bytes {
            return;
        }
        let mut bufs = self.bufs.lock();
        if bufs.len() < self.max_bufs {
            bufs.push((out.framed, out.lens));
            self.recycled.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// The `--sort-stats` `Range output buffers:` line.
    pub(crate) fn log(&self, ranges: usize) {
        log::info!(
            "Range output buffers: {} of {ranges} ranges reused a recycled buffer, {} buffers \
             recycled, {} pooled at the end (at most {} of {} MiB each)",
            self.reused.load(Ordering::Relaxed),
            self.recycled.load(Ordering::Relaxed),
            self.bufs.lock().len(),
            self.max_bufs,
            self.max_framed_bytes >> 20,
        );
    }

    /// Buffers pooled now (tests).
    #[cfg(test)]
    pub(crate) fn pooled(&self) -> usize {
        self.bufs.lock().len()
    }
}

/// What the job edge carries: one ticket per posted job.
///
/// A gather ticket carries its job (`direct`). A range-mode ticket
/// (`claim`) carries only the merge's claim table: when a worker
/// runs it, it claims the oldest unstarted job there and runs every slice of
/// it in that one dispatch (a started job is never put down behind other
/// work between slices), or, when the merge's own thread already claimed
/// everything, is spent ([`JobStep::Spent`]).
pub struct MergeTicket<O: MergeOutput> {
    /// The table a claim ticket claims from (`None`: a direct ticket).
    table: Option<Arc<ClaimTable<O>>>,
    /// The job carried or claimed, until it finishes.
    job: Option<MergeJob<O>>,
    /// The carried or claimed job's ordinal; `u64::MAX` before a claim.
    ordinal: u64,
}

impl<O: MergeOutput> MergeTicket<O> {
    /// A ticket that runs `job` itself.
    pub(crate) fn direct(job: MergeJob<O>) -> Self {
        Self { ordinal: job.ordinal(), table: None, job: Some(job) }
    }

    /// A ticket that claims from `table` when it runs.
    pub(crate) fn claim(table: &Arc<ClaimTable<O>>) -> Self {
        Self { table: Some(Arc::clone(table)), job: None, ordinal: u64::MAX }
    }
}

impl<O: MergeOutput> HeapSize for MergeTicket<O> {
    /// The descriptor itself (never zero), so the job edge's byte bound and
    /// the deadlock monitor's probe count every queued ticket.
    fn heap_size(&self) -> usize {
        std::mem::size_of::<Self>()
    }
}

impl<O: MergeOutput> Job for MergeTicket<O> {
    type Output = MergeResult<O>;

    fn ordinal(&self) -> u64 {
        self.ordinal
    }

    fn run_slice(&mut self, budget: SliceBudget) -> io::Result<JobStep<Self::Output>> {
        self.run_dispatch(budget, &|| false).map(|d| d.step)
    }

    /// Claims a job (or is spent) and runs it to completion: a job the merge
    /// waits on in order is never put down behind other work. Stops early,
    /// keeping the job, only when the pipeline is done (`stop`).
    fn run_dispatch(
        &mut self,
        budget: SliceBudget,
        stop: &dyn Fn() -> bool,
    ) -> io::Result<Dispatched<Self::Output>> {
        let spent = Dispatched { step: JobStep::Spent, slices: 0 };
        if self.job.is_none()
            && let Some(table) = &self.table
        {
            let Some(job) = table.claim(Runner::Worker) else {
                table.note_spent();
                return Ok(spent);
            };
            self.ordinal = job.ordinal();
            self.job = Some(job);
        }
        let Some(job) = self.job.as_mut() else { return Ok(spent) };
        let mut slices = 0;
        loop {
            let outcome = job.run_slice(budget)?;
            slices += 1;
            match outcome {
                JobStep::Continue => {}
                JobStep::Done(out) => {
                    self.job = None;
                    return Ok(Dispatched { step: JobStep::Done(out), slices });
                }
                JobStep::Spent => {
                    self.job = None;
                    return Ok(Dispatched { step: JobStep::Spent, slices });
                }
            }
            if stop() {
                return Ok(Dispatched { step: JobStep::Continue, slices });
            }
        }
    }
}

/// Copy the records of `batch` out of its segments, in order, into one
/// `O::Item` with `batch_serial == batch.ordinal`: into `lead` (a fused
/// range's open builder, already holding the batch's first records) when
/// given, else into a fresh builder. Each segment is a run of
/// whole `[u32 LE len][body]` records, handed to the builder as such
/// ([`MergeBatchBuilder::push_framed_records`]: each body is pushed as the
/// serial merge pushes a winner, or the span copied at once where that is
/// the builder's framing), so the item is byte-identical to the serial batch.
pub(crate) fn assemble_batch<O: MergeOutput>(
    batch: &BatchPlan,
    target_batch_count: usize,
    lead: Option<O::Builder>,
) -> io::Result<O::Item> {
    let mut builder = lead.unwrap_or_else(|| {
        let bytes_cap = usize::try_from(batch.framed_bytes).unwrap_or(usize::MAX);
        O::Builder::with_capacity(batch.ordinal, bytes_cap, target_batch_count)
    });
    for segment in &batch.segments {
        builder.push_framed_records(&segment.range.framed[segment.start..segment.end])?;
    }
    let records = builder.len() as u64;
    if records != batch.records {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "batch {} holds {records} records, its plan says {}",
                batch.ordinal, batch.records
            ),
        ));
    }
    Ok(builder.build())
}

#[cfg(test)]
mod tests {
    use std::sync::Mutex;

    use fgumi_pipeline_core::Ordered;
    use fgumi_sort::RangeStats;
    use rstest::rstest;

    use super::*;
    use crate::sort::merge::claim::ClaimPolicy;
    use crate::sort::merge::{BlockOutput, RecordBatchOutput};
    use crate::sort::range_cut::Segment;
    use crate::types::{DecompressedBlock, RecordBatch};

    /// Which [`MergeOutput`] an assembly test builds.
    #[derive(Clone, Copy, Debug)]
    enum OutKind {
        Block,
        RecordBatch,
    }

    /// The contents and ordinal of an output item, for comparing two items: a
    /// block's framed bytes as one element, a record batch's bodies one per
    /// record (so record boundaries are compared too).
    trait ItemBytes {
        fn bytes_and_ordinal(&self) -> (Vec<Vec<u8>>, u64);
    }

    impl ItemBytes for DecompressedBlock {
        fn bytes_and_ordinal(&self) -> (Vec<Vec<u8>>, u64) {
            (vec![self.bytes.clone()], self.ordinal())
        }
    }

    impl ItemBytes for RecordBatch {
        fn bytes_and_ordinal(&self) -> (Vec<Vec<u8>>, u64) {
            (self.iter_record_bytes().map(<[u8]>::to_vec).collect(), self.ordinal())
        }
    }

    /// The body of record `i` of range `r`: `len` bytes of a pattern unique to
    /// the record.
    fn body(r: usize, i: usize, len: u32) -> Vec<u8> {
        (0..len as usize).map(|k| (r * 31 + i * 7 + k).to_le_bytes()[0]).collect()
    }

    /// A merged range whose record `i` has body length `lens[i]`.
    fn make_range(r: usize, lens: &[u32]) -> Arc<RangeOut> {
        crate::sort::range_cut::range_out_for_test(lens, |i, len| body(r, i, len))
    }

    /// The segment of `range` covering records `first..first + n`.
    fn segment(range: &Arc<RangeOut>, range_idx: usize, first: usize, n: usize) -> Segment {
        let offset =
            |rec: usize| -> usize { range.lens[..rec].iter().map(|&l| 4 + l as usize).sum() };
        Segment {
            range: Arc::clone(range),
            range_idx,
            start: offset(first),
            end: offset(first + n),
            records: u32::try_from(n).unwrap(),
        }
    }

    fn check_assemble<O: MergeOutput>(
        ranges: &[Vec<u32>],
        spans: &[(usize, usize, usize)],
        ordinal: u64,
    ) where
        O::Item: ItemBytes,
    {
        let outs: Vec<_> = ranges.iter().enumerate().map(|(r, lens)| make_range(r, lens)).collect();
        let segments: Vec<_> =
            spans.iter().map(|&(r, first, n)| segment(&outs[r], r, first, n)).collect();
        let records = segments.iter().map(|s| u64::from(s.records)).sum();
        let framed_bytes = segments.iter().map(|s| (s.end - s.start) as u64).sum();

        // The serial reference: the same bodies pushed into a fresh builder.
        let mut serial = O::Builder::with_capacity(ordinal, 0, 16);
        for &(r, first, n) in spans {
            for (i, &len) in ranges[r].iter().enumerate().skip(first).take(n) {
                serial.push_record_bytes(&body(r, i, len)).unwrap();
            }
        }
        let want = serial.build();

        let mut job = MergeJob::<O>::Assemble {
            batches: vec![BatchPlan { ordinal, lead: false, segments, records, framed_bytes }],
            records,
            framed_bytes,
            job: 99,
            target_batch_count: 16,
            lead: None,
            timing: None,
            recycle: None,
            _output: PhantomData,
        };
        assert_eq!(job.ordinal(), 99, "the inbox ordinal is the job's, not the batch's");
        let JobStep::Done(MergeResult::Blocks(mut got)) = job.run_slice(MERGE_SLICE).unwrap()
        else {
            panic!("an assembly finishes in one slice with its blocks");
        };
        assert_eq!(got.len(), 1);
        let got = got.pop().unwrap();
        assert_eq!(got.bytes_and_ordinal(), want.bytes_and_ordinal());
        assert_eq!(got.ordinal(), ordinal);
    }

    #[rstest]
    #[case::one_segment(vec![vec![10, 20, 30, 40]], vec![(0, 1, 2)], 3)]
    #[case::two_segments(vec![vec![10, 20, 30], vec![5, 6]], vec![(0, 1, 2), (1, 0, 2)], 7)]
    #[case::big_record(vec![vec![8, 200_000, 8]], vec![(0, 0, 3)], 0)]
    fn assemble_matches_serial_builder(
        #[case] ranges: Vec<Vec<u32>>,
        #[case] spans: Vec<(usize, usize, usize)>,
        #[case] ordinal: u64,
        #[values(OutKind::Block, OutKind::RecordBatch)] kind: OutKind,
    ) {
        match kind {
            OutKind::Block => check_assemble::<BlockOutput>(&ranges, &spans, ordinal),
            OutKind::RecordBatch => check_assemble::<RecordBatchOutput>(&ranges, &spans, ordinal),
        }
    }

    /// A batch whose plan claims a record count its segments do not hold is
    /// corrupt, not silently emitted short or long.
    /// Columns: the plan's record count, a needle of the error.
    #[rstest]
    #[case::plan_claims_more(3, "batch 0 holds 2 records, its plan says 3")]
    #[case::plan_claims_fewer(1, "batch 0 holds 2 records, its plan says 1")]
    fn assemble_rejects_a_record_count_mismatch(#[case] claimed: u64, #[case] needle: &str) {
        let out = make_range(0, &[10, 20]);
        let mut seg = segment(&out, 0, 0, 2);
        seg.records = u32::try_from(claimed).unwrap();
        let framed_bytes = (seg.end - seg.start) as u64;
        let batch = BatchPlan {
            ordinal: 0,
            lead: false,
            segments: vec![seg],
            records: claimed,
            framed_bytes,
        };
        let Err(err) = assemble_batch::<BlockOutput>(&batch, 16, None) else {
            panic!("a mismatch is an error");
        };
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "{err}");
        assert!(err.to_string().contains(needle), "{err}");
    }

    /// A fake range that logs its calls and finishes after `records` merged
    /// records, or fails where told to.
    struct FakeRange {
        log: Arc<Mutex<Vec<String>>>,
        records: u64,
        merged: u64,
        fail_read: bool,
        fail_merge: bool,
    }

    impl RangeTask for FakeRange {
        fn read(&mut self) -> io::Result<()> {
            self.log.lock().unwrap().push("read".to_owned());
            if self.fail_read { Err(io::Error::other("fake read failed")) } else { Ok(()) }
        }

        fn merge(&mut self, max_records: u64) -> io::Result<Option<RangeOut>> {
            self.log.lock().unwrap().push(format!("merge {max_records}"));
            if self.fail_merge {
                return Err(io::Error::other("fake merge failed"));
            }
            self.merged = (self.merged + max_records).min(self.records);
            if self.merged < self.records {
                return Ok(None);
            }
            Ok(Some(RangeOut {
                framed: vec![0; 8],
                lens: vec![4],
                stats: RangeStats { records: self.records, ..RangeStats::default() },
            }))
        }
    }

    fn range_job(
        records: u64,
        fail_read: bool,
        fail_merge: bool,
    ) -> (MergeJob<BlockOutput>, Arc<Mutex<Vec<String>>>) {
        let log = Arc::new(Mutex::new(Vec::new()));
        let task = FakeRange { log: Arc::clone(&log), records, merged: 0, fail_read, fail_merge };
        (
            MergeJob::Range {
                task: Box::new(task),
                range: 5,
                job: 42,
                read: false,
                timing: None,
                buffers: None,
            },
            log,
        )
    }

    #[rstest]
    #[case::one_merge_slice(10, 16, 2, vec!["read", "merge 16"])]
    #[case::three_merge_slices(25, 10, 4, vec!["read", "merge 10", "merge 10", "merge 10"])]
    #[case::zero_budget_still_progresses(2, 0, 3, vec!["read", "merge 1", "merge 1"])]
    fn range_job_reads_then_merges_in_slices(
        #[case] records: u64,
        #[case] budget_records: u64,
        #[case] expected_slices: usize,
        #[case] expected_log: Vec<&str>,
    ) {
        let (mut job, log) = range_job(records, false, false);
        assert_eq!(job.ordinal(), 42);
        let budget = SliceBudget { bytes: 1 << 20, records: budget_records };
        let mut slices = 0;
        let out = loop {
            slices += 1;
            let step = job.run_slice(budget).unwrap();
            if slices == 1 {
                assert_eq!(*log.lock().unwrap(), ["read"], "the first slice only reads");
            }
            match step {
                JobStep::Continue => assert!(slices < 100, "the job never finished"),
                JobStep::Done(result) => break result,
                JobStep::Spent => panic!("a range job always has work"),
            }
        };
        assert_eq!(slices, expected_slices);
        assert_eq!(*log.lock().unwrap(), expected_log);
        let MergeResult::Range { range, out } = out else { panic!("a range job yields a range") };
        assert_eq!(range, 5);
        assert_eq!(out.stats.records, records);
        assert_eq!(out.lens, [4]);
    }

    /// A ticket runs its job to completion in one dispatch, unless the
    /// pipeline is done: then it stops after the slice in progress and keeps
    /// the job, which a later dispatch finishes.
    /// Columns: whether `stop` reports the pipeline done, the log after the
    /// first dispatch, whether that dispatch finished the job.
    #[rstest]
    #[case::runs_to_completion(false, vec!["read", "merge 10", "merge 10", "merge 10"], true)]
    #[case::stops_between_slices(true, vec!["read"], false)]
    fn ticket_runs_to_completion_unless_stopped(
        #[case] stopped: bool,
        #[case] expected_first_log: Vec<&str>,
        #[case] expected_done: bool,
    ) {
        let (job, log) = range_job(25, false, false);
        let mut ticket = MergeTicket::direct(job);
        let budget = SliceBudget { bytes: 1 << 20, records: 10 };
        let first = ticket.run_dispatch(budget, &|| stopped).unwrap();
        assert_eq!(*log.lock().unwrap(), expected_first_log);
        assert_eq!(first.slices, expected_first_log.len() as u64, "every slice run is reported");
        assert_eq!(matches!(first.step, JobStep::Done(_)), expected_done);
        if !expected_done {
            assert!(matches!(first.step, JobStep::Continue), "a stopped ticket keeps its job");
            let JobStep::Done(MergeResult::Range { out, .. }) =
                ticket.run_dispatch(budget, &|| false).unwrap().step
            else {
                panic!("the next dispatch finishes the kept job");
            };
            assert_eq!(out.stats.records, 25);
        }
    }

    #[rstest]
    #[case::read_error(true, false, 1, "fake read failed")]
    #[case::merge_error(false, true, 2, "fake merge failed")]
    fn range_job_error_is_reported(
        #[case] fail_read: bool,
        #[case] fail_merge: bool,
        #[case] failing_slice: usize,
        #[case] expected_message: &str,
    ) {
        let (mut job, _log) = range_job(10, fail_read, fail_merge);
        let budget = SliceBudget { bytes: 1 << 20, records: 16 };
        for _ in 1..failing_slice {
            assert!(matches!(job.run_slice(budget), Ok(JobStep::Continue)));
        }
        let Err(err) = job.run_slice(budget) else { panic!("the slice should fail") };
        assert_eq!(err.to_string(), expected_message);
    }

    /// One assembly job over the batches a range closed yields exactly the
    /// items the per-batch assemblies it replaces yield, in batch order.
    /// Columns: range body lengths, the batches as `(range, first, n)` spans.
    #[rstest]
    #[case::one_batch(vec![vec![10, 20, 30]], vec![vec![(0, 0, 3)]])]
    #[case::three_batches(
        vec![vec![10, 20, 30, 40, 50], vec![6, 7]],
        vec![vec![(0, 0, 2)], vec![(0, 2, 3), (1, 0, 1)], vec![(1, 1, 1)]]
    )]
    fn per_range_assembly_equals_per_batch(
        #[case] ranges: Vec<Vec<u32>>,
        #[case] batches: Vec<Vec<(usize, usize, usize)>>,
        #[values(OutKind::Block, OutKind::RecordBatch)] kind: OutKind,
    ) {
        fn check<O: MergeOutput>(ranges: &[Vec<u32>], batches: &[Vec<(usize, usize, usize)>])
        where
            O::Item: ItemBytes,
        {
            let outs: Vec<_> =
                ranges.iter().enumerate().map(|(r, lens)| make_range(r, lens)).collect();
            let plans = || -> Vec<BatchPlan> {
                batches
                    .iter()
                    .enumerate()
                    .map(|(i, spans)| {
                        let segments: Vec<_> = spans
                            .iter()
                            .map(|&(r, first, n)| segment(&outs[r], r, first, n))
                            .collect();
                        BatchPlan {
                            ordinal: 10 + i as u64,
                            lead: false,
                            records: segments.iter().map(|s| u64::from(s.records)).sum(),
                            framed_bytes: segments.iter().map(|s| (s.end - s.start) as u64).sum(),
                            segments,
                        }
                    })
                    .collect()
            };
            let want: Vec<_> = plans()
                .iter()
                .map(|b| assemble_batch::<O>(b, 16, None).unwrap().bytes_and_ordinal())
                .collect();
            let batches = plans();
            let records = batches.iter().map(|b| b.records).sum();
            let framed_bytes = batches.iter().map(|b| b.framed_bytes).sum();
            let mut job = MergeJob::<O>::Assemble {
                batches,
                records,
                framed_bytes,
                job: 3,
                target_batch_count: 16,
                lead: None,
                timing: None,
                recycle: None,
                _output: PhantomData,
            };
            let JobStep::Done(MergeResult::Blocks(got)) = job.run_slice(MERGE_SLICE).unwrap()
            else {
                panic!("an assembly finishes in one slice with its blocks");
            };
            let got: Vec<_> = got.iter().map(ItemBytes::bytes_and_ordinal).collect();
            assert_eq!(got, want);
        }
        match kind {
            OutKind::Block => check::<BlockOutput>(&ranges, &batches),
            OutKind::RecordBatch => check::<RecordBatchOutput>(&ranges, &batches),
        }
    }

    /// A claim ticket run by `MergeWork` claims the oldest unstarted job and
    /// runs every slice of it in one dispatch (the worker keeps nothing in
    /// progress); a ticket that finds the table empty is spent and delivers
    /// nothing. Columns: claim policy, whether the worker gets the job.
    #[rstest]
    #[case::auto(ClaimPolicy::Auto, true)]
    #[case::merge_only(ClaimPolicy::Merge, false)]
    fn worker_runs_a_claimed_range_in_one_dispatch(
        #[case] policy: ClaimPolicy,
        #[case] worker_claims: bool,
    ) {
        use fgumi_pipeline_core::StepOutcome;
        use fgumi_pipeline_core::testing::StepProbe;

        let inbox = OrderedInbox::<MergeResult<BlockOutput>>::new(true);
        let mut work = merge_work::<BlockOutput>(Arc::clone(&inbox));
        let probe = StepProbe::new(&work);
        let table = Arc::new(ClaimTable::<BlockOutput>::new(policy));
        // 100,000 records at `MERGE_SLICE`'s 16,384 per slice: 1 read + 7 merges.
        let (mut job, log) = range_job(100_000, false, false);
        if let MergeJob::Range { job: ordinal, .. } = &mut job {
            *ordinal = 0;
        }
        table.post(job);
        probe.push_input(MergeTicket::claim(&table));
        probe.push_input(MergeTicket::claim(&table));
        probe.close_input();
        assert_eq!(probe.try_run(&mut work).unwrap(), StepOutcome::Progress);
        assert_eq!(inbox.pending(), usize::from(worker_claims), "delivered in one dispatch");
        if worker_claims {
            assert_eq!(log.lock().unwrap().len(), 1 + 7, "every slice ran");
            let Some((0, MergeResult::Range { range: 5, .. })) = inbox.try_take_next() else {
                panic!("the range's result at its ordinal");
            };
        }
        // The second ticket finds nothing (claimed, or barred) and is spent.
        assert_eq!(probe.try_run(&mut work).unwrap(), StepOutcome::Progress);
        assert_eq!(inbox.pending(), 0);
        assert_eq!(probe.try_run(&mut work).unwrap(), StepOutcome::Finished);
        assert_eq!(table.len(), usize::from(!worker_claims));
    }

    /// The pool keeps at most `max_bufs` buffers and drops one whose framed
    /// capacity is past `max_framed_bytes`; only a kept buffer counts as
    /// recycled. Each case puts outputs of these framed capacities into a
    /// pool of 2 buffers of at most 1024 bytes.
    /// Columns: the capacities put, buffers pooled, buffers recycled.
    #[rstest]
    #[case::keeps_up_to_the_cap(vec![512, 1024], 2, 2)]
    #[case::keeps_at_most_max_bufs(vec![16, 16, 16], 2, 2)]
    #[case::drops_an_oversized_buffer(vec![1025], 0, 0)]
    #[case::drops_only_the_oversized_one(vec![2048, 64], 1, 1)]
    fn buf_pool_put_is_bounded(
        #[case] capacities: Vec<usize>,
        #[case] expected_pooled: usize,
        #[case] expected_recycled: u64,
    ) {
        let pool = RangeBufPool::new(2, 1024);
        for capacity in capacities {
            let framed = Vec::with_capacity(capacity);
            assert!(framed.capacity() >= capacity);
            pool.put(RangeOut { framed, lens: Vec::new(), stats: RangeStats::default() });
        }
        assert_eq!(pool.pooled(), expected_pooled);
        assert_eq!(pool.recycled.load(Ordering::Relaxed), expected_recycled);
    }

    /// The assembly job that drops a range's last segment recycles the
    /// range's buffers into the pool (one still referenced elsewhere is not),
    /// and the next range job takes them when it starts. Columns: whether
    /// another handle on the range outlives the job, buffers pooled after it.
    #[rstest]
    #[case::last_segment_recycles(false, 1)]
    #[case::shared_range_is_kept(true, 0)]
    fn assembly_recycles_a_ranges_buffers(#[case] keep_handle: bool, #[case] pooled: usize) {
        let pool = Arc::new(RangeBufPool::new(4, 1 << 20));
        let out = make_range(0, &[10, 20, 30]);
        let kept = keep_handle.then(|| Arc::clone(&out));
        let segments = vec![segment(&out, 0, 0, 2), segment(&out, 0, 2, 1)];
        drop(out);
        let (records, framed_bytes) = (3, segments.iter().map(|s| (s.end - s.start) as u64).sum());
        let mut job = MergeJob::<BlockOutput>::Assemble {
            batches: vec![BatchPlan { ordinal: 0, lead: false, segments, records, framed_bytes }],
            records,
            framed_bytes,
            job: 0,
            target_batch_count: 16,
            lead: None,
            timing: None,
            recycle: Some(Arc::clone(&pool)),
            _output: PhantomData,
        };
        assert!(matches!(job.run_slice(MERGE_SLICE), Ok(JobStep::Done(MergeResult::Blocks(_)))));
        assert_eq!(pool.pooled(), pooled);
        drop(kept);
        let (mut range, _log) = range_job(1, false, false);
        if let MergeJob::Range { buffers, .. } = &mut range {
            *buffers = Some(Arc::clone(&pool));
        }
        assert!(matches!(range.run_slice(MERGE_SLICE), Ok(JobStep::Continue)), "the read");
        assert_eq!(pool.pooled(), 0, "the range job took the pooled buffers");
    }

    #[rstest]
    #[case::empty(0, 0, 1)]
    #[case::framed_and_lens(100, 10, 140)]
    fn range_result_heap_size(
        #[case] framed_cap: usize,
        #[case] lens_cap: usize,
        #[case] expected: usize,
    ) {
        let out = RangeOut {
            framed: Vec::with_capacity(framed_cap),
            lens: Vec::with_capacity(lens_cap),
            ..RangeOut::default()
        };
        let result = MergeResult::<BlockOutput>::Range { range: 0, out: Arc::new(out) };
        assert_eq!(result.heap_size(), expected);
    }
}
