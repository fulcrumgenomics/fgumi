//! Range mode of [`SortMerge`]: the range-parallel final merge of a spilled
//! coordinate or template-coordinate sort.
//!
//! The read planner publishes `MergeMode::Range`; the merge then builds a
//! key-range plan over the spill runs' indexes, posts one range merge per key
//! range, cuts the merged ranges, in range order, into the serial merge's
//! batches ([`BatchCutter`]) and posts one assembly job per range for the
//! batches its cut closed. Results come back through the merge's inbox in
//! post order.
//!
//! Every posted job goes into the merge's [`ClaimTable`] and one ticket per
//! job onto the `MergeWork` edge. The merge's own thread and the workers
//! both claim the oldest unstarted job, so the next result the merge needs is
//! always delivered, running on a worker that finishes it within one
//! dispatch, or unstarted and within this thread's reach: the pool keeps its
//! compression-first walk (compression is the one job this thread cannot
//! do), and this thread runs whatever the pool has not reached.
//!
//! A range this thread claims once every earlier range is taken and cut (a
//! fused range) is merged straight into output batches: the cutter's carry is
//! known, so each winner goes into the open batch's builder as in the serial
//! merge's loop, with no merged-range buffer and no assembly job. While no
//! assembly is outstanding the batches it closes are streamed to the output
//! after every slice; otherwise they are emitted at a fresh ordinal posted
//! when its result is taken, behind every assembly posted before it. Its open
//! batch's builder (the *lead*) carries into the next batch to close.
//!
//! This thread posts ranges up to the window whatever the job edge does (a
//! posted job without a ticket is claimed here, or by a worker once the edge
//! drains). Downstream backpressure: while compression's queue holds at least
//! `p` blocks (see [`ClaimTable`]), compression is the bottleneck, so the
//! claim table bars workers from claiming (they compress instead) and this
//! thread keeps merging ahead, in order and taking results, within the
//! window: a block waiting for the output keeps its range live. With nothing
//! left to merge and a block waiting, it compresses that block itself. It
//! runs short slices so its blocks reach the output one or two at a time, and
//! claims assemblies before ranges.

use std::collections::VecDeque;
use std::io;
use std::sync::Arc;
use std::time::Instant;

use fgumi_pipeline_core::step::{StepCtx, StepOutcome};
use fgumi_pipeline_core::{HeldSlot, Job, JobError, JobStep, OrderedInbox, Unpushed};
use fgumi_sort::RecordSink;
use fgumi_sort::{
    CbKey32, RangeMerge, RangePlan, RangeSizing, RangeStats, RawCoordinateKey, SortMergeSlot,
    SortOrder, SpillKeyKind, TemplateKey, TemplateKey24, TertKey32,
};

use super::claim::ClaimTable;
use super::{
    MemoryChunksByKind, MergeBatchBuilder, MergeOutput, SortMerge, SortMergeState, TemplateLane,
    ensure_all_records_merged, slots_key_kind, template_lane_unknown,
};
use crate::output_compress::OutputCompressor;
use crate::sort::merge_job::{
    JobTimes, JobTiming, MERGE_SLICE, MergeJob, MergeResult, MergeTicket, RangeBufPool,
};
use crate::sort::merge_mode::{FitInputs, RangeBudget, RangeFit};
use crate::sort::range_cut::{BatchCutter, BatchPlan};
use fgumi_pipeline_core::{HeapSize, SliceBudget};

/// The range plan over `slots` (every record is in a spill run): `build_driver`'s
/// key dispatch, restricted to the orders whose spill runs carry a run index
/// (coordinate and the template lanes). A template merge's lane is `kind`,
/// the slots' common key kind (`None` when a slot carries none, which fails
/// closed as `build_driver` does).
///
/// # Errors
///
/// A queryname order (its runs are never indexed, so the planner never chooses
/// range mode for it), an unknown or non-template lane for a template merge,
/// or the plan's own validation ([`RangePlan::new`]).
fn build_range_plan(
    sort_order: SortOrder,
    slots: Vec<Arc<SortMergeSlot>>,
    kind: Option<SpillKeyKind>,
    total_records: u64,
    sizing: RangeSizing,
) -> io::Result<Arc<dyn RangeMerge>> {
    fn plan<K>(
        slots: Vec<Arc<SortMergeSlot>>,
        total_records: u64,
        sizing: RangeSizing,
    ) -> io::Result<Arc<dyn RangeMerge>>
    where
        K: fgumi_sort::RawSortKey + Default + Copy + Send + Sync + 'static,
    {
        Ok(Arc::new(RangePlan::<K>::new(slots, total_records, sizing)?))
    }
    match sort_order {
        SortOrder::Coordinate => plan::<RawCoordinateKey>(slots, total_records, sizing),
        SortOrder::Queryname(_) => Err(io::Error::other(format!(
            "SortMerge: a {sort_order:?} merge cannot run in range mode (its spill runs carry no \
             run index)"
        ))),
        SortOrder::TemplateCoordinate => {
            match TemplateLane::of(kind.ok_or_else(template_lane_unknown)?)? {
                TemplateLane::K24 => plan::<TemplateKey24>(slots, total_records, sizing),
                TemplateLane::Cb32 => plan::<CbKey32>(slots, total_records, sizing),
                TemplateLane::Tert32 => plan::<TertKey32>(slots, total_records, sizing),
                TemplateLane::K40 => plan::<TemplateKey>(slots, total_records, sizing),
            }
        }
    }
}

/// The slice the merge thread runs of its own job per dispatch: two output
/// batches' records (about 1 ms), so the blocks it builds reach the output
/// queue faster than compressors drain it, and an output refusal is retried
/// between slices. Workers keep [`MERGE_SLICE`] inside their completion loop.
fn own_slice(target_batch_count: usize) -> SliceBudget {
    SliceBudget {
        bytes: MERGE_SLICE.bytes,
        records: u64::try_from(target_batch_count.saturating_mul(2)).unwrap_or(u64::MAX).max(1),
    }
}

/// A job range mode has posted and not yet taken back, in post (= inbox
/// ordinal) order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Posted {
    /// The merge of range `.0`.
    Range(usize),
    /// The assembly of one range's batches, whose first segment is from range
    /// `first_range` and which hold `records` records.
    Assemble { first_range: usize, records: u64 },
}

/// Counters of a range-mode merge, for `--sort-stats`.
#[derive(Default)]
pub(super) struct RangeModeStats {
    /// The ranges' counters, summed.
    totals: RangeStats,
    /// Records of the largest range.
    largest_range_records: u64,
    /// Ranges the merge's own thread merged (the rest ran on `MergeWork`).
    pub(super) own_ranges: usize,
    /// Batches the merge's own thread assembled.
    pub(super) own_assemblies: usize,
    /// Ranges the merge's own thread merged straight into output batches.
    fused_ranges: usize,
    /// Batches those fused ranges closed.
    fused_batches: usize,
    /// Of those, batches streamed to the output during their range (the
    /// rest were posted with its result).
    streamed_batches: usize,
    /// Output blocks this thread compressed while its output refused them.
    compressed_blocks: u64,
    /// Output blocks the output edge accepted.
    blocks_out: u64,
    /// Diagnostic timers (`--sort-stats` only).
    diag: Option<RangeDiag>,
    /// The most ranges ever live at once (see `RangeModeState::live_floor`).
    #[cfg(test)]
    pub(super) max_live_ranges: usize,
}

/// Why the merge thread ended a dispatch with nothing to do.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Idle {
    /// The next result in order is running on a worker.
    WaitingOnWorker = 0,
    /// The next result in order is unstarted in the claim table, yet this
    /// thread did not claim it. A self-check: zero unless
    /// `FGUMI_SORT_RANGE_CLAIM=pool` bars this thread from claiming.
    WaitingOnNothing = 1,
    /// The output branch refused a block, and the window is full (so this
    /// thread could merge no further ahead).
    OutputHeld = 2,
    /// Nothing in flight to wait on; the window is full, so no range may be
    /// posted.
    WindowFull = 3,
    /// Nothing in flight to wait on; every range is posted.
    Draining = 4,
    /// A ticket is held on the job edge with nothing else to wait on. A
    /// self-check: posting does not wait on the ticket edge, so zero.
    TicketEdge = 5,
}

/// What a dispatch that made no progress saw, for [`Idle::of`].
#[derive(Clone, Copy, Debug)]
struct IdleState {
    /// A block is held on the output branch.
    output_held: bool,
    /// Ranges remain to post (the window, not the plan, stops posting).
    ranges_left: bool,
    /// The next result in order: `None` with nothing posted and not taken,
    /// else whether it is unstarted in the claim table.
    head_unstarted: Option<bool>,
    /// A ticket is held on the job edge.
    ticket_held: bool,
}

impl Idle {
    /// Why a dispatch that made no progress stopped. In such a dispatch the
    /// next result in order is not delivered (it would have been taken) and
    /// this thread has no job of its own (it would have run a slice).
    fn of(s: IdleState) -> Self {
        match s.head_unstarted {
            _ if s.output_held && s.ranges_left => Self::OutputHeld,
            Some(false) => Self::WaitingOnWorker,
            Some(true) => Self::WaitingOnNothing,
            None if s.ticket_held => Self::TicketEdge,
            None if s.ranges_left => Self::WindowFull,
            None => Self::Draining,
        }
    }
}

/// The state a dispatch ended in, charged for the time until the next.
#[derive(Clone, Copy)]
struct DispatchEnd {
    /// Why it stopped (`None`: progress).
    idle: Option<Idle>,
    /// Live ranges.
    live: usize,
    /// The output edge held a refused block.
    backlog: bool,
    /// Compression's queue held at least the gate's width of blocks.
    depth_full: bool,
}

/// Diagnostic timers of the merge thread in range mode.
struct RangeDiag {
    /// Job latencies, shared with every posted job.
    jobs: Arc<JobTimes>,
    /// Nanoseconds inside `emit_range_mode`, in total and running the
    /// thread's own range / assembly slices.
    busy_ns: u64,
    own_range_ns: u64,
    own_assemble_ns: u64,
    compress_ns: u64,
    /// Nanoseconds between dispatches, by why the last one stopped
    /// (`Idle`), and between dispatches that made progress.
    idle_ns: [u64; 6],
    between_ns: u64,
    /// Time-weighted live ranges: nanoseconds spent with `i` live
    /// ranges (the last bucket counts every count at or above it).
    live_ns: Vec<u64>,
    /// Nanoseconds spent with a positive output backlog (blocks built and
    /// refused by the output edge), with compression's queue at least the
    /// gate's width deep, and in total.
    backlog_ns: u64,
    depth_ns: u64,
    total_ns: u64,
    /// End of the last dispatch and its state.
    last_end: Option<(Instant, DispatchEnd)>,
}

impl RangeDiag {
    fn new(window: usize) -> Self {
        Self {
            jobs: Arc::new(JobTimes::default()),
            busy_ns: 0,
            own_range_ns: 0,
            own_assemble_ns: 0,
            compress_ns: 0,
            idle_ns: [0; 6],
            between_ns: 0,
            live_ns: vec![0; window + 3],
            backlog_ns: 0,
            depth_ns: 0,
            total_ns: 0,
            last_end: None,
        }
    }

    fn ns(d: std::time::Duration) -> u64 {
        u64::try_from(d.as_nanos()).unwrap_or(u64::MAX)
    }

    /// A dispatch starts now: charge the gap since the last one.
    fn begin(&mut self, now: Instant) {
        if let Some((end, state)) = self.last_end {
            let gap = Self::ns(now - end);
            match state.idle {
                Some(i) => self.idle_ns[i as usize] += gap,
                None => self.between_ns += gap,
            }
            self.charge(gap, state);
        }
    }

    /// Charge `ns` spent in `state`.
    fn charge(&mut self, ns: u64, state: DispatchEnd) {
        let b = state.live.min(self.live_ns.len() - 1);
        self.live_ns[b] += ns;
        self.total_ns += ns;
        if state.backlog {
            self.backlog_ns += ns;
        }
        if state.depth_full {
            self.depth_ns += ns;
        }
    }

    /// A dispatch that began at `start` ends now in `state` (charged for the
    /// dispatch too).
    fn end(&mut self, start: Instant, state: DispatchEnd) {
        let now = Instant::now();
        let d = Self::ns(now - start);
        self.busy_ns += d;
        self.charge(d, state);
        self.last_end = Some((now, state));
    }

    /// Share of the time spent with `ns` set, in percent.
    #[allow(clippy::cast_precision_loss, reason = "a diagnostic share")]
    fn share(&self, ns: u64) -> f64 {
        if self.total_ns == 0 { 0.0 } else { ns as f64 * 100.0 / self.total_ns as f64 }
    }

    #[allow(clippy::cast_precision_loss, reason = "diagnostic shares")]
    fn log(&self, own_ranges: usize, own_assemblies: usize, compressed: u64) {
        let s = |ns: u64| ns as f64 / 1e9;
        let idle = |i: Idle| s(self.idle_ns[i as usize]);
        let work = self.own_range_ns + self.own_assemble_ns + self.compress_ns;
        log::info!(
            "Range merge thread: busy {:.1} s (own ranges {:.1} s over {own_ranges}, own \
             assemblies {:.1} s over {own_assemblies} batches, compression {:.1} s over \
             {compressed} blocks, other {:.1} s); idle {:.1} s: waiting on worker {:.1} s, \
             output held {:.1} s, window full {:.1} s, draining {:.1} s, waiting on nothing \
             (self-check) {:.1} s, ticket edge (self-check) {:.1} s; between busy dispatches \
             {:.1} s",
            s(self.busy_ns),
            s(self.own_range_ns),
            s(self.own_assemble_ns),
            s(self.compress_ns),
            s(self.busy_ns.saturating_sub(work)),
            s(self.idle_ns.iter().sum()),
            idle(Idle::WaitingOnWorker),
            idle(Idle::OutputHeld),
            idle(Idle::WindowFull),
            idle(Idle::Draining),
            idle(Idle::WaitingOnNothing),
            idle(Idle::TicketEdge),
            s(self.between_ns),
        );
        let total: u64 = self.live_ns.iter().sum();
        let last = self.live_ns.len() - 1;
        let parts: Vec<String> = self
            .live_ns
            .iter()
            .enumerate()
            .filter(|(_, ns)| **ns > 0)
            .map(|(i, &ns)| {
                let pct = if total == 0 { 0.0 } else { ns as f64 * 100.0 / total as f64 };
                if i == last { format!(">={i}: {pct:.1}%") } else { format!("{i}: {pct:.1}%") }
            })
            .collect();
        log::info!("Range window occupancy (live ranges: share of time): {}", parts.join(", "));
        log::info!(
            "Range output backlog: blocks refused by the output edge {:.1}% of the time",
            self.share(self.backlog_ns),
        );
        self.jobs.log();
    }

    /// The `Range output depth:` line: how long compression's queue held at
    /// least the gate's width `width` of blocks, and the claims the gate
    /// refused.
    fn log_depth(&self, width: Option<u64>, gate: bool, refused: u64) {
        let gate = if gate { "on" } else { "off" };
        match width {
            Some(p) => log::info!(
                "Range output depth: compression queue >= {p} blocks {:.1}% of the time; claims \
                 refused by the gate {refused} (gate {gate}: workers claim only below {p})",
                self.share(self.depth_ns),
            ),
            None => log::info!(
                "Range output depth: not counted; claims refused by the gate {refused} (gate \
                 {gate}: workers claim only with no refused block)",
            ),
        }
    }
}

impl RangeModeStats {
    /// Count one merged range's counters.
    fn add(&mut self, s: &RangeStats) {
        self.totals.add(s);
        self.largest_range_records = self.largest_range_records.max(s.records);
    }
}

/// The state of [`SortMergeState::RangeMode`], over the output's batch
/// builder `B`.
pub(super) struct RangeModeState<B: MergeBatchBuilder> {
    /// The key ranges.
    pub(super) plan: Arc<dyn RangeMerge>,
    /// Live ranges kept (see [`Self::live_floor`]), as the planner chose.
    window: usize,
    /// Raw bytes per range the plan was cut for.
    target_range_bytes: u64,
    /// Cuts the merged ranges, in range order, into the serial batches.
    cutter: BatchCutter,
    /// The next range to post.
    next_range: usize,
    /// Ranges taken back from the inbox (and cut), in range order.
    ranges_taken: usize,
    /// Jobs posted and not yet taken back, in post order.
    posted: VecDeque<Posted>,
    /// The next job's inbox ordinal (dense over every posted job).
    next_job: u64,
    /// Tickets owed to the job edge: jobs posted to the claim table whose
    /// ticket is not pushed yet (beyond a held one).
    tickets_owed: usize,
    /// Records emitted on the output branch.
    emitted: u64,
    /// Merged ranges' output buffers, recycled into later range jobs.
    pub(super) buffers: Arc<RangeBufPool>,
    /// The open batch's builder when the batch holds a fused range's records
    /// (see the module docs); it goes to the next batch plan that closes.
    lead: Option<B>,
    /// Batches the fused range in progress closed, in order.
    fused_closed: Vec<B>,
    /// The merge thread's own job is a fused range.
    own_fused: bool,
    /// Blocks built (an assembly's, or streamed from a fused range) and not
    /// yet pushed on the output branch, in output order, each with the range
    /// it keeps live (see [`Self::live_floor`]).
    out_blocks: VecDeque<(usize, B::Item)>,
    /// The range the block held on the output branch keeps live, while one
    /// is held.
    held_range: Option<usize>,
    /// `out_blocks[..compressed_prefix]` are compressed (blocks are
    /// compressed front to back and leave from the front).
    compressed_prefix: usize,
    pub(super) stats: RangeModeStats,
}

impl<B: MergeBatchBuilder> RangeModeState<B> {
    /// Every range was taken and cut; the cutter's last partial batch was
    /// posted with the last range's assembly.
    fn cut_finished(&self) -> bool {
        self.ranges_taken == self.plan.range_count()
    }

    /// The oldest range whose data may still be held: by a block built and
    /// not yet accepted by the output (held or queued), the oldest posted
    /// assembly, the cutter's open batch, or a range posted and not yet taken
    /// (`ranges_taken`). Ranges below it are released; ranges
    /// `live_floor()..next_range` are live. A range stays live until its last
    /// block has left the merge, so merging ahead while the output refuses
    /// blocks stays within the window.
    fn live_floor(&self) -> usize {
        let posted = self.posted.iter().find_map(|p| match p {
            Posted::Assemble { first_range, .. } => Some(*first_range),
            Posted::Range(_) => None,
        });
        let queued = self.out_blocks.front().map(|(r, _)| *r);
        [self.held_range, queued, posted, self.cutter.oldest_open_range()]
            .into_iter()
            .flatten()
            .fold(self.ranges_taken, usize::min)
    }

    /// Live ranges now (see [`Self::live_floor`]).
    fn live(&self) -> usize {
        self.next_range - self.live_floor()
    }

    /// Whether another range may be posted: one remains, and either the live
    /// ranges are below `window`, or nothing at all is in flight (a posted
    /// job is in flight whether unstarted, started or delivered, and a built
    /// block until the output accepts it).
    ///
    /// The second keeps range mode live when the cutter's open batch alone
    /// pins the floor (consecutive ranges too small to close a batch), which
    /// the window rule alone would deadlock. Live ranges can then exceed
    /// `window`, but only by the open batch plus at most one range in flight
    /// (one range is posted at a time, and only with nothing else in
    /// flight); the open batch's share is bounded: past the floor range,
    /// every live range lies wholly inside that batch, so together they hold
    /// fewer than `target_batch_count` records and fewer than
    /// `output_byte_limit` bytes.
    fn may_post_range(&self) -> bool {
        self.next_range < self.plan.range_count()
            && (self.next_range - self.live_floor() < self.window || self.nothing_in_flight())
    }

    /// No job posted and not taken, and no block waiting for the output. A
    /// waiting block needs no job to leave (the output's consumer takes it),
    /// so it must not open the window's exception: merging ahead behind it
    /// would queue blocks without bound.
    fn nothing_in_flight(&self) -> bool {
        self.posted.is_empty() && self.out_blocks.is_empty() && self.held_range.is_none()
    }

    /// Blocks waiting for the output (queued, not held) that `is_raw`
    /// calls uncompressed, for tests.
    #[cfg(test)]
    pub(super) fn raw_waiting_for_test(&self, is_raw: impl Fn(&B::Item) -> bool) -> usize {
        self.out_blocks.iter().filter(|(_, b)| is_raw(b)).count()
    }

    /// `(live ranges, window, ranges posted, ranges)`, for tests.
    #[cfg(test)]
    pub(super) fn window_state_for_test(&self) -> (usize, usize, usize, usize) {
        (self.live(), self.window, self.next_range, self.plan.range_count())
    }

    /// Post the next range: count it as posted, charge its estimated reads to
    /// `inbox` until its merged output lands, and return its job.
    fn next_range_job<O: MergeOutput<Builder = B>>(
        &mut self,
        inbox: &OrderedInbox<MergeResult<O>>,
    ) -> MergeJob<O> {
        let (j, job) = (self.next_range, self.next_job);
        self.next_range += 1;
        self.next_job += 1;
        self.posted.push_back(Posted::Range(j));
        inbox.reserve(job, self.plan.estimated_bytes(j));
        #[cfg(test)]
        {
            let live = self.next_range - self.live_floor();
            self.stats.max_live_ranges = self.stats.max_live_ranges.max(live);
        }
        let timing = self.stats.diag.as_ref().map(|d| JobTiming::posted(&d.jobs));
        MergeJob::Range {
            task: Arc::clone(&self.plan).task(j),
            range: j,
            job,
            read: false,
            timing,
            buffers: Some(Arc::clone(&self.buffers)),
        }
    }

    /// Post the assembly of `batches` (those one range's cut closed): count
    /// it as posted, charge their framed bytes to `inbox` until the blocks
    /// land, and return its job. A first batch that continues the lead
    /// builder takes it along.
    fn assemble_job<O: MergeOutput<Builder = B>>(
        &mut self,
        batches: Vec<BatchPlan>,
        inbox: &OrderedInbox<MergeResult<O>>,
        target_batch_count: usize,
    ) -> MergeJob<O> {
        let lead = if batches.first().is_some_and(|b| b.lead) { self.lead.take() } else { None };
        let job = self.next_job;
        self.next_job += 1;
        let first_range = batches
            .first()
            .and_then(|b| b.segments.first())
            .map_or(self.ranges_taken, |s| s.range_idx);
        let records = batches.iter().map(|b| b.records).sum();
        let framed_bytes = batches.iter().map(|b| b.framed_bytes).sum();
        self.posted.push_back(Posted::Assemble { first_range, records });
        inbox.reserve(job, framed_bytes);
        let timing = self.stats.diag.as_ref().map(|d| JobTiming::posted(&d.jobs));
        MergeJob::Assemble {
            batches,
            records,
            framed_bytes,
            job,
            target_batch_count,
            lead,
            timing,
            recycle: Some(Arc::clone(&self.buffers)),
            _output: std::marker::PhantomData,
        }
    }
}

impl<B: MergeBatchBuilder> RangeModeState<B> {
    /// Post blocks that are already built (a fused range's) at the next
    /// ordinal, delivered at once: they leave behind every job posted
    /// before them. `range` is the fused range (the floor it holds).
    fn post_ready<O: MergeOutput<Builder = B>>(
        &mut self,
        range: usize,
        records: u64,
        blocks: Vec<O::Item>,
        inbox: &OrderedInbox<MergeResult<O>>,
    ) {
        let job = self.next_job;
        self.next_job += 1;
        self.posted.push_back(Posted::Assemble { first_range: range, records });
        let _ = inbox.deliver(job, MergeResult::Blocks(blocks));
    }

    /// Start a fused range (see [`begin_fused`]).
    fn begin_fused(&mut self, bytes_cap: usize, target_batch_count: usize) -> io::Result<()> {
        begin_fused(&mut self.cutter, &mut self.lead, bytes_cap, target_batch_count, &self.buffers)
    }
}

/// Start a fused range: copy the open batch's segments (earlier ranges'
/// records, at most one batch) into the lead builder (created for the open
/// batch's ordinal when there is none), recycling into `buffers` each range
/// whose last segment this drops.
pub(super) fn begin_fused<B: MergeBatchBuilder>(
    cutter: &mut BatchCutter,
    lead: &mut Option<B>,
    bytes_cap: usize,
    target_batch_count: usize,
    buffers: &RangeBufPool,
) -> io::Result<()> {
    let segments = cutter.take_open_segments();
    if segments.is_empty() {
        return Ok(());
    }
    let ordinal = cutter.next_ordinal();
    let lead = lead.get_or_insert_with(|| B::with_capacity(ordinal, bytes_cap, target_batch_count));
    for segment in segments {
        lead.push_framed_records(&segment.range.framed[segment.start..segment.end])?;
        if let Some(out) = Arc::into_inner(segment.range) {
            buffers.put(out);
        }
    }
    Ok(())
}

impl<B: MergeBatchBuilder> RangeModeState<B> {
    /// Take the result delivered at `ordinal`, which must match the job
    /// posted there: cut a merged range, take a fused one, or return an
    /// assembly's first range, records and blocks for the output.
    ///
    /// # Errors
    ///
    /// The result does not match the job posted at its ordinal.
    #[allow(clippy::type_complexity, reason = "a private helper's return")]
    fn take_result<O: MergeOutput<Builder = B>>(
        &mut self,
        ordinal: u64,
        result: MergeResult<O>,
        inbox: &OrderedInbox<MergeResult<O>>,
        table: &ClaimTable<O>,
        target_batch_count: usize,
    ) -> io::Result<Option<(usize, u64, Vec<O::Item>)>> {
        match (self.posted.pop_front(), result) {
            (Some(Posted::Range(j)), MergeResult::Range { range, out }) if j == range => {
                self.cut_merged_range(range, out, inbox, table, target_batch_count);
                Ok(None)
            }
            (Some(Posted::Range(j)), MergeResult::Fused { range, stats, blocks, records })
                if j == range =>
            {
                self.take_fused(range, &stats, records, blocks, inbox, table, target_batch_count);
                Ok(None)
            }
            (Some(Posted::Assemble { first_range, records }), MergeResult::Blocks(items)) => {
                Ok(Some((first_range, records, items)))
            }
            (posted, _) => Err(io::Error::other(format!(
                "SortMerge: range-mode result for job {ordinal} does not match the job posted \
                 there ({posted:?})"
            ))),
        }
    }

    /// Take merged range `range`: cut it, and post one assembly job for the
    /// batches the cut closed (with the cutter's last partial batch after the
    /// last range).
    fn cut_merged_range<O: MergeOutput<Builder = B>>(
        &mut self,
        range: usize,
        out: Arc<fgumi_sort::RangeOut>,
        inbox: &OrderedInbox<MergeResult<O>>,
        table: &ClaimTable<O>,
        target_batch_count: usize,
    ) {
        self.stats.add(&out.stats);
        self.ranges_taken += 1;
        let mut closed = VecDeque::new();
        self.cutter.cut(range, &out, &mut closed);
        // Drop this handle before posting the assembly, so the range's last
        // segment (and with it its buffers) is dropped by an assembly job,
        // which recycles them.
        drop(out);
        if self.cut_finished() {
            closed.extend(self.cutter.finish());
        }
        if !closed.is_empty() {
            table.post(self.assemble_job(closed.into(), inbox, target_batch_count));
            self.tickets_owed += 1;
        }
    }

    /// Take fused range `range`: post the blocks it closed and did not
    /// stream at a fresh ordinal (behind every assembly posted before; none
    /// when every block was streamed, `records` then 0), and after the last
    /// range post the cutter's last partial batch.
    #[allow(clippy::too_many_arguments, reason = "one call site; the parts of one result")]
    fn take_fused<O: MergeOutput<Builder = B>>(
        &mut self,
        range: usize,
        stats: &RangeStats,
        records: u64,
        blocks: Vec<O::Item>,
        inbox: &OrderedInbox<MergeResult<O>>,
        table: &ClaimTable<O>,
        target_batch_count: usize,
    ) {
        self.stats.add(stats);
        self.ranges_taken += 1;
        if !blocks.is_empty() {
            self.post_ready(range, records, blocks, inbox);
        }
        if self.cut_finished()
            && let Some(last) = self.cutter.finish()
        {
            table.post(self.assemble_job(vec![last], inbox, target_batch_count));
            self.tickets_owed += 1;
        }
    }

    /// Run one slice (of `budget`) of the merge thread's own `job`: a fused
    /// range's slice into the output builders, or any other job's slice,
    /// timed for the diagnostics. Returns the step and the records of the
    /// fused batches it streamed to `out_blocks` (see [`Self::stream_fused`]).
    ///
    /// # Errors
    ///
    /// The job's own error, wrapped with its ordinal.
    fn run_own_slice<O: MergeOutput<Builder = B>>(
        &mut self,
        job: &mut MergeJob<O>,
        budget: SliceBudget,
        bytes_cap: usize,
        target_batch_count: usize,
    ) -> io::Result<(JobStep<MergeResult<O>>, u64)>
    where
        B: MergeBatchBuilder<Item = O::Item>,
    {
        let ordinal = job.ordinal();
        let is_range = matches!(job, MergeJob::Range { .. });
        let own_range = if let MergeJob::Range { range, .. } = job { *range } else { 0 };
        let slice_start = Instant::now();
        let mut streamed = 0;
        let step = if self.own_fused {
            let mut sink = FusedSink {
                cutter: &mut self.cutter,
                lead: &mut self.lead,
                closed: &mut self.fused_closed,
                bytes_cap,
                target_batch_count,
            };
            let done =
                job.run_fused_slice(budget, &mut sink).map_err(|e| JobError::wrap(ordinal, e))?;
            streamed = self.stream_fused(own_range);
            match done {
                None => JobStep::Continue,
                Some((range, stats)) => {
                    self.own_fused = false;
                    self.stats.fused_ranges += 1;
                    self.stats.fused_batches += self.fused_closed.len();
                    let records = self.fused_closed.iter().map(|b| b.len() as u64).sum();
                    let blocks =
                        self.fused_closed.drain(..).map(MergeBatchBuilder::build).collect();
                    JobStep::Done(MergeResult::Fused { range, stats, blocks, records })
                }
            }
        } else {
            job.run_slice(budget).map_err(|e| JobError::wrap(ordinal, e))?
        };
        if let Some(d) = self.stats.diag.as_mut() {
            let ns = RangeDiag::ns(slice_start.elapsed());
            if is_range {
                d.own_range_ns += ns;
            } else {
                d.own_assemble_ns += ns;
            }
        }
        Ok((step, streamed))
    }

    /// Stream the batches fused range `range` has closed so far to
    /// `out_blocks` when every earlier batch is already there (no assembly
    /// is outstanding: every earlier range is taken, so only an assembly can
    /// hold earlier batches), counting them emitted. Returns their records;
    /// with an assembly outstanding they wait for the range's result.
    fn stream_fused(&mut self, range: usize) -> u64 {
        if self.fused_closed.is_empty()
            || self.posted.iter().any(|p| matches!(p, Posted::Assemble { .. }))
        {
            return 0;
        }
        let records: u64 = self.fused_closed.iter().map(|b| b.len() as u64).sum();
        self.stats.fused_batches += self.fused_closed.len();
        self.stats.streamed_batches += self.fused_closed.len();
        self.out_blocks.extend(self.fused_closed.drain(..).map(|b| (range, b.build())));
        self.emitted += records;
        records
    }

    /// Retry the block held on the output branch, then push the queued
    /// blocks in order; a refused block is held (with its range) and ends the
    /// flush. Stores the output backlog in `table` (the held block plus those
    /// queued; the gate's signal without a depth counter) and counts the
    /// accepted blocks in compression's queue depth. Returns whether any block was
    /// pushed and whether a block was newly refused.
    fn flush_output<O: MergeOutput<Builder = B>>(
        &mut self,
        held: &mut HeldSlot<Unpushed<B::Item>>,
        retry: impl Fn(Unpushed<B::Item>) -> Result<(), Unpushed<B::Item>>,
        push: impl Fn(B::Item) -> Result<(), Unpushed<B::Item>>,
        table: &ClaimTable<O>,
    ) -> (bool, bool)
    where
        B::Item: Send + HeapSize + 'static,
    {
        let (mut pushed, mut refused) = (false, false);
        let before = self.stats.blocks_out;
        if let Some(unpushed) = held.take() {
            match retry(unpushed) {
                Ok(()) => {
                    pushed = true;
                    self.held_range = None;
                    self.stats.blocks_out += 1;
                }
                Err(again) => held.put(again),
            }
        }
        while !held.is_held()
            && let Some((range, block)) = self.out_blocks.pop_front()
        {
            self.compressed_prefix = self.compressed_prefix.saturating_sub(1);
            match push(block) {
                Ok(()) => {
                    pushed = true;
                    self.stats.blocks_out += 1;
                }
                Err(unpushed) => {
                    held.put(unpushed);
                    self.held_range = Some(range);
                    refused = true;
                }
            }
        }
        // Counted once accepted: a compressor may pop a block before this
        // counts it, which only reads the depth low for that moment.
        table.note_output_pushed(self.stats.blocks_out - before);
        table.set_output_backlog(usize::from(held.is_held()) + self.out_blocks.len());
        (pushed, refused)
    }

    /// While a built block waits for the output (held, or queued behind the
    /// held one) and this merge compresses its output (with `c`), compress
    /// one waiting block ([`Self::compress_one`]).
    ///
    /// # Errors
    ///
    /// The compressor's own error.
    fn compress_waiting<O: MergeOutput<Builder = B>>(
        &mut self,
        held: &mut HeldSlot<Unpushed<B::Item>>,
        c: Option<&mut OutputCompressor>,
    ) -> io::Result<bool>
    where
        B: MergeBatchBuilder<Item = O::Item>,
    {
        match c {
            Some(c) if held.is_held() || !self.out_blocks.is_empty() => {
                self.compress_one::<O>(held, c)
            }
            _ => Ok(false),
        }
    }

    /// Compress the oldest output block not yet compressed — the held one,
    /// then the queued ones in order — with `c`, as the sink's compressor
    /// would (the output passes it through). Returns whether one was
    /// compressed; `false` when every waiting block already is, or the
    /// output's items cannot be compressed (only a
    /// [`CompressibleBlockOutput`](super::CompressibleBlockOutput) can).
    ///
    /// # Errors
    ///
    /// The compressor's own error.
    fn compress_one<O: MergeOutput<Builder = B>>(
        &mut self,
        held: &mut HeldSlot<Unpushed<B::Item>>,
        c: &mut OutputCompressor,
    ) -> io::Result<bool>
    where
        B: MergeBatchBuilder<Item = O::Item>,
    {
        let start = Instant::now();
        let mut done = false;
        if let Some(mut unpushed) = held.take() {
            let compressed = O::compress_in_place(unpushed.item_mut(), c);
            held.put(unpushed);
            done = compressed?;
        }
        while !done && self.compressed_prefix < self.out_blocks.len() {
            done = O::compress_in_place(&mut self.out_blocks[self.compressed_prefix].1, c)?;
            self.compressed_prefix += 1;
        }
        if done {
            self.stats.compressed_blocks += 1;
            if let Some(d) = self.stats.diag.as_mut() {
                d.compress_ns += RangeDiag::ns(start.elapsed());
            }
        }
        Ok(done)
    }
}

/// Where a fused range's records go: into the open batch's lead builder,
/// counted by the cutter, which says when a batch closes (the serial
/// merge's `batch_full` rule).
pub(super) struct FusedSink<'a, B> {
    pub(super) cutter: &'a mut BatchCutter,
    pub(super) lead: &'a mut Option<B>,
    pub(super) closed: &'a mut Vec<B>,
    pub(super) bytes_cap: usize,
    pub(super) target_batch_count: usize,
}

impl<B: MergeBatchBuilder> RecordSink for FusedSink<'_, B> {
    fn push_record(&mut self, body: &[u8]) -> io::Result<()> {
        let len = u32::try_from(body.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("record exceeds u32 BAM block_size: {}", body.len()),
            )
        })?;
        let ordinal = self.cutter.next_ordinal();
        let (bytes_cap, target) = (self.bytes_cap, self.target_batch_count);
        self.lead
            .get_or_insert_with(|| B::with_capacity(ordinal, bytes_cap, target))
            .push_record_bytes(body)?;
        if self.cutter.push_fused(len).is_some() {
            self.closed.push(self.lead.take().expect("the lead holds the closed batch"));
        }
        Ok(())
    }
}

/// Push the tickets owed to the job edge, after retrying a held one first.
/// Tickets are interchangeable (each claims whatever is oldest), so owed
/// tickets beyond the jobs still unstarted are dropped: every unstarted job
/// keeps at least one ticket on the edge, held or owed. Returns whether any
/// ticket was pushed.
fn push_tickets<O: MergeOutput>(
    held: &mut HeldSlot<Unpushed<MergeTicket<O>>>,
    owed: &mut usize,
    table: &Arc<ClaimTable<O>>,
    push: impl Fn(MergeTicket<O>) -> Result<(), Unpushed<MergeTicket<O>>>,
    retry: impl Fn(Unpushed<MergeTicket<O>>) -> Result<(), Unpushed<MergeTicket<O>>>,
) -> bool {
    let mut pushed = false;
    if let Some(unpushed) = held.take() {
        match retry(unpushed) {
            Ok(()) => pushed = true,
            Err(again) => {
                held.put(again);
                return false;
            }
        }
    }
    *owed = (*owed).min(table.len());
    while *owed > 0 {
        *owed -= 1;
        match push(MergeTicket::claim(table)) {
            Ok(()) => pushed = true,
            Err(unpushed) => {
                held.put(unpushed);
                break;
            }
        }
    }
    pushed
}

impl<O: MergeOutput> SortMerge<O> {
    /// Enter range mode over `slots`, sized by the planner's `fit`. Range mode
    /// merges spill runs only: the planner chooses it only when the final run
    /// spilled (a residual chunk makes it `Serial(ResidualInMemory)`), so
    /// `memory_chunks` must be empty.
    ///
    /// # Errors
    ///
    /// The merge cannot range-merge (a residual chunk, no pooled inbox to run
    /// its jobs, no spill slots, or `fgumi merge` inputs), or the plan fails
    /// ([`build_range_plan`]).
    pub(super) fn transition_to_range_mode(
        &mut self,
        slots: Vec<Arc<SortMergeSlot>>,
        memory_chunks: &MemoryChunksByKind,
        total_records: u64,
        fit: &RangeFit,
    ) -> io::Result<()> {
        if memory_chunks.total_len() != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "SortMerge: told to merge in range mode with {} residual chunk(s) in memory \
                     (range mode merges spill runs only)",
                    memory_chunks.total_len()
                ),
            ));
        }
        if self.inbox.is_none() || slots.is_empty() || self.merge_inputs.is_some() {
            return Err(io::Error::other(format!(
                "SortMerge: told to range-merge but cannot (pooled inbox: {}, spill slots: {}, \
                 merge inputs: {})",
                self.inbox.is_some(),
                slots.len(),
                self.merge_inputs.is_some(),
            )));
        }
        // The window and target the planner fitted to the budget
        // (`decide_merge_mode`).
        let sizing = RangeSizing::new(fit.target_range_bytes);
        #[cfg(test)]
        let sizing = self.range_sizing_override.unwrap_or(sizing);
        let kind = slots_key_kind(&slots)?;
        // Return the retired phase-1 arena's pages before the plan's samples
        // and the ranges' buffers are allocated. This step shares the sort's
        // coordination thread with `ReadBlocks`, which allocated the arena;
        // when a pool worker dropped the final chunk's last reference, mimalloc
        // frees the arena only at this thread's next collect, so
        // `ArenaPool::release`'s own collect (on the worker) cannot return it.
        // Once per sort.
        fgumi_sort::force_mi_collect();
        let plan = build_range_plan(self.sort_order, slots, kind, total_records, sizing)?;
        let cutter = BatchCutter::new(
            self.target_batch_count,
            self.output_byte_limit,
            O::Builder::FRAME_OVERHEAD_PER_RECORD,
        );
        // Who may claim posted jobs: a hidden A/B switch, not a tuning option.
        let policy = self.claim_policy;
        if policy != super::claim::ClaimPolicy::Auto {
            log::info!("FGUMI_SORT_RANGE_CLAIM={}: range-mode claims restricted", policy.name());
        }
        let table = ClaimTable::new(policy);
        let table = match &self.output_depth {
            Some(depth) => table.with_output_depth(Arc::clone(depth), self.fanout_width),
            None => table,
        };
        self.claim_table = Some(Arc::new(table));
        self.state = SortMergeState::RangeMode(Box::new(RangeModeState {
            plan,
            window: fit.window,
            target_range_bytes: sizing.target_range_bytes,
            cutter,
            next_range: 0,
            ranges_taken: 0,
            posted: VecDeque::new(),
            next_job: 0,
            tickets_owed: 0,
            emitted: 0,
            // At most one buffer per range the window keeps live, and none
            // larger than the merged output `range_fit` charges a live range
            // (`T + B`): a buffer from a range that outgrew its plan would
            // otherwise keep that excess resident for every later range.
            buffers: Arc::new(RangeBufPool::new(
                fit.window,
                usize::try_from(sizing.target_range_bytes.saturating_add(fit.boundary_bytes))
                    .unwrap_or(usize::MAX),
            )),
            lead: None,
            fused_closed: Vec::new(),
            own_fused: false,
            out_blocks: VecDeque::new(),
            held_range: None,
            compressed_prefix: 0,
            stats: RangeModeStats {
                diag: self.sort_stats.then(|| RangeDiag::new(fit.window)),
                ..RangeModeStats::default()
            },
        }));
        Ok(())
    }

    /// Log the range plan just entered: a warning when the plan as cut decodes
    /// more than the fit planned (index entries above `stride x F` raw bytes,
    /// from records longer than one frame, which the fit's `F` did not see),
    /// and, under `--sort-stats`, the `Sort merge mode: range` line.
    ///
    /// # Panics
    ///
    /// If the merge is not in range mode.
    pub(super) fn log_range_mode(
        &self,
        budget: &RangeBudget,
        inputs: &FitInputs,
        fit: &RangeFit,
        index_bytes: u64,
    ) {
        let SortMergeState::RangeMode(rm) = &self.state else {
            unreachable!("log_range_mode called outside RangeMode")
        };
        let plan_amplification = rm.plan.planned_decode_amplification();
        // A small slack keeps float rounding from warning on an exact match.
        if plan_amplification > fit.amplification() + 1e-9 {
            log::warn!(
                "Range plan decodes {plan_amplification:.3}x its run bytes, above the \
                 fit's planned {:.3}x (index entries larger than the fit's frame size); range \
                 memory may exceed its estimate",
                fit.amplification(),
            );
        }
        if self.sort_stats {
            // `planned decode amplification` is the fit's (`1 + kNF / T`);
            // `plan amplification` is what the cut plan's frame spans decode.
            log::info!(
                "Sort merge mode: range ({} ranges, window {}, {} workers, ~{} MiB per range, \
                 index stride {}, index {} MiB, mem {} MiB of {budget}, planned decode \
                 amplification {:.3}, plan amplification {plan_amplification:.3}, rho {:.2}, \
                 final run spilled)",
                rm.plan.range_count(),
                rm.window,
                self.fanout_width,
                rm.target_range_bytes >> 20,
                fit.stride,
                index_bytes >> 20,
                fit.mem >> 20,
                fit.amplification(),
                inputs.rho(),
            );
        }
    }

    /// One range-mode dispatch (see the module docs): flush the output
    /// (store its refused-block backlog) and push owed tickets, take ready
    /// results in post order (cut a merged range and post
    /// its assembly, queue assembled blocks for output), finish when
    /// everything is emitted, post ranges while the window allows, claim a
    /// job for this thread when it has none (an assembly first), and run one
    /// short slice of it. A refused output block does not stop the dispatch:
    /// this thread merges ahead within the window meanwhile, and with nothing
    /// left to merge it compresses the waiting blocks itself.
    ///
    /// # Errors
    ///
    /// A result that does not match the job posted at its ordinal, a job's
    /// own error, or the completion checks of [`Self::finish_range_mode`].
    pub(super) fn emit_range_mode(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        records_out: &mut u64,
    ) -> io::Result<StepOutcome> {
        // The dispatch is timed only with `--sort-stats` diagnostics on.
        let mut start = None;
        if let SortMergeState::RangeMode(rm) = &mut self.state
            && let Some(d) = rm.stats.diag.as_mut()
        {
            let now = Instant::now();
            d.begin(now);
            start = Some(now);
        }
        let outcome = self.emit_range_mode_inner(ctx, records_out)?;
        if let Some(start) = start
            && let SortMergeState::RangeMode(rm) = &self.state
            && rm.stats.diag.is_some()
        {
            let idle = matches!(outcome, StepOutcome::NoProgress | StepOutcome::Contention)
                .then(|| self.range_idle_reason());
            let backlog = self.held.is_held();
            let depth_full =
                self.claim_table.as_ref().and_then(|t| t.output_depth_full()).unwrap_or(false);
            let SortMergeState::RangeMode(rm) = &mut self.state else { unreachable!() };
            let live = rm.live();
            if let Some(d) = rm.stats.diag.as_mut() {
                d.end(start, DispatchEnd { idle, live, backlog, depth_full });
            }
        }
        Ok(outcome)
    }

    /// Why a range-mode dispatch that made no progress stopped.
    fn range_idle_reason(&self) -> Idle {
        let SortMergeState::RangeMode(rm) = &self.state else {
            unreachable!("range_idle_reason called outside RangeMode")
        };
        // The next ordinal in order is the oldest posted job's.
        let head = rm.next_job - rm.posted.len() as u64;
        let head_unstarted = (!rm.posted.is_empty()).then(|| {
            self.claim_table
                .as_ref()
                .and_then(|t| t.front_ordinal())
                .is_some_and(|front| front <= head)
        });
        Idle::of(IdleState {
            output_held: self.held.is_held(),
            ranges_left: rm.next_range < rm.plan.range_count(),
            head_unstarted,
            ticket_held: self.held_job.is_held(),
        })
    }

    fn emit_range_mode_inner(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        records_out: &mut u64,
    ) -> io::Result<StepOutcome> {
        let inbox = Arc::clone(self.inbox.as_ref().expect("range mode runs only with an inbox"));
        let table =
            Arc::clone(self.claim_table.as_ref().expect("range mode runs with a claim table"));
        let target_batch_count = self.target_batch_count;
        // A fused range's builders start as the serial merge's do.
        let bytes_cap = usize::try_from(self.output_byte_limit).unwrap_or(usize::MAX);
        let bytes_cap = super::INITIAL_OUTPUT_BUFFER_BYTES.min(bytes_cap);
        let SortMergeState::RangeMode(rm) = &mut self.state else {
            unreachable!("emit_range_mode called outside RangeMode")
        };
        let outputs = ctx.outputs.view();
        // While the gate bars workers, tickets stay owed: one pushed now would
        // only be refused (and this thread claims the jobs meanwhile).
        let tickets = |held: &mut HeldSlot<Unpushed<MergeTicket<O>>>, owed: &mut usize| {
            !table.workers_gated()
                && push_tickets(held, owed, &table, |t| outputs.b.push(t), |u| outputs.b.retry(u))
        };
        let mut flush = |rm: &mut RangeModeState<O::Builder>, held: &mut HeldSlot<_>| {
            let (pushed, refused) =
                rm.flush_output(held, |u| outputs.a.retry(u), |b| outputs.a.push(b), &table);
            self.dbg.output_full += u64::from(refused);
            pushed || refused
        };

        // 1. Flush the output: the held block, then the queued ones (this
        //    stores the backlog). Then re-owe tickets the gate refused and
        //    push the owed ones (none while the gate is closed).
        let mut progressed = flush(rm, &mut self.held);
        rm.tickets_owed += table.take_refused();
        progressed |= tickets(&mut self.held_job, &mut rm.tickets_owed);

        // 2. Take ready results in post order, whether or not the output is
        //    held: taking keeps `ranges_taken` advancing, so this thread's
        //    next range can still fuse. Their blocks queue behind the held
        //    one, each keeping its range live.
        let mut took = false;
        while let Some((ordinal, result)) = inbox.try_take_next() {
            took = true;
            if let Some((first_range, records, items)) =
                rm.take_result(ordinal, result, &inbox, &table, target_batch_count)?
            {
                rm.emitted += records;
                *records_out += records;
                rm.out_blocks.extend(items.into_iter().map(|b| (first_range, b)));
            }
        }
        if took {
            flush(rm, &mut self.held);
        }

        // 3. Every range cut and every batch emitted and accepted: done.
        if rm.cut_finished()
            && rm.posted.is_empty()
            && rm.out_blocks.is_empty()
            && !self.held.is_held()
        {
            return self.finish_range_mode();
        }

        // 4. Post ranges while the window allows, whatever the ticket edge
        //    does: a posted job without a ticket on the edge is still reached,
        //    by this thread's claim or by a worker once the edge drains
        //    (tickets are interchangeable). Posting before claiming keeps
        //    tickets ahead of the workers while this thread is in its slice.
        while rm.may_post_range() {
            table.post(rm.next_range_job(&inbox));
            rm.tickets_owed += 1;
            progressed = true;
        }
        progressed |= tickets(&mut self.held_job, &mut rm.tickets_owed);

        // 5. Claim this thread's own job when it has none: the oldest
        //    unstarted assembly (it produces output and releases the window),
        //    else the oldest unstarted range, the one this thread would
        //    otherwise wait on soonest.
        if self.own_range.is_none() {
            self.own_range = table.claim_prefer_assembly();
            // Every range before it is taken and cut: the cutter's carry is
            // known, so merge it straight into output batches.
            if let Some(MergeJob::Range { range, .. }) = &self.own_range
                && *range == rm.ranges_taken
            {
                rm.begin_fused(bytes_cap, target_batch_count)?;
                rm.own_fused = true;
            }
        }

        // 6. Run one short slice of this thread's own job (a range slice, or
        //    a whole assembly). Its result is delivered to the inbox at its
        //    ordinal like any worker's; a fused range's closed batches may go
        //    straight to the output.
        if let Some(mut job) = self.own_range.take() {
            let ordinal = job.ordinal();
            let budget = own_slice(target_batch_count);
            let (step, streamed) =
                rm.run_own_slice(&mut job, budget, bytes_cap, target_batch_count)?;
            if streamed > 0 {
                *records_out += streamed;
                flush(rm, &mut self.held);
            }
            match step {
                JobStep::Continue => self.own_range = Some(job),
                JobStep::Done(out) => {
                    match &out {
                        MergeResult::Range { .. } | MergeResult::Fused { .. } => {
                            rm.stats.own_ranges += 1;
                        }
                        MergeResult::Blocks(items) => rm.stats.own_assemblies += items.len(),
                        MergeResult::Block(_) => {}
                    }
                    let _ = inbox.deliver(ordinal, out);
                }
                JobStep::Spent => unreachable!("a merge job always has work to run"),
            }
            return Ok(StepOutcome::Progress);
        }
        // 7. No job of its own (the window is full or every job is claimed)
        //    and the output refuses a block: compress the oldest uncompressed
        //    waiting block, as the sink's compressor would (it passes the
        //    block through), then retry the output.
        if rm.compress_waiting::<O>(&mut self.held, self.output_compressor.as_mut())? {
            flush(rm, &mut self.held);
            return Ok(StepOutcome::Progress);
        }
        if took || progressed {
            return Ok(StepOutcome::Progress);
        }
        // 8. Nothing taken, pushed, posted, claimed, run or compressed: the next result in
        //    order is running on a worker (or barred to this thread by the
        //    claim policy), or the window is full behind a refused output
        //    block. Re-record the wait on the inbox; a held block's consumer
        //    pop wakes this thread too, so it parks on `Contention` as the
        //    serial merge does on a refused batch.
        Ok(if inbox.await_next() {
            StepOutcome::Progress
        } else if self.held.is_held() {
            StepOutcome::Contention
        } else {
            StepOutcome::NoProgress
        })
    }

    /// The `--sort-stats` lines of a finished range merge.
    fn log_range_stats(&self, rm: &RangeModeState<O::Builder>, table: &ClaimTable<O>) {
        let t = &rm.stats.totals;
        #[allow(clippy::cast_precision_loss, reason = "a diagnostic ratio")]
        let (largest_pct, amplification) = (
            if t.records == 0 {
                0.0
            } else {
                rm.stats.largest_range_records as f64 * 100.0 / t.records as f64
            },
            if t.kept_run_bytes == 0 {
                0.0
            } else {
                t.decoded_bytes as f64 / t.kept_run_bytes as f64
            },
        );
        log::info!(
            "Range merge: {} ranges ({} on the merge thread, which also assembled {} \
             batches), {} records, largest range {largest_pct:.3}% of records, decode \
             amplification {amplification:.3}, {} MiB read",
            rm.plan.range_count(),
            rm.stats.own_ranges,
            rm.stats.own_assemblies,
            t.records,
            t.read_bytes >> 20,
        );
        if let Some(d) = rm.stats.diag.as_ref() {
            d.log(rm.stats.own_ranges, rm.stats.own_assemblies, rm.stats.compressed_blocks);
            d.log_depth(table.depth_width(), table.gate_enabled(), table.refused_total());
        }
        table.log();
        log::info!(
            "Range fused: {} of the merge thread's {} ranges merged straight into {} output \
             batches ({} streamed to the output during their range)",
            rm.stats.fused_ranges,
            rm.stats.own_ranges,
            rm.stats.fused_batches,
            rm.stats.streamed_batches,
        );
        rm.buffers.log(rm.plan.range_count());
        self.log_demand_lines();
        self.log_supply_lines();
        #[allow(clippy::cast_precision_loss, reason = "a diagnostic share")]
        let pct = if rm.stats.blocks_out == 0 {
            0.0
        } else {
            rm.stats.compressed_blocks as f64 * 100.0 / rm.stats.blocks_out as f64
        };
        log::info!(
            "Range output: merge thread compressed {} of {} blocks ({pct:.1}%; merge-thread \
             compression {})",
            rm.stats.compressed_blocks,
            rm.stats.blocks_out,
            if self.output_compressor.is_some() { "on" } else { "not configured" },
        );
    }

    /// Complete range mode: check every ingested record was emitted and no
    /// job was left unclaimed, log the `--sort-stats` lines, fill the stats
    /// slot, and move to `Done`.
    ///
    /// # Errors
    ///
    /// The emitted record count differs from the records ingested or from
    /// the plan's, or a posted job was never claimed.
    fn finish_range_mode(&mut self) -> io::Result<StepOutcome> {
        let SortMergeState::RangeMode(rm) =
            std::mem::replace(&mut self.state, SortMergeState::Done)
        else {
            unreachable!("finish_range_mode called outside RangeMode")
        };
        let table = self.claim_table.take().expect("range mode runs with a claim table");
        if table.len() != 0
            || rm.lead.is_some()
            || !rm.fused_closed.is_empty()
            || table.output_backlog() != 0
        {
            return Err(io::Error::other(format!(
                "SortMerge: range merge finished with {} posted job(s) never claimed, an open \
                 lead batch: {}, {} fused batch(es) not emitted, output backlog {}",
                table.len(),
                rm.lead.is_some(),
                rm.fused_closed.len(),
                table.output_backlog(),
            )));
        }
        ensure_all_records_merged(self.processed, rm.emitted)?;
        if rm.emitted != rm.plan.total_records() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "SortMerge: range merge emitted {} records but its plan holds {}",
                    rm.emitted,
                    rm.plan.total_records()
                ),
            ));
        }
        log::info!("Sort merge complete: {} records merged", rm.emitted);
        if self.sort_stats {
            self.log_range_stats(&rm, &table);
        }
        self.publish_merge_stats(rm.emitted);
        Ok(StepOutcome::Progress)
    }
}

#[cfg(test)]
mod idle_tests {
    use rstest::rstest;

    use super::{Idle, IdleState};

    /// The idle reason of a dispatch that made no progress.
    /// Columns: output held, ranges left to post, the next result (`None`:
    /// nothing posted; `Some(unstarted)`), a ticket held, expected reason.
    #[rstest]
    #[case::output_held_window_full(true, true, Some(false), false, Idle::OutputHeld)]
    #[case::output_held_draining(true, false, Some(false), false, Idle::WaitingOnWorker)]
    #[case::running_on_a_worker(false, true, Some(false), true, Idle::WaitingOnWorker)]
    #[case::unstarted_in_the_table(false, false, Some(true), false, Idle::WaitingOnNothing)]
    #[case::ticket_edge(false, true, None, true, Idle::TicketEdge)]
    #[case::window_full(false, true, None, false, Idle::WindowFull)]
    #[case::draining(false, false, None, false, Idle::Draining)]
    fn idle_reason(
        #[case] output_held: bool,
        #[case] ranges_left: bool,
        #[case] head_unstarted: Option<bool>,
        #[case] ticket_held: bool,
        #[case] expected: Idle,
    ) {
        let state = IdleState { output_held, ranges_left, head_unstarted, ticket_held };
        assert_eq!(Idle::of(state), expected);
    }
}
