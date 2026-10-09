//! `PreadSlices` — the pool step that performs the sort's positional reads.
//!
//! Shared by both read paths, so it belongs to neither: `source` plans and
//! frames the input BAM's reads, `sort` plans and parses the spill files'.
//!
//! A serial planner (`PlanInputReads` for the input BAM, `SpillReadPlanner` for
//! the merge's spill files) decides *what* to read and emits [`ReadRequest`]s;
//! this `Parallel` step performs each one as a single positional read
//! ([`read_at_exact`]) into a recycled buffer from a [`SliceBufferPool`] and
//! emits the bytes as a [`ReadSlice`] carrying a refcounted [`SliceLease`].
//! Concurrency comes from the pool: several workers each running one read is
//! what raises the device's queue depth, with no thread of its own.
//!
//! One step type, two instances, so each path has its own stats row, phase
//! classification and cap:
//!
//! - `PreadInputSlices` ([`PreadSlices::input`]): the input BAM. Its output edge
//!   is ordinal-ordered (`ByItemOrdinal`) because the input is one stream and
//!   the framer needs file order. It runs on every pool worker except the one
//!   hosting the reader-affine planner and framer
//!   ([`PoolPlacement::ExcludeReader`]), so a blocking read never stalls them.
//!   After each read it adds the slice's bytes to the [`InputLedger`]'s
//!   `landed`, half of the read-streams ratchet's signal.
//! - `PreadSpillSlices` ([`PreadSlices::spill`]): the spill files. Its output
//!   edge is unordered: spill slices belong to k independent streams, and a
//!   global reorder would make every slot's delivery wait on the slowest
//!   outstanding read of any other slot. Each slice carries its stream and
//!   per-stream `seq`, and the owning slot reassembles its own stream.
//!
//! Worker occupancy: a 1 MiB read blocks its worker for a few milliseconds on a
//! network volume, so the planners bound the number of outstanding slices
//! rather than relying on queue limits.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

use fgumi_bam_io::pread::{PositionalSource, SliceBufferPool, SliceLease, read_at_exact};
use fgumi_pipeline_core::{
    HeapSize, HeldRetry, Ordered, PhaseCap, PoolPlacement, Unpushed, admit_input,
    held::HeldSlot,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{CounterSpec, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};
use fgumi_sort::SortMergeSlot;

/// Counter slot index for slices read.
const SLICES: usize = 0;
/// Counter slot index for bytes read.
const BYTES_READ: usize = 1;

/// Which read-ahead bucket a spill read is charged to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpillClass {
    /// A spill file the merge needs next (awaited, predicted or frontier).
    Hot,
    /// Any other live spill file.
    Cold,
}

/// Where a slice's bytes go once read.
#[derive(Clone)]
pub enum ReadTarget {
    /// The input BAM's framer.
    Input,
    /// The merge slot that owns this spill stream, and the read-ahead bucket
    /// the read is charged to.
    Spill {
        /// The slot.
        slot: Arc<SortMergeSlot>,
        /// The bucket.
        class: SpillClass,
    },
}

impl std::fmt::Debug for ReadTarget {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Input => f.write_str("Input"),
            Self::Spill { slot, class } => write!(f, "Spill({}, {class:?})", slot.file_id),
        }
    }
}

/// One positional read a planner asks for.
pub struct ReadRequest {
    /// Dense per planner; the input edge is reordered on it, the spill edge is not.
    pub ordinal: u64,
    /// `0` for the input; the spill file's id for a spill.
    pub stream: u32,
    /// Per-stream slice sequence, dense from 0: a spill slot reorders on it.
    pub seq: u32,
    /// What to read from.
    pub source: Arc<dyn PositionalSource>,
    /// Byte offset of the read.
    pub offset: u64,
    /// Bytes to read; the read fails unless exactly this many arrive.
    pub len: u32,
    /// Whether the read covers the stream's final byte.
    pub last: bool,
    /// Where the bytes go.
    pub target: ReadTarget,
}

impl HeapSize for ReadRequest {
    /// The bytes this request commits the reader to holding once it lands.
    fn heap_size(&self) -> usize {
        self.len as usize
    }
}

impl Ordered for ReadRequest {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

/// The bytes of one completed [`ReadRequest`].
pub struct ReadSlice {
    /// The request's ordinal.
    pub ordinal: u64,
    /// The request's stream.
    pub stream: u32,
    /// The request's per-stream sequence.
    pub seq: u32,
    /// The request's offset.
    pub offset: u64,
    /// Exactly `len` bytes (a short read is an error, never a short slice).
    pub bytes: SliceLease,
    /// Whether the slice covers the stream's final byte.
    pub last: bool,
    /// Where the bytes go.
    pub target: ReadTarget,
}

impl HeapSize for ReadSlice {
    /// The slice's buffer capacity: a pooled buffer reused for a shorter read
    /// keeps its allocation (up to twice the length), and that allocation is
    /// what the queue holds — the slice pool's gauges charge the same.
    fn heap_size(&self) -> usize {
        self.bytes.capacity()
    }
}

impl Ordered for ReadSlice {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

/// The input path's byte ledger, the read-streams ratchet's signal: the
/// planner adds to `issued`, `PreadInputSlices` to `landed`, the framer to
/// `framed`, and the framer sets `framer_waiting` while its head slice has not
/// arrived.
#[derive(Debug, Default)]
pub struct InputLedger {
    /// Bytes requested by the planner.
    pub issued: AtomicU64,
    /// Bytes read by `PreadInputSlices`.
    pub landed: AtomicU64,
    /// Bytes framed by the framer.
    pub framed: AtomicU64,
    /// Whether the framer's last pop found no slice.
    pub framer_waiting: AtomicBool,
}

/// Number of log2 buckets in [`RequestSizeHist`]; bucket `i` holds requests in
/// `[2^(i + 10), 2^(i + 11))` bytes (the first also holds anything smaller,
/// the last anything larger).
const HIST_BUCKETS: usize = 16;
/// log2 of the first bucket's lower bound (1 KiB).
const HIST_BASE_LOG2: u32 = 10;

/// Request-size histogram and in-flight gauge for `--sort-stats`.
#[derive(Debug)]
pub struct RequestSizeHist {
    buckets: [AtomicU64; HIST_BUCKETS],
    min: AtomicU64,
    inflight: AtomicU64,
    inflight_max: AtomicU64,
    inflight_sum: AtomicU64,
    inflight_samples: AtomicU64,
}

impl Default for RequestSizeHist {
    fn default() -> Self {
        Self {
            buckets: std::array::from_fn(|_| AtomicU64::new(0)),
            min: AtomicU64::new(u64::MAX),
            inflight: AtomicU64::new(0),
            inflight_max: AtomicU64::new(0),
            inflight_sum: AtomicU64::new(0),
            inflight_samples: AtomicU64::new(0),
        }
    }
}

impl RequestSizeHist {
    /// Count one completed request of `len` bytes.
    pub fn record(&self, len: u32) {
        let log2 = 31 - len.max(1).leading_zeros();
        let i = (log2.saturating_sub(HIST_BASE_LOG2) as usize).min(HIST_BUCKETS - 1);
        self.buckets[i].fetch_add(1, Ordering::Relaxed);
        self.min.fetch_min(u64::from(len), Ordering::Relaxed);
    }

    /// A read starts: bump the gauge and sample it.
    pub fn inflight_enter(&self) {
        let now = self.inflight.fetch_add(1, Ordering::Relaxed) + 1;
        self.inflight_max.fetch_max(now, Ordering::Relaxed);
        self.inflight_sum.fetch_add(now, Ordering::Relaxed);
        self.inflight_samples.fetch_add(1, Ordering::Relaxed);
    }

    /// A read ends (successfully or not).
    pub fn inflight_exit(&self) {
        self.inflight.fetch_sub(1, Ordering::Relaxed);
    }

    /// Requests recorded.
    #[must_use]
    pub fn count(&self) -> u64 {
        self.buckets.iter().map(|b| b.load(Ordering::Relaxed)).sum()
    }

    /// `(p50, p90, min)` request size in bytes; the percentiles are the lower
    /// bound of the log2 bucket holding them. All zero when nothing was
    /// recorded.
    #[must_use]
    pub fn p50_p90_min(&self) -> (u64, u64, u64) {
        let counts: Vec<u64> = self.buckets.iter().map(|b| b.load(Ordering::Relaxed)).collect();
        let total: u64 = counts.iter().sum();
        if total == 0 {
            return (0, 0, 0);
        }
        let at = |q_num: u64, q_den: u64| -> u64 {
            let target = (total * q_num).div_ceil(q_den).max(1);
            let mut seen = 0;
            let mut bound = 1u64 << HIST_BASE_LOG2;
            for c in &counts {
                seen += c;
                if seen >= target {
                    break;
                }
                bound <<= 1;
            }
            bound
        };
        (at(1, 2), at(9, 10), self.min.load(Ordering::Relaxed))
    }

    /// `(mean, max)` reads in flight, sampled as each read starts.
    #[must_use]
    #[allow(clippy::cast_precision_loss, reason = "a diagnostic mean")]
    pub fn inflight_mean_max(&self) -> (f64, u64) {
        let n = self.inflight_samples.load(Ordering::Relaxed);
        let mean =
            if n == 0 { 0.0 } else { self.inflight_sum.load(Ordering::Relaxed) as f64 / n as f64 };
        (mean, self.inflight_max.load(Ordering::Relaxed))
    }
}

/// `Parallel` positional-read step (see the module docs).
pub struct PreadSlices {
    name: &'static str,
    ordered: bool,
    placement: PoolPlacement,
    pool: Arc<SliceBufferPool>,
    held: HeldSlot<Unpushed<ReadSlice>>,
    output_byte_limit: u64,
    /// `Some` only under `--sort-stats`.
    hist: Option<Arc<RequestSizeHist>>,
    /// `Some` only for the input instance.
    input_ledger: Option<Arc<InputLedger>>,
    /// The phase cap shared with the phase's other pool steps; `None` unless a
    /// per-phase thread flag was given.
    cap: Option<Arc<PhaseCap>>,
}

impl PreadSlices {
    /// `PreadInputSlices`: ordered output, off the reader worker, bumps
    /// `ledger.landed` per slice.
    #[must_use]
    pub fn input(
        pool: Arc<SliceBufferPool>,
        ledger: Arc<InputLedger>,
        output_byte_limit: u64,
    ) -> Self {
        Self {
            name: "PreadInputSlices",
            ordered: true,
            placement: PoolPlacement::ExcludeReader,
            pool,
            held: HeldSlot::new(),
            output_byte_limit,
            hist: None,
            input_ledger: Some(ledger),
            cap: None,
        }
    }

    /// `PreadSpillSlices`: unordered output, on every pool worker.
    #[must_use]
    pub fn spill(pool: Arc<SliceBufferPool>, output_byte_limit: u64) -> Self {
        Self {
            name: "PreadSpillSlices",
            ordered: false,
            placement: PoolPlacement::AllWorkers,
            pool,
            held: HeldSlot::new(),
            output_byte_limit,
            hist: None,
            input_ledger: None,
            cap: None,
        }
    }

    /// Share a phase admission cap with the other pool steps of this phase
    /// (`None` = uncapped).
    #[must_use]
    pub fn with_phase_cap(mut self, cap: Option<Arc<PhaseCap>>) -> Self {
        self.cap = cap;
        self
    }

    /// Record request sizes and in-flight reads into `hist`.
    #[must_use]
    pub fn with_hist(mut self, hist: Arc<RequestSizeHist>) -> Self {
        self.hist = Some(hist);
        self
    }

    fn run_once(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        cap: Option<&PhaseCap>,
    ) -> io::Result<StepOutcome> {
        if matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld) {
            return Ok(StepOutcome::Contention);
        }
        // The permit is taken before the pop (an empty poll takes none) and held
        // across the read: the read is the capped work.
        let _permit = match admit_input(ctx.input, cap) {
            Ok(permit) => permit,
            Err(outcome) => return Ok(outcome),
        };
        let Some(req) = ctx.input.pop() else {
            return Ok(if ctx.input.is_drained() {
                StepOutcome::Finished
            } else {
                StepOutcome::NoProgress
            });
        };
        let slice = self.read(req)?;
        let len = slice.bytes.len() as u64;
        if let Err(unpushed) = ctx.outputs.push(slice) {
            self.held.put(unpushed);
        }
        ctx.counters.add(SLICES, 1);
        ctx.counters.add(BYTES_READ, len);
        Ok(StepOutcome::Progress)
    }

    /// Perform one request.
    fn read(&self, req: ReadRequest) -> io::Result<ReadSlice> {
        let mut buf = self.pool.take(req.len as usize);
        if let Some(h) = &self.hist {
            h.inflight_enter();
        }
        let result = read_at_exact(&*req.source, &mut buf, req.offset);
        if let Some(h) = &self.hist {
            h.inflight_exit();
        }
        result?;
        if let Some(h) = &self.hist {
            h.record(req.len);
        }
        if let Some(l) = &self.input_ledger {
            l.landed.fetch_add(u64::from(req.len), Ordering::Relaxed);
        }
        Ok(ReadSlice {
            ordinal: req.ordinal,
            stream: req.stream,
            seq: req.seq,
            offset: req.offset,
            bytes: self.pool.lease(buf),
            last: req.last,
            target: req.target,
        })
    }
}

impl Clone for PreadSlices {
    fn clone(&self) -> Self {
        Self {
            name: self.name,
            ordered: self.ordered,
            placement: self.placement,
            pool: Arc::clone(&self.pool),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            hist: self.hist.clone(),
            input_ledger: self.input_ledger.clone(),
            cap: self.cap.clone(),
        }
    }
}

impl Step for PreadSlices {
    type Input = ReadRequest;
    type Outputs = OrderedBytesSingle<ReadSlice>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![if self.ordered {
                BranchOrdering::ByItemOrdinal
            } else {
                BranchOrdering::None
            }],
        }
    }

    fn pool_placement(&self) -> PoolPlacement {
        self.placement
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] =
            &[CounterSpec::new("slices", "slices"), CounterSpec::new("bytes_read", "bytes")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // Moved out for the call: the permit borrows the cap while the read
        // needs `&mut self`.
        let cap = self.cap.take();
        let outcome = self.run_once(ctx, cap.as_deref());
        self.cap = cap;
        outcome
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }

    fn phase_cap(&self) -> Option<&PhaseCap> {
        self.cap.as_deref()
    }
}

#[cfg(test)]
mod tests;
