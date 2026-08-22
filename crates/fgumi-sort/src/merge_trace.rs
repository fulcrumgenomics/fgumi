//! End-to-end timing for every stage a spill block passes through.
//!
//! [`crate::merge_stalls`] narrowed the merge's stall to one state: the awaited
//! file's reorder buffer is empty and a decompression is in flight, 73% of the
//! time on a measured 780M-record merge. That says *where* the consumer is
//! blocked. It does not say why the buffer was allowed to drain in the first
//! place, which is the actual defect — the pool has 45% of its capacity idle
//! and 85 other files buffered to their caps while the one file being consumed
//! runs dry.
//!
//! Answering that needs the block's whole journey, not a snapshot at the end of
//! it:
//!
//! ```text
//!   disk ──read──▶ raw FIFO ──claim──▶ decompress ──insert──▶ reorder ──pop──▶ consumer
//!         (A)        (B)        (C)        (D)         (E)      (F)      (G)
//! ```
//!
//! - **A** `read_batch` — how long a disk read of `PHASE2_READ_BATCH` blocks takes.
//! - **B** `raw_dwell` — how long compressed bytes sat in the FIFO before a
//!   worker claimed them. Pure scheduling latency: the data was there and the
//!   pool did not pick it up.
//! - **D** `decompress` — the per-block cost itself, the irreducible part.
//! - **F** `reorder_dwell` — how long a decompressed block waited to be
//!   consumed. Near zero means the lookahead buffer is not actually buffering:
//!   blocks are being consumed as fast as they appear, so the cap is not what
//!   is limiting anything.
//!
//! And the refill cycle, which is the hungry-file question stated directly.
//! When the consumer drains a file's reorder buffer to nothing, the clock
//! starts:
//!
//! - `claim_lag` — until some worker claims a raw block for that file
//! - `insert_lag` — until a block actually lands back in the buffer
//! - `read_lag` — for the subset where the raw FIFO was empty too, until a read
//!   completes
//!
//! Splitting the empties by what the file had available at the moment it
//! drained (`raw ready` / `already decompressing` / `nothing`) says which of
//! those three is on the critical path, and `claim_lag` versus `insert_lag`
//! splits the remainder into "nobody picked it up" versus "picking it up is
//! slow".
//!
//! # Everything here is a distribution, not a mean
//!
//! Merge stalls are heavy-tailed: a mean park of 83 µs is consistent with
//! almost every park costing 83 µs, and equally consistent with 90% costing 10
//! µs and 1% costing 7 ms. Those two have different causes and different fixes.
//! Every timing below is therefore a log2 histogram, and the reports quote
//! percentiles beside the mean.
//!
//! # Cost
//!
//! One `Instant::now()` pair per *block*, not per record. Blocks are ~64 KB and
//! records ~100 bytes, so this is roughly 1/600th the rate of anything on the
//! record path, and the measured merge does about 5M blocks against 780M
//! records. The per-file depth counters are plain relaxed atomics maintained
//! next to the mutexes that already guard the collections.

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use crate::merge_stalls::AwaitedState;

// ============================================================================
// Histogram
// ============================================================================

/// Number of log2 buckets. Bucket 0 is sub-microsecond; bucket `k` covers
/// `[2^(k-1), 2^k)` µs, so the top bucket opens at ~4.2 s.
pub(crate) const HIST_BUCKETS: usize = 24;

/// Bucket index for a duration in nanoseconds.
#[must_use]
pub(crate) fn bucket_of_nanos(nanos: u64) -> usize {
    let micros = nanos / 1_000;
    if micros == 0 {
        return 0;
    }
    (micros.ilog2() as usize + 1).min(HIST_BUCKETS - 1)
}

/// Inclusive lower bound of `bucket`, in microseconds.
#[must_use]
pub(crate) fn bucket_floor_micros(bucket: usize) -> u64 {
    if bucket == 0 { 0 } else { 1u64 << (bucket - 1) }
}

/// A lock-free log2 histogram of durations, plus an exact sum and count.
///
/// The sum is kept exactly rather than estimated from buckets so the mean stays
/// trustworthy; the buckets exist for the shape, which the mean cannot show.
#[derive(Debug, Default)]
pub(crate) struct DurationHistogram {
    buckets: [AtomicU64; HIST_BUCKETS],
    count: AtomicU64,
    total_nanos: AtomicU64,
}

impl DurationHistogram {
    /// Record one observation.
    pub(crate) fn record(&self, nanos: u64) {
        self.buckets[bucket_of_nanos(nanos)].fetch_add(1, Ordering::Relaxed);
        self.count.fetch_add(1, Ordering::Relaxed);
        self.total_nanos.fetch_add(nanos, Ordering::Relaxed);
    }

    /// Record one *count* observation (a queue depth, a run length) into a
    /// histogram built for durations, scaling one unit to [`BLOCKS_TO_NANOS`]
    /// so counts land in the microsecond lane the log2 bucketing is built
    /// around. Without the scale every count below `BLOCKS_TO_NANOS` buckets to
    /// zero and both the mean and the tail read as zero.
    ///
    /// Read the result back with [`HistogramReport::summary_blocks`], never
    /// [`HistogramReport::summary`] -- the latter labels these as durations.
    pub(crate) fn record_count(&self, count: u64) {
        self.record(count.saturating_mul(BLOCKS_TO_NANOS));
    }

    /// Time `f`, recording how long it took, and return its value.
    pub(crate) fn time<R>(&self, f: impl FnOnce() -> R) -> R {
        let start = Instant::now();
        let result = f();
        self.record(elapsed_nanos(start));
        result
    }

    pub(crate) fn snapshot(&self) -> HistogramReport {
        HistogramReport {
            buckets: std::array::from_fn(|i| self.buckets[i].load(Ordering::Relaxed)),
            count: self.count.load(Ordering::Relaxed),
            total_nanos: self.total_nanos.load(Ordering::Relaxed),
        }
    }
}

/// Read-only view of a [`DurationHistogram`].
#[derive(Debug, Clone, Copy, Default)]
pub struct HistogramReport {
    /// Observations per log2 bucket.
    pub buckets: [u64; HIST_BUCKETS],
    /// Total observations.
    pub count: u64,
    /// Exact sum of all observations, in nanoseconds.
    pub total_nanos: u64,
}

impl HistogramReport {
    /// Whether anything was recorded.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.count == 0
    }

    /// Exact total, in seconds.
    #[must_use]
    #[allow(clippy::cast_precision_loss, reason = "nanosecond totals stay far below 2^52")]
    pub fn total_secs(self) -> f64 {
        self.total_nanos as f64 / 1e9
    }

    /// Exact mean, in microseconds.
    #[must_use]
    #[allow(clippy::cast_precision_loss, reason = "nanosecond totals stay far below 2^52")]
    pub fn mean_micros(self) -> f64 {
        if self.count == 0 {
            return 0.0;
        }
        self.total_nanos as f64 / self.count as f64 / 1_000.0
    }

    /// Approximate percentile, in microseconds.
    ///
    /// Returns the lower bound of the bucket the percentile falls in, so it
    /// understates by up to a factor of two. That is the right trade for a
    /// distribution spanning microseconds to seconds: the question these answer
    /// is "which order of magnitude", and an exact quantile would cost either
    /// reservoir sampling or a far larger table.
    #[must_use]
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "`fraction` is a caller-supplied percentile in [0, 1] and `count` is an observation count, so the product is non-negative and bounded by `count`"
    )]
    pub fn percentile_micros(self, fraction: f64) -> u64 {
        if self.count == 0 {
            return 0;
        }
        // Walk against the bucket total rather than `count`. `snapshot` reads
        // the buckets and the counter as separate relaxed loads while workers
        // are still recording, so `count` can be ahead of the buckets by the
        // observations that landed in between. Targeting the larger number
        // means no bucket ever reaches it, the loop falls through, and every
        // percentile reports the top bucket -- a p50 of 8us printing as ~4.2s.
        let observed: u64 = self.buckets.iter().sum();
        if observed == 0 {
            return 0;
        }
        let target = (observed as f64 * fraction.clamp(0.0, 1.0)).ceil() as u64;
        let mut cumulative = 0u64;
        for (idx, &n) in self.buckets.iter().enumerate() {
            cumulative += n;
            if cumulative >= target {
                return bucket_floor_micros(idx);
            }
        }
        bucket_floor_micros(HIST_BUCKETS - 1)
    }

    /// `count`, total, mean and the p50/p90/p99 tail for a histogram whose
    /// observations are *counts*, not durations.
    ///
    /// [`ConsumerTraceStats::record_source_run`] borrows this histogram's log2
    /// bucketing for run lengths, scaling one block to [`BLOCKS_TO_NANOS`]. That
    /// puts blocks in the microsecond lane, so `mean_micros` and
    /// `percentile_micros` already read as blocks -- but [`Self::summary`] would
    /// label them `us`, and would render `total` through `total_secs()`, which is
    /// neither blocks nor seconds. This labels every field as what it is.
    #[must_use]
    pub fn summary_blocks(self) -> String {
        format!(
            "n={} total={} blocks mean={:.1} p50={} p90={} p99={}",
            self.count,
            self.total_nanos / BLOCKS_TO_NANOS,
            self.mean_micros(),
            self.percentile_micros(0.50),
            self.percentile_micros(0.90),
            self.percentile_micros(0.99),
        )
    }

    /// `count`, mean and the p50/p90/p99 tail, formatted for one log line.
    ///
    /// Durations only. For a histogram of counts use [`Self::summary_blocks`].
    #[must_use]
    pub fn summary(self) -> String {
        format!(
            "n={} total={:.1}s mean={:.1}us p50={}us p90={}us p99={}us",
            self.count,
            self.total_secs(),
            self.mean_micros(),
            self.percentile_micros(0.50),
            self.percentile_micros(0.90),
            self.percentile_micros(0.99),
        )
    }
}

/// Nanoseconds since `start`, saturating rather than wrapping.
#[must_use]
pub(crate) fn elapsed_nanos(start: Instant) -> u64 {
    u64::try_from(start.elapsed().as_nanos()).unwrap_or(u64::MAX)
}

// ============================================================================
// Block lifecycle
// ============================================================================

/// The two *dwell* stages of a spill block's journey: time it spent queued
/// rather than time it spent being worked on.
///
/// The working stages — the disk read and the decompression — are already timed
/// by [`MergePhaseCounters`](crate::merge_phases::MergePhaseCounters), which the
/// workers record into on the production path. Duplicating them here would give
/// two sources for one number, so this struct owns only what nothing else
/// measures.
#[derive(Debug, Default)]
pub(crate) struct BlockLifecycleStats {
    /// Compressed bytes sitting in the raw FIFO, unclaimed.
    pub(crate) raw_dwell: DurationHistogram,
    /// Decompressed block sitting in the reorder buffer, waiting to be consumed.
    pub(crate) reorder_dwell: DurationHistogram,
}

impl BlockLifecycleStats {
    /// Snapshot the dwell stages, leaving the working stages at their default.
    ///
    /// `read_batch` and `decompress` come from
    /// [`MergePhaseCounters`](crate::merge_phases::MergePhaseCounters), which this
    /// struct cannot see. [`SortWorkerPool::block_lifecycle_report`] is the only
    /// caller that holds both halves, and it fills them in.
    ///
    /// [`SortWorkerPool::block_lifecycle_report`]: crate::worker_pool::SortWorkerPool
    pub(crate) fn snapshot(&self) -> BlockLifecycleReport {
        BlockLifecycleReport {
            raw_dwell: self.raw_dwell.snapshot(),
            reorder_dwell: self.reorder_dwell.snapshot(),
            ..BlockLifecycleReport::default()
        }
    }
}

/// One spill block's full journey: the dwell stages from
/// [`BlockLifecycleStats`] plus the working stages from
/// [`MergePhaseCounters`](crate::merge_phases::MergePhaseCounters).
#[derive(Debug, Clone, Copy, Default)]
pub struct BlockLifecycleReport {
    /// Disk reads.
    pub read_batch: HistogramReport,
    /// Raw FIFO dwell.
    pub raw_dwell: HistogramReport,
    /// Decompression.
    pub decompress: HistogramReport,
    /// Reorder buffer dwell.
    pub reorder_dwell: HistogramReport,
}

impl BlockLifecycleReport {
    /// Whether any block completed a stage.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.decompress.is_empty() && self.read_batch.is_empty()
    }

    /// Whether the reorder buffer is doing any real buffering.
    ///
    /// A block that is consumed almost the instant it is inserted was never
    /// lookahead — it was a block the consumer was already waiting for. When
    /// this is true, raising the per-file cap cannot help, because the cap is
    /// not what the pipeline is running into.
    #[must_use]
    pub fn reorder_is_pass_through(self) -> bool {
        !self.reorder_dwell.is_empty()
            && self.reorder_dwell.percentile_micros(0.90) <= PASS_THROUGH_MICROS
    }
}

/// p90 reorder dwell at or below which the buffer is a pass-through.
///
/// One decompression costs ~50 µs on the measured workload, so a block that
/// waits less than that to be consumed did not spend meaningful time buffered.
const PASS_THROUGH_MICROS: u64 = 50;

// ============================================================================
// Refill latency
// ============================================================================

/// What a file had available at the moment its reorder buffer drained to zero.
///
/// This is the fork in the road for the hungry-file question, and each branch
/// has a different fix: `RawReady` means the bytes were there and nothing
/// picked them up (scheduling); `Decompressing` means the pipeline was already
/// working and just is not ahead (lookahead depth); `Dry` means nothing had
/// been read yet (read concurrency).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EmptyCause {
    /// Raw blocks were queued and unclaimed.
    RawReady,
    /// A decompression was already in flight.
    Decompressing,
    /// Neither: the file had nothing anywhere.
    Dry,
}

impl EmptyCause {
    /// Number of variants, so the cause arrays indexed by `self as usize`
    /// cannot fall behind the enum. Sizing them with a literal instead lets a
    /// new variant compile and then index out of bounds inside
    /// [`RefillStats::record_empty`] -- on the consumer thread, where it kills
    /// the merge.
    pub(crate) const COUNT: usize = 3;

    /// Round-trip through a `u8`, so a file can carry the cause of its current
    /// empty in an atomic. Any value outside the enum decodes to [`Self::Dry`],
    /// which is unreachable in practice -- the only writer is
    /// [`Self::as_u8`].
    #[must_use]
    pub(crate) const fn as_u8(self) -> u8 {
        self as u8
    }

    /// Inverse of [`Self::as_u8`].
    #[must_use]
    pub(crate) const fn from_u8(value: u8) -> Self {
        match value {
            0 => Self::RawReady,
            1 => Self::Decompressing,
            _ => Self::Dry,
        }
    }

    /// Classify from a file's depths at the moment it emptied.
    #[must_use]
    pub(crate) fn classify(raw_len: usize, in_flight: usize) -> Self {
        if in_flight > 0 {
            Self::Decompressing
        } else if raw_len > 0 {
            Self::RawReady
        } else {
            Self::Dry
        }
    }
}

// Keeps `COUNT` and `from_u8` honest: the match is exhaustive, so a fourth
// variant fails to compile here instead of indexing past the end of the cause
// arrays at run time, or being folded silently into `Dry` by `from_u8`'s
// catch-all.
const _: () = {
    const fn assert_count(cause: EmptyCause) -> usize {
        match cause {
            EmptyCause::RawReady => 0,
            EmptyCause::Decompressing => 1,
            EmptyCause::Dry => EmptyCause::COUNT - 1,
        }
    }
    assert!(assert_count(EmptyCause::Dry) == 2);
};

/// How long it takes to put a block back into a file that has run dry.
#[derive(Debug, Default)]
pub(crate) struct RefillStats {
    /// Buffer emptied to a worker claiming a raw block for that file.
    pub(crate) claim_lag: DurationHistogram,
    /// Buffer emptied to a block landing back in it.
    pub(crate) insert_lag: DurationHistogram,
    /// Buffer emptied to a disk read completing, for files that were also dry.
    pub(crate) read_lag: DurationHistogram,
    /// Empties by [`EmptyCause`].
    causes: [AtomicU64; EmptyCause::COUNT],
}

impl RefillStats {
    /// Record that a file's reorder buffer drained to nothing.
    pub(crate) fn record_empty(&self, cause: EmptyCause) {
        self.causes[cause as usize].fetch_add(1, Ordering::Relaxed);
    }

    pub(crate) fn snapshot(&self) -> RefillReport {
        RefillReport {
            claim_lag: self.claim_lag.snapshot(),
            insert_lag: self.insert_lag.snapshot(),
            read_lag: self.read_lag.snapshot(),
            causes: std::array::from_fn(|i| self.causes[i].load(Ordering::Relaxed)),
        }
    }
}

/// Read-only view of [`RefillStats`].
#[derive(Debug, Clone, Copy, Default)]
pub struct RefillReport {
    /// Emptied to claimed.
    pub claim_lag: HistogramReport,
    /// Emptied to inserted.
    pub insert_lag: HistogramReport,
    /// Emptied to read completed, where the file was dry.
    pub read_lag: HistogramReport,
    /// Empties by [`EmptyCause`] as `usize`.
    pub causes: [u64; EmptyCause::COUNT],
}

impl RefillReport {
    /// Whether any file ever ran dry.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.causes.iter().all(|&n| n == 0)
    }

    /// Total empties observed.
    #[must_use]
    pub fn empties(self) -> u64 {
        self.causes.iter().sum()
    }

    /// Share of empties that hit `cause`, in `[0, 1]`.
    #[must_use]
    #[allow(clippy::cast_precision_loss, reason = "empty counts stay far below 2^52")]
    pub fn cause_share(self, cause: EmptyCause) -> f64 {
        let total = self.empties();
        if total == 0 {
            return 0.0;
        }
        self.causes[cause as usize] as f64 / total as f64
    }

    /// Share of refill latency spent waiting for a worker to claim the work,
    /// as opposed to doing it, in `[0, 1]`.
    ///
    /// The single most actionable number here. High means the pipeline is
    /// mostly waiting to *start*, so the fix is scheduling — wake a worker
    /// sooner, or apply more of them. Low means the pipeline starts promptly
    /// and the wait is the work itself, so the fix is to be further ahead when
    /// the buffer drains.
    ///
    /// The two histograms have slightly different populations — a file at end
    /// of input can empty and never be claimed again — so this is a ratio of
    /// totals rather than a paired per-cycle mean. At the scale it is read at
    /// (thousands of cycles, a handful unmatched) that does not change the
    /// reading, but it is why the counts printed beside them differ.
    #[must_use]
    pub fn claim_share(self) -> f64 {
        let insert = self.insert_lag.total_secs();
        if insert <= 0.0 {
            return 0.0;
        }
        (self.claim_lag.total_secs() / insert).min(1.0)
    }
}

// ============================================================================
// Consumer trace
// ============================================================================

/// Largest concurrent-decompression count tracked exactly; higher is clamped.
pub(crate) const MAX_TRACKED_IN_FLIGHT: usize = 8;

/// Nanoseconds one block is scaled to when a count reuses the duration
/// histogram's log2 bucketing.
///
/// 1 µs per block, so the bucketing sees blocks where it expects microseconds
/// and [`HistogramReport::summary_blocks`] divides back out. Shared by the
/// recorder and the formatter so the scale cannot drift between them.
pub(crate) const BLOCKS_TO_NANOS: u64 = 1_000;

/// What the consumer's own access pattern looks like, and what it costs.
#[derive(Debug)]
pub(crate) struct ConsumerTraceStats {
    /// Park duration, split by what the awaited file was doing.
    park_by_state: [DurationHistogram; AwaitedState::COUNT],
    /// Consecutive block pulls from one source before the merge switches.
    pub(crate) source_run_length: DurationHistogram,
    /// Concurrent decompressions on the awaited file at a park, clamped.
    in_flight_at_park: [AtomicU64; MAX_TRACKED_IN_FLIGHT + 1],
}

impl Default for ConsumerTraceStats {
    fn default() -> Self {
        Self {
            park_by_state: std::array::from_fn(|_| DurationHistogram::default()),
            source_run_length: DurationHistogram::default(),
            in_flight_at_park: std::array::from_fn(|_| AtomicU64::new(0)),
        }
    }
}

impl ConsumerTraceStats {
    /// Record a park of `nanos` that began with the awaited file in `state`,
    /// with `in_flight` decompressions running on it.
    pub(crate) fn record_park(&self, state: AwaitedState, nanos: u64, in_flight: usize) {
        self.park_by_state[state as usize].record(nanos);
        self.in_flight_at_park[in_flight.min(MAX_TRACKED_IN_FLIGHT)]
            .fetch_add(1, Ordering::Relaxed);
    }

    /// Record that the merge drew `blocks` consecutive blocks from one source.
    ///
    /// This is the test for whether the merge has a *hot file* worth steering
    /// workers toward, and it must be read weighted by **blocks, not by runs** --
    /// the two give opposite answers.
    ///
    /// By run count the merge looks perfectly interleaved: the measured median run
    /// is 1 block, so most of the time the loser tree does pull the next record
    /// from some other source. An earlier version of this comment stopped there and
    /// concluded that any fix framed as "prefetch the file the consumer needs next"
    /// was answering a question the workload does not pose.
    ///
    /// Weighted by blocks the same histogram says the opposite. On a measured
    /// 44-way merge (`n=167624 total=5306698 mean=31.7 p50=1 p90=2 p99=512`), if
    /// the bottom 90% of runs average two blocks or fewer they account for at most
    /// 302k of 5.31M blocks -- so **at least 94% of all blocks come from the top
    /// decile of runs**, which are hundreds of consecutive blocks from one source.
    /// Steering read-ahead at the source the merge is draining, and at the one it
    /// will drain next, is worth 2.8-5.2% and 2.7% of merge wall respectively.
    ///
    /// Both readings are true; only the block-weighted one predicted a result. It
    /// is the same count-versus-time trap that made the park census point at a
    /// scheduling change that then measured 0.06%. A short floor for the frequent
    /// one-block runs still has to exist -- what does not work is spreading depth
    /// uniformly, which has been measured negative four separate times.
    ///
    /// Reuses the duration histogram's bucketing for a count: run lengths are as
    /// heavy-tailed as the timings and want the same log2 treatment. One block
    /// is scaled to [`BLOCKS_TO_NANOS`] so it lands in the microsecond lane the
    /// bucketing is built around.
    ///
    /// Read the result back with [`HistogramReport::summary_blocks`], never
    /// [`HistogramReport::summary`] -- the latter labels these as durations.
    pub(crate) fn record_source_run(&self, blocks: u64) {
        self.source_run_length.record_count(blocks);
    }

    pub(crate) fn snapshot(&self) -> ConsumerTraceReport {
        ConsumerTraceReport {
            park_by_state: std::array::from_fn(|i| self.park_by_state[i].snapshot()),
            source_run_length: self.source_run_length.snapshot(),
            in_flight_at_park: std::array::from_fn(|i| {
                self.in_flight_at_park[i].load(Ordering::Relaxed)
            }),
        }
    }
}

/// Read-only view of [`ConsumerTraceStats`].
#[derive(Debug, Clone, Copy)]
pub struct ConsumerTraceReport {
    /// Park durations by [`AwaitedState`] as `usize`.
    pub park_by_state: [HistogramReport; AwaitedState::COUNT],
    /// Consecutive blocks drawn from one source (values are blocks, not µs).
    pub source_run_length: HistogramReport,
    /// Parks by concurrent decompressions on the awaited file, clamped at
    /// [`MAX_TRACKED_IN_FLIGHT`].
    pub in_flight_at_park: [u64; MAX_TRACKED_IN_FLIGHT + 1],
}

impl ConsumerTraceReport {
    /// Whether the consumer recorded nothing at all.
    ///
    /// Parks are not the only observation here: a merge that never stalls still
    /// records a source run on every block it pulls, and that distribution is
    /// what says whether there is a hot file worth steering workers toward. So
    /// both halves must be empty, or the healthy merge loses the one diagnostic
    /// it had.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.park_by_state.iter().all(|h| h.is_empty()) && self.source_run_length.is_empty()
    }

    /// Parks where exactly one worker was decompressing the awaited file.
    ///
    /// The diagnostic for a demand-blind scheduler: if the consumer is starved
    /// and only one of N workers is working on the file it needs, the pool is
    /// not applying its capacity where the merge is actually blocked.
    #[must_use]
    pub fn single_worker_parks(self) -> u64 {
        self.in_flight_at_park[1]
    }

    /// Parks with no decompression running on the awaited file at all.
    #[must_use]
    pub fn idle_file_parks(self) -> u64 {
        self.in_flight_at_park[0]
    }

    /// Parks where two or more workers were decompressing the awaited file.
    ///
    /// The complement of [`Self::idle_file_parks`] and
    /// [`Self::single_worker_parks`], and the one that distinguishes the two
    /// explanations for a starved consumer. If the pool is mostly *not* on the
    /// awaited file, it is scheduling badly and can be steered. If it is on that
    /// file several workers deep and the consumer still waits, the workers are
    /// decompressing *different* serials while the merge needs the lowest one,
    /// so the park is the head block's decompress latency and no amount of
    /// additional concurrency on that file will shorten it.
    #[must_use]
    pub fn multi_worker_parks(self) -> u64 {
        self.in_flight_at_park[2..].iter().sum()
    }

    /// Mean concurrent decompressions on the awaited file at a park.
    ///
    /// A *lower bound* once the pool saturates: the histogram clamps at
    /// [`MAX_TRACKED_IN_FLIGHT`], so every deeper park is counted as exactly
    /// that. Read it against the cap rather than as an absolute -- a mean at the
    /// clamp means "at least this deep", not "exactly this deep".
    #[must_use]
    #[allow(
        clippy::cast_precision_loss,
        reason = "park counts stay far below 2^52; depths are bounded by MAX_TRACKED_IN_FLIGHT"
    )]
    pub fn mean_in_flight(self) -> f64 {
        let parks = self.parks();
        if parks == 0 {
            return 0.0;
        }
        let weighted: u64 = self
            .in_flight_at_park
            .iter()
            .enumerate()
            .map(|(depth, count)| depth as u64 * count)
            .sum();
        weighted as f64 / parks as f64
    }

    /// Total parks counted in the in-flight histogram.
    #[must_use]
    pub fn parks(self) -> u64 {
        self.in_flight_at_park.iter().sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(0, 0)]
    #[case(999, 0)]
    #[case(1_000, 1)]
    #[case(1_999, 1)]
    #[case(2_000, 2)]
    #[case(50_000, 6)]
    // 1 ms == 1000 us, which lands in [512, 1024) -- bucket 10, not 11.
    #[case(1_000_000, 10)]
    #[case(1_048_576_000, 21)]
    fn test_bucket_of_nanos(#[case] nanos: u64, #[case] expected: usize) {
        assert_eq!(bucket_of_nanos(nanos), expected);
    }

    /// Anything longer than the top bucket's floor lands in the top bucket
    /// rather than indexing off the end.
    #[test]
    fn test_bucket_saturates() {
        assert_eq!(bucket_of_nanos(u64::MAX), HIST_BUCKETS - 1);
        assert_eq!(bucket_of_nanos(3_600 * 1_000_000_000), HIST_BUCKETS - 1);
    }

    #[rstest]
    #[case(0, 0)]
    #[case(1, 1)]
    #[case(6, 32)]
    #[case(11, 1024)]
    fn test_bucket_floor(#[case] bucket: usize, #[case] expected: u64) {
        assert_eq!(bucket_floor_micros(bucket), expected);
    }

    #[test]
    fn test_histogram_keeps_an_exact_mean_and_a_bucketed_tail() {
        let hist = DurationHistogram::default();
        for _ in 0..99 {
            hist.record(10_000); // 10 us
        }
        hist.record(10_000_000); // one 10 ms outlier
        let report = hist.snapshot();
        assert_eq!(report.count, 100);
        // The mean is exact, so the outlier moves it: (99*10 + 10000)/100 us.
        assert!((report.mean_micros() - 109.9).abs() < 0.01, "{}", report.mean_micros());
        // ...while the median is untouched by it, which is the point.
        assert_eq!(report.percentile_micros(0.50), 8);
        assert_eq!(report.percentile_micros(0.99), 8);
        assert_eq!(report.percentile_micros(1.00), 8192);
    }

    /// A queue depth is a single-digit count. Recording it raw drops every
    /// observation into bucket 0 and the report reads zero -- the merge-trace
    /// bug this guards. [`DurationHistogram::record_count`] scales one unit to
    /// [`BLOCKS_TO_NANOS`] first, so the mean and the tail read as the depths
    /// themselves.
    #[test]
    fn test_record_count_scales_small_counts_out_of_the_zero_bucket() {
        let depths = [2_u64, 3, 4, 9];

        // Raw counts collapse: every depth below `BLOCKS_TO_NANOS` buckets to
        // zero, so the mean rounds to 0.0 and the tail is a flat zero.
        let raw = DurationHistogram::default();
        for &depth in &depths {
            raw.record(depth);
        }
        let raw = raw.snapshot();
        assert_eq!(raw.count, 4);
        assert_eq!(format!("{:.1}", raw.mean_micros()), "0.0", "raw counts read as zero");
        assert_eq!(raw.percentile_micros(1.0), 0, "raw counts collapse to a zero tail");

        // Scaled through `record_count`, the same depths read back as blocks.
        let scaled = DurationHistogram::default();
        for &depth in &depths {
            scaled.record_count(depth);
        }
        let scaled = scaled.snapshot();
        assert_eq!(scaled.count, 4);
        // (2 + 3 + 4 + 9) / 4 = 4.5 blocks, reported exactly.
        assert!((scaled.mean_micros() - 4.5).abs() < 0.01, "{}", scaled.mean_micros());
        // Bucketed as if microseconds, so the tail reads as blocks (log2 floor).
        assert_eq!(scaled.percentile_micros(1.0), 8);
        assert_eq!(
            scaled.summary_blocks(),
            "n=4 total=18 blocks mean=4.5 p50=2 p90=8 p99=8",
            "the depth line must read as blocks"
        );
    }

    /// A mean alone cannot tell these two apart, and they have different
    /// causes: one is uniformly slow, the other is fast with a bad tail.
    #[test]
    fn test_percentiles_separate_distributions_with_equal_means() {
        let uniform = DurationHistogram::default();
        for _ in 0..100 {
            uniform.record(100_000); // every observation 100 us
        }
        let tailed = DurationHistogram::default();
        for _ in 0..99 {
            tailed.record(1_000); // 1 us
        }
        tailed.record(9_901_000); // one 9.9 ms outlier
        let (u, t) = (uniform.snapshot(), tailed.snapshot());
        assert!((u.mean_micros() - t.mean_micros()).abs() < 1.0, "means should match");
        assert!(t.percentile_micros(0.50) < u.percentile_micros(0.50) / 10);
    }

    /// `snapshot` reads the buckets and the counter as separate relaxed loads,
    /// so a report taken while workers are still recording can carry a `count`
    /// larger than the buckets sum to. The percentiles must still describe the
    /// observations actually present, not fall through to the top bucket --
    /// which would print a p50 of 8us as ~4.2s in the merge trace.
    #[test]
    fn test_percentiles_ignore_a_count_that_ran_ahead_of_the_buckets() {
        let hist = DurationHistogram::default();
        for _ in 0..100 {
            hist.record(10_000); // 10 us
        }
        let mut skewed = hist.snapshot();
        // Stand in for observations counted after the buckets were read.
        skewed.count += 50;

        assert_eq!(skewed.percentile_micros(0.50), 8, "p50 must stay in the populated bucket");
        assert_eq!(skewed.percentile_micros(0.99), 8);
        assert_eq!(skewed.percentile_micros(1.00), 8);
    }

    /// A report whose buckets are all empty has nothing to take a percentile
    /// of, whatever `count` claims.
    #[test]
    fn test_percentiles_are_zero_when_only_the_count_is_populated() {
        let mut report = DurationHistogram::default().snapshot();
        report.count = 7;
        assert_eq!(report.percentile_micros(0.50), 0);
    }

    #[test]
    fn test_empty_histogram_reports_zeroes() {
        let report = DurationHistogram::default().snapshot();
        assert!(report.is_empty());
        assert!((report.mean_micros() - 0.0).abs() < f64::EPSILON);
        assert_eq!(report.percentile_micros(0.99), 0);
    }

    #[rstest]
    // A decompression already running means the pipeline is working but is not
    // ahead -- a lookahead-depth problem, not a scheduling one.
    #[case(0, 1, EmptyCause::Decompressing)]
    #[case(4, 2, EmptyCause::Decompressing)]
    // Bytes on hand and nobody on them: scheduling.
    #[case(4, 0, EmptyCause::RawReady)]
    // Nothing anywhere: the read side is behind.
    #[case(0, 0, EmptyCause::Dry)]
    fn test_empty_cause(
        #[case] raw_len: usize,
        #[case] in_flight: usize,
        #[case] expected: EmptyCause,
    ) {
        assert_eq!(EmptyCause::classify(raw_len, in_flight), expected);
    }

    #[test]
    fn test_claim_share_splits_waiting_to_start_from_doing_the_work() {
        let stats = RefillStats::default();
        // 8 ms of refill latency, 6 ms of it spent before any worker started.
        stats.claim_lag.record(6_000_000);
        stats.insert_lag.record(8_000_000);
        stats.record_empty(EmptyCause::RawReady);
        let report = stats.snapshot();
        assert!((report.claim_share() - 0.75).abs() < 1e-9);
        assert!((report.cause_share(EmptyCause::RawReady) - 1.0).abs() < 1e-9);
        assert_eq!(report.empties(), 1);
    }

    /// Clock skew between two histograms must not produce a share above 1.
    #[test]
    fn test_claim_share_is_clamped() {
        let stats = RefillStats::default();
        stats.claim_lag.record(9_000_000);
        stats.insert_lag.record(8_000_000);
        assert!((stats.snapshot().claim_share() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_empty_refill_report() {
        let report = RefillStats::default().snapshot();
        assert!(report.is_empty());
        assert!((report.claim_share() - 0.0).abs() < f64::EPSILON);
        assert!((report.cause_share(EmptyCause::Dry) - 0.0).abs() < f64::EPSILON);
    }

    /// A reorder buffer whose blocks are consumed the moment they arrive is not
    /// buffering, so its cap is not the constraint -- the check that stops a
    /// cap sweep being run for the wrong reason.
    #[test]
    fn test_pass_through_reorder_buffer_is_detected() {
        let stats = BlockLifecycleStats::default();
        for _ in 0..100 {
            stats.reorder_dwell.record(2_000); // 2 us: consumed immediately
        }
        assert!(stats.snapshot().reorder_is_pass_through());

        let buffered = BlockLifecycleStats::default();
        for _ in 0..100 {
            buffered.reorder_dwell.record(5_000_000); // 5 ms: genuinely queued
        }
        assert!(!buffered.snapshot().reorder_is_pass_through());
    }

    #[test]
    fn test_empty_lifecycle_report_is_not_pass_through() {
        let report = BlockLifecycleStats::default().snapshot();
        assert!(report.is_empty());
        assert!(!report.reorder_is_pass_through(), "no data must not read as a finding");
    }

    #[test]
    fn test_consumer_trace_counts_parks_by_concurrency() {
        let stats = ConsumerTraceStats::default();
        stats.record_park(AwaitedState::Decompressing, 80_000, 1);
        stats.record_park(AwaitedState::Decompressing, 90_000, 1);
        stats.record_park(AwaitedState::RawQueued, 40_000, 0);
        // Clamped: 12 concurrent decompressions is reported as the max bucket.
        stats.record_park(AwaitedState::Decompressing, 10_000, 12);
        let report = stats.snapshot();
        assert_eq!(report.parks(), 4);
        assert_eq!(report.single_worker_parks(), 2);
        assert_eq!(report.idle_file_parks(), 1);
        assert_eq!(report.in_flight_at_park[MAX_TRACKED_IN_FLIGHT], 1);
        assert_eq!(report.park_by_state[AwaitedState::Decompressing as usize].count, 3);
    }

    /// Park time per state, not just park count.
    ///
    /// The two disagree and only one is a cost: the park-supply census had a class
    /// that was 41% of parks and 13% of park time, and reading the counts picked a
    /// fix that then measured at 0.06%. These nanoseconds were collected from the
    /// start and reported one log level down, so every experiment in that
    /// investigation had the answer and none printed it.
    #[test]
    fn test_consumer_trace_carries_park_time_per_state_not_only_counts() {
        let stats = ConsumerTraceStats::default();
        // Many cheap parks in one state, one expensive park in another: the count
        // majority and the time majority are deliberately opposite.
        for _ in 0..10 {
            stats.record_park(AwaitedState::Decompressing, 1_000, 1);
        }
        stats.record_park(AwaitedState::RawQueued, 500_000, 0);

        let report = stats.snapshot();
        let decomp = report.park_by_state[AwaitedState::Decompressing as usize];
        let raw = report.park_by_state[AwaitedState::RawQueued as usize];

        assert_eq!(decomp.count, 10, "counts say Decompressing dominates 10:1");
        assert_eq!(raw.count, 1);
        assert!(
            (decomp.total_secs() - 10e-6).abs() < 1e-12,
            "10 parks of 1us is 10us, got {}",
            decomp.total_secs()
        );
        assert!(
            (raw.total_secs() - 500e-6).abs() < 1e-12,
            "one park of 500us is 500us, got {}",
            raw.total_secs()
        );
        assert!(
            raw.total_secs() > decomp.total_secs() * 40.0,
            "time says RawQueued dominates 50:1 -- the opposite of the counts, which is \
             precisely the reading that must be available at info level"
        );
    }

    #[test]
    fn test_consumer_trace_reports_depth_beyond_one_worker() {
        // `none` and `exactly one` cannot distinguish a pool running two deep
        // from one running at the cap, yet that is the difference between
        // "scheduling headroom remains" and "decompress latency is irreducible".
        let stats = ConsumerTraceStats::default();
        stats.record_park(AwaitedState::Decompressing, 10_000, 0);
        stats.record_park(AwaitedState::Decompressing, 10_000, 1);
        stats.record_park(AwaitedState::Decompressing, 10_000, 4);
        stats.record_park(AwaitedState::Decompressing, 10_000, 4);
        let report = stats.snapshot();

        assert_eq!(report.multi_worker_parks(), 2, "two parks had 2+ workers");
        assert!(
            (report.mean_in_flight() - 2.25).abs() < 1e-9,
            "mean of 0,1,4,4 is 2.25, got {}",
            report.mean_in_flight()
        );
    }

    #[test]
    fn test_mean_in_flight_is_zero_without_parks() {
        let report = ConsumerTraceStats::default().snapshot();
        assert_eq!(report.parks(), 0);
        assert!(report.mean_in_flight().abs() < f64::EPSILON, "no parks must not divide by zero");
    }

    #[test]
    fn test_source_run_length_is_recorded_in_blocks() {
        let stats = ConsumerTraceStats::default();
        stats.record_source_run(1);
        stats.record_source_run(64);
        let report = stats.snapshot();
        assert_eq!(report.source_run_length.count, 2);
        // Bucketed as if microseconds, so the "us" figures read as blocks.
        assert_eq!(report.source_run_length.percentile_micros(1.0), 64);
    }

    /// Run lengths are counts riding in a duration histogram, so the *formatter*
    /// is what keeps them honest. `summary` would call the total seconds and the
    /// rest microseconds; `summary_blocks` must report the block count itself.
    #[test]
    fn test_source_run_summary_reports_blocks_not_durations() {
        let stats = ConsumerTraceStats::default();
        for _ in 0..3 {
            stats.record_source_run(4);
        }
        stats.record_source_run(64);
        let report = stats.snapshot().source_run_length;

        // 3x4 + 64 = 76 blocks over 4 runs, mean 19.
        assert_eq!(
            report.summary_blocks(),
            "n=4 total=76 blocks mean=19.0 p50=4 p90=64 p99=64",
            "the run-length line must read as blocks"
        );
        // What the duration formatter would have printed instead: a 76-block
        // total rendered as 0.0s, and every figure labelled `us`.
        assert!(
            report.summary().contains("total=0.0s"),
            "guard the reason summary_blocks exists: {}",
            report.summary()
        );
    }

    #[test]
    fn test_empty_consumer_trace() {
        assert!(ConsumerTraceStats::default().snapshot().is_empty());
    }

    /// A merge that never stalls still has a run-length distribution, and it is
    /// the one that answers whether there is a hot file to steer workers
    /// toward. Gating the whole consumer block on parks would drop it in
    /// exactly the healthy case it is cheapest to read.
    #[test]
    fn test_source_runs_alone_make_the_consumer_trace_non_empty() {
        let stats = ConsumerTraceStats::default();
        stats.record_source_run(4);
        assert!(
            !stats.snapshot().is_empty(),
            "a merge with source runs but no parks must still report them"
        );
    }
}
