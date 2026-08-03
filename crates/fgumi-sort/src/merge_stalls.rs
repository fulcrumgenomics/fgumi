//! Why the merge stalled — the counters that separate the candidate causes.
//!
//! [`crate::merge_phases`] answers *where the merge's time went*: how much
//! thread-busy time each component consumed, and what fraction of worker
//! capacity that used. On a spill-heavy whole-genome sort it reports a merge in
//! which the worker pool is ~55% utilized while the consumer spends ~79% of its
//! loop fetching records. Both halves are waiting, and neither number says what
//! for.
//!
//! It cannot say, because every candidate cause produces the same reading. A
//! worker that finds nothing to do sleeps and adds to one idle total, whether it
//! found nothing because
//!
//! - every file's reorder buffer was at [`PHASE2_DECOMP_CAP`] (the consumer is
//!   behind, and the pool is correctly throttled),
//! - every file's raw FIFO was at [`PHASE2_RAW_CAP`] (decompression is behind
//!   the read side, so read-ahead depth is the binding constraint),
//! - every file with data left already had a peer worker inside a disk read
//!   (waiting on I/O, which wants read concurrency rather than deeper buffers),
//!   or
//! - its `try_lock` calls lost races against other workers (contention).
//!
//! [`PHASE2_DECOMP_CAP`]: crate::worker_pool::PHASE2_DECOMP_CAP
//! [`PHASE2_RAW_CAP`]: crate::worker_pool::PHASE2_RAW_CAP
//!
//! The file-scan loop already distinguishes all four and then discards the
//! distinction. This module keeps it, plus two things nothing currently measures
//! at all:
//!
//! - **Where the consumer actually blocks.** `advance_to_next_block` parks the
//!   main thread waiting for one specific file's next block. Parks happen once
//!   per ~64 KB block rather than once per ~100-byte record, so unlike the
//!   sampled per-record timers they can be measured exactly and for free, and
//!   they isolate *waiting* from the record copy that the sampled "fetch next
//!   record" figure folds in with it. Censusing the other files at that moment
//!   is what separates head-of-line blocking (this file starved, the rest at
//!   cap) from a uniformly starved pool.
//! - **How late workers discover work.** Nothing ever unparks a pool worker: the
//!   wake path runs one way, workers to main thread. A worker that has backed
//!   off to [`MAX_BACKOFF_US`](crate::worker_pool::MAX_BACKOFF_US) therefore
//!   finds newly-available work up to a millisecond after it appears, and that
//!   delay is indistinguishable from I/O latency in every counter that exists
//!   today.
//!
//! # Cost
//!
//! Nothing here runs on a hot path. Skip reasons are tallied in a worker's stack
//! frame and flushed only on the scan that found no work — the path that is
//! about to sleep anyway. Park timing is per block. The park census is the one
//! genuinely non-trivial cost (a `try_lock` per file) and is sampled.
//!
//! # Calibration status
//!
//! [`classify_merge`](crate::merge_phases::classify_merge)'s thresholds were
//! fitted to measured runs. [`classify_stall`]'s were not — they encode the
//! hypotheses this instrumentation exists to test, and the first instrumented
//! run should either confirm or replace them. The test cases are named for the
//! hypothesis each one pins so that is not forgotten.

use std::sync::atomic::{AtomicU64, Ordering};

// ============================================================================
// Worker side: why a file-scan found no work
// ============================================================================

/// Why a worker skipped one spill file during a Phase 2 scan.
///
/// # There is deliberately no "starved" variant
///
/// A file with an empty FIFO and a reader that is neither at EOF nor held by
/// another worker is not skipped — the scan *reads* it, and the scan is
/// therefore productive. Starvation cannot appear as a skip reason; the closest
/// thing is [`ReadInProgress`](Self::ReadInProgress), which says a peer is
/// already fetching. Starvation shows up on the consumer side instead, as a park
/// while the pool is busy, which is why both halves of this module are needed to
/// answer the question.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Phase2Skip {
    /// The reader is at EOF and the raw FIFO is empty — nothing left for a
    /// worker to take, though blocks may still be in flight or waiting for the
    /// consumer to collect them.
    Drained,
    /// The reorder buffer is at `PHASE2_DECOMP_CAP` and the head raw block is
    /// not the gap-filler the consumer is stuck on. Backpressure working as
    /// intended: the consumer is behind, not the pool.
    DecompCapped,
    /// The raw FIFO is at `PHASE2_RAW_CAP`, so no more disk read-ahead is
    /// admitted. Blocks are arriving faster than they are being decompressed.
    RawFull,
    /// Another worker holds this file's reader lock, i.e. a disk read for it is
    /// already underway. Not a defect — but a scan dominated by it means
    /// workers are queued behind in-flight reads, and more read concurrency,
    /// not deeper buffers, is what would help.
    ReadInProgress,
    /// A `try_lock` on the raw FIFO or the reorder buffer lost a race. Work may
    /// well have been available; this worker could not see it.
    LockContended,
}

impl Phase2Skip {
    /// Number of variants, for fixed-size counter arrays.
    pub(crate) const COUNT: usize = 5;

    /// Every variant, in declaration order, for iteration and display.
    pub(crate) const ALL: [Self; Self::COUNT] = [
        Self::Drained,
        Self::DecompCapped,
        Self::RawFull,
        Self::ReadInProgress,
        Self::LockContended,
    ];

    /// Short label for log output.
    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::Drained => "drained",
            Self::DecompCapped => "decomp-capped",
            Self::RawFull => "raw-full",
            Self::ReadInProgress => "read-in-progress",
            Self::LockContended => "lock-contended",
        }
    }
}

/// Why the decompress half of a file scan declined to take a block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PopSkip {
    /// The raw FIFO's `try_lock` lost a race.
    RawLockContended,
    /// The raw FIFO is empty; there is nothing decompressible buffered.
    RawEmpty,
    /// The reorder buffer's `try_lock` lost a race.
    DecompLockContended,
    /// The reorder buffer is at cap and the head block is not the gap-filler.
    DecompCapped,
}

/// Why the read half of a file scan declined to fetch more blocks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ReadSkip {
    /// Another worker holds the reader lock.
    ReaderBusy,
    /// The reader has reached EOF.
    ReaderEof,
    /// The raw FIFO's `try_lock` lost a race, so its depth is unknown.
    RawLockContended,
    /// The raw FIFO is at `PHASE2_RAW_CAP`.
    RawFull,
}

/// Reduce a file's two half-scan outcomes to the one reason worth counting.
///
/// Both halves decline independently and often for related reasons (a file at
/// its decompress cap will usually also have a full raw FIFO), so the pair has
/// to be collapsed rather than counted twice. Backpressure wins over everything
/// else because it is the only outcome that says the pipeline is *working* —
/// this file has data and the pool is deliberately not taking it.
#[must_use]
pub(crate) fn combine_skip(pop: PopSkip, read: ReadSkip) -> Phase2Skip {
    match (pop, read) {
        (PopSkip::DecompCapped, _) => Phase2Skip::DecompCapped,
        (PopSkip::RawLockContended | PopSkip::DecompLockContended, _)
        | (_, ReadSkip::RawLockContended) => Phase2Skip::LockContended,
        (_, ReadSkip::ReaderBusy) => Phase2Skip::ReadInProgress,
        (_, ReadSkip::RawFull) => Phase2Skip::RawFull,
        (PopSkip::RawEmpty, ReadSkip::ReaderEof) => Phase2Skip::Drained,
    }
}

/// Skip reasons accumulated over one worker's scan of every spill file.
///
/// Lives in the worker's stack frame during the scan and is discarded if the
/// scan finds work, so a productive scan pays five increments of a local array
/// and nothing else.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Phase2ScanTally {
    counts: [u32; Phase2Skip::COUNT],
}

impl Phase2ScanTally {
    /// Note that one file was skipped for `reason`.
    #[inline]
    pub(crate) fn note(&mut self, reason: Phase2Skip) {
        self.counts[reason as usize] += 1;
    }

    /// How many files were skipped for `reason`.
    #[must_use]
    pub(crate) fn count(&self, reason: Phase2Skip) -> u32 {
        self.counts[reason as usize]
    }

    /// Files skipped for any reason.
    #[must_use]
    pub(crate) fn total(&self) -> u32 {
        self.counts.iter().sum()
    }

    /// Files that could still have yielded work — everything except the ones
    /// that are legitimately finished.
    ///
    /// The denominator for every share below. A scan of 86 files where 80 are
    /// drained and 6 are contended is a *contended* scan, not a 7%-contended
    /// one; dividing by the live files rather than all of them is what makes the
    /// verdict hold up as the merge nears completion and file after file drops
    /// out.
    #[must_use]
    pub(crate) fn live(&self) -> u32 {
        self.total() - self.count(Phase2Skip::Drained)
    }
}

/// What a fruitless file scan says about the pipeline's state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ScanVerdict {
    /// Every file is finished; the merge is winding down.
    Drained,
    /// The live files are buffered to their caps. Workers are ahead of the
    /// consumer and idling *because they are throttled*, which means raising
    /// the caps buys buffering, not throughput.
    Backpressured,
    /// Lock contention hid a material share of the live files from this worker.
    Contended,
    /// The live files already have reads underway. Workers are queued behind
    /// in-flight disk I/O, so read concurrency — not buffer depth — is the
    /// constraint.
    IoWait,
    /// No single explanation dominates.
    Mixed,
}

impl ScanVerdict {
    /// Number of variants, so the counter arrays indexed by `self as usize`
    /// cannot fall behind the enum. Sizing them with a literal instead lets a
    /// new variant compile and then index out of bounds inside
    /// [`Phase2ScanStats::record_fruitless_scan`] -- on a worker thread, where
    /// it surfaces only as a panicked sort worker.
    pub(crate) const COUNT: usize = 5;
}

// Keeps `COUNT` honest: the match is exhaustive, so a sixth variant fails to
// compile here instead of indexing past the end of the counter arrays at run
// time.
const _: () = {
    const fn assert_count(verdict: ScanVerdict) -> usize {
        match verdict {
            ScanVerdict::Drained => 0,
            ScanVerdict::Backpressured => 1,
            ScanVerdict::Contended => 2,
            ScanVerdict::IoWait => 3,
            ScanVerdict::Mixed => ScanVerdict::COUNT - 1,
        }
    }
    assert!(assert_count(ScanVerdict::Mixed) == 4);
};

/// Share of live files above which one skip reason explains the scan.
const DOMINANT_SHARE: f64 = 0.60;
/// Share of live files that must be lock-contended to call the scan contended.
///
/// Lower than [`DOMINANT_SHARE`] on purpose: contention is a defect rather than
/// a régime, and a quarter of the pool's view being blocked by `try_lock` misses
/// is already worth reporting even when another reason is more common.
const CONTENDED_SHARE: f64 = 0.25;

/// Classify one fruitless scan.
///
/// Pure so the thresholds are unit-testable without a pipeline, following
/// [`classify_merge`](crate::merge_phases::classify_merge).
#[must_use]
#[allow(
    clippy::cast_precision_loss,
    reason = "file counts are bounded by the open-file limit; f64 represents them exactly"
)]
pub(crate) fn classify_scan(tally: Phase2ScanTally) -> ScanVerdict {
    let live = tally.live();
    if live == 0 {
        return ScanVerdict::Drained;
    }
    let share = |reason: Phase2Skip| f64::from(tally.count(reason)) / f64::from(live);

    if share(Phase2Skip::LockContended) >= CONTENDED_SHARE {
        ScanVerdict::Contended
    } else if share(Phase2Skip::DecompCapped) + share(Phase2Skip::RawFull) >= DOMINANT_SHARE {
        ScanVerdict::Backpressured
    } else if share(Phase2Skip::ReadInProgress) >= DOMINANT_SHARE {
        ScanVerdict::IoWait
    } else {
        ScanVerdict::Mixed
    }
}

/// Pool-wide totals for fruitless Phase 2 scans.
///
/// Written by every worker, so atomic; `Relaxed` throughout because these are
/// diagnostics with no ordering relationship to the data they describe.
#[derive(Debug, Default)]
pub(crate) struct Phase2ScanStats {
    /// Fruitless scans, by [`ScanVerdict`].
    verdicts: [AtomicU64; ScanVerdict::COUNT],
    /// Files skipped, by [`Phase2Skip`], summed over fruitless scans.
    skips: [AtomicU64; Phase2Skip::COUNT],
    /// Number of fruitless scans.
    scans: AtomicU64,
}

impl Phase2ScanStats {
    /// Record one scan that found no work.
    pub(crate) fn record_fruitless_scan(&self, tally: Phase2ScanTally) {
        self.scans.fetch_add(1, Ordering::Relaxed);
        for reason in Phase2Skip::ALL {
            let n = u64::from(tally.count(reason));
            if n > 0 {
                self.skips[reason as usize].fetch_add(n, Ordering::Relaxed);
            }
        }
        let verdict = classify_scan(tally);
        self.verdicts[verdict as usize].fetch_add(1, Ordering::Relaxed);
    }

    pub(crate) fn snapshot(&self) -> Phase2ScanReport {
        Phase2ScanReport {
            scans: self.scans.load(Ordering::Relaxed),
            verdicts: std::array::from_fn(|i| self.verdicts[i].load(Ordering::Relaxed)),
            skips: std::array::from_fn(|i| self.skips[i].load(Ordering::Relaxed)),
        }
    }
}

/// Read-only view of [`Phase2ScanStats`], for logging.
#[derive(Debug, Clone, Copy)]
pub struct Phase2ScanReport {
    /// Number of scans that found no work.
    pub scans: u64,
    /// Fruitless scans by verdict, indexed by [`ScanVerdict`] as `usize`.
    pub verdicts: [u64; ScanVerdict::COUNT],
    /// Files skipped by reason, indexed by [`Phase2Skip`] as `usize`.
    pub skips: [u64; Phase2Skip::COUNT],
}

impl Phase2ScanReport {
    /// Whether any fruitless scan was recorded, so a sort that never merged
    /// stays silent.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.scans == 0
    }

    /// Share of fruitless scans that reached `verdict`, in `[0, 1]`.
    #[must_use]
    #[allow(
        clippy::cast_precision_loss,
        reason = "scan counts stay far below 2^52 on any real sort"
    )]
    pub(crate) fn verdict_share(self, verdict: ScanVerdict) -> f64 {
        if self.scans == 0 {
            return 0.0;
        }
        self.verdicts[verdict as usize] as f64 / self.scans as f64
    }
}

// ============================================================================
// Worker side: how late work is discovered
// ============================================================================

/// Backoff levels a worker can sleep at, in microseconds.
///
/// The worker loop doubles its backoff from `MIN_BACKOFF_US` and clamps at
/// `MAX_BACKOFF_US`, so these are exactly the reachable values.
/// `test_bucket_table_matches_the_worker_loop_backoff` pins that correspondence
/// — if the loop's constants or its doubling change, this table has to follow or
/// the histogram silently mislabels itself.
pub(crate) const WAKE_BUCKET_US: [u64; 8] = [10, 20, 40, 80, 160, 320, 640, 1000];

/// Bucket index for a sleep of `backoff_us`.
#[must_use]
pub(crate) fn wake_bucket(backoff_us: u64) -> usize {
    let last = WAKE_BUCKET_US.len() - 1;
    WAKE_BUCKET_US.iter().position(|&us| backoff_us <= us).unwrap_or(last)
}

/// How long workers sleep, and how much of that sleep precedes finding work.
///
/// The second half is the point. Idle time alone cannot distinguish a worker
/// that slept because there was nothing to do from one that slept *through* work
/// becoming available, and only the latter is a cost the pipeline pays. A
/// successful step recorded against a 1 ms bucket says the work was found no
/// sooner than a millisecond after the worker last looked.
#[derive(Debug, Default)]
pub(crate) struct WakeLatencyStats {
    /// Nanoseconds slept, by backoff level.
    idle_nanos: [AtomicU64; WAKE_BUCKET_US.len()],
    /// Sleeps, by backoff level.
    sleeps: [AtomicU64; WAKE_BUCKET_US.len()],
    /// Sleeps that were immediately followed by a successful step.
    productive_sleeps: [AtomicU64; WAKE_BUCKET_US.len()],
}

impl WakeLatencyStats {
    /// Record a sleep of `idle_nanos` at backoff level `backoff_us`.
    pub(crate) fn record_sleep(&self, backoff_us: u64, idle_nanos: u64) {
        let b = wake_bucket(backoff_us);
        self.idle_nanos[b].fetch_add(idle_nanos, Ordering::Relaxed);
        self.sleeps[b].fetch_add(1, Ordering::Relaxed);
    }

    /// Record that the sleep at `backoff_us` was followed by finding work.
    pub(crate) fn record_productive_sleep(&self, backoff_us: u64) {
        self.productive_sleeps[wake_bucket(backoff_us)].fetch_add(1, Ordering::Relaxed);
    }

    pub(crate) fn snapshot(&self) -> WakeLatencyReport {
        WakeLatencyReport {
            idle_nanos: std::array::from_fn(|i| self.idle_nanos[i].load(Ordering::Relaxed)),
            sleeps: std::array::from_fn(|i| self.sleeps[i].load(Ordering::Relaxed)),
            productive_sleeps: std::array::from_fn(|i| {
                self.productive_sleeps[i].load(Ordering::Relaxed)
            }),
        }
    }
}

impl WakeLatencyReport {
    /// This report minus `baseline`, i.e. what accumulated between them.
    ///
    /// The pool exists for the whole sort, so these counters include Phase 1's
    /// sleeps. Reporting the running total against merge wall clock overstates
    /// the merge's share -- on a measured run, 1815s of "deep-sleep
    /// worker-seconds" against a merge that only had 8 x 305s = 2442
    /// worker-seconds of capacity in the first place. Subtracting a snapshot
    /// taken at the phase boundary is what makes the figure comparable to the
    /// numbers printed beside it.
    #[must_use]
    pub fn since(self, baseline: Self) -> Self {
        Self {
            idle_nanos: std::array::from_fn(|i| {
                self.idle_nanos[i].saturating_sub(baseline.idle_nanos[i])
            }),
            sleeps: std::array::from_fn(|i| self.sleeps[i].saturating_sub(baseline.sleeps[i])),
            productive_sleeps: std::array::from_fn(|i| {
                self.productive_sleeps[i].saturating_sub(baseline.productive_sleeps[i])
            }),
        }
    }
}

/// Read-only view of [`WakeLatencyStats`], for logging.
#[derive(Debug, Clone, Copy)]
pub struct WakeLatencyReport {
    /// Nanoseconds slept, by backoff level.
    pub idle_nanos: [u64; WAKE_BUCKET_US.len()],
    /// Sleeps, by backoff level.
    pub sleeps: [u64; WAKE_BUCKET_US.len()],
    /// Sleeps immediately followed by a successful step, by backoff level.
    pub productive_sleeps: [u64; WAKE_BUCKET_US.len()],
}

impl WakeLatencyReport {
    /// Whether any sleep was recorded.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.sleeps.iter().all(|&n| n == 0)
    }

    /// Sleeps that were immediately followed by finding work.
    ///
    /// The denominator behind [`Self::deep_sleep_wake_share`]. Reported
    /// alongside it because that share is meaningless at small counts, and a
    /// pool that idles rarely produces very small counts.
    #[must_use]
    pub fn productive_sleep_count(self) -> u64 {
        self.productive_sleeps.iter().sum()
    }

    /// Estimated worker-seconds lost to discovering work late.
    ///
    /// A sleep followed by a successful step means the work was there when the
    /// worker woke and may have arrived at any point during the sleep, so half
    /// the sleep is the expected lag under a uniform arrival assumption.
    ///
    /// **This is per-worker discovery lag, not pipeline delay.** It only holds
    /// the merge up when every worker is asleep at once; with eight workers
    /// jittered against each other it usually does not. Read it as a ceiling on
    /// the cost, and compare it against the consumer's park total — if it is a
    /// rounding error next to that, scheduling latency is not the problem and
    /// the question is settled cheaply.
    #[must_use]
    #[allow(
        clippy::cast_precision_loss,
        reason = "sleep counts stay far below 2^52 on any real sort"
    )]
    pub fn estimated_discovery_lag_secs(self) -> f64 {
        self.productive_sleeps
            .iter()
            .zip(WAKE_BUCKET_US)
            .map(|(&n, us)| n as f64 * (us as f64 / 2.0) / 1e6)
            .sum()
    }

    /// Worker-seconds slept at the deepest backoff levels.
    ///
    /// The pool of idle time that discovery lag is drawn from. If this is small,
    /// [`Self::estimated_discovery_lag_secs`] cannot be large no matter how the
    /// productive sleeps fall, which makes it a quick sanity check on the
    /// estimate rather than a second measurement of it.
    #[must_use]
    #[allow(
        clippy::cast_precision_loss,
        reason = "nanosecond totals stay far below 2^52; a sort would need ~52 days of idle time to lose a nanosecond here"
    )]
    pub fn deep_sleep_idle_secs(self) -> f64 {
        self.idle_nanos[DEEP_SLEEP_FIRST_BUCKET..].iter().sum::<u64>() as f64 / 1e9
    }

    /// Share of productive sleeps that happened at the deepest backoff levels,
    /// in `[0, 1]`.
    ///
    /// High means workers are routinely finding work immediately after their
    /// longest sleep, which is the signature of a pool that is discovering work
    /// late rather than lacking it.
    #[must_use]
    #[allow(
        clippy::cast_precision_loss,
        reason = "sleep counts stay far below 2^52 on any real sort"
    )]
    pub fn deep_sleep_wake_share(self) -> f64 {
        let total: u64 = self.productive_sleeps.iter().sum();
        if total == 0 {
            return 0.0;
        }
        let deep: u64 = self.productive_sleeps[DEEP_SLEEP_FIRST_BUCKET..].iter().sum();
        deep as f64 / total as f64
    }
}

/// First bucket counted as a "deep" sleep (>= 320 µs).
const DEEP_SLEEP_FIRST_BUCKET: usize = 5;

// ============================================================================
// Consumer side: where the merge loop blocks
// ============================================================================

/// What the file the consumer is actually waiting on was doing.
///
/// The pool-wide shares below say what the *other* files were doing, which
/// turned out to answer a question nobody was asking: on a measured 780M-record
/// merge, 98% of files sat at their cap while the awaited file was starved only
/// 3% of the time. Knowing the rest of the pool is full says nothing about why
/// the one file that matters has not produced its next block.
///
/// These five states are mutually exclusive and, crucially, are observed at a
/// point where the consumer has *just* failed `try_pop_next`. So a non-empty
/// reorder buffer is proof the buffer holds the wrong serials, not a guess.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AwaitedState {
    /// The buffer holds blocks, but not the next serial, and a decompression is
    /// underway — the gap-filler is being produced right now.
    ReorderGapFilling,
    /// The buffer holds blocks, but not the next serial, and nothing is being
    /// decompressed for this file. The work exists and no worker has claimed
    /// it.
    ReorderGapStalled,
    /// Nothing buffered, but a decompression for this file is in flight.
    Decompressing,
    /// Nothing buffered and nothing in flight, but raw blocks are queued —
    /// read from disk and waiting for a worker to pick them up.
    RawQueued,
    /// Nothing anywhere and the reader is not at EOF: waiting on disk.
    Starved,
}

impl AwaitedState {
    /// Number of variants, for fixed-size counter arrays.
    pub(crate) const COUNT: usize = 5;

    /// Every variant, in declaration order.
    pub(crate) const ALL: [Self; Self::COUNT] = [
        Self::ReorderGapFilling,
        Self::ReorderGapStalled,
        Self::Decompressing,
        Self::RawQueued,
        Self::Starved,
    ];

    /// Short label for log output.
    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::ReorderGapFilling => "gap-filling",
            Self::ReorderGapStalled => "gap-stalled",
            Self::Decompressing => "decompressing",
            Self::RawQueued => "raw-queued",
            Self::Starved => "starved",
        }
    }

    /// Classify from a file's depths, observed at a point where the consumer
    /// has just failed to pop the serial it wants.
    ///
    /// That precondition is what makes a non-empty reorder buffer mean "holds
    /// the wrong serials" rather than merely "holds something", so this must
    /// only be called from the park path.
    #[must_use]
    pub(crate) fn classify(decomp_len: usize, in_flight: usize, raw_len: usize) -> Self {
        if decomp_len > 0 {
            if in_flight > 0 { Self::ReorderGapFilling } else { Self::ReorderGapStalled }
        } else if in_flight > 0 {
            Self::Decompressing
        } else if raw_len > 0 {
            Self::RawQueued
        } else {
            Self::Starved
        }
    }
}

// Keeps `COUNT` and `ALL` honest: the match is exhaustive, so a sixth variant
// fails to compile here instead of indexing past the end of `park_by_state`
// inside [`ConsumerTraceStats::record_park`](crate::merge_trace::ConsumerTraceStats)
// -- on a worker-facing path, where it surfaces only as a panic.
const _: () = {
    const fn assert_count(state: AwaitedState) -> usize {
        match state {
            AwaitedState::ReorderGapFilling => 0,
            AwaitedState::ReorderGapStalled => 1,
            AwaitedState::Decompressing => 2,
            AwaitedState::RawQueued => 3,
            AwaitedState::Starved => AwaitedState::COUNT - 1,
        }
    }
    assert!(assert_count(AwaitedState::Starved) == 4);
};

/// The other files' states at the moment the consumer parked.
///
/// Built by the caller, which owns the locks; keeping the observation out of
/// this module leaves the accounting pure and testable.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ParkCensus {
    /// Files whose reorder buffer is at `PHASE2_DECOMP_CAP`.
    pub(crate) capped: u32,
    /// Files with nothing buffered, nothing in flight, and a reader not at EOF.
    pub(crate) starved: u32,
    /// Files holding data but not yet at cap.
    pub(crate) working: u32,
    /// Files that have produced everything they ever will.
    pub(crate) drained: u32,
    /// Files whose state could not be read without blocking.
    pub(crate) contended: u32,
    /// What the awaited file itself was doing, when it could be read.
    pub(crate) awaited: Option<AwaitedState>,
}

impl ParkCensus {
    /// Files that are not finished.
    #[must_use]
    fn live(self) -> u32 {
        self.capped + self.starved + self.working + self.contended
    }
}

/// Consumer-side stall accounting.
///
/// Owned by the merge loop's `MainThreadChunkConsumer` and touched only by the
/// main thread, so no atomics: the consumer is single-threaded by construction
/// and paying for interlocked adds on its own counters would be pure waste.
#[derive(Debug)]
pub(crate) struct ConsumerStallTracker {
    /// Nanoseconds parked, per source.
    park_nanos_by_source: Vec<u64>,
    /// Parks, per source.
    parks_by_source: Vec<u64>,
    /// Calls to `advance_to_next_block`, i.e. attempts to pull a block.
    block_pulls: u64,
    /// Pulls that had to park at least once.
    stalled_pulls: u64,
    /// Total parks (a single pull may park more than once).
    parks: u64,
    /// Total nanoseconds parked.
    park_nanos: u64,
    /// Censuses taken, and the totals they observed.
    censuses: u64,
    census_capped: u64,
    census_starved: u64,
    census_contended: u64,
    census_live: u64,
    /// Awaited-file states observed, indexed by [`AwaitedState`] as `usize`.
    census_awaited: [u64; AwaitedState::COUNT],
    /// Countdown to the next census.
    census_countdown: u64,
}

/// Parks between censuses.
///
/// A census `try_lock`s every spill file, so it is the one part of this module
/// with a cost worth bounding. Sampling it is safe in a way that sampling the
/// park *timing* would not be: the shares it estimates are stable properties of
/// the pipeline's régime, whereas park durations are heavy-tailed and a sampled
/// mean of them would be badly biased. Prime, for the same aliasing reason as
/// the merge loop's 1021.
const CENSUS_INTERVAL: u64 = 61;

impl ConsumerStallTracker {
    /// Create a tracker for `num_sources` spill files.
    #[must_use]
    pub(crate) fn new(num_sources: usize) -> Self {
        Self {
            park_nanos_by_source: vec![0; num_sources],
            parks_by_source: vec![0; num_sources],
            block_pulls: 0,
            stalled_pulls: 0,
            parks: 0,
            park_nanos: 0,
            censuses: 0,
            census_capped: 0,
            census_starved: 0,
            census_contended: 0,
            census_live: 0,
            census_awaited: [0; AwaitedState::COUNT],
            census_countdown: 0,
        }
    }

    /// Restart the counters, discarding everything observed so far.
    ///
    /// The tracker is built with the consumer, which happens before the merge
    /// pulls one block per source to seed the loser tree. Those pulls are the
    /// likeliest of any to park -- every file is cold at merge start -- and they
    /// fall outside `loop_total`, which starts once the tree is built. Left in,
    /// they put `park_fraction` and `classify_stall` on a longer interval than
    /// the wall clock they are divided by, and on a short merge the K seeding
    /// parks can be most of what the verdict sees.
    ///
    /// Keeps the per-source vectors allocated; only their contents reset.
    pub(crate) fn restart(&mut self) {
        let Self {
            park_nanos_by_source,
            parks_by_source,
            block_pulls,
            stalled_pulls,
            parks,
            park_nanos,
            censuses,
            census_capped,
            census_starved,
            census_contended,
            census_live,
            census_awaited,
            census_countdown,
        } = self;
        park_nanos_by_source.fill(0);
        parks_by_source.fill(0);
        *block_pulls = 0;
        *stalled_pulls = 0;
        *parks = 0;
        *park_nanos = 0;
        *censuses = 0;
        *census_capped = 0;
        *census_starved = 0;
        *census_contended = 0;
        *census_live = 0;
        census_awaited.fill(0);
        *census_countdown = 0;
    }

    /// Note an attempt to pull the next block for a source.
    #[inline]
    pub(crate) fn note_block_pull(&mut self) {
        self.block_pulls += 1;
    }

    /// Whether the caller should census the files before this park.
    ///
    /// Consumes one step of the countdown, so call it exactly once per park.
    pub(crate) fn should_census(&mut self) -> bool {
        if self.census_countdown == 0 {
            self.census_countdown = CENSUS_INTERVAL - 1;
            true
        } else {
            self.census_countdown -= 1;
            false
        }
    }

    /// Record one park of `nanos` waiting on `source_id`.
    ///
    /// `first_park_of_pull` distinguishes a pull that stalled from the repeated
    /// parks of a single stalled pull, so `stall_rate` measures how often a
    /// block was not ready rather than how many spurious wake-ups it took.
    pub(crate) fn record_park(&mut self, source_id: usize, nanos: u64, first_park_of_pull: bool) {
        self.parks += 1;
        self.park_nanos += nanos;
        if first_park_of_pull {
            self.stalled_pulls += 1;
        }
        if let Some(slot) = self.park_nanos_by_source.get_mut(source_id) {
            *slot += nanos;
        }
        if let Some(slot) = self.parks_by_source.get_mut(source_id) {
            *slot += 1;
        }
    }

    /// Record a census taken at a park.
    pub(crate) fn record_census(&mut self, census: ParkCensus) {
        self.censuses += 1;
        self.census_capped += u64::from(census.capped);
        self.census_starved += u64::from(census.starved);
        self.census_contended += u64::from(census.contended);
        self.census_live += u64::from(census.live());
        if let Some(state) = census.awaited {
            self.census_awaited[state as usize] += 1;
        }
    }

    #[allow(
        clippy::cast_precision_loss,
        reason = "park counts and nanosecond totals stay far below 2^52"
    )]
    pub(crate) fn snapshot(&self) -> ConsumerStallReport {
        let live = self.census_live as f64;
        let share = |n: u64| if live > 0.0 { n as f64 / live } else { 0.0 };
        let top = self
            .park_nanos_by_source
            .iter()
            .enumerate()
            .max_by_key(|&(_, &nanos)| nanos)
            .map_or((0, 0), |(idx, &nanos)| (idx, nanos));
        ConsumerStallReport {
            block_pulls: self.block_pulls,
            stalled_pulls: self.stalled_pulls,
            parks: self.parks,
            park_secs: self.park_nanos as f64 / 1e9,
            capped_share: share(self.census_capped),
            starved_share: share(self.census_starved),
            contended_share: share(self.census_contended),
            awaited: AwaitedShares::from_counts(self.census_awaited),
            censuses: self.censuses,
            top_source: top.0,
            top_source_park_secs: top.1 as f64 / 1e9,
            sources_parked_on: self.parks_by_source.iter().filter(|&&n| n > 0).count(),
        }
    }
}

/// Read-only view of [`ConsumerStallTracker`], for logging.
#[derive(Debug, Clone, Copy)]
pub struct ConsumerStallReport {
    /// Attempts to pull a decompressed block.
    pub block_pulls: u64,
    /// Pulls where the block was not ready and the consumer had to park.
    pub stalled_pulls: u64,
    /// Total parks, including repeated parks within one pull.
    pub parks: u64,
    /// Total seconds the consumer spent parked. Exact, not sampled.
    pub park_secs: f64,
    /// Mean share of live files sitting at `PHASE2_DECOMP_CAP` at a park.
    pub capped_share: f64,
    /// Mean share of live files with nothing buffered at a park.
    pub starved_share: f64,
    /// Mean share of live files whose state could not be read at a park.
    pub contended_share: f64,
    /// What the awaited file was doing, as shares of censuses taken.
    pub awaited: AwaitedShares,
    /// Censuses taken.
    pub censuses: u64,
    /// Source the consumer parked on longest.
    pub top_source: usize,
    /// Seconds parked on that source.
    pub top_source_park_secs: f64,
    /// Distinct sources the consumer ever parked on.
    pub sources_parked_on: usize,
}

impl ConsumerStallReport {
    /// Whether the consumer ever pulled a block.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.block_pulls == 0
    }

    /// Share of block pulls that had to wait, in `[0, 1]`.
    #[must_use]
    #[allow(clippy::cast_precision_loss, reason = "pull counts stay far below 2^52")]
    pub fn stall_rate(self) -> f64 {
        if self.block_pulls == 0 {
            return 0.0;
        }
        self.stalled_pulls as f64 / self.block_pulls as f64
    }

    /// Parks per pull that had to wait.
    ///
    /// Ideally 1. Higher means the consumer is being woken and re-parked
    /// repeatedly within a single stall: every worker unparks the main thread on
    /// every block it publishes, for *any* file, but the consumer can only
    /// proceed when the one file it is waiting on advances. Each of those wakes
    /// costs a syscall and a re-lock of that file's reorder buffer, and buys
    /// nothing. A large ratio is its own coordination problem, distinct from the
    /// stall it is measured inside.
    #[must_use]
    #[allow(clippy::cast_precision_loss, reason = "park counts stay far below 2^52")]
    pub fn parks_per_stall(self) -> f64 {
        if self.stalled_pulls == 0 {
            return 0.0;
        }
        self.parks as f64 / self.stalled_pulls as f64
    }

    /// Share of total park time spent on the single worst source, in `[0, 1]`.
    ///
    /// Near `1 / num_sources` means the wait is spread evenly; much higher means
    /// one run is holding the merge up.
    #[must_use]
    pub fn top_source_share(self) -> f64 {
        if self.park_secs <= 0.0 { 0.0 } else { self.top_source_park_secs / self.park_secs }
    }
}

/// What the awaited file was doing, as shares of the censuses taken.
///
/// The shares sum to 1 when every census could read the awaited file.
#[derive(Debug, Default, Clone, Copy)]
pub struct AwaitedShares {
    /// Buffer holds the wrong serials; the gap-filler is being decompressed.
    pub reorder_gap_filling: f64,
    /// Buffer holds the wrong serials and nothing is being decompressed.
    pub reorder_gap_stalled: f64,
    /// Nothing buffered; a decompression is in flight.
    pub decompressing: f64,
    /// Raw blocks queued with no worker on them.
    pub raw_queued: f64,
    /// Nothing anywhere, reader not at EOF.
    pub starved: f64,
}

impl AwaitedShares {
    #[allow(
        clippy::cast_precision_loss,
        reason = "census counts stay far below 2^52 on any real sort"
    )]
    fn from_counts(counts: [u64; AwaitedState::COUNT]) -> Self {
        let total: u64 = counts.iter().sum();
        if total == 0 {
            return Self::default();
        }
        let share = |state: AwaitedState| counts[state as usize] as f64 / total as f64;
        Self {
            reorder_gap_filling: share(AwaitedState::ReorderGapFilling),
            reorder_gap_stalled: share(AwaitedState::ReorderGapStalled),
            decompressing: share(AwaitedState::Decompressing),
            raw_queued: share(AwaitedState::RawQueued),
            starved: share(AwaitedState::Starved),
        }
    }

    /// Share where the block the consumer needs exists but no worker is on it.
    ///
    /// Actionable: the pipeline has the data and is not applying capacity to
    /// it, which is a scheduling or discovery problem.
    #[must_use]
    pub fn unclaimed(self) -> f64 {
        self.reorder_gap_stalled + self.raw_queued
    }

    /// Share where a worker is already producing the needed block.
    ///
    /// Not actionable by scheduling: the wait is the per-block decompression
    /// cost, and only running further ahead or decompressing faster removes it.
    #[must_use]
    pub fn in_progress(self) -> f64 {
        self.reorder_gap_filling + self.decompressing
    }
}

/// The shape of a consumer stall — which of the candidate causes it fits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StallShape {
    /// The consumer barely parked; whatever it spends its loop on, it is not
    /// waiting for blocks.
    NotStalled,
    /// The consumer waits on a file with nothing anywhere in its pipeline —
    /// the block has not even been read from disk. The constraint is upstream
    /// of the pool: disk, or read concurrency.
    HeadOfLine,
    /// The block the consumer needs already exists, in raw or partially
    /// reordered form, and no worker is applying capacity to it. A scheduling
    /// and discovery problem, not a capacity one: workers are idle, or looking
    /// at the wrong file, or asleep.
    WorkUnclaimed,
    /// A worker is already producing the block the consumer needs. The wait is
    /// the per-block decompression cost itself, so only running further ahead
    /// (a deeper cap) or decompressing faster removes it.
    DecompressLatency,
    /// Enough of the pipeline's state was lock-contended that contention is the
    /// first thing to fix; the other shares are unreliable until it is.
    Contended,
    /// No candidate fits cleanly.
    Mixed,
}

/// Park share of the merge loop below which the consumer is not stalling.
const NOT_STALLED_PARK_FRACTION: f64 = 0.05;
/// Share of the awaited-file census one explanation must reach to be reported.
const AWAITED_DOMINANT_SHARE: f64 = 0.50;
/// Share of live files that must be unreadable to blame contention.
const CONTENDED_FILE_SHARE: f64 = 0.25;

/// Classify a consumer stall against the candidate causes.
///
/// `park_fraction` is exact park time over merge loop wall clock, which makes
/// the `NotStalled` gate trustworthy in a way the sampled fetch fraction is not.
/// `awaited` comes from the sampled park census.
///
/// # Why this classifies on the awaited file, not on the pool
///
/// It used to key off "most files at cap while the awaited file is starved",
/// on the theory that read-ahead was being misallocated across files. Measured
/// on a 780M-record merge, 98% of files were at cap and the awaited file was
/// starved 3% of the time — the pool-wide share was near-constant across arms
/// that stalled 8x differently, so it discriminated nothing. What the consumer
/// is waiting *for* does.
///
/// The useful split is whether the needed block exists yet, and if it does,
/// whether anyone is working on it. Those three cases have three different
/// fixes, and no pool-wide statistic distinguishes them.
#[must_use]
pub fn classify_stall(
    park_fraction: f64,
    contended_share: f64,
    awaited: AwaitedShares,
) -> StallShape {
    if park_fraction < NOT_STALLED_PARK_FRACTION {
        StallShape::NotStalled
    } else if contended_share >= CONTENDED_FILE_SHARE {
        StallShape::Contended
    } else if awaited.starved >= AWAITED_DOMINANT_SHARE {
        StallShape::HeadOfLine
    } else if awaited.unclaimed() >= AWAITED_DOMINANT_SHARE {
        StallShape::WorkUnclaimed
    } else if awaited.in_progress() >= AWAITED_DOMINANT_SHARE {
        StallShape::DecompressLatency
    } else {
        StallShape::Mixed
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    // ---- Phase2ScanTally / classify_scan -----------------------------------

    #[test]
    fn test_tally_counts_by_reason() {
        let mut tally = Phase2ScanTally::default();
        tally.note(Phase2Skip::ReadInProgress);
        tally.note(Phase2Skip::ReadInProgress);
        tally.note(Phase2Skip::Drained);
        assert_eq!(tally.count(Phase2Skip::ReadInProgress), 2);
        assert_eq!(tally.count(Phase2Skip::Drained), 1);
        assert_eq!(tally.count(Phase2Skip::RawFull), 0);
        assert_eq!(tally.total(), 3);
    }

    #[rstest]
    // Backpressure is the one outcome that says the pipeline is working, so it
    // wins over the read half's (correlated) complaint.
    #[case(PopSkip::DecompCapped, ReadSkip::RawFull, Phase2Skip::DecompCapped)]
    #[case(PopSkip::DecompCapped, ReadSkip::ReaderBusy, Phase2Skip::DecompCapped)]
    // A lost try_lock on either side means this worker's view is unreliable,
    // which matters more than what the other half happened to see.
    #[case(PopSkip::RawLockContended, ReadSkip::RawFull, Phase2Skip::LockContended)]
    #[case(PopSkip::DecompLockContended, ReadSkip::ReaderEof, Phase2Skip::LockContended)]
    #[case(PopSkip::RawEmpty, ReadSkip::RawLockContended, Phase2Skip::LockContended)]
    // Nothing buffered and a peer already fetching: waiting on disk, not on
    // buffer depth.
    #[case(PopSkip::RawEmpty, ReadSkip::ReaderBusy, Phase2Skip::ReadInProgress)]
    #[case(PopSkip::RawEmpty, ReadSkip::RawFull, Phase2Skip::RawFull)]
    #[case(PopSkip::RawEmpty, ReadSkip::ReaderEof, Phase2Skip::Drained)]
    fn test_combine_skip(
        #[case] pop: PopSkip,
        #[case] read: ReadSkip,
        #[case] expected: Phase2Skip,
    ) {
        assert_eq!(combine_skip(pop, read), expected);
    }

    /// Drained files are excluded from the denominator, so a scan late in the
    /// merge is judged on the files that could still have produced work.
    #[test]
    fn test_live_excludes_drained_files() {
        let mut tally = Phase2ScanTally::default();
        for _ in 0..80 {
            tally.note(Phase2Skip::Drained);
        }
        for _ in 0..6 {
            tally.note(Phase2Skip::LockContended);
        }
        assert_eq!(tally.total(), 86);
        assert_eq!(tally.live(), 6);
        assert_eq!(classify_scan(tally), ScanVerdict::Contended);
    }

    fn scan_of(pairs: &[(Phase2Skip, u32)]) -> Phase2ScanTally {
        let mut tally = Phase2ScanTally::default();
        for &(reason, n) in pairs {
            for _ in 0..n {
                tally.note(reason);
            }
        }
        tally
    }

    #[rstest]
    // The whole point of the counter: workers idle because every live file is
    // buffered to its cap. Raising the caps buys buffering, not throughput.
    #[case::all_files_at_cap(&[(Phase2Skip::DecompCapped, 86)], ScanVerdict::Backpressured)]
    // Read-side and decompress-side caps are the same régime for this purpose,
    // so they are summed rather than competing.
    #[case::caps_split_across_both_stages(
        &[(Phase2Skip::DecompCapped, 45), (Phase2Skip::RawFull, 41)],
        ScanVerdict::Backpressured
    )]
    // Workers queued behind reads their peers are already performing: the pool
    // is waiting on disk, and deeper buffers cannot help.
    #[case::everyone_waiting_on_a_peer_read(
        &[(Phase2Skip::ReadInProgress, 86)],
        ScanVerdict::IoWait
    )]
    // Contention is judged at a lower bar than the others — a quarter of the
    // pool's view being invisible is a defect even when something else is more
    // common.
    #[case::quarter_contended_beats_a_larger_io_share(
        &[(Phase2Skip::LockContended, 22), (Phase2Skip::ReadInProgress, 64)],
        ScanVerdict::Contended
    )]
    // 40/86 buffered and 46/86 waiting on a peer read: neither reaches the
    // dominance bar, and reporting either would be a guess.
    #[case::no_reason_dominates(
        &[
            (Phase2Skip::DecompCapped, 20),
            (Phase2Skip::RawFull, 20),
            (Phase2Skip::ReadInProgress, 46),
        ],
        ScanVerdict::Mixed
    )]
    #[case::merge_winding_down(&[(Phase2Skip::Drained, 86)], ScanVerdict::Drained)]
    #[case::nothing_observed(&[], ScanVerdict::Drained)]
    fn test_classify_scan(#[case] pairs: &[(Phase2Skip, u32)], #[case] expected: ScanVerdict) {
        assert_eq!(classify_scan(scan_of(pairs)), expected);
    }

    #[test]
    fn test_scan_stats_accumulate_across_scans() {
        let stats = Phase2ScanStats::default();
        stats.record_fruitless_scan(scan_of(&[(Phase2Skip::DecompCapped, 8)]));
        stats.record_fruitless_scan(scan_of(&[(Phase2Skip::DecompCapped, 8)]));
        stats.record_fruitless_scan(scan_of(&[(Phase2Skip::ReadInProgress, 8)]));
        let report = stats.snapshot();
        assert_eq!(report.scans, 3);
        assert_eq!(report.skips[Phase2Skip::DecompCapped as usize], 16);
        assert_eq!(report.skips[Phase2Skip::ReadInProgress as usize], 8);
        assert!((report.verdict_share(ScanVerdict::Backpressured) - 2.0 / 3.0).abs() < 1e-9);
        assert!((report.verdict_share(ScanVerdict::IoWait) - 1.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_empty_scan_report_is_reported_as_empty() {
        assert!(Phase2ScanStats::default().snapshot().is_empty());
    }

    // ---- Wake latency ------------------------------------------------------

    /// The bucket table must be exactly the backoff levels the worker loop can
    /// reach, or the histogram mislabels every row. Pinned against the loop's
    /// own constants so a change there fails here rather than silently.
    #[test]
    fn test_bucket_table_matches_the_worker_loop_backoff() {
        use crate::worker_pool::{MAX_BACKOFF_US, MIN_BACKOFF_US};
        let mut reachable = Vec::new();
        let mut us = MIN_BACKOFF_US;
        loop {
            reachable.push(us);
            if us == MAX_BACKOFF_US {
                break;
            }
            us = (us * 2).min(MAX_BACKOFF_US);
        }
        assert_eq!(reachable, WAKE_BUCKET_US.to_vec());
    }

    #[rstest]
    #[case(10, 0)]
    #[case(20, 1)]
    #[case(320, 5)]
    #[case(640, 6)]
    #[case(1000, 7)]
    // Jitter can push an actual sleep past its nominal level; it must not fall
    // off the end of the table.
    #[case(4000, 7)]
    fn test_wake_bucket(#[case] backoff_us: u64, #[case] expected: usize) {
        assert_eq!(wake_bucket(backoff_us), expected);
    }

    #[test]
    fn test_discovery_lag_is_half_the_sleep_that_found_work() {
        let stats = WakeLatencyStats::default();
        // 1000 sleeps at 1 ms that each found work: 1000 x 0.5 ms = 0.5 s.
        for _ in 0..1000 {
            stats.record_sleep(1000, 1_000_000);
            stats.record_productive_sleep(1000);
        }
        // Sleeps that found nothing cost the pipeline nothing and must not count.
        for _ in 0..1000 {
            stats.record_sleep(1000, 1_000_000);
        }
        let report = stats.snapshot();
        assert!((report.estimated_discovery_lag_secs() - 0.5).abs() < 1e-9);
        assert_eq!(report.sleeps[7], 2000);
    }

    #[test]
    fn test_deep_sleep_wake_share() {
        let stats = WakeLatencyStats::default();
        for _ in 0..3 {
            stats.record_productive_sleep(1000);
        }
        stats.record_productive_sleep(10);
        assert!((stats.snapshot().deep_sleep_wake_share() - 0.75).abs() < 1e-9);
    }

    /// Phase 1's sleeps must not be charged to the merge.
    #[test]
    fn test_wake_report_subtracts_a_phase_boundary_baseline() {
        let stats = WakeLatencyStats::default();
        for _ in 0..100 {
            stats.record_sleep(1000, 1_000_000);
            stats.record_productive_sleep(1000);
        }
        let at_merge_start = stats.snapshot();
        for _ in 0..10 {
            stats.record_sleep(1000, 1_000_000);
            stats.record_productive_sleep(1000);
        }
        let merge_only = stats.snapshot().since(at_merge_start);
        assert_eq!(merge_only.sleeps[7], 10, "only the merge's sleeps may be counted");
        assert!((merge_only.estimated_discovery_lag_secs() - 0.005).abs() < 1e-9);
        // The running total is still 110; `since` must not mutate the source.
        assert_eq!(stats.snapshot().sleeps[7], 110);
    }

    /// A counter that somehow went backwards must clamp rather than wrap into a
    /// nonsensical multi-century total.
    #[test]
    fn test_wake_report_since_saturates() {
        let stats = WakeLatencyStats::default();
        stats.record_sleep(1000, 1_000_000);
        let later = stats.snapshot();
        let empty = WakeLatencyStats::default().snapshot();
        assert_eq!(empty.since(later).sleeps[7], 0);
    }

    #[test]
    fn test_empty_wake_report_is_reported_as_empty() {
        let report = WakeLatencyStats::default().snapshot();
        assert!(report.is_empty());
        assert!((report.estimated_discovery_lag_secs() - 0.0).abs() < f64::EPSILON);
        assert!((report.deep_sleep_wake_share() - 0.0).abs() < f64::EPSILON);
    }

    // ---- Consumer stalls ---------------------------------------------------

    #[test]
    fn test_stall_rate_counts_pulls_not_parks() {
        let mut tracker = ConsumerStallTracker::new(4);
        for _ in 0..10 {
            tracker.note_block_pull();
        }
        // One pull that stalled and was woken spuriously twice before its block
        // arrived is one stalled pull, not three.
        tracker.record_park(0, 1_000_000, true);
        tracker.record_park(0, 1_000_000, false);
        tracker.record_park(0, 1_000_000, false);
        let report = tracker.snapshot();
        assert_eq!(report.parks, 3);
        assert_eq!(report.stalled_pulls, 1);
        assert!((report.stall_rate() - 0.1).abs() < 1e-9);
        assert!((report.park_secs - 0.003).abs() < 1e-9);
        // Two of the three parks bought nothing, which is the whole point of
        // keeping the two counts apart.
        assert!((report.parks_per_stall() - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_park_time_is_attributed_per_source() {
        let mut tracker = ConsumerStallTracker::new(3);
        tracker.note_block_pull();
        tracker.record_park(0, 1_000_000, true);
        tracker.note_block_pull();
        tracker.record_park(2, 9_000_000, true);
        let report = tracker.snapshot();
        assert_eq!(report.top_source, 2);
        assert_eq!(report.sources_parked_on, 2);
        assert!((report.top_source_share() - 0.9).abs() < 1e-9);
    }

    /// Census interval is a countdown, so the first park censuses and the next
    /// `CENSUS_INTERVAL - 1` do not.
    #[test]
    fn test_census_sampling_interval() {
        let mut tracker = ConsumerStallTracker::new(1);
        assert!(tracker.should_census(), "the first park must be censused");
        for _ in 0..CENSUS_INTERVAL - 1 {
            assert!(!tracker.should_census());
        }
        assert!(tracker.should_census(), "and again exactly CENSUS_INTERVAL parks later");

        let mut tracker = ConsumerStallTracker::new(1);
        let censused = (0..CENSUS_INTERVAL * 2).filter(|_| tracker.should_census()).count();
        assert_eq!(censused, 2);
    }

    #[test]
    fn test_census_shares_use_live_files() {
        let mut tracker = ConsumerStallTracker::new(10);
        tracker.record_census(ParkCensus {
            capped: 6,
            starved: 1,
            working: 1,
            drained: 90,
            contended: 2,
            awaited: Some(AwaitedState::Starved),
        });
        let report = tracker.snapshot();
        // 10 live files out of 100 observed: the 90 drained ones do not dilute.
        assert!((report.capped_share - 0.6).abs() < 1e-9);
        assert!((report.starved_share - 0.1).abs() < 1e-9);
        assert!((report.contended_share - 0.2).abs() < 1e-9);
        assert!((report.awaited.starved - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_awaited_states_are_grouped_by_what_can_be_done_about_them() {
        let mut tracker = ConsumerStallTracker::new(4);
        for state in [
            AwaitedState::ReorderGapStalled,
            AwaitedState::RawQueued,
            AwaitedState::ReorderGapFilling,
            AwaitedState::Decompressing,
        ] {
            tracker.record_census(ParkCensus { awaited: Some(state), ..ParkCensus::default() });
        }
        let awaited = tracker.snapshot().awaited;
        // Half the waits are on work nobody has claimed (fixable by scheduling)
        // and half on work already underway (fixable only by getting ahead).
        assert!((awaited.unclaimed() - 0.5).abs() < 1e-9);
        assert!((awaited.in_progress() - 0.5).abs() < 1e-9);
    }

    /// A census that could not read the awaited file must not be counted as an
    /// observation of it, or an unreadable file silently reads as a state.
    #[test]
    fn test_unreadable_awaited_file_is_not_counted() {
        let mut tracker = ConsumerStallTracker::new(4);
        tracker.record_census(ParkCensus { contended: 4, awaited: None, ..ParkCensus::default() });
        let report = tracker.snapshot();
        assert_eq!(report.censuses, 1);
        assert!((report.awaited.starved - 0.0).abs() < f64::EPSILON);
        assert!((report.awaited.unclaimed() - 0.0).abs() < f64::EPSILON);
    }

    /// Seeding the loser tree happens before the merge clock starts, so the
    /// stalls it incurs must not survive into the merge's report -- they would
    /// be divided by a wall time that never contained them. Those seeding pulls
    /// are also the likeliest to park, since every file is cold at merge start,
    /// so leaving them in biases the verdict in exactly one direction.
    #[test]
    fn test_restart_drops_everything_observed_before_the_merge_clock() {
        let mut tracker = ConsumerStallTracker::new(3);
        // One seeding pull per source, each parking on a cold file.
        for source in 0..3 {
            tracker.note_block_pull();
            tracker.record_park(source, 5_000_000, true);
        }
        tracker.record_census(ParkCensus {
            starved: 3,
            awaited: Some(AwaitedState::Starved),
            ..ParkCensus::default()
        });
        assert!(!tracker.snapshot().is_empty(), "the seeding stalls were recorded");

        tracker.restart();

        let report = tracker.snapshot();
        assert!(report.is_empty(), "the merge starts from nothing");
        assert_eq!(report.parks, 0);
        assert_eq!(report.block_pulls, 0);
        assert_eq!(report.stalled_pulls, 0);
        assert_eq!(report.censuses, 0);
        assert!((report.park_secs - 0.0).abs() < f64::EPSILON);
        assert!((report.top_source_share() - 0.0).abs() < f64::EPSILON);

        // And it still counts normally afterwards -- a reset that also broke
        // collection would pass every assertion above.
        tracker.note_block_pull();
        tracker.record_park(1, 2_000_000, true);
        let report = tracker.snapshot();
        assert_eq!(report.parks, 1);
        assert_eq!(report.block_pulls, 1);
        assert!((report.park_secs - 0.002).abs() < 1e-9);
    }

    #[test]
    fn test_empty_consumer_report_is_reported_as_empty() {
        let report = ConsumerStallTracker::new(4).snapshot();
        assert!(report.is_empty());
        assert!((report.stall_rate() - 0.0).abs() < f64::EPSILON);
        assert!((report.top_source_share() - 0.0).abs() < f64::EPSILON);
    }

    /// Build shares in declaration order: filling, stalled, decompressing,
    /// raw-queued, starved.
    fn awaited_shares(s: [f64; AwaitedState::COUNT]) -> AwaitedShares {
        AwaitedShares {
            reorder_gap_filling: s[0],
            reorder_gap_stalled: s[1],
            decompressing: s[2],
            raw_queued: s[3],
            starved: s[4],
        }
    }

    #[rstest]
    // The block has not been read from disk at all: nothing the pool does with
    // its buffers or its scheduler can help.
    #[case::block_not_read_yet(0.60, 0.02, [0.05, 0.05, 0.10, 0.10, 0.70], StallShape::HeadOfLine)]
    // The block exists -- as raw bytes, or as a gap with nothing being
    // decompressed -- and no worker has claimed it. Scheduling, not capacity.
    #[case::work_sitting_unclaimed(0.60, 0.02, [0.10, 0.40, 0.10, 0.35, 0.05], StallShape::WorkUnclaimed)]
    // A worker is already on it, so the wait is the decompression itself.
    #[case::already_being_produced(0.60, 0.02, [0.45, 0.05, 0.40, 0.05, 0.05], StallShape::DecompressLatency)]
    // Contention is checked first: while a quarter of the census is unreadable,
    // the other shares are estimates of an unobserved population.
    #[case::contention_masks_the_other_shares(0.60, 0.30, [0.45, 0.05, 0.40, 0.05, 0.05], StallShape::Contended)]
    // Measured on the storedout arm, where the consumer waits far less; the
    // gate has to be on exact park time, not the sampled fetch fraction.
    #[case::consumer_not_waiting(0.01, 0.30, [0.45, 0.05, 0.40, 0.05, 0.05], StallShape::NotStalled)]
    // Unclaimed and in-progress each fall short of dominance; reporting either
    // would be a guess.
    #[case::nothing_dominates(0.60, 0.02, [0.25, 0.25, 0.20, 0.20, 0.10], StallShape::Mixed)]
    fn test_classify_stall(
        #[case] park_fraction: f64,
        #[case] contended_share: f64,
        #[case] awaited: [f64; AwaitedState::COUNT],
        #[case] expected: StallShape,
    ) {
        assert_eq!(
            classify_stall(park_fraction, contended_share, awaited_shares(awaited)),
            expected
        );
    }

    /// The three groups partition the census, so a stall always has an answer
    /// even when no single state dominates.
    #[test]
    fn test_awaited_groups_partition_the_census() {
        let shares = awaited_shares([0.2, 0.3, 0.15, 0.25, 0.1]);
        let total = shares.starved + shares.unclaimed() + shares.in_progress();
        assert!((total - 1.0).abs() < 1e-9);
    }
}
