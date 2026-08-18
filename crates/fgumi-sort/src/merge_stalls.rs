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
//! - every file's raw FIFO was at its read-ahead allowance -- [`PHASE2_RAW_CAP`],
//!   or [`PHASE2_STARVING_RAW_CAP`] for the file at the drain frontier
//!   (decompression is behind the read side, so read-ahead depth is the binding
//!   constraint),
//! - every file with data left already had a peer worker inside a disk read
//!   (waiting on I/O, which wants read concurrency rather than deeper buffers),
//!   or
//! - its `try_lock` calls lost races against other workers (contention).
//!
//! [`PHASE2_DECOMP_CAP`]: crate::worker_pool::PHASE2_DECOMP_CAP
//! [`PHASE2_RAW_CAP`]: crate::worker_pool::PHASE2_RAW_CAP
//! [`PHASE2_STARVING_RAW_CAP`]: crate::worker_pool::PHASE2_STARVING_RAW_CAP
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
//! - **How late workers discover work.** One wake path runs main thread to
//!   workers: when a file's reorder buffer drains, the consumer unparks a single
//!   worker through
//!   [`SharedPipelineState::wake_one_worker`](crate::worker_pool::SharedPipelineState).
//!   Everything else a worker could pick up arrives with no wake at all, and
//!   that is what this measures: a worker that has backed off to
//!   [`MAX_BACKOFF_US`](crate::worker_pool::MAX_BACKOFF_US) finds unannounced
//!   work up to a millisecond after it appears, and that delay is
//!   indistinguishable from I/O latency in every other counter.
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
//! Both [`classify_merge`](crate::merge_phases::classify_merge)'s thresholds and
//! [`classify_stall`]'s are fitted to measured runs, and the `measured_*` test
//! cases pin them to the specific merges that produced them — so moving a
//! threshold means confronting the evidence rather than editing a number.
//!
//! Those runs cover two regimes that behave very differently: a merge of
//! *disjoint* spill runs (which arises when the input is already sorted in the
//! requested order, and which idles half the pool) and a merge of *interleaved*
//! runs (which saturates it). A classifier calibrated on only one of them
//! misreads the other — notably, "the awaited block exists and no worker is on
//! it" describes a scheduling defect in the first regime and simple saturation
//! in the second, which is why `classify_stall` takes worker utilization.

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
    /// The raw FIFO is at this file's read-ahead allowance while the reorder
    /// buffer still has room, so decompression — not the consumer — is what this
    /// file is waiting on. The allowance is `PHASE2_RAW_CAP` for every file
    /// except the one at the drain frontier, which is allowed
    /// `PHASE2_STARVING_RAW_CAP`.
    ///
    /// Reported only when the decompress half is *not* capped; when both are
    /// full the reason is [`FullyBuffered`](Self::FullyBuffered), because a full
    /// raw FIFO is the downstream consequence of workers having stopped
    /// decompressing. A zero count here therefore does **not** mean the
    /// allowance was never reached — only that it was never the root reason.
    /// Misreading that cost real time once already.
    RawFull,
    /// Both caps are full: the consumer is behind, and the read-ahead behind it
    /// has backed up to its own limit. The file is as buffered as it can be.
    FullyBuffered,
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
    pub(crate) const COUNT: usize = 6;

    /// Every variant, in declaration order, for iteration and display.
    pub(crate) const ALL: [Self; Self::COUNT] = [
        Self::Drained,
        Self::DecompCapped,
        Self::RawFull,
        Self::FullyBuffered,
        Self::ReadInProgress,
        Self::LockContended,
    ];

    /// Short label for log output.
    pub(crate) fn label(self) -> &'static str {
        match self {
            Self::Drained => "drained",
            Self::DecompCapped => "decomp-capped",
            Self::RawFull => "raw-full",
            Self::FullyBuffered => "fully-buffered",
            Self::ReadInProgress => "read-in-progress",
            Self::LockContended => "lock-contended",
        }
    }
}

// Keeps `COUNT` honest: the match is exhaustive, so a seventh variant fails to
// compile here instead of indexing past the end of the counter arrays at run
// time.
const _: () = {
    const fn assert_count(skip: Phase2Skip) -> usize {
        match skip {
            Phase2Skip::Drained => 0,
            Phase2Skip::DecompCapped => 1,
            Phase2Skip::RawFull => 2,
            Phase2Skip::FullyBuffered => 3,
            Phase2Skip::ReadInProgress => 4,
            Phase2Skip::LockContended => Phase2Skip::COUNT - 1,
        }
    }
    assert!(assert_count(Phase2Skip::LockContended) == 5);
};

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
    /// The raw FIFO is at this file's read-ahead allowance.
    RawFull,
}

/// Reduce a file's two half-scan outcomes to the one reason worth counting.
///
/// Both halves decline independently and often for causally related reasons: a
/// file at its decompress cap has workers refusing to pull raw blocks, so its
/// raw FIFO fills to its own cap shortly after. Counting both would double-count
/// one condition, and reporting the downstream one would name the symptom.
///
/// So the reasons are ordered by causal depth, not severity. `DecompCapped`
/// means the consumer is behind; `FullyBuffered` is the same thing after the
/// read-ahead behind it has also backed up. Only when the decompress side has
/// room does `RawFull` mean what it says — that decompression is the laggard.
#[must_use]
pub(crate) fn combine_skip(pop: PopSkip, read: ReadSkip) -> Phase2Skip {
    match (pop, read) {
        (PopSkip::RawLockContended | PopSkip::DecompLockContended, _)
        | (_, ReadSkip::RawLockContended) => Phase2Skip::LockContended,
        (PopSkip::DecompCapped, ReadSkip::RawFull) => Phase2Skip::FullyBuffered,
        (PopSkip::DecompCapped, _) => Phase2Skip::DecompCapped,
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
    } else if share(Phase2Skip::DecompCapped)
        + share(Phase2Skip::FullyBuffered)
        + share(Phase2Skip::RawFull)
        >= DOMINANT_SHARE
    {
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

/// Which set of wake counters a wait belongs to.
///
/// A wait is filed under the phase the worker believed it was in when it *took*
/// the wait — a value decided once, at the top of the worker's iteration, and
/// carried through to the record. That is what makes the merge's share of these
/// counters race-free: a Phase 1 wait lands in `Other` however late the worker
/// gets around to recording it.
///
/// Subtracting a baseline snapshot taken at the phase boundary cannot promise
/// that, which is why it is no longer how the merge is scoped. The boundary can
/// fall between a worker's phase load and its record, and the wait then lands
/// *after* the baseline while still belonging to the phase before it — charging
/// the merge for a wait taken before it started, and leaving a productive wait in
/// the merge window whose matching sleep was subtracted away.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum WakePhase {
    /// Every phase other than the merge, principally Phase 1's read/sort/spill.
    Other = 0,
    /// Phase 2: the k-way merge, the window [`WakeLatencyStats::snapshot`] is
    /// asked for by the merge diagnostics.
    Merge = 1,
}

impl WakePhase {
    /// Number of variants, so the counter arrays indexed by `self as usize`
    /// cannot fall behind the enum. Sizing them with a literal instead lets a
    /// new variant compile and then index out of bounds inside
    /// [`WakeLatencyStats::record_sleep`] -- on a worker thread, where it
    /// surfaces only as a panicked sort worker.
    pub(crate) const COUNT: usize = 2;
}

// Keeps `COUNT` honest: the match is exhaustive, so a third variant fails to
// compile here instead of indexing past the end of the counter arrays at run
// time.
const _: () = {
    const fn assert_count(phase: WakePhase) -> usize {
        match phase {
            WakePhase::Other => 0,
            WakePhase::Merge => WakePhase::COUNT - 1,
        }
    }
    assert!(assert_count(WakePhase::Merge) == 1);
};

/// One phase's worth of [`WakeLatencyStats`] counters.
#[derive(Debug, Default)]
struct WakeCounters {
    /// Nanoseconds slept, by backoff level.
    idle_nanos: [AtomicU64; WAKE_BUCKET_US.len()],
    /// Sleeps, by backoff level.
    sleeps: [AtomicU64; WAKE_BUCKET_US.len()],
    /// Waits that were immediately followed by a successful step, bucketed by
    /// how long they *actually* lasted.
    productive_sleeps: [AtomicU64; WAKE_BUCKET_US.len()],
    /// Observed nanoseconds of those waits.
    ///
    /// Kept alongside the count because the count alone can only be turned into
    /// a duration via the bucket's nominal width, and since waits became
    /// `park_timeout` an `unpark` can end one at any point. See
    /// [`WakeLatencyStats::record_productive_sleep`].
    productive_nanos: [AtomicU64; WAKE_BUCKET_US.len()],
}

/// How long workers sleep, and how much of that sleep precedes finding work.
///
/// The second half is the point. Idle time alone cannot distinguish a worker
/// that slept because there was nothing to do from one that slept *through* work
/// becoming available, and only the latter is a cost the pipeline pays. A
/// successful step recorded against a 1 ms bucket says the work was found no
/// sooner than a millisecond after the worker last looked.
///
/// The pool exists for the whole sort, so a single running total would include
/// Phase 1's sleeps: on a measured run, 1815s of "deep-sleep worker-seconds"
/// against a merge that only had 8 x 305s = 2442 worker-seconds of capacity in
/// the first place. Counters are therefore kept per [`WakePhase`], so the merge
/// reads its own and nothing else's.
#[derive(Debug, Default)]
pub(crate) struct WakeLatencyStats {
    /// One set of counters per [`WakePhase`], indexed by its discriminant.
    by_phase: [WakeCounters; WakePhase::COUNT],
}

impl WakeLatencyStats {
    /// Record a sleep of `idle_nanos` at backoff level `backoff_us`, taken in `phase`.
    pub(crate) fn record_sleep(&self, phase: WakePhase, backoff_us: u64, idle_nanos: u64) {
        let counters = &self.by_phase[phase as usize];
        let b = wake_bucket(backoff_us);
        counters.idle_nanos[b].fetch_add(idle_nanos, Ordering::Relaxed);
        counters.sleeps[b].fetch_add(1, Ordering::Relaxed);
    }

    /// Record that a wait of `observed_nanos`, taken in `phase`, was followed by
    /// finding work.
    ///
    /// Bucketed by the observed duration, not by the backoff level that was
    /// requested. Since waits became `park_timeout`,
    /// [`SharedPipelineState::wake_one_worker`](crate::worker_pool::SharedPipelineState)
    /// can end a 1 ms wait after a few microseconds -- and charging that to the
    /// 1 ms bucket would report discovery lag the wake had just eliminated,
    /// making the pool look worst exactly where this path works best.
    pub(crate) fn record_productive_sleep(&self, phase: WakePhase, observed_nanos: u64) {
        let counters = &self.by_phase[phase as usize];
        let b = wake_bucket(observed_nanos / 1_000);
        counters.productive_sleeps[b].fetch_add(1, Ordering::Relaxed);
        counters.productive_nanos[b].fetch_add(observed_nanos, Ordering::Relaxed);
    }

    /// The counters for `phase` alone.
    pub(crate) fn snapshot(&self, phase: WakePhase) -> WakeLatencyReport {
        let counters = &self.by_phase[phase as usize];
        WakeLatencyReport {
            idle_nanos: std::array::from_fn(|i| counters.idle_nanos[i].load(Ordering::Relaxed)),
            sleeps: std::array::from_fn(|i| counters.sleeps[i].load(Ordering::Relaxed)),
            productive_sleeps: std::array::from_fn(|i| {
                counters.productive_sleeps[i].load(Ordering::Relaxed)
            }),
            productive_nanos: std::array::from_fn(|i| {
                counters.productive_nanos[i].load(Ordering::Relaxed)
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
    /// Waits immediately followed by a successful step, by observed duration.
    pub productive_sleeps: [u64; WAKE_BUCKET_US.len()],
    /// Observed nanoseconds of those waits, by the same bucketing.
    pub productive_nanos: [u64; WAKE_BUCKET_US.len()],
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
    /// A wait followed by a successful step means the work was there when the
    /// worker woke and may have arrived at any point during it, so half the
    /// wait is the expected lag under a uniform arrival assumption.
    ///
    /// Computed from the *observed* wait, not the requested backoff level. Since
    /// waits became `park_timeout`, an `unpark` from
    /// [`SharedPipelineState::wake_one_worker`](crate::worker_pool::SharedPipelineState)
    /// can cut a 1 ms wait to microseconds; charging the nominal level would
    /// report that truncated wait as a full millisecond of lag, i.e. attribute
    /// to discovery latency exactly the latency the wake removed.
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
        reason = "nanosecond totals stay far below 2^52 on any real sort"
    )]
    pub fn estimated_discovery_lag_secs(self) -> f64 {
        self.productive_nanos.iter().sum::<u64>() as f64 / 2.0 / 1e9
    }

    /// Worker-seconds spent waiting at the deepest *requested* backoff levels.
    ///
    /// The pool of idle time most of the discovery lag is drawn from, printed
    /// beside it for scale: a lag that is a rounding error against this is not
    /// worth chasing. Not a strict bound on
    /// [`Self::estimated_discovery_lag_secs`] — that sums productive waits at
    /// every observed duration, including ones a shallow backoff produced — and
    /// the two are bucketed differently on purpose: this by the level requested,
    /// so it describes the backoff policy, and the lag by the duration observed,
    /// so an unparked wait costs what it actually cost.
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
// Consumer park decomposition
// ============================================================================

/// One consumer park, split into the stages of fetching the block it wanted.
///
/// The four fields partition the park by construction — see [`split_park`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct ParkSegments {
    /// Park spent before any worker claimed the needed block: the cost of
    /// getting *a* worker onto it, whichever one arrives first.
    pub to_claim: u64,
    /// Park spent between the claim and the block being published: the read and
    /// decompress work itself.
    pub work: u64,
    /// Park spent after the block was published, before the consumer resumed:
    /// the consumer's own wake latency.
    pub to_resume: u64,
    /// Park not attributable to fetching this block — the honesty term.
    pub unattributed: u64,
}

/// Split a park into [`ParkSegments`].
///
/// `claim` and `publish` are `None` when that stamp did not land inside this
/// park: the block may have been claimed before the consumer parked, or not yet
/// claimed when it resumed. Both are ordinary, and neither may be allowed to
/// inflate a segment — an unstamped stage contributes to `unattributed`
/// instead, so the reader can see how much of park the model fails to explain.
///
/// Stamps come from other threads and are only `Relaxed`, so they can be
/// observed out of order relative to `park_start`. Every subtraction therefore
/// saturates and every segment is clamped to what remains, which is what keeps
/// the sum exactly equal to `resume - park_start` rather than merely close.
pub(crate) fn split_park(
    park_start: u64,
    claim: Option<u64>,
    publish: Option<u64>,
    resume: u64,
) -> ParkSegments {
    let total = resume.saturating_sub(park_start);
    let mut seg = ParkSegments::default();
    let mut spent = 0u64;

    // Each stage is measured from the later of its own start and the previous
    // stage's end, so a stamp that predates the park contributes zero rather
    // than borrowing time from a neighbour.
    let mut cursor = park_start;
    if let Some(claim) = claim {
        seg.to_claim = claim.saturating_sub(cursor).min(total - spent);
        spent += seg.to_claim;
        cursor = cursor.max(claim);
    }
    if let Some(publish) = publish {
        seg.work = publish.saturating_sub(cursor).min(total - spent);
        spent += seg.work;
        cursor = cursor.max(publish);
        seg.to_resume = resume.saturating_sub(cursor).min(total - spent);
        spent += seg.to_resume;
    }
    seg.unattributed = total - spent;
    seg
}

/// Where the merge consumer's park time actually goes.
///
/// Every earlier attempt to explain park was either a share of park *events* —
/// which hides a rare-but-long cause — or a sum over workers, which overcounts
/// because the consumer waits for the *first* worker to deliver while the sum
/// counts all of them. Measured on one 16-thread merge, the worker-side sum read
/// 96.7s of critical-path lag against 97.0s of park, then a change that removed
/// 62.9s of that lag moved park by 7.8s. These counters are consumer-side and
/// additive instead, so the segments cannot exceed the park they came from.
#[derive(Debug, Default)]
pub(crate) struct ParkAttribution {
    to_claim_nanos: AtomicU64,
    work_nanos: AtomicU64,
    to_resume_nanos: AtomicU64,
    unattributed_nanos: AtomicU64,
    parks: AtomicU64,
    /// Parks where no worker claimed the block during the park, so the whole
    /// wait was for a stage that had already started or already finished.
    unclaimed_parks: AtomicU64,
    /// Blocks the awaited file had ready when the consumer resumed, summed.
    ///
    /// Divided by [`Self::parks`] this is pipeline depth on the critical path. A
    /// mean near 1 means every block the consumer needs is fetched on demand and
    /// costs a full round trip, which is a different problem from a slow round
    /// trip and has different fixes.
    ready_on_resume: AtomicU64,
}

impl ParkAttribution {
    /// Record one park, split by [`split_park`].
    pub(crate) fn record(&self, seg: ParkSegments, claimed: bool, ready_on_resume: u64) {
        self.to_claim_nanos.fetch_add(seg.to_claim, Ordering::Relaxed);
        self.work_nanos.fetch_add(seg.work, Ordering::Relaxed);
        self.to_resume_nanos.fetch_add(seg.to_resume, Ordering::Relaxed);
        self.unattributed_nanos.fetch_add(seg.unattributed, Ordering::Relaxed);
        self.parks.fetch_add(1, Ordering::Relaxed);
        if !claimed {
            self.unclaimed_parks.fetch_add(1, Ordering::Relaxed);
        }
        self.ready_on_resume.fetch_add(ready_on_resume, Ordering::Relaxed);
    }

    /// Snapshot for logging.
    pub(crate) fn snapshot(&self) -> ParkAttributionReport {
        ParkAttributionReport {
            to_claim_nanos: self.to_claim_nanos.load(Ordering::Relaxed),
            work_nanos: self.work_nanos.load(Ordering::Relaxed),
            to_resume_nanos: self.to_resume_nanos.load(Ordering::Relaxed),
            unattributed_nanos: self.unattributed_nanos.load(Ordering::Relaxed),
            parks: self.parks.load(Ordering::Relaxed),
            unclaimed_parks: self.unclaimed_parks.load(Ordering::Relaxed),
            ready_on_resume: self.ready_on_resume.load(Ordering::Relaxed),
        }
    }
}

/// Read-only view of [`ParkAttribution`].
#[derive(Debug, Clone, Copy)]
pub struct ParkAttributionReport {
    /// Park spent waiting for any worker to claim the needed block.
    pub to_claim_nanos: u64,
    /// Park spent on the read and decompress themselves.
    pub work_nanos: u64,
    /// Park spent after publication, waiting for the consumer's own wake.
    pub to_resume_nanos: u64,
    /// Park the model does not explain.
    pub unattributed_nanos: u64,
    /// Parks recorded.
    pub parks: u64,
    /// Parks in which no claim landed.
    pub unclaimed_parks: u64,
    /// Summed blocks ready on the awaited file at resume.
    pub ready_on_resume: u64,
}

impl ParkAttributionReport {
    /// Whether anything was recorded.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.parks == 0
    }

    /// Total park accounted for, which must match the exact park clock.
    #[must_use]
    pub fn total_nanos(self) -> u64 {
        self.to_claim_nanos + self.work_nanos + self.to_resume_nanos + self.unattributed_nanos
    }

    /// Mean blocks ready on the awaited file when the consumer resumed.
    ///
    /// Near 1.0 means the critical path has no pipeline depth: each block is
    /// fetched on demand, so the merge pays one full round trip per block.
    #[must_use]
    #[expect(clippy::cast_precision_loss, reason = "counts are well within f64's exact range")]
    pub fn mean_ready_on_resume(self) -> f64 {
        if self.parks == 0 {
            return 0.0;
        }
        self.ready_on_resume as f64 / self.parks as f64
    }
}

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
    /// Share of live-file observations sitting at `PHASE2_DECOMP_CAP` at a
    /// park.
    ///
    /// Pooled across censuses -- summed over the numerator and over `live`,
    /// not averaged per park -- so a census taken while more files were still
    /// live weighs proportionally more. Live count falls as sources retire,
    /// which is exactly when the shares are least interesting.
    pub capped_share: f64,
    /// Share of live-file observations with nothing buffered at a park, pooled
    /// the same way as [`Self::capped_share`].
    pub starved_share: f64,
    /// Share of live-file observations whose state could not be read at a park,
    /// pooled the same way as [`Self::capped_share`].
    pub contended_share: f64,
    /// What the awaited file was doing, as shares of the censuses that could
    /// read it -- which is at most [`Self::censuses`], not necessarily all of
    /// them. See [`AwaitedShares`].
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

/// What the awaited file was doing, as shares of the censuses that could read
/// it.
///
/// A census that found the awaited file unreadable is excluded from the
/// denominator, so these sum to 1 whenever any observation was made -- and the
/// count they are shares of can be below [`ConsumerStallReport::censuses`].
/// Multiplying a share by `censuses` therefore does not recover a count.
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
    /// Shares of `counts.iter().sum()` -- the censuses that read the awaited
    /// file, which is what `counts` tallies. Censuses that could not read it
    /// were never tallied here and so are not in the denominator.
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
    /// The pool is saturated. The consumer waits because every worker is busy,
    /// which is the expected shape of a healthy CPU-bound merge rather than a
    /// defect. Checked before the awaited-file states, because those cannot
    /// tell "nobody picked this up" from "nobody was free to pick it up" — on
    /// two measured 780M-record merges at 98% utilization the awaited file read
    /// as 79% unclaimed, which would otherwise be reported as a scheduling bug.
    PoolSaturated,
    /// The consumer waits on a file with nothing anywhere in its pipeline —
    /// the block has not even been read from disk. The constraint is upstream
    /// of the pool: disk, or read concurrency.
    HeadOfLine,
    /// The block the consumer needs already exists, in raw or partially
    /// reordered form, and no worker is applying capacity to it *while capacity
    /// is available*. A scheduling and discovery problem, not a capacity one:
    /// workers are idle, or looking at the wrong file, or asleep.
    WorkUnclaimed,
    /// A worker is already producing the block the consumer needs, so the wait
    /// is the per-block decompression cost paid serially.
    ///
    /// Deliberately does not prescribe a deeper cap. Whether running further
    /// ahead is even possible depends on the reorder dwell measured alongside
    /// this: if blocks are consumed as fast as they are inserted, the buffer is
    /// a pass-through and its cap is not what the pipeline is running into. See
    /// [`crate::merge_trace::BlockLifecycleReport::reorder_is_pass_through`].
    DecompressLatency,
    /// Enough of the pipeline's state was lock-contended that contention is the
    /// first thing to fix; the other shares are unreliable until it is.
    Contended,
    /// No candidate fits cleanly.
    Mixed,
}

/// Park share of the merge loop below which the consumer is not stalling.
const NOT_STALLED_PARK_FRACTION: f64 = 0.05;
/// Worker utilization at or above which a stall is simply saturation.
///
/// Shared with [`crate::merge_phases::classify_merge`]'s `CpuBound` threshold on
/// purpose: the two classifiers must not disagree about whether the pool is
/// busy.
const SATURATED_UTILIZATION: f64 = 0.85;
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
///
/// `utilization` is the one pool-wide figure that is still required, and it
/// gates everything else. "The block exists and no worker is on it" describes
/// both a scheduling defect and a fully-busy pool, and only utilization tells
/// them apart. Two measured merges at 98% utilization show 79% of waits on
/// unclaimed work; calling that a scheduling bug would be wrong.
#[must_use]
pub fn classify_stall(
    park_fraction: f64,
    utilization: f64,
    contended_share: f64,
    awaited: AwaitedShares,
) -> StallShape {
    if park_fraction < NOT_STALLED_PARK_FRACTION {
        StallShape::NotStalled
    } else if contended_share >= CONTENDED_FILE_SHARE {
        StallShape::Contended
    } else if utilization >= SATURATED_UTILIZATION {
        StallShape::PoolSaturated
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

// ============================================================================
// Worker recruitment: why the consumer waits for a worker it already woke
// ============================================================================

/// Why nobody had already started on the block the consumer is about to wait for.
///
/// The three cases are not degrees of one problem -- they imply different fixes,
/// which is why they are counted apart:
///
/// - `SleeperAvailable`: capacity was idle and unused. Pure coordination loss;
///   the pool should have been on it.
/// - `AllBusyCompressing`: every worker was busy and output compression was
///   queued. Priority inversion -- `get_sort_priorities` puts `Compress` ahead
///   of `Phase2FileWork` whenever the compress queue is non-empty, so the block
///   the consumer is blocked on waits behind output work.
/// - `AllBusyMerging`: every worker was busy on merge work. Genuine capacity;
///   nothing to schedule better.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ParkSupply {
    /// At least one worker was parked.
    SleeperAvailable,
    /// Nobody parked, and output compression was queued.
    AllBusyCompressing,
    /// Nobody parked, and the compress queue was empty.
    AllBusyMerging,
}

impl ParkSupply {
    /// Number of variants, so the census arrays indexed by [`Self::index`]
    /// cannot fall behind the enum. Sizing them with a literal instead lets a
    /// new variant compile and then index out of bounds inside
    /// [`ParkSupplyCensus::record`] -- on a worker-adjacent path, where it
    /// surfaces only as a panicked sort worker.
    pub(crate) const COUNT: usize = 3;

    /// Stable index for the census arrays.
    fn index(self) -> usize {
        match self {
            Self::SleeperAvailable => 0,
            Self::AllBusyCompressing => 1,
            Self::AllBusyMerging => Self::COUNT - 1,
        }
    }
}

// Keeps `COUNT` honest: the match is exhaustive, so a fourth variant fails to
// compile here instead of indexing past the end of the census arrays at run
// time.
const _: () = {
    const fn assert_count(supply: ParkSupply) -> usize {
        match supply {
            ParkSupply::SleeperAvailable => 0,
            ParkSupply::AllBusyCompressing => 1,
            ParkSupply::AllBusyMerging => ParkSupply::COUNT - 1,
        }
    }
    assert!(assert_count(ParkSupply::AllBusyMerging) == 2);
};

/// Classify the pool's state at the instant the consumer parks.
pub(crate) fn classify_park_supply(parked_workers: usize, compress_depth: usize) -> ParkSupply {
    // A sleeper outranks a queued compress deliberately: the two fixes are not
    // interchangeable. Waking an idle worker costs nothing, while reordering
    // priorities trades output throughput for merge latency, so a park with both
    // conditions true belongs in the cheaper bucket.
    if parked_workers > 0 {
        ParkSupply::SleeperAvailable
    } else if compress_depth > 0 {
        ParkSupply::AllBusyCompressing
    } else {
        ParkSupply::AllBusyMerging
    }
}

/// Park time and park counts split by [`ParkSupply`].
#[derive(Debug, Default)]
pub(crate) struct ParkSupplyCensus {
    counts: [AtomicU64; ParkSupply::COUNT],
    nanos: [AtomicU64; ParkSupply::COUNT],
}

impl ParkSupplyCensus {
    /// Attribute one park of `nanos` to `supply`.
    pub(crate) fn record(&self, supply: ParkSupply, nanos: u64) {
        let i = supply.index();
        self.counts[i].fetch_add(1, Ordering::Relaxed);
        self.nanos[i].fetch_add(nanos, Ordering::Relaxed);
    }

    /// Counts and nanoseconds per class, indexed as [`ParkSupply::index`].
    pub(crate) fn snapshot(&self) -> ParkSupplyReport {
        ParkSupplyReport {
            counts: std::array::from_fn(|i| self.counts[i].load(Ordering::Relaxed)),
            nanos: std::array::from_fn(|i| self.nanos[i].load(Ordering::Relaxed)),
        }
    }
}

/// Read-only view of [`ParkSupplyCensus`].
#[derive(Debug, Clone, Copy, Default)]
pub struct ParkSupplyReport {
    /// Parks per class: sleeper-available, all-busy-compressing, all-busy-merging.
    pub counts: [u64; ParkSupply::COUNT],
    /// Park nanoseconds per class, same order.
    pub nanos: [u64; ParkSupply::COUNT],
}

impl ParkSupplyReport {
    /// Total parks censused.
    #[must_use]
    pub fn total_parks(self) -> u64 {
        self.counts.iter().sum()
    }

    /// Total park nanoseconds censused.
    #[must_use]
    pub fn total_nanos(self) -> u64 {
        self.nanos.iter().sum()
    }
}

/// Takes a predicate rather than a slice so the wake path can read the atomics
/// directly: this runs on every wake -- millions of them -- and materializing a
/// `Vec<bool>` there would allocate on the hottest coordination path in the
/// merge. Production and tests call the same function.
pub(crate) fn first_parked_from<F>(cursor: usize, width: usize, is_parked: F) -> Option<usize>
where
    F: Fn(usize) -> bool,
{
    if width == 0 {
        return None;
    }
    (0..width).map(|offset| (cursor + offset) % width).find(|&i| is_parked(i))
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
    // Both caps full is its own state: the consumer is behind AND the read-ahead
    // behind it has backed up. Reporting this as plain `RawFull` would name the
    // symptom; reporting it as `DecompCapped` would lose that the FIFO filled.
    #[case(PopSkip::DecompCapped, ReadSkip::RawFull, Phase2Skip::FullyBuffered)]
    // Backpressure is the one outcome that says the pipeline is working, so it
    // wins over the read half's other (correlated) complaints.
    #[case(PopSkip::DecompCapped, ReadSkip::ReaderBusy, Phase2Skip::DecompCapped)]
    #[case(PopSkip::DecompCapped, ReadSkip::ReaderEof, Phase2Skip::DecompCapped)]
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
    // `FullyBuffered` is both caps at once, so it is summed with them rather
    // than falling through to `Mixed` — a guess, about a state the code can
    // name exactly.
    #[case::every_file_fully_buffered(
        &[(Phase2Skip::FullyBuffered, 86)],
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
            stats.record_sleep(WakePhase::Merge, 1000, 1_000_000);
            stats.record_productive_sleep(WakePhase::Merge, 1_000_000);
        }
        // Sleeps that found nothing cost the pipeline nothing and must not count.
        for _ in 0..1000 {
            stats.record_sleep(WakePhase::Merge, 1000, 1_000_000);
        }
        let report = stats.snapshot(WakePhase::Merge);
        assert!((report.estimated_discovery_lag_secs() - 0.5).abs() < 1e-9);
        assert_eq!(report.sleeps[7], 2000);
    }

    /// The four segments must partition the park exactly, whatever order the
    /// stamps arrive in.
    ///
    /// This is the property that makes the decomposition trustworthy where the
    /// worker-side lag sum was not: a segment cannot be inflated without another
    /// shrinking, so no single number can be read as larger than the park it
    /// came from.
    #[rstest::rstest]
    // Ordinary case: claimed then published inside the park.
    #[case(1_000, Some(1_100), Some(1_150), 1_200, (100, 50, 50, 0))]
    // Claimed but not yet published when the consumer resumed: no work segment,
    // and the remainder is unattributed rather than silently folded into work.
    #[case(1_000, Some(1_100), None, 1_200, (100, 0, 0, 100))]
    // Already claimed before the park: nothing to wait for a worker on.
    #[case(1_000, None, Some(1_050), 1_200, (0, 50, 150, 0))]
    // Neither stamp: the block was already available, so the park is not
    // attributable to fetching it at all.
    #[case(1_000, None, None, 1_200, (0, 0, 0, 200))]
    // Both stamps predate the park: the block was already published when the
    // consumer parked, so the whole wait is the consumer failing to notice --
    // charged to `to_resume`, not to a fetch that had already finished. Stamps
    // arrive from other threads under `Relaxed`, so this ordering is reachable
    // and must not yield a negative segment or a sum over the park.
    #[case(1_000, Some(900), Some(950), 1_200, (0, 0, 200, 0))]
    fn test_park_segments_always_partition_the_park(
        #[case] park_start: u64,
        #[case] claim: Option<u64>,
        #[case] publish: Option<u64>,
        #[case] resume: u64,
        #[case] want: (u64, u64, u64, u64),
    ) {
        let (want_to_claim, want_work, want_to_resume, want_unattributed) = want;
        let seg = split_park(park_start, claim, publish, resume);
        assert_eq!(seg.to_claim, want_to_claim, "to_claim");
        assert_eq!(seg.work, want_work, "work");
        assert_eq!(seg.to_resume, want_to_resume, "to_resume");
        assert_eq!(seg.unattributed, want_unattributed, "unattributed");
        assert_eq!(
            seg.to_claim + seg.work + seg.to_resume + seg.unattributed,
            resume - park_start,
            "segments must sum to the measured park"
        );
    }

    proptest::proptest! {
        /// Property: the segments partition the park for *any* arrangement of
        /// stamps, including ones that predate the park or arrive reversed.
        ///
        /// This is the whole basis for trusting the decomposition over the
        /// worker-side lag sum it replaces, and the case list cannot cover the
        /// orderings that `Relaxed` loads across threads make reachable.
        #[test]
        fn prop_park_segments_partition_the_park(
            park_start in 0u64..1_000_000,
            span in 1u64..1_000_000,
            claim_off in proptest::option::of(-500_000i64..1_500_000),
            publish_off in proptest::option::of(-500_000i64..1_500_000),
        ) {
            let resume = park_start + span;
            let stamp = |off: Option<i64>| {
                off.map(|o| u64::try_from(i64::try_from(park_start).unwrap_or(i64::MAX) + o)
                    .unwrap_or(0))
            };
            let seg = split_park(park_start, stamp(claim_off), stamp(publish_off), resume);
            proptest::prop_assert_eq!(
                seg.to_claim + seg.work + seg.to_resume + seg.unattributed,
                span,
                "segments must sum to the park exactly"
            );
        }
    }

    #[test]
    fn test_deep_sleep_wake_share() {
        let stats = WakeLatencyStats::default();
        for _ in 0..3 {
            stats.record_productive_sleep(WakePhase::Merge, 1_000_000);
        }
        stats.record_productive_sleep(WakePhase::Merge, 10_000);
        assert!((stats.snapshot(WakePhase::Merge).deep_sleep_wake_share() - 0.75).abs() < 1e-9);
    }

    /// A wait cut short by `wake_one_worker` must be charged what it cost, not
    /// what it asked for.
    ///
    /// Waits are `park_timeout`, so an unpark can end a 1 ms backoff after
    /// microseconds. Bucketing by the requested level would book each of these
    /// as a full millisecond and report ~0.5 s of discovery lag — attributing to
    /// late discovery precisely the latency the wake removed, so the metric
    /// would look worst exactly where this path works best.
    #[test]
    fn test_discovery_lag_uses_the_observed_wait_not_the_requested_backoff() {
        let stats = WakeLatencyStats::default();
        // 1000 waits that asked for 1 ms but were unparked after 5 us.
        for _ in 0..1000 {
            stats.record_sleep(WakePhase::Merge, 1000, 5_000);
            stats.record_productive_sleep(WakePhase::Merge, 5_000);
        }
        let report = stats.snapshot(WakePhase::Merge);

        // 1000 x 5us / 2 = 2.5 ms, not the 500 ms the nominal level implies.
        assert!(
            (report.estimated_discovery_lag_secs() - 0.0025).abs() < 1e-9,
            "got {}",
            report.estimated_discovery_lag_secs()
        );
        assert_eq!(report.productive_sleep_count(), 1000);
        assert!(
            (report.deep_sleep_wake_share() - 0.0).abs() < f64::EPSILON,
            "a 5us wait is not a deep sleep, whatever backoff level requested it"
        );
    }

    /// Phase 1's sleeps must not be charged to the merge.
    #[test]
    fn test_the_merge_report_holds_only_the_merges_own_waits() {
        let stats = WakeLatencyStats::default();
        for _ in 0..100 {
            stats.record_sleep(WakePhase::Other, 1000, 1_000_000);
            stats.record_productive_sleep(WakePhase::Other, 1_000_000);
        }
        for _ in 0..10 {
            stats.record_sleep(WakePhase::Merge, 1000, 1_000_000);
            stats.record_productive_sleep(WakePhase::Merge, 1_000_000);
        }
        let merge_only = stats.snapshot(WakePhase::Merge);
        assert_eq!(merge_only.sleeps[7], 10, "only the merge's sleeps may be counted");
        assert!((merge_only.estimated_discovery_lag_secs() - 0.005).abs() < 1e-9);
        // Phase 1's counters are untouched and still readable in their own right.
        assert_eq!(stats.snapshot(WakePhase::Other).sleeps[7], 100);
    }

    /// The race that a phase-boundary baseline could not close: a worker loads
    /// `PHASE1`, the main thread flips the pool to `PHASE2` (which is where the
    /// baseline used to be pinned), and only *then* does the worker record the
    /// wait it took back in Phase 1.
    ///
    /// Order is the whole test — every record here happens after the transition,
    /// which is exactly when a baseline subtraction would have swept the Phase 1
    /// wait into the merge. Filing by the phase carried from the worker's load
    /// makes the ordering irrelevant, so this holds for any interleaving rather
    /// than merely for the one a threaded test happened to produce.
    #[test]
    fn test_a_phase1_wait_recorded_after_the_merge_began_stays_out_of_the_merge() {
        let stats = WakeLatencyStats::default();

        // The merge has begun and taken one short wait of its own.
        stats.record_sleep(WakePhase::Merge, 10, 5_000);
        stats.record_productive_sleep(WakePhase::Merge, 5_000);

        // A worker that loaded PHASE1 before the flip now records its wait.
        stats.record_sleep(WakePhase::Other, 1000, 1_000_000);
        stats.record_productive_sleep(WakePhase::Other, 1_000_000);

        let merge = stats.snapshot(WakePhase::Merge);
        assert_eq!(
            merge.idle_nanos.iter().sum::<u64>(),
            5_000,
            "the Phase 1 wait must not appear in the merge's idle time"
        );
        assert_eq!(
            merge.productive_nanos.iter().sum::<u64>(),
            5_000,
            "nor in the merge's productive waits, which is what discovery lag is derived from"
        );
        // 5us / 2, not the 0.5 s the Phase 1 wait would have implied.
        assert!((merge.estimated_discovery_lag_secs() - 0.000_002_5).abs() < 1e-12);
    }

    #[test]
    fn test_empty_wake_report_is_reported_as_empty() {
        let report = WakeLatencyStats::default().snapshot(WakePhase::Merge);
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
    // Measured: the degenerate coord->coord cell (run f080bff). Half the pool
    // idle, the awaited block already being produced.
    #[case::measured_disjoint_runs(0.79, 0.55, 0.00, [0.11, 0.02, 0.73, 0.11, 0.03], StallShape::DecompressLatency)]
    // Measured: the same records from interleaved runs, on a saturated pool.
    // The awaited file reads as 79% unclaimed, which is NOT a scheduling defect
    // -- there was simply no free worker. Utilization is the only input that
    // separates this from `WorkUnclaimed`, which is why it is an argument.
    #[case::measured_interleaved_saturated(0.34, 0.98, 0.00, [0.00, 0.00, 0.08, 0.79, 0.13], StallShape::PoolSaturated)]
    #[case::measured_queryname_saturated(0.25, 0.98, 0.00, [0.00, 0.00, 0.07, 0.78, 0.15], StallShape::PoolSaturated)]
    // Same awaited shares as the two above but with capacity to spare: now the
    // unclaimed work really is a scheduling problem.
    #[case::unclaimed_with_idle_capacity(0.60, 0.40, 0.02, [0.10, 0.40, 0.10, 0.35, 0.05], StallShape::WorkUnclaimed)]
    // The block has not been read from disk at all: nothing the pool does with
    // its buffers or its scheduler can help.
    #[case::block_not_read_yet(0.60, 0.40, 0.02, [0.05, 0.05, 0.10, 0.10, 0.70], StallShape::HeadOfLine)]
    // Contention is checked before utilization: while a quarter of the census is
    // unreadable, every other share is an estimate of an unobserved population.
    #[case::contention_masks_the_other_shares(0.60, 0.98, 0.30, [0.45, 0.05, 0.40, 0.05, 0.05], StallShape::Contended)]
    // Exact park time gates the whole thing, not the sampled fetch fraction.
    #[case::consumer_not_waiting(0.01, 0.30, 0.30, [0.45, 0.05, 0.40, 0.05, 0.05], StallShape::NotStalled)]
    #[case::nothing_dominates(0.60, 0.40, 0.02, [0.25, 0.25, 0.20, 0.20, 0.10], StallShape::Mixed)]
    fn test_classify_stall(
        #[case] park_fraction: f64,
        #[case] utilization: f64,
        #[case] contended_share: f64,
        #[case] awaited: [f64; AwaitedState::COUNT],
        #[case] expected: StallShape,
    ) {
        assert_eq!(
            classify_stall(park_fraction, utilization, contended_share, awaited_shares(awaited)),
            expected
        );
    }

    /// The saturation threshold must agree with `classify_merge`'s, or the two
    /// verdicts in one log can contradict each other.
    #[test]
    fn test_saturation_threshold_agrees_with_classify_merge() {
        use crate::merge_phases::{MergeVerdict, classify_merge};
        assert_eq!(classify_merge(SATURATED_UTILIZATION, 0.9), MergeVerdict::CpuBound);
        assert_eq!(
            classify_stall(
                0.9,
                SATURATED_UTILIZATION,
                0.0,
                awaited_shares([0.0, 0.0, 0.0, 1.0, 0.0])
            ),
            StallShape::PoolSaturated
        );
    }

    /// Agreeing on the threshold is only half the contract -- both verdicts have
    /// to be computed from the *same* utilization figure, which means both
    /// loggers must divide worker busy time by the same pair.
    ///
    /// Dividing by the pool width instead of the Phase 2 active cap is not a
    /// rounding difference: on a run whose Phase 1 is wider than its Phase 2 it
    /// moves the number across `SATURATED_UTILIZATION`, so one log reports a
    /// saturated pool while the other blames scheduling for the same stall.
    #[test]
    fn test_pool_width_denominator_contradicts_the_active_cap_verdict() {
        use crate::merge_phases::{MergePhaseBreakdown, MergeVerdict, classify_merge};

        // A pool 8 threads wide, capped to 3 for Phase 2, whose workers were
        // busy 26s across a 10s merge -- 2.6x the wall clock, so the three
        // threads allowed to help were saturated.
        let breakdown = MergePhaseBreakdown {
            read: (6.0, 100),
            decompress: (20.0, 100),
            output_compress: (0.0, 1),
            spill_compress: (0.0, 0),
        };
        let (merge_wall, active_workers, pool_width) = (10.0, 3, 8);
        let awaited = awaited_shares([0.0, 0.0, 0.08, 0.79, 0.13]);

        let by_cap = breakdown.worker_utilization(merge_wall, active_workers).unwrap();
        assert!(by_cap >= SATURATED_UTILIZATION, "26s over 3 x 10s is saturated, got {by_cap}");
        assert_eq!(classify_merge(by_cap, 0.1), MergeVerdict::CpuBound);
        assert_eq!(classify_stall(0.9, by_cap, 0.0, awaited), StallShape::PoolSaturated);

        let by_width = breakdown.worker_utilization(merge_wall, pool_width).unwrap();
        assert!(by_width < SATURATED_UTILIZATION, "the wrong denominator understates it");
        assert_ne!(
            classify_stall(0.9, by_width, 0.0, awaited),
            StallShape::PoolSaturated,
            "the pool width turns the same saturated merge into a scheduling defect"
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
    // ========================================================================
    // Worker recruitment
    // ========================================================================

    /// The three classes imply three different fixes, so the boundaries matter
    /// more than the counts: a parked worker means coordination lost the block,
    /// a queued compress means priority did, and neither means capacity did.
    #[rstest]
    #[case::one_sleeper_is_coordination_loss(1, 0, ParkSupply::SleeperAvailable)]
    #[case::a_sleeper_outranks_a_full_compress_queue(3, 9, ParkSupply::SleeperAvailable)]
    #[case::nobody_free_with_compress_queued_is_priority(0, 5, ParkSupply::AllBusyCompressing)]
    #[case::nobody_free_and_nothing_queued_is_capacity(0, 0, ParkSupply::AllBusyMerging)]
    fn test_park_supply_separates_coordination_from_priority_from_capacity(
        #[case] parked_workers: usize,
        #[case] compress_depth: usize,
        #[case] expected: ParkSupply,
    ) {
        assert_eq!(classify_park_supply(parked_workers, compress_depth), expected);
    }

    /// A sleeper is reported even when compression is also backed up, because the
    /// fixes are not interchangeable: waking the sleeper costs nothing, whereas
    /// reordering priorities trades output throughput for merge latency.
    #[test]
    fn test_a_sleeper_is_not_masked_by_a_busy_compress_queue() {
        assert_eq!(classify_park_supply(1, 100), ParkSupply::SleeperAvailable);
    }

    #[test]
    fn test_census_attributes_parks_and_time_by_class() {
        let census = ParkSupplyCensus::default();
        census.record(ParkSupply::SleeperAvailable, 100);
        census.record(ParkSupply::SleeperAvailable, 50);
        census.record(ParkSupply::AllBusyCompressing, 700);
        census.record(ParkSupply::AllBusyMerging, 5);

        let report = census.snapshot();
        assert_eq!(report.counts, [2, 1, 1]);
        assert_eq!(report.nanos, [150, 700, 5]);
        assert_eq!(report.total_parks(), 4);
        assert_eq!(report.total_nanos(), 855);
    }

    #[test]
    fn test_census_starts_empty() {
        let report = ParkSupplyCensus::default().snapshot();
        assert_eq!(report.total_parks(), 0);
        assert_eq!(report.total_nanos(), 0);
    }

    #[test]
    fn test_first_parked_starts_at_the_cursor_and_wraps() {
        let parked = [false, false, true, false];
        let at = |i: usize| parked[i];
        assert_eq!(first_parked_from(0, 4, at), Some(2), "scans forward from the cursor");
        assert_eq!(first_parked_from(3, 4, at), Some(2), "and wraps around to find it");
        assert_eq!(first_parked_from(2, 4, at), Some(2), "the cursor itself counts");
    }

    /// The point of the counter this feeds: when nobody is parked a wake has
    /// nowhere useful to go, and must be recorded as wasted rather than silently
    /// spent on a running worker.
    #[test]
    fn test_first_parked_reports_none_when_every_worker_is_running() {
        assert_eq!(first_parked_from(0, 3, |_| false), None);
    }

    /// Wakes stay inside the active window, exactly as `wake_target` does -- a
    /// capped worker will not take Phase 2 work, so waking it is the same lost
    /// wake by another route.
    #[test]
    fn test_first_parked_ignores_workers_outside_the_active_limit() {
        let parked = [false, false, true, true];
        let at = |i: usize| parked[i];
        assert_eq!(
            first_parked_from(0, 2, at),
            None,
            "workers 2 and 3 are parked but capped out of Phase 2"
        );
        assert_eq!(first_parked_from(0, 3, at), Some(2));
    }

    #[test]
    fn test_first_parked_tolerates_an_empty_pool() {
        assert_eq!(first_parked_from(0, 0, |_| true), None);
        assert_eq!(first_parked_from(5, 0, |_| true), None, "a zero limit admits nobody");
    }
}
