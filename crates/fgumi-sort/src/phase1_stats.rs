//! What Phase 1's serial ingest thread waits for.
//!
//! Phase 2 has a floor line ([`crate::merge_headroom`]) because its three
//! limits -- serial consumer, worker capacity, coordination -- imply unrelated
//! fixes and are routinely confused. Phase 1 has the same three limits and had
//! none of the accounting: its report was four wall-clock spans (read, sort,
//! spill write, merge) with no way to tell a thread that is busy from one that
//! is waiting.
//!
//! That gap matters more than it did for the merge. External sampling of a
//! 16-thread whole-genome sort put Phase 1 at **60% of total wall clock with its
//! main thread 91% busy** while all 16 cores averaged 5.3 -- so the phase is
//! bound by one thread's serial CPU, and every worker-side change is pushing on
//! a wall that is not there. In-process numbers should say that without needing
//! a `/proc` sampler attached from outside.
//!
//! # The two waits
//!
//! The ingest thread blocks in exactly two places, and neither was measured:
//!
//! 1. **Waiting for a decompressed block.** [`crate::read_ahead::PooledInputStream`]
//!    parks when the next serial it needs has not arrived. Because blocks are
//!    consumed in serial order through a reorder buffer, this can fire while
//!    other blocks are ready -- head-of-line blocking -- which is a different
//!    problem from an empty queue and is counted separately here.
//! 2. **Waiting for the previous spill to finish.** `drain_pending_spill` waits
//!    on the prior chunk's write handle between the read span ending and the
//!    in-memory sort starting, so that time lands in *no* phase bucket at all
//!    and shows up only as an unexplained residual against total wall clock.
//!
//! Timing is exact rather than sampled: both waits are milliseconds against a
//! ~30 ns clock read, so the clock is 5 orders of magnitude below the quantity
//! and costs nothing to read (the same argument [`crate::merge_trace`] makes for
//! the merge's block pull, where an exact timer cost 0.15%).

use std::sync::atomic::{AtomicU64, Ordering};

/// Why the ingest thread parked waiting for a block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ParkCause {
    /// Nothing decompressed was available at all: the workers are behind.
    Starved,
    /// Blocks were buffered, but not the serial the consumer needs next. The
    /// pipeline has work in flight and the consumer still cannot use it.
    HeadOfLine,
}

/// Counters for the ingest thread's waits, shared with the pool's input stream.
#[derive(Debug, Default)]
pub(crate) struct Phase1IngestStats {
    /// Times the ingest thread parked waiting for its next block.
    parks: AtomicU64,
    /// Total nanoseconds parked waiting for a block.
    park_nanos: AtomicU64,
    /// Parks where the reorder buffer held nothing at all.
    parks_starved: AtomicU64,
    /// Parks where the reorder buffer held blocks, but not the next serial.
    parks_head_of_line: AtomicU64,
    /// Times the ingest thread waited on the previous spill to complete.
    spill_waits: AtomicU64,
    /// Total nanoseconds waiting on a previous spill.
    spill_wait_nanos: AtomicU64,
}

impl Phase1IngestStats {
    /// Record one park and what the reorder buffer looked like when it happened.
    pub(crate) fn record_park(&self, elapsed_nanos: u64, cause: ParkCause) {
        self.parks.fetch_add(1, Ordering::Relaxed);
        self.park_nanos.fetch_add(elapsed_nanos, Ordering::Relaxed);
        match cause {
            ParkCause::Starved => self.parks_starved.fetch_add(1, Ordering::Relaxed),
            ParkCause::HeadOfLine => self.parks_head_of_line.fetch_add(1, Ordering::Relaxed),
        };
    }

    /// Record one wait on the previous chunk's spill write.
    pub(crate) fn record_spill_wait(&self, elapsed_nanos: u64) {
        self.spill_waits.fetch_add(1, Ordering::Relaxed);
        self.spill_wait_nanos.fetch_add(elapsed_nanos, Ordering::Relaxed);
    }

    /// A consistent view of the counters, for reporting.
    pub(crate) fn snapshot(&self) -> Phase1IngestReport {
        Phase1IngestReport {
            parks: self.parks.load(Ordering::Relaxed),
            park_secs: secs(self.park_nanos.load(Ordering::Relaxed)),
            parks_starved: self.parks_starved.load(Ordering::Relaxed),
            parks_head_of_line: self.parks_head_of_line.load(Ordering::Relaxed),
            spill_waits: self.spill_waits.load(Ordering::Relaxed),
            spill_wait_secs: secs(self.spill_wait_nanos.load(Ordering::Relaxed)),
        }
    }
}

#[allow(clippy::cast_precision_loss, reason = "nanosecond totals stay far below 2^52")]
fn secs(nanos: u64) -> f64 {
    nanos as f64 / 1_000_000_000.0
}

/// What the ingest thread waited for, as reported.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Phase1IngestReport {
    pub(crate) parks: u64,
    pub(crate) park_secs: f64,
    pub(crate) parks_starved: u64,
    pub(crate) parks_head_of_line: u64,
    pub(crate) spill_waits: u64,
    pub(crate) spill_wait_secs: f64,
}

impl Phase1IngestReport {
    /// Mean park in microseconds, or `None` when it never parked.
    ///
    /// The mean is the discriminant between "parked rarely and long" (a supply
    /// problem) and "parked constantly and briefly" (a handoff problem), which
    /// the merge campaign found to be the difference between a fixable stall and
    /// an unfixable one.
    pub(crate) fn mean_park_micros(&self) -> Option<f64> {
        #[allow(clippy::cast_precision_loss, reason = "park counts stay far below 2^52")]
        (self.parks > 0).then(|| self.park_secs * 1_000_000.0 / self.parks as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_park_causes_are_counted_separately() {
        let stats = Phase1IngestStats::default();
        stats.record_park(1_000, ParkCause::Starved);
        stats.record_park(3_000, ParkCause::HeadOfLine);
        stats.record_park(6_000, ParkCause::HeadOfLine);

        let report = stats.snapshot();
        assert_eq!(report.parks, 3);
        assert_eq!(report.parks_starved, 1);
        assert_eq!(report.parks_head_of_line, 2);
        // Starvation and head-of-line blocking have different fixes, so a report
        // that only totalled them would not distinguish the two.
        assert_eq!(report.parks_starved + report.parks_head_of_line, report.parks);
        assert!((report.park_secs - 10e-6).abs() < 1e-12, "got {}", report.park_secs);
    }

    #[test]
    fn test_the_spill_handoff_wait_is_counted_on_its_own() {
        let stats = Phase1IngestStats::default();
        stats.record_park(2_000_000_000, ParkCause::Starved);
        stats.record_spill_wait(3_000_000_000);

        let report = stats.snapshot();
        // The spill wait sits between the read span ending and the sort starting,
        // so it is in no phase bucket. It is reported separately from the park
        // rather than summed into it: one is the pool failing to keep the ingest
        // thread fed, the other is the previous chunk's write not being done, and
        // a single "blocked" total would hide which.
        assert!((report.spill_wait_secs - 3.0).abs() < 1e-9, "got {}", report.spill_wait_secs);
        assert!((report.park_secs - 2.0).abs() < 1e-9, "got {}", report.park_secs);
        assert_eq!(report.spill_waits, 1);
    }

    #[test]
    fn test_mean_park_is_absent_rather_than_zero_when_it_never_parked() {
        let idle = Phase1IngestStats::default().snapshot();
        // A reported 0 us mean would read as "parked, instantly", which is the
        // opposite of what never parking means.
        assert_eq!(idle.mean_park_micros(), None);

        let stats = Phase1IngestStats::default();
        stats.record_park(4_000_000, ParkCause::Starved);
        stats.record_park(6_000_000, ParkCause::Starved);
        let mean = stats.snapshot().mean_park_micros().expect("parked twice");
        assert!((mean - 5_000.0).abs() < 1e-6, "got {mean}");
    }
}
