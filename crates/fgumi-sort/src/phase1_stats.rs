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

/// One record in this many is timed for the sub-phase partition.
///
/// Prime, so the sampled set cannot align with any periodic structure in the
/// input (read groups, tile boundaries, alternating mate records) and bias the
/// partition toward whichever records happen to be cheap.
pub(crate) const INGEST_SAMPLE_INTERVAL: u64 = 1021;

/// Where the ingest thread's serial CPU goes, per record.
///
/// The floor line says this thread *is* the limit -- on a 16-thread whole-genome
/// sort it is 137.2s of a 145.7s read span, against a worker-capacity floor of
/// 22.4s -- so the only question left is what the 137.2s is made of. Nothing else
/// in the phase can answer it: worker counters describe the pool, and wall-clock
/// spans describe the phase, and neither looks inside the loop.
///
/// Sampled rather than timed on every record, for the same reason
/// [`crate::merge_headroom::ConsumerSample`] is: the loop runs at ~175 ns/record
/// and an `Instant::now()` pair costs 15-35 ns on aarch64, so timing five
/// segments on every record would cost more than several of the segments it
/// measures. One record in [`INGEST_SAMPLE_INTERVAL`] is timed and scaled, and
/// the scale is reported next to the result.
///
/// Each field is **exactly one** timed region in the loop, which is what makes
/// [`Self::corrected`] valid: it subtracts one clock pair per field per sample,
/// so a field spanning two bracketed regions would be under-corrected by a
/// whole pair. That is why progress counting and spill probing are separate
/// fields rather than one "bookkeeping" bucket -- they sit at opposite ends of
/// the loop body and cannot share a bracket.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct IngestSample {
    /// Pulling the next record's bytes from the pool's decompressed stream.
    /// **Includes park time**, so it is not pure CPU -- `park_secs` measures that
    /// part exactly and separately.
    pub(crate) fetch: f64,
    /// Extracting the sort key from the record bytes.
    pub(crate) key: f64,
    /// Verifying that the lanes the chosen key variant drops are constant.
    pub(crate) verify: f64,
    /// Copying the record into the arena and appending its ref.
    pub(crate) push: f64,
    /// Counting the record toward the progress log.
    pub(crate) tick: f64,
    /// The spill probe's sample check and the memory-limit test that follows it.
    pub(crate) probe: f64,
}

impl IngestSample {
    /// Every segment multiplied by the sampling scale.
    #[must_use]
    pub(crate) fn scaled(self, scale: f64) -> Self {
        Self {
            fetch: self.fetch * scale,
            key: self.key * scale,
            verify: self.verify * scale,
            push: self.push * scale,
            tick: self.tick * scale,
            probe: self.probe * scale,
        }
    }

    /// Every segment with its own measurement overhead removed.
    ///
    /// One `Instant::now()`/`elapsed()` pair per segment per sampled record, and
    /// that pair's cost lands inside the interval it times. Clamped at zero: a
    /// segment cheaper than the clock measuring it cannot be resolved this way,
    /// and zero says so where a negative would read as a bug.
    #[must_use]
    pub(crate) fn corrected(self, samples: u64, overhead_nanos: u64) -> Self {
        if samples == 0 || overhead_nanos == 0 {
            return self;
        }
        #[expect(clippy::cast_precision_loss, reason = "sample counts stay below 2^52")]
        let per_segment = (samples * overhead_nanos) as f64 / 1e9;
        let fix = |v: f64| (v - per_segment).max(0.0);
        Self {
            fetch: fix(self.fetch),
            key: fix(self.key),
            verify: fix(self.verify),
            push: fix(self.push),
            tick: fix(self.tick),
            probe: fix(self.probe),
        }
    }

    /// The six segments summed.
    #[must_use]
    pub(crate) fn total(self) -> f64 {
        self.fetch + self.key + self.verify + self.push + self.tick + self.probe
    }
}

/// A scaled ingest sample checked against the read span it should partition.
#[derive(Debug, Clone, Copy)]
pub(crate) struct IngestPartition {
    /// Scaled, corrected per-segment seconds.
    pub(crate) segments: IngestSample,
    /// Measured read span, exact.
    pub(crate) read_secs: f64,
    /// Park time inside `segments.fetch`, measured exactly and separately.
    pub(crate) park_secs: f64,
}

impl IngestPartition {
    /// Read-span time the segments do not account for.
    ///
    /// **Signed on purpose.** A negative residual means the sample
    /// over-attributes -- clock overhead left inside the timed regions, or a
    /// sampling bias -- and the merge's first partition did exactly that,
    /// summing to 321.5s of a 189.3s loop. Only the sign made it visible; a
    /// clamped residual would have reported a tidy zero and the partition would
    /// have been believed.
    #[must_use]
    pub(crate) fn residual_secs(self) -> f64 {
        self.read_secs - self.segments.total()
    }

    /// Residual as a share of the read span, for judging whether the partition
    /// is trustworthy at all.
    #[must_use]
    pub(crate) fn residual_share(self) -> f64 {
        if self.read_secs > 0.0 { self.residual_secs() / self.read_secs } else { 0.0 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_the_residual_is_signed_so_over_attribution_is_visible() {
        // The merge's first partition summed to 321.5s of a 189.3s loop. A
        // clamped residual would have shown 0.0 and the numbers would have been
        // believed; the sign is what exposed the clock overhead inside them.
        let over = IngestPartition {
            segments: IngestSample { fetch: 200.0, ..IngestSample::default() },
            read_secs: 100.0,
            park_secs: 0.0,
        };
        assert!(over.residual_secs() < 0.0, "got {}", over.residual_secs());
        assert!((over.residual_share() + 1.0).abs() < 1e-9, "got {}", over.residual_share());
    }

    #[test]
    fn test_clock_correction_subtracts_one_pair_per_segment_per_sample() {
        // Ten samples, 20 ns per pair, five segments: each segment carries
        // 10 x 20 ns = 200 ns of clock, and each is corrected independently.
        let raw = IngestSample {
            fetch: 1e-6,
            key: 1e-6,
            verify: 1e-7,
            push: 1e-6,
            tick: 1e-6,
            probe: 1e-6,
        };
        let fixed = raw.corrected(10, 20);
        assert!((fixed.fetch - 0.8e-6).abs() < 1e-12, "got {}", fixed.fetch);
        // A segment cheaper than the clock that measured it clamps to zero rather
        // than going negative, which would read as a bug rather than as
        // "unresolvable by this method".
        assert!((fixed.verify - 0.0).abs() < 1e-12, "got {}", fixed.verify);
    }

    #[test]
    fn test_scaling_happens_before_correction_is_meaningful() {
        // Scale multiplies the sampled segments up to the whole loop; correction
        // works on the sampled scale. Applying them in the wrong order would
        // subtract one pair's cost from the *scaled* total rather than from each
        // sample, understating the correction by the scale factor.
        let raw = IngestSample { key: 2e-6, ..IngestSample::default() };
        let corrected_then_scaled = raw.corrected(10, 20).scaled(1000.0);
        let scaled_then_corrected = raw.scaled(1000.0).corrected(10, 20);
        assert!(corrected_then_scaled.key < scaled_then_corrected.key);
    }

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
