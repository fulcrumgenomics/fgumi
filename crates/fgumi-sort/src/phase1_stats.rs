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

// ============================================================================
// The input reader
// ============================================================================

/// What the input reader's exclusive thread spends its time on.
///
/// `ReadInputBlocks` is a *second* serial resource in Phase 1, and its cost is
/// not explained by the disk. Worker 0 owns the step exclusively
/// (`exclusive_step_for` returns `None` for every other worker), so one thread
/// frames the whole file: on the production cell that is 124.5s across 5,114,408
/// blocks -- **24.3 us per 8.4 KB block** -- while the same 43.1 GB at that
/// volume's measured single-stream ceiling (605 MB/s, direct I/O, caches
/// dropped) is only ~71s. The step is not disk-bound, and nothing accounts for
/// the remaining ~55s.
///
/// It cannot be split from outside the reader, because the disk wait is *inside*
/// it: `read_raw_blocks` frames blocks out of a 2 MiB `BufReader`, so whichever
/// block happens to exhaust the buffer pays for the refill and every other block
/// pays nothing. Timing the step measures framing and I/O fused together, and no
/// per-block average can separate them. [`TimedReader`] sits one layer *below*
/// the buffer and times the underlying `read()` itself, which is what makes
/// [`ReaderReport::framing_secs`] -- time inside `read_raw_blocks` that was not
/// spent in a `read()` -- mean anything.
///
/// Timed exactly rather than sampled: ~320k batches and ~21k refills against a
/// ~30 ns clock is well under 30 ms of measurement on a 124.5s step, five orders
/// of magnitude below the quantity. [`IngestSample`] has to sample only because
/// its loop runs 780 million times; this one does not.
#[derive(Debug, Default)]
pub(crate) struct ReaderStats {
    /// `read()` calls made against the file, below the 2 MiB buffer.
    refills: AtomicU64,
    /// Nanoseconds inside those calls.
    refill_nanos: AtomicU64,
    /// Bytes they returned.
    refill_bytes: AtomicU64,
    /// Successful `ReadInputBlocks` batches.
    batches: AtomicU64,
    /// Blocks those batches framed.
    blocks: AtomicU64,
    /// Nanoseconds inside `read_raw_blocks`, refill time included.
    read_raw_nanos: AtomicU64,
    /// Nanoseconds pushing framed blocks onto the decompress queue.
    dispatch_nanos: AtomicU64,
}

impl ReaderStats {
    /// Record one `read()` against the underlying file.
    pub(crate) fn record_refill(&self, elapsed_nanos: u64, bytes: usize) {
        self.refills.fetch_add(1, Ordering::Relaxed);
        self.refill_nanos.fetch_add(elapsed_nanos, Ordering::Relaxed);
        self.refill_bytes.fetch_add(bytes as u64, Ordering::Relaxed);
    }

    /// Record one framed batch: how long `read_raw_blocks` took and what it returned.
    pub(crate) fn record_batch(&self, elapsed_nanos: u64, blocks: usize) {
        self.batches.fetch_add(1, Ordering::Relaxed);
        self.blocks.fetch_add(blocks as u64, Ordering::Relaxed);
        self.read_raw_nanos.fetch_add(elapsed_nanos, Ordering::Relaxed);
    }

    /// Record one dispatch of a framed batch onto the queue.
    pub(crate) fn record_dispatch(&self, elapsed_nanos: u64) {
        self.dispatch_nanos.fetch_add(elapsed_nanos, Ordering::Relaxed);
    }

    /// A consistent view of the counters, against the step total that should
    /// contain them.
    pub(crate) fn snapshot(&self, step_secs: f64) -> ReaderReport {
        ReaderReport {
            refills: self.refills.load(Ordering::Relaxed),
            refill_secs: secs(self.refill_nanos.load(Ordering::Relaxed)),
            refill_bytes: self.refill_bytes.load(Ordering::Relaxed),
            batches: self.batches.load(Ordering::Relaxed),
            blocks: self.blocks.load(Ordering::Relaxed),
            read_raw_secs: secs(self.read_raw_nanos.load(Ordering::Relaxed)),
            dispatch_secs: secs(self.dispatch_nanos.load(Ordering::Relaxed)),
            step_secs,
        }
    }
}

/// Where the input reader's serial time went, as reported.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct ReaderReport {
    pub(crate) refills: u64,
    pub(crate) refill_secs: f64,
    pub(crate) refill_bytes: u64,
    pub(crate) batches: u64,
    pub(crate) blocks: u64,
    /// Time inside `read_raw_blocks`, **including** the refills it triggered.
    pub(crate) read_raw_secs: f64,
    pub(crate) dispatch_secs: f64,
    /// Exact `ReadInputBlocks` busy total the parts should add up to.
    pub(crate) step_secs: f64,
}

impl ReaderReport {
    /// Time inside `read_raw_blocks` that was **not** spent in a `read()`:
    /// header parse, validation, per-block allocation, and the body copy.
    ///
    /// This is the number the whole struct exists to produce. It is signed, and
    /// a small negative is a known, benign case rather than a bug: the header
    /// parse refills the buffer once before Phase 1 starts, so a few
    /// milliseconds of refill can sit outside every timed batch. A *large*
    /// negative means the reader was rebuilt or the counters were shared across
    /// runs, and should not be explained away.
    pub(crate) fn framing_secs(&self) -> f64 {
        self.read_raw_secs - self.refill_secs
    }

    /// Step time the parts do not account for -- the `try_lock`, the serial
    /// reservation, and dropping the guard.
    ///
    /// **Signed on purpose**, for the reason [`IngestPartition::residual_secs`]
    /// gives: the merge's first partition summed to 321.5s of a 189.3s loop, and
    /// only the sign made that visible.
    pub(crate) fn residual_secs(&self) -> f64 {
        self.step_secs - self.read_raw_secs - self.dispatch_secs
    }

    /// Residual as a share of the step, for judging the partition at a glance.
    pub(crate) fn residual_share(&self) -> f64 {
        if self.step_secs > 0.0 { self.residual_secs() / self.step_secs } else { 0.0 }
    }

    /// Spread `secs` over the blocks framed, in microseconds -- the unit the
    /// 24.3 us/block question is asked in.
    #[allow(clippy::cast_precision_loss, reason = "block counts stay far below 2^52")]
    pub(crate) fn per_block_micros(&self, secs: f64) -> f64 {
        if self.blocks == 0 { 0.0 } else { secs * 1_000_000.0 / self.blocks as f64 }
    }

    /// Throughput the refills actually achieved, in MB/s.
    ///
    /// The discriminant against the volume's measured ceiling: at the ceiling the
    /// refills are disk-bound and the remaining time is framing; well under it,
    /// the reader is leaving bandwidth on the floor and the fix is upstream of
    /// framing entirely.
    #[allow(clippy::cast_precision_loss, reason = "byte totals stay far below 2^52")]
    pub(crate) fn refill_mb_per_sec(&self) -> f64 {
        if self.refill_secs > 0.0 { self.refill_bytes as f64 / self.refill_secs / 1e6 } else { 0.0 }
    }
}

/// Times every `read()` made against whatever it wraps.
///
/// Belongs *below* the input `BufReader`, so it sees buffer refills rather than
/// per-block reads -- see [`ReaderStats`] for why that placement is the whole
/// point.
pub(crate) struct TimedReader<R> {
    inner: R,
    stats: std::sync::Arc<ReaderStats>,
}

impl<R> TimedReader<R> {
    pub(crate) fn new(inner: R, stats: std::sync::Arc<ReaderStats>) -> Self {
        Self { inner, stats }
    }
}

impl<R: std::io::Read> std::io::Read for TimedReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let started = std::time::Instant::now();
        let result = self.inner.read(buf);
        let elapsed = u64::try_from(started.elapsed().as_nanos()).unwrap_or(u64::MAX);
        // A failed read still cost time; counting zero bytes for it keeps the
        // throughput figure honest rather than crediting the failure with bytes.
        self.stats.record_refill(elapsed, result.as_ref().copied().unwrap_or(0));
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A reader report built from the production cell's shape: 124.5s of step
    /// time over 5,114,408 blocks, of which ~71s is disk.
    fn production_shape() -> ReaderReport {
        ReaderReport {
            refills: 20_560,
            refill_secs: 71.0,
            refill_bytes: 43_131_188_552,
            batches: 319_651,
            blocks: 5_114_408,
            read_raw_secs: 121.4,
            dispatch_secs: 3.1,
            step_secs: 124.5,
        }
    }

    #[test]
    fn test_framing_is_read_raw_time_with_the_refills_taken_out() {
        // The whole point of putting a timer below the BufReader: the refill
        // happens inside `read_raw_blocks`, so the raw span is framing and I/O
        // fused. Reporting the raw span as framing would credit userspace with
        // every second the disk spent.
        let report = production_shape();
        assert!((report.framing_secs() - 50.4).abs() < 1e-9, "got {}", report.framing_secs());
    }

    #[test]
    fn test_a_reader_report_partitions_the_step_with_a_signed_residual() {
        let report = production_shape();
        assert!((report.residual_secs() - 0.0).abs() < 1e-9, "got {}", report.residual_secs());

        // Over-attribution must be visible. If the parts claim more than the
        // step, a clamped residual reports a tidy zero and the partition gets
        // believed -- exactly how the merge's first partition summed to 321.5s
        // of a 189.3s loop and survived a reading.
        let over = ReaderReport { read_raw_secs: 200.0, ..production_shape() };
        assert!(over.residual_secs() < 0.0, "got {}", over.residual_secs());
        assert!(over.residual_share() < 0.0, "got {}", over.residual_share());
    }

    #[test]
    fn test_per_block_micros_reproduces_the_number_under_investigation() {
        // 124.5s over 5,114,408 blocks is the 24.3 us/block the campaign is
        // trying to explain; the report has to say so in that unit.
        let report = production_shape();
        let per_block = report.per_block_micros(report.step_secs);
        assert!((per_block - 24.34).abs() < 0.01, "got {per_block}");
        assert!((report.per_block_micros(0.0) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_refill_throughput_is_measured_against_the_bytes_actually_returned() {
        // 43.1 GB in 71.0s is ~607 MB/s, which is what makes "the disk is at its
        // ceiling" a claim rather than an assumption.
        let report = production_shape();
        assert!(
            (report.refill_mb_per_sec() - 607.5).abs() < 1.0,
            "got {}",
            report.refill_mb_per_sec()
        );
        assert!((ReaderReport::default().refill_mb_per_sec() - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_a_failed_read_is_timed_but_credited_no_bytes() {
        use std::io::Read;
        struct Failing;
        impl Read for Failing {
            fn read(&mut self, _buf: &mut [u8]) -> std::io::Result<usize> {
                Err(std::io::Error::other("device fell off"))
            }
        }
        let stats = std::sync::Arc::new(ReaderStats::default());
        let mut reader = TimedReader::new(Failing, std::sync::Arc::clone(&stats));
        let mut buf = [0u8; 8];
        assert!(reader.read(&mut buf).is_err());

        let report = stats.snapshot(0.0);
        assert_eq!(report.refills, 1, "the attempt still cost time and must be counted");
        assert_eq!(report.refill_bytes, 0, "a failed read must not inflate throughput");
    }

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
