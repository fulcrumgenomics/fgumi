//! Which wall the merge is against, and how much of its wall clock is
//! recoverable without doing less work.
//!
//! A spill-heavy sort's merge has exactly three limits, and every optimisation
//! either moves one or moves nothing:
//!
//! 1. **The serial consumer.** One thread owns the loser tree and emits records
//!    in a strict global order, so its own CPU time is a hard floor. Nothing in
//!    the pool can shorten it.
//! 2. **Worker capacity.** Total worker busy divided by the active thread count.
//!    Read, decompress and output-compress are all parallel, so this is the floor
//!    they collectively impose.
//! 3. **Neither.** If the loop is longer than both floors the difference is
//!    coordination -- the consumer waiting for work that exists, or for work that
//!    has not been started.
//!
//! Reporting this matters because the three imply completely different fixes and
//! are routinely confused. Measured on one cell, 8 threads against 16: at 8 the
//! consumer was 98% of the loop and the recoverable share was 1.8%, so a
//! scheduling change could not have paid however well designed; at 16 the same
//! build left 28% of the loop recoverable. Same engine, same data, opposite
//! advice -- and no way to tell from wall clock alone.

/// Which limit the merge is actually against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Binding {
    /// The serial consumer's own CPU is the largest floor.
    Consumer,
    /// Worker capacity is the largest floor.
    Workers,
    /// Both floors are well below the loop: the gap is coordination.
    Coordination,
}

impl Binding {
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Consumer => "consumer serial CPU",
            Self::Workers => "worker capacity",
            Self::Coordination => "coordination (neither floor is close)",
        }
    }
}

/// The measured merge loop against the two floors it cannot go below.
#[derive(Debug, Clone, Copy)]
pub(crate) struct MergeFloors {
    /// Measured merge loop wall clock.
    pub(crate) loop_secs: f64,
    /// The consumer thread's own CPU: `loop - park - backpressure`, exact.
    pub(crate) consumer_secs: f64,
    /// Total worker busy seconds across every phase-2 worker.
    pub(crate) worker_busy_secs: f64,
    /// Active worker threads the busy total is spread over.
    pub(crate) threads: usize,
}

impl MergeFloors {
    /// Worker capacity floor: busy time spread perfectly over the active threads.
    ///
    /// Zero threads is reported as no floor rather than a division by zero: a
    /// merge with no workers is bounded by the consumer alone.
    #[must_use]
    pub(crate) fn worker_floor_secs(&self) -> f64 {
        if self.threads == 0 {
            return 0.0;
        }
        #[expect(clippy::cast_precision_loss, reason = "thread counts are tiny")]
        let threads = self.threads as f64;
        self.worker_busy_secs / threads
    }

    /// The larger of the two floors -- the shortest this merge could run without
    /// doing less work.
    #[must_use]
    pub(crate) fn floor_secs(&self) -> f64 {
        self.consumer_secs.max(self.worker_floor_secs())
    }

    /// Which floor is closest to the measured loop.
    ///
    /// `Coordination` when neither floor is within 10% of the loop: at that point
    /// the gap is larger than either limit and naming a "binding" floor would be
    /// misleading.
    #[must_use]
    pub(crate) fn binding(&self) -> Binding {
        let floor = self.floor_secs();
        if self.loop_secs > 0.0 && floor < 0.9 * self.loop_secs {
            return Binding::Coordination;
        }
        if self.consumer_secs >= self.worker_floor_secs() {
            Binding::Consumer
        } else {
            Binding::Workers
        }
    }

    /// Wall clock that could be recovered without reducing total work.
    ///
    /// Never negative: a floor above the measured loop means the estimate is off,
    /// not that time is owed.
    #[must_use]
    pub(crate) fn recoverable_secs(&self) -> f64 {
        (self.loop_secs - self.floor_secs()).max(0.0)
    }

    /// [`Self::recoverable_secs`] as a share of the loop, 0.0 when the loop is empty.
    #[must_use]
    pub(crate) fn recoverable_share(&self) -> f64 {
        if self.loop_secs <= 0.0 {
            return 0.0;
        }
        self.recoverable_secs() / self.loop_secs
    }
}

/// Sampled seconds for each step the merge consumer takes per record, in loop
/// order.
///
/// Sampled rather than timed on every record because the consumer runs at
/// 144-239 ns/record and one `Instant::now()` costs ~20-25 ns on aarch64: timing
/// five steps on every record would add more overhead than the steps it measures.
/// One in `MERGE_SAMPLE_INTERVAL` records is timed and the result scaled, so the
/// partition holds over the sample and the scale factor is reported alongside it.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct ConsumerSample {
    /// Reading the winner and publishing the predicted next source.
    pub(crate) publish: f64,
    /// Presenting the winning record: borrowed in place, or reassembled across a
    /// block boundary.
    pub(crate) present: f64,
    /// Handing the record to the output writer.
    pub(crate) write: f64,
    /// Advancing that source, which is where the consumer parks or decompresses a
    /// block itself. **Includes wait time**, so it is not pure CPU.
    pub(crate) advance: f64,
    /// Replaying the loser tree, or removing an exhausted source.
    pub(crate) tree: f64,
}

impl ConsumerSample {
    /// Every segment multiplied by the sampling scale.
    #[must_use]
    pub(crate) fn scaled(self, scale: f64) -> Self {
        Self {
            publish: self.publish * scale,
            present: self.present * scale,
            write: self.write * scale,
            advance: self.advance * scale,
            tree: self.tree * scale,
        }
    }

    /// The five segments summed.
    #[must_use]
    pub(crate) fn total(self) -> f64 {
        self.publish + self.present + self.write + self.advance + self.tree
    }
}

/// A scaled consumer sample checked against the loop it is supposed to partition.
///
/// The check is the point. Three of these segments were already collected and
/// scaled before this type existed, and were reported without ever being summed
/// against the measured loop -- so there was no way to know whether they
/// accounted for most of it or a third of it. A partition whose residual is not
/// reported is an estimate wearing a budget's clothes.
#[derive(Debug, Clone, Copy)]
pub(crate) struct LoopPartition {
    /// Scaled per-segment seconds.
    pub(crate) segments: ConsumerSample,
    /// Measured loop wall clock, exact.
    pub(crate) loop_secs: f64,
    /// Consumer park within `segments.advance`, measured exactly and separately.
    pub(crate) park_secs: f64,
}

impl LoopPartition {
    /// Loop time the segments do not account for.
    ///
    /// Signed on purpose: a negative residual means the sample over-attributes
    /// (clock overhead inside the timed regions, or a sampling bias), and hiding
    /// that behind a clamp would hide the one number that says the partition is
    /// unsound.
    #[must_use]
    pub(crate) fn unattributed_secs(self) -> f64 {
        self.loop_secs - self.segments.total()
    }

    /// `unattributed` as a share of the loop, for judging whether the partition is
    /// trustworthy at a glance.
    #[must_use]
    pub(crate) fn unattributed_share(self) -> f64 {
        if self.loop_secs <= 0.0 {
            return 0.0;
        }
        self.unattributed_secs() / self.loop_secs
    }

    /// The part of `advance` that was the consumer working rather than waiting.
    ///
    /// `advance` covers both parking on a block and decompressing one inline, and
    /// only the second is CPU the consumer could shed. Park is measured exactly, so
    /// subtracting it separates the two. Clamped at zero: park is exact while
    /// `advance` is sampled, so at small sample counts the estimate can land below
    /// it without meaning anything.
    #[must_use]
    pub(crate) fn advance_work_secs(self) -> f64 {
        (self.segments.advance - self.park_secs).max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::{Binding, ConsumerSample, LoopPartition, MergeFloors};
    use rstest::rstest;

    /// Assertions use an epsilon rather than exact float equality: every field
    /// here is a measured duration, so the arithmetic is inherently approximate.
    const EPS: f64 = 1e-9;

    /// The measured t8 cell after the targeted-depth changes: consumer 186.5s of a
    /// 190.0s loop, workers 1412.3s over 8 threads. The consumer binds and almost
    /// nothing is recoverable -- which is what says a scheduling fix cannot pay
    /// here however well designed.
    #[test]
    fn test_a_saturated_consumer_binds_and_leaves_almost_nothing() {
        let f = MergeFloors {
            loop_secs: 190.0,
            consumer_secs: 186.5,
            worker_busy_secs: 1412.3,
            threads: 8,
        };
        assert!((f.worker_floor_secs() - 176.537_5).abs() < 1e-4);
        assert_eq!(f.binding(), Binding::Consumer);
        assert!((f.floor_secs() - 186.5).abs() < EPS, "the larger floor is the consumer's");
        assert!((f.recoverable_secs() - 3.5).abs() < 1e-9);
        assert!(f.recoverable_share() < 0.02, "under 2%: got {}", f.recoverable_share());
    }

    /// The same build at 16 threads: consumer 112.1s of a 156.0s loop, workers
    /// 1420.2s over 16. Neither floor is near the loop, so the honest verdict is
    /// coordination and 28% is recoverable.
    #[test]
    fn test_a_loop_far_above_both_floors_reads_as_coordination() {
        let f = MergeFloors {
            loop_secs: 156.0,
            consumer_secs: 112.1,
            worker_busy_secs: 1420.2,
            threads: 16,
        };
        assert!((f.worker_floor_secs() - 88.762_5).abs() < 1e-4);
        assert_eq!(f.binding(), Binding::Coordination);
        assert!((f.recoverable_secs() - 43.9).abs() < 1e-9);
        assert!((f.recoverable_share() - 0.281_41).abs() < 1e-4, "got {}", f.recoverable_share());
    }

    /// Worker capacity can be the larger floor, and must be named when it is:
    /// the fix for it (do less work per record) is nothing like the fix for a
    /// slow consumer.
    #[test]
    fn test_worker_capacity_binds_when_it_is_the_larger_floor() {
        let f = MergeFloors {
            loop_secs: 200.0,
            consumer_secs: 60.0,
            worker_busy_secs: 1560.0, // 195s over 8
            threads: 8,
        };
        assert!((f.worker_floor_secs() - 195.0).abs() < EPS);
        assert_eq!(f.binding(), Binding::Workers);
        assert!((f.recoverable_secs() - 5.0).abs() < EPS);
    }

    /// A floor above the measured loop means the inputs disagree, which must clamp
    /// to zero rather than report negative recoverable time.
    #[test]
    fn test_a_floor_above_the_loop_never_reports_negative_headroom() {
        let f = MergeFloors {
            loop_secs: 100.0,
            consumer_secs: 130.0,
            worker_busy_secs: 0.0,
            threads: 8,
        };
        assert!((f.recoverable_secs() - 0.0).abs() < EPS);
        assert!((f.recoverable_share() - 0.0).abs() < EPS);
    }

    /// Degenerate inputs must not divide by zero or panic: a merge with no workers
    /// is bounded by its consumer, and an empty loop has no share to report.
    #[rstest]
    #[case::no_workers(100.0, 90.0, 500.0, 0)]
    #[case::empty_loop(0.0, 0.0, 0.0, 8)]
    fn test_degenerate_inputs_are_finite(
        #[case] loop_secs: f64,
        #[case] consumer_secs: f64,
        #[case] worker_busy_secs: f64,
        #[case] threads: usize,
    ) {
        let f = MergeFloors { loop_secs, consumer_secs, worker_busy_secs, threads };
        assert!(f.worker_floor_secs().is_finite());
        assert!(f.floor_secs().is_finite());
        assert!(f.recoverable_secs().is_finite());
        assert!((0.0..=1.0).contains(&f.recoverable_share()));
    }

    /// The residual is the whole reason this type exists, so it must be reported
    /// and not clamped. Segments summing to less than the loop means real time is
    /// unaccounted for.
    #[test]
    fn test_an_incomplete_partition_reports_the_time_it_cannot_account_for() {
        let p = LoopPartition {
            segments: ConsumerSample {
                publish: 5.0,
                present: 10.0,
                write: 20.0,
                advance: 50.0,
                tree: 15.0,
            },
            loop_secs: 156.0,
            park_secs: 43.3,
        };
        assert!((p.segments.total() - 100.0).abs() < EPS);
        assert!((p.unattributed_secs() - 56.0).abs() < EPS, "156 - 100");
        assert!((p.unattributed_share() - 0.358_97).abs() < 1e-4);
    }

    /// Over-attribution must surface as a negative residual rather than be hidden:
    /// it means the sampled timings include their own clock overhead, or the sample
    /// is biased toward expensive records. Either way the partition is unsound and
    /// the reader has to be told.
    #[test]
    fn test_over_attribution_stays_negative_instead_of_clamping() {
        let p = LoopPartition {
            segments: ConsumerSample {
                publish: 0.0,
                present: 0.0,
                write: 0.0,
                advance: 120.0,
                tree: 0.0,
            },
            loop_secs: 100.0,
            park_secs: 0.0,
        };
        assert!(p.unattributed_secs() < 0.0, "got {}", p.unattributed_secs());
        assert!((p.unattributed_secs() + 20.0).abs() < EPS);
    }

    /// `advance` mixes waiting with working and only the working half is sheddable,
    /// so park -- which is exact -- is subtracted out.
    #[test]
    fn test_advance_separates_the_consumer_working_from_the_consumer_waiting() {
        let p = LoopPartition {
            segments: ConsumerSample {
                publish: 0.0,
                present: 0.0,
                write: 0.0,
                advance: 60.0,
                tree: 0.0,
            },
            loop_secs: 156.0,
            park_secs: 43.3,
        };
        assert!((p.advance_work_secs() - 16.7).abs() < 1e-9);
    }

    /// Park is exact and `advance` is sampled, so the estimate can fall below it
    /// without that meaning negative work.
    #[test]
    fn test_advance_work_clamps_when_the_sample_lands_below_exact_park() {
        let p = LoopPartition {
            segments: ConsumerSample {
                publish: 0.0,
                present: 0.0,
                write: 0.0,
                advance: 30.0,
                tree: 0.0,
            },
            loop_secs: 156.0,
            park_secs: 43.3,
        };
        assert!((p.advance_work_secs() - 0.0).abs() < EPS);
    }

    /// Scaling multiplies every segment by the same factor and leaves the ratios
    /// between them alone.
    #[test]
    fn test_scaling_preserves_the_shape_of_the_sample() {
        let raw =
            ConsumerSample { publish: 1.0, present: 2.0, write: 3.0, advance: 4.0, tree: 5.0 };
        let scaled = raw.scaled(1021.0);
        assert!((scaled.total() - 15.0 * 1021.0).abs() < 1e-6);
        assert!((scaled.publish / scaled.tree - raw.publish / raw.tree).abs() < 1e-9);
    }
}
