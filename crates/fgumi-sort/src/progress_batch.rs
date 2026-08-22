//! Batched progress ticks for per-record loops.

use fgumi_bam_io::progress::ProgressTracker;

/// Records counted locally before touching the tracker's atomic.
///
/// A merge emits hundreds of millions of records on ONE serial thread, and
/// `ProgressTracker::log_if_needed` does a relaxed `fetch_add` per call. On
/// aarch64 Rust defaults to outline-atomics, so that is not an inline `LDADD`
/// but a call into `__aarch64_ldadd8_relax`, which runtime-checks for LSE and
/// can fall back to an `ldxr`/`stxr` loop. Profiling the merge consumer on
/// `c7g.4xlarge` put that helper at **28% of the thread's cycles** -- the single
/// largest entry, ahead of `LoserTree::replay` at 11% and the record memcpy at
/// 8% -- on a thread whose 145 ns/record now sets the merge's wall clock.
///
/// 4096 against a 1,000,000-record logging interval means a milestone is noticed
/// at most 0.4% late, which no progress line cares about, and removes 99.98% of
/// the atomic traffic.
const TICK_BATCH: u64 = 4096;

/// Accumulates record counts and forwards them to a [`ProgressTracker`] in
/// batches.
///
/// [`Self::flush`] must be called before the tracker's final log, or the last
/// partial batch is never counted and the totals come up short.
pub(crate) struct BatchedProgress {
    pending: u64,
}

impl BatchedProgress {
    pub(crate) const fn new() -> Self {
        Self { pending: 0 }
    }

    /// Count one record, forwarding to `tracker` once a batch has accumulated.
    #[inline]
    pub(crate) fn tick(&mut self, tracker: &ProgressTracker) {
        self.pending += 1;
        if self.pending >= TICK_BATCH {
            tracker.log_if_needed(self.pending);
            self.pending = 0;
        }
    }

    /// Forward whatever is left. Idempotent.
    pub(crate) fn flush(&mut self, tracker: &ProgressTracker) {
        if self.pending > 0 {
            tracker.log_if_needed(self.pending);
            self.pending = 0;
        }
    }

    /// Records not yet forwarded. For tests.
    #[cfg(test)]
    pub(crate) const fn pending(&self) -> u64 {
        self.pending
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The whole point is that the tracker's total is unchanged: batching may
    /// notice a milestone late, but it must never lose a record.
    #[test]
    fn test_batched_ticks_preserve_the_total() {
        let tracker = ProgressTracker::new("test").with_interval(1_000_000);
        let mut batch = BatchedProgress::new();
        for _ in 0..10_000 {
            batch.tick(&tracker);
        }
        batch.flush(&tracker);
        assert_eq!(tracker.count(), 10_000, "every record must reach the tracker");
    }

    /// Without the flush the last partial batch is silently dropped, which is the
    /// one way this change can corrupt a total.
    #[test]
    fn test_unflushed_remainder_is_the_only_loss_and_flush_recovers_it() {
        let tracker = ProgressTracker::new("test").with_interval(1_000_000);
        let mut batch = BatchedProgress::new();
        for _ in 0..(TICK_BATCH + 5) {
            batch.tick(&tracker);
        }
        assert_eq!(tracker.count(), TICK_BATCH, "one full batch forwarded");
        assert_eq!(batch.pending(), 5, "remainder still held locally");

        batch.flush(&tracker);
        assert_eq!(tracker.count(), TICK_BATCH + 5);
        batch.flush(&tracker);
        assert_eq!(tracker.count(), TICK_BATCH + 5, "flush is idempotent");
    }

    #[test]
    fn test_atomic_is_touched_once_per_batch_not_once_per_record() {
        let tracker = ProgressTracker::new("test").with_interval(1_000_000);
        let mut batch = BatchedProgress::new();
        for _ in 0..(TICK_BATCH - 1) {
            batch.tick(&tracker);
        }
        assert_eq!(tracker.count(), 0, "nothing forwarded until a batch fills");
        batch.tick(&tracker);
        assert_eq!(tracker.count(), TICK_BATCH);
    }
}
