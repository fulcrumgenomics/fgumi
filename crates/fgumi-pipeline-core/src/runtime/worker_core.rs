//! `WorkerCore`: per-thread state carried by the worker loop.

use std::time::Duration;

use crate::topology::StepIdx;

const BACKOFF_INITIAL_US: u64 = 1;
const BACKOFF_MAX_US: u64 = 1_000_000; // 1 second

pub struct WorkerCore {
    /// `0..n_workers`
    pub thread_id: usize,
    /// If this worker owns an `Exclusive` step, that step's index.
    /// Used by the runtime to skip `WorkerStepEntry::Skip` placeholders
    /// for Exclusive steps owned by other workers.
    pub exclusive_owner: Option<StepIdx>,
    /// If this worker is the sole eligible dispatcher for a `sticky` step
    /// (either an `Exclusive sticky` step it owns, or a `Serial + sticky`
    /// step whose `Affinity` targets this worker), the step's index.
    /// The driver drives this step in a tight inner loop until it returns
    /// `NoProgress` / `Contention` / `Finished`, then yields to round-
    /// robin. Mirrors legacy `pipeline/base.rs:4360-4392` sticky
    /// read.
    pub sticky_owner: Option<StepIdx>,
    /// Backoff duration in microseconds. Doubled on no-progress; reset on progress.
    backoff_us: u64,
}

impl WorkerCore {
    #[must_use]
    pub fn new(
        thread_id: usize,
        exclusive_owner: Option<StepIdx>,
        sticky_owner: Option<StepIdx>,
    ) -> Self {
        Self { thread_id, exclusive_owner, sticky_owner, backoff_us: BACKOFF_INITIAL_US }
    }

    pub fn reset_backoff(&mut self) {
        self.backoff_us = BACKOFF_INITIAL_US;
    }

    pub fn sleep_backoff(&self) {
        std::thread::sleep(Duration::from_micros(self.backoff_us));
    }

    pub fn increase_backoff(&mut self) {
        self.backoff_us = self.backoff_us.saturating_mul(2).min(BACKOFF_MAX_US);
    }

    #[cfg(test)]
    fn current_backoff_us(&self) -> u64 {
        self.backoff_us
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fresh_backoff_is_initial() {
        let w = WorkerCore::new(0, None, None);
        assert_eq!(w.current_backoff_us(), BACKOFF_INITIAL_US);
    }

    #[test]
    fn backoff_doubles_then_caps() {
        let mut w = WorkerCore::new(0, None, None);
        for _ in 0..50 {
            w.increase_backoff();
        }
        assert_eq!(w.current_backoff_us(), BACKOFF_MAX_US);
    }

    #[test]
    fn reset_after_progress() {
        let mut w = WorkerCore::new(0, None, None);
        for _ in 0..5 {
            w.increase_backoff();
        }
        w.reset_backoff();
        assert_eq!(w.current_backoff_us(), BACKOFF_INITIAL_US);
    }

    #[test]
    fn worker_with_exclusive_role() {
        let w = WorkerCore::new(2, Some(StepIdx(7)), Some(StepIdx(7)));
        assert_eq!(w.thread_id, 2);
        assert_eq!(w.exclusive_owner, Some(StepIdx(7)));
        assert_eq!(w.sticky_owner, Some(StepIdx(7)));
    }

    #[test]
    fn worker_with_sticky_serial_owner_only() {
        // E.g., the Serial+Affinity::Reader source — sticky_owner set,
        // exclusive_owner None.
        let w = WorkerCore::new(0, None, Some(StepIdx(0)));
        assert_eq!(w.exclusive_owner, None);
        assert_eq!(w.sticky_owner, Some(StepIdx(0)));
    }
}
