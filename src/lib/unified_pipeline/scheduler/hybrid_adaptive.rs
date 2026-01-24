//! Hybrid Adaptive scheduler.
//!
//! Starts with fixed-priority scheduling for efficiency. When experiencing
//! consecutive failures, switches to chase-bottleneck mode for adaptability.
//! Returns to fixed-priority when stable (consecutive successes).

use super::chase_bottleneck::ChaseBottleneckScheduler;
use super::fixed_priority::FixedPriorityScheduler;
use super::{BackpressureState, Scheduler};
use crate::unified_pipeline::base::PipelineStep;

/// Operating mode for the hybrid scheduler.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    /// Use fixed-priority scheduling.
    Fixed,
    /// Use chase-bottleneck scheduling.
    Chase,
}

/// Hybrid scheduler that switches between fixed-priority and chase-bottleneck.
pub struct HybridAdaptiveScheduler {
    /// Thread ID.
    thread_id: usize,
    /// Total number of threads.
    num_threads: usize,
    /// Fixed-priority scheduler instance.
    fixed: FixedPriorityScheduler,
    /// Chase-bottleneck scheduler instance.
    chase: ChaseBottleneckScheduler,
    /// Current operating mode.
    mode: Mode,
    /// Count of consecutive failures.
    consecutive_failures: usize,
    /// Count of consecutive successes.
    consecutive_successes: usize,
    /// Threshold to switch modes.
    switch_threshold: usize,
}

impl HybridAdaptiveScheduler {
    /// Default threshold for mode switching.
    const DEFAULT_THRESHOLD: usize = 5;

    /// Create a new hybrid adaptive scheduler.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize) -> Self {
        Self {
            thread_id,
            num_threads,
            fixed: FixedPriorityScheduler::new(thread_id, num_threads),
            chase: ChaseBottleneckScheduler::new(thread_id, num_threads),
            mode: Mode::Fixed,
            consecutive_failures: 0,
            consecutive_successes: 0,
            switch_threshold: Self::DEFAULT_THRESHOLD,
        }
    }
}

impl Scheduler for HybridAdaptiveScheduler {
    fn get_priorities(&mut self, backpressure: BackpressureState) -> &[PipelineStep] {
        match self.mode {
            Mode::Fixed => self.fixed.get_priorities(backpressure),
            Mode::Chase => self.chase.get_priorities(backpressure),
        }
    }

    fn record_outcome(&mut self, step: PipelineStep, success: bool, was_contention: bool) {
        // Forward outcome to current scheduler
        match self.mode {
            Mode::Fixed => self.fixed.record_outcome(step, success, was_contention),
            Mode::Chase => self.chase.record_outcome(step, success, was_contention),
        }

        // Track consecutive outcomes for mode switching
        if success {
            self.consecutive_successes += 1;
            self.consecutive_failures = 0;

            // Return to fixed mode when stable
            if self.mode == Mode::Chase && self.consecutive_successes >= self.switch_threshold {
                self.mode = Mode::Fixed;
                self.consecutive_successes = 0;
            }
        } else {
            self.consecutive_failures += 1;
            self.consecutive_successes = 0;

            // Switch to chase mode when struggling
            if self.mode == Mode::Fixed && self.consecutive_failures >= self.switch_threshold {
                self.mode = Mode::Chase;
                self.consecutive_failures = 0;
            }
        }
    }

    fn thread_id(&self) -> usize {
        self.thread_id
    }

    fn num_threads(&self) -> usize {
        self.num_threads
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_starts_in_fixed_mode() {
        let scheduler = HybridAdaptiveScheduler::new(0, 8);
        assert_eq!(scheduler.mode, Mode::Fixed);
    }

    #[test]
    fn test_switch_to_chase_after_failures() {
        let mut scheduler = HybridAdaptiveScheduler::new(0, 8);

        // Record failures to trigger switch
        for _ in 0..5 {
            scheduler.record_outcome(PipelineStep::Read, false, false);
        }

        assert_eq!(scheduler.mode, Mode::Chase);
    }

    #[test]
    fn test_return_to_fixed_after_successes() {
        let mut scheduler = HybridAdaptiveScheduler::new(0, 8);
        scheduler.mode = Mode::Chase;

        // Record successes to trigger return
        for _ in 0..5 {
            scheduler.record_outcome(PipelineStep::Read, true, false);
        }

        assert_eq!(scheduler.mode, Mode::Fixed);
    }

    #[test]
    fn test_consecutive_counter_reset() {
        let mut scheduler = HybridAdaptiveScheduler::new(0, 8);

        // Build up failures
        for _ in 0..3 {
            scheduler.record_outcome(PipelineStep::Read, false, false);
        }
        assert_eq!(scheduler.consecutive_failures, 3);

        // Success resets failure counter
        scheduler.record_outcome(PipelineStep::Read, true, false);
        assert_eq!(scheduler.consecutive_failures, 0);
        assert_eq!(scheduler.consecutive_successes, 1);
    }
}
