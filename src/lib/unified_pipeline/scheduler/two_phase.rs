//! Two-Phase Pipeline scheduler.
//!
//! Optimizes for different pipeline phases:
//! - Startup: Chase-bottleneck to fill pipeline quickly
//! - Steady-state: Fixed-priority for efficiency
//! - Drain: Chase-bottleneck to empty pipeline
//!
//! Transitions based on `read_done` flag and success rates.

use super::chase_bottleneck::ChaseBottleneckScheduler;
use super::fixed_priority::FixedPriorityScheduler;
use super::{BackpressureState, Scheduler};
use crate::unified_pipeline::base::{ActiveSteps, PipelineStep};

/// Pipeline phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Phase {
    /// Initial phase - filling the pipeline.
    Startup,
    /// Main processing phase.
    SteadyState,
    /// Final phase - draining the pipeline.
    Drain,
}

/// Two-phase scheduler optimized for pipeline lifecycle.
pub struct TwoPhaseScheduler {
    /// Thread ID.
    thread_id: usize,
    /// Total number of threads.
    num_threads: usize,
    /// Fixed-priority scheduler.
    fixed: FixedPriorityScheduler,
    /// Chase-bottleneck scheduler.
    chase: ChaseBottleneckScheduler,
    /// Current phase.
    phase: Phase,
    /// Items processed counter.
    items_processed: u64,
    /// Threshold to exit startup phase.
    startup_threshold: u64,
    /// Active steps in the pipeline.
    active_steps: ActiveSteps,
}

impl TwoPhaseScheduler {
    /// Items to process before transitioning from startup.
    const DEFAULT_STARTUP_THRESHOLD: u64 = 1000;

    /// Create a new two-phase scheduler.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize, active_steps: ActiveSteps) -> Self {
        Self {
            thread_id,
            num_threads,
            fixed: FixedPriorityScheduler::new(thread_id, num_threads, active_steps.clone()),
            chase: ChaseBottleneckScheduler::new(thread_id, num_threads, active_steps.clone()),
            phase: Phase::Startup,
            items_processed: 0,
            startup_threshold: Self::DEFAULT_STARTUP_THRESHOLD,
            active_steps,
        }
    }
}

impl Scheduler for TwoPhaseScheduler {
    fn get_priorities(&mut self, backpressure: BackpressureState) -> &[PipelineStep] {
        // Update phase based on pipeline state
        match self.phase {
            Phase::Startup => {
                if self.items_processed >= self.startup_threshold {
                    self.phase = Phase::SteadyState;
                }
            }
            Phase::SteadyState => {
                if backpressure.read_done {
                    self.phase = Phase::Drain;
                }
            }
            Phase::Drain => {
                // Stay in drain until complete
            }
        }

        // Use appropriate scheduler
        match self.phase {
            Phase::Startup | Phase::Drain => self.chase.get_priorities(backpressure),
            Phase::SteadyState => self.fixed.get_priorities(backpressure),
        }
    }

    fn record_outcome(&mut self, step: PipelineStep, success: bool, was_contention: bool) {
        if success {
            self.items_processed += 1;
        }

        // Forward to current scheduler
        match self.phase {
            Phase::Startup | Phase::Drain => {
                self.chase.record_outcome(step, success, was_contention);
            }
            Phase::SteadyState => {
                self.fixed.record_outcome(step, success, was_contention);
            }
        }
    }

    fn thread_id(&self) -> usize {
        self.thread_id
    }

    fn num_threads(&self) -> usize {
        self.num_threads
    }

    fn active_steps(&self) -> &ActiveSteps {
        &self.active_steps
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_starts_in_startup_phase() {
        let scheduler = TwoPhaseScheduler::new(0, 8, ActiveSteps::all());
        assert_eq!(scheduler.phase, Phase::Startup);
    }

    #[test]
    fn test_transition_to_steady_state() {
        let mut scheduler = TwoPhaseScheduler::new(0, 8, ActiveSteps::all());
        let bp = BackpressureState::default();

        // Process enough items to exit startup
        for _ in 0..1000 {
            scheduler.record_outcome(PipelineStep::Process, true, false);
        }
        scheduler.get_priorities(bp);

        assert_eq!(scheduler.phase, Phase::SteadyState);
    }

    #[test]
    fn test_transition_to_drain() {
        let mut scheduler = TwoPhaseScheduler::new(0, 8, ActiveSteps::all());
        scheduler.phase = Phase::SteadyState;

        let bp = BackpressureState {
            output_high: false,
            input_low: false,
            read_done: true,
            memory_high: false,
            memory_drained: true,
        };
        scheduler.get_priorities(bp);

        assert_eq!(scheduler.phase, Phase::Drain);
    }
}
