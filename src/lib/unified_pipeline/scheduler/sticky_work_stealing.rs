//! Sticky Work-Stealing scheduler.
//!
//! Each thread has a "home" step based on its role. When the home step
//! has no work, the thread "steals" work from adjacent steps first,
//! then from any step. After success elsewhere, gradually returns home.

use super::{BackpressureState, Scheduler};
use crate::unified_pipeline::base::PipelineStep;

/// Sticky work-stealing scheduler.
pub struct StickyWorkStealingScheduler {
    /// Thread ID.
    thread_id: usize,
    /// Total number of threads.
    num_threads: usize,
    /// Home step for this thread.
    home_step: PipelineStep,
    /// Current step (may be away from home).
    current_step: PipelineStep,
    /// Steps since last work at home.
    away_counter: usize,
    /// Threshold to try returning home.
    home_return_threshold: usize,
    /// Priority buffer.
    priority_buffer: [PipelineStep; 9],
}

impl StickyWorkStealingScheduler {
    /// All pipeline steps in order.
    const STEPS: [PipelineStep; 9] = [
        PipelineStep::Read,
        PipelineStep::Decompress,
        PipelineStep::FindBoundaries,
        PipelineStep::Decode,
        PipelineStep::Group,
        PipelineStep::Process,
        PipelineStep::Serialize,
        PipelineStep::Compress,
        PipelineStep::Write,
    ];

    /// Steps to try before returning home.
    const DEFAULT_HOME_RETURN_THRESHOLD: usize = 10;

    /// Create a new sticky work-stealing scheduler.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize) -> Self {
        let home_step = Self::determine_home_step(thread_id, num_threads);

        Self {
            thread_id,
            num_threads,
            home_step,
            current_step: home_step,
            away_counter: 0,
            home_return_threshold: Self::DEFAULT_HOME_RETURN_THRESHOLD,
            priority_buffer: Self::STEPS,
        }
    }

    /// Determine home step based on thread role.
    fn determine_home_step(thread_id: usize, num_threads: usize) -> PipelineStep {
        use PipelineStep::{
            Compress, Decode, Decompress, FindBoundaries, Group, Process, Read, Serialize, Write,
        };

        if thread_id == 0 {
            Read
        } else if thread_id == num_threads - 1 && num_threads > 1 {
            Write
        } else if thread_id == 1 && num_threads > 2 {
            FindBoundaries
        } else if thread_id == num_threads - 2 && num_threads > 3 {
            Group
        } else {
            // Middle threads: spread across parallel steps
            let parallel_steps = [Decompress, Decode, Process, Serialize, Compress];
            let idx = (thread_id.saturating_sub(2)) % parallel_steps.len();
            parallel_steps[idx]
        }
    }

    /// Get the index of a step.
    fn step_index(step: PipelineStep) -> usize {
        Self::STEPS.iter().position(|&s| s == step).unwrap_or(0)
    }

    /// Build priorities: home first, then adjacent, then others.
    fn build_priorities(&mut self) {
        let current_idx = Self::step_index(self.current_step);

        let mut idx = 0;

        // 1. Current step first (sticky)
        self.priority_buffer[idx] = self.current_step;
        idx += 1;

        // 2. Home step (if different)
        if self.current_step != self.home_step {
            self.priority_buffer[idx] = self.home_step;
            idx += 1;
        }

        // 3. Adjacent to current step
        if current_idx > 0 {
            let prev = Self::STEPS[current_idx - 1];
            if prev != self.home_step {
                self.priority_buffer[idx] = prev;
                idx += 1;
            }
        }
        if current_idx < 8 {
            let next = Self::STEPS[current_idx + 1];
            if next != self.home_step {
                self.priority_buffer[idx] = next;
                idx += 1;
            }
        }

        // 4. Remaining steps
        for step in Self::STEPS {
            if !self.priority_buffer[..idx].contains(&step) {
                self.priority_buffer[idx] = step;
                idx += 1;
            }
        }
    }
}

impl Scheduler for StickyWorkStealingScheduler {
    fn get_priorities(&mut self, _backpressure: BackpressureState) -> &[PipelineStep] {
        // Periodically try to return home
        if self.away_counter >= self.home_return_threshold && self.current_step != self.home_step {
            self.current_step = self.home_step;
            self.away_counter = 0;
        }

        self.build_priorities();
        &self.priority_buffer
    }

    fn record_outcome(&mut self, step: PipelineStep, success: bool, _was_contention: bool) {
        if success {
            self.current_step = step;
            if step == self.home_step {
                self.away_counter = 0;
            }
        } else {
            self.away_counter += 1;
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
    fn test_reader_home_step() {
        let scheduler = StickyWorkStealingScheduler::new(0, 8);
        assert_eq!(scheduler.home_step, PipelineStep::Read);
    }

    #[test]
    fn test_writer_home_step() {
        let scheduler = StickyWorkStealingScheduler::new(7, 8);
        assert_eq!(scheduler.home_step, PipelineStep::Write);
    }

    #[test]
    fn test_sticky_on_success() {
        let mut scheduler = StickyWorkStealingScheduler::new(3, 8);
        scheduler.record_outcome(PipelineStep::Compress, true, false);
        assert_eq!(scheduler.current_step, PipelineStep::Compress);
    }

    #[test]
    fn test_return_home_after_threshold() {
        let mut scheduler = StickyWorkStealingScheduler::new(0, 8);
        scheduler.current_step = PipelineStep::Write; // Away from home

        // Accumulate away counter
        for _ in 0..10 {
            scheduler.record_outcome(PipelineStep::Write, false, false);
        }

        let bp = BackpressureState::default();
        scheduler.get_priorities(bp);

        assert_eq!(scheduler.current_step, PipelineStep::Read); // Returned home
    }
}
