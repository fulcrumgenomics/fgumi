//! Chase-bottleneck scheduler.
//!
//! Threads follow work: move downstream when blocked on output,
//! move upstream when starved on input, stay sticky on success.
//!
//! Algorithm:
//! 1. Try current step N first
//! 2. If output blocked → prefer downstream steps (N+1, N+2, ...)
//! 3. If input empty → prefer upstream steps (N-1, N-2, ...)
//! 4. On success → stay sticky on N
//! 5. Always return all steps to try in priority order

use super::{BackpressureState, Direction, Scheduler};
use crate::unified_pipeline::base::{ActiveSteps, PipelineStep};

/// Chase-bottleneck scheduler with dynamic step selection.
pub struct ChaseBottleneckScheduler {
    /// Thread ID.
    thread_id: usize,
    /// Total number of threads.
    #[allow(dead_code)]
    num_threads: usize,
    /// Current step this thread is focused on.
    current_step: PipelineStep,
    /// Current direction of movement.
    direction: Direction,
    /// Priority buffer for returning all steps in adaptive order.
    priority_buffer: [PipelineStep; 9],
    /// Active steps in the pipeline.
    active_steps: ActiveSteps,
}

impl ChaseBottleneckScheduler {
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

    /// Create a new chase-bottleneck scheduler.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize, active_steps: ActiveSteps) -> Self {
        let current_step = Self::initial_step(thread_id, num_threads);
        Self {
            thread_id,
            num_threads,
            current_step,
            direction: Direction::Forward,
            priority_buffer: Self::STEPS, // Will be reordered in get_priorities
            active_steps,
        }
    }

    /// Determine initial step based on thread role.
    fn initial_step(thread_id: usize, num_threads: usize) -> PipelineStep {
        use PipelineStep::{Compress, Decode, Decompress, Process, Read, Serialize, Write};
        if thread_id == 0 {
            Read // Reader starts at beginning
        } else if thread_id == num_threads - 1 && num_threads > 1 {
            Write // Writer starts at end
        } else {
            // Middle threads spread across middle steps
            let middle_steps = [Decompress, Decode, Process, Serialize, Compress];
            let idx = (thread_id - 1) % middle_steps.len();
            middle_steps[idx]
        }
    }

    /// Get the index of a step in the pipeline.
    fn step_index(step: PipelineStep) -> usize {
        Self::STEPS.iter().position(|&s| s == step).unwrap_or(0)
    }

    /// Build priority list starting from current step, expanding outward.
    /// If direction is Forward, prefer downstream steps (toward Write).
    /// If direction is Backward, prefer upstream steps (toward Read).
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap, clippy::cast_sign_loss)]
    fn build_priorities(&mut self) {
        let current_idx = Self::step_index(self.current_step);
        let mut priorities = [PipelineStep::Read; 9];
        let mut idx = 0;

        // Start with current step
        priorities[idx] = self.current_step;
        idx += 1;

        // Expand outward from current step
        // Direction determines which way we prefer to expand first
        let (first_dir, second_dir) = match self.direction {
            Direction::Forward => (1i32, -1i32),  // downstream first
            Direction::Backward => (-1i32, 1i32), // upstream first
        };

        let mut forward_offset = 1;
        let mut backward_offset = 1;

        while idx < 9 {
            // Try in preferred direction first
            let first_idx = (current_idx as i32 + first_dir * forward_offset as i32) as usize;
            if first_dir > 0 && forward_offset <= 8 - current_idx {
                if first_idx < 9 {
                    priorities[idx] = Self::STEPS[first_idx];
                    idx += 1;
                    forward_offset += 1;
                }
            } else if first_dir < 0 && backward_offset <= current_idx {
                let back_idx = current_idx - backward_offset;
                if back_idx < 9 {
                    priorities[idx] = Self::STEPS[back_idx];
                    idx += 1;
                    backward_offset += 1;
                }
            }

            if idx >= 9 {
                break;
            }

            // Then try in secondary direction
            let second_idx = (current_idx as i32 + second_dir * backward_offset as i32) as usize;
            if second_dir > 0 && forward_offset <= 8 - current_idx {
                if second_idx < 9 {
                    priorities[idx] = Self::STEPS[second_idx];
                    idx += 1;
                    forward_offset += 1;
                }
            } else if second_dir < 0 && backward_offset <= current_idx {
                let back_idx = current_idx - backward_offset;
                if back_idx < 9 {
                    priorities[idx] = Self::STEPS[back_idx];
                    idx += 1;
                    backward_offset += 1;
                }
            }

            // Fill remaining in any order
            if idx < 9 {
                let remaining_forward = 8 - current_idx - (forward_offset - 1);
                let remaining_backward = current_idx - (backward_offset - 1);
                if remaining_forward > 0 {
                    let fwd_idx = current_idx + forward_offset;
                    if fwd_idx < 9 {
                        priorities[idx] = Self::STEPS[fwd_idx];
                        idx += 1;
                        forward_offset += 1;
                    }
                }
                if idx < 9 && remaining_backward > 0 && backward_offset <= current_idx {
                    let back_idx = current_idx - backward_offset;
                    priorities[idx] = Self::STEPS[back_idx];
                    idx += 1;
                    backward_offset += 1;
                }
            }
        }

        self.priority_buffer = priorities;
    }
}

impl Scheduler for ChaseBottleneckScheduler {
    fn get_priorities(&mut self, bp: BackpressureState) -> &[PipelineStep] {
        // Apply global backpressure hints to adjust direction
        if bp.output_high {
            // Output backing up - prefer downstream (toward Write)
            self.direction = Direction::Forward;
        } else if bp.input_low && !bp.read_done {
            // Input starving - prefer upstream (toward Read)
            self.direction = Direction::Backward;
        }

        // Build priority list with current step first, expanding in preferred direction
        self.build_priorities();
        let n = self.active_steps.filter_in_place(&mut self.priority_buffer);
        &self.priority_buffer[..n]
    }

    fn record_outcome(&mut self, step: PipelineStep, success: bool, _was_contention: bool) {
        if success {
            // Stay sticky on successful step
            self.current_step = step;
        } else {
            // On failure, move in current direction to find work
            let idx = Self::step_index(self.current_step);
            match self.direction {
                Direction::Forward => {
                    // Move downstream
                    if idx < 8 {
                        self.current_step = Self::STEPS[idx + 1];
                    } else {
                        self.current_step = Self::STEPS[0]; // Wrap
                    }
                }
                Direction::Backward => {
                    // Move upstream
                    if idx > 0 {
                        self.current_step = Self::STEPS[idx - 1];
                    } else {
                        self.current_step = Self::STEPS[8]; // Wrap
                    }
                }
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

    fn all() -> ActiveSteps {
        ActiveSteps::all()
    }

    #[test]
    fn test_initial_step_reader() {
        let scheduler = ChaseBottleneckScheduler::new(0, 8, all());
        assert_eq!(scheduler.current_step, PipelineStep::Read);
    }

    #[test]
    fn test_initial_step_writer() {
        let scheduler = ChaseBottleneckScheduler::new(7, 8, all());
        assert_eq!(scheduler.current_step, PipelineStep::Write);
    }

    #[test]
    fn test_initial_step_middle() {
        let scheduler = ChaseBottleneckScheduler::new(2, 8, all());
        // Thread 2 (index 1 in middle threads) should start at Decode
        assert_eq!(scheduler.current_step, PipelineStep::Decode);
    }

    #[test]
    fn test_get_priorities_returns_all_steps() {
        let mut scheduler = ChaseBottleneckScheduler::new(3, 8, all());
        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);
        assert_eq!(priorities.len(), 9);
    }

    #[test]
    fn test_sticky_on_success() {
        let mut scheduler = ChaseBottleneckScheduler::new(3, 8, all());
        scheduler.record_outcome(PipelineStep::Compress, true, false);
        assert_eq!(scheduler.current_step, PipelineStep::Compress);
    }

    #[test]
    fn test_move_on_failure_forward() {
        let mut scheduler = ChaseBottleneckScheduler::new(3, 8, all());
        scheduler.current_step = PipelineStep::Process;
        scheduler.direction = Direction::Forward;
        scheduler.record_outcome(PipelineStep::Process, false, false);
        assert_eq!(scheduler.current_step, PipelineStep::Serialize);
    }

    #[test]
    fn test_move_on_failure_backward() {
        let mut scheduler = ChaseBottleneckScheduler::new(3, 8, all());
        scheduler.current_step = PipelineStep::Process;
        scheduler.direction = Direction::Backward;
        scheduler.record_outcome(PipelineStep::Process, false, false);
        assert_eq!(scheduler.current_step, PipelineStep::Group);
    }

    #[test]
    fn test_backpressure_output_high_sets_forward() {
        let mut scheduler = ChaseBottleneckScheduler::new(3, 8, all());
        scheduler.direction = Direction::Backward;
        let bp = BackpressureState {
            output_high: true,
            input_low: false,
            read_done: false,
            memory_high: false,
            memory_drained: true,
        };
        scheduler.get_priorities(bp);
        assert_eq!(scheduler.direction, Direction::Forward);
    }

    #[test]
    fn test_backpressure_input_low_sets_backward() {
        let mut scheduler = ChaseBottleneckScheduler::new(3, 8, all());
        scheduler.direction = Direction::Forward;
        let bp = BackpressureState {
            output_high: false,
            input_low: true,
            read_done: false,
            memory_high: false,
            memory_drained: true,
        };
        scheduler.get_priorities(bp);
        assert_eq!(scheduler.direction, Direction::Backward);
    }
}
