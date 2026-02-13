//! Balanced Chase Drain scheduler.
//!
//! A variant of balanced-chase that adds "drain mode" - when Serialize
//! fails due to output queue full, temporarily prioritize Compress to
//! drain the queue before retrying Serialize.
//!
//! The key insight is that when Q6 (serialized queue) is full, threads
//! should focus on Compress to drain it rather than immediately retrying
//! Serialize (which will fail again).
//!
//! Drain mode is backpressure-driven: once entered, it stays active until
//! the output queue drops below the high-water mark (`output_high` = false).

use super::{BackpressureState, Direction, Scheduler};
use crate::unified_pipeline::base::{ActiveSteps, PipelineStep};

/// Balanced chase scheduler with drain mode for output backpressure.
pub struct BalancedChaseDrainScheduler {
    /// Thread ID.
    thread_id: usize,
    /// Total number of threads.
    #[allow(dead_code)]
    num_threads: usize,
    /// Current preferred step.
    current_step: PipelineStep,
    /// Current direction of movement.
    direction: Direction,
    /// Priority buffer.
    priority_buffer: [PipelineStep; 9],
    /// Which exclusive step this thread is responsible for (if any).
    exclusive_role: Option<PipelineStep>,
    /// Drain mode flag - when true, prioritize Compress over Serialize.
    /// Activated when Serialize fails, deactivated when `output_high` clears.
    drain_mode: bool,
    /// Memory drain mode - when true, prioritize Process over Group.
    /// Activated when Group fails due to memory pressure, deactivated when `memory_high` clears.
    memory_drain_mode: bool,
    /// Active steps in the pipeline.
    active_steps: ActiveSteps,
}

impl BalancedChaseDrainScheduler {
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

    /// Create a new balanced chase drain scheduler.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize, active_steps: ActiveSteps) -> Self {
        let (current_step, exclusive_role) = Self::determine_role(thread_id, num_threads);

        Self {
            thread_id,
            num_threads,
            current_step,
            direction: Direction::Forward,
            priority_buffer: Self::STEPS,
            exclusive_role,
            drain_mode: false,
            memory_drain_mode: false,
            active_steps,
        }
    }

    /// Determine thread role - same as balanced-chase.
    fn determine_role(
        thread_id: usize,
        num_threads: usize,
    ) -> (PipelineStep, Option<PipelineStep>) {
        use PipelineStep::{Compress, FindBoundaries, Group, Read, Serialize, Write};

        if thread_id == 0 {
            // T0 handles Read, but starts ready to help with Compress
            (Compress, Some(Read))
        } else if thread_id == num_threads - 1 && num_threads > 1 {
            // TN-1 handles Write, but starts ready to help with Compress
            (Compress, Some(Write))
        } else if thread_id == 1 && num_threads > 2 {
            // T1 handles FindBoundaries
            (FindBoundaries, Some(FindBoundaries))
        } else if thread_id == num_threads - 2 && num_threads > 3 {
            // TN-2 handles Group
            (Group, Some(Group))
        } else {
            // Middle threads start spread across bottleneck steps
            let step = if thread_id.is_multiple_of(2) { Compress } else { Serialize };
            (step, None)
        }
    }

    /// Get the index of a step.
    fn step_index(step: PipelineStep) -> usize {
        Self::STEPS.iter().position(|&s| s == step).unwrap_or(0)
    }

    /// Build priority list with drain mode support.
    fn build_priorities(&mut self, bp: BackpressureState) {
        use PipelineStep::{
            Compress, Decode, Decompress, FindBoundaries, Group, Process, Read, Serialize, Write,
        };

        let mut priorities = Vec::with_capacity(9);

        // 1. If we have an exclusive role, always try it first
        //    (but we won't be sticky on it after success)
        if let Some(role) = self.exclusive_role {
            priorities.push(role);
        }

        // 2. Memory drain mode - prioritize Process to drain q4_groups
        if (self.memory_drain_mode || bp.memory_high) && !priorities.contains(&Process) {
            priorities.push(Process);
        }

        // 3. Bottleneck steps - with drain mode logic
        //    In drain mode OR when output is high, prioritize Compress over Serialize
        if self.drain_mode || bp.output_high {
            // DRAIN MODE: Compress first to drain Q6
            if !priorities.contains(&Compress) {
                priorities.push(Compress);
            }
            // Serialize goes last (will be added at end)

            // 4. Current step (if not already added)
            if !priorities.contains(&self.current_step) {
                priorities.push(self.current_step);
            }

            // 5. Other parallel steps based on direction
            let parallel_order: &[PipelineStep] = match self.direction {
                Direction::Forward => &[Process, Decode, Decompress],
                Direction::Backward => &[Decompress, Decode, Process],
            };

            for &step in parallel_order {
                if !priorities.contains(&step) {
                    priorities.push(step);
                }
            }

            // 6. Exclusive steps last (except our role which is first)
            let exclusive_steps = [Read, FindBoundaries, Group, Write];
            for &step in &exclusive_steps {
                if !priorities.contains(&step) {
                    priorities.push(step);
                }
            }

            // 7. Add Serialize at end if in drain mode (so it's still in the list, just deprioritized)
            if self.drain_mode && !priorities.contains(&Serialize) {
                priorities.push(Serialize);
            }
        } else {
            // Normal mode: Serialize then Compress
            if !priorities.contains(&Serialize) {
                priorities.push(Serialize);
            }
            if !priorities.contains(&Compress) {
                priorities.push(Compress);
            }

            // Current step (if not already added)
            if !priorities.contains(&self.current_step) {
                priorities.push(self.current_step);
            }

            // Other parallel steps based on direction
            let parallel_order: &[PipelineStep] = match self.direction {
                Direction::Forward => &[Process, Decode, Decompress],
                Direction::Backward => &[Decompress, Decode, Process],
            };

            for &step in parallel_order {
                if !priorities.contains(&step) {
                    priorities.push(step);
                }
            }

            // Exclusive steps last (except our role which is first)
            let exclusive_steps = [Read, FindBoundaries, Group, Write];
            for &step in &exclusive_steps {
                if !priorities.contains(&step) {
                    priorities.push(step);
                }
            }
        }

        // Copy to buffer
        for (i, &step) in priorities.iter().take(9).enumerate() {
            self.priority_buffer[i] = step;
        }
    }
}

impl Scheduler for BalancedChaseDrainScheduler {
    fn get_priorities(&mut self, bp: BackpressureState) -> &[PipelineStep] {
        // Update direction based on backpressure
        if bp.output_high || bp.memory_high {
            self.direction = Direction::Forward;
        } else if bp.input_low && !bp.read_done {
            self.direction = Direction::Backward;
        }

        // Exit drain mode when backpressure clears
        // (output_high = false means Q6 has space again)
        if self.drain_mode && !bp.output_high {
            self.drain_mode = false;
        }

        // Exit memory drain mode when memory has drained below 50% (hysteresis)
        // This prevents thrashing at the limit boundary
        if self.memory_drain_mode && bp.memory_drained {
            self.memory_drain_mode = false;
        }

        self.build_priorities(bp);
        let n = self.active_steps.filter_in_place(&mut self.priority_buffer);
        &self.priority_buffer[..n]
    }

    fn record_outcome(&mut self, step: PipelineStep, success: bool, _was_contention: bool) {
        if success {
            // Pivot logic from balanced-chase:
            // After completing exclusive work, pivot to bottleneck
            if self.exclusive_role == Some(step) {
                self.current_step = PipelineStep::Compress;
            } else {
                // Normal case: stay on successful step
                self.current_step = step;
            }
        } else {
            // DRAIN MODE TRIGGER: Serialize failure enters drain mode
            if step == PipelineStep::Serialize {
                self.drain_mode = true;
            }

            // MEMORY DRAIN MODE TRIGGER: Group failure enters memory drain mode
            // This prioritizes Process to drain q4_groups and free memory
            if step == PipelineStep::Group {
                self.memory_drain_mode = true;
            }

            // Movement logic from balanced-chase
            let idx = Self::step_index(self.current_step);

            // Bias movement toward Compress (index 7) and Serialize (index 6)
            self.current_step = match self.direction {
                Direction::Forward => {
                    if idx < 7 {
                        Self::STEPS[idx + 1]
                    } else {
                        PipelineStep::Compress
                    }
                }
                Direction::Backward => {
                    if idx > 1 && idx != 7 && idx != 6 {
                        Self::STEPS[idx - 1]
                    } else {
                        // Stay near bottleneck
                        PipelineStep::Serialize
                    }
                }
            };
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
    fn test_reader_starts_on_compress() {
        let scheduler = BalancedChaseDrainScheduler::new(0, 8, all());
        assert_eq!(scheduler.current_step, PipelineStep::Compress);
        assert_eq!(scheduler.exclusive_role, Some(PipelineStep::Read));
        assert!(!scheduler.drain_mode);
    }

    #[test]
    fn test_serialize_failure_triggers_drain_mode() {
        let mut scheduler = BalancedChaseDrainScheduler::new(3, 8, all());

        // Simulate Serialize failure
        scheduler.record_outcome(PipelineStep::Serialize, false, false);

        // Drain mode should be active
        assert!(scheduler.drain_mode);
    }

    #[test]
    fn test_drain_mode_exits_when_backpressure_clears() {
        let mut scheduler = BalancedChaseDrainScheduler::new(3, 8, all());

        // Enter drain mode
        scheduler.drain_mode = true;

        // Backpressure still high - should stay in drain mode
        let bp_high = BackpressureState {
            output_high: true,
            input_low: false,
            read_done: false,
            memory_high: false,
            memory_drained: true,
        };
        scheduler.get_priorities(bp_high);
        assert!(scheduler.drain_mode);

        // Backpressure cleared - should exit drain mode
        let bp_clear = BackpressureState::default(); // output_high = false
        scheduler.get_priorities(bp_clear);
        assert!(!scheduler.drain_mode);
    }

    #[test]
    fn test_drain_mode_prioritizes_compress() {
        let mut scheduler = BalancedChaseDrainScheduler::new(3, 8, all());

        // Enter drain mode
        scheduler.drain_mode = true;

        // Keep output_high so we don't exit drain mode
        let bp = BackpressureState {
            output_high: true,
            input_low: false,
            read_done: false,
            memory_high: false,
            memory_drained: true,
        };
        let priorities = scheduler.get_priorities(bp);

        // Compress should be first (before Serialize)
        let compress_pos = priorities.iter().position(|&s| s == PipelineStep::Compress);
        let serialize_pos = priorities.iter().position(|&s| s == PipelineStep::Serialize);

        assert!(compress_pos.unwrap() < serialize_pos.unwrap());
    }

    #[test]
    fn test_normal_mode_serializes_first() {
        let mut scheduler = BalancedChaseDrainScheduler::new(3, 8, all());

        // Not in drain mode
        assert!(!scheduler.drain_mode);

        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);

        // Serialize should be before Compress in normal mode
        let compress_pos = priorities.iter().position(|&s| s == PipelineStep::Compress);
        let serialize_pos = priorities.iter().position(|&s| s == PipelineStep::Serialize);

        assert!(serialize_pos.unwrap() < compress_pos.unwrap());
    }

    #[test]
    fn test_memory_high_prioritizes_process() {
        let mut scheduler = BalancedChaseDrainScheduler::new(3, 8, all());

        // Memory high - should prioritize Process to drain q4_groups
        let bp = BackpressureState {
            output_high: false,
            input_low: false,
            read_done: false,
            memory_high: true,
            memory_drained: false,
        };
        let priorities = scheduler.get_priorities(bp);

        // Process should be early in the priority list (after exclusive role if any)
        let process_pos = priorities.iter().position(|&s| s == PipelineStep::Process);
        assert!(process_pos.unwrap() <= 2, "Process should be prioritized when memory is high");
    }
}
