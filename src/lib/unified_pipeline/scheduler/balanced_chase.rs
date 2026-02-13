//! Balanced Chase scheduler.
//!
//! A variant of chase-bottleneck optimized for even work distribution.
//! Key insight: exclusive step specialists (T0=Read, T7=Write) should
//! aggressively help with bottleneck steps (Compress/Serialize) when
//! their exclusive work is done, rather than staying sticky.
//!
//! Design:
//! - No thread is permanently "sticky" on any step
//! - After completing exclusive work, immediately pivot to bottleneck steps
//! - Compress and Serialize are always high priority for ALL threads
//! - Minimize idle time by keeping everyone busy on the real bottleneck

use super::{BackpressureState, Direction, Scheduler};
use crate::unified_pipeline::base::{ActiveSteps, PipelineStep};

/// Balanced chase scheduler focused on even work distribution.
pub struct BalancedChaseScheduler {
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
    /// Active steps in the pipeline.
    active_steps: ActiveSteps,
}

impl BalancedChaseScheduler {
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

    /// Create a new balanced chase scheduler.
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
            active_steps,
        }
    }

    /// Determine thread role - but DON'T make them sticky.
    fn determine_role(
        thread_id: usize,
        num_threads: usize,
    ) -> (PipelineStep, Option<PipelineStep>) {
        use PipelineStep::{Compress, FindBoundaries, Group, Read, Serialize, Write};

        if thread_id == 0 {
            // T0 handles Read, but starts ready to help with Compress
            (Compress, Some(Read))
        } else if thread_id == num_threads - 1 && num_threads > 1 {
            // T7 handles Write, but starts ready to help with Compress
            (Compress, Some(Write))
        } else if thread_id == 1 && num_threads > 2 {
            // T1 handles FindBoundaries
            (FindBoundaries, Some(FindBoundaries))
        } else if thread_id == num_threads - 2 && num_threads > 3 {
            // T6 handles Group
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

    /// Build priority list with bottleneck focus.
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

        // 2. Bottleneck steps are ALWAYS high priority for everyone
        //    Order based on backpressure
        if bp.output_high {
            // Output backing up: Compress then Serialize (push data out)
            if !priorities.contains(&Compress) {
                priorities.push(Compress);
            }
            if !priorities.contains(&Serialize) {
                priorities.push(Serialize);
            }
        } else {
            // Normal: Serialize then Compress (keep pipeline flowing)
            if !priorities.contains(&Serialize) {
                priorities.push(Serialize);
            }
            if !priorities.contains(&Compress) {
                priorities.push(Compress);
            }
        }

        // 3. Current step (if not already added)
        if !priorities.contains(&self.current_step) {
            priorities.push(self.current_step);
        }

        // 4. Other parallel steps based on direction
        let parallel_order: &[PipelineStep] = match self.direction {
            Direction::Forward => &[Process, Decode, Decompress],
            Direction::Backward => &[Decompress, Decode, Process],
        };

        for &step in parallel_order {
            if !priorities.contains(&step) {
                priorities.push(step);
            }
        }

        // 5. Exclusive steps last (except our role which is first)
        //    Only add if not already present
        let exclusive_steps = [Read, FindBoundaries, Group, Write];
        for &step in &exclusive_steps {
            if !priorities.contains(&step) {
                priorities.push(step);
            }
        }

        // Copy to buffer
        for (i, &step) in priorities.iter().take(9).enumerate() {
            self.priority_buffer[i] = step;
        }
    }
}

impl Scheduler for BalancedChaseScheduler {
    fn get_priorities(&mut self, bp: BackpressureState) -> &[PipelineStep] {
        // Update direction based on backpressure
        if bp.output_high {
            self.direction = Direction::Forward;
        } else if bp.input_low && !bp.read_done {
            self.direction = Direction::Backward;
        }

        self.build_priorities(bp);
        let n = self.active_steps.filter_in_place(&mut self.priority_buffer);
        &self.priority_buffer[..n]
    }

    fn record_outcome(&mut self, step: PipelineStep, success: bool, _was_contention: bool) {
        if success {
            // Key difference: Don't stay sticky on exclusive steps!
            // After completing exclusive work, pivot to bottleneck
            if self.exclusive_role == Some(step) {
                // Just completed our exclusive role - now help with bottleneck
                self.current_step = PipelineStep::Compress;
            } else {
                // Normal case: stay on successful step
                self.current_step = step;
            }
        } else {
            // On failure, move toward bottleneck steps
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
        let scheduler = BalancedChaseScheduler::new(0, 8, all());
        assert_eq!(scheduler.current_step, PipelineStep::Compress);
        assert_eq!(scheduler.exclusive_role, Some(PipelineStep::Read));
    }

    #[test]
    fn test_writer_starts_on_compress() {
        let scheduler = BalancedChaseScheduler::new(7, 8, all());
        assert_eq!(scheduler.current_step, PipelineStep::Compress);
        assert_eq!(scheduler.exclusive_role, Some(PipelineStep::Write));
    }

    #[test]
    fn test_exclusive_role_first_in_priorities() {
        let mut scheduler = BalancedChaseScheduler::new(0, 8, all());
        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);

        // Read should be first (exclusive role)
        assert_eq!(priorities[0], PipelineStep::Read);
        // Bottleneck steps should follow
        assert!(
            priorities[1] == PipelineStep::Serialize || priorities[1] == PipelineStep::Compress
        );
    }

    #[test]
    fn test_pivot_to_compress_after_exclusive() {
        let mut scheduler = BalancedChaseScheduler::new(0, 8, all());

        // Complete Read successfully
        scheduler.record_outcome(PipelineStep::Read, true, false);

        // Should pivot to Compress, not stay on Read
        assert_eq!(scheduler.current_step, PipelineStep::Compress);
    }

    #[test]
    fn test_middle_thread_no_exclusive_role() {
        let scheduler = BalancedChaseScheduler::new(3, 8, all());
        assert!(scheduler.exclusive_role.is_none());
    }

    #[test]
    fn test_bottleneck_always_in_top_priorities() {
        let mut scheduler = BalancedChaseScheduler::new(3, 8, all());
        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);

        let compress_pos = priorities.iter().position(|&s| s == PipelineStep::Compress);
        let serialize_pos = priorities.iter().position(|&s| s == PipelineStep::Serialize);

        // Both should be in top 3
        assert!(compress_pos.unwrap() < 3);
        assert!(serialize_pos.unwrap() < 3);
    }
}
