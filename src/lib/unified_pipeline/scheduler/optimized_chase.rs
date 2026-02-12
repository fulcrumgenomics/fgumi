//! Optimized Chase-Bottleneck scheduler.
//!
//! An enhanced version of chase-bottleneck with optimizations based on profiling:
//!
//! 1. **Compress-biased prioritization**: Compress is consistently the bottleneck
//!    (~40% of work time), so it's always in the top priorities for parallel threads.
//!
//! 2. **Exclusive step avoidance**: Non-specialized threads deprioritize exclusive
//!    steps (Read, `FindBoundaries`, Group, Write) to reduce contention.
//!
//! 3. **Bottleneck stickiness**: Extra sticky on Compress/Serialize - don't wander
//!    away from these high-value steps even on temporary failures.
//!
//! 4. **Contention-aware backoff**: After failing due to contention on an exclusive
//!    step, temporarily avoid it to let the specialized thread work.

use super::{BackpressureState, Direction, Scheduler};
use crate::unified_pipeline::base::PipelineStep;

/// Optimized chase-bottleneck scheduler.
pub struct OptimizedChaseScheduler {
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
    /// Consecutive successes on bottleneck steps (for extra stickiness).
    bottleneck_streak: u8,
    /// Backoff counter for exclusive steps after contention.
    exclusive_backoff: u8,
    /// Whether this thread is specialized for an exclusive step.
    is_exclusive_specialist: bool,
    /// The exclusive step this thread specializes in (if any).
    exclusive_specialty: Option<PipelineStep>,
}

impl OptimizedChaseScheduler {
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

    /// Parallel steps that any thread can help with.
    #[allow(dead_code)]
    const PARALLEL_STEPS: [PipelineStep; 5] = [
        PipelineStep::Decompress,
        PipelineStep::Decode,
        PipelineStep::Process,
        PipelineStep::Serialize,
        PipelineStep::Compress,
    ];

    /// Exclusive steps that only one thread can execute at a time.
    const EXCLUSIVE_STEPS: [PipelineStep; 4] = [
        PipelineStep::Read,
        PipelineStep::FindBoundaries,
        PipelineStep::Group,
        PipelineStep::Write,
    ];

    /// Bottleneck steps that deserve extra stickiness.
    const BOTTLENECK_STEPS: [PipelineStep; 2] = [PipelineStep::Compress, PipelineStep::Serialize];

    /// How many failures before leaving a bottleneck step.
    const BOTTLENECK_STICKY_THRESHOLD: u8 = 3;

    /// How many cycles to avoid exclusive steps after contention.
    const EXCLUSIVE_BACKOFF_CYCLES: u8 = 5;

    /// Create a new optimized chase scheduler.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize) -> Self {
        let (current_step, is_exclusive_specialist, exclusive_specialty) =
            Self::determine_role(thread_id, num_threads);

        Self {
            thread_id,
            num_threads,
            current_step,
            direction: Direction::Forward,
            priority_buffer: Self::STEPS,
            bottleneck_streak: 0,
            exclusive_backoff: 0,
            is_exclusive_specialist,
            exclusive_specialty,
        }
    }

    /// Determine thread role and initial step.
    fn determine_role(
        thread_id: usize,
        num_threads: usize,
    ) -> (PipelineStep, bool, Option<PipelineStep>) {
        use PipelineStep::{Compress, FindBoundaries, Group, Read, Serialize, Write};

        if thread_id == 0 {
            // Thread 0: Reader specialist
            (Read, true, Some(Read))
        } else if thread_id == num_threads - 1 && num_threads > 1 {
            // Last thread: Writer specialist
            (Write, true, Some(Write))
        } else if thread_id == 1 && num_threads > 2 {
            // Thread 1: Boundary specialist
            (FindBoundaries, true, Some(FindBoundaries))
        } else if thread_id == num_threads - 2 && num_threads > 3 {
            // Second-to-last: Group specialist
            (Group, true, Some(Group))
        } else {
            // Middle threads: Start on bottleneck steps (Compress/Serialize)
            // Spread them out to maximize parallelism
            let bottleneck_idx = (thread_id - 2) % 2;
            let step = if bottleneck_idx == 0 { Compress } else { Serialize };
            (step, false, None)
        }
    }

    /// Get the index of a step in the pipeline.
    fn step_index(step: PipelineStep) -> usize {
        Self::STEPS.iter().position(|&s| s == step).unwrap_or(0)
    }

    /// Check if a step is exclusive.
    fn is_exclusive(step: PipelineStep) -> bool {
        Self::EXCLUSIVE_STEPS.contains(&step)
    }

    /// Check if a step is a bottleneck.
    fn is_bottleneck(step: PipelineStep) -> bool {
        Self::BOTTLENECK_STEPS.contains(&step)
    }

    /// Build optimized priority list.
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap, clippy::cast_sign_loss)]
    fn build_priorities(&mut self, _bp: BackpressureState) {
        let mut priorities = Vec::with_capacity(9);
        let current_idx = Self::step_index(self.current_step);

        // 1. Current step first (if not backing off from exclusive)
        let skip_current = Self::is_exclusive(self.current_step)
            && !self.is_exclusive_specialist
            && self.exclusive_backoff > 0;

        if !skip_current {
            priorities.push(self.current_step);
        }

        // 2. For exclusive specialists, always include their specialty second
        if let Some(specialty) = self.exclusive_specialty {
            if specialty != self.current_step {
                priorities.push(specialty);
            }
        }

        // 3. Bottleneck boost: Compress and Serialize always in top priorities
        //    for non-exclusive-specialist threads
        if !self.is_exclusive_specialist {
            for &step in &Self::BOTTLENECK_STEPS {
                if !priorities.contains(&step) {
                    priorities.push(step);
                }
            }
        }

        // 4. Build remaining priorities based on direction
        let (first_dir, second_dir): (i32, i32) = match self.direction {
            Direction::Forward => (1, -1),  // downstream first
            Direction::Backward => (-1, 1), // upstream first
        };

        // Expand from current position
        for distance in 1..=8 {
            // Primary direction
            let idx1 = current_idx as i32 + first_dir * distance;
            if (0..9).contains(&idx1) {
                let step = Self::STEPS[idx1 as usize];
                if !priorities.contains(&step) && self.should_include_step(step) {
                    priorities.push(step);
                }
            }

            // Secondary direction
            let idx2 = current_idx as i32 + second_dir * distance;
            if (0..9).contains(&idx2) {
                let step = Self::STEPS[idx2 as usize];
                if !priorities.contains(&step) && self.should_include_step(step) {
                    priorities.push(step);
                }
            }
        }

        // 5. Fill any remaining steps (shouldn't happen, but safety)
        for &step in &Self::STEPS {
            if !priorities.contains(&step) {
                priorities.push(step);
            }
        }

        // Copy to buffer
        for (i, &step) in priorities.iter().take(9).enumerate() {
            self.priority_buffer[i] = step;
        }
    }

    /// Check if this thread should include a step in its priorities.
    fn should_include_step(&self, step: PipelineStep) -> bool {
        if !Self::is_exclusive(step) {
            // Always include parallel steps
            return true;
        }

        // For exclusive steps:
        if self.is_exclusive_specialist {
            // Specialists include all exclusive steps (they might help others)
            true
        } else {
            // Non-specialists only include exclusive steps if not backing off
            self.exclusive_backoff == 0
        }
    }
}

impl Scheduler for OptimizedChaseScheduler {
    fn get_priorities(&mut self, bp: BackpressureState) -> &[PipelineStep] {
        // Decrement backoff counter
        if self.exclusive_backoff > 0 {
            self.exclusive_backoff -= 1;
        }

        // Apply global backpressure hints
        if bp.output_high {
            self.direction = Direction::Forward;
        } else if bp.input_low && !bp.read_done {
            self.direction = Direction::Backward;
        }

        self.build_priorities(bp);
        &self.priority_buffer
    }

    fn record_outcome(&mut self, step: PipelineStep, success: bool, was_contention: bool) {
        if success {
            // Success: stay sticky
            self.current_step = step;

            // Track bottleneck streak
            if Self::is_bottleneck(step) {
                self.bottleneck_streak = self.bottleneck_streak.saturating_add(1);
            } else {
                self.bottleneck_streak = 0;
            }

            // Reset backoff on any success
            self.exclusive_backoff = 0;
        } else {
            // Failure handling

            // Contention on exclusive step: back off if not specialist
            if was_contention && Self::is_exclusive(step) && !self.is_exclusive_specialist {
                self.exclusive_backoff = Self::EXCLUSIVE_BACKOFF_CYCLES;
            }

            // Extra stickiness for bottleneck steps
            if Self::is_bottleneck(self.current_step)
                && self.bottleneck_streak < Self::BOTTLENECK_STICKY_THRESHOLD
            {
                // Stay on bottleneck step despite failure
                self.bottleneck_streak = self.bottleneck_streak.saturating_add(1);
                return;
            }

            // Normal movement: advance in current direction
            let idx = Self::step_index(self.current_step);
            self.current_step = match self.direction {
                Direction::Forward => {
                    if idx < 8 {
                        Self::STEPS[idx + 1]
                    } else {
                        // Wrap to first parallel step, not Read
                        PipelineStep::Decompress
                    }
                }
                Direction::Backward => {
                    if idx > 1 {
                        Self::STEPS[idx - 1]
                    } else {
                        // Wrap to last parallel step, not Write
                        PipelineStep::Compress
                    }
                }
            };

            // Reset bottleneck streak when moving
            self.bottleneck_streak = 0;
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
    fn test_thread_0_is_reader_specialist() {
        let scheduler = OptimizedChaseScheduler::new(0, 8);
        assert_eq!(scheduler.current_step, PipelineStep::Read);
        assert!(scheduler.is_exclusive_specialist);
        assert_eq!(scheduler.exclusive_specialty, Some(PipelineStep::Read));
    }

    #[test]
    fn test_thread_7_is_writer_specialist() {
        let scheduler = OptimizedChaseScheduler::new(7, 8);
        assert_eq!(scheduler.current_step, PipelineStep::Write);
        assert!(scheduler.is_exclusive_specialist);
        assert_eq!(scheduler.exclusive_specialty, Some(PipelineStep::Write));
    }

    #[test]
    fn test_thread_1_is_boundary_specialist() {
        let scheduler = OptimizedChaseScheduler::new(1, 8);
        assert_eq!(scheduler.current_step, PipelineStep::FindBoundaries);
        assert!(scheduler.is_exclusive_specialist);
    }

    #[test]
    fn test_thread_6_is_group_specialist() {
        let scheduler = OptimizedChaseScheduler::new(6, 8);
        assert_eq!(scheduler.current_step, PipelineStep::Group);
        assert!(scheduler.is_exclusive_specialist);
    }

    #[test]
    fn test_middle_threads_start_on_bottleneck() {
        let scheduler = OptimizedChaseScheduler::new(2, 8);
        assert!(
            scheduler.current_step == PipelineStep::Compress
                || scheduler.current_step == PipelineStep::Serialize
        );
        assert!(!scheduler.is_exclusive_specialist);
    }

    #[test]
    fn test_bottleneck_steps_in_priorities() {
        let mut scheduler = OptimizedChaseScheduler::new(3, 8);
        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);

        // Compress and Serialize should be in top 4 for non-specialist
        let compress_pos = priorities.iter().position(|&s| s == PipelineStep::Compress);
        let serialize_pos = priorities.iter().position(|&s| s == PipelineStep::Serialize);

        assert!(compress_pos.unwrap() < 4);
        assert!(serialize_pos.unwrap() < 4);
    }

    #[test]
    fn test_exclusive_backoff_after_contention() {
        let mut scheduler = OptimizedChaseScheduler::new(3, 8); // Non-specialist
        scheduler.current_step = PipelineStep::Group;

        // Fail with contention on exclusive step
        scheduler.record_outcome(PipelineStep::Group, false, true);

        assert_eq!(scheduler.exclusive_backoff, OptimizedChaseScheduler::EXCLUSIVE_BACKOFF_CYCLES);
    }

    #[test]
    fn test_specialist_no_backoff() {
        let mut scheduler = OptimizedChaseScheduler::new(0, 8); // Reader specialist

        // Fail with contention on Read
        scheduler.record_outcome(PipelineStep::Read, false, true);

        // Specialist should not back off from their specialty
        assert_eq!(scheduler.exclusive_backoff, 0);
    }

    #[test]
    fn test_bottleneck_stickiness() {
        let mut scheduler = OptimizedChaseScheduler::new(3, 8);
        scheduler.current_step = PipelineStep::Compress;

        // First failure: stay sticky
        scheduler.record_outcome(PipelineStep::Compress, false, false);
        assert_eq!(scheduler.current_step, PipelineStep::Compress);

        // Second failure: still sticky
        scheduler.record_outcome(PipelineStep::Compress, false, false);
        assert_eq!(scheduler.current_step, PipelineStep::Compress);

        // Third failure: still sticky (at threshold)
        scheduler.record_outcome(PipelineStep::Compress, false, false);
        assert_eq!(scheduler.current_step, PipelineStep::Compress);

        // Fourth failure: now move
        scheduler.record_outcome(PipelineStep::Compress, false, false);
        assert_ne!(scheduler.current_step, PipelineStep::Compress);
    }

    #[test]
    fn test_wrap_to_parallel_steps() {
        let mut scheduler = OptimizedChaseScheduler::new(3, 8);
        scheduler.current_step = PipelineStep::Write;
        scheduler.direction = Direction::Forward;
        scheduler.bottleneck_streak = 10; // Past threshold

        scheduler.record_outcome(PipelineStep::Write, false, false);

        // Should wrap to Decompress, not Read
        assert_eq!(scheduler.current_step, PipelineStep::Decompress);
    }

    #[test]
    fn test_priorities_returns_all_steps() {
        let mut scheduler = OptimizedChaseScheduler::new(3, 8);
        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);
        assert_eq!(priorities.len(), 9);
    }
}
