//! Fixed-priority scheduler (original implementation).
//!
//! Each thread has predetermined priorities based on thread role.
//! Includes backpressure override for output queue pressure.

use super::{BackpressureState, Scheduler};
use crate::unified_pipeline::base::{ActiveSteps, PipelineStep};

/// Fixed-priority scheduler with per-thread role-based priorities.
pub struct FixedPriorityScheduler {
    /// Thread ID (0 = reader, N-1 = writer).
    thread_id: usize,
    /// Total number of threads.
    #[allow(dead_code)]
    num_threads: usize,
    /// Fixed priority order for this thread.
    priorities: [PipelineStep; 9],
    /// Backpressure priority (drain output).
    backpressure_drain: [PipelineStep; 9],
    /// Backpressure priority (fill input).
    backpressure_fill: [PipelineStep; 9],
    /// Active steps in the pipeline.
    active_steps: ActiveSteps,
    /// Scratch buffer for filtering active priorities.
    active_priorities: [PipelineStep; 9],
    /// Number of active priorities in the scratch buffer.
    num_active: usize,
}

impl FixedPriorityScheduler {
    /// Create a new fixed-priority scheduler for the given thread.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize, active_steps: ActiveSteps) -> Self {
        use PipelineStep::{
            Compress, Decode, Decompress, FindBoundaries, Group, Process, Read, Serialize, Write,
        };

        // Backpressure overrides (same for all threads)
        let backpressure_drain =
            [Write, Compress, Serialize, Process, Group, Decode, FindBoundaries, Decompress, Read];
        let backpressure_fill =
            [Read, Decompress, FindBoundaries, Decode, Group, Process, Serialize, Compress, Write];

        // Thread-specific priorities (from original ThreadScheduler)
        let is_reader = thread_id == 0;
        let is_writer = num_threads > 1 && thread_id == num_threads - 1;
        let is_boundary_early = num_threads > 2 && thread_id == 1;
        let is_boundary_late = num_threads > 3 && thread_id == num_threads - 2;

        let priorities = if is_reader {
            // Thread 0: prioritize early pipeline stages
            [Read, Decompress, FindBoundaries, Decode, Group, Process, Serialize, Compress, Write]
        } else if is_writer {
            // Thread N-1: prioritize late pipeline stages
            [Write, Compress, Serialize, Process, Group, Decode, FindBoundaries, Decompress, Read]
        } else if is_boundary_early {
            // Thread 1: boundary specialist, focus on early exclusive steps
            [FindBoundaries, Decompress, Decode, Group, Process, Serialize, Compress, Write, Read]
        } else if is_boundary_late {
            // Thread N-2: boundary specialist, focus on Group
            [Group, Decode, FindBoundaries, Process, Decompress, Serialize, Compress, Write, Read]
        } else {
            // Middle threads: focus on parallel work, rotate based on thread_id
            let rotation = (thread_id - 1) % 5;
            let mut p = [
                Decompress,
                Decode,
                Process,
                Serialize,
                Compress,
                Group,
                FindBoundaries,
                Write,
                Read,
            ];
            p[0..5].rotate_left(rotation);
            p
        };

        Self {
            thread_id,
            num_threads,
            priorities,
            backpressure_drain,
            backpressure_fill,
            active_steps,
            active_priorities: [PipelineStep::Read; 9],
            num_active: 9,
        }
    }
}

impl Scheduler for FixedPriorityScheduler {
    fn get_priorities(&mut self, bp: BackpressureState) -> &[PipelineStep] {
        let source = if bp.output_high {
            &self.backpressure_drain
        } else if bp.input_low && !bp.read_done {
            &self.backpressure_fill
        } else {
            &self.priorities
        };
        self.active_priorities = *source;
        self.num_active = self.active_steps.filter_in_place(&mut self.active_priorities);
        &self.active_priorities[..self.num_active]
    }

    fn record_outcome(&mut self, _step: PipelineStep, _success: bool, _was_contention: bool) {
        // Fixed priority scheduler doesn't learn from outcomes
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
    fn test_reader_priorities() {
        let mut scheduler = FixedPriorityScheduler::new(0, 8, all());
        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);
        assert_eq!(priorities[0], PipelineStep::Read);
    }

    #[test]
    fn test_writer_priorities() {
        let mut scheduler = FixedPriorityScheduler::new(7, 8, all());
        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);
        assert_eq!(priorities[0], PipelineStep::Write);
    }

    #[test]
    fn test_boundary_specialist_early() {
        let mut scheduler = FixedPriorityScheduler::new(1, 8, all());
        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);
        assert_eq!(priorities[0], PipelineStep::FindBoundaries);
    }

    #[test]
    fn test_boundary_specialist_late() {
        let mut scheduler = FixedPriorityScheduler::new(6, 8, all());
        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);
        assert_eq!(priorities[0], PipelineStep::Group);
    }

    #[test]
    fn test_backpressure_drain() {
        let mut scheduler = FixedPriorityScheduler::new(2, 8, all());
        let bp = BackpressureState {
            output_high: true,
            input_low: false,
            read_done: false,
            memory_high: false,
            memory_drained: true,
        };
        let priorities = scheduler.get_priorities(bp);
        assert_eq!(priorities[0], PipelineStep::Write);
    }

    #[test]
    fn test_backpressure_fill() {
        let mut scheduler = FixedPriorityScheduler::new(2, 8, all());
        let bp = BackpressureState {
            output_high: false,
            input_low: true,
            read_done: false,
            memory_high: false,
            memory_drained: true,
        };
        let priorities = scheduler.get_priorities(bp);
        assert_eq!(priorities[0], PipelineStep::Read);
    }

    #[test]
    fn test_backpressure_fill_not_when_read_done() {
        let mut scheduler = FixedPriorityScheduler::new(2, 8, all());
        let bp = BackpressureState {
            output_high: false,
            input_low: true,
            read_done: true,
            memory_high: false,
            memory_drained: true,
        };
        let priorities = scheduler.get_priorities(bp);
        // Should use normal priorities when read is done
        assert_ne!(priorities[0], PipelineStep::Read);
    }

    #[test]
    fn test_two_threads() {
        let mut s0 = FixedPriorityScheduler::new(0, 2, all());
        let mut s1 = FixedPriorityScheduler::new(1, 2, all());
        let bp = BackpressureState::default();

        assert_eq!(s0.get_priorities(bp)[0], PipelineStep::Read);
        assert_eq!(s1.get_priorities(bp)[0], PipelineStep::Write);
    }
}
