//! Upper Confidence Bound (UCB) scheduler.
//!
//! Selects steps based on upper confidence bounds that balance exploitation
//! (high success rate) with exploration (uncertainty). UCB score is:
//! `UCB_i` = `mean_i` + c * sqrt(ln(n) / `n_i`)
//!
//! Where `mean_i` is the success rate for step i, n is total attempts,
//! `n_i` is attempts for step i, and c is the exploration constant.

use super::{BackpressureState, Scheduler};
use crate::unified_pipeline::base::PipelineStep;

/// UCB scheduler with configurable exploration constant.
pub struct UCBScheduler {
    /// Thread ID.
    thread_id: usize,
    /// Total number of threads.
    num_threads: usize,
    /// Total attempts across all steps.
    total_attempts: u64,
    /// Attempts per step.
    attempts: [u64; 9],
    /// Successes per step.
    successes: [u64; 9],
    /// Exploration constant (typically sqrt(2) ≈ 1.414).
    exploration_c: f64,
    /// Priority buffer.
    priority_buffer: [PipelineStep; 9],
}

impl UCBScheduler {
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

    /// Default exploration constant (sqrt(2)).
    const DEFAULT_EXPLORATION_C: f64 = 1.414;

    /// Create a new UCB scheduler.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize) -> Self {
        Self {
            thread_id,
            num_threads,
            total_attempts: 0,
            attempts: [0; 9],
            successes: [0; 9],
            exploration_c: Self::DEFAULT_EXPLORATION_C,
            priority_buffer: Self::STEPS,
        }
    }

    /// Get the index of a step.
    fn step_index(step: PipelineStep) -> usize {
        Self::STEPS.iter().position(|&s| s == step).unwrap_or(0)
    }

    /// Calculate UCB score for a step.
    fn ucb_score(&self, step_idx: usize) -> f64 {
        let n_i = self.attempts[step_idx];

        // If never tried, return infinity to ensure exploration
        if n_i == 0 {
            return f64::INFINITY;
        }

        let mean = self.successes[step_idx] as f64 / n_i as f64;
        let ln_n = (self.total_attempts.max(1) as f64).ln();
        let exploration = self.exploration_c * (ln_n / n_i as f64).sqrt();

        mean + exploration
    }
}

impl Scheduler for UCBScheduler {
    fn get_priorities(&mut self, _backpressure: BackpressureState) -> &[PipelineStep] {
        // Calculate UCB score for each step
        let mut scores: [(f64, usize); 9] = [(0.0, 0); 9];
        for (i, score) in scores.iter_mut().enumerate() {
            *score = (self.ucb_score(i), i);
        }

        // Sort by UCB score (descending)
        scores.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

        // Build priority buffer
        for (priority, (_, step_idx)) in scores.iter().enumerate() {
            self.priority_buffer[priority] = Self::STEPS[*step_idx];
        }

        &self.priority_buffer
    }

    fn record_outcome(&mut self, step: PipelineStep, success: bool, _was_contention: bool) {
        let idx = Self::step_index(step);
        self.total_attempts += 1;
        self.attempts[idx] += 1;
        if success {
            self.successes[idx] += 1;
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
    fn test_unexplored_steps_prioritized() {
        let mut scheduler = UCBScheduler::new(0, 8);

        // Try one step many times
        for _ in 0..100 {
            scheduler.record_outcome(PipelineStep::Read, true, false);
        }

        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);

        // Unexplored steps should come first (infinite UCB)
        // Read should be later since it's well-explored
        let read_pos = priorities.iter().position(|&s| s == PipelineStep::Read).unwrap();
        assert!(read_pos > 0, "Well-explored Read should not be first");
    }

    #[test]
    fn test_high_success_rate_preferred() {
        let mut scheduler = UCBScheduler::new(0, 8);

        // Give Read high success rate
        for _ in 0..100 {
            scheduler.record_outcome(PipelineStep::Read, true, false);
        }
        // Give Write low success rate
        for _ in 0..100 {
            scheduler.record_outcome(PipelineStep::Write, false, false);
        }

        // Read should have higher UCB than Write (same exploration, better mean)
        let read_score = scheduler.ucb_score(0);
        let write_score = scheduler.ucb_score(8);
        assert!(read_score > write_score);
    }

    #[test]
    fn test_exploration_decreases_with_attempts() {
        let mut scheduler = UCBScheduler::new(0, 8);
        scheduler.total_attempts = 1000;

        // Step with few attempts has higher exploration bonus
        scheduler.attempts[0] = 10;
        scheduler.attempts[1] = 100;
        scheduler.successes[0] = 5; // 50% success
        scheduler.successes[1] = 50; // 50% success

        let score_few = scheduler.ucb_score(0);
        let score_many = scheduler.ucb_score(1);

        // Same success rate, but fewer attempts = higher exploration bonus
        assert!(score_few > score_many);
    }
}
