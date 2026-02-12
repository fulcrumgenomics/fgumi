//! Epsilon-Greedy scheduler.
//!
//! With probability (1-ε), exploits the best-known step (highest success rate).
//! With probability ε, explores by using a random order.
//! Simple but effective baseline for bandit algorithms.

use std::sync::atomic::{AtomicU64, Ordering};

use rand::RngExt;
use rand::SeedableRng;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;

use super::{BackpressureState, Scheduler};
use crate::unified_pipeline::base::PipelineStep;

/// Atomic counter for unique seeds across threads.
static SEED_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Epsilon-Greedy scheduler.
pub struct EpsilonGreedyScheduler {
    /// Thread ID.
    thread_id: usize,
    /// Total number of threads.
    num_threads: usize,
    /// Exploration probability (0.0-1.0).
    epsilon: f64,
    /// Attempts per step.
    attempts: [u64; 9],
    /// Successes per step.
    successes: [u64; 9],
    /// Random number generator.
    rng: SmallRng,
    /// Priority buffer.
    priority_buffer: [PipelineStep; 9],
}

impl EpsilonGreedyScheduler {
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

    /// Default epsilon (10% exploration).
    const DEFAULT_EPSILON: f64 = 0.1;

    /// Create a new Epsilon-Greedy scheduler.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize) -> Self {
        let seed = SEED_COUNTER
            .fetch_add(1, Ordering::Relaxed)
            .wrapping_add(thread_id as u64)
            .wrapping_mul(0x9E37_79B9_7F4A_7C15);
        Self {
            thread_id,
            num_threads,
            epsilon: Self::DEFAULT_EPSILON,
            attempts: [0; 9],
            successes: [0; 9],
            rng: SmallRng::seed_from_u64(seed),
            priority_buffer: Self::STEPS,
        }
    }

    /// Get the index of a step.
    fn step_index(step: PipelineStep) -> usize {
        Self::STEPS.iter().position(|&s| s == step).unwrap_or(0)
    }

    /// Get success rate for a step.
    #[allow(clippy::cast_precision_loss)]
    fn success_rate(&self, step_idx: usize) -> f64 {
        if self.attempts[step_idx] == 0 {
            0.5 // Neutral prior for unexplored steps
        } else {
            self.successes[step_idx] as f64 / self.attempts[step_idx] as f64
        }
    }
}

impl Scheduler for EpsilonGreedyScheduler {
    fn get_priorities(&mut self, _backpressure: BackpressureState) -> &[PipelineStep] {
        if self.rng.random::<f64>() < self.epsilon {
            // Explore: random order
            self.priority_buffer = Self::STEPS;
            self.priority_buffer.shuffle(&mut self.rng);
        } else {
            // Exploit: sort by success rate (descending)
            let mut rates: [(f64, usize); 9] = [(0.0, 0); 9];
            for (i, rate) in rates.iter_mut().enumerate() {
                *rate = (self.success_rate(i), i);
            }
            rates.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

            for (priority, (_, step_idx)) in rates.iter().enumerate() {
                self.priority_buffer[priority] = Self::STEPS[*step_idx];
            }
        }

        &self.priority_buffer
    }

    fn record_outcome(&mut self, step: PipelineStep, success: bool, _was_contention: bool) {
        let idx = Self::step_index(step);
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
    fn test_default_epsilon() {
        let scheduler = EpsilonGreedyScheduler::new(0, 8);
        assert!((scheduler.epsilon - 0.1).abs() < 0.001);
    }

    #[test]
    fn test_success_rate_calculation() {
        let mut scheduler = EpsilonGreedyScheduler::new(0, 8);
        scheduler.attempts[0] = 10;
        scheduler.successes[0] = 7;

        let rate = scheduler.success_rate(0);
        assert!((rate - 0.7).abs() < 0.001);
    }

    #[test]
    fn test_unexplored_neutral_prior() {
        let scheduler = EpsilonGreedyScheduler::new(0, 8);
        let rate = scheduler.success_rate(0);
        assert!((rate - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_priorities_returns_all_steps() {
        let mut scheduler = EpsilonGreedyScheduler::new(0, 8);
        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);
        assert_eq!(priorities.len(), 9);
    }
}
