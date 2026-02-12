//! Learned Affinity scheduler.
//!
//! Tracks success rates for each (thread, step) pair and builds a learned
//! priority order. Periodically explores to avoid local minima.
//! Self-tuning to workload characteristics over time.

use std::sync::atomic::{AtomicU64, Ordering};

use rand::Rng;
use rand::SeedableRng;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;

use super::{BackpressureState, Scheduler};
use crate::unified_pipeline::base::PipelineStep;

/// Atomic counter for unique seeds across threads.
static SEED_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Learned affinity scheduler.
pub struct LearnedAffinityScheduler {
    /// Thread ID.
    thread_id: usize,
    /// Total number of threads.
    num_threads: usize,
    /// Success counts per step.
    successes: [u64; 9],
    /// Attempt counts per step.
    attempts: [u64; 9],
    /// Exploration probability.
    exploration_rate: f64,
    /// Exploration rate decay per 1000 attempts.
    exploration_decay: f64,
    /// Minimum exploration rate.
    min_exploration_rate: f64,
    /// Total attempts (for decay calculation).
    total_attempts: u64,
    /// Random number generator.
    rng: SmallRng,
    /// Priority buffer.
    priority_buffer: [PipelineStep; 9],
}

impl LearnedAffinityScheduler {
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

    /// Initial exploration rate.
    const DEFAULT_EXPLORATION_RATE: f64 = 0.3;
    /// Decay factor per 1000 attempts.
    const DEFAULT_EXPLORATION_DECAY: f64 = 0.95;
    /// Minimum exploration rate.
    const DEFAULT_MIN_EXPLORATION: f64 = 0.05;

    /// Create a new learned affinity scheduler.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize) -> Self {
        let seed = SEED_COUNTER
            .fetch_add(1, Ordering::Relaxed)
            .wrapping_add(thread_id as u64)
            .wrapping_mul(0x9E37_79B9_7F4A_7C15);
        Self {
            thread_id,
            num_threads,
            successes: [0; 9],
            attempts: [0; 9],
            exploration_rate: Self::DEFAULT_EXPLORATION_RATE,
            exploration_decay: Self::DEFAULT_EXPLORATION_DECAY,
            min_exploration_rate: Self::DEFAULT_MIN_EXPLORATION,
            total_attempts: 0,
            rng: SmallRng::seed_from_u64(seed),
            priority_buffer: Self::STEPS,
        }
    }

    /// Get the index of a step.
    fn step_index(step: PipelineStep) -> usize {
        Self::STEPS.iter().position(|&s| s == step).unwrap_or(0)
    }

    /// Get affinity score for a step.
    #[expect(
        clippy::cast_precision_loss,
        reason = "affinity ratio doesn't need full u64 precision"
    )]
    fn affinity(&self, step_idx: usize) -> f64 {
        if self.attempts[step_idx] == 0 {
            0.5 // Neutral prior
        } else {
            self.successes[step_idx] as f64 / self.attempts[step_idx] as f64
        }
    }

    /// Get current exploration rate (decays over time).
    #[expect(
        clippy::cast_possible_truncation,
        reason = "decay_periods won't exceed i32 range in any practical run"
    )]
    fn current_exploration_rate(&self) -> f64 {
        let decay_periods = self.total_attempts / 1000;
        let decayed = self.exploration_rate * self.exploration_decay.powi(decay_periods as i32);
        decayed.max(self.min_exploration_rate)
    }
}

impl Scheduler for LearnedAffinityScheduler {
    fn get_priorities(&mut self, _backpressure: BackpressureState) -> &[PipelineStep] {
        let explore_rate = self.current_exploration_rate();

        if self.rng.random::<f64>() < explore_rate {
            // Explore: random order
            self.priority_buffer = Self::STEPS;
            self.priority_buffer.shuffle(&mut self.rng);
        } else {
            // Exploit: sort by learned affinity
            let mut affinities: [(f64, usize); 9] = [(0.0, 0); 9];
            for (i, affinity) in affinities.iter_mut().enumerate() {
                *affinity = (self.affinity(i), i);
            }
            affinities.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

            for (priority, (_, step_idx)) in affinities.iter().enumerate() {
                self.priority_buffer[priority] = Self::STEPS[*step_idx];
            }
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
    fn test_initial_exploration_rate() {
        let scheduler = LearnedAffinityScheduler::new(0, 8);
        assert!((scheduler.exploration_rate - 0.3).abs() < 0.001);
    }

    #[test]
    fn test_exploration_decay() {
        let mut scheduler = LearnedAffinityScheduler::new(0, 8);
        let initial_rate = scheduler.current_exploration_rate();

        // Simulate 2000 attempts
        scheduler.total_attempts = 2000;
        let decayed_rate = scheduler.current_exploration_rate();

        assert!(decayed_rate < initial_rate);
        assert!(decayed_rate >= scheduler.min_exploration_rate);
    }

    #[test]
    fn test_affinity_learning() {
        let mut scheduler = LearnedAffinityScheduler::new(0, 8);

        // Build high affinity for Read
        for _ in 0..100 {
            scheduler.record_outcome(PipelineStep::Read, true, false);
        }
        // Build low affinity for Write
        for _ in 0..100 {
            scheduler.record_outcome(PipelineStep::Write, false, false);
        }

        assert!(scheduler.affinity(0) > 0.9); // Read
        assert!(scheduler.affinity(8) < 0.1); // Write
    }

    #[test]
    fn test_minimum_exploration_rate() {
        let mut scheduler = LearnedAffinityScheduler::new(0, 8);
        scheduler.total_attempts = 1_000_000; // Many attempts

        let rate = scheduler.current_exploration_rate();
        assert!((rate - scheduler.min_exploration_rate).abs() < 0.001);
    }
}
