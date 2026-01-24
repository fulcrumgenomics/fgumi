//! Thompson Sampling scheduler.
//!
//! Uses Bayesian inference with Beta distributions to balance exploration
//! and exploitation. Each step maintains a Beta(α, β) distribution where
//! α = successes + 1, β = failures + 1. Steps are prioritized by sampling
//! from their distributions.

use rand::SeedableRng;
use rand::rngs::SmallRng;
use rand_distr::{Beta, Distribution};

use std::sync::atomic::{AtomicU64, Ordering};

use super::{BackpressureState, Scheduler};
use crate::unified_pipeline::base::PipelineStep;

/// Atomic counter for unique seeds across threads.
static SEED_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Thompson Sampling scheduler with Beta distribution priors.
pub struct ThompsonSamplingScheduler {
    /// Thread ID.
    thread_id: usize,
    /// Total number of threads.
    num_threads: usize,
    /// Alpha (success) counts per step - starts at 1.0 for Beta(1,1) uniform prior.
    alphas: [f64; 9],
    /// Beta (failure) counts per step - starts at 1.0.
    betas: [f64; 9],
    /// Random number generator for sampling.
    rng: SmallRng,
    /// Priority buffer for returning step order.
    priority_buffer: [PipelineStep; 9],
}

impl ThompsonSamplingScheduler {
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

    /// Create a new Thompson Sampling scheduler.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize) -> Self {
        // Use a unique seed for each scheduler instance
        let seed = SEED_COUNTER
            .fetch_add(1, Ordering::Relaxed)
            .wrapping_add(thread_id as u64)
            .wrapping_mul(0x9E37_79B9_7F4A_7C15); // Golden ratio hash
        Self {
            thread_id,
            num_threads,
            alphas: [1.0; 9], // Uniform prior Beta(1, 1)
            betas: [1.0; 9],
            rng: SmallRng::seed_from_u64(seed),
            priority_buffer: Self::STEPS,
        }
    }

    /// Get the index of a step.
    fn step_index(step: PipelineStep) -> usize {
        Self::STEPS.iter().position(|&s| s == step).unwrap_or(0)
    }

    /// Sample from Beta distribution, clamping to avoid edge cases.
    fn sample_beta(&mut self, alpha: f64, beta: f64) -> f64 {
        // Clamp to avoid numerical issues with very small/large values
        let alpha = alpha.clamp(0.001, 10000.0);
        let beta = beta.clamp(0.001, 10000.0);

        match Beta::new(alpha, beta) {
            Ok(dist) => dist.sample(&mut self.rng),
            Err(_) => alpha / (alpha + beta), // Fallback to mean
        }
    }
}

impl Scheduler for ThompsonSamplingScheduler {
    fn get_priorities(&mut self, _backpressure: BackpressureState) -> &[PipelineStep] {
        // Sample from each Beta distribution
        let mut samples: [(f64, usize); 9] = [(0.0, 0); 9];
        #[allow(clippy::needless_range_loop)]
        for i in 0..9 {
            samples[i] = (self.sample_beta(self.alphas[i], self.betas[i]), i);
        }

        // Sort by sampled value (descending - higher samples = higher priority)
        samples.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

        // Build priority buffer
        for (priority, (_, step_idx)) in samples.iter().enumerate() {
            self.priority_buffer[priority] = Self::STEPS[*step_idx];
        }

        &self.priority_buffer
    }

    fn record_outcome(&mut self, step: PipelineStep, success: bool, _was_contention: bool) {
        let idx = Self::step_index(step);
        if success {
            self.alphas[idx] += 1.0;
        } else {
            self.betas[idx] += 1.0;
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
    fn test_initial_uniform_prior() {
        let scheduler = ThompsonSamplingScheduler::new(0, 8);
        assert!((scheduler.alphas[0] - 1.0).abs() < f64::EPSILON);
        assert!((scheduler.betas[0] - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_update_on_success() {
        let mut scheduler = ThompsonSamplingScheduler::new(0, 8);
        scheduler.record_outcome(PipelineStep::Read, true, false);
        assert!((scheduler.alphas[0] - 2.0).abs() < f64::EPSILON);
        assert!((scheduler.betas[0] - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_update_on_failure() {
        let mut scheduler = ThompsonSamplingScheduler::new(0, 8);
        scheduler.record_outcome(PipelineStep::Read, false, false);
        assert!((scheduler.alphas[0] - 1.0).abs() < f64::EPSILON);
        assert!((scheduler.betas[0] - 2.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_get_priorities_returns_all_steps() {
        let mut scheduler = ThompsonSamplingScheduler::new(0, 8);
        let bp = BackpressureState::default();
        let priorities = scheduler.get_priorities(bp);
        assert_eq!(priorities.len(), 9);
    }

    #[test]
    fn test_learned_preference() {
        let mut scheduler = ThompsonSamplingScheduler::new(0, 8);
        // Heavily reward Read step
        for _ in 0..100 {
            scheduler.record_outcome(PipelineStep::Read, true, false);
        }
        // Heavily penalize Write step
        for _ in 0..100 {
            scheduler.record_outcome(PipelineStep::Write, false, false);
        }

        // Read should have high alpha, Write should have high beta
        assert!(scheduler.alphas[0] > 50.0); // Read
        assert!(scheduler.betas[8] > 50.0); // Write
    }
}
