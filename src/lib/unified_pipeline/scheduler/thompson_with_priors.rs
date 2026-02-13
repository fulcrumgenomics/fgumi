//! Thompson Sampling with thread-specific priors.
//!
//! Like Thompson Sampling, but initializes with biased priors based on
//! thread role. Reader thread (0) has strong prior for Read step,
//! writer thread (N-1) has strong prior for Write step, etc.

use std::sync::atomic::{AtomicU64, Ordering};

use rand::SeedableRng;
use rand::rngs::SmallRng;
use rand_distr::{Beta, Distribution};

use super::{BackpressureState, Scheduler};
use crate::unified_pipeline::base::{ActiveSteps, PipelineStep};

/// Atomic counter for unique seeds across threads.
static SEED_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Thompson Sampling scheduler with thread-specific biased priors.
pub struct ThompsonWithPriorsScheduler {
    /// Thread ID.
    thread_id: usize,
    /// Total number of threads.
    num_threads: usize,
    /// Alpha (success) counts per step.
    alphas: [f64; 9],
    /// Beta (failure) counts per step.
    betas: [f64; 9],
    /// Random number generator.
    rng: SmallRng,
    /// Priority buffer.
    priority_buffer: [PipelineStep; 9],
    /// Active steps in the pipeline.
    active_steps: ActiveSteps,
}

impl ThompsonWithPriorsScheduler {
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

    /// Prior strength for specialized threads.
    const STRONG_PRIOR: f64 = 20.0;
    /// Prior strength for secondary preferences.
    const MEDIUM_PRIOR: f64 = 5.0;
    /// Base prior (uniform).
    const BASE_PRIOR: f64 = 1.0;

    /// Create a new Thompson Sampling scheduler with thread-specific priors.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize, active_steps: ActiveSteps) -> Self {
        let mut alphas = [Self::BASE_PRIOR; 9];
        let betas = [Self::BASE_PRIOR; 9];

        // Assign priors based on thread role
        if thread_id == 0 {
            // Reader thread: strong prior for Read, medium for Decompress
            alphas[0] = Self::STRONG_PRIOR; // Read
            alphas[1] = Self::MEDIUM_PRIOR; // Decompress
        } else if thread_id == num_threads - 1 && num_threads > 1 {
            // Writer thread: strong prior for Write, medium for Compress
            alphas[8] = Self::STRONG_PRIOR; // Write
            alphas[7] = Self::MEDIUM_PRIOR; // Compress
        } else if thread_id == 1 && num_threads > 2 {
            // Boundary thread: strong prior for FindBoundaries
            alphas[2] = Self::STRONG_PRIOR; // FindBoundaries
            alphas[3] = Self::MEDIUM_PRIOR; // Decode
        } else if thread_id == num_threads - 2 && num_threads > 3 {
            // Group thread: strong prior for Group
            alphas[4] = Self::STRONG_PRIOR; // Group
            alphas[5] = Self::MEDIUM_PRIOR; // Process
        } else {
            // Middle threads: spread priors across parallel steps
            let parallel_steps = [1, 3, 5, 6, 7]; // Decompress, Decode, Process, Serialize, Compress
            let offset = (thread_id.saturating_sub(2)) % parallel_steps.len();
            alphas[parallel_steps[offset]] = Self::MEDIUM_PRIOR;
        }

        let seed = SEED_COUNTER
            .fetch_add(1, Ordering::Relaxed)
            .wrapping_add(thread_id as u64)
            .wrapping_mul(0x9E37_79B9_7F4A_7C15);
        Self {
            thread_id,
            num_threads,
            alphas,
            betas,
            rng: SmallRng::seed_from_u64(seed),
            priority_buffer: Self::STEPS,
            active_steps,
        }
    }

    /// Get the index of a step.
    fn step_index(step: PipelineStep) -> usize {
        Self::STEPS.iter().position(|&s| s == step).unwrap_or(0)
    }

    /// Sample from Beta distribution.
    fn sample_beta(&mut self, alpha: f64, beta: f64) -> f64 {
        let alpha = alpha.clamp(0.001, 10000.0);
        let beta = beta.clamp(0.001, 10000.0);

        match Beta::new(alpha, beta) {
            Ok(dist) => dist.sample(&mut self.rng),
            Err(_) => alpha / (alpha + beta),
        }
    }
}

impl Scheduler for ThompsonWithPriorsScheduler {
    fn get_priorities(&mut self, _backpressure: BackpressureState) -> &[PipelineStep] {
        let mut samples: [(f64, usize); 9] = [(0.0, 0); 9];
        #[allow(clippy::needless_range_loop)]
        for i in 0..9 {
            samples[i] = (self.sample_beta(self.alphas[i], self.betas[i]), i);
        }

        samples.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

        for (priority, (_, step_idx)) in samples.iter().enumerate() {
            self.priority_buffer[priority] = Self::STEPS[*step_idx];
        }

        let n = self.active_steps.filter_in_place(&mut self.priority_buffer);
        &self.priority_buffer[..n]
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

    fn active_steps(&self) -> &ActiveSteps {
        &self.active_steps
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reader_thread_prior() {
        let scheduler = ThompsonWithPriorsScheduler::new(0, 8, ActiveSteps::all());
        assert!(scheduler.alphas[0] > 10.0); // Read has strong prior
        assert!(scheduler.alphas[8] < 5.0); // Write has weak prior
    }

    #[test]
    fn test_writer_thread_prior() {
        let scheduler = ThompsonWithPriorsScheduler::new(7, 8, ActiveSteps::all());
        assert!(scheduler.alphas[8] > 10.0); // Write has strong prior
        assert!(scheduler.alphas[0] < 5.0); // Read has weak prior
    }

    #[test]
    fn test_boundary_thread_prior() {
        let scheduler = ThompsonWithPriorsScheduler::new(1, 8, ActiveSteps::all());
        assert!(scheduler.alphas[2] > 10.0); // FindBoundaries has strong prior
    }

    #[test]
    fn test_group_thread_prior() {
        let scheduler = ThompsonWithPriorsScheduler::new(6, 8, ActiveSteps::all()); // N-2 = 6
        assert!(scheduler.alphas[4] > 10.0); // Group has strong prior
    }
}
