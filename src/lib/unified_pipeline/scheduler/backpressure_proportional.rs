//! Backpressure-Proportional scheduler.
//!
//! Dynamically adjusts step weights based on queue depths using exponential
//! moving averages. When output backs up, downstream steps get higher priority.
//! When input starves, upstream steps get higher priority.

use super::{BackpressureState, Scheduler};
use crate::unified_pipeline::base::PipelineStep;

/// Backpressure-proportional scheduler with EMA smoothing.
pub struct BackpressureProportionalScheduler {
    /// Thread ID.
    thread_id: usize,
    /// Total number of threads.
    num_threads: usize,
    /// EMA weights per step (higher = higher priority).
    weights: [f64; 9],
    /// EMA smoothing factor (0.0-1.0, higher = faster adaptation).
    alpha: f64,
    /// Priority buffer.
    priority_buffer: [PipelineStep; 9],
}

impl BackpressureProportionalScheduler {
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

    /// Default EMA alpha (10% per update).
    const DEFAULT_ALPHA: f64 = 0.1;
    /// High weight value.
    const HIGH_WEIGHT: f64 = 1.0;
    /// Low weight value.
    const LOW_WEIGHT: f64 = 0.3;

    /// Create a new backpressure-proportional scheduler.
    #[must_use]
    pub fn new(thread_id: usize, num_threads: usize) -> Self {
        // Initialize with balanced weights
        let weights = [0.5; 9];

        Self {
            thread_id,
            num_threads,
            weights,
            alpha: Self::DEFAULT_ALPHA,
            priority_buffer: Self::STEPS,
        }
    }

    /// Linearly interpolate between two values.
    fn lerp(current: f64, target: f64, alpha: f64) -> f64 {
        current + alpha * (target - current)
    }
}

impl Scheduler for BackpressureProportionalScheduler {
    fn get_priorities(&mut self, backpressure: BackpressureState) -> &[PipelineStep] {
        // Adjust weights based on backpressure
        if backpressure.output_high {
            // Output backing up: increase priority for downstream steps
            for i in 0..4 {
                self.weights[i] = Self::lerp(self.weights[i], Self::LOW_WEIGHT, self.alpha);
            }
            for i in 4..9 {
                self.weights[i] = Self::lerp(self.weights[i], Self::HIGH_WEIGHT, self.alpha);
            }
        } else if backpressure.input_low && !backpressure.read_done {
            // Input starving: increase priority for upstream steps
            for i in 0..5 {
                self.weights[i] = Self::lerp(self.weights[i], Self::HIGH_WEIGHT, self.alpha);
            }
            for i in 5..9 {
                self.weights[i] = Self::lerp(self.weights[i], Self::LOW_WEIGHT, self.alpha);
            }
        }
        // If neither condition, weights gradually return to balanced

        // Sort steps by weight (descending)
        let mut weighted: [(f64, usize); 9] = [(0.0, 0); 9];
        for (i, (w, weight)) in weighted.iter_mut().zip(self.weights.iter()).enumerate() {
            *w = (*weight, i);
        }
        weighted.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

        for (priority, (_, step_idx)) in weighted.iter().enumerate() {
            self.priority_buffer[priority] = Self::STEPS[*step_idx];
        }

        &self.priority_buffer
    }

    fn record_outcome(&mut self, _step: PipelineStep, _success: bool, _was_contention: bool) {
        // This scheduler doesn't use outcome feedback - purely backpressure-driven
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
    fn test_initial_balanced_weights() {
        let scheduler = BackpressureProportionalScheduler::new(0, 8);
        for weight in scheduler.weights {
            assert!((weight - 0.5).abs() < 0.001);
        }
    }

    #[test]
    fn test_output_high_increases_downstream() {
        let mut scheduler = BackpressureProportionalScheduler::new(0, 8);
        let bp = BackpressureState {
            output_high: true,
            input_low: false,
            read_done: false,
            memory_high: false,
            memory_drained: true,
        };

        // Apply backpressure multiple times
        for _ in 0..10 {
            scheduler.get_priorities(bp);
        }

        // Downstream steps should have higher weights
        assert!(scheduler.weights[8] > scheduler.weights[0]); // Write > Read
    }

    #[test]
    fn test_input_low_increases_upstream() {
        let mut scheduler = BackpressureProportionalScheduler::new(0, 8);
        let bp = BackpressureState {
            output_high: false,
            input_low: true,
            read_done: false,
            memory_high: false,
            memory_drained: true,
        };

        // Apply backpressure multiple times
        for _ in 0..10 {
            scheduler.get_priorities(bp);
        }

        // Upstream steps should have higher weights
        assert!(scheduler.weights[0] > scheduler.weights[8]); // Read > Write
    }
}
