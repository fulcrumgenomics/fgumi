//! Periodic thread/edge/summary telemetry: config, per-worker sample binning,
//! and the four-file TSV writer. Driven by the occupancy sampler (`sampler.rs`).

use std::path::PathBuf;
use std::time::Duration;

use crate::runtime::worker_state::WorkerState;
use crate::topology::StepIdx;

/// Where + how often to emit the tick telemetry.
#[derive(Debug, Clone)]
pub struct TelemetryConfig {
    /// File stem; the four TSVs are `<stem>.ticks.{summary,edges,workers,steps}.tsv`.
    pub stem: PathBuf,
    /// Emit/sample cadence.
    pub interval: Duration,
}

/// Fraction of an emit window a worker spent in each step / non-running state.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowFractions {
    pub per_step: Vec<f32>,
    pub idle: f32,
    pub waiting: f32,
    pub parked: f32,
}

/// Accumulates point-samples for ONE worker over ONE emit window. The sampler
/// calls `record` once per tick per worker, then `fractions` at emit time and
/// `reset` to start the next window.
#[derive(Debug)]
pub struct WorkerBins {
    per_step: Vec<u64>,
    idle: u64,
    waiting: u64,
    parked: u64,
    total: u64,
}

impl WorkerBins {
    #[must_use]
    pub fn new(n_steps: usize) -> Self {
        Self { per_step: vec![0; n_steps], idle: 0, waiting: 0, parked: 0, total: 0 }
    }

    pub fn record(&mut self, state: WorkerState, step: Option<StepIdx>) {
        self.total += 1;
        match (state, step) {
            (WorkerState::Running, Some(idx)) if idx.0 < self.per_step.len() => {
                self.per_step[idx.0] += 1;
            }
            // Running with no/out-of-range step is folded into idle (defensive).
            (WorkerState::Running | WorkerState::Idle, _) => self.idle += 1,
            (WorkerState::Waiting, _) => self.waiting += 1,
            (WorkerState::Parked, _) => self.parked += 1,
        }
    }

    #[must_use]
    pub fn total(&self) -> u64 {
        self.total
    }

    #[allow(clippy::cast_precision_loss)]
    #[must_use]
    pub fn fractions(&self) -> WindowFractions {
        if self.total == 0 {
            return WindowFractions {
                per_step: vec![0.0; self.per_step.len()],
                idle: 0.0,
                waiting: 0.0,
                parked: 0.0,
            };
        }
        let t = self.total as f32;
        WindowFractions {
            per_step: self.per_step.iter().map(|&c| c as f32 / t).collect(),
            idle: self.idle as f32 / t,
            waiting: self.waiting as f32 / t,
            parked: self.parked as f32 / t,
        }
    }

    pub fn reset(&mut self) {
        self.per_step.iter_mut().for_each(|c| *c = 0);
        self.idle = 0;
        self.waiting = 0;
        self.parked = 0;
        self.total = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_fractions_sum_to_one_and_split_states() {
        let mut bins = WorkerBins::new(2); // steps 0,1
        bins.record(WorkerState::Running, Some(StepIdx(1)));
        bins.record(WorkerState::Running, Some(StepIdx(1)));
        bins.record(WorkerState::Parked, None);
        bins.record(WorkerState::Idle, None);
        assert_eq!(bins.total(), 4);
        let f = bins.fractions();
        assert_eq!(f.per_step.len(), 2);
        assert!((f.per_step[0] - 0.0).abs() < 1e-6);
        assert!((f.per_step[1] - 0.5).abs() < 1e-6);
        assert!((f.parked - 0.25).abs() < 1e-6);
        assert!((f.idle - 0.25).abs() < 1e-6);
        assert!((f.waiting - 0.0).abs() < 1e-6);
        let sum: f32 = f.per_step.iter().sum::<f32>() + f.idle + f.waiting + f.parked;
        assert!((sum - 1.0).abs() < 1e-5, "fractions sum to 1, got {sum}");
    }

    #[test]
    fn empty_window_is_all_zero() {
        let bins = WorkerBins::new(3);
        assert_eq!(bins.total(), 0);
        let f = bins.fractions();
        assert!(f.per_step.iter().all(|&x| x == 0.0));
        assert_eq!((f.idle, f.waiting, f.parked), (0.0, 0.0, 0.0));
    }

    #[test]
    fn running_with_no_step_counts_as_idle_bucket() {
        // Running but step==None (shouldn't happen, but be defensive) must not
        // index out of bounds; fold it into idle.
        let mut bins = WorkerBins::new(1);
        bins.record(WorkerState::Running, None);
        let f = bins.fractions();
        assert!((f.idle - 1.0).abs() < 1e-6);
    }
}
