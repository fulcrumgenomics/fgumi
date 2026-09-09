//! Per-OS-thread state board for scheduling telemetry.
//!
//! Each pipeline OS thread (pool worker or detached driver) owns one `AtomicU64`
//! slot holding a packed `(WorkerState, Option<StepIdx>)`. The worker loop stamps
//! it with a single `Relaxed` store per dispatch (only when telemetry is on); the
//! background sampler point-samples every slot each tick to build a thread-state
//! profile. This is statistics, not synchronization — staleness is fine.

use crate::topology::StepIdx;

/// What an OS thread is doing at the instant it is sampled.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkerState {
    /// Inside a step's `try_run` (executing `step`).
    Running,
    /// Between dispatches / round-robin walk found no work this pass.
    Idle,
    /// Holding an item it could not push (producer-side backpressure).
    Waiting,
    /// Sleeping/parked on the backoff (no work; ramped backoff).
    Parked,
}

// State occupies the top 3 bits; the step index occupies the low 32 bits.
// `Option<StepIdx>` uses a sentinel (`u32::MAX`) for `None`, so step indices are
// capped at `u32::MAX - 1` (far above any real chain length).
const STATE_SHIFT: u64 = 61;
const STEP_MASK: u64 = 0xFFFF_FFFF;
const NO_STEP: u64 = u32::MAX as u64;

#[must_use]
pub fn pack(state: WorkerState, step: Option<StepIdx>) -> u64 {
    let s = match state {
        WorkerState::Running => 0u64,
        WorkerState::Idle => 1,
        WorkerState::Waiting => 2,
        WorkerState::Parked => 3,
    };
    let step_bits = match step {
        Some(idx) => {
            debug_assert!((idx.0 as u64) < NO_STEP, "step index collides with the None sentinel");
            (idx.0 as u64) & STEP_MASK
        }
        None => NO_STEP,
    };
    (s << STATE_SHIFT) | step_bits
}

#[must_use]
pub fn unpack(bits: u64) -> (WorkerState, Option<StepIdx>) {
    let state = match bits >> STATE_SHIFT {
        0 => WorkerState::Running,
        1 => WorkerState::Idle,
        2 => WorkerState::Waiting,
        _ => WorkerState::Parked,
    };
    let step_bits = bits & STEP_MASK;
    let step = if step_bits == NO_STEP { None } else { Some(StepIdx(step_bits as usize)) };
    (state, step)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case::running_step0(WorkerState::Running, Some(StepIdx(0)))]
    #[case::running_high(WorkerState::Running, Some(StepIdx(4_000_000_000)))]
    #[case::idle_none(WorkerState::Idle, None)]
    #[case::waiting_step7(WorkerState::Waiting, Some(StepIdx(7)))]
    #[case::parked_none(WorkerState::Parked, None)]
    fn pack_round_trips(#[case] state: WorkerState, #[case] step: Option<StepIdx>) {
        assert_eq!(unpack(pack(state, step)), (state, step));
    }
}
