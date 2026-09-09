//! Per-OS-thread state board for scheduling telemetry.
//!
//! Each pipeline OS thread (pool worker or detached driver) owns one `AtomicU64`
//! slot holding a packed `(WorkerState, Option<StepIdx>)`. The worker loop stamps
//! it with a single `Relaxed` store per dispatch (only when telemetry is on); the
//! background sampler point-samples every slot each tick to build a thread-state
//! profile. This is statistics, not synchronization — staleness is fine.

use std::sync::atomic::{AtomicU64, Ordering};

use crate::topology::StepIdx;

/// What an OS thread is doing at the instant it is sampled.
///
/// # v1 semantics — which variants the worker loop actually stamps
///
/// Only two variants are stamped at runtime (`runtime/driver.rs`):
/// - `Running` is stamped in `dispatch_one_step` right before a step's
///   `try_run`, so it spans the whole dispatch **whether or not the step made
///   progress**.
/// - `Parked` is stamped on a no-work iteration when the loop backs off.
///
/// `Idle` is only this board's initial pre-run value (see
/// [`WorkerStateBoard::new`]); once work starts it folds into `Parked` and is
/// effectively never sampled. `Waiting` is NOT stamped anywhere in v1: a worker
/// holding an item it cannot push is stamped `Running`, because the held item
/// lives inside the step's output handle and the loop cannot see it. Producer
/// backpressure is instead observable in the EDGES telemetry via `d_push_rej`
/// (push rejections per tick). The two inert variants are kept for schema
/// stability — a telemetry reader must treat `f_waiting`/`f_idle` as
/// non-meaningful in v1.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkerState {
    /// Inside a step's `try_run` (executing `step`), progress or not. The only
    /// "busy" state stamped in v1.
    Running,
    /// Board init value only; folds into `Parked` at runtime (never stamped by
    /// the worker loop once work starts). Inert in v1.
    Idle,
    /// Intended: holding an item it could not push (producer-side backpressure).
    /// NOT stamped in v1 — such a worker reads as `Running`; see the enum-level
    /// note. Retained for schema stability. Backpressure shows up as EDGES
    /// `d_push_rej` instead.
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

/// One `AtomicU64` per pipeline OS thread. Sized at `Pipeline::run` start to
/// `n_pool_workers + n_detached_driver_threads`; each thread is handed a fixed
/// slot index (its "state slot") distinct from `WorkerCore::thread_id` — a
/// detached driver reuses `thread_id` 0, so keying the board on `thread_id`
/// would collide it with pool worker 0.
#[derive(Debug)]
pub struct WorkerStateBoard {
    slots: Box<[AtomicU64]>,
}

impl WorkerStateBoard {
    #[must_use]
    pub fn new(n_slots: usize) -> Self {
        let init = pack(WorkerState::Idle, None);
        let slots = (0..n_slots).map(|_| AtomicU64::new(init)).collect();
        Self { slots }
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.slots.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }

    /// Stamp `slot`'s current `(state, step)` — one `Relaxed` store. Out-of-range
    /// slots are ignored (telemetry must never panic a worker).
    pub fn stamp(&self, slot: usize, state: WorkerState, step: Option<StepIdx>) {
        if let Some(cell) = self.slots.get(slot) {
            cell.store(pack(state, step), Ordering::Relaxed);
        }
    }

    /// Point-read `slot`. Out-of-range reads report `(Idle, None)`.
    #[must_use]
    pub fn read(&self, slot: usize) -> (WorkerState, Option<StepIdx>) {
        self.slots
            .get(slot)
            .map_or((WorkerState::Idle, None), |c| unpack(c.load(Ordering::Relaxed)))
    }
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

    #[test]
    fn board_stamp_and_read_per_slot_are_independent() {
        let board = WorkerStateBoard::new(3);
        assert_eq!(board.len(), 3);
        // All slots start Idle/None.
        for slot in 0..3 {
            assert_eq!(board.read(slot), (WorkerState::Idle, None));
        }
        board.stamp(0, WorkerState::Running, Some(StepIdx(2)));
        board.stamp(2, WorkerState::Parked, None);
        assert_eq!(board.read(0), (WorkerState::Running, Some(StepIdx(2))));
        assert_eq!(board.read(1), (WorkerState::Idle, None), "untouched slot unchanged");
        assert_eq!(board.read(2), (WorkerState::Parked, None));
    }
}
