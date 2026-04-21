//! Multi-stage scheduler: priority assignment and worker affinity.
//!
//! The scheduler is stateless per iteration: given a [`BackpressureState`],
//! it returns a priority-ordered list of stage indices for a worker to try.
//!
//! Strategy: `BalancedChaseDrain` (generalized).
//! - Sticky: prefer the worker's last successful stage.
//! - Chase downstream: when a stage's output queue is high, prioritize
//!   the downstream stage to drain.
//! - Chase upstream: when a stage's input queue is low, prioritize
//!   the upstream stage to feed.
//!
//! The scheduler is live production code: each pool worker calls
//! [`priorities_for_with_role`] once per loop iteration to choose the
//! stage to attempt next.

use super::backpressure::BackpressureState;
use super::stage::Parallelism;

/// Worker role hint for priority assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum WorkerRole {
    /// No role affinity (middle worker).
    #[default]
    Generic,
    /// Prefer the head of the pipeline (feed the source-adjacent stage).
    HeadAffinity,
    /// Prefer the tail of the pipeline (drain the sink-adjacent stage).
    TailAffinity,
}

/// Assign exclusive ownership of Sequential stages to workers.
///
/// Returns a `Vec<Vec<usize>>` indexed by worker id. Each inner vec lists
/// the Sequential stage indices that worker owns. Assignment alternates
/// from the ends: Sequential[0] → worker 0, Sequential[1] → worker N-1,
/// Sequential[2] → worker 1, Sequential[3] → worker N-2, ...
///
/// # Panics
/// Panics if the number of Sequential stages exceeds `num_workers` — a
/// configuration error (each Sequential stage needs its own owner).
#[must_use]
pub fn assign_sequential_ownership(
    num_workers: usize,
    parallelism: &[Parallelism],
) -> Vec<Vec<usize>> {
    assert!(num_workers > 0, "num_workers must be >= 1");
    let sequential_indices: Vec<usize> = parallelism
        .iter()
        .enumerate()
        .filter_map(|(i, p)| (*p == Parallelism::Sequential).then_some(i))
        .collect();
    assert!(
        sequential_indices.len() <= num_workers,
        "Sequential stage count ({}) exceeds worker count ({})",
        sequential_indices.len(),
        num_workers,
    );

    let mut owned: Vec<Vec<usize>> = (0..num_workers).map(|_| Vec::new()).collect();
    for (i, stage_idx) in sequential_indices.iter().enumerate() {
        // Alternate from the ends: 0, N-1, 1, N-2, 2, N-3, ...
        let worker_id = if i % 2 == 0 { i / 2 } else { num_workers - 1 - (i / 2) };
        owned[worker_id].push(*stage_idx);
    }
    owned
}

/// Outcome of a single stage attempt, reduced to the variants the scheduler
/// uses to update its internal state. The pool's richer `StageAttemptResult`
/// is translated into this by the worker loop; variants like `Error` and
/// `Cancelled` terminate the worker before any scheduler update.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StageAttemptOutcome {
    /// The stage produced an output (or consumed an input via
    /// suppress-on-empty) — forward progress was made.
    Success,
    /// Stage's input queue was empty but not closed.
    BackpressureInput,
    /// Stage's output queue was full; the item was held for later.
    BackpressureOutput,
    /// Sequential stage was held by another worker.
    SequentialBusy,
    /// Stage's input is closed and empty — no more work possible here.
    Drained,
}

/// Drain-mode state: track which producer has an overflowing output queue
/// so its consumer can be prioritized and the producer itself deprioritized
/// until the queue drops below the hysteresis threshold.
#[derive(Debug, Clone, Copy)]
struct DrainModeState {
    /// Stage index that produced a `BackpressureOutput`; its output queue
    /// is the one we watch for hysteresis exit.
    stuck_producer: usize,
    /// Queue index to watch for drain hysteresis. Always
    /// `stuck_producer + 1` (the producer's output queue).
    stuck_queue: usize,
}

/// Per-worker scheduler state.
///
/// Owns every decision-affecting piece of state: sticky last-success, drain
/// mode, exclusive ownership of Sequential stages, and role affinity. The
/// pool drives it by calling [`Scheduler::priorities`] once per iteration to
/// pick a stage order, then [`Scheduler::record_outcome`] after each stage
/// attempt to update state.
pub struct Scheduler {
    num_stages: usize,
    role: WorkerRole,
    /// Sequential stage indices this worker exclusively owns (see
    /// `assign_sequential_ownership`). Typically 0–2 entries.
    exclusive_owner_of: Vec<usize>,
    /// Last stage this worker succeeded at. Prioritized next iteration.
    sticky: Option<usize>,
    /// Active drain-mode state, if any. Cleared by
    /// `maybe_clear_drain_mode` when the stuck queue's `memory_fill` drops
    /// below 0.5.
    drain_mode: Option<DrainModeState>,
}

impl Scheduler {
    /// Construct a scheduler for one pool worker.
    #[must_use]
    pub fn new(num_stages: usize, role: WorkerRole, exclusive_owner_of: Vec<usize>) -> Self {
        Self { num_stages, role, exclusive_owner_of, sticky: None, drain_mode: None }
    }

    /// This worker's role affinity.
    #[must_use]
    pub fn role(&self) -> WorkerRole {
        self.role
    }

    /// Whether this worker exclusively owns a given Sequential stage.
    #[must_use]
    pub fn owns(&self, stage_idx: usize) -> bool {
        self.exclusive_owner_of.contains(&stage_idx)
    }

    /// Slice of Sequential stage indices this worker exclusively owns.
    /// Used by the pool's eager-exclusive-step execution before the
    /// priority list runs.
    #[must_use]
    pub fn owned_stages(&self) -> &[usize] {
        &self.exclusive_owner_of
    }

    /// Whether this worker should attempt the given Sequential stage.
    /// Always true in drain mode (any worker helps drain); otherwise only
    /// true for the exclusive owner.
    #[must_use]
    pub fn should_attempt_sequential(&self, stage_idx: usize) -> bool {
        self.drain_mode.is_some() || self.owns(stage_idx)
    }

    /// Clear drain mode if the stuck queue's `memory_fill` has dropped
    /// below 0.5 (hysteresis threshold). Called at the top of each
    /// `priorities` invocation.
    fn maybe_clear_drain_mode(&mut self, bp: &BackpressureState) {
        if let Some(dm) = self.drain_mode {
            let cleared =
                bp.queue_summaries.get(dm.stuck_queue).is_some_and(|q| q.memory_fill < 0.5);
            if cleared {
                self.drain_mode = None;
            }
        }
    }

    /// Priority-ordered list of stage indices this worker should try this
    /// iteration.
    #[must_use]
    pub fn priorities(&mut self, bp: &BackpressureState) -> Vec<usize> {
        self.maybe_clear_drain_mode(bp);

        let mut priorities = Vec::with_capacity(self.num_stages);

        // 1. Memory-pressure drain overrides everything: prefer
        //    downstream-most first.
        if bp.memory_high {
            for i in (0..self.num_stages).rev() {
                priorities.push(i);
            }
            return priorities;
        }

        // 2. Role affinity: head/tail workers prefer their endpoints.
        match self.role {
            WorkerRole::HeadAffinity if self.num_stages > 0 => priorities.push(0),
            WorkerRole::TailAffinity if self.num_stages > 0 => {
                priorities.push(self.num_stages - 1);
            }
            _ => {}
        }

        // 3. Sticky: worker's last successful stage.
        if let Some(sticky) = self.sticky
            && sticky < self.num_stages
            && !priorities.contains(&sticky)
        {
            priorities.push(sticky);
        }

        // 4. Chase downstream: high queues -> run the consumer to drain.
        for (qid, q) in bp.queue_summaries.iter().enumerate() {
            if q.is_high() && qid < self.num_stages && !priorities.contains(&qid) {
                priorities.push(qid);
            }
        }

        // 5. Chase upstream: low non-empty queues -> run the producer to feed.
        for (qid, q) in bp.queue_summaries.iter().enumerate() {
            if q.is_low() && !q.empty && qid > 0 {
                let producer = qid - 1;
                if producer < self.num_stages && !priorities.contains(&producer) {
                    priorities.push(producer);
                }
            }
        }

        // 6. Fill in remaining stages in order.
        for i in 0..self.num_stages {
            if !priorities.contains(&i) {
                priorities.push(i);
            }
        }

        // 7. Drain-mode skip: remove the stuck producer from the list.
        //    If that leaves the list empty (single-stage pipeline edge
        //    case), append the stuck producer so the worker still has
        //    something to attempt.
        if let Some(dm) = self.drain_mode {
            priorities.retain(|&idx| idx != dm.stuck_producer);
            if priorities.is_empty() && dm.stuck_producer < self.num_stages {
                priorities.push(dm.stuck_producer);
            }
        }

        priorities
    }

    /// Priority-ordered stages with Sequential stages filtered to those
    /// this worker owns (or all, in drain mode). `parallelism` has one
    /// entry per stage matching `num_stages`.
    ///
    /// This is the method `pool.rs` calls; `priorities` is kept as the
    /// low-level primitive used by tests and the drain-mode fallback.
    #[must_use]
    pub fn priorities_filtered(
        &mut self,
        bp: &BackpressureState,
        parallelism: &[Parallelism],
    ) -> Vec<usize> {
        let base = self.priorities(bp);
        base.into_iter()
            .filter(|&idx| match parallelism.get(idx) {
                Some(Parallelism::Sequential) => self.should_attempt_sequential(idx),
                _ => true,
            })
            .collect()
    }

    /// Whether this worker's scheduler is currently in drain mode.
    #[must_use]
    pub fn is_draining(&self) -> bool {
        self.drain_mode.is_some()
    }

    /// Update scheduler state after a single stage attempt.
    ///
    /// Sticky update rules:
    /// - `Success` → `sticky = Some(stage_idx)`
    /// - `BackpressureInput`, `BackpressureOutput`, `Drained` → clear
    ///   sticky if it pointed at `stage_idx`
    /// - `SequentialBusy` → no change
    ///
    /// (Drain-mode updates are added in a later task.)
    pub fn record_outcome(
        &mut self,
        stage_idx: usize,
        outcome: StageAttemptOutcome,
        _bp: &BackpressureState,
    ) {
        match outcome {
            StageAttemptOutcome::Success => {
                self.sticky = Some(stage_idx);
            }
            StageAttemptOutcome::BackpressureInput | StageAttemptOutcome::Drained => {
                if self.sticky == Some(stage_idx) {
                    self.sticky = None;
                }
            }
            StageAttemptOutcome::BackpressureOutput => {
                if self.sticky == Some(stage_idx) {
                    self.sticky = None;
                }
                // Deprioritize this producer until its output queue
                // drops below hysteresis. `stuck_queue = stage_idx + 1`
                // is the producer's output queue by definition.
                self.drain_mode =
                    Some(DrainModeState { stuck_producer: stage_idx, stuck_queue: stage_idx + 1 });
            }
            StageAttemptOutcome::SequentialBusy => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::backpressure::QueueSummary;
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::stage::Parallelism;
    use std::sync::Arc;

    fn mk_bp(queues: Vec<QueueSummary>) -> BackpressureState {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        BackpressureState::new(queues, &tracker)
    }

    fn mk_q(id: usize, slot: f64, mem: f64) -> QueueSummary {
        QueueSummary { id, slot_fill: slot, memory_fill: mem, closed: false, empty: false }
    }

    #[test]
    fn test_sticky_affinity_first() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5)]);
        let mut s = Scheduler::new(2, WorkerRole::Generic, vec![]);
        s.sticky = Some(1);
        let prio = s.priorities(&bp);
        assert_eq!(prio[0], 1); // sticky first
    }

    #[test]
    fn test_chase_downstream_on_high_queue() {
        // 3 stages, queue 1 (between stage 0 and 1) is high. Worker without
        // affinity should prioritize stage 1 (the consumer of queue 1).
        let bp = mk_bp(vec![
            mk_q(0, 0.5, 0.5),
            mk_q(1, 0.9, 0.5), // high
            mk_q(2, 0.5, 0.5),
        ]);
        let mut s = Scheduler::new(3, WorkerRole::Generic, vec![]);
        let prio = s.priorities(&bp);
        assert_eq!(prio[0], 1); // stage 1 first
    }

    #[test]
    fn test_all_stages_covered() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5), mk_q(2, 0.5, 0.5)]);
        let mut s = Scheduler::new(3, WorkerRole::Generic, vec![]);
        let prio = s.priorities(&bp);
        assert_eq!(prio.len(), 3);
        let mut sorted = prio.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2]);
    }

    #[test]
    fn test_sticky_combined_with_chase() {
        // Sticky on stage 0, but queue 2 is high (consumer stage 2).
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5), mk_q(2, 0.9, 0.5)]);
        let mut s = Scheduler::new(3, WorkerRole::Generic, vec![]);
        s.sticky = Some(0);
        let prio = s.priorities(&bp);
        assert_eq!(prio[0], 0); // sticky
        assert_eq!(prio[1], 2); // then drain
    }

    #[test]
    fn test_drain_mode_prioritizes_downstream_when_memory_high() {
        let tracker = Arc::new(MemoryTracker::new(100));
        tracker.add(150); // over limit
        let bp = BackpressureState::new(
            vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5), mk_q(2, 0.5, 0.5)],
            &tracker,
        );
        let mut s = Scheduler::new(3, WorkerRole::Generic, vec![]);
        let prio = s.priorities(&bp);
        // Under memory pressure: prefer downstream-most stages.
        assert_eq!(prio[0], 2); // last stage first
    }

    #[test]
    fn test_role_affinity_worker_zero_prefers_stage_zero() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5), mk_q(2, 0.5, 0.5)]);
        let mut s = Scheduler::new(3, WorkerRole::HeadAffinity, vec![]);
        let prio = s.priorities(&bp);
        assert_eq!(prio[0], 0); // head-affinity worker prefers stage 0
    }

    #[test]
    fn test_role_affinity_worker_last_prefers_stage_last() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5), mk_q(2, 0.5, 0.5)]);
        let mut s = Scheduler::new(3, WorkerRole::TailAffinity, vec![]);
        let prio = s.priorities(&bp);
        assert_eq!(prio[0], 2); // tail-affinity worker prefers last stage
    }

    #[test]
    fn test_scheduler_new_defaults() {
        let s = Scheduler::new(3, WorkerRole::Generic, vec![]);
        assert_eq!(s.role(), WorkerRole::Generic);
        assert!(s.sticky.is_none());
        assert!(s.drain_mode.is_none());
        assert!(s.exclusive_owner_of.is_empty());
    }

    #[test]
    fn test_scheduler_owns_returns_true_for_owned_stage() {
        let s = Scheduler::new(4, WorkerRole::HeadAffinity, vec![1, 3]);
        assert!(s.owns(1));
        assert!(s.owns(3));
        assert!(!s.owns(0));
        assert!(!s.owns(2));
    }

    #[test]
    fn test_priorities_method_covers_all_stages() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5), mk_q(2, 0.5, 0.5)]);
        let mut sched = Scheduler::new(3, WorkerRole::Generic, vec![]);
        let prio = sched.priorities(&bp);
        assert_eq!(prio.len(), 3);
        let mut sorted = prio.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2]);
    }

    #[test]
    fn test_record_outcome_success_sets_sticky() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5)]);
        let mut s = Scheduler::new(2, WorkerRole::Generic, vec![]);
        s.record_outcome(1, StageAttemptOutcome::Success, &bp);
        assert_eq!(s.sticky, Some(1));
    }

    #[test]
    fn test_record_outcome_backpressure_input_clears_matching_sticky() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5)]);
        let mut s = Scheduler::new(2, WorkerRole::Generic, vec![]);
        s.sticky = Some(1);
        s.record_outcome(1, StageAttemptOutcome::BackpressureInput, &bp);
        assert!(s.sticky.is_none());
    }

    #[test]
    fn test_record_outcome_backpressure_input_leaves_nonmatching_sticky() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5)]);
        let mut s = Scheduler::new(2, WorkerRole::Generic, vec![]);
        s.sticky = Some(0);
        s.record_outcome(1, StageAttemptOutcome::BackpressureInput, &bp);
        assert_eq!(s.sticky, Some(0));
    }

    #[test]
    fn test_record_outcome_sequential_busy_does_not_change_sticky() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5)]);
        let mut s = Scheduler::new(2, WorkerRole::Generic, vec![]);
        s.sticky = Some(0);
        s.record_outcome(1, StageAttemptOutcome::SequentialBusy, &bp);
        assert_eq!(s.sticky, Some(0));
    }

    #[test]
    fn test_record_outcome_drained_clears_matching_sticky() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5)]);
        let mut s = Scheduler::new(2, WorkerRole::Generic, vec![]);
        s.sticky = Some(0);
        s.record_outcome(0, StageAttemptOutcome::Drained, &bp);
        assert!(s.sticky.is_none());
    }

    #[test]
    fn test_drain_mode_activates_on_backpressure_output() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.9, 0.9), mk_q(2, 0.5, 0.5)]);
        let mut s = Scheduler::new(2, WorkerRole::Generic, vec![]);
        s.record_outcome(0, StageAttemptOutcome::BackpressureOutput, &bp);
        // stuck_producer = 0, stuck_queue = 1
        assert!(s.drain_mode.is_some());
        let dm = s.drain_mode.unwrap();
        assert_eq!(dm.stuck_producer, 0);
        assert_eq!(dm.stuck_queue, 1);
    }

    #[test]
    fn test_drain_mode_deprioritizes_stuck_producer() {
        // Stage 0's output queue (queue 1) is full. After the output-push
        // failure, priorities must exclude stage 0 from the main list —
        // other stages have priority, and stage 0 only appears as the
        // last-resort fallback.
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.9, 0.9), mk_q(2, 0.5, 0.5)]);
        let mut s = Scheduler::new(3, WorkerRole::Generic, vec![]);
        s.record_outcome(0, StageAttemptOutcome::BackpressureOutput, &bp);
        let prio = s.priorities(&bp);
        // Stage 0 must NOT be first.
        assert_ne!(prio[0], 0);
        // Stage 1 (consumer of the stuck queue) should be up front because
        // queue 1 is_high() chases it.
        assert_eq!(prio[0], 1);
    }

    #[test]
    fn test_drain_mode_clears_when_queue_below_50pct() {
        let bp_high = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.9, 0.9), mk_q(2, 0.5, 0.5)]);
        let bp_low = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.4, 0.4), mk_q(2, 0.5, 0.5)]);
        let mut s = Scheduler::new(3, WorkerRole::Generic, vec![]);
        s.record_outcome(0, StageAttemptOutcome::BackpressureOutput, &bp_high);
        assert!(s.drain_mode.is_some());
        // Calling priorities with a backpressure where stuck_queue is
        // drained must clear drain mode.
        let _ = s.priorities(&bp_low);
        assert!(s.drain_mode.is_none());
    }

    #[test]
    fn test_drain_mode_stays_above_50pct() {
        let bp_high = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.9, 0.9), mk_q(2, 0.5, 0.5)]);
        let bp_mid = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.6, 0.6), mk_q(2, 0.5, 0.5)]);
        let mut s = Scheduler::new(3, WorkerRole::Generic, vec![]);
        s.record_outcome(0, StageAttemptOutcome::BackpressureOutput, &bp_high);
        let _ = s.priorities(&bp_mid);
        assert!(s.drain_mode.is_some(), "drain mode must persist above 50% fill");
    }

    #[test]
    fn test_drain_mode_fallback_appends_stuck_producer_if_list_would_be_empty() {
        // Edge case: single-stage pipeline (num_stages=1), drain mode
        // active on stage 0. Without the fallback, the list would be
        // empty. The fallback must append stage 0.
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.9, 0.9)]);
        let mut s = Scheduler::new(1, WorkerRole::Generic, vec![]);
        s.record_outcome(0, StageAttemptOutcome::BackpressureOutput, &bp);
        let prio = s.priorities(&bp);
        assert_eq!(prio, vec![0]);
    }

    #[test]
    fn test_should_attempt_sequential_owner_sees_it() {
        let s = Scheduler::new(3, WorkerRole::Generic, vec![1]);
        assert!(s.should_attempt_sequential(1));
    }

    #[test]
    fn test_should_attempt_sequential_nonowner_does_not() {
        let s = Scheduler::new(3, WorkerRole::Generic, vec![]);
        assert!(!s.should_attempt_sequential(1));
    }

    #[test]
    fn test_should_attempt_sequential_drain_mode_relaxes() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.9, 0.9), mk_q(2, 0.5, 0.5)]);
        let mut s = Scheduler::new(3, WorkerRole::Generic, vec![]);
        s.record_outcome(0, StageAttemptOutcome::BackpressureOutput, &bp);
        // Non-owner, but in drain mode — relaxed.
        assert!(s.should_attempt_sequential(1));
    }

    #[test]
    fn test_priorities_filters_sequential_for_nonowner() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5), mk_q(2, 0.5, 0.5)]);
        let mut s = Scheduler::new(3, WorkerRole::Generic, vec![]);
        let kinds = vec![Parallelism::Parallel, Parallelism::Sequential, Parallelism::Parallel];
        let prio = s.priorities_filtered(&bp, &kinds);
        assert!(!prio.contains(&1), "non-owner must not see Sequential stage 1");
        assert_eq!(prio.len(), 2);
    }

    #[test]
    fn test_priorities_includes_sequential_for_owner() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.5, 0.5), mk_q(2, 0.5, 0.5)]);
        let mut s = Scheduler::new(3, WorkerRole::Generic, vec![1]);
        let kinds = vec![Parallelism::Parallel, Parallelism::Sequential, Parallelism::Parallel];
        let prio = s.priorities_filtered(&bp, &kinds);
        assert!(prio.contains(&1));
        assert_eq!(prio.len(), 3);
    }

    #[test]
    fn test_priorities_drain_mode_allows_sequential_for_nonowner() {
        let bp = mk_bp(vec![mk_q(0, 0.5, 0.5), mk_q(1, 0.9, 0.9), mk_q(2, 0.5, 0.5)]);
        let mut s = Scheduler::new(3, WorkerRole::Generic, vec![]);
        // Stage 1 is Sequential; simulate its producer (stage 0) hitting
        // BackpressureOutput so drain mode activates.
        s.record_outcome(0, StageAttemptOutcome::BackpressureOutput, &bp);
        let kinds = vec![Parallelism::Parallel, Parallelism::Sequential, Parallelism::Parallel];
        let prio = s.priorities_filtered(&bp, &kinds);
        // In drain mode, non-owner may see Sequential stage 1.
        assert!(prio.contains(&1), "drain mode must relax ownership filter");
    }

    #[test]
    fn test_assign_sequential_ownership_zero_sequential() {
        let kinds = vec![Parallelism::Parallel, Parallelism::Parallel];
        let out = assign_sequential_ownership(4, &kinds);
        assert_eq!(out.len(), 4);
        for worker_owned in &out {
            assert!(worker_owned.is_empty());
        }
    }

    #[test]
    fn test_assign_sequential_ownership_one_sequential_goes_to_worker_zero() {
        let kinds = vec![Parallelism::Parallel, Parallelism::Sequential, Parallelism::Parallel];
        let out = assign_sequential_ownership(4, &kinds);
        assert_eq!(out[0], vec![1]);
        assert!(out[1].is_empty());
        assert!(out[2].is_empty());
        assert!(out[3].is_empty());
    }

    #[test]
    fn test_assign_sequential_ownership_two_sequential_alternates_ends() {
        let kinds = vec![
            Parallelism::Sequential,
            Parallelism::Parallel,
            Parallelism::Sequential,
            Parallelism::Parallel,
        ];
        let out = assign_sequential_ownership(4, &kinds);
        // Sequential indices are [0, 2]. First goes to worker 0, second
        // goes to worker 3 (N-1).
        assert_eq!(out[0], vec![0]);
        assert!(out[1].is_empty());
        assert!(out[2].is_empty());
        assert_eq!(out[3], vec![2]);
    }

    #[test]
    fn test_assign_sequential_ownership_four_sequential_alternates_fully() {
        let kinds = vec![
            Parallelism::Sequential,
            Parallelism::Sequential,
            Parallelism::Sequential,
            Parallelism::Sequential,
        ];
        let out = assign_sequential_ownership(4, &kinds);
        // Sequential indices [0, 1, 2, 3] -> workers [0, 3, 1, 2].
        assert_eq!(out[0], vec![0]);
        assert_eq!(out[3], vec![1]);
        assert_eq!(out[1], vec![2]);
        assert_eq!(out[2], vec![3]);
    }

    #[test]
    #[should_panic(expected = "Sequential stage count")]
    fn test_assign_sequential_ownership_panics_when_sequential_exceeds_workers() {
        let kinds = vec![Parallelism::Sequential, Parallelism::Sequential, Parallelism::Sequential];
        let _ = assign_sequential_ownership(2, &kinds);
    }
}
