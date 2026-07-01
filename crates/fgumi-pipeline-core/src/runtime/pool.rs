//! Exclusive owner assignment.

use crate::erased::ErasedStep;
use crate::signal::PipelineError;
use crate::step::StepKind;
use crate::topology::StepIdx;

/// Assign Exclusive steps to specific worker threads in chain declaration
/// order. Returns `Ok(owners)` where `owners[step_idx] == Some(worker_id)`
/// for Exclusive steps and `None` otherwise.
///
/// # Errors
///
/// Returns `PipelineError::NotEnoughThreads` if more Exclusive steps exist
/// than worker threads available.
pub fn assign_exclusive_owners(
    steps: &[Box<dyn ErasedStep>],
    n_threads: usize,
) -> Result<Vec<Option<usize>>, PipelineError> {
    let total_exclusive = steps.iter().filter(|s| s.kind() == StepKind::Exclusive).count();
    if total_exclusive > n_threads {
        return Err(PipelineError::NotEnoughThreads {
            required: total_exclusive,
            available: n_threads,
        });
    }

    let mut owners: Vec<Option<usize>> = vec![None; steps.len()];
    let mut next_owner = 0usize;
    for (idx, step) in steps.iter().enumerate() {
        if step.kind() == StepKind::Exclusive {
            owners[idx] = Some(next_owner);
            next_owner += 1;
        }
    }
    Ok(owners)
}

/// Compute each worker's sticky-driven step (the step its `WorkerCore::sticky_owner`
/// will hold). A step contributes to a worker's `sticky_owner` iff it's flagged
/// `sticky=true` in its profile AND that worker is its sole eligible dispatcher:
///   - `Exclusive sticky` step → sticky owner is the worker assigned by
///     `assign_exclusive_owners` (read from `exclusive_owners[step_idx]`).
///   - `Serial sticky` + `Affinity::Reader` → sticky owner is worker 0.
///   - `Serial sticky` + `Affinity::Writer` → sticky owner is worker `N-1`.
///   - `Serial sticky` + `Affinity::Worker(idx)` → sticky owner is worker `idx`.
///   - `Serial sticky` + `Affinity::None` → no sticky owner (every worker
///     is eligible, so no single worker can drive sticky without starving
///     the others; we silently drop the sticky hint here).
///   - `Parallel` steps are never sticky-driven (each worker has its own
///     clone — sticky drive on one would not gate others).
///
/// If two steps' sticky-ownership rules collide on the same worker, the first
/// writer to that worker's slot wins — the framework only models one
/// sticky-owned step per worker today. The Exclusive pass runs before the Serial
/// pass and each pass only fills empty slots, so an Exclusive owner beats a
/// later Serial-sticky-Affinity target on the same worker; the same first-wins
/// rule resolves Exclusive-vs-Exclusive and Serial-vs-Serial collisions too.
/// Returns `None` for workers without any sticky-owned step.
#[must_use]
pub fn assign_sticky_owners(
    steps: &[Box<dyn ErasedStep>],
    exclusive_owners: &[Option<usize>],
    n_workers: usize,
) -> Vec<Option<StepIdx>> {
    debug_assert_eq!(steps.len(), exclusive_owners.len());
    let mut sticky: Vec<Option<StepIdx>> = vec![None; n_workers];

    // First pass: Exclusive-sticky owners (highest priority).
    for (step_usize, step) in steps.iter().enumerate() {
        if step.kind() == StepKind::Exclusive && step.sticky() {
            if let Some(owner) = exclusive_owners[step_usize] {
                if owner < n_workers && sticky[owner].is_none() {
                    sticky[owner] = Some(StepIdx(step_usize));
                }
            }
        }
    }

    // Second pass: Serial-sticky-Affinity owners (only fill empty slots).
    for (step_usize, step) in steps.iter().enumerate() {
        if step.kind() == StepKind::Serial && step.sticky() {
            let target = match step.affinity() {
                crate::step::Affinity::None => continue,
                crate::step::Affinity::Reader => 0,
                crate::step::Affinity::Writer => n_workers.saturating_sub(1),
                crate::step::Affinity::Worker(idx) => idx,
            };
            if target < n_workers && sticky[target].is_none() {
                sticky[target] = Some(StepIdx(step_usize));
            }
        }
    }

    sticky
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    use rstest::rstest;

    use crate::erased::TypedStep;
    use crate::outputs::Single;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;
    use crate::step::{Affinity, Step, StepCtx, StepOutcome, StepProfile};

    fn stub_step(kind: StepKind) -> Box<dyn ErasedStep> {
        #[derive(Clone)]
        struct StubStep(StepKind);
        impl Step for StubStep {
            type Input = u32;
            type Outputs = Single<u32>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "Stub",
                    kind: self.0,
                    sticky: false,
                    output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                    branch_ordering: vec![BranchOrdering::None],
                }
            }
            fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                Ok(StepOutcome::NoProgress)
            }
            fn new_worker_copy(&self) -> Self {
                self.clone()
            }
        }
        Box::new(TypedStep::new(StubStep(kind)))
    }

    #[test]
    fn no_exclusives_returns_all_none() {
        let steps = vec![stub_step(StepKind::Parallel), stub_step(StepKind::Serial)];
        let owners = assign_exclusive_owners(&steps, 4).unwrap();
        assert_eq!(owners, vec![None, None]);
    }

    #[test]
    fn one_exclusive_assigned_to_thread_zero() {
        let steps = vec![stub_step(StepKind::Exclusive), stub_step(StepKind::Parallel)];
        let owners = assign_exclusive_owners(&steps, 4).unwrap();
        assert_eq!(owners, vec![Some(0), None]);
    }

    #[test]
    fn two_exclusives_assigned_zero_and_one() {
        let steps = vec![
            stub_step(StepKind::Exclusive),
            stub_step(StepKind::Parallel),
            stub_step(StepKind::Exclusive),
        ];
        let owners = assign_exclusive_owners(&steps, 4).unwrap();
        assert_eq!(owners, vec![Some(0), None, Some(1)]);
    }

    #[test]
    fn too_many_exclusives_returns_not_enough_threads() {
        let steps = vec![
            stub_step(StepKind::Exclusive),
            stub_step(StepKind::Exclusive),
            stub_step(StepKind::Exclusive),
        ];
        let result = assign_exclusive_owners(&steps, 2);
        assert!(matches!(
            result,
            Err(PipelineError::NotEnoughThreads { required: 3, available: 2 })
        ));
    }

    fn sticky_exclusive_step(owner_idx: usize) -> Box<dyn ErasedStep> {
        #[derive(Clone)]
        struct StickyExclusive;
        impl Step for StickyExclusive {
            type Input = u32;
            type Outputs = Single<u32>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "StickyExclusive",
                    kind: StepKind::Exclusive,
                    sticky: true,
                    output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                    branch_ordering: vec![BranchOrdering::None],
                }
            }
            fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                Ok(StepOutcome::NoProgress)
            }
            fn new_worker_copy(&self) -> Self {
                self.clone()
            }
        }
        let _ = owner_idx; // owner_idx used by caller, not embedded in step
        Box::new(TypedStep::new(StickyExclusive))
    }

    fn sticky_serial_step(affinity: crate::step::Affinity) -> Box<dyn ErasedStep> {
        struct StickySerial(crate::step::Affinity);
        impl Step for StickySerial {
            type Input = u32;
            type Outputs = Single<u32>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "StickySerial",
                    kind: StepKind::Serial,
                    sticky: true,
                    output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                    branch_ordering: vec![BranchOrdering::None],
                }
            }
            fn affinity(&self) -> crate::step::Affinity {
                self.0
            }
            fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                Ok(StepOutcome::NoProgress)
            }
        }
        Box::new(TypedStep::new(StickySerial(affinity)))
    }

    /// Assert that `sticky` holds `Some(StepIdx(0))` in exactly `expected_slot`
    /// (when `Some`) and `None` everywhere else across `n_workers` slots.
    fn assert_only_slot(
        sticky: &[Option<StepIdx>],
        n_workers: usize,
        expected_slot: Option<usize>,
    ) {
        for (slot, &got) in sticky.iter().enumerate().take(n_workers) {
            let want = if Some(slot) == expected_slot { Some(StepIdx(0)) } else { None };
            assert_eq!(got, want, "slot {slot}; expected owner slot {expected_slot:?}");
        }
    }

    #[rstest]
    #[case::in_range(2, Some(2))]
    #[case::out_of_range(10, None)]
    fn sticky_exclusive_owner_maps_to_slot(
        #[case] owner: usize,
        #[case] expected_slot: Option<usize>,
    ) {
        // A sticky-exclusive step's slot is its (in-range) owner worker; an
        // out-of-range owner is skipped, leaving every slot empty.
        let steps = vec![sticky_exclusive_step(0), stub_step(StepKind::Parallel)];
        let exclusive_owners = vec![Some(owner), None];
        let sticky = assign_sticky_owners(&steps, &exclusive_owners, 4);
        assert_only_slot(&sticky, 4, expected_slot);
    }

    #[test]
    fn sticky_exclusive_occupied_slot_not_overwritten() {
        // Two sticky-exclusive steps competing for the same worker slot.
        let steps = vec![sticky_exclusive_step(0), sticky_exclusive_step(0)];
        let exclusive_owners = vec![Some(0_usize), Some(0_usize)];
        let sticky = assign_sticky_owners(&steps, &exclusive_owners, 4);
        // First step wins; second is skipped because slot[0] is already occupied.
        assert_eq!(sticky[0], Some(StepIdx(0)));
    }

    #[rstest]
    #[case::reader(Affinity::Reader, Some(0))]
    #[case::writer(Affinity::Writer, Some(3))]
    #[case::worker_in_range(Affinity::Worker(2), Some(2))]
    #[case::worker_out_of_range(Affinity::Worker(10), None)]
    #[case::none(Affinity::None, None)]
    fn sticky_serial_affinity_maps_to_slot(
        #[case] affinity: Affinity,
        #[case] expected_slot: Option<usize>,
    ) {
        // Serial affinity resolves to a single worker slot: Reader→0,
        // Writer→last, Worker(i)→i; an out-of-range Worker index and None are
        // skipped, leaving every slot empty.
        let steps = vec![sticky_serial_step(affinity)];
        let exclusive_owners = vec![None];
        let sticky = assign_sticky_owners(&steps, &exclusive_owners, 4);
        assert_only_slot(&sticky, 4, expected_slot);
    }

    #[test]
    fn sticky_exclusive_beats_sticky_serial_on_same_slot() {
        // Exclusive pass runs first; serial pass only fills empty slots.
        let exc = sticky_exclusive_step(0);
        let ser = sticky_serial_step(crate::step::Affinity::Reader); // also targets slot 0
        let steps: Vec<Box<dyn ErasedStep>> = vec![exc, ser];
        let exclusive_owners = vec![Some(0_usize), None];
        let sticky = assign_sticky_owners(&steps, &exclusive_owners, 4);
        // Exclusive (step 0) wins slot 0; Serial (step 1) is blocked.
        assert_eq!(sticky[0], Some(StepIdx(0)));
    }
}
