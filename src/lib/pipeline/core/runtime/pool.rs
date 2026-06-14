//! Exclusive owner assignment.

use crate::pipeline::core::erased::ErasedStep;
use crate::pipeline::core::signal::PipelineError;
use crate::pipeline::core::step::StepKind;
use crate::pipeline::core::topology::StepIdx;

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
    let total_exclusive = steps.iter().filter(|s| s.profile().kind == StepKind::Exclusive).count();
    if total_exclusive > n_threads {
        return Err(PipelineError::NotEnoughThreads {
            required: total_exclusive,
            available: n_threads,
        });
    }

    let mut owners: Vec<Option<usize>> = vec![None; steps.len()];
    let mut next_owner = 0usize;
    for (idx, step) in steps.iter().enumerate() {
        if step.profile().kind == StepKind::Exclusive {
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
/// If two steps' sticky-ownership rules collide on the same worker (e.g.,
/// the worker is both an Exclusive-sticky owner AND a Serial-sticky-Affinity
/// target), the exclusive ownership wins — the framework only models one
/// sticky-owned step per worker today. Returns `None` for workers without
/// any sticky-owned step.
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
        let profile = step.profile();
        if profile.kind == StepKind::Exclusive && profile.sticky {
            if let Some(owner) = exclusive_owners[step_usize] {
                if owner < n_workers && sticky[owner].is_none() {
                    sticky[owner] = Some(StepIdx(step_usize));
                }
            }
        }
    }

    // Second pass: Serial-sticky-Affinity owners (only fill empty slots).
    for (step_usize, step) in steps.iter().enumerate() {
        let profile = step.profile();
        if profile.kind == StepKind::Serial && profile.sticky {
            let target = match step.affinity() {
                crate::pipeline::core::step::Affinity::None => continue,
                crate::pipeline::core::step::Affinity::Reader => 0,
                crate::pipeline::core::step::Affinity::Writer => n_workers.saturating_sub(1),
                crate::pipeline::core::step::Affinity::Worker(idx) => idx,
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

    use crate::pipeline::core::erased::TypedStep;
    use crate::pipeline::core::outputs::Single;
    use crate::pipeline::core::queues::QueueSpec;
    use crate::pipeline::core::reorder::BranchOrdering;
    use crate::pipeline::core::step::{Step, StepCtx, StepOutcome, StepProfile};

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
}
