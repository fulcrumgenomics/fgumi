//! `WorkerStepEntry`: per-worker per-step storage shape. Determined by
//! `StepKind` at run start; immutable thereafter.

use parking_lot::Mutex;
use std::sync::{Arc, OnceLock};

use crate::pipeline::core::erased::ErasedStep;
use crate::pipeline::core::step::{Affinity, StepKind};

/// Shared "this step has finished" latch for a `Serial` step.
///
/// A `Serial` step is one shared instance behind a `Mutex` across all workers,
/// but each worker tracks its own worklist. Without a shared latch, after one
/// worker runs the step to `Finished` the others still hold it in their
/// worklist and would re-`try_lock`/re-run an already-finished step. The
/// finishing worker sets this latch under the dispatch guard (before
/// `mark_outputs_drained`); every other worker observes it and short-circuits
/// to a synthetic `Finished` instead of re-running the step.
#[derive(Default)]
pub struct DrainGate {
    finished: OnceLock<()>,
}

impl DrainGate {
    /// True once the step has run to `StepOutcome::Finished` on some worker.
    pub fn is_finished(&self) -> bool {
        self.finished.get().is_some()
    }

    /// Latch the step as finished (idempotent). Set under the dispatch guard,
    /// before `mark_outputs_drained`, so a concurrent worker that observes the
    /// latch never re-runs the step.
    pub fn mark_finished(&self) {
        let _ = self.finished.set(());
    }
}

/// One per (worker, step) cell.
pub enum WorkerStepEntry {
    /// `Parallel` step: this worker owns its private clone via `clone_boxed`.
    /// Direct `&mut` access; no locking on `try_run` dispatch. Completion is
    /// coordinated by the per-step `StepDrainCounter` (init N) in the driver:
    /// every clone returns `Finished` when the input drains, but only the last
    /// to finish closes the shared output.
    Owned { step: Box<dyn ErasedStep> },
    /// `Serial` step: shared instance, mutex-protected. Any worker can acquire.
    /// The shared `DrainGate` finished-latch lets a worker that finishes the
    /// step stop the others from re-running it.
    Shared { step: Arc<Mutex<Box<dyn ErasedStep>>>, drain: Arc<DrainGate> },
    /// `Exclusive` step: only this worker (the owner) ever runs it.
    /// Stored locally on the owner; other workers have `Skip`.
    Exclusive { step: Box<dyn ErasedStep> },
    /// This worker doesn't run this step (it's an `Exclusive` step owned by
    /// another worker). The worker loop skips it in dispatch.
    Skip,
}

impl WorkerStepEntry {
    /// True if this worker should attempt to dispatch this step.
    #[must_use]
    pub fn is_dispatchable(&self) -> bool {
        !matches!(self, Self::Skip)
    }
}

/// Assemble per-worker step storage.
///
/// `steps`: the chain in order (consumed).
/// `exclusive_owners[step_idx] == Some(worker_id)` for each Exclusive step.
/// `n_workers`: number of worker threads.
///
/// Returns `entries[worker_id][step_idx] = WorkerStepEntry`.
///
/// # Panics
///
/// Panics if `exclusive_owners.len() != steps.len()` or if an Exclusive
/// step has no owner assignment (`assign_exclusive_owners` must run first).
#[must_use]
pub fn build_worker_storage(
    steps: Vec<Box<dyn ErasedStep>>,
    exclusive_owners: &[Option<usize>],
    n_workers: usize,
) -> Vec<Vec<WorkerStepEntry>> {
    assert_eq!(
        steps.len(),
        exclusive_owners.len(),
        "exclusive_owners length must match step count"
    );
    assert!(n_workers > 0, "build_worker_storage requires at least one worker");

    let mut entries: Vec<Vec<WorkerStepEntry>> =
        (0..n_workers).map(|_| Vec::with_capacity(steps.len())).collect();

    for (step_idx, step) in steps.into_iter().enumerate() {
        let kind = step.profile().kind;
        match kind {
            StepKind::Parallel => {
                // Each worker gets its own clone; the original goes to worker N-1.
                for entries_for_worker in entries.iter_mut().take(n_workers - 1) {
                    entries_for_worker.push(WorkerStepEntry::Owned { step: step.clone_boxed() });
                }
                entries[n_workers - 1].push(WorkerStepEntry::Owned { step });
            }
            StepKind::Serial => {
                // Snapshot the affinity hint before moving `step` into
                // the shared `Arc<Mutex<...>>` — affinity is part of the
                // static `Step` description, so calling it here is a
                // single virtual dispatch.
                let affinity = step.affinity();
                // Always-on (not `debug_assert!`): an out-of-range affinity
                // gates the Serial step out of every worker, so it is never
                // dispatched and the pipeline deadlocks. The check is a cheap
                // pure predicate, so it stays enabled in release builds too.
                assert!(
                    affinity_in_range(affinity, n_workers),
                    "Serial step affinity {affinity:?} requests a worker out of range \
                     (n_workers = {n_workers}); pipeline would deadlock"
                );
                let drain = Arc::new(DrainGate::default());
                let shared = Arc::new(Mutex::new(step));
                for (worker_id, entries_for_worker) in entries.iter_mut().enumerate() {
                    if affinity.eligible(worker_id, n_workers) {
                        entries_for_worker.push(WorkerStepEntry::Shared {
                            step: Arc::clone(&shared),
                            drain: Arc::clone(&drain),
                        });
                    } else {
                        // This worker is gated out by the Serial step's
                        // affinity hint. It will never `try_lock` the step's
                        // mutex, eliminating the `Contention` thrash that
                        // pure `Serial` exhibits at high thread counts.
                        entries_for_worker.push(WorkerStepEntry::Skip);
                    }
                }
            }
            StepKind::Exclusive => {
                let owner = exclusive_owners[step_idx].expect(
                    "Exclusive step has no owner assignment; \
                     assign_exclusive_owners must run before build_worker_storage",
                );
                // Push Skip placeholders for all workers, then replace owner's slot.
                for entries_for_worker in &mut entries {
                    entries_for_worker.push(WorkerStepEntry::Skip);
                }
                entries[owner][step_idx] = WorkerStepEntry::Exclusive { step };
            }
        }
    }

    entries
}

/// Validate that a Serial step's affinity refers to a worker that exists.
/// `Affinity::None` / `Reader` / `Writer` always resolve to a real worker
/// when `n_workers >= 1` (which the runtime's `assert!(n_workers > 0)`
/// already guarantees). `Worker(idx)` is only valid when `idx < n_workers`.
#[inline]
fn affinity_in_range(affinity: Affinity, n_workers: usize) -> bool {
    match affinity {
        Affinity::None | Affinity::Reader | Affinity::Writer => true,
        Affinity::Worker(idx) => idx < n_workers,
    }
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

    fn profile_for(name: &'static str, kind: StepKind, sticky: bool) -> StepProfile {
        StepProfile {
            name,
            kind,
            sticky,
            output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    #[derive(Clone)]
    struct ParallelStep;
    impl Step for ParallelStep {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            profile_for("Par", StepKind::Parallel, false)
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    #[derive(Clone)]
    struct SerialStep;
    impl Step for SerialStep {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            profile_for("Ser", StepKind::Serial, false)
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    #[derive(Clone)]
    struct ExclusiveStep;
    impl Step for ExclusiveStep {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            profile_for("Excl", StepKind::Exclusive, false)
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    #[test]
    fn skip_is_not_dispatchable() {
        assert!(!WorkerStepEntry::Skip.is_dispatchable());
    }

    #[test]
    fn parallel_step_yields_owned_per_worker() {
        let steps: Vec<Box<dyn ErasedStep>> = vec![Box::new(TypedStep::new(ParallelStep))];
        let owners = vec![None];
        let entries = build_worker_storage(steps, &owners, 4);
        assert_eq!(entries.len(), 4);
        for w in &entries {
            assert_eq!(w.len(), 1);
            assert!(matches!(w[0], WorkerStepEntry::Owned { .. }));
        }
    }

    #[test]
    fn serial_step_yields_shared_arc_for_all_workers() {
        let steps: Vec<Box<dyn ErasedStep>> = vec![Box::new(TypedStep::new(SerialStep))];
        let owners = vec![None];
        let entries = build_worker_storage(steps, &owners, 3);
        for w in &entries {
            assert!(matches!(w[0], WorkerStepEntry::Shared { .. }));
        }
        let arc0 = if let WorkerStepEntry::Shared { step, .. } = &entries[0][0] {
            Arc::clone(step)
        } else {
            panic!()
        };
        let arc1 = if let WorkerStepEntry::Shared { step, .. } = &entries[1][0] {
            Arc::clone(step)
        } else {
            panic!()
        };
        assert!(Arc::ptr_eq(&arc0, &arc1));
    }

    #[test]
    fn exclusive_step_owner_gets_exclusive_others_skip() {
        let steps: Vec<Box<dyn ErasedStep>> = vec![Box::new(TypedStep::new(ExclusiveStep))];
        let owners = vec![Some(2)];
        let entries = build_worker_storage(steps, &owners, 4);
        assert!(matches!(entries[0][0], WorkerStepEntry::Skip));
        assert!(matches!(entries[1][0], WorkerStepEntry::Skip));
        assert!(matches!(entries[2][0], WorkerStepEntry::Exclusive { .. }));
        assert!(matches!(entries[3][0], WorkerStepEntry::Skip));
    }

    #[test]
    fn mixed_chain_assigns_correctly() {
        let steps: Vec<Box<dyn ErasedStep>> = vec![
            Box::new(TypedStep::new(ExclusiveStep)),
            Box::new(TypedStep::new(ParallelStep)),
            Box::new(TypedStep::new(SerialStep)),
            Box::new(TypedStep::new(ExclusiveStep)),
        ];
        let owners = vec![Some(0), None, None, Some(1)];
        let entries = build_worker_storage(steps, &owners, 4);

        assert!(matches!(entries[0][0], WorkerStepEntry::Exclusive { .. }));
        assert!(matches!(entries[1][0], WorkerStepEntry::Skip));

        for w in &entries {
            assert!(matches!(w[1], WorkerStepEntry::Owned { .. }));
        }

        for w in &entries {
            assert!(matches!(w[2], WorkerStepEntry::Shared { .. }));
        }

        assert!(matches!(entries[1][3], WorkerStepEntry::Exclusive { .. }));
        assert!(matches!(entries[0][3], WorkerStepEntry::Skip));
    }
}
