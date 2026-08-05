//! `LiveSteps`: per-worker worklist of still-dispatchable steps.
//!
//! Each worker dispatches steps by walking a `Vec<StepIdx>` of the steps it
//! still has work for, in chain order. A step is **removed** from this list
//! when it returns `StepOutcome::Finished`. The worker exits once the list is
//! empty.
//!
//! This replaces the previous scheme of mutating each finished step's
//! `WorkerStepEntry` to `WorkerStepEntry::Skip` in place and re-scanning the
//! full `entries` vec every pass (`entries.iter().all(|e| !e.is_dispatchable())`
//! for the exit check, plus a per-pass `continue` over every inert `Skip`
//! slot). The worklist:
//!
//! - never re-visits a finished or build-time-excluded step,
//! - makes "done" an honest removal rather than an in-place inert variant, and
//! - keeps `StepIdx` stable as the canonical identity into the parallel
//!   `entries` / `contexts.inputs` / `contexts.outputs` / `drain_counters`
//!   arrays (the worklist holds indices *into* that stable storage; it does
//!   not renumber anything).
//!
//! Removal is **stable** (`Vec::remove`, not `swap_remove`): the chain-order
//! invariant the round-robin dispatch relies on (attempt upstream steps before
//! downstream, restart from the front on `Progress`) must survive a removal. At
//! the handful of steps in a chain the linear `remove`/`position` cost is
//! irrelevant.

use crate::runtime::storage::WorkerStepEntry;
use crate::topology::StepIdx;

/// One per worker: the steps this worker can still dispatch, in chain order.
pub struct LiveSteps {
    order: Vec<StepIdx>,
}

impl LiveSteps {
    /// Build from this worker's storage. Every dispatchable entry (i.e. not a
    /// build-time `WorkerStepEntry::Skip` placeholder for an Exclusive step
    /// owned by another worker, or a Serial step this worker's affinity gates
    /// out) enters the worklist, in `StepIdx` (chain) order.
    #[must_use]
    pub fn from_entries(entries: &[WorkerStepEntry]) -> Self {
        let order = entries
            .iter()
            .enumerate()
            .filter(|(_, entry)| entry.is_dispatchable())
            .map(|(idx, _)| StepIdx(idx))
            .collect();
        Self { order }
    }

    /// The steps this worker can still dispatch, in chain order.
    #[must_use]
    pub fn order(&self) -> &[StepIdx] {
        &self.order
    }

    /// Number of still-dispatchable steps.
    #[must_use]
    pub fn len(&self) -> usize {
        self.order.len()
    }

    /// `true` once this worker has nothing left to dispatch (loop exit).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    /// `true` if `step_idx` is still dispatchable by this worker.
    #[must_use]
    pub fn contains(&self, step_idx: StepIdx) -> bool {
        self.order.contains(&step_idx)
    }

    /// Stably remove a finished step. No-op if already absent (idempotent).
    /// Preserves the relative order of the remaining steps.
    pub fn remove(&mut self, step_idx: StepIdx) {
        if let Some(pos) = self.order.iter().position(|&s| s == step_idx) {
            self.order.remove(pos);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use parking_lot::Mutex;

    use crate::erased::{ErasedStep, TypedStep};
    use crate::outputs::Single;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;
    use crate::runtime::storage::DrainGate;
    use crate::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

    #[derive(Clone)]
    struct Nop;
    impl Step for Nop {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Nop",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::CountBounded { capacity: 1 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> std::io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    fn shared_entry() -> WorkerStepEntry {
        let step: Box<dyn ErasedStep> = Box::new(TypedStep::new(Nop));
        WorkerStepEntry::Shared {
            step: Arc::new(Mutex::new(step)),
            drain: Arc::new(DrainGate::default()),
        }
    }

    #[test]
    fn from_entries_excludes_build_time_skips() {
        let entries =
            vec![shared_entry(), WorkerStepEntry::Skip, shared_entry(), WorkerStepEntry::Skip];
        let live = LiveSteps::from_entries(&entries);
        assert_eq!(live.order(), &[StepIdx(0), StepIdx(2)]);
        assert!(!live.is_empty());
        assert_eq!(live.len(), 2);
    }

    #[test]
    fn remove_is_stable_and_idempotent() {
        let entries = vec![shared_entry(), shared_entry(), shared_entry()];
        let mut live = LiveSteps::from_entries(&entries);
        assert_eq!(live.order(), &[StepIdx(0), StepIdx(1), StepIdx(2)]);

        // Remove the middle; the survivors keep their relative order.
        live.remove(StepIdx(1));
        assert_eq!(live.order(), &[StepIdx(0), StepIdx(2)]);
        assert!(!live.contains(StepIdx(1)));

        // Removing an already-absent step is a no-op.
        live.remove(StepIdx(1));
        assert_eq!(live.order(), &[StepIdx(0), StepIdx(2)]);

        live.remove(StepIdx(0));
        live.remove(StepIdx(2));
        assert!(live.is_empty());
    }
}
