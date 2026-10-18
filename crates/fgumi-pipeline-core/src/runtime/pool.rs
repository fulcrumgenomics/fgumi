//! Exclusive owner assignment.

use std::collections::BTreeSet;

use parking_lot::Mutex;

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
///   - `Parallel` steps are never sticky-driven (each host worker has its own
///     clone — sticky drive on one would not gate others).
///
/// If two steps' sticky-ownership rules collide on the same worker, the first
/// writer to that worker's slot wins — the framework only models one
/// sticky-owned step per worker today. The Exclusive pass runs before the Serial
/// pass and each pass only fills empty slots, so an Exclusive owner beats a
/// later Serial-sticky-Affinity target on the same worker; the same first-wins
/// rule resolves Exclusive-vs-Exclusive and Serial-vs-Serial collisions too.
///
/// Each losing request is logged at debug level, naming the worker, the step
/// that kept the slot and the step that lost it, once per such triple per
/// process: the loser runs without the sticky fast path, a throughput change
/// that is otherwise invisible. Only a step that declares `sticky: true`
/// reaches this function's collision check, so a step that opts out at its
/// profile is never logged, even when it shares a worker with a sticky step
/// (the sort's `FrameBgzfBlocks` declares `sticky: false` because the planner
/// holds the reader worker's slot, and batches its pushes instead).
///
/// Returns `None` for workers without any sticky-owned step.
#[must_use]
pub fn assign_sticky_owners(
    steps: &[Box<dyn ErasedStep>],
    exclusive_owners: &[Option<usize>],
    n_workers: usize,
) -> Vec<Option<StepIdx>> {
    debug_assert_eq!(steps.len(), exclusive_owners.len());
    let mut sticky: Vec<Option<StepIdx>> = vec![None; n_workers];
    // Claim `worker`'s slot for `step_usize`; an occupied slot keeps its owner
    // and the dropped request is reported.
    let mut claim = |worker: usize, step_usize: usize| match sticky[worker] {
        None => sticky[worker] = Some(StepIdx(step_usize)),
        Some(StepIdx(kept)) => {
            log_dropped_sticky(worker, steps[kept].name(), steps[step_usize].name());
        }
    };

    // First pass: Exclusive-sticky owners (highest priority).
    for (step_usize, step) in steps.iter().enumerate() {
        if step.kind() == StepKind::Exclusive
            && step.sticky()
            && let Some(owner) = exclusive_owners[step_usize]
            && owner < n_workers
        {
            claim(owner, step_usize);
        }
    }

    // Second pass: Serial-sticky-Affinity owners (only fill empty slots).
    for (step_usize, step) in steps.iter().enumerate() {
        // Resolve through the shared `Affinity::target_worker` helper rather
        // than re-matching the variants here: `Affinity::eligible` gates which
        // worker may dispatch the step, and a local copy of that mapping could
        // drift and hand a sticky owner to a worker that `Skip`s the step.
        if step.kind() == StepKind::Serial
            && step.sticky()
            && let Some(target) = step.affinity().target_worker(n_workers)
            && target < n_workers
        {
            claim(target, step_usize);
        }
    }

    sticky
}

/// Log (at debug level) that `dropped`'s sticky request lost `worker`'s slot
/// to `kept`, once per `(worker, kept, dropped)` per process: a chain is built
/// once per run, but a command may build several with the same steps, and a
/// repeated line adds nothing. The worker is part of the key, so the same pair
/// colliding on another worker (a pool of another size puts a writer-affine
/// step on another last worker) is logged with its own worker index.
///
/// Debug, not warn: collisions are routine in production chains (two writer
/// sinks, e.g. `filter --rejects`, both target the last worker; a reader and a
/// writer share the only worker of a one-worker pool), and the user can do
/// nothing about any of them. The line is for whoever is diagnosing a step's
/// dispatch cadence.
fn log_dropped_sticky(worker: usize, kept: &'static str, dropped: &'static str) {
    static LOGGED: Mutex<BTreeSet<(usize, &'static str, &'static str)>> =
        Mutex::new(BTreeSet::new());
    let first = LOGGED.lock().insert((worker, kept, dropped));
    if first {
        log::debug!(
            "Sticky dispatch for step `{dropped}` was dropped: worker {worker} already drives \
             sticky step `{kept}`, and a worker drives at most one sticky step. `{dropped}` \
             runs without the sticky fast path."
        );
    }
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

    fn sticky_exclusive_step() -> Box<dyn ErasedStep> {
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
        Box::new(TypedStep::new(StickyExclusive))
    }

    fn sticky_serial_step(affinity: crate::step::Affinity) -> Box<dyn ErasedStep> {
        named_sticky_serial_step("StickySerial", affinity)
    }

    fn named_sticky_serial_step(
        name: &'static str,
        affinity: crate::step::Affinity,
    ) -> Box<dyn ErasedStep> {
        struct StickySerial(&'static str, crate::step::Affinity);
        impl Step for StickySerial {
            type Input = u32;
            type Outputs = Single<u32>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: self.0,
                    kind: StepKind::Serial,
                    sticky: true,
                    output_queues: vec![QueueSpec::CountBounded { capacity: 4 }],
                    branch_ordering: vec![BranchOrdering::None],
                }
            }
            fn affinity(&self) -> crate::step::Affinity {
                self.1
            }
            fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                Ok(StepOutcome::NoProgress)
            }
        }
        Box::new(TypedStep::new(StickySerial(name, affinity)))
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
        let steps = vec![sticky_exclusive_step(), stub_step(StepKind::Parallel)];
        let exclusive_owners = vec![Some(owner), None];
        let sticky = assign_sticky_owners(&steps, &exclusive_owners, 4);
        assert_only_slot(&sticky, 4, expected_slot);
    }

    #[test]
    fn sticky_exclusive_occupied_slot_not_overwritten() {
        // Two sticky-exclusive steps competing for the same worker slot.
        let steps = vec![sticky_exclusive_step(), sticky_exclusive_step()];
        let exclusive_owners = vec![Some(0_usize), Some(0_usize)];
        let sticky = assign_sticky_owners(&steps, &exclusive_owners, 4);
        // First step wins; second is skipped because slot[0] is already occupied.
        // The loser must not land in any other slot either.
        assert_only_slot(&sticky, 4, Some(0));
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
        let exc = sticky_exclusive_step();
        let ser = sticky_serial_step(crate::step::Affinity::Reader); // also targets slot 0
        let steps: Vec<Box<dyn ErasedStep>> = vec![exc, ser];
        let exclusive_owners = vec![Some(0_usize), None];
        let sticky = assign_sticky_owners(&steps, &exclusive_owners, 4);
        // Exclusive (step 0) wins slot 0; Serial (step 1) is blocked — and the
        // blocked serial step must not land in any other slot either.
        assert_only_slot(&sticky, 4, Some(0));
    }

    /// Captured records from this module that mention `needle`.
    fn records_mentioning(needle: &str) -> Vec<crate::test_support::CapturedRecord> {
        crate::test_support::captured_logs(module_path!().trim_end_matches("::tests"), needle)
    }

    #[test]
    fn dropped_sticky_request_is_logged_once_per_worker_and_pair() {
        // A synthetic collision: three sticky steps target the reader worker,
        // and the first keeps its slot.
        crate::test_support::install_log_capture();
        let reader_steps = vec![
            named_sticky_serial_step("DropLogKept", Affinity::Reader),
            stub_step(StepKind::Parallel),
            named_sticky_serial_step("DropLogFirst", Affinity::Reader),
            named_sticky_serial_step("DropLogSecond", Affinity::Reader),
        ];
        let sticky = assign_sticky_owners(&reader_steps, &[None, None, None, None], 4);
        assert_only_slot(&sticky, 4, Some(0));
        // A second chain with the same collisions (a command may build several)
        // repeats nothing.
        let _ = assign_sticky_owners(&reader_steps, &[None, None, None, None], 4);

        // Two writer-affine sticky steps collide on the last worker: worker 3
        // of four (twice), then worker 1 of two, the same pair on a new worker.
        let writer_steps = vec![
            named_sticky_serial_step("DropLogSinkA", Affinity::Writer),
            named_sticky_serial_step("DropLogSinkB", Affinity::Writer),
        ];
        let _ = assign_sticky_owners(&writer_steps, &[None, None], 4);
        let _ = assign_sticky_owners(&writer_steps, &[None, None], 4);
        let _ = assign_sticky_owners(&writer_steps, &[None, None], 2);

        let records = records_mentioning("DropLog");
        let expected = [
            (0, "DropLogKept", "DropLogFirst"),
            (0, "DropLogKept", "DropLogSecond"),
            (3, "DropLogSinkA", "DropLogSinkB"),
            (1, "DropLogSinkA", "DropLogSinkB"),
        ];
        assert_eq!(
            records.len(),
            expected.len(),
            "one record per (worker, kept, dropped): {records:?}"
        );
        for (worker, kept, dropped) in expected {
            let hits: Vec<_> = records
                .iter()
                .filter(|r| {
                    r.message.contains(&format!("step `{dropped}` was dropped"))
                        && r.message.contains(&format!("worker {worker} already"))
                        && r.message.contains(&format!("sticky step `{kept}`"))
                })
                .collect();
            assert_eq!(hits.len(), 1, "({worker}, {kept}, {dropped}): {records:?}");
            assert_eq!(hits[0].level, log::Level::Debug, "routine in production chains, so debug");
        }
    }

    #[test]
    fn sticky_requests_on_distinct_workers_are_not_logged() {
        crate::test_support::install_log_capture();
        let steps = vec![
            named_sticky_serial_step("NoLogReader", Affinity::Reader),
            named_sticky_serial_step("NoLogWriter", Affinity::Writer),
        ];
        let sticky = assign_sticky_owners(&steps, &[None, None], 4);
        assert_eq!(sticky, vec![Some(StepIdx(0)), None, None, Some(StepIdx(1))]);
        assert!(records_mentioning("NoLog").is_empty());
    }
}
