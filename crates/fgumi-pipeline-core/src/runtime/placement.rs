//! Where each `Parallel` step's clones run. One function
//! ([`plan_parallel_hosts`]) decides it; the drain-counter init
//! ([`drain_counter_inits`]), worker storage, driver-group membership, the wake
//! plan ([`wake_driver_of`] and the fallback masks) and `dag_at` all read its
//! result, so the driver's "counter init == clone count" invariant holds by
//! construction.

use crate::erased::ErasedStep;
use crate::step::{Affinity, PoolPlacement, StepKind};
use crate::topology::{ChainGraph, StepIdx};

/// Index of a dedicated driver thread (one per `DetachedDriverGroup`, in the
/// order `extract_detached_steps` returns them). Defined here, where placement
/// names the driver hosting a clone; the wake plan re-exports it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DriverIdx(pub usize);

/// Where a `Parallel` step's clones live: pool workers by id and/or one driver.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ParallelHosts {
    /// Pool workers holding a clone, ascending.
    pub workers: Vec<usize>,
    /// A driver hosting the single clone (only when `workers` is empty).
    pub driver: Option<DriverIdx>,
}

impl ParallelHosts {
    /// Total clones — the step's `StepDrainCounter` init.
    #[must_use]
    pub fn clone_count(&self) -> usize {
        self.workers.len() + usize::from(self.driver.is_some())
    }
}

/// The pool worker hosting an `Affinity::Reader` `Serial` step, if the chain has
/// one. Only `Reader` counts: a `Writer` or `Worker(k)` pin is not the read
/// path the hint keeps free.
#[must_use]
pub fn reader_worker(steps: &[Box<dyn ErasedStep>], n_workers: usize) -> Option<usize> {
    steps
        .iter()
        .filter(|s| s.kind() == StepKind::Serial && s.affinity() == Affinity::Reader)
        .find_map(|s| s.affinity().target_worker(n_workers))
        .filter(|&w| w < n_workers)
}

/// The placement rule for one `Parallel` step (see [`PoolPlacement`]).
///
/// `ExcludeReader` with a reader drops the reader's worker; if that leaves no
/// worker, the single clone goes to `producer_driver` when the step's producer
/// runs on a driver, else to worker 0. Every other combination is every worker.
/// For a multi-input step, "the producer" is the one the pipeline's host planner
/// passes: its lowest-indexed producer.
#[must_use]
pub fn parallel_hosts(
    placement: PoolPlacement,
    reader: Option<usize>,
    n_workers: usize,
    producer_driver: Option<DriverIdx>,
) -> ParallelHosts {
    match (placement, reader) {
        (PoolPlacement::ExcludeReader, Some(r)) => {
            let workers: Vec<usize> = (0..n_workers).filter(|&w| w != r).collect();
            if !workers.is_empty() {
                ParallelHosts { workers, driver: None }
            } else if let Some(d) = producer_driver {
                ParallelHosts { workers: Vec::new(), driver: Some(d) }
            } else {
                ParallelHosts { workers: vec![0], driver: None }
            }
        }
        _ => ParallelHosts { workers: (0..n_workers).collect(), driver: None },
    }
}

/// Hosts for every step, in `StepIdx` order. `driver_of[i]` is the driver of a
/// Detached step `i` (`driver_index_of`). Non-Parallel steps get
/// `ParallelHosts::default()`.
///
/// A step's producer is [`ChainGraph::first_producer_into`]: the
/// lowest-indexed step wired into it. A single-input step has exactly one; for
/// a multi-input (`Step2` / `StepK`) step only that one decides whether its
/// clone can be hosted on a driver — another Detached producer of the same
/// step does not.
///
/// Called only by `Pipeline::run` and `Pipeline::dag_at`; every consumer of
/// placement reads this one result rather than re-deriving it.
///
/// # Panics
///
/// Panics if `driver_of.len() != steps.len()`.
#[must_use]
pub fn plan_parallel_hosts(
    steps: &[Box<dyn ErasedStep>],
    graph: &ChainGraph,
    driver_of: &[Option<DriverIdx>],
    n_workers: usize,
) -> Vec<ParallelHosts> {
    assert_eq!(driver_of.len(), steps.len(), "driver_of must cover every step");
    let reader = reader_worker(steps, n_workers);
    steps
        .iter()
        .enumerate()
        .map(|(i, s)| {
            if s.kind() != StepKind::Parallel {
                return ParallelHosts::default();
            }
            let producer_driver =
                graph.first_producer_into(StepIdx(i)).and_then(|p| driver_of[p.0]);
            parallel_hosts(s.pool_placement(), reader, n_workers, producer_driver)
        })
        .collect()
}

/// Drain-counter init per step: the Parallel clone count from `hosts`, else 1.
/// (A Parallel *source*'s clones each return `Finished` once their shared
/// counter is exhausted, so init N still matches.) Reads the cached `kind()`,
/// so it is safe over extracted placeholders.
///
/// # Panics
///
/// Panics if `hosts` does not cover every step.
#[must_use]
pub(crate) fn drain_counter_inits(
    steps: &[Box<dyn ErasedStep>],
    hosts: &[ParallelHosts],
) -> Vec<usize> {
    assert_eq!(steps.len(), hosts.len(), "hosts must cover every step");
    steps
        .iter()
        .zip(hosts)
        .map(|(s, h)| match s.kind() {
            StepKind::Parallel => h.clone_count(),
            // Serial/Exclusive have one shared instance, and a `Detached` step
            // runs on a single dedicated thread: exactly one finisher closes its
            // output edge.
            StepKind::Serial | StepKind::Exclusive | StepKind::Detached => 1,
        })
        .collect()
}

/// The `driver_of` the wake plan routes by: every Detached step's driver, plus
/// every hosted Parallel step mapped to the driver hosting it. `run` and
/// `dag_at` both build `detached_driver_of` with `driver_index_of` and call
/// this, so the two renderings cannot disagree.
#[must_use]
pub(crate) fn wake_driver_of(
    detached_driver_of: &[Option<DriverIdx>],
    hosts: &[ParallelHosts],
) -> Vec<Option<DriverIdx>> {
    detached_driver_of.iter().zip(hosts).map(|(d, h)| h.driver.or(*d)).collect()
}

/// Check the assumption the wake plan's event-count path rests on: a worker a
/// `Parallel` step's placement leaves out holds no clone of it, so it must not
/// wait on the shared event-count (`notify_one` could pick it, and the wake
/// would be spent on a thread that cannot run the step). Every worker outside a
/// restricted host set must therefore be pinned — a `Serial` step's affinity
/// target or an `Exclusive` owner (`pinned_worker`, from `wake_inputs`), which
/// idles on its timer. `ExcludeReader` leaves out only the reader's worker,
/// which a `Serial` `Reader` step pins, so this holds; a new placement that
/// drops an unpinned worker trips it at run start.
///
/// # Panics
///
/// Panics if a worker outside some `Parallel` step's host workers is not
/// pinned.
pub(crate) fn assert_excluded_workers_are_pinned(
    hosts: &[ParallelHosts],
    pinned_worker: &[Option<usize>],
    n_workers: usize,
) {
    for (step, h) in hosts.iter().enumerate() {
        if h.clone_count() == 0 || h.workers.len() >= n_workers {
            continue; // not a Parallel step, or on every worker
        }
        for w in (0..n_workers).filter(|w| !h.workers.contains(w)) {
            assert!(
                pinned_worker.contains(&Some(w)),
                "step {step}'s placement leaves out worker {w}, which is not pinned; an \
                 unpinned worker waits on the event-count, where a wake for the step could \
                 reach it"
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io;

    use rstest::rstest;

    use super::*;
    use crate::builder::PipelineBuilder;
    use crate::outputs::Single;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;
    use crate::step::{Step, StepCtx, StepOutcome, StepProfile};

    /// Serial source with a configurable affinity.
    struct AffinitySource(Affinity);
    impl Step for AffinitySource {
        type Input = ();
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "AffinitySource",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![QueueSpec::Unbounded],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn affinity(&self) -> Affinity {
            self.0
        }
        fn try_run(&mut self, _: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::Finished)
        }
    }

    /// Parallel pass-through with a configurable placement.
    #[derive(Clone)]
    struct PlacedPar(PoolPlacement);
    impl Step for PlacedPar {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "PlacedPar",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![QueueSpec::Unbounded],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn pool_placement(&self) -> PoolPlacement {
            self.0
        }
        fn try_run(&mut self, _: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::Finished)
        }
        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    /// Detached pass-through (never run).
    struct DetachedPass;
    impl Step for DetachedPass {
        type Input = u32;
        type Outputs = Single<u32>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "DetachedPass",
                kind: StepKind::Detached,
                sticky: false,
                output_queues: vec![QueueSpec::Unbounded],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, _: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::Finished)
        }
    }

    /// Exclusive sink (never run).
    struct Sink;
    impl Step for Sink {
        type Input = u32;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Sink",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, _: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::Finished)
        }
    }

    #[rstest]
    #[case::reader(Affinity::Reader, Some(0))]
    #[case::writer_is_not_a_reader(Affinity::Writer, None)]
    #[case::pinned_worker_is_not_a_reader(Affinity::Worker(1), None)]
    #[case::none(Affinity::None, None)]
    fn reader_worker_finds_only_affinity_reader_serial_steps(
        #[case] affinity: Affinity,
        #[case] want: Option<usize>,
    ) {
        let b = PipelineBuilder::new();
        b.chain(AffinitySource(affinity))
            .chain(PlacedPar(PoolPlacement::AllWorkers))
            .chain(Sink)
            .into_sink_marker();
        let p = b.build().unwrap();
        assert_eq!(reader_worker(&p.steps, 4), want);
    }

    /// The producer lookup uses `ChainGraph::first_producer_into` and the
    /// producer's `driver_of` entry: hosted only when that producer is on a driver.
    #[test]
    fn plan_hosts_on_the_producers_driver_at_one_worker() {
        let b = PipelineBuilder::new();
        b.chain(AffinitySource(Affinity::Reader))
            .chain(DetachedPass)
            .chain(PlacedPar(PoolPlacement::ExcludeReader))
            .chain(Sink)
            .into_sink_marker();
        let p = b.build().unwrap();
        let driver_of = vec![None, Some(DriverIdx(0)), None, None];
        let at1 = plan_parallel_hosts(&p.steps, &p.graph, &driver_of, 1);
        assert_eq!(at1[2], ParallelHosts { workers: vec![], driver: Some(DriverIdx(0)) });
        let at4 = plan_parallel_hosts(&p.steps, &p.graph, &driver_of, 4);
        assert_eq!(at4[2].workers, vec![1, 2, 3]);
        for i in [0, 1, 3] {
            assert_eq!(at4[i], ParallelHosts::default(), "non-Parallel step {i}");
        }
        assert_eq!(drain_counter_inits(&p.steps, &at1), vec![1, 1, 1, 1]);
        assert_eq!(drain_counter_inits(&p.steps, &at4), vec![1, 1, 3, 1]);
        assert_eq!(
            wake_driver_of(&driver_of, &at1),
            vec![None, Some(DriverIdx(0)), Some(DriverIdx(0)), None],
            "the hosted step routes as its driver"
        );
        assert_eq!(wake_driver_of(&driver_of, &at4), driver_of);
    }

    /// Excluding the pinned reader worker passes the event-count check; at one
    /// worker the hosted clone leaves out worker 0, which the reader pins.
    #[rstest]
    #[case::excluded_reader(2, vec![1], None, vec![Some(0), None])]
    #[case::hosted_at_one(1, vec![], Some(DriverIdx(0)), vec![Some(0), None])]
    #[case::every_worker(2, vec![0, 1], None, vec![None, None])]
    fn excluded_workers_that_are_pinned_pass(
        #[case] n: usize,
        #[case] workers: Vec<usize>,
        #[case] driver: Option<DriverIdx>,
        #[case] pinned: Vec<Option<usize>>,
    ) {
        let hosts = [ParallelHosts::default(), ParallelHosts { workers, driver }];
        assert_excluded_workers_are_pinned(&hosts, &pinned, n);
    }

    /// A placement that leaves out an unpinned worker trips the check: that
    /// worker waits on the event-count, where `notify_one` could pick it.
    #[test]
    #[should_panic(expected = "leaves out worker 0, which is not pinned")]
    fn an_excluded_unpinned_worker_panics() {
        let hosts = [ParallelHosts::default(), ParallelHosts { workers: vec![1], driver: None }];
        assert_excluded_workers_are_pinned(&hosts, &[None, None], 2);
    }

    #[rstest]
    #[case::all_workers_ignores_reader(
        PoolPlacement::AllWorkers,
        Some(0),
        4,
        None,
        vec![0, 1, 2, 3],
        None
    )]
    #[case::exclude_reader_at_four(PoolPlacement::ExcludeReader, Some(0), 4, None, vec![1, 2, 3], None)]
    #[case::exclude_reader_without_reader_is_noop(
        PoolPlacement::ExcludeReader,
        None,
        4,
        None,
        vec![0, 1, 2, 3],
        None
    )]
    #[case::exclude_reader_at_two(PoolPlacement::ExcludeReader, Some(0), 2, None, vec![1], None)]
    #[case::exclude_reader_at_one_hosts_on_driver(
        PoolPlacement::ExcludeReader,
        Some(0),
        1,
        Some(DriverIdx(0)),
        vec![],
        Some(DriverIdx(0))
    )]
    #[case::exclude_reader_at_one_without_driver_runs_on_zero(
        PoolPlacement::ExcludeReader,
        Some(0),
        1,
        None,
        vec![0],
        None
    )]
    #[case::all_workers_at_one_never_hosts(
        PoolPlacement::AllWorkers,
        Some(0),
        1,
        Some(DriverIdx(0)),
        vec![0],
        None
    )]
    fn parallel_hosts_table(
        #[case] placement: PoolPlacement,
        #[case] reader: Option<usize>,
        #[case] n: usize,
        #[case] producer_driver: Option<DriverIdx>,
        #[case] workers: Vec<usize>,
        #[case] driver: Option<DriverIdx>,
    ) {
        let h = parallel_hosts(placement, reader, n, producer_driver);
        assert_eq!(h.workers, workers);
        assert_eq!(h.driver, driver);
        assert!(h.clone_count() >= 1, "a Parallel step always has at least one clone");
        assert_eq!(h.clone_count(), h.workers.len() + usize::from(h.driver.is_some()));
    }

    /// The excluded worker is the reader's index, whatever it is (here 2).
    #[test]
    fn exclusion_is_by_reader_index() {
        assert_eq!(
            parallel_hosts(PoolPlacement::ExcludeReader, Some(2), 4, None).workers,
            vec![0, 1, 3]
        );
    }
}
