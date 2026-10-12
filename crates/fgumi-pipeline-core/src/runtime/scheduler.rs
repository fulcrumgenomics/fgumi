//! Pluggable per-worker dispatch-order policy for the round-robin pool driver.
//!
//! The worker loop (`run_worker_loop`) walks
//! each worker's *live* steps once per pass and runs the first that makes
//! progress. A [`Scheduler`] decides the ORDER of that walk — the only thing it
//! controls; it never changes which steps exist, the sticky source/sink
//! fast-path, or the Serial/Exclusive contention rules.
//!
//! Three policies ship:
//!
//! - [`ChainOrderScheduler`] (the default) — walk **upstream-first** (chain
//!   order). A worker attempts the earliest-in-chain step with work, favouring
//!   production. This is the historical behaviour; every command keeps it unless
//!   it opts out, so existing pipelines are byte-for-byte unaffected.
//! - [`DrainFirstScheduler`] — walk **downstream-first** (reverse chain order).
//!   A worker attempts the deepest step with work first, favouring *draining*
//!   buffered work before producing more. Combined with skip-on-Serial-contention
//!   this self-balances: for a Serial step fed by an N-way Parallel producer, one
//!   worker grabs the Serial drain (mutex) while the rest find it contended, skip,
//!   and fall through to the producer — so the drain overlaps production instead
//!   of starving behind it on the shared pool. A sticky step is **not** exempt
//!   from the walk: it takes its bounded burst on the sticky fast-path first and
//!   is then still visited in the walk itself (under `Reverse`, last rather than
//!   first). Only its `Progress` priority restart is suppressed, so the walk
//!   continues past it to the steps that drain its output — see
//!   `super::driver::round_robin_dispatch`.
//! - [`RefillDrainScheduler`] — the chain's own direction (`Forward` or
//!   `Reverse`), except while a stage has raised one of the chain's refill
//!   hints: then the steps feeding that stage walk upstream-first and the rest
//!   downstream-first. A chain may carry several hints (one per stage).
//!
//! This mirrors main's `Scheduler`-trait design (`BalancedChaseDrainScheduler`
//! among others), but generically over an arbitrary step list rather than a
//! fixed set of named BAM stages: the only lever exposed here is walk direction,
//! which is all the generic driver needs to express drain-first scheduling.

/// The order in which a worker attempts its live steps in one round-robin pass.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WalkDirection {
    /// Chain order (upstream → downstream): favour production.
    Forward,
    /// Reverse chain order (downstream → upstream): favour draining.
    Reverse,
    /// Chain order through the given step, then reverse chain order for the
    /// steps after it: refill a starved stage's input first, while everything
    /// downstream of it keeps draining before it produces more.
    RefillThenReverse(crate::topology::StepIdx),
}

/// A per-worker dispatch-order policy. Selected per pipeline via
/// [`PipelineConfig::with_scheduler`](crate::builder::PipelineConfig::with_scheduler)
/// and shared across all workers (the shipped policies are stateless).
pub trait Scheduler: Send + Sync + std::fmt::Debug {
    /// Direction to walk this worker's live steps this pass. Called once per
    /// round-robin pass, so it must be cheap.
    fn walk(&self) -> WalkDirection;

    /// Human-readable name for diagnostics / `--pipeline-stats`.
    fn name(&self) -> &'static str;

    /// Called once per run with the run's byte-bounded queues, after they exist
    /// and before any worker dispatches. A policy that watches a queue resolves
    /// its handle here. The shipped stateless policies ignore it.
    fn bind(&self, _queues: &[crate::runtime::contexts::RegisteredQueue]) {}
}

/// Default upstream-first (chain-order) walk. Preserves the historical dispatch
/// behaviour for every pipeline that does not opt into a different policy.
#[derive(Debug, Default, Clone, Copy)]
pub struct ChainOrderScheduler;

impl Scheduler for ChainOrderScheduler {
    #[inline]
    fn walk(&self) -> WalkDirection {
        WalkDirection::Forward
    }
    fn name(&self) -> &'static str {
        "chain-order"
    }
}

/// Downstream-first (reverse chain-order) walk — drain buffered work before
/// producing more. Opt-in per pipeline (e.g. the sort chain, to overlap the
/// serial boundary/key scan with the parallel inflate instead of starving it on
/// the shared pool).
#[derive(Debug, Default, Clone, Copy)]
pub struct DrainFirstScheduler;

impl Scheduler for DrainFirstScheduler {
    #[inline]
    fn walk(&self) -> WalkDirection {
        WalkDirection::Reverse
    }
    fn name(&self) -> &'static str {
        "drain-first"
    }
}

/// One stage's refill hint as the scheduler sees it: while `signal` is raised
/// and the `refill_branch` output queue of `refill_through` holds less than
/// `cap_bytes`, the steps up to and including `refill_through` walk
/// upstream-first.
#[derive(Debug, Clone)]
pub struct RefillSource {
    /// Raised by the starved stage while it wants its input refilled.
    pub signal: std::sync::Arc<std::sync::atomic::AtomicBool>,
    /// The deepest step of the refill walk (the stage's feeding producer).
    pub refill_through: crate::topology::StepIdx,
    /// The output branch of `refill_through` whose depth caps the read-ahead.
    pub refill_branch: crate::topology::BranchIdx,
    /// Read-ahead cap on that branch's queue, in bytes (`u64::MAX` = uncapped).
    pub cap_bytes: u64,
}

/// A [`RefillSource`] plus its watched queue, resolved in [`Scheduler::bind`].
struct BoundRefillSource {
    source: RefillSource,
    queue: std::sync::OnceLock<std::sync::Arc<dyn crate::queues::BoundedQueueHandle>>,
}

impl BoundRefillSource {
    /// Whether this source's refill walk applies now.
    fn refilling(&self) -> bool {
        self.source.signal.load(std::sync::atomic::Ordering::Relaxed)
            && self.queue.get().is_none_or(|q| q.current_bytes() < self.source.cap_bytes)
    }

    /// Resolve the watched output branch's queue. If it is not byte-bounded
    /// the cap cannot be read and the refill walk follows the signal alone.
    ///
    /// Binding a second run's (different) queue is a contract violation — the
    /// `OnceLock` would keep watching the first run's queue — so it trips a
    /// `debug_assert!` and logs a warning in release.
    fn bind(&self, queues: &[crate::runtime::contexts::RegisteredQueue]) {
        let RefillSource { refill_through, refill_branch, .. } = self.source;
        if let Some(queue) =
            queues.iter().find(|q| q.producer_step == refill_through && q.branch == refill_branch)
        {
            if let Err(handle) = self.queue.set(std::sync::Arc::clone(&queue.handle)) {
                let same_queue = self.queue.get().is_some_and(|bound| {
                    std::ptr::addr_eq(
                        std::sync::Arc::as_ptr(bound),
                        std::sync::Arc::as_ptr(&handle),
                    )
                });
                debug_assert!(
                    same_queue,
                    "RefillDrainScheduler bound to a second run's queue; one instance serves \
                     one Pipeline::run"
                );
                if !same_queue {
                    log::warn!(
                        "refill-drain: scheduler reused across pipeline runs; its read-ahead \
                         cap keeps watching the first run's queue"
                    );
                }
            }
        } else {
            log::debug!(
                "refill-drain: no byte-bounded queue for step {refill_through:?} branch \
                 {refill_branch:?}; read-ahead is uncapped and follows the refill signal alone",
            );
        }
    }
}

impl std::fmt::Debug for BoundRefillSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RefillSource")
            .field("refill_through", &self.source.refill_through)
            .field("refill_branch", &self.source.refill_branch)
            .field("cap_bytes", &self.source.cap_bytes)
            .field("bound", &self.queue.get().is_some())
            .finish_non_exhaustive()
    }
}

/// The chain's own walk direction, except while one or more stages have
/// raised a *refill* signal: then the steps up to and including the deepest
/// raised stage's feeding step walk upstream-first, and the rest walk
/// downstream-first.
///
/// [`DrainFirstScheduler`] visits the source and early steps last, so while a
/// heavy downstream step always has work, nothing upstream runs and the next
/// batch of input is never read ahead; the stage later waits on just-in-time
/// decoding. A stage that needs input raises its signal, and while it is
/// raised its feeding steps (`refill_through` and everything before it) come
/// first in the walk. Only those steps move: turning the whole walk
/// upstream-first would also let the stage's own heavy consumers run ahead of
/// the steps draining older work, which grows the working set and costs CPU.
/// Each signal is read once per worker pass with a relaxed load.
///
/// A chain may carry several hints — an in-process aligner's (input starvation
/// during ingest) and a sort merge's (spill starvation in its merge phase) —
/// one [`RefillSource`] each. Their signals are not raised at the same time in
/// practice; if they are, the deepest raised `refill_through` wins, since its
/// walk covers the shallower one's steps too. While no signal is raised the
/// walk is `base`: `Reverse` for a chain that chose drain-first dispatch,
/// `Forward` for one that kept chain order, so the hints change nothing until
/// a stage starves.
///
/// Read-ahead is capped per source: a source's refill walk applies only while
/// the watched output branch of its `refill_through` (the queue feeding the
/// stage) holds less than its `cap_bytes`, so the stage gets about one batch of
/// work ready rather than every upstream queue filling to its limit. Naming the
/// branch keeps the cap on the right queue when `refill_through` fans out (e.g.
/// a rejects branch).
///
/// One instance serves one [`Pipeline::run`](crate::builder::Pipeline::run): it
/// names that pipeline's steps and branches, and binds that run's queues once.
/// Build a fresh one per run rather than reusing a cloned `PipelineConfig`.
pub struct RefillDrainScheduler {
    base: WalkDirection,
    sources: Box<[BoundRefillSource]>,
}

impl std::fmt::Debug for RefillDrainScheduler {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RefillDrainScheduler")
            .field("base", &self.base)
            .field("sources", &self.sources)
            .finish()
    }
}

impl RefillDrainScheduler {
    /// A scheduler that walks `base` while no source is refilling, and
    /// `RefillThenReverse(deepest refilling source's refill_through)` while
    /// some are.
    ///
    /// # Panics
    ///
    /// Panics if `base` is not [`WalkDirection::Forward`] or
    /// [`WalkDirection::Reverse`] (the base is a chain's own direction, never a
    /// refill walk).
    #[must_use]
    pub fn new(base: WalkDirection, sources: Vec<RefillSource>) -> Self {
        assert!(
            matches!(base, WalkDirection::Forward | WalkDirection::Reverse),
            "RefillDrainScheduler base must be Forward or Reverse, got {base:?}"
        );
        let sources = sources
            .into_iter()
            .map(|source| BoundRefillSource { source, queue: std::sync::OnceLock::new() })
            .collect();
        Self { base, sources }
    }
}

impl Scheduler for RefillDrainScheduler {
    #[inline]
    fn walk(&self) -> WalkDirection {
        self.sources
            .iter()
            .filter(|s| s.refilling())
            .map(|s| s.source.refill_through)
            .max_by_key(|step| step.0)
            .map_or(self.base, WalkDirection::RefillThenReverse)
    }

    fn name(&self) -> &'static str {
        "refill-drain"
    }

    /// Resolve each source's watched output branch's queue (see
    /// `BoundRefillSource::bind`).
    fn bind(&self, queues: &[crate::runtime::contexts::RegisteredQueue]) {
        for source in &self.sources {
            source.bind(queues);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_forward() {
        assert_eq!(ChainOrderScheduler.walk(), WalkDirection::Forward);
        assert_eq!(ChainOrderScheduler.name(), "chain-order");
    }

    #[test]
    fn drain_first_is_reverse() {
        assert_eq!(DrainFirstScheduler.walk(), WalkDirection::Reverse);
        assert_eq!(DrainFirstScheduler.name(), "drain-first");
    }

    // The position→step arithmetic these policies drive is tested where it
    // lives, in driver.rs: `walk_visits_live_steps_in_order` pins every walk's
    // visit order, and `driver_round_robins_all_live_before_parking` /
    // `refill_walk_reaches_the_consumer_of_a_holding_step` run it end-to-end
    // against the real dispatcher.

    fn source(
        through: usize,
        signal: &std::sync::Arc<std::sync::atomic::AtomicBool>,
    ) -> RefillSource {
        RefillSource {
            signal: std::sync::Arc::clone(signal),
            refill_through: crate::topology::StepIdx(through),
            refill_branch: crate::topology::BranchIdx(0),
            cap_bytes: u64::MAX,
        }
    }

    /// Each hint's signal moves the walk on its own; the deepest raised
    /// `refill_through` wins; with none raised the base direction applies.
    #[rstest::rstest]
    #[case::reverse_base(WalkDirection::Reverse)]
    #[case::forward_base(WalkDirection::Forward)]
    fn refill_drain_follows_each_hint_and_falls_back_to_its_base(#[case] base: WalkDirection) {
        use crate::topology::StepIdx;
        use std::sync::Arc;
        use std::sync::atomic::{AtomicBool, Ordering};
        let align = Arc::new(AtomicBool::new(false));
        let merge = Arc::new(AtomicBool::new(false));
        // Listed deepest-first, so a "first raised wins" walk is caught too.
        let s = RefillDrainScheduler::new(base, vec![source(7, &merge), source(3, &align)]);
        assert_eq!(s.walk(), base);
        align.store(true, Ordering::Relaxed);
        assert_eq!(s.walk(), WalkDirection::RefillThenReverse(StepIdx(3)));
        align.store(false, Ordering::Relaxed);
        merge.store(true, Ordering::Relaxed);
        assert_eq!(s.walk(), WalkDirection::RefillThenReverse(StepIdx(7)));
        align.store(true, Ordering::Relaxed);
        assert_eq!(s.walk(), WalkDirection::RefillThenReverse(StepIdx(7)), "deepest raised wins");
        merge.store(false, Ordering::Relaxed);
        assert_eq!(s.walk(), WalkDirection::RefillThenReverse(StepIdx(3)));
        align.store(false, Ordering::Relaxed);
        assert_eq!(s.walk(), base);
    }

    /// The deepest-wins rule also holds when the deeper hint is listed last.
    #[test]
    fn refill_drain_deepest_wins_regardless_of_list_order() {
        use std::sync::Arc;
        use std::sync::atomic::AtomicBool;
        let raised = Arc::new(AtomicBool::new(true));
        let s = RefillDrainScheduler::new(
            WalkDirection::Reverse,
            vec![source(3, &raised), source(7, &raised)],
        );
        assert_eq!(s.walk(), WalkDirection::RefillThenReverse(crate::topology::StepIdx(7)));
    }

    /// A refill walk is not a base direction.
    #[test]
    #[should_panic(expected = "base must be Forward or Reverse")]
    fn refill_drain_rejects_a_refill_walk_as_its_base() {
        let _ = RefillDrainScheduler::new(
            WalkDirection::RefillThenReverse(crate::topology::StepIdx(1)),
            Vec::new(),
        );
    }

    #[test]
    fn refill_drain_walks_upstream_only_while_the_signal_is_raised() {
        use std::sync::Arc;
        use std::sync::atomic::{AtomicBool, Ordering};
        let refill = Arc::new(AtomicBool::new(false));
        let through = crate::topology::StepIdx(3);
        let scheduler = RefillDrainScheduler::new(
            WalkDirection::Reverse,
            vec![RefillSource {
                signal: Arc::clone(&refill),
                refill_through: through,
                refill_branch: crate::topology::BranchIdx(0),
                cap_bytes: 100,
            }],
        );
        assert_eq!(scheduler.walk(), WalkDirection::Reverse);
        refill.store(true, Ordering::Relaxed);
        assert_eq!(scheduler.walk(), WalkDirection::RefillThenReverse(through));
        refill.store(false, Ordering::Relaxed);
        assert_eq!(scheduler.walk(), WalkDirection::Reverse);
    }

    /// A byte-bounded queue stand-in whose depth the test sets directly.
    struct FakeQueue(std::sync::atomic::AtomicU64);

    impl crate::queues::BoundedQueueHandle for FakeQueue {
        fn current_bytes(&self) -> u64 {
            self.0.load(std::sync::atomic::Ordering::Relaxed)
        }
        fn limit_bytes(&self) -> u64 {
            u64::MAX
        }
        fn set_limit_bytes(&self, _new_limit: u64) {}
        fn enable_tracking(&self, _t: crate::queues::EdgeTracking, _: crate::queues::Sealed) {}
        fn pushed_total(&self, _: crate::queues::Sealed) -> u64 {
            0
        }
        fn take_holders(&self, _f: &mut dyn FnMut(usize), _: crate::queues::Sealed) -> usize {
            0
        }
        fn is_tracked(&self, _: crate::queues::Sealed) -> bool {
            false
        }
    }

    fn registered(
        producer: usize,
        branch: usize,
        queue: &std::sync::Arc<FakeQueue>,
    ) -> crate::runtime::contexts::RegisteredQueue {
        crate::runtime::contexts::RegisteredQueue {
            producer_step_name: "producer",
            producer_step: crate::topology::StepIdx(producer),
            branch: crate::topology::BranchIdx(branch),
            handle: std::sync::Arc::clone(queue) as _,
            reorder_cap: None,
        }
    }

    /// Without a byte-bounded queue on the watched branch, `bind` leaves the
    /// scheduler unbound: the read-ahead cap cannot be read, so the refill walk
    /// follows the signal alone however full other queues are. `Debug` reports
    /// the binding state.
    #[test]
    fn refill_drain_without_a_watched_queue_follows_the_signal_alone() {
        use std::sync::Arc;
        use std::sync::atomic::{AtomicBool, AtomicU64};
        let through = crate::topology::StepIdx(3);
        let scheduler = RefillDrainScheduler::new(
            WalkDirection::Reverse,
            vec![RefillSource {
                signal: Arc::new(AtomicBool::new(true)),
                refill_through: through,
                refill_branch: crate::topology::BranchIdx(0),
                cap_bytes: 100,
            }],
        );
        assert_eq!(scheduler.name(), "refill-drain");
        let full = Arc::new(FakeQueue(AtomicU64::new(1_000)));
        scheduler.bind(&[registered(2, 0, &full), registered(3, 1, &full)]);

        assert!(scheduler.sources[0].queue.get().is_none(), "no queue on the watched branch");
        assert_eq!(scheduler.walk(), WalkDirection::RefillThenReverse(through));
        let debug = format!("{scheduler:?}");
        assert!(debug.contains("bound: false"), "{debug}");

        let ours = Arc::new(FakeQueue(AtomicU64::new(0)));
        scheduler.bind(&[registered(3, 0, &ours)]);
        let debug = format!("{scheduler:?}");
        assert!(debug.contains("bound: true"), "{debug}");
    }

    /// Re-binding the same queue is harmless (idempotent), but binding a
    /// second run's different queue violates the one-run contract and trips
    /// the debug assertion rather than silently watching the stale queue.
    #[test]
    fn refill_drain_rebind_to_the_same_queue_is_idempotent() {
        use std::sync::Arc;
        use std::sync::atomic::{AtomicBool, AtomicU64};
        let scheduler = RefillDrainScheduler::new(
            WalkDirection::Reverse,
            vec![RefillSource {
                signal: Arc::new(AtomicBool::new(true)),
                refill_through: crate::topology::StepIdx(0),
                refill_branch: crate::topology::BranchIdx(0),
                cap_bytes: 100,
            }],
        );
        let ours = Arc::new(FakeQueue(AtomicU64::new(0)));
        scheduler.bind(&[registered(0, 0, &ours)]);
        scheduler.bind(&[registered(0, 0, &ours)]);
        assert!(
            scheduler.sources[0].queue.get().is_some(),
            "still bound after an idempotent rebind"
        );
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "bound to a second run's queue")]
    fn refill_drain_rebind_to_another_queue_panics_in_debug() {
        use std::sync::Arc;
        use std::sync::atomic::{AtomicBool, AtomicU64};
        let scheduler = RefillDrainScheduler::new(
            WalkDirection::Reverse,
            vec![RefillSource {
                signal: Arc::new(AtomicBool::new(true)),
                refill_through: crate::topology::StepIdx(0),
                refill_branch: crate::topology::BranchIdx(0),
                cap_bytes: 100,
            }],
        );
        let first = Arc::new(FakeQueue(AtomicU64::new(0)));
        let second = Arc::new(FakeQueue(AtomicU64::new(0)));
        scheduler.bind(&[registered(0, 0, &first)]);
        scheduler.bind(&[registered(0, 0, &second)]);
    }

    /// Once bound to the watched output branch, the refill walk stops at the
    /// cap even while the signal is raised, and resumes below it. Queues from
    /// another step, or from the refill step's other branch (a fan-out such as
    /// a rejects output), are ignored.
    #[test]
    fn refill_drain_stops_reading_ahead_at_the_cap() {
        use std::sync::Arc;
        use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
        let through = crate::topology::StepIdx(3);
        let scheduler = RefillDrainScheduler::new(
            WalkDirection::Reverse,
            vec![RefillSource {
                signal: Arc::new(AtomicBool::new(true)),
                refill_through: through,
                refill_branch: crate::topology::BranchIdx(1),
                cap_bytes: 100,
            }],
        );
        let ours = Arc::new(FakeQueue(AtomicU64::new(0)));
        let other = Arc::new(FakeQueue(AtomicU64::new(1_000)));
        scheduler.bind(&[
            registered(2, 1, &other),
            registered(3, 0, &other),
            registered(3, 1, &ours),
        ]);

        assert_eq!(scheduler.walk(), WalkDirection::RefillThenReverse(through));
        ours.0.store(100, Ordering::Relaxed);
        assert_eq!(scheduler.walk(), WalkDirection::Reverse, "at the cap");
        ours.0.store(99, Ordering::Relaxed);
        assert_eq!(scheduler.walk(), WalkDirection::RefillThenReverse(through));
    }
}
