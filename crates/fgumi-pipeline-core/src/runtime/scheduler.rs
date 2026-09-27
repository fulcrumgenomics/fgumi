//! Pluggable per-worker dispatch-order policy for the round-robin pool driver.
//!
//! The worker loop ([`run_worker_loop`](crate::runtime::run_worker_loop)) walks
//! each worker's *live* steps once per pass and runs the first that makes
//! progress. A [`Scheduler`] decides the ORDER of that walk — the only thing it
//! controls; it never changes which steps exist, the sticky source/sink
//! fast-path, or the Serial/Exclusive contention rules.
//!
//! Two policies ship:
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

/// Downstream-first dispatch, except while a stage has raised a shared
/// *refill* signal: then the steps up to and including that stage walk
/// upstream-first, and the rest stay downstream-first.
///
/// [`DrainFirstScheduler`] visits the source and early steps last, so while a
/// heavy downstream step always has work, nothing upstream runs and the next
/// batch of input is never read ahead; the stage later waits on just-in-time
/// decoding. A stage that needs input raises the signal, and while it is raised
/// its feeding steps (`refill_through` and everything before it) come first in
/// the walk. Only those steps move: turning the whole walk upstream-first would
/// also let the stage's own heavy consumers run ahead of the steps draining
/// older work, which grows the working set and costs CPU. The signal is read
/// once per worker pass with a relaxed load.
///
/// Read-ahead is capped: the refill walk applies only while the watched output
/// branch of `refill_through` (the queue feeding the stage) holds less than
/// `cap_bytes`, so the stage gets about one batch of work ready rather than
/// every upstream queue filling to its limit. Naming the branch keeps the cap
/// on the right queue when `refill_through` fans out (e.g. a rejects branch).
///
/// One instance serves one [`Pipeline::run`](crate::builder::Pipeline::run): it
/// names that pipeline's step and branch, and binds that run's queue once.
/// Build a fresh one per run rather than reusing a cloned `PipelineConfig`.
pub struct RefillDrainScheduler {
    refill: std::sync::Arc<std::sync::atomic::AtomicBool>,
    refill_through: crate::topology::StepIdx,
    refill_branch: crate::topology::BranchIdx,
    cap_bytes: u64,
    /// The watched output queue, resolved in [`Scheduler::bind`].
    queue: std::sync::OnceLock<std::sync::Arc<dyn crate::queues::BoundedQueueHandle>>,
}

impl std::fmt::Debug for RefillDrainScheduler {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RefillDrainScheduler")
            .field("refill_through", &self.refill_through)
            .field("refill_branch", &self.refill_branch)
            .field("cap_bytes", &self.cap_bytes)
            .field("bound", &self.queue.get().is_some())
            .finish_non_exhaustive()
    }
}

impl RefillDrainScheduler {
    /// A scheduler that walks `refill_through` and the steps before it
    /// upstream-first while `refill` is set and `refill_through`'s output
    /// branch `refill_branch` holds less than `cap_bytes`.
    #[must_use]
    pub fn new(
        refill: std::sync::Arc<std::sync::atomic::AtomicBool>,
        refill_through: crate::topology::StepIdx,
        refill_branch: crate::topology::BranchIdx,
        cap_bytes: u64,
    ) -> Self {
        Self { refill, refill_through, refill_branch, cap_bytes, queue: std::sync::OnceLock::new() }
    }

    /// Whether the refill walk applies now.
    fn refilling(&self) -> bool {
        self.refill.load(std::sync::atomic::Ordering::Relaxed)
            && self.queue.get().is_none_or(|q| q.current_bytes() < self.cap_bytes)
    }
}

impl Scheduler for RefillDrainScheduler {
    #[inline]
    fn walk(&self) -> WalkDirection {
        if self.refilling() {
            WalkDirection::RefillThenReverse(self.refill_through)
        } else {
            WalkDirection::Reverse
        }
    }

    fn name(&self) -> &'static str {
        "refill-drain"
    }

    /// Resolve the watched output branch's queue. If it is not byte-bounded
    /// the cap cannot be read and the refill walk follows the signal alone.
    ///
    /// Binding a second run's (different) queue is a contract violation — the
    /// `OnceLock` would keep watching the first run's queue — so it trips a
    /// `debug_assert!` and logs a warning in release.
    fn bind(&self, queues: &[crate::runtime::contexts::RegisteredQueue]) {
        if let Some(queue) = queues
            .iter()
            .find(|q| q.producer_step == self.refill_through && q.branch == self.refill_branch)
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
                "refill-drain: no byte-bounded queue for step {:?} branch {:?}; \
                 read-ahead is uncapped and follows the refill signal alone",
                self.refill_through,
                self.refill_branch,
            );
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

    #[test]
    fn refill_drain_walks_upstream_only_while_the_signal_is_raised() {
        use std::sync::Arc;
        use std::sync::atomic::{AtomicBool, Ordering};
        let refill = Arc::new(AtomicBool::new(false));
        let through = crate::topology::StepIdx(3);
        let scheduler = RefillDrainScheduler::new(
            Arc::clone(&refill),
            through,
            crate::topology::BranchIdx(0),
            100,
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
            Arc::new(AtomicBool::new(true)),
            through,
            crate::topology::BranchIdx(0),
            100,
        );
        assert_eq!(scheduler.name(), "refill-drain");
        let full = Arc::new(FakeQueue(AtomicU64::new(1_000)));
        scheduler.bind(&[registered(2, 0, &full), registered(3, 1, &full)]);

        assert!(scheduler.queue.get().is_none(), "no queue on the watched branch");
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
            Arc::new(AtomicBool::new(true)),
            crate::topology::StepIdx(0),
            crate::topology::BranchIdx(0),
            100,
        );
        let ours = Arc::new(FakeQueue(AtomicU64::new(0)));
        scheduler.bind(&[registered(0, 0, &ours)]);
        scheduler.bind(&[registered(0, 0, &ours)]);
        assert!(scheduler.queue.get().is_some(), "still bound after an idempotent rebind");
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "bound to a second run's queue")]
    fn refill_drain_rebind_to_another_queue_panics_in_debug() {
        use std::sync::Arc;
        use std::sync::atomic::{AtomicBool, AtomicU64};
        let scheduler = RefillDrainScheduler::new(
            Arc::new(AtomicBool::new(true)),
            crate::topology::StepIdx(0),
            crate::topology::BranchIdx(0),
            100,
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
            Arc::new(AtomicBool::new(true)),
            through,
            crate::topology::BranchIdx(1),
            100,
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
