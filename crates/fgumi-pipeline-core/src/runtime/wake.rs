//! Directed wake routing (`WakePlan`): every scheduled chain's wakes.
//!
//! A blanket wake — every `Progress` notifying one parked pool worker
//! (`notify_one`) whatever consumes what was produced — wakes the wrong
//! thread when the consumer is not a pool step, and never reaches a dedicated
//! driver thread. So every scheduled chain routes its wakes through a plan.
//!
//! A `WakePlan` is a per-chain table, built once in `Pipeline::run`, that says
//! for each `(step, output branch)` *who* consumes it and therefore who to
//! wake: nobody (same thread), the pool event-count, a driver thread, or one
//! pinned pool worker. Drivers and workers are woken with
//! `std::thread::Thread::unpark`, and only when armed for their timer park
//! (the parking bit, below).
//!
//! A forward wake is push-gated: a producer's branch with a byte-bounded
//! transport wakes its consumer only when the dispatch actually pushed to it
//! (the transport's push counter advanced; an ordered branch's stash insert
//! counts), so progress that only held, accumulated or retried wakes nobody.
//!
//! A held item wakes its holder: every bounded edge the plan tracks records
//! the refusing thread's wake slot (`runtime::wake_slot`), and the pop (or, on
//! a byte-bounded edge, a limit raise) that makes room unparks exactly those
//! threads. On an ordered branch the release is the in-order pop, and the
//! stage's fast path re-checks `next_serial_buffered` after the latch; an
//! unbounded edge refuses only at an ordered branch's stash cap.
//!
//! A flushed held retry (`BranchOutputHandle::retry` succeeded) gets the same
//! forward wake as a `Progress`, whatever the dispatch then reports, and a
//! same-thread consumer gets a re-poll (`on_flushed`).
//!
//! Every scheduled chain gets a plan, with or without a `Detached` step: every
//! pool worker and driver registers itself, and a `Pool` wake that finds no
//! event-count waiter reaches an armed direct-parked worker, so a chain whose
//! workers are all pinned still runs on wakes. The fused single-thread path
//! uses [`WakePlan::standalone`], which routes nothing.
//!
//! A directed wake (a `Worker`/`Driver` forward target, a holder or
//! cap-release wake) unparks only a thread armed for its timer park — one bit
//! per wake slot (`parking`). A running thread costs the deliverer a load, not
//! an `unpark`: every timer park first sets its bit, fences and re-polls, so a
//! wake that finds the bit clear is seen by that re-poll. A holder wake reads
//! the bit behind a second fence, after its claim (`WakePlan::deliver_holder`).
//! A plan with one thread neither arms nor gates: its lone thread is unparked
//! by every directed wake, since the queue-memory rebalancer, outside the
//! plan, still delivers there.
//!
//! Lost-wakeup arguments, in short: the pool path is the event-count's own
//! fence pair; the unpark path is closed by the std park token (`unpark` before
//! `park` makes the park return at once) plus two `SeqCst` fence pairs with the
//! delivery — registration's (`ThreadSlots`) and the park arm's (`parking`:
//! the thread sets its bit, fences and re-polls; the deliverer fences and
//! reads the bit; a holder take fences once more between its claim and that
//! read); the reverse path is a Dekker
//! pair — the producer latches its slot, fences and re-checks for room
//! (`HolderSet::latch_and_recheck`, one implementation for every transport);
//! the consumer pops, fences and takes the latches; the direct-park fallback is
//! the same shape over the armed bits (`DirectParked`: the worker arms, fences
//! and re-polls; the producer publishes, fences inside `notify_one` and claims
//! a bit). The fences of these protocols live in `runtime::wake_slot` and, for
//! the producer half of the pool and direct-park paths, in
//! `PoolEventCount::notify_one` (`event_count.rs`); `tests/loom_wake.rs`
//! model-checks both as they run here.
//!
//! The plan is crate-internal: nothing outside the crate builds, reads or
//! drives one, except the loom models, through a `cfg(loom)` re-export.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::admission::PhaseCap;
use crate::padded::Padded;
use crate::queues::{BoundedQueueHandle, HolderQueue, SEALED};
use crate::runtime::contexts::{RegisteredHolderOnlyQueue, RegisteredQueue};
use crate::runtime::event_count::{NotifyOutcome, PoolEventCount};
use crate::runtime::placement::ParallelHosts;
use crate::runtime::stats::WakeCounts;
use crate::runtime::wake_slot::{
    DirectParked, ParkThread, ThreadSlots, arm_fence, current_slot, delivery_fence, set_cap_owed,
    take_cap_owed, take_cap_recorded,
};
use crate::step::StepKind;
use crate::topology::{BranchIdx, ChainGraph, StepIdx};

pub use crate::runtime::placement::DriverIdx;

/// Who a forward wake goes to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WakeTarget {
    /// Nobody: the consumer runs on the thread that just produced, or the
    /// branch has no consumer. (A step's `Finished` is still broadcast.)
    None,
    /// The pool event-count (`notify_one`), falling back to one direct-parked
    /// worker when no event-count waiter exists. A pinned worker is never an
    /// event-count waiter, but it arms the fallback before its timer park, so a
    /// chain whose workers are all pinned still receives `Pool` work.
    Pool,
    /// `Thread::unpark` the driver hosting a `Detached` consumer.
    Driver(DriverIdx),
    /// `Thread::unpark` pool worker `w`: the consumer is pinned to `w`, or
    /// `n_threads == 1` (the lone pool worker is `Worker(0)`).
    Worker(usize),
}

/// Size of the per-step push-counter snapshot. A typed step has at most
/// `outputs::MAX_ARITY` branches, so every gateable branch of a producer fits.
pub const MAX_GATED_BRANCHES: usize = 4;
const _: () = assert!(MAX_GATED_BRANCHES >= crate::outputs::MAX_ARITY);

/// Push counters of a step's gated branches before a dispatch.
#[derive(Debug, Default, Clone, Copy)]
pub struct PushSnapshot([u64; MAX_GATED_BRANCHES]);

struct BranchWake {
    target: WakeTarget,
    /// `Some` ⇒ wake only if `pushed_total()` advanced during the dispatch (a
    /// reorder-stash insert counts). Set for every producer's byte-bounded
    /// branch, whatever the producer's kind; `None` on a count-bounded or
    /// unbounded transport, which has no push counter.
    gate: Option<Arc<dyn BoundedQueueHandle>>,
    /// The branch has a consumer and it runs on this step's thread (`target`
    /// is then `None`). Read only by `on_flushed`.
    same_thread: bool,
    /// The branch's consumer, for its phase cap (`WakePlan::caps`, stored
    /// once per step). A `Pool` wake for it claims a cap-parked worker only
    /// while that cap admits (`PhaseCap::has_free_permit`): when it is full,
    /// the worker would only be refused again, and the cap's release wakes the
    /// threads recorded on it.
    consumer: Option<StepIdx>,
    /// `Some` ⇒ a `Pool` target whose consumer has clones on only some workers
    /// (`ParallelHosts::workers`): one bit per host worker. The direct-park and
    /// cap-parked fallbacks may claim only a worker in it — any other worker
    /// holds no clone of the consumer, so waking it would leave the work for a
    /// timer.
    fallback_mask: Option<Box<[u64]>>,
}

struct ReverseWake {
    queue: HolderQueue,
    producer: StepIdx,
}

struct StepWake {
    forward: Box<[BranchWake]>,
    reverse: Box<[ReverseWake]>,
    /// Precomputed: any gated branch (the dispatcher snapshots only then).
    gated: bool,
    /// Precomputed: a `Driver`/`Worker` target or a reverse edge exists, so the
    /// unpark/holder fence is needed. A step whose only target is `Pool` skips
    /// it (`notify_one` has its own).
    needs_fence: bool,
    /// This step's own thread (`Driver(d)` or `Worker(w)`) when some forward
    /// branch is `same_thread`; `None` otherwise. Read only by `on_flushed`.
    own: WakeTarget,
}

/// The bounded edges a plan may track, by kind.
#[derive(Clone, Copy, Default)]
pub struct WakeEdges<'a> {
    /// Gated (every producer's branch on one of these edges, ordered branches
    /// included: a reorder-stash insert counts as a push) and reverse-woken.
    pub byte_bounded: &'a [RegisteredQueue],
    /// Reverse-woken only: there is no push counter to gate on.
    pub(crate) holder_only: &'a [RegisteredHolderOnlyQueue],
}

impl WakeEdges<'_> {
    /// No edges (a plan that routes forward wakes only, as `Pipeline::dag_at`).
    pub const NONE: Self = Self { byte_bounded: &[], holder_only: &[] };
}

#[cfg(test)]
impl<'a> WakeEdges<'a> {
    /// Byte-bounded edges only.
    #[must_use]
    pub fn byte_bounded(queues: &'a [RegisteredQueue]) -> Self {
        Self { byte_bounded: queues, holder_only: &[] }
    }
}

/// The reverse-wake edges of `step`: every input edge of `step` in either
/// registry whose producer runs on another thread, whatever that producer's
/// role (a same-thread holder runs the consumer on its next pass).
/// Byte-bounded inputs first, then holder-only, each in registry order.
fn reverse_edges<'a>(
    graph: &'a ChainGraph,
    edges: WakeEdges<'a>,
    step: StepIdx,
    same_thread: &'a dyn Fn(StepIdx, StepIdx) -> bool,
) -> impl Iterator<Item = (StepIdx, BranchIdx, HolderQueue)> + 'a {
    let byte = edges
        .byte_bounded
        .iter()
        .map(|q| (q.producer_step, q.branch, HolderQueue::Bytes(Arc::clone(&q.handle))));
    let holder_only = edges
        .holder_only
        .iter()
        .map(|q| (q.producer_step, q.branch, HolderQueue::HolderOnly(Arc::clone(&q.handle))));
    byte.chain(holder_only)
        .filter(move |(p, b, _)| graph.consumer(*p, *b) == Some(step) && !same_thread(*p, step))
}

/// Why a `request_worker` did or did not wake anyone.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PoolRequest {
    /// A parked pool worker other than the requester was woken: an
    /// event-count waiter (ungated by the ceiling), else a worker that armed
    /// for a timer park, else one parked only because a phase cap refused it,
    /// or — at one thread — worker 0 was unparked. A cap-parked
    /// worker is woken whatever its caps' occupancy (a request names no
    /// step), so if the requested work sits behind a full cap that wake is
    /// wasted: the worker is refused again and re-records on that cap.
    Woken,
    /// No other worker to wake: besides the requester itself, no event-count
    /// waiter, no worker armed for a timer park, and no worker parked only
    /// because a phase cap refused it. A worker idling on its timer without
    /// arming (a pinned one between arms) is not reachable
    /// here; it finds the work on its next pass.
    AllAwake,
    /// An earlier request woke an event-count waiter, and no worker has left
    /// the event-count since: suppressed. At most one such request is in flight
    /// per pipeline. The direct-parked fallback is not rate-limited this way:
    /// each request claims (and so clears) one armed worker's bit, so it cannot
    /// wake the same worker twice.
    Pending,
    /// No pool event-count and no registered worker 0 (the fused path's
    /// standalone plan).
    Unavailable,
}

/// Step-side view of the pool, handed to every `try_run` via the context.
/// `Copy`; three words.
///
/// Call [`Self::request_worker`] only when this `try_run` is about to return
/// `NoProgress`/`Contention` because it is waiting on work that a **pool** step
/// produces and that reaches it **off any queue edge** (a decompressed block
/// landing in a side-channel slot). For queue-edge inputs the wake plan already
/// woke the consumer when the item was pushed. A step that calls it anyway is
/// visible in `pool_requests_{woken,all_awake,pending}` and costs one fence
/// and two loads while the pool is awake. The requester is never the worker a
/// request targets: a worker armed on the event-count (its re-poll) is not
/// counted as a waiter, and its own armed timer-park bit is skipped. An armed
/// requester whose request does go to another event-count waiter still sees
/// the generation that notify bumps: its own pending wait returns at once, so
/// it makes one more pass instead of parking, and it acknowledges the request
/// (clearing [`PoolRequest::Pending`]) whether or not the target has run yet.
#[derive(Clone, Copy)]
pub struct PoolHandle<'a> {
    plan: &'a WakePlan,
    stats: Option<&'a crate::runtime::stats::PipelineStats>,
    step: StepIdx,
}

impl<'a> PoolHandle<'a> {
    /// A handle for `step` over `plan`, counting outcomes into `stats` when present.
    #[must_use]
    pub(crate) fn new(
        plan: &'a WakePlan,
        stats: Option<&'a crate::runtime::stats::PipelineStats>,
        step: StepIdx,
    ) -> Self {
        Self { plan, stats, step }
    }

    /// Ask for ONE parked pool worker, other than the caller, to run a pass
    /// now. Never blocks, takes no queue lock and is safe from any thread. An
    /// event-count wake is idempotent until a worker leaves the event-count
    /// ([`PoolRequest::Pending`]); a direct-park wake claims one armed worker.
    /// Cost when nobody is parked: one `SeqCst` fence, the waiter-count load
    /// and one load per 64 workers of the direct-park bits.
    #[inline]
    #[must_use = "the outcome says whether a worker was woken"]
    pub fn request_worker(&self) -> PoolRequest {
        let r = self.plan.request_pool_worker();
        if let Some(stats) = self.stats {
            stats.record_pool_request(self.step, r);
        }
        r
    }

    /// Pool workers currently armed or parked on the event-count (diagnostics
    /// only); 0 without a pool.
    #[inline]
    #[must_use]
    pub fn parked_workers(&self) -> usize {
        self.plan.pool.as_ref().map_or(0, |p| p.waiters())
    }
}

/// The per-chain wake routing table. See the module docs.
pub struct WakePlan {
    /// By `StepIdx`; empty in the standalone plan.
    steps: Box<[StepWake]>,
    /// By `DriverIdx`.
    drivers: ThreadSlots,
    /// By worker id; every worker registers itself before its first pass.
    workers: ThreadSlots,
    /// The pool workers armed for a timer park.
    direct_parked: DirectParked,
    /// The pool workers parked on their timer only because a phase cap
    /// refused them: claimed after `direct_parked`, and only for a `Pool` wake
    /// whose consumer is uncapped or capped by a cap with a free permit.
    cap_parked: DirectParked,
    /// One bit per wake slot (workers `0..n`, then drivers): the threads
    /// armed for a timer park. A directed wake (a `Worker`/`Driver` forward
    /// target, a holder or cap-release wake through [`Self::deliver_slot`])
    /// unparks only a thread whose bit is set; a running thread costs the
    /// deliverer one shared load instead of an `unpark`. Sound because every
    /// timer park arms first: set the bit, fence, re-poll ([`Self::arm_direct`],
    /// [`Self::arm_cap_parked`], [`Self::arm_driver`]), and every deliverer
    /// fences after its publish and before it reads the bit, so a delivery
    /// either finds the bit or the re-poll finds the work. Cleared by the owner
    /// after it wakes. Not used in a plan with one thread ([`Self::gates`]).
    parking: DirectParked,
    /// Whether directed wakes are gated on `parking`: the plan has more than
    /// one thread (workers plus drivers). A plan with one thread never arms —
    /// its own steps deliver nothing across threads — and delivers every
    /// directed wake unconditionally instead, because threads outside the plan
    /// still deliver there: the queue-memory rebalancer's holder wakes
    /// ([`Self::wake_holders`]) reach the lone worker of a one-thread run.
    gates: bool,
    /// By `StepIdx`: the phase cap the step reports, if any. Empty in the
    /// standalone plan.
    caps: Box<[Option<Arc<PhaseCap>>]>,
    /// Each cap in `caps` once, in a fixed order: a cap's position here is its
    /// run index (`PhaseCap::set_plan_index`), the bit a worker's
    /// thread-local record masks use for it, which is what lets a pass start
    /// tell which wakes the thread owes ([`Self::start_pass`]). The dedup is
    /// what makes the index well defined for a cap several steps share.
    /// Empty in the standalone plan.
    distinct_caps: Box<[Arc<PhaseCap>]>,
    /// `None` when `n_threads == 1` (no event-count exists).
    pool: Option<Arc<PoolEventCount>>,
    /// Edges `build` gates or reverse-wakes, by `(producer, branch)`.
    tracked: Box<[(StepIdx, BranchIdx)]>,
    /// One unacknowledged `request_worker` in flight per pipeline; cleared by
    /// any pool worker's return from its event-count wait. Its own cache line:
    /// it is written by requesters and by every waking worker.
    request_pending: Padded<AtomicBool>,
}

/// Each distinct cap among `caps` once, in step order, each told its index
/// in the list (`PhaseCap::set_plan_index`): the bit a worker's record masks
/// use for it.
///
/// # Panics
///
/// Panics if `caps` holds more than `MAX_PHASE_CAPS` distinct caps.
fn index_caps(caps: &[Option<Arc<PhaseCap>>]) -> Box<[Arc<PhaseCap>]> {
    let distinct = crate::admission::distinct_caps(caps.iter().flatten().map(|c| &**c));
    // `PipelineBuilder::build` refuses such a pipeline; this is the backstop
    // for a plan built any other way.
    assert!(
        distinct.len() <= crate::admission::MAX_PHASE_CAPS,
        "WakePlan::build: a run may use at most {} distinct phase caps",
        crate::admission::MAX_PHASE_CAPS
    );
    distinct
        .into_iter()
        .enumerate()
        .map(|(i, cap)| {
            cap.set_plan_index(i);
            cap.arc().expect("the plan holds the cap")
        })
        .collect()
}

/// One bit per host worker of `consumer` when its clones live on only some of
/// the `n_threads` workers (`ParallelHosts::workers`); `None` when it runs on
/// every worker (or `hosts` does not cover it).
fn host_mask(hosts: &[ParallelHosts], consumer: StepIdx, n_threads: usize) -> Option<Box<[u64]>> {
    let h = hosts.get(consumer.0)?;
    if h.workers.is_empty() || h.workers.len() >= n_threads {
        return None;
    }
    Some(DirectParked::mask_for(&h.workers, n_threads))
}

impl WakePlan {
    /// The fused single-thread path's plan: no pool, no drivers, no workers.
    /// Every forward target is `None`, `on_progress` delivers nothing,
    /// registration is a no-op and `request_pool_worker` is `Unavailable`.
    /// A fresh `Arc` per run; it holds nothing.
    #[must_use]
    pub fn standalone() -> Arc<Self> {
        Arc::new(Self {
            steps: Box::new([]),
            drivers: ThreadSlots::new(0),
            workers: ThreadSlots::new(0),
            direct_parked: DirectParked::new(0),
            cap_parked: DirectParked::new(0),
            parking: DirectParked::new(0),
            gates: false,
            caps: Box::new([]),
            distinct_caps: Box::new([]),
            pool: None,
            tracked: Box::new([]),
            request_pending: Padded(AtomicBool::new(false)),
        })
    }

    /// Build the routing table from the chain. Pure over its inputs.
    /// `pinned_worker[i]` is the one pool worker that runs step `i` (Serial
    /// affinity target, range-checked, or Exclusive owner); `driver_of[i]` the
    /// driver hosting a Detached step; `caps[i]` the phase cap step `i`
    /// reports (a step past the end of a shorter slice, `&[]` included, is
    /// uncapped); `edges` the bounded transports, by kind.
    ///
    /// Gates: every producer's byte-bounded branch, whatever the producer's
    /// kind, wakes its consumer only if the dispatch pushed; on an ordered
    /// branch a reorder-stash insert counts as a push
    /// (`ItemQueue::note_stash_accepted`). A holder-only (count-bounded) edge,
    /// an unbounded edge and a zero-output Parallel step's self-cascade have no
    /// push counter and stay ungated: they wake on every `Progress`. Reverse
    /// edges: every input
    /// edge in either registry whose producer runs on another thread, whatever
    /// that producer's role.
    ///
    /// # Panics
    ///
    /// Panics if `kinds`, `pinned_worker` and `driver_of` are not all
    /// `graph.n_steps()` long (a `Pipeline::run` bug), or if `caps` holds more
    /// than [`MAX_PHASE_CAPS`](crate::MAX_PHASE_CAPS) distinct caps.
    ///
    /// Tests and the loom models only: `Pipeline::run` and `dag_at` place
    /// Parallel clones and call [`Self::build_with_hosts`].
    #[cfg(any(test, loom))]
    #[allow(clippy::too_many_arguments)] // the chain's per-step inputs, one slice each
    #[must_use]
    pub fn build(
        graph: &ChainGraph,
        kinds: &[StepKind],
        pinned_worker: &[Option<usize>],
        driver_of: &[Option<DriverIdx>],
        caps: &[Option<Arc<PhaseCap>>],
        edges: WakeEdges<'_>,
        pool: Option<Arc<PoolEventCount>>,
        n_threads: usize,
    ) -> Arc<Self> {
        Self::build_with_hosts(
            graph,
            kinds,
            pinned_worker,
            driver_of,
            caps,
            edges,
            pool,
            n_threads,
            &[],
        )
    }

    /// `build` (same inputs) plus where each `Parallel` step's clones live
    /// (`plan_parallel_hosts`; empty ⇒ every Parallel step on every worker). A
    /// `Pool` forward branch whose consumer has a restricted worker set records
    /// that set, and the direct-park and cap-parked fallbacks claim only a
    /// worker in it. The event-count path needs no restriction because every
    /// worker a placement leaves out is pinned, so never an event-count waiter;
    /// `run` and `dag_at` check that with
    /// `placement::assert_excluded_workers_are_pinned` before building the plan.
    ///
    /// # Panics
    ///
    /// As `build`.
    #[must_use]
    #[allow(clippy::too_many_arguments)] // `build`'s inputs plus the placement
    pub fn build_with_hosts(
        graph: &ChainGraph,
        kinds: &[StepKind],
        pinned_worker: &[Option<usize>],
        driver_of: &[Option<DriverIdx>],
        caps: &[Option<Arc<PhaseCap>>],
        edges: WakeEdges<'_>,
        pool: Option<Arc<PoolEventCount>>,
        n_threads: usize,
        hosts: &[ParallelHosts],
    ) -> Arc<Self> {
        let n = graph.n_steps();
        assert!(
            kinds.len() == n && pinned_worker.len() == n && driver_of.len() == n,
            "WakePlan::build: per-step inputs must be {n} long"
        );
        let n_drivers = driver_of.iter().flatten().map(|d| d.0 + 1).max().unwrap_or(0);

        // Who to wake so that `step` runs. `Pool` resolves to `Worker(0)` at one
        // thread (no event-count exists; the lone worker parks on its timer).
        let thread_of = |step: StepIdx| -> WakeTarget {
            if let Some(d) = driver_of[step.0] {
                WakeTarget::Driver(d)
            } else if let Some(w) = pinned_worker[step.0].filter(|&w| w < n_threads) {
                WakeTarget::Worker(w)
            } else if n_threads == 1 {
                WakeTarget::Worker(0)
            } else {
                WakeTarget::Pool
            }
        };
        // Same OS thread ⇒ no wake needed (the next pass visits it).
        let same_thread = |a: StepIdx, b: StepIdx| -> bool {
            match (thread_of(a), thread_of(b)) {
                (WakeTarget::Driver(x), WakeTarget::Driver(y)) => x == y,
                (WakeTarget::Worker(x), WakeTarget::Worker(y)) => x == y,
                _ => false,
            }
        };
        let queue_for = |step: StepIdx, branch: BranchIdx| -> Option<&RegisteredQueue> {
            edges.byte_bounded.iter().find(|q| q.producer_step == step && q.branch == branch)
        };
        debug_assert!(
            edges.holder_only.iter().all(|h| queue_for(h.producer_step, h.branch).is_none()),
            "WakePlan::build: an edge is in both registries"
        );

        let mask_of = |c: StepIdx| host_mask(hosts, c, n_threads);

        let mut tracked = Vec::new();
        let steps = (0..n)
            .map(|i| {
                let step = StepIdx(i);
                let mut forward: Box<[BranchWake]> = (0..graph.branch_count(step))
                    .map(|b| {
                        let branch = BranchIdx(b);
                        let consumer = graph.consumer(step, branch);
                        let same = matches!(consumer, Some(c) if same_thread(step, c));
                        let target = match consumer {
                            None => WakeTarget::None,
                            Some(_) if same => WakeTarget::None,
                            Some(c) => thread_of(c),
                        };
                        // Gate every producer whose branch has a byte-bounded
                        // transport. A `Progress` that pushed nothing (a held
                        // retry, an accumulation step, a Parallel clone that lost
                        // the race to pop) wakes nobody; a concurrent clone's push
                        // can only turn "not pushed" into "pushed" — a spurious,
                        // safe wake. An ordered branch's stash insert counts as a
                        // push (`ItemQueue::note_stash_accepted`). Non-byte-bounded edges
                        // have no counter and stay ungated.
                        let gate = queue_for(step, branch).map(|q| Arc::clone(&q.handle));
                        if gate.is_some() {
                            tracked.push((step, branch));
                        }
                        let fallback_mask =
                            consumer.filter(|_| target == WakeTarget::Pool).and_then(mask_of);
                        BranchWake { target, gate, same_thread: same, consumer, fallback_mask }
                    })
                    .collect();
                if forward.is_empty() && kinds[i] == StepKind::Parallel && n_threads > 1 {
                    forward = Box::new([Self::self_cascade(step, mask_of(step))]);
                }
                // One list from both registries, built before `needs_fence` so a
                // holder-only reverse edge brings the consumer's fence with it.
                let reverse: Box<[ReverseWake]> = reverse_edges(graph, edges, step, &same_thread)
                    .map(|(producer, branch, queue)| {
                        tracked.push((producer, branch));
                        ReverseWake { queue, producer }
                    })
                    .collect();
                let gated = forward.iter().any(|b| b.gate.is_some());
                let needs_fence = !reverse.is_empty()
                    || forward
                        .iter()
                        .any(|b| matches!(b.target, WakeTarget::Driver(_) | WakeTarget::Worker(_)));
                let own = forward.iter().any(|b| b.same_thread).then(|| thread_of(step));
                let own = own.unwrap_or(WakeTarget::None);
                StepWake { forward, reverse, gated, needs_fence, own }
            })
            .collect();
        tracked.sort_by_key(|&(s, b)| (s.0, b.0));
        tracked.dedup();
        let caps: Box<[Option<Arc<PhaseCap>>]> =
            (0..n).map(|i| caps.get(i).cloned().flatten()).collect();
        let distinct_caps = index_caps(&caps);

        Arc::new(Self {
            steps,
            drivers: ThreadSlots::new(n_drivers),
            workers: ThreadSlots::new(n_threads),
            direct_parked: DirectParked::new(n_threads),
            cap_parked: DirectParked::new(n_threads),
            parking: DirectParked::new(n_threads + n_drivers),
            gates: n_threads + n_drivers > 1,
            caps,
            distinct_caps,
            pool,
            tracked: tracked.into_boxed_slice(),
            request_pending: Padded(AtomicBool::new(false)),
        })
    }

    /// The implicit forward branch of a `Parallel` step with no output edge.
    /// Such a step publishes through a side channel (the sort's slot FIFOs);
    /// its `Progress` wakes one pool worker (a `notify_one` cascade) so a peer clone
    /// can claim the work it exposed. Ungated — there is no transport to gate
    /// on — and built only with a peer to wake: at one thread the lone worker
    /// is the one making the progress. A `Pool` target needs no fence. The
    /// woken peer runs `step` itself, so `step` is the consumer whose phase cap
    /// gates the cap-parked fallback.
    fn self_cascade(step: StepIdx, fallback_mask: Option<Box<[u64]>>) -> BranchWake {
        BranchWake {
            target: WakeTarget::Pool,
            gate: None,
            same_thread: false,
            consumer: Some(step),
            fallback_mask,
        }
    }

    /// The pool event-count (`None` at one thread).
    #[must_use]
    pub fn pool(&self) -> Option<&PoolEventCount> {
        self.pool.as_deref()
    }

    /// Number of dedicated driver threads the plan routes to.
    #[cfg(test)]
    #[must_use]
    pub fn n_drivers(&self) -> usize {
        self.drivers.len()
    }

    /// Wake slots: pool workers `0..n_threads`, then drivers.
    #[cfg(test)]
    #[must_use]
    pub fn n_slots(&self) -> usize {
        self.workers.len() + self.drivers.len()
    }

    /// Forward branches in the built table for `step` (test support): a
    /// zero-output Parallel step's implicit cascade counts as one.
    #[cfg(test)]
    pub(crate) fn forward_branch_count_for_test(&self, step: StepIdx) -> usize {
        self.steps.get(step.0).map_or(0, |s| s.forward.len())
    }

    /// The forward target of `(step, branch)`; `None` for a step or branch the
    /// plan does not know (e.g. any step of the standalone plan).
    #[cfg(test)]
    #[must_use]
    pub fn forward_target(&self, step: StepIdx, branch: BranchIdx) -> WakeTarget {
        self.steps
            .get(step.0)
            .and_then(|s| s.forward.get(branch.0))
            .map_or(WakeTarget::None, |b| b.target)
    }

    /// Every pool worker some forward branch wakes with a direct `unpark`
    /// (`WakeTarget::Worker`). `Pipeline::run` checks each idles on its timer,
    /// where that unpark reaches it.
    pub fn direct_worker_targets(&self) -> impl Iterator<Item = usize> + '_ {
        self.steps.iter().flat_map(|s| s.forward.iter()).filter_map(|b| match b.target {
            WakeTarget::Worker(w) => Some(w),
            WakeTarget::None | WakeTarget::Pool | WakeTarget::Driver(_) => None,
        })
    }

    /// Whether any branch of `step` is push-gated (so `dispatch_one_step` must snapshot).
    #[inline]
    #[must_use]
    pub fn is_gated(&self, step: StepIdx) -> bool {
        self.steps.get(step.0).is_some_and(|s| s.gated)
    }

    /// How many branches of `step` are push-gated.
    #[cfg(test)]
    #[must_use]
    pub fn gated_branches(&self, step: StepIdx) -> usize {
        self.steps.get(step.0).map_or(0, |s| s.forward.iter().filter(|b| b.gate.is_some()).count())
    }

    /// The producers of `step`'s reverse-wake edges: byte-bounded inputs first,
    /// then holder-only, each in registry order.
    #[cfg(test)]
    #[must_use]
    pub fn reverse_producers(&self, step: StepIdx) -> Vec<StepIdx> {
        self.steps
            .get(step.0)
            .map_or_else(Vec::new, |s| s.reverse.iter().map(|r| r.producer).collect())
    }

    /// Whether `Pipeline::run` must install a holder set (and push counting) on `q`.
    #[must_use]
    pub fn tracks(&self, q: &RegisteredQueue) -> bool {
        self.is_tracked_key(q.producer_step, q.branch)
    }

    /// Whether `Pipeline::run` must install a holder set on the holder-only edge `q`.
    #[must_use]
    pub(crate) fn tracks_holder_only(&self, q: &RegisteredHolderOnlyQueue) -> bool {
        self.is_tracked_key(q.producer_step, q.branch)
    }

    fn is_tracked_key(&self, step: StepIdx, branch: BranchIdx) -> bool {
        self.tracked.binary_search_by_key(&(step.0, branch.0), |&(s, b)| (s.0, b.0)).is_ok()
    }

    /// Whether every reverse edge's queue has a holder set installed. A reverse
    /// edge on an untracked queue would take nothing forever and strand its
    /// holder on its timer; `Pipeline::run` asserts this right after it
    /// installs tracking.
    #[must_use]
    pub fn reverse_edges_are_tracked(&self) -> bool {
        self.steps.iter().all(|s| s.reverse.iter().all(|r| r.queue.is_tracked()))
    }

    /// Whether every gated branch's queue has a holder set installed (the push
    /// counter advances only then). A gate on an untracked queue would read a
    /// frozen counter, hold every forward wake off and strand the consumer on
    /// its timer; `Pipeline::run` asserts this right after it installs
    /// tracking.
    #[must_use]
    pub fn gated_edges_are_tracked(&self) -> bool {
        self.steps
            .iter()
            .all(|s| s.forward.iter().filter_map(|b| b.gate.as_ref()).all(|q| q.is_tracked(SEALED)))
    }

    /// Register driver `d`'s OS thread. Called by the driver thread itself as
    /// its first statement, before its first pass. The registration's `SeqCst`
    /// fence pairs with the delivery fence in `on_progress`: a producer's push
    /// either finds this handle or is found by the first pass
    /// ([`ThreadSlots`]). (A registration made by another thread after `spawn`
    /// has no such ordering against the driver's first pass and parks.)
    pub fn register_driver(&self, d: DriverIdx, t: ParkThread) {
        self.drivers.register(d.0, t);
    }

    /// Register pool worker `w`'s OS thread — called by that thread itself before
    /// its first pass, for the reason given on [`Self::register_driver`]; the
    /// one-worker path registers the caller thread, which runs worker 0. Every
    /// worker of every chain registers, pinned or not: a holder idles on its
    /// timer and must be unparkable. Ignored for an out-of-range `w` (and so on
    /// the standalone plan). A registered thread can receive an `unpark` it
    /// does not need; the std park token then makes a later `park_timeout`
    /// return early, which every park site tolerates by re-checking its
    /// condition.
    pub fn register_worker(&self, w: usize, t: ParkThread) {
        self.workers.register(w, t);
    }

    /// Whether pool worker `w`'s thread is registered.
    #[cfg(test)]
    #[must_use]
    pub fn is_worker_registered(&self, w: usize) -> bool {
        self.workers.is_registered(w)
    }

    /// Snapshot the gated branches' push counters before a dispatch (a counter
    /// advances on every push and on every reorder-stash insert). One `Relaxed`
    /// load per gated branch; call only when `is_gated`.
    #[inline]
    pub fn snapshot_pushes(&self, step: StepIdx, out: &mut PushSnapshot) {
        if let Some(s) = self.steps.get(step.0) {
            for (slot, b) in out.0.iter_mut().zip(s.forward.iter()) {
                if let Some(g) = &b.gate {
                    *slot = g.pushed_total(SEALED);
                }
            }
        }
    }

    /// Deliver the wakes a `Progress` of `step` implies: each forward branch's
    /// target (skipped on a gated branch that pushed nothing since `before`),
    /// then the holders recorded on every reverse edge.
    #[inline]
    #[must_use]
    pub(crate) fn on_progress(&self, step: StepIdx, before: &PushSnapshot) -> WakeCounts {
        let mut r = WakeCounts::default();
        let Some(s) = self.steps.get(step.0) else {
            // A step the plan does not know (the standalone plan): nothing to wake.
            return r;
        };
        let _ = self.forward(s, before, &mut r);
        for rev in &s.reverse {
            let n = rev.queue.take_holders(&mut |slot| {
                let _ = self.deliver_holder(slot);
            });
            r.reverse = r.reverse.saturating_add(u8::try_from(n).unwrap_or(u8::MAX));
        }
        r
    }

    /// [`Self::on_progress`] with no gated branch, for the loom models in
    /// `tests/loom_wake.rs`, which drive the real delivery path.
    #[cfg(loom)]
    #[doc(hidden)]
    pub fn on_progress_ungated(&self, step: StepIdx) {
        let _ = self.on_progress(step, &PushSnapshot::default());
    }

    /// A dispatch of `step` flushed a held item (`BranchOutputHandle::retry`
    /// returned `Ok`) and reported `NoProgress`, `Contention` or `Capped`.
    /// Deliver the forward wakes that push implies, exactly as `on_progress`
    /// does (the same fence, gates and targets: one shared body). Then, if a
    /// branch whose consumer runs on this thread passed its gate, unpark this
    /// thread itself: without a `Progress` the thread may idle without another
    /// pass, and the park token makes its next timer park return at once. Every
    /// thread that can host a same-thread consumer (a driver, a pinned worker,
    /// the lone worker) idles only through `park_timeout`. No reverse take: a
    /// dispatch that popped reports `Progress` (the `StepOutcome::Progress`
    /// contract), so this one made no room. A step the plan does not know (the
    /// standalone plan): an empty report, nothing delivered.
    #[inline]
    #[must_use]
    pub(crate) fn on_flushed(&self, step: StepIdx, before: &PushSnapshot) -> WakeCounts {
        let mut r = WakeCounts::default();
        let Some(s) = self.steps.get(step.0) else {
            return r; // the standalone plan: no per-step table
        };
        if self.forward(s, before, &mut r) {
            // `Driver`/`Worker`: an unpark of this very thread, which is
            // running, so its parking bit is clear — unconditional.
            let woke = match s.own {
                WakeTarget::Worker(w) => self.workers.unpark(w),
                WakeTarget::Driver(d) => self.drivers.unpark(d.0),
                WakeTarget::None | WakeTarget::Pool => false,
            };
            if woke {
                r.unparked = r.unparked.saturating_add(1);
            }
        }
        r
    }

    /// The forward half of `on_progress` and `on_flushed`: the delivery fence
    /// when an unpark or a reverse take follows, then each forward branch's
    /// target unless its gate saw no push since `before`. Returns whether a
    /// delivered branch's consumer runs on this thread (`on_flushed` re-polls
    /// it).
    #[inline]
    fn forward(&self, s: &StepWake, before: &PushSnapshot, r: &mut WakeCounts) -> bool {
        if s.needs_fence {
            // Producer half of the unpark-path fence rule and the consumer half
            // of the holder protocol (`wake_slot::delivery_fence`).
            delivery_fence();
        }
        let mut same_thread = false;
        for (b, wake) in s.forward.iter().enumerate() {
            if let Some(g) = &wake.gate
                && g.pushed_total(SEALED) == before.0[b]
            {
                r.gated_off = r.gated_off.saturating_add(1);
                continue;
            }
            same_thread |= wake.same_thread;
            self.deliver(wake.target, wake.consumer, wake.fallback_mask.as_deref(), r);
        }
        same_thread
    }

    /// A capacity release that is not a pop (the queue-memory rebalancer raised
    /// `q`'s limit): claim `q`'s recorded holders and deliver each through
    /// [`Self::deliver_holder`] (an `unpark` only for one armed for its timer
    /// park, or the lone thread of a one-thread plan; a running holder
    /// re-polls on its own). Fences first. Returns the records claimed.
    pub fn wake_holders(&self, q: &dyn BoundedQueueHandle) -> usize {
        delivery_fence();
        q.take_holders(
            &mut |slot| {
                let _ = self.deliver_holder(slot);
            },
            SEALED,
        )
    }

    /// A `Finished` closed an edge: broadcast to the pool and unpark every
    /// registered thread (~one call per step per run).
    pub fn on_finished(&self, _step: StepIdx) {
        if let Some(p) = &self.pool {
            let _ = p.notify_all();
        }
        self.unpark_all();
    }

    /// Unpark every registered driver and worker. Used by `on_finished` and by
    /// the cancel/error path (`PipelineSignal`) so a parked thread observes
    /// `is_done()` without waiting out its timer.
    pub fn unpark_all(&self) {
        self.drivers.unpark_all();
        self.workers.unpark_all();
    }

    /// Worker `w` is about to park on its timer: set its parking bit, so a
    /// directed wake reaches it, and — when an event-count exists — its
    /// direct-parked bit, so a `Pool` wake that finds no event-count waiter
    /// does too; then fence once. Returns `true` when armed: the caller must
    /// then re-poll once before parking (the fence orders the arm before that
    /// re-poll), and that one re-poll serves both bits. `false` in a plan with
    /// one thread and no event-count (the lone worker of a run with no driver,
    /// or the standalone plan): such a plan does not gate its directed wakes
    /// (they unpark the lone worker unconditionally), so there is nothing to
    /// arm and no extra pass.
    #[must_use]
    pub fn arm_direct(&self, w: usize) -> bool {
        let parking = self.mark_parking(w);
        let direct = self.pool.is_some() && self.direct_parked.mark(w);
        Self::fence_arm(parking || direct)
    }

    /// Worker `w` is running again.
    pub fn disarm_direct(&self, w: usize) {
        self.direct_parked.disarm(w);
        self.parking.disarm(w);
    }

    /// Driver `d` is about to park on its timer: set its parking bit and
    /// fence. Same protocol as [`Self::arm_direct`]: `true` means armed, and
    /// the caller re-polls once before parking. `false` only in a plan with
    /// one thread (the standalone plan of a test driving a driver loop by
    /// hand), which gates nothing.
    ///
    /// # Panics
    ///
    /// If the plan gates its wakes and has no driver `d`: such a driver would
    /// park with no bit of its own, so every directed wake would skip it.
    #[must_use]
    pub fn arm_driver(&self, d: DriverIdx) -> bool {
        if !self.gates {
            return false;
        }
        assert!(
            d.0 < self.drivers.len(),
            "arm_driver: driver {} is not in this plan ({} drivers)",
            d.0,
            self.drivers.len()
        );
        let parking = self.mark_parking(self.workers.len() + d.0);
        Self::fence_arm(parking)
    }

    /// Driver `d` is running again.
    pub fn disarm_driver(&self, d: DriverIdx) {
        self.parking.disarm(self.workers.len() + d.0);
    }

    /// Set wake slot `slot`'s parking bit, unless the plan has one thread
    /// (which does not gate on the bit).
    #[inline]
    fn mark_parking(&self, slot: usize) -> bool {
        self.gates && self.parking.mark(slot)
    }

    /// The one fence of a park arm, when it set any bit; returns `armed`.
    #[inline]
    fn fence_arm(armed: bool) -> bool {
        if armed {
            arm_fence();
        }
        armed
    }

    /// Whether wake slot `slot` is armed for a timer park (tests).
    #[cfg(test)]
    pub(crate) fn is_parking(&self, slot: usize) -> bool {
        self.parking.is_set(slot)
    }

    /// [`Self::arm_direct`]'s parking bit without its fence: the negative
    /// control of the loom parking model, which must lose a wake.
    #[cfg(loom)]
    #[doc(hidden)]
    pub fn mark_parking_unfenced(&self, w: usize) {
        let _ = self.mark_parking(w);
    }

    /// Worker `w`, refused by a phase cap and holding nothing, is about to park
    /// on its timer: make it reachable by a `Pool` wake that finds no
    /// event-count waiter and no direct-parked worker, for a consumer that is
    /// uncapped or whose cap has a free permit.
    /// Same protocol as [`Self::arm_direct`]: `true` means armed, and the
    /// caller re-polls once before parking.
    #[must_use]
    pub fn arm_cap_parked(&self, w: usize) -> bool {
        let parking = self.mark_parking(w);
        let cap = self.pool.is_some() && self.cap_parked.mark(w);
        Self::fence_arm(parking || cap)
    }

    /// The phase cap `step` reports, if any (`None` in the standalone plan).
    #[inline]
    fn cap_of(&self, step: StepIdx) -> Option<&PhaseCap> {
        self.caps.get(step.0).and_then(Option::as_deref)
    }

    /// Whether a worker woken for `consumer` could be admitted there: it is
    /// uncapped (or there is no consumer), or its cap has a free permit.
    #[inline]
    fn admits(&self, consumer: Option<StepIdx>) -> bool {
        consumer.and_then(|c| self.cap_of(c)).is_none_or(PhaseCap::has_free_permit)
    }

    /// A cap-parked worker's re-poll is about to visit `step`: whether to skip
    /// it because its phase cap is full right now (no free permit, or a
    /// whole-cap reservation pending). A skip is recorded on that cap as a
    /// refusal would be (`PhaseCap::record_skip`: record, fence, re-check), so
    /// the release that frees a permit wakes this thread; a re-check that sees
    /// a free permit polls the step instead. `false` for an uncapped step and
    /// in the standalone plan.
    #[inline]
    #[must_use]
    pub(crate) fn skip_full_cap(&self, step: StepIdx) -> bool {
        self.cap_of(step).is_some_and(PhaseCap::record_skip_current)
    }

    /// The calling thread starts a dispatch pass: clear its own refusal
    /// record on every cap, so a record covers only the pass that made it
    /// (`PhaseCap::try_acquire`). Per cap one `Relaxed` load, plus a
    /// `fetch_and` only when the bit is set; nothing without caps.
    ///
    /// A cap the thread recorded on since its last pass start, whose bit is
    /// now clear, had the record claimed by a release: the thread was handed
    /// that cap's wake. It owes the wake until it takes a permit there in
    /// this pass; [`Self::end_pass`] forwards it otherwise.
    #[inline]
    pub(crate) fn start_pass(&self) {
        if self.distinct_caps.is_empty() {
            return;
        }
        let slot = current_slot();
        let recorded = take_cap_recorded();
        let mut owed = 0u64;
        for (i, cap) in self.distinct_caps.iter().enumerate() {
            let still_set = cap.clear_own_record(slot);
            if recorded & (1 << i) != 0 && !still_set {
                owed |= 1 << i;
            }
        }
        set_cap_owed(owed);
    }

    /// The calling thread ends a dispatch pass: forward every cap wake it
    /// owes and did not spend (`PhaseCap::forward_wake`: if a permit is
    /// still free, wake the next recorded thread). Returns at once in a plan
    /// without caps, and after one thread-local read when nothing is owed.
    #[inline]
    pub(crate) fn end_pass(&self) {
        if self.distinct_caps.is_empty() {
            return;
        }
        let mut owed = take_cap_owed();
        while owed != 0 {
            let i = owed.trailing_zeros() as usize;
            owed &= owed - 1;
            if let Some(cap) = self.distinct_caps.get(i) {
                cap.forward_wake();
            }
        }
    }

    /// Whether worker `w` is armed as cap-parked (tests).
    #[cfg(test)]
    pub(crate) fn is_cap_parked(&self, w: usize) -> bool {
        self.cap_parked.is_armed(w)
    }

    /// Worker `w` is running again.
    pub fn disarm_cap_parked(&self, w: usize) {
        self.cap_parked.disarm(w);
        self.parking.disarm(w);
    }

    /// Deliver one forward wake. `fallback_mask` restricts a `Pool` wake's
    /// direct-park and cap-parked fallbacks to the consumer's host workers
    /// (`None`: any worker).
    ///
    /// Only a `Pool` target goes through the event-count, whose concurrency
    /// ceiling may suppress it; reverse wakes are unparks of the holder's slot
    /// (`deliver_slot`), never `notify_one`, so the ceiling cannot suppress a
    /// wake meant for a specific holder.
    #[inline]
    fn deliver(
        &self,
        target: WakeTarget,
        consumer: Option<StepIdx>,
        fallback_mask: Option<&[u64]>,
        r: &mut WakeCounts,
    ) {
        match target {
            WakeTarget::None => {}
            WakeTarget::Pool => {
                if let Some(p) = &self.pool {
                    let o = p.notify_one();
                    Self::count_notify(o, r);
                    // `notify_one` fenced before reading `waiters`; the same fence
                    // orders these reads of the armed bits and of the consumer
                    // cap's occupancy. A wake for a consumer whose cap is full
                    // never takes a cap-parked worker: it would be refused
                    // again, and the release that frees a permit wakes the
                    // threads recorded on that cap. The cap's lines are read
                    // only when some worker is cap-parked.
                    if o == NotifyOutcome::NoWaiters
                        && (self.unpark_one_direct_parked(fallback_mask)
                            || (!self.cap_parked.is_empty()
                                && self.admits(consumer)
                                && self.unpark_one_cap_parked(fallback_mask)))
                    {
                        r.fallback = r.fallback.saturating_add(1);
                    }
                }
            }
            WakeTarget::Driver(_) | WakeTarget::Worker(_) => {
                if self.unpark_if_parking(target) {
                    r.unparked = r.unparked.saturating_add(1);
                }
            }
        }
    }

    /// Unpark `target` (a `Driver` or `Worker`) only if it is armed for a
    /// timer park (`parking`); `false` when it is running or unregistered.
    /// The caller has fenced after its publish: a thread whose bit is clear
    /// here sets it, fences and re-polls before it parks, and that re-poll
    /// sees the publish. A plan with one thread does not gate
    /// ([`Self::gates`]): its lone thread is unparked whenever registered.
    #[inline]
    fn unpark_if_parking(&self, target: WakeTarget) -> bool {
        let armed = |slot: usize| !self.gates || self.parking.is_set(slot);
        match target {
            WakeTarget::Worker(w) => armed(w) && self.workers.unpark(w),
            WakeTarget::Driver(d) => armed(self.workers.len() + d.0) && self.drivers.unpark(d.0),
            WakeTarget::None | WakeTarget::Pool => false,
        }
    }

    /// Unpark the holder in wake slot `slot` (worker, driver, or — for a thread
    /// without a slot — everyone). A worker or driver is unparked only if it
    /// is armed for a timer park (in a plan with one thread, whenever it is
    /// registered): a running holder's next pass, or the re-poll it runs
    /// before it parks, retries on its own. The caller has fenced. Returns
    /// whether an `unpark` was issued (`true` for the broadcast).
    #[inline]
    pub fn deliver_slot(&self, slot: usize) -> bool {
        let nw = self.workers.len();
        if slot < nw {
            self.unpark_if_parking(WakeTarget::Worker(slot))
        } else if slot - nw < self.drivers.len() {
            self.unpark_if_parking(WakeTarget::Driver(DriverIdx(slot - nw)))
        } else {
            self.deliver_anonymous();
            true
        }
    }

    /// Deliver a holder wake for `slot`, whose record the caller has just
    /// claimed from an edge's holder set (after its own fence): fence again,
    /// then [`Self::deliver_slot`].
    ///
    /// The second fence is what makes the parking gate sound for holder
    /// claims. A claim need not come with room for the holder (an in-order pop
    /// that leaves the next ordinal buffered, a pop that frees too little), and
    /// it can claim a record the holder wrote after the caller's first fence —
    /// one its re-poll made after it armed. Reading the parking bit unfenced
    /// there can miss that arm: the holder is skipped, its record is gone, and
    /// its re-check saw no room, so it parks with nothing left to wake it.
    /// With the fence, either the claimed record was written after the
    /// holder's arm (the arm's fence then orders the bit before this read, and
    /// the holder is unparked) or every record its re-poll writes after the
    /// arm lands after this claim and stands for the next take. A cap release
    /// needs no such fence: it always frees a permit, and a thread whose
    /// claimed cap record went unwoken owes that wake and forwards it
    /// (`WakePlan::start_pass`). `tests/loom_wake.rs`'s ordered fast-path
    /// model loses the wake without it.
    #[inline]
    pub fn deliver_holder(&self, slot: usize) -> bool {
        delivery_fence();
        self.deliver_slot(slot)
    }

    /// Wake a holder whose slot is unknown: every thread, pool waiters
    /// included.
    pub(crate) fn deliver_anonymous(&self) {
        if let Some(p) = &self.pool {
            let _ = p.notify_all();
        }
        self.unpark_all();
    }

    /// Claim one armed worker — within `mask` when given — and unpark it,
    /// trying the next armed bit when a claimed one has no registered thread.
    fn unpark_one_direct_parked(&self, mask: Option<&[u64]>) -> bool {
        self.direct_parked.claim_one(mask, &mut |w| self.workers.unpark(w))
    }

    /// Claim one cap-parked worker — within `mask` when given — and unpark it
    /// (see [`Self::arm_cap_parked`]).
    fn unpark_one_cap_parked(&self, mask: Option<&[u64]>) -> bool {
        self.cap_parked.claim_one(mask, &mut |w| self.workers.unpark(w))
    }

    #[inline]
    fn count_notify(o: NotifyOutcome, r: &mut WakeCounts) {
        match o {
            NotifyOutcome::Woken => r.notified = r.notified.saturating_add(1),
            NotifyOutcome::NoWaiters => r.no_waiters = r.no_waiters.saturating_add(1),
            NotifyOutcome::Suppressed => r.suppressed = r.suppressed.saturating_add(1),
        }
    }

    /// Semantics of [`PoolHandle::request_worker`]: at one thread, unpark the
    /// lone worker (unless it is the caller); otherwise wake one parked worker
    /// other than the caller: an event-count waiter (ungated by the ceiling;
    /// `Pending` while an earlier such request is unacknowledged) or, when
    /// there is none, one direct-parked worker, else one cap-parked worker
    /// (see [`PoolRequest::Woken`] for when that wake is wasted). `AllAwake`
    /// only when none exists.
    #[inline]
    #[must_use]
    pub(crate) fn request_pool_worker(&self) -> PoolRequest {
        let caller = crate::runtime::wake_slot::current_slot();
        let Some(pool) = &self.pool else {
            // One-thread run: the lone worker is `Worker(0)`, and a request
            // from worker 0 itself has no one else to wake. The standalone plan
            // has no worker at all: `Unavailable`.
            if caller == Some(0) && self.workers.len() > 0 {
                return PoolRequest::AllAwake;
            }
            return if self.workers.unpark(0) {
                PoolRequest::Woken
            } else {
                PoolRequest::Unavailable
            };
        };
        // Producer half of both fence pairs: it orders the caller's publish (the
        // off-edge work the requester is waiting on) before the census below and
        // before the armed-bit read, pairing with `prepare_wait`'s fence and with
        // `arm_direct`'s.
        delivery_fence();
        // A caller armed on the event-count (a step in its re-poll) is one of
        // the waiters; it must not count as the worker to wake.
        let own = usize::from(crate::runtime::wake_slot::ec_armed());
        if pool.waiters() <= own {
            return self.request_direct_parked(caller);
        }
        if self.request_pending.0.swap(true, Ordering::AcqRel) {
            return PoolRequest::Pending;
        }
        match pool.notify_one_demand() {
            NotifyOutcome::Woken => PoolRequest::Woken,
            // The waiter left between the census and the notify: nothing
            // pending, and a timer-parked worker may still be there.
            NotifyOutcome::NoWaiters | NotifyOutcome::Suppressed => {
                self.request_pending.0.store(false, Ordering::Release);
                self.request_direct_parked(caller)
            }
        }
    }

    /// No other event-count waiter: unpark one worker that armed for a timer
    /// park (a pinned worker, a holder), else one parked only because a phase
    /// cap refused it, as a `Pool` forward wake for an uncapped consumer does
    /// (the work requested comes off any queue edge, so no consumer cap is
    /// known; a cap-parked worker woken for work behind a full cap is refused
    /// again, a wasted wake), skipping the caller's own bit. The claim clears that worker's
    /// bit, so a repeated request cannot wake it twice; it is not a pending
    /// event-count request (that worker never returns from an event-count
    /// wait), so `request_pending` is not set. The caller fenced after its
    /// publish.
    #[inline]
    fn request_direct_parked(&self, caller: Option<usize>) -> PoolRequest {
        let unpark = &mut |w| self.workers.unpark(w);
        let woke = self.direct_parked.claim_one_except(caller, unpark)
            || self.cap_parked.claim_one_except(caller, unpark);
        if woke { PoolRequest::Woken } else { PoolRequest::AllAwake }
    }

    /// Worker side: called after **every** exit from the event-count (a return
    /// from its wait, a productive arm→wait re-poll, or a re-poll that was
    /// refused and leaves for the timer park), so one request is outstanding
    /// until some worker wakes for any reason — bounding requests at ≤ one per
    /// worker wake. A `Relaxed` load, plus a `Release` store only when a
    /// request is pending, so waking workers do not all write the line when
    /// no request exists.
    #[inline]
    pub(crate) fn acknowledge_request(&self) {
        if self.request_pending.0.load(Ordering::Relaxed) {
            self.request_pending.0.store(false, Ordering::Release);
        }
    }

    /// One `wake:` line per step for `Pipeline::dag()`, so a chain's routing
    /// can be reviewed without running it.
    #[must_use]
    pub fn dag_lines(&self, graph: &ChainGraph) -> Vec<String> {
        let n = graph.n_steps();
        let fmt_target = |t: WakeTarget| match t {
            WakeTarget::None => "None".to_string(),
            WakeTarget::Pool => "Pool".to_string(),
            WakeTarget::Driver(d) => format!("Driver({})", d.0),
            WakeTarget::Worker(w) => format!("Worker({w})"),
        };
        (0..n)
            .map(|i| {
                let step = StepIdx(i);
                let s = &self.steps[i];
                let mut line = String::from("      wake: ");
                if s.forward.is_empty() {
                    line.push_str("(sink)");
                } else if graph.branch_count(step) == 0 {
                    // A zero-output Parallel step's implicit self-cascade.
                    let parts: Vec<String> =
                        s.forward.iter().map(|w| fmt_target(w.target)).collect();
                    line.push_str("(sink) ");
                    line.push_str(&parts.join(", "));
                } else {
                    let parts: Vec<String> = s
                        .forward
                        .iter()
                        .enumerate()
                        .map(|(b, w)| {
                            let consumer = graph
                                .consumer(step, BranchIdx(b))
                                .map_or("<unwired>", |c| graph.step_name(c));
                            let gated = if w.gate.is_some() { " [gated]" } else { "" };
                            format!(".{b} → {consumer} {}{gated}", fmt_target(w.target))
                        })
                        .collect();
                    line.push_str(&parts.join(", "));
                }
                if !s.reverse.is_empty() {
                    line.push_str("; reverse: ");
                    let producers: Vec<String> = s
                        .reverse
                        .iter()
                        .map(|r| format!("{} (holder)", graph.step_name(r.producer)))
                        .collect();
                    line.push_str(&producers.join(", "));
                }
                line
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use std::time::Duration;

    use rstest::rstest;

    use super::*;
    use crate::queues::{BoundedQueueHandle, ByteBoundedQueue, EdgeTracking, ItemQueue};
    use crate::runtime::contexts::{RegisteredHolderOnlyQueue, RegisteredQueue};
    use crate::runtime::event_count::{PoolEventCount, WaitOutcome};
    use crate::runtime::wake_slot::{HolderSet, SlotGuard};
    use crate::step::StepKind;
    use crate::topology::{BranchIdx, ChainGraph, StepIdx};

    /// 8-byte items so a 1-byte limit admits exactly one.
    #[derive(Debug, PartialEq)]
    struct Sized8;
    impl crate::item::HeapSize for Sized8 {
        fn heap_size(&self) -> usize {
            8
        }
    }

    fn q(producer: usize, branch: usize) -> RegisteredQueue {
        let handle: Arc<dyn BoundedQueueHandle> = Arc::new(ByteBoundedQueue::<u32>::new(1 << 20));
        RegisteredQueue {
            producer_step_name: "q",
            producer_step: StepIdx(producer),
            branch: BranchIdx(branch),
            handle,
            reorder_cap: None,
        }
    }

    /// A one-item tracked queue plus its typed handle, for reverse-wake units.
    fn one_item_tracked(
        producer: usize,
        n_slots: usize,
    ) -> (RegisteredQueue, Arc<ByteBoundedQueue<Sized8>>) {
        let typed = Arc::new(ByteBoundedQueue::<Sized8>::new(1));
        typed.enable_tracking(
            EdgeTracking { clock: None, holders: Some(HolderSet::new(n_slots)) },
            crate::queues::SEALED,
        );
        let handle: Arc<dyn BoundedQueueHandle> = Arc::clone(&typed) as _;
        let rq = RegisteredQueue {
            producer_step_name: "q",
            producer_step: StepIdx(producer),
            branch: BranchIdx(0),
            handle,
            reorder_cap: None,
        };
        (rq, typed)
    }

    /// The two transports that can refuse a direct push.
    #[derive(Clone, Copy, Debug)]
    enum Bounded {
        Byte,
        Count,
    }

    /// A one-item edge `(producer, 0)` with a holder set installed, as the
    /// registry entry the plan reads (byte-bounded or holder-only) plus the
    /// transport.
    struct OneItemEdge {
        byte: Vec<RegisteredQueue>,
        holder_only: Vec<RegisteredHolderOnlyQueue>,
        queue: Arc<dyn ItemQueue<Sized8>>,
    }

    impl OneItemEdge {
        fn tracked(kind: Bounded, producer: usize, n_slots: usize) -> Self {
            use crate::queues::{CountBoundedQueue, HolderOnlyHandle};
            match kind {
                Bounded::Byte => {
                    let (rq, typed) = one_item_tracked(producer, n_slots);
                    Self { byte: vec![rq], holder_only: vec![], queue: typed }
                }
                Bounded::Count => {
                    let typed = Arc::new(CountBoundedQueue::<Sized8>::new(1));
                    typed.enable_holder_tracking(n_slots);
                    let rq = RegisteredHolderOnlyQueue {
                        producer_step: StepIdx(producer),
                        branch: BranchIdx(0),
                        handle: Arc::clone(&typed) as _,
                    };
                    Self { byte: vec![], holder_only: vec![rq], queue: typed }
                }
            }
        }

        fn edges(&self) -> WakeEdges<'_> {
            WakeEdges { byte_bounded: &self.byte, holder_only: &self.holder_only }
        }
    }

    /// A holder-only registry entry for `(producer, 0)` over a fresh
    /// one-item count-bounded queue (untracked until the plan says so).
    fn count_edge(producer: usize) -> RegisteredHolderOnlyQueue {
        RegisteredHolderOnlyQueue {
            producer_step: StepIdx(producer),
            branch: BranchIdx(0),
            handle: Arc::new(crate::queues::CountBoundedQueue::<u32>::new(1)),
        }
    }

    type Topology = (
        ChainGraph,
        Vec<StepKind>,
        Vec<Option<usize>>,
        Vec<Option<DriverIdx>>,
        Vec<RegisteredQueue>,
    );

    /// A synthetic 11-step topology shaped like an earlier sort chain
    /// (template-coordinate, t16, `ReadBgzfBlocks` front, no merge-supply
    /// planner): two Shared driver groups (coord = {1,3,4,8}, io = {6,10}),
    /// `ReadBgzfBlocks` pinned to w0, every edge byte-bounded. It mirrors no
    /// production chain today and nothing checks it against one; it is kept
    /// as a mixed-kind fixture whose expectations are its own. One Detached
    /// `ArenaSort` step stands for the arena front's serial sort work on the
    /// coordination driver.
    /// [`native_input_sort_topology`] models the current chain.
    pub(super) fn synthetic_sort_topology() -> Topology {
        use StepKind::{Detached as D, Parallel as P, Serial as S};
        let names = [
            "ReadBgzfBlocks",
            "ReadBlocks",
            "InflateToArena",
            "ArenaSort",
            "SpillGather",
            "SpillBlockCompress",
            "SpillWrite",
            "SortSpillDecompress",
            "SortMerge",
            "BgzfCompress",
            "WriteBgzfFile",
        ];
        let kinds = vec![S, D, P, D, D, P, D, P, D, P, D];
        let mut g = ChainGraph::new();
        let idx: Vec<StepIdx> = names
            .iter()
            .enumerate()
            .map(|(i, n)| g.register_step(n, usize::from(i != 10)))
            .collect();
        for i in 0..10 {
            g.wire(idx[i], BranchIdx(0), idx[i + 1]);
        }
        let coord = Some(DriverIdx(0));
        let io = Some(DriverIdx(1));
        let driver_of = vec![None, coord, None, coord, coord, None, io, None, coord, None, io];
        let pinned = vec![Some(0), None, None, None, None, None, None, None, None, None, None];
        let queues = (0..10).map(|i| q(i, 0)).collect();
        (g, kinds, pinned, driver_of, queues)
    }

    /// `Src(Serial w0) → Det(Detached d0) → Par(Parallel, no output edge)`.
    fn three_step_with_parallel_sink() -> Topology {
        let mut g = ChainGraph::new();
        let src = g.register_step("Src", 1);
        let det = g.register_step("Det", 1);
        let par = g.register_step("Par", 0);
        g.wire(src, BranchIdx(0), det);
        g.wire(det, BranchIdx(0), par);
        let kinds = vec![StepKind::Serial, StepKind::Detached, StepKind::Parallel];
        let pinned = vec![Some(0), None, None];
        let driver_of = vec![None, Some(DriverIdx(0)), None];
        let queues = (0..2).map(|i| q(i, 0)).collect();
        (g, kinds, pinned, driver_of, queues)
    }

    /// A zero-output Parallel step gets one implicit ungated `Pool` target (read
    /// from the built table, not from `forward_target`'s `Pool` default), and its
    /// `Progress` reaches a real armed waiter.
    #[test]
    fn zero_output_parallel_step_cascades_to_the_pool() {
        let (g, kinds, pinned, driver_of, queues) = three_step_with_parallel_sink();
        let pool = Arc::new(PoolEventCount::new(4));
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            Some(Arc::clone(&pool)),
            4,
        );
        assert_eq!(plan.forward_branch_count_for_test(StepIdx(2)), 1, "one implicit branch");
        assert!(!plan.is_gated(StepIdx(2)), "the implicit branch is ungated");
        assert!(plan.dag_lines(&g)[2].contains("(sink) Pool"), "{:?}", plan.dag_lines(&g));
        let (tx, rx) = std::sync::mpsc::channel();
        let p2 = Arc::clone(&pool);
        let waiter = std::thread::spawn(move || {
            let key = p2.prepare_wait();
            tx.send(()).unwrap();
            p2.wait(key, Duration::from_secs(5))
        });
        rx.recv().unwrap();
        let _ = plan.on_progress(StepIdx(2), &PushSnapshot::default());
        assert!(
            matches!(waiter.join().unwrap(), WaitOutcome::Woken | WaitOutcome::Notified),
            "the Progress must wake the armed waiter before its 5 s timer"
        );
    }

    /// Only Parallel sinks get the implicit target: a Serial/Exclusive sink is
    /// woken by its producers, never by itself; and at one thread there is no
    /// peer to wake.
    #[rstest]
    #[case::serial(StepKind::Serial, 4)]
    #[case::exclusive(StepKind::Exclusive, 4)]
    #[case::parallel_one_thread(StepKind::Parallel, 1)]
    fn zero_output_cascade_only_for_parallel_sinks(#[case] kind: StepKind, #[case] threads: usize) {
        let (g, mut kinds, pinned, driver_of, queues) = three_step_with_parallel_sink();
        kinds[2] = kind;
        let pool = (threads > 1).then(|| Arc::new(PoolEventCount::new(threads)));
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            pool,
            threads,
        );
        assert_eq!(plan.forward_branch_count_for_test(StepIdx(2)), 0);
        assert_eq!(plan.dag_lines(&g)[2], "      wake: (sink); reverse: Det (holder)");
    }

    #[test]
    fn wake_plan_targets_for_the_sort_topology() {
        use WakeTarget::{Driver, None as NoneT, Pool};
        let (g, kinds, pinned, driver_of, queues) = synthetic_sort_topology();
        let pool = Arc::new(PoolEventCount::new(16));
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            Some(pool),
            16,
        );
        assert_eq!(plan.n_drivers(), 2);
        assert_eq!(plan.n_slots(), 18);
        let coord = Driver(DriverIdx(0));
        let io = Driver(DriverIdx(1));
        // (step, forward, gated, reverse producers)
        let expected: [(usize, WakeTarget, bool, &[usize]); 11] = [
            (0, coord, true, &[]),
            (1, Pool, true, &[0]),
            (2, coord, true, &[1]),
            (3, NoneT, true, &[2]),
            (4, Pool, true, &[]), // from ArenaSort: same thread → omitted
            (5, io, true, &[4]),
            (6, Pool, true, &[5]),
            (7, coord, true, &[6]),
            (8, Pool, true, &[7]),
            (9, io, true, &[8]),
            (10, NoneT, false, &[9]), // sink: no forward branch
        ];
        for (i, fwd, gated, rev) in expected {
            if i != 10 {
                assert_eq!(
                    plan.forward_target(StepIdx(i), BranchIdx(0)),
                    fwd,
                    "forward of step {i}"
                );
            }
            assert_eq!(plan.is_gated(StepIdx(i)), gated, "gate of step {i}");
            let want: Vec<StepIdx> = rev.iter().map(|&p| StepIdx(p)).collect();
            assert_eq!(plan.reverse_producers(StepIdx(i)), want, "reverse of step {i}");
        }
        // Every edge is tracked: every producer with a byte-bounded output is
        // gated, whatever its kind — the same-thread ArenaSort → SpillGather edge
        // has no reverse wake, but ArenaSort is a gated producer, so it is counted.
        for (i, rq) in queues.iter().enumerate() {
            assert!(plan.tracks(rq), "tracking of edge {i}");
        }
    }

    #[rstest]
    #[case::terminal_branch_unwired(false, WakeTarget::None)]
    #[case::fan_out_to_pool(true, WakeTarget::Pool)]
    fn fan_out_branches_resolve_independently(
        #[case] wire_rejects: bool,
        #[case] expected: WakeTarget,
    ) {
        // Src(Exclusive w0) → Det(2 branches: .0 → Par, .1 → Rejects?) ; Par → Sink(Serial, None)
        let mut g = ChainGraph::new();
        let s = g.register_step("Src", 1);
        let d = g.register_step("Det", 2);
        let p = g.register_step("Par", 1);
        let k = g.register_step("Sink", 0);
        let r = g.register_step("Rejects", 0);
        g.wire(s, BranchIdx(0), d);
        g.wire(d, BranchIdx(0), p);
        g.wire(p, BranchIdx(0), k);
        if wire_rejects {
            g.wire(d, BranchIdx(1), r);
        }
        let kinds = [
            StepKind::Exclusive,
            StepKind::Detached,
            StepKind::Parallel,
            StepKind::Serial,
            StepKind::Parallel,
        ];
        let pinned = [Some(0), None, None, None, None];
        let driver_of = [None, Some(DriverIdx(0)), None, None, None];
        let queues = vec![q(0, 0), q(1, 0), q(1, 1), q(2, 0)];
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            Some(Arc::new(PoolEventCount::new(4))),
            4,
        );
        assert_eq!(plan.forward_target(StepIdx(1), BranchIdx(0)), WakeTarget::Pool);
        assert_eq!(plan.forward_target(StepIdx(1), BranchIdx(1)), expected);
        assert_eq!(
            plan.gated_branches(StepIdx(1)),
            2,
            "both byte-bounded branches of a Detached producer are gated"
        );
    }

    #[test]
    fn step_k_consumer_has_one_reverse_edge_per_input() {
        // A(Det, driver 0) ─┐
        //                   ├→ K(Serial, None) → Sink
        // B(Par)          ──┘
        let mut g = ChainGraph::new();
        let det = g.register_step("A", 1);
        let par = g.register_step("B", 1);
        let two_in = g.register_step_with_input_arity("K", 1, 2);
        let sink = g.register_step("Sink", 0);
        g.wire_to_slot(det, BranchIdx(0), two_in, 0);
        g.wire_to_slot(par, BranchIdx(0), two_in, 1);
        g.wire(two_in, BranchIdx(0), sink);
        let kinds = [StepKind::Detached, StepKind::Parallel, StepKind::Serial, StepKind::Exclusive];
        let pinned = [None, None, None, Some(0)];
        let driver_of = [Some(DriverIdx(0)), None, None, None];
        let queues = vec![q(0, 0), q(1, 0), q(2, 0)];
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        assert_eq!(plan.reverse_producers(StepIdx(2)), vec![StepIdx(0), StepIdx(1)]);
        assert_eq!(
            plan.forward_target(StepIdx(2), BranchIdx(0)),
            WakeTarget::Worker(0),
            "Exclusive consumer → its owner"
        );
    }

    #[test]
    fn single_worker_resolves_pool_to_worker_zero() {
        let (g, kinds, pinned, driver_of, queues) = synthetic_sort_topology();
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            None,
            1,
        );
        assert_eq!(plan.forward_target(StepIdx(1), BranchIdx(0)), WakeTarget::Worker(0));
        // ReadBgzfBlocks (w0) → ReadBlocks (driver) still crosses threads at t1.
        assert_eq!(plan.reverse_producers(StepIdx(1)), vec![StepIdx(0)]);
        assert!(plan.pool().is_none());
    }

    /// Every typed step has at most `MAX_ARITY` outputs, so the gate snapshot
    /// never overflows. Pinned at compile time (`const _` in `wake.rs`); this
    /// test documents it at the value level.
    #[test]
    #[allow(clippy::assertions_on_constants)]
    fn gate_snapshot_covers_every_typed_arity() {
        assert!(MAX_GATED_BRANCHES >= crate::outputs::MAX_ARITY);
    }

    /// Unit gate, observed through the event-count: a gated Detached producer's
    /// `Progress` with no push wakes nobody; with a push it wakes the pool once.
    #[test]
    fn gate_wakes_only_when_the_dispatch_pushed() {
        let mut g = ChainGraph::new();
        let d = g.register_step("Det", 1);
        let p = g.register_step("Par", 0);
        g.wire(d, BranchIdx(0), p);
        let kinds = [StepKind::Detached, StepKind::Parallel];
        let typed = Arc::new(ByteBoundedQueue::<u32>::new(1 << 20));
        let handle: Arc<dyn BoundedQueueHandle> = Arc::clone(&typed) as _;
        let rq = RegisteredQueue {
            producer_step_name: "Det",
            producer_step: d,
            branch: BranchIdx(0),
            handle,
            reorder_cap: None,
        };
        let pool = Arc::new(PoolEventCount::new(2));
        let plan = WakePlan::build(
            &g,
            &kinds,
            &[None, None],
            &[Some(DriverIdx(0)), None],
            &[],
            WakeEdges::byte_bounded(std::slice::from_ref(&rq)),
            Some(Arc::clone(&pool)),
            2,
        );
        assert!(plan.tracks(&rq));
        typed.enable_tracking(
            EdgeTracking { clock: None, holders: Some(HolderSet::new(plan.n_slots())) },
            crate::queues::SEALED,
        );

        let mut snap = PushSnapshot::default();
        plan.snapshot_pushes(d, &mut snap);
        let key = pool.prepare_wait();
        let r = plan.on_progress(d, &snap);
        assert_eq!(pool.wait(key, Duration::ZERO), WaitOutcome::TimedOut, "no push → no wake");
        assert_eq!(r.gated_off, 1);

        plan.snapshot_pushes(d, &mut snap);
        typed.try_push(7).unwrap();
        let key = pool.prepare_wait();
        let _ = plan.on_progress(d, &snap);
        assert_eq!(pool.wait(key, Duration::ZERO), WaitOutcome::Woken, "a push → one pool wake");
    }

    /// The reverse wake goes to the recorded holder, by direct unpark, and only
    /// when a rejection was recorded. The consumer is a sink (no forward branch),
    /// so nothing but the reverse edge can wake anyone. The holder is a real
    /// parked thread with a 10 s park; it must return under a 5 s bound.
    #[rstest]
    #[case::byte_bounded(Bounded::Byte)]
    #[case::count_bounded(Bounded::Count)]
    fn reverse_wake_unparks_the_recorded_holder(#[case] kind: Bounded) {
        // P (Parallel, pool) → C (Detached sink on driver 0); t2 → slots 0,1 workers, 2 driver.
        let mut g = ChainGraph::new();
        let p = g.register_step("P", 1);
        let c = g.register_step("C", 0);
        g.wire(p, BranchIdx(0), c);
        let kinds = [StepKind::Parallel, StepKind::Detached];
        let driver_of = [None, Some(DriverIdx(0))];
        let edge = OneItemEdge::tracked(kind, 0, 3);
        let pool = Arc::new(PoolEventCount::new(2));
        let plan = WakePlan::build(
            &g,
            &kinds,
            &[None, None],
            &driver_of,
            &[],
            edge.edges(),
            Some(Arc::clone(&pool)),
            2,
        );
        assert_eq!(plan.reverse_producers(c), vec![p]);

        edge.queue.try_push(Sized8).unwrap(); // full
        let (ready_tx, ready_rx) = std::sync::mpsc::channel::<()>();
        let queue_h = Arc::clone(&edge.queue);
        let plan_h = Arc::clone(&plan);
        let holder = std::thread::spawn(move || {
            let _slot = SlotGuard::enter(1); // worker 1 holds
            assert!(queue_h.try_push(Sized8).is_err(), "rejected → latched as slot 1");
            assert!(plan_h.arm_direct(1), "armed for its timer park");
            ready_tx.send(()).unwrap();
            let t = std::time::Instant::now();
            std::thread::park_timeout(Duration::from_secs(10));
            t.elapsed()
        });
        plan.register_worker(1, holder.thread().clone());
        ready_rx.recv().unwrap();
        // First prove "no rejection → no wake" on a fresh edge.
        let fresh = OneItemEdge::tracked(kind, 0, 3);
        let plan2 = WakePlan::build(
            &g,
            &kinds,
            &[None, None],
            &driver_of,
            &[],
            fresh.edges(),
            Some(Arc::clone(&pool)),
            2,
        );
        assert_eq!(
            plan2.on_progress(c, &PushSnapshot::default()).reverse,
            0,
            "nothing latched → no reverse wake"
        );

        assert_eq!(edge.queue.try_pop(), Some(Sized8)); // the consumer's pop
        let r = plan.on_progress(c, &PushSnapshot::default());
        assert_eq!(r.reverse, 1);
        let waited = holder.join().unwrap();
        assert!(
            waited < Duration::from_secs(5),
            "the holder must be unparked, not timed out: {waited:?}"
        );
    }

    /// Which pool wake `pool_wake_falls_back_to_a_direct_parked_worker` issues.
    #[derive(Clone, Copy, Debug)]
    enum PoolWake {
        Progress,
        RequestWorker,
    }

    /// A pool wake that finds no event-count waiter unparks one direct-parked
    /// worker, whether it is a `Progress`'s forward wake or a step's explicit
    /// `request_worker()`: a timer-parked worker is parked, so "every worker
    /// is awake" would be false. The worker is a real thread parked for 10 s.
    #[rstest]
    #[case::progress(PoolWake::Progress)]
    #[case::request_worker(PoolWake::RequestWorker)]
    fn pool_wake_falls_back_to_a_direct_parked_worker(#[case] via: PoolWake) {
        let mut g = ChainGraph::new();
        let d = g.register_step("Det", 1);
        let p = g.register_step("Par", 0);
        g.wire(d, BranchIdx(0), p);
        let pool = Arc::new(PoolEventCount::new(2));
        let plan = WakePlan::build(
            &g,
            &[StepKind::Detached, StepKind::Parallel],
            &[None, None],
            &[Some(DriverIdx(0)), None],
            &[],
            WakeEdges::NONE,
            Some(Arc::clone(&pool)),
            2,
        );
        let (armed_tx, armed_rx) = std::sync::mpsc::channel::<()>();
        let plan_w = Arc::clone(&plan);
        let worker = std::thread::spawn(move || {
            assert!(plan_w.arm_direct(0));
            armed_tx.send(()).unwrap();
            let t = std::time::Instant::now();
            std::thread::park_timeout(Duration::from_secs(10));
            plan_w.disarm_direct(0);
            t.elapsed()
        });
        plan.register_worker(0, worker.thread().clone());
        armed_rx.recv().unwrap();
        let handle = PoolHandle::new(&plan, None, d);
        match via {
            PoolWake::Progress => {
                // Ungated (no queue): wakes Pool.
                let r = plan.on_progress(d, &PushSnapshot::default());
                assert_eq!((r.no_waiters, r.fallback), (1, 1));
            }
            PoolWake::RequestWorker => {
                assert_eq!(handle.request_worker(), PoolRequest::Woken, "the armed worker");
            }
        }
        assert!(worker.join().unwrap() < Duration::from_secs(5));
        // With nobody armed, the fallback finds nobody.
        match via {
            PoolWake::Progress => {
                assert_eq!(plan.on_progress(d, &PushSnapshot::default()).fallback, 0);
            }
            PoolWake::RequestWorker => {
                assert_eq!(handle.request_worker(), PoolRequest::AllAwake);
            }
        }
    }

    /// `host_mask` restricts only a consumer whose clones live on some of the
    /// workers, and its bits span words past 64 workers.
    #[test]
    fn host_mask_sets_the_host_workers_bits_across_words() {
        use crate::runtime::placement::ParallelHosts;
        let all_but = |skip: &[usize]| ParallelHosts {
            workers: (0..130).filter(|w| !skip.contains(w)).collect(),
            driver: None,
        };
        let hosts = [
            all_but(&[0, 65]),
            all_but(&[]),
            ParallelHosts { workers: vec![], driver: Some(DriverIdx(0)) },
        ];
        let mask = host_mask(&hosts, StepIdx(0), 130).expect("two workers left out");
        assert_eq!(&*mask, &[!1u64, !(1u64 << 1), 0b11]);
        assert_eq!(host_mask(&hosts, StepIdx(1), 130), None, "every worker: unrestricted");
        assert_eq!(host_mask(&hosts, StepIdx(2), 130), None, "hosted on a driver: no pool wake");
        assert_eq!(host_mask(&hosts, StepIdx(3), 130), None, "not covered");
    }

    /// The direct-park fallback for a `Pool` wake claims only a worker that
    /// holds a clone of the consumer: with `Par` placed on worker 1 alone, an
    /// armed worker 0 is passed over (and stays parked), and an armed worker 1
    /// is unparked — a real thread parked for 10 s, returning in < 5 s.
    #[test]
    fn pool_fallback_skips_workers_without_a_clone() {
        use crate::runtime::placement::ParallelHosts;
        let mut g = ChainGraph::new();
        let d = g.register_step("Det", 1);
        let p = g.register_step("Par", 0);
        g.wire(d, BranchIdx(0), p);
        let pool = Arc::new(PoolEventCount::new(2));
        let plan = WakePlan::build_with_hosts(
            &g,
            &[StepKind::Detached, StepKind::Parallel],
            &[None, None],
            &[Some(DriverIdx(0)), None],
            &[],
            WakeEdges::NONE,
            Some(Arc::clone(&pool)),
            2,
            &[ParallelHosts::default(), ParallelHosts { workers: vec![1], driver: None }],
        );
        assert_eq!(plan.forward_target(d, BranchIdx(0)), WakeTarget::Pool);
        let (w0, w0_parked) = parked_thread();
        plan.register_worker(0, w0.thread().clone());
        assert!(plan.arm_direct(0));
        w0_parked.recv().unwrap();
        let r = plan.on_progress(d, &PushSnapshot::default());
        assert_eq!((r.no_waiters, r.fallback), (1, 0), "worker 0 holds no clone of Par");

        let (w1, w1_parked) = parked_thread();
        plan.register_worker(1, w1.thread().clone());
        assert!(plan.arm_direct(1));
        w1_parked.recv().unwrap();
        let r = plan.on_progress(d, &PushSnapshot::default());
        assert_eq!(r.fallback, 1, "worker 1 hosts Par");
        assert!(w1.join().unwrap() < Duration::from_secs(5));
        // Worker 0 was passed over by both wakes: still parked.
        std::thread::sleep(Duration::from_millis(100));
        assert!(!w0.is_finished(), "a wake reached worker 0, which holds no clone of Par");
        // Release worker 0, which no wake above reached.
        plan.disarm_direct(0);
        plan.unpark_all();
        let _ = w0.join().unwrap();
    }

    /// The cap-parked fallback honours the same host mask as the direct-park
    /// one: a cap-parked worker holding no clone of the consumer is never
    /// claimed for it (the wake would be spent on a thread that cannot run it).
    #[test]
    fn cap_parked_fallback_skips_workers_without_a_clone() {
        use crate::runtime::placement::ParallelHosts;
        let mut g = ChainGraph::new();
        let d = g.register_step("Det", 1);
        let p = g.register_step("Par", 0);
        g.wire(d, BranchIdx(0), p);
        let pool = Arc::new(PoolEventCount::new(2));
        let plan = WakePlan::build_with_hosts(
            &g,
            &[StepKind::Detached, StepKind::Parallel],
            &[None, None],
            &[Some(DriverIdx(0)), None],
            &[],
            WakeEdges::NONE,
            Some(Arc::clone(&pool)),
            2,
            &[ParallelHosts::default(), ParallelHosts { workers: vec![1], driver: None }],
        );
        let (w0, w0_parked) = parked_thread();
        plan.register_worker(0, w0.thread().clone());
        assert!(plan.arm_cap_parked(0));
        w0_parked.recv().unwrap();
        let r = plan.on_progress(d, &PushSnapshot::default());
        assert_eq!((r.no_waiters, r.fallback), (1, 0), "worker 0 holds no clone of Par");
        assert!(plan.is_cap_parked(0), "worker 0's cap-parked bit was not claimed");
        let (w1, w1_parked) = parked_thread();
        plan.register_worker(1, w1.thread().clone());
        assert!(plan.arm_cap_parked(1));
        w1_parked.recv().unwrap();
        let r = plan.on_progress(d, &PushSnapshot::default());
        assert_eq!(r.fallback, 1, "worker 1 hosts Par");
        assert!(w1.join().unwrap() < Duration::from_secs(5));
        std::thread::sleep(Duration::from_millis(100));
        assert!(!w0.is_finished(), "a wake reached worker 0, which holds no clone of Par");
        plan.disarm_cap_parked(0);
        plan.unpark_all();
        let _ = w0.join().unwrap();
    }

    #[test]
    fn registration_is_once_and_mode_gated() {
        let (g, kinds, pinned, driver_of, queues) = synthetic_sort_topology();
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            Some(Arc::new(PoolEventCount::new(16))),
            16,
        );
        assert!(!plan.is_worker_registered(0));
        plan.register_worker(0, std::thread::current());
        assert!(plan.is_worker_registered(0));
        plan.register_worker(0, std::thread::current()); // second set is ignored, no panic
        plan.register_driver(DriverIdx(1), std::thread::current());
        // Out-of-range registrations are ignored rather than panicking: the plan
        // sizes `workers` to `n_threads`, so a caller bug here must not take down a run.
        plan.register_worker(99, std::thread::current());
        assert!(!plan.is_worker_registered(99));
    }

    /// `Par(Parallel) → Ser(Serial, pinned to worker 1)` at two workers, plus
    /// `Feed(Parallel) → Det(Detached, driver 0) → Par2(Parallel)` so the plan
    /// has a driver: `Par`'s branch targets `Worker(1)` and `Feed`'s targets
    /// `Driver(0)`.
    fn worker_and_driver_plan() -> Arc<WakePlan> {
        let mut g = ChainGraph::new();
        let par = g.register_step("Par", 1);
        let ser = g.register_step("Ser", 0);
        let feed = g.register_step("Feed", 1);
        let det = g.register_step("Det", 1);
        let par2 = g.register_step("Par2", 0);
        g.wire(par, BranchIdx(0), ser);
        g.wire(feed, BranchIdx(0), det);
        g.wire(det, BranchIdx(0), par2);
        WakePlan::build(
            &g,
            &[
                StepKind::Parallel,
                StepKind::Serial,
                StepKind::Parallel,
                StepKind::Detached,
                StepKind::Parallel,
            ],
            &[None, Some(1), None, None, None],
            &[None, None, None, Some(DriverIdx(0)), None],
            &[],
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        )
    }

    /// Drain this thread's park token, so the next check starts clean.
    fn drain_token() {
        std::thread::park_timeout(Duration::ZERO);
    }

    /// A park of this thread ended before its 10 s timeout: an unpark reached it.
    fn was_unparked() -> bool {
        let t = std::time::Instant::now();
        std::thread::park_timeout(Duration::from_secs(10));
        t.elapsed() < Duration::from_secs(5)
    }

    /// A directed wake unparks only a thread armed for its timer park. While
    /// worker 1 runs (its parking bit clear) a `Worker(1)` forward wake and a
    /// holder wake through `deliver_slot` issue no unpark; once it arms they
    /// do; disarming clears the bit again. The same for driver 0, whose bit is
    /// its wake slot `n_workers + 0`, through a `Driver(0)` forward wake and
    /// `deliver_slot(2)`. The skips are asserted through the returned counts,
    /// not through a park's timing.
    #[test]
    fn a_directed_wake_unparks_only_an_armed_thread() {
        let plan = worker_and_driver_plan();
        assert_eq!(plan.forward_target(StepIdx(0), BranchIdx(0)), WakeTarget::Worker(1));
        assert_eq!(plan.forward_target(StepIdx(2), BranchIdx(0)), WakeTarget::Driver(DriverIdx(0)));
        plan.register_worker(1, std::thread::current());
        plan.register_driver(DriverIdx(0), std::thread::current());
        drain_token();

        // Worker 1.
        assert!(!plan.is_parking(1));
        assert_eq!(plan.on_progress(StepIdx(0), &PushSnapshot::default()).unparked, 0);
        assert!(!plan.deliver_slot(1), "a running worker was unparked");
        assert!(plan.arm_direct(1));
        assert!(plan.is_parking(1));
        assert_eq!(plan.on_progress(StepIdx(0), &PushSnapshot::default()).unparked, 1);
        assert!(was_unparked(), "the armed worker was not unparked");
        assert!(plan.deliver_slot(1), "deliver_slot skipped an armed worker");
        drain_token();
        plan.disarm_direct(1);
        assert!(!plan.is_parking(1));
        assert!(!plan.deliver_slot(1), "a disarmed worker was unparked");

        assert!(plan.arm_cap_parked(1));
        assert!(plan.is_parking(1), "a cap-parked worker is armed for directed wakes too");
        assert!(plan.deliver_slot(1), "a cap-release wake skipped a cap-parked worker");
        drain_token();
        plan.disarm_cap_parked(1);
        assert!(!plan.is_parking(1));

        // Driver 0, wake slot 2.
        assert!(!plan.is_parking(2));
        assert_eq!(plan.on_progress(StepIdx(2), &PushSnapshot::default()).unparked, 0);
        assert!(!plan.deliver_slot(2), "a running driver was unparked");
        assert!(plan.arm_driver(DriverIdx(0)));
        assert!(plan.is_parking(2), "driver 0 is wake slot n_workers + 0");
        assert!(!plan.is_parking(0) && !plan.is_parking(1));
        assert_eq!(plan.on_progress(StepIdx(2), &PushSnapshot::default()).unparked, 1);
        assert!(was_unparked(), "the armed driver was not unparked");
        assert!(plan.deliver_slot(2), "deliver_slot skipped an armed driver");
        drain_token();
        plan.disarm_driver(DriverIdx(0));
        assert!(!plan.is_parking(2));
        assert!(!plan.deliver_slot(2), "a disarmed driver was unparked");
    }

    /// A driver index the plan does not have is loud: a driver parked with no
    /// bit of its own would be skipped by every directed wake, silently.
    #[test]
    #[should_panic(expected = "arm_driver: driver 1 is not in this plan (1 drivers)")]
    fn arming_a_driver_the_plan_lacks_panics() {
        let _ = worker_and_driver_plan().arm_driver(DriverIdx(1));
    }

    /// A plan with one thread arms nothing (the arm and its re-poll would be
    /// pure cost) and gates nothing: a directed wake to its lone worker is an
    /// unconditional unpark, since a thread outside the plan (the queue-memory
    /// rebalancer) still delivers there. One driver beside the lone worker
    /// makes the arm, and the gate, real. The standalone plan arms nothing.
    #[test]
    fn a_one_thread_plan_arms_nothing_and_gates_nothing() {
        let one_step = |driver: bool| {
            let mut g = ChainGraph::new();
            let a = g.register_step("A", 1);
            let b = g.register_step("B", 0);
            g.wire(a, BranchIdx(0), b);
            let kinds =
                [if driver { StepKind::Detached } else { StepKind::Serial }, StepKind::Serial];
            let driver_of = [driver.then_some(DriverIdx(0)), None];
            WakePlan::build(&g, &kinds, &[None, None], &driver_of, &[], WakeEdges::NONE, None, 1)
        };
        let lone = one_step(false);
        assert!(!lone.arm_direct(0) && !lone.arm_cap_parked(0));
        assert!(!lone.is_parking(0));
        lone.register_worker(0, std::thread::current());
        drain_token();
        assert!(lone.deliver_slot(0), "the unarmed lone worker is still unparked");
        assert!(was_unparked());

        let with_driver = one_step(true);
        with_driver.register_worker(0, std::thread::current());
        assert!(!with_driver.deliver_slot(0), "a plan with two threads gates the lone worker");
        assert!(with_driver.arm_direct(0) && with_driver.is_parking(0));
        assert!(with_driver.deliver_slot(0));
        drain_token();
        assert!(with_driver.arm_driver(DriverIdx(0)) && with_driver.is_parking(1));
        assert!(!WakePlan::standalone().arm_direct(0));
    }

    #[test]
    fn dag_lines_render_one_line_per_step() {
        let (g, kinds, pinned, driver_of, queues) = synthetic_sort_topology();
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            Some(Arc::new(PoolEventCount::new(16))),
            16,
        );
        let lines = plan.dag_lines(&g);
        assert_eq!(lines.len(), 11);
        assert_eq!(lines[0], "      wake: .0 → ReadBlocks Driver(0) [gated]");
        assert_eq!(
            lines[1],
            "      wake: .0 → InflateToArena Pool [gated]; reverse: ReadBgzfBlocks (holder)"
        );
        assert_eq!(
            lines[2],
            "      wake: .0 → ArenaSort Driver(0) [gated]; reverse: ReadBlocks (holder)"
        );
        assert_eq!(lines[10], "      wake: (sink); reverse: BgzfCompress (holder)");
    }

    /// Arm a pool waiter on its own thread; returns once `prepare_wait` ran.
    fn armed_pool_waiter(pool: &Arc<PoolEventCount>) -> std::thread::JoinHandle<WaitOutcome> {
        let (tx, rx) = std::sync::mpsc::channel::<()>();
        let p = Arc::clone(pool);
        let h = std::thread::spawn(move || {
            let key = p.prepare_wait();
            tx.send(()).unwrap();
            p.wait(key, Duration::from_secs(30))
        });
        rx.recv().unwrap();
        h
    }

    /// The four outcomes and the one-pending-request rule.
    #[test]
    fn request_worker_semantics() {
        use crate::runtime::event_count::NotifyOutcome;
        // No pool, no worker (the standalone plan): Unavailable.
        let standalone = WakePlan::standalone();
        assert_eq!(
            PoolHandle::new(&standalone, None, StepIdx(0)).request_worker(),
            PoolRequest::Unavailable
        );

        // Pool with nobody armed: AllAwake, and the generation does not move.
        let pool = Arc::new(PoolEventCount::new(2));
        let plan = super::tests_support::plan_with_pool(Arc::clone(&pool), 2);
        let h = PoolHandle::new(&plan, None, StepIdx(0));
        assert_eq!(h.parked_workers(), 0);
        assert_eq!(h.request_worker(), PoolRequest::AllAwake);
        let key = pool.prepare_wait();
        assert_eq!(
            pool.wait(key, Duration::ZERO),
            WaitOutcome::TimedOut,
            "AllAwake must not bump the generation"
        );

        // One armed worker: Woken, and the waiter is released (not timed out).
        let waiter = armed_pool_waiter(&pool);
        assert_eq!(h.parked_workers(), 1);
        assert_eq!(h.request_worker(), PoolRequest::Woken);
        assert!(matches!(waiter.join().unwrap(), WaitOutcome::Woken | WaitOutcome::Notified));

        // A second request before any acknowledgement is Pending; after
        // acknowledge → Woken again.
        let waiter = armed_pool_waiter(&pool);
        assert_eq!(h.request_worker(), PoolRequest::Pending, "first request still unacknowledged");
        plan.acknowledge_request();
        assert_eq!(h.request_worker(), PoolRequest::Woken);
        assert!(matches!(waiter.join().unwrap(), WaitOutcome::Woken | WaitOutcome::Notified));

        // Ceiling 1 with one worker awake: `notify_one` is Suppressed, the demand
        // path is ungated → Woken.
        pool.set_ceiling(1);
        plan.acknowledge_request();
        let waiter = armed_pool_waiter(&pool);
        assert_eq!(pool.notify_one(), NotifyOutcome::Suppressed);
        assert_eq!(h.request_worker(), PoolRequest::Woken);
        assert!(matches!(waiter.join().unwrap(), WaitOutcome::Woken | WaitOutcome::Notified));
    }

    /// One worker: no pool, worker 0 registered → the request
    /// unparks it and reports Woken; unregistered → Unavailable.
    #[test]
    fn request_worker_at_one_thread_unparks_worker_zero() {
        let (g, kinds, pinned, driver_of, queues) = synthetic_sort_topology();
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            None,
            1,
        );
        let h = PoolHandle::new(&plan, None, StepIdx(8));
        assert_eq!(h.request_worker(), PoolRequest::Unavailable);
        plan.register_worker(0, std::thread::current());
        assert_eq!(h.request_worker(), PoolRequest::Woken);
        // The token is pending on this thread: a park returns at once.
        let t = std::time::Instant::now();
        std::thread::park_timeout(Duration::from_secs(5));
        assert!(t.elapsed() < Duration::from_millis(500));
        // From worker 0 itself there is nobody else to wake: no self-unpark.
        let _slot = SlotGuard::enter(0);
        assert_eq!(
            h.request_worker(),
            PoolRequest::AllAwake,
            "the lone worker does not wake itself"
        );
        let t = std::time::Instant::now();
        std::thread::park_timeout(Duration::from_millis(200));
        assert!(t.elapsed() >= Duration::from_millis(150), "no token was left on this thread");
    }

    /// The request `SortMerge` makes when it stalls (`PoolHandle::request_worker`)
    /// reaches a worker parked on its timer (a pinned worker or a holder, armed
    /// direct-parked) when no worker waits on the event-count: the timer-parked
    /// worker is unparked at once rather than left to its timer.
    #[test]
    fn a_merge_request_reaches_a_timer_parked_worker() {
        let (plan, _g) = tests_support::directed_plan_with_workers(2);
        plan.register_worker(1, std::thread::current());
        assert!(plan.arm_direct(1));
        let h = PoolHandle::new(&plan, None, StepIdx(0));
        assert_eq!(h.request_worker(), PoolRequest::Woken);
        let t = std::time::Instant::now();
        std::thread::park_timeout(Duration::from_secs(10));
        assert!(t.elapsed() < Duration::from_secs(5), "worker 1 (this thread) was unparked");
        assert_eq!(h.request_worker(), PoolRequest::AllAwake, "its bit was claimed");
    }

    /// With no event-count waiter and no direct-parked worker, a request wakes a
    /// worker parked only because a phase cap refused it.
    #[test]
    fn a_request_reaches_a_cap_parked_worker() {
        let (plan, _g) = tests_support::directed_plan_with_workers(2);
        plan.register_worker(1, std::thread::current());
        assert!(plan.arm_cap_parked(1));
        assert_eq!(plan.request_pool_worker(), PoolRequest::Woken);
        let t = std::time::Instant::now();
        std::thread::park_timeout(Duration::from_secs(10));
        assert!(t.elapsed() < Duration::from_secs(5), "worker 1 (this thread) was unparked");
        assert_eq!(plan.request_pool_worker(), PoolRequest::AllAwake, "its bit was claimed");
    }

    /// A worker's request never wakes the worker itself: armed on the
    /// event-count (a step in its re-poll) it is not counted as a waiter, and
    /// its own armed timer-park bit is skipped. With nobody else parked the
    /// request is `AllAwake`, and the requester's own wait is not cut short.
    #[test]
    fn a_request_never_wakes_the_requester() {
        let (plan, _g) = tests_support::directed_plan_with_workers(2);
        let pool = plan.pool().expect("two workers: an event-count");
        let _slot = SlotGuard::enter(1);
        plan.register_worker(1, std::thread::current());
        // Armed on the event-count, as in the worker loop's re-poll.
        let key = pool.prepare_wait();
        crate::runtime::wake_slot::set_ec_armed(true);
        assert_eq!(
            plan.request_pool_worker(),
            PoolRequest::AllAwake,
            "the only waiter is the caller"
        );
        crate::runtime::wake_slot::set_ec_armed(false);
        assert_eq!(pool.wait(key, Duration::from_millis(50)), WaitOutcome::TimedOut, "not woken");
        // Armed for a timer park, as in the worker loop's direct-park re-poll.
        assert!(plan.arm_direct(1));
        assert_eq!(plan.request_pool_worker(), PoolRequest::AllAwake, "its own bit is skipped");
        let t = std::time::Instant::now();
        std::thread::park_timeout(Duration::from_millis(200));
        assert!(t.elapsed() >= Duration::from_millis(150), "no self-unpark token");
        plan.disarm_direct(1);
    }

    /// Stats-on: each outcome lands in its per-step counter.
    #[test]
    fn request_worker_records_per_step_outcomes() {
        let stats = crate::runtime::stats::PipelineStats::new(vec!["A", "B"]);
        let standalone = WakePlan::standalone();
        let h = PoolHandle::new(&standalone, Some(&stats), StepIdx(1));
        assert_eq!(h.request_worker(), PoolRequest::Unavailable);
        let pool = Arc::new(PoolEventCount::new(2));
        let plan = super::tests_support::plan_with_pool(Arc::clone(&pool), 2);
        let h = PoolHandle::new(&plan, Some(&stats), StepIdx(1));
        assert_eq!(h.request_worker(), PoolRequest::AllAwake);
        // An armed waiter: the first request wakes it, a second one before any
        // acknowledgement is pending.
        let waiter = armed_pool_waiter(&pool);
        assert_eq!(h.request_worker(), PoolRequest::Woken);
        assert_eq!(h.request_worker(), PoolRequest::Pending);
        assert_eq!(h.request_worker(), PoolRequest::Pending);
        assert!(matches!(waiter.join().unwrap(), WaitOutcome::Woken | WaitOutcome::Notified));
        let s = stats.snapshot().steps[1].1;
        assert_eq!(
            (s.pool_requests_woken, s.pool_requests_all_awake, s.pool_requests_pending),
            (1, 1, 2),
            "each outcome in its own counter; Unavailable is not counted"
        );
        assert_eq!(stats.snapshot().steps[0].1.pool_requests_woken, 0, "step A untouched");
    }

    /// A thread that parks for up to 10 s after telling the test it is about to,
    /// and returns how long it was parked.
    fn parked_thread() -> (std::thread::JoinHandle<Duration>, std::sync::mpsc::Receiver<()>) {
        let (tx, rx) = std::sync::mpsc::channel::<()>();
        let h = std::thread::spawn(move || {
            tx.send(()).unwrap();
            let t = std::time::Instant::now();
            std::thread::park_timeout(Duration::from_secs(10));
            t.elapsed()
        });
        (h, rx)
    }

    /// `on_finished` from a step whose forward target is `None` (same-driver
    /// consumer) still unparks every registered thread — two drivers and a
    /// worker, each parked for 10 s.
    #[test]
    fn on_finished_unparks_every_registered_thread() {
        // Src(Det d0) → A(Det d0) → Sink(Det d1) at t2: Src's forward is None.
        let mut g = ChainGraph::new();
        let src = g.register_step("Src", 1);
        let a = g.register_step("A", 1);
        let k = g.register_step("Sink", 0);
        g.wire(src, BranchIdx(0), a);
        g.wire(a, BranchIdx(0), k);
        let kinds = [StepKind::Detached; 3];
        let driver_of = [Some(DriverIdx(0)), Some(DriverIdx(0)), Some(DriverIdx(1))];
        let plan = WakePlan::build(
            &g,
            &kinds,
            &[None; 3],
            &driver_of,
            &[],
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        assert_eq!(plan.forward_target(src, BranchIdx(0)), WakeTarget::None);
        let threads: Vec<_> = (0..3).map(|_| parked_thread()).collect();
        plan.register_driver(DriverIdx(0), threads[0].0.thread().clone());
        plan.register_driver(DriverIdx(1), threads[1].0.thread().clone());
        plan.register_worker(0, threads[2].0.thread().clone());
        for (_, ready) in &threads {
            ready.recv().unwrap();
        }
        plan.on_finished(src);
        for (h, _) in threads {
            assert!(
                h.join().unwrap() < Duration::from_secs(5),
                "every registered thread must be unparked"
            );
        }
    }

    /// A smoke test of the registration race on real threads: a third thread
    /// registers the target while the producer runs `on_progress`; the
    /// target's first pass is sequenced after its registration (as in
    /// `Pipeline::run`, where a thread is registered before it can run a pass on
    /// pushed work). A lost unpark leaves the target parked for 10 s and the
    /// watchdog fires. Not fence coverage: the channel that sequences the pass
    /// after the registration hides a missing fence; the loom registration
    /// models in `tests/loom_wake.rs` are what check the fences.
    #[test]
    fn stress_no_lost_unpark() {
        super::e2e::run_under_watchdog("stress unpark", 30, || {
            let (g, kinds, pinned, driver_of, queues) = synthetic_sort_topology();
            for _ in 0..2_000 {
                // Fresh plan per iteration: `OnceLock` registers once.
                let plan = WakePlan::build(
                    &g,
                    &kinds,
                    &pinned,
                    &driver_of,
                    &[],
                    WakeEdges::byte_bounded(&queues),
                    None,
                    1,
                );
                let work = Arc::new(std::sync::atomic::AtomicBool::new(false));
                let (registered_tx, registered_rx) = std::sync::mpsc::channel::<()>();
                let work_t = Arc::clone(&work);
                let plan_t = Arc::clone(&plan);
                let target = std::thread::spawn(move || {
                    registered_rx.recv().unwrap(); // first pass after registration
                    // The driver arms for its timer park before that pass.
                    assert!(plan_t.arm_driver(DriverIdx(0)));
                    // `Relaxed`, as a queue publish is: only the registration
                    // and delivery fences may order this load against the store.
                    if !work_t.load(std::sync::atomic::Ordering::Relaxed) {
                        std::thread::park_timeout(Duration::from_secs(10));
                    }
                });
                let plan_r = Arc::clone(&plan);
                let handle = target.thread().clone();
                // Release the registrar and the producer together, so the
                // registration races the publish + `on_progress` below.
                let start = Arc::new(std::sync::Barrier::new(2));
                let start_r = Arc::clone(&start);
                let registrar = std::thread::spawn(move || {
                    start_r.wait();
                    plan_r.register_driver(DriverIdx(0), handle);
                    registered_tx.send(()).unwrap();
                });
                start.wait();
                // `Relaxed` publish: a `SeqCst` store would order itself against
                // the registration and let the test pass without `on_progress`'s
                // fence.
                work.store(true, std::sync::atomic::Ordering::Relaxed);
                // InflateToArena (2) → ArenaSort on Driver(0). The
                // branch is push-gated; a snapshot no counter can equal stands
                // for a dispatch that pushed.
                let pushed = PushSnapshot([u64::MAX; MAX_GATED_BRANCHES]);
                let _ = plan.on_progress(StepIdx(2), &pushed);
                registrar.join().unwrap();
                target.join().unwrap();
            }
        });
    }

    /// A reorder-cap fixture: the transport, its registry entry (in one of the
    /// two lists), how many ordinals fill it, and whether a drain moves them
    /// into the stash.
    type StashFixture = (
        Arc<dyn ItemQueue<crate::reorder::Sequenced<Sized8>>>,
        Vec<RegisteredQueue>,
        Vec<RegisteredHolderOnlyQueue>,
        u64,
        bool,
    );

    /// The transports an ordered branch can sit on.
    #[derive(Clone, Copy, Debug)]
    enum Stash {
        Byte,
        Count,
        Unbounded,
    }

    /// The reorder-stash cap is a rejection too, on every transport: the holder
    /// is recorded under the stage's lock and woken by the consumer's in-order
    /// pop. The stash cap is 8 bytes. Byte (limit 9) and count (capacity 2):
    /// ordinals 1 and 2 fill the transport and 3 goes to the stash. Unbounded:
    /// ordinals 1–3 go to the transport, and an in-order pop with 0 missing
    /// moves them into the stash. Ordinal 4 then hits the cap and is held. The
    /// stash inserts record nothing, so the cap refusal is the only holder. The
    /// test thread pushes ordinal 0 (exempt from the cap) as slot 2, an
    /// unregistered driver slot.
    #[rstest]
    #[case::byte_bounded(Stash::Byte)]
    #[case::count_bounded(Stash::Count)]
    #[case::unbounded(Stash::Unbounded)]
    fn reorder_cap_rejection_wakes_the_holder_on_the_in_order_pop(#[case] kind: Stash) {
        use std::sync::mpsc::channel;

        use crate::queues::{CountBoundedQueue, UnboundedQueue};
        use crate::reorder::{ReorderStage, Sequenced};
        use crate::runtime::wake_slot::rejected_pending;
        let mut g = ChainGraph::new();
        let p = g.register_step("P", 1);
        let c = g.register_step("C", 0);
        g.wire(p, BranchIdx(0), c);
        let (queue, byte, holder_only, fill, drain): StashFixture = match kind {
            Stash::Byte => {
                let t = Arc::new(ByteBoundedQueue::<Sequenced<Sized8>>::new(9));
                let rq = RegisteredQueue {
                    producer_step_name: "P",
                    producer_step: p,
                    branch: BranchIdx(0),
                    handle: Arc::clone(&t) as _,
                    reorder_cap: None,
                };
                (t, vec![rq], vec![], 3, false)
            }
            Stash::Count => {
                let t = Arc::new(CountBoundedQueue::<Sequenced<Sized8>>::new(2));
                let rq = RegisteredHolderOnlyQueue {
                    producer_step: p,
                    branch: BranchIdx(0),
                    handle: Arc::clone(&t) as _,
                };
                (t, vec![], vec![rq], 3, false)
            }
            Stash::Unbounded => {
                let t = Arc::new(UnboundedQueue::<Sequenced<Sized8>>::new());
                let rq = RegisteredHolderOnlyQueue {
                    producer_step: p,
                    branch: BranchIdx(0),
                    handle: Arc::clone(&t) as _,
                };
                (t, vec![], vec![rq], 3, true)
            }
        };
        // t2 + one driver: slots 0, 1 are workers, slot 2 is driver 0 (never registered here).
        let plan = WakePlan::build(
            &g,
            &[StepKind::Parallel, StepKind::Detached],
            &[None, None],
            &[None, Some(DriverIdx(0))],
            &[],
            WakeEdges { byte_bounded: &byte, holder_only: &holder_only },
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        // Install tracking where the plan says, as `Pipeline::run` does.
        for rq in &byte {
            assert!(plan.tracks(rq));
            rq.handle.enable_tracking(
                EdgeTracking { clock: None, holders: Some(HolderSet::new(plan.n_slots())) },
                crate::queues::SEALED,
            );
        }
        for rq in &holder_only {
            assert!(plan.tracks_holder_only(rq));
            rq.handle.enable_holder_tracking(plan.n_slots());
        }
        assert!(plan.reverse_edges_are_tracked());
        let stage = Arc::new(ReorderStage::with_max_overflow_bytes(queue, 8));

        let (ready_tx, ready_rx) = channel::<()>();
        let stage_h = Arc::clone(&stage);
        let plan_h = Arc::clone(&plan);
        let holder = std::thread::spawn(move || {
            let _slot = SlotGuard::enter(1);
            let _ = crate::runtime::wake_slot::take_rejected();
            for ordinal in 1..=fill {
                stage_h.try_push(ordinal, Sized8).unwrap();
            }
            if drain {
                assert_eq!(
                    stage_h.try_pop_in_order(),
                    None,
                    "0 is missing: the drain only stashes"
                );
            }
            assert!(!rejected_pending(), "the stash inserts mark nothing");
            assert!(stage_h.try_push(4, Sized8).is_err(), "stash at its cap: ordinal 4 is held");
            assert!(rejected_pending(), "a stash-cap hold marks the thread");
            assert!(plan_h.arm_direct(1), "armed for its timer park");
            ready_tx.send(()).unwrap();
            let t = std::time::Instant::now();
            std::thread::park_timeout(Duration::from_secs(10));
            t.elapsed()
        });
        plan.register_worker(1, holder.thread().clone());
        ready_rx.recv().unwrap();
        {
            let _slot = SlotGuard::enter(2);
            stage.try_push(0, Sized8).unwrap(); // next_serial: exempt from the cap
        }
        assert_eq!(stage.try_pop_in_order(), Some(Sized8));
        let r = plan.on_progress(c, &PushSnapshot::default());
        assert_eq!(r.reverse, 1, "slot 1, from the cap refusal only");
        assert!(
            holder.join().unwrap() < Duration::from_secs(5),
            "the stash-cap holder must be unparked"
        );
    }

    /// Holder-only edges get a reverse wake for every cross-thread producer,
    /// whatever its role, and are never gated. Every edge is count-bounded:
    /// `Src(Exclusive w0) → Par(Parallel) → Det(Detached d0) → Det2(Detached d0)
    /// → Sink(Serial)` at t2. `Det2 ← Det` is same-driver, so it has none.
    #[test]
    fn holder_only_edges_get_reverse_wakes_for_every_cross_thread_producer() {
        let mut g = ChainGraph::new();
        let src = g.register_step("Src", 1);
        let par = g.register_step("Par", 1);
        let det = g.register_step("Det", 1);
        let det2 = g.register_step("Det2", 1);
        let sink = g.register_step("Sink", 0);
        for (a, b) in [(src, par), (par, det), (det, det2), (det2, sink)] {
            g.wire(a, BranchIdx(0), b);
        }
        let kinds = [
            StepKind::Exclusive,
            StepKind::Parallel,
            StepKind::Detached,
            StepKind::Detached,
            StepKind::Serial,
        ];
        let d0 = Some(DriverIdx(0));
        let holder_only: Vec<RegisteredHolderOnlyQueue> = (0..4).map(count_edge).collect();
        let plan = WakePlan::build(
            &g,
            &kinds,
            &[Some(0), None, None, None, None],
            &[None, None, d0, d0, None],
            &[],
            WakeEdges { byte_bounded: &[], holder_only: &holder_only },
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        assert_eq!(plan.reverse_producers(par), vec![src]);
        assert_eq!(plan.reverse_producers(det), vec![par]);
        assert!(plan.reverse_producers(det2).is_empty(), "same driver: no reverse edge");
        assert_eq!(plan.reverse_producers(sink), vec![det2]);
        let tracked: Vec<bool> = holder_only.iter().map(|q| plan.tracks_holder_only(q)).collect();
        assert_eq!(tracked, vec![true, true, false, true]);
        for step in [src, par, det, det2, sink] {
            assert!(!plan.is_gated(step), "a holder-only edge is never gated ({step:?})");
        }
        assert!(plan.dag_lines(&g)[par.0].ends_with("; reverse: Src (holder)"));
        assert!(plan.steps[sink.0].needs_fence, "Sink's only reason is its holder-only input");
        assert!(plan.steps[det.0].needs_fence, "Det's forward target is None (same driver)");
    }

    /// An ordered byte-bounded branch of a Detached producer is gated like a
    /// direct one: its must-accept stash insert counts as a push
    /// (`ItemQueue::note_stash_accepted`), so the gate cannot hide it. The reverse edge
    /// and its tracking are kept either way.
    /// `A(Exclusive w0)` and `B(Exclusive w0)` feed the two input slots of
    /// `M(Detached d0)`, then `M → C(Parallel)`, at t2.
    #[rstest]
    #[case::direct(false)]
    #[case::ordered(true)]
    fn ordered_branch_of_a_detached_producer_is_gated(#[case] ordered: bool) {
        use crate::reorder::{ReorderCapHandle, ReorderStage, Sequenced};
        let mut g = ChainGraph::new();
        let a = g.register_step("A", 1);
        let b = g.register_step("B", 1);
        let m = g.register_step_with_input_arity("M", 1, 2);
        let c = g.register_step("C", 0);
        g.wire_to_slot(a, BranchIdx(0), m, 0);
        g.wire_to_slot(b, BranchIdx(0), m, 1);
        g.wire(m, BranchIdx(0), c);
        let transport = Arc::new(ByteBoundedQueue::<Sequenced<Sized8>>::new(1 << 20));
        let reorder_cap = ordered.then(|| {
            Arc::new(ReorderStage::with_max_overflow_bytes(Arc::clone(&transport) as _, 8))
                as Arc<dyn ReorderCapHandle>
        });
        let rq = RegisteredQueue {
            producer_step_name: "M",
            producer_step: m,
            branch: BranchIdx(0),
            handle: Arc::clone(&transport) as _,
            reorder_cap,
        };
        let plan = WakePlan::build(
            &g,
            &[StepKind::Exclusive, StepKind::Exclusive, StepKind::Detached, StepKind::Parallel],
            &[Some(0), Some(0), None, None],
            &[None, None, Some(DriverIdx(0)), None],
            &[],
            WakeEdges::byte_bounded(std::slice::from_ref(&rq)),
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        assert!(plan.is_gated(m), "ordered: {ordered}");
        assert_eq!(plan.gated_branches(m), 1);
        assert_eq!(plan.forward_target(m, BranchIdx(0)), WakeTarget::Pool);
        assert_eq!(plan.reverse_producers(c), vec![m]);
        assert!(plan.tracks(&rq), "the reverse edge keeps the edge tracked");
    }

    /// The producer's forward target in `on_flushed_delivers_the_push_s_forward_wake`.
    #[derive(Clone, Copy, Debug)]
    enum Fwd {
        Pool,
        Driver,
        Worker,
        Terminal,
        GatedPushed,
        GatedNotPushed,
        GatedPoolPushed,
        GatedPoolNotPushed,
    }

    /// `on_flushed` delivers the flushed push's forward wake to every target
    /// kind, through the gate, when the producer has no reverse edge and no
    /// same-thread consumer. The test thread is registered as P's own thread in
    /// `pool` and `terminal`, so a self-unpark there would add an `unparked`.
    #[rstest]
    #[case::pool(Fwd::Pool, WakeCounts { no_waiters: 1, ..WakeCounts::default() })]
    #[case::driver(Fwd::Driver, WakeCounts { unparked: 1, ..WakeCounts::default() })]
    #[case::worker(Fwd::Worker, WakeCounts { unparked: 1, ..WakeCounts::default() })]
    #[case::terminal(Fwd::Terminal, WakeCounts::default())]
    #[case::gated_pushed(Fwd::GatedPushed, WakeCounts { no_waiters: 1, ..WakeCounts::default() })]
    #[case::gated_not_pushed(Fwd::GatedNotPushed, WakeCounts { gated_off: 1, ..WakeCounts::default() })]
    #[case::gated_pool_pushed(Fwd::GatedPoolPushed, WakeCounts { no_waiters: 1, ..WakeCounts::default() })]
    #[case::gated_pool_not_pushed(Fwd::GatedPoolNotPushed, WakeCounts { gated_off: 1, ..WakeCounts::default() })]
    fn on_flushed_delivers_the_push_s_forward_wake(#[case] fwd: Fwd, #[case] expected: WakeCounts) {
        use StepKind::{Detached as D, Exclusive as E, Parallel as P, Serial as S};
        let mut g = ChainGraph::new();
        let p = g.register_step("P", 1);
        let c = g.register_step("C", usize::from(!matches!(fwd, Fwd::Driver)));
        if !matches!(fwd, Fwd::Terminal) {
            g.wire(p, BranchIdx(0), c);
        }
        let d0 = Some(DriverIdx(0));
        let gated = matches!(
            fwd,
            Fwd::GatedPushed | Fwd::GatedNotPushed | Fwd::GatedPoolPushed | Fwd::GatedPoolNotPushed
        );
        // (kinds, pinned, driver_of) per case.
        let (kinds, pinned, driver_of): (Vec<StepKind>, Vec<Option<usize>>, Vec<_>) = match fwd {
            Fwd::Pool => (vec![E, P, D], vec![Some(0), None, None], vec![None, None, d0]),
            Fwd::Driver => (vec![P, D], vec![None, None], vec![None, d0]),
            Fwd::Worker => (vec![P, S, D], vec![None, Some(1), None], vec![None, None, d0]),
            Fwd::Terminal => (vec![E, D], vec![Some(0), None], vec![None, d0]),
            Fwd::GatedPushed | Fwd::GatedNotPushed => {
                (vec![D, P], vec![None, None], vec![d0, None])
            }
            // A pool producer's byte-bounded branch is push-gated too.
            Fwd::GatedPoolPushed | Fwd::GatedPoolNotPushed => {
                (vec![P, P, D], vec![None, None, None], vec![None, None, d0])
            }
        };
        if kinds.len() == 3 {
            let t = g.register_step("Tail", 0);
            g.wire(c, BranchIdx(0), t);
        }
        let typed = Arc::new(ByteBoundedQueue::<u32>::new(1 << 20));
        let rq = RegisteredQueue {
            producer_step_name: "P",
            producer_step: p,
            branch: BranchIdx(0),
            handle: Arc::clone(&typed) as _,
            reorder_cap: None,
        };
        let registry = if gated { vec![rq] } else { vec![] };
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&registry),
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        assert_eq!(plan.is_gated(p), gated);
        if gated {
            typed.enable_tracking(
                EdgeTracking { clock: None, holders: Some(HolderSet::new(plan.n_slots())) },
                crate::queues::SEALED,
            );
        }
        let me = std::thread::current();
        // The consumer's thread is registered and, for a directed target,
        // armed for its timer park (a running target is not unparked).
        match fwd {
            Fwd::Pool | Fwd::Terminal => plan.register_worker(0, me),
            Fwd::Driver => {
                plan.register_driver(DriverIdx(0), me);
                assert!(plan.arm_driver(DriverIdx(0)));
            }
            Fwd::Worker => {
                plan.register_worker(1, me);
                assert!(plan.arm_direct(1));
            }
            Fwd::GatedPushed
            | Fwd::GatedNotPushed
            | Fwd::GatedPoolPushed
            | Fwd::GatedPoolNotPushed => {}
        }
        let mut snap = PushSnapshot::default();
        plan.snapshot_pushes(p, &mut snap);
        if matches!(fwd, Fwd::GatedPushed | Fwd::GatedPoolPushed) {
            typed.try_push(7).unwrap();
        }
        assert_eq!(plan.on_flushed(p, &snap), expected, "the flushed push's forward wake");
        std::thread::park_timeout(Duration::ZERO); // consume any token left on this thread
    }

    /// A flushed retry on a gated ordered branch whose item the reorder stage
    /// stashed (the transport is full) is judged by the push counter like a
    /// push: the stash insert advances it (`ItemQueue::note_stash_accepted`), so the
    /// forward wake is delivered, not gated off. `P(Detached d0) → C(Parallel)`
    /// at t2.
    #[test]
    fn on_flushed_counts_a_stash_insert_on_an_ordered_branch_as_a_push() {
        use crate::reorder::{ReorderStage, Sequenced};
        let mut g = ChainGraph::new();
        let p = g.register_step("P", 1);
        let c = g.register_step("C", 0);
        g.wire(p, BranchIdx(0), c);
        let transport = Arc::new(ByteBoundedQueue::<Sequenced<Sized8>>::new(8));
        let stage = Arc::new(ReorderStage::new(Arc::clone(&transport) as _));
        let rq = RegisteredQueue {
            producer_step_name: "P",
            producer_step: p,
            branch: BranchIdx(0),
            handle: Arc::clone(&transport) as _,
            reorder_cap: Some(Arc::clone(&stage) as _),
        };
        let plan = WakePlan::build(
            &g,
            &[StepKind::Detached, StepKind::Parallel],
            &[None, None],
            &[Some(DriverIdx(0)), None],
            &[],
            WakeEdges::byte_bounded(std::slice::from_ref(&rq)),
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        assert!(plan.is_gated(p));
        transport.enable_tracking(
            EdgeTracking { clock: None, holders: Some(HolderSet::new(plan.n_slots())) },
            crate::queues::SEALED,
        );
        stage.try_push(1, Sized8).unwrap(); // out of order: fills the transport
        let mut snap = PushSnapshot::default();
        plan.snapshot_pushes(p, &mut snap);
        stage.try_push(0, Sized8).unwrap(); // next_serial, transport full → stashed
        assert_eq!(
            plan.on_flushed(p, &snap),
            WakeCounts { no_waiters: 1, ..WakeCounts::default() },
            "the stash insert is a push"
        );
        let _ = crate::runtime::wake_slot::take_rejected();
    }

    /// `on_flushed` takes no holders: an idle-outcome dispatch made no room.
    /// `P(Parallel) →(one-item, tracked) C(Detached d0) →(count) Sink(Detached
    /// d0)`. C's only forward branch is same-driver, so its flushed wake is the
    /// self-unpark.
    #[test]
    fn on_flushed_takes_no_holders() {
        let mut g = ChainGraph::new();
        let p = g.register_step("P", 1);
        let c = g.register_step("C", 1);
        let sink = g.register_step("Sink", 0);
        g.wire(p, BranchIdx(0), c);
        g.wire(c, BranchIdx(0), sink);
        let d0 = Some(DriverIdx(0));
        let (rq, typed) = one_item_tracked(0, 3);
        let holder_only = vec![count_edge(1)];
        let plan = WakePlan::build(
            &g,
            &[StepKind::Parallel, StepKind::Detached, StepKind::Detached],
            &[None, None, None],
            &[None, d0, d0],
            &[],
            WakeEdges { byte_bounded: std::slice::from_ref(&rq), holder_only: &holder_only },
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        assert!(!plan.tracks_holder_only(&holder_only[0]), "C → Sink is same-driver");
        plan.register_driver(DriverIdx(0), std::thread::current());
        typed.try_push(Sized8).unwrap(); // full
        {
            let _slot = SlotGuard::enter(1);
            assert!(typed.try_push(Sized8).is_err(), "slot 1 recorded on (P, 0)");
        }
        let _ = crate::runtime::wake_slot::take_rejected();
        let r = plan.on_flushed(c, &PushSnapshot::default());
        assert_eq!(r.reverse, 0, "a flush takes no holders");
        assert_eq!(r.unparked, 1, "C's self-unpark (Sink shares its driver)");
        std::thread::park_timeout(Duration::ZERO);
        let r = plan.on_progress(c, &PushSnapshot::default());
        assert_eq!(r.reverse, 1, "the holder bit was still there for the pop's take");
        std::thread::park_timeout(Duration::ZERO);
    }

    /// Where a same-thread consumer runs in `on_flushed_repolls_a_same_thread_consumer`.
    #[derive(Clone, Copy, Debug)]
    enum SameThread {
        Driver,
        PinnedWorker,
        LoneWorker,
        DriverGatedNotPushed,
    }

    /// A flushed dispatch whose branch's consumer runs on the producer's own
    /// thread unparks that thread, so its next timer park returns at once and
    /// the next pass visits the consumer. A same-thread `Progress` delivers
    /// nothing (it relies on `did_work`). A gated branch that did not push gets
    /// no re-poll.
    #[rstest]
    #[case::same_driver(SameThread::Driver)]
    #[case::same_pinned_worker(SameThread::PinnedWorker)]
    #[case::lone_worker(SameThread::LoneWorker)]
    #[case::same_driver_gated_not_pushed(SameThread::DriverGatedNotPushed)]
    fn on_flushed_repolls_a_same_thread_consumer(#[case] case: SameThread) {
        use StepKind::{Detached as D, Exclusive as E, Parallel as P, Serial as S};
        let mut g = ChainGraph::new();
        let p = g.register_step("P", 1);
        let c = g.register_step("C", 1);
        let tail = g.register_step("Tail", 0);
        g.wire(p, BranchIdx(0), c);
        g.wire(c, BranchIdx(0), tail);
        let d0 = Some(DriverIdx(0));
        let d1 = Some(DriverIdx(1));
        let (kinds, pinned, driver_of, n_threads) = match case {
            SameThread::Driver | SameThread::DriverGatedNotPushed => {
                (vec![D, D, D], vec![None, None, None], vec![d0, d0, d1], 2)
            }
            SameThread::PinnedWorker => {
                (vec![E, S, D], vec![Some(0), Some(0), None], vec![None, None, d0], 2)
            }
            SameThread::LoneWorker => {
                (vec![P, P, D], vec![None, None, None], vec![None, None, d0], 1)
            }
        };
        let gated = matches!(case, SameThread::DriverGatedNotPushed);
        let (rq, _typed) = one_item_tracked(0, 4);
        let registry = if gated { vec![rq] } else { vec![] };
        let pool = (n_threads > 1).then(|| Arc::new(PoolEventCount::new(n_threads)));
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&registry),
            pool,
            n_threads,
        );
        assert_eq!(plan.forward_target(p, BranchIdx(0)), WakeTarget::None, "same thread");
        let me = std::thread::current();
        match case {
            SameThread::Driver | SameThread::DriverGatedNotPushed => {
                plan.register_driver(DriverIdx(0), me);
            }
            SameThread::PinnedWorker | SameThread::LoneWorker => plan.register_worker(0, me),
        }
        let mut snap = PushSnapshot::default();
        plan.snapshot_pushes(p, &mut snap);
        if gated {
            let r = plan.on_flushed(p, &snap);
            assert_eq!(r, WakeCounts { gated_off: 1, ..WakeCounts::default() });
            std::thread::park_timeout(Duration::ZERO);
            return;
        }
        assert_eq!(
            plan.on_progress(p, &snap),
            WakeCounts::default(),
            "Progress relies on did_work"
        );
        assert_eq!(plan.on_flushed(p, &snap).unparked, 1, "the self-unpark");
        let t = std::time::Instant::now();
        std::thread::park_timeout(Duration::from_secs(10));
        assert!(t.elapsed() < Duration::from_secs(5), "the token ends the next park at once");
    }

    /// A notify the concurrency ceiling keeps from a parked waiter is counted
    /// as `suppressed`, not as finding no waiter.
    #[test]
    fn a_ceiling_suppressed_notify_is_counted_apart() {
        let ec = Arc::new(PoolEventCount::new(2));
        ec.set_ceiling(1);
        let plan = super::tests_support::plan_with_pool(Arc::clone(&ec), 2);
        let key = ec.prepare_wait(); // one waiter: awake 1 >= ceiling 1
        let r = plan.on_progress(StepIdx(0), &PushSnapshot::default());
        ec.cancel_wait(key);
        assert_eq!(r, WakeCounts { suppressed: 1, ..WakeCounts::default() });
    }

    /// A `Pool` wake claims a cap-parked worker only when the consumer's cap has
    /// a free permit (or the consumer is uncapped): with the cap full, the bit
    /// stays armed and nobody is unparked (the release that frees a permit
    /// wakes the threads recorded on the cap). Once that cap has room, the
    /// next wake claims the still-armed worker. `Det (driver 0) → Par` at two
    /// workers; worker 1 is cap-parked.
    #[rstest]
    #[case::uncapped(None, 1)]
    #[case::cap_with_room(Some(false), 1)]
    #[case::full_cap(Some(true), 0)]
    fn a_cap_parked_claim_needs_the_consumer_s_cap_to_have_room(
        #[case] cap_full: Option<bool>,
        #[case] fallback: u8,
    ) {
        let mut g = ChainGraph::new();
        let d = g.register_step("Det", 1);
        let p = g.register_step("Par", 0);
        g.wire(d, BranchIdx(0), p);
        let cap = cap_full.map(|_| crate::PhaseCap::new("t", 1));
        let mut held = cap.as_ref().filter(|_| cap_full == Some(true)).map(|c| c.try_acquire_as(9));
        let plan = WakePlan::build(
            &g,
            &[StepKind::Detached, StepKind::Parallel],
            &[None, None],
            &[Some(DriverIdx(0)), None],
            &[None, cap.clone()],
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        plan.register_worker(1, std::thread::current());
        assert!(plan.arm_cap_parked(1));
        let r = plan.on_progress(d, &PushSnapshot::default());
        assert_eq!((r.no_waiters, r.fallback), (1, fallback), "{r:?}");
        if held.take().is_some() {
            let r = plan.on_progress(d, &PushSnapshot::default());
            assert_eq!(r.fallback, 1, "the cap has room now: the armed worker is claimed");
        }
        plan.disarm_cap_parked(1);
        std::thread::park_timeout(Duration::ZERO); // consume any token
    }

    /// `Det (driver 0) → Par` at two workers, `Par` capped by `cap`.
    fn plan_with_capped_consumer(cap: Option<Arc<crate::PhaseCap>>) -> (Arc<WakePlan>, StepIdx) {
        let mut g = ChainGraph::new();
        let d = g.register_step("Det", 1);
        let p = g.register_step("Par", 0);
        g.wire(d, BranchIdx(0), p);
        let plan = WakePlan::build(
            &g,
            &[StepKind::Detached, StepKind::Parallel],
            &[None, None],
            &[Some(DriverIdx(0)), None],
            &[None, cap],
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        (plan, p)
    }

    /// The occupancy rule a cap-parked re-poll skips by: an uncapped step and
    /// a cap with room are polled; a full cap, or one with room but a
    /// whole-cap reservation pending, is skipped.
    #[rstest]
    #[case::uncapped(None, false)]
    #[case::cap_with_room(Some(0), false)]
    #[case::full_cap(Some(2), true)]
    #[case::room_but_reserved(Some(1), true)]
    fn a_re_poll_skips_exactly_a_full_cap(#[case] held: Option<usize>, #[case] skip: bool) {
        let cap = held.map(|_| crate::PhaseCap::new("t", 2));
        let (plan, p) = plan_with_capped_consumer(cap.clone());
        let permits: Vec<_> = cap
            .iter()
            .flat_map(|c| (0..held.unwrap_or(0)).map(|i| c.try_acquire_as(10 + i)))
            .collect();
        let signal = crate::PipelineSignal::new();
        let reserver = cap.as_ref().filter(|_| held == Some(1)).map(|c| {
            c.bind_signal(&signal);
            let c = Arc::clone(c);
            let h = std::thread::spawn(move || c.acquire_whole().is_some());
            let deadline = std::time::Instant::now() + Duration::from_secs(10);
            while cap.as_ref().is_some_and(|c| c.pending_reservations() == 0) {
                assert!(std::time::Instant::now() < deadline, "reservation never registered");
                std::thread::yield_now();
            }
            h
        });
        assert_eq!(plan.skip_full_cap(p), skip);
        drop(permits);
        if let Some(h) = reserver {
            assert!(h.join().expect("reserver"), "the reservation is granted once the cap empties");
        }
    }

    /// A cap-parked worker that skips a full cap is recorded on it: the release
    /// that frees a permit wakes it.
    #[test]
    fn a_skipped_full_cap_wakes_the_skipping_worker_on_release() {
        let y = crate::PhaseCap::new("y", 1);
        let woken = crate::admission::bind_recording_waker(&y);
        let (plan, p) = plan_with_capped_consumer(Some(Arc::clone(&y)));
        let held = y.try_acquire_as(9).expect("an idle cap admits");
        {
            let _slot = SlotGuard::enter(1);
            assert!(plan.skip_full_cap(p), "the cap is full");
        }
        drop(held);
        assert_eq!(*woken.lock(), vec![Some(1)]);
    }

    /// A wake handed to a thread that does not spend it is forwarded. Slot 1
    /// (this thread) is refused, then slot 2; the holder's release claims slot
    /// 1 (the lower). Slot 1's next pass starts, owes that wake, and takes no
    /// permit, so its end forwards it to slot 2 while the permit is free.
    /// `admit_first`: the pass takes a permit (and releases it, waking slot 2
    /// itself), which pays the wake: with slot 3 also recorded, the end of the
    /// pass forwards nothing more.
    #[rstest]
    #[case::unspent(false, vec![Some(1), Some(2)])]
    #[case::admitted_first(true, vec![Some(1), Some(2)])]
    fn an_unspent_cap_wake_is_forwarded(
        #[case] admit_first: bool,
        #[case] expected: Vec<Option<usize>>,
    ) {
        let y = crate::PhaseCap::new("y", 1);
        let woken = crate::admission::bind_recording_waker(&y);
        let (plan, _p) = plan_with_capped_consumer(Some(Arc::clone(&y)));
        let held = y.try_acquire_as(9).expect("an idle cap admits");
        let _slot = SlotGuard::enter(1);
        plan.start_pass();
        assert!(y.try_acquire().is_none(), "slot 1 is refused");
        assert!(y.try_acquire_as(2).is_none(), "slot 2 is refused");
        if admit_first {
            assert!(y.try_acquire_as(3).is_none(), "slot 3 is refused");
        }
        plan.end_pass();
        drop(held);
        assert_eq!(*woken.lock(), vec![Some(1)], "the release claimed slot 1");
        plan.start_pass();
        if admit_first {
            drop(y.try_acquire().expect("the freed permit")); // its release wakes slot 2
        }
        plan.end_pass();
        assert_eq!(*woken.lock(), expected);
    }

    /// The cap-parked re-poll's skip, through the calling thread's own slot:
    /// a release that claims the skip's fresh bit before its re-check (the
    /// hook) hands the thread its wake, and the thread polls instead of
    /// skipping. The skip notes the record, so the next pass start owes the
    /// wake, and with no permit taken the pass end forwards it to slot 2.
    #[test]
    fn a_skip_whose_fresh_record_is_claimed_is_owed_and_forwarded() {
        let y: &'static Arc<crate::PhaseCap> = Box::leak(Box::new(crate::PhaseCap::new("y", 1)));
        let woken = crate::admission::bind_recording_waker(y);
        let (plan, p) = plan_with_capped_consumer(Some(Arc::clone(y)));
        let held = y.try_acquire_as(9).expect("an idle cap admits");
        assert!(y.try_acquire_as(2).is_none(), "slot 2 waits on the cap");
        let _slot = SlotGuard::enter(1);
        plan.start_pass();
        let release = std::cell::Cell::new(Some(held));
        let release: Box<dyn FnOnce()> = Box::new(move || drop(release.take()));
        crate::queues::test_hooks::after_record(release);
        assert!(!plan.skip_full_cap(p), "the freed permit: poll, do not skip");
        plan.end_pass();
        assert_eq!(*woken.lock(), vec![Some(1)], "the release claimed slot 1's fresh bit");
        plan.start_pass();
        plan.end_pass();
        assert_eq!(*woken.lock(), vec![Some(1), Some(2)], "slot 1 forwarded the unspent wake");
    }

    /// Only a wake the thread was handed is forwarded. Slots 2 and 3 are
    /// recorded and the release claims slot 2, so a permit is free and slot 3
    /// still waits; a pass on a thread that recorded nothing
    /// (`never_recorded`, slot 1), or only on the shared anonymous bit, which
    /// it cannot tell apart from a claim (`anonymous`, a slot past
    /// `REFUSED_SLOTS`), owes
    /// nothing, and its pass end forwards nothing.
    #[rstest]
    #[case::never_recorded(1)]
    #[case::anonymous(crate::admission::REFUSED_SLOTS + 45)]
    fn a_pass_with_no_claimed_record_forwards_nothing(#[case] slot: usize) {
        let y = crate::PhaseCap::new("y", 1);
        let woken = crate::admission::bind_recording_waker(&y);
        let (plan, _p) = plan_with_capped_consumer(Some(Arc::clone(&y)));
        let held = y.try_acquire_as(9).expect("an idle cap admits");
        let _slot = SlotGuard::enter(slot);
        plan.start_pass();
        if slot >= crate::admission::REFUSED_SLOTS {
            assert!(y.try_acquire().is_none(), "recorded on the anonymous bit");
        }
        assert!(y.try_acquire_as(2).is_none() && y.try_acquire_as(3).is_none());
        plan.end_pass();
        drop(held);
        assert_eq!(*woken.lock(), vec![Some(2)], "the release claimed slot 2");
        plan.start_pass();
        plan.end_pass();
        assert_eq!(*woken.lock(), vec![Some(2)], "slot 3 is not woken by a forward");
    }

    /// A plan over a Detached source and `n_caps` steps, each capped
    /// by its own cap; returns the caps in step order.
    fn plan_with_own_caps(n_caps: usize) -> Vec<Arc<crate::PhaseCap>> {
        let n = n_caps + 1;
        let mut g = ChainGraph::new();
        let mut prev = g.register_step("Det", 1);
        for _ in 1..n {
            let s = g.register_step("S", 1);
            g.wire(prev, BranchIdx(0), s);
            prev = s;
        }
        let mut kinds = vec![StepKind::Parallel; n];
        kinds[0] = StepKind::Detached;
        let mut driver_of = vec![None; n];
        driver_of[0] = Some(DriverIdx(0));
        let caps: Vec<_> = (0..n).map(|i| (i > 0).then(|| crate::PhaseCap::new("c", 1))).collect();
        let _ = WakePlan::build(
            &g,
            &kinds,
            &vec![None; n],
            &driver_of,
            &caps,
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        caps.into_iter().flatten().collect()
    }

    /// The plan takes exactly `MAX_PHASE_CAPS` distinct caps, the last at run
    /// index 63.
    #[test]
    fn a_plan_takes_max_phase_caps() {
        let caps = plan_with_own_caps(crate::MAX_PHASE_CAPS);
        assert_eq!(caps.last().and_then(|c| c.plan_bit()), Some(63));
    }

    /// The plan itself refuses more than `MAX_PHASE_CAPS` distinct caps, in
    /// every build profile (`PipelineBuilder::build` refuses them first).
    #[test]
    #[should_panic(expected = "a run may use at most 64 distinct phase caps")]
    fn a_plan_with_too_many_caps_panics() {
        let _ = plan_with_own_caps(crate::MAX_PHASE_CAPS + 1);
    }

    /// A cap has one run index at a time: a plan's build sets it, and a later
    /// plan sets its own.
    #[test]
    fn a_reused_cap_takes_each_run_s_index() {
        let cap = crate::PhaseCap::new("c", 1);
        let other = crate::PhaseCap::new("o", 1);
        let _first = plan_with_capped_consumer(Some(Arc::clone(&cap)));
        cap.bind_waker(Arc::new(|_: Option<usize>| {}));
        assert_eq!(cap.plan_bit(), Some(0), "the first run");
        let mut g = ChainGraph::new();
        let d = g.register_step("Det", 1);
        let a = g.register_step("A", 1);
        let b = g.register_step("B", 0);
        g.wire(d, BranchIdx(0), a);
        g.wire(a, BranchIdx(0), b);
        let _second = WakePlan::build(
            &g,
            &[StepKind::Detached, StepKind::Parallel, StepKind::Parallel],
            &[None, None, None],
            &[Some(DriverIdx(0)), None, None],
            &[None, Some(other), Some(Arc::clone(&cap))],
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        cap.bind_waker(Arc::new(|_: Option<usize>| {}));
        assert_eq!(cap.plan_bit(), Some(1), "the second run");
    }

    /// A pass start clears the thread's own refusal record on every cap the
    /// plan knows, and nothing else: another slot's bit stays, as does the
    /// record on a cap no step reports.
    #[test]
    fn a_pass_start_clears_only_the_thread_s_own_records() {
        let tracked = |name| {
            let c = crate::PhaseCap::new(name, 1);
            c.bind_waker(Arc::new(|_: Option<usize>| {}));
            c
        };
        let (cap_a, cap_b, other) = (tracked("a"), tracked("b"), tracked("other"));
        let mut graph = ChainGraph::new();
        let det = graph.register_step("Det", 1);
        let step_a = graph.register_step("A", 1);
        let step_b = graph.register_step("B", 0);
        graph.wire(det, BranchIdx(0), step_a);
        graph.wire(step_a, BranchIdx(0), step_b);
        let plan = WakePlan::build(
            &graph,
            &[StepKind::Detached, StepKind::Parallel, StepKind::Parallel],
            &[None, None, None],
            &[Some(DriverIdx(0)), None, None],
            &[None, Some(Arc::clone(&cap_a)), Some(Arc::clone(&cap_b))],
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(2))),
            2,
        );
        let held: Vec<_> =
            [&cap_a, &cap_b, &other].map(|c| c.try_acquire_as(9)).into_iter().collect();
        for c in [&cap_a, &cap_b, &other] {
            assert!(c.try_acquire_as(1).is_none() && c.try_acquire_as(0).is_none());
        }
        {
            let _slot = SlotGuard::enter(1);
            plan.start_pass();
        }
        assert_eq!(
            [&cap_a, &cap_b, &other].map(|c| (c.is_recorded(0), c.is_recorded(1))),
            [(true, false), (true, false), (true, true)]
        );
        drop(held);
    }

    /// Park `who` on its own thread for up to 10 s; returns once it has
    /// registered (through `register`), and the handle yields how long its
    /// park lasted.
    fn parked(
        register: impl FnOnce(std::thread::Thread) + Send + 'static,
    ) -> std::thread::JoinHandle<Duration> {
        let (tx, rx) = std::sync::mpsc::channel();
        let h = std::thread::spawn(move || {
            register(std::thread::current());
            tx.send(()).unwrap();
            let t = std::time::Instant::now();
            std::thread::park_timeout(Duration::from_secs(10));
            t.elapsed()
        });
        rx.recv().unwrap();
        h
    }

    /// A holder recorded in the anonymous slot (a slot outside the plan, or a
    /// thread with none) is woken by a broadcast: every registered driver and
    /// worker is unparked, and the event-count is notified.
    #[test]
    fn an_anonymous_holder_wake_is_a_broadcast() {
        let (plan, _graph) = tests_support::directed_plan_with_workers(2);
        let driver = {
            let plan = Arc::clone(&plan);
            parked(move |t| plan.register_driver(DriverIdx(0), t))
        };
        let worker = {
            let plan = Arc::clone(&plan);
            parked(move |t| plan.register_worker(1, t))
        };
        let ec = plan.pool().expect("two workers: an event-count");
        let k = ec.prepare_wait();
        assert!(plan.deliver_slot(plan.n_slots())); // the anonymous bit
        assert_eq!(ec.wait(k, Duration::ZERO), WaitOutcome::Woken, "the event-count is notified");
        for (who, h) in [("driver", driver), ("worker", worker)] {
            assert!(h.join().unwrap() < Duration::from_secs(5), "the {who} was unparked");
        }
    }

    /// The direct-park fallback reaches a worker past the first 64: worker 64
    /// arms (word 1, bit 0), and a `Pool` wake that finds no event-count waiter
    /// unparks exactly that registered thread.
    #[test]
    fn the_pool_fallback_reaches_a_worker_past_the_first_word() {
        let (plan, _graph) = tests_support::directed_plan_with_workers(65);
        let worker = {
            let plan = Arc::clone(&plan);
            parked(move |t| {
                plan.register_worker(64, t);
                assert!(plan.arm_direct(64));
            })
        };
        let r = plan.on_progress(StepIdx(0), &PushSnapshot::default());
        assert_eq!((r.no_waiters, r.fallback), (1, 1), "no waiter, so the armed worker: {r:?}");
        assert!(worker.join().unwrap() < Duration::from_secs(5), "worker 64 was unparked");
    }

    type ParIntoPar = (
        ChainGraph,
        Vec<StepKind>,
        Vec<Option<usize>>,
        Vec<Option<DriverIdx>>,
        Arc<ByteBoundedQueue<Sized8>>,
        Vec<RegisteredQueue>,
    );

    /// `Par(P) → Cons(P) → Tap(Detached, driver 0)`, n = 4. The Par→Cons edge is
    /// a one-item tracked queue (`one_item_tracked`); Cons→Tap is an ordinary
    /// `q`. Returns the plan inputs and the typed Par→Cons handle.
    fn parallel_into_parallel() -> ParIntoPar {
        let mut g = ChainGraph::new();
        let par = g.register_step("Par", 1);
        let cons = g.register_step("Cons", 1);
        let tap = g.register_step("Tap", 0);
        g.wire(par, BranchIdx(0), cons);
        g.wire(cons, BranchIdx(0), tap);
        let kinds = vec![StepKind::Parallel, StepKind::Parallel, StepKind::Detached];
        let pinned = vec![None, None, None];
        let driver_of = vec![None, None, Some(DriverIdx(0))];
        // 4 workers + 1 driver + the anonymous slot.
        let (rq0, typed) = one_item_tracked(0, 4 + 1 + 1);
        let queues = vec![rq0, q(1, 0)];
        (g, kinds, pinned, driver_of, typed, queues)
    }

    /// Every merge of `a` and `b` that preserves each sequence's own order.
    fn interleavings<'a>(a: &[&'a str], b: &[&'a str]) -> Vec<Vec<&'a str>> {
        if a.is_empty() {
            return vec![b.to_vec()];
        }
        if b.is_empty() {
            return vec![a.to_vec()];
        }
        let mut out = Vec::new();
        for mut rest in interleavings(&a[1..], b) {
            rest.insert(0, a[0]);
            out.push(rest);
        }
        for mut rest in interleavings(a, &b[1..]) {
            rest.insert(0, b[0]);
            out.push(rest);
        }
        out
    }

    /// Asserted at the top of every gate test: a fixture whose producer is not
    /// gated would make the gate tests pass for any implementation.
    fn assert_gated(plan: &WakePlan, step: StepIdx) {
        assert!(plan.is_gated(step), "fixture's producer must be gated");
    }

    /// Step names of [`native_input_sort_topology`], in `StepIdx` order.
    const NATIVE_INPUT_SORT_NAMES: [&str; 17] = [
        "PlanInputReads",
        "PreadInputSlices",
        "FrameBgzfBlocks",
        "ReadBlocks",
        "InflateToArena",
        "FindTemplateBoundaries",
        "ExtractTemplateKeys",
        "SortTemplateRun",
        "SpillGather",
        "SpillBlockCompress",
        "SpillWrite",
        "SpillReadPlanner",
        "PreadSpillSlices",
        "SortSpillDecompress",
        "SortMerge",
        "BgzfCompress",
        "WriteBgzfFile",
    ];

    /// A hand-written model of `fgumi sort`'s template-coordinate chain with
    /// native input, at t16: 17 steps, two driver groups (coord, io), the
    /// planner and framer pinned to the reader worker, every edge
    /// byte-bounded. Its step names, kinds, branch counts and per-branch
    /// gating mirror the table fgumi's `sort_chain_gating_matches_the_table`
    /// checks against the real chain (fgumi keeps its own copy, so this crate
    /// exports nothing chain-specific); its driver groups, pins and wiring are
    /// checked by nothing but the tests below.
    pub(super) fn native_input_sort_topology() -> Topology {
        use StepKind::{Detached as D, Parallel as P, Serial as S};
        let kinds = vec![S, P, S, D, P, D, P, D, D, P, D, S, P, P, D, P, D];
        let branches = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 0, 1, 1, 0];
        let mut g = ChainGraph::new();
        let idx: Vec<StepIdx> = NATIVE_INPUT_SORT_NAMES
            .iter()
            .zip(branches)
            .map(|(n, b)| g.register_step(n, b))
            .collect();
        for i in 0..11 {
            g.wire(idx[i], BranchIdx(0), idx[i + 1]);
        }
        g.wire(idx[11], BranchIdx(0), idx[12]);
        g.wire(idx[12], BranchIdx(0), idx[13]);
        g.wire(idx[11], BranchIdx(1), idx[14]);
        g.wire(idx[14], BranchIdx(0), idx[15]);
        g.wire(idx[15], BranchIdx(0), idx[16]);
        let (coord, io) = (Some(DriverIdx(0)), Some(DriverIdx(1)));
        let driver_of = vec![
            None, None, None, coord, None, coord, None, coord, coord, None, io, None, None, None,
            coord, None, io,
        ];
        let mut pinned = vec![None; 17];
        pinned[0] = Some(0);
        pinned[2] = Some(0);
        let mut queues: Vec<_> = (0..=11).map(|s| q(s, 0)).collect();
        queues.extend([q(11, 1), q(12, 0), q(14, 0), q(15, 0)]);
        (g, kinds, pinned, driver_of, queues)
    }

    /// Every producer with a byte-bounded output is gated, whatever its kind;
    /// the zero-output Parallel sink keeps its implicit ungated cascade.
    #[test]
    fn every_producer_with_a_byte_bounded_output_is_gated() {
        let (g, kinds, pinned, driver_of, queues) = native_input_sort_topology();
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            Some(Arc::new(PoolEventCount::new(16))),
            16,
        );
        // Every branch of a step with outputs is gated; `SortSpillDecompress`
        // (a zero-output Parallel step) and the `WriteBgzfFile` sink have none.
        let want = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 0, 1, 1, 0];
        for (i, want) in want.iter().enumerate() {
            assert_eq!(
                plan.gated_branches(StepIdx(i)),
                *want,
                "{} ({:?})",
                NATIVE_INPUT_SORT_NAMES[i],
                kinds[i]
            );
        }
        assert_eq!(
            plan.forward_branch_count_for_test(StepIdx(13)),
            1,
            "SortSpillDecompress: implicit cascade"
        );
        for rq in &queues {
            assert!(plan.tracks(rq), "every gated edge is tracked");
        }
    }

    /// Gate oracle at the boundary: an armed event-count key.
    #[test]
    fn gated_pool_producer_does_not_notify_without_a_push() {
        let (g, kinds, pinned, driver_of, typed, queues) = parallel_into_parallel();
        let pool = Arc::new(PoolEventCount::new(4));
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            Some(Arc::clone(&pool)),
            4,
        );
        assert_gated(&plan, StepIdx(0));
        let mut snap = PushSnapshot::default();
        plan.snapshot_pushes(StepIdx(0), &mut snap);
        let key = pool.prepare_wait();
        let r = plan.on_progress(StepIdx(0), &snap);
        assert_eq!(pool.wait(key, Duration::ZERO), WaitOutcome::TimedOut, "no push → no notify");
        assert_eq!(r.gated_off, 1);
        typed.try_push(Sized8).unwrap();
        let key = pool.prepare_wait();
        let _ = plan.on_progress(StepIdx(0), &snap);
        assert_eq!(pool.wait(key, Duration::ZERO), WaitOutcome::Woken, "a push → one notify");
    }

    /// Two clones A (pushes) and B (does not) of one producer, and a consumer
    /// that arms at every possible point, then re-polls and checks
    /// `wait(key, 0)`. In every one of the 10 producer orders × 6 arm points,
    /// the re-poll saw the item or the key reads `Woken`.
    #[test]
    fn concurrent_clone_orders_never_lose_a_wake() {
        let orders = interleavings(&["A.snap", "A.push", "A.prog"], &["B.snap", "B.prog"]);
        assert_eq!(orders.len(), 10);
        for order in &orders {
            for arm_at in 0..=order.len() {
                let (g, kinds, pinned, driver_of, typed, queues) = parallel_into_parallel();
                let pool = Arc::new(PoolEventCount::new(4));
                let plan = WakePlan::build(
                    &g,
                    &kinds,
                    &pinned,
                    &driver_of,
                    &[],
                    WakeEdges::byte_bounded(&queues),
                    Some(Arc::clone(&pool)),
                    4,
                );
                assert_gated(&plan, StepIdx(0));
                let (mut a, mut b) = (PushSnapshot::default(), PushSnapshot::default());
                let (mut key, mut saw_item) = (None, false);
                for (i, ev) in order.iter().enumerate() {
                    if i == arm_at {
                        key = Some(pool.prepare_wait());
                        saw_item = !typed.is_empty();
                    }
                    match *ev {
                        "A.snap" => plan.snapshot_pushes(StepIdx(0), &mut a),
                        "B.snap" => plan.snapshot_pushes(StepIdx(0), &mut b),
                        "A.push" => typed.try_push(Sized8).unwrap(),
                        "A.prog" => {
                            let _ = plan.on_progress(StepIdx(0), &a);
                        }
                        "B.prog" => {
                            let _ = plan.on_progress(StepIdx(0), &b);
                        }
                        _ => unreachable!(),
                    }
                }
                if key.is_none() {
                    key = Some(pool.prepare_wait());
                    saw_item = !typed.is_empty();
                }
                let woken = pool.wait(key.unwrap(), Duration::ZERO) == WaitOutcome::Woken;
                assert!(
                    saw_item || woken,
                    "order {order:?}, consumer armed at {arm_at}: wake lost"
                );
            }
        }
    }

    type DirectedCountBounded =
        (ChainGraph, Vec<StepKind>, Vec<Option<usize>>, Vec<Option<DriverIdx>>, Vec<ParallelHosts>);

    /// `Src(Exclusive w0) → Par(P) → Sink(Serial, None) → Tap(Detached)`; every
    /// edge `CountBounded`, so no byte-bounded registry entry exists: a routed
    /// plan with no gates.
    fn directed_count_bounded_chain(n: usize, hosts_exclude_reader: bool) -> DirectedCountBounded {
        let mut g = ChainGraph::new();
        let src = g.register_step("Src", 1);
        let par = g.register_step("Par", 1);
        let sink = g.register_step("Sink", 1);
        let tap = g.register_step("Tap", 0);
        g.wire(src, BranchIdx(0), par);
        g.wire(par, BranchIdx(0), sink);
        g.wire(sink, BranchIdx(0), tap);
        let kinds =
            vec![StepKind::Exclusive, StepKind::Parallel, StepKind::Serial, StepKind::Detached];
        let pinned = vec![Some(0), None, None, None];
        let driver_of = vec![None, None, None, Some(DriverIdx(0))];
        let workers: Vec<usize> =
            if hosts_exclude_reader { (1..n).collect() } else { (0..n).collect() };
        let hosts = vec![
            ParallelHosts::default(),
            ParallelHosts { workers, driver: None },
            ParallelHosts::default(),
            ParallelHosts::default(),
        ];
        (g, kinds, pinned, driver_of, hosts)
    }

    /// Under the widened gate (every producer kind), a count-bounded branch
    /// stays ungated — it has no push counter — whatever the producer's role,
    /// a Pool producer's `Progress` on it still wakes its consumer, and every
    /// cross-thread producer gets a reverse wake on it (count-bounded holder
    /// wakes).
    #[test]
    fn count_bounded_branch_stays_ungated_with_a_reverse_edge() {
        let (g, kinds, pinned, driver_of, hosts) = directed_count_bounded_chain(4, false);
        let holder_only: Vec<RegisteredHolderOnlyQueue> = (0..3).map(count_edge).collect();
        let pool = Arc::new(PoolEventCount::new(4));
        let plan = WakePlan::build_with_hosts(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges { byte_bounded: &[], holder_only: &holder_only },
            Some(Arc::clone(&pool)),
            4,
            &hosts,
        );
        for p in 0..3 {
            assert!(!plan.is_gated(StepIdx(p)), "a count-bounded branch is never gated");
        }
        assert_eq!(plan.reverse_producers(StepIdx(1)), vec![StepIdx(0)]);
        assert_eq!(plan.reverse_producers(StepIdx(2)), vec![StepIdx(1)]);
        assert_eq!(plan.reverse_producers(StepIdx(3)), vec![StepIdx(2)]);
        let key = pool.prepare_wait();
        let _ = plan.on_progress(StepIdx(1), &PushSnapshot::default());
        assert_eq!(pool.wait(key, Duration::ZERO), WaitOutcome::Woken, "Par's ungated Pool wake");
    }

    /// The stashed `next_serial` reaches a consumer armed before the push: an
    /// ordered branch is gated, and the gate sees the stash insert as a push.
    #[test]
    fn stashed_must_accept_push_wakes_the_consumer() {
        use crate::reorder::{ReorderStage, Sequenced};

        let mut g = ChainGraph::new();
        let par = g.register_step("Par", 1);
        let cons = g.register_step("Cons", 1);
        let tap = g.register_step("Tap", 0);
        g.wire(par, BranchIdx(0), cons);
        g.wire(cons, BranchIdx(0), tap);
        let transport = Arc::new(ByteBoundedQueue::<Sequenced<Sized8>>::new(8));
        let stage = Arc::new(ReorderStage::new(Arc::clone(&transport) as _));
        let rq = RegisteredQueue {
            producer_step_name: "Par",
            producer_step: par,
            branch: BranchIdx(0),
            handle: Arc::clone(&transport) as _,
            reorder_cap: Some(Arc::clone(&stage) as _),
        };
        let pool = Arc::new(PoolEventCount::new(4));
        let plan = WakePlan::build(
            &g,
            &[StepKind::Parallel, StepKind::Parallel, StepKind::Detached],
            &[None, None, None],
            &[None, None, Some(DriverIdx(0))],
            &[],
            WakeEdges::byte_bounded(&[rq, q(1, 0)]),
            Some(Arc::clone(&pool)),
            4,
        );
        assert_gated(&plan, par);
        transport.enable_tracking(
            EdgeTracking { clock: None, holders: Some(HolderSet::new(plan.n_slots())) },
            crate::queues::SEALED,
        );
        stage.try_push(1, Sized8).unwrap(); // fills the transport (out of order)
        let mut snap = PushSnapshot::default();
        plan.snapshot_pushes(par, &mut snap);
        let key = pool.prepare_wait();
        stage.try_push(0, Sized8).unwrap(); // next_serial → stashed
        let _ = plan.on_progress(par, &snap);
        assert_eq!(pool.wait(key, Duration::ZERO), WaitOutcome::Woken, "stashed push gated off");
    }

    type PlainChain = (
        ChainGraph,
        Vec<StepKind>,
        Vec<Option<usize>>,
        Vec<Option<DriverIdx>>,
        Vec<Arc<ByteBoundedQueue<Sized8>>>,
        Vec<RegisteredQueue>,
    );

    /// `Src(Serial, Reader w0) → Par(P) → Grouper(Serial, Worker(1)) → Par2(P)
    /// → Sink(Serial, Writer w3)`, n = 4, no Detached step. Every edge is a
    /// one-item `Sized8` queue; returns the typed handles (index = producer step).
    fn plain_chain_with_grouper() -> PlainChain {
        let mut g = ChainGraph::new();
        let names = ["Src", "Par", "Grouper", "Par2", "Sink"];
        let idx: Vec<StepIdx> = names
            .iter()
            .enumerate()
            .map(|(i, n)| g.register_step(n, usize::from(i != 4)))
            .collect();
        for i in 0..4 {
            g.wire(idx[i], BranchIdx(0), idx[i + 1]);
        }
        let kinds = vec![
            StepKind::Serial,
            StepKind::Parallel,
            StepKind::Serial,
            StepKind::Parallel,
            StepKind::Serial,
        ];
        let pinned = vec![Some(0), None, Some(1), None, Some(3)];
        let driver_of = vec![None; 5];
        let mut typed = Vec::new();
        let mut queues = Vec::new();
        for (p, name) in names.iter().take(4).enumerate() {
            let t = Arc::new(ByteBoundedQueue::<Sized8>::new(1));
            queues.push(RegisteredQueue {
                producer_step_name: name,
                producer_step: StepIdx(p),
                branch: BranchIdx(0),
                handle: Arc::clone(&t) as _,
                reorder_cap: None,
            });
            typed.push(t);
        }
        (g, kinds, pinned, driver_of, typed, queues)
    }

    /// Every chain gets tables: the plain-chain shape (no Detached step) routes
    /// pinned consumers to their worker and unpinned ones to the pool, and
    /// tracks every byte-bounded edge. A rendering check of `build`'s rules;
    /// delivery is `plain_chain_progress_unparks_the_pinned_consumer`.
    #[test]
    fn plain_chain_builds_tables() {
        let (g, kinds, pinned, driver_of, _typed, queues) = plain_chain_with_grouper();
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            Some(Arc::new(PoolEventCount::new(4))),
            4,
        );
        assert_eq!(plan.forward_target(StepIdx(0), BranchIdx(0)), WakeTarget::Pool);
        assert_eq!(plan.forward_target(StepIdx(1), BranchIdx(0)), WakeTarget::Worker(1));
        assert_eq!(plan.forward_target(StepIdx(3), BranchIdx(0)), WakeTarget::Worker(3));
        assert!(
            queues.iter().all(|rq| plan.tracks(rq)),
            "every byte-bounded edge is tracked (S2 gate + holder sets)"
        );
        assert!(
            plan.dag_lines(&g).iter().all(|l| !l.contains("Legacy")),
            "{:?}",
            plan.dag_lines(&g)
        );
    }

    /// The fused path's plan: no tables, no pool, nothing to request.
    #[test]
    fn standalone_plan_has_no_tables_and_requests_are_unavailable() {
        let plan = WakePlan::standalone();
        assert_eq!(plan.forward_target(StepIdx(0), BranchIdx(0)), WakeTarget::None);
        assert_eq!(plan.gated_branches(StepIdx(0)), 0);
        assert!(plan.reverse_producers(StepIdx(0)).is_empty());
        assert_eq!(plan.n_slots(), 0);
        assert_eq!(
            PoolHandle::new(&plan, None, StepIdx(0)).request_worker(),
            PoolRequest::Unavailable
        );
        let _slot = SlotGuard::enter(0);
        assert_eq!(
            PoolHandle::new(&plan, None, StepIdx(0)).request_worker(),
            PoolRequest::Unavailable,
            "even from slot 0: the standalone plan has no worker"
        );
    }

    /// A delivery to a pinned consumer of a plain chain reaches a real parked
    /// thread: the helper registers as worker 1 from its own thread, parks on a
    /// 10 s timer, and must return under a 5 s watchdog once the Par step's
    /// progress (one push onto `Par → Grouper`) is delivered.
    #[test]
    fn plain_chain_progress_unparks_the_pinned_consumer() {
        let (g, kinds, pinned, driver_of, typed, queues) = plain_chain_with_grouper();
        let plan = WakePlan::build(
            &g,
            &kinds,
            &pinned,
            &driver_of,
            &[],
            WakeEdges::byte_bounded(&queues),
            Some(Arc::new(PoolEventCount::new(4))),
            4,
        );
        // Mirror `Pipeline::run`: a holder set on every edge the plan tracks.
        for rq in &queues {
            if plan.tracks(rq) {
                rq.handle.enable_tracking(
                    EdgeTracking { clock: None, holders: Some(HolderSet::new(plan.n_slots())) },
                    crate::queues::SEALED,
                );
            }
        }
        assert!(plan.is_gated(StepIdx(1)), "fixture: Par → Grouper(Worker(1)) must be gated");
        let (ready_tx, ready_rx) = std::sync::mpsc::channel();
        let (done_tx, done_rx) = std::sync::mpsc::channel();
        let p = Arc::clone(&plan);
        std::thread::spawn(move || {
            p.register_worker(1, std::thread::current());
            assert!(p.arm_direct(1), "armed for its timer park");
            ready_tx.send(()).unwrap();
            let t0 = std::time::Instant::now();
            std::thread::park_timeout(Duration::from_secs(10));
            done_tx.send(t0.elapsed()).unwrap();
        });
        ready_rx.recv().unwrap();
        // Let it park; an earlier unpark is still caught (the token persists).
        std::thread::sleep(Duration::from_millis(20));
        let mut snap = PushSnapshot::default();
        plan.snapshot_pushes(StepIdx(1), &mut snap);
        typed[1].try_push(Sized8).unwrap();
        let r = plan.on_progress(StepIdx(1), &snap);
        assert_eq!(r.unparked, 1, "{r:?}");
        let woke = done_rx.recv_timeout(Duration::from_secs(5)).expect("worker 1 not unparked");
        assert!(woke < Duration::from_secs(5));
    }
}

/// Fixtures shared with other modules' tests.
#[cfg(test)]
pub(crate) mod tests_support {
    use std::sync::Arc;

    use super::{DriverIdx, WakeEdges, WakePlan};
    use crate::runtime::event_count::PoolEventCount;
    use crate::step::StepKind;
    use crate::topology::{BranchIdx, ChainGraph};

    /// A pool handle over a process-wide standalone plan (no pool), for tests
    /// that build an `ErasedStepCtx` by hand.
    pub(crate) fn standalone_pool_handle() -> super::PoolHandle<'static> {
        static PLAN: std::sync::OnceLock<Arc<WakePlan>> = std::sync::OnceLock::new();
        super::PoolHandle::new(
            PLAN.get_or_init(WakePlan::standalone),
            None,
            crate::topology::StepIdx(0),
        )
    }

    /// A real two-step plan over `pool` (built for `n` workers) for tests
    /// that dispatch step 0 by hand: `Par(Parallel) → Par2(Parallel)`, no edge
    /// registered, so step 0's one branch targets `Pool`, is ungated and has
    /// no reverse edge — every `on_progress(StepIdx(0), ..)` is exactly one
    /// `notify_one`.
    pub(crate) fn plan_with_pool(pool: Arc<PoolEventCount>, n: usize) -> Arc<WakePlan> {
        let mut g = ChainGraph::new();
        let a = g.register_step("Par", 1);
        let b = g.register_step("Par2", 0);
        g.wire(a, BranchIdx(0), b);
        WakePlan::build(
            &g,
            &[StepKind::Parallel, StepKind::Parallel],
            &[None, None],
            &[None, None],
            &[],
            WakeEdges::NONE,
            Some(pool),
            n,
        )
    }

    /// The plan of a plain chain (no Detached step) at `n` threads:
    /// `A(Serial) → B(Serial)`, a pool event-count only when `n > 1`. At one
    /// thread it is the lone worker's plan, with nothing else in it.
    pub(crate) fn plain_plan(n: usize) -> Arc<WakePlan> {
        let mut g = ChainGraph::new();
        let a = g.register_step("A", 1);
        let b = g.register_step("B", 0);
        g.wire(a, BranchIdx(0), b);
        WakePlan::build(
            &g,
            &[StepKind::Serial, StepKind::Serial],
            &[None, None],
            &[None, None],
            &[],
            WakeEdges::NONE,
            (n > 1).then(|| Arc::new(PoolEventCount::new(n))),
            n,
        )
    }

    /// A plan for `Det (driver 0) → Par` with `n` pool workers and a
    /// pool event-count.
    pub(crate) fn directed_plan_with_workers(n: usize) -> (Arc<WakePlan>, ChainGraph) {
        let mut g = ChainGraph::new();
        let d = g.register_step("Det", 1);
        let p = g.register_step("Par", 0);
        g.wire(d, BranchIdx(0), p);
        let plan = WakePlan::build(
            &g,
            &[StepKind::Detached, StepKind::Parallel],
            &[None, None],
            &[Some(DriverIdx(0)), None],
            &[],
            WakeEdges::NONE,
            Some(Arc::new(PoolEventCount::new(n))),
            n,
        );
        (plan, g)
    }
}

/// End-to-end liveness of the directed wakes. Each test gives **only the thread
/// whose wake it checks** a 10 s timer (`PipelineConfig::test_backoff`), so that
/// wake — and nothing else — can finish the run under the watchdog.
#[cfg(test)]
mod e2e {
    use std::io;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering as AO};
    use std::time::{Duration, Instant};

    use rstest::rstest;

    use crate::builder::{Pipeline, PipelineConfig};
    use crate::handles::Unpushed;
    use crate::outputs::Single;
    use crate::queues::QueueSpec;
    use crate::reorder::BranchOrdering;
    use crate::runtime::driver::test_hooks;
    use crate::runtime::event_count::WaitOutcome;
    use crate::runtime::stats::StatsSnapshot;
    use crate::runtime::wake_slot;
    use crate::runtime::worker_core::{TestBackoff, TestBackoffTarget as T};
    use crate::step::{Affinity, DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile};

    const TEN_S: u64 = 10_000_000;
    const EDGE: QueueSpec = QueueSpec::ByteBounded { limit_bytes: 1 << 20 };
    const ONE_BLOB: QueueSpec = QueueSpec::ByteBounded { limit_bytes: 64 };
    const ONE_COUNT: QueueSpec = QueueSpec::CountBounded { capacity: 1 };

    fn slow(target: T) -> Vec<TestBackoff> {
        vec![TestBackoff { target, us: TEN_S }]
    }

    /// 64-byte items so a 64-byte-limit edge holds exactly one.
    #[derive(Clone, Copy, Debug)]
    struct Blob(#[allow(dead_code)] u32);
    impl crate::item::HeapSize for Blob {
        fn heap_size(&self) -> usize {
            64
        }
    }

    /// Run `f` on its own thread; fail if it has not returned in `secs`.
    pub(super) fn run_under_watchdog(context: &str, secs: u64, f: impl FnOnce() + Send + 'static) {
        let (tx, rx) = std::sync::mpsc::channel::<()>();
        let handle = std::thread::spawn(move || {
            f();
            let _ = tx.send(());
        });
        match rx.recv_timeout(Duration::from_secs(secs)) {
            Ok(()) => handle.join().expect("run thread panicked"),
            Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => {
                std::panic::resume_unwind(handle.join().expect_err("closure panicked"));
            }
            Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
                panic!("{context}: WEDGED — only a 10 s timer would have finished it")
            }
        }
    }

    /// Source of `n` items. `delay` (from its first call) holds the first item
    /// back so every consumer parks before any item exists; `kind`/`group`
    /// choose where it runs (Exclusive on a pool worker, or Detached on a driver).
    struct Source {
        n: u32,
        next: u32,
        delay: Duration,
        started: Option<Instant>,
        kind: StepKind,
        group: DetachedGroup,
        affinity: Affinity,
        held: Option<Blob>,
        /// When set, `Finished` waits until the sink has seen every item, so a
        /// `Finished` broadcast cannot stand in for the wakes under test.
        finish_after_sink: Option<Arc<AtomicU64>>,
        /// The output edge (`EDGE` unless `with_edge`).
        edge: QueueSpec,
        /// When set, every push and re-push is counted (see [`RefusalCounts`]).
        counts: Option<RefusalCounts>,
    }
    impl Source {
        fn exclusive(n: u32, delay_ms: u64) -> Self {
            Self {
                n,
                next: 0,
                delay: Duration::from_millis(delay_ms),
                started: None,
                kind: StepKind::Exclusive,
                group: DetachedGroup::PerStep,
                affinity: Affinity::None,
                held: None,
                finish_after_sink: None,
                edge: EDGE,
                counts: None,
            }
        }
        fn with_edge(self, edge: QueueSpec) -> Self {
            Self { edge, ..self }
        }
        fn with_counts(self, counts: &RefusalCounts) -> Self {
            Self { counts: Some(counts.clone()), ..self }
        }
        fn push(&self, ctx: &StepCtx<'_, Self>, b: Blob) -> Result<(), Unpushed<Blob>> {
            match &self.counts {
                Some(c) => c.count(Result::is_err, || ctx.outputs.push(b)),
                None => ctx.outputs.push(b),
            }
        }
        /// A Serial source pinned by `affinity` (`Reader` → worker 0).
        fn serial_with_affinity(n: u32, delay_ms: u64, affinity: Affinity) -> Self {
            Self { kind: StepKind::Serial, affinity, ..Self::exclusive(n, delay_ms) }
        }
        /// Like [`Self::exclusive`], but `Finished` only once `seen` reaches `n`.
        fn exclusive_finishing_last(n: u32, delay_ms: u64, seen: &Arc<AtomicU64>) -> Self {
            Self { finish_after_sink: Some(Arc::clone(seen)), ..Self::exclusive(n, delay_ms) }
        }
        fn detached(n: u32, delay_ms: u64, group: &'static str) -> Self {
            Self {
                kind: StepKind::Detached,
                group: DetachedGroup::Shared(group),
                ..Self::exclusive(n, delay_ms)
            }
        }
    }
    impl Step for Source {
        type Input = ();
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Source",
                kind: self.kind,
                sticky: false,
                output_queues: vec![self.edge],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn affinity(&self) -> Affinity {
            self.affinity
        }
        fn detached_group(&self) -> DetachedGroup {
            self.group
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let started = *self.started.get_or_insert_with(Instant::now);
            if started.elapsed() < self.delay {
                std::thread::sleep(Duration::from_millis(1));
                return Ok(StepOutcome::NoProgress);
            }
            if let Some(b) = self.held.take() {
                return match self.push(ctx, b) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(u) => {
                        self.held = Some(u.into_item());
                        Ok(StepOutcome::NoProgress)
                    }
                };
            }
            if self.next == self.n {
                if let Some(seen) = &self.finish_after_sink
                    && seen.load(AO::Relaxed) < u64::from(self.n)
                {
                    std::thread::sleep(Duration::from_millis(1));
                    return Ok(StepOutcome::NoProgress);
                }
                return Ok(StepOutcome::Finished);
            }
            let b = Blob(self.next);
            self.next += 1;
            match self.push(ctx, b) {
                Ok(()) => Ok(StepOutcome::Progress),
                Err(u) => {
                    self.held = Some(u.into_item());
                    Ok(StepOutcome::Progress)
                }
            }
        }
    }

    /// Shared across every worker copy of the step that carries it: how many
    /// of its pushes and held retries were refused, and how many of those
    /// refusals marked the thread as holding (`note_held`), so a holder can be
    /// shown to idle where a direct unpark reaches it.
    #[derive(Clone, Default)]
    struct RefusalCounts {
        refusals: Arc<AtomicU64>,
        latched: Arc<AtomicU64>,
    }
    impl RefusalCounts {
        fn count<R>(&self, refused: impl Fn(&R) -> bool, push: impl FnOnce() -> R) -> R {
            // REJECTED stays set until the thread idles, so a peek after a later
            // refusal in the same pass would read an earlier refusal's flag. Take
            // it first and restore it after, so the peek sees only what this push
            // did. REJECTED is read only at idle, never inside a dispatch, so the
            // take and restore cannot change the loop's behaviour.
            let before = wake_slot::take_rejected();
            let r = push();
            if refused(&r) {
                self.refusals.fetch_add(1, AO::Relaxed);
                if wake_slot::rejected_pending() {
                    self.latched.fetch_add(1, AO::Relaxed);
                }
            }
            if before {
                wake_slot::note_held();
            }
            r
        }
        fn refusals(&self) -> u64 {
            self.refusals.load(AO::Relaxed)
        }
        fn latched(&self) -> u64 {
            self.latched.load(AO::Relaxed)
        }
    }

    /// Mirrors `Process::try_run` (`src/lib/pipeline/steps/process.rs`) line
    /// for line: the held retry first, `Contention` on a refused retry, then the
    /// input pop even after a flushed retry (an empty, undrained input then
    /// returns `NoProgress`, whose forward wake `on_flushed` delivers), and
    /// `Progress` on a new hold. Its output is one count-bounded slot.
    struct ProcessShaped {
        counts: RefusalCounts,
        held: crate::held::HeldSlot<Unpushed<Blob>>,
    }
    impl Step for ProcessShaped {
        type Input = Blob;
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ProcessShaped",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![ONE_COUNT],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(unpushed) = self.held.take() {
                match self.counts.count(Result::is_err, || ctx.outputs.retry(unpushed)) {
                    Ok(()) => {}
                    Err(again) => {
                        self.held.put(again);
                        return Ok(StepOutcome::Contention);
                    }
                }
            }
            let Some(item) = ctx.input.pop() else {
                if ctx.input.is_drained() {
                    return Ok(StepOutcome::Finished);
                }
                return Ok(StepOutcome::NoProgress);
            };
            match self.counts.count(Result::is_err, || ctx.outputs.push(item)) {
                Ok(()) => Ok(StepOutcome::Progress),
                Err(unpushed) => {
                    self.held.put(unpushed);
                    Ok(StepOutcome::Progress)
                }
            }
        }
        fn new_worker_copy(&self) -> Self {
            Self { counts: self.counts.clone(), held: crate::held::HeldSlot::new() }
        }
    }

    /// Pass-through with configurable kind/affinity/group/edge and an optional per-item sleep.
    #[derive(Clone)]
    struct Pass {
        name: &'static str,
        kind: StepKind,
        affinity: Affinity,
        group: DetachedGroup,
        edge: QueueSpec,
        work: Duration,
        held: Option<Blob>,
    }
    impl Pass {
        fn new(name: &'static str, kind: StepKind) -> Self {
            Self {
                name,
                kind,
                affinity: Affinity::None,
                group: DetachedGroup::PerStep,
                edge: EDGE,
                work: Duration::ZERO,
                held: None,
            }
        }
    }
    impl Step for Pass {
        type Input = Blob;
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: self.name,
                kind: self.kind,
                sticky: false,
                output_queues: vec![self.edge],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn affinity(&self) -> Affinity {
            self.affinity
        }
        fn detached_group(&self) -> DetachedGroup {
            self.group
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(b) = self.held.take() {
                return match ctx.outputs.push(b) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(u) => {
                        self.held = Some(u.into_item());
                        Ok(StepOutcome::NoProgress)
                    }
                };
            }
            match ctx.input.pop() {
                Some(b) => {
                    if !self.work.is_zero() {
                        std::thread::sleep(self.work);
                    }
                    match ctx.outputs.push(b) {
                        Ok(()) => Ok(StepOutcome::Progress),
                        Err(u) => {
                            self.held = Some(u.into_item());
                            Ok(StepOutcome::Progress)
                        }
                    }
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
        fn new_worker_copy(&self) -> Self {
            Self { held: None, ..self.clone() }
        }
    }

    /// Counting sink of configurable kind/affinity/group, with an optional
    /// per-item sleep and an optional count of its empty, undrained polls.
    struct Sink {
        kind: StepKind,
        affinity: Affinity,
        group: DetachedGroup,
        seen: Arc<AtomicU64>,
        work: Duration,
        empty_polls: Option<Arc<AtomicU64>>,
    }
    impl Sink {
        fn new(kind: StepKind, seen: &Arc<AtomicU64>) -> Self {
            Self {
                kind,
                affinity: Affinity::None,
                group: DetachedGroup::PerStep,
                seen: Arc::clone(seen),
                work: Duration::ZERO,
                empty_polls: None,
            }
        }
        fn counting_empty_polls(self, polls: &Arc<AtomicU64>) -> Self {
            Self { empty_polls: Some(Arc::clone(polls)), ..self }
        }
    }
    impl Step for Sink {
        type Input = Blob;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Sink",
                kind: self.kind,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn affinity(&self) -> Affinity {
            self.affinity
        }
        fn detached_group(&self) -> DetachedGroup {
            self.group
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(_) => {
                    if !self.work.is_zero() {
                        std::thread::sleep(self.work);
                    }
                    self.seen.fetch_add(1, AO::Relaxed);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => {
                    if let Some(c) = &self.empty_polls {
                        c.fetch_add(1, AO::Relaxed);
                    }
                    Ok(StepOutcome::NoProgress)
                }
            }
        }
    }

    fn run(
        threads: usize,
        test_backoff: Vec<TestBackoff>,
        build: impl FnOnce(&crate::builder::PipelineBuilder),
    ) -> StatsSnapshot {
        let builder = Pipeline::builder();
        build(&builder);
        let pipeline = builder.build().expect("build");
        let stats = pipeline.stats();
        pipeline
            .run(
                PipelineConfig { threads, test_backoff, ..Default::default() }
                    .with_stats(Arc::clone(&stats)),
            )
            .expect("run");
        stats.snapshot()
    }

    /// Parallel producer on the pool → Detached consumer on a driver with a
    /// 10 s timer (only that driver). 1,000 items complete under a 5 s watchdog
    /// only because every push unparks the driver: the source finishes only
    /// after the sink has seen every item, so no `Finished` broadcast can
    /// deliver them instead.
    #[test]
    fn directed_driver_wake_is_not_the_timer() {
        let seen = Arc::new(AtomicU64::new(0));
        let s2 = Arc::clone(&seen);
        run_under_watchdog("driver unpark", 5, move || {
            let _ = run(2, slow(T::Driver(0)), |b| {
                b.chain(Source::exclusive_finishing_last(1_000, 100, &s2))
                    .chain(Pass::new("Par", StepKind::Parallel))
                    .chain(Pass::new("Det", StepKind::Detached))
                    .chain(Sink::new(StepKind::Serial, &s2))
                    .into_sink_marker();
            });
        });
        assert_eq!(seen.load(AO::Relaxed), 1_000);
    }

    /// Parallel → `Serial` (`Affinity::Worker(1)`) → a Serial sink, no Detached
    /// step; only worker 1 has the 10 s timer. Worker 0 is the Source's owner,
    /// so both workers are pinned: Par's work reaches them through the
    /// direct-park fallback and Ser's through `Worker(1)` unparks. The source
    /// finishes only after the sink has seen every item, so no `Finished`
    /// broadcast can deliver them instead.
    #[test]
    fn pinned_worker_wake_is_not_the_timer() {
        let seen = Arc::new(AtomicU64::new(0));
        let s2 = Arc::clone(&seen);
        run_under_watchdog("pinned unpark", 5, move || {
            let mut ser = Pass::new("Ser", StepKind::Serial);
            ser.affinity = Affinity::Worker(1);
            let _ = run(2, slow(T::Worker(1)), |b| {
                b.chain(Source::exclusive_finishing_last(1_000, 100, &s2))
                    .chain(Pass::new("Par", StepKind::Parallel))
                    .chain(ser)
                    .chain(Sink::new(StepKind::Serial, &s2))
                    .into_sink_marker();
            });
        });
        assert_eq!(seen.load(AO::Relaxed), 1_000);
    }

    /// End-to-end liveness of the push gate (the gate itself is pinned by the
    /// unit test `gate_wakes_only_when_the_dispatch_pushed`): a Detached
    /// producer that reports `Progress` 10,000 times while pushing 3 items
    /// still delivers all 3.
    struct ChattyDetached {
        calls: u32,
        pushed: u32,
        /// The sink's count: `Finished` waits for all three items, so no
        /// `Finished` broadcast can deliver them in place of the gated wakes.
        seen: Arc<AtomicU64>,
    }
    impl Step for ChattyDetached {
        type Input = ();
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Chatty",
                kind: StepKind::Detached,
                sticky: false,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            self.calls += 1;
            if self.calls > 10_000 && self.seen.load(AO::Relaxed) == 3 {
                return Ok(StepOutcome::Finished);
            }
            if self.calls.is_multiple_of(3_000) && self.pushed < 3 {
                ctx.outputs
                    .push(Blob(self.pushed))
                    .map_err(|_| io::Error::other("edge has room"))?;
                self.pushed += 1;
            }
            Ok(StepOutcome::Progress)
        }
    }

    /// Every pool worker has a 10 s timer and the producer finishes only after
    /// the sink has counted every item, so under the 5 s watchdog the items
    /// reach the pool only through the gated branch's wakes: a gate that hid a
    /// push wedges the run.
    #[test]
    fn gated_detached_producer_still_delivers_every_push() {
        let seen = Arc::new(AtomicU64::new(0));
        let s2 = Arc::clone(&seen);
        run_under_watchdog("gate", 5, move || {
            let _ = run(2, slow(T::AllWorkers), |b| {
                b.chain(ChattyDetached { calls: 0, pushed: 0, seen: Arc::clone(&s2) })
                    .chain(Pass::new("Par", StepKind::Parallel))
                    .chain(Sink::new(StepKind::Serial, &s2))
                    .into_sink_marker();
            });
        });
        assert_eq!(seen.load(AO::Relaxed), 3);
    }

    /// One worker. The source is Detached and finishes early (its 1 MiB edge
    /// takes all 50 items), the slow consumer and the sink are on two more
    /// drivers, so the lone worker runs only `Par` — and once the source has
    /// finished, the only thing that can wake it while it holds an item is the
    /// slow consumer's pop (the reverse edge). Worker 0 alone has the 10 s timer.
    #[test]
    fn held_item_reverse_wake_unparks_the_lone_worker() {
        let seen = Arc::new(AtomicU64::new(0));
        let s2 = Arc::clone(&seen);
        run_under_watchdog("reverse wake t1", 8, move || {
            let mut par = Pass::new("Par", StepKind::Parallel);
            par.edge = ONE_BLOB;
            let mut slow_step = Pass::new("Slow", StepKind::Detached);
            slow_step.group = DetachedGroup::Shared("slow");
            slow_step.work = Duration::from_millis(2);
            let mut sink = Sink::new(StepKind::Detached, &s2);
            sink.group = DetachedGroup::Shared("sink");
            let _ = run(1, slow(T::Worker(0)), |b| {
                b.chain(Source::detached(50, 0, "src"))
                    .chain(par)
                    .chain(slow_step)
                    .chain(sink)
                    .into_sink_marker();
            });
        });
        assert_eq!(seen.load(AO::Relaxed), 50);
    }

    /// Two unpinned pool workers, both on 10 s timers, run only `Par`; the
    /// source (finished early), the slow consumer and the sink are on drivers.
    /// A worker whose push into `Slow`'s one-item edge was rejected is holding
    /// an item, and its only wake is the reverse `unpark` from `Slow`'s pop — so
    /// it must idle on its timer park, which an `unpark` ends, rather than on
    /// the event-count, which an `unpark` cannot reach.
    #[test]
    fn held_item_on_an_unpinned_worker_idles_where_unpark_reaches_it() {
        let seen = Arc::new(AtomicU64::new(0));
        let s2 = Arc::clone(&seen);
        run_under_watchdog("unpinned holder t2", 8, move || {
            let mut par = Pass::new("Par", StepKind::Parallel);
            par.edge = ONE_BLOB;
            let mut slow_step = Pass::new("Slow", StepKind::Detached);
            slow_step.group = DetachedGroup::Shared("slow");
            slow_step.work = Duration::from_millis(2);
            let mut sink = Sink::new(StepKind::Detached, &s2);
            sink.group = DetachedGroup::Shared("sink");
            let _ = run(2, slow(T::AllWorkers), |b| {
                b.chain(Source::detached(50, 0, "src"))
                    .chain(par)
                    .chain(slow_step)
                    .chain(sink)
                    .into_sink_marker();
            });
        });
        assert_eq!(seen.load(AO::Relaxed), 50);
    }

    /// Every pool worker pinned (`Source` Exclusive on w0, `Sink` Exclusive on
    /// w1) with a Parallel step between them, at `--threads 2`, both workers on
    /// 10 s timers. Items held by `Par` on either worker are released only by
    /// the holder-directed reverse wake from `Slow`'s pops; `Par`'s input
    /// reaches a parked worker only through the direct-park fallback. Neither
    /// worker may ever wait on the event-count.
    #[test]
    fn held_item_on_a_pinned_worker_is_woken_directly_at_t2() {
        let seen = Arc::new(AtomicU64::new(0));
        let s2 = Arc::clone(&seen);
        let out = Arc::new(parking_lot::Mutex::new(None));
        let out2 = Arc::clone(&out);
        run_under_watchdog("pinned holder t2", 8, move || {
            let mut par = Pass::new("Par", StepKind::Parallel);
            par.edge = ONE_BLOB;
            let mut slow_step = Pass::new("Slow", StepKind::Detached);
            slow_step.work = Duration::from_millis(1);
            let snap = run(2, slow(T::AllWorkers), |b| {
                b.chain(Source::exclusive(200, 0))
                    .chain(par)
                    .chain(slow_step)
                    .chain(Sink::new(StepKind::Exclusive, &s2))
                    .into_sink_marker();
            });
            *out2.lock() = Some(snap);
        });
        assert_eq!(seen.load(AO::Relaxed), 200);
        let snap = out.lock().take().unwrap();
        for w in [0usize, 1] {
            let ec = snap.worker_waits.iter().find(|x| x.0 == w).map_or(0, |x| x.1 + x.2 + x.3);
            assert_eq!(ec, 0, "pinned worker {w} must never wait on the event-count: {snap:?}");
        }
    }

    /// The finishing step (`Source`, driver "x") and its consumer (`A`, driver
    /// "x") share a driver, so `Source`'s forward target is `None`; the sink's
    /// driver "y" has the 10 s timer and is parked with an empty input until
    /// the edges close. (The unit test `on_finished_unparks_every_registered_thread`
    /// is the discriminating check that `Finished` is a broadcast.)
    #[test]
    fn finished_reaches_a_parked_driver() {
        let seen = Arc::new(AtomicU64::new(0));
        let s2 = Arc::clone(&seen);
        run_under_watchdog("finished broadcast", 5, move || {
            let mut a = Pass::new("A", StepKind::Detached);
            a.group = DetachedGroup::Shared("x");
            let mut sink = Sink::new(StepKind::Detached, &s2);
            sink.group = DetachedGroup::Shared("y");
            let _ = run(1, slow(T::Driver(1)), |b| {
                b.chain(Source::detached(0, 200, "x")).chain(a).chain(sink).into_sink_marker();
            });
        });
        assert_eq!(seen.load(AO::Relaxed), 0);
    }

    /// Cancel while the drivers are parked on 10 s timers returns promptly at
    /// one worker (no event-count exists; only the bound plan can unpark the
    /// driver) and at two. The source is on a pool worker with the normal ramp
    /// and never emits within the test.
    #[rstest]
    #[case::one_worker(1)]
    #[case::two_workers(2)]
    fn cancel_unparks_parked_drivers(#[case] threads: usize) {
        run_under_watchdog("cancel", 5, move || {
            let seen = Arc::new(AtomicU64::new(0));
            let builder = Pipeline::builder();
            builder
                .chain(Source::exclusive(1, 60_000))
                .chain(Pass::new("Det", StepKind::Detached))
                .chain(Sink::new(StepKind::Serial, &seen))
                .into_sink_marker();
            let pipeline = builder.build().expect("build");
            let cancel = pipeline.cancel_handle();
            std::thread::spawn(move || {
                std::thread::sleep(Duration::from_millis(200));
                cancel.cancel();
            });
            let t = Instant::now();
            let r = pipeline.run(PipelineConfig {
                threads,
                test_backoff: slow(T::AllDrivers),
                ..Default::default()
            });
            assert!(matches!(r, Err(crate::signal::PipelineError::Cancelled)), "{r:?}");
            assert!(
                t.elapsed() < Duration::from_secs(1),
                "cancel must unpark the driver: {:?}",
                t.elapsed()
            );
        });
    }

    /// A `Parallel` step whose work is gated by a cap: it takes its permit with
    /// [`crate::admit_input`], and holds it for `work` per item. The first item
    /// it admits waits until the cap has refused someone, then (after `settle`)
    /// sets `busy`; the last item clears it.
    #[derive(Clone)]
    struct CappedStep {
        cap: Arc<crate::PhaseCap>,
        work: Duration,
        settle: Duration,
        first: Arc<std::sync::atomic::AtomicBool>,
        busy: Arc<std::sync::atomic::AtomicBool>,
        done: Arc<AtomicU64>,
        total: u64,
    }
    impl Step for CappedStep {
        type Input = Blob;
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Capped",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let _permit = match crate::admit_input(ctx.input, Some(&self.cap)) {
                Ok(p) => p,
                Err(outcome) => return Ok(outcome),
            };
            let Some(b) = ctx.input.pop() else {
                return Ok(StepOutcome::NoProgress);
            };
            if self.first.swap(false, AO::Relaxed) && !self.settle.is_zero() {
                let deadline = Instant::now() + Duration::from_secs(5);
                while self.cap.refused() == 0 && Instant::now() < deadline {
                    std::thread::sleep(Duration::from_millis(1));
                }
                std::thread::sleep(self.settle);
                self.busy.store(true, AO::Relaxed);
            }
            if !self.work.is_zero() {
                std::thread::sleep(self.work);
            }
            assert!(ctx.outputs.push(b).is_ok(), "the output edge holds every item");
            if self.done.fetch_add(1, AO::Relaxed) + 1 == self.total {
                self.busy.store(false, AO::Relaxed);
            }
            Ok(StepOutcome::Progress)
        }
        fn phase_cap(&self) -> Option<&crate::PhaseCap> {
            Some(&self.cap)
        }
        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    /// A pass-through that, while `busy` is set and its input is empty, keeps
    /// reporting `Progress` — an uncapped step that always has work, and
    /// that the chain-order walk reaches before the capped step. It finishes
    /// only once the capped step has done all `total` items, so it cannot drop
    /// out of a worker's walk before the scenario starts.
    #[derive(Clone)]
    struct AlwaysBusy {
        busy: Arc<std::sync::atomic::AtomicBool>,
        done: Arc<AtomicU64>,
        total: u64,
    }
    impl Step for AlwaysBusy {
        type Input = Blob;
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "AlwaysBusy",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(b) = ctx.input.pop() {
                assert!(ctx.outputs.push(b).is_ok(), "the output edge holds every item");
                return Ok(StepOutcome::Progress);
            }
            if self.busy.load(AO::Relaxed) {
                std::thread::sleep(Duration::from_micros(20));
                return Ok(StepOutcome::Progress);
            }
            if ctx.input.is_drained() && self.done.load(AO::Relaxed) == self.total {
                Ok(StepOutcome::Finished)
            } else {
                Ok(StepOutcome::NoProgress)
            }
        }
        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    /// A cap release wakes the refused worker, and that worker spends the permit
    /// on the step that refused it. Two pool workers, both on 10 s timers, a
    /// cap of 1, two items. The worker that admits the first item waits until
    /// the other has been refused and parked, then makes the uncapped step —
    /// which the chain-order walk visits first — report work forever, and
    /// releases. From then on the only way the second item completes is the
    /// woken worker polling the capped step before its walk: the releasing
    /// worker spins on the busy step, the busy step ends only when the second
    /// item is done, and the woken worker's walk restarts at the busy step on
    /// every `Progress`. Without the release wake the refused worker sleeps
    /// its 10 s timer; without the resume hint it never reaches the step.
    #[test]
    fn cap_release_wakes_the_refused_worker_into_the_capped_step() {
        let cap = crate::PhaseCap::new("test-cap", 1);
        let done = Arc::new(AtomicU64::new(0));
        let seen = Arc::new(AtomicU64::new(0));
        let (cap2, done2, seen2) = (Arc::clone(&cap), Arc::clone(&done), Arc::clone(&seen));
        run_under_watchdog("cap-release resume", 5, move || {
            let busy = Arc::new(std::sync::atomic::AtomicBool::new(false));
            let always_busy =
                AlwaysBusy { busy: Arc::clone(&busy), done: Arc::clone(&done2), total: 2 };
            let capped = CappedStep {
                cap: cap2,
                work: Duration::ZERO,
                settle: Duration::from_millis(20),
                first: Arc::new(std::sync::atomic::AtomicBool::new(true)),
                busy,
                done: done2,
                total: 2,
            };
            let _ = run(2, slow(T::AllWorkers), |b| {
                b.chain(Source::detached(2, 0, "src"))
                    .chain(always_busy)
                    .chain(capped)
                    .chain(Sink::new(StepKind::Serial, &seen2))
                    .into_sink_marker();
            });
        });
        assert_eq!(done.load(AO::Relaxed), 2);
        assert_eq!(seen.load(AO::Relaxed), 2);
        assert!(cap.refused() >= 1, "the scenario needs a refusal");
    }

    /// A Detached step that takes a permit of `cap` for every item
    /// ([`crate::admit_input`]) and passes it on. Before its first try at its
    /// second item it waits (bounded) until another thread holds the cap, so
    /// that admission is refused.
    struct DetachedAdmitted {
        cap: Arc<crate::PhaseCap>,
        admitted: u32,
        waited: bool,
    }
    impl Step for DetachedAdmitted {
        type Input = Blob;
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "DetachedAdmitted",
                kind: StepKind::Detached,
                sticky: false,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if self.admitted == 1 && !self.waited && !ctx.input.is_empty() {
                self.waited = true;
                let deadline = Instant::now() + Duration::from_secs(5);
                while self.cap.active() == 0 && Instant::now() < deadline {
                    std::thread::sleep(Duration::from_millis(1));
                }
            }
            let _permit = match crate::admit_input(ctx.input, Some(&self.cap)) {
                Ok(p) => p,
                Err(outcome) => return Ok(outcome),
            };
            match ctx.input.pop() {
                Some(b) => {
                    self.admitted += 1;
                    assert!(ctx.outputs.push(b).is_ok(), "the output edge holds every item");
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
        fn phase_cap(&self) -> Option<&crate::PhaseCap> {
            Some(&self.cap)
        }
    }

    /// A cap release wakes a refused driver: the driver arms its parking bit
    /// for a cap-only park (`arm_driver`), and the release's
    /// `deliver_slot(n_workers + 0)` finds it. A cap of 1 shared by a Detached
    /// step (driver 0, 10 s timer) and the Parallel step after it. The driver
    /// passes item 1 on and, before item 2, waits until the pool step holds
    /// the cap, so its admission is refused and it parks. The pool step waits
    /// for that refusal, lets the driver park, and releases. Nothing else
    /// reaches the driver: the source pushed both items before, finishes only
    /// after the sink has seen them (no `Finished` broadcast), and the pool
    /// step's own wakes go to the sink. A driver that parked unarmed would be
    /// skipped by the release and sleep its 10 s timer.
    #[test]
    fn cap_release_wakes_a_refused_driver() {
        let cap = crate::PhaseCap::new("test-cap", 1);
        let done = Arc::new(AtomicU64::new(0));
        let seen = Arc::new(AtomicU64::new(0));
        let (cap2, done2, seen2) = (Arc::clone(&cap), Arc::clone(&done), Arc::clone(&seen));
        let snap = Arc::new(parking_lot::Mutex::new(None));
        let snap2 = Arc::clone(&snap);
        run_under_watchdog("cap-release driver", 5, move || {
            let capped = CappedStep {
                cap: Arc::clone(&cap2),
                work: Duration::ZERO,
                settle: Duration::from_millis(20),
                first: Arc::new(std::sync::atomic::AtomicBool::new(true)),
                busy: Arc::new(std::sync::atomic::AtomicBool::new(false)),
                done: done2,
                total: 2,
            };
            *snap2.lock() = Some(run(2, slow(T::Driver(0)), |b| {
                b.chain(Source::exclusive_finishing_last(2, 0, &seen2))
                    .chain(DetachedAdmitted { cap: cap2, admitted: 0, waited: false })
                    .chain(capped)
                    .chain(Sink::new(StepKind::Serial, &seen2))
                    .into_sink_marker();
            }));
        });
        assert_eq!(done.load(AO::Relaxed), 2);
        assert_eq!(seen.load(AO::Relaxed), 2);
        assert!(cap.refused() >= 1, "the scenario needs the driver refused");
        let snap = snap.lock().take().expect("snapshot");
        assert!(row(&snap, "DetachedAdmitted").capped_count >= 1, "the driver was refused");
        assert!(
            snap.detached_park_outcomes
                .iter()
                .any(|&(step, unparked, _)| step == 1 && unparked > 0),
            "the refused driver parked and was unparked: {:?}",
            snap.detached_park_outcomes
        );
    }

    /// A Detached step that consumes every item and pushes nothing: its
    /// `Progress` is push-gated, so it wakes no pool worker.
    struct Swallow;
    impl Step for Swallow {
        type Input = Blob;
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Swallow",
                kind: StepKind::Detached,
                sticky: false,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn detached_group(&self) -> DetachedGroup {
            DetachedGroup::Shared("swallow")
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(_) => Ok(StepOutcome::Progress),
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// A refused worker idles where the cap release reaches it — a timer park
    /// that `unpark` ends — and not on the event-count, which a release does
    /// not notify. Same two-worker, 10 s, cap-of-1 scenario as above, but with
    /// no pool notify left that could stand in for the release: the capped
    /// step feeds a Detached step that pushes nothing (push-gated, so no pool
    /// wake), and the busy step comes after it under a downstream-first walk,
    /// so its `Progress` goes to a driver. A refused worker parked on the
    /// event-count would sleep its 10 s timer.
    #[test]
    fn a_refused_worker_parks_where_the_release_reaches_it() {
        let cap = crate::PhaseCap::new("test-cap", 1);
        let done = Arc::new(AtomicU64::new(0));
        let seen = Arc::new(AtomicU64::new(0));
        let (cap2, done2, seen2) = (Arc::clone(&cap), Arc::clone(&done), Arc::clone(&seen));
        run_under_watchdog("cap-release park", 5, move || {
            let busy = Arc::new(std::sync::atomic::AtomicBool::new(false));
            let always_busy =
                AlwaysBusy { busy: Arc::clone(&busy), done: Arc::clone(&done2), total: 2 };
            let capped = CappedStep {
                cap: cap2,
                work: Duration::ZERO,
                settle: Duration::from_millis(20),
                first: Arc::new(std::sync::atomic::AtomicBool::new(true)),
                busy,
                done: done2,
                total: 2,
            };
            let builder = Pipeline::builder();
            builder
                .chain(Source::detached(2, 0, "src"))
                .chain(capped)
                .chain(Swallow)
                .chain(always_busy)
                .chain(Sink::new(StepKind::Detached, &seen2))
                .into_sink_marker();
            let pipeline = builder.build().expect("build");
            pipeline
                .run(
                    PipelineConfig {
                        threads: 2,
                        test_backoff: slow(T::AllWorkers),
                        ..Default::default()
                    }
                    .with_scheduler(Arc::new(crate::runtime::DrainFirstScheduler)),
                )
                .expect("run");
        });
        assert_eq!(done.load(AO::Relaxed), 2);
        assert!(cap.refused() >= 1, "the scenario needs a refusal");
    }

    /// A Detached step whose phase cap no pool step reports.
    struct DetachedCapped(Arc<crate::PhaseCap>);
    impl Step for DetachedCapped {
        type Input = Blob;
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "DetachedCapped",
                kind: StepKind::Detached,
                sticky: false,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(b) => {
                    assert!(ctx.outputs.push(b).is_ok(), "the output edge holds every item");
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
        fn phase_cap(&self) -> Option<&crate::PhaseCap> {
            Some(&self.0)
        }
    }

    /// A cap reported only by a Detached step is bound to the run's waker too:
    /// `Pipeline::run` binds every step's cap from the list it takes before the
    /// Detached steps are extracted (their placeholders report none). What is
    /// asserted is the binding; that a bound cap records and wakes a refused
    /// worker is `cap_release_wakes_the_refused_worker_into_the_capped_step`'s.
    #[test]
    fn a_cap_only_a_detached_step_reports_gets_the_release_waker() {
        let cap = crate::PhaseCap::new("detached-only", 1);
        let seen = Arc::new(AtomicU64::new(0));
        let (cap2, seen2) = (Arc::clone(&cap), Arc::clone(&seen));
        let _ = run(2, Vec::new(), move |b| {
            b.chain(Source::detached(3, 0, "src"))
                .chain(DetachedCapped(cap2))
                .chain(Sink::new(StepKind::Serial, &seen2))
                .into_sink_marker();
        });
        assert_eq!(seen.load(AO::Relaxed), 3);
        assert!(cap.is_tracking(), "the run bound a waker to the Detached step's cap");
    }

    /// Emits one item, then a second 300 ms later (once every pool worker is
    /// parked), then waits for `finish` before reporting `Finished`, so no
    /// `Finished` broadcast can stand in for the wake under test.
    struct TwoPhaseSource {
        sent: u32,
        first_at: Option<Instant>,
        finish: Arc<std::sync::atomic::AtomicBool>,
    }
    impl Step for TwoPhaseSource {
        type Input = ();
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "TwoPhaseSource",
                kind: StepKind::Detached,
                sticky: false,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let ready = match self.sent {
                0 => true,
                1 => self.first_at.is_some_and(|t| t.elapsed() >= Duration::from_millis(300)),
                _ => {
                    return Ok(if self.finish.load(AO::Relaxed) {
                        StepOutcome::Finished
                    } else {
                        StepOutcome::NoProgress
                    });
                }
            };
            if !ready {
                return Ok(StepOutcome::NoProgress);
            }
            assert!(ctx.outputs.push(Blob(self.sent)).is_ok(), "the edge holds both items");
            self.sent += 1;
            self.first_at.get_or_insert_with(Instant::now);
            Ok(StepOutcome::Progress)
        }
    }

    /// A Parallel step that counts the items it passes on, optionally admitted
    /// by its own phase cap.
    #[derive(Clone)]
    struct CountingPass(Arc<AtomicU64>, Option<Arc<crate::PhaseCap>>);
    impl Step for CountingPass {
        type Input = Blob;
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Uncapped",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let _permit = match crate::admit_input(ctx.input, self.1.as_deref()) {
                Ok(p) => p,
                Err(outcome) => return Ok(outcome),
            };
            match ctx.input.pop() {
                Some(b) => {
                    assert!(ctx.outputs.push(b).is_ok(), "the edge holds every item");
                    self.0.fetch_add(1, AO::Relaxed);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
        fn phase_cap(&self) -> Option<&crate::PhaseCap> {
            self.1.as_deref()
        }
        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    /// A `Pool` wake for a step that is uncapped, or capped by another cap Y
    /// with a free permit, reaches a worker parked only because cap X refused
    /// it. Two pool workers on 10 s timers; X's one permit is held by the test
    /// for the whole scenario. The first item passes the first step and stalls
    /// at the X-capped one, so both workers are refused by X and park; the
    /// second item, 300 ms later, is for the first step: with no event-count
    /// waiter and no direct-parked worker, only the cap-parked claim can wake a
    /// worker for it before its 10 s timer.
    #[rstest]
    #[case::uncapped(false)]
    #[case::capped_by_another_cap_with_room(true)]
    fn an_item_for_an_admissible_step_wakes_a_cap_parked_worker(#[case] other_cap: bool) {
        let cap = crate::PhaseCap::new("test-cap", 1);
        let passed = Arc::new(AtomicU64::new(0));
        let seen = Arc::new(AtomicU64::new(0));
        let finish = Arc::new(std::sync::atomic::AtomicBool::new(false));
        let held = cap.try_acquire().expect("the test holds the only permit");
        let run_thread = {
            let (cap, passed, seen, finish) =
                (Arc::clone(&cap), Arc::clone(&passed), Arc::clone(&seen), Arc::clone(&finish));
            std::thread::spawn(move || {
                let capped = CappedStep {
                    cap,
                    work: Duration::ZERO,
                    settle: Duration::ZERO,
                    first: Arc::new(std::sync::atomic::AtomicBool::new(false)),
                    busy: Arc::new(std::sync::atomic::AtomicBool::new(false)),
                    done: Arc::new(AtomicU64::new(0)),
                    total: 2,
                };
                let y = other_cap.then(|| crate::PhaseCap::new("cap-y", 2));
                let _ = run(2, slow(T::AllWorkers), |b| {
                    b.chain(TwoPhaseSource { sent: 0, first_at: None, finish })
                        .chain(CountingPass(passed, y))
                        .chain(capped)
                        .chain(Sink::new(StepKind::Detached, &seen))
                        .into_sink_marker();
                });
            })
        };
        let deadline = Instant::now() + Duration::from_secs(5);
        while passed.load(AO::Relaxed) < 2 && Instant::now() < deadline {
            std::thread::sleep(Duration::from_millis(5));
        }
        let reached = passed.load(AO::Relaxed);
        drop(held); // the cap's release wakes the refused workers for the capped step
        finish.store(true, AO::Relaxed);
        run_watchdog_join(run_thread, 30);
        assert_eq!(reached, 2, "the uncapped item waited for a 10 s timer: no worker was woken");
        assert_eq!(seen.load(AO::Relaxed), 2);
    }

    /// Join `h`, failing if it has not finished in `secs`.
    fn run_watchdog_join(h: std::thread::JoinHandle<()>, secs: u64) {
        let deadline = Instant::now() + Duration::from_secs(secs);
        while !h.is_finished() {
            assert!(Instant::now() < deadline, "the run did not finish");
            std::thread::sleep(Duration::from_millis(10));
        }
        h.join().expect("run thread");
    }

    /// Refused workers park instead of spinning: with a cap of 1 on a step that
    /// holds its permit ~200 µs per item, four pool workers refuse about 4.0
    /// times per item (804–808 for 200 items over 20 runs, with refused workers
    /// armed as cap-parked; 786–801 before that; 798–811 over 60
    /// runs, idle and fully loaded, before a refusal admitted by its own
    /// re-check withdrew its record): each release wakes exactly one refused
    /// worker, whose resume poll and pass both meet the releaser re-admitting
    /// first (permits are not handed off, so a free permit is never left
    /// waiting on a thread being scheduled), plus one refusal from each of the
    /// two other workers' passes. Since every item is one release, this is also
    /// the refusals-per-release mean. The bound, 4.8 per item (960), leaves ~19%
    /// headroom above the measured ~808 for an oversubscribed runner (each
    /// extra timer expiry or pass adds refusals) and sits ~4% below the nearest
    /// wrong idle path, the pool's direct-park fallback (1003–1502 over 40
    /// runs: a pool wake lands on the refused worker on top of its cap wake);
    /// the event-count (3710) and a timer ramp from the floor (1908) are far
    /// above. That margin is thin, so the cap-parked claim rule is pinned
    /// deterministically by `a_cap_parked_claim_needs_the_consumer_s_cap_to_have_room`,
    /// not by this bound.
    #[test]
    fn refused_workers_park_instead_of_spinning() {
        const ITEMS: u64 = 200;
        let cap = crate::PhaseCap::new("test-cap", 1);
        let seen = Arc::new(AtomicU64::new(0));
        let (cap2, seen2) = (Arc::clone(&cap), Arc::clone(&seen));
        run_under_watchdog("cap spin", 30, move || {
            let capped = CappedStep {
                cap: cap2,
                work: Duration::from_micros(200),
                settle: Duration::ZERO,
                first: Arc::new(std::sync::atomic::AtomicBool::new(false)),
                busy: Arc::new(std::sync::atomic::AtomicBool::new(false)),
                done: Arc::new(AtomicU64::new(0)),
                total: ITEMS,
            };
            let _ = run(4, Vec::new(), |b| {
                b.chain(Source::detached(u32::try_from(ITEMS).unwrap(), 0, "src"))
                    .chain(capped)
                    .chain(Sink::new(StepKind::Serial, &seen2))
                    .into_sink_marker();
            });
        });
        assert_eq!(seen.load(AO::Relaxed), ITEMS);
        assert!(
            cap.refused() * 5 <= 24 * ITEMS,
            "refused workers must park, not spin: {} refusals for {ITEMS} items",
            cap.refused()
        );
    }

    /// One step's row of a stats snapshot, by name.
    fn row(snap: &StatsSnapshot, name: &str) -> crate::runtime::stats::StepStatsSnapshot {
        snap.steps.iter().find(|(n, _)| *n == name).map(|(_, s)| *s).expect("step row")
    }

    /// Every pool worker and every driver on a 10 s timer.
    fn all_slow() -> Vec<TestBackoff> {
        vec![
            TestBackoff { target: T::AllWorkers, us: TEN_S },
            TestBackoff { target: T::AllDrivers, us: TEN_S },
        ]
    }

    /// Who holds an item on the one-slot count-bounded edge.
    #[derive(Clone, Copy, Debug)]
    enum Holder {
        ParallelPool,
        ExclusivePinned,
        DetachedDriver,
    }

    /// A producer refused by a full one-slot count-bounded edge is woken by the
    /// consumer's pop, whatever its role: unpinned pool workers (a
    /// `Process`-shaped step into a Detached sink), a pinned worker, or a
    /// driver. Every worker and driver has a 10 s timer and the slow step
    /// takes 200 µs per item, so holds keep happening, and only the pop's
    /// reverse wake can release them under the watchdog. Every refusal must
    /// mark its thread (`latched == refusals`), or an unpinned holder would
    /// idle on the event-count, which an unpark cannot reach.
    #[rstest]
    #[case::parallel_pool_holder(Holder::ParallelPool)]
    #[case::exclusive_pinned_holder(Holder::ExclusivePinned)]
    #[case::detached_driver_holder(Holder::DetachedDriver)]
    fn count_bounded_holder_is_woken_by_the_pop(#[case] holder: Holder) {
        const N: u32 = 200;
        let work = Duration::from_micros(200);
        let seen = Arc::new(AtomicU64::new(0));
        let counts = RefusalCounts::default();
        let out = Arc::new(parking_lot::Mutex::new(None));
        let (s2, c2, out2) = (Arc::clone(&seen), counts.clone(), Arc::clone(&out));
        run_under_watchdog("count-bounded holder", 5, move || {
            let snap = run(2, all_slow(), |b| match holder {
                Holder::ParallelPool => {
                    let sink = Sink { work, ..Sink::new(StepKind::Detached, &s2) };
                    b.chain(Source::detached(N, 0, "src"))
                        .chain(ProcessShaped { counts: c2, held: crate::held::HeldSlot::new() })
                        .chain(sink)
                        .into_sink_marker();
                }
                Holder::ExclusivePinned => {
                    let sink = Sink { work, ..Sink::new(StepKind::Detached, &s2) };
                    b.chain(Source::exclusive(N, 0).with_edge(ONE_COUNT).with_counts(&c2))
                        .chain(sink)
                        .into_sink_marker();
                }
                Holder::DetachedDriver => {
                    let mut pass = Pass::new("Pass", StepKind::Parallel);
                    pass.work = work;
                    b.chain(Source::detached(N, 0, "src").with_edge(ONE_COUNT).with_counts(&c2))
                        .chain(pass)
                        .chain(Sink::new(StepKind::Detached, &s2))
                        .into_sink_marker();
                }
            });
            *out2.lock() = Some(snap);
        });
        let snap = out.lock().take().expect("snapshot");
        assert_eq!(seen.load(AO::Relaxed), u64::from(N));
        assert!(counts.refusals() > 0, "the one-slot edge must refuse, or the test proves nothing");
        assert_eq!(counts.latched(), counts.refusals(), "every refusal marks its thread");
        let consumer = match holder {
            Holder::ParallelPool | Holder::ExclusivePinned => "Sink",
            Holder::DetachedDriver => "Pass",
        };
        assert!(row(&snap, consumer).reverse_wakes > 0, "the consumer's pops woke the holder");
    }

    /// Who holds an item on the one-slot count-bounded edge of a plain chain.
    #[derive(Clone, Copy, Debug)]
    enum PlainHolder {
        ExclusivePinned,
        ParallelPool,
    }

    /// The plain-chain twin of `count_bounded_holder_is_woken_by_the_pop`: no
    /// Detached step anywhere, the consumer a Serial sink pinned to worker 1.
    /// A producer refused by the full one-slot count-bounded edge — the
    /// Exclusive source on worker 0, or `ProcessShaped` clones on the pool —
    /// is woken by the sink's pop through the edge's holder set. Every worker
    /// has a 10 s timer and the sink takes 200 µs per item, so holds keep
    /// happening, and only the pop's reverse wake can release them under the
    /// watchdog. Every refusal must mark its thread.
    #[rstest]
    #[case::exclusive_pinned_holder(PlainHolder::ExclusivePinned)]
    #[case::parallel_pool_holder(PlainHolder::ParallelPool)]
    fn plain_count_bounded_holder_is_woken_by_the_pop(#[case] holder: PlainHolder) {
        const N: u32 = 200;
        let work = Duration::from_micros(200);
        let seen = Arc::new(AtomicU64::new(0));
        let counts = RefusalCounts::default();
        let out = Arc::new(parking_lot::Mutex::new(None));
        let (s2, c2, out2) = (Arc::clone(&seen), counts.clone(), Arc::clone(&out));
        run_under_watchdog("plain count-bounded holder", 5, move || {
            let sink =
                Sink { work, affinity: Affinity::Worker(1), ..Sink::new(StepKind::Serial, &s2) };
            let snap = match holder {
                PlainHolder::ExclusivePinned => run(2, all_slow(), |b| {
                    b.chain(Source::exclusive(N, 0).with_edge(ONE_COUNT).with_counts(&c2))
                        .chain(sink)
                        .into_sink_marker();
                }),
                PlainHolder::ParallelPool => run(3, all_slow(), |b| {
                    b.chain(Source::exclusive(N, 0))
                        .chain(ProcessShaped { counts: c2, held: crate::held::HeldSlot::new() })
                        .chain(sink)
                        .into_sink_marker();
                }),
            };
            *out2.lock() = Some(snap);
        });
        let snap = out.lock().take().expect("snapshot");
        assert_eq!(seen.load(AO::Relaxed), u64::from(N));
        assert!(counts.refusals() > 0, "the one-slot edge must refuse, or the test proves nothing");
        assert_eq!(counts.latched(), counts.refusals(), "every refusal marks its thread");
        assert!(row(&snap, "Sink").reverse_wakes > 0, "the sink's pops woke the holder");
    }

    /// Emits `rounds` inputs, one per round: input k only once the sink has
    /// seen 2k items and `settle` has passed since, so the sink is parked.
    /// `Finished` only once the sink has seen 2 × rounds, so no `Finished`
    /// broadcast can stand in for the wake under test. Its 1 MiB edge never
    /// refuses.
    struct Feeder {
        rounds: u32,
        next: u32,
        seen: Arc<AtomicU64>,
        start_delay: Duration,
        settle: Duration,
        started: Option<Instant>,
        ready_at: Option<Instant>,
    }
    impl Feeder {
        fn new(seen: &Arc<AtomicU64>) -> Self {
            Self {
                rounds: 8,
                next: 0,
                seen: Arc::clone(seen),
                start_delay: Duration::from_millis(100),
                settle: Duration::from_millis(20),
                started: None,
                ready_at: None,
            }
        }
    }

    /// The feeder's idle tick: a short sleep, then `NoProgress`.
    fn feeder_waits() -> StepOutcome {
        std::thread::sleep(Duration::from_millis(1));
        StepOutcome::NoProgress
    }
    impl Step for Feeder {
        type Input = ();
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "Feeder",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let started = *self.started.get_or_insert_with(Instant::now);
            if started.elapsed() < self.start_delay {
                return Ok(feeder_waits());
            }
            let seen = self.seen.load(AO::Relaxed);
            if self.next == self.rounds {
                return Ok(if seen < 2 * u64::from(self.rounds) {
                    feeder_waits()
                } else {
                    StepOutcome::Finished
                });
            }
            if seen < 2 * u64::from(self.next) {
                self.ready_at = None;
                return Ok(feeder_waits());
            }
            if self.ready_at.get_or_insert_with(Instant::now).elapsed() < self.settle {
                return Ok(feeder_waits());
            }
            self.ready_at = None;
            ctx.outputs.push(Blob(self.next)).map_err(|_| io::Error::other("the edge has room"))?;
            self.next += 1;
            Ok(StepOutcome::Progress)
        }
    }

    /// When `FanOutTwo` may retry its held y. Both gates open only once the
    /// sink has taken x.
    #[derive(Clone)]
    enum RetryGate {
        /// `d` after `FanOutTwo` first sees that the sink took x (the
        /// cross-thread tests).
        After(Duration),
        /// Once the sink has polled empty at least once since `FanOutTwo` held
        /// y (the same-driver test), so the flush lands in a pass with no other
        /// progress. The counter is the sink's `empty_polls`, snapshotted at the
        /// hold. Under the driver's reverse walk this is already implied by
        /// `took_x`; the gate pins the precondition explicitly.
        SinkPolledEmpty(Arc<AtomicU64>),
    }

    /// Shared across `FanOutTwo`'s worker copies.
    #[derive(Clone, Default)]
    struct FlushCounts {
        /// Dispatches that flushed y, found the input empty and returned
        /// `after_flush`.
        flushed_idle: Arc<AtomicU64>,
        /// Rounds whose y was pushed at once (the sink popped x in between), so
        /// they could not exercise the path.
        vacuous: Arc<AtomicU64>,
    }

    /// Process-shaped producer: one input → two outputs (x = 2k, y = 2k + 1)
    /// over a one-item edge, so y is always held. The held retry waits for
    /// `gate`; it then flushes and falls through to an empty, undrained input
    /// (the feeder is waiting), and the dispatch returns `after_flush`.
    struct FanOutTwo {
        kind: StepKind,
        affinity: Affinity,
        group: DetachedGroup,
        edge: QueueSpec,
        gate: RetryGate,
        after_flush: StepOutcome,
        seen: Arc<AtomicU64>,
        counts: FlushCounts,
        held: Option<crate::handles::Unpushed<Blob>>,
        seen_at_hold: u64,
        took_x_at: Option<Instant>,
        empties_at_hold: u64,
    }
    impl FanOutTwo {
        fn new(
            kind: StepKind,
            edge: QueueSpec,
            gate: RetryGate,
            after_flush: StepOutcome,
            seen: &Arc<AtomicU64>,
            counts: &FlushCounts,
        ) -> Self {
            Self {
                kind,
                affinity: Affinity::None,
                group: DetachedGroup::PerStep,
                edge,
                gate,
                after_flush,
                seen: Arc::clone(seen),
                counts: counts.clone(),
                held: None,
                seen_at_hold: 0,
                took_x_at: None,
                empties_at_hold: 0,
            }
        }
    }
    impl Step for FanOutTwo {
        type Input = Blob;
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "FanOutTwo",
                kind: self.kind,
                sticky: false,
                output_queues: vec![self.edge],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn detached_group(&self) -> DetachedGroup {
            self.group
        }
        fn affinity(&self) -> Affinity {
            self.affinity
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let mut flushed = false;
            if let Some(u) = self.held.take() {
                let took_x = self.seen.load(AO::Relaxed) > self.seen_at_hold;
                let open = took_x
                    && match &self.gate {
                        RetryGate::After(d) => {
                            self.took_x_at.get_or_insert_with(Instant::now).elapsed() >= *d
                        }
                        // The base is the snapshot taken at the hold.
                        RetryGate::SinkPolledEmpty(e) => e.load(AO::Relaxed) > self.empties_at_hold,
                    };
                if !open {
                    self.held = Some(u);
                    return Ok(StepOutcome::Contention); // no retry call, so no flag
                }
                match ctx.outputs.retry(u) {
                    Ok(()) => {
                        flushed = true;
                        self.took_x_at = None;
                    }
                    Err(again) => {
                        self.held = Some(again);
                        return Ok(StepOutcome::Contention);
                    }
                }
            }
            let Some(Blob(k)) = ctx.input.pop() else {
                if ctx.input.is_drained() {
                    return Ok(StepOutcome::Finished);
                }
                if flushed {
                    self.counts.flushed_idle.fetch_add(1, AO::Relaxed);
                    return Ok(self.after_flush);
                }
                return Ok(StepOutcome::NoProgress);
            };
            // Snapshot before x is pushed: a pop of x at any later point (even
            // between y's refusal and the hold) then reads as `took_x`.
            let seen_before_x = self.seen.load(AO::Relaxed);
            ctx.outputs
                .push(Blob(2 * k))
                .map_err(|_| io::Error::other("x: the edge is empty at a round's start"))?;
            match ctx.outputs.push(Blob(2 * k + 1)) {
                Err(u) => {
                    self.held = Some(u);
                    self.seen_at_hold = seen_before_x;
                    if let RetryGate::SinkPolledEmpty(e) = &self.gate {
                        self.empties_at_hold = e.load(AO::Relaxed);
                    }
                }
                Ok(()) => {
                    self.counts.vacuous.fetch_add(1, AO::Relaxed);
                }
            }
            Ok(StepOutcome::Progress)
        }
        fn new_worker_copy(&self) -> Self {
            // The configuration and the shared counters, with a fresh hold.
            Self {
                group: self.group,
                affinity: self.affinity,
                ..Self::new(
                    self.kind,
                    self.edge,
                    self.gate.clone(),
                    self.after_flush,
                    &self.seen,
                    &self.counts,
                )
            }
        }
    }

    /// The cross-thread cases of `flushed_retry_wakes_a_parked_detached_consumer`.
    #[derive(Clone, Copy, Debug)]
    enum FlushCase {
        /// A Parallel producer on `edge`, returning `after` after the flush.
        Pool(QueueSpec, StepOutcome),
        /// A Detached producer on its own driver over a gated one-blob edge.
        DetachedProducerGated,
    }

    /// A flushed held retry that falls through to an empty input wakes its
    /// parked Detached consumer, whatever the dispatch then reports. Only the
    /// sink's driver has the 10 s timer. After the sink pops x it polls empty
    /// and parks; the flush of y happens ≥ 50 ms later on another thread; the
    /// feeder waits for the sink, so nothing else reports `Progress`; a driver
    /// is never reached by the direct-park fallback and no holder record names
    /// it. So every round needs the flushed wake. `detached_producer_gated` is
    /// the behavioural check that a flushed push counts as a push for the
    /// Detached-producer gate.
    #[rstest]
    #[case::one_blob_no_progress(FlushCase::Pool(ONE_BLOB, StepOutcome::NoProgress))]
    #[case::one_blob_contention(FlushCase::Pool(ONE_BLOB, StepOutcome::Contention))]
    #[case::one_blob_capped(FlushCase::Pool(ONE_BLOB, StepOutcome::Capped))]
    #[case::one_count_no_progress(FlushCase::Pool(ONE_COUNT, StepOutcome::NoProgress))]
    #[case::one_count_contention(FlushCase::Pool(ONE_COUNT, StepOutcome::Contention))]
    #[case::one_count_capped(FlushCase::Pool(ONE_COUNT, StepOutcome::Capped))]
    #[case::detached_producer_gated(FlushCase::DetachedProducerGated)]
    fn flushed_retry_wakes_a_parked_detached_consumer(#[case] case: FlushCase) {
        let seen = Arc::new(AtomicU64::new(0));
        let counts = FlushCounts::default();
        let out = Arc::new(parking_lot::Mutex::new(None));
        let (s2, c2, out2) = (Arc::clone(&seen), counts.clone(), Arc::clone(&out));
        run_under_watchdog("flushed retry wake", 5, move || {
            let gate = RetryGate::After(Duration::from_millis(50));
            let (producer, sink_driver) = match case {
                FlushCase::Pool(edge, after) => {
                    (FanOutTwo::new(StepKind::Parallel, edge, gate, after, &s2, &c2), 0)
                }
                FlushCase::DetachedProducerGated => {
                    let p = FanOutTwo::new(
                        StepKind::Detached,
                        ONE_BLOB,
                        gate,
                        StepOutcome::NoProgress,
                        &s2,
                        &c2,
                    );
                    // Its own driver, which is driver 0 by chain order.
                    (FanOutTwo { group: DetachedGroup::Shared("p"), ..p }, 1)
                }
            };
            let snap = run(2, slow(T::Driver(sink_driver)), |b| {
                b.chain(Feeder::new(&s2))
                    .chain(producer)
                    .chain(Sink::new(StepKind::Detached, &s2))
                    .into_sink_marker();
            });
            *out2.lock() = Some(snap);
        });
        let snap = out.lock().take().expect("snapshot");
        assert_eq!(seen.load(AO::Relaxed), 16);
        let flushed_idle = counts.flushed_idle.load(AO::Relaxed);
        let vacuous = counts.vacuous.load(AO::Relaxed);
        assert!(flushed_idle >= 1, "at least one round must exercise the flushed path");
        assert_eq!(flushed_idle + vacuous, 8, "every round is flushed or vacuous");
        // At most one unpark of the sink's driver per Progress and per flushed
        // idle dispatch. Not exactly one: a delivery that finds the sink
        // running (its parking bit clear) is skipped and not counted, and
        // whether the sink has armed by then is scheduling. That every round's
        // flushed wake reached it is the watchdog's to show: the sink's only
        // other way out of its park is the 10 s timer.
        let f = row(&snap, "FanOutTwo");
        assert!(
            f.unparks_issued <= f.progress_count + flushed_idle,
            "more unparks than deliveries: {} > {} + {flushed_idle}",
            f.unparks_issued,
            f.progress_count
        );
        if matches!(case, FlushCase::DetachedProducerGated) {
            assert_eq!(f.gated_off, 0, "a flushed push counts as a push for the gate");
        }
    }

    /// The plain-chain twin of `flushed_retry_wakes_a_parked_detached_consumer`'s
    /// Pool cases: no Detached step, the sink a Serial step pinned to worker 1,
    /// which alone has the 10 s timer. The producer is a Serial step pinned to
    /// worker 0 (beside the Exclusive feeder, also on worker 0), so its held y
    /// is never on the sink's slow worker. After the sink pops x it polls empty
    /// and parks; the flush of y happens ≥ 50 ms later on worker 0 and the
    /// dispatch reports `NoProgress`; the feeder waits for the sink, so nothing
    /// else reports `Progress`. Every round therefore needs the flushed push's
    /// forward wake (a `Worker(1)` unpark).
    #[rstest]
    #[case::one_blob(ONE_BLOB)]
    #[case::one_count(ONE_COUNT)]
    fn plain_flushed_retry_wakes_its_consumer(#[case] edge: QueueSpec) {
        let seen = Arc::new(AtomicU64::new(0));
        let counts = FlushCounts::default();
        let out = Arc::new(parking_lot::Mutex::new(None));
        let (s2, c2, out2) = (Arc::clone(&seen), counts.clone(), Arc::clone(&out));
        run_under_watchdog("plain flushed retry wake", 5, move || {
            let gate = RetryGate::After(Duration::from_millis(50));
            let producer = FanOutTwo {
                affinity: Affinity::Worker(0),
                ..FanOutTwo::new(StepKind::Serial, edge, gate, StepOutcome::NoProgress, &s2, &c2)
            };
            let sink = Sink { affinity: Affinity::Worker(1), ..Sink::new(StepKind::Serial, &s2) };
            let snap = run(2, slow(T::Worker(1)), |b| {
                b.chain(Feeder::new(&s2)).chain(producer).chain(sink).into_sink_marker();
            });
            *out2.lock() = Some(snap);
        });
        let snap = out.lock().take().expect("snapshot");
        assert_eq!(seen.load(AO::Relaxed), 16);
        let flushed_idle = counts.flushed_idle.load(AO::Relaxed);
        let vacuous = counts.vacuous.load(AO::Relaxed);
        assert!(flushed_idle >= 1, "at least one round must exercise the flushed path");
        assert_eq!(flushed_idle + vacuous, 8, "every round is flushed or vacuous");
        // At most one unpark of the sink's worker per Progress and per flushed
        // idle dispatch. Not exactly one: a delivery that finds the sink
        // running (its parking bit clear) is skipped and not counted, and
        // whether the sink has armed by then is scheduling. That every round's
        // flushed wake reached it is the watchdog's to show: the sink's only
        // other way out of its park is the 10 s timer.
        let f = row(&snap, "FanOutTwo");
        assert!(
            f.unparks_issued <= f.progress_count + flushed_idle,
            "more unparks than deliveries: {} > {} + {flushed_idle}",
            f.unparks_issued,
            f.progress_count
        );
    }

    /// The reverse-walk flush-first hole: producer and consumer share one
    /// driver, which walks the sink first and has no sticky owner, so the
    /// sink's pop of x ends the pass and the producer's flush of y lands in a
    /// pass where nothing makes progress. Only the self-unpark ends the
    /// driver's 10 s timer park before the watchdog.
    #[rstest]
    #[case::one_blob(ONE_BLOB)]
    #[case::one_count(ONE_COUNT)]
    fn flushed_retry_repolls_a_same_driver_consumer(#[case] edge: QueueSpec) {
        let seen = Arc::new(AtomicU64::new(0));
        let empties = Arc::new(AtomicU64::new(0));
        let counts = FlushCounts::default();
        let out = Arc::new(parking_lot::Mutex::new(None));
        let (s2, e2, c2, out2) =
            (Arc::clone(&seen), Arc::clone(&empties), counts.clone(), Arc::clone(&out));
        run_under_watchdog("same-driver flushed retry", 5, move || {
            let producer = FanOutTwo {
                group: DetachedGroup::Shared("g"),
                ..FanOutTwo::new(
                    StepKind::Detached,
                    edge,
                    RetryGate::SinkPolledEmpty(Arc::clone(&e2)),
                    StepOutcome::NoProgress,
                    &s2,
                    &c2,
                )
            };
            let sink = Sink {
                group: DetachedGroup::Shared("g"),
                ..Sink::new(StepKind::Detached, &s2).counting_empty_polls(&e2)
            };
            let snap = run(2, slow(T::Driver(0)), |b| {
                b.chain(Feeder::new(&s2)).chain(producer).chain(sink).into_sink_marker();
            });
            *out2.lock() = Some(snap);
        });
        let snap = out.lock().take().expect("snapshot");
        assert_eq!(seen.load(AO::Relaxed), 16);
        assert_eq!(
            counts.flushed_idle.load(AO::Relaxed),
            8,
            "every round flushes into an idle pass"
        );
        assert_eq!(counts.vacuous.load(AO::Relaxed), 0, "one thread: x is never popped in between");
        assert_eq!(
            row(&snap, "FanOutTwo").unparks_issued,
            8,
            "the self-unparks; a same-thread Progress delivers nothing"
        );
    }

    /// Records each Blob's value (in arrival order) and passes it on; a fixed
    /// per-item delay makes the downstream edge fill.
    #[derive(Clone)]
    struct SlowRecordingPass {
        seen: Arc<parking_lot::Mutex<Vec<u32>>>,
        delay: Duration,
        affinity: Affinity,
        held: Option<Blob>,
    }
    impl SlowRecordingPass {
        fn serial(seen: &Arc<parking_lot::Mutex<Vec<u32>>>, delay: Duration) -> Self {
            Self { seen: Arc::clone(seen), delay, affinity: Affinity::None, held: None }
        }
    }
    impl Step for SlowRecordingPass {
        type Input = Blob;
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "SlowRecordingPass",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn affinity(&self) -> Affinity {
            self.affinity
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(b) = self.held.take() {
                return match ctx.outputs.push(b) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(u) => {
                        self.held = Some(u.into_item());
                        Ok(StepOutcome::NoProgress)
                    }
                };
            }
            match ctx.input.pop() {
                Some(b) => {
                    std::thread::sleep(self.delay);
                    self.seen.lock().push(b.0);
                    match ctx.outputs.push(b) {
                        Ok(()) => Ok(StepOutcome::Progress),
                        Err(u) => {
                            self.held = Some(u.into_item());
                            Ok(StepOutcome::Progress)
                        }
                    }
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
        fn new_worker_copy(&self) -> Self {
            Self { held: None, ..self.clone() }
        }
    }

    /// `Source(Exclusive w0) → Holder(Parallel, one-item edge) →
    /// SlowRecordingPass(Serial, w0, 1 ms) → Det(Detached) → Sink(Serial, w0)` at
    /// t2, every worker and driver on a 10 s timer. Worker 1 runs only its
    /// Holder clone, and no forward wake targets it once the source has pushed
    /// everything (every consumer downstream of the holder is on worker 0 or a
    /// driver). Its retries are push-gated, so when its clone holds an item the
    /// slow pass's pop must unpark it through the edge's holder set. The source
    /// finishes only after the sink has seen every item, so no `Finished`
    /// broadcast can release a held item instead; the sink shares worker 0 with
    /// the source, so the pass that counts the last item lets the source finish.
    /// (Worker 0 also pays the source's 1 ms wait per pass, ~4.5 ms per item, so
    /// 300 items keep the run well inside the watchdog.) Preconditions, from the
    /// run's own stats: the Holder held items (299 of 300 when measured) and the
    /// slow pass delivered reverse wakes, so the run went through the mechanism
    /// the test is named for.
    #[test]
    fn directed_held_item_resolves_by_reverse_wake() {
        let seen = Arc::new(parking_lot::Mutex::new(Vec::new()));
        let s2 = Arc::clone(&seen);
        let snap = Arc::new(parking_lot::Mutex::new(None));
        let sn2 = Arc::clone(&snap);
        let counted = Arc::new(AtomicU64::new(0));
        let c2 = Arc::clone(&counted);
        run_under_watchdog("directed held item", 5, move || {
            let mut holder = Pass::new("Holder", StepKind::Parallel);
            holder.edge = ONE_BLOB;
            let mut slow = SlowRecordingPass::serial(&s2, Duration::from_millis(1));
            slow.affinity = Affinity::Worker(0);
            let mut sink = Sink::new(StepKind::Serial, &c2);
            sink.affinity = Affinity::Worker(0);
            let stats = run(
                2,
                vec![
                    TestBackoff { target: T::AllWorkers, us: TEN_S },
                    TestBackoff { target: T::AllDrivers, us: TEN_S },
                ],
                |b| {
                    b.chain(Source::exclusive_finishing_last(300, 0, &c2))
                        .chain(holder)
                        .chain(slow)
                        .chain(Pass::new("Det", StepKind::Detached))
                        .chain(sink)
                        .into_sink_marker();
                },
            );
            *sn2.lock() = Some(stats);
        });
        let mut got = seen.lock().clone();
        got.sort_unstable();
        assert_eq!(got, (0..300).collect::<Vec<u32>>());
        assert_eq!(counted.load(AO::Relaxed), 300);
        let snap = snap.lock().take().expect("snapshot");
        let holds = row(&snap, "Holder").holds;
        assert!(holds >= 100, "fixture must hold items on the one-item edge: holds={holds}");
        let reverse = row(&snap, "SlowRecordingPass").reverse_wakes;
        assert!(reverse > 0, "the slow pass's pops must deliver reverse wakes: {reverse}");
    }

    /// A Blob's value is its ordinal: `Source` emits `Blob(0..n)` densely, so a
    /// `ByItemOrdinal` edge delivers in value order, and a refused item keeps its
    /// ordinal on retry.
    impl crate::item::Ordered for Blob {
        fn ordinal(&self) -> u64 {
            u64::from(self.0)
        }
    }

    /// Ordered Parallel pass-through: `OrderedBytesSingle<Blob>` with
    /// `ByItemOrdinal`. A refused item is held as `Unpushed` and re-sent with
    /// `retry`.
    struct OrderedPass {
        name: &'static str,
        edge: QueueSpec,
        held: Option<crate::handles::Unpushed<Blob>>,
    }
    impl OrderedPass {
        fn new(name: &'static str, edge: QueueSpec) -> Self {
            Self { name, edge, held: None }
        }
    }
    impl Step for OrderedPass {
        type Input = Blob;
        type Outputs = crate::outputs::OrderedBytesSingle<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: self.name,
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![self.edge],
                branch_ordering: vec![BranchOrdering::ByItemOrdinal],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(u) = self.held.take() {
                return match ctx.outputs.retry(u) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(u) => {
                        self.held = Some(u);
                        Ok(StepOutcome::NoProgress)
                    }
                };
            }
            match ctx.input.pop() {
                Some(b) => match ctx.outputs.push(b) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(u) => {
                        self.held = Some(u);
                        Ok(StepOutcome::Progress)
                    }
                },
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
        fn new_worker_copy(&self) -> Self {
            Self::new(self.name, self.edge)
        }
    }

    /// `OrderedPass` behind a one-item transport and a 1-byte reorder-stash cap
    /// into a slow Serial consumer, Detached tail; every thread on a 10 s timer.
    /// With one item in the transport and the next stashed, every further push
    /// is a cap rejection even when items arrive in order. Cap rejections latch
    /// the holder; in-order pops unpark it. `holds >= 100` proves the fixture
    /// reached the cap (precondition): with the cap applied ~496 of 500 pushes
    /// hold, without it at most a handful of fast-path refusals do (0-1 when
    /// measured). The oracle is the in-order values inside the watchdog.
    #[test]
    fn stash_cap_holder_is_woken_in_a_directed_chain() {
        let seen = Arc::new(parking_lot::Mutex::new(Vec::new()));
        let s2 = Arc::clone(&seen);
        let snap = Arc::new(parking_lot::Mutex::new(None));
        let sn2 = Arc::clone(&snap);
        let counted = Arc::new(AtomicU64::new(0));
        let c2 = Arc::clone(&counted);
        run_under_watchdog("stash-cap holder", 5, move || {
            let builder = Pipeline::builder();
            builder
                .chain(Source::exclusive(500, 0))
                .chain(OrderedPass::new("OrderedPass", ONE_BLOB))
                .chain(SlowRecordingPass::serial(&s2, Duration::from_micros(200)))
                .chain(Sink::new(StepKind::Detached, &c2))
                .into_sink_marker();
            let pipeline = builder.build().expect("build");
            let stats = pipeline.stats();
            pipeline
                .run(
                    PipelineConfig {
                        threads: 3,
                        test_backoff: vec![
                            TestBackoff { target: T::AllWorkers, us: TEN_S },
                            TestBackoff { target: T::AllDrivers, us: TEN_S },
                        ],
                        test_reorder_overflow_bytes: Some(1),
                        ..Default::default()
                    }
                    .with_stats(Arc::clone(&stats)),
                )
                .expect("run");
            *sn2.lock() = Some(stats.snapshot());
        });
        assert_eq!(*seen.lock(), (0..500).collect::<Vec<u32>>());
        let s = snap.lock().take().unwrap();
        let holds = row(&s, "OrderedPass").holds;
        assert!(holds >= 100, "fixture never hit the stash cap: holds={holds}");
    }

    // ───── Plain chains (no Detached step) ─────

    /// A Serial sink of configurable affinity that records each value in
    /// arrival order (chains must end in a zero-output step).
    struct RecordingSink {
        seen: Arc<parking_lot::Mutex<Vec<u32>>>,
        affinity: Affinity,
        delay: Duration,
    }
    impl Step for RecordingSink {
        type Input = Blob;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "RecordingSink",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn affinity(&self) -> Affinity {
            self.affinity
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(b) => {
                    if !self.delay.is_zero() {
                        std::thread::sleep(self.delay);
                    }
                    self.seen.lock().push(b.0);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// A Parallel pass-through that works only on wake slot `w`: other clones
    /// return `NoProgress` while input exists and `Finished` once it is
    /// drained, so the Parallel drain counter still reaches zero.
    #[derive(Clone)]
    struct OnlyOnWorker {
        name: &'static str,
        w: usize,
        edge: QueueSpec,
        held: Option<Blob>,
    }
    impl Step for OnlyOnWorker {
        type Input = Blob;
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: self.name,
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![self.edge],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if crate::runtime::wake_slot::current_slot() != Some(self.w) {
                return Ok(if ctx.input.is_drained() {
                    StepOutcome::Finished
                } else {
                    StepOutcome::NoProgress
                });
            }
            if let Some(b) = self.held.take() {
                return match ctx.outputs.push(b) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(u) => {
                        self.held = Some(u.into_item());
                        Ok(StepOutcome::NoProgress)
                    }
                };
            }
            match ctx.input.pop() {
                Some(b) => match ctx.outputs.push(b) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(u) => {
                        self.held = Some(u.into_item());
                        Ok(StepOutcome::Progress)
                    }
                },
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
        fn new_worker_copy(&self) -> Self {
            Self { held: None, ..self.clone() }
        }
    }

    /// A Serial source pinned to the reader that never produces.
    struct NeverSource;
    impl Step for NeverSource {
        type Input = ();
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "NeverSource",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn affinity(&self) -> Affinity {
            Affinity::Reader
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::NoProgress)
        }
    }

    /// `Reader → Par(ordered) → Writer` at t2 with no Detached step: both
    /// workers are pinned, so no event-count waiter exists. Forward `Pool`
    /// wakes must reach an armed direct-parked worker, and a Par clone holding
    /// on the full one-item writer edge must be unparked by the writer's pop.
    /// The edge is ordered by item (`OrderedPass`, `ByItemOrdinal`), so a held
    /// item also blocks the writer's in-order pop — the worst case — and value
    /// order is delivery order, so the order assertion is exact. A 1-byte
    /// reorder-stash cap makes the ordered edge refuse (an unbounded stash
    /// would take every out-of-order push, and nothing would hold); `holds`
    /// and the writer's `reverse_wakes` confirm it did. Every worker is on a
    /// 10 s timer.
    #[test]
    fn parallel_middle_at_t2_runs_on_unparks_not_timers() {
        let seen = Arc::new(parking_lot::Mutex::new(Vec::new()));
        let s2 = Arc::clone(&seen);
        let snap = Arc::new(parking_lot::Mutex::new(None));
        let sn2 = Arc::clone(&snap);
        run_under_watchdog("reader → par → writer at t2", 5, move || {
            let builder = Pipeline::builder();
            builder
                .chain(Source::serial_with_affinity(500, 0, Affinity::Reader))
                .chain(OrderedPass::new("Par", ONE_BLOB))
                .chain(RecordingSink {
                    seen: s2,
                    affinity: Affinity::Writer,
                    delay: Duration::ZERO,
                })
                .into_sink_marker();
            let pipeline = builder.build().expect("build");
            let stats = pipeline.stats();
            pipeline
                .run(
                    PipelineConfig {
                        threads: 2,
                        test_backoff: slow(T::AllWorkers),
                        test_reorder_overflow_bytes: Some(1),
                        ..Default::default()
                    }
                    .with_stats(Arc::clone(&stats)),
                )
                .expect("run");
            *sn2.lock() = Some(stats.snapshot());
        });
        assert_eq!(*seen.lock(), (0..500).collect::<Vec<u32>>(), "values in order, none lost");
        let snap = snap.lock().take().expect("snapshot");
        // The mechanism the test is named for ran (measured: ~475 holds and
        // ~760 reverse claims of 500 items).
        let (holds, reverse) = (row(&snap, "Par").holds, row(&snap, "RecordingSink").reverse_wakes);
        assert!(holds >= 100, "the Par clones must hold on the one-item edge: holds={holds}");
        assert!(reverse >= 100, "the writer's pops must claim the holders: reverse={reverse}");
    }

    /// At t4 only worker 0's Par clone (the reader's worker) produces, and it
    /// holds on the full one-item writer edge; only an unpark of w0 itself can
    /// make progress — a reverse wake sent to the pool would reach a clone that
    /// does nothing. One producing clone ⇒ arrival order is exact.
    #[test]
    fn held_clone_on_the_reader_worker_is_unparked_by_the_pop() {
        let seen = Arc::new(parking_lot::Mutex::new(Vec::new()));
        let s2 = Arc::clone(&seen);
        let snap = Arc::new(parking_lot::Mutex::new(None));
        let sn2 = Arc::clone(&snap);
        run_under_watchdog("held clone on w0", 5, move || {
            let par = OnlyOnWorker { name: "ParW0", w: 0, edge: ONE_BLOB, held: None };
            *sn2.lock() = Some(run(4, slow(T::AllWorkers), |b| {
                b.chain(Source::serial_with_affinity(300, 0, Affinity::Reader))
                    .chain(par)
                    .chain(RecordingSink {
                        seen: s2,
                        affinity: Affinity::Writer,
                        delay: Duration::from_millis(1),
                    })
                    .into_sink_marker();
            }));
        });
        assert_eq!(*seen.lock(), (0..300).collect::<Vec<u32>>());
        let snap = snap.lock().take().expect("snapshot");
        // The mechanism the test is named for ran (measured: ~298 holds and
        // ~299 reverse claims of 300 items).
        let (holds, reverse) =
            (row(&snap, "ParW0").holds, row(&snap, "RecordingSink").reverse_wakes);
        assert!(holds >= 100, "w0's clone must hold on the one-item edge: holds={holds}");
        assert!(reverse >= 100, "the writer's pops must claim the holder: reverse={reverse}");
    }

    /// `directed_held_item_resolves_by_reverse_wake`'s twin for a plain chain
    /// (no Detached step) at t3: a Parallel holder on a one-item edge into a
    /// slow Serial sink, every worker on a 10 s timer. As in the twin, `holds`
    /// and the sink's `reverse_wakes` confirm the run went through the
    /// mechanism; the source's one early `Finished` broadcast cannot stand in
    /// for the ~1,000 holds after it.
    #[test]
    fn plain_held_item_resolves_by_reverse_wake() {
        let seen = Arc::new(parking_lot::Mutex::new(Vec::new()));
        let s2 = Arc::clone(&seen);
        let snap = Arc::new(parking_lot::Mutex::new(None));
        let sn2 = Arc::clone(&snap);
        run_under_watchdog("plain held item", 5, move || {
            let mut holder = Pass::new("Holder", StepKind::Parallel);
            holder.edge = ONE_BLOB;
            *sn2.lock() = Some(run(3, slow(T::AllWorkers), |b| {
                b.chain(Source::exclusive(1_000, 0))
                    .chain(holder)
                    .chain(RecordingSink {
                        seen: s2,
                        affinity: Affinity::None,
                        delay: Duration::from_millis(1),
                    })
                    .into_sink_marker();
            }));
        });
        let mut got = seen.lock().clone();
        got.sort_unstable(); // one Parallel stage, unordered edge
        assert_eq!(got, (0..1_000).collect::<Vec<u32>>());
        let snap = snap.lock().take().expect("snapshot");
        // The mechanism the test is named for ran (measured: ~998 holds and
        // ~2,990 reverse claims of 1,000 items).
        let (holds, reverse) =
            (row(&snap, "Holder").holds, row(&snap, "RecordingSink").reverse_wakes);
        assert!(holds >= 100, "the Holder clones must hold on the one-item edge: holds={holds}");
        assert!(reverse >= 100, "the sink's pops must claim the holders: reverse={reverse}");
    }

    /// Cancel unparks a pinned worker parked on its 10 s timer in a plain
    /// chain, on the scheduled path. At t1 `Summary` instrumentation declines
    /// fusion, and non-fusion is observed afterwards: per-worker stats rows
    /// exist only on the scheduled path.
    #[rstest]
    #[case::t1_scheduled(1)]
    #[case::t2(2)]
    fn cancel_unparks_parked_pinned_workers_in_plain_chains(#[case] threads: usize) {
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let builder = Pipeline::builder();
            let seen = Arc::new(parking_lot::Mutex::new(Vec::new()));
            builder
                .chain(NeverSource)
                .chain(RecordingSink { seen, affinity: Affinity::Writer, delay: Duration::ZERO })
                .into_sink_marker();
            let pipeline = builder.build().expect("build");
            let stats = pipeline.stats();
            let cancel = pipeline.cancel_handle();
            std::thread::spawn(move || {
                std::thread::sleep(Duration::from_millis(50));
                cancel.cancel();
            });
            let t0 = Instant::now();
            let res = pipeline.run(
                PipelineConfig {
                    threads,
                    instrumentation: crate::builder::InstrumentationLevel::Summary,
                    test_backoff: slow(T::AllWorkers),
                    ..Default::default()
                }
                .with_stats(Arc::clone(&stats)),
            );
            let _ = tx.send((res, t0.elapsed(), stats.snapshot()));
        });
        let (res, elapsed, snap) = rx
            .recv_timeout(Duration::from_secs(5))
            .expect("WEDGED: cancel did not unpark the pinned worker");
        assert!(matches!(res, Err(crate::signal::PipelineError::Cancelled)), "{res:?}");
        assert!(elapsed < Duration::from_secs(1), "cancel took {elapsed:?}");
        assert!(!snap.workers.is_empty(), "the run must be scheduled (per-worker rows), not fused");
    }

    // ───── Event-count idle policy: one pass after a timed-out wait ─────

    /// A Parallel step that ignores its (empty) queue input and instead takes
    /// items from a test queue, on wake slot `w` only. A hook injects an item
    /// with no other worker's help, and worker `w` is the only thread that can
    /// deliver it. Non-target clones never touch the queue's contents
    /// (`is_empty` does not consume), so they cannot steal the item. `closed` is
    /// set by the hook after its last push; every clone finishes once its input
    /// is drained, `closed` is set and the queue is empty, so the Parallel drain
    /// counter reaches zero.
    #[derive(Clone)]
    struct ChannelPar {
        w: usize,
        rx: Arc<crossbeam_queue::SegQueue<Blob>>,
        closed: Arc<AtomicBool>,
        held: Option<Blob>,
    }
    impl ChannelPar {
        fn finished(&self, drained: bool) -> bool {
            drained && self.held.is_none() && self.closed.load(AO::SeqCst) && self.rx.is_empty()
        }
    }
    impl Step for ChannelPar {
        type Input = Blob;
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ChannelPar",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            let _ = ctx.input.pop(); // the source sends nothing; keep the input drained
            if wake_slot::current_slot() != Some(self.w) {
                return Ok(if self.finished(ctx.input.is_drained()) {
                    StepOutcome::Finished
                } else {
                    StepOutcome::NoProgress
                });
            }
            let item = match self.held.take() {
                Some(b) => Some(b),
                None => self.rx.pop(),
            };
            match item {
                Some(b) => match ctx.outputs.push(b) {
                    Ok(()) => Ok(StepOutcome::Progress),
                    Err(u) => {
                        self.held = Some(u.into_item());
                        Ok(StepOutcome::NoProgress)
                    }
                },
                None if self.finished(ctx.input.is_drained()) => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
        fn new_worker_copy(&self) -> Self {
            Self { held: None, ..self.clone() }
        }
    }

    /// A source that finishes at once (its edge carries nothing). A sticky
    /// Reader, so it pins worker 0 (`ParkedSleep`), leaving worker 1 as the
    /// only `EventCount` waiter.
    struct EmptySource;
    impl Step for EmptySource {
        type Input = ();
        type Outputs = Single<Blob>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "EmptySource",
                kind: StepKind::Serial,
                sticky: true,
                output_queues: vec![EDGE],
                branch_ordering: vec![BranchOrdering::None],
            }
        }
        fn affinity(&self) -> Affinity {
            Affinity::Reader
        }
        fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            Ok(StepOutcome::Finished)
        }
    }

    /// A Parallel sink that consumes on any worker (unpinned, so worker 1 stays
    /// on the `EventCount` policy) and counts items.
    #[derive(Clone)]
    struct ParallelCountingSink {
        seen: Arc<AtomicU64>,
    }
    impl Step for ParallelCountingSink {
        type Input = Blob;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ParallelCountingSink",
                kind: StepKind::Parallel,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(_) => {
                    self.seen.fetch_add(1, AO::Relaxed);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
        fn new_worker_copy(&self) -> Self {
            self.clone()
        }
    }

    /// Like [`run`], with this run's own worker-loop hooks
    /// (`PipelineConfig::test_hooks`): only this run's threads see them.
    fn run_with_hooks(
        threads: usize,
        test_backoff: Vec<TestBackoff>,
        hooks: &Arc<test_hooks::TestHooks>,
        build: impl FnOnce(&crate::builder::PipelineBuilder),
    ) -> StatsSnapshot {
        let builder = Pipeline::builder();
        build(&builder);
        let pipeline = builder.build().expect("build");
        let stats = pipeline.stats();
        pipeline
            .run(
                PipelineConfig {
                    threads,
                    test_backoff,
                    test_hooks: Some(Arc::clone(hooks)),
                    ..Default::default()
                }
                .with_stats(Arc::clone(&stats)),
            )
            .expect("run");
        stats.snapshot()
    }

    /// `EmptySource` (w0, pinned) → `ChannelPar` (w1 only) →
    /// `ParallelCountingSink`, 2 threads, with `hooks`. Worker 0 has a 10 s
    /// timer; worker 1's starts at `w1_us`.
    fn run_channel_chain(
        rx: Arc<crossbeam_queue::SegQueue<Blob>>,
        closed: Arc<AtomicBool>,
        w1_us: u64,
        seen: &Arc<AtomicU64>,
        hooks: &Arc<test_hooks::TestHooks>,
    ) {
        let seen = Arc::clone(seen);
        let backoff = vec![
            TestBackoff { target: T::Worker(0), us: TEN_S },
            TestBackoff { target: T::Worker(1), us: w1_us },
        ];
        let _ = run_with_hooks(2, backoff, hooks, move |b| {
            b.chain(EmptySource)
                .chain(ChannelPar { w: 1, rx, closed, held: None })
                .chain(ParallelCountingSink { seen })
                .into_sink_marker();
        });
    }

    /// Worker 1's idle iterations from `hooks`: `(passes, outcome)`.
    fn worker1_iterations(hooks: &test_hooks::TestHooks) -> Vec<(u64, Option<WaitOutcome>)> {
        hooks.iterations.lock().iter().filter(|(w, ..)| *w == 1).map(|&(_, p, o)| (p, o)).collect()
    }

    /// Worker 1 is the only `EventCount` waiter (worker 0 is pinned by
    /// `EmptySource`). Its override is 1 ms, so idle iterations time out
    /// deterministically (nothing can wake it). On its 5th armed wait, the
    /// hook — running after `prepare_wait` and before `wait` — pushes an item
    /// and calls `ec.notify_one()`, so that wait's outcome is deterministically
    /// `Woken`. Expected pass counts: 2 for the first idle iteration, 1 after a
    /// `TimedOut`, 2 after anything else. The hooks are this run's own, so the
    /// test holds under plain `cargo test` too.
    #[test]
    fn timed_out_wait_skips_the_pre_arm_pass() {
        let tx = Arc::new(crossbeam_queue::SegQueue::new());
        let rx = Arc::clone(&tx);
        let closed = Arc::new(AtomicBool::new(false));
        let waits = AtomicU64::new(0);
        let c2 = Arc::clone(&closed);
        let hooks = Arc::new(test_hooks::TestHooks::with_before_wait(Box::new(move |w, ec| {
            if w == 1 && waits.fetch_add(1, AO::SeqCst) == 4 {
                tx.push(Blob(1));
                c2.store(true, AO::SeqCst);
                let _ = ec.notify_one(); // after prepare_wait, before wait → that wait returns Woken
            }
        })));
        let seen = Arc::new(AtomicU64::new(0));
        run_under_watchdog("pass counts", 5, {
            let seen = Arc::clone(&seen);
            let hooks = Arc::clone(&hooks);
            move || run_channel_chain(rx, closed, 1_000, &seen, &hooks)
        });
        assert_eq!(seen.load(AO::Relaxed), 1);
        let w1 = worker1_iterations(&hooks);
        assert!(w1.len() >= 5, "worker 1 waited at least five times: {w1:?}");
        assert_eq!(w1[0].0, 2, "first idle iteration has two passes: {w1:?}");
        for pair in w1.windows(2) {
            let expected = if pair[0].1 == Some(WaitOutcome::TimedOut) { 1 } else { 2 };
            assert_eq!(pair[1].0, expected, "{w1:?}");
        }
        // EmptySource's Finished (on_finished → notify_all) can wake w1's first
        // wait, so only "some timeout before the hook" is asserted.
        assert!(
            w1[1..4].iter().any(|(_, o)| *o == Some(WaitOutcome::TimedOut)),
            "a timed-out wait precedes the hook: {w1:?}"
        );
        assert_eq!(
            w1[4].1,
            Some(WaitOutcome::Woken),
            "fifth wait sees the hook's notify (generation moved before wait): {w1:?}"
        );
    }

    /// After a `TimedOut`, the pre-arm pass is skipped but the post-arm re-poll
    /// is not. At the loop top after the first timeout, the hook pushes one item
    /// (no waiter is armed, so nothing notifies), closes the queue, and raises
    /// worker 1's deadline to 10 s. Worker 0 is on a 10 s timer and cannot
    /// deliver (`ChannelPar` runs only on worker 1). Only worker 1's post-arm
    /// re-poll can find the item: the iteration right after the hook must end
    /// in that re-poll (`cancel_wait`, outcome `None`) after one pass. A loop
    /// that skipped the re-poll waits 10 s and trips the 5 s watchdog, or, if a
    /// stray notify cut that wait short, ends the iteration with a wait outcome.
    #[test]
    fn skip_after_timeout_keeps_the_post_arm_repoll() {
        let tx = Arc::new(crossbeam_queue::SegQueue::new());
        let rx = Arc::clone(&tx);
        let closed = Arc::new(AtomicBool::new(false));
        let c2 = Arc::clone(&closed);
        let hooks = Arc::new(test_hooks::TestHooks::default());
        *hooks.at_loop_top_after_timeout.lock() = Some(Box::new(move |w, core| {
            assert_eq!(w, 1);
            tx.push(Blob(7)); // no waiter is armed, and nothing notifies
            c2.store(true, AO::SeqCst);
            core.set_override_for_test(TEN_S);
        }));
        let seen = Arc::new(AtomicU64::new(0));
        run_under_watchdog("post-arm re-poll after timeout", 5, {
            let seen = Arc::clone(&seen);
            let hooks = Arc::clone(&hooks);
            move || run_channel_chain(rx, closed, 1_000, &seen, &hooks)
        });
        assert_eq!(seen.load(AO::Relaxed), 1);
        assert!(hooks.at_loop_top_after_timeout.lock().is_none(), "the hook ran");
        let w1 = worker1_iterations(&hooks);
        let first_timeout = w1
            .iter()
            .position(|(_, o)| *o == Some(WaitOutcome::TimedOut))
            .unwrap_or_else(|| panic!("worker 1 timed out before the hook: {w1:?}"));
        assert_eq!(
            w1.get(first_timeout + 1).copied(),
            Some((1, None)),
            "the iteration after the hook found the item in its own re-poll: {w1:?}"
        );
    }

    /// A regression guard for the stash-cap reverse wake under the event-count
    /// idle policy. An ordered Parallel pass behind a one-item stash cap
    /// (`test_reorder_overflow_bytes`), into a slow recording sink; every worker
    /// on a 10 s timer, so only the reverse wake can finish. `holds > 0` is the
    /// precondition; the oracle is the in-order values inside the watchdog.
    #[test]
    fn stash_cap_hold_is_woken_by_its_reverse_wake() {
        let seen = Arc::new(parking_lot::Mutex::new(Vec::new()));
        let s2 = Arc::clone(&seen);
        let snap = Arc::new(parking_lot::Mutex::new(None));
        let sn2 = Arc::clone(&snap);
        run_under_watchdog("stash-cap hold", 5, move || {
            let builder = Pipeline::builder();
            builder
                .chain(Source::exclusive(400, 0))
                .chain(OrderedPass::new("OrderedPar", ONE_BLOB))
                .chain(RecordingSink {
                    seen: s2,
                    affinity: Affinity::None,
                    delay: Duration::from_micros(300),
                })
                .into_sink_marker();
            let pipeline = builder.build().expect("build");
            let stats = pipeline.stats();
            pipeline
                .run(
                    PipelineConfig {
                        threads: 4,
                        test_backoff: vec![TestBackoff { target: T::AllWorkers, us: TEN_S }],
                        test_reorder_overflow_bytes: Some(1),
                        ..Default::default()
                    }
                    .with_stats(Arc::clone(&stats)),
                )
                .expect("run");
            *sn2.lock() = Some(stats.snapshot());
        });
        assert_eq!(*seen.lock(), (0..400).collect::<Vec<u32>>());
        let s = snap.lock().take().unwrap();
        assert!(row(&s, "OrderedPar").holds > 0, "fixture never hit the stash cap");
    }

    /// The same chain with NO timer override (the 200µs floor in force) and an
    /// immediate sink: items held at the reorder stash cap are released by the
    /// directed holder wake. The stash-cap refusal latches its holder
    /// (`note_rejected_holder`, which also marks the thread for this test) and
    /// the in-order pop that shrinks the stash unparks it. Oracle: the run's
    /// hooks count the timer parks of exactly those latched holders that an
    /// `unpark` ended before their deadline. Without the latch no thread is
    /// ever marked, so the count is exactly 0; with it, it is positive and
    /// grows with the holds. Nothing else feeds it — not the pinned source's
    /// unparks, nor a transport-full hold, nor a holder's timer — so the oracle
    /// does not depend on how far scheduling delays a wake. `holds > 0` is the
    /// precondition; 1 000 items keep the run short under a loaded host (the
    /// 10 s watchdog is the hang guard, not a latency bound). Hold latency (p99 / mean / max) is printed, not asserted:
    /// it tracks the consumer's scheduling, so it is measured on a quiet host.
    #[test]
    fn stash_cap_holds_are_released_by_the_directed_wake() {
        let snap = Arc::new(parking_lot::Mutex::new(None));
        let sn2 = Arc::clone(&snap);
        let seen = Arc::new(parking_lot::Mutex::new(Vec::new()));
        let s2 = Arc::clone(&seen);
        let hooks = Arc::new(test_hooks::TestHooks::default());
        let h2 = Arc::clone(&hooks);
        run_under_watchdog("stash-cap holds at the floor", 10, move || {
            let builder = Pipeline::builder();
            builder
                .chain(Source::exclusive(1_000, 0))
                .chain(OrderedPass::new("OrderedPar", ONE_BLOB))
                .chain(RecordingSink { seen: s2, affinity: Affinity::None, delay: Duration::ZERO })
                .into_sink_marker();
            let pipeline = builder.build().expect("build");
            let stats = pipeline.stats();
            pipeline
                .run(
                    PipelineConfig {
                        threads: 4,
                        test_reorder_overflow_bytes: Some(1),
                        test_hooks: Some(h2),
                        ..Default::default()
                    }
                    .with_stats(Arc::clone(&stats)),
                )
                .expect("run");
            *sn2.lock() = Some(stats.snapshot());
        });
        assert_eq!(*seen.lock(), (0..1_000).collect::<Vec<u32>>());
        let s = snap.lock().take().unwrap();
        let step = row(&s, "OrderedPar");
        assert!(step.holds > 0, "fixture never hit the stash cap");
        let parks = hooks.stash_holder_parks.load(AO::Relaxed);
        let unparked = hooks.stash_holder_unparks.load(AO::Relaxed);
        eprintln!(
            "stash-cap holds={} held_retries={} p99={} ns mean={} ns max={} ns; \
             stash-cap holder parks={parks} unparked={unparked}; worker_waits={:?}",
            step.holds,
            step.held_retries,
            step.held_wait_quantile_ns(0.99),
            step.held_wait_ns / step.holds,
            step.held_wait_max_ns,
            s.worker_waits,
        );
        assert!(
            unparked > 0,
            "no stash-cap holder was unparked by the in-order pop: {parks} parks for {} holds",
            step.holds
        );
    }

    // ───── Idle policies by role ─────

    /// Runs `Src(Serial, Reader) → Par → [Tap(Detached) →] Sink(Serial,
    /// Writer)` with a 20 ms first-item delay, so every pool worker and driver
    /// goes idle before any item exists (a worker that never idles would leave
    /// fewer records than workers). 64 items keep the run short. Returns each
    /// thread's first idle deadline, from this run's own hooks.
    fn run_and_record_first_deadlines(n: usize, with_detached: bool) -> Vec<(String, Duration)> {
        let seen = Arc::new(AtomicU64::new(0));
        let hooks = Arc::new(test_hooks::TestHooks::default());
        let builder = Pipeline::builder();
        let mut sink = Sink::new(StepKind::Serial, &seen);
        sink.affinity = Affinity::Writer;
        let par = builder
            .chain(Source::serial_with_affinity(64, 20, Affinity::Reader))
            .chain(Pass::new("Par", StepKind::Parallel));
        if with_detached {
            par.chain(Pass::new("Tap", StepKind::Detached)).chain(sink).into_sink_marker();
        } else {
            par.chain(sink).into_sink_marker();
        }
        builder
            .build()
            .expect("build")
            .run(PipelineConfig {
                threads: n,
                test_hooks: Some(Arc::clone(&hooks)),
                ..Default::default()
            })
            .expect("run");
        assert_eq!(seen.load(AO::Relaxed), 64);
        std::mem::take(&mut *hooks.first_deadlines.lock())
    }

    /// Behavioural oracle: the first idle deadline each thread computes, before
    /// its first episode is counted. An unpinned worker (`EventCount`) records
    /// `jittered_wait_micros(200, w, 0)` ∈ [150, 250]µs, a pinned or lone
    /// worker (`ParkedSleep`) 20µs, a driver (`Park`) 10µs. No `test_backoff`
    /// overrides are set.
    #[rstest]
    #[case::four_workers(4, false)]
    #[case::one_worker_with_a_driver(1, true)]
    fn policies_are_assigned_by_role(#[case] n: usize, #[case] with_detached: bool) {
        // n = 4: w0 (Reader) and w3 (Writer) pinned; w1, w2 not.
        // n = 1: the lone worker; the Detached step keeps the chain from fusing.
        assert!(n <= test_hooks::MAX_TEST_WORKERS);
        let firsts = run_and_record_first_deadlines(n, with_detached);
        for (name, d) in &firsts {
            match name.as_str() {
                "fgumi-worker-1" | "fgumi-worker-2" => {
                    assert!((150..=250).contains(&d.as_micros()), "{name}: {d:?}");
                }
                "fgumi-worker-0" | "fgumi-worker-3" => {
                    assert_eq!(*d, Duration::from_micros(20), "{name}");
                }
                n if n.starts_with("fgumi-driver") => {
                    assert_eq!(*d, Duration::from_micros(10), "{n}");
                }
                other => panic!("unexpected thread {other}"),
            }
        }
        let workers = firsts.iter().filter(|(n, _)| n.starts_with("fgumi-worker")).count();
        assert_eq!(workers, n, "exactly one record per pool worker: {firsts:?}");
        assert_eq!(firsts.iter().any(|(n, _)| n.starts_with("fgumi-driver")), with_detached);
    }
}
