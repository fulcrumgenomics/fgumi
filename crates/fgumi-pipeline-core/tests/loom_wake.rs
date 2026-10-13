//! Loom models of the wake protocols in `runtime::wake`.
//!
//! Every model drives the REAL protocol types through the crate's `cfg(loom)`
//! re-exports. Under `--cfg loom` their words, flags and fences are loom's,
//! and so is the thread handle a registry slot holds (`runtime::wake_slot`):
//!
//! - holders: `HolderSet` (`latch_as`, `take`) — the record → fence → re-check
//!   every transport runs (`HolderSet::latch_and_recheck`), and the consumer's
//!   `delivery_fence` before its take;
//! - direct park: `DirectParked` (`arm`, `claim_one`) and the real
//!   `PoolEventCount::notify_one`, whose fence is loom's too;
//! - registration: `ThreadSlots` (`register`, `unpark`) and the producer's
//!   `delivery_fence` — the very calls `WakePlan::register_*`, `on_progress`
//!   and `deliver` make;
//! - the parking gate: `WakePlan`'s arms (`arm_direct`, `arm_driver`,
//!   `arm_cap_parked`) and its deliveries (`on_progress`, `deliver_holder` for
//!   a holder take, `deliver_slot` for a cap release, bound as the cap's waker
//!   exactly as `Pipeline::run` binds it), raced against the real `PhaseCap`
//!   release, record, skip and forward.
//!
//! What is not driven is the queues: their `crossbeam_queue` storage is not
//! loom-instrumented, so the room predicate each producer passes and the work
//! a re-poll looks for are loom atomics written line for line with the code;
//! the deterministic hook tests in `queues.rs` and `reorder.rs` pin each
//! production predicate.
//!
//! ```text
//! Holder.        producer: refused -> HolderSet::latch_as(slot, has_room)
//!                          (record -> fence(SeqCst) -> has_room())
//!                          -> room ? retry push : hold
//!                consumer: pop -> fence(SeqCst) -> HolderSet::take -> unpark each bit
//!                rebalancer: limit.store(higher) -> fence(SeqCst) -> take
//! Ordered stage. producer: fast-path refusal -> latch_as(slot, || !next_serial_buffered)
//!                consumer: in-order pop stores the flag -> fence(SeqCst) -> take
//! Direct park.   worker:   WakePlan::arm_direct (fetch_or -> fence(SeqCst))
//!                          -> re-poll (load work) -> park
//!                producer: work.store(true) -> notify_one (fence(SeqCst)
//!                          -> load ec waiters, 0 here)
//!                          -> DirectParked::claim_one (fetch_and claims) -> unpark
//! Registration.  thread:   ThreadSlots::register (slot store -> fence(SeqCst))
//!                          -> first pass (load work) -> arm -> re-poll -> park
//!                producer: work.store(true) -> delivery_fence (on_progress)
//!                          -> ThreadSlots::unpark (load slot, unpark if set)
//! Parking.       thread:   WakePlan::arm_direct / arm_driver (parking bit
//!                          fetch_or -> fence(SeqCst)) -> re-poll -> park
//!                producer: work.store(true) -> delivery_fence (on_progress)
//!                          -> parking bit set ? unpark : skip
//!                holder take: delivery_fence -> claim -> fence(SeqCst)
//!                          (deliver_holder) -> parking bit set ? unpark : skip
//! Cap release.   worker:   refused (record -> fence -> re-check) -> pass end
//!                          -> arm_cap_parked -> pass start (claimed ? owe)
//!                          -> record_skip (record -> fence -> re-check) -> park
//!                release:  free -> fence(SeqCst) -> claim -> deliver_slot
//!                          (parking bit set ? unpark : skip); an owed wake is
//!                          forwarded (fence -> claim -> deliver_slot)
//! ```
//!
//! `unpark`/`park` are modelled by the park token in the single-wake holder
//! models: they record "the waker unparked the holder" (an `AtomicBool`), since
//! the std token makes an unpark before the park return at once, and one
//! release can only owe one wake. The ordered fast-path, direct-park and
//! registration models park and unpark real loom threads, so a lost wake —
//! including a later round's, after the holder re-latched — is a loom deadlock.
//! `wake::tests::stress_no_lost_unpark` and
//! `event_count::tests::stress_no_lost_wakeup` are std smoke tests, not fence
//! coverage.
//!
//! Negative controls are revert-checks, since every fence and record order is
//! production code with no switch: deleting the fence in `HolderSet::latch_slot`
//! (or evaluating the predicate before the record), in `delivery_fence`, in
//! `wake_slot::arm_fence`, or in `PoolEventCount::notify_one`, makes a model
//! below fail. `ThreadSlots::register`'s fence is now covered by the arm fence
//! every timer park issues after the registration (a directed wake reaches only
//! an armed thread); the spawner-side registration controls still fail, since
//! there the registration is not ordered before the arm at all. The parking
//! gate also has in-suite controls, each a deadlock: an unfenced arm (worker
//! and holder), a holder take with no fence between its claim and the bit read,
//! a cap-refused worker that parks unarmed, and a skipped cap wake that is not
//! forwarded.
//!
//! Run: `RUSTFLAGS="--cfg loom" cargo test -p fgumi-pipeline-core --test loom_wake --release`
#![cfg(loom)]

use fgumi_pipeline_core::runtime::PoolEventCount;
use fgumi_pipeline_core::{
    BranchIdx, ChainGraph, DriverIdx, HolderSet, PhaseCap, StepIdx, StepKind, WakeEdges, WakePlan,
    current_thread, delivery_fence,
};
use loom::sync::Arc;
use loom::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering, fence};
use loom::thread;

const ITEM: u64 = 8;

/// The consumer half of every holder model: after the release, the real
/// delivery fence and take; a non-empty take unparks the holder (sets its
/// token).
fn fence_and_take(holders: &HolderSet, token: &AtomicBool) {
    delivery_fence();
    if holders.take(&mut |_| {}) != 0 {
        token.store(true, Ordering::Release); // unpark the holder
    }
}

/// A byte-bounded edge with room for one item, its real holder set and the
/// producer's park token.
struct Edge {
    bytes: AtomicU64,
    limit: AtomicU64,
    holders: HolderSet,
    token: AtomicBool,
}

impl Edge {
    fn full() -> Self {
        Self::full_for(1)
    }

    /// A full edge whose holder set has `n_slots` slots.
    fn full_for(n_slots: usize) -> Self {
        Self {
            bytes: AtomicU64::new(ITEM),
            limit: AtomicU64::new(ITEM),
            holders: HolderSet::new(n_slots),
            token: AtomicBool::new(false),
        }
    }

    fn has_room(&self) -> bool {
        self.bytes.load(Ordering::Relaxed) < self.limit.load(Ordering::Relaxed)
    }

    /// One `try_push` with the latch-and-recheck loop bounded to two retries
    /// (loom needs a bounded model; one retry already covers a pop or a raise).
    fn try_push(&self) -> bool {
        for _ in 0..2 {
            if self.has_room() {
                self.bytes.fetch_add(ITEM, Ordering::Relaxed);
                return true;
            }
            if !self.holders.latch_as(0, || self.has_room()) {
                return false; // hold
            }
        }
        false
    }

    /// The consumer's pop-side `Progress`, or the rebalancer's raise.
    fn release_and_take(&self, by_pop: bool) {
        if by_pop {
            self.bytes.fetch_sub(ITEM, Ordering::Relaxed);
        } else {
            self.limit.store(4 * ITEM, Ordering::Relaxed);
        }
        fence_and_take(&self.holders, &self.token);
    }
}

fn held_item_model(by_pop: bool) {
    loom::model(move || {
        let edge = Arc::new(Edge::full());
        let producer = {
            let edge = Arc::clone(&edge);
            thread::spawn(move || edge.try_push())
        };
        let releaser = {
            let edge = Arc::clone(&edge);
            thread::spawn(move || edge.release_and_take(by_pop))
        };
        let pushed = producer.join().unwrap();
        releaser.join().unwrap();
        // Either the push landed (it saw the room), or the holder was unparked.
        // The bad state — item held, room present, no token — is unreachable.
        let woke = edge.token.load(Ordering::Acquire);
        assert!(pushed || woke, "held item with room and no wake");
    });
}

#[test]
fn held_item_is_pushed_or_its_holder_woken_after_a_pop() {
    held_item_model(true);
}

#[test]
fn held_item_is_pushed_or_its_holder_woken_after_a_limit_raise() {
    held_item_model(false);
}

/// The same real `HolderSet`, with the caller re-checking for room BEFORE the
/// latch instead of inside it. This model can fail, and does: it shows that a
/// pass of the models above is evidence, and why the predicate must be
/// evaluated after the record and the fence (which `latch_slot` owns).
#[test]
#[should_panic(expected = "held item with room and no wake")]
fn a_recheck_before_the_latch_loses_the_wake() {
    loom::model(|| {
        let edge = Arc::new(Edge::full());
        let producer = {
            let edge = Arc::clone(&edge);
            thread::spawn(move || {
                if edge.has_room() {
                    edge.bytes.fetch_add(ITEM, Ordering::Relaxed);
                    return true;
                }
                let room = edge.has_room();
                if edge.holders.latch_as(0, || room) {
                    edge.bytes.fetch_add(ITEM, Ordering::Relaxed);
                    return true;
                }
                false
            })
        };
        let releaser = {
            let edge = Arc::clone(&edge);
            thread::spawn(move || edge.release_and_take(true))
        };
        let pushed = producer.join().unwrap();
        releaser.join().unwrap();
        assert!(pushed || edge.token.load(Ordering::Acquire), "held item with room and no wake");
    });
}

/// A count-bounded edge of capacity 1: a refused producer latches with the
/// room predicate `CountBoundedQueue` passes (`!is_full()`; a loom length
/// stands in for crossbeam's `head`/`tail`, whose loads are `SeqCst`), and
/// retries the push on room. The consumer pops, fences and takes.
#[test]
fn count_bounded_held_item_is_pushed_or_its_holder_woken() {
    loom::model(|| {
        let len = Arc::new(AtomicUsize::new(1)); // full
        let holders = Arc::new(HolderSet::new(1));
        let token = Arc::new(AtomicBool::new(false));
        let producer = {
            let (len, holders) = (Arc::clone(&len), Arc::clone(&holders));
            thread::spawn(move || {
                for _ in 0..2 {
                    if len.compare_exchange(0, 1, Ordering::SeqCst, Ordering::SeqCst).is_ok() {
                        return true;
                    }
                    if !holders.latch_as(0, || len.load(Ordering::SeqCst) < 1) {
                        return false; // hold
                    }
                }
                false
            })
        };
        let consumer = {
            let (len, holders, token) =
                (Arc::clone(&len), Arc::clone(&holders), Arc::clone(&token));
            thread::spawn(move || {
                len.fetch_sub(1, Ordering::SeqCst); // the pop
                fence_and_take(&holders, &token);
            })
        };
        let pushed = producer.join().unwrap();
        consumer.join().unwrap();
        assert!(pushed || token.load(Ordering::Acquire), "held item with room and no wake");
    });
}

/// The ordered fast path: the refused producer latches with the stage term
/// (`next_serial` no longer buffered) and retries if the flag went false. The
/// consumer's in-order pop stores the flag under the stage lock, then its
/// `Progress` fences and takes, delivering each taken holder through the real
/// `WakePlan::deliver_holder` (fence, then `deliver_slot`, which unparks it
/// only if its parking bit is set).
///
/// The producer is a real loom thread, worker 1 of the plain two-worker plan:
/// it registers its own handle, then loops as the worker loop does — latch;
/// while held, arm for its timer park (`arm_direct`), re-poll (the retry
/// re-latches), park, disarm — so a wake that finds the flag still set sends
/// it back to re-latch. `two_pops`: the first pop leaves the next ordinal
/// buffered (stores `true`), so a holder it wakes re-latches and must be woken
/// again by the second pop. A lost wake in either round leaves the producer
/// parked: a loom deadlock. `second_wake = false` drops the second pop's take —
/// the negative control below. `claim_fence = false` delivers through
/// `deliver_slot` with no fence after the claim — the second negative control.
/// (`std::sync::Arc`: loom aborts on a loom `Arc` dropped while it unwinds a
/// deadlocked model.)
fn ordered_fast_path_model(two_pops: bool, second_wake: bool, claim_fence: bool) {
    loom::model(move || {
        let flag = std::sync::Arc::new(AtomicBool::new(true)); // next_serial buffered
        let holders = std::sync::Arc::new(HolderSet::new(2));
        let plan = plain_two_step_plan();
        let producer = {
            let (flag, holders, plan) = (
                std::sync::Arc::clone(&flag),
                std::sync::Arc::clone(&holders),
                std::sync::Arc::clone(&plan),
            );
            thread::spawn(move || {
                plan.register_worker(1, current_thread());
                let retry = || holders.latch_as(1, || !flag.load(Ordering::Acquire));
                while !retry() {
                    assert!(plan.arm_direct(1), "a two-thread plan arms");
                    let pushed = retry(); // the re-poll retries the held item
                    if !pushed {
                        thread::park(); // held: wait for a take to unpark us
                    }
                    plan.disarm_direct(1);
                    if pushed {
                        break;
                    }
                }
            })
        };
        let consumer = {
            let (flag, holders, plan) = (
                std::sync::Arc::clone(&flag),
                std::sync::Arc::clone(&holders),
                std::sync::Arc::clone(&plan),
            );
            thread::spawn(move || {
                let fence_take_unpark = || {
                    delivery_fence();
                    holders.take(&mut |s| {
                        let _ =
                            if claim_fence { plan.deliver_holder(s) } else { plan.deliver_slot(s) };
                    });
                };
                if two_pops {
                    flag.store(true, Ordering::Release); // pop n, with n+1 buffered
                    fence_take_unpark();
                }
                flag.store(false, Ordering::Release); // pop of the buffered next_serial
                if second_wake {
                    fence_take_unpark();
                }
            })
        };
        consumer.join().unwrap();
        producer.join().unwrap();
    });
}

#[test]
fn ordered_fast_path_hold_is_retried_or_woken() {
    ordered_fast_path_model(false, true, true);
}

#[test]
fn ordered_fast_path_hold_is_woken_by_the_next_pop() {
    ordered_fast_path_model(true, true, true);
}

/// Negative control: with the second pop's take dropped, a holder woken by the
/// first pop re-latches while `next_serial` is still buffered, parks, and
/// nothing wakes it again. The model must catch that as a deadlock, so a pass
/// of the model above is evidence that the second round's wake is checked.
#[test]
#[should_panic(expected = "deadlock")]
fn ordered_fast_path_without_the_second_wake_deadlocks() {
    ordered_fast_path_model(true, false, true);
}

/// A two-step plan `P → C`: `C` runs on driver 0 (`consumer_on_driver`),
/// or `P` runs on driver 0 and `C` on the pool (two workers, so a `Pool` wake
/// is an event-count notify with the direct-park fallback).
fn two_step_plan(consumer_on_driver: bool) -> std::sync::Arc<WakePlan> {
    two_step_plan_capped(consumer_on_driver, None)
}

/// [`two_step_plan`] with `C` capped by `cap`.
fn two_step_plan_capped(
    consumer_on_driver: bool,
    cap: Option<std::sync::Arc<PhaseCap>>,
) -> std::sync::Arc<WakePlan> {
    let mut g = ChainGraph::new();
    let p = g.register_step("P", 1);
    let c = g.register_step("C", 0);
    g.wire(p, BranchIdx(0), c);
    let d0 = Some(DriverIdx(0));
    let (kinds, driver_of, pool) = if consumer_on_driver {
        ([StepKind::Parallel, StepKind::Detached], [None, d0], None)
    } else {
        ([StepKind::Detached, StepKind::Parallel], [d0, None], Some(PoolEventCount::new(2)))
    };
    WakePlan::build(
        &g,
        &kinds,
        &[None, None],
        &driver_of,
        &[None, cap],
        WakeEdges::NONE,
        pool.map(Into::into),
        2,
    )
}

/// A `Pool` wake that finds no event-count waiter reaches a worker that armed
/// for a timer park — or that worker's re-poll sees the work. The real plan
/// throughout: the worker registers and arms (`register_worker`,
/// `arm_direct`), re-polls, and parks (a loom park, no timer) if it found
/// nothing; the producer publishes and runs `on_progress` (`notify_one`'s
/// fence, no waiter, then the claim and unpark). A lost wake is a loom
/// deadlock. (`std::sync::Arc`: loom aborts on a loom `Arc` dropped while it
/// unwinds a deadlocked model.)
/// `cap_parked`: the worker is parked only because a phase cap refused it, so
/// it arms the cap-parked set (`arm_cap_parked`), which a `Pool` wake claims
/// after the direct-parked set when the consumer is uncapped or its cap has a
/// free permit. `consumer_cap`: `C` is capped by a (loom) cap of two with one
/// permit held, so the claim reads the cap's occupancy under the same fence.
fn pool_wake_model(cap_parked: bool, consumer_cap: bool) {
    loom::model(move || {
        let cap = consumer_cap.then(|| PhaseCap::new("c", 2));
        let held = cap.as_ref().map(|c| c.try_acquire_as(9).expect("room"));
        let plan = two_step_plan_capped(false, cap.clone());
        let work = std::sync::Arc::new(AtomicBool::new(false));
        let worker = {
            let (plan, work) = (std::sync::Arc::clone(&plan), std::sync::Arc::clone(&work));
            thread::spawn(move || {
                plan.register_worker(1, current_thread());
                let armed = if cap_parked { plan.arm_cap_parked(1) } else { plan.arm_direct(1) };
                assert!(armed, "a plan with a pool arms");
                if !work.load(Ordering::Relaxed) {
                    thread::park(); // the re-poll found nothing
                }
                if cap_parked {
                    plan.disarm_cap_parked(1);
                } else {
                    plan.disarm_direct(1);
                }
            })
        };
        let producer = {
            let (plan, work) = (std::sync::Arc::clone(&plan), std::sync::Arc::clone(&work));
            thread::spawn(move || {
                work.store(true, Ordering::Relaxed); // publish
                plan.on_progress_ungated(StepIdx(0));
            })
        };
        producer.join().unwrap();
        worker.join().unwrap();
        drop(held);
    });
}

#[test]
fn pool_work_reaches_a_direct_parked_worker() {
    pool_wake_model(false, false);
}

#[test]
fn uncapped_pool_work_reaches_a_cap_parked_worker() {
    pool_wake_model(true, false);
}

#[test]
fn pool_work_for_a_cap_with_room_reaches_a_cap_parked_worker() {
    pool_wake_model(true, true);
}

/// A driver registering its handle, raced against a producer that publishes
/// work and runs `on_progress` for the driver's consumer — the real
/// `register_driver` (slot store, then fence) and the real `on_progress` (its
/// delivery fence, then the slot read and unpark). The driver parks (a loom
/// park, no timer) when its first pass finds nothing, so a lost wake is a loom
/// deadlock.
///
/// `self_register = true` is `Pipeline::run`: the spawned thread registers
/// itself before its first pass. `false` is a spawner registering the handle
/// after `spawn`, unordered against the thread's first pass — the version
/// that loses the wake.
fn registration_model(self_register: bool) {
    // `std::sync::Arc`, not loom's: when the negative control deadlocks, loom
    // unwinds the parked target, and dropping a loom `Arc` during that cleanup
    // aborts the test binary instead of failing the model.
    loom::model(move || {
        let plan = two_step_plan(true);
        let work = std::sync::Arc::new(AtomicBool::new(false));
        let target = {
            let (plan, work) = (std::sync::Arc::clone(&plan), std::sync::Arc::clone(&work));
            thread::spawn(move || {
                if self_register {
                    plan.register_driver(DriverIdx(0), current_thread());
                }
                // First pass found nothing: arm, re-poll, park (no timer).
                if !work.load(Ordering::Relaxed) {
                    assert!(plan.arm_driver(DriverIdx(0)));
                    if !work.load(Ordering::Relaxed) {
                        thread::park();
                    }
                    plan.disarm_driver(DriverIdx(0));
                }
            })
        };
        // The spawner registers the target from its own thread, unordered
        // against both the target's first pass and the producer.
        let spawner = (!self_register).then(|| {
            let (plan, handle) = (std::sync::Arc::clone(&plan), target.thread().clone());
            thread::spawn(move || plan.register_driver(DriverIdx(0), handle))
        });
        let producer = {
            let (plan, work) = (std::sync::Arc::clone(&plan), std::sync::Arc::clone(&work));
            thread::spawn(move || {
                work.store(true, Ordering::Relaxed); // publish
                plan.on_progress_ungated(StepIdx(0));
            })
        };
        producer.join().unwrap();
        if let Some(s) = spawner {
            s.join().unwrap();
        }
        target.join().unwrap();
    });
}

#[test]
fn a_self_registered_thread_sees_the_work_or_is_unparked() {
    registration_model(true);
}

/// Negative control: registration from the spawner after `spawn` loses the
/// wake when the thread's first pass and the producer's slot read both run
/// before it: the thread parks and nothing unparks it.
#[test]
#[should_panic(expected = "deadlock")]
fn spawner_side_registration_loses_the_wake() {
    registration_model(false);
}

/// The push gate under a racing consumer. Two clones of one producer: A pushes
/// (`item = true`, `pushed += 1`) and B does not. Each snapshots `pushed`
/// before its dispatch and, after it, notifies only if `pushed` moved — behind
/// `notify_one`'s fence (`producer_fence`) and its read of the waiter count.
/// The consumer arms (`waiters += 1`), fences (`prepare_wait`), re-polls `item`.
/// In every interleaving the consumer saw the item or a notify set `token`; a
/// clone whose stale snapshot equals the counter can only skip a wake that the
/// pushing clone delivers itself.
fn push_gate_model(producer_fence: bool) {
    loom::model(move || {
        let pushed = Arc::new(AtomicU64::new(0));
        let item = Arc::new(AtomicBool::new(false));
        let waiters = Arc::new(AtomicU64::new(0));
        let token = Arc::new(AtomicBool::new(false));
        let producer = |push: bool| {
            let (pushed, item, waiters, token) =
                (pushed.clone(), item.clone(), waiters.clone(), token.clone());
            thread::spawn(move || {
                let snap = pushed.load(Ordering::Relaxed);
                if push {
                    item.store(true, Ordering::Relaxed);
                    pushed.fetch_add(1, Ordering::Relaxed);
                }
                if producer_fence {
                    fence(Ordering::SeqCst); // notify_one's fence
                }
                if pushed.load(Ordering::Relaxed) != snap && waiters.load(Ordering::Relaxed) > 0 {
                    token.store(true, Ordering::Relaxed);
                }
            })
        };
        let a = producer(true);
        let b = producer(false);
        waiters.fetch_add(1, Ordering::Relaxed);
        fence(Ordering::SeqCst);
        let saw = item.load(Ordering::Relaxed);
        a.join().unwrap();
        b.join().unwrap();
        assert!(saw || token.load(Ordering::Relaxed), "lost wake");
    });
}

#[test]
fn loom_two_producers_one_consumer_gate() {
    push_gate_model(true);
}

/// Negative control: without the producer's fence the gate can read a waiter
/// count that predates the consumer's arm while the consumer's re-poll misses
/// the item, so the model is known to be able to fail.
#[test]
#[should_panic(expected = "lost wake")]
fn the_unfenced_push_gate_loses_the_wake() {
    push_gate_model(false);
}

/// Negative control: a take that reads the parking bit right after its claim,
/// with no fence between them, loses the wake. The first pop's take claims the
/// record the holder's re-poll wrote after it armed, reads the bit clear (no
/// fence orders the arm before that read), and skips it; that pop left the next
/// ordinal buffered, so the holder's re-check sees no room and it parks; the
/// second pop finds no record.
#[test]
#[should_panic(expected = "deadlock")]
fn ordered_fast_path_without_the_claim_fence_deadlocks() {
    ordered_fast_path_model(true, true, false);
}

/// The plan of a chain with no Detached step — the shape that used to build
/// an empty-table plan: `Par(Parallel) → Ser(Serial, pinned to worker 1)` at
/// two workers, so `Par`'s branch targets `Worker(1)`.
fn plain_two_step_plan() -> std::sync::Arc<WakePlan> {
    let mut g = ChainGraph::new();
    let p = g.register_step("Par", 1);
    let s = g.register_step("Ser", 0);
    g.wire(p, BranchIdx(0), s);
    WakePlan::build(
        &g,
        &[StepKind::Parallel, StepKind::Serial],
        &[None, Some(1)],
        &[None, None],
        &[],
        WakeEdges::NONE,
        Some(PoolEventCount::new(2).into()),
        2,
    )
}

/// A pinned worker of a plain chain (no Detached step) registering itself,
/// raced against a producer that publishes work and runs `on_progress` for
/// the worker's consumer — the registrations that take effect now that every
/// chain routes wakes (a plain chain's `register_worker` used to be ignored).
/// The real plan throughout: `register_worker` (slot store, then fence) and
/// `on_progress` (its delivery fence, then the slot read and unpark). The
/// worker parks (a loom park, no timer) when its first pass finds nothing, so
/// a lost wake is a loom deadlock. `self_register = false` registers the
/// worker from the spawner after `spawn`, unordered against its first pass —
/// the version that loses the wake.
fn plain_chain_registration_model(self_register: bool) {
    // `std::sync::Arc`: see `registration_model`.
    loom::model(move || {
        let plan = plain_two_step_plan();
        let work = std::sync::Arc::new(AtomicBool::new(false));
        let worker = {
            let (plan, work) = (std::sync::Arc::clone(&plan), std::sync::Arc::clone(&work));
            thread::spawn(move || {
                if self_register {
                    plan.register_worker(1, current_thread());
                }
                // First pass found nothing: arm, re-poll, park (no timer).
                if !work.load(Ordering::Relaxed) {
                    assert!(plan.arm_direct(1));
                    if !work.load(Ordering::Relaxed) {
                        thread::park();
                    }
                    plan.disarm_direct(1);
                }
            })
        };
        let spawner = (!self_register).then(|| {
            let (plan, handle) = (std::sync::Arc::clone(&plan), worker.thread().clone());
            thread::spawn(move || plan.register_worker(1, handle))
        });
        let producer = {
            let (plan, work) = (std::sync::Arc::clone(&plan), std::sync::Arc::clone(&work));
            thread::spawn(move || {
                work.store(true, Ordering::Relaxed); // publish
                plan.on_progress_ungated(StepIdx(0));
            })
        };
        producer.join().unwrap();
        if let Some(s) = spawner {
            s.join().unwrap();
        }
        worker.join().unwrap();
    });
}

#[test]
fn plain_chain_worker_registration_races_a_push() {
    plain_chain_registration_model(true);
}

/// Negative control: registration from the spawner after `spawn` loses the
/// wake for the plain chain too.
#[test]
#[should_panic(expected = "deadlock")]
fn plain_chain_spawner_side_registration_loses_the_wake() {
    plain_chain_registration_model(false);
}

/// Unpark only a parked thread: a directed wake reads the target's parking bit
/// and skips a thread that is running. The worker of the plain plan finds
/// nothing, arms for its timer park (`arm_direct`: the parking bit, then the
/// arm fence), re-polls, and parks (a loom park, no timer) if the re-poll found
/// nothing; the producer publishes and runs `on_progress`, whose `Worker(1)`
/// delivery fences and unparks worker 1 only if its bit is set. A lost wake is
/// a loom deadlock. `fenced = false` sets the bit with no fence
/// (`mark_parking_unfenced`): the delivery can read the bit clear while the
/// re-poll misses the work.
fn parking_model(fenced: bool) {
    // `std::sync::Arc`: see `registration_model`.
    loom::model(move || {
        let plan = plain_two_step_plan();
        let work = std::sync::Arc::new(AtomicBool::new(false));
        let worker = {
            let (plan, work) = (std::sync::Arc::clone(&plan), std::sync::Arc::clone(&work));
            thread::spawn(move || {
                plan.register_worker(1, current_thread());
                if fenced {
                    assert!(plan.arm_direct(1), "a two-thread plan arms");
                } else {
                    plan.mark_parking_unfenced(1);
                }
                if !work.load(Ordering::Relaxed) {
                    thread::park(); // the re-poll found nothing
                }
                plan.disarm_direct(1);
            })
        };
        let producer = {
            let (plan, work) = (std::sync::Arc::clone(&plan), std::sync::Arc::clone(&work));
            thread::spawn(move || {
                work.store(true, Ordering::Relaxed); // publish
                plan.on_progress_ungated(StepIdx(0));
            })
        };
        producer.join().unwrap();
        worker.join().unwrap();
    });
}

#[test]
fn a_directed_wake_reaches_an_armed_thread_or_its_re_poll_sees_the_work() {
    parking_model(true);
}

/// Negative control: an arm with no fence lets the delivery read the parking
/// bit clear (and skip the unpark) while the re-poll misses the work.
#[test]
#[should_panic(expected = "deadlock")]
fn an_unfenced_park_arm_loses_the_wake() {
    parking_model(false);
}

/// The holder protocol composed with the parking bits: a refused producer
/// latches its slot (`HolderSet::latch_as`: record, fence, re-check), and only
/// when it holds does it arm for its timer park (`arm_direct`), re-check the
/// room (the re-poll's retry) and park. The consumer pops, fences
/// (`delivery_fence`) and takes the holders through the real
/// `WakePlan::deliver_holder`, which fences again and unparks slot 1 only if
/// its parking bit is set. Either the holder retries into the room or it is unparked; a lost wake
/// is a loom deadlock. `fenced = false` sets the parking bit with no fence
/// (`mark_parking_unfenced`): the latch's own fence precedes the arm, so the
/// take can read the bit clear while the re-poll's retry misses the room.
/// (`std::sync::Arc`: see `registration_model`.)
fn gated_holder_model(fenced: bool) {
    loom::model(move || {
        let plan = plain_two_step_plan();
        let edge = std::sync::Arc::new(Edge::full_for(2));
        let holder = {
            let (plan, edge) = (std::sync::Arc::clone(&plan), std::sync::Arc::clone(&edge));
            thread::spawn(move || {
                plan.register_worker(1, current_thread());
                if edge.has_room() || edge.holders.latch_as(1, || edge.has_room()) {
                    return; // pushed or retried
                }
                if fenced {
                    assert!(plan.arm_direct(1));
                } else {
                    plan.mark_parking_unfenced(1);
                }
                if !edge.has_room() {
                    thread::park(); // the re-poll's retry was refused again
                }
                plan.disarm_direct(1);
            })
        };
        let consumer = {
            let (plan, edge) = (std::sync::Arc::clone(&plan), std::sync::Arc::clone(&edge));
            thread::spawn(move || {
                edge.bytes.fetch_sub(ITEM, Ordering::Relaxed); // the pop
                delivery_fence();
                edge.holders.take(&mut |slot| {
                    let _ = plan.deliver_holder(slot);
                });
            })
        };
        consumer.join().unwrap();
        holder.join().unwrap();
    });
}

#[test]
fn a_holder_armed_for_its_park_is_unparked_or_retries() {
    gated_holder_model(true);
}

/// Negative control: a holder that sets its parking bit with no fence after
/// its latch loses the wake.
#[test]
#[should_panic(expected = "deadlock")]
fn an_unfenced_holder_park_arm_loses_the_wake() {
    gated_holder_model(false);
}

/// Bind `cap`'s release waker to `plan` exactly as `Pipeline::run` does: a
/// recorded slot is delivered through the real, parking-gated
/// `WakePlan::deliver_slot` (a `Weak`, as in the run, so the plan and the cap
/// do not keep each other alive). The models record no anonymous slot.
fn bind_plan_waker(cap: &PhaseCap, plan: &std::sync::Arc<WakePlan>) {
    let plan = std::sync::Arc::downgrade(plan);
    cap.bind_waker_for_test(std::sync::Arc::new(move |slot: Option<usize>| {
        if let (Some(plan), Some(slot)) = (plan.upgrade(), slot) {
            let _ = plan.deliver_slot(slot);
        }
    }));
}

/// Worker 1 of the capped plan taking a permit of `cap` the way the worker loop
/// does, with the plan's real arm and the cap's real record protocol (the
/// `_as` variants stand in for the thread-local wake slot and record mask):
/// a main pass (pass start: clear its own record, owing the wake if a release
/// claimed it; poll: `try_acquire`, which records on a refusal, fences and
/// re-checks; pass end: forward an owed wake it did not spend); then, refused,
/// the timer park — arm (`arm_cap_parked`, unless `arm` is off), a re-poll
/// pass that skips the full cap through `record_skip` (record, fence,
/// re-check) or polls it, park (a loom park, no timer), disarm — and around
/// again.
fn capped_worker(
    cap: std::sync::Arc<PhaseCap>,
    plan: std::sync::Arc<WakePlan>,
    arm: bool,
) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        plan.register_worker(1, current_thread());
        let mut recorded = false; // the thread's record mask bit for `cap`
        // One pass over the capped step; whether it took (and spent) a permit.
        let pass = |recorded: &mut bool, skip_full: bool| -> bool {
            let owed = std::mem::take(recorded) && !cap.clear_own_record_as(1);
            let poll = if skip_full {
                let o = cap.record_skip_as(1);
                *recorded |= o.noted;
                !o.skip
            } else {
                true
            };
            let got = poll
                && match cap.try_acquire_as(1) {
                    Some(permit) => {
                        drop(permit);
                        true
                    }
                    None => {
                        *recorded = true;
                        false
                    }
                };
            if owed && !got {
                cap.forward_wake_for_test();
            }
            got
        };
        loop {
            if pass(&mut recorded, false) {
                return;
            }
            if arm {
                assert!(plan.arm_cap_parked(1), "a plan with a pool arms");
            }
            let got = pass(&mut recorded, true);
            if !got {
                thread::park(); // the re-poll found the cap still full
            }
            if arm {
                plan.disarm_cap_parked(1);
            }
            if got {
                return;
            }
        }
    })
}

/// A cap release reaches a refused worker through the parking gate: the
/// release (`PhaseCap::release`: free, fence, claim a record) delivers the
/// claimed slot through the plan's real `deliver_slot`, which unparks worker 1
/// only if its parking bit is set. Worker 1 is refused by the one permit a
/// holder keeps, then follows the loop's real sequence (`capped_worker`).
/// Either it takes the freed permit on a poll or the release unparks it; a
/// lost wake is a loom deadlock. `arm = false` parks with no parking bit (the
/// shape of a thread whose park skips the arm): the release claims its record
/// and skips it.
fn cap_release_gated_model(arm: bool) {
    let mut builder = loom::model::Builder::new();
    builder.preemption_bound = Some(3);
    builder.check(move || {
        let cap = PhaseCap::new("c", 1);
        let plan = two_step_plan_capped(false, Some(std::sync::Arc::clone(&cap)));
        bind_plan_waker(&cap, &plan);
        let held = cap.try_acquire_as(9).expect("an idle cap admits the holder");
        let worker = capped_worker(std::sync::Arc::clone(&cap), std::sync::Arc::clone(&plan), arm);
        drop(held);
        worker.join().unwrap();
        assert_eq!(cap.active(), 0, "every permit released");
    });
}

#[test]
fn a_cap_release_reaches_an_armed_refused_worker_or_its_poll() {
    cap_release_gated_model(true);
}

/// Negative control: a refused worker that parks without arming is skipped by
/// the gated release.
#[test]
#[should_panic(expected = "deadlock")]
fn a_cap_release_skips_a_worker_parked_unarmed() {
    cap_release_gated_model(false);
}

/// A forwarded cap wake through the parking gate. Worker 0 holds a record from
/// an earlier pass and is running (never armed); worker 1 is refused and
/// parks armed (`capped_worker`). The holder releases the one permit. If the
/// release claims worker 0's stale record, `deliver_slot(0)` skips it (its bit
/// is clear) and worker 0's next pass start finds the record claimed: it owes
/// the wake, takes no permit, and forwards it (`forward_wake`: fence, claim
/// the next record, deliver — again through the gate, to worker 1). Whichever
/// record the release claims, worker 1 gets the permit; a lost wake is a loom
/// deadlock. `forward = false` drops worker 0's forward.
fn forwarded_cap_wake_gated_model(forward: bool) {
    let mut builder = loom::model::Builder::new();
    builder.preemption_bound = Some(3);
    builder.check(move || {
        let cap = PhaseCap::new("c", 1);
        let plan = two_step_plan_capped(false, Some(std::sync::Arc::clone(&cap)));
        bind_plan_waker(&cap, &plan);
        let held = cap.try_acquire_as(9).expect("an idle cap admits the holder");
        assert!(cap.try_acquire_as(0).is_none(), "worker 0's record from an earlier pass");
        let zero = {
            let (cap, plan) = (std::sync::Arc::clone(&cap), std::sync::Arc::clone(&plan));
            thread::spawn(move || {
                plan.register_worker(0, current_thread());
                // Its next pass start; that pass polls nothing on the cap.
                if !cap.clear_own_record_as(0) && forward {
                    cap.forward_wake_for_test();
                }
            })
        };
        let one = capped_worker(std::sync::Arc::clone(&cap), std::sync::Arc::clone(&plan), true);
        drop(held);
        zero.join().unwrap();
        one.join().unwrap();
        assert_eq!(cap.active(), 0, "every permit released");
    });
}

#[test]
fn a_cap_wake_skipped_for_a_running_worker_is_forwarded() {
    forwarded_cap_wake_gated_model(true);
}

/// Negative control: without worker 0's forward, a release that claims its
/// stale record (and skips it: worker 0 is not armed) leaves worker 1 parked.
#[test]
#[should_panic(expected = "deadlock")]
fn a_skipped_cap_wake_without_the_forward_is_lost() {
    forwarded_cap_wake_gated_model(false);
}
