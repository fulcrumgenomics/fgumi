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
//!   and `deliver` make.
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
//! Direct park.   worker:   DirectParked::arm (fetch_or -> fence(SeqCst))
//!                          -> re-poll (load work) -> park
//!                producer: work.store(true) -> notify_one (fence(SeqCst)
//!                          -> load ec waiters, 0 here)
//!                          -> DirectParked::claim_one (fetch_and claims) -> unpark
//! Registration.  thread:   ThreadSlots::register (slot store -> fence(SeqCst))
//!                          -> first pass (load work) -> park
//!                producer: work.store(true) -> delivery_fence (on_progress)
//!                          -> ThreadSlots::unpark (load slot, unpark if set)
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
//! `DirectParked::arm`, in `ThreadSlots::register`, or in
//! `PoolEventCount::notify_one`, makes a model below fail.
//!
//! Run: `RUSTFLAGS="--cfg loom" cargo test -p fgumi-pipeline-core --test loom_wake --release`
#![cfg(loom)]

use fgumi_pipeline_core::runtime::PoolEventCount;
use fgumi_pipeline_core::{
    BranchIdx, ChainGraph, DriverIdx, HolderSet, PhaseCap, StepIdx, StepKind, ThreadSlots,
    WakeEdges, WakePlan, current_thread, delivery_fence,
};
use loom::sync::Arc;
use loom::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
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
        Self {
            bytes: AtomicU64::new(ITEM),
            limit: AtomicU64::new(ITEM),
            holders: HolderSet::new(1),
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
/// `Progress` fences and takes, unparking each taken holder through its
/// registered handle (`ThreadSlots`).
///
/// The producer is a real loom thread: it registers its own handle, then loops
/// latch → park while held, so a wake that finds the flag still set sends it
/// back to re-latch. `two_pops`: the first pop leaves the next ordinal buffered
/// (stores `true`), so a holder it wakes re-latches and must be woken again by
/// the second pop. A lost wake in either round leaves the producer parked: a
/// loom deadlock. `second_wake = false` drops the second pop's take — the
/// negative control below. (`std::sync::Arc`: loom aborts on a loom `Arc`
/// dropped while it unwinds a deadlocked model.)
fn ordered_fast_path_model(two_pops: bool, second_wake: bool) {
    loom::model(move || {
        let flag = std::sync::Arc::new(AtomicBool::new(true)); // next_serial buffered
        let holders = std::sync::Arc::new(HolderSet::new(1));
        let slots = std::sync::Arc::new(ThreadSlots::new(1));
        let producer = {
            let (flag, holders, slots) = (
                std::sync::Arc::clone(&flag),
                std::sync::Arc::clone(&holders),
                std::sync::Arc::clone(&slots),
            );
            thread::spawn(move || {
                slots.register(0, current_thread());
                while !holders.latch_as(0, || !flag.load(Ordering::Acquire)) {
                    thread::park(); // held: wait for a take to unpark us, then re-latch
                }
            })
        };
        let consumer = {
            let (flag, holders, slots) = (
                std::sync::Arc::clone(&flag),
                std::sync::Arc::clone(&holders),
                std::sync::Arc::clone(&slots),
            );
            thread::spawn(move || {
                let fence_take_unpark = || {
                    delivery_fence();
                    holders.take(&mut |s| {
                        slots.unpark(s);
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
    ordered_fast_path_model(false, true);
}

#[test]
fn ordered_fast_path_hold_is_woken_by_the_next_pop() {
    ordered_fast_path_model(true, true);
}

/// Negative control: with the second pop's take dropped, a holder woken by the
/// first pop re-latches while `next_serial` is still buffered, parks, and
/// nothing wakes it again. The model must catch that as a deadlock, so a pass
/// of the model above is evidence that the second round's wake is checked.
#[test]
#[should_panic(expected = "deadlock")]
fn ordered_fast_path_without_the_second_wake_deadlocks() {
    ordered_fast_path_model(true, false);
}

/// A two-step Directed plan `P → C`: `C` runs on driver 0 (`consumer_on_driver`),
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
                assert!(armed, "a Directed plan with a pool arms");
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
                if !work.load(Ordering::Relaxed) {
                    thread::park(); // first pass found nothing: park, no timer
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
