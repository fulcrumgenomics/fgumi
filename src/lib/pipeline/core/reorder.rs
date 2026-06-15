//! Reorder operator on top of an `ItemQueue<Sequenced<T>>`.
//!
//! `ReorderStage<T>` is layered between a producer's transport queue and
//! the consumer's `InputHandle<T>`. The producer pushes items wrapped in
//! a framework-managed `Sequenced<T> { ordinal, item }`. The reorder
//! stage buffers items until their ordinal equals `next_serial`, then
//! releases them in order.
//!
//! Smart backpressure (deadlock avoidance):
//!   - Until `next_serial` is *observable in the reorder buffer*, the stage
//!     MUST accept everything (refusing could deadlock — the producer of
//!     `next_serial` may be one of the backpressured producers).
//!   - Once `next_serial` is in the buffer, the consumer can make progress
//!     and producers MAY be rejected by the underlying transport's normal
//!     backpressure.
//!
//! "Observable in the reorder buffer" is the conservative test: items
//! already in transport but not yet pulled into the buffer count as
//! *not* observable. This means producers over-accept when `next_serial` is
//! sitting in transport waiting to be pulled — but it can never deadlock,
//! and the over-acceptance window closes the next time a consumer calls
//! `try_pop_in_order` (which drains transport into the buffer).
//!
//! Storage layering:
//!   - "In flight" items live in the underlying `ItemQueue<Sequenced<T>>`
//!     transport (`CountBounded` / `Unbounded` — see PR 1 caveat below).
//!   - "Stashed" items (must-accept overflow when transport rejected, or
//!     items pulled but not yet at their turn) live in `state.buffer:
//!     HashMap<u64, T>`.
//!   - The transport's backpressure budget covers only in-flight items.
//!     The overflow stash has its own byte cap (`max_overflow_bytes`). The
//!     framework sizes that cap thread-awarely from the per-edge transport
//!     budget (see `apply_initial_queue_budget` / `set_max_overflow_bytes`),
//!     clamped to a fixed ceiling — so at low thread counts the stash stays
//!     small (a streaming footprint) and at high thread counts it keeps the
//!     prior lookahead headroom. `next_serial` is always exempt from the cap,
//!     so the stash bound is purely a memory/throughput knob, never a
//!     liveness constraint (any cap ≥ 0 is deadlock-free).
//!
//! `Sequenced<T>` impls `HeapSize` (see below), so `BranchOrdering::ByOrdinal`
//! / `ByItemOrdinal` compose with `QueueSpec::ByteBounded` — the canonical BAM
//! step output shape (`build_branch_ordered_bytes`).

use parking_lot::Mutex;
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

use super::item::{HeapSize, Ordered};
use super::queues::ItemQueue;

/// Default cap on a `ReorderStage`'s must-accept overflow buffer, in
/// bytes. Mirrors legacy `BACKPRESSURE_THRESHOLD_BYTES / 2` (`base.rs:704`):
/// the legacy pipeline gates non-`next_seq` reorder pushes at half the
/// 512 MB threshold (= 256 MB). We use the same value per branch so a
/// multi-stage pipeline with four ordered edges peaks around 1 GB of
/// reorder overflow under heterogeneous load — far below the unbounded
/// growth Task #29 hit (23+ GB on a 53M-record group workload). Items
/// at `next_serial` are exempt to preserve liveness.
pub const DEFAULT_REORDER_OVERFLOW_BYTES: u64 = 256 * 1024 * 1024;

/// Per-branch ordering directive in `StepProfile::branch_ordering`.
///
/// Three modes; pick based on what the consumer needs:
///
/// - **`None`** — FIFO by arrival. Cheapest. Consumer sees items in the
///   order workers happened to push them, which under multi-producer
///   concurrency is non-deterministic.
///
/// - **`ByOrdinal`** — producer-allocated ordinals via a per-branch
///   `AtomicU64` counter. The framework wraps each pushed item in a
///   `Sequenced<T>` and inserts a `ReorderStage` in front of the consumer.
///   Imposes a total order at the producer's emission point but **does
///   not preserve any pre-existing global ordering** — under multi-
///   producer Parallel concurrency, the ordinal a worker gets bears no
///   relation to the order of the input it processed. Useful for
///   single-producer steps (sources) and for cases where any deterministic
///   total order suffices.
///
/// - **`ByItemOrdinal`** — items carry their own serial via the [`Ordered`]
///   trait. The framework reads `item.ordinal()` instead of allocating
///   one. Preserves global ordering across multi-step Parallel transforms
///   when each step propagates the input's serial onto its outputs (the
///   canonical pattern for BAM pipelines: every batch carries
///   `batch_serial: u64` from its read order, and every transform
///   preserves that serial on its output items).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BranchOrdering {
    /// FIFO by arrival.
    None,
    /// Producer-allocated ordinal via a per-branch `AtomicU64` counter.
    ByOrdinal,
    /// Items carry their own ordinal via the `Ordered` trait. Requires
    /// `T: Ordered` at branch construction time.
    ByItemOrdinal,
}

/// Framework-internal wrapper carrying a producer-assigned ordinal.
/// Step authors never see this type.
pub struct Sequenced<T> {
    pub ordinal: u64,
    pub item: T,
}

impl<T> Ordered for Sequenced<T> {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

/// `Sequenced<T>` forwards `HeapSize` to the inner item, allowing byte-
/// bounded queues to wrap ordered items. The `ordinal` field is stack-
/// only and contributes nothing to the heap budget.
impl<T: HeapSize> HeapSize for Sequenced<T> {
    fn heap_size(&self) -> usize {
        self.item.heap_size()
    }
}

/// Reorder operator. Wraps an `Arc<dyn ItemQueue<Sequenced<T>>>` and
/// presents a `try_pop_in_order() -> Option<T>` surface plus an
/// ordinal-tagged `try_push(ordinal, item)`.
pub struct ReorderStage<T: Send + HeapSize + 'static> {
    transport: Arc<dyn ItemQueue<Sequenced<T>>>,
    state: Mutex<ReorderState<T>>,
    /// Cached "is `next_serial` currently in the buffer?" snapshot,
    /// updated under the state lock by writers and read **without** the
    /// lock by `try_push`'s fast path. Mirrors legacy
    /// `OrderedQueue::has_next` (`queue.rs:54`).
    ///
    /// Semantics:
    /// - `false` → `next_serial` NOT in buffer → producers MUST overflow
    ///   into the buffer if transport rejects (the slow path takes the
    ///   state lock and re-checks under it; the I1 fix from Phase 1).
    /// - `true`  → `next_serial` IS in buffer → consumer can drain →
    ///   producers CAN apply backpressure. Fast path: a lock-free
    ///   `transport.try_push` whose rejection surfaces as `Err`.
    ///
    /// Stale-`true` is harmless (the slow path is always correct).
    /// Stale-`false` is also harmless (just a missed fast-path
    /// opportunity; producer takes the slow path which still does the
    /// right thing).
    next_serial_buffered: AtomicBool,
    /// Sticky drained-observed flag for `is_drained()` callers that want
    /// to short-circuit on subsequent calls. Set lazily by `is_drained`
    /// when it observes the drained-and-empty condition.
    drained_observed: AtomicBool,
    /// Byte cap on the must-accept overflow stash: once the buffered bytes
    /// reach the cap, must-accept rejects new pushes for ordinals other than
    /// `next_serial` (the producer then holds via its `held_slot`).
    /// `next_serial` is always exempt — its producer cannot be backpressured
    /// by this cap (see Task #29) — so any cap value is deadlock-free.
    /// `u64::MAX` (`REORDER_OVERFLOW_UNBOUNDED`) means unbounded.
    ///
    /// Set once at construction (`with_max_overflow_bytes`, the no-budget
    /// fallback) and re-set once before workers spawn by the budget pass
    /// (`set_max_overflow_bytes` via `apply_initial_queue_budget`) to a
    /// thread-aware value from the per-edge transport budget. Read only on the
    /// must-accept slow path (under the `state` lock), so the atomic costs
    /// nothing on the lock-free fast path. Mirrors legacy
    /// `ReorderBufferState::can_proceed` (`base.rs:766-781`).
    max_overflow_bytes: AtomicU64,
}

/// Unbounded sentinel for [`ReorderStage::max_overflow_bytes`].
const REORDER_OVERFLOW_UNBOUNDED: u64 = u64::MAX;

/// Type-erased setter for a [`ReorderStage`]'s overflow cap, so the runtime
/// can size it without naming the branch's item type `T`. The framework
/// collects one per ordered byte-bounded branch and sets the cap in
/// `apply_initial_queue_budget` (before any worker spawns). Kept separate
/// from `BoundedQueueHandle` (the transport-resize handle) so existing
/// transport-handle impls are untouched.
pub trait ReorderCapHandle: Send + Sync {
    /// Set the overflow byte cap (`u64::MAX` = unbounded).
    fn set_max_overflow_bytes(&self, bytes: u64);
}

impl<T: Send + HeapSize + 'static> ReorderCapHandle for ReorderStage<T> {
    fn set_max_overflow_bytes(&self, bytes: u64) {
        // Set once before workers spawn; `Relaxed` is sufficient because the
        // worker pool's spawn establishes the happens-before edge, and the
        // value is only read on the must-accept slow path under `state.lock()`.
        self.max_overflow_bytes.store(bytes, Ordering::Relaxed);
    }
}

#[cfg(test)]
impl<T: Send + HeapSize + 'static> ReorderStage<T> {
    /// Read the current overflow cap (`u64::MAX` = unbounded). Test-only —
    /// lets the budget-wiring test assert `apply_initial_queue_budget` set it.
    pub(crate) fn current_max_overflow_bytes(&self) -> u64 {
        self.max_overflow_bytes.load(Ordering::Relaxed)
    }

    /// Read the current stash (overflow buffer) byte count. Test-only — lets
    /// the cap-enforcement tests observe the stash size without the production
    /// peak-tracking diagnostics.
    pub(crate) fn current_buffer_bytes(&self) -> u64 {
        self.state.lock().buffer_bytes
    }
}

struct ReorderState<T> {
    /// Items pulled from transport but not yet at their turn, indexed by
    /// ordinal. Each entry caches the item's `heap_size()` measured at
    /// insert time so we don't pay an O(item) walk on every pop or
    /// transport→buffer drain. Mirrors legacy `ReorderBuffer<T>`
    /// (`fgumi-bam-io/src/reorder.rs:50`) which stores `(T, usize)`.
    /// The size is consumed by `buffer_bytes` accounting; for items
    /// whose `heap_size()` is O(records) (e.g. position groups) this
    /// caching avoids 3× the `heap_size` cost per item flowing through.
    buffer: HashMap<u64, (T, u64)>,
    /// Tracked heap bytes of items currently in `buffer`. Updated on
    /// insert + on `try_pop_in_order` drain. Used to enforce
    /// `max_overflow_bytes`. Mirrors legacy
    /// `ReorderBufferState::heap_bytes` (`base.rs:746`).
    buffer_bytes: u64,
    /// Next ordinal we'll release.
    next_serial: u64,
}

impl<T: Send + HeapSize + 'static> ReorderStage<T> {
    #[must_use]
    pub fn new(transport: Arc<dyn ItemQueue<Sequenced<T>>>) -> Self {
        Self {
            transport,
            state: Mutex::new(ReorderState {
                buffer: HashMap::new(),
                buffer_bytes: 0,
                next_serial: 0,
            }),
            next_serial_buffered: AtomicBool::new(false),
            drained_observed: AtomicBool::new(false),
            max_overflow_bytes: AtomicU64::new(REORDER_OVERFLOW_UNBOUNDED),
        }
    }

    /// Variant that caps the must-accept overflow buffer to `max_bytes`
    /// of accumulated heap (computed via `T::heap_size()` per item).
    /// `next_serial` is exempt (always accepted to preserve liveness);
    /// other ordinals are rejected back to the producer when the buffer
    /// is at the byte cap, so producers hold via their `held_slot` and
    /// retry. Without a cap, heterogeneous-size workloads (e.g. one
    /// large position group from a busy locus while others crank
    /// through small groups) can OOM via must-accept overflow.
    ///
    /// Mirrors legacy `ReorderBufferState::can_proceed` semantics
    /// (`base.rs:766-781`): the legacy gates non-`next_seq` pushes on
    /// `heap_bytes < memory_limit / 2`. We collapse the halving into
    /// the caller-supplied cap so the hot-path check is one load + one
    /// compare.
    ///
    /// In production the framework constructs the stage with a fallback cap
    /// here and then RE-SIZES it thread-awarely via [`set_max_overflow_bytes`]
    /// in `apply_initial_queue_budget` (from the same per-edge budget as the
    /// transport queue). So this constructor's value is the no-budget fallback;
    /// the live cap tracks the transport budget.
    ///
    /// [`set_max_overflow_bytes`]: ReorderCapHandle::set_max_overflow_bytes
    #[must_use]
    pub fn with_max_overflow_bytes(
        transport: Arc<dyn ItemQueue<Sequenced<T>>>,
        max_bytes: u64,
    ) -> Self {
        Self {
            transport,
            state: Mutex::new(ReorderState {
                buffer: HashMap::new(),
                buffer_bytes: 0,
                next_serial: 0,
            }),
            next_serial_buffered: AtomicBool::new(false),
            drained_observed: AtomicBool::new(false),
            max_overflow_bytes: AtomicU64::new(max_bytes),
        }
    }

    /// Producer-side push. The framework allocates `ordinal` from a
    /// per-branch `AtomicU64` counter (see `handles.rs`).
    ///
    /// Backpressure semantics:
    ///   - If `next_serial` is not yet in the buffer, MUST accept (overflows
    ///     into the buffer if transport rejects). This prevents the deadlock
    ///     where the producer of `next_serial` is itself backpressured.
    ///   - If `next_serial` is in the buffer, the transport's normal
    ///     `try_push` rules apply; rejection surfaces as `Err((ordinal, item))`.
    ///
    /// Performance: when `next_serial` is observed in the buffer (the
    /// steady-state case once the consumer is keeping up), the push is
    /// **lock-free** — only the transport's atomic `try_push` runs. The
    /// state mutex is acquired only when overflow into the buffer might
    /// be needed.
    ///
    /// # Errors
    ///
    /// Returns `Err((ordinal, item))` when the transport rejected and we
    /// were not in must-accept mode (i.e., `next_serial` is already
    /// observable, so the consumer can drain).
    pub fn try_push(&self, ordinal: u64, item: T) -> Result<(), (u64, T)> {
        // Lock-free fast path. If next_serial is in the buffer (consumer
        // can drain), the producer is in pure-backpressure mode: a
        // transport push that succeeds returns Ok; a rejection returns
        // Err. No overflow into the buffer is possible, so we don't need
        // the state lock at all.
        //
        // Liveness: returning `Err` here (instead of falling back to the
        // slow path) is safe because `try_pop_in_order` re-derives
        // `next_serial_buffered` under the state lock, and the round-robin
        // worker driver guarantees a `try_pop_in_order` runs between
        // producer retries — so a stale-`true` flag here is corrected on
        // the next consumer poll and the producer is re-dispatched.
        if self.next_serial_buffered.load(Ordering::Acquire) {
            let seq = Sequenced { ordinal, item };
            return match self.transport.try_push(seq) {
                Ok(()) => Ok(()),
                Err(seq) => Err((seq.ordinal, seq.item)),
            };
        }

        // Slow path: must_accept may apply. Hold the state lock across
        // the must_accept check AND the transport push so a concurrent
        // consumer can't drain `next_serial` between the two — which
        // would let producers keep over-accepting into the overflow
        // buffer indefinitely (I1 from the Phase 1 review).
        //
        // Transport pushes are non-blocking, so holding the lock briefly
        // is safe (we don't risk priority inversion against blocking I/O).
        let mut state = self.state.lock();
        let must_accept = !state.buffer.contains_key(&state.next_serial);
        let seq = Sequenced { ordinal, item };

        if must_accept {
            // Apply the byte-aware overflow cap *before* attempting the
            // transport push, not only when transport is full. The consumer's
            // `try_pop_in_order` drain loop relocates the entire transport into
            // the stash while `next_serial` is absent, so a cap that only fires
            // on transport-full never binds — the consumer keeps transport
            // non-full — and the stash grows without bound (#330 zipper OOM:
            // 7.7 GB of reorder stash vs 466 MB of transport). Gating the push
            // on `buffer_bytes >= cap` regardless of transport room bounds the
            // stash to roughly `cap + one transport-worth`.
            //
            // Liveness preserved: `next_serial` is exempt (always accepted),
            // so the producer of `next_serial` can never be backpressured by
            // this cap and the consumer can always make progress. Any cap value
            // is deadlock-free (see `concurrent_tiny_cap_drains_all_in_order`,
            // which proves this with a 1-byte cap).
            let landed_next = ordinal == state.next_serial;
            let cap = self.max_overflow_bytes.load(Ordering::Relaxed);
            if !landed_next && cap != REORDER_OVERFLOW_UNBOUNDED && state.buffer_bytes >= cap {
                let Sequenced { ordinal, item } = seq;
                return Err((ordinal, item));
            }
            match self.transport.try_push(seq) {
                Ok(()) => Ok(()),
                Err(seq) => {
                    // Transport full: overflow into the stash. The cap was
                    // already checked above; `next_serial` is exempt either
                    // way, so the must-accept liveness guarantee holds.
                    let item_bytes = seq.item.heap_size() as u64;
                    state.buffer.insert(seq.ordinal, (seq.item, item_bytes));
                    state.buffer_bytes = state.buffer_bytes.saturating_add(item_bytes);
                    if landed_next {
                        self.next_serial_buffered.store(true, Ordering::Release);
                    }
                    Ok(())
                }
            }
        } else {
            // Cache was stale (false) but next_serial actually IS in the
            // buffer. Update the cache so subsequent producers take the
            // fast path, then apply backpressure.
            self.next_serial_buffered.store(true, Ordering::Release);
            match self.transport.try_push(seq) {
                Ok(()) => Ok(()),
                Err(seq) => Err((seq.ordinal, seq.item)),
            }
        }
    }

    /// Consumer-side pop. Returns `Some(T)` if `next_serial` is available,
    /// `None` if still waiting for it.
    pub fn try_pop_in_order(&self) -> Option<T> {
        let mut state = self.state.lock();
        let next = state.next_serial;

        // Drain transport into buffer until we see next_serial or transport
        // is empty. This is the only place transport → buffer movement
        // happens; it's the visibility synchronizer between producer pushes
        // and consumer reads. Each move costs one `heap_size()` call,
        // cached alongside the item so subsequent pop / cap accounting
        // doesn't pay it again.
        if !state.buffer.contains_key(&next) {
            while let Some(seq) = self.transport.try_pop() {
                let bytes = seq.item.heap_size() as u64;
                state.buffer.insert(seq.ordinal, (seq.item, bytes));
                state.buffer_bytes = state.buffer_bytes.saturating_add(bytes);
                if state.buffer.contains_key(&next) {
                    break;
                }
            }
        }

        if let Some((item, item_bytes)) = state.buffer.remove(&next) {
            state.buffer_bytes = state.buffer_bytes.saturating_sub(item_bytes);
            let new_next = next + 1;
            state.next_serial = new_next;
            // Update the cache: is the NEW next_serial in the buffer?
            // Producers reading post-update see the right state.
            let new_buffered = state.buffer.contains_key(&new_next);
            self.next_serial_buffered.store(new_buffered, Ordering::Release);
            Some(item)
        } else {
            None
        }
    }

    /// True iff transport is drained, transport is empty, and the reorder
    /// buffer is empty. Sticky: once observed true, stays true.
    pub fn is_drained(&self) -> bool {
        if self.drained_observed.load(Ordering::Acquire) {
            return true;
        }
        if !self.transport.is_drained() {
            return false;
        }
        if !self.transport.is_empty() {
            return false;
        }
        let buf_empty = self.state.lock().buffer.is_empty();
        if buf_empty {
            self.drained_observed.store(true, Ordering::Release);
        }
        buf_empty
    }

    /// Producer-side: mark transport drained. Once drained-and-empty is
    /// observed by a consumer, `is_drained()` returns true.
    pub fn mark_drained(&self) {
        self.transport.mark_drained();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::core::queues::CountBoundedQueue;

    fn make_stage(transport_capacity: usize) -> ReorderStage<u32> {
        let q: Arc<dyn ItemQueue<Sequenced<u32>>> =
            Arc::new(CountBoundedQueue::<Sequenced<u32>>::new(transport_capacity));
        ReorderStage::new(q)
    }

    #[test]
    fn pops_in_serial_order_regardless_of_push_order() {
        let s = make_stage(8);
        s.try_push(2, 200).unwrap();
        s.try_push(0, 100).unwrap();
        s.try_push(1, 150).unwrap();
        assert_eq!(s.try_pop_in_order(), Some(100));
        assert_eq!(s.try_pop_in_order(), Some(150));
        assert_eq!(s.try_pop_in_order(), Some(200));
        assert_eq!(s.try_pop_in_order(), None);
    }

    #[test]
    fn waits_for_next_serial() {
        let s = make_stage(8);
        s.try_push(1, 100).unwrap();
        s.try_push(2, 200).unwrap();
        assert_eq!(s.try_pop_in_order(), None);
        s.try_push(0, 0).unwrap();
        assert_eq!(s.try_pop_in_order(), Some(0));
        assert_eq!(s.try_pop_in_order(), Some(100));
        assert_eq!(s.try_pop_in_order(), Some(200));
    }

    #[test]
    fn must_accept_overflows_when_transport_full() {
        // Transport capacity 1. Push 4 items in arrival order; all must be
        // accepted because we're waiting for serial 0.
        let s = make_stage(1);
        s.try_push(3, 30).unwrap();
        s.try_push(2, 20).unwrap();
        s.try_push(1, 10).unwrap();
        s.try_push(0, 0).unwrap();
        // Now drain in order.
        assert_eq!(s.try_pop_in_order(), Some(0));
        assert_eq!(s.try_pop_in_order(), Some(10));
        assert_eq!(s.try_pop_in_order(), Some(20));
        assert_eq!(s.try_pop_in_order(), Some(30));
    }

    #[test]
    fn max_overflow_bytes_caps_must_accept_buffer() {
        // Use `Heavy(u32)` carrying a pretend heap size of 100 bytes per
        // item. Transport capacity 1, byte cap 200. Buffer accepts items
        // until `buffer_bytes >= 200`, then rejects (except for
        // `next_serial`, which always gets in).
        #[derive(Debug, PartialEq)]
        struct Heavy(u32);
        impl HeapSize for Heavy {
            fn heap_size(&self) -> usize {
                100
            }
        }
        let q: Arc<dyn ItemQueue<Sequenced<Heavy>>> =
            Arc::new(CountBoundedQueue::<Sequenced<Heavy>>::new(1));
        let s = ReorderStage::with_max_overflow_bytes(q, 200);

        // Fill transport (capacity 1).
        s.try_push(3, Heavy(30)).unwrap();
        // Overflow into buffer: ordinal 4 (100 B). buffer_bytes = 100,
        // < cap=200, accepted.
        s.try_push(4, Heavy(40)).unwrap();
        // Ordinal 5 lands in buffer too (buffer_bytes = 100 at the cap
        // check, before insert; insert lifts it to 200).
        s.try_push(5, Heavy(50)).unwrap();
        // Pushing ordinal 6: buffer_bytes = 200 >= cap. Reject. Producer's
        // job to retry via held_slot.
        assert_eq!(s.try_push(6, Heavy(60)), Err((6, Heavy(60))));
        // But the next_serial=0 ALWAYS gets in (liveness exemption).
        s.try_push(0, Heavy(0)).unwrap();
        // Drain to advance next_serial.
        assert_eq!(s.try_pop_in_order(), Some(Heavy(0)));
        // next_serial is now 1; ordinal 1 IS the new next_serial → must
        // always be accepted regardless of byte cap.
        s.try_push(1, Heavy(10)).unwrap();
        assert_eq!(s.try_pop_in_order(), Some(Heavy(10)));
        s.try_push(2, Heavy(20)).unwrap();
        assert_eq!(s.try_pop_in_order(), Some(Heavy(20)));
        assert_eq!(s.try_pop_in_order(), Some(Heavy(30)));
        assert_eq!(s.try_pop_in_order(), Some(Heavy(40)));
        assert_eq!(s.try_pop_in_order(), Some(Heavy(50)));
        // After buffer drains, we can finally push the rejected 6.
        s.try_push(6, Heavy(60)).unwrap();
        assert_eq!(s.try_pop_in_order(), Some(Heavy(60)));
    }

    #[test]
    fn drain_into_stash_respects_cap_on_success_path() {
        // Regression for the zipper reorder-stash OOM (issue #330). When
        // `next_serial` is withheld by a lagging producer, the consumer's
        // `try_pop_in_order` drain loop relocates the *entire* transport into
        // the stash hunting for the absent serial. If the must-accept push
        // only consults the byte cap on the transport-FULL path, the cap never
        // fires — the consumer keeps transport non-full — and the stash grows
        // without bound (measured: 7.7 GB on a 60M-read zipper run). The cap
        // must gate non-`next_serial` pushes regardless of transport room.
        #[derive(Debug, PartialEq)]
        struct Heavy(u32);
        impl HeapSize for Heavy {
            fn heap_size(&self) -> usize {
                100
            }
        }

        // Transport capacity 1 item (≤100 B in flight); stash byte cap 200 B.
        let q: Arc<dyn ItemQueue<Sequenced<Heavy>>> =
            Arc::new(CountBoundedQueue::<Sequenced<Heavy>>::new(1));
        let s = ReorderStage::with_max_overflow_bytes(q, 200);

        // Serial 0 is never pushed (the lagging worker). Push ordinals 1..=50
        // in order, polling the consumer after each push — the poll returns
        // None (serial 0 absent) but drains transport into the stash as a side
        // effect, which is the relocation that bypassed the cap. Track the peak
        // stash size after each drain.
        let mut rejected = 0u32;
        let mut peak_stash = 0u64;
        for ord in 1u32..=50 {
            if s.try_push(u64::from(ord), Heavy(ord)).is_err() {
                rejected += 1;
            }
            assert_eq!(s.try_pop_in_order(), None, "serial 0 absent → no item releases");
            peak_stash = peak_stash.max(s.current_buffer_bytes());
        }

        // The cap must fire (without the fix, the success path bypasses it
        // entirely and `rejected` stays 0).
        assert!(rejected > 0, "stash cap never fired: drain-into-stash bypassed the byte cap");
        // Peak stash is bounded by the cap (200 B) plus at most one
        // transport-worth (100 B) of overshoot.
        assert!(peak_stash <= 300, "stash grew past cap + one transport-worth: {peak_stash} B");
    }

    #[test]
    fn setter_tightens_cap_and_tiny_cap_preserves_liveness() {
        #[derive(Debug, PartialEq)]
        struct Heavy(u32);
        impl HeapSize for Heavy {
            fn heap_size(&self) -> usize {
                100
            }
        }

        // (a) The `ReorderCapHandle` setter tightens an initially-unbounded
        //     stage's cap (this is how `apply_initial_queue_budget` sizes it).
        let q: Arc<dyn ItemQueue<Sequenced<Heavy>>> =
            Arc::new(CountBoundedQueue::<Sequenced<Heavy>>::new(1));
        let s = ReorderStage::new(q);
        s.set_max_overflow_bytes(150);
        s.try_push(3, Heavy(30)).unwrap(); // -> transport (cap 1)
        s.try_push(4, Heavy(40)).unwrap(); // overflow buffer: 0 < 150
        s.try_push(5, Heavy(50)).unwrap(); // overflow buffer: 100 < 150
        // buffer_bytes is now 200 >= 150 → a non-next push is rejected,
        // proving the setter's value took effect.
        assert_eq!(s.try_push(6, Heavy(60)), Err((6, Heavy(60))), "tightened cap rejects");

        // (b) Liveness for any cap (§10.3): a near-zero (1-byte) cap, far below
        //     any item, still drains every item in serial order with none lost.
        //     The `next_serial` exemption guarantees progress; the round bound
        //     catches a livelock blow-up (not just a hard deadlock).
        let q2: Arc<dyn ItemQueue<Sequenced<Heavy>>> =
            Arc::new(CountBoundedQueue::<Sequenced<Heavy>>::new(1));
        let s2 = ReorderStage::new(q2);
        s2.set_max_overflow_bytes(1);
        let n = 16u32;
        let mut pending: Vec<(u64, Heavy)> =
            (0..n).rev().map(|i| (u64::from(i), Heavy(i))).collect();
        let mut out: Vec<u32> = Vec::new();
        let mut rounds = 0u32;
        while out.len() < n as usize {
            rounds += 1;
            assert!(rounds <= 4 * n, "tiny cap must converge (no livelock blow-up)");
            // Retry every pending push (held-slot semantics), keeping rejects.
            let mut still = Vec::new();
            for (ord, item) in pending.drain(..) {
                if let Err(rej) = s2.try_push(ord, item) {
                    still.push(rej);
                }
            }
            pending = still;
            while let Some(Heavy(v)) = s2.try_pop_in_order() {
                out.push(v);
            }
        }
        assert_eq!(out, (0..n).collect::<Vec<_>>(), "all drain in serial order, none lost");
    }

    #[test]
    fn concurrent_tiny_cap_drains_all_in_order() {
        // The load-bearing liveness test (§10.3/§10.5): at t>1, deadlock-freedom
        // is "not proven by the next_serial exemption alone" — the worry is all
        // producers wedged on rejected later-ordinal pushes while next_serial is
        // unpushed. Here N producer threads concurrently push a disjoint,
        // interleaved set of ordinals into one stage with a 1-byte cap (far
        // below any item), retrying rejects (held-slot semantics); a consumer
        // drains. All items must emerge in serial order with none lost. Progress
        // is guaranteed because the producer owning `next_serial` always has it
        // as its current (ascending) push, and `next_serial` is cap-exempt.
        use std::thread;

        // Heavy carries its ordinal as `usize` (identity) with a fixed heap
        // size, so the byte cap binds without any narrowing casts in the test.
        #[derive(Debug, PartialEq)]
        struct Heavy(usize);
        impl HeapSize for Heavy {
            fn heap_size(&self) -> usize {
                100
            }
        }

        let n_threads = 4usize;
        let per = 64usize;
        let total = n_threads * per;
        let q: Arc<dyn ItemQueue<Sequenced<Heavy>>> =
            Arc::new(CountBoundedQueue::<Sequenced<Heavy>>::new(4));
        let stage = Arc::new(ReorderStage::<Heavy>::new(q));
        stage.set_max_overflow_bytes(1); // 1 byte ≪ any item

        let producers: Vec<_> = (0..n_threads)
            .map(|t| {
                let s = Arc::clone(&stage);
                thread::spawn(move || {
                    // This thread owns ordinals {t, t+N, t+2N, ...}, pushed in
                    // ascending order (each retried until accepted).
                    let mine: Vec<usize> = (0..per).map(|k| t + k * n_threads).collect();
                    let mut i = 0;
                    while i < mine.len() {
                        let ord = mine[i];
                        match s.try_push(ord as u64, Heavy(ord)) {
                            Ok(()) => i += 1,
                            Err(_) => thread::yield_now(), // held-slot retry
                        }
                    }
                })
            })
            .collect();

        let mut out: Vec<usize> = Vec::with_capacity(total);
        while out.len() < total {
            match stage.try_pop_in_order() {
                Some(Heavy(v)) => out.push(v),
                None => thread::yield_now(),
            }
        }
        for p in producers {
            p.join().unwrap();
        }
        assert_eq!(out, (0..total).collect::<Vec<_>>(), "all drain in serial order, none lost");
    }

    #[test]
    fn drained_propagates() {
        let s = make_stage(4);
        s.try_push(0, 0).unwrap();
        s.try_push(1, 1).unwrap();
        s.mark_drained();
        // Not drained while items remain (need at least one pop to hydrate
        // the buffer/transport visibility test).
        assert!(!s.is_drained());
        let _ = s.try_pop_in_order().unwrap();
        let _ = s.try_pop_in_order().unwrap();
        assert!(s.is_drained());
    }

    #[test]
    fn backpressure_applies_after_next_serial_buffered() {
        // Drive the stage into a state where the reorder buffer holds
        // `next_serial`, then verify producers see backpressure rejection.
        //
        // Step 1: push ordinal 1 with next_serial=0 still pending.
        //   transport: [(1,1)], buf: {}, next_serial=0.
        // Step 2: a consumer call can't return anything (still waiting for 0),
        //   but it drains transport into buf as a side effect.
        //   transport: [], buf: {1}, next_serial=0.
        // Step 3: push ordinal 0 — must_accept; goes to transport.
        //   transport: [(0,0)], buf: {1}, next_serial=0.
        // Step 4: pop returns 0; advances next_serial=1; buf still has 1.
        //   transport: [], buf: {1}, next_serial=1.
        // Step 5: now buf has next_serial=1, so backpressure path is active.
        //   Pushes 2, 3 fill transport (capacity 2). Push 4 rejects.
        let s = make_stage(2);
        s.try_push(1, 100).unwrap();
        assert_eq!(s.try_pop_in_order(), None); // drains transport into buf as side effect
        s.try_push(0, 0).unwrap();
        assert_eq!(s.try_pop_in_order(), Some(0));
        // buf now contains next_serial=1; backpressure path active.
        s.try_push(2, 200).unwrap();
        s.try_push(3, 300).unwrap();
        // Transport full, buf has next_serial=1: producer 4 must be rejected.
        let result = s.try_push(4, 400);
        assert!(result.is_err(), "expected backpressure rejection");
        let (ord, item) = result.unwrap_err();
        assert_eq!((ord, item), (4, 400));
    }
}
