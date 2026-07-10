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
//!     AHashMap<u64, T>` (ahash — the ordinals are trivial monotonic `u64`
//!     keys, so the default `SipHash` buys nothing on this per-item path).
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

use ahash::AHashMap;
use parking_lot::Mutex;
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
///
/// Step authors never see this type; it is `pub` only because it appears in the
/// public `ItemQueue<Sequenced<T>>` trait bounds that ordered queues (and their
/// tests) instantiate, so callers may construct it directly when wiring queues.
pub struct Sequenced<T> {
    /// Producer-assigned monotonic ordinal used to restore emission order.
    pub ordinal: u64,
    /// The wrapped payload carried alongside its `ordinal`.
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
    /// Producer-side push metrics for an **ordered** edge, recorded at THIS
    /// boundary rather than on the transport queue. `try_push` can turn a
    /// full-transport `Err` into a must-accept stash `Ok`, so recording on the
    /// transport would miscount a stashed (accepted) item as a rejection. Here we
    /// record `record_push` on every accepted push (transport OR stash) and
    /// `record_reject` only on a true `Err`, keeping `pushed_items` /
    /// `push_rejections` accurate under producer skew. `None` when the edge is
    /// not instrumented; the pop side is recorded separately by the input handle.
    push_metrics: Option<Arc<crate::runtime::metrics::EdgeMetrics>>,
    /// Whether `push_metrics` records each push's `heap_size` (byte-bounded edge)
    /// or `0` (count/unbounded edge) — mirroring the underlying queue's own byte
    /// accounting. Set from the branch's queue kind at construction.
    record_item_bytes: bool,
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

    /// Bytes currently held in the must-accept overflow stash. The deadlock
    /// monitor sums this across branches (with the transport queues) to tell a
    /// real wedge (work stuck) from upstream starvation (everything empty); the
    /// cap-enforcement tests use it to observe stash growth.
    fn current_buffer_bytes(&self) -> u64;
}

impl<T: Send + HeapSize + 'static> ReorderCapHandle for ReorderStage<T> {
    fn set_max_overflow_bytes(&self, bytes: u64) {
        // Set once before workers spawn; `Relaxed` is sufficient because the
        // worker pool's spawn establishes the happens-before edge, and the
        // value is only read on the must-accept slow path under `state.lock()`.
        self.max_overflow_bytes.store(bytes, Ordering::Relaxed);
    }

    fn current_buffer_bytes(&self) -> u64 {
        self.state.lock().buffer_bytes
    }
}

#[cfg(test)]
impl<T: Send + HeapSize + 'static> ReorderStage<T> {
    /// Read the current overflow cap (`u64::MAX` = unbounded). Test-only —
    /// lets the budget-wiring test assert `apply_initial_queue_budget` set it.
    pub(crate) fn current_max_overflow_bytes(&self) -> u64 {
        self.max_overflow_bytes.load(Ordering::Relaxed)
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
    buffer: AHashMap<u64, (T, u64)>,
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
                buffer: AHashMap::new(),
                buffer_bytes: 0,
                next_serial: 0,
            }),
            next_serial_buffered: AtomicBool::new(false),
            drained_observed: AtomicBool::new(false),
            max_overflow_bytes: AtomicU64::new(REORDER_OVERFLOW_UNBOUNDED),
            push_metrics: None,
            record_item_bytes: false,
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
                buffer: AHashMap::new(),
                buffer_bytes: 0,
                next_serial: 0,
            }),
            next_serial_buffered: AtomicBool::new(false),
            drained_observed: AtomicBool::new(false),
            max_overflow_bytes: AtomicU64::new(max_bytes),
            push_metrics: None,
            record_item_bytes: false,
        }
    }

    /// Attach producer-side push metrics recorded at the `ReorderStage` boundary
    /// (see [`push_metrics`](Self::push_metrics)). `record_item_bytes` is `true`
    /// for a byte-bounded edge (record each push's `heap_size`) and `false` for a
    /// count/unbounded edge (record `0`), matching the underlying queue's byte
    /// accounting. Called once by the ordered-branch builder.
    #[must_use]
    pub fn with_push_metrics(
        mut self,
        push_metrics: Option<Arc<crate::runtime::metrics::EdgeMetrics>>,
        record_item_bytes: bool,
    ) -> Self {
        self.push_metrics = push_metrics;
        self.record_item_bytes = record_item_bytes;
        self
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
        // Record push-side metrics at THIS boundary (see `push_metrics`). The
        // must-accept path turns a full-transport `Err` into a stash `Ok`, so the
        // transport queue can't tell a stashed (accepted) push from a reject —
        // only the final `Result` here can. `heap_size()` is read before `item`
        // moves into the inner push. Every `Ok` (transport OR stash) is a push;
        // every `Err` (backpressure or stash-cap held) is a reject.
        //
        // Gate the `heap_size()` call on metrics being present: on the default
        // instrumentation-off path (`push_metrics == None`) the byte figure is
        // never recorded, so computing it would be pure hot-path overhead.
        let bytes = if self.push_metrics.is_some() && self.record_item_bytes {
            item.heap_size() as u64
        } else {
            0
        };
        let result = self.try_push_inner(ordinal, item);
        if let Some(m) = &self.push_metrics {
            match &result {
                Ok(()) => m.record_push(bytes),
                Err(_) => m.record_reject(),
            }
        }
        result
    }

    /// Inner push: the transport / must-accept-stash decision, without metrics.
    /// See [`try_push`](Self::try_push) for the public contract; metrics are
    /// recorded there so a stashed push is not miscounted as a rejection.
    fn try_push_inner(&self, ordinal: u64, item: T) -> Result<(), (u64, T)> {
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
                    // `usize → u64` is a lossless widen on every supported
                    // (≤64-bit) target, so this cast never truncates.
                    let item_bytes = seq.item.heap_size() as u64;
                    // A duplicate ordinal would silently drop the buffered item and
                    // leak its bytes into `buffer_bytes`. `ByOrdinal` serials are
                    // unique by construction; a hit here means a `ByItemOrdinal`
                    // upstream emitted two items with the same serial (a step bug).
                    // Fail loud in release too (like the `next_serial` overflow guard
                    // below): silently dropping a buffered record is a data-integrity bug.
                    assert!(
                        !state.buffer.contains_key(&seq.ordinal),
                        "duplicate reorder ordinal {} — ByItemOrdinal upstream serials must be unique",
                        seq.ordinal
                    );
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
    ///
    /// # Panics
    ///
    /// Panics if the in-order ordinal counter would overflow `u64` (i.e.
    /// `next_serial == u64::MAX`). Branch ordinals start at 0 and increment, so
    /// this is unreachable on any real workload (~1.8e19 items on one edge); the
    /// guard exists only to fail loudly rather than silently wrap and wait
    /// forever for ordinal 0.
    pub fn try_pop_in_order(&self) -> Option<T> {
        self.try_pop_in_order_reporting_blocked().0
    }

    /// Like [`try_pop_in_order`](Self::try_pop_in_order), but also reports whether
    /// a `None` result was *reorder-blocked* — later ordinals are buffered while
    /// the stage waits for an earlier one — rather than genuinely drained. The
    /// pop and the reorder-blocked check are computed under the SAME state lock,
    /// so a shared (`Parallel`) consumer cannot observe a torn `(item, blocked)`
    /// pair that would skew the `pop_empties` starvation metric. The flag is
    /// always `false` when an item is returned.
    ///
    /// # Panics
    ///
    /// Panics if the in-order ordinal counter would overflow `u64` (i.e.
    /// `next_serial == u64::MAX`). Branch ordinals start at 0 and increment, so
    /// this is unreachable on any real workload (~1.8e19 items on one edge); the
    /// guard exists only to fail loudly rather than silently wrap and wait
    /// forever for ordinal 0.
    pub fn try_pop_in_order_reporting_blocked(&self) -> (Option<T>, bool) {
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
                // `usize → u64` is a lossless widen on every supported
                // (≤64-bit) target, so this cast never truncates.
                let bytes = seq.item.heap_size() as u64;
                // See the drain-path insert above: a duplicate ordinal here would
                // silently drop the buffered item and corrupt `buffer_bytes`.
                // Fail loud in release too (like the `next_serial` overflow guard):
                // silently dropping a buffered record is a data-integrity bug.
                assert!(
                    !state.buffer.contains_key(&seq.ordinal),
                    "duplicate reorder ordinal {} — ByItemOrdinal upstream serials must be unique",
                    seq.ordinal
                );
                state.buffer.insert(seq.ordinal, (seq.item, bytes));
                state.buffer_bytes = state.buffer_bytes.saturating_add(bytes);
                if state.buffer.contains_key(&next) {
                    break;
                }
            }
        }

        if let Some((item, item_bytes)) = state.buffer.remove(&next) {
            state.buffer_bytes = state.buffer_bytes.saturating_sub(item_bytes);
            // Advance the in-order cursor. `ByOrdinal` ordinals come from an
            // `AtomicU64` allocator and `ByItemOrdinal` from upstream item
            // ordinals; both start at 0 and increment, so `u64::MAX` is
            // unreachable on any real workload (~1.8e19 items on one edge).
            // `checked_add` makes that invariant explicit: a wrap here would
            // silently wait forever for ordinal 0, so we fail loudly instead.
            let new_next = next.checked_add(1).expect("reorder ordinal overflow (next_serial)");
            state.next_serial = new_next;
            // Update the cache: is the NEW next_serial in the buffer?
            // Producers reading post-update see the right state.
            let new_buffered = state.buffer.contains_key(&new_next);
            self.next_serial_buffered.store(new_buffered, Ordering::Release);
            (Some(item), false)
        } else {
            // No in-order item. Reorder-blocked (backlog, not starvation) iff the
            // buffer still holds out-of-order items awaiting an earlier ordinal —
            // `next` was just confirmed absent, so any remaining entry is a later
            // ordinal. Computed here under the same lock as the pop above.
            let reorder_blocked = !state.buffer.is_empty();
            (None, reorder_blocked)
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
    use crate::queues::CountBoundedQueue;

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

    /// Q1 (audit D4): a `ByItemOrdinal` upstream supplies each item's serial, so
    /// a buggy step could emit two items with the same ordinal. Inserting the
    /// second would silently drop the first and leak its bytes into
    /// `buffer_bytes`; the always-on `assert!` must turn that into a loud failure
    /// (in release builds too).
    #[test]
    #[should_panic(expected = "duplicate reorder ordinal")]
    fn duplicate_ordinal_trips_assert() {
        let s = make_stage(8);
        s.try_push(1, 100).unwrap();
        s.try_push(1, 200).unwrap(); // duplicate ordinal (simulated upstream bug)
        s.try_push(0, 0).unwrap();
        // Draining moves both ordinal-1 items from the transport into the
        // in-order buffer; the second insert hits the duplicate guard.
        let _ = s.try_pop_in_order();
    }

    /// Sibling of `duplicate_ordinal_trips_assert`, covering the OTHER always-on
    /// duplicate-ordinal guard: the must-accept *stash-insert* path in
    /// `try_push_inner` (transport full → overflow into `state.buffer`), not the
    /// transport-drain path in `try_pop_in_order_reporting_blocked`.
    ///
    /// With transport capacity 1 and `next_serial` (0) still absent, every push
    /// is must-accept: the first fills the transport, and each later push finds
    /// the transport full and overflows into the stash. The second and third
    /// pushes both carry ordinal 1, so the second stashes it and the third trips
    /// the stash-insert `assert!` — no `try_pop_in_order` runs, so the drain-path
    /// guard is never reached. This test fails if that stash-insert guard is
    /// removed.
    #[test]
    #[should_panic(expected = "duplicate reorder ordinal")]
    fn duplicate_ordinal_trips_stash_assert() {
        let s = make_stage(1); // capacity 1 so the transport fills and later items stash
        s.try_push(1, 100).unwrap(); // ordinal 1 -> transport (fills capacity-1 transport)
        s.try_push(1, 200).unwrap(); // transport full -> ordinal 1 stashed into buffer
        // Transport still full and ordinal 1 already in the stash: the overflow
        // insert hits the duplicate guard on the push path.
        let _ = s.try_push(1, 300);
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
    fn reporting_pop_flags_reorder_blocked_vs_drained() {
        // The reporting pop returns the reorder-blocked flag from the same locked
        // path as the pop itself, so a shared consumer sees a consistent pair.
        let s = make_stage(8);
        // Buffer a later ordinal while ordinal 0 is absent → reorder-blocked.
        s.try_push(1, 100).unwrap();
        assert_eq!(s.try_pop_in_order_reporting_blocked(), (None, true), "blocked, not drained");
        // Ordinal 0 arrives; the pop yields it and is not blocked.
        s.try_push(0, 0).unwrap();
        assert_eq!(s.try_pop_in_order_reporting_blocked(), (Some(0), false));
        assert_eq!(s.try_pop_in_order_reporting_blocked(), (Some(100), false));
        // Genuinely drained now → None and NOT reorder-blocked (a true empty pop).
        assert_eq!(s.try_pop_in_order_reporting_blocked(), (None, false), "drained, not blocked");
    }

    #[test]
    fn push_metrics_count_stashed_item_as_push_not_reject() {
        use crate::runtime::metrics::EdgeMetrics;
        // Transport capacity 1, next_serial (0) absent → must-accept. The first
        // push fills the transport; the second must-accept-overflows into the
        // stash, returning Ok. Recording at the ReorderStage boundary must count
        // BOTH as pushes and NEITHER as a rejection — the transport's internal
        // `Err` on the stashed push is not a real backpressure event. (Recording
        // on the transport, as before, miscounted the stashed item as a reject.)
        let m = EdgeMetrics::new();
        let q: Arc<dyn ItemQueue<Sequenced<u32>>> = Arc::new(CountBoundedQueue::new(1));
        let stage = ReorderStage::new(q).with_push_metrics(Some(Arc::clone(&m)), false);
        assert!(stage.try_push(1, 10).is_ok(), "ordinal 1 into transport");
        assert!(stage.try_push(2, 20).is_ok(), "ordinal 2 stashed (transport full)");
        let s = m.snapshot();
        assert_eq!(s.pushed_items, 2, "both accepted pushes counted");
        assert_eq!(s.push_rejections, 0, "a stashed push is not a rejection");
        assert_eq!(s.pushed_bytes, 0, "count-bounded edge records 0 push bytes");
    }

    #[test]
    fn push_metrics_record_item_bytes_when_byte_bounded() {
        use crate::runtime::metrics::EdgeMetrics;
        #[derive(Debug)]
        struct Heavy;
        impl HeapSize for Heavy {
            fn heap_size(&self) -> usize {
                100
            }
        }
        // With `record_item_bytes = true` (byte-bounded edge), each accepted push
        // records its `heap_size` — so `pushed_bytes` tracks the bare `T`, the
        // same size the pop side records.
        let m = EdgeMetrics::new();
        let q: Arc<dyn ItemQueue<Sequenced<Heavy>>> = Arc::new(CountBoundedQueue::new(8));
        let stage = ReorderStage::new(q).with_push_metrics(Some(Arc::clone(&m)), true);
        assert!(stage.try_push(0, Heavy).is_ok());
        assert!(stage.try_push(1, Heavy).is_ok());
        let s = m.snapshot();
        assert_eq!(s.pushed_items, 2);
        assert_eq!(s.pushed_bytes, 200, "byte-bounded edge records heap_size per push");
    }

    #[test]
    fn push_without_metrics_skips_heap_size() {
        use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
        // On the default instrumentation-off path (`push_metrics == None`) the
        // push wrapper must NOT compute `heap_size` — that byte figure is only
        // needed to record a push, so computing it would be pure hot-path cost.
        #[derive(Debug)]
        struct Counted(Arc<AtomicUsize>);
        impl HeapSize for Counted {
            fn heap_size(&self) -> usize {
                self.0.fetch_add(1, AtomicOrdering::Relaxed);
                0
            }
        }
        let calls = Arc::new(AtomicUsize::new(0));
        let q: Arc<dyn ItemQueue<Sequenced<Counted>>> = Arc::new(CountBoundedQueue::new(8));
        // Byte-bounded edge (`record_item_bytes = true`) but no push metrics.
        let stage = ReorderStage::new(q).with_push_metrics(None, true);
        // ordinal 0 == next_serial → transport accepts (no stash), so the only
        // `heap_size` call would be the wrapper's — which the gate must skip.
        stage.try_push(0, Counted(Arc::clone(&calls))).unwrap();
        assert_eq!(
            calls.load(AtomicOrdering::Relaxed),
            0,
            "heap_size must not be computed when push metrics are disabled"
        );
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

        // Deadline guard: this test proves deadlock-freedom, so a liveness
        // regression must *fail fast* rather than hang the whole `nextest` run
        // until a global harness timeout (if any) fires. 30s is generous —
        // 256 tiny items drain in milliseconds when healthy — so it only trips
        // on a genuine wedge, not on a slow CI host. Mirrors the
        // bounded-convergence guard the single-threaded sibling uses.
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(30);
        let mut out: Vec<usize> = Vec::with_capacity(total);
        while out.len() < total {
            assert!(
                std::time::Instant::now() < deadline,
                "concurrent reorder drain did not complete within 30s ({} of {total} drained) — \
                 likely a deadlock/livelock regression",
                out.len(),
            );
            match stage.try_pop_in_order() {
                Some(Heavy(v)) => out.push(v),
                None => thread::yield_now(),
            }
        }
        // Bound producer completion with the same deadline. The drain loop above
        // only proves the *consumer* made progress; if a regression let
        // `out.len()` reach `total` (e.g. a double-emit) while a producer is
        // still stuck in its `try_push` retry loop, an unbounded `join()` would
        // hang the run forever. Poll `is_finished()` against the deadline so a
        // stuck producer fails fast instead.
        for p in producers {
            while !p.is_finished() {
                assert!(
                    std::time::Instant::now() < deadline,
                    "a producer thread did not finish within 30s — likely a stuck try_push \
                     retry loop (liveness regression)",
                );
                thread::yield_now();
            }
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
