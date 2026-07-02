//! Transport-layer queue trait + three concrete impls.
//!
//! Concerns: pure transport (push/pop, drained signal). **Not** ordering —
//! see `core/reorder.rs` for the `ReorderStage<T>` operator that adds
//! ordinal-based reordering on top of any `ItemQueue<T>`. **Not** memory
//! bookkeeping at the trait level — `ByteBoundedQueue<T: HeapSize>` is a
//! concrete impl that knows about heap size, but the trait surface is
//! type-uniform.
//!
//! Backpressure is expressed as `try_push -> Result<(), T>`: `Err(item)`
//! returns the rejected item back to the producer (which holds it in a
//! `HeldSlot<T>` and re-pushes on the next worker iteration). No blocking,
//! no awaiting — pure non-blocking surface.
//!
//! Drained-signal protocol:
//!   - Producer (output side) calls `mark_drained()` exactly once when the
//!     producing step returns `StepOutcome::Finished` (counter-gated for
//!     `Parallel` so only the last clone closes the shared queue). Subsequent
//!     `try_push` calls panic in debug builds (a contract violation: producer
//!     pushed after declaring done).
//!   - Consumer (input side) checks `is_drained() && is_empty()` to detect
//!     end-of-stream. Once both are true, no further items will arrive.

use crossbeam_queue::{ArrayQueue, SegQueue};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering, fence};
use std::sync::{Arc, Condvar, Mutex};
use std::time::Instant;

use super::item::HeapSize;
use super::runtime::metrics::EdgeMetrics;

/// Elapsed ns since `since`, saturating (park durations never realistically
/// overflow `u64` ns ≈ 584 years).
fn park_ns(since: Instant) -> u64 {
    u64::try_from(since.elapsed().as_nanos()).unwrap_or(u64::MAX)
}

// ─────────────────────────────────────────────────────────────────────────────
// Blocking facet (StepKind::Detached bridge)
// ─────────────────────────────────────────────────────────────────────────────

/// Park/notify rendezvous bolted onto a lock-free bounded queue so a
/// [`StepKind::Detached`](crate::step::StepKind::Detached) thread can *block*
/// on an empty input or a full output instead of cooperatively re-polling like
/// a pool worker.
///
/// **Opt-in and hot-path-free.** A notifier exists only on a queue constructed
/// `*_blocking`; non-blocking queues hold `blocking: None` and their
/// `try_push`/`try_pop` skip every notifier touch via one predictable nullable
/// branch — no atomic, no fence, no lock. Pool↔pool edges therefore pay
/// nothing.
///
/// **Lost-wakeup protocol (M5).** Over a lock-free `ArrayQueue`, a naive
/// "publish then notify-if-waiter" races: a producer can publish and read
/// `waiters == 0` in the same window a consumer reads "empty" and decides to
/// park, losing the wakeup. We close it with a Dekker-style `SeqCst` store-load
/// fence pair. The **producer** (`note_item_available`, called after a push or
/// `mark_drained`) publishes, then `fence(SeqCst)`, then loads the waiter count
/// and notifies only if `> 0`. The **consumer** (`blocking_pop`), under the
/// condvar mutex, registers as a waiter, then `fence(SeqCst)`, then re-checks
/// `try_pop`/`is_drained`, and parks only if still empty-and-open. The two
/// `SeqCst` fences guarantee at least one side observes the other, so a parked
/// consumer always has a pending notify (exercised by the high-iteration
/// `blocking_pop_stress_no_lost_wakeup` test; a formal `loom` model is added at
/// the lever-2 gate). The producer-parks-on-full side
/// (`note_space_available` / `blocking_push`) is symmetric.
struct BlockingNotifier {
    mutex: Mutex<()>,
    /// Signalled when an item is published or the queue is drained.
    item_cv: Condvar,
    /// Signalled when a pop frees space.
    space_cv: Condvar,
    /// Consumers currently parked (or about to park) on `item_cv`.
    item_waiters: AtomicUsize,
    /// Producers currently parked (or about to park) on `space_cv`.
    space_waiters: AtomicUsize,
    /// Deep-mode (`--pipeline-trace deep`) park-time accounting: total ns a
    /// consumer spent blocked on an empty queue (starvation latency) and a
    /// producer spent blocked on a full queue (backpressure latency). Measured
    /// only around the actual `wait` calls, so the no-waiter fast path pays
    /// nothing. The exact starvation/backpressure *time* for a Detached edge,
    /// where the derived (occupancy/throughput) estimate is weakest.
    park_on_empty_ns: AtomicU64,
    park_on_full_ns: AtomicU64,
}

impl BlockingNotifier {
    fn new() -> Self {
        Self {
            mutex: Mutex::new(()),
            item_cv: Condvar::new(),
            space_cv: Condvar::new(),
            item_waiters: AtomicUsize::new(0),
            space_waiters: AtomicUsize::new(0),
            park_on_empty_ns: AtomicU64::new(0),
            park_on_full_ns: AtomicU64::new(0),
        }
    }

    /// `(park_on_empty_ns, park_on_full_ns)` — exact consumer-starvation and
    /// producer-backpressure time accumulated across all `blocking_pop` /
    /// `blocking_push` parks. Surfaced into the edge report once Detached edges
    /// are wired (lever 2).
    #[cfg(test)]
    fn park_times(&self) -> (u64, u64) {
        (
            self.park_on_empty_ns.load(Ordering::Relaxed),
            self.park_on_full_ns.load(Ordering::Relaxed),
        )
    }

    /// Producer side: call right after publishing an item (or marking drained).
    /// Wakes a consumer parked in `blocking_pop`. `fence(SeqCst)` before the
    /// waiter load is the producer half of the Dekker pair.
    fn note_item_available(&self) {
        fence(Ordering::SeqCst);
        if self.item_waiters.load(Ordering::SeqCst) > 0 {
            let _g = self.mutex.lock().expect("blocking notifier mutex poisoned");
            self.item_cv.notify_all();
        }
    }

    /// Consumer side: call right after a pop frees space. Wakes a producer
    /// parked in `blocking_push`.
    fn note_space_available(&self) {
        fence(Ordering::SeqCst);
        if self.space_waiters.load(Ordering::SeqCst) > 0 {
            let _g = self.mutex.lock().expect("blocking notifier mutex poisoned");
            self.space_cv.notify_all();
        }
    }

    /// Block until `try_pop` yields an item, or return `None` once the queue is
    /// drained and empty. `try_pop`/`is_empty`/`is_drained` close over the
    /// owning queue's non-blocking primitives.
    fn blocking_pop<T>(
        &self,
        try_pop: impl Fn() -> Option<T>,
        is_empty: impl Fn() -> bool,
        is_drained: impl Fn() -> bool,
    ) -> Option<T> {
        loop {
            if let Some(item) = try_pop() {
                return Some(item);
            }
            // Empty on the fast path. Register as a waiter and re-check under
            // the fence before committing to `wait` (Dekker consumer half).
            let guard = self.mutex.lock().expect("blocking notifier mutex poisoned");
            self.item_waiters.fetch_add(1, Ordering::SeqCst);
            fence(Ordering::SeqCst);
            if !is_empty() {
                // A producer published between the fast-path `try_pop` and now;
                // don't park — loop back and pop it.
                self.item_waiters.fetch_sub(1, Ordering::SeqCst);
                continue;
            }
            if is_drained() {
                // Drained and empty → end of stream.
                self.item_waiters.fetch_sub(1, Ordering::SeqCst);
                return None;
            }
            // Genuinely empty and open: park. `wait` atomically releases the
            // mutex; a producer's `note_item_available` (which locks the same
            // mutex when `item_waiters > 0`) wakes us. Time the park (deep-mode
            // starvation latency); only the actual wait pays the `Instant` cost.
            let parked = Instant::now();
            let _unused = self.item_cv.wait(guard).expect("blocking notifier mutex poisoned");
            self.park_on_empty_ns.fetch_add(park_ns(parked), Ordering::Relaxed);
            self.item_waiters.fetch_sub(1, Ordering::SeqCst);
        }
    }

    /// Block until `try_push` accepts `item` (parks on a full queue, woken by a
    /// consumer's pop). Assumes the consumer keeps draining (true for the sort
    /// merge→writer edge: the writer pops until the merge marks drained).
    fn blocking_push<T>(&self, item: T, try_push: impl Fn(T) -> Result<(), T>) {
        let mut pending = item;
        loop {
            match try_push(pending) {
                Ok(()) => return,
                Err(returned) => pending = returned,
            }
            // Full. Register as a space waiter and retry under the fence before
            // parking (Dekker producer-parks half).
            let guard = self.mutex.lock().expect("blocking notifier mutex poisoned");
            self.space_waiters.fetch_add(1, Ordering::SeqCst);
            fence(Ordering::SeqCst);
            match try_push(pending) {
                Ok(()) => {
                    self.space_waiters.fetch_sub(1, Ordering::SeqCst);
                    return;
                }
                Err(returned) => pending = returned,
            }
            let parked = Instant::now();
            let _unused = self.space_cv.wait(guard).expect("blocking notifier mutex poisoned");
            self.park_on_full_ns.fetch_add(park_ns(parked), Ordering::Relaxed);
            self.space_waiters.fetch_sub(1, Ordering::SeqCst);
        }
    }
}

/// Blocking facet for [`StepKind::Detached`](crate::step::StepKind::Detached)
/// threads, layered on top of [`ItemQueue`]. Only the two bounded queues
/// implement it, and only when constructed `*_blocking`; a Detached step's
/// dedicated thread uses these instead of the pool's cooperative re-poll.
pub trait BlockingQueue<T: Send + 'static>: ItemQueue<T> {
    /// Block until an item is available, or `None` once drained and empty.
    fn blocking_pop(&self) -> Option<T>;

    /// Block until `item` is accepted (parks on a full queue).
    fn blocking_push(&self, item: T);
}

/// Transport-layer queue trait. Type-uniform across queue impls: the
/// `try_push` surface accepts any `T` regardless of whether the impl uses
/// item-count or memory bookkeeping internally.
///
/// `Send + Sync`: queues are shared between worker threads via `Arc`.
pub trait ItemQueue<T: Send + 'static>: Send + Sync {
    /// Non-blocking push. `Err(item)` returns the rejected item to the
    /// caller; the framework holds it in a `HeldSlot<T>` and retries.
    ///
    /// # Errors
    ///
    /// Returns `Err(item)` when the queue is at its backpressure limit
    /// (item-count or byte-budget, depending on the impl).
    fn try_push(&self, item: T) -> Result<(), T>;

    /// Non-blocking pop. `None` means the queue is currently empty (which
    /// is *not* the same as drained — combine with `is_drained()`).
    fn try_pop(&self) -> Option<T>;

    /// True when no items are currently buffered. May race with concurrent
    /// pushes; consumers that need a quiescent check combine with
    /// `is_drained()`.
    fn is_empty(&self) -> bool;

    /// Mark the queue drained (producer-side: "I'm done pushing"). Idempotent.
    /// In debug builds a `try_push` after `mark_drained` panics.
    fn mark_drained(&self);

    /// True if `mark_drained` has been called.
    fn is_drained(&self) -> bool;
}

/// One entry per output branch in `StepProfile::output_queues`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QueueSpec {
    /// Item-count bounded. `try_push` rejects when `len() >= capacity`.
    /// Best for fixed-size items (parsed records, compressed BGZF blocks
    /// of known size, etc).
    CountBounded { capacity: usize },
    /// Memory-bounded. `try_push` rejects when adding the item would push
    /// the running byte counter past `limit_bytes`. Requires `T: HeapSize`.
    /// Best for variable-size BAM batches and FASTQ batches.
    ByteBounded { limit_bytes: u64 },
    /// No backpressure. `try_push` always succeeds. Use only when the
    /// branch is naturally rate-limited upstream (e.g., a header-once
    /// emit on pipeline start).
    Unbounded,
}

// ─────────────────────────────────────────────────────────────────────────────
// CountBoundedQueue
// ─────────────────────────────────────────────────────────────────────────────

/// Item-count bounded transport. Backed by `crossbeam_queue::ArrayQueue<T>`.
pub struct CountBoundedQueue<T: Send + 'static> {
    inner: ArrayQueue<T>,
    drained: AtomicBool,
    /// `Some` only on a `new_blocking` queue feeding/feeding-from a Detached
    /// step; `None` keeps `try_push`/`try_pop` notifier-free for pool edges.
    blocking: Option<BlockingNotifier>,
    /// `Some` only on an instrumented edge (`--pipeline-trace`); `None` keeps the
    /// hot path metric-free. Producer-push counts are recorded here; consumer-pop
    /// counts are recorded at the `BranchInputHandle` (see `handles.rs`).
    metrics: Option<Arc<EdgeMetrics>>,
}

impl<T: Send + 'static> CountBoundedQueue<T> {
    /// Construct a count-bounded transport with the given capacity.
    ///
    /// # Panics
    ///
    /// Panics if `capacity == 0` (a zero-capacity queue would always reject).
    #[must_use]
    pub fn new(capacity: usize) -> Self {
        Self::build(capacity, None, None)
    }

    /// Like [`new`](Self::new) but with the blocking facet enabled, for an edge
    /// adjacent to a [`StepKind::Detached`](crate::step::StepKind::Detached)
    /// thread. The non-blocking `try_*` surface is unchanged; only the extra
    /// notifier touches (a `SeqCst` fence + waiter-count load, then a lock+notify
    /// only when a thread is actually parked) are added.
    ///
    /// # Panics
    ///
    /// Panics if `capacity == 0` (a zero-capacity queue would always reject).
    #[must_use]
    pub fn new_blocking(capacity: usize) -> Self {
        Self::build(capacity, Some(BlockingNotifier::new()), None)
    }

    /// Like [`new`](Self::new) but recording producer-push metrics into `metrics`
    /// (an instrumented edge). The non-blocking `try_*` surface is unchanged.
    ///
    /// # Panics
    ///
    /// Panics if `capacity == 0`.
    #[must_use]
    pub fn new_instrumented(capacity: usize, metrics: Arc<EdgeMetrics>) -> Self {
        Self::build(capacity, None, Some(metrics))
    }

    /// [`new`](Self::new) when `metrics` is `None`, [`new_instrumented`](Self::new_instrumented)
    /// when `Some`. Lets branch builders thread an optional metrics handle uniformly.
    ///
    /// # Panics
    ///
    /// Panics if `capacity == 0`.
    #[must_use]
    pub fn maybe_instrumented(capacity: usize, metrics: Option<Arc<EdgeMetrics>>) -> Self {
        Self::build(capacity, None, metrics)
    }

    fn build(
        capacity: usize,
        blocking: Option<BlockingNotifier>,
        metrics: Option<Arc<EdgeMetrics>>,
    ) -> Self {
        assert!(capacity > 0, "CountBoundedQueue capacity must be > 0");
        Self {
            inner: ArrayQueue::new(capacity),
            drained: AtomicBool::new(false),
            blocking,
            metrics,
        }
    }

    /// Occupancy fraction `0.0..=1.0` (`len / capacity`). For the sampler.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn depth_fraction(&self) -> f32 {
        self.inner.len() as f32 / self.inner.capacity() as f32
    }

    /// Number of items currently buffered. For the sampler.
    #[must_use]
    pub fn occupancy_len(&self) -> usize {
        self.inner.len()
    }
}

impl<T: Send + 'static> ItemQueue<T> for CountBoundedQueue<T> {
    fn try_push(&self, item: T) -> Result<(), T> {
        debug_assert!(
            !self.drained.load(Ordering::Acquire),
            "try_push after mark_drained — producer contract violation"
        );
        if let Err(item) = self.inner.push(item) {
            if let Some(m) = &self.metrics {
                m.record_reject();
            }
            return Err(item);
        }
        if let Some(m) = &self.metrics {
            m.record_push(0); // count-bounded: items only, no byte size
        }
        if let Some(n) = &self.blocking {
            n.note_item_available();
        }
        Ok(())
    }

    fn try_pop(&self) -> Option<T> {
        let item = self.inner.pop()?;
        if let Some(n) = &self.blocking {
            n.note_space_available();
        }
        Some(item)
    }

    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn mark_drained(&self) {
        self.drained.store(true, Ordering::Release);
        if let Some(n) = &self.blocking {
            n.note_item_available();
        }
    }

    fn is_drained(&self) -> bool {
        self.drained.load(Ordering::Acquire)
    }
}

impl<T: Send + 'static> BlockingQueue<T> for CountBoundedQueue<T> {
    fn blocking_pop(&self) -> Option<T> {
        let n = self.blocking.as_ref().expect("blocking_pop on a non-blocking CountBoundedQueue");
        // Pop via `try_pop` so a producer parked in `blocking_push` is woken.
        n.blocking_pop(|| self.try_pop(), || self.inner.is_empty(), || self.is_drained())
    }

    fn blocking_push(&self, item: T) {
        let n = self.blocking.as_ref().expect("blocking_push on a non-blocking CountBoundedQueue");
        // Route the push through `ItemQueue::try_push` so the item-available
        // notify still fires for a consumer parked in `blocking_pop`.
        n.blocking_push(item, |it| self.try_push(it));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ByteBoundedQueue
// ─────────────────────────────────────────────────────────────────────────────

/// Backing slot capacity for `ByteBoundedQueue`. Since the queue's real
/// gate is the byte budget, this just needs to be large enough that the
/// count never matters for any sane workload. 1024 slots is well past
/// the working set of any single pipeline edge — even for the smallest
/// items the byte cap (default 4 MiB) imposes a tighter bound.
///
/// Sized in pages of `crossbeam_queue::ArrayQueue` storage (one
/// pre-allocated slot array, no per-push allocation). Mirrors the
/// `ArrayQueue::new(queue_capacity)` strategy the legacy pipeline used — it
/// also used a fixed-capacity `ArrayQueue` everywhere for the same
/// reason: `SegQueue` allocates segments on demand under load, and
/// the resulting allocator churn shows up as `mi_*` overhead in
/// profiles (≈260 samples vs legacy on CODEC 8M).
const BYTE_BOUNDED_QUEUE_SLOT_CAPACITY: usize = 1024;

/// Memory-bounded transport. Backed by
/// `crossbeam_queue::ArrayQueue<(T, u64)>` plus an atomic byte counter.
/// `try_push` rejects when the running byte counter has already reached
/// `limit_bytes`. Requires `T: HeapSize`.
///
/// ## Concurrency / ordering
///
/// The check-then-add is two atomics, so two concurrent pushes can both
/// observe `cur < limit` and both succeed, yielding a small overshoot.
/// The next push will see the overshoot and reject; the budget is
/// enforced as "approximate within one item's worth per producer." This
/// trade-off avoids a CAS loop and is fine for backpressure semantics.
///
/// `current_bytes` is **only** a backpressure heuristic — it's not used
/// to synchronize handoff of the items themselves. The handoff is the
/// `ArrayQueue`'s job; `ArrayQueue`'s internal atomics provide the
/// happens-before relationship between `inner.push` and `inner.pop`. We
/// therefore use `Relaxed` ordering on every `current_bytes` access:
/// the worst case is a slightly stale reading of the budget, never an
/// observability violation on the items.
///
/// ## Cached size at push
///
/// The size is stored alongside the item in the inner queue
/// (`ArrayQueue<(T, u64)>`) so `try_pop` doesn't need to recompute
/// `T::heap_size()` for the budget update. For types whose
/// `heap_size()` is O(items inside) (e.g. `BatchedRawPositionGroups`,
/// `OrderedRawPositionGroup`) this avoids recomputing a O(group)
/// walk on every pop. Mirrors the legacy `ReorderBuffer<T>`'s
/// cached-size storage strategy (`(T, usize)` there; `(T, u64)` here,
/// matching `inner`'s `ArrayQueue<(T, u64)>` above).
pub struct ByteBoundedQueue<T: Send + HeapSize + 'static> {
    inner: ArrayQueue<(T, u64)>,
    current_bytes: AtomicU64,
    /// Mutable byte-budget cap. The rebalancer (when enabled via
    /// `PipelineConfig::queue_memory_total`) updates this atomic at
    /// runtime to shift budget across queues based on observed
    /// fullness. Producers read it on every `try_push`; the
    /// `Relaxed` ordering matches `current_bytes` (this is a
    /// best-effort backpressure heuristic, not a correctness gate).
    limit_bytes: AtomicU64,
    drained: AtomicBool,
    /// Per-instance one-shot guard so the "slot cap hit before byte budget"
    /// warning (see `try_push`) is emitted at most once *per queue*, not once
    /// per process. A process-global flag would silence the warning for every
    /// later queue (e.g. a second `runall` stage, or many pipelines in one
    /// long-lived host / test harness) after the first occurrence. The hot-path
    /// cost is a single relaxed swap after the first hit.
    slot_cap_warned: AtomicBool,
    /// `Some` only on a `new_blocking` queue adjacent to a Detached step; `None`
    /// keeps `try_push`/`try_pop` notifier-free for pool↔pool edges.
    blocking: Option<BlockingNotifier>,
    /// `Some` only on an instrumented edge; producer-push (items + bytes) and
    /// rejections are recorded here. Consumer-pop is recorded at the
    /// `BranchInputHandle` (see `handles.rs`).
    metrics: Option<Arc<EdgeMetrics>>,
}

impl<T: Send + HeapSize + 'static> ByteBoundedQueue<T> {
    /// Construct a byte-bounded transport with the given memory limit.
    ///
    /// # Panics
    ///
    /// Panics if `limit_bytes == 0` (a zero-budget queue would always reject).
    #[must_use]
    pub fn new(limit_bytes: u64) -> Self {
        Self::build(limit_bytes, None, None)
    }

    /// Like [`new`](Self::new) but with the blocking facet enabled, for an edge
    /// adjacent to a [`StepKind::Detached`](crate::step::StepKind::Detached)
    /// thread. The non-blocking `try_*` surface and byte-budget semantics are
    /// unchanged.
    ///
    /// # Panics
    ///
    /// Panics if `limit_bytes == 0` (a zero-budget queue would always reject).
    #[must_use]
    pub fn new_blocking(limit_bytes: u64) -> Self {
        Self::build(limit_bytes, Some(BlockingNotifier::new()), None)
    }

    /// Like [`new`](Self::new) but recording producer-push metrics (items + bytes
    /// + rejections) into `metrics`. Byte-budget semantics unchanged.
    ///
    /// # Panics
    ///
    /// Panics if `limit_bytes == 0`.
    #[must_use]
    pub fn new_instrumented(limit_bytes: u64, metrics: Arc<EdgeMetrics>) -> Self {
        Self::build(limit_bytes, None, Some(metrics))
    }

    /// [`new`](Self::new) when `metrics` is `None`, [`new_instrumented`](Self::new_instrumented)
    /// when `Some`.
    ///
    /// # Panics
    ///
    /// Panics if `limit_bytes == 0`.
    #[must_use]
    pub fn maybe_instrumented(limit_bytes: u64, metrics: Option<Arc<EdgeMetrics>>) -> Self {
        Self::build(limit_bytes, None, metrics)
    }

    fn build(
        limit_bytes: u64,
        blocking: Option<BlockingNotifier>,
        metrics: Option<Arc<EdgeMetrics>>,
    ) -> Self {
        assert!(limit_bytes > 0, "ByteBoundedQueue limit_bytes must be > 0");
        Self {
            inner: ArrayQueue::new(BYTE_BOUNDED_QUEUE_SLOT_CAPACITY),
            current_bytes: AtomicU64::new(0),
            limit_bytes: AtomicU64::new(limit_bytes),
            drained: AtomicBool::new(false),
            slot_cap_warned: AtomicBool::new(false),
            blocking,
            metrics,
        }
    }

    /// Occupancy fraction `0.0..=1.0` (`current_bytes / limit_bytes`, clamped).
    /// For the sampler.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn depth_fraction(&self) -> f32 {
        let limit = self.limit_bytes.load(Ordering::Relaxed);
        if limit == 0 {
            return 0.0;
        }
        (self.current_bytes.load(Ordering::Relaxed) as f32 / limit as f32).clamp(0.0, 1.0)
    }

    /// Number of items currently buffered. For the sampler.
    #[must_use]
    pub fn occupancy_len(&self) -> usize {
        self.inner.len()
    }

    /// Best-effort, stale-tolerant `Relaxed` read of the running byte
    /// counter. Used by the rebalancer as a budget heuristic, not as a
    /// correctness gate — it may lag a concurrent `try_push`/`try_pop`.
    #[must_use]
    pub fn current_bytes(&self) -> u64 {
        self.current_bytes.load(Ordering::Relaxed)
    }

    /// Best-effort, stale-tolerant `Relaxed` read of the byte-budget cap.
    /// A concurrent `set_limit_bytes` (rebalancer) may not yet be visible;
    /// callers use this as a heuristic, never as a correctness gate.
    #[must_use]
    pub fn limit_bytes(&self) -> u64 {
        self.limit_bytes.load(Ordering::Relaxed)
    }

    /// Update the byte-budget cap. Called by the rebalancer when
    /// reallocating budget across queues. Concurrent `try_push`es
    /// see the new cap on their next read; transient overshoot
    /// (pushes already in flight that read the old cap) is
    /// self-correcting.
    pub fn set_limit_bytes(&self, new_limit: u64) {
        self.limit_bytes.store(new_limit, Ordering::Relaxed);
    }
}

/// Type-erased handle for a byte-bounded queue. The pipeline
/// rebalancer iterates over registered handles to read fullness
/// (`current_bytes / limit_bytes`) and reallocate budget across
/// queues by calling `set_limit_bytes`. The trait deliberately
/// does not surface the queue's item type or its `ItemQueue`
/// methods — rebalancing only needs the byte counters.
pub trait BoundedQueueHandle: Send + Sync {
    /// Bytes currently held in the queue.
    fn current_bytes(&self) -> u64;
    /// Current byte-budget cap. May change between calls if a
    /// rebalancer is active.
    fn limit_bytes(&self) -> u64;
    /// Update the byte-budget cap. Concurrent producers see the
    /// new value on their next push.
    fn set_limit_bytes(&self, new_limit: u64);
}

impl<T: Send + HeapSize + 'static> BoundedQueueHandle for ByteBoundedQueue<T> {
    fn current_bytes(&self) -> u64 {
        self.current_bytes()
    }
    fn limit_bytes(&self) -> u64 {
        self.limit_bytes()
    }
    fn set_limit_bytes(&self, new_limit: u64) {
        self.set_limit_bytes(new_limit);
    }
}

impl<T: Send + HeapSize + 'static> ItemQueue<T> for ByteBoundedQueue<T> {
    fn try_push(&self, item: T) -> Result<(), T> {
        debug_assert!(
            !self.drained.load(Ordering::Acquire),
            "try_push after mark_drained — producer contract violation"
        );
        // Like the legacy `ReorderBufferState::can_proceed`, this
        // gates on `heap_bytes < limit` — accept if currently *under*
        // budget, regardless of incoming item size. Per-item-larger-than
        // -limit is a real case (busy-locus position-group batches can
        // be tens of MB while the queue limit is 4 MiB), so a strict
        // `cur + size <= limit` would deadlock the producer.
        //
        // Once `cur` reaches `limit_bytes`, subsequent pushes reject
        // until a consumer drains. Transient overshoot under concurrent
        // pushes is self-correcting on the next round.
        let cur = self.current_bytes.load(Ordering::Relaxed);
        let limit = self.limit_bytes.load(Ordering::Relaxed);
        if cur >= limit {
            if let Some(m) = &self.metrics {
                m.record_reject();
            }
            return Err(item);
        }
        let size = item.heap_size() as u64;
        // Reserve bytes before pushing so a concurrent consumer cannot pop and
        // decrement the counter before we add our share, which would cause the
        // counter to underflow and create permanent false backpressure.
        self.current_bytes.fetch_add(size, Ordering::Relaxed);
        // ArrayQueue::push returns Err((item, size)) on full; roll back the
        // reservation and return the item to the caller for retry. (In practice
        // the slot cap should never be hit before the byte budget triggers a
        // reject above, but defend against it anyway.)
        match self.inner.push((item, size)) {
            Ok(()) => {
                if let Some(m) = &self.metrics {
                    m.record_push(size);
                }
                if let Some(n) = &self.blocking {
                    n.note_item_available();
                }
                Ok(())
            }
            Err((item, _size)) => {
                // Roll back the byte reservation — the item never entered the queue.
                self.current_bytes.fetch_sub(size, Ordering::Relaxed);
                if let Some(m) = &self.metrics {
                    m.record_reject();
                }
                // The fixed 1024-slot backing was hit before the byte budget.
                // This degrades byte-backpressure into a hard count cap for
                // small items (heap_size ≲ limit/1024) — correctness is
                // preserved (the producer retries) but throughput silently
                // suffers. Surface it once so it is observable rather than a
                // silent foot-gun; near-zero cost after the first hit.
                if !self.slot_cap_warned.swap(true, Ordering::Relaxed) {
                    log::warn!(
                        "ByteBoundedQueue hit its {BYTE_BOUNDED_QUEUE_SLOT_CAPACITY}-slot count \
                         cap before the byte budget; small items are degrading byte-backpressure \
                         into a count cap (throughput, not correctness, is affected)."
                    );
                }
                Err(item)
            }
        }
    }

    fn try_pop(&self) -> Option<T> {
        let (item, size) = self.inner.pop()?;
        self.current_bytes.fetch_sub(size, Ordering::Relaxed);
        if let Some(n) = &self.blocking {
            n.note_space_available();
        }
        Some(item)
    }

    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn mark_drained(&self) {
        self.drained.store(true, Ordering::Release);
        if let Some(n) = &self.blocking {
            n.note_item_available();
        }
    }

    fn is_drained(&self) -> bool {
        self.drained.load(Ordering::Acquire)
    }
}

impl<T: Send + HeapSize + 'static> BlockingQueue<T> for ByteBoundedQueue<T> {
    fn blocking_pop(&self) -> Option<T> {
        let n = self.blocking.as_ref().expect("blocking_pop on a non-blocking ByteBoundedQueue");
        n.blocking_pop(|| self.try_pop(), || self.inner.is_empty(), || self.is_drained())
    }

    fn blocking_push(&self, item: T) {
        let n = self.blocking.as_ref().expect("blocking_push on a non-blocking ByteBoundedQueue");
        n.blocking_push(item, |it| self.try_push(it));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// UnboundedQueue
// ─────────────────────────────────────────────────────────────────────────────

/// Unbounded transport. `try_push` always succeeds. Backed by `SegQueue<T>`.
pub struct UnboundedQueue<T: Send + 'static> {
    inner: SegQueue<T>,
    drained: AtomicBool,
    /// `Some` only on an instrumented edge; producer-push items are recorded
    /// here (unbounded → never rejects, no byte tracking). Consumer-pop is at the
    /// `BranchInputHandle`.
    metrics: Option<Arc<EdgeMetrics>>,
}

impl<T: Send + 'static> UnboundedQueue<T> {
    #[must_use]
    pub fn new() -> Self {
        Self { inner: SegQueue::new(), drained: AtomicBool::new(false), metrics: None }
    }

    /// Like [`new`](Self::new) but recording producer-push item counts into
    /// `metrics`. Unbounded edges have no byte budget and never reject; depth is
    /// reported as raw length only.
    #[must_use]
    pub fn new_instrumented(metrics: Arc<EdgeMetrics>) -> Self {
        Self { inner: SegQueue::new(), drained: AtomicBool::new(false), metrics: Some(metrics) }
    }

    /// [`new`](Self::new) when `metrics` is `None`, [`new_instrumented`](Self::new_instrumented)
    /// when `Some`.
    #[must_use]
    pub fn maybe_instrumented(metrics: Option<Arc<EdgeMetrics>>) -> Self {
        Self { inner: SegQueue::new(), drained: AtomicBool::new(false), metrics }
    }

    /// Number of items currently buffered. For the sampler. (Unbounded edges
    /// have no meaningful depth *fraction*; the sampler reports raw length.)
    #[must_use]
    pub fn occupancy_len(&self) -> usize {
        self.inner.len()
    }
}

impl<T: Send + 'static> Default for UnboundedQueue<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Send + 'static> ItemQueue<T> for UnboundedQueue<T> {
    fn try_push(&self, item: T) -> Result<(), T> {
        debug_assert!(
            !self.drained.load(Ordering::Acquire),
            "try_push after mark_drained — producer contract violation"
        );
        self.inner.push(item);
        if let Some(m) = &self.metrics {
            m.record_push(0);
        }
        Ok(())
    }

    fn try_pop(&self) -> Option<T> {
        self.inner.pop()
    }

    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn mark_drained(&self) {
        self.drained.store(true, Ordering::Release);
    }

    fn is_drained(&self) -> bool {
        self.drained.load(Ordering::Acquire)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn count_bounded_round_trip() {
        let q: Arc<dyn ItemQueue<u32>> = Arc::new(CountBoundedQueue::new(2));
        assert!(q.try_push(1).is_ok());
        assert!(q.try_push(2).is_ok());
        assert_eq!(q.try_push(3), Err(3));
        assert_eq!(q.try_pop(), Some(1));
        assert_eq!(q.try_pop(), Some(2));
        assert_eq!(q.try_pop(), None);
    }

    #[test]
    fn count_bounded_drain_signal() {
        let q = CountBoundedQueue::<u32>::new(4);
        assert!(!q.is_drained());
        q.mark_drained();
        assert!(q.is_drained());
    }

    #[derive(Debug)]
    struct Heavy(Vec<u8>);
    impl HeapSize for Heavy {
        fn heap_size(&self) -> usize {
            self.0.len()
        }
    }

    #[test]
    fn byte_bounded_slot_cap_reject_is_observable() {
        // Tiny (0-byte heap) items with a huge byte limit: the byte budget is
        // never reached, so the fixed slot backing becomes the binding cap.
        // The first SLOT_CAPACITY pushes succeed; the next rejects on the slot
        // cap even though current_bytes is far below the limit. Regression for
        // F02 — this path silently degraded byte-backpressure into a count cap;
        // it is now warn-once observable, and this pins the reject behaviour.
        let q = ByteBoundedQueue::<Heavy>::new(1_000_000);
        for i in 0..BYTE_BOUNDED_QUEUE_SLOT_CAPACITY {
            assert!(q.try_push(Heavy(Vec::new())).is_ok(), "push {i} within slot cap");
        }
        assert_eq!(q.current_bytes(), 0, "0-byte items leave the byte budget unused");
        assert!(
            q.try_push(Heavy(Vec::new())).is_err(),
            "push #{} must reject on the slot cap, not the byte budget",
            BYTE_BOUNDED_QUEUE_SLOT_CAPACITY + 1
        );
    }

    #[test]
    fn slot_cap_warn_flag_is_per_instance_not_process_global() {
        // The "slot cap hit before byte budget" warn-once guard lives on the
        // queue instance, so a second queue (e.g. a later runall stage, or a new
        // pipeline in a long-lived host) still warns on its own first hit — the
        // signal is not silenced process-wide by an earlier queue.
        let fill_to_slot_cap = |q: &ByteBoundedQueue<Heavy>| {
            for _ in 0..BYTE_BOUNDED_QUEUE_SLOT_CAPACITY {
                q.try_push(Heavy(Vec::new())).expect("push within slot cap");
            }
            // This push trips the slot cap and (first time) sets the flag.
            assert!(q.try_push(Heavy(Vec::new())).is_err(), "push must reject on slot cap");
        };

        let q1 = ByteBoundedQueue::<Heavy>::new(1_000_000);
        assert!(!q1.slot_cap_warned.load(Ordering::Relaxed));
        fill_to_slot_cap(&q1);
        assert!(q1.slot_cap_warned.load(Ordering::Relaxed), "first queue must warn on its hit");

        // A fresh queue starts un-warned even though q1 already warned, so it
        // will warn on its own first hit (per-queue, not process-global).
        let q2 = ByteBoundedQueue::<Heavy>::new(1_000_000);
        assert!(
            !q2.slot_cap_warned.load(Ordering::Relaxed),
            "a second queue must NOT inherit the first queue's warned state"
        );
        fill_to_slot_cap(&q2);
        assert!(
            q2.slot_cap_warned.load(Ordering::Relaxed),
            "second queue must warn on its own hit"
        );
    }

    #[test]
    fn byte_bounded_respects_limit() {
        let q = ByteBoundedQueue::<Heavy>::new(100);
        // Empty queue accepts even an oversized item (legacy semantics:
        // gate on `cur < limit`, not `cur + size <= limit`). This is the
        // fix for the per-item-larger-than-limit deadlock.
        assert!(q.try_push(Heavy(vec![0; 200])).is_ok());
        assert_eq!(q.current_bytes(), 200);
        // Now `cur >= limit`, all subsequent pushes reject regardless
        // of size.
        let rejected = q.try_push(Heavy(vec![0; 1]));
        assert!(rejected.is_err(), "queue at/over budget should reject");
        assert_eq!(q.current_bytes(), 200);
        // After a pop drops `cur` below limit, pushes succeed again.
        let _ = q.try_pop().unwrap();
        assert_eq!(q.current_bytes(), 0);
        assert!(q.try_push(Heavy(vec![0; 50])).is_ok());
        assert_eq!(q.current_bytes(), 50);
    }

    #[test]
    fn byte_bounded_oversized_first_push_succeeds() {
        // Regression: previously a single push larger than `limit_bytes`
        // would always reject (`0 + size > limit`), deadlocking
        // producers that emit oversized batches (e.g. busy-locus
        // position-group batches). With the legacy `cur < limit`
        // semantics, the oversized push goes through.
        let q = ByteBoundedQueue::<Heavy>::new(100);
        assert!(q.try_push(Heavy(vec![0; 1024])).is_ok());
    }

    #[test]
    fn byte_bounded_decrements_on_pop() {
        let q = ByteBoundedQueue::<Heavy>::new(1000);
        q.try_push(Heavy(vec![0; 200])).unwrap();
        assert_eq!(q.current_bytes(), 200);
        let _ = q.try_pop().unwrap();
        assert_eq!(q.current_bytes(), 0);
    }

    #[test]
    fn unbounded_never_rejects() {
        let q = UnboundedQueue::<u32>::new();
        for i in 0..1024 {
            assert!(q.try_push(i).is_ok());
        }
    }

    // ── Blocking facet (L2.2) ────────────────────────────────────────────────

    use std::time::Duration;

    #[test]
    fn non_blocking_queue_has_no_notifier() {
        // The pool hot path must stay notifier-free: a plain `new` queue holds
        // `blocking: None`, so `try_push`/`try_pop` take the nullable-branch
        // fast path (no fence, no atomic, no lock).
        assert!(CountBoundedQueue::<u32>::new(4).blocking.is_none());
        assert!(ByteBoundedQueue::<Heavy>::new(100).blocking.is_none());
        assert!(CountBoundedQueue::<u32>::new_blocking(4).blocking.is_some());
        assert!(ByteBoundedQueue::<Heavy>::new_blocking(100).blocking.is_some());
    }

    #[test]
    fn blocking_pop_wakes_on_push() {
        let q: Arc<CountBoundedQueue<u32>> = Arc::new(CountBoundedQueue::new_blocking(4));
        let qc = Arc::clone(&q);
        let consumer = std::thread::spawn(move || qc.blocking_pop());
        // Give the consumer time to park on the empty queue.
        std::thread::sleep(Duration::from_millis(50));
        q.try_push(99).unwrap();
        assert_eq!(consumer.join().unwrap(), Some(99), "consumer wakes with the pushed item");
    }

    #[test]
    fn blocking_pop_returns_none_on_drain() {
        // The shutdown signal: a parked consumer wakes with `None` when the
        // producer marks the queue drained without pushing.
        let q: Arc<CountBoundedQueue<u32>> = Arc::new(CountBoundedQueue::new_blocking(4));
        let qc = Arc::clone(&q);
        let consumer = std::thread::spawn(move || qc.blocking_pop());
        std::thread::sleep(Duration::from_millis(50));
        q.mark_drained();
        assert_eq!(consumer.join().unwrap(), None, "drain wakes the parked consumer with None");
    }

    #[test]
    fn blocking_pop_drains_remaining_then_none() {
        // Items already queued at drain time are returned before `None`.
        let q: Arc<CountBoundedQueue<u32>> = Arc::new(CountBoundedQueue::new_blocking(4));
        q.try_push(1).unwrap();
        q.try_push(2).unwrap();
        q.mark_drained();
        assert_eq!(q.blocking_pop(), Some(1));
        assert_eq!(q.blocking_pop(), Some(2));
        assert_eq!(q.blocking_pop(), None);
    }

    #[test]
    fn blocking_push_wakes_on_pop() {
        // Output backpressure: a producer parked on a full queue wakes when a
        // consumer frees a slot.
        let q: Arc<CountBoundedQueue<u32>> = Arc::new(CountBoundedQueue::new_blocking(1));
        q.try_push(1).unwrap(); // queue now full (capacity 1)
        let qp = Arc::clone(&q);
        let producer = std::thread::spawn(move || qp.blocking_push(2));
        std::thread::sleep(Duration::from_millis(50));
        assert_eq!(q.try_pop(), Some(1), "consumer frees the slot");
        producer.join().unwrap();
        assert_eq!(q.try_pop(), Some(2), "producer's parked item landed");
    }

    #[test]
    fn detached_edge_records_park_time() {
        // Deep mode: a consumer that parks on an empty queue accrues
        // park_on_empty_ns ≈ the wait duration; the producer side stays 0.
        let q: Arc<CountBoundedQueue<u32>> = Arc::new(CountBoundedQueue::new_blocking(4));
        let qc = Arc::clone(&q);
        let consumer = std::thread::spawn(move || qc.blocking_pop());
        std::thread::sleep(Duration::from_millis(20));
        q.try_push(7).unwrap();
        assert_eq!(consumer.join().unwrap(), Some(7));
        let (empty_ns, full_ns) = q.blocking.as_ref().unwrap().park_times();
        assert!(empty_ns >= 10_000_000, "parked ~20ms on empty, got {empty_ns} ns");
        assert_eq!(full_ns, 0, "no producer parked on full");
    }

    #[test]
    fn blocking_facet_works_on_byte_bounded_queue() {
        let q: Arc<ByteBoundedQueue<Heavy>> = Arc::new(ByteBoundedQueue::new_blocking(1000));
        let qc = Arc::clone(&q);
        let consumer = std::thread::spawn(move || qc.blocking_pop().map(|h| h.0.len()));
        std::thread::sleep(Duration::from_millis(50));
        q.try_push(Heavy(vec![0; 42])).unwrap();
        assert_eq!(consumer.join().unwrap(), Some(42));
    }

    // ── Per-edge metrics (L2-instrumentation Task 2) ─────────────────────────

    #[test]
    fn instrumented_queue_counts_push_and_reject() {
        let m = EdgeMetrics::new();
        let q = CountBoundedQueue::<u32>::new_instrumented(1, Arc::clone(&m));
        assert!(q.try_push(1).is_ok());
        assert_eq!(q.try_push(2), Err(2)); // full (cap 1) → reject
        let s = m.snapshot();
        assert_eq!(s.pushed_items, 1, "one successful push");
        assert_eq!(s.push_rejections, 1, "one rejection");
        // Producer-push only at this layer; pop is counted at the input handle.
        assert_eq!(s.popped_items, 0);
    }

    #[test]
    fn byte_bounded_instrumented_push_bytes_and_depth() {
        let m = EdgeMetrics::new();
        let q = ByteBoundedQueue::<Heavy>::new_instrumented(1000, Arc::clone(&m));
        q.try_push(Heavy(vec![0; 200])).unwrap();
        let s = m.snapshot();
        assert_eq!(s.pushed_items, 1);
        assert_eq!(s.pushed_bytes, 200);
        assert!((q.depth_fraction() - 0.2).abs() < 0.01, "200/1000 ≈ 0.2");
        assert_eq!(q.occupancy_len(), 1);
    }

    #[test]
    fn byte_bounded_instrumented_counts_reject() {
        // The byte-budget reject path increments push_rejections (distinct from
        // the CountBounded slot reject above). Fill to budget, then a push rejects.
        let m = EdgeMetrics::new();
        let q = ByteBoundedQueue::<Heavy>::new_instrumented(100, Arc::clone(&m));
        q.try_push(Heavy(vec![0; 200])).unwrap(); // accepted (cur<limit on empty), now over budget
        assert!(q.try_push(Heavy(vec![0; 1])).is_err(), "at/over budget rejects");
        let s = m.snapshot();
        assert_eq!(s.pushed_items, 1);
        assert_eq!(s.push_rejections, 1, "byte-budget reject counted");
    }

    #[test]
    fn non_instrumented_queue_has_no_metrics() {
        // Hot-path guard: plain `new`/`new_blocking` queues are metric-free
        // (one nullable branch, no atomics). Mirrors non_blocking_queue_has_no_notifier.
        assert!(CountBoundedQueue::<u32>::new(4).metrics.is_none());
        assert!(CountBoundedQueue::<u32>::new_blocking(4).metrics.is_none());
        assert!(ByteBoundedQueue::<Heavy>::new(100).metrics.is_none());
        assert!(UnboundedQueue::<u32>::new().metrics.is_none());
        // And instrumented constructors do attach metrics.
        assert!(
            CountBoundedQueue::<u32>::new_instrumented(4, EdgeMetrics::new()).metrics.is_some()
        );
        assert!(UnboundedQueue::<u32>::new_instrumented(EdgeMetrics::new()).metrics.is_some());
    }

    /// Hammer the exact lost-wakeup window (M5): each round the producer pushes
    /// with no delay so its push frequently lands while the consumer is between
    /// its fast-path `try_pop` and committing to `wait`. The Dekker `SeqCst`
    /// fence pair must make every push observable or every park notified; a lost
    /// wakeup would hang one round, tripped by the watchdog. (A formal `loom`
    /// model of this interleaving is added at the lever-2 gate alongside the
    /// merge bridge.)
    #[test]
    fn blocking_pop_stress_no_lost_wakeup() {
        use std::sync::atomic::AtomicBool;
        let done = Arc::new(AtomicBool::new(false));
        // Watchdog: a lost wakeup parks a consumer forever; fail loudly instead.
        {
            let done = Arc::clone(&done);
            std::thread::spawn(move || {
                for _ in 0..200 {
                    std::thread::sleep(Duration::from_millis(50));
                    if done.load(Ordering::SeqCst) {
                        return;
                    }
                }
                // 10s elapsed without completion → a round is wedged.
                eprintln!("blocking_pop_stress_no_lost_wakeup: WEDGED (lost wakeup)");
                std::process::abort();
            });
        }
        for round in 0..5_000u32 {
            let q: Arc<CountBoundedQueue<u32>> = Arc::new(CountBoundedQueue::new_blocking(2));
            let qc = Arc::clone(&q);
            let consumer = std::thread::spawn(move || qc.blocking_pop());
            // No sleep: push races the consumer's park decision.
            q.try_push(round).unwrap();
            assert_eq!(consumer.join().unwrap(), Some(round), "round {round} lost its item");
        }
        done.store(true, Ordering::SeqCst);
    }
}
