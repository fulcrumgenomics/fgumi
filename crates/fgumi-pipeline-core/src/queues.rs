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
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

use super::item::HeapSize;
use super::runtime::metrics::EdgeMetrics;

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
        Self::build(capacity, None)
    }

    /// Like [`new`](Self::new) but recording producer-push metrics into `metrics`
    /// (an instrumented edge). The non-blocking `try_*` surface is unchanged.
    ///
    /// # Panics
    ///
    /// Panics if `capacity == 0`.
    #[must_use]
    pub fn new_instrumented(capacity: usize, metrics: Arc<EdgeMetrics>) -> Self {
        Self::build(capacity, Some(metrics))
    }

    /// [`new`](Self::new) when `metrics` is `None`, [`new_instrumented`](Self::new_instrumented)
    /// when `Some`. Lets branch builders thread an optional metrics handle uniformly.
    ///
    /// # Panics
    ///
    /// Panics if `capacity == 0`.
    #[must_use]
    pub fn maybe_instrumented(capacity: usize, metrics: Option<Arc<EdgeMetrics>>) -> Self {
        Self::build(capacity, metrics)
    }

    fn build(capacity: usize, metrics: Option<Arc<EdgeMetrics>>) -> Self {
        assert!(capacity > 0, "CountBoundedQueue capacity must be > 0");
        Self { inner: ArrayQueue::new(capacity), drained: AtomicBool::new(false), metrics }
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
        Ok(())
    }

    fn try_pop(&self) -> Option<T> {
        let item = self.inner.pop()?;
        Some(item)
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
        Self::build(limit_bytes, None)
    }

    /// Like [`new`](Self::new) but recording producer-push metrics (items + bytes
    /// + rejections) into `metrics`. Byte-budget semantics unchanged.
    ///
    /// # Panics
    ///
    /// Panics if `limit_bytes == 0`.
    #[must_use]
    pub fn new_instrumented(limit_bytes: u64, metrics: Arc<EdgeMetrics>) -> Self {
        Self::build(limit_bytes, Some(metrics))
    }

    /// [`new`](Self::new) when `metrics` is `None`, [`new_instrumented`](Self::new_instrumented)
    /// when `Some`.
    ///
    /// # Panics
    ///
    /// Panics if `limit_bytes == 0`.
    #[must_use]
    pub fn maybe_instrumented(limit_bytes: u64, metrics: Option<Arc<EdgeMetrics>>) -> Self {
        Self::build(limit_bytes, metrics)
    }

    fn build(limit_bytes: u64, metrics: Option<Arc<EdgeMetrics>>) -> Self {
        assert!(limit_bytes > 0, "ByteBoundedQueue limit_bytes must be > 0");
        Self {
            inner: ArrayQueue::new(BYTE_BOUNDED_QUEUE_SLOT_CAPACITY),
            current_bytes: AtomicU64::new(0),
            limit_bytes: AtomicU64::new(limit_bytes),
            drained: AtomicBool::new(false),
            slot_cap_warned: AtomicBool::new(false),
            metrics,
        }
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
        Some(item)
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
        // Hot-path guard: a plain `new` queue is metric-free (one nullable
        // branch, no atomics).
        assert!(CountBoundedQueue::<u32>::new(4).metrics.is_none());
        assert!(ByteBoundedQueue::<Heavy>::new(100).metrics.is_none());
        assert!(UnboundedQueue::<u32>::new().metrics.is_none());
        // And instrumented constructors do attach metrics.
        assert!(
            CountBoundedQueue::<u32>::new_instrumented(4, EdgeMetrics::new()).metrics.is_some()
        );
        assert!(UnboundedQueue::<u32>::new_instrumented(EdgeMetrics::new()).metrics.is_some());
    }
}
