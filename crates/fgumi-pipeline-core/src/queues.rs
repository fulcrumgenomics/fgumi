//! Transport-layer queue trait + three concrete impls.
//!
//! Concerns: pure transport (push/pop, drained signal). **Not** ordering —
//! see [`crate::reorder`] for the `ReorderStage<T>` operator that adds
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
//!     `try_push` calls panic — in every build, not just debug (a contract
//!     violation: producer pushed after declaring done, and the item would be
//!     silently lost). See `assert_not_drained`.
//!   - Consumer (input side) checks `is_drained() && is_empty()` to detect
//!     end-of-stream. Once both are true, no further items will arrive.

use crossbeam_queue::{ArrayQueue, SegQueue};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, OnceLock};

use super::item::HeapSize;
use super::runtime::metrics::EdgeMetrics;
use super::runtime::wake_slot::{HolderSet, note_held};

/// Transport-layer queue trait. Type-uniform across queue impls: the
/// `try_push` surface accepts any `T` regardless of whether the impl uses
/// item-count or memory bookkeeping internally.
///
/// `Send + Sync`: queues are shared between worker threads via `Arc`.
///
/// # Sealed
///
/// The trait cannot be implemented outside this crate: its wake-protocol
/// methods take a crate-private token. Use the queues this crate provides.
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
    /// A `try_push` after `mark_drained` panics — see `assert_not_drained`.
    fn mark_drained(&self);

    /// True if `mark_drained` has been called.
    fn is_drained(&self) -> bool;

    // The methods below are the protocol between a transport and a layer
    // stacked on it (the `ReorderStage`), which pushes without the transport
    // deciding on its own that a refusal is final. Each takes a [`Sealed`]
    // token, which only this crate can name or build, and none has a provided
    // body: a provided default either records nothing or latches where it must
    // not, so the compiler makes every transport state its own behaviour. The
    // token also seals the trait: code outside this crate cannot implement it.

    /// One push attempt with `try_push`'s success bookkeeping (push counter,
    /// hold-clock end, push metric) and none of its refusal bookkeeping: no
    /// holder record, no hold stamp, no reject count. The layer above decides
    /// whether a refusal is final (the reorder stage's must-accept arm turns it
    /// into a stash insert, which is an accepted push).
    ///
    /// # Errors
    ///
    /// Returns `Err(item)` when the transport is at its limit.
    #[doc(hidden)]
    fn try_push_unlatched(&self, item: T, _: Sealed) -> Result<(), T>;

    /// After a `try_push_unlatched` refusal the caller would hold. On a tracked
    /// edge: record the calling thread, fence, and return whether the transport
    /// now has room or `also()` holds, in which case the caller retries. On an
    /// untracked edge: `false` at once, nothing recorded. `also` is the layer's
    /// own release condition, which the transport cannot see (the reorder
    /// stage's `next_serial` no longer buffered).
    #[doc(hidden)]
    fn latch_holder_and_recheck(&self, also: &dyn Fn() -> bool, _: Sealed) -> bool;

    /// The caller is returning that refusal: mark the thread as holding if the
    /// edge is tracked, time the hold if it has a hold clock, count the reject
    /// if it is instrumented.
    #[doc(hidden)]
    fn on_final_reject(&self, _: Sealed);

    /// A push refused by [`Self::try_push_unlatched`] was accepted by the layer
    /// above instead (the reorder stash): if the calling producer was holding
    /// an item on this edge, its hold ends now, as on a successful push.
    #[doc(hidden)]
    fn note_stash_accepted(&self, _: Sealed);

    /// Record the calling thread as holding an item that a cap layered above
    /// this transport refused (the reorder stash cap). On a tracked edge the
    /// holder is recorded with no fence (the layer's lock orders it against the
    /// consumer's in-order pop, which takes the same lock first) and the thread
    /// is marked as holding; with a hold clock the hold is timed.
    #[doc(hidden)]
    fn note_rejected_holder(&self, _: Sealed);
}

/// The token that seals the transport-protocol methods of [`ItemQueue`] and
/// [`BoundedQueueHandle`]: a type in a private module whose only value is
/// [`SEALED`], so code outside this crate can neither call nor override a
/// method that takes one.
mod sealed {
    #[derive(Clone, Copy, Debug)]
    pub struct Sealed(pub(super) ());
}
pub(crate) use sealed::Sealed;

/// The one [`Sealed`] value.
pub(crate) const SEALED: Sealed = Sealed(());

/// Per-edge bookkeeping that `Pipeline::run` installs once on the edges it
/// tracks. Absent (`OnceLock` empty) on every other edge, so an untracked push
/// pays one predictable branch.
pub struct EdgeTracking {
    /// Hold timing for `--pipeline-stats` (`None` when stats are off).
    pub(crate) clock: Option<crate::runtime::wake_slot::HoldClock>,
    /// The threads holding an item this edge refused, and the push counting
    /// behind the wake plan's push gate (`None` on edges no wake plan routes).
    pub(crate) holders: Option<crate::runtime::wake_slot::HolderSet>,
}

/// Panic if `try_push` is called after `mark_drained`.
///
/// The consumer treats a drained queue as closed, so an item pushed afterwards
/// may never be popped: the item is silently lost and the loss surfaces (if at
/// all) as a short output far from its cause. That makes this a framework
/// contract violation rather than a recoverable condition, so — like
/// `BranchOutputHandle::retry`'s `Ordered` + `ordinal = None` arm — it fails
/// loudly in **every** build. It was `debug_assert!`-only, which left release
/// builds performing exactly the silent push the message warns about.
///
/// `Relaxed` is sufficient here and is the cheaper load on a per-item path.
/// `drained` is monotonic — its only write anywhere is `store(true, Release)` in
/// `mark_drained` — so a `Relaxed` load can return a stale `false` (a missed
/// detection when the producer races the close on another thread) but never a
/// spurious `true`. It cannot panic a correct program.
#[inline]
fn assert_not_drained(drained: &AtomicBool, queue_kind: &'static str) {
    assert!(
        !drained.load(Ordering::Relaxed),
        "{queue_kind}::try_push after mark_drained — producer contract violation"
    );
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
///
/// Tracking: in a Directed chain, `Pipeline::run` installs a holder set
/// (`HolderOnlyHandle::enable_holder_tracking`) when the wake plan
/// reverse-wakes this edge. A refused push then records the refusing thread,
/// fences and re-checks for room, and the consumer's pop-side `Progress` takes
/// and unparks it. There is no hold clock and no push counter: a count-bounded
/// edge is never push-gated (see `HolderOnlyHandle`).
pub struct CountBoundedQueue<T: Send + 'static> {
    inner: ArrayQueue<T>,
    drained: AtomicBool,
    /// `Some` only on an instrumented edge (`--pipeline-trace`); `None` keeps the
    /// hot path metric-free. Producer-push counts are recorded here; consumer-pop
    /// counts are recorded at the `BranchInputHandle` (see `handles.rs`).
    metrics: Option<Arc<EdgeMetrics>>,
    /// The threads this edge refused, installed by `Pipeline::run` on an edge
    /// the wake plan reverse-wakes; empty otherwise (Legacy, same-thread, fused).
    holders: OnceLock<HolderSet>,
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
        Self {
            inner: ArrayQueue::new(capacity),
            drained: AtomicBool::new(false),
            metrics,
            holders: OnceLock::new(),
        }
    }

    /// The one push site behind `try_push` (`latch`) and `try_push_unlatched`:
    /// one attempt with the success-path bookkeeping, and on refusal either
    /// `push_refused` (`latch`) or `Err(item)`. A single caller of the
    /// `ArrayQueue` push lets the compiler inline that push here, and both
    /// trait methods reduce to a tail call; with a copy per entry point (or a
    /// const-generic body inlined into each) the push stays out of line, a
    /// call per push (measured on `pipeline_dispatch`).
    #[inline(never)]
    fn push_with(&self, item: T, latch: bool) -> Result<(), T> {
        assert_not_drained(&self.drained, "CountBoundedQueue");
        match self.inner.push(item) {
            Ok(()) => {
                if let Some(m) = &self.metrics {
                    m.record_push(0); // count-bounded: items only, no byte size
                }
                Ok(())
            }
            Err(item) if latch => self.push_refused(item),
            Err(item) => Err(item),
        }
    }

    /// `try_push`'s refusal path: latch → fence → re-check, retrying while
    /// room appears, else book the final reject. Out of line and cold so the
    /// success path stays a few instructions: inlined, the latch's thread-local
    /// slot lookup is hoisted to the top of the push and paid on every push.
    /// A retry happens only when a pop landed between the refusal and the
    /// latch; each one that loses again means another producer pushed, so the
    /// loop is lock-free.
    #[cold]
    #[inline(never)]
    fn push_refused(&self, mut item: T) -> Result<(), T> {
        loop {
            if !self.latch_holder_and_recheck(&|| false, SEALED) {
                self.on_final_reject(SEALED);
                return Err(item);
            }
            match self.push_with(item, false) {
                Ok(()) => return Ok(()),
                Err(back) => item = back,
            }
        }
    }
}

impl<T: Send + 'static> ItemQueue<T> for CountBoundedQueue<T> {
    fn try_push(&self, item: T) -> Result<(), T> {
        self.push_with(item, true)
    }

    fn try_push_unlatched(&self, item: T, _: Sealed) -> Result<(), T> {
        self.push_with(item, false)
    }

    fn latch_holder_and_recheck(&self, also: &dyn Fn() -> bool, _: Sealed) -> bool {
        // A pop frees a slot with a SeqCst CAS on the ArrayQueue's `head`, and
        // `is_full` loads `tail` then `head` (SeqCst, crossbeam-queue 0.3.13), so
        // after the latch's fence the re-check sees any pop whose consumer fence
        // came first.
        self.holders.get().is_some_and(|h| h.latch_and_recheck(|| !self.inner.is_full() || also()))
    }

    fn on_final_reject(&self, _: Sealed) {
        if self.holders.get().is_some() {
            note_held();
        }
        if let Some(m) = &self.metrics {
            m.record_reject();
        }
    }

    /// No hold clock on a count-bounded edge: nothing to end.
    fn note_stash_accepted(&self, _: Sealed) {}

    fn note_rejected_holder(&self, _: Sealed) {
        if let Some(h) = self.holders.get() {
            h.record_current();
            note_held();
        }
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
    /// Hold timing and holder set installed by `Pipeline::run` on the edges
    /// it tracks; empty otherwise.
    tracking: OnceLock<EdgeTracking>,
    /// Successful pushes on an edge with a holder set, bumped `Relaxed` on the
    /// success path. Read by `WakePlan` as the push gate: a Detached producer's
    /// `Progress` wakes its consumer only if this advanced during the dispatch.
    /// Same struct as `current_bytes`, which that path already RMWs. Not
    /// counted on untracked edges.
    pushed_total: AtomicU64,
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
            tracking: OnceLock::new(),
            pushed_total: AtomicU64::new(0),
        }
    }

    /// The one push site behind `try_push` (`latch`) and `try_push_unlatched`:
    /// one attempt with the success-path bookkeeping (push counter on a
    /// tracked edge, hold-clock end, push metric), and on refusal either
    /// `push_refused` (`latch`) or `Err(item)`. A push is refused at either
    /// point: the byte budget is reached, or the fixed slot array is full (the
    /// byte reservation is rolled back and the slot-cap warning fires once). A
    /// single caller of the `ArrayQueue` push lets the compiler inline that
    /// push here, and both trait methods reduce to a tail call; with a copy per
    /// entry point (or a const-generic body inlined into each) the push stays
    /// out of line, a call per push (measured on `pipeline_dispatch`).
    #[inline(never)]
    fn push_with(&self, item: T, latch: bool) -> Result<(), T> {
        assert_not_drained(&self.drained, "ByteBoundedQueue");
        // Like the legacy `ReorderBufferState::can_proceed`, this gates on
        // `heap_bytes < limit` — accept if currently *under* budget, regardless
        // of incoming item size. Per-item-larger-than-limit is a real case
        // (busy-locus position-group batches can be tens of MB while the queue
        // limit is 4 MiB), so a strict `cur + size <= limit` would deadlock the
        // producer.
        //
        // Once `cur` reaches `limit_bytes`, subsequent pushes reject until a
        // consumer drains. Transient overshoot under concurrent pushes is
        // self-correcting on the next round.
        let refused = 'push: {
            let cur = self.current_bytes.load(Ordering::Relaxed);
            let limit = self.limit_bytes.load(Ordering::Relaxed);
            if cur >= limit {
                break 'push item;
            }
            let size = item.heap_size() as u64;
            // Reserve bytes before pushing so a concurrent consumer cannot pop
            // and decrement the counter before we add our share, which would
            // cause the counter to underflow and create permanent false
            // backpressure.
            self.current_bytes.fetch_add(size, Ordering::Relaxed);
            // ArrayQueue::push returns Err((item, size)) on full; roll back the
            // reservation and return the item to the caller for retry. (In
            // practice the slot cap should never be hit before the byte budget
            // triggers a reject above, but defend against it anyway.)
            match self.inner.push((item, size)) {
                Ok(()) => {
                    if let Some(t) = self.tracking.get() {
                        if t.holders.is_some() {
                            self.pushed_total.fetch_add(1, Ordering::Relaxed);
                        }
                        if let Some(c) = &t.clock {
                            c.on_push();
                        }
                    }
                    if let Some(m) = &self.metrics {
                        m.record_push(size);
                    }
                    return Ok(());
                }
                Err((rejected, _size)) => {
                    // Roll back the byte reservation — the item never entered
                    // the queue.
                    self.current_bytes.fetch_sub(size, Ordering::Relaxed);
                    self.warn_slot_cap_once();
                    rejected
                }
            }
        };
        if latch { self.push_refused(refused) } else { Err(refused) }
    }

    /// The fixed 1024-slot backing was hit before the byte budget. This
    /// degrades byte-backpressure into a hard count cap for small items
    /// (`heap_size` ≲ limit/1024) — correctness is preserved (the producer
    /// retries) but throughput silently suffers. Surface it once so it is
    /// observable rather than a silent foot-gun; near-zero cost after the first
    /// hit. Cold and out of line, so the logging code stays off the push path.
    #[cold]
    #[inline(never)]
    fn warn_slot_cap_once(&self) {
        if !self.slot_cap_warned.swap(true, Ordering::Relaxed) {
            log::warn!(
                "ByteBoundedQueue hit its {BYTE_BOUNDED_QUEUE_SLOT_CAPACITY}-slot count cap \
                 before the byte budget; small items are degrading byte-backpressure into a \
                 count cap (throughput, not correctness, is affected)."
            );
        }
    }

    /// `try_push`'s refusal path: latch → fence → re-check, retrying while
    /// room appears, else book the final reject. Out of line and cold so the
    /// success path stays short: inlined, the latch's thread-local slot lookup
    /// is hoisted to the top of the push and paid on every push. A retry
    /// happens only when a pop or a limit raise landed between the refusal and
    /// the latch; each one that loses again means another producer pushed, so
    /// the loop is lock-free.
    #[cold]
    #[inline(never)]
    fn push_refused(&self, mut item: T) -> Result<(), T> {
        loop {
            if !self.latch_holder_and_recheck(&|| false, SEALED) {
                self.on_final_reject(SEALED);
                return Err(item);
            }
            match self.push_with(item, false) {
                Ok(()) => return Ok(()),
                Err(back) => item = back,
            }
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
        // Floor at 1. `try_push` rejects when `current_bytes >= limit_bytes`, so a
        // limit of 0 rejects unconditionally — even on an empty edge — and wedges
        // the producer permanently. `new` asserts `limit_bytes > 0` for exactly
        // this reason; without a floor here that invariant could be undone after
        // construction, which is the one case the constructor cannot guard.
        //
        // Clamped rather than asserted: this runs on a live pipeline (the budget
        // pass and the rebalancer), where degrading to a 1-byte limit still makes
        // progress — `try_push` admits an item whenever `current_bytes` is under
        // the limit, regardless of item size — while a panic would take down a
        // running pipeline over a recoverable arithmetic slip. Every current caller
        // already applies its own positive per-queue floor.
        self.limit_bytes.store(new_limit.max(1), Ordering::Relaxed);
    }
}

/// An edge whose holders the wake plan may track but that has no push counter
/// and no byte budget: a [`CountBoundedQueue`], or the [`UnboundedQueue`]
/// behind a `ReorderStage` (whose stash cap can refuse). With no push counter
/// such an edge cannot be push-gated, by type: a gate built on a counter that
/// does not count would drop every forward wake. Crate-internal: only
/// `Pipeline::run` and the wake plan use it.
pub(crate) trait HolderOnlyHandle: Send + Sync {
    /// Install a holder set sized for `n_slots` pipeline threads (once; a second
    /// call is ignored). Takes a count rather than a set so the sizing rule
    /// (`n_threads + n_detached`) stays in `Pipeline::run`.
    fn enable_holder_tracking(&self, n_slots: usize);
    /// Take (and clear) the recorded holders, calling `f(slot)` per holder;
    /// returns how many. 0 on an edge without a holder set.
    fn take_holders(&self, f: &mut dyn FnMut(usize)) -> usize;
    /// Whether a holder set is installed on this edge.
    fn is_tracked(&self) -> bool;
}

/// The consumer-side view of an edge that records the threads it refused: the
/// wake plan's reverse edge takes them after the consumer's pop and unparks
/// each one. One variant per registry an edge can sit in.
#[derive(Clone)]
pub(crate) enum HolderQueue {
    /// A byte-bounded edge (which may also be push-gated).
    Bytes(Arc<dyn BoundedQueueHandle>),
    /// A count-bounded edge, or the unbounded transport behind a reorder stage.
    HolderOnly(Arc<dyn HolderOnlyHandle>),
}

impl HolderQueue {
    /// Take (and clear) the recorded holders, calling `f(slot)` per holder;
    /// returns how many. 0 on an edge without a holder set.
    #[inline]
    pub(crate) fn take_holders(&self, f: &mut dyn FnMut(usize)) -> usize {
        match self {
            Self::Bytes(q) => q.take_holders(f, SEALED),
            Self::HolderOnly(q) => q.take_holders(f),
        }
    }

    /// Whether a holder set is installed on this edge.
    pub(crate) fn is_tracked(&self) -> bool {
        match self {
            Self::Bytes(q) => q.is_tracked(SEALED),
            Self::HolderOnly(q) => q.is_tracked(),
        }
    }
}

/// Type-erased handle for a byte-bounded queue. The pipeline
/// rebalancer iterates over registered handles to read fullness
/// (`current_bytes / limit_bytes`) and reallocate budget across
/// queues by calling `set_limit_bytes`. The trait deliberately
/// does not surface the queue's item type or its `ItemQueue`
/// methods — rebalancing only needs the byte counters.
///
/// # Sealed
///
/// The trait cannot be implemented outside this crate: its wake-protocol
/// methods take a crate-private token (see [`ItemQueue`]), so outside this
/// crate they can be neither called nor implemented.
pub trait BoundedQueueHandle: Send + Sync {
    /// Bytes currently held in the queue.
    fn current_bytes(&self) -> u64;
    /// Current byte-budget cap. May change between calls if a
    /// rebalancer is active.
    fn limit_bytes(&self) -> u64;
    /// Update the byte-budget cap. Concurrent producers see the
    /// new value on their next push.
    fn set_limit_bytes(&self, new_limit: u64);
    /// Install this edge's tracking (once; a second call is ignored).
    #[doc(hidden)]
    fn enable_tracking(&self, t: EdgeTracking, _: Sealed);
    /// Successful pushes so far (`Relaxed`). Advances only on edges with a
    /// holder set.
    #[doc(hidden)]
    fn pushed_total(&self, _: Sealed) -> u64;
    /// Take (and clear) the threads this edge refused, calling `f(slot)` per
    /// holder; returns how many. 0 on an edge without a holder set.
    #[doc(hidden)]
    fn take_holders(&self, f: &mut dyn FnMut(usize), _: Sealed) -> usize;
    /// Whether a holder set is installed on this edge.
    #[doc(hidden)]
    fn is_tracked(&self, _: Sealed) -> bool;
    /// How many hold stamps the edge's hold clock keeps (tests): `None` with
    /// no clock.
    #[cfg(test)]
    fn hold_clock_stamps(&self) -> Option<usize> {
        None
    }
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
    fn enable_tracking(&self, t: EdgeTracking, _: Sealed) {
        let _ = self.tracking.set(t);
    }
    fn pushed_total(&self, _: Sealed) -> u64 {
        self.pushed_total.load(Ordering::Relaxed)
    }
    fn take_holders(&self, f: &mut dyn FnMut(usize), _: Sealed) -> usize {
        self.tracking.get().and_then(|t| t.holders.as_ref()).map_or(0, |h| h.take(f))
    }
    fn is_tracked(&self, _: Sealed) -> bool {
        self.tracking.get().is_some_and(|t| t.holders.is_some())
    }
    #[cfg(test)]
    fn hold_clock_stamps(&self) -> Option<usize> {
        self.tracking
            .get()
            .and_then(|t| t.clock.as_ref())
            .map(crate::runtime::wake_slot::HoldClock::stamps)
    }
}

impl<T: Send + HeapSize + 'static> ItemQueue<T> for ByteBoundedQueue<T> {
    fn try_push(&self, item: T) -> Result<(), T> {
        self.push_with(item, true)
    }

    fn try_push_unlatched(&self, item: T, _: Sealed) -> Result<(), T> {
        self.push_with(item, false)
    }

    fn latch_holder_and_recheck(&self, also: &dyn Fn() -> bool, _: Sealed) -> bool {
        let Some(holders) = self.tracking.get().and_then(|t| t.holders.as_ref()) else {
            return false;
        };
        // Both loads `Relaxed`: after the latch's fence they see any pop or limit
        // raise whose fence came first. The parentheses matter: without them a
        // reading of `a && (b || c)` would drop the layer's term whenever the
        // edge is over budget.
        holders.latch_and_recheck(|| {
            (self.current_bytes.load(Ordering::Relaxed) < self.limit_bytes.load(Ordering::Relaxed)
                && !self.inner.is_full())
                || also()
        })
    }

    fn on_final_reject(&self, _: Sealed) {
        if let Some(t) = self.tracking.get() {
            if t.holders.is_some() {
                note_held();
            }
            if let Some(c) = &t.clock {
                c.on_reject();
            }
        }
        if let Some(m) = &self.metrics {
            m.record_reject();
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

    fn note_stash_accepted(&self, _: Sealed) {
        if let Some(t) = self.tracking.get()
            && let Some(c) = &t.clock
        {
            c.on_push();
        }
    }

    /// No fence and no re-check, unlike the transport's own reject path: the
    /// reorder cap runs under `ReorderStage`'s state lock, and the consumer's
    /// in-order pop takes the same lock before its holder take, so lock order
    /// already decides whether the take sees this latch or the producer's
    /// retry sees the room.
    fn note_rejected_holder(&self, _: Sealed) {
        if let Some(t) = self.tracking.get() {
            if let Some(h) = &t.holders {
                h.record_current();
                note_held();
            }
            if let Some(c) = &t.clock {
                c.on_reject();
            }
        }
    }
}

impl<T: Send + 'static> HolderOnlyHandle for CountBoundedQueue<T> {
    fn enable_holder_tracking(&self, n_slots: usize) {
        let _ = self.holders.set(HolderSet::new(n_slots));
    }
    fn take_holders(&self, f: &mut dyn FnMut(usize)) -> usize {
        self.holders.get().map_or(0, |h| h.take(f))
    }
    fn is_tracked(&self) -> bool {
        self.holders.get().is_some()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// UnboundedQueue
// ─────────────────────────────────────────────────────────────────────────────

/// Unbounded transport. `try_push` always succeeds. Backed by `SegQueue<T>`.
///
/// Tracking: the transport itself never refuses, but an ordered branch puts a
/// `ReorderStage` with a stash cap over it, and that cap does refuse. On such a
/// branch `Pipeline::run` installs a holder set when the wake plan
/// reverse-wakes the edge, and the cap's refusal records the holder through
/// [`ItemQueue::note_rejected_holder`]. A direct unbounded edge is never
/// tracked.
pub struct UnboundedQueue<T: Send + 'static> {
    inner: SegQueue<T>,
    drained: AtomicBool,
    /// `Some` only on an instrumented edge; producer-push items are recorded
    /// here (unbounded → never rejects, no byte tracking). Consumer-pop is at the
    /// `BranchInputHandle`.
    metrics: Option<Arc<EdgeMetrics>>,
    /// The threads the reorder stash cap above this transport refused; used
    /// only by `note_rejected_holder`. Empty unless the wake plan tracks the
    /// branch.
    holders: OnceLock<HolderSet>,
}

impl<T: Send + 'static> UnboundedQueue<T> {
    #[must_use]
    pub fn new() -> Self {
        Self::maybe_instrumented(None)
    }

    /// Like [`new`](Self::new) but recording producer-push item counts into
    /// `metrics`. Unbounded edges have no byte budget and never reject; depth is
    /// reported as raw length only.
    #[must_use]
    pub fn new_instrumented(metrics: Arc<EdgeMetrics>) -> Self {
        Self::maybe_instrumented(Some(metrics))
    }

    /// [`new`](Self::new) when `metrics` is `None`, [`new_instrumented`](Self::new_instrumented)
    /// when `Some`.
    #[must_use]
    pub fn maybe_instrumented(metrics: Option<Arc<EdgeMetrics>>) -> Self {
        Self {
            inner: SegQueue::new(),
            drained: AtomicBool::new(false),
            metrics,
            holders: OnceLock::new(),
        }
    }
}

impl<T: Send + 'static> Default for UnboundedQueue<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Send + 'static> ItemQueue<T> for UnboundedQueue<T> {
    fn try_push(&self, item: T) -> Result<(), T> {
        assert_not_drained(&self.drained, "UnboundedQueue");
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

    /// The transport never refuses, so this is `try_push`.
    fn try_push_unlatched(&self, item: T, _: Sealed) -> Result<(), T> {
        self.try_push(item)
    }

    /// Never reached (the transport never refuses); records nothing.
    fn latch_holder_and_recheck(&self, _also: &dyn Fn() -> bool, _: Sealed) -> bool {
        false
    }

    /// Never reached (the transport never refuses); nothing to book.
    fn on_final_reject(&self, _: Sealed) {}

    /// No hold clock on an unbounded edge: nothing to end.
    fn note_stash_accepted(&self, _: Sealed) {}

    /// The reorder stash cap above this transport refused: record the holder
    /// under the stage's lock and mark the thread, as on the bounded transports.
    fn note_rejected_holder(&self, _: Sealed) {
        if let Some(h) = self.holders.get() {
            h.record_current();
            note_held();
        }
    }
}

impl<T: Send + 'static> HolderOnlyHandle for UnboundedQueue<T> {
    fn enable_holder_tracking(&self, n_slots: usize) {
        let _ = self.holders.set(HolderSet::new(n_slots));
    }
    fn take_holders(&self, f: &mut dyn FnMut(usize)) -> usize {
        self.holders.get().map_or(0, |h| h.take(f))
    }
    fn is_tracked(&self) -> bool {
        self.holders.get().is_some()
    }
}

/// Test-only interleaving hook: lets a test run code on the producer thread
/// between a refusal and its holder latch.
#[cfg(test)]
pub(crate) mod test_hooks {
    use std::cell::RefCell;
    thread_local! {
        static BEFORE_LATCH: RefCell<Option<Box<dyn FnOnce()>>> = const { RefCell::new(None) };
    }
    /// Run `f` once, on this thread, between a refusal and its holder latch. It
    /// runs inside `HolderSet::latch_slot`, so it fires for every transport and
    /// for the reorder stage's fast-path latch. It may run under a
    /// `ReorderStage` lock (the stale-cache arm calls the latching `try_push`
    /// there, and that mutex is not reentrant): a hook must not touch that
    /// stage there.
    pub(crate) fn before_latch(f: impl FnOnce() + 'static) {
        BEFORE_LATCH.with(|h| *h.borrow_mut() = Some(Box::new(f)));
    }
    pub(crate) fn run_before_latch() {
        if let Some(f) = BEFORE_LATCH.with(|h| h.borrow_mut().take()) {
            f();
        }
    }
    thread_local! {
        static AFTER_RECORD: RefCell<Option<Box<dyn FnOnce()>>> = const { RefCell::new(None) };
    }
    /// Run `f` once, on this thread, between a latch's record (and fence) and
    /// its re-check: where a release that claims the fresh bit lands.
    pub(crate) fn after_record(f: impl FnOnce() + 'static) {
        AFTER_RECORD.with(|h| *h.borrow_mut() = Some(Box::new(f)));
    }
    pub(crate) fn run_after_record() {
        if let Some(f) = AFTER_RECORD.with(|h| h.borrow_mut().take()) {
            f();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::sync::Arc;

    /// A `try_push` after `mark_drained` must panic on every transport impl, in
    /// every build. It was `debug_assert!`-only, so a release build pushed the
    /// item into a queue the consumer had already closed — a silent loss that
    /// surfaces only as a short output. `#[values]` covers all three impls so a
    /// new transport that forgets the guard is caught by the same table.
    #[rstest]
    #[case::count_bounded(Arc::new(CountBoundedQueue::new(2)) as Arc<dyn ItemQueue<u32>>)]
    #[case::byte_bounded(Arc::new(ByteBoundedQueue::new(1024)) as Arc<dyn ItemQueue<u32>>)]
    #[case::unbounded(Arc::new(UnboundedQueue::new()) as Arc<dyn ItemQueue<u32>>)]
    #[should_panic(expected = "try_push after mark_drained — producer contract violation")]
    fn try_push_after_mark_drained_panics(#[case] q: Arc<dyn ItemQueue<u32>>) {
        q.mark_drained();
        let _ = q.try_push(1);
    }

    /// The guard must not fire before `mark_drained` — a plain push on a fresh
    /// queue still succeeds on every impl.
    #[rstest]
    #[case::count_bounded(Arc::new(CountBoundedQueue::new(2)) as Arc<dyn ItemQueue<u32>>)]
    #[case::byte_bounded(Arc::new(ByteBoundedQueue::new(1024)) as Arc<dyn ItemQueue<u32>>)]
    #[case::unbounded(Arc::new(UnboundedQueue::new()) as Arc<dyn ItemQueue<u32>>)]
    fn try_push_before_mark_drained_succeeds(#[case] q: Arc<dyn ItemQueue<u32>>) {
        assert!(q.try_push(1).is_ok(), "an undrained queue must still accept a push");
        assert_eq!(q.try_pop(), Some(1));
    }

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

    /// The slot-cap reject path reserves `size` bytes *before* pushing and rolls
    /// the reservation back when `ArrayQueue::push` reports full. The sibling test
    /// above fills with 0-byte items, so that `fetch_sub` runs with `size == 0`
    /// and a leak or double-subtract is invisible. Fill with nonzero items under a
    /// limit large enough that the slot cap still binds, and pin the byte counter
    /// across the rejected push.
    #[test]
    fn byte_bounded_slot_cap_reject_rolls_back_reserved_bytes() {
        const ITEM_BYTES: usize = 8;
        // Large enough that 1024 * 8 bytes never reaches it, so the reject below
        // is the slot cap and not the byte budget.
        let q = ByteBoundedQueue::<Heavy>::new(1_000_000);
        for i in 0..BYTE_BOUNDED_QUEUE_SLOT_CAPACITY {
            assert!(q.try_push(Heavy(vec![0; ITEM_BYTES])).is_ok(), "push {i} within slot cap");
        }
        let before = q.current_bytes();
        assert_eq!(
            before,
            (BYTE_BOUNDED_QUEUE_SLOT_CAPACITY * ITEM_BYTES) as u64,
            "every admitted item's bytes are accounted"
        );
        assert!(before < 1_000_000, "the byte budget must not be the binding cap here");
        assert!(
            q.try_push(Heavy(vec![0; ITEM_BYTES])).is_err(),
            "push #{} must reject on the slot cap",
            BYTE_BOUNDED_QUEUE_SLOT_CAPACITY + 1
        );
        assert_eq!(
            q.current_bytes(),
            before,
            "a slot-cap reject must roll its reservation back exactly — leaking bytes here \
             would create permanent false backpressure"
        );
    }

    /// A 0 limit makes `try_push` reject unconditionally (`current_bytes >= 0` is
    /// always true), wedging the producer forever. `new` asserts against it, so
    /// the setter must not be able to reintroduce it after construction. Clamping
    /// to 1 keeps the edge alive: `try_push` admits an item whenever
    /// `current_bytes` is *under* the limit, whatever the item's size.
    #[test]
    fn set_limit_bytes_clamps_zero_to_one_so_the_edge_still_admits() {
        let q = ByteBoundedQueue::<Heavy>::new(4096);
        q.set_limit_bytes(0);
        assert_eq!(q.limit_bytes(), 1, "a 0 limit is floored to 1, never stored as 0");
        assert!(
            q.try_push(Heavy(vec![0; 64])).is_ok(),
            "an empty edge must still admit one item — a 0 limit would reject forever"
        );
        // Now over the 1-byte limit, so the next push rejects: still a real bound,
        // not a silent promotion to unbounded.
        assert!(q.try_push(Heavy(vec![0; 64])).is_err(), "the clamped limit still applies");
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
        assert!(
            ByteBoundedQueue::<Heavy>::new_instrumented(100, EdgeMetrics::new()).metrics.is_some()
        );

        // `maybe_instrumented` is the constructor the branch builders actually
        // call, and it was the only one with no coverage: one that dropped a
        // `Some(metrics)` would silently produce an edge that reports nothing
        // under `--pipeline-trace`, with every other test still green. Both
        // directions, all three impls.
        assert!(CountBoundedQueue::<u32>::maybe_instrumented(4, None).metrics.is_none());
        assert!(ByteBoundedQueue::<Heavy>::maybe_instrumented(100, None).metrics.is_none());
        assert!(UnboundedQueue::<u32>::maybe_instrumented(None).metrics.is_none());
        assert!(
            CountBoundedQueue::<u32>::maybe_instrumented(4, Some(EdgeMetrics::new()))
                .metrics
                .is_some()
        );
        assert!(
            ByteBoundedQueue::<Heavy>::maybe_instrumented(100, Some(EdgeMetrics::new()))
                .metrics
                .is_some()
        );
        assert!(
            UnboundedQueue::<u32>::maybe_instrumented(Some(EdgeMetrics::new())).metrics.is_some()
        );
    }

    /// 8-byte items so a 1-byte limit admits exactly one.
    #[derive(Debug, PartialEq)]
    struct Sized8;
    impl HeapSize for Sized8 {
        fn heap_size(&self) -> usize {
            8
        }
    }

    /// A hold is timed from the producer thread's first rejection to its next
    /// successful push on the same edge; every rejection in between is a retry.
    /// The interleaving is forced with channels, and the wait is bounded below
    /// by a sleep that happens strictly inside the hold.
    #[test]
    fn hold_clock_times_first_reject_to_next_push() {
        use std::sync::mpsc::channel;
        use std::time::Duration;

        use crate::runtime::stats::PipelineStats;
        use crate::runtime::wake_slot::{HoldClock, SlotGuard};
        use crate::topology::StepIdx;
        let stats = Arc::new(PipelineStats::new(vec!["P", "C"]));
        let q = Arc::new(ByteBoundedQueue::<Sized8>::new(1));
        q.enable_tracking(
            EdgeTracking {
                clock: Some(HoldClock::per_thread(StepIdx(0), Arc::clone(&stats), 2)),
                holders: None,
            },
            SEALED,
        );
        let (rejected_tx, rejected_rx) = channel::<()>();
        let (popped_tx, popped_rx) = channel::<()>();
        let qp = Arc::clone(&q);
        let producer = std::thread::spawn(move || {
            let _slot = SlotGuard::enter(1);
            qp.try_push(Sized8).unwrap();
            assert!(qp.try_push(Sized8).is_err(), "full → first rejection stamps the hold");
            assert!(qp.try_push(Sized8).is_err(), "second rejection is a retry");
            rejected_tx.send(()).unwrap();
            popped_rx.recv().unwrap();
            qp.try_push(Sized8).unwrap(); // releases the hold
        });
        rejected_rx.recv().unwrap();
        std::thread::sleep(Duration::from_millis(5)); // strictly inside the hold
        assert_eq!(q.try_pop(), Some(Sized8));
        popped_tx.send(()).unwrap();
        producer.join().unwrap();
        let s = stats.snapshot().steps[0].1;
        assert_eq!((s.holds, s.held_retries), (1, 1));
        assert!(s.held_wait_ns >= 5_000_000, "hold spans the 5 ms sleep: {}", s.held_wait_ns);
        assert_eq!(stats.snapshot().steps[1].1.holds, 0, "the consumer step is untouched");
    }

    /// Without tracking (every edge when stats are off) the reject path is today's.
    #[test]
    fn untracked_edge_records_nothing_and_rejects_as_before() {
        let q = ByteBoundedQueue::<Sized8>::new(1);
        q.try_push(Sized8).unwrap();
        assert!(q.try_push(Sized8).is_err());
        assert!(q.tracking.get().is_none());
    }

    /// A non-`Parallel` producer's held item may be retried by another worker
    /// (a `Serial` step runs on whichever worker wins its lock): with the
    /// per-step clock, worker 2's flush of the item worker 1 was refused ends
    /// the hold, and worker 1's next unrelated push books nothing more. The
    /// per-thread clock, a `Parallel` producer's, keeps the two apart.
    #[rstest]
    #[case::per_step(true, (1, 0))]
    #[case::per_thread(false, (0, 0))]
    fn a_flush_by_another_worker_ends_a_per_step_hold(
        #[case] per_step: bool,
        #[case] after_flush: (u64, u64),
    ) {
        use std::time::Duration;

        use crate::runtime::stats::PipelineStats;
        use crate::runtime::wake_slot::{HoldClock, SlotGuard};
        use crate::topology::StepIdx;
        fn on_worker(
            q: &Arc<ByteBoundedQueue<Sized8>>,
            slot: usize,
            f: impl FnOnce(&ByteBoundedQueue<Sized8>) + Send + 'static,
        ) {
            let q = Arc::clone(q);
            std::thread::spawn(move || {
                let _slot = SlotGuard::enter(slot);
                f(&q);
            })
            .join()
            .unwrap();
        }
        let stats = Arc::new(PipelineStats::new(vec!["P", "C"]));
        let q = Arc::new(ByteBoundedQueue::<Sized8>::new(1));
        let clock = if per_step {
            HoldClock::per_step(StepIdx(0), Arc::clone(&stats))
        } else {
            HoldClock::per_thread(StepIdx(0), Arc::clone(&stats), 3)
        };
        q.enable_tracking(EdgeTracking { clock: Some(clock), holders: None }, SEALED);
        q.try_push(Sized8).unwrap();
        on_worker(&q, 1, |q| assert!(q.try_push(Sized8).is_err(), "worker 1 is refused"));
        assert_eq!(q.try_pop(), Some(Sized8));
        on_worker(&q, 2, |q| q.try_push(Sized8).unwrap()); // worker 2 flushes the item
        let s = stats.snapshot().steps[0].1;
        assert_eq!((s.holds, s.held_retries), after_flush);
        assert_eq!(q.try_pop(), Some(Sized8));
        std::thread::sleep(Duration::from_millis(500));
        on_worker(&q, 1, |q| q.try_push(Sized8).unwrap()); // unrelated push
        let s = stats.snapshot().steps[0].1;
        if per_step {
            assert_eq!(s.holds, 1, "the hold was booked once, at the flush");
            assert!(s.held_wait_ns < 250_000_000, "not stretched to the unrelated push");
        } else {
            assert_eq!(s.holds, 1, "per thread: worker 1's own next push ends its hold");
        }
    }

    /// A tracked byte-bounded edge still under its byte budget refuses once its
    /// fixed slot array is full. The re-check after the latch must count that
    /// as no room (`!inner.is_full()`), or the refusal loop retries a push that
    /// can never land, forever. Run under a watchdog so that spin fails the test
    /// rather than hanging it; the refusal must also record and mark its holder.
    #[test]
    fn a_slot_cap_refusal_on_a_tracked_edge_under_budget_is_a_hold() {
        use std::sync::mpsc::channel;
        use std::time::Duration;

        use crate::runtime::wake_slot::{SlotGuard, take_rejected};
        let q = tracked(1 << 30, 4);
        for i in 0..BYTE_BOUNDED_QUEUE_SLOT_CAPACITY {
            assert!(q.try_push(Sized8).is_ok(), "push {i} within the slot cap");
        }
        let (tx, rx) = channel();
        let qp = Arc::clone(&q);
        std::thread::spawn(move || {
            let _slot = SlotGuard::enter(2);
            let _ = take_rejected();
            let refused = qp.try_push(Sized8).is_err();
            let _ = tx.send((refused, take_rejected()));
        });
        let (refused, marked) = rx
            .recv_timeout(Duration::from_secs(10))
            .expect("the refusal must return, not spin on a retry that cannot land");
        assert!(refused, "the slot cap refuses although the byte budget has room");
        assert!(marked, "...and the refusing thread is marked as holding");
        assert_eq!(holders_of(&q), vec![2], "...and recorded");
    }

    fn tracked(limit: u64, n_slots: usize) -> Arc<ByteBoundedQueue<Sized8>> {
        use crate::runtime::wake_slot::HolderSet;
        let q = Arc::new(ByteBoundedQueue::<Sized8>::new(limit));
        q.enable_tracking(
            EdgeTracking { clock: None, holders: Some(HolderSet::new(n_slots)) },
            SEALED,
        );
        q
    }

    fn holders_of(q: &ByteBoundedQueue<Sized8>) -> Vec<usize> {
        let mut v = Vec::new();
        q.take_holders(&mut |s| v.push(s), SEALED);
        v
    }

    /// `pushed_total` advances once per successful push on a tracked edge and
    /// never on a reject; a rejection records the rejecting thread's slot, and
    /// one take returns and clears it.
    #[test]
    fn tracked_edge_counts_pushes_and_records_the_holder_slot() {
        use crate::runtime::wake_slot::SlotGuard;
        let q = tracked(1, 4);
        let _slot = SlotGuard::enter(2);
        assert!(q.try_push(Sized8).is_ok());
        assert_eq!(q.pushed_total(SEALED), 1);
        assert!(q.try_push(Sized8).is_err(), "cur (8) >= limit (1) rejects");
        assert_eq!(q.pushed_total(SEALED), 1, "a reject does not count as a push");
        assert_eq!(holders_of(&q), vec![2], "the rejection names slot 2");
        assert!(holders_of(&q).is_empty(), "...and the take cleared it");
        assert!(crate::runtime::wake_slot::take_rejected(), "the thread knows it was rejected");
        assert!(!crate::runtime::wake_slot::take_rejected(), "...once");
    }

    /// Untracked (the default): no counting, no holders — today's behaviour.
    #[test]
    fn untracked_edge_counts_nothing() {
        use crate::runtime::wake_slot::SlotGuard;
        let q = ByteBoundedQueue::<Sized8>::new(1);
        let _slot = SlotGuard::enter(0);
        q.try_push(Sized8).unwrap();
        assert!(q.try_push(Sized8).is_err());
        assert_eq!(q.pushed_total(SEALED), 0);
        assert_eq!(holders_of(&q), Vec::<usize>::new());
        assert!(!q.is_tracked(SEALED));
        assert!(!crate::runtime::wake_slot::take_rejected(), "an untracked reject latches nothing");
    }

    /// The two transports that can refuse a direct push.
    #[derive(Clone, Copy, Debug)]
    enum Kind {
        Byte,
        Count,
    }

    /// A transport that admits exactly one `Sized8`, seen as a queue and as a
    /// holder edge; `bytes` is the concrete byte-bounded queue, for a limit
    /// raise.
    struct OneItem {
        q: Arc<dyn ItemQueue<Sized8>>,
        edge: HolderQueue,
        bytes: Option<Arc<ByteBoundedQueue<Sized8>>>,
    }

    impl OneItem {
        /// Tracked with a `HolderSet::new(n)` when `slots` is `Some(n)`;
        /// untracked otherwise.
        fn new(kind: Kind, slots: Option<usize>) -> Self {
            use crate::runtime::wake_slot::HolderSet;
            match kind {
                Kind::Byte => {
                    let b = Arc::new(ByteBoundedQueue::<Sized8>::new(1));
                    if let Some(n) = slots {
                        b.enable_tracking(
                            EdgeTracking { clock: None, holders: Some(HolderSet::new(n)) },
                            SEALED,
                        );
                    }
                    Self {
                        q: Arc::clone(&b) as _,
                        edge: HolderQueue::Bytes(Arc::clone(&b) as _),
                        bytes: Some(b),
                    }
                }
                Kind::Count => {
                    let c = Arc::new(CountBoundedQueue::<Sized8>::new(1));
                    if let Some(n) = slots {
                        c.enable_holder_tracking(n);
                    }
                    Self { q: Arc::clone(&c) as _, edge: HolderQueue::HolderOnly(c), bytes: None }
                }
            }
        }

        fn take(&self) -> Vec<usize> {
            let mut v = Vec::new();
            self.edge.take_holders(&mut |s| v.push(s));
            v
        }
    }

    /// A final refusal on a tracked edge records the refusing thread's slot,
    /// and marks the thread as holding, exactly once; one take clears it.
    #[rstest]
    #[case::byte_bounded(Kind::Byte)]
    #[case::count_bounded(Kind::Count)]
    fn tracked_reject_records_the_holder_and_marks_the_thread(#[case] kind: Kind) {
        use crate::runtime::wake_slot::{SlotGuard, take_rejected};
        let _ = take_rejected();
        let e = OneItem::new(kind, Some(4));
        let _slot = SlotGuard::enter(2);
        e.q.try_push(Sized8).unwrap();
        assert!(e.q.try_push(Sized8).is_err(), "full");
        assert_eq!(e.take(), vec![2], "the refusal names slot 2");
        assert!(e.take().is_empty(), "...and the take cleared it");
        assert!(take_rejected(), "the thread knows it is holding");
        assert!(!take_rejected(), "...once");
    }

    /// Untracked (Legacy, same-thread, fused): a refusal records nothing and
    /// marks nothing.
    #[rstest]
    #[case::byte_bounded(Kind::Byte)]
    #[case::count_bounded(Kind::Count)]
    fn untracked_reject_records_nothing(#[case] kind: Kind) {
        use crate::runtime::wake_slot::{SlotGuard, take_rejected};
        let _ = take_rejected();
        let e = OneItem::new(kind, None);
        let _slot = SlotGuard::enter(0);
        e.q.try_push(Sized8).unwrap();
        assert!(e.q.try_push(Sized8).is_err());
        assert_eq!(e.edge.take_holders(&mut |_| {}), 0);
        assert!(!e.edge.is_tracked());
        assert!(!take_rejected(), "an untracked refusal marks nothing");
    }

    /// A thread with no wake slot is recorded in the anonymous slot.
    #[rstest]
    #[case::byte_bounded(Kind::Byte)]
    #[case::count_bounded(Kind::Count)]
    fn slotless_thread_is_recorded_as_anonymous(#[case] kind: Kind) {
        let e = OneItem::new(kind, Some(3));
        e.q.try_push(Sized8).unwrap();
        assert!(e.q.try_push(Sized8).is_err());
        assert_eq!(e.take(), vec![3]);
    }

    /// How room appears between the refusal and the latch.
    #[derive(Clone, Copy, Debug)]
    enum Release {
        Pop,
        LimitRaise,
    }

    /// The lost-wake interleaving: the consumer's pop (and its holder take), or
    /// a limit raise, lands between the producer's refusal and its latch. With
    /// the latch → fence → re-check protocol the producer sees the room and
    /// retries instead of holding, and the retry leaves the thread unmarked.
    /// The hook runs on the producer thread at exactly that point, so the
    /// interleaving is forced, not hoped for.
    #[rstest]
    #[case::room_from_a_pop(Kind::Byte, Release::Pop)]
    #[case::room_from_a_limit_raise(Kind::Byte, Release::LimitRaise)]
    #[case::count_bounded_pop(Kind::Count, Release::Pop)]
    fn room_between_check_and_latch_is_retried_not_held(
        #[case] kind: Kind,
        #[case] release: Release,
    ) {
        use crate::runtime::wake_slot::{SlotGuard, take_rejected};
        let _ = take_rejected();
        let e = OneItem::new(kind, Some(2));
        let _slot = SlotGuard::enter(1);
        e.q.try_push(Sized8).unwrap(); // full
        let (q2, edge2, bytes2) = (Arc::clone(&e.q), e.edge.clone(), e.bytes.clone());
        let taken = Arc::new(std::sync::atomic::AtomicUsize::new(usize::MAX));
        let taken2 = Arc::clone(&taken);
        test_hooks::before_latch(move || {
            match release {
                Release::Pop => assert_eq!(q2.try_pop(), Some(Sized8)),
                Release::LimitRaise => {
                    bytes2
                        .expect("a limit raise needs the byte transport")
                        .set_limit_bytes(1 << 20);
                }
            }
            let n = edge2.take_holders(&mut |_| {});
            taken2.store(n, std::sync::atomic::Ordering::SeqCst);
        });
        assert!(
            e.q.try_push(Sized8).is_ok(),
            "room appeared before the latch: the push must be retried"
        );
        assert_eq!(
            taken.load(std::sync::atomic::Ordering::SeqCst),
            0,
            "the take ran before the latch and saw no holder"
        );
        assert_eq!(
            e.take(),
            vec![1],
            "the latch stays set until a consumer-side take (one spurious wake)"
        );
        assert!(!take_rejected(), "a retried refusal is not a hold");
    }

    /// The stage protocol's unlatched push neither records nor marks: a layer
    /// above decides whether the refusal is final.
    #[rstest]
    #[case::byte_bounded(Kind::Byte)]
    #[case::count_bounded(Kind::Count)]
    fn unlatched_refusal_records_nothing(#[case] kind: Kind) {
        use crate::runtime::wake_slot::{SlotGuard, take_rejected};
        let _ = take_rejected();
        let e = OneItem::new(kind, Some(2));
        let _slot = SlotGuard::enter(1);
        e.q.try_push_unlatched(Sized8, SEALED).unwrap();
        assert!(e.q.try_push_unlatched(Sized8, SEALED).is_err(), "full");
        assert_eq!(e.edge.take_holders(&mut |_| {}), 0, "nothing recorded");
        assert!(!take_rejected(), "nothing marked");
    }

    /// A Legacy byte-bounded edge under `--pipeline-stats` has a hold clock and
    /// no holder set: its refusal books the hold but must not mark the thread,
    /// or the Legacy idle branch would leave the event-count for a timer park
    /// that nothing unparks.
    #[test]
    fn clock_only_reject_leaves_the_thread_unmarked() {
        use crate::runtime::stats::PipelineStats;
        use crate::runtime::wake_slot::{HoldClock, SlotGuard, take_rejected};
        use crate::topology::StepIdx;
        let _ = take_rejected();
        let stats = Arc::new(PipelineStats::new(vec!["P", "C"]));
        let q = ByteBoundedQueue::<Sized8>::new(1);
        q.enable_tracking(
            EdgeTracking {
                clock: Some(HoldClock::per_thread(StepIdx(0), Arc::clone(&stats), 2)),
                holders: None,
            },
            SEALED,
        );
        let _slot = SlotGuard::enter(1);
        q.try_push(Sized8).unwrap();
        assert!(q.try_push(Sized8).is_err());
        assert_eq!(q.try_pop(), Some(Sized8));
        q.try_push(Sized8).unwrap(); // ends the hold
        assert_eq!(stats.snapshot().steps[0].1.holds, 1, "the clock booked the hold");
        assert!(!take_rejected(), "a clock-only refusal leaves the thread unmarked");
    }
}
