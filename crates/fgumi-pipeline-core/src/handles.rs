//! Framework-side input/output handles wired to the layered queue stack.
//!
//! Per-branch architecture:
//!
//! ```text
//!     producer step
//!         │ outputs.push(item)
//!         ▼
//!     BranchOutputHandle<T>
//!         │  - allocates ordinal (if BranchOrdering::ByOrdinal)
//!         │  - heap-size accounting happens inside ByteBoundedQueue
//!         │  - forwards to either ReorderStage or transport directly
//!         ▼
//!  ┌────────────────────────────────┐
//!  │  ReorderStage<T> (optional)    │  ← only if BranchOrdering::ByOrdinal
//!  └────────────────────────────────┘
//!         ▼
//!  ┌──────────────────────────────────────┐
//!  │  Arc<dyn ItemQueue<Wrapped>>         │  ← Wrapped = Sequenced<T> if ordered, else T
//!  └──────────────────────────────────────┘
//!         ▼
//!     BranchInputHandle<T>
//!         │ pop() -> Option<T>
//!         ▼
//!     consumer step
//! ```
//!
//! No code in this module surfaces serials or heap sizes to step authors —
//! those are framework-managed. Step authors see `outputs.push(item)` and
//! `input.pop() -> Option<T>`.
//!
//! `QueueSpec::ByteBounded` requires `T: HeapSize`. `Single<T>` (and the
//! ordered-bytes paths) support `ByteBounded` directly: `build_single_queues`
//! bounds `T: HeapSize` and routes `ByteBounded` through
//! `build_branch_byte_aware`. The only remaining panics are for tuple branches
//! built via `build_branch` — which has the `HeapSize` bound but is the
//! non-byte-aware path, so it never constructs a byte-bounded queue and panics
//! on `ByteBounded` (the byte-aware `build_branch_byte_aware` must be used
//! instead) — and for `ByItemOrdinal` declared without `Ordered`.

use std::any::Any;
use std::marker::PhantomData;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

use super::item::{HeapSize, Ordered};
use super::outputs::Single;
use super::queues::{ByteBoundedQueue, CountBoundedQueue, ItemQueue, QueueSpec, UnboundedQueue};
use super::reorder::{
    BranchOrdering, DEFAULT_REORDER_OVERFLOW_BYTES, ReorderCapHandle, ReorderStage, Sequenced,
};
use super::step::{InputHandle, OutputHandles, OutputsViewAny};

// ─────────────────────────────────────────────────────────────────────────────
// BranchOutputHandle<T> — one per output branch on the producer side.
// ─────────────────────────────────────────────────────────────────────────────

/// Producer-side handle for one output branch. Step authors call `push(item)`.
/// The framework manages ordinal assignment and routes to the right destination
/// (transport or reorder stage).
pub struct BranchOutputHandle<T: Send + HeapSize + 'static> {
    inner: BranchOutputInner<T>,
}

enum BranchOutputInner<T: Send + HeapSize + 'static> {
    /// `BranchOrdering::None`: push directly to transport.
    Direct(Arc<dyn ItemQueue<T>>),
    /// `BranchOrdering::ByOrdinal` or `ByItemOrdinal`: push through
    /// `ReorderStage`. The stage internally wraps in `Sequenced<T>` and
    /// forwards to its transport. `OrdinalSource` decides where the
    /// ordinal comes from on each push.
    Ordered { stage: Arc<ReorderStage<T>>, ordinal_source: OrdinalSource<T> },
}

/// Where each push's ordinal comes from on an `Ordered` branch.
enum OrdinalSource<T: Send + HeapSize + 'static> {
    /// `BranchOrdering::ByOrdinal`: producer allocates a fresh ordinal via
    /// this counter on every push.
    Allocated(Arc<AtomicU64>),
    /// `BranchOrdering::ByItemOrdinal`: items carry their own ordinal. The
    /// function pointer is `|item: &T| item.ordinal()`, monomorphized at
    /// branch construction time (requires `T: Ordered`).
    ItemSerial(fn(&T) -> u64),
}

impl<T: Send + HeapSize + 'static> OrdinalSource<T> {
    #[inline]
    fn next(&self, item: &T) -> u64 {
        match self {
            Self::Allocated(c) => c.fetch_add(1, AtomicOrdering::AcqRel),
            Self::ItemSerial(f) => f(item),
        }
    }
}

/// An item the underlying queue rejected due to backpressure.
///
/// `Unpushed<T>` carries the rejected item back to the producer for retry,
/// **plus an opaque token** that — for ordered (`ByOrdinal`) branches —
/// preserves the ordinal allocated on the original push. Retrying via
/// [`BranchOutputHandle::retry`] reuses that ordinal so the consumer's
/// `ReorderStage` doesn't stall waiting for a sequence number that was
/// burned on a rejected attempt.
///
/// Step authors store `HeldSlot<Unpushed<T>>` (not `HeldSlot<T>`) when their
/// step pushes to an ordered branch. Direct branches don't carry an ordinal,
/// so `Unpushed<T>` is just a typed wrapper around the item there.
pub struct Unpushed<T: Send + HeapSize + 'static> {
    item: T,
    /// Pre-allocated ordinal from a rejected `Ordered` push. `None` for
    /// `Direct` (FIFO) branches.
    ordinal: Option<u64>,
}

/// Outcome of [`OutputHandles::retry_held`] — a step's drain-time held-slot
/// retry. Lets the caller distinguish "nothing was held" from "flushed
/// successfully" (the two differ for a single-final-batch flusher: the former
/// must still build its batch, the latter is done).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HeldRetry {
    /// The held slot was empty — nothing to retry.
    WasEmpty,
    /// A held item was retried and accepted; the slot is now empty.
    Flushed,
    /// A held item was retried, still rejected by backpressure, and put back
    /// in the slot. The caller must yield (return `NoProgress`/`Contention`
    /// and retry on the next dispatch).
    StillHeld,
}

impl<T: Send + HeapSize + 'static> std::fmt::Debug for Unpushed<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Unpushed").field("ordinal", &self.ordinal).field("item", &"<...>").finish()
    }
}

impl<T: Send + HeapSize + 'static> Unpushed<T> {
    /// Borrow the rejected item (read-only).
    pub fn item(&self) -> &T {
        &self.item
    }

    /// Discard the framework token and recover the item, abandoning the
    /// pre-allocated ordinal (if any). After calling this, **do not** push
    /// the item back into an ordered branch — the missing ordinal will
    /// stall the consumer's `ReorderStage`. Use only when abandoning the
    /// item entirely (e.g., during cancellation).
    pub fn into_item(self) -> T {
        self.item
    }
}

impl<T: Send + HeapSize + 'static> BranchOutputHandle<T> {
    /// Try to push a fresh item.
    ///
    /// # Errors
    ///
    /// Returns `Err(Unpushed<T>)` when the underlying queue (or reorder
    /// stage, for `ByOrdinal` branches) applied backpressure. The producer
    /// must hand the `Unpushed<T>` to a subsequent [`Self::retry`] call —
    /// **not** call `push` with the recovered item — so any pre-allocated
    /// ordinal is preserved.
    pub fn push(&self, item: T) -> Result<(), Unpushed<T>> {
        match &self.inner {
            BranchOutputInner::Direct(q) => {
                q.try_push(item).map_err(|item| Unpushed { item, ordinal: None })
            }
            BranchOutputInner::Ordered { stage, ordinal_source } => {
                let ord = ordinal_source.next(&item);
                stage
                    .try_push(ord, item)
                    .map_err(|(ord, item)| Unpushed { item, ordinal: Some(ord) })
            }
        }
    }

    /// Retry a previously-rejected push. Reuses the original ordinal for
    /// ordered branches so the consumer's `ReorderStage` sees a contiguous
    /// sequence.
    ///
    /// # Errors
    ///
    /// Returns `Err(Unpushed<T>)` if the queue still applied backpressure;
    /// the caller should re-hold and retry on the next worker iteration.
    ///
    /// # Panics
    ///
    /// Panics in debug builds if the `Unpushed<T>` came from a different
    /// branch kind than this handle (a framework invariant violation).
    pub fn retry(&self, unpushed: Unpushed<T>) -> Result<(), Unpushed<T>> {
        let Unpushed { item, ordinal } = unpushed;
        match (&self.inner, ordinal) {
            (BranchOutputInner::Direct(q), None) => {
                q.try_push(item).map_err(|item| Unpushed { item, ordinal: None })
            }
            (BranchOutputInner::Ordered { stage, .. }, Some(ord)) => stage
                .try_push(ord, item)
                .map_err(|(ord, item)| Unpushed { item, ordinal: Some(ord) }),
            (BranchOutputInner::Direct(q), Some(_)) => {
                debug_assert!(false, "Unpushed::ordinal=Some on a Direct branch");
                q.try_push(item).map_err(|item| Unpushed { item, ordinal: None })
            }
            (BranchOutputInner::Ordered { .. }, None) => {
                debug_assert!(false, "Unpushed::ordinal=None on an Ordered branch");
                self.push(item)
            }
        }
    }

    /// Mark this branch drained (producer-side). Idempotent.
    pub fn mark_drained(&self) {
        match &self.inner {
            BranchOutputInner::Direct(q) => q.mark_drained(),
            BranchOutputInner::Ordered { stage, .. } => stage.mark_drained(),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BranchInputHandle<T> — implements InputHandle<T>.
// ─────────────────────────────────────────────────────────────────────────────

/// Consumer-side input handle for one branch. Implements `InputHandle<T>`.
pub struct BranchInputHandle<T: Send + HeapSize + 'static> {
    inner: BranchInputInner<T>,
    /// `Some` only on an instrumented edge. The **consumer-pop** side of the
    /// edge's `EdgeMetrics` is recorded here (not at the transport queue),
    /// because an ordered edge pops through a `ReorderStage` that drains the
    /// transport in bulk — `pop` here is the real consumer pop. Shares the same
    /// `Arc<EdgeMetrics>` as the producer transport (push side).
    metrics: Option<Arc<crate::runtime::metrics::EdgeMetrics>>,
}

enum BranchInputInner<T: Send + HeapSize + 'static> {
    Direct(Arc<dyn ItemQueue<T>>),
    Ordered(Arc<ReorderStage<T>>),
    /// A permanently-empty, permanently-drained handle that owns no transport.
    /// Used for source steps' dummy input (their input is implicitly drained
    /// from t=0; the worker loop never pops from it), avoiding a per-source
    /// `SegQueue` allocation just to report `is_drained() == true`.
    AlwaysDrained(PhantomData<fn() -> T>),
}

impl<T: Send + HeapSize + 'static> BranchInputHandle<T> {
    /// Construct a zero-state input handle that is always empty and always
    /// drained, owning no backing queue. Used for source steps, whose input is
    /// implicitly drained from the start and never popped.
    #[must_use]
    pub fn always_drained() -> Self {
        Self { inner: BranchInputInner::AlwaysDrained(PhantomData), metrics: None }
    }
}

impl<T: Send + HeapSize + 'static> InputHandle<T> for BranchInputHandle<T> {
    fn pop(&self) -> Option<T> {
        // A `None` from an `Ordered` edge is *reorder-blocked* (not starved) while
        // the reorder stage still holds out-of-order items waiting for an earlier
        // ordinal — counting that as an empty pop inflates `pop_empties` and can
        // misclassify a backlogged edge as starved. The `reorder_blocked` flag is
        // reported by the SAME locked pop, so shared (`Parallel`) consumers can't
        // observe a torn `(item, blocked)` pair. Direct / always-drained branches
        // have no reorder buffer, so a `None` there is always a true empty pop.
        let (item, reorder_blocked) = match &self.inner {
            BranchInputInner::Direct(q) => (q.try_pop(), false),
            BranchInputInner::Ordered(stage) => stage.try_pop_in_order_reporting_blocked(),
            BranchInputInner::AlwaysDrained(_) => (None, false),
        };
        if let Some(m) = &self.metrics {
            if let Some(it) = &item {
                // usize→u64 is lossless on every target we build for.
                #[allow(clippy::cast_possible_truncation)]
                m.record_pop(it.heap_size() as u64);
            } else if !reorder_blocked {
                m.record_empty();
            }
        }
        item
    }

    fn is_drained(&self) -> bool {
        match &self.inner {
            BranchInputInner::Direct(q) => q.is_drained() && q.is_empty(),
            BranchInputInner::Ordered(stage) => stage.is_drained(),
            BranchInputInner::AlwaysDrained(_) => true,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Branch construction (one branch at a time, dispatching on QueueSpec + Ordering).
// ─────────────────────────────────────────────────────────────────────────────

/// The budget-resize handles for one byte-bounded branch: the transport
/// queue's limit setter, plus (for an ordered branch) the reorder stage's
/// overflow-cap setter. The framework's budget pass
/// (`apply_initial_queue_budget`) sets both — the transport limit and the
/// reorder cap — from the same per-edge budget, so they stay in lockstep
/// (single source of truth for `per_queue`).
pub(crate) struct BranchBudgetHandles {
    pub(crate) transport: Arc<dyn crate::queues::BoundedQueueHandle>,
    /// `Some` for an ordered (`ByOrdinal` / `ByItemOrdinal`) byte-bounded
    /// branch — its reorder overflow stash is sized from the same budget.
    /// `None` for a direct (unordered) byte-bounded branch (no reorder stage).
    pub(crate) reorder_cap: Option<Arc<dyn ReorderCapHandle>>,
}

/// One end-to-end branch: an output handle and an input handle wired
/// to the same underlying transport (and optional reorder stage).
///
/// `bounded_queue_handle` is `Some` iff the branch's transport is a
/// `ByteBoundedQueue`; it bundles the transport-limit and reorder-cap
/// setters (see [`BranchBudgetHandles`]). The pipeline-builder collects
/// these (across every branch in the chain) into a registry that the
/// budget pass + optional queue-memory rebalancer use to set/redistribute
/// budget. `None` for `CountBounded` / `Unbounded` branches.
pub(crate) struct Branch<T: Send + HeapSize + 'static> {
    pub(crate) output: BranchOutputHandle<T>,
    pub(crate) input: BranchInputHandle<T>,
    pub(crate) bounded_queue_handle: Option<BranchBudgetHandles>,
    /// `Some` on an instrumented edge — the shared `EdgeMetrics` (also held by
    /// the transport queue for push counts and the input handle for pop counts).
    /// Collected into the edge registry by `contexts.rs` Pass 1.5.
    pub(crate) metrics: Option<Arc<crate::runtime::metrics::EdgeMetrics>>,
}

/// Build one branch from a queue spec + ordering directive, where `T` does
/// not need to impl `HeapSize` or `Ordered`.
///
/// # Panics
///
/// Panics on `QueueSpec::ByteBounded` (requires `T: HeapSize`; use the
/// byte-aware build path) or on `BranchOrdering::ByItemOrdinal` (requires
/// `T: Ordered`; use `build_branch_ordered`).
/// Mint a per-edge [`EdgeMetrics`] when instrumentation is on, else `None`.
/// Called once per branch in the `build_branch*` constructors; the same handle
/// is shared by the transport (push counts) and the input handle (pop counts).
fn edge_metrics(
    level: crate::builder::InstrumentationLevel,
) -> Option<Arc<crate::runtime::metrics::EdgeMetrics>> {
    level.is_on().then(crate::runtime::metrics::EdgeMetrics::new)
}

pub(crate) fn build_branch<T: Send + HeapSize + 'static>(
    spec: QueueSpec,
    ordering: BranchOrdering,
    level: crate::builder::InstrumentationLevel,
) -> Branch<T> {
    match (spec, ordering) {
        (QueueSpec::CountBounded { capacity }, BranchOrdering::None) => {
            let m = edge_metrics(level);
            let q: Arc<dyn ItemQueue<T>> =
                Arc::new(CountBoundedQueue::<T>::maybe_instrumented(capacity, m.clone()));
            direct_branch(q, None, m)
        }
        (QueueSpec::CountBounded { capacity }, BranchOrdering::ByOrdinal) => {
            let m = edge_metrics(level);
            let transport: Arc<dyn ItemQueue<Sequenced<T>>> =
                Arc::new(CountBoundedQueue::<Sequenced<T>>::maybe_instrumented(
                    // Ordered transport is NOT instrumented: push/reject are recorded
                    // at the ReorderStage boundary (a stash turns a full-transport
                    // `Err` into an accepted `Ok`). Pop side is on the input handle.
                    capacity, None,
                ));
            ordered_branch(
                transport,
                OrdinalSource::Allocated(Arc::new(AtomicU64::new(0))),
                Some(DEFAULT_REORDER_OVERFLOW_BYTES),
                None,
                m,
            )
        }
        (QueueSpec::Unbounded, BranchOrdering::None) => {
            let m = edge_metrics(level);
            let q: Arc<dyn ItemQueue<T>> =
                Arc::new(UnboundedQueue::<T>::maybe_instrumented(m.clone()));
            direct_branch(q, None, m)
        }
        (QueueSpec::Unbounded, BranchOrdering::ByOrdinal) => {
            let m = edge_metrics(level);
            let transport: Arc<dyn ItemQueue<Sequenced<T>>> =
                // Ordered transport is NOT instrumented — push/reject recorded at
                // the ReorderStage boundary (see the count-bounded note above).
                Arc::new(UnboundedQueue::<Sequenced<T>>::maybe_instrumented(None));
            ordered_branch(
                transport,
                OrdinalSource::Allocated(Arc::new(AtomicU64::new(0))),
                Some(DEFAULT_REORDER_OVERFLOW_BYTES),
                None,
                m,
            )
        }
        (_, BranchOrdering::ByItemOrdinal) => {
            panic!(
                "BranchOrdering::ByItemOrdinal requires `T: Ordered` — \
                 use `build_branch_ordered::<T: Ordered>` instead. The plain \
                 `build_branch::<T>` path doesn't have the trait bound to \
                 read `item.ordinal()` at push time."
            );
        }
        (QueueSpec::ByteBounded { .. }, _) => {
            panic!(
                "QueueSpec::ByteBounded requires `T: HeapSize` — \
                 use the byte-aware build path (`build_branch_byte_aware` or \
                 `build_branch_ordered_bytes`)."
            );
        }
    }
}

/// Build one branch where `T: Ordered` (no `HeapSize` requirement).
/// Supports `BranchOrdering::ByItemOrdinal` (uses `item.ordinal()` for the
/// reorder stage's serial). Falls through to `build_branch::<T>` for non-
/// `ByItemOrdinal` orderings.
///
/// # Panics
///
/// Panics on `QueueSpec::ByteBounded` (requires `T: HeapSize`; use the
/// `_ordered_bytes` build path).
pub(crate) fn build_branch_ordered<T: Send + HeapSize + Ordered + 'static>(
    spec: QueueSpec,
    ordering: BranchOrdering,
    level: crate::builder::InstrumentationLevel,
) -> Branch<T> {
    match (spec, ordering) {
        (QueueSpec::CountBounded { capacity }, BranchOrdering::ByItemOrdinal) => {
            let m = edge_metrics(level);
            let transport: Arc<dyn ItemQueue<Sequenced<T>>> =
                Arc::new(CountBoundedQueue::<Sequenced<T>>::maybe_instrumented(
                    // Ordered transport is NOT instrumented: push/reject are recorded
                    // at the ReorderStage boundary (a stash turns a full-transport
                    // `Err` into an accepted `Ok`). Pop side is on the input handle.
                    capacity, None,
                ));
            ordered_branch(
                transport,
                OrdinalSource::ItemSerial(|item: &T| item.ordinal()),
                Some(DEFAULT_REORDER_OVERFLOW_BYTES),
                None,
                m,
            )
        }
        (QueueSpec::Unbounded, BranchOrdering::ByItemOrdinal) => {
            let m = edge_metrics(level);
            let transport: Arc<dyn ItemQueue<Sequenced<T>>> =
                // Ordered transport is NOT instrumented — push/reject recorded at
                // the ReorderStage boundary (see the count-bounded note above).
                Arc::new(UnboundedQueue::<Sequenced<T>>::maybe_instrumented(None));
            ordered_branch(
                transport,
                OrdinalSource::ItemSerial(|item: &T| item.ordinal()),
                Some(DEFAULT_REORDER_OVERFLOW_BYTES),
                None,
                m,
            )
        }
        (QueueSpec::ByteBounded { .. }, _) => {
            panic!(
                "QueueSpec::ByteBounded requires `T: HeapSize` — \
                 `build_branch_ordered` only bounds `T: Ordered`. Use \
                 `build_branch_ordered_bytes::<T: HeapSize + Ordered>` for the \
                 combined case."
            );
        }
        (other, ord) => build_branch::<T>(other, ord, level),
    }
}

/// Build one branch where `T: HeapSize + Ordered` (the canonical BAM step
/// case). Supports all `QueueSpec` × `BranchOrdering` combinations.
pub(crate) fn build_branch_ordered_bytes<T: Send + HeapSize + Ordered + 'static>(
    spec: QueueSpec,
    ordering: BranchOrdering,
    level: crate::builder::InstrumentationLevel,
) -> Branch<T> {
    use crate::queues::BoundedQueueHandle;

    match (spec, ordering) {
        (QueueSpec::ByteBounded { limit_bytes }, BranchOrdering::None) => {
            let m = edge_metrics(level);
            let q = Arc::new(ByteBoundedQueue::<T>::maybe_instrumented(limit_bytes, m.clone()));
            let handle: Arc<dyn BoundedQueueHandle> = Arc::clone(&q) as Arc<dyn BoundedQueueHandle>;
            let q_dyn: Arc<dyn ItemQueue<T>> = q;
            direct_branch(q_dyn, Some(handle), m)
        }
        (QueueSpec::ByteBounded { limit_bytes }, BranchOrdering::ByOrdinal) => {
            let m = edge_metrics(level);
            let transport_concrete = Arc::new(
                // Ordered transport is NOT instrumented — see the note on the
                // count-bounded ordered transport above; push/reject are recorded
                // at the ReorderStage boundary. The byte `handle` (occupancy /
                // budget resize) is derived from this same queue and is unaffected.
                ByteBoundedQueue::<Sequenced<T>>::maybe_instrumented(limit_bytes, None),
            );
            let handle: Arc<dyn BoundedQueueHandle> =
                Arc::clone(&transport_concrete) as Arc<dyn BoundedQueueHandle>;
            let transport: Arc<dyn ItemQueue<Sequenced<T>>> = transport_concrete;
            // Cap the must-accept overflow buffer so one worker grinding on a
            // large ordinal can't let every later ordinal overflow unbounded
            // (#29). This `DEFAULT_REORDER_OVERFLOW_BYTES` is the no-budget
            // FALLBACK: when `queue_memory_total` is set (production chains),
            // `apply_initial_queue_budget` re-sizes this cap thread-awarely
            // from the per-edge transport budget via the registered
            // `ReorderCapHandle` — so at low thread counts the stash is small
            // and at high thread counts it keeps this 256 MiB ceiling.
            ordered_branch(
                transport,
                OrdinalSource::Allocated(Arc::new(AtomicU64::new(0))),
                Some(DEFAULT_REORDER_OVERFLOW_BYTES),
                Some(handle),
                m,
            )
        }
        (QueueSpec::ByteBounded { limit_bytes }, BranchOrdering::ByItemOrdinal) => {
            let m = edge_metrics(level);
            let transport_concrete = Arc::new(
                // Ordered transport is NOT instrumented — see the note on the
                // count-bounded ordered transport above; push/reject are recorded
                // at the ReorderStage boundary. The byte `handle` (occupancy /
                // budget resize) is derived from this same queue and is unaffected.
                ByteBoundedQueue::<Sequenced<T>>::maybe_instrumented(limit_bytes, None),
            );
            let handle: Arc<dyn BoundedQueueHandle> =
                Arc::clone(&transport_concrete) as Arc<dyn BoundedQueueHandle>;
            let transport: Arc<dyn ItemQueue<Sequenced<T>>> = transport_concrete;
            ordered_branch(
                transport,
                OrdinalSource::ItemSerial(|item: &T| item.ordinal()),
                Some(DEFAULT_REORDER_OVERFLOW_BYTES),
                Some(handle),
                m,
            )
        }
        (other, ord) => build_branch_ordered::<T>(other, ord, level),
    }
}

/// Build one branch where `T: HeapSize`. Supports all three queue specs
/// for non-`ByItemOrdinal` orderings. `ByItemOrdinal` requires `T: Ordered`
/// — use `build_branch_ordered_bytes` for the combined case.
///
/// Currently exercised only by tests; production code uses
/// `build_branch_ordered_bytes` (BAM steps need both `HeapSize` and
/// `Ordered`). Kept for the rare case where a step needs byte-bounded
/// outputs without ordering.
#[allow(dead_code)]
pub(crate) fn build_branch_byte_aware<T: Send + HeapSize + 'static>(
    spec: QueueSpec,
    ordering: BranchOrdering,
    level: crate::builder::InstrumentationLevel,
) -> Branch<T> {
    use crate::queues::BoundedQueueHandle;

    match (spec, ordering) {
        (QueueSpec::ByteBounded { limit_bytes }, BranchOrdering::None) => {
            let m = edge_metrics(level);
            let q = Arc::new(ByteBoundedQueue::<T>::maybe_instrumented(limit_bytes, m.clone()));
            let handle: Arc<dyn BoundedQueueHandle> = Arc::clone(&q) as Arc<dyn BoundedQueueHandle>;
            let q_dyn: Arc<dyn ItemQueue<T>> = q;
            direct_branch(q_dyn, Some(handle), m)
        }
        (QueueSpec::ByteBounded { limit_bytes }, BranchOrdering::ByOrdinal) => {
            let m = edge_metrics(level);
            let transport_concrete = Arc::new(
                // Ordered transport is NOT instrumented — see the note on the
                // count-bounded ordered transport above; push/reject are recorded
                // at the ReorderStage boundary. The byte `handle` (occupancy /
                // budget resize) is derived from this same queue and is unaffected.
                ByteBoundedQueue::<Sequenced<T>>::maybe_instrumented(limit_bytes, None),
            );
            let handle: Arc<dyn BoundedQueueHandle> =
                Arc::clone(&transport_concrete) as Arc<dyn BoundedQueueHandle>;
            let transport: Arc<dyn ItemQueue<Sequenced<T>>> = transport_concrete;
            ordered_branch(
                transport,
                OrdinalSource::Allocated(Arc::new(AtomicU64::new(0))),
                Some(DEFAULT_REORDER_OVERFLOW_BYTES),
                Some(handle),
                m,
            )
        }
        (QueueSpec::ByteBounded { .. }, BranchOrdering::ByItemOrdinal) => {
            panic!(
                "ByteBounded + ByItemOrdinal requires `T: HeapSize + Ordered` — \
                 use `build_branch_ordered_bytes` instead."
            );
        }
        (other, ord) => build_branch::<T>(other, ord, level),
    }
}

fn direct_branch<T: Send + HeapSize + 'static>(
    q: Arc<dyn ItemQueue<T>>,
    transport_handle: Option<Arc<dyn crate::queues::BoundedQueueHandle>>,
    metrics: Option<Arc<crate::runtime::metrics::EdgeMetrics>>,
) -> Branch<T> {
    // A direct (unordered) branch has no reorder stage, so no reorder cap.
    let bounded_queue_handle =
        transport_handle.map(|transport| BranchBudgetHandles { transport, reorder_cap: None });
    Branch {
        output: BranchOutputHandle { inner: BranchOutputInner::Direct(Arc::clone(&q)) },
        input: BranchInputHandle { inner: BranchInputInner::Direct(q), metrics: metrics.clone() },
        bounded_queue_handle,
        metrics,
    }
}

fn ordered_branch<T: Send + HeapSize + 'static>(
    transport: Arc<dyn ItemQueue<Sequenced<T>>>,
    ordinal_source: OrdinalSource<T>,
    max_overflow_bytes: Option<u64>,
    transport_handle: Option<Arc<dyn crate::queues::BoundedQueueHandle>>,
    metrics: Option<Arc<crate::runtime::metrics::EdgeMetrics>>,
) -> Branch<T> {
    // A byte-bounded ordered edge carries a `transport_handle`; count/unbounded
    // ones don't. Record the push side's byte size only on the byte-bounded ones
    // (matching the queue's own byte accounting). Push/reject are recorded at the
    // ReorderStage boundary — NOT on the transport — because a must-accept stash
    // turns a full-transport `Err` into an accepted `Ok` (see `push_metrics`).
    let record_item_bytes = transport_handle.is_some();
    let stage = Arc::new(
        match max_overflow_bytes {
            Some(cap) => ReorderStage::<T>::with_max_overflow_bytes(transport, cap),
            None => ReorderStage::<T>::new(transport),
        }
        .with_push_metrics(metrics.clone(), record_item_bytes),
    );
    // Pair the transport-limit setter with this stage's overflow-cap setter so
    // the budget pass sizes both from one per-edge budget. Cloned (and coerced
    // to the trait object) before the stage is moved into the input handle.
    let reorder_cap: Arc<dyn ReorderCapHandle> = stage.clone();
    let bounded_queue_handle = transport_handle
        .map(|transport| BranchBudgetHandles { transport, reorder_cap: Some(reorder_cap) });
    Branch {
        output: BranchOutputHandle {
            inner: BranchOutputInner::Ordered { stage: Arc::clone(&stage), ordinal_source },
        },
        input: BranchInputHandle {
            inner: BranchInputInner::Ordered(stage),
            metrics: metrics.clone(),
        },
        bounded_queue_handle,
        metrics,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OutputQueueSet — per-step collection of typed input handles
// ─────────────────────────────────────────────────────────────────────────────

/// Type-erased per-branch input handles, owned by the producer step's chain
/// entry until the consumer claims them. The framework moves one entry from
/// here into each consumer's `TypedStep<S>` at chain-build time via
/// `take_typed_input`.
pub struct OutputQueueSet {
    pub(crate) branches: Vec<BranchEntry>,
}

pub(crate) struct BranchEntry {
    /// `Box<dyn Any>` carrying a `BranchInputHandle<T>` for the branch's `T`.
    /// Set to `Box::new(())` after `take_typed_input`.
    pub(crate) input_handle: Box<dyn Any + Send + Sync>,
    /// `Some` iff this branch's transport is a `ByteBoundedQueue`; bundles
    /// the transport-limit setter and (for an ordered branch) the reorder
    /// overflow-cap setter (see [`BranchBudgetHandles`]). The pipeline
    /// builder collects these into a registry so the budget pass +
    /// optional queue-memory rebalancer can set/reallocate budget.
    pub(crate) bounded_queue_handle: Option<BranchBudgetHandles>,
    /// `Some` iff this edge is instrumented (`--pipeline-trace`); the shared
    /// `EdgeMetrics` (also held by the transport for push counts and the input
    /// handle for pop counts). Collected into the edge registry by `contexts.rs`
    /// Pass 1.5. Cleared (`None`) by `take_typed_input`'s placeholder.
    pub(crate) metrics: Option<Arc<crate::runtime::metrics::EdgeMetrics>>,
}

impl OutputQueueSet {
    pub(crate) fn new(branches: Vec<BranchEntry>) -> Self {
        Self { branches }
    }

    #[must_use]
    pub fn n_branches(&self) -> usize {
        self.branches.len()
    }

    /// Take ownership of branch `i`'s typed input handle. Called by the
    /// consumer step's chain-build code once.
    ///
    /// # Panics
    ///
    /// Panics if `branch_idx >= n_branches()`, if the branch was already
    /// taken, or if the requested type `T` doesn't match the producer's
    /// declared branch type (build-time invariant enforced by `Chain::chain`).
    pub fn take_typed_input<T: Send + HeapSize + 'static>(
        &mut self,
        branch_idx: usize,
    ) -> BranchInputHandle<T> {
        let entry = std::mem::replace(
            &mut self.branches[branch_idx],
            BranchEntry { input_handle: Box::new(()), bounded_queue_handle: None, metrics: None },
        );
        let handle: Box<BranchInputHandle<T>> =
            entry.input_handle.downcast::<BranchInputHandle<T>>().unwrap_or_else(|_| {
                panic!("OutputQueueSet branch type mismatch at index {branch_idx}")
            });
        *handle
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Two-input handle wrapper — held inside ChainContexts.inputs[step_idx] for
// every `Step2` consumer. The `TypedStep2<S>` adapter downcasts the
// type-erased box back to `&TwoInputHandles<S::InputA, S::InputB>` per
// dispatch and exposes the per-branch refs as `ctx.a` / `ctx.b`.
//
// The wrapper itself is plain owned data; `ChainContexts` builds it once
// at chain-build time by pulling each branch's `BranchInputHandle<T>` out
// of the matching upstream's `OutputQueueSet`.
// ─────────────────────────────────────────────────────────────────────────────

/// Two per-branch input handles for a [`crate::step::Step2`]
/// consumer, paired by branch slot (`a` = consumer's input slot 0,
/// `b` = slot 1). Boxed type-erased into `ChainContexts.inputs[step_idx]`
/// at chain-build time; the
/// [`crate::erased::TypedStep2`] adapter
/// downcasts and lends per-branch references into [`crate::step::StepCtx2`].
pub struct TwoInputHandles<A, B>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
{
    pub(crate) a: BranchInputHandle<A>,
    pub(crate) b: BranchInputHandle<B>,
}

impl<A, B> TwoInputHandles<A, B>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
{
    /// Construct from two per-branch handles. Called by the
    /// chain-context builder.
    pub(crate) fn new(a: BranchInputHandle<A>, b: BranchInputHandle<B>) -> Self {
        Self { a, b }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-arity typed views — held inside OutputsViewAny.inner
// ─────────────────────────────────────────────────────────────────────────────

/// View for `Single<T>` outputs. Internal: stored inside `OutputsViewAny`
/// and accessed via `OutputHandles<Single<T>>::push` / `retry`.
pub(crate) struct SingleOutputsView<T: Send + HeapSize + 'static> {
    pub(crate) primary: BranchOutputHandle<T>,
}

impl<T: Send + HeapSize + 'static> SingleOutputsView<T> {
    pub fn mark_all_drained(&self) {
        self.primary.mark_drained();
    }
}

/// View for `OrderedBytesSingle<T>` outputs. Same shape as `SingleOutputsView`
/// but the contained `BranchOutputHandle<T>` is constructed via the
/// `_ordered_bytes` build path (which supports `T: HeapSize + Ordered` and
/// dispatches `BranchOrdering::ByItemOrdinal` correctly).
pub(crate) struct OrderedBytesSingleOutputsView<T: Send + HeapSize + Ordered + 'static> {
    pub(crate) primary: BranchOutputHandle<T>,
}

impl<T: Send + HeapSize + Ordered + 'static> OrderedBytesSingleOutputsView<T> {
    pub fn mark_all_drained(&self) {
        self.primary.mark_drained();
    }
}

/// View for `(A, B)` tuple outputs. Internal: see `SingleOutputsView`.
pub(crate) struct Tuple2OutputsView<A: Send + HeapSize + 'static, B: Send + HeapSize + 'static> {
    pub(crate) a: BranchOutputHandle<A>,
    pub(crate) b: BranchOutputHandle<B>,
}

impl<A: Send + HeapSize + 'static, B: Send + HeapSize + 'static> Tuple2OutputsView<A, B> {
    pub fn mark_all_drained(&self) {
        self.a.mark_drained();
        self.b.mark_drained();
    }
}

/// View for `(A, B, C)` tuple outputs. Internal: see `SingleOutputsView`.
pub(crate) struct Tuple3OutputsView<A, B, C>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
{
    pub(crate) a: BranchOutputHandle<A>,
    pub(crate) b: BranchOutputHandle<B>,
    pub(crate) c: BranchOutputHandle<C>,
}

impl<A, B, C> Tuple3OutputsView<A, B, C>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
{
    pub fn mark_all_drained(&self) {
        self.a.mark_drained();
        self.b.mark_drained();
        self.c.mark_drained();
    }
}

/// View for `(A, B, C, D)` tuple outputs. Internal: see `SingleOutputsView`.
pub(crate) struct Tuple4OutputsView<A, B, C, D>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
    D: Send + HeapSize + 'static,
{
    pub(crate) a: BranchOutputHandle<A>,
    pub(crate) b: BranchOutputHandle<B>,
    pub(crate) c: BranchOutputHandle<C>,
    pub(crate) d: BranchOutputHandle<D>,
}

impl<A, B, C, D> Tuple4OutputsView<A, B, C, D>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
    D: Send + HeapSize + 'static,
{
    pub fn mark_all_drained(&self) {
        self.a.mark_drained();
        self.b.mark_drained();
        self.c.mark_drained();
        self.d.mark_drained();
    }
}

/// View for `()` (sink) outputs. Internal: see `SingleOutputsView`.
/// No fields and no methods — sinks have no outputs to track or drain.
pub(crate) struct UnitOutputsView;

// ─────────────────────────────────────────────────────────────────────────────
// Typed accessors on OutputHandles<O>
// ─────────────────────────────────────────────────────────────────────────────

impl<T: Send + HeapSize + 'static> OutputHandles<Single<T>> {
    /// Push a fresh item to the (single) output.
    ///
    /// # Errors
    ///
    /// Returns `Err(Unpushed<T>)` when backpressure rejected the push. Hand
    /// the rejected item to [`Self::retry`] on the next iteration so any
    /// pre-allocated ordinal is preserved.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `SingleOutputsView<T>` (a framework invariant violation).
    #[inline]
    pub fn push(&self, item: T) -> Result<(), Unpushed<T>> {
        let view = self
            .inner
            .inner
            .downcast_ref::<SingleOutputsView<T>>()
            .expect("Single<T> outputs view downcast failed");
        view.primary.push(item)
    }

    /// Retry a previously-rejected push.
    ///
    /// # Errors
    ///
    /// Returns `Err(Unpushed<T>)` if backpressure still rejected the push.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `SingleOutputsView<T>`.
    #[inline]
    pub fn retry(&self, unpushed: Unpushed<T>) -> Result<(), Unpushed<T>> {
        let view = self
            .inner
            .inner
            .downcast_ref::<SingleOutputsView<T>>()
            .expect("Single<T> outputs view downcast failed");
        view.primary.retry(unpushed)
    }

    /// Retry a step's held output slot, re-holding on backpressure.
    ///
    /// The single canonical definition of the held-slot re-hold invariant for
    /// the `Single<T>` output shape — mirror of the `OrderedBytesSingle<T>`
    /// variant. If the slot holds an
    /// item, [`retry`](Self::retry) it; on rejection put it **back** in the slot
    /// (so it is never dropped) and report [`HeldRetry::StillHeld`]. Used by step
    /// flush-first preambles so each step doesn't re-implement the
    /// take/retry/put-back dance (a copy that forgot the put-back would silently
    /// drop a final batch). **Never spins** — the caller maps `StillHeld` to a
    /// yield (`NoProgress`/`Contention`) and retries on the next dispatch.
    #[inline]
    pub fn retry_held(&self, held: &mut crate::held::HeldSlot<Unpushed<T>>) -> HeldRetry {
        match held.take() {
            None => HeldRetry::WasEmpty,
            Some(unpushed) => match self.retry(unpushed) {
                Ok(()) => HeldRetry::Flushed,
                Err(again) => {
                    held.put(again);
                    HeldRetry::StillHeld
                }
            },
        }
    }
}

impl<T: Send + HeapSize + Ordered + 'static> OutputHandles<crate::outputs::OrderedBytesSingle<T>> {
    /// Push a fresh item to the heap-aware ordered output.
    ///
    /// # Errors
    ///
    /// Returns `Err(Unpushed<T>)` when backpressure (count or byte) rejected
    /// the push. Hand the rejected item to [`Self::retry`] on the next
    /// iteration so the pre-allocated ordinal (if any) is preserved.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `OrderedBytesSingleOutputsView<T>` (a framework invariant violation).
    #[inline]
    pub fn push(&self, item: T) -> Result<(), Unpushed<T>> {
        let view = self
            .inner
            .inner
            .downcast_ref::<OrderedBytesSingleOutputsView<T>>()
            .expect("OrderedBytesSingle<T> outputs view downcast failed");
        view.primary.push(item)
    }

    /// Retry a previously-rejected push on the heap-aware ordered output.
    ///
    /// # Errors
    ///
    /// Returns `Err(Unpushed<T>)` if backpressure still rejected the push.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `OrderedBytesSingleOutputsView<T>`.
    #[inline]
    pub fn retry(&self, unpushed: Unpushed<T>) -> Result<(), Unpushed<T>> {
        let view = self
            .inner
            .inner
            .downcast_ref::<OrderedBytesSingleOutputsView<T>>()
            .expect("OrderedBytesSingle<T> outputs view downcast failed");
        view.primary.retry(unpushed)
    }

    /// Retry a step's held output slot, re-holding on backpressure.
    ///
    /// The single canonical definition of the held-slot re-hold invariant: if
    /// the slot holds an item, [`retry`](Self::retry) it; on rejection put it
    /// **back** in the slot (so it is never dropped) and report
    /// [`HeldRetry::StillHeld`]. Used by step flush-first preambles so each
    /// step doesn't re-implement the take/retry/put-back dance (a copy that
    /// forgot the put-back would silently drop a final batch). **Never spins** —
    /// the caller maps `StillHeld` to a yield (`NoProgress`/`Contention`) and
    /// retries on the next dispatch.
    #[inline]
    pub fn retry_held(&self, held: &mut crate::held::HeldSlot<Unpushed<T>>) -> HeldRetry {
        match held.take() {
            None => HeldRetry::WasEmpty,
            Some(unpushed) => match self.retry(unpushed) {
                Ok(()) => HeldRetry::Flushed,
                Err(again) => {
                    held.put(again);
                    HeldRetry::StillHeld
                }
            },
        }
    }
}

impl<A: Send + HeapSize + 'static, B: Send + HeapSize + 'static> OutputHandles<(A, B)> {
    /// Borrow the typed per-branch view for a 2-tuple output.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `Tuple2OutputsView<A, B>` (a framework invariant violation).
    #[must_use]
    #[inline]
    pub fn view(&self) -> Tuple2View<'_, A, B> {
        let v = self
            .inner
            .inner
            .downcast_ref::<Tuple2OutputsView<A, B>>()
            .expect("Tuple-2 outputs view downcast failed");
        Tuple2View { a: &v.a, b: &v.b }
    }
}

impl<A, B> OutputHandles<crate::outputs::OrderedBytesTuple2<A, B>>
where
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
{
    /// Borrow the typed per-branch view for an ordered + byte-bounded
    /// 2-tuple output. Same view type as plain `(A, B)` (the view
    /// only carries `BranchOutputHandle`s; the ordering is encoded in
    /// the queue's transport, not the handle).
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `Tuple2OutputsView<A, B>` (a framework invariant violation).
    #[must_use]
    #[inline]
    pub fn view(&self) -> Tuple2View<'_, A, B> {
        let v = self
            .inner
            .inner
            .downcast_ref::<Tuple2OutputsView<A, B>>()
            .expect("OrderedBytesTuple2 outputs view downcast failed");
        Tuple2View { a: &v.a, b: &v.b }
    }
}

impl<A, B, C> OutputHandles<crate::outputs::OrderedBytesTuple3<A, B, C>>
where
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    C: Send + HeapSize + Ordered + 'static,
{
    /// Borrow the typed per-branch view for an ordered + byte-bounded 3-tuple
    /// output. Same view type as plain `(A, B, C)` — the view carries only
    /// `BranchOutputHandle`s; the ordering is encoded in the queue transport.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `Tuple3OutputsView<A, B, C>` (a framework invariant violation).
    #[must_use]
    #[inline]
    pub fn view(&self) -> Tuple3View<'_, A, B, C> {
        let v = self
            .inner
            .inner
            .downcast_ref::<Tuple3OutputsView<A, B, C>>()
            .expect("OrderedBytesTuple3 outputs view downcast failed");
        Tuple3View { a: &v.a, b: &v.b, c: &v.c }
    }
}

pub struct Tuple2View<'a, A: Send + HeapSize + 'static, B: Send + HeapSize + 'static> {
    pub a: &'a BranchOutputHandle<A>,
    pub b: &'a BranchOutputHandle<B>,
}

impl<A, B, C> OutputHandles<(A, B, C)>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
{
    /// Borrow the typed per-branch view for a 3-tuple output.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `Tuple3OutputsView<A, B, C>` (a framework invariant violation).
    #[must_use]
    #[inline]
    pub fn view(&self) -> Tuple3View<'_, A, B, C> {
        let v = self
            .inner
            .inner
            .downcast_ref::<Tuple3OutputsView<A, B, C>>()
            .expect("Tuple-3 outputs view downcast failed");
        Tuple3View { a: &v.a, b: &v.b, c: &v.c }
    }
}

pub struct Tuple3View<'a, A, B, C>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
{
    pub a: &'a BranchOutputHandle<A>,
    pub b: &'a BranchOutputHandle<B>,
    pub c: &'a BranchOutputHandle<C>,
}

impl<A, B, C, D> OutputHandles<(A, B, C, D)>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
    D: Send + HeapSize + 'static,
{
    /// Borrow the typed per-branch view for a 4-tuple output.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `Tuple4OutputsView<A, B, C, D>` (a framework invariant violation).
    #[must_use]
    #[inline]
    pub fn view(&self) -> Tuple4View<'_, A, B, C, D> {
        let v = self
            .inner
            .inner
            .downcast_ref::<Tuple4OutputsView<A, B, C, D>>()
            .expect("Tuple-4 outputs view downcast failed");
        Tuple4View { a: &v.a, b: &v.b, c: &v.c, d: &v.d }
    }
}

pub struct Tuple4View<'a, A, B, C, D>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
    D: Send + HeapSize + 'static,
{
    pub a: &'a BranchOutputHandle<A>,
    pub b: &'a BranchOutputHandle<B>,
    pub c: &'a BranchOutputHandle<C>,
    pub d: &'a BranchOutputHandle<D>,
}

impl OutputHandles<()> {
    /// Sinks have no output. Method exists for API symmetry.
    pub fn noop(&self) {}
}

// ─────────────────────────────────────────────────────────────────────────────
// Typed mark_all_drained on OutputHandles<O>
//
// The framework's driver calls this through `TypedStep<S>` when a step returns
// `StepOutcome::Finished` (counter-gated for `Parallel` so only the last clone
// closes the shared output). Each impl downcasts to the right per-arity view
// and forwards.
// ─────────────────────────────────────────────────────────────────────────────

impl<T: Send + HeapSize + 'static> OutputHandles<Single<T>> {
    /// Mark all output branches drained.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `SingleOutputsView<T>` (a framework invariant violation).
    pub fn mark_all_drained(&self) {
        let view = self
            .inner
            .inner
            .downcast_ref::<SingleOutputsView<T>>()
            .expect("Single<T> outputs view downcast failed");
        view.mark_all_drained();
    }
}

impl<T: Send + HeapSize + Ordered + 'static> OutputHandles<crate::outputs::OrderedBytesSingle<T>> {
    /// Mark all output branches drained.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `OrderedBytesSingleOutputsView<T>` (a framework invariant violation).
    pub fn mark_all_drained(&self) {
        let view = self
            .inner
            .inner
            .downcast_ref::<OrderedBytesSingleOutputsView<T>>()
            .expect("OrderedBytesSingle<T> outputs view downcast failed");
        view.mark_all_drained();
    }
}

impl<A: Send + HeapSize + 'static, B: Send + HeapSize + 'static> OutputHandles<(A, B)> {
    /// Mark all output branches drained.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `Tuple2OutputsView<A, B>` (a framework invariant violation).
    pub fn mark_all_drained(&self) {
        let view = self
            .inner
            .inner
            .downcast_ref::<Tuple2OutputsView<A, B>>()
            .expect("Tuple-2 outputs view downcast failed");
        view.mark_all_drained();
    }
}

impl<A, B> OutputHandles<crate::outputs::OrderedBytesTuple2<A, B>>
where
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
{
    /// Mark all output branches drained.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `Tuple2OutputsView<A, B>` (the view shape is shared with plain
    /// `(A, B)`).
    pub fn mark_all_drained(&self) {
        let view = self
            .inner
            .inner
            .downcast_ref::<Tuple2OutputsView<A, B>>()
            .expect("OrderedBytesTuple2 outputs view downcast failed");
        view.mark_all_drained();
    }
}

impl<A, B, C> OutputHandles<crate::outputs::OrderedBytesTuple3<A, B, C>>
where
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    C: Send + HeapSize + Ordered + 'static,
{
    /// Mark all output branches drained.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `Tuple3OutputsView<A, B, C>` (the view shape is shared with plain
    /// `(A, B, C)`).
    pub fn mark_all_drained(&self) {
        let view = self
            .inner
            .inner
            .downcast_ref::<Tuple3OutputsView<A, B, C>>()
            .expect("OrderedBytesTuple3 outputs view downcast failed");
        view.mark_all_drained();
    }
}

impl<A, B, C> OutputHandles<(A, B, C)>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
{
    /// Mark all output branches drained.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `Tuple3OutputsView<A, B, C>` (a framework invariant violation).
    pub fn mark_all_drained(&self) {
        let view = self
            .inner
            .inner
            .downcast_ref::<Tuple3OutputsView<A, B, C>>()
            .expect("Tuple-3 outputs view downcast failed");
        view.mark_all_drained();
    }
}

impl<A, B, C, D> OutputHandles<(A, B, C, D)>
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
    D: Send + HeapSize + 'static,
{
    /// Mark all output branches drained.
    ///
    /// # Panics
    ///
    /// Panics if the type-erased outputs view doesn't downcast to
    /// `Tuple4OutputsView<A, B, C, D>` (a framework invariant violation).
    pub fn mark_all_drained(&self) {
        let view = self
            .inner
            .inner
            .downcast_ref::<Tuple4OutputsView<A, B, C, D>>()
            .expect("Tuple-4 outputs view downcast failed");
        view.mark_all_drained();
    }
}

impl OutputHandles<()> {
    /// Mark all output branches drained. Sinks have no outputs; no-op.
    pub fn mark_all_drained(&self) {
        // No-op for unit outputs.
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// build_*_queues per arity — entry points called from outputs.rs
// ─────────────────────────────────────────────────────────────────────────────

pub(crate) fn build_single_queues<T: Send + HeapSize + 'static>(
    specs: &[QueueSpec],
    ordering: &[BranchOrdering],
    level: crate::builder::InstrumentationLevel,
) -> (OutputQueueSet, OutputsViewAny) {
    assert_eq!(specs.len(), 1, "Single<T>::build_queues requires 1 spec");
    assert_eq!(ordering.len(), 1, "Single<T>::build_queues requires 1 ordering");

    // `Single<T>` bounds `T: HeapSize`, so use the byte-aware build path: it
    // honors `QueueSpec::ByteBounded` (documented as supported for `Single<T>`
    // outputs without item-carried serials) and delegates every non-byte spec
    // straight back to `build_branch::<T>`, so count/unbounded paths are
    // unchanged.
    let branch = build_branch_byte_aware::<T>(specs[0], ordering[0], level);
    let view = SingleOutputsView { primary: branch.output };
    let outputs_view = OutputsViewAny { inner: Box::new(view) };
    let queue_set = OutputQueueSet::new(vec![BranchEntry {
        input_handle: Box::new(branch.input),
        bounded_queue_handle: branch.bounded_queue_handle,
        metrics: branch.metrics,
    }]);
    (queue_set, outputs_view)
}

/// Build queues for `OrderedBytesSingle<T>` outputs (`T: HeapSize + Ordered`).
/// Supports every `QueueSpec` × `BranchOrdering` combination.
pub(crate) fn build_single_queues_ordered_bytes<T: Send + HeapSize + Ordered + 'static>(
    specs: &[QueueSpec],
    ordering: &[BranchOrdering],
    level: crate::builder::InstrumentationLevel,
) -> (OutputQueueSet, OutputsViewAny) {
    assert_eq!(specs.len(), 1, "OrderedBytesSingle<T>::build_queues requires 1 spec");
    assert_eq!(ordering.len(), 1, "OrderedBytesSingle<T>::build_queues requires 1 ordering");

    let branch = build_branch_ordered_bytes::<T>(specs[0], ordering[0], level);
    let view = OrderedBytesSingleOutputsView { primary: branch.output };
    let outputs_view = OutputsViewAny { inner: Box::new(view) };
    let queue_set = OutputQueueSet::new(vec![BranchEntry {
        input_handle: Box::new(branch.input),
        bounded_queue_handle: branch.bounded_queue_handle,
        metrics: branch.metrics,
    }]);
    (queue_set, outputs_view)
}

pub(crate) fn build_tuple2_queues<A, B>(
    specs: &[QueueSpec],
    ordering: &[BranchOrdering],
    level: crate::builder::InstrumentationLevel,
) -> (OutputQueueSet, OutputsViewAny)
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
{
    assert_eq!(specs.len(), 2, "(A, B)::build_queues requires 2 specs");
    assert_eq!(ordering.len(), 2, "(A, B)::build_queues requires 2 orderings");

    let ba = build_branch::<A>(specs[0], ordering[0], level);
    let bb = build_branch::<B>(specs[1], ordering[1], level);
    let view = Tuple2OutputsView { a: ba.output, b: bb.output };
    let outputs_view = OutputsViewAny { inner: Box::new(view) };
    let queue_set = OutputQueueSet::new(vec![
        BranchEntry {
            input_handle: Box::new(ba.input),
            bounded_queue_handle: ba.bounded_queue_handle,
            metrics: ba.metrics,
        },
        BranchEntry {
            input_handle: Box::new(bb.input),
            bounded_queue_handle: bb.bounded_queue_handle,
            metrics: bb.metrics,
        },
    ]);
    (queue_set, outputs_view)
}

/// Build queues for `OrderedBytesTuple2<A, B>` outputs. Both branches
/// support the full `QueueSpec × BranchOrdering` cross-product because
/// `A: Ordered + HeapSize` and `B: Ordered + HeapSize`. The
/// `Tuple2OutputsView` carries only `BranchOutputHandle`s (no ordering
/// metadata), so the view type is the same as for plain `(A, B)`.
pub(crate) fn build_tuple2_queues_ordered_bytes<A, B>(
    specs: &[QueueSpec],
    ordering: &[BranchOrdering],
    level: crate::builder::InstrumentationLevel,
) -> (OutputQueueSet, OutputsViewAny)
where
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
{
    assert_eq!(specs.len(), 2, "OrderedBytesTuple2<A, B>::build_queues requires 2 specs");
    assert_eq!(ordering.len(), 2, "OrderedBytesTuple2<A, B>::build_queues requires 2 orderings");

    let ba = build_branch_ordered_bytes::<A>(specs[0], ordering[0], level);
    let bb = build_branch_ordered_bytes::<B>(specs[1], ordering[1], level);
    let view = Tuple2OutputsView { a: ba.output, b: bb.output };
    let outputs_view = OutputsViewAny { inner: Box::new(view) };
    let queue_set = OutputQueueSet::new(vec![
        BranchEntry {
            input_handle: Box::new(ba.input),
            bounded_queue_handle: ba.bounded_queue_handle,
            metrics: ba.metrics,
        },
        BranchEntry {
            input_handle: Box::new(bb.input),
            bounded_queue_handle: bb.bounded_queue_handle,
            metrics: bb.metrics,
        },
    ]);
    (queue_set, outputs_view)
}

pub(crate) fn build_tuple3_queues<A, B, C>(
    specs: &[QueueSpec],
    ordering: &[BranchOrdering],
    level: crate::builder::InstrumentationLevel,
) -> (OutputQueueSet, OutputsViewAny)
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
{
    assert_eq!(specs.len(), 3, "(A, B, C)::build_queues requires 3 specs");
    assert_eq!(ordering.len(), 3, "(A, B, C)::build_queues requires 3 orderings");

    let ba = build_branch::<A>(specs[0], ordering[0], level);
    let bb = build_branch::<B>(specs[1], ordering[1], level);
    let bc = build_branch::<C>(specs[2], ordering[2], level);
    let view = Tuple3OutputsView { a: ba.output, b: bb.output, c: bc.output };
    let outputs_view = OutputsViewAny { inner: Box::new(view) };
    let queue_set = OutputQueueSet::new(vec![
        BranchEntry {
            input_handle: Box::new(ba.input),
            bounded_queue_handle: ba.bounded_queue_handle,
            metrics: ba.metrics,
        },
        BranchEntry {
            input_handle: Box::new(bb.input),
            bounded_queue_handle: bb.bounded_queue_handle,
            metrics: bb.metrics,
        },
        BranchEntry {
            input_handle: Box::new(bc.input),
            bounded_queue_handle: bc.bounded_queue_handle,
            metrics: bc.metrics,
        },
    ]);
    (queue_set, outputs_view)
}

/// Build queues for `OrderedBytesTuple3<A, B, C>` outputs. All three branches
/// support the full `QueueSpec × BranchOrdering` cross-product because each is
/// `Ordered + HeapSize`. The `Tuple3OutputsView` carries only
/// `BranchOutputHandle`s (no ordering metadata), so the view type is the same
/// as for plain `(A, B, C)`.
pub(crate) fn build_tuple3_queues_ordered_bytes<A, B, C>(
    specs: &[QueueSpec],
    ordering: &[BranchOrdering],
    level: crate::builder::InstrumentationLevel,
) -> (OutputQueueSet, OutputsViewAny)
where
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    C: Send + HeapSize + Ordered + 'static,
{
    assert_eq!(specs.len(), 3, "OrderedBytesTuple3<A, B, C>::build_queues requires 3 specs");
    assert_eq!(ordering.len(), 3, "OrderedBytesTuple3<A, B, C>::build_queues requires 3 orderings");

    let ba = build_branch_ordered_bytes::<A>(specs[0], ordering[0], level);
    let bb = build_branch_ordered_bytes::<B>(specs[1], ordering[1], level);
    let bc = build_branch_ordered_bytes::<C>(specs[2], ordering[2], level);
    let view = Tuple3OutputsView { a: ba.output, b: bb.output, c: bc.output };
    let outputs_view = OutputsViewAny { inner: Box::new(view) };
    let queue_set = OutputQueueSet::new(vec![
        BranchEntry {
            input_handle: Box::new(ba.input),
            bounded_queue_handle: ba.bounded_queue_handle,
            metrics: ba.metrics,
        },
        BranchEntry {
            input_handle: Box::new(bb.input),
            bounded_queue_handle: bb.bounded_queue_handle,
            metrics: bb.metrics,
        },
        BranchEntry {
            input_handle: Box::new(bc.input),
            bounded_queue_handle: bc.bounded_queue_handle,
            metrics: bc.metrics,
        },
    ]);
    (queue_set, outputs_view)
}

pub(crate) fn build_tuple4_queues<A, B, C, D>(
    specs: &[QueueSpec],
    ordering: &[BranchOrdering],
    level: crate::builder::InstrumentationLevel,
) -> (OutputQueueSet, OutputsViewAny)
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
    D: Send + HeapSize + 'static,
{
    assert_eq!(specs.len(), 4, "(A, B, C, D)::build_queues requires 4 specs");
    assert_eq!(ordering.len(), 4, "(A, B, C, D)::build_queues requires 4 orderings");

    let ba = build_branch::<A>(specs[0], ordering[0], level);
    let bb = build_branch::<B>(specs[1], ordering[1], level);
    let bc = build_branch::<C>(specs[2], ordering[2], level);
    let bd = build_branch::<D>(specs[3], ordering[3], level);
    let view = Tuple4OutputsView { a: ba.output, b: bb.output, c: bc.output, d: bd.output };
    let outputs_view = OutputsViewAny { inner: Box::new(view) };
    let queue_set = OutputQueueSet::new(vec![
        BranchEntry {
            input_handle: Box::new(ba.input),
            bounded_queue_handle: ba.bounded_queue_handle,
            metrics: ba.metrics,
        },
        BranchEntry {
            input_handle: Box::new(bb.input),
            bounded_queue_handle: bb.bounded_queue_handle,
            metrics: bb.metrics,
        },
        BranchEntry {
            input_handle: Box::new(bc.input),
            bounded_queue_handle: bc.bounded_queue_handle,
            metrics: bc.metrics,
        },
        BranchEntry {
            input_handle: Box::new(bd.input),
            bounded_queue_handle: bd.bounded_queue_handle,
            metrics: bd.metrics,
        },
    ]);
    (queue_set, outputs_view)
}

pub(crate) fn build_unit_queues(
    specs: &[QueueSpec],
    ordering: &[BranchOrdering],
    _level: crate::builder::InstrumentationLevel,
) -> (OutputQueueSet, OutputsViewAny) {
    assert_eq!(specs.len(), 0, "()::build_queues requires 0 specs");
    assert_eq!(ordering.len(), 0, "()::build_queues requires 0 orderings");
    let outputs_view = OutputsViewAny { inner: Box::new(UnitOutputsView) };
    let queue_set = OutputQueueSet::new(Vec::new());
    (queue_set, outputs_view)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod handle_tests {
    use super::*;

    #[test]
    fn always_drained_handle_is_empty_and_drained() {
        // The zero-state source input handle owns no transport: it pops nothing
        // and reports drained from the start, so a source step sees its
        // (implicit) input as immediately end-of-stream.
        let h = BranchInputHandle::<()>::always_drained();
        assert_eq!(h.pop(), None, "always-drained handle yields no items");
        assert!(h.is_drained(), "always-drained handle reports drained");
        // Idempotent: still drained, still empty after repeated reads.
        assert_eq!(h.pop(), None);
        assert!(h.is_drained());
    }

    #[test]
    fn build_branch_mints_metrics_iff_on() {
        use crate::builder::InstrumentationLevel as L;
        let on = build_branch::<u32>(
            QueueSpec::CountBounded { capacity: 4 },
            BranchOrdering::None,
            L::Summary,
        );
        assert!(on.metrics.is_some(), "level on → edge metrics minted");
        let off = build_branch::<u32>(
            QueueSpec::CountBounded { capacity: 4 },
            BranchOrdering::None,
            L::Off,
        );
        assert!(off.metrics.is_none(), "level off → no metrics (hot path metric-free)");
    }

    #[test]
    fn metrics_shared_between_transport_push_and_input_handle_pop() {
        // Producer-push (transport) and consumer-pop (input handle) share one
        // EdgeMetrics: push lands via the queue, pop/empty via the input handle.
        use crate::builder::InstrumentationLevel as L;
        let b = build_branch::<u32>(
            QueueSpec::CountBounded { capacity: 4 },
            BranchOrdering::None,
            L::Summary,
        );
        let m = b.metrics.clone().expect("metrics present");
        b.output.push(7).unwrap();
        assert_eq!(b.input.pop(), Some(7));
        assert_eq!(b.input.pop(), None); // empty
        let s = m.snapshot();
        assert_eq!(s.pushed_items, 1, "producer push counted at the transport");
        assert_eq!(s.popped_items, 1, "consumer pop counted at the input handle");
        assert_eq!(s.pop_empties, 1, "empty pop counted at the input handle");
    }

    #[test]
    fn ordered_reorder_blocked_pop_is_not_counted_empty() {
        // Regression: an ordered edge's `try_pop_in_order` returns `None` both
        // when the edge is genuinely starved AND when it is only reorder-blocked
        // (later ordinals are buffered while it waits for the next in-order
        // ordinal). Counting the reorder-blocked case as an empty pop inflates
        // `pop_empties` and can misclassify a backlogged edge as starved, so
        // the pop path must skip `record_empty()` while the reorder buffer holds
        // out-of-order items.
        use crate::runtime::metrics::EdgeMetrics;

        let transport: Arc<dyn ItemQueue<Sequenced<u32>>> =
            Arc::new(CountBoundedQueue::<Sequenced<u32>>::new(8));
        let stage = Arc::new(ReorderStage::new(transport));
        // Buffer ordinal 1 while ordinal 0 is still absent: the next in-order
        // pop is blocked, not starved.
        stage.try_push(1, 100).unwrap();

        let m = EdgeMetrics::new();
        let handle = BranchInputHandle {
            inner: BranchInputInner::Ordered(stage.clone()),
            metrics: Some(m.clone()),
        };

        // Reorder-blocked: yields nothing yet, but there is buffered work.
        assert_eq!(handle.pop(), None, "next ordinal (0) absent → no in-order item");
        assert_eq!(
            m.snapshot().pop_empties,
            0,
            "a reorder-blocked pop is backlog, not starvation — must not count as empty"
        );

        // Now the missing ordinal arrives and both items drain in order; still
        // no empty pops recorded.
        stage.try_push(0, 50).unwrap();
        assert_eq!(handle.pop(), Some(50));
        assert_eq!(handle.pop(), Some(100));
        assert_eq!(m.snapshot().pop_empties, 0, "in-order drains record no empties");

        // Genuinely drained now: this pop IS a starved/empty pop.
        assert_eq!(handle.pop(), None);
        assert_eq!(
            m.snapshot().pop_empties,
            1,
            "an empty pop with no buffered work is a true empty pop"
        );
    }

    #[test]
    fn count_bounded_fifo_round_trip() {
        let b = build_branch::<u32>(
            QueueSpec::CountBounded { capacity: 4 },
            BranchOrdering::None,
            crate::builder::InstrumentationLevel::Off,
        );
        b.output.push(1).unwrap();
        b.output.push(2).unwrap();
        assert_eq!(b.input.pop(), Some(1));
        assert_eq!(b.input.pop(), Some(2));
        assert_eq!(b.input.pop(), None);
    }

    #[test]
    fn count_bounded_ordered_emits_in_order() {
        let b = build_branch::<u32>(
            QueueSpec::CountBounded { capacity: 8 },
            BranchOrdering::ByOrdinal,
            crate::builder::InstrumentationLevel::Off,
        );
        // Single producer pushes ordinals 0,1,2 in arrival order.
        b.output.push(100).unwrap();
        b.output.push(200).unwrap();
        b.output.push(300).unwrap();
        assert_eq!(b.input.pop(), Some(100));
        assert_eq!(b.input.pop(), Some(200));
        assert_eq!(b.input.pop(), Some(300));
    }

    #[test]
    fn drained_signal_propagates() {
        let b = build_branch::<u32>(
            QueueSpec::CountBounded { capacity: 4 },
            BranchOrdering::None,
            crate::builder::InstrumentationLevel::Off,
        );
        b.output.push(1).unwrap();
        b.output.mark_drained();
        // Marker set but item still buffered: not drained.
        assert!(!b.input.is_drained());
        assert_eq!(b.input.pop(), Some(1));
        assert!(b.input.is_drained());
    }

    #[test]
    fn unbounded_branch_works() {
        let b = build_branch::<u32>(
            QueueSpec::Unbounded,
            BranchOrdering::None,
            crate::builder::InstrumentationLevel::Off,
        );
        for i in 0..1024 {
            b.output.push(i).unwrap();
        }
        for i in 0..1024 {
            assert_eq!(b.input.pop(), Some(i));
        }
    }

    #[derive(Debug)]
    struct Bytes(Vec<u8>);
    impl HeapSize for Bytes {
        fn heap_size(&self) -> usize {
            self.0.len()
        }
    }

    #[test]
    fn byte_bounded_branch_works() {
        let b = build_branch_byte_aware::<Bytes>(
            QueueSpec::ByteBounded { limit_bytes: 1000 },
            BranchOrdering::None,
            crate::builder::InstrumentationLevel::Off,
        );
        b.output.push(Bytes(vec![0; 500])).unwrap();
        let popped = b.input.pop().unwrap();
        assert_eq!(popped.0.len(), 500);
    }

    #[test]
    #[should_panic(expected = "ByteBounded requires `T: HeapSize`")]
    fn byte_bounded_panics_in_non_heap_aware_builder() {
        let _ = build_branch::<u32>(
            QueueSpec::ByteBounded { limit_bytes: 1000 },
            BranchOrdering::None,
            crate::builder::InstrumentationLevel::Off,
        );
    }

    #[test]
    fn byte_bounded_single_builds_via_user_facing_build_queues() {
        // Regression: a user step declaring `QueueSpec::ByteBounded` in its
        // `StepProfile::output_queues` with a `Single<T>` output reaches
        // `StepOutputs::build_queues` → `build_single_queues`. `Single<T>`
        // bounds `T: HeapSize`, so this builds a byte-bounded queue (the
        // documented "byte-bounded without item-carried serials" shape)
        // rather than panicking.
        use crate::outputs::{Single, StepOutputs};
        use crate::step::OutputHandles;
        let (mut queue_set, outputs_view) = <Single<u32> as StepOutputs>::build_queues(
            &[QueueSpec::ByteBounded { limit_bytes: 1000 }],
            &[BranchOrdering::None],
            crate::builder::InstrumentationLevel::Off,
        );
        let outputs: OutputHandles<Single<u32>> = OutputHandles::new(outputs_view);
        outputs.push(7).unwrap();
        let input = queue_set.take_typed_input::<u32>(0);
        assert_eq!(input.pop(), Some(7));
    }

    #[test]
    fn single_retry_held_drains_then_reholds() {
        // Exercises `OutputHandles<Single<T>>::retry_held`: empty → WasEmpty,
        // a rejected push re-held → StillHeld while the queue is full → Flushed
        // once a slot frees, with the held item never dropped.
        use crate::held::HeldSlot;
        use crate::outputs::{Single, StepOutputs};
        use crate::step::OutputHandles;
        let (mut queue_set, outputs_view) = <Single<u32> as StepOutputs>::build_queues(
            &[QueueSpec::CountBounded { capacity: 1 }],
            &[BranchOrdering::None],
            crate::builder::InstrumentationLevel::Off,
        );
        let outputs: OutputHandles<Single<u32>> = OutputHandles::new(outputs_view);
        let mut held: HeldSlot<Unpushed<u32>> = HeldSlot::new();

        // Empty slot → WasEmpty (queue untouched).
        assert!(matches!(outputs.retry_held(&mut held), HeldRetry::WasEmpty));

        // Fill the capacity-1 queue; the next push is rejected and re-held.
        outputs.push(1).unwrap();
        let rejected = outputs.push(2).expect_err("capacity-1 queue must reject the 2nd push");
        held.put(rejected);

        // Queue still full → StillHeld; the item is put back, not dropped.
        assert!(matches!(outputs.retry_held(&mut held), HeldRetry::StillHeld));
        assert!(held.is_held());

        // Drain one item; the held push now flushes.
        let input = queue_set.take_typed_input::<u32>(0);
        assert_eq!(input.pop(), Some(1));
        assert!(matches!(outputs.retry_held(&mut held), HeldRetry::Flushed));
        assert!(!held.is_held());
        assert_eq!(input.pop(), Some(2));
    }

    #[test]
    fn byte_bounded_plus_byordinal_works_after_phase3_amendment() {
        // Phase 3 amendment 2 lifted the deferred ByteBounded + ByOrdinal
        // path: `Sequenced<T>: HeapSize` is now impl'd, so a byte-bounded
        // transport can wrap `Sequenced<T>` and the `Allocated` ordinal
        // source provides per-branch ordinals.
        let b = build_branch_byte_aware::<Bytes>(
            QueueSpec::ByteBounded { limit_bytes: 1000 },
            BranchOrdering::ByOrdinal,
            crate::builder::InstrumentationLevel::Off,
        );
        b.output.push(Bytes(vec![0; 100])).unwrap();
        let popped = b.input.pop().unwrap();
        assert_eq!(popped.0.len(), 100);
    }

    #[test]
    #[should_panic(expected = "BranchOrdering::ByItemOrdinal requires `T: Ordered`")]
    fn by_item_ordinal_panics_in_non_ordered_builder() {
        let _ = build_branch::<u32>(
            QueueSpec::CountBounded { capacity: 4 },
            BranchOrdering::ByItemOrdinal,
            crate::builder::InstrumentationLevel::Off,
        );
    }

    /// Trivial `Ordered`-impl test fixture for the by-item-ordinal builders.
    #[derive(Debug, PartialEq, Eq)]
    struct OrdItem {
        ord: u64,
        v: u32,
    }
    impl crate::item::Ordered for OrdItem {
        fn ordinal(&self) -> u64 {
            self.ord
        }
    }

    #[test]
    fn by_item_ordinal_uses_item_serial() {
        let b = build_branch_ordered::<OrdItem>(
            QueueSpec::CountBounded { capacity: 8 },
            BranchOrdering::ByItemOrdinal,
            crate::builder::InstrumentationLevel::Off,
        );
        // Push out of order: ordinals 2, 0, 1.
        b.output.push(OrdItem { ord: 2, v: 200 }).unwrap();
        b.output.push(OrdItem { ord: 0, v: 0 }).unwrap();
        b.output.push(OrdItem { ord: 1, v: 100 }).unwrap();
        // Consumer sees them in item-ordinal order.
        assert_eq!(b.input.pop().unwrap().v, 0);
        assert_eq!(b.input.pop().unwrap().v, 100);
        assert_eq!(b.input.pop().unwrap().v, 200);
    }

    impl HeapSize for OrdItem {
        fn heap_size(&self) -> usize {
            std::mem::size_of::<Self>()
        }
    }

    #[test]
    fn ordered_bytes_supports_byte_bounded_with_item_serial() {
        let b = build_branch_ordered_bytes::<OrdItem>(
            QueueSpec::ByteBounded { limit_bytes: 4096 },
            BranchOrdering::ByItemOrdinal,
            crate::builder::InstrumentationLevel::Off,
        );
        b.output.push(OrdItem { ord: 1, v: 1 }).unwrap();
        b.output.push(OrdItem { ord: 0, v: 0 }).unwrap();
        assert_eq!(b.input.pop().unwrap().v, 0);
        assert_eq!(b.input.pop().unwrap().v, 1);
    }

    #[test]
    fn ordered_branch_preserves_ordinal_across_retry() {
        // Regression for C1 (ordinal-burn): with capacity 2, push three items
        // through an ordered branch, draining one at a time. Without the
        // Unpushed<T> retry, the third push would burn an ordinal on
        // backpressure and stall the reorder stage.
        let b = build_branch::<u32>(
            QueueSpec::CountBounded { capacity: 2 },
            BranchOrdering::ByOrdinal,
            crate::builder::InstrumentationLevel::Off,
        );

        // Pump 5 items through the branch, draining sequentially. Capacity 2
        // means every push after the second hits backpressure once.
        let mut received = Vec::new();
        for n in 0..5_u32 {
            // Try to push. Retry on backpressure until accepted, draining the
            // input handle in between to make progress.
            let mut held: Option<Unpushed<u32>> = None;
            let mut fresh: Option<u32> = Some(n);

            loop {
                if let Some(unpushed) = held.take() {
                    match b.output.retry(unpushed) {
                        Ok(()) => {}
                        Err(again) => {
                            held = Some(again);
                        }
                    }
                }
                if held.is_none()
                    && let Some(item) = fresh.take()
                {
                    match b.output.push(item) {
                        Ok(()) => {}
                        Err(unpushed) => {
                            held = Some(unpushed);
                        }
                    }
                }
                if held.is_none() && fresh.is_none() {
                    break;
                }
                // Drain any items the consumer has ready, freeing transport space.
                while let Some(v) = b.input.pop() {
                    received.push(v);
                }
            }
        }
        // Final drain — drain everything still buffered.
        while let Some(v) = b.input.pop() {
            received.push(v);
        }
        assert_eq!(received, vec![0, 1, 2, 3, 4], "ordinals preserved across retries");
    }

    #[test]
    fn output_queue_set_take_typed_input() {
        let (mut set, _view) = build_single_queues::<u32>(
            &[QueueSpec::CountBounded { capacity: 4 }],
            &[BranchOrdering::None],
            crate::builder::InstrumentationLevel::Off,
        );
        let _input: BranchInputHandle<u32> = set.take_typed_input::<u32>(0);
    }
}

#[cfg(test)]
mod build_queues_tests {
    use super::*;
    use crate::outputs::{Single, StepOutputs};
    use crate::step::OutputHandles;

    fn count_specs(arity: usize, capacity: usize) -> (Vec<QueueSpec>, Vec<BranchOrdering>) {
        (vec![QueueSpec::CountBounded { capacity }; arity], vec![BranchOrdering::None; arity])
    }

    #[test]
    fn single_round_trip() {
        let (specs, ordering) = count_specs(1, 4);
        let (mut queue_set, outputs_view) = <Single<u32> as StepOutputs>::build_queues(
            &specs,
            &ordering,
            crate::builder::InstrumentationLevel::Off,
        );
        let outputs: OutputHandles<Single<u32>> = OutputHandles::new(outputs_view);
        outputs.push(7).unwrap();

        let input = queue_set.take_typed_input::<u32>(0);
        assert_eq!(input.pop(), Some(7));
    }

    #[test]
    fn tuple_2_round_trip() {
        let (specs, ordering) = count_specs(2, 2);
        let (mut queue_set, outputs_view) = <(u32, String) as StepOutputs>::build_queues(
            &specs,
            &ordering,
            crate::builder::InstrumentationLevel::Off,
        );
        let outputs: OutputHandles<(u32, String)> = OutputHandles::new(outputs_view);

        let v = outputs.view();
        v.a.push(10).unwrap();
        v.b.push("hello".to_string()).unwrap();

        let in_a = queue_set.take_typed_input::<u32>(0);
        let in_b = queue_set.take_typed_input::<String>(1);
        assert_eq!(in_a.pop(), Some(10));
        assert_eq!(in_b.pop(), Some("hello".to_string()));
    }

    #[test]
    fn tuple_3_round_trip() {
        let (specs, ordering) = count_specs(3, 2);
        let (mut queue_set, outputs_view) = <(u32, u64, String) as StepOutputs>::build_queues(
            &specs,
            &ordering,
            crate::builder::InstrumentationLevel::Off,
        );
        let outputs: OutputHandles<(u32, u64, String)> = OutputHandles::new(outputs_view);

        let v = outputs.view();
        v.a.push(1).unwrap();
        v.b.push(2).unwrap();
        v.c.push("three".to_string()).unwrap();

        assert_eq!(queue_set.take_typed_input::<u32>(0).pop(), Some(1));
        assert_eq!(queue_set.take_typed_input::<u64>(1).pop(), Some(2));
        assert_eq!(queue_set.take_typed_input::<String>(2).pop(), Some("three".to_string()));
    }

    #[test]
    fn tuple_4_round_trip() {
        let (specs, ordering) = count_specs(4, 2);
        let (mut queue_set, outputs_view) =
            <(u32, u64, String, Vec<u8>) as StepOutputs>::build_queues(
                &specs,
                &ordering,
                crate::builder::InstrumentationLevel::Off,
            );
        let outputs: OutputHandles<(u32, u64, String, Vec<u8>)> = OutputHandles::new(outputs_view);

        let v = outputs.view();
        v.a.push(1).unwrap();
        v.b.push(2).unwrap();
        v.c.push("three".to_string()).unwrap();
        v.d.push(vec![4u8, 5, 6]).unwrap();

        assert_eq!(queue_set.take_typed_input::<u32>(0).pop(), Some(1));
        assert_eq!(queue_set.take_typed_input::<u64>(1).pop(), Some(2));
        assert_eq!(queue_set.take_typed_input::<String>(2).pop(), Some("three".to_string()));
        assert_eq!(queue_set.take_typed_input::<Vec<u8>>(3).pop(), Some(vec![4u8, 5, 6]));
    }

    #[test]
    fn unit_build_queues_yields_empty_set() {
        let (queue_set, outputs_view) =
            <() as StepOutputs>::build_queues(&[], &[], crate::builder::InstrumentationLevel::Off);
        let outputs: OutputHandles<()> = OutputHandles::new(outputs_view);
        outputs.noop();
        assert_eq!(queue_set.n_branches(), 0);
    }

    #[test]
    fn mark_all_drained_propagates_to_input_handle() {
        let (specs, ordering) = count_specs(1, 4);
        let (mut queue_set, outputs_view) = <Single<u32> as StepOutputs>::build_queues(
            &specs,
            &ordering,
            crate::builder::InstrumentationLevel::Off,
        );
        let outputs: OutputHandles<Single<u32>> = OutputHandles::new(outputs_view);

        let input = queue_set.take_typed_input::<u32>(0);
        assert!(!input.is_drained());
        outputs.mark_all_drained();
        assert!(input.is_drained());
    }

    #[test]
    fn ordered_branch_preserves_emission_order_through_typed_path() {
        let specs = vec![QueueSpec::CountBounded { capacity: 8 }];
        let ordering = vec![BranchOrdering::ByOrdinal];
        let (mut queue_set, outputs_view) = <Single<u32> as StepOutputs>::build_queues(
            &specs,
            &ordering,
            crate::builder::InstrumentationLevel::Off,
        );
        let outputs: OutputHandles<Single<u32>> = OutputHandles::new(outputs_view);

        outputs.push(10).unwrap();
        outputs.push(20).unwrap();
        outputs.push(30).unwrap();

        let input = queue_set.take_typed_input::<u32>(0);
        assert_eq!(input.pop(), Some(10));
        assert_eq!(input.pop(), Some(20));
        assert_eq!(input.pop(), Some(30));
    }
}
