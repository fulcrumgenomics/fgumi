//! Bounded, memory-tracked stage queue.
//!
//! Wraps `crossbeam_queue::ArrayQueue` with:
//! - Slot-based bound (`ArrayQueue` capacity).
//! - Memory-based bound (atomic byte counter + configured limit).
//! - Closed flag for graceful shutdown propagation.
//!
//! When `push` fails due to either slot or memory limit being hit, the item
//! is returned to the caller. Callers (workers) store the rejected item in
//! a held-item slot (see `worker.rs`) and retry on a later iteration.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use crossbeam_queue::ArrayQueue;

use super::backpressure::QueueSummary;
use super::deadlock::QueueProgress;
use super::memory::MemoryTracker;
use super::stage::SequencedItem;

/// Bounded, memory-tracked queue between two pipeline stages.
pub struct StageQueue<T> {
    inner: Arc<ArrayQueue<SequencedItem<T>>>,
    queue_bytes: Arc<AtomicUsize>,
    queue_limit_bytes: usize,
    global_memory: Arc<MemoryTracker>,
    closed: Arc<AtomicBool>,
    /// Human-readable queue name for logs/diagnostics. Stored as `String`
    /// so callers can embed stage indices or other runtime metadata
    /// (e.g. `"stage_3_out"`).
    name: String,
    /// Optional progress counters for deadlock detection. When present,
    /// every successful push and pop is reported to the watchdog.
    progress: Option<Arc<QueueProgress>>,
}

impl<T> StageQueue<T> {
    /// Construct a new queue with `capacity` slots and per-queue memory limit.
    pub fn new(
        name: impl Into<String>,
        capacity: usize,
        queue_limit_bytes: usize,
        global_memory: Arc<MemoryTracker>,
    ) -> Self {
        Self {
            inner: Arc::new(ArrayQueue::new(capacity)),
            queue_bytes: Arc::new(AtomicUsize::new(0)),
            queue_limit_bytes,
            global_memory,
            closed: Arc::new(AtomicBool::new(false)),
            name: name.into(),
            progress: None,
        }
    }

    /// Construct a new queue with an attached `QueueProgress` counter.
    ///
    /// The counter is incremented on every successful push/pop. The pipeline
    /// driver hands the same counter to the deadlock watchdog so it can
    /// observe real queue traffic rather than a placeholder signal.
    pub fn with_progress(
        name: impl Into<String>,
        capacity: usize,
        queue_limit_bytes: usize,
        global_memory: Arc<MemoryTracker>,
        progress: Arc<QueueProgress>,
    ) -> Self {
        Self {
            inner: Arc::new(ArrayQueue::new(capacity)),
            queue_bytes: Arc::new(AtomicUsize::new(0)),
            queue_limit_bytes,
            global_memory,
            closed: Arc::new(AtomicBool::new(false)),
            name: name.into(),
            progress: Some(progress),
        }
    }

    /// Attempt to push an item. Returns `Err(item)` if the queue is full
    /// (slot, per-queue memory, or global memory limit).
    ///
    /// Memory reservations go through [`MemoryTracker::try_add`] so the global
    /// limit is enforced atomically — concurrent producers across all queues
    /// race on the same CAS. A failed reservation returns the item to the
    /// caller unchanged.
    ///
    /// # Errors
    ///
    /// Returns the item if the queue has reached any of: the slot capacity,
    /// the per-queue memory limit, or the global memory limit.
    pub fn push(&self, item: SequencedItem<T>) -> Result<(), SequencedItem<T>> {
        let bytes = item.memory_estimate;
        let current = self.queue_bytes.load(Ordering::Relaxed);
        // Use `saturating_add` so a transient racy overcount (pop updates
        // `queue_bytes` after popping, which can briefly wrap near zero
        // under heavy concurrent pool traffic) doesn't panic on overflow
        // in debug builds.
        if current.saturating_add(bytes) > self.queue_limit_bytes {
            return Err(item);
        }
        // Reserve global bytes BEFORE touching the underlying queue. If we
        // fail here, nothing has been mutated and the item is returned as-is.
        if !self.global_memory.try_add(bytes) {
            return Err(item);
        }
        match self.inner.push(item) {
            Ok(()) => {
                self.queue_bytes.fetch_add(bytes, Ordering::Relaxed);
                if let Some(progress) = &self.progress {
                    progress.note_push();
                }
                Ok(())
            }
            Err(item) => {
                // ArrayQueue rejected the push (slot full). Roll back the
                // global reservation so the counter stays accurate.
                self.global_memory.sub(bytes);
                Err(item)
            }
        }
    }

    /// Push an item, retrying on a full queue until success or cancellation.
    ///
    /// This is the correct pattern for sources and forwarder threads that do
    /// not have a higher-level loop that already re-checks cancellation each
    /// iteration. Call this instead of hand-rolling `loop { match push ... }`.
    ///
    /// Returns `Ok(())` if the push succeeded, or `Err(item)` if `cancel`
    /// fired before the push could complete (the caller owns the item back).
    ///
    /// # Errors
    ///
    /// Returns the item back to the caller if cancellation was observed
    /// before the queue accepted the item.
    pub fn push_until_cancelled(
        &self,
        mut item: SequencedItem<T>,
        cancel: &super::cancel::CancelToken,
    ) -> Result<(), SequencedItem<T>> {
        let mut backoff = super::backoff::Backoff::new();
        loop {
            if cancel.is_cancelled() {
                return Err(item);
            }
            match self.push(item) {
                Ok(()) => return Ok(()),
                Err(returned) => {
                    item = returned;
                    backoff.snooze();
                }
            }
        }
    }

    /// Attempt to pop an item. Returns `None` if the queue is empty.
    #[must_use]
    pub fn pop(&self) -> Option<SequencedItem<T>> {
        let item = self.inner.pop()?;
        let bytes = item.memory_estimate;
        // Saturating CAS loop: under heavy concurrency `push` may not have
        // fetch-added its bytes yet by the time we observe the popped item,
        // so a naive `fetch_sub` could wrap around. We keep the visible
        // counter at `>= 0` so racing pushes don't panic on overflow.
        let mut current = self.queue_bytes.load(Ordering::Relaxed);
        loop {
            let next = current.saturating_sub(bytes);
            match self.queue_bytes.compare_exchange_weak(
                current,
                next,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(observed) => current = observed,
            }
        }
        self.global_memory.sub(bytes);
        if let Some(progress) = &self.progress {
            progress.note_pop();
        }
        Some(item)
    }

    /// Number of items currently queued.
    #[must_use]
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Whether the queue is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Total capacity (slot count).
    #[must_use]
    pub fn capacity(&self) -> usize {
        self.inner.capacity()
    }

    /// Fill ratio by slot count (0.0 to 1.0).
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn slot_fill_ratio(&self) -> f64 {
        self.inner.len() as f64 / self.inner.capacity() as f64
    }

    /// Fill ratio by memory (0.0 to 1.0+).
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn memory_fill_ratio(&self) -> f64 {
        self.queue_bytes.load(Ordering::Relaxed) as f64 / self.queue_limit_bytes as f64
    }

    /// Mark the queue as closed (no more items will be pushed).
    pub fn close(&self) {
        self.closed.store(true, Ordering::Release);
    }

    /// Whether the queue has been closed.
    #[must_use]
    pub fn is_closed(&self) -> bool {
        self.closed.load(Ordering::Acquire)
    }

    /// Whether the queue is closed and empty (i.e., fully drained).
    #[must_use]
    pub fn is_drained(&self) -> bool {
        self.is_closed() && self.is_empty()
    }

    /// Stage queue name (for logs/diagnostics).
    #[must_use]
    pub fn queue_name(&self) -> &str {
        &self.name
    }

    /// Snapshot this queue as a [`QueueSummary`] for the scheduler / backpressure
    /// sampling. `id` is assigned by the caller (typically the queue index in
    /// the pipeline's queue vector).
    #[must_use]
    pub fn summarize(&self, id: usize) -> QueueSummary {
        QueueSummary {
            id,
            slot_fill: self.slot_fill_ratio(),
            memory_fill: self.memory_fill_ratio(),
            closed: self.is_closed(),
            empty: self.is_empty(),
        }
    }
}

impl<T> Clone for StageQueue<T> {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
            queue_bytes: self.queue_bytes.clone(),
            queue_limit_bytes: self.queue_limit_bytes,
            global_memory: self.global_memory.clone(),
            closed: self.closed.clone(),
            name: self.name.clone(),
            progress: self.progress.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_tracker() -> Arc<MemoryTracker> {
        Arc::new(MemoryTracker::new(1_000_000))
    }

    #[test]
    fn test_push_pop_fifo() {
        let q: StageQueue<u32> = StageQueue::new("test", 4, 1_000, make_tracker());
        q.push(SequencedItem::new(0, 10, 4)).unwrap();
        q.push(SequencedItem::new(1, 20, 4)).unwrap();
        assert_eq!(q.pop().unwrap().item, 10);
        assert_eq!(q.pop().unwrap().item, 20);
        assert!(q.pop().is_none());
    }

    #[test]
    fn test_push_until_cancelled_succeeds_when_space_available() {
        use super::super::cancel::CancelToken;
        let q: StageQueue<u32> = StageQueue::new("test", 4, 1_000, make_tracker());
        let cancel = CancelToken::new();
        q.push_until_cancelled(SequencedItem::new(0, 42, 4), &cancel).unwrap();
        assert_eq!(q.pop().unwrap().item, 42);
    }

    #[test]
    fn test_push_until_cancelled_returns_item_when_cancelled() {
        use super::super::cancel::CancelToken;
        let q: StageQueue<u32> = StageQueue::new("test", 1, 1_000, make_tracker());
        let cancel = CancelToken::new();
        // Fill the queue so the next push must retry.
        q.push(SequencedItem::new(0, 1, 4)).unwrap();
        // Cancel immediately so the retry loop exits on its first check.
        cancel.cancel();
        let result = q.push_until_cancelled(SequencedItem::new(1, 99, 4), &cancel);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err().item, 99);
    }

    #[test]
    fn test_push_until_cancelled_exits_on_late_cancellation() {
        use super::super::cancel::CancelToken;
        use std::sync::Arc;
        let q: Arc<StageQueue<u32>> = Arc::new(StageQueue::new("test", 1, 1_000, make_tracker()));
        let cancel = CancelToken::new();

        q.push(SequencedItem::new(0, 1, 4)).unwrap();

        let cancel_for_timer = cancel.clone();
        std::thread::spawn(move || {
            std::thread::sleep(std::time::Duration::from_millis(50));
            cancel_for_timer.cancel();
        });

        let result = q.push_until_cancelled(SequencedItem::new(1, 99, 4), &cancel);
        assert!(result.is_err(), "expected retry loop to exit on cancel");
    }

    #[test]
    fn test_push_returns_item_on_slot_full() {
        let q: StageQueue<u32> = StageQueue::new("test", 2, 1_000, make_tracker());
        q.push(SequencedItem::new(0, 1, 4)).unwrap();
        q.push(SequencedItem::new(1, 2, 4)).unwrap();
        let rejected = q.push(SequencedItem::new(2, 3, 4)).unwrap_err();
        assert_eq!(rejected.item, 3);
    }

    #[test]
    fn test_push_returns_item_on_memory_full() {
        let q: StageQueue<u32> = StageQueue::new("test", 100, 10, make_tracker());
        q.push(SequencedItem::new(0, 1, 4)).unwrap();
        q.push(SequencedItem::new(1, 2, 4)).unwrap();
        // 4 + 4 = 8 used, next push of 4 bytes would push total to 12 > 10.
        let rejected = q.push(SequencedItem::new(2, 3, 4)).unwrap_err();
        assert_eq!(rejected.item, 3);
    }

    #[test]
    fn test_push_returns_item_on_global_memory_full() {
        // Tight global tracker (10 bytes) shared between two queues with
        // generous per-queue limits. A push that would exceed the global
        // limit must fail without mutating the queue or the tracker.
        let tracker = Arc::new(MemoryTracker::new(10));
        let q1: StageQueue<u32> = StageQueue::new("q1", 100, 1_000, tracker.clone());
        let q2: StageQueue<u32> = StageQueue::new("q2", 100, 1_000, tracker.clone());

        q1.push(SequencedItem::new(0, 1, 4)).unwrap();
        q2.push(SequencedItem::new(1, 2, 4)).unwrap();
        // Global is at 8; a 4-byte push would push to 12 > 10.
        let rejected = q1.push(SequencedItem::new(2, 3, 4)).unwrap_err();
        assert_eq!(rejected.item, 3);
        assert_eq!(tracker.current(), 8);
    }

    #[test]
    fn test_memory_tracking_through_push_pop() {
        let tracker = make_tracker();
        let q: StageQueue<u32> = StageQueue::new("test", 100, 1_000, tracker.clone());
        q.push(SequencedItem::new(0, 1, 100)).unwrap();
        assert_eq!(tracker.current(), 100);
        q.push(SequencedItem::new(1, 2, 200)).unwrap();
        assert_eq!(tracker.current(), 300);
        let _ = q.pop();
        assert_eq!(tracker.current(), 200);
        let _ = q.pop();
        assert_eq!(tracker.current(), 0);
    }

    #[test]
    fn test_close_and_drained() {
        let q: StageQueue<u32> = StageQueue::new("test", 4, 1_000, make_tracker());
        q.push(SequencedItem::new(0, 1, 4)).unwrap();
        assert!(!q.is_drained());
        q.close();
        assert!(q.is_closed());
        assert!(!q.is_drained()); // still has items
        let _ = q.pop();
        assert!(q.is_drained());
    }

    #[test]
    fn test_clone_shares_state() {
        let q1: StageQueue<u32> = StageQueue::new("test", 4, 1_000, make_tracker());
        let q2 = q1.clone();
        q1.push(SequencedItem::new(0, 42, 4)).unwrap();
        assert_eq!(q2.pop().unwrap().item, 42);
    }

    #[test]
    fn test_fill_ratios() {
        let q: StageQueue<u32> = StageQueue::new("test", 4, 1_000, make_tracker());
        q.push(SequencedItem::new(0, 1, 100)).unwrap();
        q.push(SequencedItem::new(1, 2, 200)).unwrap();
        assert!((q.slot_fill_ratio() - 0.5).abs() < 1e-9); // 2 of 4 slots
        assert!((q.memory_fill_ratio() - 0.3).abs() < 1e-9); // 300 of 1000 bytes
    }

    #[test]
    fn test_summarize_empty_queue() {
        let q: StageQueue<u32> = StageQueue::new("test", 4, 1_000, make_tracker());
        let s = q.summarize(7);
        assert_eq!(s.id, 7);
        assert!((s.slot_fill - 0.0).abs() < 1e-9);
        assert!((s.memory_fill - 0.0).abs() < 1e-9);
        assert!(!s.closed);
        assert!(s.empty);
    }

    #[test]
    fn test_summarize_partially_filled_queue() {
        let q: StageQueue<u32> = StageQueue::new("test", 4, 1_000, make_tracker());
        q.push(SequencedItem::new(0, 1, 250)).unwrap();
        let s = q.summarize(3);
        assert_eq!(s.id, 3);
        assert!((s.slot_fill - 0.25).abs() < 1e-9); // 1/4 slots
        assert!((s.memory_fill - 0.25).abs() < 1e-9); // 250/1000 bytes
        assert!(!s.closed);
        assert!(!s.empty);
    }

    #[test]
    fn test_summarize_closed_queue() {
        let q: StageQueue<u32> = StageQueue::new("test", 4, 1_000, make_tracker());
        q.close();
        let s = q.summarize(0);
        assert!(s.closed);
        assert!(s.empty);
    }

    /// Memory-accounting audit (task M3.8).
    ///
    /// Constructs a queue whose payload is a heap-allocated `Vec<u8>` of a
    /// known size, pushes N batches with a `memory_estimate` equal to the
    /// vector's heap footprint (`capacity()`), and proves that:
    ///
    /// 1. The global tracker's `current` and `peak` both report exactly
    ///    `N * batch.capacity()` after N pushes — i.e. the accounting matches
    ///    the heap bytes the caller declared were allocated.
    /// 2. One more push past the configured global limit is rejected without
    ///    mutating any counter (backpressure kicks in).
    /// 3. Popping releases the bytes back to the tracker, so a steady-state
    ///    push/pop loop never climbs above the limit.
    ///
    /// This guards against regressions where the `memory_estimate` would be
    /// silently swapped for `size_of::<T>()` (which would under-count heap
    /// payloads by 1000x for the `Vec<u8>` case) or where `try_add` would
    /// stop enforcing the global ceiling.
    #[test]
    fn test_accounting_matches_heap_bytes_and_enforces_global_limit() {
        const BATCH_BYTES: usize = 1024 * 1024; // 1 MiB
        const N: usize = 4;

        // Global limit is exactly N * BATCH_BYTES, so the (N+1)-th push must
        // be rejected by the global gate (not the per-queue gate).
        let tracker = Arc::new(MemoryTracker::new(N * BATCH_BYTES));
        // Per-queue limit is generous — we want to prove the GLOBAL gate fires.
        // Queue capacity is N + 1 slots so the slot gate doesn't trip first.
        let q: StageQueue<Vec<u8>> =
            StageQueue::new("audit", N + 1, 10 * N * BATCH_BYTES, tracker.clone());

        // A freshly-allocated `Vec<u8>` of length L has capacity >= L. We
        // advertise `capacity()` (what a real stage would report). Zero the
        // buffer so any allocator that defers page-fault-in is forced to
        // back the allocation.
        let make_batch = || {
            let mut v = vec![0u8; BATCH_BYTES];
            v.shrink_to_fit();
            v
        };

        // Sanity: each batch's declared memory matches its heap capacity.
        let probe = make_batch();
        assert_eq!(probe.capacity(), BATCH_BYTES);
        drop(probe);

        // 1. Accounting matches N * batch.capacity() after N pushes.
        for i in 0..N {
            let payload = make_batch();
            let bytes = payload.capacity();
            assert_eq!(bytes, BATCH_BYTES, "batch capacity must be stable");
            q.push(SequencedItem::new(i as u64, payload, bytes))
                .unwrap_or_else(|_| panic!("push {i} unexpectedly failed"));
            assert_eq!(
                tracker.current(),
                (i + 1) * BATCH_BYTES,
                "after push {i} tracker must equal cumulative heap bytes",
            );
        }
        assert_eq!(tracker.peak(), N * BATCH_BYTES, "peak tracks cumulative heap usage");

        // 2. One more push trips the global limit and is rejected cleanly.
        let overflow = make_batch();
        let overflow_bytes = overflow.capacity();
        let rejected = q
            .push(SequencedItem::new(N as u64, overflow, overflow_bytes))
            .expect_err("push past limit must be rejected");
        assert_eq!(
            rejected.item.capacity(),
            BATCH_BYTES,
            "rejected payload is returned to caller untouched",
        );
        assert_eq!(
            tracker.current(),
            N * BATCH_BYTES,
            "rejected push must not mutate the global tracker",
        );
        assert_eq!(tracker.peak(), N * BATCH_BYTES, "rejected push must not mutate peak");

        // 3. Popping an item releases its bytes and frees space for exactly
        //    one more push. Steady-state push/pop never exceeds the limit.
        let popped = q.pop().expect("queue is non-empty");
        assert_eq!(popped.item.capacity(), BATCH_BYTES);
        assert_eq!(
            tracker.current(),
            (N - 1) * BATCH_BYTES,
            "pop must release exactly the bytes that were reserved",
        );

        let replacement = make_batch();
        let replacement_bytes = replacement.capacity();
        q.push(SequencedItem::new(N as u64 + 1, replacement, replacement_bytes))
            .expect("after a pop, there is room for exactly one more push");
        assert_eq!(tracker.current(), N * BATCH_BYTES);
        assert_eq!(tracker.peak(), N * BATCH_BYTES, "steady-state loop never exceeds the limit");
    }
}
