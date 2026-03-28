//! Bounded, memory-aware work queues for connecting Zone 3 pipeline stages.
//!
//! Uses `crossbeam_queue::ArrayQueue` for lock-free, non-blocking push/pop.
//! Workers that find a queue full or empty simply try another stage (work-stealing).

use std::sync::{
    Arc,
    atomic::{AtomicBool, AtomicUsize, Ordering},
};

use crossbeam_queue::ArrayQueue;

/// A lock-free bounded queue with memory tracking for backpressure.
///
/// Uses `crossbeam_queue::ArrayQueue` — never blocks. `push` returns the item
/// on failure (full), `try_pop` returns `None` on empty. Workers that find a
/// queue full/empty simply try another stage (work-stealing).
///
/// Multiple clones of a `WorkQueue` share the same underlying queue, memory counter,
/// and completion flag, allowing work-stealing schedulers to have many producers and
/// consumers on the same queue.
///
/// The queue does not automatically decrement memory when items are popped; callers must
/// invoke [`WorkQueue::release_memory`] after consuming and processing each item.
pub struct WorkQueue<T> {
    inner: Arc<ArrayQueue<T>>,
    memory_bytes: Arc<AtomicUsize>,
    memory_limit: usize,
    name: &'static str,
    /// Explicit "input exhausted" flag — set by the producer when no more items will arrive.
    closed: Arc<AtomicBool>,
    /// Number of items that have been popped but whose results have not yet been pushed
    /// to the downstream queue. This prevents premature close cascades: a queue is only
    /// truly done when `closed && is_empty() && in_flight == 0`.
    in_flight: Arc<AtomicUsize>,
}

impl<T> WorkQueue<T> {
    /// Create a new work queue with the given capacity and memory limit.
    ///
    /// # Arguments
    ///
    /// * `capacity` — maximum number of items that may be buffered in the queue.
    /// * `memory_limit` — soft memory ceiling in bytes; checked by [`WorkQueue::memory_available`].
    /// * `name` — a static label used for logging and diagnostics.
    #[must_use]
    pub fn new(capacity: usize, memory_limit: usize, name: &'static str) -> Self {
        Self {
            inner: Arc::new(ArrayQueue::new(capacity)),
            memory_bytes: Arc::new(AtomicUsize::new(0)),
            memory_limit,
            name,
            closed: Arc::new(AtomicBool::new(false)),
            in_flight: Arc::new(AtomicUsize::new(0)),
        }
    }

    /// Attempt a non-blocking push.
    ///
    /// Returns `None` if the item was queued successfully, or `Some(item)` if the
    /// queue was full (the item is returned to the caller for holding).
    /// Memory is only incremented on success.
    pub fn push(&self, item: T, size_bytes: usize) -> Option<T> {
        match self.inner.push(item) {
            Ok(()) => {
                self.memory_bytes.fetch_add(size_bytes, Ordering::Relaxed);
                None
            }
            Err(returned) => Some(returned),
        }
    }

    /// Push an item onto the queue, spinning with backoff if full.
    ///
    /// `size_bytes` is added to the shared memory counter on success so that
    /// any concurrent reader of [`WorkQueue::memory_bytes`] sees the item's footprint
    /// as soon as it enters the queue.
    ///
    /// This spins with exponential backoff when the queue is full: the first few
    /// attempts yield the thread, then subsequent attempts sleep for 1 µs … 128 µs.
    /// The output queues are sized large enough (4096 slots) that persistent fullness
    /// is rare — this naturally applies backpressure as downstream stages drain.
    ///
    /// Use this for producer threads that must block when the queue is full (e.g. the
    /// position batcher feeding the pipeline entry point).
    ///
    /// # Errors
    ///
    /// Returns an error if the queue has been closed (no more items should be pushed).
    pub fn push_blocking(&self, item: T, size_bytes: usize) -> anyhow::Result<()> {
        let mut current = item;
        let mut attempts: u32 = 0;
        loop {
            match self.inner.push(current) {
                Ok(()) => {
                    self.memory_bytes.fetch_add(size_bytes, Ordering::Relaxed);
                    return Ok(());
                }
                Err(returned) => {
                    if self.closed.load(Ordering::Acquire) {
                        anyhow::bail!("queue '{}' is closed", self.name);
                    }
                    current = returned;
                    attempts += 1;
                    if attempts < 4 {
                        std::thread::yield_now();
                    } else {
                        // Exponential backoff: 1, 2, 4, 8, … up to 128 µs.
                        let shift = (attempts - 4).min(7);
                        std::thread::sleep(std::time::Duration::from_micros(1 << shift));
                    }
                }
            }
        }
    }

    /// Blocking pop. Spins with `thread::yield_now()` until an item is available
    /// or the queue is closed, empty, and no items are in flight.
    ///
    /// Returns `None` when the queue is fully done (no more items will arrive).
    ///
    /// On success, increments the in-flight counter. The caller **must** call
    /// [`WorkQueue::complete_in_flight`] after the popped item's result has been
    /// pushed to the downstream queue (or otherwise fully handled).
    #[must_use]
    pub fn pop(&self) -> Option<T> {
        loop {
            if let Some(item) = self.inner.pop() {
                self.in_flight.fetch_add(1, Ordering::AcqRel);
                return Some(item);
            }
            if self.is_closed_and_empty() {
                return None;
            }
            std::thread::yield_now();
        }
    }

    /// Non-blocking pop. Returns `None` if the queue is currently empty.
    ///
    /// On success, increments the in-flight counter. The caller **must** call
    /// [`WorkQueue::complete_in_flight`] after the popped item's result has been
    /// pushed to the downstream queue (or otherwise fully handled).
    #[must_use]
    pub fn try_pop(&self) -> Option<T> {
        let item = self.inner.pop()?;
        self.in_flight.fetch_add(1, Ordering::AcqRel);
        Some(item)
    }

    /// Decrement the memory counter by `size_bytes`.
    ///
    /// Callers must invoke this after fully processing each popped item so that the
    /// memory accounting stays accurate and [`WorkQueue::memory_available`] remains
    /// meaningful for backpressure decisions.
    pub fn release_memory(&self, size_bytes: usize) {
        self.memory_bytes.fetch_sub(size_bytes, Ordering::Relaxed);
    }

    /// Return the current approximate memory usage tracked by this queue, in bytes.
    #[must_use]
    pub fn memory_bytes(&self) -> usize {
        self.memory_bytes.load(Ordering::Relaxed)
    }

    /// Return `true` if the tracked memory usage is below the configured limit.
    ///
    /// This is a soft check; the queue does not enforce the limit — callers should
    /// consult this before deciding whether to push additional work.
    #[must_use]
    pub fn memory_available(&self) -> bool {
        self.memory_bytes() < self.memory_limit
    }

    /// Return the queue's slot capacity.
    #[must_use]
    pub fn capacity(&self) -> usize {
        self.inner.capacity()
    }

    /// Return the fill ratio (0.0 to 1.0) for backpressure decisions.
    #[expect(
        clippy::cast_precision_loss,
        reason = "queue lengths are small enough that f32 precision is sufficient"
    )]
    #[must_use]
    pub fn fill_ratio(&self) -> f32 {
        self.len() as f32 / self.capacity() as f32
    }

    /// Return the ratio of current memory usage to the configured memory limit.
    ///
    /// A value >= 1.0 means memory usage meets or exceeds the limit.
    #[expect(
        clippy::cast_precision_loss,
        reason = "memory values are approximate for backpressure decisions"
    )]
    #[must_use]
    pub fn memory_ratio(&self) -> f32 {
        self.memory_bytes() as f32 / self.memory_limit as f32
    }

    /// Return `true` if the underlying queue currently holds no items.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Return the number of items currently buffered in the queue.
    #[must_use]
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Return the static name assigned to this queue.
    #[must_use]
    pub fn name(&self) -> &'static str {
        self.name
    }

    /// Signal that one in-flight item has been fully processed and its result
    /// pushed to the downstream queue.
    ///
    /// Must be called exactly once per successful [`try_pop`](WorkQueue::try_pop) or
    /// [`pop`](WorkQueue::pop) call, after the downstream push is complete.
    pub fn complete_in_flight(&self) {
        self.in_flight.fetch_sub(1, Ordering::AcqRel);
    }

    /// Mark the queue as closed. No more items should be pushed after this call.
    ///
    /// Existing items can still be popped. This is used to signal that the upstream
    /// producer has finished and no further work will arrive.
    pub fn close(&self) {
        self.closed.store(true, Ordering::Release);
    }

    /// Returns `true` if the queue has been closed, is currently empty, **and** no
    /// items are in flight (being processed by workers between pop and downstream push).
    ///
    /// This three-way check prevents premature close cascades: a queue is only truly
    /// done when all popped items have been fully processed and their results pushed
    /// downstream.
    #[must_use]
    pub fn is_closed_and_empty(&self) -> bool {
        self.closed.load(Ordering::Acquire)
            && self.inner.is_empty()
            && self.in_flight.load(Ordering::Acquire) == 0
    }
}

impl<T> Clone for WorkQueue<T> {
    /// Clone the queue, sharing the same underlying queue and memory counter.
    ///
    /// All clones participate in the same queue; a push on any clone is
    /// visible to a pop on any other clone.
    fn clone(&self) -> Self {
        Self {
            inner: Arc::clone(&self.inner),
            memory_bytes: Arc::clone(&self.memory_bytes),
            memory_limit: self.memory_limit,
            name: self.name,
            closed: Arc::clone(&self.closed),
            in_flight: Arc::clone(&self.in_flight),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_push_pop() {
        let q: WorkQueue<u32> = WorkQueue::new(4, 1024, "test");
        q.push_blocking(42, 4).unwrap();
        let item = q.pop().unwrap();
        assert_eq!(item, 42);
    }

    #[test]
    fn test_memory_tracking() {
        let q: WorkQueue<u32> = WorkQueue::new(8, usize::MAX, "mem-test");

        q.push_blocking(1, 100).unwrap();
        q.push_blocking(2, 200).unwrap();
        assert_eq!(q.memory_bytes(), 300);

        // Pop and release the first item.
        let _ = q.pop().unwrap();
        q.release_memory(100);
        assert_eq!(q.memory_bytes(), 200);

        // Pop and release the second item.
        let _ = q.pop().unwrap();
        q.release_memory(200);
        assert_eq!(q.memory_bytes(), 0);
    }

    #[test]
    fn test_try_pop_empty() {
        let q: WorkQueue<u32> = WorkQueue::new(4, 1024, "empty-test");
        assert!(q.try_pop().is_none());
    }

    #[test]
    fn test_bounded_capacity() {
        let q: WorkQueue<u32> = WorkQueue::new(2, usize::MAX, "bounded-test");

        // Fill to capacity.
        assert!(q.push(1, 0).is_none());
        assert!(q.push(2, 0).is_none());

        // Queue is full — push should return Some(item).
        let result = q.push(3, 0);
        assert_eq!(result, Some(3), "expected Some(3) when queue is full");

        // Drain and confirm items arrive in order.
        assert_eq!(q.pop().unwrap(), 1);
        assert_eq!(q.pop().unwrap(), 2);
    }

    #[test]
    fn test_memory_available() {
        let q: WorkQueue<u32> = WorkQueue::new(8, 500, "avail-test");
        assert!(q.memory_available());

        q.push_blocking(1, 400).unwrap();
        assert!(q.memory_available()); // 400 < 500

        q.push_blocking(2, 101).unwrap();
        assert!(!q.memory_available()); // 501 >= 500

        let _ = q.pop().unwrap();
        q.release_memory(400);
        assert!(q.memory_available()); // back to 101 < 500
    }

    #[test]
    fn test_clone_shares_channel() {
        let q1: WorkQueue<u32> = WorkQueue::new(4, 1024, "clone-test");
        let q2 = q1.clone();

        q1.push_blocking(99, 8).unwrap();
        assert_eq!(q2.memory_bytes(), 8);

        let item = q2.pop().unwrap();
        assert_eq!(item, 99);

        q2.release_memory(8);
        assert_eq!(q1.memory_bytes(), 0);
    }

    #[test]
    fn test_close_and_is_closed_and_empty() {
        let q: WorkQueue<u32> = WorkQueue::new(8, usize::MAX, "close-test");

        // Not closed initially.
        assert!(!q.is_closed_and_empty());

        // Closed but not empty.
        q.push_blocking(1, 4).unwrap();
        q.close();
        assert!(!q.is_closed_and_empty());

        // Closed and empty after draining, but still in-flight.
        let _ = q.pop().unwrap();
        assert!(!q.is_closed_and_empty()); // in_flight == 1

        // Truly done after completing in-flight.
        q.complete_in_flight();
        assert!(q.is_closed_and_empty());
    }

    #[test]
    fn test_close_shared_across_clones() {
        let q1: WorkQueue<u32> = WorkQueue::new(8, usize::MAX, "clone-close-test");
        let q2 = q1.clone();

        q1.close();
        assert!(q2.is_closed_and_empty());
    }

    #[test]
    fn test_len_and_is_empty() {
        let q: WorkQueue<u32> = WorkQueue::new(8, usize::MAX, "len-test");
        assert!(q.is_empty());
        assert_eq!(q.len(), 0);

        q.push_blocking(1, 0).unwrap();
        q.push_blocking(2, 0).unwrap();
        assert!(!q.is_empty());
        assert_eq!(q.len(), 2);

        let _ = q.pop();
        assert_eq!(q.len(), 1);
    }

    #[test]
    fn test_pop_returns_none_when_closed_and_empty() {
        let q: WorkQueue<u32> = WorkQueue::new(8, usize::MAX, "pop-close-test");
        q.close();
        assert!(q.pop().is_none());
    }

    #[test]
    fn test_pop_drains_before_returning_none() {
        let q: WorkQueue<u32> = WorkQueue::new(8, usize::MAX, "pop-drain-test");
        q.push_blocking(1, 0).unwrap();
        q.push_blocking(2, 0).unwrap();
        q.close();

        assert_eq!(q.pop().unwrap(), 1);
        q.complete_in_flight();
        assert_eq!(q.pop().unwrap(), 2);
        q.complete_in_flight();
        assert!(q.pop().is_none());
    }
}
