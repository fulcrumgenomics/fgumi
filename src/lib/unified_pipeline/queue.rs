//! Memory-bounded queue types for pipeline flow control.
//!
//! This module provides queues that enforce memory limits rather than item counts,
//! enabling precise control over pipeline memory usage.
//!
//! # Key Types
//!
//! - [`MemoryBoundedQueue`]: A FIFO queue bounded by memory usage
//! - [`OrderedQueue`]: A reorder buffer with smart backpressure to prevent deadlock
//! - [`QueueStats`]: Statistics collected per queue for dynamic rebalancing

use crossbeam_queue::ArrayQueue;
use parking_lot::Mutex;
use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

/// Memory statistics for a queue
#[derive(Debug, Clone)]
pub struct QueueMemoryStats {
    pub current_bytes: u64,
    pub peak_bytes: u64,
    pub limit_bytes: u64,
}

/// Statistics collected per queue for rebalancing decisions.
#[derive(Debug, Clone, Default)]
pub struct QueueStats {
    /// Average memory usage over the epoch.
    pub avg_bytes: u64,
    /// Peak memory usage during the epoch.
    pub peak_bytes: u64,
    /// Total time blocked waiting to push (milliseconds).
    pub time_blocked_ms: u64,
}

/// A queue bounded by memory usage rather than item count.
///
/// The caller provides the heap size when pushing items. The queue tracks
/// total memory and refuses items when over the limit.
///
/// # Thread Safety
///
/// This queue is thread-safe and lock-free for the common path. Multiple
/// producers and consumers can operate concurrently.
///
/// # Example
///
/// ```ignore
/// let queue = MemoryBoundedQueue::new(1024, 1_000_000); // 1024 items max, 1MB limit
/// let item = vec![0u8; 1000];
/// let heap_size = item.capacity();
/// queue.push(item, heap_size)?;
/// let (item, size) = queue.pop().unwrap();
/// ```
pub struct MemoryBoundedQueue<T> {
    inner: ArrayQueue<(T, usize)>,
    current_bytes: AtomicU64,
    limit_bytes: AtomicU64, // Atomic to allow dynamic adjustment

    // Stats for rebalancing
    peak_bytes: AtomicU64,
    samples_sum: AtomicU64,
    samples_count: AtomicU64,
    blocked_ns: AtomicU64,
}

impl<T> MemoryBoundedQueue<T> {
    /// Create a new memory-bounded queue.
    ///
    /// # Arguments
    ///
    /// * `capacity` - Maximum number of items (for `ArrayQueue` backing store)
    /// * `limit_bytes` - Memory limit in bytes
    #[must_use]
    pub fn new(capacity: usize, limit_bytes: u64) -> Self {
        Self {
            inner: ArrayQueue::new(capacity),
            current_bytes: AtomicU64::new(0),
            limit_bytes: AtomicU64::new(limit_bytes),
            peak_bytes: AtomicU64::new(0),
            samples_sum: AtomicU64::new(0),
            samples_count: AtomicU64::new(0),
            blocked_ns: AtomicU64::new(0),
        }
    }

    /// Check if the queue can accept an item of the given size.
    ///
    /// This is a hint only - the actual push may still fail due to races.
    pub fn can_accept(&self, heap_size: usize) -> bool {
        let limit = self.limit_bytes.load(Ordering::Acquire);
        let current = self.current_bytes.load(Ordering::Acquire);
        current + heap_size as u64 <= limit
    }

    /// Push an item with its pre-computed heap size.
    ///
    /// Returns `Err((item, heap_size))` if over memory limit or queue is full.
    ///
    /// # Errors
    ///
    /// Returns the item and heap size if the queue is full or over its memory limit.
    pub fn push(&self, item: T, heap_size: usize) -> Result<(), (T, usize)> {
        let limit = self.limit_bytes.load(Ordering::Acquire);
        let current = self.current_bytes.load(Ordering::Acquire);

        if current + heap_size as u64 > limit {
            return Err((item, heap_size));
        }

        match self.inner.push((item, heap_size)) {
            Ok(()) => {
                let new_current = self.current_bytes.fetch_add(heap_size as u64, Ordering::AcqRel)
                    + heap_size as u64;

                // Update peak using CAS loop
                let mut peak = self.peak_bytes.load(Ordering::Relaxed);
                while new_current > peak {
                    match self.peak_bytes.compare_exchange_weak(
                        peak,
                        new_current,
                        Ordering::Relaxed,
                        Ordering::Relaxed,
                    ) {
                        Ok(_) => break,
                        Err(p) => peak = p,
                    }
                }

                Ok(())
            }
            Err((item, size)) => Err((item, size)),
        }
    }

    /// Pop an item, returning the item and its heap size.
    pub fn pop(&self) -> Option<(T, usize)> {
        self.inner.pop().map(|(item, size)| {
            self.current_bytes.fetch_sub(size as u64, Ordering::AcqRel);
            (item, size)
        })
    }

    /// Current memory usage in bytes.
    pub fn current_bytes(&self) -> u64 {
        self.current_bytes.load(Ordering::Acquire)
    }

    /// Get queue memory statistics for reporting.
    pub fn memory_stats(&self) -> QueueMemoryStats {
        QueueMemoryStats {
            current_bytes: self.current_bytes(),
            peak_bytes: self.peak_bytes.load(Ordering::Relaxed),
            limit_bytes: self.limit_bytes.load(Ordering::Relaxed),
        }
    }

    /// Check if queue is at or over memory limit.
    pub fn is_full(&self) -> bool {
        self.current_bytes() >= self.limit_bytes.load(Ordering::Acquire)
    }

    /// Update the memory limit (for dynamic rebalancing).
    pub fn set_limit(&self, new_limit: u64) {
        self.limit_bytes.store(new_limit, Ordering::Release);
    }

    /// Get current limit.
    pub fn limit_bytes(&self) -> u64 {
        self.limit_bytes.load(Ordering::Acquire)
    }

    /// Record a sample for stats (call periodically).
    pub fn record_sample(&self) {
        let current = self.current_bytes.load(Ordering::Relaxed);
        self.samples_sum.fetch_add(current, Ordering::Relaxed);
        self.samples_count.fetch_add(1, Ordering::Relaxed);
    }

    /// Record blocked time in nanoseconds.
    pub fn record_blocked(&self, ns: u64) {
        self.blocked_ns.fetch_add(ns, Ordering::Relaxed);
    }

    /// Collect and reset stats for this epoch.
    pub fn collect_stats(&self) -> QueueStats {
        let peak = self.peak_bytes.swap(0, Ordering::Relaxed);
        let sum = self.samples_sum.swap(0, Ordering::Relaxed);
        let count = self.samples_count.swap(0, Ordering::Relaxed);
        let blocked = self.blocked_ns.swap(0, Ordering::Relaxed);

        QueueStats {
            avg_bytes: if count > 0 { sum / count } else { 0 },
            peak_bytes: peak,
            time_blocked_ms: blocked / 1_000_000,
        }
    }

    /// Number of items in the queue.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Check if queue is empty.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

/// A reorder buffer that outputs items in serial order.
///
/// Uses smart backpressure to prevent deadlock:
/// - When waiting for `next_seq`: MUST accept items (refusing would deadlock)
/// - When we have `next_seq`: CAN refuse items (consumer can drain)
///
/// # Deadlock Prevention
///
/// The key insight is that if we're waiting for serial N, we must accept
/// serials N+1, N+2, etc. because serial N might be produced by another
/// thread that's blocked trying to push to this queue. Only when we have
/// serial N available can we safely apply backpressure, because the consumer
/// can make progress by draining serial N.
///
/// # Example
///
/// ```ignore
/// let queue = OrderedQueue::new(1_000_000); // 1MB limit
/// queue.insert(2, item2, 1000)?; // Accepted (waiting for 0)
/// queue.insert(0, item0, 1000)?; // Accepted (now has 0)
/// let (item, size) = queue.try_pop_next().unwrap(); // Returns item0
/// queue.insert(1, item1, 1000)?; // Accepted (now has 1)
/// ```
pub struct OrderedQueue<T> {
    inner: Mutex<OrderedQueueInner<T>>,
    current_bytes: AtomicU64,
    limit_bytes: AtomicU64,
    next_seq: AtomicU64,  // Cached for lock-free checks
    has_next: AtomicBool, // Cached: do we have next_seq?

    // Stats
    peak_bytes: AtomicU64,
    samples_sum: AtomicU64,
    samples_count: AtomicU64,
    blocked_ns: AtomicU64,
}

struct OrderedQueueInner<T> {
    buffer: HashMap<u64, (T, usize)>,
    next_seq: u64,
}

impl<T> OrderedQueue<T> {
    /// Create a new ordered queue with the given memory limit.
    #[must_use]
    pub fn new(limit_bytes: u64) -> Self {
        Self {
            inner: Mutex::new(OrderedQueueInner { buffer: HashMap::new(), next_seq: 0 }),
            current_bytes: AtomicU64::new(0),
            limit_bytes: AtomicU64::new(limit_bytes),
            next_seq: AtomicU64::new(0),
            has_next: AtomicBool::new(false),
            peak_bytes: AtomicU64::new(0),
            samples_sum: AtomicU64::new(0),
            samples_count: AtomicU64::new(0),
            blocked_ns: AtomicU64::new(0),
        }
    }

    /// Check if we can accept an item (lock-free fast path).
    ///
    /// Returns true if:
    /// - We don't have `next_seq` (must accept to make progress), OR
    /// - We're under the memory limit
    pub fn can_accept(&self, heap_size: usize) -> bool {
        // If we don't have next_seq, we MUST accept (deadlock avoidance)
        if !self.has_next.load(Ordering::Acquire) {
            return true;
        }

        // We have next_seq, so consumer can drain. Apply backpressure.
        let current = self.current_bytes.load(Ordering::Acquire);
        let limit = self.limit_bytes.load(Ordering::Acquire);
        current + heap_size as u64 <= limit
    }

    /// Insert an item into the reorder buffer.
    ///
    /// Acceptance rule:
    /// - If we do NOT have `next_seq`: ACCEPT (must accumulate for progress)
    /// - If we DO have `next_seq`: only accept if under memory limit
    ///
    /// Returns `Err((item, heap_size))` if rejected due to backpressure.
    ///
    /// # Errors
    ///
    /// Returns the item and heap size if rejected due to memory backpressure.
    pub fn insert(&self, serial: u64, item: T, heap_size: usize) -> Result<(), (T, usize)> {
        let mut inner = self.inner.lock();

        let has_next = inner.buffer.contains_key(&inner.next_seq);

        if has_next {
            // Consumer can drain - apply backpressure
            let current = self.current_bytes.load(Ordering::Acquire);
            let limit = self.limit_bytes.load(Ordering::Acquire);
            if current + heap_size as u64 > limit {
                return Err((item, heap_size));
            }
        }
        // else: must accept, we need to accumulate until next_seq arrives

        inner.buffer.insert(serial, (item, heap_size));
        let new_current =
            self.current_bytes.fetch_add(heap_size as u64, Ordering::AcqRel) + heap_size as u64;

        // Update cached state
        let new_has_next = inner.buffer.contains_key(&inner.next_seq);
        self.has_next.store(new_has_next, Ordering::Release);

        // Update peak using CAS loop
        let mut peak = self.peak_bytes.load(Ordering::Relaxed);
        while new_current > peak {
            match self.peak_bytes.compare_exchange_weak(
                peak,
                new_current,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(p) => peak = p,
            }
        }

        Ok(())
    }

    /// Try to pop the next item in serial order.
    ///
    /// Returns `Some((item, heap_size))` if `next_seq` is available.
    pub fn try_pop_next(&self) -> Option<(T, usize)> {
        let mut inner = self.inner.lock();

        let next = inner.next_seq;
        if let Some((item, heap_size)) = inner.buffer.remove(&next) {
            inner.next_seq += 1;
            self.current_bytes.fetch_sub(heap_size as u64, Ordering::AcqRel);

            // Update cached state
            self.next_seq.store(inner.next_seq, Ordering::Release);
            let new_has_next = inner.buffer.contains_key(&inner.next_seq);
            self.has_next.store(new_has_next, Ordering::Release);

            Some((item, heap_size))
        } else {
            None
        }
    }

    /// Get the next expected serial number.
    pub fn next_seq(&self) -> u64 {
        self.next_seq.load(Ordering::Acquire)
    }

    /// Check if we have the next expected serial (can make progress).
    pub fn can_pop(&self) -> bool {
        self.has_next.load(Ordering::Acquire)
    }

    /// Current memory usage in bytes.
    pub fn current_bytes(&self) -> u64 {
        self.current_bytes.load(Ordering::Acquire)
    }

    /// Update the memory limit (for dynamic rebalancing).
    pub fn set_limit(&self, new_limit: u64) {
        self.limit_bytes.store(new_limit, Ordering::Release);
    }

    /// Get current limit.
    pub fn limit_bytes(&self) -> u64 {
        self.limit_bytes.load(Ordering::Acquire)
    }

    /// Number of items in the buffer.
    pub fn len(&self) -> usize {
        self.inner.lock().buffer.len()
    }

    /// Check if buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.inner.lock().buffer.is_empty()
    }

    /// Record a sample for stats.
    pub fn record_sample(&self) {
        let current = self.current_bytes.load(Ordering::Relaxed);
        self.samples_sum.fetch_add(current, Ordering::Relaxed);
        self.samples_count.fetch_add(1, Ordering::Relaxed);
    }

    /// Record blocked time in nanoseconds.
    pub fn record_blocked(&self, ns: u64) {
        self.blocked_ns.fetch_add(ns, Ordering::Relaxed);
    }

    /// Collect and reset stats.
    pub fn collect_stats(&self) -> QueueStats {
        let peak = self.peak_bytes.swap(0, Ordering::Relaxed);
        let sum = self.samples_sum.swap(0, Ordering::Relaxed);
        let count = self.samples_count.swap(0, Ordering::Relaxed);
        let blocked = self.blocked_ns.swap(0, Ordering::Relaxed);

        QueueStats {
            avg_bytes: if count > 0 { sum / count } else { 0 },
            peak_bytes: peak,
            time_blocked_ms: blocked / 1_000_000,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_memory_bounded_queue_basic() {
        let queue: MemoryBoundedQueue<Vec<u8>> = MemoryBoundedQueue::new(100, 1000);

        // Push items within limit
        let item1 = vec![0u8; 100];
        assert!(queue.push(item1, 100).is_ok());
        assert_eq!(queue.current_bytes(), 100);

        let item2 = vec![0u8; 200];
        assert!(queue.push(item2, 200).is_ok());
        assert_eq!(queue.current_bytes(), 300);

        // Pop and verify accounting
        let (_, size) = queue.pop().unwrap();
        assert_eq!(size, 100);
        assert_eq!(queue.current_bytes(), 200);
    }

    #[test]
    fn test_memory_bounded_queue_limit() {
        let queue: MemoryBoundedQueue<Vec<u8>> = MemoryBoundedQueue::new(100, 500);

        // Fill to limit
        for _ in 0..5 {
            let item = vec![0u8; 100];
            assert!(queue.push(item, 100).is_ok());
        }
        assert_eq!(queue.current_bytes(), 500);

        // Should reject
        let item = vec![0u8; 100];
        assert!(queue.push(item, 100).is_err());
    }

    #[test]
    fn test_ordered_queue_basic() {
        let queue: OrderedQueue<u32> = OrderedQueue::new(1000);

        // Insert out of order
        assert!(queue.insert(2, 200, 10).is_ok());
        assert!(queue.insert(0, 100, 10).is_ok());
        assert!(queue.insert(1, 150, 10).is_ok());

        // Pop in order
        let (val, _) = queue.try_pop_next().unwrap();
        assert_eq!(val, 100);
        let (val, _) = queue.try_pop_next().unwrap();
        assert_eq!(val, 150);
        let (val, _) = queue.try_pop_next().unwrap();
        assert_eq!(val, 200);

        assert!(queue.try_pop_next().is_none());
    }

    #[test]
    fn test_ordered_queue_backpressure_when_has_next() {
        let queue: OrderedQueue<u32> = OrderedQueue::new(100);

        // Insert serial 0 - now we have next_seq
        assert!(queue.insert(0, 100, 50).is_ok());
        assert!(queue.can_pop());

        // We have next_seq, so backpressure applies
        // Try to insert something that would exceed limit
        assert!(queue.insert(1, 200, 60).is_err());

        // But we can still insert if under limit
        assert!(queue.insert(1, 200, 40).is_ok());
    }

    #[test]
    fn test_ordered_queue_must_accept_when_waiting() {
        let queue: OrderedQueue<u32> = OrderedQueue::new(100);

        // No next_seq yet - must accept even if over limit
        assert!(queue.insert(5, 500, 200).is_ok()); // Way over 100 byte limit
        assert!(!queue.can_pop()); // Still waiting for serial 0

        // Still must accept because we don't have serial 0
        assert!(queue.insert(3, 300, 200).is_ok());
        assert!(queue.insert(1, 100, 200).is_ok());

        // Now insert serial 0
        assert!(queue.insert(0, 0, 10).is_ok());
        assert!(queue.can_pop()); // Now we have next_seq

        // NOW backpressure applies - should reject
        assert!(queue.insert(2, 200, 200).is_err());
    }

    #[test]
    fn test_stats_collection() {
        let queue: MemoryBoundedQueue<u32> = MemoryBoundedQueue::new(100, 1000);

        queue.push(1, 100).unwrap();
        queue.record_sample(); // Sample at 100 bytes
        queue.push(2, 200).unwrap();
        queue.record_sample(); // Sample at 300 bytes (100 + 200)
        queue.record_blocked(5_000_000); // 5ms

        let stats = queue.collect_stats();
        assert_eq!(stats.avg_bytes, 200); // (100 + 300) / 2
        assert_eq!(stats.peak_bytes, 300);
        assert_eq!(stats.time_blocked_ms, 5);

        // Stats should be reset
        let stats2 = queue.collect_stats();
        assert_eq!(stats2.avg_bytes, 0);
        assert_eq!(stats2.peak_bytes, 0);
        assert_eq!(stats2.time_blocked_ms, 0);
    }

    // =========================================================================
    // Backpressure Tests (Phase 4 of the plan)
    // =========================================================================

    /// Test that backpressure prevents memory exhaustion.
    ///
    /// This verifies the queue properly rejects items when memory limit is reached,
    /// preventing OOM conditions in the pipeline.
    #[test]
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    fn test_backpressure_prevents_oom() {
        // Create queue with small memory limit (1KB)
        let queue: MemoryBoundedQueue<Vec<u8>> = MemoryBoundedQueue::new(16, 1024);

        // Push items until we hit the limit
        let mut pushed = 0;
        let mut rejected = 0;
        let item_size = 200; // 200 bytes per item

        for i in 0..100 {
            let item = vec![i as u8; item_size];
            match queue.push(item, item_size) {
                Ok(()) => pushed += 1,
                Err(_) => rejected += 1,
            }
        }

        // We should have pushed some items
        assert!(pushed > 0, "Should accept some items");

        // We should have rejected items due to backpressure
        assert!(rejected > 0, "Should reject items when over memory limit");

        // Total memory should not exceed limit by much
        // (one item may slightly exceed the limit)
        let max_expected: u64 = 1024 + item_size as u64;
        assert!(
            queue.current_bytes() <= max_expected,
            "Memory should stay bounded: {} <= {}",
            queue.current_bytes(),
            max_expected
        );

        // Pop items and verify memory is freed
        while let Some((item, _)) = queue.pop() {
            drop(item);
        }
        assert_eq!(queue.current_bytes(), 0, "Queue should be empty");
    }

    /// Test that `OrderedQueue` backpressure respects memory limits.
    #[test]
    #[allow(clippy::cast_possible_truncation)]
    fn test_ordered_queue_backpressure_memory_bound() {
        // Small limit: 500 bytes
        let queue: OrderedQueue<Vec<u8>> = OrderedQueue::new(500);

        // Insert serial 0 first so backpressure can apply
        assert!(queue.insert(0, vec![0u8; 100], 100).is_ok());

        // Now try to insert items
        let mut pushed = 0;
        let mut rejected = 0;

        for i in 1..20 {
            let item = vec![i as u8; 100];
            match queue.insert(i, item, 100) {
                Ok(()) => pushed += 1,
                Err(_) => rejected += 1,
            }
        }

        // Should have some accepted and some rejected
        assert!(pushed > 0, "Should accept some items");
        assert!(rejected > 0, "Should reject items when over limit");

        // Verify we can drain the queue
        let mut count = 0;
        while queue.try_pop_next().is_some() {
            count += 1;
        }
        assert!(count > 0, "Should pop the items we inserted");
    }

    /// Test that backpressure allows recovery when items are consumed.
    #[test]
    fn test_backpressure_recovery() {
        let queue: MemoryBoundedQueue<u32> = MemoryBoundedQueue::new(8, 100);

        // Fill the queue
        for i in 0..5 {
            queue.push(i, 20).unwrap();
        }

        // Should reject new items
        assert!(queue.push(100, 20).is_err());

        // Pop some items to free space
        queue.pop();
        queue.pop();

        // Should now accept items again
        assert!(queue.push(100, 20).is_ok());
    }
}
