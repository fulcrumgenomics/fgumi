//! Reordering buffer for out-of-order batch completion.
//!
//! This module provides a buffer that accepts items tagged with sequence numbers
//! and releases them in sequential order. It's used by the ordered work-stealing
//! pipeline to maintain output order when batches complete out of order.
//!
//! # Example
//!
//! ```
//! use fgumi_lib::reorder_buffer::ReorderBuffer;
//!
//! let mut buffer: ReorderBuffer<String> = ReorderBuffer::new();
//!
//! // Insert items out of order
//! buffer.insert(2, "third".to_string());
//! buffer.insert(0, "first".to_string());
//! buffer.insert(1, "second".to_string());
//!
//! // Pop in sequence order
//! assert_eq!(buffer.try_pop_next(), Some("first".to_string()));
//! assert_eq!(buffer.try_pop_next(), Some("second".to_string()));
//! assert_eq!(buffer.try_pop_next(), Some("third".to_string()));
//! assert_eq!(buffer.try_pop_next(), None);
//! ```

use std::collections::VecDeque;

use crate::unified_pipeline::MemoryEstimate;

/// A buffer that releases items in sequential order.
///
/// Items can be inserted with any sequence number, but they are only
/// released when all prior sequence numbers have been released.
/// This enables ordered output from parallel processing where
/// items may complete out of order.
///
/// Uses a sparse `VecDeque` for O(1) insert and pop operations.
///
/// # Memory Tracking
///
/// The buffer optionally tracks heap memory usage when using the
/// `insert_with_size`/`try_pop_next_with_size` methods. This enables
/// memory-bounded backpressure without expensive per-item size calculations.
#[derive(Debug)]
pub struct ReorderBuffer<T> {
    /// Sparse buffer: index (seq - `base_seq`) maps to `Option<(T, usize)>` where
    /// usize is the pre-computed heap size (0 if not tracked).
    buffer: VecDeque<Option<(T, usize)>>,
    /// The sequence number corresponding to buffer[0].
    base_seq: u64,
    /// Next sequence number to release.
    next_seq: u64,
    /// Number of items currently stored.
    count: usize,
    /// Total tracked heap memory in bytes.
    heap_bytes: u64,
    /// Whether the next sequential item is ready to pop (optimization to avoid
    /// repeated checks). Updated on insert/pop.
    can_pop: bool,
}

impl<T> ReorderBuffer<T> {
    /// Create a new reorder buffer.
    #[must_use]
    pub fn new() -> Self {
        Self {
            buffer: VecDeque::new(),
            base_seq: 0,
            next_seq: 0,
            count: 0,
            heap_bytes: 0,
            can_pop: false,
        }
    }

    /// Insert an item with a sequence number (without memory tracking).
    ///
    /// Items can be inserted in any order. They will be released
    /// in sequence order via `try_pop_next()` or `drain_ready()`.
    ///
    /// # Arguments
    ///
    /// * `seq` - The sequence number for this item (0-indexed, monotonically increasing).
    /// * `item` - The item to buffer.
    ///
    /// # Panics
    ///
    /// Panics in debug mode if an item with the same sequence number is already buffered,
    /// or if the sequence number is before the current base.
    pub fn insert(&mut self, seq: u64, item: T) {
        self.insert_with_size(seq, item, 0);
    }

    /// Insert an item with a sequence number and pre-computed heap size.
    ///
    /// This variant tracks memory usage for backpressure control.
    /// The size should be computed once when creating the item.
    ///
    /// # Arguments
    ///
    /// * `seq` - The sequence number for this item (0-indexed, monotonically increasing).
    /// * `item` - The item to buffer.
    /// * `heap_size` - Pre-computed heap size of the item in bytes.
    ///
    /// # Panics
    ///
    /// Panics in debug mode if an item with the same sequence number is already buffered,
    /// or if the sequence number is before the current base.
    #[allow(clippy::cast_possible_truncation)]
    pub fn insert_with_size(&mut self, seq: u64, item: T, heap_size: usize) {
        debug_assert!(
            seq >= self.base_seq,
            "Sequence number {seq} is before base {}",
            self.base_seq
        );

        let index = (seq - self.base_seq) as usize;

        // Extend buffer with None entries if needed
        while self.buffer.len() <= index {
            self.buffer.push_back(None);
        }

        debug_assert!(self.buffer[index].is_none(), "Duplicate sequence number: {seq}");
        self.buffer[index] = Some((item, heap_size));
        self.count += 1;
        self.heap_bytes += heap_size as u64;

        // Update can_pop: if we inserted at index 0, we can now pop
        if index == 0 {
            self.can_pop = true;
        }
    }

    /// Pop the next sequential item if available (without returning size).
    ///
    /// Returns `Some(item)` if the item with `next_seq` is buffered,
    /// otherwise returns `None`.
    ///
    /// This advances the internal sequence counter, so subsequent calls
    /// will return the next item in sequence.
    #[must_use]
    pub fn try_pop_next(&mut self) -> Option<T> {
        self.try_pop_next_with_size().map(|(item, _size)| item)
    }

    /// Pop the next sequential item if available, returning the item and its tracked size.
    ///
    /// Returns `Some((item, heap_size))` if the item with `next_seq` is buffered,
    /// otherwise returns `None`.
    ///
    /// # Panics
    ///
    /// Panics if internal state is inconsistent (the front item was indicated as poppable but is missing).
    #[must_use]
    pub fn try_pop_next_with_size(&mut self) -> Option<(T, usize)> {
        // Check if front item is present (handles empty buffer and gaps)
        if !self.can_pop {
            return None;
        }

        // Pop the front item
        let (item, size) = self.buffer.pop_front().unwrap().unwrap();
        self.base_seq += 1;
        self.next_seq += 1;
        self.count -= 1;
        self.heap_bytes = self.heap_bytes.saturating_sub(size as u64);

        // Update can_pop for next item
        self.can_pop = self.buffer.front().is_some_and(Option::is_some);

        Some((item, size))
    }

    /// Drain all consecutive ready items starting from the current sequence.
    ///
    /// Returns an iterator that yields items in sequence order, stopping
    /// when it reaches a gap in the sequence.
    ///
    /// # Example
    ///
    /// ```
    /// use fgumi_lib::reorder_buffer::ReorderBuffer;
    ///
    /// let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();
    /// buffer.insert(0, 10);
    /// buffer.insert(1, 20);
    /// buffer.insert(3, 40);  // Gap at 2
    ///
    /// let ready: Vec<_> = buffer.drain_ready().collect();
    /// assert_eq!(ready, vec![10, 20]);  // Stops at gap
    /// ```
    pub fn drain_ready(&mut self) -> DrainReady<'_, T> {
        DrainReady { buffer: self }
    }

    /// Check if the buffer is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Get the number of items currently stored in the buffer.
    #[must_use]
    pub fn len(&self) -> usize {
        self.count
    }

    /// Get the next expected sequence number.
    #[must_use]
    pub fn next_seq(&self) -> u64 {
        self.next_seq
    }

    /// Check if the next sequential item is ready to pop.
    ///
    /// This is O(1) and can be called frequently for backpressure decisions.
    #[must_use]
    pub fn can_pop(&self) -> bool {
        self.can_pop
    }

    /// Get the tracked heap memory in bytes.
    ///
    /// This returns the sum of sizes passed to `insert_with_size`.
    /// It's O(1) and suitable for frequent backpressure checks.
    #[must_use]
    pub fn heap_bytes(&self) -> u64 {
        self.heap_bytes
    }

    /// Check if this reorder buffer would accept a new item with the given serial.
    ///
    /// This implements deadlock-free admission control:
    /// - Always accept if serial == `next_seq` (the item we're waiting for)
    /// - Always accept if we can't pop (we're stuck, need more items)
    /// - Reject if over the memory limit AND we can make progress
    ///
    /// # Arguments
    ///
    /// * `serial` - The sequence number of the item to potentially insert
    /// * `memory_limit` - The memory limit in bytes (0 = unlimited)
    ///
    /// # Returns
    ///
    /// `true` if the item should be accepted, `false` if backpressure should be applied.
    #[must_use]
    pub fn would_accept(&self, serial: u64, memory_limit: u64) -> bool {
        // No limit = always accept
        if memory_limit == 0 {
            return true;
        }

        // Always accept the item we need next (avoids deadlock)
        if serial == self.next_seq {
            return true;
        }

        // If we can't make progress, accept everything (avoids deadlock)
        if !self.can_pop {
            return true;
        }

        // We can make progress and this isn't the needed item.
        // Apply backpressure if over limit.
        self.heap_bytes < memory_limit
    }

    /// Calculate total heap memory used by all items in the buffer.
    ///
    /// This iterates all items and sums their heap sizes. Only use for
    /// profiling/debugging as it's O(n) in the number of buffered items.
    ///
    /// Note: This calculates fresh from items using `MemoryEstimate`, not from
    /// tracked sizes. Use `heap_bytes()` for the tracked value.
    #[must_use]
    pub fn total_heap_size(&self) -> usize
    where
        T: MemoryEstimate,
    {
        self.buffer
            .iter()
            .filter_map(|opt| opt.as_ref())
            .map(|(item, _size)| item.estimate_heap_size())
            .sum()
    }

    /// Set the starting sequence numbers (for testing or reset scenarios).
    #[cfg(test)]
    pub fn set_base_seq(&mut self, base: u64) {
        self.base_seq = base;
        self.next_seq = base;
    }
}

impl<T> Default for ReorderBuffer<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Iterator that drains consecutive ready items from a `ReorderBuffer`.
pub struct DrainReady<'a, T> {
    buffer: &'a mut ReorderBuffer<T>,
}

impl<T> Iterator for DrainReady<'_, T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.buffer.try_pop_next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_in_order_insertion() {
        let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();

        buffer.insert(0, 100);
        buffer.insert(1, 200);
        buffer.insert(2, 300);

        assert_eq!(buffer.try_pop_next(), Some(100));
        assert_eq!(buffer.try_pop_next(), Some(200));
        assert_eq!(buffer.try_pop_next(), Some(300));
        assert_eq!(buffer.try_pop_next(), None);
    }

    #[test]
    fn test_out_of_order_insertion() {
        let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();

        // Insert in reverse order
        buffer.insert(2, 300);
        buffer.insert(1, 200);
        buffer.insert(0, 100);

        // Should still pop in sequence order
        assert_eq!(buffer.try_pop_next(), Some(100));
        assert_eq!(buffer.try_pop_next(), Some(200));
        assert_eq!(buffer.try_pop_next(), Some(300));
        assert_eq!(buffer.try_pop_next(), None);
    }

    #[test]
    fn test_gap_blocks_progress() {
        let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();

        buffer.insert(0, 100);
        buffer.insert(2, 300); // Gap at 1

        assert_eq!(buffer.try_pop_next(), Some(100));
        assert_eq!(buffer.try_pop_next(), None); // Blocked on missing 1

        // Fill the gap
        buffer.insert(1, 200);

        assert_eq!(buffer.try_pop_next(), Some(200));
        assert_eq!(buffer.try_pop_next(), Some(300));
        assert_eq!(buffer.try_pop_next(), None);
    }

    #[test]
    fn test_drain_ready() {
        let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();

        buffer.insert(0, 100);
        buffer.insert(1, 200);
        buffer.insert(3, 400); // Gap at 2

        let ready: Vec<_> = buffer.drain_ready().collect();
        assert_eq!(ready, vec![100, 200]);

        // Buffer should have 3 still pending
        assert!(!buffer.is_empty());
        assert_eq!(buffer.next_seq(), 2);

        // Fill gap and drain again
        buffer.insert(2, 300);
        let more: Vec<_> = buffer.drain_ready().collect();
        assert_eq!(more, vec![300, 400]);
        assert!(buffer.is_empty());
    }

    #[test]
    fn test_sparse_insertion() {
        let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();

        // Insert with gaps
        buffer.insert(0, 100);
        buffer.insert(5, 600);
        buffer.insert(2, 300);

        assert_eq!(buffer.try_pop_next(), Some(100));
        assert_eq!(buffer.try_pop_next(), None); // Gap at 1

        buffer.insert(1, 200);
        assert_eq!(buffer.try_pop_next(), Some(200));
        assert_eq!(buffer.try_pop_next(), Some(300));
        assert_eq!(buffer.try_pop_next(), None); // Gap at 3, 4

        buffer.insert(3, 400);
        buffer.insert(4, 500);

        let rest: Vec<_> = buffer.drain_ready().collect();
        assert_eq!(rest, vec![400, 500, 600]);
    }

    #[test]
    fn test_large_sequence_numbers() {
        let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();

        // Start from a large sequence number
        let start = 1_000_000u64;
        buffer.set_base_seq(start);

        buffer.insert(start, 100);
        buffer.insert(start + 1, 200);

        assert_eq!(buffer.try_pop_next(), Some(100));
        assert_eq!(buffer.try_pop_next(), Some(200));
    }

    #[test]
    fn test_memory_tracking() {
        let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();

        // Insert with sizes
        buffer.insert_with_size(0, 100, 1000);
        buffer.insert_with_size(1, 200, 2000);
        buffer.insert_with_size(2, 300, 3000);

        assert_eq!(buffer.heap_bytes(), 6000);
        assert_eq!(buffer.len(), 3);

        // Pop and verify memory decreases
        let (val, size) = buffer.try_pop_next_with_size().unwrap();
        assert_eq!(val, 100);
        assert_eq!(size, 1000);
        assert_eq!(buffer.heap_bytes(), 5000);

        let (val, size) = buffer.try_pop_next_with_size().unwrap();
        assert_eq!(val, 200);
        assert_eq!(size, 2000);
        assert_eq!(buffer.heap_bytes(), 3000);
    }

    #[test]
    fn test_can_pop_tracking() {
        let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();

        // Empty buffer - can't pop
        assert!(!buffer.can_pop());

        // Insert at seq 1 (gap at 0) - still can't pop
        buffer.insert(1, 200);
        assert!(!buffer.can_pop());

        // Insert at seq 0 - now can pop
        buffer.insert(0, 100);
        assert!(buffer.can_pop());

        // Pop - should still be able to pop (seq 1 is ready)
        assert_eq!(buffer.try_pop_next(), Some(100));
        assert!(buffer.can_pop());

        // Pop again - now can't pop (empty)
        assert_eq!(buffer.try_pop_next(), Some(200));
        assert!(!buffer.can_pop());
    }

    #[test]
    fn test_would_accept_no_limit() {
        let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();
        buffer.insert_with_size(0, 100, 1_000_000_000); // 1GB

        // With no limit (0), always accept
        assert!(buffer.would_accept(1, 0));
        assert!(buffer.would_accept(999, 0));
    }

    #[test]
    fn test_would_accept_needed_serial() {
        let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();
        buffer.insert_with_size(1, 200, 1_000_000_000); // 1GB, creates gap at 0

        // We need serial 0, so always accept it even if over limit
        assert!(buffer.would_accept(0, 100)); // limit is 100 bytes, we have 1GB
    }

    #[test]
    fn test_would_accept_stuck() {
        let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();
        buffer.insert_with_size(2, 300, 1_000_000_000); // Gap at 0, 1

        // Buffer is stuck (can't pop), so accept everything
        assert!(!buffer.can_pop());
        assert!(buffer.would_accept(3, 100)); // Not the needed serial, but we're stuck
        assert!(buffer.would_accept(10, 100)); // Still stuck
    }

    #[test]
    fn test_would_accept_over_limit() {
        let mut buffer: ReorderBuffer<i32> = ReorderBuffer::new();
        buffer.insert_with_size(0, 100, 1000);

        // Can pop and over limit - reject non-needed serial
        assert!(buffer.can_pop());
        assert_eq!(buffer.heap_bytes(), 1000);

        // Serial 1 is not the needed serial (0 is), and we're at 1000 bytes
        // With limit 500, we're over limit, should reject
        assert!(!buffer.would_accept(1, 500));

        // But serial 0 (the next_seq) should always be accepted
        // Actually next_seq is 0, and we have item at 0, so next needed is after pop
        // Let me think: next_seq=0, item at 0 exists, so would_accept(0, ...) checks serial==next_seq
        // which is true, so it returns true. But we already have 0! That's fine, the check
        // is just for admission control before insertion.

        // With limit 2000, we're under limit, should accept
        assert!(buffer.would_accept(1, 2000));
    }
}
