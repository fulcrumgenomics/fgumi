//! Reorder buffer for reassembling out-of-order items into sequence order.

use std::collections::VecDeque;

/// A buffer that reassembles out-of-order items into sequence order.
///
/// Items are inserted with sequence numbers. [`try_next`](ReorderBuffer::try_next) returns items
/// in sequence order as they become available.
///
/// Internally uses a `VecDeque<Option<T>>` indexed by offset from `next_expected`,
/// giving O(1) insert and O(1) pop with better cache locality than a `HashMap`.
pub struct ReorderBuffer<T> {
    pending: VecDeque<Option<T>>,
    next_expected: u64,
}

impl<T> ReorderBuffer<T> {
    /// Create a new empty reorder buffer starting at sequence number 0.
    #[must_use]
    pub fn new() -> Self {
        Self { pending: VecDeque::new(), next_expected: 0 }
    }

    /// Insert an item at the given sequence number. Items may be inserted in any order.
    #[expect(
        clippy::cast_possible_truncation,
        reason = "reorder buffer offset from next_expected always fits in usize in practice"
    )]
    pub fn insert(&mut self, seq: u64, item: T) {
        let idx = (seq - self.next_expected) as usize;
        // Extend with None slots if the sequence number is beyond the current length.
        if idx >= self.pending.len() {
            self.pending.resize_with(idx + 1, || None);
        }
        self.pending[idx] = Some(item);
    }

    /// Get the next item in sequence order, if available.
    ///
    /// Returns `Some(item)` if the item with sequence number `next_expected` is present,
    /// advancing `next_expected` by one. Returns `None` otherwise.
    pub fn try_next(&mut self) -> Option<T> {
        if self.pending.front().is_some_and(Option::is_some) {
            self.next_expected += 1;
            self.pending.pop_front().flatten()
        } else {
            None
        }
    }

    /// Drain all consecutive available items starting from `next_expected`.
    ///
    /// Returns a `Vec` of items in sequence order. Stops as soon as there is a gap.
    pub fn drain_ready(&mut self) -> Vec<T> {
        let mut result = Vec::new();
        while let Some(item) = self.try_next() {
            result.push(item);
        }
        result
    }

    /// Number of items buffered but not yet emitted.
    #[must_use]
    pub fn pending_count(&self) -> usize {
        self.pending.iter().filter(|slot| slot.is_some()).count()
    }

    /// The next expected sequence number.
    #[must_use]
    pub fn next_expected(&self) -> u64 {
        self.next_expected
    }

    /// Returns `true` if no items are buffered.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.pending.iter().all(Option::is_none)
    }
}

impl<T> Default for ReorderBuffer<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_in_order_passthrough() {
        let mut buf: ReorderBuffer<u32> = ReorderBuffer::new();
        buf.insert(0, 10);
        buf.insert(1, 20);
        buf.insert(2, 30);

        assert_eq!(buf.try_next(), Some(10));
        assert_eq!(buf.try_next(), Some(20));
        assert_eq!(buf.try_next(), Some(30));
        assert_eq!(buf.try_next(), None);
        assert!(buf.is_empty());
    }

    #[test]
    fn test_out_of_order() {
        let mut buf: ReorderBuffer<u32> = ReorderBuffer::new();

        // Insert out of order: 2, then 0, then 1
        buf.insert(2, 30);
        assert_eq!(buf.try_next(), None); // seq 0 not yet available

        buf.insert(0, 10);
        assert_eq!(buf.try_next(), Some(10)); // seq 0 now available

        // seq 1 not yet available
        assert_eq!(buf.try_next(), None);

        buf.insert(1, 20);
        let ready = buf.drain_ready();
        assert_eq!(ready, vec![20, 30]);
        assert!(buf.is_empty());
    }

    #[test]
    fn test_drain_ready() {
        let mut buf: ReorderBuffer<u32> = ReorderBuffer::new();
        buf.insert(0, 10);
        buf.insert(1, 20);
        buf.insert(2, 30);
        buf.insert(5, 60); // gap at 3 and 4

        let ready = buf.drain_ready();
        assert_eq!(ready, vec![10, 20, 30]);
        assert_eq!(buf.pending_count(), 1);
        assert_eq!(buf.next_expected(), 3);
    }

    #[test]
    fn test_empty() {
        let mut buf: ReorderBuffer<u32> = ReorderBuffer::new();
        assert_eq!(buf.try_next(), None);
        assert!(buf.is_empty());
        assert_eq!(buf.pending_count(), 0);
    }
}
