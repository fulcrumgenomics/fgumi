//! Sink-side reorder buffer: emits items in ascending sequence-number order.
//!
//! Intermediate pipeline queues are unordered — items travel via the shortest
//! path through a worker-stealing scheduler, so they arrive at the sink in
//! arbitrary order. Any sink that needs ordered output can drop a
//! [`ReorderBuffer`] in front of its downstream writer: push every arriving
//! item, then drain in-order items whenever the next sequence number is
//! available. Gap-then-fill sequences are handled automatically.
//!
//! The engine does NOT force this buffer on every pipeline: some pipelines
//! (e.g., metrics collection that doesn't care about order) prefer the
//! unordered fast path. Ordering is opt-in at the sink.
//!
//! Sequence numbers are expected to be dense starting from `0` (the default
//! contract produced by [`crate::runall::engine::source::Source`] implementations).

use std::collections::BTreeMap;

/// Buffers items keyed by sequence number; emits them in ascending order.
///
/// Generic over the payload type `T`. The buffer stores payloads, not
/// `SequencedItem<T>`, because the sequence number is the key.
pub struct ReorderBuffer<T> {
    /// Next sequence number we are ready to emit.
    next_expected: u64,
    /// Out-of-order items awaiting their turn.
    pending: BTreeMap<u64, T>,
}

impl<T> ReorderBuffer<T> {
    /// Create an empty buffer expecting sequence number `0` first.
    #[must_use]
    pub fn new() -> Self {
        Self { next_expected: 0, pending: BTreeMap::new() }
    }

    /// Create an empty buffer expecting sequence number `start` first.
    ///
    /// Useful for testing or for resumed pipelines where the first sequence
    /// number is not `0`.
    #[must_use]
    pub fn with_start(start: u64) -> Self {
        Self { next_expected: start, pending: BTreeMap::new() }
    }

    /// Insert an item with its sequence number.
    ///
    /// The item is stored until [`ReorderBuffer::pop_ready`] or
    /// [`ReorderBuffer::drain_ready`] retrieves it. Inserting a sequence
    /// number less than the current `next_expected` is silently ignored
    /// (duplicates from retries, out-of-contract producers).
    pub fn push(&mut self, seq: u64, item: T) {
        if seq < self.next_expected {
            return;
        }
        self.pending.insert(seq, item);
    }

    /// Pop the next in-order item if it is available.
    ///
    /// Returns `None` if the buffer is empty or if the head item is ahead of
    /// `next_expected` (still waiting on a gap to fill).
    pub fn pop_ready(&mut self) -> Option<T> {
        if let Some(item) = self.pending.remove(&self.next_expected) {
            self.next_expected += 1;
            Some(item)
        } else {
            None
        }
    }

    /// Drain every contiguous in-order item currently available.
    ///
    /// Returns items in ascending sequence order. Stops at the first gap
    /// (i.e., when `next_expected` is not present in the buffer).
    pub fn drain_ready(&mut self) -> Vec<T> {
        let mut out = Vec::new();
        while let Some(item) = self.pop_ready() {
            out.push(item);
        }
        out
    }

    /// Number of items currently buffered (in-order or out-of-order).
    #[must_use]
    pub fn len(&self) -> usize {
        self.pending.len()
    }

    /// Whether the buffer holds any items.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.pending.is_empty()
    }

    /// Next sequence number the buffer will emit if available.
    #[must_use]
    pub fn next_expected(&self) -> u64 {
        self.next_expected
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
    fn test_new_starts_at_zero() {
        let buf: ReorderBuffer<u32> = ReorderBuffer::new();
        assert_eq!(buf.next_expected(), 0);
        assert!(buf.is_empty());
    }

    #[test]
    fn test_push_pop_in_order() {
        let mut buf: ReorderBuffer<u32> = ReorderBuffer::new();
        buf.push(0, 100);
        buf.push(1, 200);
        buf.push(2, 300);
        assert_eq!(buf.pop_ready(), Some(100));
        assert_eq!(buf.pop_ready(), Some(200));
        assert_eq!(buf.pop_ready(), Some(300));
        assert_eq!(buf.pop_ready(), None);
    }

    #[test]
    fn test_push_pop_out_of_order() {
        let mut buf: ReorderBuffer<u32> = ReorderBuffer::new();
        buf.push(2, 300);
        buf.push(0, 100);
        buf.push(1, 200);
        assert_eq!(buf.pop_ready(), Some(100));
        assert_eq!(buf.pop_ready(), Some(200));
        assert_eq!(buf.pop_ready(), Some(300));
        assert!(buf.is_empty());
    }

    #[test]
    fn test_gap_then_fill() {
        let mut buf: ReorderBuffer<u32> = ReorderBuffer::new();
        buf.push(0, 10);
        buf.push(2, 30);
        // 1 is missing; only 10 is ready.
        assert_eq!(buf.pop_ready(), Some(10));
        assert_eq!(buf.pop_ready(), None);
        assert_eq!(buf.len(), 1);
        // Fill the gap.
        buf.push(1, 20);
        assert_eq!(buf.pop_ready(), Some(20));
        assert_eq!(buf.pop_ready(), Some(30));
        assert!(buf.is_empty());
    }

    #[test]
    fn test_drain_ready_stops_at_gap() {
        let mut buf: ReorderBuffer<u32> = ReorderBuffer::new();
        buf.push(0, 10);
        buf.push(1, 20);
        buf.push(3, 40); // gap at 2
        let drained = buf.drain_ready();
        assert_eq!(drained, vec![10, 20]);
        assert_eq!(buf.len(), 1);
        assert_eq!(buf.next_expected(), 2);
    }

    #[test]
    fn test_drain_all_contiguous() {
        let mut buf: ReorderBuffer<u32> = ReorderBuffer::new();
        buf.push(2, 30);
        buf.push(0, 10);
        buf.push(3, 40);
        buf.push(1, 20);
        let drained = buf.drain_ready();
        assert_eq!(drained, vec![10, 20, 30, 40]);
        assert!(buf.is_empty());
    }

    #[test]
    fn test_with_start_offset() {
        let mut buf: ReorderBuffer<u32> = ReorderBuffer::with_start(100);
        buf.push(100, 1);
        buf.push(101, 2);
        assert_eq!(buf.drain_ready(), vec![1, 2]);
    }

    #[test]
    fn test_late_push_below_expected_ignored() {
        let mut buf: ReorderBuffer<u32> = ReorderBuffer::new();
        buf.push(0, 10);
        let _ = buf.pop_ready(); // advances next_expected to 1
        buf.push(0, 999); // should be ignored (already emitted)
        assert!(buf.is_empty());
        assert_eq!(buf.next_expected(), 1);
    }
}
