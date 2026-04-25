//! `crossbeam_queue::ArrayQueue`-compatible queue that pops items in the
//! order of their pipeline serial, not their arrival order.
//!
//! The BAM pipeline's `Q5` (process step output) historically used
//! `ArrayQueue<(u64, Vec<P>)>` so that serialize workers popped in completion
//! (arrival) order. That made any global counter advanced inside `serialize_fn`
//! — notably the `MoleculeId` counter in `fgumi group` and `fgumi dedup` —
//! non-deterministic across runs because the parallel process step finishes
//! batches in an arbitrary order.
//!
//! [`SerialOrderedArrayQueue`] keeps the existing ArrayQueue-like surface
//! (`push`, `pop`, `is_full`, `is_empty`, `len`) so it can be dropped in for
//! `Q5` with no other call-site changes, but internally it is a
//! `BTreeMap` keyed by serial:
//!
//! - `push((serial, item))` inserts the entry into the map (until full).
//! - `pop()` returns the entry with the smallest serial **only when that
//!   serial is the next expected one**, advancing an internal cursor on
//!   success. If the smallest serial currently buffered is greater than the
//!   cursor (because earlier serials are still in the process step), `pop`
//!   returns `None` and the worker is free to do other work.
//!
//! Because the pipeline's group step assigns serials monotonically with no
//! gaps, the cursor always eventually catches up — the queue can stall
//! serialize workers but it cannot deadlock the pipeline. The deadlock-
//! avoidance logic in [`SerialOrderedArrayQueue::push`] also lets producers
//! exceed the configured capacity when the cursor's expected serial has not
//! yet been buffered, so the missing serial can always reach the queue.

use parking_lot::Mutex;
use std::collections::BTreeMap;

/// Capacity-bounded, serial-ordered FIFO with an `ArrayQueue`-compatible API.
#[derive(Debug)]
pub struct SerialOrderedArrayQueue<T> {
    inner: Mutex<Inner<T>>,
    capacity: usize,
}

#[derive(Debug)]
struct Inner<T> {
    buf: BTreeMap<u64, T>,
    /// Next serial that is allowed to be popped. The pipeline only ever
    /// pushes monotonically increasing serials with no gaps starting from 0,
    /// so the cursor is always equal to "smallest serial that has not yet
    /// been popped or skipped".
    cursor: u64,
}

impl<T> SerialOrderedArrayQueue<T> {
    /// Create an empty queue with the given capacity.
    #[must_use]
    pub fn new(capacity: usize) -> Self {
        Self {
            inner: Mutex::new(Inner { buf: BTreeMap::new(), cursor: 0 }),
            capacity: capacity.max(1),
        }
    }

    /// Insert an `(serial, item)` pair. Returns the input back as `Err` if
    /// the queue is at capacity AND the cursor's expected serial is already
    /// buffered (so the consumer can drain).
    ///
    /// **Deadlock avoidance:** if the cursor's expected serial is **not** yet
    /// buffered, the consumer cannot make progress and the producers must be
    /// allowed to push the missing serial through. In that case the push
    /// always succeeds, even if the buffer is over capacity. This mirrors
    /// the same safety net used by [`crate::unified_pipeline::queue::OrderedQueue`]
    /// — without it, a transient burst of out-of-order completions can fill
    /// the queue and trap the worker holding the next-expected serial in its
    /// `held_processed` slot, deadlocking the pipeline.
    ///
    /// # Errors
    ///
    /// Returns the input pair as `Err` when the queue is at capacity and
    /// the consumer can drain (so backpressure should kick in at the
    /// caller).
    pub fn push(&self, item: (u64, T)) -> Result<(), (u64, T)> {
        let mut g = self.inner.lock();
        // Defensive: cursor is monotonic, items with serial < cursor have
        // already been popped and must not appear again. Convert this into a
        // panic in debug builds — silently dropping the item would mask a
        // pipeline ordering bug.
        debug_assert!(
            item.0 >= g.cursor,
            "SerialOrderedArrayQueue: stale serial {} pushed (cursor={})",
            item.0,
            g.cursor,
        );
        let consumer_can_drain = g
            .buf
            .first_key_value()
            .is_some_and(|(&k, _)| k == g.cursor);
        if g.buf.len() >= self.capacity && consumer_can_drain {
            return Err(item);
        }
        // Always accept when we don't have the cursor's expected serial yet,
        // even if over capacity, so the missing serial can arrive and unblock
        // the consumer.
        g.buf.insert(item.0, item.1);
        Ok(())
    }

    /// Pop the entry with the smallest serial **iff** that serial is the
    /// next expected one. Returns `None` if the queue is empty or if the
    /// smallest buffered serial is ahead of the cursor (i.e. some earlier
    /// batch is still in flight).
    pub fn pop(&self) -> Option<(u64, T)> {
        let mut g = self.inner.lock();
        let cursor = g.cursor;
        // Peek: only pop if the smallest key matches the cursor.
        let take = matches!(g.buf.first_key_value(), Some((&k, _)) if k == cursor);
        if !take {
            return None;
        }
        let (serial, item) = g.buf.pop_first()?;
        debug_assert_eq!(serial, cursor);
        g.cursor += 1;
        Some((serial, item))
    }

    /// True when the queue is at capacity AND the consumer can drain (i.e.
    /// the cursor's expected serial is buffered). When the cursor's serial
    /// is **not** yet buffered, the queue advertises that it has space —
    /// otherwise producer backpressure would block all workers from pushing
    /// the missing serial through, deadlocking the pipeline. See [`push`]
    /// for the matching deadlock-avoidance logic.
    pub fn is_full(&self) -> bool {
        let g = self.inner.lock();
        let consumer_can_drain = g
            .buf
            .first_key_value()
            .is_some_and(|(&k, _)| k == g.cursor);
        g.buf.len() >= self.capacity && consumer_can_drain
    }

    /// True when the queue is empty.
    pub fn is_empty(&self) -> bool {
        let g = self.inner.lock();
        g.buf.is_empty()
    }

    /// Current number of buffered items.
    pub fn len(&self) -> usize {
        let g = self.inner.lock();
        g.buf.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use std::thread;

    #[test]
    fn pops_in_serial_order_regardless_of_push_order() {
        let q: SerialOrderedArrayQueue<&'static str> = SerialOrderedArrayQueue::new(8);
        q.push((3, "three")).unwrap();
        q.push((1, "one")).unwrap();
        q.push((0, "zero")).unwrap();
        q.push((2, "two")).unwrap();

        assert_eq!(q.pop(), Some((0, "zero")));
        assert_eq!(q.pop(), Some((1, "one")));
        assert_eq!(q.pop(), Some((2, "two")));
        assert_eq!(q.pop(), Some((3, "three")));
        assert_eq!(q.pop(), None);
    }

    #[test]
    fn pop_returns_none_when_smallest_is_ahead_of_cursor() {
        let q: SerialOrderedArrayQueue<u64> = SerialOrderedArrayQueue::new(8);
        // Push out-of-order serials, but never push 0 yet.
        q.push((1, 100)).unwrap();
        q.push((2, 200)).unwrap();
        // Cursor is 0 but the smallest buffered serial is 1, so pop must
        // wait — returning None means "do other work".
        assert!(q.pop().is_none());
        assert_eq!(q.len(), 2);
        // Once the missing serial 0 arrives, pops drain in order.
        q.push((0, 10)).unwrap();
        assert_eq!(q.pop(), Some((0, 10)));
        assert_eq!(q.pop(), Some((1, 100)));
        assert_eq!(q.pop(), Some((2, 200)));
    }

    #[test]
    fn push_returns_err_when_full_and_consumer_can_drain() {
        let q: SerialOrderedArrayQueue<u8> = SerialOrderedArrayQueue::new(2);
        q.push((0, 0)).unwrap();
        q.push((1, 1)).unwrap();
        assert!(q.is_full());
        let pushed_back = q.push((2, 2));
        assert_eq!(pushed_back, Err((2, 2)));
    }

    #[test]
    fn push_succeeds_over_capacity_when_consumer_is_starved() {
        let q: SerialOrderedArrayQueue<u8> = SerialOrderedArrayQueue::new(2);
        // Cursor is 0 but only out-of-order serials are buffered.
        q.push((1, 1)).unwrap();
        q.push((2, 2)).unwrap();
        assert!(!q.is_full(), "must not advertise full while consumer is starved");
        // Going over capacity is allowed when the cursor's serial is missing
        // — otherwise the producer holding serial 0 deadlocks waiting for
        // queue space that the consumer cannot create.
        q.push((3, 3)).unwrap();
        assert_eq!(q.len(), 3);
        // Once serial 0 arrives the consumer can drain.
        q.push((0, 0)).unwrap();
        assert_eq!(q.pop(), Some((0, 0)));
        assert_eq!(q.pop(), Some((1, 1)));
        assert_eq!(q.pop(), Some((2, 2)));
        assert_eq!(q.pop(), Some((3, 3)));
    }

    #[test]
    fn many_concurrent_producers_pop_in_serial_order() {
        const N: u64 = 64;
        let q: Arc<SerialOrderedArrayQueue<u64>> = Arc::new(SerialOrderedArrayQueue::new(256));
        let mut handles = Vec::new();
        for s in 0..N {
            let q = Arc::clone(&q);
            handles.push(thread::spawn(move || {
                // Spread the pushes apart so they finish in arbitrary order.
                thread::yield_now();
                q.push((s, s * 10)).expect("push");
            }));
        }
        for h in handles {
            h.join().expect("producer");
        }
        for s in 0..N {
            // Pop must always return the next expected serial.
            assert_eq!(q.pop(), Some((s, s * 10)), "expected to pop serial {s} next");
        }
        assert!(q.pop().is_none());
    }
}
