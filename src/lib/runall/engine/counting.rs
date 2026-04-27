//! Counting adapters for source and sink queues.
//!
//! Used by the pipeline driver to catch silent data loss: if the source
//! produces N items and the sink consumes M where N != M, the driver treats
//! that as a pipeline bug. Without this check, an off-by-one in a stage,
//! a swallowed error, or a dropped item in a scheduler would return Ok from
//! the driver even though some items were lost.
//!
//! The adapters wrap an existing `OutputQueue<T>` / `InputQueue<T>` and
//! increment an `Arc<AtomicU64>` on every successful push / pop. No
//! `Source` / `Sink` trait changes are required — the driver installs the
//! adapters between the underlying queue and the caller's trait object.

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use super::cancel::CancelToken;
use super::sink::InputQueue;
use super::source::OutputQueue;
use super::stage::SequencedItem;

/// Wraps an `OutputQueue<T>` and increments a counter on every successful push.
pub struct CountingOutputQueue<T> {
    inner: Box<dyn OutputQueue<T>>,
    produced: Arc<AtomicU64>,
}

impl<T: Send + 'static> CountingOutputQueue<T> {
    #[must_use]
    pub fn new(inner: Box<dyn OutputQueue<T>>, produced: Arc<AtomicU64>) -> Self {
        Self { inner, produced }
    }
}

impl<T: Send + 'static> OutputQueue<T> for CountingOutputQueue<T> {
    fn push(&self, item: SequencedItem<T>) -> Result<(), SequencedItem<T>> {
        match self.inner.push(item) {
            Ok(()) => {
                self.produced.fetch_add(1, Ordering::Relaxed);
                Ok(())
            }
            Err(e) => Err(e),
        }
    }

    fn push_until_cancelled(
        &self,
        item: SequencedItem<T>,
        cancel: &CancelToken,
    ) -> Result<(), SequencedItem<T>> {
        match self.inner.push_until_cancelled(item, cancel) {
            Ok(()) => {
                self.produced.fetch_add(1, Ordering::Relaxed);
                Ok(())
            }
            Err(e) => Err(e),
        }
    }

    fn close(&self) {
        self.inner.close();
    }

    fn is_closed(&self) -> bool {
        self.inner.is_closed()
    }
}

/// Wraps an `InputQueue<T>` and increments a counter on every successful pop.
pub struct CountingInputQueue<T> {
    inner: Box<dyn InputQueue<T>>,
    consumed: Arc<AtomicU64>,
}

impl<T: Send + 'static> CountingInputQueue<T> {
    #[must_use]
    pub fn new(inner: Box<dyn InputQueue<T>>, consumed: Arc<AtomicU64>) -> Self {
        Self { inner, consumed }
    }
}

impl<T: Send + 'static> InputQueue<T> for CountingInputQueue<T> {
    fn pop(&self) -> Option<SequencedItem<T>> {
        let item = self.inner.pop()?;
        self.consumed.fetch_add(1, Ordering::Relaxed);
        Some(item)
    }

    fn is_drained(&self) -> bool {
        self.inner.is_drained()
    }

    fn is_closed(&self) -> bool {
        self.inner.is_closed()
    }

    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::queue::StageQueue;

    #[test]
    fn test_counting_output_queue_increments_on_push() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let q: Arc<StageQueue<u64>> = Arc::new(StageQueue::new("test", 10, 1000, tracker));
        let produced = Arc::new(AtomicU64::new(0));
        let counting = CountingOutputQueue::new(
            Box::new(q.clone()) as Box<dyn OutputQueue<u64>>,
            produced.clone(),
        );

        counting.push(SequencedItem::new(0, 42, 8)).unwrap();
        counting.push(SequencedItem::new(1, 43, 8)).unwrap();
        assert_eq!(produced.load(Ordering::Relaxed), 2);
    }

    #[test]
    fn test_counting_output_queue_does_not_increment_on_failed_push() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        // Capacity 1 so the second push fails.
        let q: Arc<StageQueue<u64>> = Arc::new(StageQueue::new("test", 1, 1000, tracker));
        let produced = Arc::new(AtomicU64::new(0));
        let counting = CountingOutputQueue::new(
            Box::new(q.clone()) as Box<dyn OutputQueue<u64>>,
            produced.clone(),
        );

        counting.push(SequencedItem::new(0, 42, 8)).unwrap();
        let result = counting.push(SequencedItem::new(1, 43, 8));
        assert!(result.is_err());
        assert_eq!(produced.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_counting_input_queue_increments_on_pop() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let q: Arc<StageQueue<u64>> = Arc::new(StageQueue::new("test", 10, 1000, tracker));
        q.push(SequencedItem::new(0, 42, 8)).unwrap();
        q.push(SequencedItem::new(1, 43, 8)).unwrap();

        let consumed = Arc::new(AtomicU64::new(0));
        let counting = CountingInputQueue::new(
            Box::new(q.clone()) as Box<dyn InputQueue<u64>>,
            consumed.clone(),
        );

        assert!(counting.pop().is_some());
        assert!(counting.pop().is_some());
        assert!(counting.pop().is_none());
        assert_eq!(consumed.load(Ordering::Relaxed), 2);
    }
}
