//! Type-erasing wrappers around [`StageQueue`] for the multi-stage driver.
//!
//! The multi-stage pipeline driver carries items between stages as
//! `Box<dyn Any + Send>` so stages with heterogeneous `Input`/`Output` types
//! can share a single queue type. The user-facing [`Source`] and [`Sink`]
//! traits, however, are written against a concrete `T`.
//!
//! `ErasingOutputQueue<T>` and `UnerasingInputQueue<T>` bridge the two:
//!
//! - `ErasingOutputQueue<T>` implements [`OutputQueue<T>`] by boxing each
//!   pushed item into `Box<dyn Any + Send>` before forwarding it to the
//!   underlying erased queue. On a failed push, it unboxes the returned
//!   item and hands `T` back to the caller.
//! - `UnerasingInputQueue<T>` implements [`InputQueue<T>`] by downcasting
//!   each popped item from `Box<dyn Any + Send>` back to `T`.
//!
//! This replaces an older architecture that used two additional typed
//! intermediate queues plus two forwarder threads to convert between
//! typed and erased representations. The wrappers eliminate those extra
//! queues and threads — items box exactly once at the source boundary
//! and unbox exactly once at the sink boundary, and the global
//! [`MemoryTracker`] sees each item in exactly one queue at a time.
//!
//! [`Source`]: super::source::Source
//! [`Sink`]: super::sink::Sink
//! [`OutputQueue<T>`]: super::source::OutputQueue
//! [`InputQueue<T>`]: super::sink::InputQueue
//! [`MemoryTracker`]: super::memory::MemoryTracker

use std::any::Any;
use std::marker::PhantomData;
use std::sync::Arc;

use super::cancel::CancelToken;
use super::queue::StageQueue;
use super::sink::InputQueue;
use super::source::OutputQueue;
use super::stage::SequencedItem;

/// `OutputQueue<T>` adapter that boxes each item as `Box<dyn Any + Send>`
/// and forwards to an erased [`StageQueue`].
pub struct ErasingOutputQueue<T> {
    inner: Arc<StageQueue<Box<dyn Any + Send>>>,
    _phantom: PhantomData<fn(T)>,
}

impl<T> ErasingOutputQueue<T> {
    /// Wrap an erased queue so sources can push typed `T` values into it.
    #[must_use]
    pub fn new(inner: Arc<StageQueue<Box<dyn Any + Send>>>) -> Self {
        Self { inner, _phantom: PhantomData }
    }
}

impl<T: Send + 'static> OutputQueue<T> for ErasingOutputQueue<T> {
    fn push(&self, item: SequencedItem<T>) -> Result<(), SequencedItem<T>> {
        let SequencedItem { seq, item, memory_estimate } = item;
        let boxed: Box<dyn Any + Send> = Box::new(item);
        let erased = SequencedItem::new(seq, boxed, memory_estimate);
        match StageQueue::push(&self.inner, erased) {
            Ok(()) => Ok(()),
            Err(returned) => Err(unerase_item::<T>(returned)),
        }
    }

    fn push_until_cancelled(
        &self,
        item: SequencedItem<T>,
        cancel: &CancelToken,
    ) -> Result<(), SequencedItem<T>> {
        let SequencedItem { seq, item, memory_estimate } = item;
        let boxed: Box<dyn Any + Send> = Box::new(item);
        let erased = SequencedItem::new(seq, boxed, memory_estimate);
        match StageQueue::push_until_cancelled(&self.inner, erased, cancel) {
            Ok(()) => Ok(()),
            Err(returned) => Err(unerase_item::<T>(returned)),
        }
    }

    fn close(&self) {
        StageQueue::close(&self.inner);
    }

    fn is_closed(&self) -> bool {
        StageQueue::is_closed(&self.inner)
    }
}

/// `InputQueue<T>` adapter that downcasts each popped
/// `Box<dyn Any + Send>` back to `T`.
pub struct UnerasingInputQueue<T> {
    inner: Arc<StageQueue<Box<dyn Any + Send>>>,
    _phantom: PhantomData<fn() -> T>,
}

impl<T> UnerasingInputQueue<T> {
    /// Wrap an erased queue so sinks can pop typed `T` values from it.
    #[must_use]
    pub fn new(inner: Arc<StageQueue<Box<dyn Any + Send>>>) -> Self {
        Self { inner, _phantom: PhantomData }
    }
}

impl<T: Send + 'static> InputQueue<T> for UnerasingInputQueue<T> {
    fn pop(&self) -> Option<SequencedItem<T>> {
        let erased = StageQueue::pop(&self.inner)?;
        Some(unerase_item::<T>(erased))
    }

    fn is_drained(&self) -> bool {
        StageQueue::is_drained(&self.inner)
    }

    fn is_closed(&self) -> bool {
        StageQueue::is_closed(&self.inner)
    }

    fn is_empty(&self) -> bool {
        StageQueue::is_empty(&self.inner)
    }
}

/// Downcast a `SequencedItem<Box<dyn Any + Send>>` back to `SequencedItem<T>`.
///
/// The multi-stage driver guarantees by construction that items flowing
/// into the sink's input queue are of type `T` (the sink's input type).
/// A downcast failure therefore indicates a programming error — the
/// pipeline was wired with mismatched types — and we panic with a
/// diagnostic that names the expected type.
fn unerase_item<T: Send + 'static>(item: SequencedItem<Box<dyn Any + Send>>) -> SequencedItem<T> {
    let SequencedItem { seq, item, memory_estimate } = item;
    let boxed: Box<T> = item.downcast::<T>().unwrap_or_else(|_| {
        panic!(
            "pipeline: erased queue item did not downcast to expected type `{}`; \
             this indicates mis-wired pipeline stages (driver invariant violation)",
            std::any::type_name::<T>(),
        )
    });
    SequencedItem::new(seq, *boxed, memory_estimate)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;

    /// Round-trip a typed item through `ErasingOutputQueue` +
    /// `UnerasingInputQueue`. The item value, sequence number, and
    /// memory estimate must all be preserved.
    #[test]
    fn test_erasing_unerasing_roundtrip_preserves_item_seq_and_memory() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let erased: Arc<StageQueue<Box<dyn Any + Send>>> =
            Arc::new(StageQueue::new("test", 16, 10_000, tracker));

        let out: ErasingOutputQueue<String> = ErasingOutputQueue::new(erased.clone());
        let inp: UnerasingInputQueue<String> = UnerasingInputQueue::new(erased.clone());

        // Push a few items with distinct seq/value/memory triples.
        out.push(SequencedItem::new(7, "alpha".to_string(), 123)).unwrap();
        out.push(SequencedItem::new(42, "beta".to_string(), 456)).unwrap();
        out.push(SequencedItem::new(100, "gamma".to_string(), 789)).unwrap();
        out.close();
        assert!(out.is_closed());
        assert!(inp.is_closed());

        let a = inp.pop().expect("first item");
        assert_eq!(a.seq, 7);
        assert_eq!(a.item, "alpha");
        assert_eq!(a.memory_estimate, 123);

        let b = inp.pop().expect("second item");
        assert_eq!(b.seq, 42);
        assert_eq!(b.item, "beta");
        assert_eq!(b.memory_estimate, 456);

        let c = inp.pop().expect("third item");
        assert_eq!(c.seq, 100);
        assert_eq!(c.item, "gamma");
        assert_eq!(c.memory_estimate, 789);

        assert!(inp.is_empty());
        assert!(inp.is_drained());
        assert!(inp.pop().is_none());
    }

    /// When the underlying erased queue rejects a push (e.g., full), the
    /// caller must get the original typed `T` back — not a boxed value.
    #[test]
    fn test_erasing_push_returns_typed_item_on_failure() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        // Capacity 1: second push must fail and return the item.
        let erased: Arc<StageQueue<Box<dyn Any + Send>>> =
            Arc::new(StageQueue::new("test", 1, 10_000, tracker));
        let out: ErasingOutputQueue<u64> = ErasingOutputQueue::new(erased);

        out.push(SequencedItem::new(0, 111, 8)).unwrap();
        let err = out.push(SequencedItem::new(1, 222, 8)).unwrap_err();
        assert_eq!(err.seq, 1);
        assert_eq!(err.item, 222);
        assert_eq!(err.memory_estimate, 8);
    }
}
