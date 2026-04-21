//! Source trait: produces the initial items fed into the pipeline.
//!
//! A Source has its own internal threads (it's the head of the pipeline).
//! The engine calls [`Source::run`] which blocks until the source is exhausted,
//! pushing items into the pipeline's first queue along the way.
//!
//! # Why Source and Sink are separate traits
//!
//! A `Source` is strictly a producer: it pushes into an `output` queue and
//! has no input queue. A `Sink` is strictly a consumer: it pops from an
//! `input` queue and has no output queue. Because their shapes differ at
//! the engine boundary — direction of data flow, ownership of the queue,
//! lifecycle of `close()` — unifying them into a single trait would
//! obscure those asymmetries at the type level.
//!
//! Keeping them separate lets the compiler enforce the data-flow direction
//! (you cannot pass a Source where a Sink is expected and vice versa) and
//! lets each trait document the contract most relevant to its role.

pub mod bam_batched;
pub use bam_batched::BamBatchedSource;

pub mod fastq;
pub use fastq::{FastqFileRead, FastqFormat, detect_fastq_format};

pub mod queue;
pub use queue::QueueSource;

use std::sync::Arc;

use super::cancel::CancelToken;
use super::queue::StageQueue;
use super::stage::SequencedItem;

/// Output queue abstraction passed to [`Source::run`].
///
/// Two concrete implementations exist:
///
/// - Direct typed case: `Arc<StageQueue<T>>` — used by the single-stage
///   driver and by tests, with zero boxing overhead.
/// - Erased case: `ErasingOutputQueue<T>` — used by the multi-stage driver,
///   which stores inter-stage items as `Box<dyn Any + Send>`. The wrapper
///   boxes items on push so sources can work with `T` directly.
///
/// The trait exists so `Source` implementations are written once against a
/// single API and the driver picks the appropriate backing queue.
pub trait OutputQueue<T>: Send + Sync {
    /// Attempt to push an item. Returns `Err(item)` on transient failure
    /// (e.g., queue full); callers should retry.
    ///
    /// # Errors
    ///
    /// Returns the item if the queue could not accept it right now.
    fn push(&self, item: SequencedItem<T>) -> Result<(), SequencedItem<T>>;

    /// Push an item, retrying until success or cancellation.
    ///
    /// # Errors
    ///
    /// Returns the item back to the caller if cancellation was observed
    /// before the queue accepted the item.
    fn push_until_cancelled(
        &self,
        item: SequencedItem<T>,
        cancel: &CancelToken,
    ) -> Result<(), SequencedItem<T>>;

    /// Mark the queue as closed (no more items will be pushed).
    fn close(&self);

    /// Whether the queue has been closed.
    fn is_closed(&self) -> bool;
}

impl<T: Send + 'static> OutputQueue<T> for Arc<StageQueue<T>> {
    fn push(&self, item: SequencedItem<T>) -> Result<(), SequencedItem<T>> {
        StageQueue::push(self, item)
    }

    fn push_until_cancelled(
        &self,
        item: SequencedItem<T>,
        cancel: &CancelToken,
    ) -> Result<(), SequencedItem<T>> {
        StageQueue::push_until_cancelled(self, item, cancel)
    }

    fn close(&self) {
        StageQueue::close(self);
    }

    fn is_closed(&self) -> bool {
        StageQueue::is_closed(self)
    }
}

/// Produces items at the pipeline head.
pub trait Source: Send {
    type Output: Send + 'static;

    /// Run the source to completion.
    ///
    /// Pushes items (with monotonically increasing sequence numbers starting at 0)
    /// into `output`. Must respect backpressure by retrying pushes that return
    /// `Err`. Must check `cancel` periodically and exit promptly if set.
    /// Calls `output.close()` when no more items will be produced.
    ///
    /// # Errors
    ///
    /// Returns an error if the source fails to run to completion.
    fn run(
        self: Box<Self>,
        output: Box<dyn OutputQueue<Self::Output>>,
        cancel: CancelToken,
    ) -> anyhow::Result<()>;

    fn name(&self) -> &'static str;
}

/// Type-erased source for storage in the pipeline builder.
pub trait AnySource: Send {
    fn name(&self) -> &'static str;
}

impl<S: Source + ?Sized> AnySource for S {
    fn name(&self) -> &'static str {
        Source::name(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;

    struct CountingSource {
        count: u64,
    }

    impl Source for CountingSource {
        type Output = u64;

        fn run(
            self: Box<Self>,
            output: Box<dyn OutputQueue<Self::Output>>,
            cancel: CancelToken,
        ) -> anyhow::Result<()> {
            for seq in 0..self.count {
                if cancel.is_cancelled() {
                    break;
                }
                let item = SequencedItem::new(seq, seq, 8);
                if output.push_until_cancelled(item, &cancel).is_err() {
                    break;
                }
            }
            output.close();
            Ok(())
        }

        fn name(&self) -> &'static str {
            "CountingSource"
        }
    }

    #[test]
    fn test_counting_source_produces_sequence() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let queue: Arc<StageQueue<u64>> = Arc::new(StageQueue::new("test", 100, 10_000, tracker));
        let source = Box::new(CountingSource { count: 5 });

        let queue_for_source = queue.clone();
        let handle = std::thread::spawn(move || {
            source.run(Box::new(queue_for_source), CancelToken::new()).unwrap();
        });
        handle.join().unwrap();

        let mut collected = vec![];
        while let Some(item) = queue.pop() {
            collected.push((item.seq, item.item));
        }
        assert_eq!(collected, vec![(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]);
        assert!(queue.is_drained());
    }

    #[test]
    fn test_source_respects_cancellation() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let queue: Arc<StageQueue<u64>> = Arc::new(StageQueue::new("test", 100, 10_000, tracker));
        let source = Box::new(CountingSource { count: 1_000_000 });
        let cancel = CancelToken::new();
        cancel.cancel();

        let queue_for_source = queue.clone();
        source.run(Box::new(queue_for_source), cancel).unwrap();

        // Source should exit quickly; queue drained/closed.
        assert!(queue.is_closed());
    }
}
