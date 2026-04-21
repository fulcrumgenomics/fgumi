//! Sink trait: consumes items at the pipeline tail.
//!
//! A Sink has its own internal threads (it's the tail of the pipeline).
//! The engine calls [`Sink::run`] which blocks until the source is drained.
//! The sink is responsible for reordering items by sequence number before
//! final output.
//!
//! # Why Source and Sink are separate traits
//!
//! A `Sink` is strictly a consumer: it pops from an `input` queue and has
//! no output queue. A `Source` is strictly a producer: it pushes into an
//! `output` queue and has no input queue. Because their shapes differ at
//! the engine boundary — direction of data flow, ownership of the queue,
//! lifecycle of `close()` — unifying them into a single trait would
//! obscure those asymmetries at the type level.
//!
//! Keeping them separate lets the compiler enforce the data-flow direction
//! (you cannot pass a Sink where a Source is expected and vice versa) and
//! lets each trait document the contract most relevant to its role.

pub mod bam;
pub use bam::{BamSink, DeferredHeader};

pub mod count;
pub use count::CountSink;

pub mod bam_file_write;
pub use bam_file_write::BamFileWrite;

pub mod fastq_file;
pub use fastq_file::FastqFileSink;

use std::sync::Arc;

use super::cancel::CancelToken;
use super::queue::StageQueue;
use super::stage::SequencedItem;

/// Input queue abstraction passed to [`Sink::run`].
///
/// Two concrete implementations exist:
///
/// - Direct typed case: `Arc<StageQueue<T>>` — used by the single-stage
///   driver and by tests, with zero boxing overhead.
/// - Erased case: `UnerasingInputQueue<T>` — used by the multi-stage driver,
///   which stores inter-stage items as `Box<dyn Any + Send>`. The wrapper
///   downcasts on pop so sinks can work with `T` directly.
///
/// The trait exists so `Sink` implementations are written once against a
/// single API and the driver picks the appropriate backing queue.
pub trait InputQueue<T>: Send {
    /// Pop an item. Returns `None` if the queue is empty (check
    /// [`InputQueue::is_drained`] to distinguish "drained" from "empty").
    fn pop(&self) -> Option<SequencedItem<T>>;

    /// Whether the queue is closed and empty (i.e., fully drained).
    fn is_drained(&self) -> bool;

    /// Whether the queue has been closed.
    fn is_closed(&self) -> bool;

    /// Whether the queue is currently empty.
    fn is_empty(&self) -> bool;
}

impl<T: Send + 'static> InputQueue<T> for Arc<StageQueue<T>> {
    fn pop(&self) -> Option<SequencedItem<T>> {
        StageQueue::pop(self)
    }

    fn is_drained(&self) -> bool {
        StageQueue::is_drained(self)
    }

    fn is_closed(&self) -> bool {
        StageQueue::is_closed(self)
    }

    fn is_empty(&self) -> bool {
        StageQueue::is_empty(self)
    }
}

/// Consumes items at the pipeline tail.
pub trait Sink: Send {
    type Input: Send + 'static;

    /// Run the sink to completion.
    ///
    /// Pops items from `input` until the queue is drained. Reorders by
    /// sequence number (for deterministic output). Must check `cancel`
    /// periodically and exit promptly if set.
    ///
    /// # Errors
    ///
    /// Returns an error if the sink fails to run to completion.
    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<Self::Input>>,
        cancel: CancelToken,
    ) -> anyhow::Result<()>;

    fn name(&self) -> &'static str;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;
    use std::sync::Mutex;

    struct RecordingSink {
        out: Arc<Mutex<Vec<u64>>>,
    }

    impl Sink for RecordingSink {
        type Input = u64;

        fn run(
            self: Box<Self>,
            input: Box<dyn InputQueue<Self::Input>>,
            cancel: CancelToken,
        ) -> anyhow::Result<()> {
            let mut buffer: Vec<(u64, u64)> = Vec::new();
            loop {
                if cancel.is_cancelled() {
                    break;
                }
                if let Some(item) = input.pop() {
                    buffer.push((item.seq, item.item));
                } else if input.is_drained() {
                    break;
                } else {
                    std::thread::yield_now();
                }
            }
            buffer.sort_by_key(|(seq, _)| *seq);
            let mut out = self.out.lock().unwrap();
            for (_, v) in buffer {
                out.push(v);
            }
            Ok(())
        }

        fn name(&self) -> &'static str {
            "RecordingSink"
        }
    }

    #[test]
    fn test_sink_reorders_by_sequence() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let queue: Arc<StageQueue<u64>> = Arc::new(StageQueue::new("test", 100, 10_000, tracker));

        // Push out of order.
        queue.push(SequencedItem::new(2, 200, 8)).unwrap();
        queue.push(SequencedItem::new(0, 100, 8)).unwrap();
        queue.push(SequencedItem::new(1, 150, 8)).unwrap();
        queue.close();

        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(RecordingSink { out: out.clone() });
        sink.run(Box::new(queue), CancelToken::new()).unwrap();

        assert_eq!(*out.lock().unwrap(), vec![100, 150, 200]);
    }
}
