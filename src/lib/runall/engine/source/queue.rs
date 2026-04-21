//! `QueueSource<T>`: an adapter that turns an existing typed [`StageQueue<T>`]
//! into a [`Source<Output = T>`].
//!
//! Purpose: bridge a hand-rolled upstream stage (for example a Barrier-style
//! sort or a Special-style position batcher) into a [`run_multi_stage`] chain
//! without requiring the upstream stage to implement the [`Source`] trait.
//!
//! The upstream producer writes into the shared `Arc<StageQueue<T>>` and
//! calls `close()` when it is done. The `QueueSource` polls `pop()` on the
//! same queue and forwards items to the pipeline's erased output queue.
//!
//! [`run_multi_stage`]: super::super::driver::run_multi_stage

use std::sync::Arc;

use anyhow::Result;

use super::super::backoff::Backoff;
use super::super::cancel::CancelToken;
use super::super::queue::StageQueue;
use super::{OutputQueue, Source};

/// A [`Source`] that drains an existing typed [`StageQueue<T>`].
///
/// The upstream producer owns the queue and must call `close()` on it when
/// it has no more items to emit; otherwise this source will never observe
/// drain and will block indefinitely.
pub struct QueueSource<T: Send + 'static> {
    /// Pre-wired input queue fed by the upstream (hand-rolled) producer.
    input: Arc<StageQueue<T>>,
}

impl<T: Send + 'static> QueueSource<T> {
    /// Wrap a typed queue so it can be used as a pipeline source.
    #[must_use]
    pub fn new(input: Arc<StageQueue<T>>) -> Self {
        Self { input }
    }
}

impl<T: Send + 'static> Source for QueueSource<T> {
    type Output = T;

    fn run(
        self: Box<Self>,
        output: Box<dyn OutputQueue<Self::Output>>,
        cancel: CancelToken,
    ) -> Result<()> {
        let mut backoff = Backoff::new();
        loop {
            if cancel.is_cancelled() {
                break;
            }
            if let Some(item) = self.input.pop() {
                if output.push_until_cancelled(item, &cancel).is_err() {
                    break;
                }
                backoff.reset();
            } else if self.input.is_drained() {
                break;
            } else {
                backoff.snooze();
            }
        }
        output.close();
        Ok(())
    }

    fn name(&self) -> &'static str {
        "QueueSource"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::stage::SequencedItem;

    /// Pushes a small batch into the upstream queue, closes it, and asserts
    /// the adapter forwards every item to the downstream queue and closes it.
    #[test]
    fn test_queue_source_drains_queue_and_pushes_to_output() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let upstream: Arc<StageQueue<u64>> =
            Arc::new(StageQueue::new("upstream", 16, 1_000_000, tracker.clone()));
        let downstream: Arc<StageQueue<u64>> =
            Arc::new(StageQueue::new("downstream", 16, 1_000_000, tracker));

        for i in 0..5u64 {
            upstream.push(SequencedItem::new(i, i, 8)).unwrap();
        }
        upstream.close();

        let source = Box::new(QueueSource::new(upstream));
        let downstream_for_source = downstream.clone();
        source
            .run(Box::new(downstream_for_source) as Box<dyn OutputQueue<u64>>, CancelToken::new())
            .unwrap();

        assert!(downstream.is_closed());
        let mut collected = Vec::new();
        while let Some(item) = downstream.pop() {
            collected.push(item.item);
        }
        assert_eq!(collected, vec![0, 1, 2, 3, 4]);
    }

    /// A pre-cancelled token should make the adapter exit promptly and close
    /// the downstream queue without attempting to drain the upstream.
    #[test]
    fn test_queue_source_exits_on_cancellation() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let upstream: Arc<StageQueue<u64>> =
            Arc::new(StageQueue::new("upstream", 16, 1_000_000, tracker.clone()));
        let downstream: Arc<StageQueue<u64>> =
            Arc::new(StageQueue::new("downstream", 16, 1_000_000, tracker));

        let cancel = CancelToken::new();
        cancel.cancel();

        let source = Box::new(QueueSource::new(upstream));
        source.run(Box::new(downstream.clone()) as Box<dyn OutputQueue<u64>>, cancel).unwrap();

        assert!(downstream.is_closed());
    }
}
