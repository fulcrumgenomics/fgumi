//! Per-thread worker loop.
//!
//! A worker repeatedly tries to pop an input, process it via a Stage, and push
//! the output. If push fails (queue full), the produced item is stashed in a
//! held-item slot and retried on the next iteration before doing new work.
//! Workers exit when the input queue is drained (closed + empty). The pipeline
//! driver is responsible for closing the output queue after all workers finish.

use std::sync::Arc;

use super::backoff::Backoff;
use super::cancel::CancelToken;
use super::queue::StageQueue;
use super::stage::{SequencedItem, Stage};

/// Execute one worker thread for a single Stage (single-stage worker loop).
///
/// Used for early tests and for stages with only one worker. Multi-stage
/// scheduling uses [`crate::runall::engine::scheduler`] in a later task.
///
/// # Panics
///
/// Panics if the stage's `process` callback is called more than once per item.
///
/// # Errors
///
/// Returns an error if the stage's `process` function fails.
#[allow(clippy::needless_pass_by_value)]
pub fn run_worker<S: Stage>(
    mut stage: S,
    input: Arc<StageQueue<S::Input>>,
    output: Arc<StageQueue<S::Output>>,
    cancel: CancelToken,
) -> anyhow::Result<()> {
    let mut held: Option<SequencedItem<S::Output>> = None;
    let mut backoff = Backoff::new();

    loop {
        if cancel.is_cancelled() {
            return Ok(());
        }

        // Try to drain any held item first (non-blocking). If push fails, we
        // do NOT tight-spin retrying only the held slot: that pattern starves
        // the pipeline when every worker holds a blocked item and nobody pops
        // new input. Instead, we skip to the held.is_some() check below and,
        // if still full, back off before re-attempting. See v1 issues
        // #251/#253 for the history of this fix.
        if let Some(item) = held.take() {
            if let Err(returned) = output.push(item) {
                held = Some(returned);
            }
        }

        // If the held slot is still full we cannot take new work (would need
        // a second held slot). Back off, then re-loop so the held-push
        // retry and cancellation check happen again.
        if held.is_some() {
            backoff.snooze();
            continue;
        }

        // Try to pop an input.
        let Some(input_item) = input.pop() else {
            if input.is_drained() {
                return Ok(());
            }
            backoff.snooze();
            continue;
        };

        // Made forward progress — reset backoff so the next stall starts fresh.
        backoff.reset();

        // Process.
        let seq = input_item.seq;
        let mut captured: Option<S::Output> = None;
        stage.process(input_item.item, &mut |v| {
            assert!(captured.is_none(), "stage emitted more than one output");
            captured = Some(v);
        })?;
        if let Some(output_value) = captured {
            let mem = stage.output_memory_estimate(&output_value);
            let output_item = SequencedItem::new(seq, output_value, mem);
            // Try push; stash on failure.
            match output.push(output_item) {
                Ok(()) => {}
                Err(returned) => held = Some(returned),
            }
        }
        // Else: suppress-on-empty; just continue the loop.
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::stage::Parallelism;

    struct DoubleStage;
    impl Stage for DoubleStage {
        type Input = u32;
        type Output = u32;
        fn process(
            &mut self,
            input: u32,
            output: &mut dyn FnMut(Self::Output),
        ) -> anyhow::Result<()> {
            output(input * 2);
            Ok(())
        }
        fn parallelism(&self) -> Parallelism {
            Parallelism::Parallel
        }
        fn output_memory_estimate(&self, _output: &u32) -> usize {
            std::mem::size_of::<u32>()
        }
        fn name(&self) -> &'static str {
            "DoubleStage"
        }
    }

    #[test]
    fn test_worker_processes_all_items() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let input: Arc<StageQueue<u32>> =
            Arc::new(StageQueue::new("in", 100, 10_000, tracker.clone()));
        let output: Arc<StageQueue<u32>> = Arc::new(StageQueue::new("out", 100, 10_000, tracker));

        for i in 0u64..10u64 {
            input.push(SequencedItem::new(i, u32::try_from(i).unwrap(), 4)).unwrap();
        }
        input.close();

        run_worker(DoubleStage, input, output.clone(), CancelToken::new()).unwrap();

        let mut results: Vec<_> =
            std::iter::from_fn(|| output.pop().map(|s| (s.seq, s.item))).collect();
        results.sort_by_key(|(s, _)| *s);
        let expected: Vec<_> = (0u64..10u64).map(|i| (i, u32::try_from(i * 2).unwrap())).collect();
        assert_eq!(results, expected);
    }

    #[test]
    fn test_worker_exits_on_cancellation() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let input: Arc<StageQueue<u32>> =
            Arc::new(StageQueue::new("in", 100, 10_000, tracker.clone()));
        let output: Arc<StageQueue<u32>> = Arc::new(StageQueue::new("out", 100, 10_000, tracker));

        // Don't close input — only cancellation should exit the worker.
        let cancel = CancelToken::new();
        cancel.cancel();

        run_worker(DoubleStage, Arc::clone(&input), Arc::clone(&output), cancel).unwrap();
    }

    #[test]
    fn test_worker_held_item_on_full_output() {
        // Output queue has only 2 slots; push 5 items, verify all propagate.
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let input: Arc<StageQueue<u32>> =
            Arc::new(StageQueue::new("in", 100, 10_000, tracker.clone()));
        let output: Arc<StageQueue<u32>> = Arc::new(StageQueue::new("out", 2, 1_000, tracker));

        for i in 0u64..5u64 {
            input.push(SequencedItem::new(i, u32::try_from(i).unwrap(), 4)).unwrap();
        }
        input.close();

        // Start worker in a thread.
        let output_for_worker = output.clone();
        let worker_handle = std::thread::spawn(move || {
            run_worker(DoubleStage, input, output_for_worker, CancelToken::new())
        });

        // Drain output in this thread.
        let mut results = vec![];
        loop {
            if let Some(item) = output.pop() {
                results.push(item.item);
                if results.len() == 5 {
                    break;
                }
            } else if output.is_drained() {
                break;
            } else {
                std::thread::yield_now();
            }
        }

        worker_handle.join().unwrap().unwrap();
        results.sort_unstable();
        assert_eq!(results, vec![0, 2, 4, 6, 8]);
    }
}
