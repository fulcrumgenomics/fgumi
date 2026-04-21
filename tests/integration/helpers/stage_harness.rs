//! Synchronous in-process harness for driving a pipeline stage with fixed
//! inputs and collecting its outputs.
//!
//! Use this to test a `Stage` or `SpecialStage` without standing up the full
//! pipeline engine. The harness calls the stage's per-item entry point (or
//! drives the stage's `run` loop, for `SpecialStage`) with a never-triggered
//! `CancelToken`, and returns the collected outputs in the order they were
//! emitted.

use std::sync::Arc;

use fgumi_lib::runall::engine::CancelToken;
use fgumi_lib::runall::engine::memory::MemoryTracker;
use fgumi_lib::runall::engine::queue::StageQueue;
use fgumi_lib::runall::engine::special_stage::SpecialStage;
use fgumi_lib::runall::engine::stage::{SequencedItem, Stage};

/// Run `stage` over `inputs` in order, collecting all emitted outputs.
///
/// Calls `Stage::process` for each input in order. The stage's per-item
/// callback appends to an output vector, which is returned at the end.
///
/// The `Stage` trait has no separate finalize/flush method — each `process`
/// call is the sole opportunity to emit outputs — so the harness simply
/// returns the accumulated outputs once every input has been processed.
///
/// # Panics
///
/// Panics if `Stage::process` returns an error.
pub fn run_stage<I, O, S>(mut stage: S, inputs: Vec<I>) -> Vec<O>
where
    I: Send + 'static,
    O: Send + 'static,
    S: Stage<Input = I, Output = O>,
{
    let mut outputs: Vec<O> = Vec::new();
    for input in inputs {
        stage
            .process(input, &mut |out| outputs.push(out))
            .expect("Stage::process returned an error in run_stage");
    }
    outputs
}

/// Run a `SpecialStage` over `inputs` synchronously, collecting all outputs.
///
/// The harness builds typed input/output queues backed by `StageQueue`,
/// pre-fills the input queue with `inputs` (assigning monotonic sequence
/// numbers and a memory estimate of zero), closes the input queue, and
/// then invokes `SpecialStage::run` on the current thread with a
/// never-triggered `CancelToken`. Outputs are drained from the output queue
/// in the order the stage produced them.
///
/// # Panics
///
/// Panics if pushing to the input queue fails or if `SpecialStage::run`
/// returns an error.
pub fn run_special_stage<I, O, S>(stage: S, inputs: Vec<I>) -> Vec<O>
where
    I: Send + 'static,
    O: Send + 'static,
    S: SpecialStage<Input = I, Output = O> + 'static,
{
    // Generous limits so the harness never blocks on backpressure. Tests
    // drive a fixed, small number of items so a large queue capacity and
    // byte budget are fine.
    let capacity = inputs.len().max(1).saturating_mul(4).max(16);
    let tracker = Arc::new(MemoryTracker::new(usize::MAX / 4));
    let input_q: Arc<StageQueue<I>> =
        Arc::new(StageQueue::new("harness_in", capacity, usize::MAX / 4, tracker.clone()));
    let output_q: Arc<StageQueue<O>> =
        Arc::new(StageQueue::new("harness_out", capacity, usize::MAX / 4, tracker));

    for (seq, item) in inputs.into_iter().enumerate() {
        let seq = u64::try_from(seq).expect("input index fits in u64");
        input_q
            .push(SequencedItem::new(seq, item, 0))
            .map_err(|_| "failed to push input into harness queue")
            .expect("harness input push");
    }
    input_q.close();

    let boxed = Box::new(stage);
    boxed
        .run(Box::new(input_q), Box::new(output_q.clone()), CancelToken::new())
        .expect("SpecialStage::run returned an error in run_special_stage");

    let mut outputs: Vec<O> = Vec::new();
    while let Some(item) = output_q.pop() {
        outputs.push(item.item);
    }
    outputs
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::sync::atomic::{AtomicU64, Ordering};

    use anyhow::Result;

    use fgumi_lib::runall::engine::fastq_types::{
        FastqTemplateV2, OwnedFastqRecord, TemplateBatch,
    };
    use fgumi_lib::runall::engine::sink::InputQueue;
    use fgumi_lib::runall::engine::source::OutputQueue;
    use fgumi_lib::runall::engine::stages::count::CountStage;

    fn template(name: &[u8]) -> FastqTemplateV2 {
        FastqTemplateV2 {
            records: vec![OwnedFastqRecord {
                name: name.to_vec(),
                sequence: b"A".to_vec(),
                quality: b"I".to_vec(),
            }],
            name: name.to_vec(),
        }
    }

    fn batch(templates: Vec<FastqTemplateV2>) -> TemplateBatch {
        TemplateBatch { templates, ordinal: 0 }
    }

    #[test]
    fn test_run_stage_drives_count_stage() {
        // `CountStage` emits `()` for every input batch and accumulates the
        // total template count into a shared atomic counter. The harness
        // should deliver every batch to `process` and surface each `()` as
        // an output item.
        let counter = Arc::new(AtomicU64::new(0));
        let stage = CountStage::new(counter.clone());

        let inputs = vec![
            batch(vec![template(b"r1"), template(b"r2"), template(b"r3")]),
            batch(vec![template(b"r4"), template(b"r5")]),
            batch(vec![]),
        ];

        let outputs = run_stage(stage, inputs);

        // One `()` emitted per input batch — three in, three out.
        assert_eq!(outputs.len(), 3);
        assert_eq!(counter.load(Ordering::Relaxed), 5);
    }

    /// Trivial `SpecialStage` that forwards every input to the output queue
    /// unchanged. Lets the smoke test exercise `run_special_stage` without
    /// depending on any concrete production stage's construction.
    struct IdentitySpecialStage;

    impl SpecialStage for IdentitySpecialStage {
        type Input = u64;
        type Output = u64;

        fn run(
            self: Box<Self>,
            input: Box<dyn InputQueue<u64>>,
            output: Box<dyn OutputQueue<u64>>,
            cancel: CancelToken,
        ) -> Result<()> {
            loop {
                if cancel.is_cancelled() {
                    break;
                }
                if let Some(item) = input.pop() {
                    let _ =
                        output.push(SequencedItem::new(item.seq, item.item, item.memory_estimate));
                } else if input.is_drained() {
                    break;
                } else {
                    std::thread::yield_now();
                }
            }
            output.close();
            Ok(())
        }

        fn name(&self) -> &'static str {
            "IdentitySpecialStage"
        }
    }

    #[test]
    fn test_run_special_stage_passthrough() {
        let outputs = run_special_stage(IdentitySpecialStage, vec![10u64, 20, 30]);
        assert_eq!(outputs, vec![10, 20, 30]);
    }
}
