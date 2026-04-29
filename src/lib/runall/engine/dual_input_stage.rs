//! [`DualInputStage`]: a pool-driven stage that consumes from two upstream
//! queues and emits to one downstream queue.
//!
//! This is a generalisation of [`super::stage::Stage`] for pipelines that
//! need fan-in: e.g. zipper-merging an unmapped stream with a mapped
//! stream, dedup with two ordering streams, or duplex consensus that
//! re-joins per-strand paths.
//!
//! The `Stage` trait is enough for almost every pipeline stage; reach for
//! `DualInputStage` only when the stage genuinely needs two upstream
//! inputs that are produced by different upstream chains. A single-input
//! `Stage` whose `Input` is an `Either<A, B>` — fed by a single upstream
//! merge queue — is usually a worse fit because both producers have to
//! agree on a tag-and-wrap protocol; the dual-input model lets each
//! upstream remain typed against its own output.
//!
//! # Parallelism
//!
//! `DualInputStage::parallelism` returns the same [`Parallelism`] enum as
//! `Stage`. The most common case is `Parallelism::Sequential` because
//! coordinated cross-input state (e.g. a name-based zipper that buffers
//! one stream while waiting for the other) cannot trivially shard across
//! workers. `Parallel` is allowed but only when the stage's two inputs
//! contain enough information per item that no cross-item state is
//! needed.
//!
//! # Drain semantics
//!
//! The driver calls [`DualInputStage::input_a_drained`] exactly once when
//! input A's upstream queue is closed and observed empty, and the
//! analogous [`DualInputStage::input_b_drained`] for input B. These hooks
//! let the stage flush any per-input buffer state. After both drain
//! callbacks have fired, the stage's run is complete and its output queue
//! is closed.

use crate::runall::engine::stage::Parallelism;

/// A pool-driven pipeline stage with two input queues and one output
/// queue.
///
/// See module-level documentation for when to choose this over [`Stage`].
///
/// [`Stage`]: super::stage::Stage
pub trait DualInputStage: Send {
    /// Item type consumed from input A.
    type InputA: Send + 'static;
    /// Item type consumed from input B.
    type InputB: Send + 'static;
    /// Item type emitted on the single output queue.
    type Output: Send + 'static;

    /// Process one item popped from input A. May emit zero or one output
    /// items via `output(...)`.
    ///
    /// # Errors
    /// Returns an error if processing fails. The driver propagates the
    /// error to the caller and cancels sibling workers.
    fn process_a(
        &mut self,
        input: Self::InputA,
        output: &mut dyn FnMut(Self::Output),
    ) -> anyhow::Result<()>;

    /// Process one item popped from input B. Same emission contract as
    /// [`process_a`](Self::process_a).
    ///
    /// # Errors
    /// Returns an error if processing fails.
    fn process_b(
        &mut self,
        input: Self::InputB,
        output: &mut dyn FnMut(Self::Output),
    ) -> anyhow::Result<()>;

    /// Called exactly once when input A's upstream is observed drained
    /// (closed and empty). The stage may flush any per-A buffer state and
    /// emit additional outputs.
    ///
    /// Default: no-op.
    ///
    /// # Errors
    /// Returns an error if flushing fails.
    fn input_a_drained(
        &mut self,
        _output: &mut dyn FnMut(Self::Output),
    ) -> anyhow::Result<()> {
        Ok(())
    }

    /// Same as [`input_a_drained`](Self::input_a_drained) for input B.
    ///
    /// # Errors
    /// Returns an error if flushing fails.
    fn input_b_drained(
        &mut self,
        _output: &mut dyn FnMut(Self::Output),
    ) -> anyhow::Result<()> {
        Ok(())
    }

    /// Parallelism constraint, mirroring `Stage::parallelism`. Most
    /// `DualInputStage` implementations return `Sequential` because
    /// matching/joining across two streams requires cross-item state.
    fn parallelism(&self) -> Parallelism;

    /// Estimated heap+stack bytes for a single output item. Same
    /// contract as `Stage::output_memory_estimate`. Used for memory-
    /// bounded queue accounting.
    fn output_memory_estimate(&self, output: &Self::Output) -> usize;

    /// Human-readable stage name for logs and diagnostics.
    fn name(&self) -> &'static str;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A minimal `DualInputStage` for trait-shape tests: pairs items from
    /// A and B in arrival order, emitting `(a, b)` tuples whenever both
    /// streams have produced at least one item that hasn't been paired
    /// yet. On drain, any unmatched items are silently dropped.
    struct PairStage {
        pending_a: std::collections::VecDeque<u32>,
        pending_b: std::collections::VecDeque<u32>,
    }

    impl PairStage {
        fn new() -> Self {
            Self { pending_a: Default::default(), pending_b: Default::default() }
        }
    }

    impl DualInputStage for PairStage {
        type InputA = u32;
        type InputB = u32;
        type Output = (u32, u32);

        fn process_a(
            &mut self,
            input: Self::InputA,
            output: &mut dyn FnMut(Self::Output),
        ) -> anyhow::Result<()> {
            if let Some(b) = self.pending_b.pop_front() {
                output((input, b));
            } else {
                self.pending_a.push_back(input);
            }
            Ok(())
        }

        fn process_b(
            &mut self,
            input: Self::InputB,
            output: &mut dyn FnMut(Self::Output),
        ) -> anyhow::Result<()> {
            if let Some(a) = self.pending_a.pop_front() {
                output((a, input));
            } else {
                self.pending_b.push_back(input);
            }
            Ok(())
        }

        fn parallelism(&self) -> Parallelism {
            Parallelism::Sequential
        }
        fn output_memory_estimate(&self, _: &Self::Output) -> usize {
            std::mem::size_of::<(u32, u32)>()
        }
        fn name(&self) -> &'static str {
            "PairStage"
        }
    }

    /// Run a sequence of (input, value) operations through the stage and
    /// return the emitted output items. `input` is `'a'` or `'b'` to
    /// route to the corresponding `process_*` method. Avoids holding a
    /// closure-captured `&mut emitted` across calls.
    fn drive(stage: &mut PairStage, ops: &[(char, u32)]) -> Vec<(u32, u32)> {
        let mut emitted: Vec<(u32, u32)> = Vec::new();
        for &(which, v) in ops {
            let mut sink = |out: (u32, u32)| emitted.push(out);
            match which {
                'a' => stage.process_a(v, &mut sink).unwrap(),
                'b' => stage.process_b(v, &mut sink).unwrap(),
                _ => panic!("unknown input"),
            }
        }
        emitted
    }

    #[test]
    fn pair_stage_emits_when_both_streams_have_data() {
        // A1 buffers, B1 pairs with A1, B2 buffers, A2 pairs with B2.
        let mut stage = PairStage::new();
        let emitted = drive(&mut stage, &[('a', 1), ('b', 10), ('b', 20), ('a', 2)]);
        assert_eq!(emitted, vec![(1, 10), (2, 20)]);
    }

    #[test]
    fn pair_stage_drain_callbacks_are_optional() {
        // The default `input_a_drained` / `input_b_drained` are no-ops;
        // a stage with leftover buffered items doesn't error on drain.
        let mut stage = PairStage::new();
        let _ = drive(&mut stage, &[('a', 1), ('a', 2)]);
        let mut emitted: Vec<(u32, u32)> = Vec::new();
        {
            let mut sink = |out: (u32, u32)| emitted.push(out);
            stage.input_a_drained(&mut sink).unwrap();
            stage.input_b_drained(&mut sink).unwrap();
        }
        // pending_a still holds (1, 2) but that's expected for this
        // minimal stage; no output is emitted.
        assert!(emitted.is_empty());
    }
}
