//! `ChainGraph`: per-step queue + branch-consumer tracking. Used by
//! `PipelineBuilder::build()` to assert all-outputs-wired and by the
//! runtime to construct queue topology.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StepIdx(pub usize);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BranchIdx(pub usize);

#[derive(Debug, Default)]
pub struct ChainGraph {
    /// `consumers[(producer, branch)] = Some(consumer)` if wired.
    consumers: Vec<Option<StepIdx>>,
    /// `consumer_input_slots[(producer, branch)] = Some(slot)` —
    /// which input branch of the consumer this edge feeds.
    /// Defaults to `0` for single-input consumers; multi-input
    /// consumers (`Step2` / future `StepN`) record their per-input
    /// slot explicitly so [`crate::runtime::contexts`]
    /// can build the right [`TwoInputHandles`] wrapper.
    consumer_input_slots: Vec<Option<usize>>,
    /// `branch_count[step]` — number of output branches for that step.
    branch_counts: Vec<usize>,
    /// `input_arities[step]` — number of input branches for that
    /// step (1 for single-input [`Step`] impls, 2 for [`Step2`],
    /// N for future `StepN`). Defaults to 1 via [`Self::register_step`];
    /// multi-input steps register via
    /// [`Self::register_step_with_input_arity`].
    input_arities: Vec<usize>,
    /// Step names for error messages.
    step_names: Vec<&'static str>,
}

impl ChainGraph {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a step with a single input branch (`input_arity = 1`).
    /// Convenience shorthand for [`Self::register_step_with_input_arity`].
    pub fn register_step(&mut self, name: &'static str, branch_count: usize) -> StepIdx {
        self.register_step_with_input_arity(name, branch_count, 1)
    }

    /// Register a step with an explicit input arity. Multi-input
    /// steps (`Step2` and future `StepN`) pass `input_arity > 1`;
    /// sources pass `input_arity = 0` (their input is implicit).
    pub fn register_step_with_input_arity(
        &mut self,
        name: &'static str,
        branch_count: usize,
        input_arity: usize,
    ) -> StepIdx {
        let idx = StepIdx(self.branch_counts.len());
        self.branch_counts.push(branch_count);
        self.input_arities.push(input_arity);
        self.step_names.push(name);
        self.consumers.resize(self.consumers.len() + branch_count, None);
        self.consumer_input_slots.resize(self.consumer_input_slots.len() + branch_count, None);
        idx
    }

    /// Wire a (producer, branch) → consumer link, into consumer's
    /// input branch 0. Convenience shorthand for
    /// [`Self::wire_to_slot`] used by single-input consumers.
    ///
    /// # Panics
    ///
    /// Panics if the (producer, branch) is already wired (the type system +
    /// `Chain` move semantics ensure single-consumer in practice; this is a
    /// defensive check).
    pub fn wire(&mut self, producer: StepIdx, branch: BranchIdx, consumer: StepIdx) {
        self.wire_to_slot(producer, branch, consumer, 0);
    }

    /// Wire a (producer, branch) → (consumer, consumer's input slot)
    /// link. Multi-input consumers (`Step2` / future `StepN`) record
    /// the consumer's input-slot index so
    /// [`crate::runtime::contexts`] can
    /// route each edge to the right per-branch input handle.
    ///
    /// # Panics
    ///
    /// Panics if the (producer, branch) is already wired (defensive), or if
    /// `producer`, `branch`, `consumer`, or `consumer_input_slot` are out of
    /// range.
    pub fn wire_to_slot(
        &mut self,
        producer: StepIdx,
        branch: BranchIdx,
        consumer: StepIdx,
        consumer_input_slot: usize,
    ) {
        assert!(
            consumer.0 < self.input_arities.len(),
            "consumer StepIdx({}) out of range (graph has {} steps)",
            consumer.0,
            self.input_arities.len()
        );
        let consumer_arity = self.input_arities[consumer.0];
        assert!(
            consumer_input_slot < consumer_arity,
            "consumer_input_slot {consumer_input_slot} out of range for step '{}' \
             with input_arity {consumer_arity}",
            self.step_names[consumer.0]
        );
        let slot = self.consumer_slot_index(producer, branch);
        assert!(
            self.consumers[slot].is_none(),
            "branch already wired: {:?} branch {:?} → {:?}",
            producer,
            branch,
            self.consumers[slot]
        );
        self.consumers[slot] = Some(consumer);
        self.consumer_input_slots[slot] = Some(consumer_input_slot);
    }

    /// Returns the consumer-input-slot this (producer, branch) edge
    /// feeds, if wired. Single-input consumers always return
    /// `Some(0)`; multi-input consumers return `Some(0)` or
    /// `Some(1)` depending on which input branch the edge feeds.
    #[must_use]
    pub fn consumer_input_slot(&self, producer: StepIdx, branch: BranchIdx) -> Option<usize> {
        let slot = self.consumer_slot_index(producer, branch);
        self.consumer_input_slots[slot]
    }

    /// Returns the input arity of a step (1 for single-input
    /// [`Step`] impls, 2 for [`Step2`], etc.). Sources have arity 0.
    #[must_use]
    pub fn input_arity(&self, step: StepIdx) -> usize {
        self.input_arities[step.0]
    }

    /// Returns the first unwired branch as `(producer, branch_idx, branch_name)`,
    /// or `None` if every output branch is wired.
    #[must_use]
    pub fn first_unwired(&self) -> Option<(StepIdx, BranchIdx, &'static str)> {
        for (producer_usize, &branch_count) in self.branch_counts.iter().enumerate() {
            for branch in 0..branch_count {
                let producer = StepIdx(producer_usize);
                let slot = self.consumer_slot_index(producer, BranchIdx(branch));
                if self.consumers[slot].is_none() {
                    let name = match branch {
                        0 => "0",
                        1 => "1",
                        2 => "2",
                        3 => "3",
                        _ => "?",
                    };
                    return Some((producer, BranchIdx(branch), name));
                }
            }
        }
        None
    }

    #[must_use]
    pub fn consumer(&self, producer: StepIdx, branch: BranchIdx) -> Option<StepIdx> {
        let slot = self.consumer_slot_index(producer, branch);
        self.consumers[slot]
    }

    #[must_use]
    pub fn step_name(&self, step: StepIdx) -> &'static str {
        self.step_names[step.0]
    }

    #[must_use]
    pub fn n_steps(&self) -> usize {
        self.branch_counts.len()
    }

    #[must_use]
    pub fn branch_count(&self, step: StepIdx) -> usize {
        self.branch_counts[step.0]
    }

    fn consumer_slot_index(&self, producer: StepIdx, branch: BranchIdx) -> usize {
        assert!(
            producer.0 < self.branch_counts.len(),
            "producer StepIdx({}) out of range (graph has {} steps)",
            producer.0,
            self.branch_counts.len()
        );
        let branch_count = self.branch_counts[producer.0];
        assert!(
            branch.0 < branch_count,
            "branch BranchIdx({}) out of range for producer '{}' with {} branches",
            branch.0,
            self.step_names[producer.0],
            branch_count
        );
        let mut offset = 0;
        for &count in &self.branch_counts[..producer.0] {
            offset += count;
        }
        offset + branch.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fresh_graph_is_empty() {
        let g = ChainGraph::new();
        assert_eq!(g.n_steps(), 0);
        assert!(g.first_unwired().is_none());
    }

    #[test]
    fn sink_with_no_branches_is_wired() {
        let mut g = ChainGraph::new();
        g.register_step("Sink", 0);
        assert!(g.first_unwired().is_none());
    }

    #[test]
    fn unwired_single_output_step_is_detected() {
        let mut g = ChainGraph::new();
        g.register_step("Source", 1);
        let (producer, branch, name) = g.first_unwired().unwrap();
        assert_eq!(producer, StepIdx(0));
        assert_eq!(branch, BranchIdx(0));
        assert_eq!(name, "0");
    }

    #[test]
    fn wired_chain_is_clean() {
        let mut g = ChainGraph::new();
        let src = g.register_step("Source", 1);
        let sink = g.register_step("Sink", 0);
        g.wire(src, BranchIdx(0), sink);
        assert!(g.first_unwired().is_none());
        assert_eq!(g.consumer(src, BranchIdx(0)), Some(sink));
    }

    #[test]
    fn unwired_fanout_branch_is_detected() {
        let mut g = ChainGraph::new();
        let src = g.register_step("Source", 1);
        let mid = g.register_step("FanOut", 2);
        let sink = g.register_step("Sink", 0);
        g.wire(src, BranchIdx(0), mid);
        g.wire(mid, BranchIdx(0), sink);
        let (producer, branch, _) = g.first_unwired().unwrap();
        assert_eq!(producer, StepIdx(1));
        assert_eq!(branch, BranchIdx(1));
    }

    #[test]
    #[should_panic(expected = "already wired")]
    fn double_wire_panics() {
        let mut g = ChainGraph::new();
        let src = g.register_step("Source", 1);
        let sink1 = g.register_step("Sink1", 0);
        let sink2 = g.register_step("Sink2", 0);
        g.wire(src, BranchIdx(0), sink1);
        g.wire(src, BranchIdx(0), sink2); // panics
    }

    #[test]
    #[should_panic(expected = "out of range")]
    fn wire_into_zero_arity_source_panics() {
        // Sources register with `input_arity = 0` (their input is implicit), so
        // even slot 0 must be rejected — there is no valid input branch to wire.
        let mut g = ChainGraph::new();
        let producer = g.register_step("Producer", 1);
        let source = g.register_step_with_input_arity("Source", 1, 0);
        g.wire(producer, BranchIdx(0), source); // slot 0, but arity 0 → panics
    }

    #[test]
    #[should_panic(expected = "consumer StepIdx(5) out of range")]
    fn wire_to_out_of_range_consumer_panics() {
        // A consumer index past the registered steps must produce a
        // deterministic range error rather than an opaque out-of-bounds panic.
        let mut g = ChainGraph::new();
        let producer = g.register_step("Producer", 1);
        g.wire(producer, BranchIdx(0), StepIdx(5)); // no step 5 registered
    }
}
