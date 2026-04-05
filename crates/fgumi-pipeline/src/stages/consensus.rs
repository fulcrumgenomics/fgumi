//! Consensus calling pipeline stage.
//!
//! This module wraps `ConsensusCaller` implementations behind a `PipelineStage` interface.
//! Because `consensus_reads` takes `&mut self`, a single shared caller would serialize all
//! consensus work across worker threads. Instead, a **caller factory** is used: each worker
//! thread lazily creates its own `ConsensusCaller` instance the first time it processes an
//! item, stored in a `thread_local!` slot for the lifetime of the pipeline.
//!
//! This approach is sound because the scheduler uses dedicated `std::thread` workers (not
//! rayon tasks), so each thread's local storage is stable and one-to-one with workers.

use std::cell::RefCell;
use std::sync::Arc;

use anyhow::Result;
use fgumi_consensus::caller::{ConsensusCaller, ConsensusOutput};

use crate::stage::{PipelineStage, with_thread_local};
use crate::stages::group_assign::MiGroup;

// ============================================================================
// CallerFactory
// ============================================================================

/// A factory that produces a new, independent `ConsensusCaller` for each worker thread.
///
/// The factory is called at most once per worker thread. The resulting caller is stored
/// in thread-local storage and reused for all subsequent work items on that thread.
pub type CallerFactory = Arc<dyn Fn() -> Box<dyn ConsensusCaller> + Send + Sync>;

// ============================================================================
// Thread-local caller storage
// ============================================================================

thread_local! {
    /// Per-worker `ConsensusCaller` instance. Initialized on first use by `ConsensusStage`.
    static THREAD_CALLER: RefCell<Option<Box<dyn ConsensusCaller>>> = const { RefCell::new(None) };
}

// ============================================================================
// ConsensusStage
// ============================================================================

/// Pipeline stage that runs consensus calling on each [`MiGroup`].
///
/// Worker threads share a single `ConsensusStage` instance, but each thread maintains its
/// own `ConsensusCaller` in thread-local storage, created on demand by `caller_factory`.
/// This gives full parallelism without lock contention.
pub struct ConsensusStage {
    /// Factory invoked once per worker thread to create a per-thread caller.
    caller_factory: CallerFactory,
}

impl ConsensusStage {
    /// Create a new `ConsensusStage` with the given caller factory.
    ///
    /// # Arguments
    ///
    /// * `caller_factory` - A factory function that returns a fresh `ConsensusCaller`.
    ///   It will be called at most once per worker thread.
    pub fn new(caller_factory: CallerFactory) -> Self {
        Self { caller_factory }
    }
}

impl PipelineStage for ConsensusStage {
    type Input = Vec<MiGroup>;
    type Output = ConsensusOutput;

    /// Call consensus on all MI groups in the batch and merge their outputs.
    ///
    /// Lazily initializes a thread-local `ConsensusCaller` on first invocation within
    /// each worker thread, then processes each [`MiGroup`] in the batch, merging all
    /// [`ConsensusOutput`]s into a single buffer via [`ConsensusOutput::merge`].
    ///
    /// # Errors
    ///
    /// Returns any error propagated from the underlying `ConsensusCaller`.
    fn process(&self, input: Self::Input) -> Result<Self::Output> {
        let factory = &self.caller_factory;
        with_thread_local(
            &THREAD_CALLER,
            || factory(),
            |caller| {
                let mut merged = ConsensusOutput::new();
                for group in input {
                    let output = caller.consensus_reads(group.records)?;
                    merged.merge(output);
                }
                Ok(merged)
            },
        )
    }

    /// Estimate memory usage as the number of bytes in the serialized consensus output.
    fn output_memory_estimate(&self, output: &Self::Output) -> usize {
        output.data.len()
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_consensus::caller::ConsensusCallingStats;

    /// A no-op `ConsensusCaller` that always returns an empty `ConsensusOutput`.
    struct MockCaller {
        stats: ConsensusCallingStats,
    }

    impl MockCaller {
        fn new() -> Self {
            Self { stats: ConsensusCallingStats::new() }
        }
    }

    impl ConsensusCaller for MockCaller {
        fn consensus_reads(&mut self, records: Vec<Vec<u8>>) -> Result<ConsensusOutput> {
            self.stats.record_input(records.len());
            Ok(ConsensusOutput::new())
        }

        fn total_reads(&self) -> usize {
            self.stats.total_reads
        }

        fn total_filtered(&self) -> usize {
            self.stats.filtered_reads
        }

        fn consensus_reads_constructed(&self) -> usize {
            self.stats.consensus_reads
        }

        fn statistics(&self) -> ConsensusCallingStats {
            self.stats.clone()
        }

        fn log_statistics(&self) {}
    }

    #[test]
    fn test_consensus_stage_process_returns_empty_output_from_mock() {
        let factory: CallerFactory = Arc::new(|| Box::new(MockCaller::new()));
        let stage = ConsensusStage::new(factory);

        let group = MiGroup { records: vec![vec![0u8; 50], vec![0u8; 50]], mi: 1, byte_size: 100 };

        let output = stage.process(vec![group]).expect("process should succeed");
        assert_eq!(output.count, 0, "mock returns empty output");
        assert!(output.data.is_empty(), "mock returns no data bytes");
    }

    #[test]
    fn test_consensus_stage_output_memory_estimate() {
        let factory: CallerFactory = Arc::new(|| Box::new(MockCaller::new()));
        let stage = ConsensusStage::new(factory);

        let output = ConsensusOutput { data: vec![0u8; 128], count: 1 };
        assert_eq!(stage.output_memory_estimate(&output), 128);
    }

    #[test]
    fn test_consensus_stage_empty_batch() {
        let factory: CallerFactory = Arc::new(|| Box::new(MockCaller::new()));
        let stage = ConsensusStage::new(factory);

        let output = stage.process(vec![]).expect("empty batch should succeed");
        assert_eq!(output.count, 0);
        assert!(output.data.is_empty());
    }

    #[test]
    fn test_consensus_stage_multiple_mi_groups_merged() {
        let factory: CallerFactory = Arc::new(|| Box::new(MockCaller::new()));
        let stage = ConsensusStage::new(factory);

        let groups = vec![
            MiGroup { records: vec![vec![0u8; 50]], mi: 0, byte_size: 50 },
            MiGroup { records: vec![vec![0u8; 30], vec![0u8; 30]], mi: 1, byte_size: 60 },
            MiGroup { records: vec![], mi: 2, byte_size: 0 },
        ];

        let output = stage.process(groups).expect("batch should succeed");
        // MockCaller returns empty ConsensusOutput for each group, so merged is also empty.
        assert_eq!(output.count, 0);
        assert!(output.data.is_empty());
    }
}
