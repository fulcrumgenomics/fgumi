//! Consensus calling stage: per-thread caller from factory.
//!
//! Parallel pool stage. Delegates to `fgumi_consensus`'s
//! [`ConsensusCaller`] trait. The caller is constructed per worker
//! thread via a factory closure (avoiding locking). Optionally owns a
//! per-worker [`OverlappingBasesConsensusCaller`] used to apply
//! overlapping-bases consensus correction to paired reads *before*
//! each MI group is handed to the caller — matching the standalone
//! `fgumi simplex` / `fgumi duplex` flow (`commands/simplex.rs::execute`).
//!
//! ## Input
//!
//! [`MiGroupBatch`] — the MI groups from one position group, sorted
//! ascending by MI (emitted by
//! [`crate::runall::engine::stages::group_assign::GroupAssignStage`]).
//!
//! ## Output
//!
//! [`ConsensusOutput`] — framed BAM-record bytes (each record prefixed
//! with a 4-byte little-endian `block_size`) representing the
//! consensus R1/R2 pair for each MI group, concatenated in input MI
//! order. Per-batch `count` is the sum of `ConsensusOutput::count`
//! reported by the caller for each MI group.
//!
//! ## Ordering guarantees
//!
//! Per-batch: consensus output follows MI order as presented by the
//! input batch (ascending MI). Cross-batch: no ordering is imposed
//! here — downstream reassembly must be by input batch `ordinal`. Note
//! that a [`ConsensusOutput`] does not itself carry an `ordinal` field,
//! so this stage is expected to be followed immediately by a
//! sink/coalesce that preserves the driver's batch sequence.
//!
//! ## Memory model
//!
//! 1:1 per input batch. Consensus output is usually smaller than input
//! (one R1/R2 template per MI group regardless of input family size).
//! The per-group scratch `Vec<Vec<u8>>` is freshly allocated from
//! concat-byte `MiGroup.data`; see D6 in the type-unification design
//! spec for why a further zero-copy refactor was deferred.
//!
//! ## Determinism
//!
//! Deterministic per batch: `caller.clear()` is invoked before every
//! MI group so the caller's RNG cursor and stats do not leak between
//! groups, matching the standalone `fgumi simplex` / `fgumi duplex`
//! invariant. Byte-identical output requires byte-identical input MI
//! ordering and record ordering within each group — both guaranteed
//! by `GroupAssignStage`.

use std::sync::Arc;

use anyhow::Result;
use fgumi_consensus::caller::{ConsensusCaller, ConsensusOutput};
use fgumi_consensus::overlapping::{
    AgreementStrategy, DisagreementStrategy, OverlappingBasesConsensusCaller,
    apply_overlapping_consensus,
};

use crate::runall::engine::grouping_types::{MiGroupBatch, iter_length_prefixed};
use crate::runall::engine::stage::{Parallelism, Stage};

/// Factory that produces a fresh per-worker [`ConsensusCaller`].
pub type CallerFactory = Arc<dyn Fn() -> Box<dyn ConsensusCaller> + Send + Sync>;

/// Consensus stage wrapping a single [`ConsensusCaller`] instance plus an
/// optional overlapping-bases pre-processor.
pub struct ConsensusStage {
    caller: Box<dyn ConsensusCaller>,
    /// Overlapping-bases correction applied in place before each MI group's
    /// consensus call.
    ///
    /// - `Some` — applied before each MI group, matching the standalone
    ///   simplex/duplex default (`--consensus-call-overlapping-bases true`).
    /// - `None` — skipped entirely. Selected when the user passes
    ///   `--consensus-call-overlapping-bases false` (honored by the `runall`
    ///   path via [`ConsensusStage::with_overlapping`]) or when the stage is
    ///   built via the bare [`ConsensusStage::new`] constructor.
    overlapping: Option<OverlappingBasesConsensusCaller>,
}

impl ConsensusStage {
    /// Construct for a single worker (calls the factory once). Overlapping
    /// consensus disabled — use [`ConsensusStage::with_overlapping`] to
    /// enable it.
    #[must_use]
    pub fn new(factory: &CallerFactory) -> Self {
        Self { caller: factory(), overlapping: None }
    }

    /// Construct with overlapping-bases pre-processing enabled. Each worker
    /// owns its own [`OverlappingBasesConsensusCaller`]; state carries across
    /// MI groups within a worker the same way it does in standalone simplex.
    #[must_use]
    pub fn with_overlapping(factory: &CallerFactory, overlapping: bool) -> Self {
        let overlapping = if overlapping {
            Some(OverlappingBasesConsensusCaller::new(
                AgreementStrategy::Consensus,
                DisagreementStrategy::Consensus,
            ))
        } else {
            None
        };
        Self { caller: factory(), overlapping }
    }
}

impl Stage for ConsensusStage {
    type Input = MiGroupBatch;
    type Output = ConsensusOutput;

    #[tracing::instrument(name = "consensus", skip_all)]
    fn process(&mut self, input: Self::Input, out: &mut dyn FnMut(Self::Output)) -> Result<()> {
        // Consensus output for a batch typically sits at a few tens of KB;
        // pre-size to skip the early reallocation ramp.
        let mut combined = ConsensusOutput { data: Vec::with_capacity(16 * 1024), count: 0 };
        for mi_group in input.groups {
            // Reset per-group transient state (stats, rejects buffer, RNG
            // cursor) before each molecule — matches the
            // `caller.clear()` call `fgumi simplex` makes between every MI
            // group. Without this, the RNG cursor from group N-1 carries into
            // group N; downsampling draws (and stats accumulation) become
            // input-order-dependent on previous groups, producing
            // consensus output that differs from the standalone command.
            self.caller.clear();

            // Convert concat-byte MiGroup.data → Vec<RawRecord> for the
            // ConsensusCaller trait.
            let mut records: Vec<fgumi_raw_bam::RawRecord> = iter_length_prefixed(&mi_group.data)
                .map(|r| r.map(|bytes| fgumi_raw_bam::RawRecord::from(bytes.to_vec())))
                .collect::<Result<Vec<_>>>()?;

            // Overlapping-bases correction on paired reads (R1/R2 overlap
            // regions), applied in place exactly as standalone simplex does
            // in `commands/simplex.rs::execute` — this is the pre-processing
            // step the previous audit missed.
            if let Some(ref mut oc) = self.overlapping {
                apply_overlapping_consensus(&mut records, oc)?;
            }

            let co = self.caller.consensus_reads(records)?;
            combined.data.extend(co.data);
            combined.count += co.count;
        }
        out(combined);
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, output: &Self::Output) -> usize {
        output.data.len()
    }

    fn name(&self) -> &'static str {
        "Consensus"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_consensus::caller::ConsensusCallingStats;

    /// A dummy caller that returns empty output — we just validate the plumbing.
    struct NoopCaller;

    impl ConsensusCaller for NoopCaller {
        fn consensus_reads(
            &mut self,
            _records: Vec<fgumi_raw_bam::RawRecord>,
        ) -> fgumi_consensus::error::Result<ConsensusOutput> {
            Ok(ConsensusOutput { data: vec![], count: 0 })
        }
        fn total_reads(&self) -> usize {
            0
        }
        fn total_filtered(&self) -> usize {
            0
        }
        fn consensus_reads_constructed(&self) -> usize {
            0
        }
        fn statistics(&self) -> ConsensusCallingStats {
            ConsensusCallingStats::new()
        }
        fn log_statistics(&self) {}
        fn clear(&mut self) {}
    }

    #[test]
    fn test_consensus_empty_input() {
        let factory: CallerFactory = Arc::new(|| Box::new(NoopCaller));
        let mut stage = ConsensusStage::new(&factory);
        let mut captured = None;
        stage
            .process(MiGroupBatch { groups: vec![], ordinal: 0, position_key: (0, 0) }, &mut |v| {
                captured = Some(v);
            })
            .unwrap();
        let out = captured.expect("stage must emit");
        assert_eq!(out.count, 0);
        assert_eq!(stage.name(), "Consensus");
        assert_eq!(stage.parallelism(), Parallelism::Parallel);
    }
}
