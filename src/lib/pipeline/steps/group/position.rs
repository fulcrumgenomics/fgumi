//! `GroupByPosition` step for the new typed-step pipeline framework.
//!
//! Mirrors the legacy pipeline's position-based grouping (`fgumi group`):
//! consumes `DecodedRecordBatch` items in template-coordinate order,
//! accumulates records by position via `RecordPositionGrouper`, and
//! emits **batches** of completed `RawPositionGroup`s as a single
//! `BatchedRawPositionGroups` per output ordinal.
//!
//! ## Why batched (not one-group-per-item)
//!
//! Legacy's chain operates on `Vec<P>` at every Serial step
//! (`mi_assign`, `serialize`, ...) so each Serial step pays one mutex
//! acquisition per *batch*, not per *group*. The new pipeline's
//! initial implementation emitted one group per item; profile diff at
//! 53M records showed `mi_assign`-style per-group serial bottlenecks
//! collapsed CPU from 4 cores to ~1 effective core mid-run.
//!
//! Batching here matches legacy's data-flow shape: the Group step
//! aggregates many small completions into one batch, every downstream
//! Serial step amortizes its lock cost over the whole batch, and the
//! Parallel Process / Serialize steps each get one fat unit of work
//! per dispatch.
//!
//! ## Ordinal allocation
//!
//! Each emitted *batch* gets a fresh `batch_serial: u64` from a per-step
//! monotonic counter. Ordinals form a contiguous `0, 1, 2, …` sequence
//! so the downstream `ReorderStage` can release batches in
//! input-record order — required for the Serial `MiAssign` step to
//! produce deterministic MI numbering across runs (matches legacy
//! `docs/design/deterministic-mi-numbering.md`).
//!
//! ## Batch sizing
//!
//! Two thresholds, whichever fires first:
//!   - `target_batch_count`: defaults to 64 groups. Mirrors legacy's
//!     `template_batch_size: 500` adjusted for position-group
//!     granularity (one position group ≈ 5–10 templates on typical
//!     workloads).
//!   - `target_batch_bytes`: derived from the producer's
//!     `output_byte_limit` so a full batch fits comfortably in the
//!     queue's transport budget.

use std::io;

use crate::grouper::{ProcessedPositionGroup, RawPositionGroup, RecordPositionGrouper};
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::DecodedRecordBatch;
use fgumi_bam_io::Grouper;
use fgumi_bam_io::MemoryEstimate;

/// Max input batches consumed per `try_run` invocation. Amortizes the
/// Serial mutex acquisition; mirrors `GroupBam`'s `MAX_BATCHES_PER_LOCK`.
const MAX_BATCHES_PER_LOCK: usize = 8;

/// Default target batch count. Mirrors legacy's `template_batch_size: 500`
/// adjusted for position-group granularity. Position grouping aggregates
/// many records into one group; legacy's 500 templates per batch maps
/// to roughly 50–100 position groups, but small caps reduce tail
/// latency on busy loci. 64 is a balance: fewer mutex acquisitions
/// downstream, but small enough that one busy-locus group dominating a
/// batch isn't catastrophic.
pub const DEFAULT_TARGET_BATCH_COUNT: usize = 64;

/// A batch of `RawPositionGroup`s carrying its monotonic ordinal.
/// Replaces the per-group `OrderedRawPositionGroup` emitted by earlier
/// versions of `GroupByPosition`. Downstream Process / `MiAssign` /
/// Serialize steps operate on the whole batch per dispatch, amortizing
/// per-item costs the way legacy does on `Vec<P>`.
#[derive(Debug)]
pub struct BatchedRawPositionGroups {
    pub batch_serial: u64,
    pub groups: Vec<RawPositionGroup>,
}

impl BatchedRawPositionGroups {
    #[must_use]
    pub fn new(batch_serial: u64, groups: Vec<RawPositionGroup>) -> Self {
        Self { batch_serial, groups }
    }
}

impl HeapSize for BatchedRawPositionGroups {
    fn heap_size(&self) -> usize {
        self.groups.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
            + self.groups.capacity() * std::mem::size_of::<RawPositionGroup>()
    }
}

impl Ordered for BatchedRawPositionGroups {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

/// A batch of `ProcessedPositionGroup`s carrying its monotonic ordinal.
/// Same shape as `BatchedRawPositionGroups`; the `process_step`
/// transforms one into the other one-batch-at-a-time.
#[derive(Debug)]
pub struct BatchedProcessedPositionGroups {
    pub batch_serial: u64,
    pub groups: Vec<ProcessedPositionGroup>,
}

impl BatchedProcessedPositionGroups {
    #[must_use]
    pub fn new(batch_serial: u64, groups: Vec<ProcessedPositionGroup>) -> Self {
        Self { batch_serial, groups }
    }
}

impl HeapSize for BatchedProcessedPositionGroups {
    fn heap_size(&self) -> usize {
        self.groups.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
            + self.groups.capacity() * std::mem::size_of::<ProcessedPositionGroup>()
    }
}

impl Ordered for BatchedProcessedPositionGroups {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

/// `Serial + ByItemOrdinal` position grouper. Wraps a
/// `RecordPositionGrouper` and emits batches of completed
/// `RawPositionGroup`s.
///
/// State is held behind the framework's per-step mutex (the runtime
/// stores a `Serial` step inside `Arc<Mutex<...>>` and acquires
/// per-`try_run`); the inner `RecordPositionGrouper` mutates freely
/// while the lock is held. Held-slot retry pattern matches `GroupBam`:
/// when an emit is rejected by downstream backpressure, the partial
/// state stays buffered until the held slot drains.
pub struct GroupByPosition {
    grouper: RecordPositionGrouper,
    /// Self-managed monotonic ordinal — each emitted *batch* gets the
    /// next value. Required for `BranchOrdering::ByItemOrdinal`'s
    /// `ReorderStage` to see contiguous serials.
    next_ordinal: u64,
    /// Accumulator: completed groups produced by `add_records` waiting
    /// to be packaged into a batch and pushed downstream. Drained when
    /// it reaches `target_batch_count`, on the input-drained completion
    /// path in `try_run`, or via the held-slot path when a prior push was
    /// rejected.
    accumulator: Vec<RawPositionGroup>,
    /// Held output slot when downstream rejected the most recent push.
    held: HeldSlot<Unpushed<BatchedRawPositionGroups>>,
    /// Set once the final-flush path has called the (non-idempotent)
    /// `grouper.finish()`; guards against a second call across the
    /// multi-pass completion drain.
    finalized: bool,
    target_batch_count: usize,
    output_byte_limit: u64,
    name: &'static str,
}

impl GroupByPosition {
    /// Construct a `GroupByPosition` step. `output_byte_limit` controls
    /// the byte-bounded queue capacity for emitted batches; sized via
    /// `BamPipelineTuning::per_step_byte_limit` in production.
    /// Uses `DEFAULT_TARGET_BATCH_COUNT` for the per-batch group count.
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self::with_target_batch_count(output_byte_limit, DEFAULT_TARGET_BATCH_COUNT)
    }

    /// Construct with a custom target batch count.
    #[must_use]
    pub fn with_target_batch_count(output_byte_limit: u64, target_batch_count: usize) -> Self {
        Self::with_grouper(RecordPositionGrouper::new(), output_byte_limit, target_batch_count)
    }

    /// Construct with `RecordPositionGrouper::with_secondary_supplementary()`
    /// — used by `fgumi dedup`, which must place secondary/supplementary
    /// records into the same position group as their adjacent primary so
    /// the duplicate flag propagates uniformly across split alignments.
    /// The standard `new` constructor uses the default grouper, which
    /// excludes secondary/supplementary records by position (matching
    /// `fgumi group` semantics).
    #[must_use]
    pub fn with_secondary_supplementary(output_byte_limit: u64) -> Self {
        Self::with_grouper(
            RecordPositionGrouper::with_secondary_supplementary(),
            output_byte_limit,
            DEFAULT_TARGET_BATCH_COUNT,
        )
    }

    fn with_grouper(
        grouper: RecordPositionGrouper,
        output_byte_limit: u64,
        target_batch_count: usize,
    ) -> Self {
        Self {
            grouper,
            next_ordinal: 0,
            accumulator: Vec::with_capacity(target_batch_count),
            held: HeldSlot::new(),
            finalized: false,
            target_batch_count: target_batch_count.max(1),
            output_byte_limit,
            name: "GroupByPosition",
        }
    }
}

impl Step for GroupByPosition {
    type Input = DecodedRecordBatch;
    type Outputs = OrderedBytesSingle<BatchedRawPositionGroups>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain held slot first.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // 2. If the accumulator is full enough, emit a batch. Don't
        // pull more input until the batch lands — keeps `accumulator`
        // bounded to ~`target_batch_count` groups.
        if self.accumulator.len() >= self.target_batch_count {
            return Ok(self.emit_batch(ctx));
        }

        // 3. Process up to `MAX_BATCHES_PER_LOCK` input batches per
        // call to amortize the Serial mutex acquisition. Stop early
        // once the accumulator hits the target so memory stays
        // bounded.
        let mut did_work = false;
        for _ in 0..MAX_BATCHES_PER_LOCK {
            let Some(batch) = ctx.input.pop() else { break };
            did_work = true;
            let records = batch.into_records();
            let groups = self.grouper.add_records(records)?;
            self.accumulator.extend(groups);
            if self.accumulator.len() >= self.target_batch_count {
                break;
            }
        }

        // 4. If we filled the accumulator, emit. Otherwise return
        // Progress (we did work) or NoProgress (input was empty).
        if self.accumulator.len() >= self.target_batch_count {
            return Ok(self.emit_batch(ctx));
        }
        if did_work {
            return Ok(StepOutcome::Progress);
        }

        // 5. No input this call. If upstream is drained, flush the inner
        // grouper's final partial group once (guarded by `finalized` —
        // `grouper.finish()` is not idempotent), emit any leftover groups as
        // a final batch, and report `Finished` once nothing remains. `held` is
        // empty here (step 1 returned `Contention` otherwise); `emit_batch`
        // parks a bounced final push in `held` for step 1 to retry next pass.
        if ctx.input.is_drained() {
            if !self.finalized {
                self.finalized = true;
                if let Some(final_group) = self.grouper.finish()? {
                    self.accumulator.push(final_group);
                }
            }
            if !self.accumulator.is_empty() {
                return Ok(self.emit_batch(ctx));
            }
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

impl GroupByPosition {
    /// Package the accumulator into a `BatchedRawPositionGroups` and
    /// emit. On rejection, hold; the next `try_run` will retry via the
    /// held slot. Always returns `Progress` because we either emitted
    /// or queued a batch for retry.
    fn emit_batch(&mut self, ctx: &mut StepCtx<'_, Self>) -> StepOutcome {
        let serial = self.next_ordinal;
        self.next_ordinal += 1;
        let groups =
            std::mem::replace(&mut self.accumulator, Vec::with_capacity(self.target_batch_count));
        let out = BatchedRawPositionGroups::new(serial, groups);
        if let Err(unpushed) = ctx.outputs.push(out) {
            self.held.put(unpushed);
        }
        StepOutcome::Progress
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_advertises_serial_byordinal() {
        let s = GroupByPosition::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "GroupByPosition");
        assert_eq!(p.kind, StepKind::Serial);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn batched_raw_groups_carries_ordinal() {
        let groups = vec![RawPositionGroup {
            group_key: fgumi_bam_io::GroupKey::default(),
            records: Vec::new(),
        }];
        let wrapped = BatchedRawPositionGroups::new(42, groups);
        assert_eq!(wrapped.ordinal(), 42);
        assert_eq!(wrapped.groups.len(), 1);
    }
}
