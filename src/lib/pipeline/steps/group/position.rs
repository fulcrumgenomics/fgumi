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
//! Legacy's chain operates on `Vec<P>` at every `Serial` step
//! (`mi_assign`, `serialize`, ...) so each `Serial` step pays one mutex
//! acquisition per *batch*, not per *group*. The new pipeline's
//! initial implementation emitted one group per item; profile diff at
//! 53M records showed `mi_assign`-style per-group serial bottlenecks
//! collapsed CPU from 4 cores to ~1 effective core mid-run.
//!
//! Batching here matches legacy's data-flow shape: the `Group` step
//! aggregates many small completions into one batch, every downstream
//! `Serial` step amortizes its lock cost over the whole batch, and the
//! `Parallel` `Process` / `Serialize` steps each get one fat unit of work
//! per dispatch.
//!
//! ## Ordinal allocation
//!
//! Each emitted *batch* gets a fresh `batch_serial: u64` from a per-step
//! monotonic counter. Ordinals form a contiguous `0, 1, 2, …` sequence
//! so the downstream `ReorderStage` can release batches in
//! input-record order — required for the `Serial` `MiAssign` step to
//! produce deterministic MI numbering across runs (matching legacy
//! fgbio's per-strategy numbering behavior).
//!
//! ## Batch sizing
//!
//! One threshold: `target_batch_count`, defaulting to 64 groups. Mirrors
//! legacy's `template_batch_size: 500` adjusted for position-group granularity
//! (one position group ≈ 5–10 templates on typical workloads). It is a hard cap
//! on the emitted batch, not just a trigger — `emit_batch` drains at most that
//! many groups, because `add_records` can complete more in one call than the
//! target and the accumulator is only checked after extending.
//!
//! Byte-based sizing is deliberately *not* a second threshold here. An earlier
//! revision of this doc described a `target_batch_bytes` derived from the
//! producer's `output_byte_limit`, but no such field was ever implemented — the
//! transport budget is enforced by the byte-bounded output queue itself, which
//! rejects a push that would exceed it and parks the batch in `held`. Adding a
//! second in-step threshold would duplicate that accounting rather than
//! reinforce it.

use std::io;

use crate::grouper::{ProcessedPositionGroup, RawPositionGroup, RecordPositionGrouper};
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{CounterSpec, Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::DecodedRecordBatch;
use fgumi_bam_io::Grouper;
use fgumi_bam_io::MemoryEstimate;

/// Max input batches consumed per `try_run` invocation. Amortizes the
/// `Serial` mutex acquisition; mirrors `GroupBam`'s `MAX_BATCHES_PER_LOCK`.
const MAX_BATCHES_PER_LOCK: usize = 8;

/// Counter slot index: records consumed this call.
///
/// NOTE: a `molecules` counter (distinct UMI/MI groups) is deliberately NOT
/// wired here. `GroupByPosition` only forms *position* groups
/// (`RawPositionGroup`) — the UMI-adjacency split into final MI/molecule
/// groups happens downstream in `MiAssign`, which this step has no visibility
/// into. Counting emitted `RawPositionGroup`s here would misrepresent them as
/// molecules when a single position group can still split into several MI
/// groups later. See T-BW2 report for the follow-up.
const RECORDS: usize = 0;

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
/// versions of `GroupByPosition`. Downstream `Process` / `MiAssign` /
/// `Serialize` steps operate on the whole batch per dispatch, amortizing
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

    /// Enable strict template-coordinate sort-order verification (`--verify`) on
    /// the inner grouper. `header` supplies the read-group -> library-ordinal
    /// mapping used to build the template-coordinate keys. Composes with either
    /// constructor (`new` / `with_secondary_supplementary`).
    #[must_use]
    pub fn verifying(mut self, header: &noodles::sam::Header) -> Self {
        self.grouper.enable_verify(header);
        self
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

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] = &[CounterSpec::new("records", "records")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // Batch total for this call, accumulated across every input batch
        // consumed in the loop below; bumped exactly once regardless of which
        // of `try_run_inner`'s several exit paths fires.
        let mut records_this_call: u64 = 0;
        let outcome = self.try_run_inner(ctx, &mut records_this_call);
        ctx.counters.add(RECORDS, records_this_call);
        outcome
    }
}

impl GroupByPosition {
    /// Body of `Step::try_run`, factored out so the counter bump in the trait
    /// method stays a single call regardless of which early-return path below
    /// fires. `records_this_call` accumulates the batch total (records
    /// consumed) for this call; the caller bumps `ctx.counters` after this
    /// returns.
    fn try_run_inner(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        records_this_call: &mut u64,
    ) -> io::Result<StepOutcome> {
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
            *records_this_call += records.len() as u64;
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
    /// Take at most `target_batch_count` groups off the accumulator, leaving
    /// any excess for the next batch.
    ///
    /// Split out from `emit_batch` so the cap is testable without standing up a
    /// pipeline: `emit_batch` needs a `StepCtx`, and the property worth pinning
    /// — that one oversized `add_records` result becomes several bounded batches
    /// rather than one unbounded one — is entirely in this decision.
    fn drain_one_batch(&mut self) -> Vec<RawPositionGroup> {
        let take = self.accumulator.len().min(self.target_batch_count);
        self.accumulator.drain(..take).collect()
    }

    /// Package at most `target_batch_count` accumulated groups into one batch.
    ///
    /// The cap is why this drains rather than taking the whole accumulator: a
    /// single `add_records` call can complete more groups than the target — the
    /// count is only checked *after* extending — so handing the accumulator over
    /// wholesale would emit a batch of unbounded size, defeating the very
    /// threshold that exists to bound it. Anything above the cap stays in the
    /// accumulator, and the `>=` checks in `try_run` re-enter here until it
    /// falls below the target.
    ///
    /// Leftovers are also what makes a rejected push safe: the bounced batch
    /// parks in `held` while the groups that did not fit in it remain queued.
    fn emit_batch(&mut self, ctx: &mut StepCtx<'_, Self>) -> StepOutcome {
        let serial = self.next_ordinal;
        self.next_ordinal += 1;
        let groups = self.drain_one_batch();
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

    /// One `add_records` call can complete more groups than `target_batch_count`
    /// — the count is only checked after extending the accumulator — so the cap
    /// has to be applied when the batch is packaged, not when it is triggered.
    ///
    /// Without it, a burst of completions is emitted as a single unbounded
    /// batch, which is exactly the memory bound `target_batch_count` exists to
    /// provide. Pins that the excess is preserved rather than dropped, and that
    /// repeated draining converges.
    #[test]
    fn draining_caps_each_batch_and_preserves_the_excess() {
        const TARGET: usize = 64;
        let mut step = GroupByPosition::with_target_batch_count(1 << 20, TARGET);

        // 150 completed groups, as one oversized `add_records` result would leave.
        step.accumulator = (0..150)
            .map(|_| RawPositionGroup {
                group_key: fgumi_bam_io::GroupKey::default(),
                records: Vec::new(),
            })
            .collect();

        let first = step.drain_one_batch();
        assert_eq!(first.len(), TARGET, "a batch must never exceed the target");
        assert_eq!(step.accumulator.len(), 86, "the excess stays queued, not dropped");

        let second = step.drain_one_batch();
        assert_eq!(second.len(), TARGET);
        assert_eq!(step.accumulator.len(), 22);

        // The tail comes out as a short final batch and the accumulator empties,
        // so the repeated `>=` checks in `try_run` terminate.
        let third = step.drain_one_batch();
        assert_eq!(third.len(), 22, "the remainder is emitted as a short batch");
        assert!(step.accumulator.is_empty());
        assert!(step.drain_one_batch().is_empty(), "draining an empty accumulator is a no-op");
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

    // ---------------------------------------------------------------------
    // Domain counter (T-BW2): records, driven through a real pipeline
    // ---------------------------------------------------------------------

    /// Minimal single-end mapped primary record, distinguished only by
    /// `qname`. Mirrors `pipeline::steps::chain_tests::single_end_record`.
    fn single_end_record(qname: &[u8]) -> fgumi_raw_bam::RawRecord {
        let mut b = fgumi_raw_bam::SamBuilder::new();
        b.read_name(qname)
            .flags(0)
            .ref_id(0)
            .pos(100)
            .cigar_ops(&[4u32 << 4])
            .sequence(b"ACGT")
            .qualities(&[30u8; 4]);
        b.build()
    }

    /// A single-end position key at `pos`. Mirrors
    /// `pipeline::steps::chain_tests::position_key_at`.
    fn position_key_at(pos: i32) -> fgumi_bam_io::GroupKey {
        fgumi_bam_io::GroupKey { ref_id1: 0, pos1: pos, strand1: 0, ..Default::default() }
    }

    /// Wrap raw records (each paired with its pre-computed `GroupKey`) into a
    /// `DecodedRecordBatch`. Mirrors
    /// `pipeline::steps::chain_tests::decoded_batch_with_keys`.
    fn decoded_batch_with_keys(
        batch_serial: u64,
        records: Vec<(fgumi_raw_bam::RawRecord, fgumi_bam_io::GroupKey)>,
    ) -> DecodedRecordBatch {
        DecodedRecordBatch::new(
            batch_serial,
            records
                .into_iter()
                .map(|(raw, key)| fgumi_bam_io::DecodedRecord::from_raw_bytes(raw, key))
                .collect(),
        )
    }

    /// `Exclusive` source draining a `Vec<DecodedRecordBatch>`, one batch per
    /// `try_run`.
    struct ReplaySource {
        items: std::collections::VecDeque<DecodedRecordBatch>,
        held: HeldSlot<Unpushed<DecodedRecordBatch>>,
    }
    impl Step for ReplaySource {
        type Input = ();
        type Outputs = OrderedBytesSingle<DecodedRecordBatch>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ReplaySource",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 20 }],
                branch_ordering: vec![BranchOrdering::ByItemOrdinal],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(unpushed) = self.held.take() {
                match ctx.outputs.retry(unpushed) {
                    Ok(()) => {}
                    Err(again) => {
                        self.held.put(again);
                        return Ok(StepOutcome::Progress);
                    }
                }
            }
            let Some(item) = self.items.pop_front() else { return Ok(StepOutcome::Finished) };
            if let Err(unpushed) = ctx.outputs.push(item) {
                self.held.put(unpushed);
            }
            Ok(StepOutcome::Progress)
        }
    }

    /// Serial sink that sleeps briefly per received `BatchedRawPositionGroups`
    /// before recording it — mirroring
    /// `fgumi_pipeline_core::tests::CountingSink`'s per-item `thread::sleep`.
    /// `GroupByPosition` is the step under test and sits upstream of this
    /// sink, so it can burst through all its `try_run` calls well within the
    /// first sampler tick (see
    /// `fgumi_pipeline_io::source::read_bam::tests::try_run_bumps_blocks_and_bytes_read_counters`
    /// for the same pattern applied to another upstream-of-the-sink counter);
    /// throttling the terminal sink keeps the run open long enough for the
    /// sampler to observe `GroupByPosition`'s counter at its frozen final
    /// value.
    struct ThrottledGroupSink {
        received: std::sync::Arc<std::sync::Mutex<Vec<BatchedRawPositionGroups>>>,
    }
    impl Step for ThrottledGroupSink {
        type Input = BatchedRawPositionGroups;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ThrottledGroupSink",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(batch) => {
                    std::thread::sleep(std::time::Duration::from_micros(300));
                    self.received.lock().unwrap().push(batch);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// Drives `ReplaySource -> GroupByPosition -> ThrottledGroupSink` with
    /// telemetry enabled and asserts the `records` counter lands in the
    /// telemetry files with the exact expected value: `N_RECORDS` single-end
    /// mapped records across `N_POSITIONS` contiguous positions, split into
    /// many small input batches (some straddling a position boundary) so
    /// `GroupByPosition` sees many `try_run` calls. A small
    /// `target_batch_count` forces several output batches, giving the
    /// throttled sink several sleep opportunities.
    #[test]
    fn try_run_bumps_records_counter() {
        use crate::pipeline::core::builder::{InstrumentationLevel, Pipeline, PipelineConfig};
        use crate::pipeline::core::runtime::telemetry::TelemetryConfig;
        use std::time::Duration;

        const N_POSITIONS: usize = 50;
        const RECORDS_PER_POSITION: usize = 4;
        const N_RECORDS: usize = N_POSITIONS * RECORDS_PER_POSITION;
        const RECORDS_PER_INPUT_BATCH: usize = 7; // deliberately not position-aligned

        let records: Vec<(fgumi_raw_bam::RawRecord, fgumi_bam_io::GroupKey)> = (0..N_POSITIONS)
            .flat_map(|p| {
                let key = position_key_at(i32::try_from(100 + p * 10).expect("pos fits i32"));
                (0..RECORDS_PER_POSITION)
                    .map(move |r| (single_end_record(format!("p{p}_{r}").as_bytes()), key))
            })
            .collect();
        assert_eq!(records.len(), N_RECORDS);

        let batches: Vec<DecodedRecordBatch> = records
            .chunks(RECORDS_PER_INPUT_BATCH)
            .enumerate()
            .map(|(i, chunk)| decoded_batch_with_keys(i as u64, chunk.to_vec()))
            .collect();
        assert!(batches.len() >= 20, "expected many small input batches, got {}", batches.len());

        let source = ReplaySource { items: batches.into(), held: HeldSlot::new() };
        let step = GroupByPosition::with_target_batch_count(1024 * 1024, 5);
        let received = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let sink = ThrottledGroupSink { received: std::sync::Arc::clone(&received) };

        let dir = std::env::temp_dir()
            .join(format!("fgumi-tbw2-group-by-position-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let telemetry_stem = dir.join("run");

        let builder = Pipeline::builder();
        builder.chain(source).chain(step).chain(sink).into_sink_marker();
        let pipeline = builder.build().expect("pipeline builds");
        pipeline
            .run(PipelineConfig {
                threads: 1,
                instrumentation: InstrumentationLevel::Summary,
                telemetry: Some(TelemetryConfig {
                    stem: telemetry_stem.clone(),
                    interval: Duration::from_millis(1),
                }),
                ..Default::default()
            })
            .expect("pipeline runs to completion");

        // Ground truth, independent of the sampled telemetry file: every
        // consumed record must land in exactly one emitted group.
        let collected = std::mem::take(&mut *received.lock().unwrap());
        let total_out_records: usize =
            collected.iter().flat_map(|b| &b.groups).map(|g| g.records.len()).sum();
        assert_eq!(total_out_records, N_RECORDS, "every record reached the sink exactly once");

        // `GroupByPosition` is step index 1 (source=0, step=1, sink=2).
        let names = std::fs::read_to_string(dir.join("run.ticks.counter_names.tsv")).unwrap();
        let name_rows: Vec<&str> = names.lines().skip(1).filter(|l| l.starts_with("1\t")).collect();
        assert_eq!(
            name_rows,
            vec!["1\t0\trecords\trecords"],
            "GroupByPosition declares exactly one records counter"
        );

        let counters = std::fs::read_to_string(dir.join("run.ticks.counters.tsv")).unwrap();
        let last_value = counters
            .lines()
            .skip(1)
            .filter(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                f[2] == "1" && f[3] == "0"
            })
            .last()
            .map(|l| l.split('\t').nth(5).unwrap().parse::<u64>().unwrap())
            .expect("records counter recorded at least once");
        assert_eq!(last_value, N_RECORDS as u64, "records counter must equal the true total");

        std::fs::remove_dir_all(&dir).ok();
    }
}
