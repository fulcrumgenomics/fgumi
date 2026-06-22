//! `SortMerge` — third step of the runall-sort three-step chain.

use std::collections::HashMap;
use std::io;
use std::sync::Arc;

use fgumi_raw_bam::RawRecord;
use fgumi_sort::{
    MergeDriver, MergeDriverDyn, MergeStep, QuerynameComparator, RawCoordinateKey, RawQuerynameKey,
    RawQuerynameLexKey, SortMergeSlot, SortOrder, TemplateKey,
};

use crate::sort::protocol::{MemoryChunkErased, SortPhase2Event};
use crate::types::{RecordBatch, RecordBatchBuilder};
use fgumi_pipeline_core::{
    Unpushed,
    held::HeldSlot,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// Default output batch size: 1024 records per emitted `RecordBatch`.
pub const DEFAULT_TARGET_BATCH_COUNT: usize = 1024;

/// Max output batches emitted per `try_run` invocation in `Merging`.
const MAX_DRAIN_BATCHES_PER_LOCK: usize = 8;

#[derive(Default)]
struct MemoryChunksByKind {
    coordinate: Vec<Vec<(RawCoordinateKey, RawRecord)>>,
    queryname_lex: Vec<Vec<(RawQuerynameLexKey, RawRecord)>>,
    queryname_natural: Vec<Vec<(RawQuerynameKey, RawRecord)>>,
    template_coordinate: Vec<Vec<(TemplateKey, RawRecord)>>,
}

impl MemoryChunksByKind {
    fn push(&mut self, chunk: MemoryChunkErased) {
        match chunk {
            MemoryChunkErased::Coordinate(v) => self.coordinate.push(v),
            MemoryChunkErased::QuerynameLex(v) => self.queryname_lex.push(v),
            MemoryChunkErased::QuerynameNatural(v) => self.queryname_natural.push(v),
            MemoryChunkErased::TemplateCoordinate(v) => self.template_coordinate.push(v),
        }
    }

    fn total_len(&self) -> usize {
        self.coordinate.len()
            + self.queryname_lex.len()
            + self.queryname_natural.len()
            + self.template_coordinate.len()
    }
}

fn build_driver(
    sort_order: SortOrder,
    slots: Vec<Arc<SortMergeSlot>>,
    chunks: MemoryChunksByKind,
    total_records: u64,
) -> Box<dyn MergeDriverDyn + Send> {
    match sort_order {
        SortOrder::Coordinate => Box::new(MergeDriver::<RawCoordinateKey>::from_slots(
            slots,
            chunks.coordinate,
            total_records,
        )),
        SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
            Box::new(MergeDriver::<RawQuerynameLexKey>::from_slots(
                slots,
                chunks.queryname_lex,
                total_records,
            ))
        }
        SortOrder::Queryname(QuerynameComparator::Natural) => {
            Box::new(MergeDriver::<RawQuerynameKey>::from_slots(
                slots,
                chunks.queryname_natural,
                total_records,
            ))
        }
        SortOrder::TemplateCoordinate => Box::new(MergeDriver::<TemplateKey>::from_slots(
            slots,
            chunks.template_coordinate,
            total_records,
        )),
    }
}

enum NextBatch {
    Batch(RecordBatch),
    Stalled(Option<RecordBatch>),
    Done(Option<RecordBatch>, u64),
}

enum SortMergeState {
    WaitingForSetup {
        slots: Vec<Arc<SortMergeSlot>>,
        slot_index: HashMap<u32, usize>,
        memory_chunks: MemoryChunksByKind,
        total_records: u64,
        expected_slot_count: Option<u32>,
        expected_memory_chunk_count: Option<u32>,
    },
    Merging {
        driver: Box<dyn MergeDriverDyn + Send>,
        builder: RecordBatchBuilder,
        next_ordinal: u64,
    },
    Done,
}

fn absorb_phase2_event(
    event: SortPhase2Event,
    slots: &mut Vec<Arc<SortMergeSlot>>,
    slot_index: &mut HashMap<u32, usize>,
    memory_chunks: &mut MemoryChunksByKind,
    total_records: &mut u64,
    expected_slot_count: &mut Option<u32>,
    expected_memory_chunk_count: &mut Option<u32>,
) {
    match event {
        SortPhase2Event::SpillReady { slot, path: _, records_ingested_so_far } => {
            if let std::collections::hash_map::Entry::Vacant(e) = slot_index.entry(slot.file_id) {
                e.insert(slots.len());
                slots.push(slot);
            }
            *total_records = (*total_records).max(records_ingested_so_far);
        }
        SortPhase2Event::MemoryChunk { chunk, records_ingested_so_far } => {
            let inner = Arc::try_unwrap(chunk).unwrap_or_else(|arc| match arc.as_ref() {
                MemoryChunkErased::Coordinate(v) => MemoryChunkErased::Coordinate(v.clone()),
                MemoryChunkErased::QuerynameLex(v) => MemoryChunkErased::QuerynameLex(v.clone()),
                MemoryChunkErased::QuerynameNatural(v) => {
                    MemoryChunkErased::QuerynameNatural(v.clone())
                }
                MemoryChunkErased::TemplateCoordinate(v) => {
                    MemoryChunkErased::TemplateCoordinate(v.clone())
                }
            });
            memory_chunks.push(inner);
            *total_records = (*total_records).max(records_ingested_so_far);
        }
        SortPhase2Event::AllAnnounced {
            slot_count,
            memory_chunk_count,
            total_records: ar_total,
        } => {
            if expected_slot_count.is_some() {
                log::warn!(
                    "SortMerge: duplicate AllAnnounced — prior {expected_slot_count:?}/\
                     {expected_memory_chunk_count:?}, new {slot_count}/{memory_chunk_count}; \
                     using the new values",
                );
            }
            *expected_slot_count = Some(slot_count);
            *expected_memory_chunk_count = Some(memory_chunk_count);
            *total_records = (*total_records).max(ar_total);
        }
    }
}

fn slot_set_complete(
    slots_len: usize,
    memory_chunks_total_len: usize,
    expected_slot_count: Option<u32>,
    expected_memory_chunk_count: Option<u32>,
) -> bool {
    matches!(
        (expected_slot_count, expected_memory_chunk_count),
        (Some(want_slots), Some(want_chunks))
            if u32::try_from(slots_len).unwrap_or(u32::MAX) == want_slots
                && u32::try_from(memory_chunks_total_len).unwrap_or(u32::MAX) == want_chunks
    )
}

/// `Serial + ByItemOrdinal` middle-of-chain merge.
pub struct SortMerge {
    state: SortMergeState,
    held: HeldSlot<Unpushed<RecordBatch>>,
    sort_order: SortOrder,
    target_batch_count: usize,
    output_byte_limit: u64,
    affinity: Affinity,
}

impl SortMerge {
    /// Build a `SortMerge` step with default batch size.
    #[must_use]
    pub fn new(sort_order: SortOrder, output_byte_limit: u64) -> Self {
        Self::with_target_batch_count(sort_order, output_byte_limit, DEFAULT_TARGET_BATCH_COUNT)
    }

    /// Build a `SortMerge` step with a custom output batch size.
    #[must_use]
    pub fn with_target_batch_count(
        sort_order: SortOrder,
        output_byte_limit: u64,
        target_batch_count: usize,
    ) -> Self {
        Self {
            state: SortMergeState::WaitingForSetup {
                slots: Vec::new(),
                slot_index: HashMap::new(),
                memory_chunks: MemoryChunksByKind::default(),
                total_records: 0,
                expected_slot_count: None,
                expected_memory_chunk_count: None,
            },
            held: HeldSlot::new(),
            sort_order,
            target_batch_count: target_batch_count.max(1),
            output_byte_limit,
            affinity: Affinity::None,
        }
    }

    /// Override the affinity hint.
    #[must_use]
    pub fn with_affinity(mut self, affinity: Affinity) -> Self {
        self.affinity = affinity;
        self
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        let Some(unpushed) = self.held.take() else {
            return true;
        };
        match ctx.outputs.retry(unpushed) {
            Ok(()) => true,
            Err(again) => {
                self.held.put(again);
                false
            }
        }
    }

    /// Drains every currently-available input event into the setup state and
    /// returns the number absorbed. The drain is intentionally unbounded — the
    /// upstream queue is byte-bounded, so memory is gated on the producer side.
    ///
    /// # Panics
    ///
    /// Panics if `self.state` is not `WaitingForSetup`.
    fn absorb_events_into_setup(&mut self, ctx: &mut StepCtx<'_, Self>) -> usize {
        let SortMergeState::WaitingForSetup {
            slots,
            slot_index,
            memory_chunks,
            total_records,
            expected_slot_count,
            expected_memory_chunk_count,
        } = &mut self.state
        else {
            unreachable!("absorb_events_into_setup called outside WaitingForSetup state");
        };
        let mut absorbed = 0usize;
        while let Some(event) = ctx.input.pop() {
            absorb_phase2_event(
                event,
                slots,
                slot_index,
                memory_chunks,
                total_records,
                expected_slot_count,
                expected_memory_chunk_count,
            );
            absorbed += 1;
        }
        absorbed
    }

    fn is_ready_to_merge(&self) -> bool {
        let SortMergeState::WaitingForSetup {
            slots,
            memory_chunks,
            expected_slot_count,
            expected_memory_chunk_count,
            ..
        } = &self.state
        else {
            return false;
        };
        slot_set_complete(
            slots.len(),
            memory_chunks.total_len(),
            *expected_slot_count,
            *expected_memory_chunk_count,
        )
    }

    /// # Panics
    ///
    /// Panics if `self.state` is not `Merging`.
    fn next_batch(&mut self) -> io::Result<NextBatch> {
        let target = self.target_batch_count;
        let byte_limit = self.output_byte_limit;
        let bytes_cap = usize::try_from(byte_limit).unwrap_or(usize::MAX);
        let SortMergeState::Merging { driver, builder, next_ordinal } = &mut self.state else {
            unreachable!("next_batch called outside Merging state");
        };

        let flush = |builder: &mut RecordBatchBuilder, next_ordinal: &mut u64| {
            *next_ordinal += 1;
            let next_builder = RecordBatchBuilder::with_capacity(*next_ordinal, bytes_cap, target);
            std::mem::replace(builder, next_builder).build()
        };
        let flush_partial = |builder: &mut RecordBatchBuilder, next_ordinal: &mut u64| {
            if builder.is_empty() { None } else { Some(flush(builder, next_ordinal)) }
        };

        loop {
            match driver
                .try_step()
                .map_err(|e| io::Error::other(format!("SortMerge: merge step failed: {e:#}")))?
            {
                MergeStep::Produced(bytes) => {
                    builder.push_record_bytes(bytes);
                    let count_full = builder.len() >= target;
                    let bytes_full = (builder.total_bytes() as u64) >= byte_limit;
                    if count_full || bytes_full {
                        return Ok(NextBatch::Batch(flush(builder, next_ordinal)));
                    }
                }
                MergeStep::Stalled => {
                    return Ok(NextBatch::Stalled(flush_partial(builder, next_ordinal)));
                }
                MergeStep::Done => {
                    return Ok(NextBatch::Done(
                        flush_partial(builder, next_ordinal),
                        driver.records_merged(),
                    ));
                }
            }
        }
    }

    /// # Panics
    ///
    /// Panics if `self.state` is not `Merging`.
    fn emit_batches_cooperative(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let mut delivered = 0usize;
        loop {
            match self.next_batch()? {
                NextBatch::Batch(batch) => {
                    if let Err(unpushed) = ctx.outputs.push(batch) {
                        self.held.put(unpushed);
                        return Ok(StepOutcome::Progress);
                    }
                    delivered += 1;
                    if delivered >= MAX_DRAIN_BATCHES_PER_LOCK {
                        return Ok(StepOutcome::Progress);
                    }
                }
                NextBatch::Stalled(partial) => {
                    if let Some(batch) = partial {
                        if let Err(unpushed) = ctx.outputs.push(batch) {
                            self.held.put(unpushed);
                            return Ok(StepOutcome::Progress);
                        }
                        delivered += 1;
                    }
                    return Ok(if delivered > 0 {
                        StepOutcome::Progress
                    } else {
                        StepOutcome::Contention
                    });
                }
                NextBatch::Done(partial, merged) => {
                    if let Some(batch) = partial {
                        if let Err(unpushed) = ctx.outputs.push(batch) {
                            self.held.put(unpushed);
                            return Ok(StepOutcome::Progress);
                        }
                        delivered += 1;
                    }
                    log::info!("Sort merge complete: {merged} records merged");
                    self.state = SortMergeState::Done;
                    return Ok(if delivered > 0 {
                        StepOutcome::Progress
                    } else {
                        StepOutcome::NoProgress
                    });
                }
            }
        }
    }

    fn transition_to_merging(&mut self) {
        if !matches!(&self.state, SortMergeState::WaitingForSetup { .. }) {
            return;
        }
        let SortMergeState::WaitingForSetup {
            mut slots,
            slot_index: _,
            memory_chunks,
            total_records,
            expected_slot_count: _,
            expected_memory_chunk_count: _,
        } = std::mem::replace(&mut self.state, SortMergeState::Done)
        else {
            unreachable!("just matched WaitingForSetup")
        };
        slots.sort_by_key(|s| s.file_id);
        let driver = build_driver(self.sort_order, slots, memory_chunks, total_records);
        let bytes_cap = usize::try_from(self.output_byte_limit).unwrap_or(usize::MAX);
        let builder = RecordBatchBuilder::with_capacity(0, bytes_cap, self.target_batch_count);
        self.state = SortMergeState::Merging { driver, builder, next_ordinal: 0 };
    }
}

impl Step for SortMerge {
    type Input = SortPhase2Event;
    type Outputs = OrderedBytesSingle<RecordBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SortMerge",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn affinity(&self) -> Affinity {
        self.affinity
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        if matches!(&self.state, SortMergeState::WaitingForSetup { .. }) {
            // Drain the input queue unbounded: the setup absorb is cheap (it just
            // moves `Arc`s/`Vec`s into the setup state) and the upstream queue is
            // already byte-bounded, so backpressure belongs on the producer, not
            // on a consumer-side drain cap. A cap here would cycle a full upstream
            // queue through repeated partial drains and add producer contention.
            let absorbed = self.absorb_events_into_setup(ctx);
            if !self.is_ready_to_merge() {
                if absorbed > 0 {
                    return Ok(StepOutcome::Progress);
                }
                if !ctx.input.is_drained() {
                    return Ok(StepOutcome::NoProgress);
                }
                // Input is drained but setup never completed. Fail closed for any
                // setup that saw payload or a (mismatched) `AllAnnounced`, so an
                // incomplete setup can never silently merge a partial result. The
                // only legitimate drained-but-not-ready case is a wholly empty
                // input (no slots, no chunks, no announcement), which merges to an
                // empty output.
                let SortMergeState::WaitingForSetup {
                    slots,
                    memory_chunks,
                    expected_slot_count,
                    expected_memory_chunk_count,
                    ..
                } = &self.state
                else {
                    unreachable!("state matched WaitingForSetup above");
                };
                let saw_payload = !slots.is_empty() || memory_chunks.total_len() > 0;
                let saw_expectations =
                    expected_slot_count.is_some() || expected_memory_chunk_count.is_some();
                if saw_payload || saw_expectations {
                    return Err(io::Error::other(format!(
                        "SortMerge: setup incomplete at input drain \
                         (slots={}, chunks={}, expected_slots={expected_slot_count:?}, \
                         expected_chunks={expected_memory_chunk_count:?})",
                        slots.len(),
                        memory_chunks.total_len(),
                    )));
                }
            }
            self.transition_to_merging();
        }

        match &self.state {
            SortMergeState::Merging { .. } => self.emit_batches_cooperative(ctx),
            SortMergeState::Done => Ok(StepOutcome::Finished),
            SortMergeState::WaitingForSetup { .. } => {
                unreachable!("Phase 1 must have left state non-WaitingForSetup")
            }
        }
    }
}

#[cfg(test)]
mod tests;
