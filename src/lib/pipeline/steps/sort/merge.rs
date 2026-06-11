//! `SortMerge` — third step of the runall-sort three-step chain.
//!
//! ```text
//! [RecordBatch] → SortAndSpill ──SortPhase1Event──> SortSpillDecompress ──SortPhase2Event──> SortMerge → [RecordBatch]
//!                   (Serial)        (Parallel)                                                    (Serial)
//! ```
//!
//! Consumes `SortPhase2Event`s from the parallel decompress stage, drives
//! a non-blocking `MergeDriverDyn` over the slot table + memory chunks, and
//! emits sorted `RecordBatch`es downstream. The decompressed BGZF block
//! bytes flow out-of-band on each `SortMergeSlot`'s bounded queue (pushed by
//! `SortSpillDecompress`, popped by the merge driver) — only the slot
//! handles and counts ride the typed `SortPhase2Event` chain.
//!
//! See `docs/design/sort-step-split.md` for the locked design and
//! `docs/design/sort-merge-nonblocking.md` for the cooperative, resumable
//! merge consumer (issue #330 BUG #3).
//!
//! ## State machine
//!
//! 1. **`WaitingForSetup`** — collect upstream events. `SpillReady` events
//!    contribute their `slot` to the slot table (deduplicated by
//!    `file_id`); `MemoryChunk` events accumulate as typed memory chunks;
//!    `AllAnnounced` carries the final counts. We `max()`
//!    `records_ingested_so_far` from every event so `total_records` is
//!    order-independent. Transitions to `Merging` once the announced
//!    slot/chunk counts have all arrived (or, as a fallback, on input
//!    drain).
//!
//! 2. **`Merging`** — build a typed `MergeDriver<K>` via
//!    [`fgumi_sort::MergeDriver::from_slots`], type-erase to
//!    [`fgumi_sort::MergeDriverDyn`], and pump
//!    [`MergeDriverDyn::try_step`]. `try_step` is **non-blocking**: when a
//!    slot's decompressed queue is momentarily empty (and not at EOF) it
//!    returns [`fgumi_sort::MergeStep::Stalled`], and this step flushes any
//!    ready batch and returns `Contention` — yielding the worker so it can
//!    round-robin to `SortSpillDecompress` to refill the slot, then resume
//!    the merge on a later dispatch. Each produced record is appended to a
//!    [`RecordBatchBuilder`]; batches flush downstream when full
//!    (`target_batch_count` or `output_byte_limit`).
//!
//! 3. **`Done`** — terminal.
//!
//! ## Why the merge must never block
//!
//! `MergeDriver` runs on whichever framework worker holds `SortMerge`'s
//! Serial-step lock. A blocking wait there would park that worker (and hold
//! the lock) until a sibling worker refilled the slot — which a single-worker
//! configuration can never do, deadlocking the chain (issue #330 BUG #3,
//! review finding M1). The non-blocking `try_step` keeps the merge
//! cooperative: the same worker that stalls can refill the slot itself. Lazy
//! priming inside `from_slots` extends the same property to driver
//! construction, so building the `LoserTree` over not-yet-decompressed slots
//! never blocks either.

use std::collections::HashMap;
use std::io;
use std::sync::Arc;

use fgumi_raw_bam::RawRecord;
use fgumi_sort::{
    MergeDriver, MergeDriverDyn, MergeStep, QuerynameComparator, RawCoordinateKey, RawQuerynameKey,
    RawQuerynameLexKey, SortMergeSlot, SortOrder, TemplateKey,
};

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::sort::protocol::{MemoryChunkErased, SortPhase2Event};
use crate::pipeline::steps::types::{RecordBatch, RecordBatchBuilder};

/// Default output batch size: 1024 records per emitted `RecordBatch`.
/// Matches the legacy `Sort::DEFAULT_TARGET_BATCH_COUNT`.
pub const DEFAULT_TARGET_BATCH_COUNT: usize = 1024;

/// Max input events drained per `try_run` invocation in the
/// `WaitingForSetup` state. Bounds the per-call work so other Serial
/// steps in the chain make progress.
const MAX_EVENTS_PER_LOCK: usize = 16;

/// Max output batches emitted per `try_run` invocation in the `Merging`
/// state. Same intent as `MAX_EVENTS_PER_LOCK`.
const MAX_DRAIN_BATCHES_PER_LOCK: usize = 8;

/// Type-erased memory-chunk accumulator. Per-variant `Vec`s grow as
/// `MemoryChunk` events arrive in `WaitingForSetup`; on transition to
/// `Merging` exactly one variant is non-empty (matching `sort_order`)
/// and we hand its contents to the typed `MergeDriver::from_slots`.
///
/// We don't gate which variant fills purely by `sort_order` — the
/// upstream `SortAndSpill` step is single-typed too — but the four
/// fields exist so the dispatch is statically typed.
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

    /// Total non-empty memory chunks across all four variants.
    fn total_len(&self) -> usize {
        self.coordinate.len()
            + self.queryname_lex.len()
            + self.queryname_natural.len()
            + self.template_coordinate.len()
    }
}

/// Build a typed merge driver from the collected slots + memory chunks
/// for the configured sort order. Construction is infallible and always
/// yields a driver; an all-empty input is reported by the driver's first
/// `try_step` returning [`MergeStep::Done`].
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

/// Result of one [`SortMerge::next_batch`] pump.
enum NextBatch {
    /// A full output batch is ready to push; more may follow immediately.
    Batch(RecordBatch),
    /// The merge driver stalled on a not-yet-ready slot. Carries any
    /// partially-filled batch to flush before yielding.
    Stalled(Option<RecordBatch>),
    /// The merge is exhausted. Carries any final partial batch to flush, plus
    /// the driver's total `records_merged()` — read in `next_batch` while the
    /// driver is still in scope — for the completion log.
    Done(Option<RecordBatch>, u64),
}

/// Three-phase state machine for [`SortMerge`].
enum SortMergeState {
    WaitingForSetup {
        /// Slots received via `SpillReady` events. Sorted by `file_id`
        /// before passing to `MergeDriver::from_slots` so the
        /// `LoserTree` tie-break for equal sort keys is deterministic
        /// and matches the legacy chunk-files order.
        slots: Vec<Arc<SortMergeSlot>>,
        /// `slot.file_id -> index into slots` so deduplication is
        /// O(1) per `SpillReady`.
        slot_index: HashMap<u32, usize>,
        /// Typed accumulators for the four sort-key variants.
        memory_chunks: MemoryChunksByKind,
        /// Running `max(records_ingested_so_far)` across events.
        total_records: u64,
        /// Expected slot count from `AllAnnounced`. Set once when
        /// the sentinel arrives. The transition gate compares this
        /// against `slots.len()`; when both match (and
        /// `expected_memory_chunk_count` matches), we can transition
        /// to `Merging` without waiting for `ctx.input.is_drained()`.
        expected_slot_count: Option<u32>,
        /// Expected non-empty memory-chunk count from
        /// `AllAnnounced`. Compared against `memory_chunks.total_len()`.
        expected_memory_chunk_count: Option<u32>,
    },
    /// Active merge state. Boxed because `dyn MergeDriverDyn` is a wide
    /// pointer; `RecordBatchBuilder` carries a `Vec<u8>` allocation.
    Merging {
        driver: Box<dyn MergeDriverDyn + Send>,
        builder: RecordBatchBuilder,
        next_ordinal: u64,
    },
    Done,
}

/// Pop one `SortPhase2Event` and route it into the in-flight setup
/// accumulators. Factored out of `try_run`'s `WaitingForSetup` arm so
/// the latter stays under clippy's `too_many_lines` threshold.
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
            let inner = Arc::try_unwrap(chunk).unwrap_or_else(|arc| {
                // Fallback: clone. Only hit if a future change adds a
                // second `Arc` holder.
                match arc.as_ref() {
                    MemoryChunkErased::Coordinate(v) => MemoryChunkErased::Coordinate(v.clone()),
                    MemoryChunkErased::QuerynameLex(v) => {
                        MemoryChunkErased::QuerynameLex(v.clone())
                    }
                    MemoryChunkErased::QuerynameNatural(v) => {
                        MemoryChunkErased::QuerynameNatural(v.clone())
                    }
                    MemoryChunkErased::TemplateCoordinate(v) => {
                        MemoryChunkErased::TemplateCoordinate(v.clone())
                    }
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

/// Returns `true` if the slot set is complete: `AllAnnounced` has
/// arrived and the actual slot + memory-chunk counts match the
/// expected counts. Used to gate the early `WaitingForSetup → Merging`
/// transition without waiting for `ctx.input.is_drained()`.
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

/// `Serial + ByItemOrdinal` middle-of-chain merge. Input =
/// `SortPhase2Event`, Output = sorted `RecordBatch`. Drives
/// [`fgumi_sort::MergeDriver`] via pausable `peek` / `advance` so the
/// merge phase runs interleaved with the framework's round-robin
/// dispatch — no dedicated OS thread.
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

    /// Override the affinity hint. Not load-bearing for correctness — the
    /// pausable driver releases the Serial lock between bounded batches.
    #[must_use]
    pub fn with_affinity(mut self, affinity: Affinity) -> Self {
        self.affinity = affinity;
        self
    }

    /// Try to deliver a previously-held batch. Returns `true` if the held
    /// slot is now empty, `false` if it's still held (caller signals
    /// `Contention`).
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

    /// If currently in `WaitingForSetup`, sort the accumulated slots,
    /// build the merge driver, and transition to `Merging` (or
    /// directly to `Done` for empty input — `build_driver` returns
    /// `Ok(None)` when there's nothing to merge). No-op when already
    /// in `Merging` or `Done`. Called from `try_run` once the setup
    /// predicate is satisfied or the input is drained.
    ///
    /// `caller_label` is interpolated into the error message on
    /// `build_driver` failure so post-mortems can attribute the
    /// failure to the correct call site.
    /// `WaitingForSetup` Phase 1: pop events from upstream and absorb
    /// each into the setup accumulators. `cap` bounds the per-call
    /// absorb count; `try_run` passes `Some(MAX_EVENTS_PER_LOCK)` to keep
    /// each dispatch bounded. Returns the count of events absorbed so
    /// `try_run` can short-circuit on `did_work`.
    ///
    /// # Panics
    ///
    /// Panics if `self.state` is not `WaitingForSetup`. Callers must
    /// check the state before invoking.
    fn absorb_events_into_setup(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        cap: Option<usize>,
    ) -> usize {
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
        while cap.is_none_or(|c| absorbed < c) {
            let Some(event) = ctx.input.pop() else { break };
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

    /// True when `self.state` is `WaitingForSetup` AND the
    /// `AllAnnounced` predicates match. False in any other state.
    /// Used by `try_run` to early-transition `WaitingForSetup ->
    /// Merging` without waiting for upstream's drained flag.
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

    /// Pump the merge driver into the current batch builder via the
    /// non-blocking [`MergeDriverDyn::try_step`] until the builder is
    /// full (a `Batch` is ready), the driver stalls on a not-yet-ready
    /// slot (`Stalled`), or the merge is exhausted (`Done`). For the
    /// latter two, any partially-filled builder is flushed alongside so
    /// the caller can push it before yielding / finishing.
    ///
    /// # Panics
    ///
    /// Panics if `self.state` is not `Merging`. Callers must check
    /// the state before invoking.
    fn next_batch(&mut self) -> io::Result<NextBatch> {
        let target = self.target_batch_count;
        let byte_limit = self.output_byte_limit;
        let bytes_cap = usize::try_from(byte_limit).unwrap_or(usize::MAX);
        let SortMergeState::Merging { driver, builder, next_ordinal } = &mut self.state else {
            unreachable!("next_batch called outside Merging state");
        };

        // Flush the current builder into a `RecordBatch`, advancing the
        // ordinal and installing a fresh builder.
        let flush = |builder: &mut RecordBatchBuilder, next_ordinal: &mut u64| {
            *next_ordinal += 1;
            let next_builder = RecordBatchBuilder::with_capacity(*next_ordinal, bytes_cap, target);
            std::mem::replace(builder, next_builder).build()
        };
        // Flush only if the builder holds at least one record.
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

    /// Cooperative merge emit. Produce up to `MAX_DRAIN_BATCHES_PER_LOCK`
    /// batches via `next_batch`, push each downstream, stash the first
    /// rejection in `held` and return `Progress`. Transitions state to
    /// `Done` when the driver is exhausted.
    ///
    /// Called from `try_run` only.
    ///
    /// # Panics
    ///
    /// Panics if `self.state` is not `Merging`.
    fn emit_batches_cooperative(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let mut delivered = 0usize;
        // Push a batch downstream; on backpressure stash it in `held` and
        // signal the caller to return `Progress` (so `flush_held` retries it
        // next dispatch). Returns `false` if the push was held.
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
                    // A slot's queue is momentarily empty. Yield so this
                    // worker round-robins to `SortSpillDecompress` to refill
                    // it, then resume the merge on a later dispatch. Either
                    // outcome yields: `Progress` (we delivered batches this
                    // call) triggers the framework's priority-restart at step
                    // 0, which subsumes the `Contention` retry — so reporting
                    // `Progress` when `delivered > 0` is the honest, stronger
                    // signal, and `Contention` covers the made-no-progress case.
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
                            // Stay in `Merging`; next dispatch flushes `held`,
                            // then `next_batch` returns `Done(None)` → `Done`.
                            return Ok(StepOutcome::Progress);
                        }
                        delivered += 1;
                    }
                    // Surface the merged-record count carried up from `next_batch`
                    // (where the driver was in scope). The streaming sort's
                    // analogue of standalone sort's Records-written summary.
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

    /// If currently `WaitingForSetup`, build the merge driver and move to
    /// `Merging`. No-op in any other state. Infallible: `build_driver` always
    /// yields a driver, and empty input is reported by the driver's first
    /// `try_step` returning `Done` (handled by the `Merging` emit paths), so
    /// there is no separate empty-input or error branch.
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
        // Sort slots by `file_id` so `MergeDriver::from_slots` source
        // order is deterministic and matches the legacy chunk-files
        // order. `SpillReady` events arrive in `SortAndSpill`'s
        // emission order (which IS `file_id` order today), but
        // defensively sort to avoid relying on that.
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

        // Phase 1: absorb events while in WaitingForSetup. Transition
        // to Merging when AllAnnounced predicates match (early-
        // transition gate — avoids waiting for upstream's drained flag
        // which would require SortSpillDecompress workers to Skip
        // themselves first, blocking further decompress work) OR
        // when upstream is drained with zero slots/chunks (empty-
        // input fall-back).
        if matches!(&self.state, SortMergeState::WaitingForSetup { .. }) {
            let absorbed = self.absorb_events_into_setup(ctx, Some(MAX_EVENTS_PER_LOCK));
            if !self.is_ready_to_merge() {
                if absorbed > 0 {
                    return Ok(StepOutcome::Progress);
                }
                if !ctx.input.is_drained() {
                    return Ok(StepOutcome::NoProgress);
                }
            }
            self.transition_to_merging();
        }

        // Phase 2: cooperative emit if Merging; Done reports `Finished`.
        // `Done` means the merge driver is exhausted and every batch has
        // landed — completion is keyed on the merge state, NOT on
        // `ctx.input.is_drained()` (the input edge drains early, before the
        // merge finishes, because the records come from spilled chunks /
        // sibling decompress workers; finishing on `is_drained` would
        // truncate the sorted output). Empty input reaches `Merging`, then
        // `emit_batches_cooperative`'s first `next_batch` returns `Done(None)`
        // → state `Done` → `Finished` on the next pass. WaitingForSetup is
        // unreachable — Phase 1 either returned early or transitioned us out.
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
