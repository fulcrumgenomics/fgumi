//! `SortMerge` — third step of the runall-sort three-step chain.

use std::collections::HashMap;
use std::io;
use std::sync::Arc;

use fgumi_sort::{
    CbKey32, InMemoryChunk, MemorySources, MergeDriver, MergeDriverDyn, MergeStep,
    QuerynameComparator, RawCoordinateKey, RawQuerynameKey, RawQuerynameLexKey, SortMergeSlot,
    SortOrder, TemplateKey, TemplateKey24, TemplateMemChunk, TertKey32,
};

use crate::sort::protocol::{MemoryChunkErased, SortPhase2Event};
use crate::types::{DecompressedBlock, RecordBatch, RecordBatchBuilder};
use fgumi_pipeline_core::{
    HeapSize, Ordered, Unpushed,
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

/// Initial reservation for an output-batch byte buffer, before any batch has
/// been emitted to size the next one from. Kept modest on purpose: most batches
/// fill on the record-count cap well below the output-queue byte budget, so
/// reserving the full budget for every buffer chronically over-allocates (and
/// inflates the byte-bounded queue's capacity-based accounting). Buffers grow
/// on demand via `extend_from_slice`, so under-reserving only costs a few
/// startup reallocations.
const INITIAL_OUTPUT_BUFFER_BYTES: usize = 64 * 1024;

// ─────────────────────────────────────────────────────────────────────────────
// MergeOutput — the framing strategy the merge accumulates winners into.
// ─────────────────────────────────────────────────────────────────────────────

/// The output-framing strategy `SortMerge` accumulates merged winner records
/// into. The merge state machine, `LoserTree` driver, source ordering and
/// tie-break are identical for every strategy; only the per-record framing and
/// the emitted item type differ.
///
/// Two implementations exist:
///
/// - [`RecordBatchOutput`] (the default) — accumulates raw record bodies into a
///   [`RecordBatch`] (flat backing buffer + per-record `(start, end)` ranges).
///   This is the **intermediate** sort output, consumed by `DecodeFromRecords`
///   downstream in a fused `--start-from sort` chain.
/// - [`BlockOutput`] — frames each record as `[u32 LE block_size][body]`
///   directly into a [`DecompressedBlock`], byte-for-byte identical to the
///   `SerializeRecordBatch` step it replaces (lever 1). This is the
///   **terminal** standalone-sort output, wired straight to `BgzfCompress`,
///   folding the former `SortMerge → SerializeRecordBatch → BgzfCompress`
///   triple into `SortMerge → BgzfCompress` (one fewer pool step, one fewer
///   reorder stage, one fewer memcpy per record).
pub trait MergeOutput: Send + 'static {
    /// The emitted batch item type.
    type Item: Send + HeapSize + Ordered + 'static;
    /// The per-batch accumulator.
    type Builder: MergeBatchBuilder<Item = Self::Item>;
}

/// A per-batch accumulator for a [`MergeOutput`] strategy. Mirrors the
/// [`RecordBatchBuilder`] surface the merge loop already drives, so the merge
/// state machine is strategy-agnostic.
pub trait MergeBatchBuilder: Send + 'static {
    /// The finalized batch item this builder produces.
    type Item;

    /// Create a builder for batch `batch_serial`, reserving `bytes_cap` bytes of
    /// payload and room for `records_cap` records.
    fn with_capacity(batch_serial: u64, bytes_cap: usize, records_cap: usize) -> Self;

    /// Append one merged winner record's raw BAM body.
    ///
    /// # Errors
    ///
    /// Returns an error if the record cannot be framed (e.g. a body whose
    /// length does not fit the strategy's length prefix).
    fn push_record_bytes(&mut self, body: &[u8]) -> io::Result<()>;

    /// Number of records appended so far.
    fn len(&self) -> usize;

    /// `true` iff no records have been appended.
    fn is_empty(&self) -> bool;

    /// Total payload bytes accumulated so far (used to size the next buffer and
    /// to enforce the per-batch byte cap).
    fn total_bytes(&self) -> usize;

    /// Finalize and produce the batch item, consuming the builder.
    fn build(self) -> Self::Item;
}

/// Intermediate-sort output: raw record bodies into a [`RecordBatch`].
pub struct RecordBatchOutput;

impl MergeOutput for RecordBatchOutput {
    type Item = RecordBatch;
    type Builder = RecordBatchBuilder;
}

impl MergeBatchBuilder for RecordBatchBuilder {
    type Item = RecordBatch;

    fn with_capacity(batch_serial: u64, bytes_cap: usize, records_cap: usize) -> Self {
        RecordBatchBuilder::with_capacity(batch_serial, bytes_cap, records_cap)
    }

    fn push_record_bytes(&mut self, body: &[u8]) -> io::Result<()> {
        RecordBatchBuilder::push_record_bytes(self, body);
        Ok(())
    }

    fn len(&self) -> usize {
        RecordBatchBuilder::len(self)
    }

    fn is_empty(&self) -> bool {
        RecordBatchBuilder::is_empty(self)
    }

    fn total_bytes(&self) -> usize {
        RecordBatchBuilder::total_bytes(self)
    }

    fn build(self) -> RecordBatch {
        RecordBatchBuilder::build(self)
    }
}

/// Terminal standalone-sort output: each record framed as
/// `[u32 LE block_size][body]` directly into a [`DecompressedBlock`], ready for
/// `BgzfCompress`. This is byte-for-byte identical to the framing the former
/// `SerializeRecordBatch` step produced (lever 1).
pub struct BlockOutput;

impl MergeOutput for BlockOutput {
    type Item = DecompressedBlock;
    type Builder = BlockBuilder;
}

/// Accumulates merged winner records as BAM on-disk framing
/// (`[u32 LE block_size][body]` per record) into a single byte buffer that
/// becomes a [`DecompressedBlock`]. This is the canonical BAM record layout;
/// the `fgumi` crate's `serialize::frame_record_into` is the sibling
/// implementation (a separate crate, so the two cannot share code) and the two
/// MUST stay byte-for-byte in sync — each has a layout test pinning it.
pub struct BlockBuilder {
    batch_serial: u64,
    bytes: Vec<u8>,
    /// Record count — tracked separately because the framed byte buffer mixes
    /// length prefixes with bodies, so it cannot be recovered from `bytes`.
    records: usize,
}

impl MergeBatchBuilder for BlockBuilder {
    type Item = DecompressedBlock;

    fn with_capacity(batch_serial: u64, bytes_cap: usize, _records_cap: usize) -> Self {
        // `_records_cap` sizes the `RecordBatch` ranges vector; the framed-block
        // builder has no separate per-record allocation to reserve.
        Self { batch_serial, bytes: Vec::with_capacity(bytes_cap), records: 0 }
    }

    fn push_record_bytes(&mut self, body: &[u8]) -> io::Result<()> {
        // `[u32 LE block_size][body]`, byte-identical to
        // `SerializeRecordBatch::frame_record_into`.
        let block_size = u32::try_from(body.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("record exceeds u32 BAM block_size: {}", body.len()),
            )
        })?;
        self.bytes.extend_from_slice(&block_size.to_le_bytes());
        self.bytes.extend_from_slice(body);
        self.records += 1;
        Ok(())
    }

    fn len(&self) -> usize {
        self.records
    }

    fn is_empty(&self) -> bool {
        self.records == 0
    }

    fn total_bytes(&self) -> usize {
        self.bytes.len()
    }

    fn build(self) -> DecompressedBlock {
        DecompressedBlock { batch_serial: self.batch_serial, bytes: self.bytes }
    }
}

/// Same-variant collector for template-coordinate residual chunks.
///
/// A single sort chooses its `--key-types` narrowed lane variant exactly once
/// (globally, on the first record) and reuses it for every run, so all template
/// chunks in one merge share one arm. The first push fixes the arm; subsequent
/// pushes assert-match it (a variant change mid-sort is impossible by
/// construction and would be a bug).
#[derive(Default)]
enum TemplateChunks {
    /// No template chunks pushed yet — the variant is not yet known.
    #[default]
    Empty,
    /// 24-byte core-only lane.
    K24(Vec<InMemoryChunk<TemplateKey24>>),
    /// 32-byte lane carrying `cb_hash`.
    Cb32(Vec<InMemoryChunk<CbKey32>>),
    /// 32-byte lane carrying the tertiary word.
    Tert32(Vec<InMemoryChunk<TertKey32>>),
    /// Full 40-byte key (all lanes) — the legacy owned path and full variant.
    K40(Vec<InMemoryChunk<TemplateKey>>),
}

impl TemplateChunks {
    fn push(&mut self, chunk: TemplateMemChunk) {
        match chunk {
            TemplateMemChunk::K24(c) => match self {
                Self::Empty => *self = Self::K24(vec![c]),
                Self::K24(v) => v.push(c),
                _ => unreachable!("template chunk variant changed mid-sort (variant is global)"),
            },
            TemplateMemChunk::Cb32(c) => match self {
                Self::Empty => *self = Self::Cb32(vec![c]),
                Self::Cb32(v) => v.push(c),
                _ => unreachable!("template chunk variant changed mid-sort (variant is global)"),
            },
            TemplateMemChunk::Tert32(c) => match self {
                Self::Empty => *self = Self::Tert32(vec![c]),
                Self::Tert32(v) => v.push(c),
                _ => unreachable!("template chunk variant changed mid-sort (variant is global)"),
            },
            TemplateMemChunk::K40(c) => match self {
                Self::Empty => *self = Self::K40(vec![c]),
                Self::K40(v) => v.push(c),
                _ => unreachable!("template chunk variant changed mid-sort (variant is global)"),
            },
        }
    }

    fn len(&self) -> usize {
        match self {
            Self::Empty => 0,
            Self::K24(v) => v.len(),
            Self::Cb32(v) => v.len(),
            Self::Tert32(v) => v.len(),
            Self::K40(v) => v.len(),
        }
    }

    /// Pop the sole chunk (caller guarantees exactly one) and re-erase it.
    fn pop_single(self) -> TemplateMemChunk {
        match self {
            Self::K24(mut v) => TemplateMemChunk::K24(v.pop().expect("one chunk")),
            Self::Cb32(mut v) => TemplateMemChunk::Cb32(v.pop().expect("one chunk")),
            Self::Tert32(mut v) => TemplateMemChunk::Tert32(v.pop().expect("one chunk")),
            Self::K40(mut v) => TemplateMemChunk::K40(v.pop().expect("one chunk")),
            Self::Empty => unreachable!("pop_single called with no chunk"),
        }
    }
}

#[derive(Default)]
struct MemoryChunksByKind {
    coordinate: Vec<InMemoryChunk<RawCoordinateKey>>,
    queryname_lex: Vec<InMemoryChunk<RawQuerynameLexKey>>,
    queryname_natural: Vec<InMemoryChunk<RawQuerynameKey>>,
    template_coordinate: TemplateChunks,
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

    /// Consume the single chunk held across all kinds, re-erased.
    ///
    /// # Panics
    ///
    /// Panics if `total_len() != 1` (the single-source fast path's precondition).
    fn into_single(mut self) -> MemoryChunkErased {
        debug_assert_eq!(self.total_len(), 1, "into_single requires exactly one chunk");
        if let Some(c) = self.coordinate.pop() {
            MemoryChunkErased::Coordinate(c)
        } else if let Some(c) = self.queryname_lex.pop() {
            MemoryChunkErased::QuerynameLex(c)
        } else if let Some(c) = self.queryname_natural.pop() {
            MemoryChunkErased::QuerynameNatural(c)
        } else if self.template_coordinate.len() == 1 {
            MemoryChunkErased::TemplateCoordinate(
                std::mem::take(&mut self.template_coordinate).pop_single(),
            )
        } else {
            unreachable!("into_single called with no chunk")
        }
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
            MemorySources::Shared(chunks.coordinate),
            total_records,
        )),
        SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
            Box::new(MergeDriver::<RawQuerynameLexKey>::from_slots(
                slots,
                MemorySources::Shared(chunks.queryname_lex),
                total_records,
            ))
        }
        SortOrder::Queryname(QuerynameComparator::Natural) => {
            Box::new(MergeDriver::<RawQuerynameKey>::from_slots(
                slots,
                MemorySources::Shared(chunks.queryname_natural),
                total_records,
            ))
        }
        SortOrder::TemplateCoordinate => match chunks.template_coordinate {
            // Empty (no residual chunks) only for empty input, where no spill
            // files of any other variant exist either — so any K is safe.
            TemplateChunks::Empty => Box::new(MergeDriver::<TemplateKey>::from_slots(
                slots,
                MemorySources::Shared(Vec::new()),
                total_records,
            )),
            TemplateChunks::K24(v) => Box::new(MergeDriver::<TemplateKey24>::from_slots(
                slots,
                MemorySources::Shared(v),
                total_records,
            )),
            TemplateChunks::Cb32(v) => Box::new(MergeDriver::<CbKey32>::from_slots(
                slots,
                MemorySources::Shared(v),
                total_records,
            )),
            TemplateChunks::Tert32(v) => Box::new(MergeDriver::<TertKey32>::from_slots(
                slots,
                MemorySources::Shared(v),
                total_records,
            )),
            TemplateChunks::K40(v) => Box::new(MergeDriver::<TemplateKey>::from_slots(
                slots,
                MemorySources::Shared(v),
                total_records,
            )),
        },
    }
}

enum NextBatch<I> {
    Batch(I),
    Stalled(Option<I>),
    Done(Option<I>, u64),
}

enum SortMergeState<B> {
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
        builder: B,
        next_ordinal: u64,
    },
    /// Single-source fast path: 0 spill slots and exactly one in-memory chunk, so
    /// the chunk is already globally sorted and no k-way merge is needed. Gather
    /// its records in order straight into output blocks (the dominant in-memory
    /// cost — a single-threaded k = 1 loser-tree walk — is pure overhead here).
    FastPath {
        chunk: MemoryChunkErased,
        /// Index of the next record to gather.
        cursor: usize,
        /// Total record count (`chunk.len()`), cached to avoid re-dispatching.
        total: usize,
        builder: B,
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
) -> io::Result<()> {
    match event {
        SortPhase2Event::SpillReady { slot, path: _, records_ingested_so_far } => {
            if let std::collections::hash_map::Entry::Vacant(e) = slot_index.entry(slot.file_id) {
                e.insert(slots.len());
                slots.push(slot);
            }
            *total_records = (*total_records).max(records_ingested_so_far);
        }
        SortPhase2Event::MemoryChunk { chunk, records_ingested_so_far } => {
            // The MemoryChunk `Arc` is created uniquely in `SortAndSpill` and only
            // ever *moved* (never cloned) through `SortSpillDecompress` to here, so
            // it must be uniquely owned at this single consumer. Fail closed on a
            // shared `Arc` rather than silently deep-cloning a potentially large
            // record vector on the merge setup path.
            let inner = Arc::try_unwrap(chunk).map_err(|_| {
                io::Error::other(
                    "SortMerge: MemoryChunk Arc unexpectedly shared at the merge consumer \
                     (protocol invariant: memory chunks are moved, never cloned)",
                )
            })?;
            memory_chunks.push(inner);
            *total_records = (*total_records).max(records_ingested_so_far);
        }
        SortPhase2Event::AllAnnounced {
            slot_count,
            memory_chunk_count,
            total_records: ar_total,
        } => {
            // The Phase-2 protocol emits exactly one `AllAnnounced` (the last event
            // from `SortAndSpill`). A second one is a protocol violation; fail
            // closed rather than overwrite the prior expectations and risk masking
            // the bug behind a silently-different completion target.
            if expected_slot_count.is_some() || expected_memory_chunk_count.is_some() {
                return Err(io::Error::other(format!(
                    "SortMerge: duplicate AllAnnounced — prior {expected_slot_count:?}/\
                     {expected_memory_chunk_count:?}, new {slot_count}/{memory_chunk_count}; \
                     the Phase-2 protocol emits exactly one AllAnnounced",
                )));
            }
            *expected_slot_count = Some(slot_count);
            *expected_memory_chunk_count = Some(memory_chunk_count);
            *total_records = (*total_records).max(ar_total);
        }
    }
    Ok(())
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

/// `Detached + ByItemOrdinal` terminal merge: the final of the three sort
/// steps, producing the sorted output stream consumed by the sink.
///
/// Generic over the output-framing strategy `O` (see [`MergeOutput`]):
/// [`RecordBatchOutput`] (the default) emits [`RecordBatch`] for a fused
/// intermediate sort, and [`BlockOutput`] frames records directly into
/// [`DecompressedBlock`]s for the standalone-sort terminal so the chain can
/// wire `SortMerge → BgzfCompress` with no intervening serialize step
/// (lever 1). The merge core — `LoserTree` driver, source ordering, tie-break,
/// cooperative `try_run` body — is identical for both.
pub struct SortMerge<O: MergeOutput = RecordBatchOutput> {
    state: SortMergeState<O::Builder>,
    held: HeldSlot<Unpushed<O::Item>>,
    sort_order: SortOrder,
    target_batch_count: usize,
    output_byte_limit: u64,
    affinity: Affinity,
    /// Optional sink for the end-of-run sort summary, filled when the merge
    /// reaches `Done`. The standalone-sort summary finalize hook reads it to
    /// log records processed/written and the spill-chunk count.
    stats_slot: Option<Arc<std::sync::Mutex<Option<fgumi_sort::SortStats>>>>,
    /// Total records ingested, captured at the merge transition (the summary's
    /// "records processed").
    processed: u64,
    /// Number of spill files, captured at the merge transition (the summary's
    /// "temporary chunks"). Zero for a fully in-memory sort.
    chunk_count: usize,
    /// INSTRUMENTATION (lever-2 merge-stall diagnosis; `RUST_LOG=info` at Done).
    /// `SortMerge` runs on a single dedicated `Detached` thread (one instance,
    /// never `new_worker_copy`'d), so plain `&mut self` counters are sound — no
    /// atomics needed.
    dbg: MergeDiag,
}

/// Lever-2 diagnostic counters: is the serial merge starved on decompress
/// (`input-empty`/`stalls`) or blocked on the downstream writer
/// (`output_full`), and how much does its worker spin (`contention`)?
#[derive(Default, Clone, Copy)]
struct MergeDiag {
    /// Merge-loop passes that ended `Stalled` — the winning source's next block
    /// was not yet decompressed (INPUT-STARVED: the lever-2 hypothesis).
    stalls: u64,
    /// `ctx.outputs.push` returned `Err` — downstream (compress/write) full
    /// (OUTPUT-BACKPRESSURE).
    output_full: u64,
    /// `try_run` returned `Contention` — the merge worker had nothing to do this
    /// dispatch and spun/yielded (pure under-utilization).
    contention: u64,
    /// `try_run` calls that delivered ≥1 batch (PROGRESS dispatches).
    progress_dispatches: u64,
}

impl<O: MergeOutput> SortMerge<O> {
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
            stats_slot: None,
            processed: 0,
            chunk_count: 0,
            dbg: MergeDiag::default(),
        }
    }

    /// Override the affinity hint.
    #[must_use]
    pub fn with_affinity(mut self, affinity: Affinity) -> Self {
        self.affinity = affinity;
        self
    }

    /// Attach a slot to receive the end-of-run [`fgumi_sort::SortStats`] when
    /// the merge completes (records processed/written + spill-chunk count). Used
    /// by the standalone-sort summary finalize hook; runall leaves it unset.
    #[must_use]
    pub fn with_stats_slot(
        mut self,
        slot: Arc<std::sync::Mutex<Option<fgumi_sort::SortStats>>>,
    ) -> Self {
        self.stats_slot = Some(slot);
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
    ///
    /// # Errors
    ///
    /// Returns an error on a Phase-2 protocol violation (a duplicate
    /// `AllAnnounced`, or a `MemoryChunk` whose `Arc` is unexpectedly shared).
    fn absorb_events_into_setup(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<usize> {
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
            )?;
            absorbed += 1;
        }
        Ok(absorbed)
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
    fn next_batch(&mut self) -> io::Result<NextBatch<O::Item>> {
        let target = self.target_batch_count;
        let byte_limit = self.output_byte_limit;
        let bytes_cap = usize::try_from(byte_limit).unwrap_or(usize::MAX);
        let SortMergeState::Merging { driver, builder, next_ordinal } = &mut self.state else {
            unreachable!("next_batch called outside Merging state");
        };

        let buffer_floor = INITIAL_OUTPUT_BUFFER_BYTES.min(bytes_cap);
        let flush = |builder: &mut O::Builder, next_ordinal: &mut u64| {
            *next_ordinal += 1;
            // Size the next buffer to the batch we just filled, clamped to
            // `[buffer_floor, bytes_cap]`. Count-bound batches stay small; a
            // byte-bound batch carries ~`bytes_cap` forward. This avoids
            // reserving the full byte budget for every (typically count-bound)
            // batch — see `INITIAL_OUTPUT_BUFFER_BYTES`.
            let hint = builder.total_bytes().clamp(buffer_floor, bytes_cap);
            let next_builder = O::Builder::with_capacity(*next_ordinal, hint, target);
            std::mem::replace(builder, next_builder).build()
        };
        let flush_partial = |builder: &mut O::Builder, next_ordinal: &mut u64| {
            if builder.is_empty() { None } else { Some(flush(builder, next_ordinal)) }
        };

        loop {
            match driver
                .try_step()
                .map_err(|e| io::Error::other(format!("SortMerge: merge step failed: {e:#}")))?
            {
                MergeStep::Produced(bytes) => {
                    builder.push_record_bytes(bytes)?;
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
                        self.dbg.output_full += 1;
                        self.held.put(unpushed);
                        return Ok(StepOutcome::Progress);
                    }
                    delivered += 1;
                    if delivered >= MAX_DRAIN_BATCHES_PER_LOCK {
                        self.dbg.progress_dispatches += 1;
                        return Ok(StepOutcome::Progress);
                    }
                }
                NextBatch::Stalled(partial) => {
                    // INPUT-STARVED: the driver couldn't advance because the
                    // winning source's next block isn't decompressed yet.
                    self.dbg.stalls += 1;
                    if let Some(batch) = partial {
                        if let Err(unpushed) = ctx.outputs.push(batch) {
                            self.dbg.output_full += 1;
                            self.held.put(unpushed);
                            return Ok(StepOutcome::Progress);
                        }
                        delivered += 1;
                    }
                    return Ok(if delivered > 0 {
                        self.dbg.progress_dispatches += 1;
                        StepOutcome::Progress
                    } else {
                        // Pure under-utilization: this dispatch did nothing.
                        self.dbg.contention += 1;
                        StepOutcome::Contention
                    });
                }
                NextBatch::Done(partial, merged) => {
                    if let Some(batch) = partial {
                        if let Err(unpushed) = ctx.outputs.push(batch) {
                            self.dbg.output_full += 1;
                            self.held.put(unpushed);
                            return Ok(StepOutcome::Progress);
                        }
                        delivered += 1;
                    }
                    log::info!("Sort merge complete: {merged} records merged");
                    // INSTRUMENTATION (lever-2): is the serial merge starved on
                    // decompress (stalls/contention high) or blocked on the
                    // writer (output_full high)? `stalls` counts merge-loop
                    // passes that ended input-starved; `contention` counts
                    // dispatches that produced nothing (pure idle spin);
                    // `output_full` counts downstream-backpressure events.
                    let d = self.dbg;
                    log::info!(
                        "Sort merge diag: stalls={} contention={} output_full={} \
                         progress_dispatches={} ({} records, {} sources)",
                        d.stalls,
                        d.contention,
                        d.output_full,
                        d.progress_dispatches,
                        merged,
                        self.chunk_count,
                    );
                    if let Some(slot) = &self.stats_slot {
                        *slot.lock().expect("sort stats slot poisoned") =
                            Some(fgumi_sort::SortStats {
                                total_records: self.processed,
                                output_records: merged,
                                chunks_written: self.chunk_count,
                            });
                    }
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

    /// Build the next output batch for the single-source fast path: gather records
    /// from the sorted chunk into the builder until the count/byte cap, or `Done`
    /// when the chunk is exhausted. Mirrors [`next_batch`](Self::next_batch)'s
    /// framing and buffer-sizing exactly, so the output is byte-identical to a
    /// (k = 1) loser-tree merge of the same chunk — only the record SOURCE differs
    /// (a direct cursor instead of `driver.try_step()`).
    ///
    /// # Panics
    ///
    /// Panics if `self.state` is not `FastPath`.
    fn next_fast_batch(&mut self) -> io::Result<NextBatch<O::Item>> {
        let target = self.target_batch_count;
        let byte_limit = self.output_byte_limit;
        let bytes_cap = usize::try_from(byte_limit).unwrap_or(usize::MAX);
        let buffer_floor = INITIAL_OUTPUT_BUFFER_BYTES.min(bytes_cap);
        let SortMergeState::FastPath { chunk, cursor, total, builder, next_ordinal } =
            &mut self.state
        else {
            unreachable!("next_fast_batch called outside FastPath state");
        };

        let flush = |builder: &mut O::Builder, next_ordinal: &mut u64| {
            *next_ordinal += 1;
            let hint = builder.total_bytes().clamp(buffer_floor, bytes_cap);
            let next_builder = O::Builder::with_capacity(*next_ordinal, hint, target);
            std::mem::replace(builder, next_builder).build()
        };
        let flush_partial = |builder: &mut O::Builder, next_ordinal: &mut u64| {
            if builder.is_empty() { None } else { Some(flush(builder, next_ordinal)) }
        };

        loop {
            if *cursor >= *total {
                return Ok(NextBatch::Done(flush_partial(builder, next_ordinal), *total as u64));
            }
            builder.push_record_bytes(chunk.record_bytes(*cursor))?;
            *cursor += 1;
            let count_full = builder.len() >= target;
            let bytes_full = (builder.total_bytes() as u64) >= byte_limit;
            if count_full || bytes_full {
                return Ok(NextBatch::Batch(flush(builder, next_ordinal)));
            }
        }
    }

    /// Cooperative emit loop for the single-source fast path. Mirrors
    /// [`emit_batches_cooperative`](Self::emit_batches_cooperative) but never
    /// `Stalled` (every record is already in memory).
    ///
    /// # Panics
    ///
    /// Panics if `self.state` is not `FastPath`.
    fn emit_fast_batches(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let mut delivered = 0usize;
        loop {
            match self.next_fast_batch()? {
                NextBatch::Batch(batch) => {
                    if let Err(unpushed) = ctx.outputs.push(batch) {
                        self.dbg.output_full += 1;
                        self.held.put(unpushed);
                        return Ok(StepOutcome::Progress);
                    }
                    delivered += 1;
                    if delivered >= MAX_DRAIN_BATCHES_PER_LOCK {
                        self.dbg.progress_dispatches += 1;
                        return Ok(StepOutcome::Progress);
                    }
                }
                NextBatch::Stalled(_) => unreachable!("FastPath never stalls (all in memory)"),
                NextBatch::Done(partial, count) => {
                    if let Some(batch) = partial {
                        if let Err(unpushed) = ctx.outputs.push(batch) {
                            self.dbg.output_full += 1;
                            self.held.put(unpushed);
                            return Ok(StepOutcome::Progress);
                        }
                        delivered += 1;
                    }
                    log::info!(
                        "Sort in-memory fast path complete: {count} records (single source, \
                         no merge)"
                    );
                    if let Some(slot) = &self.stats_slot {
                        *slot.lock().expect("sort stats slot poisoned") =
                            Some(fgumi_sort::SortStats {
                                total_records: self.processed,
                                output_records: count,
                                chunks_written: 0,
                            });
                    }
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
        // Capture summary inputs before `slots` is consumed by the driver:
        // total records ingested and the spill-file count.
        self.processed = total_records;
        self.chunk_count = slots.len();
        let bytes_cap_for_init = usize::try_from(self.output_byte_limit).unwrap_or(usize::MAX);
        let initial_bytes_for_init = INITIAL_OUTPUT_BUFFER_BYTES.min(bytes_cap_for_init);
        // FAST PATH: zero spill slots + exactly one in-memory chunk → the chunk is
        // already globally sorted, so skip the (k = 1) loser-tree merge and gather
        // it directly. This is the in-memory regime's dominant cost.
        if slots.is_empty() && memory_chunks.total_len() == 1 {
            let chunk = memory_chunks.into_single();
            let total = chunk.len();
            let builder =
                O::Builder::with_capacity(0, initial_bytes_for_init, self.target_batch_count);
            self.state =
                SortMergeState::FastPath { chunk, cursor: 0, total, builder, next_ordinal: 0 };
            return;
        }
        let driver = build_driver(self.sort_order, slots, memory_chunks, total_records);
        let bytes_cap = usize::try_from(self.output_byte_limit).unwrap_or(usize::MAX);
        // Seed the first buffer modestly; subsequent buffers are sized from the
        // prior batch's actual byte length (see `next_batch`).
        let initial_bytes = INITIAL_OUTPUT_BUFFER_BYTES.min(bytes_cap);
        let builder = O::Builder::with_capacity(0, initial_bytes, self.target_batch_count);
        self.state = SortMergeState::Merging { driver, builder, next_ordinal: 0 };
    }
}

impl<O: MergeOutput> Step for SortMerge<O> {
    type Input = SortPhase2Event;
    type Outputs = OrderedBytesSingle<O::Item>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SortMerge",
            // Lever 2: the merge runs on its own dedicated OS thread (off the
            // work-stealing pool), mirroring the legacy sort's main-thread
            // merge in the "N + 2" model. Its cooperative `try_run` body is
            // UNCHANGED — `run_detached_step` drives it with a private
            // blocking-style loop, parking on `Contention`/`NoProgress`
            // (winner-slot momentarily empty / output full) instead of the
            // pool re-dispatching it. `Detached` collapses the declared
            // `ByItemOrdinal` output to `None` exactly as `Serial` did (see
            // `effective_branch_orderings`), so the output transport — a direct
            // byte-bounded queue, no reorder stage — is byte-for-byte identical
            // to the pre-lever-2 Serial merge; the LoserTree core, source order
            // (`slots.sort_by_key(file_id)` + residual last), and tie-break are
            // untouched. SortMerge is only ever built by the sort chain's
            // `add_sort`, so this is sort-chain-only by construction.
            kind: StepKind::Detached,
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
            let absorbed = self.absorb_events_into_setup(ctx)?;
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
            SortMergeState::FastPath { .. } => self.emit_fast_batches(ctx),
            SortMergeState::Done => Ok(StepOutcome::Finished),
            SortMergeState::WaitingForSetup { .. } => {
                unreachable!("Phase 1 must have left state non-WaitingForSetup")
            }
        }
    }
}

#[cfg(test)]
mod tests;
