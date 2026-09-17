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
    HeapSize, HeldRetry, Ordered, Unpushed,
    held::HeldSlot,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{CounterSpec, DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// Default output batch size: 1024 records per emitted `RecordBatch`.
pub const DEFAULT_TARGET_BATCH_COUNT: usize = 1024;

/// Max output batches emitted per `try_run` invocation in `Merging`.
const MAX_DRAIN_BATCHES_PER_LOCK: usize = 8;

/// Minimum record count for the fast path to fan its gather across the pool.
/// Below this the serial gather is faster: the rayon pool build, block plan, and
/// per-window join overhead would dominate a small chunk (which the serial path
/// clears in well under a millisecond). ~1 output block per worker at the
/// default 1024-record batch, so the fan-out always has real work to spread.
const FAST_PATH_PARALLEL_MIN_RECORDS: usize = 64 * 1024;

/// Output blocks gathered per parallel window per worker. The parallel fast path
/// gathers `fast_path_threads * this` blocks per burst, pushes them in order,
/// then repeats — bounding peak extra memory to that many uncompressed blocks
/// (vs. materialising the whole output) while keeping every worker fed and the
/// `try_run` body cooperative.
const FAST_PATH_BLOCKS_PER_WORKER_WINDOW: usize = 4;

/// `SortMerge` counter slot index: records merged/gathered this call.
const RECORDS: usize = 0;

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

    /// Bytes this builder's [`total_bytes`](Self::total_bytes) grows by, PER
    /// record, beyond the record body itself — i.e. the per-record framing
    /// overhead the byte-cap accounting must include. [`BlockBuilder`] frames
    /// each record as `[u32 LE block_size][body]`, so it is `4`;
    /// [`RecordBatchBuilder`] stores bare bodies in a flat backing buffer, so it
    /// is `0`.
    ///
    /// The parallel fast-path block planner (`plan_fast_path_blocks`) uses this
    /// to reproduce the serial byte-cap split byte-for-byte from the per-record
    /// lengths alone — WITHOUT constructing a builder or touching record
    /// bytes — so a parallel gather emits exactly the same block boundaries (and
    /// therefore the same output) as the serial `next_fast_batch` loop.
    const FRAME_OVERHEAD_PER_RECORD: usize;

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

    // Bare bodies into a flat backing buffer — no per-record framing prefix.
    const FRAME_OVERHEAD_PER_RECORD: usize = 0;

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

    // `[u32 LE block_size][body]` framing adds a 4-byte length prefix per record.
    const FRAME_OVERHEAD_PER_RECORD: usize = 4;

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
    /// Name the narrowed-lane variant for diagnostics.
    fn variant_name(&self) -> &'static str {
        match self {
            Self::Empty => "empty",
            Self::K24(_) => "K24",
            Self::Cb32(_) => "Cb32",
            Self::Tert32(_) => "Tert32",
            Self::K40(_) => "K40",
        }
    }

    /// Accumulate one template chunk, which must keep the lane variant fixed.
    ///
    /// The `--key-types` narrowed-lane variant is chosen once per sort and is
    /// global to the run, so every template chunk reaching the merge must carry
    /// the same one. A variant change means the phase-1 producer and the merge
    /// consumer disagree about the key width, and merging on would compare keys
    /// of different layouts and emit silently mis-ordered output.
    ///
    /// # Errors
    ///
    /// Returns `InvalidData` if `chunk`'s variant differs from the accumulated
    /// one. This is the same fail-closed treatment the sibling protocol
    /// violations get (`ensure_single_lane`, `build_driver`) rather than a
    /// panic, so a corrupt stream aborts the sort with a diagnosable error.
    fn push(&mut self, chunk: TemplateMemChunk) -> io::Result<()> {
        /// Build the mismatch error, naming both variants.
        fn mismatch(found: &str, have: &str) -> io::Error {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "SortMerge: template chunk variant changed mid-sort \
                     (accumulated {have}, got {found}); the --key-types lane \
                     variant is global to a sort and must not change"
                ),
            )
        }
        match chunk {
            TemplateMemChunk::K24(c) => match self {
                Self::Empty => *self = Self::K24(vec![c]),
                Self::K24(v) => v.push(c),
                other => return Err(mismatch("K24", other.variant_name())),
            },
            TemplateMemChunk::Cb32(c) => match self {
                Self::Empty => *self = Self::Cb32(vec![c]),
                Self::Cb32(v) => v.push(c),
                other => return Err(mismatch("Cb32", other.variant_name())),
            },
            TemplateMemChunk::Tert32(c) => match self {
                Self::Empty => *self = Self::Tert32(vec![c]),
                Self::Tert32(v) => v.push(c),
                other => return Err(mismatch("Tert32", other.variant_name())),
            },
            TemplateMemChunk::K40(c) => match self {
                Self::Empty => *self = Self::K40(vec![c]),
                Self::K40(v) => v.push(c),
                other => return Err(mismatch("K40", other.variant_name())),
            },
        }
        Ok(())
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
    /// Accumulate one erased chunk into its per-order bucket.
    ///
    /// # Errors
    ///
    /// Propagates the template lane-variant mismatch from
    /// [`TemplateChunks::push`]; the other orders are infallible.
    fn push(&mut self, chunk: MemoryChunkErased) -> io::Result<()> {
        match chunk {
            MemoryChunkErased::Coordinate(v) => self.coordinate.push(v),
            MemoryChunkErased::QuerynameLex(v) => self.queryname_lex.push(v),
            MemoryChunkErased::QuerynameNatural(v) => self.queryname_natural.push(v),
            MemoryChunkErased::TemplateCoordinate(v) => {
                return self.template_coordinate.push(v);
            }
        }
        Ok(())
    }

    fn total_len(&self) -> usize {
        self.coordinate.len()
            + self.queryname_lex.len()
            + self.queryname_natural.len()
            + self.template_coordinate.len()
    }

    /// Fail closed if any lane other than the one `sort_order` selects holds a
    /// chunk. `build_driver` (and the single-chunk fast path) consume only the
    /// selected lane, so a chunk in another lane — a Phase-2 protocol violation
    /// emitting the wrong `MemoryChunkErased` variant — would be silently dropped
    /// even though `total_len()` counted it toward setup completion. Reject it
    /// rather than merge a partial result.
    fn ensure_single_lane(&self, sort_order: SortOrder) -> io::Result<()> {
        let selected_len = match sort_order {
            SortOrder::Coordinate => self.coordinate.len(),
            SortOrder::Queryname(QuerynameComparator::Lexicographic) => self.queryname_lex.len(),
            SortOrder::Queryname(QuerynameComparator::Natural) => self.queryname_natural.len(),
            SortOrder::TemplateCoordinate => self.template_coordinate.len(),
        };
        let stray = self.total_len() - selected_len;
        if stray > 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "SortMerge: {stray} residual memory chunk(s) in a lane not matching the \
                     {sort_order:?} sort order — Phase-2 emitted a mismatched MemoryChunkErased \
                     variant; failing closed rather than silently dropping records",
                ),
            ));
        }
        Ok(())
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
) -> io::Result<Box<dyn MergeDriverDyn + Send>> {
    Ok(match sort_order {
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
            // `Empty` means no residual chunk identified the `--key-types` lane.
            // For valid input this only happens with empty input (no spill files
            // either) — Phase-1's deferred seal always emits a variant-tagged
            // residual otherwise. So `Empty` WITH spill slots can only arise from
            // the documented "defensive/unreachable" no-residual finalize branch
            // (a seal-logic regression). Defaulting to `TemplateKey` (K40) there
            // would decode narrow (K24/Cb32/Tert32) spill files at the wrong key
            // width and silently corrupt output, so fail closed instead of
            // guessing the width. With no slots, any K is safe (nothing to merge).
            TemplateChunks::Empty => {
                if !slots.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "SortMerge: template-coordinate spill slots present but no residual \
                         chunk to identify the key-types lane — refusing to guess the key \
                         width (would mis-decode narrow spill files). This indicates a \
                         Phase-1 seal-logic regression.",
                    ));
                }
                Box::new(MergeDriver::<TemplateKey>::from_slots(
                    slots,
                    MemorySources::Shared(Vec::new()),
                    total_records,
                ))
            }
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
    })
}

enum NextBatch<I> {
    Batch(I),
    Stalled(Option<I>),
    Done(Option<I>, u64),
}

/// Plan the output-block boundaries for a parallel fast-path gather of `chunk`,
/// reproducing the serial [`SortMerge::next_fast_batch`] count/byte-cap split
/// byte-for-byte — but from the per-record lengths alone, without constructing a
/// builder or slicing any record body (a cheap, cache-friendly scan over the
/// chunk's contiguous `len` index).
///
/// Matches the serial loop's framing exactly: each record is "pushed" (the
/// running record count and byte total advance by `1` and
/// `len + O::Builder::FRAME_OVERHEAD_PER_RECORD` respectively), THEN the caps are
/// checked — `count >= target_batch_count` OR `bytes >= output_byte_limit` closes
/// the block after that record. So an oversized single record (whose framed size
/// alone meets the byte cap) closes its own block immediately, identical to the
/// serial path. The returned ranges are contiguous and cover `0..chunk.len()`;
/// an empty chunk yields no blocks (the serial path emits nothing then too).
///
/// `byte_limit == 0` would make every record trip the byte cap (one record per
/// block) — the same degenerate the serial `>=` comparison produces — so no
/// special-casing is needed; the two stay in lockstep.
fn plan_fast_path_blocks<O: MergeOutput>(
    chunk: &MemoryChunkErased,
    target_batch_count: usize,
    output_byte_limit: u64,
) -> Vec<FastBlock> {
    let total = chunk.len();
    let mut blocks = Vec::new();
    let mut start = 0usize;
    let mut count = 0usize;
    let mut bytes: u64 = 0;
    for i in 0..total {
        let framed = u64::from(chunk.record_len(i)) + O::Builder::FRAME_OVERHEAD_PER_RECORD as u64;
        count += 1;
        bytes = bytes.saturating_add(framed);
        let count_full = count >= target_batch_count;
        let bytes_full = bytes >= output_byte_limit;
        if count_full || bytes_full {
            blocks.push(FastBlock { start, end: i + 1 });
            start = i + 1;
            count = 0;
            bytes = 0;
        }
    }
    // Trailing partial block (the serial path's `flush_partial` on `Done`).
    if start < total {
        blocks.push(FastBlock { start, end: total });
    }
    blocks
}

/// Gather one planned [`FastBlock`] into a finished output item, framing each
/// record body `[start, end)` through a fresh `O::Builder` seeded with `ordinal`
/// as its `batch_serial`. Pure/read-only over `chunk` (`record_bytes` is a
/// shared arena slice), so many blocks can be gathered concurrently.
///
/// Byte-identical to the serial `next_fast_batch` for the same record range: the
/// builder frames records the same way and the block boundaries were planned to
/// match (see [`plan_fast_path_blocks`]).
fn gather_fast_block<O: MergeOutput>(
    chunk: &MemoryChunkErased,
    block: FastBlock,
    ordinal: u64,
    initial_bytes: usize,
    target_batch_count: usize,
) -> io::Result<O::Item> {
    let mut builder = O::Builder::with_capacity(ordinal, initial_bytes, target_batch_count);
    for i in block.start..block.end {
        builder.push_record_bytes(chunk.record_bytes(i))?;
    }
    Ok(builder.build())
}

/// Total records covered by a block plan (the sum of each block's range width).
fn blocks_total_records(blocks: &[FastBlock]) -> u64 {
    blocks.iter().map(|b| (b.end - b.start) as u64).sum()
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
    /// Parallel single-source fast path: same precondition as [`SortMergeState::FastPath`], but
    /// the gather is fanned across a bounded rayon pool. The output-block
    /// boundaries are precomputed up front (`blocks`) from the per-record
    /// lengths ALONE — reproducing the serial [`SortMergeState::FastPath`] count/byte-cap split
    /// byte-for-byte — so each block can be gathered independently and the
    /// emitted stream is identical to the serial path (same records, same block
    /// boundaries, same dense ordinals). The step wraps `chunk` in an `Arc` so
    /// the rayon closures can share it (`MemoryChunkErased` is `Send + Sync` and
    /// `record_bytes` is a read-only arena slice).
    FastPathParallel {
        chunk: Arc<MemoryChunkErased>,
        /// Precomputed output-block record ranges `[start, end)`, contiguous and
        /// covering `0..total`; block `i` emits with ordinal `i` (dense).
        blocks: Vec<FastBlock>,
        /// Index of the next block to gather+emit (in ascending, i.e. ordinal,
        /// order — the detached output edge has no reorder stage).
        next_block: usize,
    },
    Done,
}

/// One planned output block for the parallel fast path: the half-open record
/// range `[start, end)` into the sorted chunk. The block's ordinal is its index
/// in the plan vector (dense, `0..blocks.len()`), matching what the serial
/// [`SortMergeState::FastPath`] mints via `next_ordinal`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct FastBlock {
    start: usize,
    end: usize,
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
            // The MemoryChunk `Arc` is created uniquely in the Phase-1 producer
            // (`CompressSpill` / `SpillWrite`) and only ever *moved* (never cloned)
            // through `SortSpillDecompress` to here, so
            // it must be uniquely owned at this single consumer. Fail closed on a
            // shared `Arc` rather than silently deep-cloning a potentially large
            // record vector on the merge setup path.
            let inner = Arc::try_unwrap(chunk).map_err(|_| {
                io::Error::other(
                    "SortMerge: MemoryChunk Arc unexpectedly shared at the merge consumer \
                     (protocol invariant: memory chunks are moved, never cloned)",
                )
            })?;
            memory_chunks.push(inner)?;
            *total_records = (*total_records).max(records_ingested_so_far);
        }
        SortPhase2Event::AllAnnounced {
            slot_count,
            memory_chunk_count,
            total_records: ar_total,
        } => {
            // The Phase-2 protocol emits exactly one `AllAnnounced` (the last event
            // from the Phase-1 producer, `CompressSpill` / `SpillWrite`). A second
            // one is a protocol violation; fail
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
    /// Optional sink for the end-of-run sort summary, filled when the merge
    /// reaches `Done`. The standalone-sort summary finalize hook reads it to
    /// log records processed/written and the spill-chunk count.
    stats_slot: Option<Arc<parking_lot::Mutex<Option<fgumi_sort::SortStats>>>>,
    /// Whether to log the `--sort-stats` merge-loop performance diagnostic
    /// (`Sort merge diag: ...`: stalls/contention/backpressure counters). Off
    /// by default -- it is instrumentation for performance investigations, not
    /// something a normal run should show. See [`Self::with_sort_stats`].
    sort_stats: bool,
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
    /// Worker budget for the in-memory fast-path parallel gather (see
    /// [`Self::fast_path_pool`]). `1` (the default) keeps the historical
    /// single-threaded serial gather; `> 1` fans the gather across a bounded,
    /// step-owned rayon pool. Wired from the sort chain's phase-2 (`merge`)
    /// thread budget in `add_sort`.
    ///
    /// This affects ONLY the fast path (single already-sorted in-memory chunk).
    /// The k-way merge path is untouched — it is a genuinely serial loser-tree
    /// walk whose winner order cannot be produced out of order.
    fast_path_threads: usize,
    /// Minimum record count for the fast path to fan across the pool
    /// ([`FAST_PATH_PARALLEL_MIN_RECORDS`] in production). A field, not the bare
    /// const, only so tests can force the parallel path on a small chunk to
    /// assert byte-for-byte parity with the serial gather without materialising
    /// 64 Ki records.
    fast_path_min_records: usize,
    /// Bounded rayon pool for the fast-path parallel gather, built once on first
    /// use and sized to [`fast_path_threads`](Self::fast_path_threads).
    ///
    /// Owned by the step (not the pipeline work-stealing pool, which is not
    /// reachable from a step body) and bounded to the sort's thread budget so
    /// the gather never oversubscribes past `--threads` — the same discipline
    /// the chunk-sort front uses (`CoordinateChunkSorter::rayon_pool`). Only
    /// built when `fast_path_threads > 1`, and only for a fast-path run.
    fast_path_pool: Option<rayon::ThreadPool>,
    /// Blocks gathered by `emit_fast_batches_parallel` but not yet pushed
    /// downstream (the queue rejected a push mid-window). Drained IN ORDER — with
    /// the `held` slot holding at most the single most-recently-rejected item —
    /// at the top of the next `try_run` before any new gather, preserving the
    /// dense ascending-ordinal emission the detached output edge requires.
    fast_pending: std::collections::VecDeque<O::Item>,
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
            stats_slot: None,
            sort_stats: false,
            processed: 0,
            chunk_count: 0,
            dbg: MergeDiag::default(),
            fast_path_threads: 1,
            fast_path_min_records: FAST_PATH_PARALLEL_MIN_RECORDS,
            fast_path_pool: None,
            fast_pending: std::collections::VecDeque::new(),
        }
    }

    /// Set the worker budget for the in-memory fast-path parallel gather.
    ///
    /// `1` (the default) preserves the single-threaded serial gather. `> 1` fans
    /// the gather of the single already-sorted in-memory chunk across a bounded
    /// step-owned rayon pool (built lazily, capped at `threads`), which removes
    /// the flat ~serial gather cost that otherwise floors the in-memory sort's
    /// wall clock while the pipeline pool sits idle. Clamped to `>= 1`.
    ///
    /// Only the fast path is affected; the k-way merge path is untouched.
    #[must_use]
    pub fn with_fast_path_threads(mut self, threads: usize) -> Self {
        self.fast_path_threads = threads.max(1);
        self
    }

    /// Test-only: lower the record-count threshold that gates the parallel
    /// fast-path gather, so a small synthetic chunk exercises the parallel path
    /// (parity vs. the serial gather) without materialising
    /// [`FAST_PATH_PARALLEL_MIN_RECORDS`] records.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_fast_path_min_records(mut self, min: usize) -> Self {
        self.fast_path_min_records = min;
        self
    }

    /// Attach a slot to receive the end-of-run [`fgumi_sort::SortStats`] when
    /// the merge completes (records processed/written + spill-chunk count). Used
    /// by the standalone-sort summary finalize hook; runall leaves it unset.
    #[must_use]
    pub fn with_stats_slot(
        mut self,
        slot: Arc<parking_lot::Mutex<Option<fgumi_sort::SortStats>>>,
    ) -> Self {
        self.stats_slot = Some(slot);
        self
    }

    /// Enable or disable the `--sort-stats` merge-loop performance diagnostic
    /// (`Sort merge diag: ...`). Off by default.
    #[must_use]
    pub fn with_sort_stats(mut self, enabled: bool) -> Self {
        self.sort_stats = enabled;
        self
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        // `true` once the slot is clear (was empty, or the held event flushed);
        // `false` while it's still held under backpressure. Uses the canonical
        // re-hold helper so the put-back-on-reject invariant lives in one place.
        !matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld)
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

    /// `records_out` accumulates one per [`MergeStep::Produced`] — the caller
    /// (`emit_batches_cooperative`, in turn `try_run`) sums this across every
    /// `next_batch` call in a dispatch and bumps the `records` counter ONCE per
    /// `try_run` with the batch total, rather than per record here.
    ///
    /// # Panics
    ///
    /// Panics if `self.state` is not `Merging`.
    fn next_batch(&mut self, records_out: &mut u64) -> io::Result<NextBatch<O::Item>> {
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
                    *records_out += 1;
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

    /// `records_out` accumulates the batch total of records merged this call;
    /// see [`Self::next_batch`].
    ///
    /// # Panics
    ///
    /// Panics if `self.state` is not `Merging`.
    fn emit_batches_cooperative(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        records_out: &mut u64,
    ) -> io::Result<StepOutcome> {
        let mut delivered = 0usize;
        loop {
            match self.next_batch(records_out)? {
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
                    //
                    // Gated on `--sort-stats` (`self.sort_stats`): this is
                    // performance-investigation instrumentation, not something a
                    // normal run should print.
                    if self.sort_stats {
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
                    }
                    if let Some(slot) = &self.stats_slot {
                        *slot.lock() = Some(fgumi_sort::SortStats {
                            total_records: self.processed,
                            output_records: merged,
                            runs_written: self.chunk_count,
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
    /// `records_out` accumulates one per record gathered from the chunk; see
    /// [`Self::next_batch`] for how the caller uses this.
    ///
    /// # Panics
    ///
    /// Panics if `self.state` is not `FastPath`.
    fn next_fast_batch(&mut self, records_out: &mut u64) -> io::Result<NextBatch<O::Item>> {
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
            *records_out += 1;
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
    /// `records_out` accumulates the batch total of records gathered this call;
    /// see [`Self::next_batch`].
    ///
    /// # Panics
    ///
    /// Panics if `self.state` is not `FastPath`.
    fn emit_fast_batches(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        records_out: &mut u64,
    ) -> io::Result<StepOutcome> {
        let mut delivered = 0usize;
        loop {
            match self.next_fast_batch(records_out)? {
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
                    // Shared with the parallel `finish_fast_path` so the message,
                    // diag, and `SortStats` shape cannot drift between the two.
                    self.finalize_fast_path_done(count);
                    return Ok(if delivered > 0 {
                        StepOutcome::Progress
                    } else {
                        StepOutcome::NoProgress
                    });
                }
            }
        }
    }

    /// Build the step-owned bounded rayon pool for the parallel fast-path
    /// gather, if not already built. Sized to `fast_path_threads` and capped
    /// so the gather never oversubscribes past `--threads`. Idempotent.
    ///
    /// On a `ThreadPoolBuilder` failure this leaves `fast_path_pool` `None`; the
    /// gather then runs on rayon's global pool via a plain `par_iter` (still
    /// correct, just not thread-budget-bounded), which never happens in practice
    /// (the builder only fails on an OS thread-spawn error).
    fn ensure_fast_path_pool(&mut self) {
        if self.fast_path_pool.is_some() {
            return;
        }
        match rayon::ThreadPoolBuilder::new()
            .num_threads(self.fast_path_threads)
            // Name the workers so they are attributable in the repo's profiling
            // workflow (tricorder / samply), matching the other bounded sort
            // pools (`build_sort_rayon_pool`, the arena seal pools).
            .thread_name(|i| format!("fast-path-gather-{i}"))
            .build()
        {
            Ok(pool) => self.fast_path_pool = Some(pool),
            Err(e) => log::warn!(
                "SortMerge: failed to build fast-path rayon pool ({e}); \
                 falling back to the global pool for the parallel gather"
            ),
        }
    }

    /// Cooperative emit loop for the PARALLEL single-source fast path.
    ///
    /// Gathers the next window of up to
    /// `fast_path_threads * FAST_PATH_BLOCKS_PER_WORKER_WINDOW` planned blocks
    /// concurrently (each block frames its own record range into a finished
    /// output item — read-only over the shared `chunk`), then pushes them in
    /// ascending block/ordinal order through the same held/backpressure idiom as
    /// the serial path. Blocks MUST leave in order: `SortMerge` is `Detached`, so
    /// its output edge has no reorder stage (the by-ordinal reassembly is
    /// downstream at `BgzfCompress`), and the reorder there is a dense
    /// single-cursor stream.
    ///
    /// Emitting a bounded window per burst (rather than the whole plan at once)
    /// keeps peak extra memory to at most `fast_path_threads *
    /// FAST_PATH_BLOCKS_PER_WORKER_WINDOW` uncompressed blocks in flight (a
    /// function of config, not input size — it scales with the thread budget, so
    /// it is tens of blocks on a high-core box, not "a handful") and preserves
    /// output backpressure. The gathered-but-not-yet-pushed blocks are stashed in
    /// `self.held` one at a time via the existing `Unpushed` slot; a full
    /// downstream queue returns `Progress` and the remaining gathered blocks are
    /// re-pushed on the next dispatch before any new gather.
    ///
    /// `records_out` accumulates one per record emitted this call.
    fn emit_fast_batches_parallel(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        records_out: &mut u64,
    ) -> io::Result<StepOutcome> {
        // 1. Re-push any blocks gathered on a prior dispatch but rejected under
        //    backpressure. `flush_held` (called at the top of `try_run`) has
        //    already cleared the single `held` item; drain the rest in order.
        if !self.drain_fast_pending(ctx) {
            return Ok(StepOutcome::Progress); // still backpressured
        }

        let total_blocks = match &self.state {
            SortMergeState::FastPathParallel { blocks, .. } => blocks.len(),
            _ => unreachable!("emit_fast_batches_parallel called outside FastPathParallel state"),
        };
        let next_block = match &self.state {
            SortMergeState::FastPathParallel { next_block, .. } => *next_block,
            _ => unreachable!(),
        };

        // 2. Nothing left to gather → finalize (this also covers an empty chunk,
        //    whose plan has zero blocks).
        if next_block >= total_blocks {
            let count = match &self.state {
                SortMergeState::FastPathParallel { blocks, .. } => blocks_total_records(blocks),
                _ => unreachable!(),
            };
            return Ok(self.finish_fast_path(count));
        }

        // 3. Gather the next window of blocks in parallel.
        let bytes_cap = usize::try_from(self.output_byte_limit).unwrap_or(usize::MAX);
        let initial_bytes = INITIAL_OUTPUT_BUFFER_BYTES.min(bytes_cap);
        let target = self.target_batch_count;
        let window =
            self.fast_path_threads.saturating_mul(FAST_PATH_BLOCKS_PER_WORKER_WINDOW).max(1);
        let window_end = (next_block + window).min(total_blocks);
        let base_ordinal = next_block as u64;

        let items: Vec<O::Item> = {
            use rayon::prelude::*;
            let (chunk, window_blocks) = match &self.state {
                SortMergeState::FastPathParallel { chunk, blocks, .. } => {
                    (Arc::clone(chunk), blocks[next_block..window_end].to_vec())
                }
                _ => unreachable!(),
            };
            let run = || -> io::Result<Vec<O::Item>> {
                window_blocks
                    .par_iter()
                    .enumerate()
                    .map(|(offset, &block)| {
                        gather_fast_block::<O>(
                            &chunk,
                            block,
                            base_ordinal + offset as u64,
                            initial_bytes,
                            target,
                        )
                    })
                    .collect()
            };
            match &self.fast_path_pool {
                Some(pool) => pool.install(run)?,
                None => run()?,
            }
        };

        let window_records: u64 = match &self.state {
            SortMergeState::FastPathParallel { blocks, .. } => {
                blocks_total_records(&blocks[next_block..window_end])
            }
            _ => unreachable!(),
        };
        *records_out += window_records;

        // Advance past the gathered window, then push the gathered items in order.
        if let SortMergeState::FastPathParallel { next_block, .. } = &mut self.state {
            *next_block = window_end;
        }
        self.fast_pending.extend(items);
        if !self.drain_fast_pending(ctx) {
            return Ok(StepOutcome::Progress); // downstream full mid-window
        }

        // 4. Window pushed cleanly; finalize if that was the last one.
        if window_end >= total_blocks {
            let count = match &self.state {
                SortMergeState::FastPathParallel { blocks, .. } => blocks_total_records(blocks),
                _ => unreachable!(),
            };
            return Ok(self.finish_fast_path(count));
        }
        Ok(StepOutcome::Progress)
    }

    /// Push every buffered `fast_pending` block in order, honouring the held-slot
    /// backpressure idiom. Returns `true` once the buffer is fully drained,
    /// `false` if a push was rejected (the offending item is put in `self.held`
    /// and the remainder stays in `fast_pending` for the next dispatch).
    fn drain_fast_pending(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        while let Some(item) = self.fast_pending.pop_front() {
            if let Err(unpushed) = ctx.outputs.push(item) {
                self.dbg.output_full += 1;
                self.held.put(unpushed);
                return false;
            }
        }
        true
    }

    /// Emit the fast-path completion log + `--sort-stats` diagnostic, populate
    /// the stats slot, and transition to `Done`. Shared by BOTH fast-path
    /// finalizers — the serial `emit_fast_batches` `Done` arm and the parallel
    /// [`finish_fast_path`](Self::finish_fast_path) — so the completion message,
    /// the diag prefix, and the `SortStats` shape cannot drift between them. The
    /// `StepOutcome` return is deliberately NOT part of this helper: the serial
    /// arm returns `Progress`/`NoProgress` off its own `delivered` counter, so it
    /// stays at each call site.
    fn finalize_fast_path_done(&mut self, count: u64) {
        log::info!("Sort in-memory fast path complete: {count} records (single source, no merge)");
        // `--sort-stats` (`self.sort_stats`): the merge-loop diagnostic in
        // `emit_batches_cooperative` never runs on this path (there is no k-way
        // merge to diagnose -- the single already-sorted chunk is gathered
        // directly), so say so explicitly rather than staying silent, which is
        // indistinguishable from the flag being ignored. Deliberately a different
        // prefix from "Sort merge diag:" (the k-way-merge counters emitted by
        // `emit_batches_cooperative`) so the two are never mistaken for one
        // another in a log or a test.
        if self.sort_stats {
            log::info!("Sort fast-path diag: in-memory fast path taken, no k-way merge occurred");
        }
        if let Some(slot) = &self.stats_slot {
            *slot.lock() = Some(fgumi_sort::SortStats {
                total_records: self.processed,
                output_records: count,
                runs_written: 0,
            });
        }
        self.state = SortMergeState::Done;
    }

    /// Finalize the fast path (parallel variant): run the shared
    /// [`finalize_fast_path_done`](Self::finalize_fast_path_done) work and return
    /// `Progress`. Unlike the serial arm (which returns `NoProgress` when it
    /// delivered nothing this call), the parallel path only reaches here after
    /// its window drained cleanly, so a plain `Progress` is correct.
    fn finish_fast_path(&mut self, count: u64) -> StepOutcome {
        self.finalize_fast_path_done(count);
        StepOutcome::Progress
    }

    fn transition_to_merging(&mut self) -> io::Result<()> {
        if !matches!(&self.state, SortMergeState::WaitingForSetup { .. }) {
            return Ok(());
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
        // Fail closed before consuming only the selected lane (fast path or
        // build_driver): a chunk stranded in a non-selected lane would otherwise
        // be dropped silently.
        memory_chunks.ensure_single_lane(self.sort_order)?;
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
            // Parallel gather when a worker budget was requested AND there is
            // enough work to amortise the fan-out (a tiny chunk is faster serial:
            // the pool build + block-plan + join overhead would dominate). Below
            // the threshold, or with a single thread, fall through to the serial
            // gather — byte-identical, just cheaper for small inputs.
            if self.fast_path_threads > 1 && total >= self.fast_path_min_records {
                let blocks = plan_fast_path_blocks::<O>(
                    &chunk,
                    self.target_batch_count,
                    self.output_byte_limit,
                );
                // A plan that collapses to a single block (or none) has no
                // parallelism to exploit, so the pool build + rayon fan-out would
                // be pure overhead. Fall through to the serial gather, which
                // emits byte-identically. (`plan_fast_path_blocks` only borrows
                // `chunk`, so it is still owned here.) Unreachable under the
                // production defaults — `target_batch_count` (1024) caps every
                // block, so `total >= FAST_PATH_PARALLEL_MIN_RECORDS` (64Ki)
                // always yields >= 64 blocks — but guards a degenerate plan if
                // those knobs ever become independently tunable.
                if blocks.len() > 1 {
                    self.ensure_fast_path_pool();
                    self.state = SortMergeState::FastPathParallel {
                        chunk: Arc::new(chunk),
                        blocks,
                        next_block: 0,
                    };
                    return Ok(());
                }
            }
            let builder =
                O::Builder::with_capacity(0, initial_bytes_for_init, self.target_batch_count);
            self.state =
                SortMergeState::FastPath { chunk, cursor: 0, total, builder, next_ordinal: 0 };
            return Ok(());
        }
        let driver = build_driver(self.sort_order, slots, memory_chunks, total_records)?;
        let bytes_cap = usize::try_from(self.output_byte_limit).unwrap_or(usize::MAX);
        // Seed the first buffer modestly; subsequent buffers are sized from the
        // prior batch's actual byte length (see `next_batch`).
        let initial_bytes = INITIAL_OUTPUT_BUFFER_BYTES.min(bytes_cap);
        let builder = O::Builder::with_capacity(0, initial_bytes, self.target_batch_count);
        self.state = SortMergeState::Merging { driver, builder, next_ordinal: 0 };
        Ok(())
    }
}

impl<O: MergeOutput> Step for SortMerge<O> {
    type Input = SortPhase2Event;
    type Outputs = OrderedBytesSingle<O::Item>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SortMerge",
            // The merge runs off the work-stealing pool, on the sort's shared
            // COORDINATION driver thread (N+2) — the same thread that ran the
            // phase-1 admit/sort/frame steps, which have Finished and left the
            // driver's live set by the time the phase-2 merge runs, so the merge
            // effectively gets a dedicated thread in phase 2 (mirrors main's main
            // thread). Its cooperative `try_run` body is UNCHANGED —
            // `run_detached_driver` drives it with the same `run_worker_loop` the
            // pool uses (Park backoff), parking on `Contention`/`NoProgress`
            // (winner-slot momentarily empty / output full) instead of the pool
            // re-dispatching it. `Detached` collapses the declared `ByItemOrdinal`
            // output to `None` exactly as `Serial` would (see
            // `effective_branch_orderings`), so the output transport — a direct
            // byte-bounded queue, no reorder stage — is byte-for-byte identical;
            // the LoserTree core, source order (`slots.sort_by_key(file_id)` +
            // residual last), and tie-break are untouched. SortMerge is only ever
            // built by the sort chain's `add_sort`, so this is sort-chain-only.
            kind: StepKind::Detached,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn detached_group(&self) -> DetachedGroup {
        // Co-located with the phase-1 coordination steps on ONE driver thread —
        // phase 1 and phase 2 are temporally disjoint, so this is the true N+2.
        DetachedGroup::Shared(crate::sort::SORT_COORD_GROUP)
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] = &[CounterSpec::new("records", "records")];
        SPECS
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
            self.transition_to_merging()?;
        }

        // Batch total for this call: one bump per `try_run`, accumulated across
        // however many `next_batch`/`next_fast_batch` calls the emit loop below
        // makes (never a per-record bump).
        let mut records_this_call: u64 = 0;
        let outcome = match &self.state {
            SortMergeState::Merging { .. } => {
                self.emit_batches_cooperative(ctx, &mut records_this_call)
            }
            SortMergeState::FastPath { .. } => self.emit_fast_batches(ctx, &mut records_this_call),
            SortMergeState::FastPathParallel { .. } => {
                self.emit_fast_batches_parallel(ctx, &mut records_this_call)
            }
            SortMergeState::Done => Ok(StepOutcome::Finished),
            SortMergeState::WaitingForSetup { .. } => {
                unreachable!("Phase 1 must have left state non-WaitingForSetup")
            }
        };
        ctx.counters.add(RECORDS, records_this_call);
        outcome
    }
}

#[cfg(test)]
mod tests;
