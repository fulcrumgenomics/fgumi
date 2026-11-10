//! `SortMerge` — third step of the runall-sort three-step chain.

use std::collections::HashMap;
use std::io;
use std::sync::{Arc, OnceLock};

#[cfg(test)]
use fgumi_sort::RangeSizing;
use fgumi_sort::{
    CbKey32, InMemoryChunk, MergeDriver, MergeDriverDyn, MergeStep, QnClassTable, QnKey, QnLex,
    QnMergeKey, QnNatural, QuerynameComparator, RECORD_LEN_BYTES, RawCoordinateKey, SortMergeSlot,
    SortOrder, SpillKeyKind, TemplateKey, TemplateKey24, TemplateMemChunk, TertKey32,
};

use crate::output_compress::{OutputBlock, OutputCompressor, OutputDepth};
use crate::sort::merge_job::{GatherSizing, MergeJob, MergeResult, MergeTicket};
use crate::sort::merge_mode::{MergeMode, SerialReason};
use crate::sort::protocol::{MemoryChunkErased, SortPhase2Event};
use crate::types::{DecompressedBlock, RecordBatch, RecordBatchBuilder};
use fgumi_pipeline_core::{
    HeapSize, Ordered, OrderedInbox, SideChannel, Unpushed,
    held::HeldSlot,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{CounterSpec, DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

pub(crate) mod claim;

pub use claim::ClaimPolicy;
mod range_mode;
use range_mode::RangeModeState;

/// Default output batch size: 1024 records per emitted `RecordBatch`.
pub const DEFAULT_TARGET_BATCH_COUNT: usize = 1024;

/// Max output batches emitted per `try_run` invocation in `Merging`.
const MAX_DRAIN_BATCHES_PER_LOCK: usize = 8;

/// Minimum record count for the fast path to fan its gather across the pool.
/// Below this the serial gather is faster: the block plan and the per-block job
/// round trip would dominate a small chunk (which the serial path
/// clears in well under a millisecond). ~1 output block per worker at the
/// default 1024-record batch, so the fan-out always has real work to spread.
const FAST_PATH_PARALLEL_MIN_RECORDS: usize = 64 * 1024;

/// Gather jobs in flight per worker. Gather mode keeps at most
/// `fanout_width * this` jobs posted but not yet taken back, bounding peak
/// extra memory to that many uncompressed blocks (vs. materialising the whole
/// output) while keeping every worker fed.
const FAST_PATH_BLOCKS_PER_WORKER_WINDOW: usize = 4;

/// Every record ingested must come out of the merge: a slot that finalized a
/// clean EOF early (a truncated or misparsed spill run) would otherwise end the
/// merge short with a clean exit.
fn ensure_all_records_merged(ingested: u64, merged: u64) -> io::Result<()> {
    if merged == ingested {
        return Ok(());
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "sort lost records: read {ingested} but merged {merged} (differ by {})",
            ingested.abs_diff(merged)
        ),
    ))
}

/// `fgumi merge`'s counterpart of [`ensure_all_records_merged`]: its inputs'
/// record counts are not known up front (the merge itself reads them), so
/// instead every input must have been read to its end, with no error, when the
/// merge reports done. The driver only retires a source at its clean EOF, so
/// this fails only if a source was retired early: the merge would otherwise
/// end short with a clean exit. `inputs` names each slot by `file_id`.
fn ensure_merge_inputs_read(slots: &[Arc<SortMergeSlot>], inputs: &[String]) -> io::Result<()> {
    for slot in slots {
        let name = usize::try_from(slot.file_id)
            .ok()
            .and_then(|i| inputs.get(i))
            .map_or_else(|| format!("slot {}", slot.file_id), Clone::clone);
        if slot.has_error() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("merge failed reading input {name}"),
            ));
        }
        if !slot.is_drained() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("merge lost records: it ended before input {name} was read to its end"),
            ));
        }
    }
    Ok(())
}

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

    /// Compress `item` in place with `c` (the range merge's thread does this
    /// to blocks its output edge refuses); `false` when the item cannot be or
    /// already is compressed. Only [`CompressibleBlockOutput`] compresses;
    /// every other output is left as it is.
    ///
    /// # Errors
    ///
    /// The compressor's own error.
    fn compress_in_place(_item: &mut Self::Item, _c: &mut OutputCompressor) -> io::Result<bool> {
        Ok(false)
    }
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

    /// Append every record of `framed`, whole `[u32 LE len][body]` records
    /// back to back, as [`Self::push_record_bytes`] would push each body;
    /// returns the record count. The default pushes each body; a builder whose
    /// own framing is that layout copies the span at once.
    ///
    /// # Errors
    ///
    /// `framed` ends inside a record, or a body cannot be framed.
    fn push_framed_records(&mut self, framed: &[u8]) -> io::Result<u64> {
        walk_framed_records(framed, |body| self.push_record_bytes(body))
    }
}

/// Call `f` on each body of `framed` (whole `[u32 LE len][body]` records back
/// to back) and return the record count.
///
/// # Errors
///
/// `framed` ends inside a record or its length prefix, or `f` fails.
fn walk_framed_records(
    framed: &[u8],
    mut f: impl FnMut(&[u8]) -> io::Result<()>,
) -> io::Result<u64> {
    let mut span = framed;
    let mut records = 0u64;
    while let Some((prefix, rest)) = span.split_first_chunk::<RECORD_LEN_BYTES>() {
        let len = u32::from_le_bytes(*prefix) as usize;
        let (body, rest) = rest.split_at_checked(len).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "batch segment ends inside a record")
        })?;
        f(body)?;
        records += 1;
        span = rest;
    }
    if !span.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "batch segment ends inside a record prefix",
        ));
    }
    Ok(records)
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

/// Terminal standalone-sort output whose blocks the merge may compress
/// itself: [`BlockOutput`]'s framing, emitted as [`OutputBlock::Raw`], which
/// the range merge's thread turns into [`OutputBlock::Compressed`] (with the
/// output's settings, see [`SortMerge::with_output_compression`]) while its
/// output edge refuses them. Wired to `BgzfCompress<OutputBlock>`, which
/// compresses the raw blocks and passes the compressed ones through.
pub struct CompressibleBlockOutput;

impl MergeOutput for CompressibleBlockOutput {
    type Item = OutputBlock;
    type Builder = CompressibleBlockBuilder;

    fn compress_in_place(item: &mut OutputBlock, c: &mut OutputCompressor) -> io::Result<bool> {
        item.compress_in_place(c)
    }
}

/// [`BlockBuilder`] emitting an [`OutputBlock::Raw`].
pub struct CompressibleBlockBuilder(BlockBuilder);

impl MergeBatchBuilder for CompressibleBlockBuilder {
    type Item = OutputBlock;

    const FRAME_OVERHEAD_PER_RECORD: usize = BlockBuilder::FRAME_OVERHEAD_PER_RECORD;

    fn with_capacity(batch_serial: u64, bytes_cap: usize, records_cap: usize) -> Self {
        Self(BlockBuilder::with_capacity(batch_serial, bytes_cap, records_cap))
    }

    fn push_record_bytes(&mut self, body: &[u8]) -> io::Result<()> {
        self.0.push_record_bytes(body)
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn total_bytes(&self) -> usize {
        self.0.total_bytes()
    }

    fn build(self) -> OutputBlock {
        OutputBlock::Raw(self.0.build())
    }

    fn push_framed_records(&mut self, framed: &[u8]) -> io::Result<u64> {
        self.0.push_framed_records(framed)
    }
}

/// The settings of the output's BGZF compression, for a merge that may
/// compress its own blocks ([`SortMerge::with_output_compression`]): the
/// same two values the sink's `BgzfCompress` is built with.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OutputCompression {
    /// The compression level.
    pub level: u32,
    /// Whether the sink indexes the BAM (each block carries its manifest).
    pub index_bam: bool,
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
    const FRAME_OVERHEAD_PER_RECORD: usize = RECORD_LEN_BYTES;

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

    /// The span is already this builder's framing: check it holds whole
    /// records, then copy it at once.
    fn push_framed_records(&mut self, framed: &[u8]) -> io::Result<u64> {
        let records = walk_framed_records(framed, |_| Ok(()))?;
        self.bytes.extend_from_slice(framed);
        self.records += usize::try_from(records).expect("records of an in-memory span fit usize");
        Ok(records)
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

/// A template-coordinate merge's `--key-types` lane, as its spill slots' key
/// kind names it: the one dispatch from a [`SpillKeyKind`] to a template key
/// type, shared by the serial merge's collector and the range plan.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TemplateLane {
    /// [`TemplateKey24`].
    K24,
    /// [`CbKey32`].
    Cb32,
    /// [`TertKey32`].
    Tert32,
    /// [`TemplateKey`].
    K40,
}

impl TemplateLane {
    /// The lane `kind` names.
    ///
    /// # Errors
    ///
    /// `InvalidData` if `kind` is not a template-coordinate kind.
    fn of(kind: SpillKeyKind) -> io::Result<Self> {
        match kind {
            SpillKeyKind::TemplateK24 => Ok(Self::K24),
            SpillKeyKind::TemplateCb32 => Ok(Self::Cb32),
            SpillKeyKind::TemplateTert32 => Ok(Self::Tert32),
            SpillKeyKind::TemplateK40 => Ok(Self::K40),
            other => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("SortMerge: a template-coordinate merge got spill slots keyed {other:?}"),
            )),
        }
    }
}

impl TemplateChunks {
    /// The empty collector of the lane `kind` names: a merge with no residual
    /// chunk still knows its lane from its spill slots.
    ///
    /// # Errors
    ///
    /// `InvalidData` if `kind` is not a template-coordinate kind.
    fn empty_for(kind: SpillKeyKind) -> io::Result<Self> {
        Ok(match TemplateLane::of(kind)? {
            TemplateLane::K24 => Self::K24(Vec::new()),
            TemplateLane::Cb32 => Self::Cb32(Vec::new()),
            TemplateLane::Tert32 => Self::Tert32(Vec::new()),
            TemplateLane::K40 => Self::K40(Vec::new()),
        })
    }

    /// Resolve an `Empty` collector over `slots` to the lane their key kind
    /// names; any other collector (or no slots) is returned unchanged.
    ///
    /// # Errors
    ///
    /// `InvalidData` if the slots disagree on their kind, name a non-template
    /// kind, or any slot carries none (`template_lane_unknown`).
    fn resolve_with_slots(self, slots: &[Arc<SortMergeSlot>]) -> io::Result<Self> {
        if !matches!(self, Self::Empty) || slots.is_empty() {
            return Ok(self);
        }
        match slots_key_kind(slots)? {
            Some(kind) => Self::empty_for(kind),
            None => Err(template_lane_unknown()),
        }
    }

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
    queryname_lex: Vec<InMemoryChunk<QnKey>>,
    queryname_natural: Vec<InMemoryChunk<QnKey>>,
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

/// How `SortMerge` merges the inputs of `fgumi merge` rather than a sort's
/// spill runs ([`SortMerge::with_merge_inputs`]): every slot is a body-keyed BAM
/// input, merged on the full key with each input's order checked.
#[derive(Clone)]
pub struct MergeInputs {
    /// Each input's name, by `file_id` (the order check's error names it).
    pub inputs: Vec<String>,
    /// The order's `fgumi sort --order` name, for the order check's error.
    pub order_name: &'static str,
    /// The full template key's context (template-coordinate only).
    pub template_ctx: Option<Arc<fgumi_sort::TemplateKeyContext>>,
    /// A rewrite of every record before its key is extracted (the read-tag
    /// renames of `fgumi merge`), if any input needs one.
    pub rewrite: Option<fgumi_sort::RecordRewrite>,
}

/// A driver over `fgumi merge`'s inputs (see [`MergeInputs`]): the kind
/// dispatch ([`fgumi_sort::merge_driver_for_kind`]) over the order's merge key
/// (template-coordinate on the full key, with its context), then the input
/// order check and the read-tag rewrite.
fn build_merge_inputs_driver(
    sort_order: SortOrder,
    slots: Vec<Arc<SortMergeSlot>>,
    demand: Option<&Arc<fgumi_sort::MergeDemand>>,
    class_table: Option<&Arc<QnClassTable>>,
    mode: &MergeInputs,
) -> io::Result<Box<dyn MergeDriverDyn + Send>> {
    use fgumi_sort::{MergeCtxErased, SpillKeyKind};
    let table =
        || MergeCtxErased::Queryname(class_table.map_or_else(QnClassTable::new, Arc::clone));
    let (kind, ctx) = match sort_order {
        SortOrder::Coordinate => (SpillKeyKind::Coordinate, MergeCtxErased::None),
        SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
            (SpillKeyKind::QuerynameLex, table())
        }
        SortOrder::Queryname(QuerynameComparator::Natural) => {
            (SpillKeyKind::QuerynameNatural, table())
        }
        SortOrder::TemplateCoordinate => {
            let ctx = mode.template_ctx.clone().ok_or_else(|| {
                io::Error::other(
                    "SortMerge: a template-coordinate merge of BAM inputs needs a key context",
                )
            })?;
            (SpillKeyKind::TemplateK40, MergeCtxErased::Template(ctx))
        }
    };
    let demand = demand.map_or_else(|| Arc::new(fgumi_sort::MergeDemand::new()), Arc::clone);
    let mut driver = fgumi_sort::merge_driver_for_kind(kind, slots, ctx, demand);
    driver.set_input_order_check(mode.inputs.clone(), mode.order_name);
    if let Some(rewrite) = &mode.rewrite {
        driver.set_record_rewrite(rewrite);
    }
    Ok(driver)
}

/// The key kind every slot of `slots` carries: `None` when there are no slots
/// or any slot has no kind attached.
///
/// # Errors
///
/// `InvalidData` if two slots carry different kinds.
fn slots_key_kind(slots: &[Arc<SortMergeSlot>]) -> io::Result<Option<SpillKeyKind>> {
    let mut found: Option<SpillKeyKind> = None;
    for slot in slots {
        let Some(kind) = slot.key_kind() else { return Ok(None) };
        match found {
            Some(prev) if prev != kind => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "SortMerge: spill slots carry different key kinds ({prev:?}, {kind:?})"
                    ),
                ));
            }
            _ => found = Some(kind),
        }
    }
    Ok(found)
}

fn build_driver(
    sort_order: SortOrder,
    slots: Vec<Arc<SortMergeSlot>>,
    chunks: MemoryChunksByKind,
    total_records: u64,
    demand: Option<&Arc<fgumi_sort::MergeDemand>>,
    class_table: Option<&Arc<QnClassTable>>,
) -> io::Result<Box<dyn MergeDriverDyn + Send>> {
    let table = || class_table.map_or_else(QnClassTable::new, Arc::clone);
    Ok(match sort_order {
        SortOrder::Coordinate => Box::new(attach(
            MergeDriver::<RawCoordinateKey>::from_slots(
                slots,
                chunks.coordinate,
                total_records,
                (),
            ),
            demand,
        )),
        SortOrder::Queryname(QuerynameComparator::Lexicographic) => Box::new(attach(
            MergeDriver::<QnMergeKey<QnLex>, QnKey>::from_slots(
                slots,
                chunks.queryname_lex,
                total_records,
                table(),
            ),
            demand,
        )),
        SortOrder::Queryname(QuerynameComparator::Natural) => Box::new(attach(
            MergeDriver::<QnMergeKey<QnNatural>, QnKey>::from_slots(
                slots,
                chunks.queryname_natural,
                total_records,
                table(),
            ),
            demand,
        )),
        SortOrder::TemplateCoordinate => {
            match chunks.template_coordinate.resolve_with_slots(&slots)? {
                // `Empty` here means no slots (a slot-bearing `Empty` was resolved
                // from the slots' key kind, or failed closed, above): empty input,
                // so any K is safe (nothing to merge).
                TemplateChunks::Empty => Box::new(attach(
                    MergeDriver::<TemplateKey>::from_slots(slots, Vec::new(), total_records, None),
                    demand,
                )),
                TemplateChunks::K24(v) => Box::new(attach(
                    MergeDriver::<TemplateKey24>::from_slots(slots, v, total_records, ()),
                    demand,
                )),
                TemplateChunks::Cb32(v) => Box::new(attach(
                    MergeDriver::<CbKey32>::from_slots(slots, v, total_records, ()),
                    demand,
                )),
                TemplateChunks::Tert32(v) => Box::new(attach(
                    MergeDriver::<TertKey32>::from_slots(slots, v, total_records, ()),
                    demand,
                )),
                TemplateChunks::K40(v) => Box::new(attach(
                    MergeDriver::<TemplateKey>::from_slots(slots, v, total_records, None),
                    demand,
                )),
            }
        }
    })
}

/// The fail-closed error for template-coordinate spill slots with no residual
/// chunk to name their `--key-types` lane and no key kind attached to a slot
/// (see `build_driver`).
fn template_lane_unknown() -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        "SortMerge: template-coordinate spill slots present but no residual \
         chunk to identify the key-types lane — refusing to guess the key \
         width (would mis-decode narrow spill files). This indicates a \
         Phase-1 seal-logic regression.",
    )
}

/// Hand the sort's merge demand to `driver`, which publishes its predictions
/// and frontier through it.
fn attach<K: fgumi_sort::RawSortKey + Default + Send + 'static, R: fgumi_sort::ChunkRef<K>>(
    driver: MergeDriver<K, R>,
    demand: Option<&Arc<fgumi_sort::MergeDemand>>,
) -> MergeDriver<K, R> {
    match demand {
        Some(d) => driver.with_demand(Arc::clone(d)),
        None => driver,
    }
}

enum NextBatch<I> {
    Batch(I),
    /// The driver stalled on a slot. Nothing is emitted: the partial batch
    /// stays in the builder until it fills or the merge is `Done`.
    Stalled,
    Done(Option<I>, u64),
}

/// Ask the pool for one parked worker, mirroring the outcome for the merge
/// demand's counters (fgumi-sort does not depend on pipeline-core).
fn merge_pool_request(pool: &fgumi_pipeline_core::PoolHandle<'_>) -> fgumi_sort::MergePoolRequest {
    use fgumi_pipeline_core::PoolRequest;
    use fgumi_sort::MergePoolRequest as M;
    match pool.request_worker() {
        PoolRequest::Woken => M::Woken,
        PoolRequest::AllAwake => M::AllAwake,
        PoolRequest::Pending => M::Pending,
        PoolRequest::Unavailable => M::Unavailable,
    }
}

/// The serial merge's batch-close rule, checked after each record is pushed:
/// a batch closes once it holds `target_batch_count` records or
/// `output_byte_limit` bytes (record bytes plus the builder's per-record frame
/// overhead). The serial merge, the fast-path block plan and the range
/// merge's batch cutter all close batches through this one predicate, which is
/// what keeps their batch boundaries identical.
pub(crate) fn batch_full(
    count: usize,
    bytes: u64,
    target_batch_count: usize,
    output_byte_limit: u64,
) -> bool {
    count >= target_batch_count || bytes >= output_byte_limit
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
pub(crate) fn plan_fast_path_blocks<O: MergeOutput>(
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
        if batch_full(count, bytes, target_batch_count, output_byte_limit) {
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
pub(crate) fn gather_fast_block<O: MergeOutput>(
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

enum SortMergeState<B: MergeBatchBuilder> {
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
    /// the chunk is already globally sorted and no k-way merge is needed. (A
    /// front splits its final run into pieces only after a spill, so a sort that
    /// fits in memory always arrives here.) Gather
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
    /// Pooled single-source fast path: same precondition as
    /// [`SortMergeState::FastPath`], but each output block is gathered by a
    /// `MergeWork` pool worker. The block boundaries are precomputed (`blocks`)
    /// from the per-record lengths ALONE — reproducing the serial
    /// [`SortMergeState::FastPath`] count/byte-cap split byte-for-byte — so each
    /// block is gathered independently and the emitted stream is identical to
    /// the serial path (same records, same block boundaries, same dense
    /// ordinals). The merge posts one [`MergeJob::Gather`] per block on its
    /// second branch and takes the results back from its inbox in ordinal order
    /// (the detached output edge has no reorder stage). `chunk` is shared with
    /// the jobs through an `Arc` (`MemoryChunkErased` is `Send + Sync` and
    /// `record_bytes` is a read-only arena slice).
    GatherMode {
        chunk: Arc<MemoryChunkErased>,
        /// Precomputed output-block record ranges `[start, end)`, contiguous and
        /// covering `0..total`; block `i` emits with ordinal `i` (dense).
        blocks: Vec<FastBlock>,
        /// Index of the next block to post as a job.
        next_to_post: usize,
        /// Ordinal of the next result to take from the inbox and emit.
        next_to_take: u64,
    },
    /// Range-parallel final merge of a spilled coordinate or template sort
    /// (the read planner published [`MergeMode::Range`]): `MergeWork` merges
    /// each key range from the spill files (range mode spills the final run,
    /// so nothing is merged from memory), this step cuts
    /// the merged ranges into the serial merge's batches, and `MergeWork`
    /// assembles each batch. Boxed: it is much larger than the other states.
    RangeMode(Box<RangeModeState<B>>),
    Done,
}

/// One planned output block for the pooled fast path: the half-open record
/// range `[start, end)` into the sorted chunk. The block's ordinal is its index
/// in the plan vector (dense, `0..blocks.len()`), matching what the serial
/// `SortMergeState::FastPath` mints via `next_ordinal`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FastBlock {
    pub(crate) start: usize,
    pub(crate) end: usize,
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
            if let std::collections::hash_map::Entry::Vacant(e) = slot_index.entry(slot.slot_id()) {
                e.insert(slots.len());
                slots.push(slot);
            }
            *total_records = (*total_records).max(records_ingested_so_far);
        }
        SortPhase2Event::MemoryChunk { chunk, records_ingested_so_far } => {
            // The MemoryChunk `Arc` is created uniquely in the Phase-1 producer
            // (`SpillWrite`) and only ever *moved* (never cloned)
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
            // from the Phase-1 producer, `SpillWrite`). A second
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
    /// Spill-run accounting from `SpillWrite`, so the summary's `runs_written`
    /// counts the runs the spill phase wrote rather than the (possibly fewer,
    /// after consolidation) runs this merge read.
    spill_stats: Option<Arc<crate::sort::SpillRunStats>>,
    /// Whether to log the `--sort-stats` merge-loop performance diagnostic
    /// (`Sort merge diag: ...`: stalls/contention/backpressure counters, then
    /// the merge demand's `MergeDemandSnapshot::log_lines` when a demand is
    /// attached, then the `MergeDriver`'s `Merge output: records borrowed B /
    /// assembled A` split, printed whenever the k-way merge ran and never on
    /// the single-chunk fast path). Off by default -- it is instrumentation for
    /// performance investigations, not something a normal run should show. See
    /// [`Self::with_sort_stats`].
    sort_stats: bool,
    /// The merge supply's end-of-merge lines (`--sort-stats` only).
    supply_diagnostics: Option<super::supply_ledger::SupplyDiagnostics>,
    /// Total records ingested, captured at the merge transition (the summary's
    /// "records processed").
    processed: u64,
    /// Spill runs read by the final merge (after any consolidation), captured at
    /// the merge transition. Zero for a fully in-memory sort.
    chunk_count: usize,
    /// INSTRUMENTATION (lever-2 merge-stall diagnosis; `RUST_LOG=info` at Done).
    /// `SortMerge` runs on a single dedicated `Detached` thread (one instance,
    /// never `new_worker_copy`'d), so plain `&mut self` counters are sound — no
    /// atomics needed.
    dbg: MergeDiag,
    /// Pool workers the in-memory fast path fans its gather across (see
    /// [`Self::with_fanout_width`]); sets the gather window. `1` (the default)
    /// keeps the serial gather; `> 1` (with an inbox) posts gather jobs to
    /// `MergeWork`. Wired from the sort chain's phase-2 (`merge`) thread budget
    /// in `add_sort`.
    ///
    /// Range mode's window of live key ranges comes from the planner's
    /// decision ([`MergeMode::Range`]), not from this. The serial k-way merge
    /// path is untouched — it is a genuinely serial loser-tree walk whose
    /// winner order cannot be produced out of order.
    fanout_width: usize,
    /// Minimum record count for the fast path to fan across the pool
    /// ([`FAST_PATH_PARALLEL_MIN_RECORDS`] in production). A field, not the bare
    /// const, only so tests can force the parallel path on a small chunk to
    /// assert byte-for-byte parity with the serial gather without materialising
    /// 64 Ki records.
    fast_path_min_records: usize,
    /// The inbox `MergeWork` delivers gather results to (see
    /// [`Self::with_inbox`]). `None` keeps the fast path serial whatever the
    /// fan-out width: nothing would run the posted jobs.
    inbox: Option<Arc<OrderedInbox<MergeResult<O>>>>,
    /// The inbox's owner waker was bound (on the first dispatch, which is the
    /// first time this step sees its pool handle).
    inbox_bound: bool,
    /// A ticket the job branch refused, re-posted before any new one.
    held_job: HeldSlot<Unpushed<MergeTicket<O>>>,
    /// Range mode: the job this thread claimed and runs itself, between
    /// dispatches (see `emit_range_mode`).
    own_range: Option<MergeJob<O>>,
    /// Range mode: the jobs posted and not yet started, which this thread and
    /// the `MergeWork` workers claim from (see `claim::ClaimTable`).
    claim_table: Option<Arc<claim::ClaimTable<O>>>,
    /// Who may claim range mode's posted jobs (see [`Self::with_claim_policy`]).
    claim_policy: claim::ClaimPolicy,
    /// The compressor range mode's thread compresses refused output blocks
    /// with (see [`Self::with_output_compression`]); `None` unless set, and
    /// unused by every output but [`CompressibleBlockOutput`].
    output_compressor: Option<OutputCompressor>,
    /// Compression's queue depth on this merge's output edge, shared with
    /// the sink's `BgzfCompress` (see [`Self::with_output_depth`]).
    output_depth: Option<Arc<OutputDepth>>,
    /// Test-only override of [`Self::job_edge_limit`], so a test can make the
    /// job branch refuse a post within the window.
    #[cfg(test)]
    job_edge_limit_override: Option<u64>,
    /// The sort's merge-wide demand, shared with `SortSpillDecompress`: on a
    /// stall the merge declares the slot it waits on, and the delivery to that
    /// slot unparks this step's driver thread. `None` when the merge runs
    /// without a spill supply that notifies (unit tests that drive it alone).
    demand: Option<Arc<fgumi_sort::MergeDemand>>,
    /// This dispatch's stall registered on the merge demand (its `await_slot`
    /// found nothing). `try_run` settles it once the dispatch's outcome is
    /// known: a park (`Contention`) is recorded as one, and any other outcome
    /// withdraws the registration, because the merge runs again before it
    /// parks and re-registers then.
    awaiting: bool,
    /// The sort's queryname class table: the queryname merge keys compare
    /// through the windows phase 1 built against it. `None` merges against an
    /// empty table (class 0: correct, without the window shortcut).
    class_table: Option<Arc<QnClassTable>>,
    /// Merge `fgumi merge`'s BAM inputs, not spill runs (see [`MergeInputs`]).
    merge_inputs: Option<MergeInputs>,
    /// The merge inputs' slots, kept beside the driver for the completion
    /// check (`ensure_merge_inputs_read`); empty for a sort.
    merge_input_slots: Vec<Arc<SortMergeSlot>>,
    /// The current stall episode (one episode spans every `Stalled` until the
    /// next `Produced` or `Done`).
    episode: StallEpisode,
    /// The merge supply's ledger: self-serve claims are booked to the
    /// consumer. Set with the demand by [`Self::with_spill_supply`], so the
    /// merge cannot book to a ledger the planner does not read; `None` in unit
    /// tests that drive the merge alone.
    ledger: Option<Arc<super::supply_ledger::SupplyLedger>>,
    /// Decompressor for self-serve, built on the first self-served block.
    decomp: Option<fgumi_sort::SpillBlockDecompressor>,
    /// The per-slot reorder-window byte budget self-serve claims under (the
    /// decompress step's).
    window_budget: u64,
    /// The final merge's mode, published by the read planner at
    /// `AllAnnounced` (before this merge leaves setup). Set with the demand by
    /// [`Self::with_spill_supply`]; `None`, or an unset cell, merges serially.
    merge_mode: Option<Arc<OnceLock<MergeMode>>>,
    /// Test-only override of the range plan's sizing, so a small fixture cuts
    /// many ranges.
    #[cfg(test)]
    range_sizing_override: Option<RangeSizing>,
}

/// The merge's own supply path: on a stall it claims, decompresses and
/// publishes the awaited slot's stashed blocks itself, through the same
/// [`fgumi_sort::ClaimedBlock::decompress_and_publish`] a `SortSpillDecompress`
/// worker uses (same claim admission, same publish), and never reads the spill
/// file.
pub(crate) struct SelfServe<'a> {
    pub(crate) decomp: &'a mut Option<fgumi_sort::SpillBlockDecompressor>,
    pub(crate) ledger: Option<&'a super::supply_ledger::SupplyLedger>,
    pub(crate) window_budget: u64,
}

impl SelfServe<'_> {
    /// Serve `slot` until it can make progress or nothing more is claimable or
    /// drainable. Returns the blocks served and whether the slot now has a
    /// block or EOF.
    ///
    /// The drain-only call before giving up moves a block that waits in
    /// `reorder` behind a full FIFO (the consumer has just emptied it) without
    /// any worker.
    ///
    /// # Errors
    /// A decompress failure: the slot is marked failed first, and the error
    /// names the slot as the merge's own failed-slot error does.
    pub(crate) fn run(&mut self, slot: &fgumi_sort::SortMergeSlot) -> io::Result<(u64, bool)> {
        let mut served = 0u64;
        loop {
            if let Some(block) = slot.bp_claim_raw(self.window_budget) {
                if let Some(l) = self.ledger {
                    l.sub_stash(block.charge, true);
                }
                let decomp =
                    self.decomp.get_or_insert_with(fgumi_sort::SpillBlockDecompressor::new);
                block.decompress_and_publish(decomp).map_err(|e| {
                    io::Error::new(
                        e.kind(),
                        format!(
                            "SortMerge: spill decompression error on slot {}: {e}",
                            slot.file_id
                        ),
                    )
                })?;
                served += 1;
            } else if !slot.bp_drain_and_finalize() {
                return Ok((served, slot.has_block_or_eof()));
            }
            if slot.has_block_or_eof() {
                return Ok((served, true));
            }
        }
    }
}

/// Bookkeeping for one stall episode of the merge: when it began, and whether
/// its single pool-worker request has been made.
#[derive(Default)]
struct StallEpisode {
    /// When the episode began; `None` while the merge is not stalled.
    started: Option<std::time::Instant>,
    /// The episode's pool-worker request was made (at most one per episode).
    requested: bool,
    /// That request found a sleeping worker: one parked on the event-count
    /// when it asked, or one it woke from a timer park
    /// (`MergeDemandStats::record_pool_request`).
    had_sleeper: bool,
}

impl StallEpisode {
    /// Open an episode on the first stall that will park, classifying the
    /// awaited slot. A re-stall within the episode only refreshes `awaited`.
    fn begin(&mut self, demand: &fgumi_sort::MergeDemand, slot: &SortMergeSlot) {
        if self.started.is_none() {
            self.started = Some(std::time::Instant::now());
            demand.stats().record_stall(slot.awaited_state());
        }
        demand.set_starved(true);
    }

    /// Close the episode, if one is open: book its duration (and, when its
    /// request found a sleeping worker, the sleeper time), lower `starved`,
    /// and re-arm the once-per-episode request.
    fn end(&mut self, demand: Option<&fgumi_sort::MergeDemand>) {
        let Some(t0) = self.started.take() else { return };
        if let Some(d) = demand {
            let ns = fgumi_pipeline_core::runtime::elapsed_ns(t0);
            d.stats().record_stall_ns(ns);
            if self.requested && self.had_sleeper {
                d.stats().record_stall_with_sleeper_ns(ns);
            }
            d.set_starved(false);
        }
        self.requested = false;
        self.had_sleeper = false;
    }
}

/// Lever-2 diagnostic counters: is the serial merge starved on decompress
/// (`input-empty`/`stalls`) or blocked on the downstream writer
/// (`output_full`), and how often does a dispatch deliver nothing (`contention`)?
#[derive(Default, Clone, Copy)]
struct MergeDiag {
    /// Merge-loop passes that ended `Stalled` — the winning source's next block
    /// was not yet decompressed (INPUT-STARVED: the lever-2 hypothesis) — and
    /// that `MergeDemand::await_slot`'s re-check did not resolve. A stall whose
    /// re-check finds the block or EOF already landed keeps merging and is not
    /// counted; before the merge demand every stall was, so this count is not
    /// comparable with a pre-`MergeDemand` run's `stalls=`.
    stalls: u64,
    /// `ctx.outputs.push` returned `Err` — downstream (compress/write) full
    /// (OUTPUT-BACKPRESSURE).
    output_full: u64,
    /// `try_run` returned `Contention` — the merge had nothing to deliver this
    /// dispatch, so its driver parked until the awaited slot's delivery woke it
    /// (or its idle timer expired).
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
            spill_stats: None,
            sort_stats: false,
            supply_diagnostics: None,
            processed: 0,
            chunk_count: 0,
            dbg: MergeDiag::default(),
            fanout_width: 1,
            fast_path_min_records: FAST_PATH_PARALLEL_MIN_RECORDS,
            inbox: None,
            inbox_bound: false,
            held_job: HeldSlot::new(),
            own_range: None,
            claim_table: None,
            claim_policy: claim::ClaimPolicy::Auto,
            output_compressor: None,
            output_depth: None,
            #[cfg(test)]
            job_edge_limit_override: None,
            demand: None,
            awaiting: false,
            class_table: None,
            merge_inputs: None,
            merge_input_slots: Vec::new(),
            episode: StallEpisode::default(),
            ledger: None,
            decomp: None,
            window_budget: super::spill_decompress::reorder_window_budget(output_byte_limit),
            merge_mode: None,
            #[cfg(test)]
            range_sizing_override: None,
        }
    }

    /// Merge `fgumi merge`'s BAM inputs instead of a sort's spill runs (see
    /// [`MergeInputs`]).
    #[must_use]
    pub fn with_merge_inputs(mut self, mode: MergeInputs) -> Self {
        self.merge_inputs = Some(mode);
        self
    }

    /// Join the merge supply: share its merge demand (a stalled merge parks
    /// until the block it awaits lands) and book self-serve claims in its
    /// ledger — one value, so the merge cannot pair one supply's demand with
    /// another's ledger.
    #[must_use]
    pub fn with_spill_supply(mut self, supply: &super::supply_ledger::SpillSupply) -> Self {
        self.demand = Some(Arc::clone(&supply.demand));
        self.ledger = Some(Arc::clone(&supply.ledger));
        self.merge_mode = Some(Arc::clone(&supply.merge_mode));
        self
    }

    /// Test-only: size the range plan with `sizing` instead of from the
    /// published range budget.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_range_sizing_for_test(mut self, sizing: RangeSizing) -> Self {
        self.range_sizing_override = Some(sizing);
        self
    }

    /// Run range mode under claim `policy` (default [`ClaimPolicy::Auto`]):
    /// the chain builder passes the hidden diagnostic
    /// `FGUMI_SORT_RANGE_CLAIM`'s value, resolved with
    /// [`ClaimPolicy::from_env_value`].
    #[must_use]
    pub fn with_claim_policy(mut self, policy: ClaimPolicy) -> Self {
        self.claim_policy = policy;
        self
    }

    /// Let range mode's thread compress the output blocks its output edge
    /// refuses, with `settings` (the sink compressor's level and indexing),
    /// when it has nothing else to do. Only a [`CompressibleBlockOutput`]
    /// merge compresses; for any other output this has no effect.
    #[must_use]
    pub fn with_output_compression(mut self, settings: OutputCompression) -> Self {
        self.output_compressor = Some(OutputCompressor::new(settings.level, settings.index_bam));
        self
    }

    /// Count this merge's output in `depth`, which the sink's compressor
    /// pops from: range mode then gates the pool's merge-side claims on
    /// compression's queue depth (workers claim only while fewer blocks than
    /// the phase-2 width wait). Without it the gate reads the merge's own
    /// refused-block backlog.
    #[must_use]
    pub fn with_output_depth(mut self, depth: Arc<OutputDepth>) -> Self {
        self.output_depth = Some(depth);
        self
    }

    /// The class table this merge compares queryname keys through (test
    /// support: wiring checks compare it by pointer with the sort's).
    #[cfg(any(test, feature = "test-utils"))]
    #[must_use]
    pub fn class_table_for_test(&self) -> Option<&Arc<QnClassTable>> {
        self.class_table.as_ref()
    }

    /// Merge queryname orders against the sort's class table, so spilled and
    /// residual keys compare through the same class windows phase 1 used.
    #[must_use]
    pub fn with_class_table(mut self, table: Arc<QnClassTable>) -> Self {
        self.class_table = Some(table);
        self
    }

    /// Share the sort's [`fgumi_sort::MergeDemand`] with the spill supply, so a
    /// stalled merge parks until the block it awaits lands instead of until its
    /// driver's idle timer.
    #[must_use]
    pub fn with_merge_demand(mut self, demand: Arc<fgumi_sort::MergeDemand>) -> Self {
        self.demand = Some(demand);
        self
    }

    /// Fan the in-memory fast-path gather across `width` pool workers, and
    /// return this merge together with the `MergeWork` step that must be wired
    /// on its branch 1 (`(merge_tail.0, BranchIdx(1))`). The helper creates the
    /// inbox and hands both ends their half, so a merge built through it and the
    /// `MergeWork` it returns share one inbox. That pairing holds only while the
    /// caller wires this `MergeWork` (not another, such as the idle one
    /// [`serial_merge_output`](crate::sort::merge_job::serial_merge_output)
    /// adds) on branch 1. `cap` is the phase-2 cap the `MergeWork` step admits
    /// under (`None` leaves it uncapped). See `Self::with_fanout_width`.
    #[must_use]
    pub fn with_pooled_gather(
        self,
        width: usize,
        cap: Option<Arc<fgumi_pipeline_core::PhaseCap>>,
    ) -> (Self, fgumi_pipeline_core::FanoutWorker<MergeTicket<O>>) {
        let inbox = OrderedInbox::new(true);
        let work = crate::sort::merge_job::merge_work(Arc::clone(&inbox))
            .with_step_cap(cap.map(fgumi_pipeline_core::StepCap::Fixed));
        (self.with_inbox(inbox).with_fanout_width(width), work)
    }

    /// Take pooled gather results from `inbox`. Crate-private: the caller must
    /// wire a `MergeWork` step (a [`fgumi_pipeline_core::FanoutWorker`] over
    /// [`MergeTicket`] delivering to this same inbox) on this step's branch 1, and
    /// a mismatch strands every gather result in an orphan inbox. Production
    /// wiring goes through [`Self::with_pooled_gather`], which returns the
    /// merge and the matching `MergeWork` together; the pairing is guaranteed
    /// only through it. With an inbox and a fan-out width above 1, a large
    /// single-chunk fast path posts its gather jobs there and waits for their
    /// results.
    #[must_use]
    pub(crate) fn with_inbox(mut self, inbox: Arc<OrderedInbox<MergeResult<O>>>) -> Self {
        self.inbox = Some(inbox);
        self
    }

    /// Set how many pool workers the in-memory fast-path gather fans across.
    ///
    /// `1` (the default) preserves the single-threaded serial gather. `> 1`
    /// (with [`Self::with_inbox`]) gathers the single already-sorted in-memory
    /// chunk's blocks on the pipeline pool, at most `width *
    /// FAST_PATH_BLOCKS_PER_WORKER_WINDOW` jobs ahead of the blocks emitted,
    /// which removes the flat ~serial gather cost that otherwise floors the
    /// in-memory sort's wall clock. Clamped to `>= 1`.
    ///
    /// Range mode takes its window of live key ranges from the planner's
    /// decision ([`MergeMode::Range`]), not from this; the serial k-way merge
    /// path is untouched. Crate-private: a width does nothing
    /// without the inbox, so callers outside the crate set both through
    /// [`Self::with_pooled_gather`].
    #[must_use]
    pub(crate) fn with_fanout_width(mut self, width: usize) -> Self {
        self.fanout_width = width.max(1);
        self
    }

    /// Gather jobs kept in flight (posted, not yet taken) in gather mode.
    fn gather_window(&self) -> usize {
        self.fanout_width.saturating_mul(FAST_PATH_BLOCKS_PER_WORKER_WINDOW).max(1)
    }

    /// The byte bound of the job branch: room for twice the gather window of
    /// tickets, so a gather post within the window is never refused. Range
    /// mode does not wait on it: tickets it cannot push stay owed while it
    /// keeps posting and claiming jobs itself.
    fn job_edge_limit(&self) -> u64 {
        #[cfg(test)]
        if let Some(limit) = self.job_edge_limit_override {
            return limit;
        }
        let jobs = self.gather_window().saturating_mul(2);
        (jobs.saturating_mul(std::mem::size_of::<MergeTicket<O>>())) as u64
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

    /// Test-only: shrink the job branch's byte bound to `limit`, so a post
    /// within the gather window is refused and the held-job path runs.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_job_edge_limit(mut self, limit: u64) -> Self {
        self.job_edge_limit_override = Some(limit);
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

    /// Report `SortStats::runs_written` from `SpillWrite`'s accounting (runs the
    /// spill phase wrote) instead of this merge's source count, which is smaller
    /// whenever `--max-temp-files` consolidated runs.
    #[must_use]
    pub fn with_spill_stats(mut self, spill_stats: Arc<crate::sort::SpillRunStats>) -> Self {
        self.spill_stats = Some(spill_stats);
        self
    }

    /// Enable or disable the `--sort-stats` merge-loop performance diagnostic
    /// (`Sort merge diag: ...`, the merge-demand lines, and the k-way merge's
    /// `Merge output:` borrowed/assembled line). Off by default.
    #[must_use]
    pub fn with_sort_stats(mut self, enabled: bool) -> Self {
        self.sort_stats = enabled;
        self
    }

    /// Log the merge supply's lines at the end of the merge (built only under
    /// `--sort-stats`).
    #[must_use]
    pub fn with_supply_diagnostics(
        mut self,
        d: Option<super::supply_ledger::SupplyDiagnostics>,
    ) -> Self {
        self.supply_diagnostics = d;
        self
    }

    /// Settle this dispatch's merge-demand registration (`awaiting`) once its
    /// outcome is known. A `Contention` parks the driver until the awaited
    /// slot's delivery (or its idle timer), so it is a parking registration.
    /// Any other outcome runs the merge again before it parks — after a stall
    /// in a dispatch that delivered a batch — so the registration is
    /// withdrawn: a delivery must not unpark a merge that is not waiting on
    /// it, and the next stall re-registers (and re-checks) before any park.
    fn settle_await(&mut self, outcome: &io::Result<StepOutcome>) {
        if !std::mem::take(&mut self.awaiting) {
            return;
        }
        if let Some(d) = &self.demand {
            if matches!(outcome, Ok(StepOutcome::Contention)) {
                d.stats().record_park();
            } else {
                d.clear_awaited();
            }
        }
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        // `true` once the slot is clear (was empty, or the held event flushed);
        // `false` while it's still held under backpressure (put back, never
        // dropped). Branch 0 is the merged output.
        match self.held.take() {
            None => true,
            Some(unpushed) => match ctx.outputs.view().a.retry(unpushed) {
                Ok(()) => true,
                Err(again) => {
                    self.held.put(again);
                    false
                }
            },
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
        let demand = self.demand.as_deref();
        let awaiting = &mut self.awaiting;
        let episode = &mut self.episode;
        let mut serve = SelfServe {
            decomp: &mut self.decomp,
            ledger: self.ledger.as_deref(),
            window_budget: self.window_budget,
        };
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
                    episode.end(demand);
                    if batch_full(builder.len(), builder.total_bytes() as u64, target, byte_limit) {
                        return Ok(NextBatch::Batch(flush(builder, next_ordinal)));
                    }
                }
                MergeStep::Stalled => {
                    if let (Some(d), Some(slot)) = (demand, driver.stalled_slot()) {
                        // Declare the slot awaited, then re-check it (the wake
                        // protocol in `fgumi_sort::MergeDemand`): a block or EOF
                        // that already landed means keep merging; otherwise the
                        // delivery to this slot unparks the driver.
                        if d.await_slot(slot) {
                            continue;
                        }
                        // Before parking, decompress the slot's already-read
                        // blocks on this thread; the slot stays awaited, so a
                        // worker's delivery meanwhile still wakes the driver.
                        let (served, ready) = serve.run(slot)?;
                        d.stats().record_self_serve(served, ready);
                        if ready {
                            d.clear_awaited();
                            continue;
                        }
                        *awaiting = true;
                        episode.begin(d, slot);
                    }
                    // The partial batch stays in the builder: a batch closes
                    // only when it is full (count or byte cap) or at `Done`.
                    // The state persists in `Merging`, so the next call keeps
                    // filling it.
                    return Ok(NextBatch::Stalled);
                }
                MergeStep::Done => {
                    // The last episode (typically a run-end EOF stall) ends here.
                    episode.end(demand);
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
                    if let Err(unpushed) = ctx.outputs.view().a.push(batch) {
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
                NextBatch::Stalled => {
                    // INPUT-STARVED: the driver couldn't advance because the
                    // winning source's next block isn't decompressed yet. The
                    // partial batch stays in the builder.
                    self.dbg.stalls += 1;
                    return Ok(if delivered > 0 {
                        self.dbg.progress_dispatches += 1;
                        StepOutcome::Progress
                    } else {
                        // Nothing delivered this dispatch: the driver parks
                        // until the awaited slot's delivery unparks it (the
                        // merge demand), with its idle timer only as the bound.
                        // The awaited block comes from pool work that reaches
                        // the merge off any queue edge, so nothing else wakes a
                        // parked pool worker for it: ask for one, once per stall
                        // episode (the pool keeps at most one request pending).
                        if let Some(d) = &self.demand
                            && self.episode.started.is_some()
                            && !self.episode.requested
                        {
                            self.episode.requested = true;
                            let parked = ctx.pool.parked_workers();
                            let outcome = merge_pool_request(&ctx.pool);
                            self.episode.had_sleeper =
                                d.stats().record_pool_request(outcome, parked);
                        }
                        self.dbg.contention += 1;
                        StepOutcome::Contention
                    });
                }
                NextBatch::Done(partial, merged) => {
                    // Before the last batch leaves: a merge of BAM inputs has
                    // no ingested count to check against (its records are read
                    // by the merge itself), so it checks every input was read
                    // to its end instead.
                    match &self.merge_inputs {
                        Some(mode) => {
                            ensure_merge_inputs_read(&self.merge_input_slots, &mode.inputs)?;
                        }
                        None => ensure_all_records_merged(self.processed, merged)?,
                    }
                    if let Some(batch) = partial
                        && let Err(unpushed) = ctx.outputs.view().a.push(batch)
                    {
                        self.dbg.output_full += 1;
                        self.held.put(unpushed);
                        return Ok(StepOutcome::Progress);
                    }
                    self.merge_input_slots.clear();
                    log::info!("Sort merge complete: {merged} records merged");
                    // INSTRUMENTATION (lever-2): is the serial merge starved on
                    // decompress (stalls/contention high) or blocked on the
                    // writer (output_full high)? `stalls` counts merge-loop
                    // passes that ended input-starved with no block landed by
                    // the await re-check (see `MergeDiag::stalls`); `contention` counts
                    // dispatches that produced nothing (the driver then parks);
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
                        self.log_demand_lines();
                        // This function runs only in `Merging` (the fast path has
                        // its own `Done`), and the state leaves it only below.
                        let SortMergeState::Merging { driver, .. } = &self.state else {
                            unreachable!("emit_batches_cooperative runs only while Merging")
                        };
                        let (borrowed, assembled) = driver.record_counts();
                        log::info!(
                            "Merge output: records borrowed {borrowed} / assembled {assembled}"
                        );
                        self.log_supply_lines();
                    }
                    self.publish_merge_stats(merged);
                    self.state = SortMergeState::Done;
                    // Reaching `Done` is progress even when this dispatch
                    // delivered nothing (the run-end stall: the last source's
                    // EOF woke the merge, which then had nothing left to emit).
                    // `NoProgress` would park the driver on its idle timer
                    // before the next dispatch reports `Finished`.
                    return Ok(StepOutcome::Progress);
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
                    if let Err(unpushed) = ctx.outputs.view().a.push(batch) {
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
                NextBatch::Stalled => unreachable!("FastPath never stalls (all in memory)"),
                NextBatch::Done(partial, count) => {
                    if let Some(batch) = partial {
                        if let Err(unpushed) = ctx.outputs.view().a.push(batch) {
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

    /// Cooperative emit loop for the POOLED single-source fast path
    /// ([`SortMergeState::GatherMode`]). One dispatch:
    ///
    /// 1. re-posts a job the job branch refused (the output branch's held block
    ///    was already flushed at the top of `try_run`);
    /// 2. takes every result that is ready, in ordinal order, and pushes it on
    ///    branch 0 (a refused push is held and ends the dispatch);
    /// 3. finishes once the last block was taken;
    /// 4. posts ONE gather job while fewer than the window are in flight (one
    ///    per dispatch: each post wakes one more pool worker, and a `Progress`
    ///    dispatch runs again at once);
    /// 5. with nothing taken, asks the inbox afresh whether the next result is
    ///    there (`await_next`), which records the wait; `false` with nothing
    ///    posted parks the driver (`NoProgress`) until that delivery wakes it.
    ///
    /// Step 5 runs on every dispatch that took nothing, never trusting an
    /// earlier dispatch's answer: a delivery clears the owner's waiting flag
    /// even when its wake was skipped because this thread was running, so only
    /// a fresh `await_next` re-records the wait before a park.
    ///
    /// `records_out` accumulates one per record emitted this call.
    fn emit_gather_mode(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        records_out: &mut u64,
    ) -> StepOutcome {
        let inbox = Arc::clone(self.inbox.as_ref().expect("gather mode runs only with an inbox"));
        let window = self.gather_window() as u64;
        let mut progressed = false;

        // 1. A job the job branch refused goes first, so jobs stay in order.
        if let Some(unpushed) = self.held_job.take() {
            match ctx.outputs.view().b.retry(unpushed) {
                Ok(()) => progressed = true,
                Err(again) => self.held_job.put(again),
            }
        }

        // 2. Take ready results in ordinal order.
        let mut took = false;
        while let Some((ordinal, result)) = inbox.try_take_next() {
            took = true;
            let SortMergeState::GatherMode { blocks, next_to_take, .. } = &mut self.state else {
                unreachable!("emit_gather_mode called outside GatherMode")
            };
            debug_assert_eq!(ordinal, *next_to_take, "the inbox yields ordinals in order");
            let block = blocks[usize::try_from(ordinal).expect("ordinal indexes the plan")];
            *next_to_take += 1;
            *records_out += (block.end - block.start) as u64;
            let MergeResult::Block(item) = result else {
                unreachable!("GatherMode posts only gather jobs, which yield blocks")
            };
            if let Err(unpushed) = ctx.outputs.view().a.push(item) {
                self.dbg.output_full += 1;
                self.held.put(unpushed);
                return StepOutcome::Progress;
            }
        }

        let SortMergeState::GatherMode { chunk, blocks, next_to_post, next_to_take } =
            &mut self.state
        else {
            unreachable!("emit_gather_mode called outside GatherMode")
        };

        // 3. Every block taken and pushed: done.
        if *next_to_take == blocks.len() as u64 {
            let count = blocks_total_records(blocks);
            return self.finish_fast_path(count);
        }

        // 4. Post one job while below the window.
        if !self.held_job.is_held()
            && *next_to_post < blocks.len()
            && (*next_to_post as u64) < *next_to_take + window
        {
            let ordinal = *next_to_post as u64;
            let bytes_cap = usize::try_from(self.output_byte_limit).unwrap_or(usize::MAX);
            let sizing = GatherSizing {
                initial_bytes: INITIAL_OUTPUT_BUFFER_BYTES.min(bytes_cap),
                target_batch_count: self.target_batch_count,
            };
            // Charge the block's builder reservation until the result lands
            // (delivery replaces it by the block's real heap size).
            inbox.reserve(ordinal, sizing.initial_bytes as u64);
            let job = MergeJob::Gather {
                chunk: Arc::clone(chunk),
                block: blocks[*next_to_post],
                ordinal,
                sizing,
                _output: std::marker::PhantomData,
            };
            *next_to_post += 1;
            if let Err(unpushed) = ctx.outputs.view().b.push(MergeTicket::direct(job)) {
                self.held_job.put(unpushed);
            }
            progressed = true;
        }

        if took || progressed {
            return StepOutcome::Progress;
        }
        // 5. Nothing taken or posted: re-record the wait on this dispatch.
        if inbox.await_next() { StepOutcome::Progress } else { StepOutcome::NoProgress }
    }

    /// Log the merge demand's `--sort-stats` lines, if a demand is wired. Part
    /// of the completion tail the serial and range merges share.
    fn log_demand_lines(&self) {
        if let Some(d) = &self.demand {
            for line in d.snapshot().log_lines() {
                log::info!("{line}");
            }
        }
    }

    /// Log the spill supply's `--sort-stats` lines, if they are collected.
    fn log_supply_lines(&self) {
        if let Some(s) = &self.supply_diagnostics {
            for line in s.lines() {
                log::info!("{line}");
            }
        }
    }

    /// Fill the stats slot for a k-way (serial or range) merge that emitted
    /// `merged` records.
    fn publish_merge_stats(&self, merged: u64) {
        if let Some(slot) = &self.stats_slot {
            *slot.lock() = Some(fgumi_sort::SortStats {
                total_records: self.processed,
                output_records: merged,
                runs_written: self.spill_stats.as_ref().map_or(self.chunk_count, |s| {
                    usize::try_from(s.runs_written()).unwrap_or(usize::MAX)
                }),
            });
        }
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
    /// delivered nothing this call), the pooled path only reaches here in a
    /// dispatch that took and pushed the last block, so a plain `Progress` is correct.
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
        let mode = self.merge_mode.as_ref().and_then(|m| m.get()).cloned();
        if let Some(MergeMode::Range { budget, inputs, fit, index_bytes }) = mode {
            self.transition_to_range_mode(slots, &memory_chunks, total_records, &fit)?;
            self.log_range_mode(&budget, &inputs, &fit, index_bytes);
            return Ok(());
        }
        if self.sort_stats {
            // Unset only when the planner failed (the run is ending) or no
            // planner joined this merge (`fgumi merge`, unit tests).
            let reason = match mode {
                Some(MergeMode::Serial(reason)) => reason,
                _ => SerialReason::NoPolicy,
            };
            log::info!("Sort merge mode: serial ({reason})");
        }
        let bytes_cap_for_init = usize::try_from(self.output_byte_limit).unwrap_or(usize::MAX);
        let initial_bytes_for_init = INITIAL_OUTPUT_BUFFER_BYTES.min(bytes_cap_for_init);
        // FAST PATH: zero spill slots + exactly one in-memory chunk → the chunk is
        // already globally sorted, so skip the (k = 1) loser-tree merge and gather
        // it directly. This is the in-memory regime's dominant cost.
        if slots.is_empty() && memory_chunks.total_len() == 1 {
            let chunk = memory_chunks.into_single();
            let total = chunk.len();
            // Pooled gather when a fan-out width was requested, a `MergeWork`
            // inbox is wired, AND there is enough work to amortise the fan-out
            // (a tiny chunk is faster serial: the block plan and the job round
            // trips would dominate). Otherwise fall through to the serial
            // gather — byte-identical, just cheaper for small inputs.
            if self.inbox.is_some() && self.fanout_width > 1 && total >= self.fast_path_min_records
            {
                let blocks = plan_fast_path_blocks::<O>(
                    &chunk,
                    self.target_batch_count,
                    self.output_byte_limit,
                );
                // A plan that collapses to a single block (or none) has no
                // parallelism to exploit, so the job round trip would be pure
                // overhead. Fall through to the serial gather, which
                // emits byte-identically. (`plan_fast_path_blocks` only borrows
                // `chunk`, so it is still owned here.) Unreachable under the
                // production defaults — `target_batch_count` (1024) caps every
                // block, so `total >= FAST_PATH_PARALLEL_MIN_RECORDS` (64Ki)
                // always yields >= 64 blocks — but guards a degenerate plan if
                // those knobs ever become independently tunable.
                if blocks.len() > 1 {
                    self.state = SortMergeState::GatherMode {
                        chunk: Arc::new(chunk),
                        blocks,
                        next_to_post: 0,
                        next_to_take: 0,
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
        let driver = match &self.merge_inputs {
            Some(mode) => build_merge_inputs_driver(
                self.sort_order,
                {
                    self.merge_input_slots.clone_from(&slots);
                    slots
                },
                self.demand.as_ref(),
                self.class_table.as_ref(),
                mode,
            )?,
            None => build_driver(
                self.sort_order,
                slots,
                memory_chunks,
                total_records,
                self.demand.as_ref(),
                self.class_table.as_ref(),
            )?,
        };
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
    /// Branch 0: the merged output (declared `ByItemOrdinal`, collapsed to a
    /// direct edge for this single-producer step). Branch 1: gather, range and
    /// assembly jobs for the `MergeWork` pool step; empty unless the fast path
    /// fans out or the merge runs in range mode.
    type Outputs = (O::Item, MergeTicket<O>);

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
            // residual last), and tie-break are untouched. SortMerge is built by
            // the sort chain's `add_sort` and by `fgumi merge`'s `build_merge`
            // (with merge inputs); both place it on this detached driver thread.
            kind: StepKind::Detached,
            sticky: false,
            output_queues: vec![
                QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit },
                // Fixed: a window of descriptors, not a share of the queue
                // budget, so the data edges keep their full share.
                QueueSpec::FixedByteBounded { limit_bytes: self.job_edge_limit() },
            ],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal, BranchOrdering::None],
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

    fn side_channels(&self) -> Vec<Arc<dyn SideChannel>> {
        self.inbox.iter().map(|inbox| Arc::clone(inbox) as Arc<dyn SideChannel>).collect()
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if !self.inbox_bound {
            if let Some(inbox) = &self.inbox {
                inbox.bind_owner(ctx.pool.waker());
            }
            self.inbox_bound = true;
        }
        // Range mode handles its own held output: it keeps merging ahead (in
        // order, within its window) while the output refuses a block.
        if !matches!(self.state, SortMergeState::RangeMode(_)) && !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        // Whether this dispatch popped input (the setup absorb): if so it must
        // report `Progress` whatever the emit below reports, so the plan's
        // reverse wake reaches an upstream producer holding an item for the slot
        // the pop freed (the `StepOutcome::Progress` contract).
        let mut took_input = false;
        if matches!(&self.state, SortMergeState::WaitingForSetup { .. }) {
            // Drain the input queue unbounded: the setup absorb is cheap (it just
            // moves `Arc`s/`Vec`s into the setup state) and the upstream queue is
            // already byte-bounded, so backpressure belongs on the producer, not
            // on a consumer-side drain cap. A cap here would cycle a full upstream
            // queue through repeated partial drains and add producer contention.
            let absorbed = self.absorb_events_into_setup(ctx)?;
            took_input = absorbed > 0;
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
            SortMergeState::GatherMode { .. } => {
                Ok(self.emit_gather_mode(ctx, &mut records_this_call))
            }
            SortMergeState::RangeMode(_) => self.emit_range_mode(ctx, &mut records_this_call),
            SortMergeState::Done => Ok(StepOutcome::Finished),
            SortMergeState::WaitingForSetup { .. } => {
                unreachable!("Phase 1 must have left state non-WaitingForSetup")
            }
        };
        ctx.counters.add(RECORDS, records_this_call);
        let outcome = match outcome {
            Ok(StepOutcome::NoProgress | StepOutcome::Contention | StepOutcome::Capped)
                if took_input =>
            {
                Ok(StepOutcome::Progress)
            }
            other => other,
        };
        self.settle_await(&outcome);
        outcome
    }
}

#[cfg(test)]
mod tests;
