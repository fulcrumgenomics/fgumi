//! `SortBuffer` — first step of the P6 Phase-1 split (`SortBuffer` →
//! `CompressSpill` → `SortSpillDecompress` → `SortMerge`).
//!
//! `SortBuffer` (`Serial`) ingests `RecordBatch`es into an in-memory arena via
//! the order-erased `ChunkSorter` (each variant drives the same per-order
//! [`ArenaSortStrategy`] as the block-input `FindBoundariesAndSort` front),
//! sorts each filled arena, and emits the sorted records as a [`SortChunkEvent`]
//! — **without touching disk**. Mid-stream spill
//! chunks (`Spill { seq, .. }`) and the final residual (`Residual`) flow to the
//! `Parallel` `CompressSpill` step, which compresses + writes spills and passes
//! the residual through. This replaces the monolithic `SortAndSpill`, which
//! drove a private `SortWorkerPool` to compress inline.
//!
//! Spill chunks are emitted **as the buffer fills** (not accumulated to the
//! end), and `try_run` drains the staged-events queue before popping the next
//! input batch, so staged chunks never accumulate *across* batches. Within a
//! single batch, `ingest_one_batch` seals one chunk each time the arena reaches
//! `memory_limit`. In production this fires at most once per batch — each input
//! `RecordBatch` is block-bounded (`ParseBamRecords` emits one batch per
//! decompressed BGZF block, orders of magnitude below the 512 MB default
//! `memory_limit`) — so peak memory stays at ~one spill chunk plus the live
//! buffer. If `memory_limit` is configured far below the batch size (as some
//! tests do to force spilling), a single batch may seal several chunks into
//! `pending` in one `ingest_one_batch` call; this only raises transient peak
//! memory — every sealed chunk is still emitted in order and no records are
//! dropped.
//!
//! All four sort orders (coordinate, template-coordinate, queryname lex +
//! natural) route through this step via the `ChunkSorter` order enum.

use std::collections::VecDeque;
use std::io;
use std::sync::Arc;

use anyhow::{Result, anyhow};
use fgumi_bam_io::ProgressTracker;
use fgumi_sort::{
    PooledSegmentedBuf, QuerynameComparator, RawExternalSorter, RawQuerynameKey,
    RawQuerynameLexKey, SegmentedBuf, SortOrder, TemplateArenaAccumulator,
};
use noodles::sam::Header;

use crate::sort::protocol::{MemoryChunkErased, SortChunkEvent};
use crate::sort::{ArenaSortStrategy, CoordinateStrategy, QuerynameStrategy, TemplateStrategy};
use crate::types::RecordBatch;
use fgumi_pipeline_core::{
    HeldRetry, Unpushed,
    held::HeldSlot,
    outputs::Single,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// Max staged events flushed to the output per `try_run` invocation.
const MAX_EVENTS_PER_LOCK: usize = 8;

/// Per-record memory overhead added to the arena byte count so the seal
/// threshold accounts for the strategy's `(key, offset, len)` ref alongside the
/// record body. Conservative; the exact value only shifts spill boundaries, not
/// the merged output (coordinate/template are globally stable; queryname's tie
/// order is unspecified).
const PER_RECORD_REF_OVERHEAD: usize = 64;

/// Record-input arena accumulator: copies each pushed record into a growing
/// [`SegmentedBuf`] and drives an [`ArenaSortStrategy`] over the arena refs, so
/// the record-input (SAM / fused) path uses the SAME per-order sort engine as the
/// block-input arena front ([`FindBoundariesAndSort`](crate::sort::FindBoundariesAndSort)).
/// At seal the filled arena is wrapped in an `Arc` and handed to the strategy (no
/// further record copies), and a fresh arena starts the next run.
struct ArenaAccum<S: ArenaSortStrategy> {
    strategy: S,
    arena: SegmentedBuf,
    memory_limit: usize,
    memory_used: usize,
    total_records: u64,
    sort_threads: usize,
}

impl<S: ArenaSortStrategy> ArenaAccum<S> {
    fn new(strategy: S, memory_limit: usize, sort_threads: usize) -> Self {
        Self {
            strategy,
            arena: SegmentedBuf::new(),
            memory_limit,
            memory_used: 0,
            total_records: 0,
            sort_threads,
        }
    }

    /// Copy one record into the arena and accumulate its sort ref. Returns `true`
    /// once the run's byte budget is reached.
    #[allow(clippy::cast_possible_truncation)] // offset/len fit usize on all supported (64-bit) targets
    fn push(&mut self, bam_bytes: &[u8]) -> Result<bool> {
        let len =
            u32::try_from(bam_bytes.len()).map_err(|_| anyhow!("record length exceeds u32"))?;
        let offset = self.arena.extend_from_slice(bam_bytes) as u64;
        let body = self.arena.slice(offset as usize, len as usize);
        // The strategy reads `body` to extract the sort key and stores only
        // `(key, offset, len)`; it does not retain `body`, so the arena may grow
        // (and this slice's borrow end) freely afterwards.
        self.strategy.push_record(body, offset, len).map_err(|e| anyhow!("{e:#}"))?;
        self.memory_used += bam_bytes.len() + PER_RECORD_REF_OVERHEAD;
        self.total_records += 1;
        Ok(self.memory_used >= self.memory_limit)
    }

    /// Seal the current run: wrap the filled arena in an `Arc`, sort + materialize
    /// via the strategy, and reset for the next run. Empty if nothing was pushed
    /// since the last seal.
    fn take_sorted_chunk(&mut self) -> MemoryChunkErased {
        let arena = std::mem::replace(&mut self.arena, SegmentedBuf::new());
        self.memory_used = 0;
        let arc = Arc::new(PooledSegmentedBuf::unpooled(arena));
        self.strategy.seal(arc, self.sort_threads)
    }

    fn total_records(&self) -> u64 {
        self.total_records
    }
}

/// Order-erased record-input arena sorter. Each variant pairs an [`ArenaAccum`]
/// with the concrete [`ArenaSortStrategy`] for its order, so `SortBuffer` drives
/// the same per-order sort engine as the block-input `FindBoundariesAndSort`.
enum ChunkSorter {
    Coordinate(ArenaAccum<CoordinateStrategy>),
    Template(ArenaAccum<TemplateStrategy>),
    QuerynameLex(ArenaAccum<QuerynameStrategy<RawQuerynameLexKey>>),
    QuerynameNatural(ArenaAccum<QuerynameStrategy<RawQuerynameKey>>),
}

impl ChunkSorter {
    /// Build the arena sorter matching `sorter.sort_order()`, provisioning each
    /// order's strategy exactly as the block-input arena front does.
    #[allow(clippy::needless_pass_by_value)] // by-value keeps the caller's move-in; only read here
    fn from_sorter(sorter: RawExternalSorter, header: &Header) -> Result<Self> {
        let memory_limit = sorter.memory_limit_bytes();
        // Phase-1 count honors the `--sort-threads` override (falls back to
        // `--threads`); `num_threads()` would drop the override silently.
        let sort_threads = sorter.phase1_threads();
        Ok(match sorter.sort_order() {
            SortOrder::Coordinate => {
                let n_ref = u32::try_from(header.reference_sequences().len())
                    .map_err(|_| anyhow!("reference sequence count overflows u32"))?;
                Self::Coordinate(ArenaAccum::new(
                    CoordinateStrategy::new(n_ref),
                    memory_limit,
                    sort_threads,
                ))
            }
            SortOrder::TemplateCoordinate => {
                let acc = TemplateArenaAccumulator::from_header(
                    header,
                    sorter.cell_tag_value(),
                    sorter.key_types_spec(),
                );
                Self::Template(ArenaAccum::new(
                    TemplateStrategy::new(acc),
                    memory_limit,
                    sort_threads,
                ))
            }
            SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
                Self::QuerynameLex(ArenaAccum::new(
                    QuerynameStrategy::new(MemoryChunkErased::QuerynameLex),
                    memory_limit,
                    sort_threads,
                ))
            }
            SortOrder::Queryname(QuerynameComparator::Natural) => {
                Self::QuerynameNatural(ArenaAccum::new(
                    QuerynameStrategy::new(MemoryChunkErased::QuerynameNatural),
                    memory_limit,
                    sort_threads,
                ))
            }
        })
    }

    /// Push one record; `true` once the run hits the memory limit.
    fn push(&mut self, bam_bytes: &[u8]) -> Result<bool> {
        match self {
            Self::Coordinate(s) => s.push(bam_bytes),
            Self::Template(s) => s.push(bam_bytes),
            Self::QuerynameLex(s) => s.push(bam_bytes),
            Self::QuerynameNatural(s) => s.push(bam_bytes),
        }
    }

    /// Seal + materialize the current run into one erased chunk (empty if nothing
    /// was pushed), resetting for the next run.
    fn take_sorted_chunk(&mut self) -> MemoryChunkErased {
        match self {
            Self::Coordinate(s) => s.take_sorted_chunk(),
            Self::Template(s) => s.take_sorted_chunk(),
            Self::QuerynameLex(s) => s.take_sorted_chunk(),
            Self::QuerynameNatural(s) => s.take_sorted_chunk(),
        }
    }

    /// Take the final residual as zero-or-one erased chunk. Every order seals one
    /// globally-sorted chunk per run (coordinate/template are stable; queryname's
    /// tie order is unspecified), so the residual is at most one chunk — the
    /// legacy multi-chunk `par_chunks_mut` split is no longer needed.
    ///
    /// An empty residual is normally dropped. The one EXCEPTION is
    /// template-coordinate when there were prior spills: the empty chunk still
    /// carries the chosen `--key-types` narrow-lane variant, which is the only
    /// signal `SortMerge`'s `build_driver` has to pick the spill files' key width.
    /// Without it, a run whose records all spilled (empty final residual) would
    /// leave the merge to default to the full 40-byte key and misread the
    /// narrow-key spills. Coordinate and queryname have a fixed key type, so their
    /// empty residual is dropped as before; with no spills there is nothing to
    /// disambiguate, so the fast path (single non-empty chunk) is preserved.
    fn take_residual_chunks(&mut self, had_spills: bool) -> Vec<MemoryChunkErased> {
        residual_chunks_for(self.take_sorted_chunk(), had_spills)
    }

    fn total_records(&self) -> u64 {
        match self {
            Self::Coordinate(s) => s.total_records(),
            Self::Template(s) => s.total_records(),
            Self::QuerynameLex(s) => s.total_records(),
            Self::QuerynameNatural(s) => s.total_records(),
        }
    }
}

/// Decide whether a run's final (residual) chunk is emitted.
///
/// Split out of `ChunkSorter::take_residual_chunks` so the rule can be tested
/// without constructing a whole accumulator: the interesting behaviour is a pure
/// function of the chunk and whether the run spilled.
///
/// An empty residual is normally dropped. The one EXCEPTION is
/// template-coordinate when there were prior spills: the empty chunk is the only
/// carrier of the chosen `--key-types` narrow-lane variant, which is the only
/// signal `SortMerge`'s `build_driver` has to pick the spill files' key width.
/// Dropping it leaves the merge defaulting to the full 40-byte key and misreading
/// the narrow-key spills — wrong output, not a crash.
fn residual_chunks_for(chunk: MemoryChunkErased, had_spills: bool) -> Vec<MemoryChunkErased> {
    let keep_empty_for_variant =
        had_spills && matches!(chunk, MemoryChunkErased::TemplateCoordinate(_));
    if chunk.is_empty() && !keep_empty_for_variant { Vec::new() } else { vec![chunk] }
}

/// Push every record in `batch` into `sorter`, staging a sealed `Spill` chunk
/// into `pending` (and bumping `next_seq`) each time the arena fills. Returns
/// the number of records traversed and the first push failure, if any.
///
/// The sorter is BORROWED, so a failure cannot leave the caller's sorter slot
/// empty — a state that would be indistinguishable from "finalized". The error
/// is returned rather than propagated with `?` so the caller can record it and
/// still account for the records it already consumed.
///
/// Ingest stops at the first failing record: the pushes that preceded it are
/// retained (they are already in the arena), but continuing would ingest records
/// *after* a rejected one and quietly sort a subset of the input.
fn ingest_batch_records(
    sorter: &mut ChunkSorter,
    batch: &RecordBatch,
    pending: &mut VecDeque<SortChunkEvent>,
    next_seq: &mut u32,
) -> (u64, Option<String>) {
    let mut batch_records = 0u64;
    for record in batch.iter_record_bytes() {
        batch_records += 1;
        let buffer_full = match sorter.push(record) {
            Ok(full) => full,
            Err(e) => return (batch_records, Some(format!("SortBuffer: push failed: {e:#}"))),
        };
        if buffer_full {
            // Seal the filled arena and stage it. We keep draining the rest
            // of the batch rather than stopping early — breaking here would
            // strand the batch's remaining records. Staging multiple chunks
            // per batch is correct (each is emitted in order downstream); in
            // production a block-bounded `RecordBatch` is far below
            // `memory_limit`, so this fires at most once per batch (see the
            // module docs).
            let chunk = sorter.take_sorted_chunk();
            if !chunk.is_empty() {
                pending.push_back(SortChunkEvent::Spill {
                    seq: *next_seq,
                    chunk,
                    records_ingested_so_far: sorter.total_records(),
                });
                *next_seq += 1;
            }
        }
    }
    (batch_records, None)
}

/// `Serial` step that buffers, sorts, and emits sorted chunks for `CompressSpill`.
pub struct SortBuffer {
    /// In-memory buffering sorter. `Some` while ingesting; `None` after the
    /// residual has been taken (finalized).
    sorter: Option<ChunkSorter>,
    /// Sorted chunks awaiting output (spill chunks during ingest; the residual +
    /// `AllAnnounced` after finalize). Drained before the next batch is ingested.
    pending: VecDeque<SortChunkEvent>,
    /// Monotonic spill index. Each `Spill` event's `seq` (= the eventual slot
    /// `file_id`) makes the merge tie-break order independent of which
    /// `CompressSpill` worker writes the file. Also the final `slot_count`.
    next_seq: u32,
    /// First ingest failure, if any. Once set the step is poisoned and every
    /// later `try_run` re-raises instead of ingesting or finalizing — a step
    /// that failed must never go on to report `Finished`, which downstream
    /// would read as a complete sort (see `try_run`). Holds the message rather
    /// than the `io::Error` because `io::Error` is not `Clone`.
    failed: Option<String>,
    held: HeldSlot<Unpushed<SortChunkEvent>>,
    output_byte_limit: u64,
    affinity: Affinity,
    /// Read+ingest progress, logged every 1M records (mirrors the legacy
    /// "Read records" tracker) so the ingest rate over time is visible under
    /// `RUST_LOG=info` — distinguishes a slow read path from a stalled one.
    ingest_progress: ProgressTracker,
}

impl SortBuffer {
    /// Build a `SortBuffer` from a configured `RawExternalSorter` (any of the
    /// four sort orders) and the output `Header`.
    ///
    /// `output_byte_limit` byte-bounds the output event queue (its chunk-bearing
    /// variants retain sorted records, so the queue budgets on bytes, not count).
    ///
    /// # Errors
    ///
    /// Returns an error if the header's reference-sequence count does not fit in
    /// a `u32` (the coordinate key's reference field). That conversion is the
    /// only fallible step: the template path's `TemplateArenaAccumulator::from_header`
    /// is infallible here.
    pub fn from_sorter(
        sorter: RawExternalSorter,
        header: &Header,
        output_byte_limit: u64,
    ) -> Result<Self> {
        let chunk_sorter = ChunkSorter::from_sorter(sorter, header)?;
        Ok(Self {
            sorter: Some(chunk_sorter),
            pending: VecDeque::new(),
            next_seq: 0,
            failed: None,
            held: HeldSlot::new(),
            output_byte_limit,
            affinity: Affinity::None,
            ingest_progress: ProgressTracker::new("Sort ingest records").with_interval(1_000_000),
        })
    }

    /// Override the affinity hint.
    #[must_use]
    pub fn with_affinity(mut self, affinity: Affinity) -> Self {
        self.affinity = affinity;
        self
    }

    /// Re-raise the first ingest failure once the step has been poisoned.
    ///
    /// # Errors
    ///
    /// Returns an error naming the original failure if `failed` is set.
    fn check_not_failed(&self) -> io::Result<()> {
        match &self.failed {
            None => Ok(()),
            Some(message) => Err(io::Error::other(format!(
                "SortBuffer: refusing to continue after an earlier failure: {message}"
            ))),
        }
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        // `true` once the slot is clear (was empty, or the held event flushed);
        // `false` while it's still held under backpressure. Uses the canonical
        // re-hold helper so the put-back-on-reject invariant lives in one place.
        !matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld)
    }

    /// Push up to `MAX_EVENTS_PER_LOCK` staged events to the output, parking the
    /// first that can't be pushed in `held`. Caller guarantees `held` is empty.
    fn emit_pending(&mut self, ctx: &mut StepCtx<'_, Self>) -> StepOutcome {
        let mut emitted = 0usize;
        while emitted < MAX_EVENTS_PER_LOCK {
            let Some(event) = self.pending.pop_front() else { break };
            if let Err(unpushed) = ctx.outputs.push(event) {
                self.held.put(unpushed);
                return StepOutcome::Progress;
            }
            emitted += 1;
        }
        if emitted > 0 { StepOutcome::Progress } else { StepOutcome::NoProgress }
    }

    /// Pop one input batch (if any) and push its records into the sorter, staging
    /// a `Spill` chunk into `pending` whenever the buffer fills. Returns `true`
    /// if a batch was consumed.
    ///
    /// # Errors
    ///
    /// Returns the first `ChunkSorter::push` failure, after poisoning the step so
    /// no later `try_run` can finalize (see the `failed` field).
    ///
    /// # Panics
    ///
    /// Panics if called after finalize (`self.sorter` is `None`).
    fn ingest_one_batch(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<bool> {
        let Some(batch) = ctx.input.pop() else {
            return Ok(false);
        };
        // Borrow the sorter rather than taking it: an error path that left
        // `self.sorter` as `None` would be indistinguishable from "finalized".
        // `pending` / `next_seq` are passed alongside as disjoint borrows.
        let sorter = self.sorter.as_mut().expect("ingest_one_batch after finalize");
        let (batch_records, push_error) =
            ingest_batch_records(sorter, &batch, &mut self.pending, &mut self.next_seq);
        // Log ingest progress (every 1M records) so the read+ingest rate over
        // wall time is visible — a steady rate means the read path is the cost;
        // a bursty/stalling rate means downstream backpressure.
        self.ingest_progress.log_if_needed(batch_records);
        if let Some(message) = push_error {
            self.failed = Some(message.clone());
            return Err(io::Error::other(message));
        }
        Ok(true)
    }

    /// Take the residual chunk and enqueue it (if non-empty) followed by the
    /// terminal `AllAnnounced`. Consumes the sorter (`self.sorter` becomes
    /// `None`), freeing the buffer + rayon pool.
    ///
    /// # Panics
    ///
    /// Panics if called twice (`self.sorter` already `None`).
    fn finalize(&mut self) {
        let mut sorter = self.sorter.take().expect("finalize called twice");
        let total_records = sorter.total_records();
        // `had_spills` keeps an EMPTY template-coordinate residual alive: it is
        // the only carrier of the `--key-types` narrowed-lane variant, without
        // which `SortMerge` falls back to the full 40-byte key and misreads the
        // spills. See `take_residual_chunks` for the authoritative rule.
        let residual_chunks = sorter.take_residual_chunks(self.next_seq > 0);
        let memory_chunk_count =
            u32::try_from(residual_chunks.len()).expect("residual chunk count fits u32");
        for chunk in residual_chunks {
            self.pending.push_back(SortChunkEvent::Residual {
                chunk,
                records_ingested_so_far: total_records,
            });
        }
        self.pending.push_back(SortChunkEvent::AllAnnounced {
            slot_count: self.next_seq,
            memory_chunk_count,
            total_records,
        });
        // `sorter` dropped here (releases the buffer + private rayon pool).
    }
}

impl Step for SortBuffer {
    type Input = RecordBatch;
    type Outputs = Single<SortChunkEvent>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SortBuffer",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn affinity(&self) -> Affinity {
        self.affinity
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // Fail closed. Today both runtimes stop dispatching a step that returned
        // `Err` (`record_error` marks the pipeline done), so this is unreachable
        // — but a scheduler that re-dispatched instead would otherwise fall
        // through to `finalize()` and publish a residual + `AllAnnounced` for a
        // sort that never ingested the rest of its input: a truncated result
        // that is structurally indistinguishable from a complete one.
        self.check_not_failed()?;

        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        // Drain staged events before ingesting more — keeps peak memory bounded
        // to ~one spill chunk plus the live buffer.
        if !self.pending.is_empty() {
            return Ok(self.emit_pending(ctx));
        }

        if self.sorter.is_some() {
            if self.ingest_one_batch(ctx)? {
                // Emit anything the batch just staged.
                if !self.pending.is_empty() {
                    return Ok(self.emit_pending(ctx));
                }
                return Ok(StepOutcome::Progress);
            }
            // No batch available right now.
            if !ctx.input.is_drained() {
                return Ok(StepOutcome::NoProgress);
            }
            // Input fully drained — produce the residual + AllAnnounced.
            self.finalize();
            return Ok(self.emit_pending(ctx));
        }

        // Finalized and all events drained.
        Ok(StepOutcome::Finished)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sort::tests::record_with_mi;
    use crate::types::RecordBatchBuilder;
    use fgumi_sort::{
        InMemoryChunk, KeyTypesSpec, RawCoordinateKey, TemplateKey24, TemplateMemChunk,
    };
    use rstest::rstest;

    fn coordinate(payloads: Vec<Vec<u8>>) -> MemoryChunkErased {
        MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(
            payloads.into_iter().map(|b| (RawCoordinateKey::default(), b)).collect(),
        ))
    }

    fn template(payloads: Vec<Vec<u8>>) -> MemoryChunkErased {
        MemoryChunkErased::TemplateCoordinate(TemplateMemChunk::K24(
            InMemoryChunk::from_owned_records(
                payloads.into_iter().map(|b| (TemplateKey24::default(), b)).collect(),
            ),
        ))
    }

    /// The empty template-coordinate residual must survive `had_spills`, because
    /// it is the only carrier of the narrowed-lane variant. Every other empty
    /// residual is dropped, and every non-empty residual is kept.
    #[rstest]
    #[case::empty_coordinate_no_spills(coordinate(vec![]), false, 0)]
    #[case::empty_coordinate_with_spills(coordinate(vec![]), true, 0)]
    #[case::empty_template_no_spills(template(vec![]), false, 0)]
    // The exception: kept purely to carry the key-width variant to the merge.
    #[case::empty_template_with_spills(template(vec![]), true, 1)]
    #[case::nonempty_coordinate_no_spills(coordinate(vec![vec![1u8; 8]]), false, 1)]
    #[case::nonempty_coordinate_with_spills(coordinate(vec![vec![1u8; 8]]), true, 1)]
    #[case::nonempty_template_no_spills(template(vec![vec![1u8; 8]]), false, 1)]
    #[case::nonempty_template_with_spills(template(vec![vec![1u8; 8]]), true, 1)]
    fn residual_chunk_emission_rule(
        #[case] chunk: MemoryChunkErased,
        #[case] had_spills: bool,
        #[case] expected_len: usize,
    ) {
        assert_eq!(residual_chunks_for(chunk, had_spills).len(), expected_len);
    }

    // ── Ingest failure handling ─────────────────────────────────────────────

    fn batch_of(records: &[Vec<u8>]) -> RecordBatch {
        let total: usize = records.iter().map(Vec::len).sum();
        let mut builder = RecordBatchBuilder::with_capacity(0, total, records.len());
        for record in records {
            builder.push_record_bytes(record);
        }
        builder.build()
    }

    /// Template-coordinate sorter with every optional lane dropped
    /// (`--key-types none`), so the first record fixes the narrowed variant and
    /// any later record whose MI differs is rejected — the one reachable
    /// `ChunkSorter::push` failure.
    fn template_sorter_dropping_mi(memory_limit: usize) -> ChunkSorter {
        let sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .memory_limit(memory_limit)
            .threads(1)
            .key_types(KeyTypesSpec::None);
        ChunkSorter::from_sorter(sorter, &Header::default()).expect("build template chunk sorter")
    }

    /// A push failure must stop ingest at the offending record and leave the
    /// caller's sorter intact. `SortBuffer` distinguishes "ingesting" from
    /// "finalized" solely by `sorter.is_some()`, so a failure that consumed the
    /// sorter would let a later `try_run` report `Finished` — a truncated sort
    /// that looks complete.
    #[test]
    fn ingest_batch_records_stops_at_the_failure_and_keeps_the_sorter() {
        let mut sorter = template_sorter_dropping_mi(256 * 1024 * 1024);
        let mut pending = VecDeque::new();
        let mut next_seq = 0u32;
        let records = vec![
            record_with_mi(10, b"r1", 1),
            record_with_mi(20, b"r2", 2), // differing MI — rejected
            record_with_mi(30, b"r3", 1),
        ];

        let (traversed, error) =
            ingest_batch_records(&mut sorter, &batch_of(&records), &mut pending, &mut next_seq);

        let message = error.expect("a differing MI under --key-types none must be rejected");
        assert!(message.starts_with("SortBuffer: push failed"), "unexpected message: {message}");
        assert_eq!(traversed, 2, "ingest stops at the offending record, not after the batch");
        assert_eq!(sorter.total_records(), 1, "only the accepted record reached the arena");
        // Borrowed, never moved out: the sorter is still usable afterwards.
        assert_eq!(sorter.take_sorted_chunk().len(), 1);
        assert!(pending.is_empty());
        assert_eq!(next_seq, 0);
    }

    /// The clean path: each time the arena reaches `memory_limit` the sealed
    /// chunk is staged and `next_seq` advances. A 1-byte limit seals per record.
    #[test]
    fn ingest_batch_records_stages_a_spill_chunk_each_time_the_arena_fills() {
        let mut sorter = template_sorter_dropping_mi(1);
        let mut pending = VecDeque::new();
        let mut next_seq = 0u32;
        let records = vec![record_with_mi(10, b"r1", 1), record_with_mi(20, b"r2", 1)];

        let (traversed, error) =
            ingest_batch_records(&mut sorter, &batch_of(&records), &mut pending, &mut next_seq);

        assert!(error.is_none(), "constant MI is not a dropped-lane violation");
        assert_eq!(traversed, 2);
        assert_eq!(pending.len(), 2, "one staged spill chunk per seal");
        assert_eq!(next_seq, 2);
    }

    /// A poisoned step re-raises on every later `try_run` instead of falling
    /// through to `finalize()`, and names the original failure so the re-raise
    /// is not mistaken for a second, unrelated error.
    #[rstest]
    #[case::clean(None, None)]
    #[case::poisoned(
        Some("SortBuffer: push failed: dropped lane MI"),
        Some(
            "refusing to continue after an earlier failure: SortBuffer: push failed: dropped lane MI"
        )
    )]
    fn check_not_failed_re_raises_the_original_failure(
        #[case] failed: Option<&str>,
        #[case] expected_message: Option<&str>,
    ) {
        let mut step = SortBuffer::from_sorter(
            RawExternalSorter::new(SortOrder::Coordinate).memory_limit(1 << 20).threads(1),
            &Header::default(),
            1 << 20,
        )
        .expect("build SortBuffer");
        step.failed = failed.map(str::to_string);

        match expected_message {
            None => step.check_not_failed().expect("a step that never failed continues"),
            Some(expected) => {
                let err = step.check_not_failed().expect_err("a poisoned step must re-raise");
                assert!(err.to_string().contains(expected), "unexpected error: {err}");
            }
        }
    }
}
