//! The split arena front: boundary scan, parallel key extraction, run sort.
//!
//! Every sort order's block-input front is the same three steps, generic over
//! an [`ArenaKeyExtractor`]:
//!
//! ```text
//! InflateToArena ─InflatedBlock─▶ FindRecordBoundaries<X>  (Detached, sort-coord)
//!                                   cuts a KeyBatch every KEY_BATCH_RECORDS records
//!                ─KeyBatch<X>───▶ ExtractKeys<X>          (Parallel, ExcludeReader, phase-1 cap)
//!                                   one ref per record, lowest violation per batch
//!                ─RefBatch<X>───▶ SortRun<X>              (Detached, sort-coord, ordinal order)
//!                                   appends; at a run's end sorts and emits Spill / Residual
//! ```
//!
//! `FindRecordBoundaries` owns the arena cursor — header skip, the straddler
//! carry, truncation checks — and never reads a record body except the sort's
//! first, which provisions the extractor's context. It cuts a batch every
//! [`KEY_BATCH_RECORDS`] records and pushes it at once, so extraction runs on
//! the pool while the scan continues. `SortRun` receives the batches in ordinal
//! (input) order, appends their refs and, on the run's last batch, sorts the
//! run on the order's bounded pool and emits exactly the events the run's
//! chunk needs. The run sort is phase-1 work off the pool: with a phase-1 cap
//! (an explicit `--sort-threads`) it first reserves the whole cap
//! (`reserve_for_run_sort` in `arena_ingest`).
//!
//! **Concurrent arena reads.** A worker reads `[byte_start, byte_end)` of a run's
//! arena while inflate workers may still write later slots of the same run. The
//! batch's bytes are all below the scan cursor, which only advances over slots
//! whose `InflatedBlock` token has arrived — the inflate worker's write
//! happens-before the token's push, the scan's pop, the batch's push and the
//! extractor's pop. Inflate slots are disjoint and at or above `FRONT_REGION`;
//! the only write below it is the carry, made before any batch of its run is
//! cut. Bytes below the cursor are never written again during the run.

use std::collections::VecDeque;
use std::io;
use std::sync::Arc;
use std::time::Instant;

use fgumi_pipeline_core::admission::admit_input;
use fgumi_pipeline_core::held::HeldSlot;
use fgumi_pipeline_core::item::{HeapSize, Ordered};
use fgumi_pipeline_core::outputs::{OrderedBytesSingle, Single};
use fgumi_pipeline_core::queues::QueueSpec;
use fgumi_pipeline_core::reorder::BranchOrdering;
use fgumi_pipeline_core::step::{
    CounterSpec, DetachedGroup, PoolPlacement, Step, StepCtx, StepKind, StepOutcome, StepProfile,
};
use fgumi_pipeline_core::{HeldRetry, PhaseCap, Unpushed};
use fgumi_raw_bam::MIN_BAM_RECORD_LEN;
use fgumi_sort::prefetch::{KEY_PREFETCH_DISTANCE, prefetch_read_l1};
use fgumi_sort::{
    BoundedSortPool, ExtractedRefs, InMemoryChunk, KeyViolation, PooledSegmentedBuf, QnBatchLcp,
    QnComparator, QnKey, QnLex, QnNatural, QnRef, QnRefBuilder, QnSortState, RecordRef,
    TemplateArenaRefs, TemplateKeyContext, TemplateKeySeed, TemplateKeyVariant, TemplateMemChunk,
    coordinate_chunk_from_refs, extract_coordinate_key_inline, queryname_chunk_from_qn_refs,
    radix_sorts_in_parallel, seal_template_refs,
};
use log::{debug, info};

use super::arena_ingest::{FRONT_REGION, InflatedBlock, reserve_for_run_sort};
use super::run_ledger::{RUN0_MEAN_RECORD_BYTES, RunLedger, run_cap_for};
use crate::boundaries::bam_header_len;
use crate::sort::protocol::{MemoryChunkErased, SortChunkEvent};

/// Records per key batch: the unit `FindRecordBoundaries` cuts and one
/// `ExtractKeys` call extracts.
pub const KEY_BATCH_RECORDS: usize = 4 * 1024;

/// Order-specific half of the split front: what a key batch extracts, how a
/// run's refs accumulate and how a run seals into a chunk.
pub trait ArenaKeyExtractor: Send + Sync + 'static {
    /// Names of the three steps, in chain order: boundaries, extraction, run sort.
    const STEP_NAMES: [&'static str; 3];
    /// Per-sort configuration known before the first record.
    type Seed: Clone + Send + 'static;
    /// Per-sort extraction context, provisioned from the sort's first record.
    type Ctx: Clone + Send + Sync + 'static;
    /// One batch's refs, in record order.
    type Refs: Send + 'static;
    /// One run's accumulated refs.
    type Accum: Send + 'static;
    /// Per-`ExtractKeys`-clone mutable state (`()` when none).
    type Scratch: Default + Clone + Send + 'static;
    /// The sort's context, from its first record's body.
    fn provision(seed: &Self::Seed, first_body: &[u8]) -> Self::Ctx;
    /// Bytes of one record's ref, when known before the first record (every
    /// order but template order under `--key-types auto`).
    fn known_ref_width(seed: &Self::Seed) -> Option<usize>;
    /// Bytes of one record's ref.
    fn ref_width(ctx: &Self::Ctx) -> usize;
    /// Re-walk `[byte_start, byte_end)` of `arena` as exactly `n`
    /// `[u32 len][body]` frames and extract one ref per record, with the
    /// lowest-indexed violation of the batch.
    ///
    /// # Errors
    ///
    /// `InvalidData` if the frames do not match the batch the scan cut.
    #[allow(clippy::too_many_arguments)]
    fn extract_batch(
        ctx: &Self::Ctx,
        scratch: &mut Self::Scratch,
        arena: &PooledSegmentedBuf,
        byte_start: u64,
        byte_end: u64,
        n: u32,
        first_index: u64,
    ) -> io::Result<(Self::Refs, Option<KeyViolation>)>;
    /// Refs in a batch.
    fn refs_len(refs: &Self::Refs) -> usize;
    /// Bytes a batch's refs occupy.
    fn refs_byte_len(refs: &Self::Refs) -> usize;
    /// An empty run accumulator for `ctx`.
    fn accum_for(ctx: &Self::Ctx) -> Self::Accum;
    /// Reserve room for `n` more refs.
    fn reserve(acc: &mut Self::Accum, n: usize);
    /// Append a batch's refs, in order.
    fn append(acc: &mut Self::Accum, refs: Self::Refs);
    /// Whether sealing `acc` at `sort_threads` sorts on more than one thread —
    /// what decides whether the seal reserves the phase-1 cap.
    fn sorts_in_parallel(acc: &Self::Accum, sort_threads: usize) -> bool;
    /// Sort this run's refs into a chunk over `arena`, at `sort_threads` (the
    /// phase-1 width) on the order's shared bounded pool, leaving `acc` empty
    /// for the next run. When phase 1 is capped, the caller holds the whole cap
    /// for the duration of this call.
    fn seal(
        acc: &mut Self::Accum,
        ctx: &Self::Ctx,
        arena: Arc<PooledSegmentedBuf>,
        sort_threads: usize,
        pool: &BoundedSortPool,
    ) -> MemoryChunkErased;
    /// The chunk a run sort that saw no record at all seals to.
    fn empty_chunk() -> MemoryChunkErased;
    /// The bounded sort pool's thread-name prefix.
    fn pool_name() -> &'static str;
}

/// A run's end, carried by its last key batch.
#[derive(Debug, Clone, Copy)]
pub struct RunEnd {
    /// Whether the run spills (`Spill`) or is the final residual.
    pub seals_to_spill: bool,
    /// Records scanned through the end of this run.
    pub records_ingested_so_far: u64,
    /// The bytes of the run's complete records (length prefixes included):
    /// neither run 0's BAM header nor a trailing partial record carried into
    /// the next run, so `run_bytes / records` is the run's exact mean.
    pub run_bytes: u64,
    /// When the scan finished the run.
    pub scan_end: Instant,
}

/// `FindRecordBoundaries` → `ExtractKeys`: a contiguous range of whole records
/// of one run.
pub struct KeyBatch<X: ArenaKeyExtractor> {
    /// The run's arena, kept alive while the batch is in flight.
    pub arena: Arc<PooledSegmentedBuf>,
    /// The sort's context; `None` only for a batch of a sort that has seen no
    /// record (then `n_records == 0`).
    pub ctx: Option<X::Ctx>,
    /// Dense across the whole sort, empty run-final batches included.
    pub ordinal: u64,
    /// The run this batch belongs to.
    pub run_seq: u32,
    /// Run-relative index of the batch's first record.
    pub first_record_index: u64,
    /// Records in the batch; 0 only for an empty run-final batch.
    pub n_records: u32,
    /// Arena offset of the first record's length prefix.
    pub byte_start: u64,
    /// Arena offset one past the last record's body.
    pub byte_end: u64,
    /// Set on a run's last batch.
    pub run_end: Option<RunEnd>,
}

impl<X: ArenaKeyExtractor> HeapSize for KeyBatch<X> {
    /// The ref bytes the batch commits the extractor to produce.
    fn heap_size(&self) -> usize {
        self.ctx.as_ref().map_or(0, |ctx| self.n_records as usize * X::ref_width(ctx))
    }
}

impl<X: ArenaKeyExtractor> Ordered for KeyBatch<X> {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

/// `ExtractKeys` → `SortRun`: a batch's refs.
pub struct RefBatch<X: ArenaKeyExtractor> {
    /// The run's arena.
    pub arena: Arc<PooledSegmentedBuf>,
    /// The sort's context (`None` as for [`KeyBatch::ctx`]).
    pub ctx: Option<X::Ctx>,
    /// The batch's ordinal.
    pub ordinal: u64,
    /// The run this batch belongs to.
    pub run_seq: u32,
    /// Run-relative index of the batch's first record.
    pub first_record_index: u64,
    /// The batch's refs, in record order (`None` for an empty batch).
    pub refs: Option<X::Refs>,
    /// The batch's lowest-indexed key violation.
    pub violation: Option<KeyViolation>,
    /// When its keys were extracted.
    pub keyed_at: Instant,
    /// Set on a run's last batch.
    pub run_end: Option<RunEnd>,
}

impl<X: ArenaKeyExtractor> RefBatch<X> {
    /// Refs in the batch.
    fn len(&self) -> usize {
        self.refs.as_ref().map_or(0, X::refs_len)
    }
}

impl<X: ArenaKeyExtractor> HeapSize for RefBatch<X> {
    fn heap_size(&self) -> usize {
        self.refs.as_ref().map_or(0, X::refs_byte_len)
    }
}

impl<X: ArenaKeyExtractor> Ordered for RefBatch<X> {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

// ============================================================================
// FindRecordBoundaries
// ============================================================================

/// `Detached` step on the `sort-coord` driver that scans each run's arena span
/// for record boundaries and cuts them into [`KeyBatch`]es. It owns the run's
/// arena cursor: the BAM header skip on run 0, the straddler carry between
/// runs, contiguity and truncation checks. It reads no record body except the
/// sort's first, which provisions the extractor's context.
pub struct FindRecordBoundaries<X: ArenaKeyExtractor> {
    seed: X::Seed,
    ctx: Option<X::Ctx>,
    batch_records: usize,
    output_byte_limit: u64,
    held: HeldSlot<Unpushed<KeyBatch<X>>>,
    pending: VecDeque<KeyBatch<X>>,
    /// The current run's arena (set on its first block, cleared at its seal).
    arena: Option<Arc<PooledSegmentedBuf>>,
    /// Arena offset where the current run's scan began.
    scan_start: u64,
    /// One past the last inflated byte of the current run.
    run_end: u64,
    /// Arena offset of the next unparsed record.
    scan_cursor: u64,
    /// Arena offset of the current run's first record (past run 0's BAM header;
    /// at the carried straddler's first byte when there is one).
    run_first_record: u64,
    /// `true` once run 0's BAM header has been skipped (always for runs > 0).
    header_skipped: bool,
    /// The previous run's trailing partial record, written into this run's
    /// front region on its first block.
    carry: Vec<u8>,
    /// Arena offset of the current batch's first record.
    batch_start: u64,
    /// Run-relative index of the current batch's first record.
    batch_first_index: u64,
    /// Records in the current batch so far.
    batch_n: u32,
    /// Records scanned across all runs.
    total_records: u64,
    next_ordinal: u64,
    /// The next run's sequence number.
    next_seq: u32,
    finalized: bool,
    /// Receives the ref width when the sort's first record provisions it.
    ledger: Option<Arc<RunLedger>>,
    /// Test probe: may hold the scan after batch 0 until it is extracted.
    #[cfg(test)]
    probe: Option<Arc<SplitProbe>>,
}

impl<X: ArenaKeyExtractor> FindRecordBoundaries<X> {
    /// A boundary step for a sort seeded with `seed`, its output queue bounded
    /// to `output_byte_limit` bytes.
    #[must_use]
    pub fn new(seed: X::Seed, output_byte_limit: u64) -> Self {
        Self {
            seed,
            ctx: None,
            batch_records: KEY_BATCH_RECORDS,
            output_byte_limit,
            held: HeldSlot::new(),
            pending: VecDeque::new(),
            arena: None,
            scan_start: 0,
            run_end: 0,
            scan_cursor: 0,
            run_first_record: 0,
            header_skipped: false,
            carry: Vec::new(),
            batch_start: 0,
            batch_first_index: 0,
            batch_n: 0,
            total_records: 0,
            next_ordinal: 0,
            next_seq: 0,
            finalized: false,
            ledger: None,
            #[cfg(test)]
            probe: None,
        }
    }

    /// Publish the provisioned ref width to the sort's run ledger, which sizes
    /// the next run's arena cap.
    #[must_use]
    pub fn with_run_ledger(mut self, ledger: Arc<RunLedger>) -> Self {
        self.ledger = Some(ledger);
        self
    }

    /// Cut batches of `n` records instead of [`KEY_BATCH_RECORDS`] (tests).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_batch_records(mut self, n: usize) -> Self {
        self.batch_records = n.max(1);
        self
    }

    /// Report to `probe`, and hold the scan after the sort's first batch until
    /// it is extracted when the probe asks (tests).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_probe(mut self, probe: Arc<SplitProbe>) -> Self {
        self.probe = Some(probe);
        self
    }

    /// Cut the current batch (possibly empty when `run_end` is set) and stage it.
    fn cut(&mut self, arena: &Arc<PooledSegmentedBuf>, run_end: Option<RunEnd>) {
        let ordinal = self.next_ordinal;
        self.next_ordinal += 1;
        let batch = KeyBatch {
            arena: Arc::clone(arena),
            ctx: self.ctx.clone(),
            ordinal,
            run_seq: self.next_seq,
            first_record_index: self.batch_first_index,
            n_records: self.batch_n,
            byte_start: self.batch_start,
            byte_end: self.scan_cursor,
            run_end,
        };
        self.batch_first_index += u64::from(self.batch_n);
        self.batch_start = self.scan_cursor;
        self.batch_n = 0;
        self.pending.push_back(batch);
    }

    /// Ingest one inflated block: extend the run's span and scan every record it
    /// completes; seal the run on its last block. Returns the records found.
    ///
    /// # Errors
    ///
    /// A carry longer than `FRONT_REGION`, a frame shorter than a BAM record, a
    /// short header (or one longer than a spilled run), or a truncated final
    /// record.
    pub(crate) fn ingest_block(&mut self, block: &InflatedBlock) -> io::Result<u64> {
        let block_end = block.offset + u64::from(block.len);
        if self.arena.is_none() {
            let scan_start = if self.carry.is_empty() {
                if self.next_seq == 0 {
                    // Run 0: the BAM header starts at the run's first byte.
                    block.offset
                } else {
                    FRONT_REGION as u64
                }
            } else {
                let l = self.carry.len();
                if l > FRONT_REGION {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "{}: carry ({l} bytes) exceeds FRONT_REGION ({FRONT_REGION}); record \
                             exceeds max straddler size — raise --max-memory",
                            X::STEP_NAMES[0]
                        ),
                    ));
                }
                let write_offset = FRONT_REGION - l;
                // SAFETY: `[FRONT_REGION - l, FRONT_REGION)` is the front region
                // `ReadBlocks::ensure_arena` reserved (`grow_uninit` of the whole
                // segment) before any block slot was handed out. Soundness rests on
                // DISJOINTNESS: this range lies strictly below `FRONT_REGION`, while
                // every inflate block slot lies at or above it, so the `&mut [u8]`
                // aliases no slot an `InflateToArena` worker may still be writing —
                // and no key batch of this run is cut before this write, so no
                // extractor reads it concurrently. The segment was frozen with
                // `reserve_full_capacity`, so no grow can move these bytes. Written
                // once per straddler per run; `u8` has no validity invariant.
                #[allow(unsafe_code)]
                let dst = unsafe { block.arena.slice_mut(write_offset, l) };
                dst.copy_from_slice(&self.carry);
                self.carry.clear();
                write_offset as u64
            };
            self.arena = Some(Arc::clone(&block.arena));
            self.scan_start = scan_start;
            self.scan_cursor = scan_start;
            self.batch_start = scan_start;
            self.batch_first_index = 0;
            self.batch_n = 0;
            self.run_first_record = scan_start;
            self.header_skipped = self.next_seq != 0;
        } else {
            debug_assert_eq!(
                block.offset,
                self.run_end,
                "{}: non-contiguous arena block (expected offset {}, got {})",
                X::STEP_NAMES[0],
                self.run_end,
                block.offset
            );
        }
        self.run_end = block_end;
        let arena = Arc::clone(&block.arena);
        let found = self.scan_available(&arena)?;
        if block.is_last_of_run {
            self.seal_run(block.seals_to_spill, block.run_seq)?;
        }
        Ok(found)
    }

    /// Scan every record fully present in `[scan_cursor, run_end)`, cutting a
    /// batch every `batch_records` records, and park the cursor at the first
    /// incomplete record. Returns the records found.
    fn scan_available(&mut self, arena: &Arc<PooledSegmentedBuf>) -> io::Result<u64> {
        let scan_start = self.scan_start;
        let start = usize::try_from(scan_start).expect("scan_start fits usize");
        let end = usize::try_from(self.run_end).expect("run_end fits usize");
        let span = arena.slice(start, end.checked_sub(start).expect("run_end >= scan_start"));

        if !self.header_skipped {
            match bam_header_len(span)? {
                Some(h) => {
                    self.scan_cursor = scan_start + h as u64;
                    self.batch_start = self.scan_cursor;
                    self.run_first_record = self.scan_cursor;
                    self.header_skipped = true;
                }
                None => return Ok(0),
            }
        }

        let mut cur =
            usize::try_from(self.scan_cursor - scan_start).expect("cursor offset fits usize");
        let mut found: u64 = 0;
        loop {
            if cur + 4 > span.len() {
                break;
            }
            let bs = u32::from_le_bytes([span[cur], span[cur + 1], span[cur + 2], span[cur + 3]])
                as usize;
            if cur + 4 + bs > span.len() {
                break;
            }
            // A frame shorter than the 32 fixed BAM bytes cannot be a record, and
            // every order's key extractor (and `provision`) indexes those fields.
            if bs < MIN_BAM_RECORD_LEN {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "{}: block_size {bs} is shorter than a BAM record \
                         ({MIN_BAM_RECORD_LEN} bytes)",
                        X::STEP_NAMES[0]
                    ),
                ));
            }
            if let Some(ahead) = span.get(cur + KEY_PREFETCH_DISTANCE) {
                prefetch_read_l1(ahead);
            }
            if self.ctx.is_none() {
                let ctx = X::provision(&self.seed, &span[cur + 4..cur + 4 + bs]);
                if let Some(ledger) = &self.ledger {
                    ledger.publish_width(X::ref_width(&ctx));
                }
                self.ctx = Some(ctx);
            }
            cur += 4 + bs;
            found += 1;
            self.batch_n += 1;
            if self.batch_n as usize == self.batch_records {
                self.scan_cursor = scan_start + cur as u64;
                self.cut(arena, None);
            }
        }
        self.scan_cursor = scan_start + cur as u64;
        self.total_records += found;
        Ok(found)
    }

    /// Seal the current run: its trailing partial record becomes the carry (a
    /// spill) or is a truncation error (the residual); then cut the run's last
    /// batch, carrying the run's end.
    fn seal_run(&mut self, seals_to_spill: bool, run_seq: u32) -> io::Result<()> {
        let arena = self.arena.take().expect("seal_run called with no arena");
        if !self.header_skipped {
            // A spilled run sealed mid-header: the header outgrew one run's arena
            // cap. The final run reaching here is genuinely short input.
            let msg = if seals_to_spill {
                format!(
                    "{}: the BAM header does not fit in one run's memory budget ({} bytes \
                     read without reaching its end) — raise --max-memory",
                    X::STEP_NAMES[0],
                    self.run_end - self.scan_start
                )
            } else {
                format!("{}: input too short to contain a valid BAM header", X::STEP_NAMES[0])
            };
            return Err(io::Error::new(io::ErrorKind::InvalidData, msg));
        }
        let cursor = usize::try_from(self.scan_cursor).expect("scan_cursor fits usize");
        let end = usize::try_from(self.run_end).expect("run_end fits usize");
        let tail_len = end.checked_sub(cursor).expect("run_end >= scan_cursor");
        if seals_to_spill {
            if tail_len > 0 {
                if tail_len > FRONT_REGION {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "{}: carry ({tail_len} bytes) exceeds FRONT_REGION \
                             ({FRONT_REGION}); record exceeds max straddler size — raise \
                             --max-memory",
                            X::STEP_NAMES[0]
                        ),
                    ));
                }
                self.carry.extend_from_slice(arena.slice(cursor, tail_len));
            }
        } else if tail_len > 0 || !self.carry.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "truncated BAM: final record incomplete at end of input",
            ));
        }
        debug_assert_eq!(
            run_seq,
            self.next_seq,
            "ReadBlocks run_seq desynced from {}'s run counter",
            X::STEP_NAMES[0]
        );
        let run_end = RunEnd {
            seals_to_spill,
            records_ingested_so_far: self.total_records,
            run_bytes: self.scan_cursor.saturating_sub(self.run_first_record),
            scan_end: Instant::now(),
        };
        self.cut(&arena, Some(run_end));
        self.next_seq += 1;
        if !seals_to_spill {
            self.finalized = true;
        }
        self.scan_start = 0;
        self.run_end = 0;
        self.scan_cursor = 0;
        self.header_skipped = false;
        Ok(())
    }

    /// Input drained: seal an open run as the residual. No run at all (empty
    /// input) emits nothing; spills with no residual (impossible under
    /// `ReadBlocks`' deferred seal) fail on a non-empty carry and otherwise
    /// leave `SortRun` to announce.
    fn finalize(&mut self) -> io::Result<()> {
        if self.finalized {
            return Ok(());
        }
        if self.arena.is_some() {
            let seq = self.next_seq;
            return self.seal_run(false, seq);
        }
        if self.next_seq > 0 && !self.carry.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "truncated BAM: final record incomplete at end of input (no-residual path)",
            ));
        }
        self.finalized = true;
        Ok(())
    }

    fn emit_pending(&mut self, ctx: &mut StepCtx<'_, Self>) -> StepOutcome {
        let Some(batch) = self.pending.pop_front() else { return StepOutcome::NoProgress };
        if let Err(unpushed) = ctx.outputs.push(batch) {
            self.held.put(unpushed);
        }
        StepOutcome::Progress
    }

    /// Whether the test probe holds the scan.
    #[cfg(test)]
    fn latched(&self) -> bool {
        self.probe
            .as_ref()
            .is_some_and(|p| p.hold_scan && self.next_ordinal > 0 && !p.extracted_first())
    }
}

impl<X: ArenaKeyExtractor> Step for FindRecordBoundaries<X> {
    type Input = InflatedBlock;
    type Outputs = OrderedBytesSingle<KeyBatch<X>>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: X::STEP_NAMES[0],
            kind: StepKind::Detached,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn detached_group(&self) -> DetachedGroup {
        DetachedGroup::Shared(crate::sort::SORT_COORD_GROUP)
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] =
            &[CounterSpec::new("records", "records"), CounterSpec::new("batches", "batches")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld) {
            return Ok(StepOutcome::Contention);
        }
        if !self.pending.is_empty() {
            ctx.counters.add(1, 1);
            return Ok(self.emit_pending(ctx));
        }
        if self.finalized {
            return Ok(StepOutcome::Finished);
        }
        #[cfg(test)]
        if self.latched() {
            return Ok(StepOutcome::NoProgress);
        }
        if let Some(block) = ctx.input.pop() {
            let found = self.ingest_block(&block)?;
            ctx.counters.add(0, found);
            if !self.pending.is_empty() {
                ctx.counters.add(1, 1);
                return Ok(self.emit_pending(ctx));
            }
            return Ok(StepOutcome::Progress);
        }
        if !ctx.input.is_drained() {
            return Ok(StepOutcome::NoProgress);
        }
        self.finalize()?;
        if !self.pending.is_empty() {
            ctx.counters.add(1, 1);
            return Ok(self.emit_pending(ctx));
        }
        Ok(StepOutcome::Finished)
    }
}

// ============================================================================
// ExtractKeys
// ============================================================================

/// `Parallel` step: extracts one batch's refs on a pool worker. Placed off the
/// reader worker by default; under a phase-1 cap it is admitted like every
/// phase-1 pool step, before it pops.
pub struct ExtractKeys<X: ArenaKeyExtractor> {
    cap: Option<Arc<PhaseCap>>,
    scratch: X::Scratch,
    held: HeldSlot<Unpushed<RefBatch<X>>>,
    output_byte_limit: u64,
    placement: PoolPlacement,
    #[cfg(test)]
    probe: Option<Arc<SplitProbe>>,
}

impl<X: ArenaKeyExtractor> ExtractKeys<X> {
    /// An extraction step, its output queue bounded to `output_byte_limit`
    /// bytes.
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self {
            cap: None,
            scratch: X::Scratch::default(),
            held: HeldSlot::new(),
            output_byte_limit,
            placement: PoolPlacement::ExcludeReader,
            #[cfg(test)]
            probe: None,
        }
    }

    /// Share the phase-1 admission cap with the other phase-1 pool steps.
    /// `None` leaves the step uncapped.
    #[must_use]
    pub fn with_phase_cap(mut self, cap: Option<Arc<PhaseCap>>) -> Self {
        self.cap = cap;
        self
    }

    /// Override the pool placement (default [`PoolPlacement::ExcludeReader`]).
    #[must_use]
    pub fn with_placement(mut self, placement: PoolPlacement) -> Self {
        self.placement = placement;
        self
    }

    /// Report each extraction to `probe` (tests).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_probe(mut self, probe: Arc<SplitProbe>) -> Self {
        self.probe = Some(probe);
        self
    }

    fn run_once(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        cap: Option<&PhaseCap>,
    ) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take()
            && let Err(again) = ctx.outputs.retry(unpushed)
        {
            self.held.put(again);
            return Ok(StepOutcome::Contention);
        }
        // Admission before the pop: an empty poll takes no permit; a refusal
        // leaves the batch queued for a clone that may run.
        let _permit = match admit_input(ctx.input, cap) {
            Ok(permit) => permit,
            Err(outcome) => return Ok(outcome),
        };
        let Some(batch) = ctx.input.pop() else {
            return Ok(if ctx.input.is_drained() {
                StepOutcome::Finished
            } else {
                StepOutcome::NoProgress
            });
        };
        let (refs, violation) = match &batch.ctx {
            Some(c) if batch.n_records > 0 => {
                let (refs, violation) = X::extract_batch(
                    c,
                    &mut self.scratch,
                    &batch.arena,
                    batch.byte_start,
                    batch.byte_end,
                    batch.n_records,
                    batch.first_record_index,
                )?;
                (Some(refs), violation)
            }
            _ => (None, None),
        };
        #[cfg(test)]
        if let Some(probe) = &self.probe {
            probe.extracted(batch.ordinal);
        }
        ctx.counters.add(0, u64::from(batch.n_records));
        let out = RefBatch {
            arena: batch.arena,
            ctx: batch.ctx,
            ordinal: batch.ordinal,
            run_seq: batch.run_seq,
            first_record_index: batch.first_record_index,
            refs,
            violation,
            keyed_at: Instant::now(),
            run_end: batch.run_end,
        };
        if let Err(unpushed) = ctx.outputs.push(out) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

impl<X: ArenaKeyExtractor> Step for ExtractKeys<X> {
    type Input = KeyBatch<X>;
    type Outputs = OrderedBytesSingle<RefBatch<X>>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: X::STEP_NAMES[1],
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn pool_placement(&self) -> PoolPlacement {
        self.placement
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] = &[CounterSpec::new("records", "records")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // Moved out for the call: the permit borrows the cap while the body
        // needs `&mut self`. A move, so no reference-count traffic.
        let cap = self.cap.take();
        let outcome = self.run_once(ctx, cap.as_deref());
        self.cap = cap;
        outcome
    }

    fn new_worker_copy(&self) -> Self {
        Self {
            cap: self.cap.clone(),
            scratch: self.scratch.clone(),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            placement: self.placement,
            #[cfg(test)]
            probe: self.probe.clone(),
        }
    }

    fn phase_cap(&self) -> Option<&PhaseCap> {
        self.cap.as_deref()
    }
}

// ============================================================================
// SortRun
// ============================================================================

/// How much of a sort's key extraction overlapped the boundary scan.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct KeyOverlapCensus {
    /// Records whose batch was keyed before its run's scan ended.
    pub overlapped_records: u64,
    /// Records keyed after their run's scan ended.
    pub tail_records: u64,
    /// Σ over runs of the time from the scan's end to the run's last batch.
    pub tail_wait_ns: u64,
    /// Batches.
    pub batches: u64,
    /// Records.
    pub records: u64,
    /// Runs.
    pub runs: u64,
}

impl KeyOverlapCensus {
    /// The census lines, `phase_width` being the resolved phase-1 width.
    #[must_use]
    #[allow(clippy::cast_precision_loss)] // display values
    pub fn lines(&self, phase_width: usize) -> [String; 2] {
        let pct = if self.records == 0 {
            100.0
        } else {
            100.0 * self.overlapped_records as f64 / self.records as f64
        };
        [
            format!(
                "Deferred key extraction: {pct:.1}% overlapped ({} of {} records; {} keyed after \
                 the scan ended)",
                self.overlapped_records, self.records, self.tail_records
            ),
            format!(
                "  tail wait (exact): {:.1}s over {} runs; {} batches of <={KEY_BATCH_RECORDS} \
                 records; phase-1 width {phase_width}",
                self.tail_wait_ns as f64 / 1e9,
                self.runs,
                self.batches
            ),
        ]
    }
}

/// `Detached` step on the `sort-coord` driver: appends each run's batches, in
/// ordinal order, and at the run's last batch sorts the run and emits its
/// chunk — `Spill` for a sealed mid-stream run, `Residual` then `AllAnnounced`
/// for the final run.
pub struct SortRun<X: ArenaKeyExtractor> {
    sort_threads: usize,
    output_byte_limit: u64,
    sort_stats: bool,
    pool: BoundedSortPool,
    cap: Option<Arc<PhaseCap>>,
    ctx: Option<X::Ctx>,
    acc: Option<X::Accum>,
    run_arena: Option<Arc<PooledSegmentedBuf>>,
    held: HeldSlot<Unpushed<SortChunkEvent>>,
    pending: VecDeque<SortChunkEvent>,
    next_seq: u32,
    spilled_run_count: u32,
    total_records: u64,
    finalized: bool,
    census: KeyOverlapCensus,
    run_batches: Vec<(u64, Instant)>,
    census_logged: bool,
    /// The sort's run ledger and in-memory budget: each run's totals are
    /// published to it, and each run's refs are reserved from it.
    budget: Option<(Arc<RunLedger>, usize)>,
    #[cfg(test)]
    probe: Option<Arc<SplitProbe>>,
}

impl<X: ArenaKeyExtractor> SortRun<X> {
    /// A run sorter. `sort_threads` is the resolved phase-1 width
    /// (`PhaseThreads::phase1`, from `--threads`/`--sort-threads`): it sizes the
    /// order's bounded sort pool. With a phase-1 cap (`with_phase_cap`, set only
    /// for an explicit `--sort-threads`), each parallel run sort first reserves
    /// the whole cap.
    #[must_use]
    pub fn new(sort_threads: usize, output_byte_limit: u64) -> Self {
        Self {
            sort_threads: sort_threads.max(1),
            output_byte_limit,
            sort_stats: false,
            pool: BoundedSortPool::new(X::pool_name()),
            cap: None,
            ctx: None,
            acc: None,
            run_arena: None,
            held: HeldSlot::new(),
            pending: VecDeque::new(),
            next_seq: 0,
            spilled_run_count: 0,
            total_records: 0,
            finalized: false,
            census: KeyOverlapCensus::default(),
            run_batches: Vec::new(),
            census_logged: false,
            budget: None,
            #[cfg(test)]
            probe: None,
        }
    }

    /// Publish each run's exact totals to the sort's run ledger (before the
    /// run's chunk is emitted, so before its arena can return to `ReadBlocks`),
    /// and reserve each run's refs for the arena cap the ledger implies under
    /// `memory_limit`.
    #[must_use]
    pub fn with_run_ledger(mut self, ledger: Arc<RunLedger>, memory_limit: usize) -> Self {
        self.budget = Some((ledger, memory_limit));
        self
    }

    /// Publish the final census to `probe` (tests).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_probe(mut self, probe: Arc<SplitProbe>) -> Self {
        self.probe = Some(probe);
        self
    }

    /// Reserve the whole phase-1 cap around each parallel run sort. `None`
    /// leaves the sort unreserved.
    #[must_use]
    pub fn with_phase_cap(mut self, cap: Option<Arc<PhaseCap>>) -> Self {
        self.cap = cap;
        self
    }

    /// Log the overlap census at `info` (`--sort-stats`) instead of `debug`.
    #[must_use]
    pub fn with_sort_stats(mut self, on: bool) -> Self {
        self.sort_stats = on;
        self
    }

    /// The census so far.
    #[must_use]
    pub fn census(&self) -> KeyOverlapCensus {
        self.census
    }

    /// Test-only: the bounded sort pool's thread count, `None` until a seal
    /// builds it.
    #[cfg(test)]
    pub(crate) fn sort_pool_threads(&self) -> Option<usize> {
        self.pool.built_threads()
    }

    /// Append one batch; seal the run on its last batch.
    fn accept(&mut self, batch: RefBatch<X>) -> io::Result<()> {
        debug_assert_eq!(batch.run_seq, self.next_seq, "batch of an unexpected run");
        if let Some(v) = batch.violation {
            // Batches arrive in input order and each carries its lowest-indexed
            // violation, so this is the input's first violation.
            return Err(v.into_io_error());
        }
        let n = batch.len();
        if self.ctx.is_none() {
            self.ctx.clone_from(&batch.ctx);
        }
        if self.run_arena.is_none() {
            self.run_arena = Some(Arc::clone(&batch.arena));
            if let (Some(ctx), Some((ledger, memory_limit))) = (&self.ctx, &self.budget) {
                // The records this run's arena cap holds at the previous run's mean
                // size. Every order's seal moved the last run's refs into its chunk,
                // so this sizes the run's fresh buffer.
                let prior = ledger.snapshot();
                let mean = prior.mean_record_bytes.unwrap_or(RUN0_MEAN_RECORD_BYTES).max(1);
                let records = run_cap_for(*memory_limit, prior) as u64 / mean;
                let acc = self.acc.get_or_insert_with(|| X::accum_for(ctx));
                X::reserve(acc, usize::try_from(records).unwrap_or(usize::MAX));
            }
        }
        if let Some(refs) = batch.refs {
            let ctx = self.ctx.as_ref().expect("a batch with refs has a context");
            X::append(self.acc.get_or_insert_with(|| X::accum_for(ctx)), refs);
        }
        self.total_records += n as u64;
        self.census.batches += 1;
        self.census.records += n as u64;
        self.run_batches.push((n as u64, batch.keyed_at));
        if let Some(end) = batch.run_end {
            self.seal_run(end)?;
        }
        Ok(())
    }

    /// Sort the run and stage its events.
    fn seal_run(&mut self, end: RunEnd) -> io::Result<()> {
        // The run's last batch has arrived: the tail wait ends here, before the
        // run sort, which the census must not count as extraction lag.
        let last_batch_at = Instant::now();
        let arena = self.run_arena.take().expect("a run with a last batch has an arena");
        if let Some((ledger, _)) = &self.budget {
            let records: u64 = self.run_batches.iter().map(|(n, _)| n).sum();
            ledger.publish_seal(end.run_bytes, records);
        }
        let chunk = match self.ctx.as_ref() {
            Some(ctx) => {
                let acc = self.acc.get_or_insert_with(|| X::accum_for(ctx));
                let parallel = X::sorts_in_parallel(acc, self.sort_threads);
                // With a phase-1 cap, a parallel sort holds the whole cap: new pool
                // admissions are refused while it waits, and it sorts at the full
                // phase width.
                let _whole = reserve_for_run_sort(self.cap.as_deref(), parallel)?;
                X::seal(acc, ctx, arena, self.sort_threads, &self.pool)
            }
            None => X::empty_chunk(),
        };
        for (n, keyed_at) in self.run_batches.drain(..) {
            if keyed_at < end.scan_end {
                self.census.overlapped_records += n;
            } else {
                self.census.tail_records += n;
            }
        }
        self.census.tail_wait_ns +=
            u64::try_from(last_batch_at.saturating_duration_since(end.scan_end).as_nanos())
                .unwrap_or(0);
        self.census.runs += 1;
        let seq = self.next_seq;
        self.next_seq += 1;
        if end.seals_to_spill {
            self.pending.push_back(SortChunkEvent::Spill {
                seq,
                chunk,
                records_ingested_so_far: end.records_ingested_so_far,
            });
            self.spilled_run_count += 1;
        } else {
            self.pending.push_back(SortChunkEvent::Residual {
                chunk,
                records_ingested_so_far: end.records_ingested_so_far,
            });
            self.pending.push_back(SortChunkEvent::AllAnnounced {
                slot_count: self.spilled_run_count,
                memory_chunk_count: 1,
                total_records: end.records_ingested_so_far,
            });
            self.finalized = true;
        }
        Ok(())
    }

    fn log_census(&mut self) {
        if self.census_logged || self.census.runs == 0 {
            return;
        }
        self.census_logged = true;
        #[cfg(test)]
        if let Some(probe) = &self.probe {
            *probe.census.lock() = Some(self.census);
        }
        for line in self.census.lines(self.sort_threads) {
            if self.sort_stats {
                info!("{line}");
            } else {
                debug!("{line}");
            }
        }
    }

    fn emit_pending(&mut self, ctx: &mut StepCtx<'_, Self>) -> StepOutcome {
        let Some(event) = self.pending.pop_front() else { return StepOutcome::NoProgress };
        if let Err(unpushed) = ctx.outputs.push(event) {
            self.held.put(unpushed);
        }
        StepOutcome::Progress
    }
}

impl<X: ArenaKeyExtractor> Step for SortRun<X> {
    type Input = RefBatch<X>;
    type Outputs = Single<SortChunkEvent>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: X::STEP_NAMES[2],
            kind: StepKind::Detached,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn detached_group(&self) -> DetachedGroup {
        DetachedGroup::Shared(crate::sort::SORT_COORD_GROUP)
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] = &[CounterSpec::new("records", "records")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld) {
            return Ok(StepOutcome::Contention);
        }
        if !self.pending.is_empty() {
            return Ok(self.emit_pending(ctx));
        }
        if self.finalized {
            self.log_census();
            return Ok(StepOutcome::Finished);
        }
        if let Some(batch) = ctx.input.pop() {
            ctx.counters.add(0, batch.len() as u64);
            self.accept(batch)?;
            if !self.pending.is_empty() {
                return Ok(self.emit_pending(ctx));
            }
            return Ok(StepOutcome::Progress);
        }
        if !ctx.input.is_drained() {
            return Ok(StepOutcome::NoProgress);
        }
        // Drained without a residual: announce the spills (unreachable under
        // the boundary step's always-residual final run, kept fail-safe).
        if self.spilled_run_count > 0 {
            self.pending.push_back(SortChunkEvent::AllAnnounced {
                slot_count: self.spilled_run_count,
                memory_chunk_count: 0,
                total_records: self.total_records,
            });
            self.finalized = true;
            return Ok(self.emit_pending(ctx));
        }
        self.log_census();
        Ok(StepOutcome::Finished)
    }

    fn phase_cap(&self) -> Option<&PhaseCap> {
        self.cap.as_deref()
    }
}

// ============================================================================
// Template-coordinate
// ============================================================================

/// Template-coordinate keys: the lane variant is chosen from the sort's first
/// record; refs are the variant's narrow lane.
pub struct TemplateKeys;

impl ArenaKeyExtractor for TemplateKeys {
    const STEP_NAMES: [&'static str; 3] =
        ["FindTemplateBoundaries", "ExtractTemplateKeys", "SortTemplateRun"];
    type Seed = TemplateKeySeed;
    type Ctx = Arc<TemplateKeyContext>;
    type Refs = ExtractedRefs;
    type Accum = TemplateArenaRefs;
    type Scratch = ();

    fn provision(seed: &Self::Seed, first_body: &[u8]) -> Self::Ctx {
        seed.provision(first_body)
    }

    fn known_ref_width(seed: &Self::Seed) -> Option<usize> {
        seed.known_variant().map(TemplateKeyVariant::ref_width)
    }

    fn ref_width(ctx: &Self::Ctx) -> usize {
        ctx.variant.ref_width()
    }

    fn extract_batch(
        ctx: &Self::Ctx,
        _scratch: &mut (),
        arena: &PooledSegmentedBuf,
        byte_start: u64,
        byte_end: u64,
        n: u32,
        first_index: u64,
    ) -> io::Result<(Self::Refs, Option<KeyViolation>)> {
        ctx.extract_batch_erased(arena, byte_start, byte_end, n, first_index)
    }

    fn refs_len(refs: &Self::Refs) -> usize {
        refs.len()
    }

    fn refs_byte_len(refs: &Self::Refs) -> usize {
        refs.byte_len()
    }

    fn accum_for(ctx: &Self::Ctx) -> Self::Accum {
        TemplateArenaRefs::for_variant(ctx.variant)
    }

    fn reserve(acc: &mut Self::Accum, n: usize) {
        acc.reserve(n);
    }

    fn append(acc: &mut Self::Accum, refs: Self::Refs) {
        acc.append(refs);
    }

    fn sorts_in_parallel(acc: &Self::Accum, sort_threads: usize) -> bool {
        radix_sorts_in_parallel(acc.len(), sort_threads)
    }

    fn seal(
        acc: &mut Self::Accum,
        _ctx: &Self::Ctx,
        arena: Arc<PooledSegmentedBuf>,
        sort_threads: usize,
        pool: &BoundedSortPool,
    ) -> MemoryChunkErased {
        MemoryChunkErased::TemplateCoordinate(seal_template_refs(
            acc,
            arena,
            sort_threads,
            pool.bounded(sort_threads),
        ))
    }

    fn empty_chunk() -> MemoryChunkErased {
        MemoryChunkErased::TemplateCoordinate(TemplateMemChunk::K40(InMemoryChunk::default()))
    }

    fn pool_name() -> &'static str {
        "tmpl-sort"
    }
}

/// The template-coordinate boundary step.
pub type FindTemplateBoundaries = FindRecordBoundaries<TemplateKeys>;
/// The template-coordinate extraction step.
pub type ExtractTemplateKeys = ExtractKeys<TemplateKeys>;
/// The template-coordinate run sort.
pub type SortTemplateRun = SortRun<TemplateKeys>;
/// A template-coordinate key batch.
pub type TemplateKeyBatch = KeyBatch<TemplateKeys>;
/// A template-coordinate ref batch.
pub type TemplateRefBatch = RefBatch<TemplateKeys>;

// ============================================================================
// Shared frame walk
// ============================================================================

/// An `InvalidData` error for a key batch whose frames do not match the scan's
/// cut.
fn malformed_batch(first_index: u64, what: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("key batch at record {first_index}: {what}"))
}

/// Re-walk `[byte_start, byte_end)` of `arena` as exactly `n` `[u32 len][body]`
/// frames, calling `f(body, body_offset, len)` per record in order, with an L1
/// prefetch ahead of the walk.
///
/// # Errors
///
/// `InvalidData` if the span holds fewer or more than `n` whole frames, or a
/// frame is shorter than a BAM record.
#[inline]
fn for_each_frame(
    arena: &PooledSegmentedBuf,
    byte_start: u64,
    byte_end: u64,
    n: u32,
    first_index: u64,
    mut f: impl FnMut(&[u8], u64, u32),
) -> io::Result<()> {
    let start = usize::try_from(byte_start)
        .map_err(|_| malformed_batch(first_index, "start overflows usize"))?;
    let end = usize::try_from(byte_end)
        .map_err(|_| malformed_batch(first_index, "end overflows usize"))?;
    if end < start {
        return Err(malformed_batch(first_index, "end precedes start"));
    }
    let span = arena.slice(start, end - start);
    let mut cur = 0usize;
    for _ in 0..n {
        if cur + 4 > span.len() {
            return Err(malformed_batch(first_index, "fewer frames than n_records"));
        }
        let bs = u32::from_le_bytes([span[cur], span[cur + 1], span[cur + 2], span[cur + 3]]);
        if (bs as usize) < MIN_BAM_RECORD_LEN {
            return Err(malformed_batch(first_index, "frame shorter than a BAM record"));
        }
        let body_start = cur + 4;
        let body_end = body_start + bs as usize;
        if body_end > span.len() {
            return Err(malformed_batch(first_index, "frame runs past byte_end"));
        }
        if let Some(ahead) = span.get(cur + KEY_PREFETCH_DISTANCE) {
            prefetch_read_l1(ahead);
        }
        f(&span[body_start..body_end], byte_start + body_start as u64, bs);
        cur = body_end;
    }
    if cur != span.len() {
        return Err(malformed_batch(first_index, "bytes remain after the last frame"));
    }
    Ok(())
}

// ============================================================================
// Coordinate
// ============================================================================

/// The coordinate key's per-sort context: the header's reference count, which
/// places unmapped records last.
#[derive(Debug, Clone, Copy)]
pub struct CoordinateCtx {
    /// BAM header reference-sequence count.
    pub n_ref: u32,
}

/// Coordinate keys: one 24-byte `RecordRef` (a packed `u64` key, the body's
/// arena offset and length) per record, radix-sorted at seal.
///
/// The radix spawns exactly `sort_threads` threads itself, so the bounded pool
/// the trait hands `seal` is unused for coordinate order; the trait stays
/// uniform across orders.
pub struct CoordinateKeys;

impl ArenaKeyExtractor for CoordinateKeys {
    const STEP_NAMES: [&'static str; 3] =
        ["FindCoordinateBoundaries", "ExtractCoordinateKeys", "SortCoordinateRun"];
    /// The header's reference count.
    type Seed = u32;
    type Ctx = CoordinateCtx;
    type Refs = Vec<RecordRef>;
    type Accum = Vec<RecordRef>;
    type Scratch = ();

    fn provision(seed: &Self::Seed, _first_body: &[u8]) -> Self::Ctx {
        CoordinateCtx { n_ref: *seed }
    }

    fn known_ref_width(_seed: &Self::Seed) -> Option<usize> {
        Some(size_of::<RecordRef>())
    }

    fn ref_width(_ctx: &Self::Ctx) -> usize {
        size_of::<RecordRef>()
    }

    fn extract_batch(
        ctx: &Self::Ctx,
        _scratch: &mut (),
        arena: &PooledSegmentedBuf,
        byte_start: u64,
        byte_end: u64,
        n: u32,
        first_index: u64,
    ) -> io::Result<(Self::Refs, Option<KeyViolation>)> {
        let mut refs: Vec<RecordRef> = Vec::with_capacity(n as usize);
        for_each_frame(arena, byte_start, byte_end, n, first_index, |body, offset, len| {
            refs.push(RecordRef::new(extract_coordinate_key_inline(body, ctx.n_ref), offset, len));
        })?;
        Ok((refs, None))
    }

    fn refs_len(refs: &Self::Refs) -> usize {
        refs.len()
    }

    fn refs_byte_len(refs: &Self::Refs) -> usize {
        refs.len() * size_of::<RecordRef>()
    }

    fn accum_for(_ctx: &Self::Ctx) -> Self::Accum {
        Vec::new()
    }

    fn reserve(acc: &mut Self::Accum, n: usize) {
        acc.reserve(n);
    }

    fn append(acc: &mut Self::Accum, mut refs: Self::Refs) {
        acc.append(&mut refs);
    }

    fn sorts_in_parallel(acc: &Self::Accum, sort_threads: usize) -> bool {
        radix_sorts_in_parallel(acc.len(), sort_threads)
    }

    fn seal(
        acc: &mut Self::Accum,
        _ctx: &Self::Ctx,
        arena: Arc<PooledSegmentedBuf>,
        sort_threads: usize,
        _pool: &BoundedSortPool,
    ) -> MemoryChunkErased {
        // The refs move into the chunk (converted in place), leaving `acc`
        // empty with no capacity.
        MemoryChunkErased::Coordinate(coordinate_chunk_from_refs(
            arena,
            std::mem::take(acc),
            sort_threads,
        ))
    }

    fn empty_chunk() -> MemoryChunkErased {
        MemoryChunkErased::Coordinate(InMemoryChunk::default())
    }

    fn pool_name() -> &'static str {
        "coord-sort"
    }
}

/// The coordinate boundary step.
pub type FindCoordinateBoundaries = FindRecordBoundaries<CoordinateKeys>;
/// The coordinate extraction step.
pub type ExtractCoordinateKeys = ExtractKeys<CoordinateKeys>;
/// The coordinate run sort.
pub type SortCoordinateRun = SortRun<CoordinateKeys>;

// ============================================================================
// Queryname
// ============================================================================

/// A queryname comparator that knows the erased chunk its runs seal into.
pub trait QnChunkOrder: QnComparator + Send + Sync + 'static {
    /// Erase a sorted queryname chunk into this order's variant.
    fn wrap(chunk: InMemoryChunk<QnKey>) -> MemoryChunkErased;
}

impl QnChunkOrder for QnLex {
    fn wrap(chunk: InMemoryChunk<QnKey>) -> MemoryChunkErased {
        MemoryChunkErased::QuerynameLex(chunk)
    }
}

impl QnChunkOrder for QnNatural {
    fn wrap(chunk: InMemoryChunk<QnKey>) -> MemoryChunkErased {
        MemoryChunkErased::QuerynameNatural(chunk)
    }
}

/// Queryname keys: 40-byte refs over the arena (a 24-byte name window behind
/// a learned prefix class; see `fgumi_sort::queryname_ref`). Each extraction
/// clone carries its own class guess (`QnRefBuilder`) across the batches it
/// extracts; every batch is classified, and re-windowed when it learns a new
/// class, before it leaves the extractor. The order is `C`'s.
pub struct QuerynameKeys<C>(std::marker::PhantomData<C>);

impl<C: QnChunkOrder> ArenaKeyExtractor for QuerynameKeys<C> {
    const STEP_NAMES: [&'static str; 3] =
        ["FindQuerynameBoundaries", "ExtractQuerynameKeys", "SortQuerynameRun"];
    /// The sort's class table and counters, shared with its merge.
    type Seed = QnSortState;
    type Ctx = QnSortState;
    type Refs = Vec<QnRef>;
    type Accum = Vec<QnRef>;
    type Scratch = QnRefBuilder;

    fn provision(seed: &Self::Seed, _first_body: &[u8]) -> Self::Ctx {
        seed.clone()
    }

    fn known_ref_width(_seed: &Self::Seed) -> Option<usize> {
        Some(size_of::<QnRef>())
    }

    fn ref_width(_ctx: &Self::Ctx) -> usize {
        size_of::<QnRef>()
    }

    fn extract_batch(
        ctx: &Self::Ctx,
        builder: &mut QnRefBuilder,
        arena: &PooledSegmentedBuf,
        byte_start: u64,
        byte_end: u64,
        n: u32,
        first_index: u64,
    ) -> io::Result<(Self::Refs, Option<KeyViolation>)> {
        let table = &ctx.table;
        let mut lcp = QnBatchLcp::default();
        let mut refs: Vec<QnRef> = Vec::with_capacity(n as usize);
        for_each_frame(arena, byte_start, byte_end, n, first_index, |body, offset, len| {
            refs.push(builder.build(table, &mut lcp, body, offset, len));
        })?;
        if !refs.is_empty() {
            let reclassified = builder.reclassify(table, arena, &mut refs, &lcp);
            let class0 =
                if reclassified { 0 } else { refs.iter().filter(|r| r.0.class == 0).count() };
            ctx.stats.note_batch(refs.len(), class0, reclassified);
        }
        Ok((refs, None))
    }

    fn refs_len(refs: &Self::Refs) -> usize {
        refs.len()
    }

    fn refs_byte_len(refs: &Self::Refs) -> usize {
        refs.len() * size_of::<QnRef>()
    }

    fn accum_for(_ctx: &Self::Ctx) -> Self::Accum {
        Vec::new()
    }

    fn reserve(acc: &mut Self::Accum, n: usize) {
        acc.reserve(n);
    }

    fn append(acc: &mut Self::Accum, mut refs: Self::Refs) {
        acc.append(&mut refs);
    }

    /// The comparator sort always runs on the bounded pool, so it is parallel
    /// whenever that pool has more than one thread.
    fn sorts_in_parallel(_acc: &Self::Accum, sort_threads: usize) -> bool {
        sort_threads > 1
    }

    /// The run's refs move into the chunk; nothing is copied.
    fn seal(
        acc: &mut Self::Accum,
        ctx: &Self::Ctx,
        arena: Arc<PooledSegmentedBuf>,
        sort_threads: usize,
        pool: &BoundedSortPool,
    ) -> MemoryChunkErased {
        let refs = std::mem::take(acc);
        pool.bounded(sort_threads).install(move || {
            C::wrap(queryname_chunk_from_qn_refs::<C>(&ctx.table, arena, refs, ctx.compare_stats()))
        })
    }

    fn empty_chunk() -> MemoryChunkErased {
        C::wrap(InMemoryChunk::default())
    }

    fn pool_name() -> &'static str {
        "qname-sort"
    }
}

/// The queryname boundary step, in comparator `C`'s order.
pub type FindQuerynameBoundaries<C> = FindRecordBoundaries<QuerynameKeys<C>>;
/// The queryname extraction step.
pub type ExtractQuerynameKeys<C> = ExtractKeys<QuerynameKeys<C>>;
/// The queryname run sort.
pub type SortQuerynameRun<C> = SortRun<QuerynameKeys<C>>;

/// Test probe shared by a split front's steps: can hold the scan after the
/// sort's first batch until that batch is extracted (proving extraction runs
/// while the scan goes on), and records the threads extraction ran on and the
/// run sort's final census.
#[cfg(test)]
#[derive(Default)]
pub(crate) struct SplitProbe {
    /// Hold the scan after batch 0 until it is extracted.
    pub(crate) hold_scan: bool,
    extracted_first: std::sync::atomic::AtomicBool,
    /// Names of the threads extraction ran on.
    pub(crate) extract_threads: parking_lot::Mutex<std::collections::BTreeSet<String>>,
    /// The run sort's census at its end.
    pub(crate) census: parking_lot::Mutex<Option<KeyOverlapCensus>>,
}

#[cfg(test)]
impl SplitProbe {
    /// A probe that holds the scan after batch 0 until it is extracted.
    pub(crate) fn holding_scan() -> Self {
        Self { hold_scan: true, ..Self::default() }
    }

    fn extracted(&self, ordinal: u64) {
        if ordinal == 0 {
            self.extracted_first.store(true, std::sync::atomic::Ordering::Release);
        }
        let name = std::thread::current().name().unwrap_or("").to_string();
        self.extract_threads.lock().insert(name);
    }

    fn extracted_first(&self) -> bool {
        self.extracted_first.load(std::sync::atomic::Ordering::Acquire)
    }
}

#[cfg(test)]
mod tests;
