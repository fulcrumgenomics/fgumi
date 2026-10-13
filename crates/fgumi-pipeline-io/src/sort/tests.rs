//! Tests for the sort chains.
//!
//! Two phase-1 fronts feed the same spill and merge steps: the record-input
//! front (`SortBuffer`) and the block-input arena front (`ReadBlocks` →
//! `InflateToArena` → the order's arena steps, see [`production_arena_front`]).
//! Both spill through `SpillGather` → `SpillBlockCompress` → `SpillWrite` and
//! merge through `SpillReadPlanner` → `PreadSpillSlices` →
//! `SortSpillDecompress` → `SortMerge`.
//!
//! The output contract every parity test checks is [`reference_sort`]: the
//! in-memory stable sort of the input records by the order's key. Equal keys
//! keep input order — at every thread count, across spill boundaries and on
//! both fronts — so every order, queryname included, is checked for byte
//! identity.

use rand::RngExt as _;
use std::io;
use std::sync::Arc;

use anyhow::Result;
use fgumi_raw_bam::RawRecord;
use fgumi_raw_bam::testutil::make_bam_bytes;
use fgumi_sort::{KeyTypesSpec, QuerynameComparator, SortOrder, SpillCodec};
use noodles::sam::Header;
use parking_lot::Mutex;
use rstest::rstest;

use super::*;
use crate::sort::protocol::SortChunkEvent;
use crate::sort::testing::{
    N_TEST_REFS, RecordBytes, assert_input_order_within_equal_keys, chance,
    equal_key_group_spans_a_boundary, fixture_rng, max_equal_key_group, reference_sort,
    synthetic_header, tie_heavy_records, xi_of,
};
use crate::types::RecordBatch;
use fgumi_pipeline_core::{
    HeldRetry, Unpushed,
    builder::{Chain, Pipeline, PipelineConfig},
    held::HeldSlot,
    outputs::{OrderedBytesSingle, Single},
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

// ── In-memory source / sink test steps ──────────────────────────────────────

/// `Exclusive` source that drains a `Vec<RecordBatch>` one batch per `try_run` call.
struct VecSource {
    batches: Vec<RecordBatch>,
    held: HeldSlot<Unpushed<RecordBatch>>,
    output_byte_limit: u64,
}

impl VecSource {
    fn new(mut batches: Vec<RecordBatch>, output_byte_limit: u64) -> Self {
        batches.reverse();
        Self { batches, held: HeldSlot::new(), output_byte_limit }
    }
}

impl Step for VecSource {
    type Input = ();
    type Outputs = OrderedBytesSingle<RecordBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "VecSource",
            kind: StepKind::Exclusive,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
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
        let Some(batch) = self.batches.pop() else {
            return Ok(StepOutcome::Finished);
        };
        match ctx.outputs.push(batch) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }
}

/// Sink that appends every received batch into a shared `Vec<RecordBatch>`.
struct VecSink {
    received: Arc<Mutex<Vec<RecordBatch>>>,
    kind: StepKind,
    /// Records received so far, for a test feeder that gates on sink progress.
    progress: Option<Arc<std::sync::atomic::AtomicUsize>>,
}

impl Step for VecSink {
    type Input = RecordBatch;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "VecSink",
            kind: self.kind,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(batch) => {
                let n = batch.len();
                self.received.lock().push(batch);
                if let Some(p) = &self.progress {
                    p.fetch_add(n, std::sync::atomic::Ordering::Release);
                }
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// Pass-through between a phase-1 front and `SpillGather` that records every
/// `Spill`'s `records_ingested_so_far`: the run boundaries, as input record
/// counts, in run order.
struct SpillBoundaryTap {
    boundaries: Arc<Mutex<Vec<u64>>>,
    held: HeldSlot<Unpushed<SortChunkEvent>>,
    output_byte_limit: u64,
}

impl Step for SpillBoundaryTap {
    type Input = SortChunkEvent;
    type Outputs = Single<SortChunkEvent>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SpillBoundaryTap",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld) {
            return Ok(StepOutcome::Contention);
        }
        let Some(event) = ctx.input.pop() else {
            return Ok(if ctx.input.is_drained() {
                StepOutcome::Finished
            } else {
                StepOutcome::NoProgress
            });
        };
        if let SortChunkEvent::Spill { records_ingested_so_far, .. } = &event {
            self.boundaries.lock().push(*records_ingested_so_far);
        }
        if let Err(unpushed) = ctx.outputs.push(event) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

/// Spill codec and compression level a test chain writes its spill files
/// with. They shape only the temporary files, never the merged records.
type SpillFormat = (SpillCodec, u32);

/// The spill format the parity tests use unless they vary it.
const ZSTD_1: SpillFormat = (SpillCodec::Zstd, 1);

/// A [`SortBufferConfig`] with no cell tag and auto-detected key types.
fn buffer_config(
    sort_order: SortOrder,
    memory_limit: usize,
    sort_threads: usize,
) -> SortBufferConfig {
    SortBufferConfig {
        sort_order,
        memory_limit,
        sort_threads,
        cell_tag: None,
        key_types: KeyTypesSpec::Auto,
    }
}

// ── Synthetic-record helpers ────────────────────────────────────────────────

fn synthesize_records(n: usize, seed: u64) -> (Header, Vec<RawRecord>) {
    synthesize_sized_records(n, seed, 0)
}

fn pack_batches(records: &[RawRecord], batch_size: usize) -> Vec<RecordBatch> {
    use crate::types::RecordBatchBuilder;
    records
        .chunks(batch_size)
        .enumerate()
        .map(|(i, chunk)| {
            let total: usize = chunk.iter().map(RawRecord::len).sum();
            let mut b = RecordBatchBuilder::with_capacity(i as u64, total, chunk.len());
            for r in chunk {
                b.push_record_bytes(r.as_ref());
            }
            b.build()
        })
        .collect()
}

// ── Chain harness ───────────────────────────────────────────────────────────

/// What one sort chain produced.
#[derive(Debug, Default)]
struct ChainRun {
    /// The merged record bytes, in output order.
    records: RecordBytes,
    /// Spill files the merge saw (`SortStats::runs_written`): the runs the
    /// `SpillGather` run-former formed.
    runs: usize,
    /// Spill runs handed to the final merge (`SpillRunStats::merge_sources`).
    merge_sources: u64,
    /// `records_ingested_so_far` of every `Spill` the front emitted, in order.
    spill_boundaries: Vec<u64>,
    /// Queryname classes the sort's shared class table holds after the run
    /// (the front learned them; the merge rebuilt its slot keys against them).
    qn_classes: usize,
}

/// The merge's output end: `SortMerge<RecordBatchOutput>` into a [`VecSink`],
/// or the terminal `SortMerge<BlockOutput>` into a [`BlockSink`].
#[derive(Clone, Copy, Debug)]
enum MergeSink {
    Batches(StepKind),
    Blocks(StepKind),
}

/// The merge-supply knobs the parity tests vary; production resolves them
/// from `--max-memory`, the slot count and the read-stream policy.
#[derive(Clone, Copy, Debug)]
struct SupplyTuning {
    /// `(hot_fill, cold_fill)` replacing the resolved budget's fill sizes.
    fills: Option<(u64, u64)>,
    /// Slices per fill (`--read-streams`).
    read_streams: usize,
    /// `(max_inflight_slices, cold_inflight_bytes)`, forced small.
    inflight_caps: Option<(u32, u64)>,
}

impl Default for SupplyTuning {
    fn default() -> Self {
        Self { fills: None, read_streams: 1, inflight_caps: None }
    }
}

/// The merge supply as `add_sort` builds it — `SpillReadPlanner` (fan-out:
/// read requests on branch 0, phase events on branch 1) → `PreadSpillSlices`
/// → `SortSpillDecompress` — sharing one `MergeDemand` and ledger with the
/// merge the caller wires to the planner's branch 1.
fn spill_supply(
    output_byte_limit: u64,
    threads: usize,
    tuning: SupplyTuning,
) -> (SpillReadPlanner, crate::pread::PreadSlices, SortSpillDecompress, Arc<fgumi_sort::MergeDemand>)
{
    let supply = SpillSupply::new();
    let demand = Arc::clone(&supply.demand);
    let slices = fgumi_bam_io::pread::SliceBufferPool::new(2 * threads + 8);
    let mut planner = SpillReadPlanner::new(
        768 << 20,
        threads,
        threads,
        output_byte_limit,
        &supply,
        Arc::clone(&slices),
    )
    .with_read_streams(fgumi_bam_io::pread::ReadStreamsPolicy::fixed(tuning.read_streams));
    if let Some((hot, cold)) = tuning.fills {
        planner = planner.with_fills_for_test(hot, cold);
    }
    if let Some((slices, cold_bytes)) = tuning.inflight_caps {
        planner = planner.with_inflight_caps_for_test(slices, cold_bytes);
    }
    let pread = crate::pread::PreadSlices::spill(slices, output_byte_limit);
    let decompress = SortSpillDecompress::new(output_byte_limit, &supply);
    (planner, pread, decompress, demand)
}

/// Optional phase caps for the pool steps of a test chain: phase 1 on
/// `InflateToArena`, the front's capped steps and `SpillBlockCompress`;
/// phase 2 on `PreadSpillSlices` and `SortSpillDecompress`.
#[derive(Clone, Default)]
struct ChainCaps {
    phase1: Option<Arc<fgumi_pipeline_core::PhaseCap>>,
    phase2: Option<Arc<fgumi_pipeline_core::PhaseCap>>,
    /// Records per key batch on the split front (`None`: the production size).
    batch_records: Option<usize>,
    /// Probe on the split front's steps.
    probe: Option<Arc<crate::sort::arena_split::SplitProbe>>,
    /// Make the block source a `Serial` reader-affinity step (worker 0), as
    /// production's input reader is.
    reader_source: bool,
}

/// What a test chain's taps and sinks collect while it runs; read back by
/// [`ChainTaps::collect`] once the pipeline has finished.
struct ChainTaps {
    class_table: Arc<fgumi_sort::QnClassTable>,
    boundaries: Arc<Mutex<Vec<u64>>>,
    spill_stats: Arc<SpillRunStats>,
    merge_stats: Arc<Mutex<Option<fgumi_sort::SortStats>>>,
    batches: Arc<Mutex<Vec<RecordBatch>>>,
    blocks: Arc<Mutex<Vec<crate::types::DecompressedBlock>>>,
}

impl ChainTaps {
    fn collect(self) -> ChainRun {
        let mut records = Vec::new();
        for batch in std::mem::take(&mut *self.batches.lock()) {
            records.extend(batch.iter_record_bytes().map(<[u8]>::to_vec));
        }
        let mut blocks = std::mem::take(&mut *self.blocks.lock());
        // `Detached` collapses ordering to None, so blocks arrive in merge-emit
        // order on the single sink; sort by serial defensively in case a future
        // change makes the sink parallel.
        blocks.sort_by_key(|b| b.batch_serial);
        for block in &blocks {
            records.extend(unframe_block_records(&block.bytes));
        }
        let runs = self.merge_stats.lock().take().map_or(0, |s| s.runs_written);
        let spill_boundaries = std::mem::take(&mut *self.boundaries.lock());
        ChainRun {
            records,
            runs,
            merge_sources: self.spill_stats.merge_sources(),
            spill_boundaries,
            qn_classes: self.class_table.len(),
        }
    }
}

/// Everything below a phase-1 front, shared by both fronts: the boundary tap,
/// the spill split (`SpillGather` → `SpillBlockCompress` → `SpillWrite`), the
/// merge supply and the merge into `sink`. The merge holds `class_table`, the
/// one the front's queryname state windowed every key against, as production
/// wires it.
#[allow(clippy::too_many_arguments)]
fn attach_spill_and_merge(
    front: Chain<'_, Single<SortChunkEvent>>,
    class_table: &Arc<fgumi_sort::QnClassTable>,
    sort_order: SortOrder,
    output_byte_limit: u64,
    threads: usize,
    sink: MergeSink,
    supply: SupplyTuning,
    (codec, compression): SpillFormat,
    caps: &ChainCaps,
) -> Result<ChainTaps> {
    use fgumi_sort::TmpDirAllocator;

    let dir = tempfile::TempDir::new()?;
    let alloc =
        TmpDirAllocator::with_probe(vec![dir.path().to_path_buf()], Box::new(|_| Ok(u64::MAX)), 0)?;
    let temp_dirs = Arc::new(vec![dir]);
    let boundaries = Arc::new(Mutex::new(Vec::new()));
    let spill_stats = Arc::new(SpillRunStats::default());
    let merge_stats: Arc<Mutex<Option<fgumi_sort::SortStats>>> = Arc::new(Mutex::new(None));
    let batches: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let blocks: Arc<Mutex<Vec<crate::types::DecompressedBlock>>> = Arc::new(Mutex::new(Vec::new()));

    let tap = SpillBoundaryTap {
        boundaries: Arc::clone(&boundaries),
        held: HeldSlot::new(),
        output_byte_limit,
    };
    let (planner, pread, decompress, demand) = spill_supply(output_byte_limit, threads, supply);
    let fan = front
        .chain(tap)
        .chain(SpillGather::new(output_byte_limit))
        .chain(
            SpillBlockCompress::new(codec, compression, output_byte_limit)
                .with_step_cap(caps.phase1.clone().map(fgumi_pipeline_core::StepCap::from)),
        )
        .chain(
            SpillWrite::new(Arc::new(Mutex::new(alloc)), codec, output_byte_limit, temp_dirs)
                .with_stats(Arc::clone(&spill_stats)),
        )
        .chain(planner)
        .into_multi();
    fan.b2.chain(crate::sort::NoConsolidation).into_sink_marker();
    fan.b0
        .chain(pread.with_step_cap(caps.phase2.clone().map(fgumi_pipeline_core::StepCap::from)))
        .chain(
            decompress.with_step_cap(caps.phase2.clone().map(fgumi_pipeline_core::StepCap::from)),
        )
        .into_sink_marker();
    match sink {
        MergeSink::Batches(kind) => fan
            .b1
            .chain(
                SortMerge::<RecordBatchOutput>::with_target_batch_count(
                    sort_order,
                    output_byte_limit,
                    256,
                )
                .with_stats_slot(Arc::clone(&merge_stats))
                .with_merge_demand(demand)
                .with_class_table(Arc::clone(class_table)),
            )
            .chain(VecSink { received: Arc::clone(&batches), kind, progress: None })
            .into_sink_marker(),
        MergeSink::Blocks(kind) => fan
            .b1
            .chain(
                SortMerge::<BlockOutput>::with_target_batch_count(
                    sort_order,
                    output_byte_limit,
                    256,
                )
                .with_stats_slot(Arc::clone(&merge_stats))
                .with_merge_demand(demand)
                .with_class_table(Arc::clone(class_table)),
            )
            .chain(BlockSink { received: Arc::clone(&blocks), kind })
            .into_sink_marker(),
    }
    Ok(ChainTaps {
        class_table: Arc::clone(class_table),
        boundaries,
        spill_stats,
        merge_stats,
        batches,
        blocks,
    })
}

/// Drive the record-input front: `VecSource` → `SortBuffer` → the spill and
/// merge steps of [`attach_spill_and_merge`].
#[allow(clippy::too_many_arguments)]
fn drive_buffer_chain(
    config: SortBufferConfig,
    header: &Header,
    batches: Vec<RecordBatch>,
    output_byte_limit: u64,
    threads: usize,
    sink: MergeSink,
    supply: SupplyTuning,
    spill: SpillFormat,
) -> Result<ChainRun> {
    let qn_keys = fgumi_sort::QnSortState::new(false);
    let builder = Pipeline::builder();
    let front = builder.chain(VecSource::new(batches, output_byte_limit)).chain(
        SortBuffer::new(config, header, output_byte_limit)?.with_queryname_keys(qn_keys.clone()),
    );
    let taps = attach_spill_and_merge(
        front,
        &qn_keys.table,
        config.sort_order,
        output_byte_limit,
        threads,
        sink,
        supply,
        spill,
        &ChainCaps::default(),
    )?;
    builder.build()?.run(PipelineConfig { threads, ..Default::default() })?;
    Ok(taps.collect())
}

/// Drive the production sort chain (`VecSource` → `SortBuffer` →
/// `SpillGather` → `SpillBlockCompress` → `SpillWrite` → merge supply →
/// `SortMerge` → `VecSink`) and collect the merged record bytes.
fn drive_sort_pipeline(
    config: SortBufferConfig,
    header: &Header,
    batches: Vec<RecordBatch>,
    output_byte_limit: u64,
    threads: usize,
    sink_kind: StepKind,
) -> Result<Vec<Vec<u8>>> {
    drive_sort_pipeline_tuned(
        config,
        header,
        batches,
        output_byte_limit,
        threads,
        sink_kind,
        SupplyTuning::default(),
        ZSTD_1,
    )
}

/// [`drive_sort_pipeline`] with the merge-supply knobs and the spill format
/// chosen by the caller.
#[allow(clippy::too_many_arguments)]
fn drive_sort_pipeline_tuned(
    config: SortBufferConfig,
    header: &Header,
    batches: Vec<RecordBatch>,
    output_byte_limit: u64,
    threads: usize,
    sink_kind: StepKind,
    supply: SupplyTuning,
    spill: SpillFormat,
) -> Result<Vec<Vec<u8>>> {
    Ok(drive_buffer_chain(
        config,
        header,
        batches,
        output_byte_limit,
        threads,
        MergeSink::Batches(sink_kind),
        supply,
        spill,
    )?
    .records)
}

/// Terminal-path sink: collects the [`DecompressedBlock`]s emitted by
/// `SortMerge<BlockOutput>` (the framed-bytes terminal output, lever 1).
struct BlockSink {
    received: Arc<Mutex<Vec<crate::types::DecompressedBlock>>>,
    kind: StepKind,
}

impl Step for BlockSink {
    type Input = crate::types::DecompressedBlock;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "BlockSink",
            kind: self.kind,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(block) => {
                self.received.lock().push(block);
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// Parse a `[u32 LE block_size][body]`-framed byte buffer (the terminal-sort
/// `BlockOutput` framing, identical to the former `SerializeRecordBatch`) into
/// its record bodies. Mirrors what `BgzfDecompress → FindBamBoundaries` does
/// downstream; kept local to the test so parity is checked against an
/// independent re-implementation of the layout.
fn unframe_block_records(bytes: &[u8]) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < bytes.len() {
        let len = u32::from_le_bytes(bytes[i..i + 4].try_into().unwrap()) as usize;
        i += 4;
        out.push(bytes[i..i + len].to_vec());
        i += len;
    }
    out
}

/// Drive the production sort chain with the **terminal** merge
/// (`SortMerge<BlockOutput>` → `BlockSink`) and recover the merged record
/// bodies by un-framing each `DecompressedBlock`. Used to prove the framed
/// terminal output carries byte-identical records to the `RecordBatchOutput`
/// path and the reference sort (lever 1 parity gate).
fn drive_sort_block_pipeline(
    config: SortBufferConfig,
    header: &Header,
    batches: Vec<RecordBatch>,
    output_byte_limit: u64,
    threads: usize,
    sink_kind: StepKind,
) -> Result<Vec<Vec<u8>>> {
    Ok(drive_buffer_chain(
        config,
        header,
        batches,
        output_byte_limit,
        threads,
        MergeSink::Blocks(sink_kind),
        SupplyTuning::default(),
        ZSTD_1,
    )?
    .records)
}

// ── Fronts ──────────────────────────────────────────────────────────────────

/// The phase-1 front a test drives.
#[derive(Clone, Copy, Debug)]
enum Front {
    /// `ReadBlocks` → `InflateToArena` → the order's arena steps
    /// ([`production_arena_front`]).
    Arena,
    /// `SortBuffer`, the record-input front.
    SortBuffer,
}

/// The arena steps of an order's block-input front.
#[derive(Clone, Copy, Debug)]
enum ArenaFront {
    /// The split front (`FindRecordBoundaries` → `ExtractKeys` → `SortRun`)
    /// over the order's `ArenaKeyExtractor`.
    Split,
}

/// The arena front each order ships with at this commit. When an order moves
/// to another front, this match is the one place that changes, so every
/// arena-front test gates the front that ships.
fn production_arena_front(_order: SortOrder) -> ArenaFront {
    ArenaFront::Split
}

/// The sort's run ledger, with the ref width `config`'s order knows up front.
fn run_ledger_for(config: SortBufferConfig, header: &Header) -> Arc<RunLedger> {
    let n_ref = u32::try_from(header.reference_sequences().len()).expect("n_ref fits u32");
    RunLedger::new(match config.sort_order {
        SortOrder::Coordinate => CoordinateKeys::known_ref_width(&n_ref),
        SortOrder::Queryname(_) => QuerynameKeys::<fgumi_sort::QnLex>::known_ref_width(
            &fgumi_sort::QnSortState::new(false),
        ),
        SortOrder::TemplateCoordinate => TemplateKeys::known_ref_width(
            &fgumi_sort::TemplateKeySeed::from_header(header, config.cell_tag, config.key_types),
        ),
    })
}

/// The split front's three steps for extractor `X`, with `caps`' test knobs.
fn split_steps<'b, X: ArenaKeyExtractor>(
    head: Chain<'b, OrderedBytesSingle<InflatedBlock>>,
    seed: X::Seed,
    sort_threads: usize,
    output_byte_limit: u64,
    (ledger, memory_limit): (&Arc<RunLedger>, usize),
    caps: &ChainCaps,
) -> Chain<'b, Single<SortChunkEvent>> {
    let mut boundaries =
        FindRecordBoundaries::<X>::new(seed, output_byte_limit).with_run_ledger(Arc::clone(ledger));
    let mut extract = ExtractKeys::<X>::new(output_byte_limit).with_phase_cap(caps.phase1.clone());
    let mut run = SortRun::<X>::new(sort_threads, output_byte_limit)
        .with_phase_cap(caps.phase1.clone())
        .with_run_ledger(Arc::clone(ledger), memory_limit);
    if let Some(n) = caps.batch_records {
        boundaries = boundaries.with_batch_records(n);
    }
    if let Some(probe) = &caps.probe {
        boundaries = boundaries.with_probe(Arc::clone(probe));
        extract = extract.with_probe(Arc::clone(probe));
        run = run.with_probe(Arc::clone(probe));
    }
    head.chain(boundaries).chain(extract).chain(run)
}

/// Append `front`'s steps for `config.sort_order` to an arena chain head
/// (`… → InflateToArena`), returning the front's chunk-event tail.
#[allow(clippy::too_many_arguments)]
fn arena_front<'b>(
    head: Chain<'b, OrderedBytesSingle<InflatedBlock>>,
    front: ArenaFront,
    config: SortBufferConfig,
    header: &Header,
    output_byte_limit: u64,
    caps: &ChainCaps,
    ledger: &Arc<RunLedger>,
    qn_keys: &fgumi_sort::QnSortState,
) -> Chain<'b, Single<SortChunkEvent>> {
    let n_ref = u32::try_from(header.reference_sequences().len()).expect("n_ref fits u32");
    let SortBufferConfig { sort_order, sort_threads, cell_tag, key_types, memory_limit } = config;
    let budget = (ledger, memory_limit);
    match (front, sort_order) {
        (ArenaFront::Split, SortOrder::TemplateCoordinate) => split_steps::<TemplateKeys>(
            head,
            fgumi_sort::TemplateKeySeed::from_header(header, cell_tag, key_types),
            sort_threads,
            output_byte_limit,
            budget,
            caps,
        ),
        (ArenaFront::Split, SortOrder::Queryname(QuerynameComparator::Lexicographic)) => {
            split_steps::<QuerynameKeys<fgumi_sort::QnLex>>(
                head,
                qn_keys.clone(),
                sort_threads,
                output_byte_limit,
                budget,
                caps,
            )
        }
        (ArenaFront::Split, SortOrder::Queryname(QuerynameComparator::Natural)) => {
            split_steps::<QuerynameKeys<fgumi_sort::QnNatural>>(
                head,
                qn_keys.clone(),
                sort_threads,
                output_byte_limit,
                budget,
                caps,
            )
        }
        (ArenaFront::Split, SortOrder::Coordinate) => split_steps::<CoordinateKeys>(
            head,
            n_ref,
            sort_threads,
            output_byte_limit,
            budget,
            caps,
        ),
    }
}

/// Drive the block-input front `BgzfBlockSource` → `ReadBlocks` →
/// `InflateToArena` → `front`'s steps for `config.sort_order`, then the spill
/// and merge steps of [`attach_spill_and_merge`]. `ReadBlocks` seals a run at
/// `config.memory_limit`; the run sort uses `config.sort_threads`.
#[allow(clippy::too_many_arguments)]
fn drive_arena_chain(
    front: ArenaFront,
    config: SortBufferConfig,
    header: &Header,
    blocks: Vec<crate::types::BgzfBlock>,
    output_byte_limit: u64,
    threads: usize,
    sink: MergeSink,
    caps: &ChainCaps,
) -> Result<ChainRun> {
    let ledger = run_ledger_for(config, header);
    let qn_keys = fgumi_sort::QnSortState::new(false);
    let builder = Pipeline::builder();
    let head = builder
        .chain(BgzfBlockSource::new(blocks, output_byte_limit).reader(caps.reader_source))
        .chain(ReadBlocks::new(config.memory_limit, output_byte_limit, Arc::clone(&ledger)))
        .chain(InflateToArena::new(output_byte_limit).with_phase_cap(caps.phase1.clone()));
    let taps = attach_spill_and_merge(
        arena_front(head, front, config, header, output_byte_limit, caps, &ledger, &qn_keys),
        &qn_keys.table,
        config.sort_order,
        output_byte_limit,
        threads,
        sink,
        SupplyTuning::default(),
        ZSTD_1,
        caps,
    )?;
    builder.build()?.run(PipelineConfig { threads, ..Default::default() })?;
    Ok(taps.collect())
}

/// Drive the block-input front alone (`BgzfBlockSource` → `ReadBlocks` →
/// `InflateToArena` → `front`'s steps) into a [`ChunkEventSink`], returning
/// every emitted chunk's records, its key kind, and the `AllAnnounced` counts.
fn drive_arena_chunks(
    front: ArenaFront,
    config: SortBufferConfig,
    header: &Header,
    blocks: Vec<crate::types::BgzfBlock>,
    threads: usize,
    caps: &ChainCaps,
) -> Result<ArenaFrontOutput> {
    const OUTPUT_BYTE_LIMIT: u64 = 4 * 1024 * 1024;
    let received: Arc<Mutex<ArenaFrontOutput>> = Arc::new(Mutex::new(ArenaFrontOutput::default()));
    let ledger = run_ledger_for(config, header);
    let builder = Pipeline::builder();
    let head = builder
        .chain(BgzfBlockSource::new(blocks, OUTPUT_BYTE_LIMIT).reader(caps.reader_source))
        .chain(ReadBlocks::new(config.memory_limit, OUTPUT_BYTE_LIMIT, Arc::clone(&ledger)))
        .chain(InflateToArena::new(OUTPUT_BYTE_LIMIT).with_phase_cap(caps.phase1.clone()));
    let qn_keys = fgumi_sort::QnSortState::new(false);
    arena_front(head, front, config, header, OUTPUT_BYTE_LIMIT, caps, &ledger, &qn_keys)
        .chain(ChunkEventSink { received: Arc::clone(&received) })
        .into_sink_marker();
    builder.build()?.run(PipelineConfig { threads, ..Default::default() })?;
    Ok(std::mem::take(&mut *received.lock()))
}

/// Drives the PRODUCTION front of `front`'s kind for `config.sort_order`:
/// `Front::Arena` dispatches per order through [`production_arena_front`].
/// The merge sink is `Serial`, which fits every thread count.
fn drive_front(
    front: Front,
    config: SortBufferConfig,
    header: &Header,
    records: &[RawRecord],
    pipeline_threads: usize,
) -> Result<ChainRun> {
    const OUTPUT_BYTE_LIMIT: u64 = 4 * 1024 * 1024;
    let sink = MergeSink::Batches(StepKind::Serial);
    match front {
        Front::SortBuffer => drive_buffer_chain(
            config,
            header,
            pack_batches(records, 256),
            OUTPUT_BYTE_LIMIT,
            pipeline_threads,
            sink,
            SupplyTuning::default(),
            ZSTD_1,
        ),
        Front::Arena => {
            let n_ref = u32::try_from(header.reference_sequences().len()).expect("n_ref fits u32");
            drive_arena_chain(
                production_arena_front(config.sort_order),
                config,
                header,
                bgzf_blocks_for(records, n_ref, 4096),
                OUTPUT_BYTE_LIMIT,
                pipeline_threads,
                sink,
                &ChainCaps::default(),
            )
        }
    }
}

// ── The output contract ─────────────────────────────────────────────────────

/// Every order's output is the stable sort of its input by key — equal keys
/// keep input order — at every thread count and across spill boundaries, on
/// both fronts.
#[rstest]
fn ties_across_run_boundaries_keep_input_order(
    #[values(
        SortOrder::Coordinate,
        SortOrder::TemplateCoordinate,
        SortOrder::Queryname(QuerynameComparator::Lexicographic),
        SortOrder::Queryname(QuerynameComparator::Natural)
    )]
    order: SortOrder,
    #[values(1usize, 4)] threads: usize,
    #[values(Front::Arena, Front::SortBuffer)] front: Front,
) {
    let (header, records) = tie_heavy_records(6_000, 11);
    let expected = reference_sort(order, &header, &records, None);
    assert_ne!(
        records.iter().map(|r| r.as_ref().to_vec()).collect::<Vec<_>>(),
        expected,
        "{order:?}: the fixture must not already be sorted"
    );
    let config = buffer_config(order, 96 * 1024, threads);
    let run = drive_front(front, config, &header, &records, threads).unwrap();
    assert!(
        run.merge_sources >= 3,
        "{order:?}: the fixture must merge ≥ 3 sources, got {}",
        run.merge_sources
    );
    assert!(
        equal_key_group_spans_a_boundary(order, &header, &expected, &run.spill_boundaries),
        "{order:?}: no equal-key group straddles a spill boundary — the test cannot catch a \
         cross-run reorder"
    );
    if matches!(order, SortOrder::Queryname(_)) {
        assert!(
            run.qn_classes >= 2,
            "{order:?}: the front learned no class, so the merge's shared-table compares go \
             unexercised"
        );
    }
    assert_record_parity(&format!("{order:?} {front:?} t{threads}"), &run.records, &expected);
    assert_input_order_within_equal_keys(order, &header, &run.records);
}

/// The oracle itself is stable, observed through the unique `XI` tags, and the
/// fixture really has ties in every order.
#[rstest]
fn reference_sort_is_stable(
    #[values(
        SortOrder::Coordinate,
        SortOrder::TemplateCoordinate,
        SortOrder::Queryname(QuerynameComparator::Lexicographic),
        SortOrder::Queryname(QuerynameComparator::Natural)
    )]
    order: SortOrder,
) {
    let (header, records) = tie_heavy_records(500, 3);
    let out = reference_sort(order, &header, &records, None);
    assert_input_order_within_equal_keys(order, &header, &out);
    assert!(max_equal_key_group(order, &header, &out) >= 2, "{order:?}: the fixture has no ties");
}

// ── Record-input chain parity (SortBuffer → spill split → supply → merge) ────

/// The record-input chain must produce the reference sort, byte for byte, for
/// EVERY sort order, across the in-memory and multi-spill regimes at 2 and 4
/// threads. It exercises `SortBuffer`'s streaming-spill emission, the spill
/// split's writes and slot opens (`file_id == seq`), and the single-residual
/// fast path all the way through the real `SortSpillDecompress` + `SortMerge`
/// — for coordinate, template, and both queryname comparators.
#[rstest]
#[case::coord_inmem_t2(SortOrder::Coordinate, 5_000, 256 * 1024 * 1024, 2)]
#[case::coord_inmem_t4(SortOrder::Coordinate, 5_000, 256 * 1024 * 1024, 4)]
#[case::coord_spill_t2(SortOrder::Coordinate, 20_000, 256 * 1024, 2)]
#[case::coord_spill_t4(SortOrder::Coordinate, 20_000, 256 * 1024, 4)]
#[case::template_inmem(SortOrder::TemplateCoordinate, 5_000, 256 * 1024 * 1024, 2)]
#[case::template_spill(SortOrder::TemplateCoordinate, 20_000, 256 * 1024, 2)]
#[case::template_spill_t4(SortOrder::TemplateCoordinate, 20_000, 256 * 1024, 4)]
// Many-spill regime (~20+ spill runs at a 128 KiB buffer): exercises the deeper
// per-slot FIFO (cap 32) and the emptiest-first refill order across many slots, where
// the merge must still emit byte-identical output to the oracle.
#[case::coord_manyspill_t4(SortOrder::Coordinate, 60_000, 128 * 1024, 4)]
#[case::template_manyspill_t4(SortOrder::TemplateCoordinate, 60_000, 128 * 1024, 4)]
#[case::qname_lex_inmem(SortOrder::Queryname(QuerynameComparator::Lexicographic), 5_000, 256 * 1024 * 1024, 2)]
#[case::qname_lex_spill(SortOrder::Queryname(QuerynameComparator::Lexicographic), 20_000, 256 * 1024, 2)]
#[case::qname_nat_inmem(SortOrder::Queryname(QuerynameComparator::Natural), 5_000, 256 * 1024 * 1024, 2)]
#[case::qname_nat_spill(SortOrder::Queryname(QuerynameComparator::Natural), 20_000, 256 * 1024, 2)]
fn sort_buffer_chain_matches_reference_all_orders(
    #[case] sort_order: SortOrder,
    #[case] n: usize,
    #[case] memory_limit: usize,
    #[case] threads: usize,
) {
    let (header, records) = synthesize_records(n, 0x5A17_C0DE);
    let sorter = buffer_config(sort_order, memory_limit, threads);
    let new_out = drive_sort_pipeline_tuned(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        // Run the pipeline with exactly the case's thread count — production
        // sort uses `num_threads` workers (no `max(3)` floor), so the parity
        // test must too.
        threads,
        StepKind::Exclusive,
        SupplyTuning::default(),
        ZSTD_1,
    )
    .expect("buffer pipeline drives to completion");

    let expected = reference_sort(sort_order, &header, &records, None);
    assert_record_parity(&format!("{sort_order:?}"), &new_out, &expected);
}

/// Raw BAM record carrying `MI:i:<mi>` (aux tag bytes + `i` type + i32 LE).
/// Shared with `sort_buffer`'s unit tests, which drive the same dropped-lane
/// rejection directly against `ingest_batch_records`.
pub(super) fn record_with_mi(pos: i32, name: &[u8], mi: i32) -> Vec<u8> {
    let mut aux = vec![b'M', b'I', b'i'];
    aux.extend_from_slice(&mi.to_le_bytes());
    make_bam_bytes(0, pos, 0, name, &[], 40, -1, -1, &aux)
}

/// An ingest failure must fail the whole pipeline, never produce a partial sort.
/// `SortBuffer` tracks "still ingesting" as `sorter.is_some()`, so a failure that
/// left that slot empty would be indistinguishable from "finalized" and could
/// surface as a clean `Finished` — a truncated result with a valid structure.
/// Drives the one reachable `ChunkSorter::push` failure: under
/// `--key-types none` the first record fixes the narrowed lane set, so a later
/// record with a differing MI is rejected.
#[test]
fn sort_buffer_chain_fails_the_pipeline_on_a_dropped_lane_violation() {
    let records: Vec<RawRecord> = (0..8i32)
        .map(|i| {
            // Record 0 fixes the variant with MI 1; record 4 violates it.
            let mi = if i == 4 { 2 } else { 1 };
            RawRecord::from(record_with_mi(100 - i, format!("rd_{i}").as_bytes(), mi))
        })
        .collect();

    let sorter = SortBufferConfig {
        key_types: KeyTypesSpec::None,
        ..buffer_config(SortOrder::TemplateCoordinate, 256 * 1024 * 1024, 1)
    };
    // Serial sink (production's actual `WriteBgzfFile` kind), not the parity
    // matrix's `Exclusive`: a sort chain declares a Detached step (the merge), so
    // at `--threads 1` it no longer fuses (see
    // `should_fuse_single_thread`/`has_detached_step`) and runs on the scheduled
    // path. `VecSource` above is itself `StepKind::Exclusive`, so an `Exclusive`
    // sink here would make TWO exclusive steps, and `assign_exclusive_owners`
    // fails when `total_exclusive (2) > n_threads (1)` — before the ingest
    // failure this test targets is even reached. A `Serial` sink leaves exactly
    // one exclusive step (the source), which fits one thread, and also matches
    // production's `WriteBgzfFile` kind — keeping the dropped-lane propagation
    // assertion meaningful at one thread. (A *lone* Exclusive step is fine at
    // `--threads 1`; only two or more exceed the budget — see `fused.rs`.)
    let err = drive_sort_pipeline_tuned(
        sorter,
        &Header::default(),
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        1,
        StepKind::Serial,
        SupplyTuning::default(),
        ZSTD_1,
    )
    .expect_err("a dropped-lane violation must fail the pipeline, not truncate the sort");

    let message = format!("{err:#}");
    assert!(
        message.contains("SortBuffer: push failed"),
        "the ingest failure must reach the caller: {message}"
    );
}

/// Regression guard for the `SortBuffer` peak-memory invariant: when a single
/// input `RecordBatch` is larger than `memory_limit`, `ingest_one_batch` seals
/// several chunks in one call (staging multiple `Spill` events into `pending`
/// before `emit_pending` runs). That transiently exceeds the "~one spill chunk"
/// production bound, but it MUST stay correct — every record emitted in sorted
/// order, none stranded by the loop. Packs all records into ONE oversized batch
/// against a small `memory_limit` (the case the 256-record-per-batch parity
/// matrix never hits) and asserts byte-for-byte parity with the reference sort
/// (coordinate is stable, so equal keys keep input order deterministically).
#[test]
fn sort_buffer_single_oversized_batch_seals_multiple_chunks_without_dropping() {
    let sort_order = SortOrder::Coordinate;
    let n = 20_000;
    // Far below the single batch's byte size, so the one batch seals many chunks.
    let memory_limit = 128 * 1024;
    let threads = 2;
    let (header, records) = synthesize_records(n, 0x0B16_BA7C);

    // One batch holding every record — forces multiple seals per ingest call.
    let batches = pack_batches(&records, records.len());
    assert_eq!(batches.len(), 1, "test must drive a single oversized batch");

    let sorter = buffer_config(sort_order, memory_limit, threads);
    let new_out = drive_sort_pipeline_tuned(
        sorter,
        &header,
        batches,
        4 * 1024 * 1024,
        threads,
        StepKind::Exclusive,
        SupplyTuning::default(),
        ZSTD_1,
    )
    .expect("buffer pipeline drives to completion");

    let expected = reference_sort(sort_order, &header, &records, None);

    assert_eq!(new_out.len(), records.len(), "every record must survive — none dropped");
    assert_eq!(new_out.len(), expected.len(), "record count matches oracle");
    for (i, (got, want)) in new_out.iter().zip(expected.iter()).enumerate() {
        assert_eq!(got, want, "record {i} bytes differ from oracle");
    }
}

/// Lever 1 parity gate: the terminal `SortMerge<BlockOutput>` path (framed
/// `DecompressedBlock`s, wired straight to `BgzfCompress`) must carry
/// byte-identical records to the reference sort for every sort order, across the
/// in-memory and spill regimes. The framed bytes are `[u32 LE len][body]` per
/// record (identical to the removed `SerializeRecordBatch`); un-framing them
/// recovers the same record bodies the `RecordBatchOutput` path emits. This is
/// the gate proving the merge-side framing did not change the output bytes.
#[rstest]
#[case::coord_inmem(SortOrder::Coordinate, 5_000, 256 * 1024 * 1024, 2)]
#[case::coord_spill(SortOrder::Coordinate, 20_000, 256 * 1024, 4)]
#[case::template_inmem(SortOrder::TemplateCoordinate, 5_000, 256 * 1024 * 1024, 2)]
#[case::template_spill(SortOrder::TemplateCoordinate, 20_000, 256 * 1024, 4)]
#[case::qname_lex_spill(SortOrder::Queryname(QuerynameComparator::Lexicographic), 20_000, 256 * 1024, 2)]
#[case::qname_nat_spill(SortOrder::Queryname(QuerynameComparator::Natural), 20_000, 256 * 1024, 2)]
fn sort_merge_block_output_matches_reference(
    #[case] sort_order: SortOrder,
    #[case] n: usize,
    #[case] memory_limit: usize,
    #[case] threads: usize,
) {
    let (header, records) = synthesize_records(n, 0x5A17_C0DE);
    let make_sorter = || buffer_config(sort_order, memory_limit, threads);

    // Terminal framed-block path.
    let block_out = drive_sort_block_pipeline(
        make_sorter(),
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        threads,
        StepKind::Exclusive,
    )
    .expect("block pipeline drives to completion");

    // Intermediate RecordBatch path — the framed bytes must un-frame to exactly
    // the records this path emits (cross-check the two SortMerge framings agree).
    let batch_out = drive_sort_pipeline_tuned(
        make_sorter(),
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        threads,
        StepKind::Exclusive,
        SupplyTuning::default(),
        ZSTD_1,
    )
    .expect("buffer pipeline drives to completion");

    let expected = reference_sort(sort_order, &header, &records, None);

    assert_record_parity(&format!("{sort_order:?} block path"), &block_out, &expected);
    assert_record_parity(&format!("{sort_order:?} RecordBatch path"), &batch_out, &expected);
}

/// L2.6: a `StepKind::Detached` SINK — the `WriteBgzfFile` analogue, i.e. the
/// detached-thread runtime driving a pure consumer — yields byte-identical
/// merged output to the reference sort, exactly like the pool-scheduled sink.
/// The chain's `SortMerge` is already `Detached`, so this drives the full chain
/// through `pipeline.run` with TWO off-pool threads (merge + sink) over a
/// multi-spill coordinate workload, pinning that `run_detached_driver` preserves
/// the output bytes for both the producing (merge) and consuming (sink) roles.
#[test]
fn detached_sink_chain_matches_reference_coordinate() {
    let memory_limit = 64 * 1024; // small → forces many real spill files
    let threads = 4;
    let (header, records) = synthesize_records(20_000, 0xD17A_C4ED);
    let sorter = buffer_config(SortOrder::Coordinate, memory_limit, threads);
    let detached_out = drive_sort_pipeline_tuned(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        threads,
        StepKind::Detached,
        SupplyTuning::default(),
        ZSTD_1,
    )
    .expect("detached-sink buffer pipeline drives to completion");

    let expected = reference_sort(SortOrder::Coordinate, &header, &records, None);

    assert_eq!(detached_out.len(), expected.len(), "record count mismatch");
    for (i, (got, want)) in detached_out.iter().zip(expected.iter()).enumerate() {
        assert_eq!(got, want, "detached-sink record {i} bytes differ from oracle");
    }
}

/// The buffer chain must hold coordinate parity across the merge supply's
/// knobs — fills small enough that frames straddle slices, four slices per
/// fill, and in-flight caps forced to bind — with a multi-spill workload that
/// forces real spill files. Guards against any spill-format / slot-ordering
/// drift between the spill split and the supply.
#[rstest]
#[case::default(SupplyTuning::default())]
#[case::fills_64k(SupplyTuning { fills: Some((64 << 10, 64 << 10)), ..SupplyTuning::default() })]
#[case::streams_4(SupplyTuning { read_streams: 4, ..SupplyTuning::default() })]
#[case::caps_bind(SupplyTuning {
    fills: Some((4 << 10, 4 << 10)),
    inflight_caps: Some((1, 4 << 10)),
    ..SupplyTuning::default()
})]
fn sort_buffer_chain_coordinate_matches_reference_across_supply_tunings(
    #[case] supply: SupplyTuning,
) {
    let (header, records) = synthesize_records(20_000, 0xBADD_CAFE);
    let memory_limit = 256 * 1024;
    let threads = 2;
    let sorter = buffer_config(SortOrder::Coordinate, memory_limit, threads);
    let new_out = drive_sort_pipeline_tuned(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        // Exactly the case's thread count (see the sibling test): production
        // runs the pipeline with `num_threads`, no `max(3)` floor.
        threads,
        StepKind::Exclusive,
        supply,
        ZSTD_1,
    )
    .expect("buffer pipeline drives to completion");

    let expected = reference_sort(SortOrder::Coordinate, &header, &records, None);

    assert_eq!(new_out.len(), expected.len(), "record count mismatch");
    for (i, (got, want)) in new_out.iter().zip(expected.iter()).enumerate() {
        assert_eq!(got, want, "record {i} bytes differ");
    }
}

/// Equal-key stability across spill boundaries — the output-identity-critical
/// tie-break the single-residual design depends on. Every record shares one
/// coordinate (tid 0, pos 0), so a stable sort must emit them in input order;
/// the buffer chain must preserve that across multiple spill chunks (tie-broken
/// by `file_id == seq`) and the residual. Pinned both directly (output == input
/// order) and against the reference sort.
#[test]
fn sort_buffer_chain_preserves_equal_key_input_order_across_spills() {
    let header = Header::default();
    // All at tid 0, pos 0 → identical coordinate key; names encode input order.
    // 6_000 records against a tiny memory limit forces several spill chunks.
    let records: Vec<RawRecord> = (0..6_000u32)
        .map(|i| {
            let name = format!("r{i:06}");
            RawRecord::from(make_bam_bytes(0, 0, 0, name.as_bytes(), &[], 80, -1, -1, &[]))
        })
        .collect();
    let input_bytes: RecordBytes = records.iter().map(|r| r.as_ref().to_vec()).collect();
    let memory_limit = 256 * 1024;

    let sorter = buffer_config(SortOrder::Coordinate, memory_limit, 2);
    let out = drive_sort_pipeline_tuned(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        2,
        StepKind::Exclusive,
        SupplyTuning::default(),
        ZSTD_1,
    )
    .expect("buffer pipeline drives to completion");

    assert_eq!(out.len(), input_bytes.len(), "record count mismatch");
    assert_eq!(out, input_bytes, "equal-key records must preserve input order across spills");

    let expected = reference_sort(SortOrder::Coordinate, &header, &records, None);
    assert_eq!(out, expected, "equal-key order must also match the reference sort");
}

// ── Merge-supply fill parity (fill size × codec × read streams) ──────────────

/// The streaming sort must produce byte-identical output whatever the fill
/// size — 1-byte and 7/3-byte fills (every frame, and every BGZF block's
/// 18-byte header, straddles slices), a frame-sized fill, and the production
/// sizes — for both spill codecs and one or four slices per fill, validated
/// against the reference sort. Tiny fills run a smaller fixture (each byte is
/// a read); every fixture spills to several files.
#[rstest]
fn three_step_chain_fill_matrix_matches_reference(
    #[values((1, 1), (7, 3), (64 << 10, 64 << 10), (4 << 20, 1 << 20))] fills: (u64, u64),
    #[values(SpillCodec::Zstd, SpillCodec::Bgzf)] codec: SpillCodec,
    #[values(1, 4)] read_streams: usize,
) {
    let sort_order = SortOrder::Coordinate;
    let threads = 4;
    let tiny = fills.0 < 1024;
    let (n_records, memory_limit) = if tiny { (2_000, 32 * 1024) } else { (30_000, 256 * 1024) };
    let (header, records) = synthesize_sized_records(n_records, 0x5EED_1234, 120);
    let sorter = buffer_config(sort_order, memory_limit, 2);
    let supply = SupplyTuning { fills: Some(fills), read_streams, inflight_caps: None };
    let new_out = drive_sort_pipeline_tuned(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        threads,
        StepKind::Exclusive,
        supply,
        (codec, 1),
    )
    .expect("pipeline drives to completion");
    let expected = reference_sort(sort_order, &header, &records, None);
    assert_eq!(new_out.len(), expected.len(), "record count mismatch ({supply:?}, {codec:?})");
    assert_eq!(new_out, expected, "sorted bytes differ ({supply:?}, {codec:?})");
}

/// Block-parallel decompression over BGZF spill files (the non-default codec)
/// must also match the reference sort — the matrix/soak/proptest exercise the
/// default zstd spills, so this confirms the block-parallel reorder path is
/// codec-agnostic. (Output records are codec-independent: the spill codec only
/// affects temp files, not the sorted output.)
#[test]
fn block_parallel_bgzf_spill_matches_reference() {
    let sort_order = SortOrder::Coordinate;
    let (header, records) = synthesize_sized_records(30_000, 0x5EED_BEEF, 120);
    let memory_limit = 256 * 1024;
    let sorter = buffer_config(sort_order, memory_limit, 2);
    let new_out = drive_sort_pipeline_tuned(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        8,
        StepKind::Exclusive,
        SupplyTuning { read_streams: 4, ..SupplyTuning::default() },
        (SpillCodec::Bgzf, 1),
    )
    .expect("pipeline drives to completion");
    let expected = reference_sort(sort_order, &header, &records, None);
    assert_eq!(new_out.len(), expected.len(), "record count mismatch (bgzf block-parallel)");
    assert_eq!(new_out, expected, "sorted bytes differ (bgzf block-parallel)");
}

// ── Production spill split + run extension (SpillGather run-former) ────────────
//
// The tests below exercise the run-forming `SpillGather` end-to-end:
// contiguous chunks coalesce into one run, and the merge sees one slot per run.

/// Drive the production spill split and return `(merged record bytes, run count)`.
/// The run count is the `SortMerge` stats `runs_written`, which equals the number
/// of spill files (`SpillReady`s) the merge saw — i.e. the number of runs the
/// `SpillGather` run-former formed. Mirrors `ChainBuilder::add_sort`'s wiring.
fn drive_spill_split_pipeline(
    config: SortBufferConfig,
    header: &Header,
    batches: Vec<RecordBatch>,
    output_byte_limit: u64,
    threads: usize,
) -> Result<(Vec<Vec<u8>>, usize)> {
    let run = drive_buffer_chain(
        config,
        header,
        batches,
        output_byte_limit,
        threads,
        MergeSink::Batches(StepKind::Serial),
        SupplyTuning::default(),
        ZSTD_1,
    )?;
    Ok((run.records, run.runs))
}

/// Records sorted into `order`, obtained by running them through the oracle and
/// parsing its (byte-identical) sorted output back into `RawRecord`s. Feeding
/// this back in is the already-sorted-input case the run-former must coalesce.
fn presorted_records(order: SortOrder, header: &Header, records: &[RawRecord]) -> Vec<RawRecord> {
    reference_sort(order, header, records, None).into_iter().map(RawRecord::from).collect()
}

/// An already-sorted input, re-sorted under a small memory budget through the
/// production split, must coalesce every contiguous chunk into ONE run — and
/// stay byte-identical to the oracle. This is the t16 fix's headline behaviour,
/// exercised across all four sort orders via the shared `SpillGather` run-former.
#[rstest]
#[case::coordinate(SortOrder::Coordinate)]
#[case::queryname_lex(SortOrder::Queryname(QuerynameComparator::Lexicographic))]
#[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
#[case::template(SortOrder::TemplateCoordinate)]
fn presorted_input_coalesces_to_one_run(#[case] order: SortOrder) {
    let (header, records) = synthesize_sized_records(20_000, 0xC0FF_EE00, 40);
    let sorted = presorted_records(order, &header, &records);

    let sorter = buffer_config(order, 256 * 1024, 2); // small ⇒ many contiguous chunks to coalesce
    let (out, runs) =
        drive_spill_split_pipeline(sorter, &header, pack_batches(&sorted, 256), 4 * 1024 * 1024, 4)
            .expect("spill-split pipeline drives to completion");

    let expected = reference_sort(order, &header, &sorted, None);
    assert_eq!(out, expected, "{order:?}: coalesced output diverged from the oracle");
    assert_eq!(
        runs, 1,
        "{order:?}: contiguous chunks must coalesce into a single run (SpillReady==#runs)"
    );
}

/// A reverse-ordered input cannot coalesce: every chunk's minimum precedes the
/// open run's maximum, so the run-former opens a fresh run per chunk. Output must
/// still be byte-identical to the oracle — extension is a layout optimisation,
/// never a reordering.
#[test]
fn reverse_ordered_input_forms_many_runs_but_stays_byte_parity() {
    let order = SortOrder::Coordinate;
    let (header, records) = synthesize_sized_records(20_000, 0xBEEF_0001, 40);
    let mut reversed = presorted_records(order, &header, &records);
    reversed.reverse();

    let sorter = buffer_config(order, 256 * 1024, 2);
    let (out, runs) = drive_spill_split_pipeline(
        sorter,
        &header,
        pack_batches(&reversed, 256),
        4 * 1024 * 1024,
        4,
    )
    .expect("spill-split pipeline drives to completion");

    let expected = reference_sort(order, &header, &reversed, None);
    assert_eq!(out, expected, "reverse-ordered sort diverged from the oracle");
    assert!(runs >= 2, "reverse-ordered chunks must not coalesce (got {runs} run(s))");
}

/// A fully shuffled input: the run-former coalesces the runs it can and opens new
/// ones where it cannot, and the merged output must still match the oracle
/// byte-for-byte. This is the partial-coalescing correctness case.
#[rstest]
#[case::coordinate(SortOrder::Coordinate)]
#[case::queryname_lex(SortOrder::Queryname(QuerynameComparator::Lexicographic))]
#[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
#[case::template(SortOrder::TemplateCoordinate)]
fn shuffled_input_matches_the_oracle_through_the_split(#[case] order: SortOrder) {
    let (header, records) = synthesize_sized_records(20_000, 0x5EED_1234, 40);

    let sorter = buffer_config(order, 256 * 1024, 2);
    let (out, runs) = drive_spill_split_pipeline(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        4,
    )
    .expect("spill-split pipeline drives to completion");

    let expected = reference_sort(order, &header, &records, None);
    assert_eq!(out, expected, "{order:?}: shuffled-input sort diverged from the oracle");
    // Shuffled input under a tiny budget must open more than one run (partial
    // coalescing: a new run wherever a chunk cannot extend the previous one).
    assert!(runs >= 2, "{order:?}: shuffled input must form >=2 runs (got {runs})");
}

/// An in-memory sort (memory budget above the whole input) never spills, so the
/// run-former forms zero runs — the residual fast path stays intact under the
/// production split.
#[test]
fn in_memory_input_forms_no_runs() {
    let order = SortOrder::Coordinate;
    let (header, records) = synthesize_sized_records(5_000, 0xA110_0C00, 40);

    let sorter = buffer_config(order, 256 * 1024 * 1024, 2); // above the whole input ⇒ no spill
    let (out, runs) = drive_spill_split_pipeline(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        4,
    )
    .expect("spill-split pipeline drives to completion");

    let expected = reference_sort(order, &header, &records, None);
    assert_eq!(out, expected, "in-memory sort diverged from the oracle");
    assert_eq!(runs, 0, "a budget-sized run never spills");
}

/// The single-thread scheduled path the fusion-policy change (t1 fix) newly
/// activates: at `threads == 1` this chain declares Detached steps (`SpillGather`
/// and `SortMerge`), so it does NOT fuse and runs on the scheduled path (one pool
/// worker plus the detached driver threads). Pin that path's spilling output
/// against the oracle end-to-end through the production spill split — the
/// buffer-chain parity matrix only covers threads 2/4, and its harness cannot
/// take a `threads == 1` case (its `VecSource` + `VecSink` are both Exclusive).
#[rstest]
#[case::coordinate(SortOrder::Coordinate)]
#[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
fn sorts_at_one_thread_through_the_scheduled_path(#[case] order: SortOrder) {
    let (header, records) = synthesize_sized_records(8_000, 0x0057_0001, 40);

    let sorter = buffer_config(order, 256 * 1024, 1); // small ⇒ spills, exercising the merge at t1
    let (out, runs) = drive_spill_split_pipeline(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        1,
    )
    .expect("spill-split pipeline drives to completion at one thread");

    let expected = reference_sort(order, &header, &records, None);
    assert_eq!(out, expected, "{order:?}: single-thread scheduled sort diverged from the oracle");
    // The small budget must actually spill; otherwise this degrades to an
    // in-memory sort and no longer exercises the merge path it names.
    assert!(runs >= 1, "{order:?}: the small budget must spill (got {runs} run(s))");
}

/// Out-of-order supply — slices land out of order (four per fill), workers
/// decompress one file's blocks concurrently, and small fills make frames
/// straddle slices — must reassemble byte-identically to the production-size
/// single-stream supply. Property test over record counts, fill sizes and a
/// high pipeline-thread count; a straggler finalizing EOF while another
/// worker holds an in-flight block must not truncate the output (record count
/// is asserted equal).
#[cfg(test)]
// Soak/matrix/proptest suites: multi-minute, so gated off the default test
// target and run on the nightly `cargo ci-test-stress` job instead.
#[cfg(feature = "stress-tests")]
mod proptests {
    use super::*;
    use proptest::prelude::*;
    proptest! {
        #![proptest_config(ProptestConfig { cases: 24, ..ProptestConfig::default() })]
        #[test]
        fn small_fill_supply_matches_production_fills(
            n_records in 2_000usize..18_000,
            fill in 512u64..=65_536,
            read_streams in 1usize..=4,
            seed in any::<u64>(),
        ) {
            let sort_order = SortOrder::Coordinate;
            let (header, records) = synthesize_sized_records(n_records, seed, 100);
            // Force spilling so slots (and the reorder path) are exercised.
            let memory_limit = 256 * 1024;
            let pipeline_threads = 6;
            let make_sorter = || buffer_config(sort_order, memory_limit, 2);
            let production = drive_sort_pipeline_tuned(
                make_sorter(),
                &header,
                pack_batches(&records, 256),
                4 * 1024 * 1024,
                pipeline_threads,
                StepKind::Exclusive,
                SupplyTuning::default(),
                ZSTD_1,
            ).expect("production-fill pipeline");
            let small = drive_sort_pipeline_tuned(
                make_sorter(),
                &header,
                pack_batches(&records, 256),
                4 * 1024 * 1024,
                pipeline_threads,
                StepKind::Exclusive,
                SupplyTuning { fills: Some((fill, fill)), read_streams, inflight_caps: None },
                ZSTD_1,
            ).expect("small-fill pipeline");
            prop_assert_eq!(small.len(), records.len(), "no truncation");
            prop_assert_eq!(small, production, "small fills diverge from production fills");
        }
    }
}
/// Soak-iteration count, collapsed to a minimal pass under coverage
/// instrumentation.
///
/// `cargo llvm-cov` (the coverage CI job) sets `cfg(coverage)`. Under it these
/// soak tests run only to exercise the ~200 lines nothing else reaches, and a
/// couple of iterations covers every one of them — so we drop to `min` there and
/// skip the sustained repetition whose sole purpose is hunting rare races. The
/// FULL soak (its actual race-hunting value) still runs on the nightly
/// `cargo ci-test-stress` job, which does NOT set `cfg(coverage)`, so nothing is
/// lost on the path where the soak matters; only the coverage job — which repeats
/// the same covered lines dozens of times for no extra coverage — gets faster.
#[cfg(feature = "stress-tests")]
fn soak_iterations(full: usize, min: usize) -> usize {
    if cfg!(coverage) { min } else { full }
}

#[cfg(feature = "stress-tests")]
/// Maximum-contention soak for the merge supply. Drives the path repeatedly under the most
/// adversarial settings the knobs allow — many spill files, a tiny reorder
/// window (so stragglers continuously hit the reorder-window backpressure and
/// the drain-only path), 4 KiB fills in four slices with the in-flight caps at
/// two slices (maximum per-slice churn and the most frequent `reader_eof` /
/// `in_flight` transitions), and far more
/// pipeline worker threads (12) than sorter threads (so many workers race to
/// decompress one file's blocks concurrently and finalize out of order).
///
/// Each iteration uses a fresh random seed and is checked for *byte-identity*
/// against the reference sort — so a lost, duplicated,
/// or reordered block (the truncation class the `reader_eof`/`in_flight`
/// protocol guards against) fails the assertion. Each iteration runs under a
/// per-iteration wall-clock watchdog (via `run_watchdogged_parity`): a livelock
/// (e.g. `queue_eof` never finalizing) trips the timeout and fails the test
/// instead of hanging CI.
///
/// This complements the loom model (exhaustive but tiny) and the proptest
/// (random sizes, moderate threads) by hammering the REAL pipeline under
/// sustained high contention for many iterations.
#[test]
fn block_parallel_high_contention_soak_matches_reference() {
    use std::time::Duration;

    const RECORDS_PER_ITER: usize = 15_000;
    const PIPELINE_THREADS: usize = 12;
    const SORTER_THREADS: usize = 2;
    // Tiny per-thread sort memory ⇒ many spill files ⇒ many slot readers.
    const MEMORY_LIMIT: usize = 128 * 1024;
    // Tiny output/reorder-window budget ⇒ the block-parallel reorder window is
    // ~1 block, so `bp_reorder_admits` backpressures aggressively and workers
    // are repeatedly forced through the Phase-B drain-only path.
    const OUTPUT_BYTE_LIMIT: u64 = 64 * 1024;
    const SUPPLY: SupplyTuning = SupplyTuning {
        fills: Some((4 << 10, 4 << 10)),
        read_streams: 4,
        inflight_caps: Some((2, 8 << 10)),
    };
    // Per-iteration watchdog: a livelock in any single iteration fails fast.
    const WATCHDOG: Duration = Duration::from_secs(60);

    let iterations = soak_iterations(40, 2);
    let sort_order = SortOrder::Coordinate;
    for iter in 0..iterations {
        let seed = 0xA5A5_0000_u64 ^ (iter as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        let (header, records) = synthesize_sized_records(RECORDS_PER_ITER, seed, 110);
        run_watchdogged_parity(
            &format!("hc-soak-i{iter}"),
            sort_order,
            header,
            records,
            MEMORY_LIMIT,
            SORTER_THREADS,
            PIPELINE_THREADS,
            OUTPUT_BYTE_LIMIT,
            SUPPLY,
            WATCHDOG,
        );
    }
}

#[cfg(feature = "stress-tests")]
/// Maximum-churn soak over the `SortBuffer` chain: a *tiny* memory limit seals a
/// run after only a handful of records, so the chain cycles
/// seal → compress → spill → decompress → merge as fast as it can while eight
/// pipeline threads contend for the `Serial` steps. A wedge (a step parking
/// while a downstream holds the chunk it is waiting on) trips the per-iteration
/// watchdog and fails fast instead of hanging, and byte-parity against the
/// reference sort proves nothing is lost, duplicated, or reordered under churn.
///
/// NOTE: this drives `drive_sort_pipeline_tuned` (the record-input `SortBuffer`
/// chain), NOT the block-input arena front. The front's capacity-1 arena cycle —
/// `ReadBlocks` acquire+admit and the split front's run sort seal+free on the
/// coordination driver with `InflateToArena` on the pool in between — is covered
/// by `arena_front_chain_seals_multiple_runs_without_losing_records`, which seals
/// several runs through the real runtime and so exercises acquire/seal/free
/// across runs. That test is not a *soak*: there is no watchdogged high-churn
/// coverage of the arena front yet.
#[test]
fn sort_buffer_chain_tight_memory_soak_no_deadlock() {
    use std::time::Duration;

    const RECORDS_PER_ITER: usize = 12_000;
    const PIPELINE_THREADS: usize = 8;
    const SORTER_THREADS: usize = 4;
    // Very tight sort memory ⇒ `SortBuffer` seals after only a handful of
    // records ⇒ maximal seal/spill/merge churn across the pipeline threads.
    const MEMORY_LIMIT: usize = 24 * 1024;
    const OUTPUT_BYTE_LIMIT: u64 = 48 * 1024;
    const WATCHDOG: Duration = Duration::from_secs(60);

    let iterations = soak_iterations(24, 2);
    for iter in 0..iterations {
        let seed = 0xC0FF_EE00_u64 ^ (iter as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        let (header, records) = synthesize_sized_records(RECORDS_PER_ITER, seed, 130);
        run_watchdogged_parity(
            &format!("tight-memory-soak-i{iter}"),
            SortOrder::Coordinate,
            header,
            records,
            MEMORY_LIMIT,
            SORTER_THREADS,
            PIPELINE_THREADS,
            OUTPUT_BYTE_LIMIT,
            SupplyTuning { fills: Some((4 << 10, 4 << 10)), ..SupplyTuning::default() },
            WATCHDOG,
        );
    }
}

/// Assert two record-byte streams are identical, reporting only the FIRST
/// mismatch (index + lengths + 16-byte prefixes). A blanket `assert_eq!` on the
/// two `Vec<Vec<u8>>` would dump tens of thousands of binary records into CI
/// logs on divergence; this keeps a failure readable while still catching any
/// lost / duplicated / reordered / corrupted record.
fn assert_record_parity(label: &str, actual: &[Vec<u8>], expected: &[Vec<u8>]) {
    if let Some((idx, (a, b))) = actual.iter().zip(expected).enumerate().find(|(_, (a, b))| a != b)
    {
        panic!(
            "{label}: output diverges from the reference sort at record {idx}: \
             actual_len={}, expected_len={}, actual_prefix={:?}, expected_prefix={:?}",
            a.len(),
            b.len(),
            &a[..a.len().min(16)],
            &b[..b.len().min(16)],
        );
    }
    // No record differs within the common prefix; a length delta is the only
    // remaining divergence (truncation or duplication).
    assert_eq!(
        actual.len(),
        expected.len(),
        "{label}: record count mismatch (truncation/duplication)"
    );
}

/// Run one watchdog-guarded streaming sort and assert byte-identity against the
/// reference sort. Owns `header`/`records` so the worker thread can take them.
///
/// Both the reference sort AND the streaming pipeline run *inside* the worker, so
/// the `recv_timeout` watchdog covers both — a stall in either (the path under
/// test or, defensively, the reference sorter) fails the test fast instead of
/// hanging the test process. The merge sink is `Serial` (every watchdog'd parity
/// caller drives the streaming three-step chain). Panics (failing the test) on
/// divergence, pipeline/oracle error, watchdog timeout (livelock), or a worker
/// panic — never hangs.
#[allow(clippy::too_many_arguments)]
fn run_watchdogged_parity(
    label: &str,
    sort_order: SortOrder,
    header: Header,
    records: Vec<RawRecord>,
    memory_limit: usize,
    sorter_threads: usize,
    pipeline_threads: usize,
    output_byte_limit: u64,
    tuning: SupplyTuning,
    watchdog: std::time::Duration,
) {
    use std::sync::mpsc;

    let (tx, rx) = mpsc::channel();
    let worker = std::thread::Builder::new()
        .name(label.to_string())
        .spawn(move || {
            // (pipeline_out, expected) — both computed under the watchdog.
            let result = (|| -> Result<(RecordBytes, RecordBytes)> {
                let expected = reference_sort(sort_order, &header, &records, None);
                let batches = pack_batches(&records, 256);
                let sorter = buffer_config(sort_order, memory_limit, sorter_threads);
                let out = drive_sort_pipeline_tuned(
                    sorter,
                    &header,
                    batches,
                    output_byte_limit,
                    pipeline_threads,
                    StepKind::Serial,
                    tuning,
                    ZSTD_1,
                )?;
                Ok((out, expected))
            })();
            let _ = tx.send(result);
        })
        .expect("spawn soak worker");

    match rx.recv_timeout(watchdog) {
        Ok(Ok((out, expected))) => {
            assert_record_parity(label, &out, &expected);
            worker.join().expect("soak worker panicked");
        }
        Ok(Err(e)) => panic!("{label}: pipeline/oracle errored: {e:#}"),
        Err(mpsc::RecvTimeoutError::Timeout) => {
            panic!("{label}: DEADLOCK/LIVELOCK — sort did not complete within {watchdog:?}")
        }
        Err(mpsc::RecvTimeoutError::Disconnected) => {
            panic!("{label}: soak worker dropped its sender (panicked?)")
        }
    }
}

#[cfg(feature = "stress-tests")]
/// Spill-pressure regimes for the block-parallel soak matrix. Both force real
/// spill files (so slot readers and the Phase-2 reorder path run); they differ
/// in how the spilled data is shaped across files.
#[derive(Clone, Copy, Debug)]
enum SoakRegime {
    /// Larger-than-budget workload: total spilled bytes vastly exceed the
    /// in-memory sort budget, producing MANY small spill files (tens). This is
    /// the mandatory larger-than-RAM run — the external merge over many slot
    /// readers, with frequent cross-file `reader_eof`/`in_flight` transitions,
    /// is the case the truncation protocol must survive.
    ManySmallFiles,
    /// A handful of LARGE spill files: the budget admits a big batch, so only a
    /// few (but > 1) files spill, each with many blocks. Stresses long per-file
    /// block runs and the per-slot reorder window rather than cross-file churn.
    FewLargeFiles,
}

#[cfg(feature = "stress-tests")]
struct SoakParams {
    records: usize,
    seq_len: usize,
    memory_limit: usize,
    output_byte_limit: u64,
    /// Fill size of the supply.
    fill: u64,
}

#[cfg(feature = "stress-tests")]
impl SoakRegime {
    /// Discriminant folded into the per-case seed so each regime sorts a
    /// distinct record set.
    fn seed_salt(self) -> u64 {
        match self {
            SoakRegime::ManySmallFiles => 0x1111_1111_1111_1111,
            SoakRegime::FewLargeFiles => 0x2222_2222_2222_2222,
        }
    }

    fn params(self) -> SoakParams {
        match self {
            // ~40k records ≈ 10 MB spilled into many (~100) small files at a
            // 96 KiB budget, with a tiny reorder window and 4 KiB fills (max
            // per-slice churn and the most `reader_eof`/`in_flight` events).
            SoakRegime::ManySmallFiles => SoakParams {
                records: 40_000,
                seq_len: 150,
                memory_limit: 96 * 1024,
                output_byte_limit: 128 * 1024,
                fill: 4 << 10,
            },
            // ~12k × 150B ≈ 1.8 MB spilled into ~4 large files at a 512 KiB
            // budget, with a roomy window and production-size fills.
            SoakRegime::FewLargeFiles => SoakParams {
                records: 12_000,
                seq_len: 150,
                memory_limit: 512 * 1024,
                output_byte_limit: 4 * 1024 * 1024,
                fill: 1 << 20,
            },
        }
    }
}

#[cfg(feature = "stress-tests")]
/// External-watchdog soak MATRIX for the merge supply. Crosses pipeline-thread
/// count × in-flight caps (production, or forced to one slice and one fill of
/// cold bytes) × spill regime, so the reorder/in-flight/EOF protocol is
/// hammered across {1, 2, 8} workers, {many small, few large} spill-file
/// shapes, and both a free-running and a cap-bound planner. Each (case × iteration) runs under a wall-clock
/// watchdog and is checked byte-for-byte against the reference sort, so a
/// livelock fails fast and any lost / duplicated / reordered block is caught.
///
/// `rstest` generates the full cross product (3 × 2 × 2 = 12 cases); nextest runs
/// them as independent parallel tests. The mandatory larger-than-budget run is
/// `ManySmallFiles` (~100 spill files at a 96 KiB budget). This is the real P5
/// hardening gate the OFF-default flip rests on; it complements the loom model
/// (exhaustive but tiny), the proptest (random sizes), and the single-corner
/// high-contention soak (12 threads, 1-block window).
#[rstest]
fn block_parallel_soak_matrix_matches_reference(
    #[values(1, 2, 8)] pipeline_threads: usize,
    #[values(false, true)] caps_bind: bool,
    #[values(SoakRegime::ManySmallFiles, SoakRegime::FewLargeFiles)] regime: SoakRegime,
) {
    use std::time::Duration;

    const SORTER_THREADS: usize = 2;
    const WATCHDOG: Duration = Duration::from_secs(120);

    let iterations = soak_iterations(4, 1);
    let sort_order = SortOrder::Coordinate;
    let SoakParams { records, seq_len, memory_limit, output_byte_limit, fill } = regime.params();
    let supply = SupplyTuning {
        fills: Some((4 * fill, fill)),
        read_streams: 4,
        inflight_caps: caps_bind.then_some((1, fill)),
    };

    for iter in 0..iterations {
        // Distinct seed per (regime, threads, granularity, iter).
        let seed = 0x50A4_0000_u64
            .wrapping_add((iter as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15))
            .wrapping_add((pipeline_threads as u64) << 40)
            .wrapping_add(u64::from(caps_bind) << 32)
            ^ regime.seed_salt();
        let (header, recs) = synthesize_sized_records(records, seed, seq_len);
        let label =
            format!("soak-{regime:?}-t{pipeline_threads}-caps{caps_bind}-fill{fill}-i{iter}");
        run_watchdogged_parity(
            &label,
            sort_order,
            header,
            recs,
            memory_limit,
            SORTER_THREADS,
            pipeline_threads,
            output_byte_limit,
            supply,
            WATCHDOG,
        );
    }
}

fn synthesize_sized_records(n: usize, seed: u64, seq_len: usize) -> (Header, Vec<RawRecord>) {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(1);
    let mut next_u32 = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        #[allow(clippy::cast_possible_truncation)]
        let v = state as u32;
        v
    };

    // Reference sequences must exist for the mapped records below to carry a
    // meaningful coordinate key.
    let header = synthetic_header();
    let mut records = Vec::with_capacity(n);
    for i in 0..n {
        let name = format!("rd_{}", next_u32() % 100_000);
        let pos: i32 = (next_u32() % 1_000_000).cast_signed();
        let is_paired = i % 2 == 0;
        // Most records are MAPPED across a handful of references. With
        // `tid = -1` everywhere, `extract_coordinate_key_inline` returns
        // `RawCoordinateKey::unmapped()` (`u64::MAX`) for every record, `pos` is
        // never read, and the coordinate cases degenerate into one equal-key
        // bucket — they would pass even if the key comparison were broken. Every
        // eighth record stays unmapped so the equal-key/tie path is still covered.
        let unmapped = i % 8 == 0;
        #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
        let tid: i32 = if unmapped { -1 } else { (i % N_TEST_REFS) as i32 };
        let flags: u16 = if unmapped { 0x4 } else { 0 } | if is_paired { 0x1 | 0x8 } else { 0 };
        let bytes = make_bam_bytes(tid, pos, flags, name.as_bytes(), &[], seq_len, -1, -1, &[]);
        records.push(RawRecord::from(bytes));
    }
    (header, records)
}

#[rstest]
#[case::t1(1)]
#[case::t2(2)]
#[case::t4(4)]
#[case::t8(8)]
fn three_step_chain_large_spill_completes(#[case] pipeline_threads: usize) {
    use std::time::Duration;

    let (header, records) = synthesize_sized_records(60_000, 0xBADD_CAFE, 200);
    // Production supply tuning under a per-case watchdog: the regression this pins is a deadlock at high
    // pipeline-thread counts, so the watchdog converts a hang into a failure.
    run_watchdogged_parity(
        &format!("large-spill-t{pipeline_threads}"),
        SortOrder::Coordinate,
        header,
        records,
        1024 * 1024, // sort memory_limit
        2,           // sorter_threads
        pipeline_threads,
        256 * 1024 * 1024, // output queue limit
        SupplyTuning::default(),
        Duration::from_secs(90),
    );
}

#[test]
fn three_step_chain_empty_input_drains_cleanly() {
    let header = Header::default();
    let sorter = buffer_config(SortOrder::Coordinate, 1024 * 1024, 1);
    let out = drive_sort_pipeline(sorter, &header, Vec::new(), 64 * 1024, 3, StepKind::Exclusive)
        .expect("empty pipeline");
    assert!(out.is_empty());
}

// ── SortMerge fail-closed regression tests ──────────────────────────────────

/// `Exclusive` source that drains a `Vec<SortPhase2Event>` one event per
/// `try_run`, feeding `SortMerge` directly. Used to drive the merge into a
/// drained-but-incomplete-setup state without standing up the spill machinery.
struct Phase2EventSource {
    events: Vec<crate::sort::protocol::SortPhase2Event>,
    held: HeldSlot<Unpushed<crate::sort::protocol::SortPhase2Event>>,
    output_byte_limit: u64,
}

impl Phase2EventSource {
    fn new(
        mut events: Vec<crate::sort::protocol::SortPhase2Event>,
        output_byte_limit: u64,
    ) -> Self {
        events.reverse();
        Self { events, held: HeldSlot::new(), output_byte_limit }
    }
}

impl Step for Phase2EventSource {
    type Input = ();
    type Outputs = fgumi_pipeline_core::outputs::Single<crate::sort::protocol::SortPhase2Event>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "Phase2EventSource",
            kind: StepKind::Exclusive,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
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
        let Some(event) = self.events.pop() else {
            return Ok(StepOutcome::Finished);
        };
        match ctx.outputs.push(event) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }
}

fn run_merge_over_events(
    events: Vec<crate::sort::protocol::SortPhase2Event>,
) -> Result<Vec<Vec<u8>>> {
    // Delegates to `collect_merge_batches` so the merge chain is wired in exactly
    // one place; this helper only flattens the batches into per-record bytes.
    let batches = collect_merge_batches(events, 1 << 20, 256)?;
    let mut out = Vec::new();
    for batch in batches {
        for bytes in batch.iter_record_bytes() {
            out.push(bytes.to_vec());
        }
    }
    Ok(out)
}

/// A wholly empty event stream (no payload, no `AllAnnounced`) is the one
/// legitimate drained-but-not-ready case: it merges to an empty output.
#[test]
fn test_sort_merge_empty_input_merges_to_empty() {
    let out = run_merge_over_events(Vec::new()).expect("empty merge should succeed");
    assert!(out.is_empty(), "expected no records, got {}", out.len());
}

/// An `AllAnnounced` that promises a slot which never arrives leaves the setup
/// incomplete when the input drains; `SortMerge` must fail closed rather than
/// silently merge a partial result.
#[test]
fn test_sort_merge_fails_closed_on_incomplete_setup() {
    let events = vec![crate::sort::protocol::SortPhase2Event::AllAnnounced {
        slot_count: 1,
        memory_chunk_count: 0,
        total_records: 0,
    }];
    let err = run_merge_over_events(events).expect_err("incomplete setup must error");
    let msg = err.to_string();
    assert!(msg.contains("setup incomplete at input drain"), "unexpected error message: {msg}");
}

// ── SortMerge output-buffer sizing + duplicate-AllAnnounced regression ───────

/// Drive `SortMerge` over `events` and return the *batches* it emits (not the
/// flattened records), so tests can inspect per-batch buffer capacity.
fn collect_merge_batches(
    events: Vec<crate::sort::protocol::SortPhase2Event>,
    output_byte_limit: u64,
    target_batch_count: usize,
) -> Result<Vec<RecordBatch>> {
    Ok(collect_merge_batches_with(
        events,
        output_byte_limit,
        target_batch_count,
        MergeRunOptions::default(),
    )?
    .batches)
}

/// How [`collect_merge_batches_with`] runs the merge: the merge demand handed
/// to `SortMerge`, the pool size, per-thread test timers, a sink-side record
/// counter a feeder can gate on, an optional pool-side supply between the
/// source and the merge, and the sink's kind.
#[derive(Default)]
struct MergeRunOptions {
    demand: Option<Arc<fgumi_sort::MergeDemand>>,
    /// The supply ledger handed to `SortMerge` (self-serve claims).
    ledger: Option<Arc<crate::sort::supply_ledger::SupplyLedger>>,
    /// `0` means one thread.
    threads: usize,
    test_backoff: Vec<fgumi_pipeline_core::runtime::TestBackoff>,
    progress: Option<Arc<std::sync::atomic::AtomicUsize>>,
    /// A pool step that delivers the slots' blocks (see [`FeedOnDemand`]).
    feeder: Option<FeedOnDemand>,
    /// Hold the setup's `AllAnnounced` back until a pool worker is parked
    /// (see [`HoldUntilPoolParked`]), so the merge's first stall finds one.
    hold_until_parked: bool,
    /// `None` means `Serial` (a pool step).
    sink_kind: Option<StepKind>,
}

/// What [`collect_merge_batches_with`] observed.
struct MergeRun {
    batches: Vec<RecordBatch>,
    /// The run's per-step stats, by step name.
    stats: fgumi_pipeline_core::runtime::StatsSnapshot,
}

impl MergeRun {
    /// The `SortMerge` step's stats.
    fn merge_stats(&self) -> &fgumi_pipeline_core::runtime::StepStatsSnapshot {
        &self.stats.steps.iter().find(|(name, _)| *name == "SortMerge").expect("SortMerge").1
    }

    /// Positions of every merged record, in output order.
    fn positions(&self) -> Vec<i32> {
        self.batches
            .iter()
            .flat_map(|b| b.iter_record_bytes().map(fgumi_raw_bam::pos).collect::<Vec<_>>())
            .collect()
    }
}

/// [`collect_merge_batches`] with a merge demand, a pool size, per-thread test
/// timers, a sink progress counter and an optional pool-side supply.
fn collect_merge_batches_with(
    events: Vec<crate::sort::protocol::SortPhase2Event>,
    output_byte_limit: u64,
    target_batch_count: usize,
    opts: MergeRunOptions,
) -> Result<MergeRun> {
    let received: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let source = Phase2EventSource::new(events, output_byte_limit);
    let mut merge = SortMerge::<RecordBatchOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        output_byte_limit,
        target_batch_count,
    );
    match (opts.demand, opts.ledger) {
        (Some(demand), Some(ledger)) => {
            merge = merge.with_spill_supply(&SpillSupply {
                demand,
                ledger,
                retire_epoch: Arc::default(),
            });
        }
        (Some(demand), None) => merge = merge.with_merge_demand(demand),
        (None, Some(_)) => panic!("a supply ledger needs the supply's merge demand"),
        (None, None) => {}
    }
    let sink = VecSink {
        received: Arc::clone(&received),
        kind: opts.sink_kind.unwrap_or(StepKind::Serial),
        progress: opts.progress,
    };

    let builder = Pipeline::builder();
    let chain = builder.chain(source);
    let chain = match opts.feeder {
        Some(feeder) if opts.hold_until_parked => {
            chain.chain(feeder).chain(HoldUntilPoolParked::new(output_byte_limit)).chain(merge)
        }
        Some(feeder) => chain.chain(feeder).chain(merge),
        None => chain.chain(merge),
    };
    chain.chain(sink).into_sink_marker();
    let pipeline = builder.build()?;
    let stats = pipeline.stats();
    pipeline.run(PipelineConfig {
        threads: opts.threads.max(1),
        test_backoff: opts.test_backoff,
        stats: Some(Arc::clone(&stats)),
        ..Default::default()
    })?;

    Ok(MergeRun { batches: std::mem::take(&mut *received.lock()), stats: stats.snapshot() })
}

/// One record per block for the gated feeds: a coordinate-sortable BAM record
/// at `pos` on reference 0, framed as a spill block (the coordinate key is
/// embedded in the record, so the frame is `[u32 len][record]`).
fn gated_block(pos: usize) -> Vec<u8> {
    let pos = i32::try_from(pos).expect("small test positions");
    let record = make_bam_bytes(0, pos, 0, format!("r{pos:04}").as_bytes(), &[], 4, -1, -1, &[]);
    let mut block = Vec::new();
    fgumi_sort::frame_keyed_record_into(
        &mut block,
        &fgumi_sort::RawCoordinateKey::default(),
        &record,
    )
    .expect("frame record");
    block
}

/// Two empty, open `SortMergeSlot`s over tempfiles (file ids 0 and 1).
fn two_open_slots() -> Vec<Arc<fgumi_sort::SortMergeSlot>> {
    (0..2u32)
        .map(|file_id| Arc::new(fgumi_sort::SortMergeSlot::for_test(file_id, SpillCodec::Bgzf)))
        .collect()
}

/// The phase-2 setup events announcing `slots` (no memory chunks) for a merge
/// of `total` records.
fn announce_slots(
    slots: &[Arc<fgumi_sort::SortMergeSlot>],
    total: usize,
) -> Vec<crate::sort::protocol::SortPhase2Event> {
    use crate::sort::protocol::SortPhase2Event;
    let mut events: Vec<SortPhase2Event> = slots
        .iter()
        .map(|slot| SortPhase2Event::SpillReady {
            slot: Arc::clone(slot),
            path: std::path::PathBuf::new(),
            records_ingested_so_far: total as u64,
        })
        .collect();
    events.push(SortPhase2Event::AllAnnounced {
        slot_count: u32::try_from(slots.len()).expect("few slots"),
        memory_chunk_count: 0,
        total_records: total as u64,
    });
    events
}

/// Two `SortMergeSlot`s over tempfiles, announced through phase-2 events, and a
/// feeder thread that pushes one single-record block at a time into alternating
/// slots (record `j` at position `j` into slot `j % 2`, so the merge switches
/// sources every record), then marks both slots EOF. Record `j` is pushed only
/// once the sink holds `j - 1` records, and `delay` after that: the merge
/// needs both slots' heads to emit, so every push lands while the merge is
/// parked or about to park on the slot it needs. The EOFs land only once the
/// sink holds every record but the last (and `delay` after that): the last
/// record needs the other slot's EOF to win, so the run ends on a run-end EOF
/// stall that only the EOF's notify can end. Each delivery is the producer half of the wake
/// protocol without a `SortSpillDecompress`: push under `decompressed`,
/// release, then `notify_delivered`.
fn gated_two_slot_feed(
    n_blocks_per_slot: usize,
    demand: &Arc<fgumi_sort::MergeDemand>,
    progress: &Arc<std::sync::atomic::AtomicUsize>,
    delay: std::time::Duration,
) -> (Vec<crate::sort::protocol::SortPhase2Event>, std::thread::JoinHandle<()>) {
    use std::sync::atomic::Ordering;
    let slots = two_open_slots();
    let total = 2 * n_blocks_per_slot;
    let events = announce_slots(&slots, total);
    let demand = Arc::clone(demand);
    let progress = Arc::clone(progress);
    let feeder = std::thread::spawn(move || {
        for j in 0..total {
            while progress.load(Ordering::Acquire) < j.saturating_sub(1) {
                std::thread::sleep(std::time::Duration::from_micros(200));
            }
            std::thread::sleep(delay);
            let slot = &slots[j % 2];
            slot.push_decompressed_for_test(gated_block(j));
            demand.notify_delivered(slot.file_id);
        }
        while progress.load(Ordering::Acquire) < total - 1 {
            std::thread::sleep(std::time::Duration::from_micros(200));
        }
        std::thread::sleep(delay);
        for slot in &slots {
            slot.set_queue_eof_for_test();
            demand.notify_delivered(slot.file_id);
        }
    });
    (events, feeder)
}

/// The state [`FeedOnDemand`]'s clones share: the two slots, the blocks still
/// to deliver (record `j` to slot `j % 2`, in key order), and the merge demand.
struct FeedShared {
    slots: Vec<Arc<fgumi_sort::SortMergeSlot>>,
    blocks: Mutex<std::collections::VecDeque<(usize, Vec<u8>)>>,
    eof_sent: std::sync::atomic::AtomicBool,
    demand: Arc<fgumi_sort::MergeDemand>,
    starved: Arc<std::sync::atomic::AtomicBool>,
}

impl FeedShared {
    /// While the merge is starved, deliver the next block — or, once every
    /// block is out, both slots' EOF — the producer half of the wake protocol.
    /// `true` iff something was delivered.
    fn feed_if_starved(&self) -> bool {
        use std::sync::atomic::Ordering;
        if !self.starved.load(Ordering::Relaxed) || self.eof_sent.load(Ordering::Acquire) {
            return false;
        }
        let mut blocks = self.blocks.lock();
        if let Some((slot_idx, block)) = blocks.pop_front() {
            let slot = &self.slots[slot_idx];
            slot.push_decompressed_for_test(block);
            self.demand.notify_delivered(slot.file_id);
            return true;
        }
        for slot in &self.slots {
            slot.set_queue_eof_for_test();
            self.demand.notify_delivered(slot.file_id);
        }
        self.eof_sent.store(true, Ordering::Release);
        true
    }
}

/// A pool (Parallel) step standing in for the spill supply: it forwards the
/// phase-2 setup events and, on a pass that finds the merge starved, delivers
/// the next prepared block. It runs only when a pool worker passes over it, so
/// a merge parked on an empty slot is fed only once a pool worker is awake.
struct FeedOnDemand {
    shared: Arc<FeedShared>,
    held: HeldSlot<Unpushed<crate::sort::protocol::SortPhase2Event>>,
    output_byte_limit: u64,
}

impl FeedOnDemand {
    /// The setup events for two open slots holding `n_blocks_per_slot`
    /// single-record blocks each, and the step that delivers them.
    fn new(
        n_blocks_per_slot: usize,
        demand: &Arc<fgumi_sort::MergeDemand>,
        output_byte_limit: u64,
    ) -> (Vec<crate::sort::protocol::SortPhase2Event>, Self) {
        let slots = two_open_slots();
        let total = 2 * n_blocks_per_slot;
        let events = announce_slots(&slots, total);
        let blocks = (0..total).map(|j| (j % 2, gated_block(j))).collect();
        let shared = Arc::new(FeedShared {
            slots,
            blocks: Mutex::new(blocks),
            eof_sent: std::sync::atomic::AtomicBool::new(false),
            demand: Arc::clone(demand),
            starved: demand.starved_signal(),
        });
        (events, Self { shared, held: HeldSlot::new(), output_byte_limit })
    }
}

impl Step for FeedOnDemand {
    type Input = crate::sort::protocol::SortPhase2Event;
    type Outputs = fgumi_pipeline_core::outputs::Single<crate::sort::protocol::SortPhase2Event>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "FeedOnDemand",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take()
            && let Err(again) = ctx.outputs.retry(unpushed)
        {
            self.held.put(again);
            return Ok(StepOutcome::Contention);
        }
        if let Some(event) = ctx.input.pop() {
            if let Err(unpushed) = ctx.outputs.push(event) {
                self.held.put(unpushed);
            }
            return Ok(StepOutcome::Progress);
        }
        if self.shared.feed_if_starved() {
            return Ok(StepOutcome::Progress);
        }
        if ctx.input.is_drained() && self.shared.eof_sent.load(std::sync::atomic::Ordering::Acquire)
        {
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }

    fn new_worker_copy(&self) -> Self {
        Self {
            shared: Arc::clone(&self.shared),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
        }
    }
}

/// A Detached pass-through that holds the setup's `AllAnnounced` back until a
/// pool worker is parked on the pool's event-count (bounded at 5 s), so the
/// merge's first stall happens while every pool worker that could feed it is
/// asleep — only the merge's `request_worker` can then wake one.
struct HoldUntilPoolParked {
    held: HeldSlot<Unpushed<crate::sort::protocol::SortPhase2Event>>,
    output_byte_limit: u64,
}

impl HoldUntilPoolParked {
    fn new(output_byte_limit: u64) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit }
    }
}

impl Step for HoldUntilPoolParked {
    type Input = crate::sort::protocol::SortPhase2Event;
    type Outputs = fgumi_pipeline_core::outputs::Single<crate::sort::protocol::SortPhase2Event>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "HoldUntilPoolParked",
            kind: StepKind::Detached,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        use crate::sort::protocol::SortPhase2Event;
        if let Some(unpushed) = self.held.take()
            && let Err(again) = ctx.outputs.retry(unpushed)
        {
            self.held.put(again);
            return Ok(StepOutcome::Contention);
        }
        let Some(event) = ctx.input.pop() else {
            return Ok(if ctx.input.is_drained() {
                StepOutcome::Finished
            } else {
                StepOutcome::NoProgress
            });
        };
        if matches!(event, SortPhase2Event::AllAnnounced { .. }) {
            let deadline = std::time::Instant::now() + std::time::Duration::from_secs(5);
            while ctx.pool.parked_workers() == 0 && std::time::Instant::now() < deadline {
                std::thread::sleep(std::time::Duration::from_millis(1));
            }
        }
        if let Err(unpushed) = ctx.outputs.push(event) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

/// A merge parked on an empty slot gets its block from pool work that reaches
/// it off any queue edge: at threads 2, the merge's driver, the sink's driver
/// and both pool workers all have 10 s park timers, and the only supply is
/// [`FeedOnDemand`], which runs only when a pool worker passes over it. The
/// sink is Detached, so no pool worker is woken by the merge's output either:
/// after the merge stalls, nothing reaches a parked worker except the merge's
/// `request_worker`. Finishing under the watchdog proves the request woke a
/// worker; `pool_requests_woken` (recorded by pipeline-core's `PoolHandle`, not
/// by the merge) confirms it, and the first stall is held until a worker is
/// parked so at least one request must wake one.
#[test]
fn stalled_merge_requests_a_worker_that_feeds_it() {
    use fgumi_pipeline_core::runtime::{TestBackoff, TestBackoffTarget};
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let (events, feeder) = FeedOnDemand::new(10, &demand, 1 << 20);
    let opts = MergeRunOptions {
        demand: Some(Arc::clone(&demand)),
        threads: 2,
        test_backoff: vec![
            TestBackoff { target: TestBackoffTarget::AllDrivers, us: 10_000_000 },
            TestBackoff { target: TestBackoffTarget::AllWorkers, us: 10_000_000 },
        ],
        feeder: Some(feeder),
        hold_until_parked: true,
        sink_kind: Some(StepKind::Detached),
        ..Default::default()
    };
    let (done_tx, done_rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = done_tx.send(collect_merge_batches_with(events, 1 << 20, 1, opts));
    });
    let run = done_rx
        .recv_timeout(std::time::Duration::from_secs(15))
        .expect("WEDGED: the starved merge's pool-worker request woke nobody")
        .unwrap();
    assert_eq!(run.positions(), (0..20).collect::<Vec<i32>>(), "every record, in key order");
    assert!(
        run.merge_stats().pool_requests_woken > 0,
        "a request must have woken a parked worker: {:?}",
        run.merge_stats()
    );
    assert!(
        !demand.starved_signal().load(std::sync::atomic::Ordering::Relaxed),
        "starved is lowered when the merge finishes"
    );
}

/// At one thread a sort chain is never fused (it has Detached steps), and the
/// pool's lone worker 0 is the only thread that runs [`FeedOnDemand`]: with
/// worker 0 and every driver on 10 s timers, the run finishes under the
/// watchdog only if the merge's request unparks worker 0.
#[test]
fn request_worker_at_one_thread_unparks_worker_zero() {
    use fgumi_pipeline_core::runtime::{TestBackoff, TestBackoffTarget};
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let (events, feeder) = FeedOnDemand::new(10, &demand, 1 << 20);
    let shared = Arc::clone(&feeder.shared);
    let progress = Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let opts = MergeRunOptions {
        demand: Some(Arc::clone(&demand)),
        progress: Some(Arc::clone(&progress)),
        threads: 1,
        test_backoff: vec![
            TestBackoff { target: TestBackoffTarget::AllDrivers, us: 10_000_000 },
            TestBackoff { target: TestBackoffTarget::Worker(0), us: 10_000_000 },
        ],
        feeder: Some(feeder),
        sink_kind: Some(StepKind::Detached),
        ..Default::default()
    };
    let (done_tx, done_rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = done_tx.send(collect_merge_batches_with(events, 1 << 20, 1, opts));
    });
    let run = done_rx
        .recv_timeout(std::time::Duration::from_secs(10))
        .unwrap_or_else(|_| {
            panic!(
                "WEDGED: the merge's request did not unpark worker 0 (sink records {}, blocks \
                 left {}, eof sent {}, starved {}, awaited {:?}, fifo {:?}, {:?})",
                progress.load(std::sync::atomic::Ordering::Relaxed),
                shared.blocks.lock().len(),
                shared.eof_sent.load(std::sync::atomic::Ordering::Relaxed),
                shared.starved.load(std::sync::atomic::Ordering::Relaxed),
                demand.awaited(),
                shared.slots.iter().map(|s| s.fifo_len()).collect::<Vec<_>>(),
                demand.snapshot()
            )
        })
        .unwrap();
    assert_eq!(run.positions(), (0..20).collect::<Vec<i32>>(), "every record, in key order");
    assert!(run.merge_stats().pool_requests_woken > 0, "{:?}", run.merge_stats());
    // Worker 0 was asleep on its timer, which no event-count census sees: a
    // request that woke it still found a sleeping worker.
    let s = demand.snapshot();
    assert!(s.pool_woken > 0 && s.requests_with_sleeper >= s.pool_woken, "{s:?}");
}

/// One pool request per stall episode, never per re-stall: the merge's driver
/// wakes on a 1 ms timer while the gated feeder waits 20 ms before each push,
/// so every episode re-stalls many times before its block lands. Requests
/// (counted by pipeline-core for the `SortMerge` step) stay at or below the
/// stall episodes (counted by the merge demand); a request per re-stall would
/// exceed them many times over.
#[test]
fn stalled_merge_requests_once_per_episode() {
    use fgumi_pipeline_core::runtime::{TestBackoff, TestBackoffTarget};
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let progress = Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let (events, feeder) =
        gated_two_slot_feed(10, &demand, &progress, std::time::Duration::from_millis(20));
    let opts = MergeRunOptions {
        demand: Some(Arc::clone(&demand)),
        threads: 1,
        test_backoff: vec![TestBackoff { target: TestBackoffTarget::AllDrivers, us: 1_000 }],
        progress: Some(Arc::clone(&progress)),
        ..Default::default()
    };
    let run = collect_merge_batches_with(events, 1 << 20, 1, opts).unwrap();
    feeder.join().unwrap();
    assert_eq!(run.positions(), (0..20).collect::<Vec<i32>>());
    let m = run.merge_stats();
    let requests = m.pool_requests_woken + m.pool_requests_all_awake + m.pool_requests_pending;
    let episodes = demand.snapshot().stall_episodes;
    assert!(episodes >= 2, "the fixture must stall more than once ({episodes})");
    assert!(requests >= 1, "a parking stall must request a worker: {m:?}");
    assert!(requests <= episodes, "at most one request per episode: {requests} > {episodes}");
}

/// The run-end stall: the merge has emitted every record and waits only on a
/// slot's EOF. That EOF's notify must let it finish without another park —
/// its driver parks for 10 s, so a merge that reached `Done` and then waited
/// out its timer before reporting `Finished` fails the 5 s watchdog.
#[test]
fn run_end_eof_wake_finishes_the_merge_without_a_timer_park() {
    use fgumi_pipeline_core::runtime::{TestBackoff, TestBackoffTarget};
    use std::sync::atomic::Ordering;
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let progress = Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let slots = two_open_slots();
    // Slot 0: record 0, then EOF. Slot 1: record 1, left open.
    slots[0].push_decompressed_for_test(gated_block(0));
    slots[0].set_queue_eof_for_test();
    slots[1].push_decompressed_for_test(gated_block(1));
    let events = announce_slots(&slots, 2);
    let feeder = {
        let (slot, demand, progress) =
            (Arc::clone(&slots[1]), Arc::clone(&demand), Arc::clone(&progress));
        std::thread::spawn(move || {
            // Both records are out; the merge is stalled (or about to stall)
            // refilling slot 1, with nothing left to emit.
            while progress.load(Ordering::Acquire) < 2 {
                std::thread::sleep(std::time::Duration::from_micros(200));
            }
            slot.set_queue_eof_for_test();
            demand.notify_delivered(slot.file_id);
        })
    };
    let opts = MergeRunOptions {
        demand: Some(Arc::clone(&demand)),
        threads: 1,
        test_backoff: vec![TestBackoff { target: TestBackoffTarget::AllDrivers, us: 10_000_000 }],
        progress: Some(Arc::clone(&progress)),
        ..Default::default()
    };
    let (done_tx, done_rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = done_tx.send(collect_merge_batches_with(events, 1 << 20, 1, opts));
    });
    let run = done_rx
        .recv_timeout(std::time::Duration::from_secs(5))
        .expect("WEDGED: the merge parked on its timer after its last wake")
        .unwrap();
    feeder.join().unwrap();
    assert_eq!(run.positions(), vec![0, 1]);
}

/// The merge's stalls end on the awaited delivery, not on the driver's timer:
/// the merge's driver parks for 10 s, so each of the ~100 stalls ends in time
/// only if the awaited delivery's `notify_delivered` unparks it — including the
/// run-end stall, whose EOF must wake the merge and let it finish without
/// another park. A missing wake fails the 10 s watchdog.
#[test]
fn stalled_merge_is_woken_by_the_awaited_delivery() {
    use fgumi_pipeline_core::runtime::{TestBackoff, TestBackoffTarget};
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let progress = Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let (events, feeder) = gated_two_slot_feed(50, &demand, &progress, std::time::Duration::ZERO);
    let opts = MergeRunOptions {
        demand: Some(Arc::clone(&demand)),
        threads: 1,
        test_backoff: vec![TestBackoff { target: TestBackoffTarget::AllDrivers, us: 10_000_000 }],
        progress: Some(Arc::clone(&progress)),
        ..Default::default()
    };
    let (done_tx, done_rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = done_tx.send(collect_merge_batches_with(events, 1 << 20, 1, opts));
    });
    let out = done_rx
        .recv_timeout(std::time::Duration::from_secs(10))
        .expect("WEDGED: a stall was not woken by its delivery")
        .unwrap();
    feeder.join().unwrap();
    assert_eq!(out.positions(), (0..100).collect::<Vec<i32>>(), "every record, in key order");
    let s = demand.snapshot();
    assert!(s.stall_episodes > 0, "the gated feed must make the merge stall: {s:?}");
}

/// Two hand-driven dispatches over two slots: slot 0 holds record 0 and slot 1
/// record 1, both left open, so the merge emits record 0 and stalls on slot 0.
/// One-record batches, so record 0 closes a full batch before the stall.
fn stalled_after_one_record(
    demand: &Arc<fgumi_sort::MergeDemand>,
) -> (
    SortMerge<RecordBatchOutput>,
    fgumi_pipeline_core::testing::StepProbe<SortMerge<RecordBatchOutput>>,
) {
    use crate::sort::protocol::SortPhase2Event;
    let step =
        SortMerge::<RecordBatchOutput>::with_target_batch_count(SortOrder::Coordinate, 1 << 20, 1)
            .with_merge_demand(Arc::clone(demand));
    let probe = fgumi_pipeline_core::testing::StepProbe::new(&step);
    for file_id in 0..2u32 {
        let slot = Arc::new(fgumi_sort::SortMergeSlot::for_test(file_id, SpillCodec::Bgzf));
        slot.push_decompressed_for_test(gated_block(file_id as usize));
        probe.push_input(SortPhase2Event::SpillReady {
            slot,
            path: std::path::PathBuf::new(),
            records_ingested_so_far: 2,
        });
    }
    probe.push_input(SortPhase2Event::AllAnnounced {
        slot_count: 2,
        memory_chunk_count: 0,
        total_records: 2,
    });
    (step, probe)
}

/// A merge that ends short of the records ingested fails, not exits clean: two
/// slots each deliver their one record and a clean EOF, but three records were
/// announced (a run that finalized early would look exactly like this).
#[rstest]
#[case::exact(2, true)]
#[case::short(3, false)]
fn a_merge_that_ends_short_of_the_ingested_records_fails(#[case] announced: u64, #[case] ok: bool) {
    use crate::sort::protocol::SortPhase2Event;
    let mut step = SortMerge::<RecordBatchOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        1 << 20,
        256,
    );
    let probe = fgumi_pipeline_core::testing::StepProbe::new(&step);
    for file_id in 0..2u32 {
        let slot = Arc::new(fgumi_sort::SortMergeSlot::for_test(file_id, SpillCodec::Bgzf));
        slot.push_decompressed_for_test(gated_block(file_id as usize));
        slot.set_queue_eof_for_test();
        probe.push_input(SortPhase2Event::SpillReady {
            slot,
            path: std::path::PathBuf::new(),
            records_ingested_so_far: 2,
        });
    }
    probe.push_input(SortPhase2Event::AllAnnounced {
        slot_count: 2,
        memory_chunk_count: 0,
        total_records: announced,
    });
    probe.close_input();
    let mut result = Ok(StepOutcome::Progress);
    for _ in 0..100 {
        result = probe.try_run(&mut step);
        if !matches!(result, Ok(StepOutcome::Progress | StepOutcome::Contention)) {
            break;
        }
    }
    match (ok, result) {
        (true, Ok(_)) => {}
        (false, Err(e)) => {
            assert!(e.to_string().contains("sort lost records: read 3 but merged 2"), "{e}");
        }
        (_, other) => panic!("announced {announced}: unexpected {other:?}"),
    }
}

/// A registration counts as a parking registration only when the merge parks:
/// a stall in a dispatch that delivered a batch runs again at once, so it
/// withdraws its registration (a delivery must not unpark a merge that is not
/// waiting); the stall after it delivers nothing, returns `Contention`, and
/// stays registered.
#[test]
fn only_a_registration_the_merge_parks_on_counts_as_parking() {
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let (mut step, probe) = stalled_after_one_record(&demand);
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
    let s = demand.snapshot();
    assert_eq!(
        (s.registrations, s.parking_registrations),
        (1, 0),
        "record 0 went out as a full batch before the stall, so the merge did not park: {s:?}"
    );
    assert_eq!(demand.awaited(), None, "the registration was withdrawn");
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Contention);
    let s = demand.snapshot();
    assert_eq!(
        (s.registrations, s.parking_registrations),
        (2, 1),
        "an empty builder: nothing delivered, so the merge parks: {s:?}"
    );
    assert_eq!(demand.awaited(), Some(0), "a parked merge stays registered on slot 0");
}

/// Two BGZF spill files whose blocks each hold one coordinate record (record
/// `j` at position `j` in slot `j % 2`, so the merge switches slots every
/// record), opened with the production opener and ingested whole into their
/// raw stashes: every block is read and parsed, none decompressed. With
/// `corrupt`, slot 1's first block carries a bad CRC.
fn stashed_two_slots(
    dir: &std::path::Path,
    per_slot: usize,
    corrupt: bool,
) -> Vec<Arc<fgumi_sort::SortMergeSlot>> {
    use std::io::Write;
    let pool = fgumi_bam_io::pread::SliceBufferPool::new(0);
    (0..2usize)
        .map(|id| {
            let path = dir.join(format!("run{id}.spill"));
            let mut file = std::fs::File::create(&path).expect("create spill file");
            file.write_all(fgumi_sort::spill_magic(SpillCodec::Bgzf)).unwrap();
            let mut compressor =
                fgumi_sort::SpillBlockCompressor::new(SpillCodec::Bgzf, 1).unwrap();
            for k in 0..per_slot {
                let mut block = compressor.compress_block(&gated_block(2 * k + id)).unwrap();
                if corrupt && id == 1 && k == 0 {
                    let crc = block.len() - 8;
                    block[crc] ^= 0xFF;
                }
                file.write_all(&block).unwrap();
            }
            file.write_all(fgumi_sort::spill_trailer(SpillCodec::Bgzf)).unwrap();
            drop(file);
            let slot = fgumi_sort::open_spill_slot(&path, u32::try_from(id).unwrap())
                .expect("open spill slot");
            let start = slot.body_start();
            let mut body = vec![0u8; usize::try_from(slot.len() - start).unwrap()];
            fgumi_bam_io::pread::read_at_exact(&**slot.source(), &mut body, start).unwrap();
            slot.bp_note_issued(body.len() as u64);
            let out = slot.bp_ingest_slice(0, pool.lease(body), true).expect("ingest");
            assert_eq!(out.frames, per_slot, "every block stashed");
            slot
        })
        .collect()
}

/// Run `SortMerge` alone (no decompress step: nothing but the merge itself can
/// decompress a stashed block) over `slots` under a 30 s watchdog.
fn merge_alone_over_stashes(
    slots: &[Arc<fgumi_sort::SortMergeSlot>],
    total: usize,
    demand: &Arc<fgumi_sort::MergeDemand>,
    ledger: &Arc<crate::sort::supply_ledger::SupplyLedger>,
) -> Result<MergeRun> {
    use fgumi_pipeline_core::runtime::{TestBackoff, TestBackoffTarget};
    let events = announce_slots(slots, total);
    let opts = MergeRunOptions {
        demand: Some(Arc::clone(demand)),
        ledger: Some(Arc::clone(ledger)),
        threads: 2,
        test_backoff: vec![TestBackoff { target: TestBackoffTarget::AllDrivers, us: 10_000_000 }],
        ..Default::default()
    };
    let (done_tx, done_rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = done_tx.send(collect_merge_batches_with(events, 1 << 20, 4, opts));
    });
    done_rx
        .recv_timeout(std::time::Duration::from_secs(30))
        .expect("WEDGED: the merge parked on a slot whose blocks it could decompress itself")
}

/// With the merge stalled on a slot whose stash holds blocks and no
/// decompress worker in the chain, the merge serves itself and finishes:
/// every record in key order, every block decompressed by the merge, and the
/// claims booked to the consumer.
#[test]
fn stalled_merge_self_serves_from_the_stash() {
    let tmp = tempfile::tempdir().unwrap();
    let slots = stashed_two_slots(tmp.path(), 40, false);
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let ledger = Arc::new(crate::sort::supply_ledger::SupplyLedger::default());
    let run = merge_alone_over_stashes(&slots, 80, &demand, &ledger).unwrap();
    assert_eq!(run.positions(), (0..80).collect::<Vec<i32>>(), "every record, in key order");
    let s = demand.snapshot();
    assert_eq!(s.self_served_blocks, 80, "{s:?}");
    assert!(s.self_served_episodes >= 1, "{s:?}");
    assert_eq!(ledger.claims(), (0, 80), "every claim by the consumer");
    // The `--sort-stats` line reports it.
    let line = s
        .log_lines()
        .into_iter()
        .find(|l| l.starts_with("Consumer served itself: "))
        .expect("the self-serve line");
    assert!(
        line.ends_with(&format!(
            "{} parks avoided by decompressing already-read blocks inline (80 blocks \
             decompressed by the merge)",
            s.self_served_episodes
        )),
        "{line}"
    );
    // A stall the merge served itself never parks: its registration is
    // withdrawn, not counted as a parking registration.
    assert!(s.registrations >= s.self_served_episodes, "{s:?}");
    assert_eq!(s.parking_registrations, 0, "no stall parked: {s:?}");
    assert_eq!(demand.awaited(), None, "every self-served registration withdrawn");
}

/// A self-serve claim releases from the ledger's stash tally exactly the
/// charge the ingest booked, so a merge served entirely by the consumer leaves
/// the tally at zero and its peak at the ingests' charge (the slots' resident
/// bytes), never above it.
#[test]
fn self_serve_returns_the_ledger_stash_to_zero() {
    let tmp = tempfile::tempdir().unwrap();
    let slots = stashed_two_slots(tmp.path(), 40, false);
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let ledger = Arc::new(crate::sort::supply_ledger::SupplyLedger::default());
    // Book each slot's ingest as `SortSpillDecompress` does.
    let mut bound = 0u64;
    for slot in &slots {
        let charged = slot.stash_bytes();
        bound += charged;
        ledger.add_stash(i64::try_from(charged).unwrap(), 40);
    }
    let run = merge_alone_over_stashes(&slots, 80, &demand, &ledger).unwrap();
    assert_eq!(run.positions(), (0..80).collect::<Vec<i32>>());
    assert_eq!(ledger.claims(), (0, 80), "every claim by the consumer");
    assert_eq!(ledger.stash_bytes(), 0, "every charge released");
    assert_eq!(ledger.stash_peak_bytes(), bound, "the peak is the booked charge");
}

/// Self-serve never reads the spill file: with every slot's file truncated to
/// zero and unlinked after its bytes were stashed, the merge still completes
/// from the stash alone.
#[test]
fn self_serve_never_reads() {
    let tmp = tempfile::tempdir().unwrap();
    let slots = stashed_two_slots(tmp.path(), 40, false);
    for id in 0..2 {
        let path = tmp.path().join(format!("run{id}.spill"));
        std::fs::OpenOptions::new().write(true).open(&path).unwrap().set_len(0).unwrap();
        std::fs::remove_file(&path).unwrap();
    }
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let ledger = Arc::new(crate::sort::supply_ledger::SupplyLedger::default());
    let run = merge_alone_over_stashes(&slots, 80, &demand, &ledger).unwrap();
    assert_eq!(run.positions(), (0..80).collect::<Vec<i32>>());
}

/// Self-serve runs the drain-only path: slot 0's last block was claimed (the
/// front, over a full FIFO) and decompressed by a worker, so it waits in
/// `reorder` with the stash empty. Once the merge pops the FIFO and stalls on
/// slot 0, no claim remains; only the drain moves the block and finalizes the
/// slot, without any worker.
#[test]
fn self_serve_drains_a_block_waiting_behind_a_full_fifo() {
    let tmp = tempfile::tempdir().unwrap();
    let slots = stashed_two_slots(tmp.path(), 2, false);
    let s0 = &slots[0];
    s0.set_fifo_cap(1);
    let mut dec = fgumi_sort::SpillBlockDecompressor::new();
    for _ in 0..2 {
        let b = s0.bp_claim_raw(u64::MAX).expect("block 0, then the front over the cap");
        b.decompress_and_publish(&mut dec).unwrap();
    }
    assert_eq!((s0.fifo_len(), s0.reorder_len_relaxed()), (1, 1), "block 1 waits in reorder");
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let ledger = Arc::new(crate::sort::supply_ledger::SupplyLedger::default());
    let run = merge_alone_over_stashes(&slots, 4, &demand, &ledger).unwrap();
    assert_eq!(run.positions(), vec![0, 1, 2, 3]);
}

/// A self-serve decompress error marks the slot failed and fails the merge
/// with the slot's error text.
#[test]
fn self_serve_decompress_error_fails_closed() {
    let tmp = tempfile::tempdir().unwrap();
    let slots = stashed_two_slots(tmp.path(), 4, true);
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let ledger = Arc::new(crate::sort::supply_ledger::SupplyLedger::default());
    let Err(err) = merge_alone_over_stashes(&slots, 8, &demand, &ledger) else {
        panic!("a corrupt stashed block must fail the merge");
    };
    let msg = format!("{err:#}");
    assert!(msg.contains("spill decompression error on slot 1"), "{msg}");
    assert!(slots[1].has_error(), "the slot is marked failed");
}

/// Like [`collect_merge_batches`] but drives the terminal `SortMerge<BlockOutput>`
/// fast path with an explicit `fast_path_threads` and a lowered parallel gate
/// (`min_records`), returning each emitted block's `(batch_serial, framed bytes)`
/// so a test can assert the parallel gather is byte-for-byte identical to the
/// serial one. `threads == 1` (or a `min_records` above the record count) takes
/// the serial gather; `threads > 1` with a low `min_records` takes the parallel
/// gather.
fn collect_fastpath_blocks(
    events: Vec<crate::sort::protocol::SortPhase2Event>,
    output_byte_limit: u64,
    target_batch_count: usize,
    fast_path_threads: usize,
    min_records: usize,
) -> Result<Vec<(u64, Vec<u8>)>> {
    use crate::sort::merge::BlockOutput;
    use crate::types::DecompressedBlock;

    let received: Arc<Mutex<Vec<DecompressedBlock>>> = Arc::new(Mutex::new(Vec::new()));
    let source = Phase2EventSource::new(events, output_byte_limit);
    let merge = SortMerge::<BlockOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        output_byte_limit,
        target_batch_count,
    )
    .with_fast_path_threads(fast_path_threads)
    .with_fast_path_min_records(min_records);
    let sink = BlockSink { received: Arc::clone(&received), kind: StepKind::Serial };

    let builder = Pipeline::builder();
    builder.chain(source).chain(merge).chain(sink).into_sink_marker();
    let pipeline = builder.build()?;
    pipeline.run(PipelineConfig { threads: 4, ..Default::default() })?;

    let blocks = std::mem::take(&mut *received.lock());
    Ok(blocks.into_iter().map(|b| (b.batch_serial, b.bytes)).collect())
}

/// The terminal fast-path parallel gather must emit byte-for-byte the same
/// blocks (same dense ordinals, same framed bytes, same boundaries) as the
/// serial gather — the parallelism is a pure throughput change, never an output
/// change. Exercised across the two block-closing regimes (count-cap and
/// byte-cap) and a boundary case where the last block is partial.
///
/// The `backpressure` case deliberately drives `emit_fast_batches_parallel`'s
/// reject/resume path: `SortMerge` wires `output_byte_limit` as BOTH the output
/// queue's `ByteBounded` limit AND the block-closing byte cap, so a small
/// `byte_limit` makes the queue fill after ~1 block while a full
/// `fast_path_threads * FAST_PATH_BLOCKS_PER_WORKER_WINDOW`-block window is
/// mid-drain — forcing `ctx.outputs.push` to reject deterministically and the
/// `held` + `fast_pending` resume to carry the remainder to the next dispatch,
/// across many windows. If the byte-identical / dense-ordinal assertions still
/// hold under sustained mid-window rejection, the resume path preserves order.
/// (`byte_capped` also trips this incidentally; `backpressure` pins it.)
#[rstest]
#[case::count_capped(2000, 300, 256, 1 << 30)] // many count-capped blocks + partial tail
#[case::byte_capped(2000, 120, 1 << 20, 8 * 1024)] // byte cap trips before the count cap
#[case::few_blocks(100, 64, 40, 1 << 30)] // small multi-block plan (40/40/20)
#[case::uneven_tail(1023, 200, 512, 1 << 30)] // last block short of the count cap
#[case::backpressure(5000, 128, 1 << 20, 4 * 1024)] // tiny queue: sustained mid-window reject/resume
fn test_fastpath_parallel_matches_serial(
    #[case] n_records: usize,
    #[case] rec_size: usize,
    #[case] target_batch_count: usize,
    #[case] byte_limit: u64,
) {
    let (_header, records) = synthesize_sized_records(n_records, 42, rec_size);
    let make_events = || {
        vec![
            coordinate_memory_chunk_event(records.clone()),
            crate::sort::protocol::SortPhase2Event::AllAnnounced {
                slot_count: 0,
                memory_chunk_count: 1,
                total_records: n_records as u64,
            },
        ]
    };

    // Serial gather: threads=1.
    let serial = collect_fastpath_blocks(make_events(), byte_limit, target_batch_count, 1, 0)
        .expect("serial fast path");
    // Parallel gather: threads=4, min_records=0 forces the parallel path.
    let parallel = collect_fastpath_blocks(make_events(), byte_limit, target_batch_count, 4, 0)
        .expect("parallel fast path");

    assert_eq!(
        serial.len(),
        parallel.len(),
        "block count differs (serial {} vs parallel {})",
        serial.len(),
        parallel.len()
    );
    for (i, (s, p)) in serial.iter().zip(parallel.iter()).enumerate() {
        assert_eq!(s.0, p.0, "block {i} ordinal differs: serial {} vs parallel {}", s.0, p.0);
        assert_eq!(s.0, i as u64, "block {i} ordinal must be dense (== index)");
        assert_eq!(
            s.1,
            p.1,
            "block {i} framed bytes differ (serial {} bytes vs parallel {} bytes)",
            s.1.len(),
            p.1.len()
        );
    }
}

/// A `Parallel` pass-through for gathered blocks that takes a permit of `cap`
/// for each block it forwards — the role `BgzfCompress` plays beside the
/// fast-path gather in a capped standalone sort.
#[derive(Clone)]
struct CappedBlockPassThrough {
    cap: Arc<fgumi_pipeline_core::PhaseCap>,
}

impl Step for CappedBlockPassThrough {
    type Input = crate::types::DecompressedBlock;
    type Outputs = fgumi_pipeline_core::outputs::Single<crate::types::DecompressedBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "CappedBlockPassThrough",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::Unbounded],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let Some(_permit) = self.cap.try_acquire() else {
            return Ok(StepOutcome::Capped);
        };
        match ctx.input.pop() {
            Some(block) => {
                assert!(ctx.outputs.push(block).is_ok(), "unbounded queue accepts");
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }

    fn step_cap(&self) -> Option<fgumi_pipeline_core::StepCapRef<'_>> {
        Some(fgumi_pipeline_core::StepCapRef::Fixed(&self.cap))
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

/// The capped parallel fast-path gather takes the WHOLE phase-2 cap for each
/// burst, completes beside a capped pool step sharing that cap, and is
/// byte-identical to the uncapped serial gather. The pool has 2 workers, so
/// pool steps alone can never hold more than 2 of the cap's 4 permits: a
/// high-water mark of 4 is reachable only through `acquire_whole`, so dropping
/// it from the gather fails this test.
#[test]
fn capped_parallel_fast_path_gather_holds_the_whole_phase2_cap() {
    use crate::sort::merge::BlockOutput;
    use fgumi_pipeline_core::PhaseCap;

    const GATHER_THREADS: usize = 4;
    const POOL_THREADS: usize = 2;
    let (_header, records) = synthesize_sized_records(4000, 0x0FA5_7CA9, 120);
    let events = || {
        vec![
            coordinate_memory_chunk_event(records.clone()),
            crate::sort::protocol::SortPhase2Event::AllAnnounced {
                slot_count: 0,
                memory_chunk_count: 1,
                total_records: records.len() as u64,
            },
        ]
    };
    let serial =
        collect_fastpath_blocks(events(), 1 << 30, 128, 1, 0).expect("uncapped serial gather");

    let cap = PhaseCap::new("sort-phase2", GATHER_THREADS);
    let received: Arc<Mutex<Vec<crate::types::DecompressedBlock>>> =
        Arc::new(Mutex::new(Vec::new()));
    let merge =
        SortMerge::<BlockOutput>::with_target_batch_count(SortOrder::Coordinate, 1 << 30, 128)
            .with_fast_path_threads(GATHER_THREADS)
            .with_fast_path_min_records(0)
            .with_fast_path_cap(Some(Arc::clone(&cap)));
    let builder = Pipeline::builder();
    builder
        .chain(Phase2EventSource::new(events(), 1 << 30))
        .chain(merge)
        .chain(CappedBlockPassThrough { cap: Arc::clone(&cap) })
        .chain(BlockSink { received: Arc::clone(&received), kind: StepKind::Serial })
        .into_sink_marker();
    builder
        .build()
        .expect("build")
        .run(PipelineConfig { threads: POOL_THREADS, ..Default::default() })
        .expect("the capped gather completes beside the capped pool step");

    let mut capped: Vec<(u64, Vec<u8>)> = std::mem::take(&mut *received.lock())
        .into_iter()
        .map(|b| (b.batch_serial, b.bytes))
        .collect();
    capped.sort_by_key(|(ordinal, _)| *ordinal);
    assert!(serial.len() > 1, "fixture must plan several blocks");
    assert_eq!(capped, serial, "the capped parallel gather is byte-identical");
    assert_eq!(cap.peak(), GATHER_THREADS, "only the whole-cap gather reaches the full cap");
    assert!(cap.peak() > POOL_THREADS);
    assert_eq!(cap.active(), 0, "every permit released");
}

/// Like [`collect_fastpath_blocks`] but for the default `SortMerge<RecordBatchOutput>`
/// (the intermediate-sort output feeding group/consensus), whose byte-cap plan
/// relies on `FRAME_OVERHEAD_PER_RECORD = 0` rather than `BlockOutput`'s `4`.
/// Returns each emitted batch's `(batch_serial, per-record bodies)` so a test
/// can assert the parallel gather is byte-for-byte identical to the serial one —
/// pinning both the dense ordinal and every block boundary (via the per-record
/// body vector, not just the concatenated bytes).
fn collect_fastpath_batches(
    events: Vec<crate::sort::protocol::SortPhase2Event>,
    output_byte_limit: u64,
    target_batch_count: usize,
    fast_path_threads: usize,
    min_records: usize,
) -> Result<Vec<(u64, Vec<Vec<u8>>)>> {
    let received: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let source = Phase2EventSource::new(events, output_byte_limit);
    let merge = SortMerge::<RecordBatchOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        output_byte_limit,
        target_batch_count,
    )
    .with_fast_path_threads(fast_path_threads)
    .with_fast_path_min_records(min_records);
    let sink = VecSink { received: Arc::clone(&received), kind: StepKind::Serial, progress: None };

    let builder = Pipeline::builder();
    builder.chain(source).chain(merge).chain(sink).into_sink_marker();
    let pipeline = builder.build()?;
    pipeline.run(PipelineConfig { threads: 4, ..Default::default() })?;

    let batches = std::mem::take(&mut *received.lock());
    Ok(batches
        .into_iter()
        .map(|b| (b.batch_serial(), b.iter_record_bytes().map(<[u8]>::to_vec).collect()))
        .collect())
}

/// Same parity contract as [`test_fastpath_parallel_matches_serial`] but for the
/// `RecordBatchOutput` builder (`FRAME_OVERHEAD_PER_RECORD = 0`), covering the
/// 0-overhead byte-cap plan the intermediate-sort branch now drives in parallel.
/// The `backpressure` case pins the reject/resume path for this builder too —
/// see the mechanism note on [`test_fastpath_parallel_matches_serial`].
#[rstest]
#[case::count_capped(2000, 300, 256, 1 << 30)] // many count-capped batches + partial tail
#[case::byte_capped(2000, 120, 1 << 20, 8 * 1024)] // byte cap trips before the count cap
#[case::few_batches(100, 64, 40, 1 << 30)] // small multi-block plan (40/40/20)
#[case::uneven_tail(1023, 200, 512, 1 << 30)] // last batch short of the count cap
#[case::backpressure(5000, 128, 1 << 20, 4 * 1024)] // tiny queue: sustained mid-window reject/resume
fn test_fastpath_parallel_matches_serial_record_batches(
    #[case] n_records: usize,
    #[case] rec_size: usize,
    #[case] target_batch_count: usize,
    #[case] byte_limit: u64,
) {
    let (_header, records) = synthesize_sized_records(n_records, 42, rec_size);
    let make_events = || {
        vec![
            coordinate_memory_chunk_event(records.clone()),
            crate::sort::protocol::SortPhase2Event::AllAnnounced {
                slot_count: 0,
                memory_chunk_count: 1,
                total_records: n_records as u64,
            },
        ]
    };

    // Serial gather: threads=1.
    let serial = collect_fastpath_batches(make_events(), byte_limit, target_batch_count, 1, 0)
        .expect("serial fast path");
    // Parallel gather: threads=4, min_records=0 forces the parallel path.
    let parallel = collect_fastpath_batches(make_events(), byte_limit, target_batch_count, 4, 0)
        .expect("parallel fast path");

    assert_eq!(
        serial.len(),
        parallel.len(),
        "batch count differs (serial {} vs parallel {})",
        serial.len(),
        parallel.len()
    );
    for (i, (s, p)) in serial.iter().zip(parallel.iter()).enumerate() {
        assert_eq!(s.0, p.0, "batch {i} ordinal differs: serial {} vs parallel {}", s.0, p.0);
        assert_eq!(s.0, i as u64, "batch {i} ordinal must be dense (== index)");
        assert_eq!(
            s.1,
            p.1,
            "batch {i} record bodies/boundaries differ (serial {} records vs parallel {})",
            s.1.len(),
            p.1.len()
        );
    }
}

/// Wrap `records` as a single coordinate-sorted in-memory chunk event. All keys
/// are `default()` (equal) — order does not matter for the buffer-sizing and
/// duplicate-announcement assertions, only that the chunk merges cleanly.
fn coordinate_memory_chunk_event(
    records: Vec<RawRecord>,
) -> crate::sort::protocol::SortPhase2Event {
    use crate::sort::protocol::MemoryChunkErased;
    let total = records.len() as u64;
    let chunk = fgumi_sort::InMemoryChunk::from_owned_records(
        records
            .into_iter()
            .map(|r| (fgumi_sort::RawCoordinateKey::default(), r.into_inner()))
            .collect(),
    );
    crate::sort::protocol::SortPhase2Event::MemoryChunk {
        chunk: Arc::new(MemoryChunkErased::Coordinate(chunk)),
        records_ingested_so_far: total,
    }
}

/// Count-bound output batches must not each reserve the full output-queue byte
/// budget. Drive a many-small-record merge that emits several count-capped
/// batches and assert their total resident capacity stays well under one byte
/// budget — it would be ~`num_batches * byte_limit` if every buffer reserved
/// the full budget (the pre-fix behavior).
#[test]
fn test_sort_merge_does_not_over_reserve_output_buffers() {
    use fgumi_pipeline_core::item::HeapSize;

    let (_header, records) = synthesize_records(600, 7);
    let byte_limit: u64 = 1 << 20; // 1 MiB
    let target = 256;
    let events = vec![
        coordinate_memory_chunk_event(records),
        crate::sort::protocol::SortPhase2Event::AllAnnounced {
            slot_count: 0,
            memory_chunk_count: 1,
            total_records: 600,
        },
    ];
    let batches = collect_merge_batches(events, byte_limit, target).expect("merge should succeed");

    let emitted: usize = batches.iter().map(|b| b.iter_record_bytes().count()).sum();
    assert_eq!(emitted, 600, "all records must be emitted");
    assert!(
        batches.len() >= 2,
        "workload must span multiple count-bound batches, got {}",
        batches.len()
    );
    let total_heap: usize = batches.iter().map(HeapSize::heap_size).sum();
    assert!(
        (total_heap as u64) < byte_limit,
        "output buffers over-reserved: {total_heap} bytes across {} batches \
         (would exceed one byte budget if each reserved the full {byte_limit})",
        batches.len(),
    );
}

/// The Phase-2 protocol emits exactly one `AllAnnounced`. A second one is a
/// protocol violation; `SortMerge` must fail closed rather than overwrite its
/// completion expectations. The first announcement over-promises (2 chunks) so
/// setup never completes and the duplicate is still absorbed in setup.
#[test]
fn test_sort_merge_fails_closed_on_duplicate_all_announced() {
    let (_header, records) = synthesize_records(1, 1);
    let byte_limit: u64 = 1 << 20;
    let events = vec![
        coordinate_memory_chunk_event(records),
        crate::sort::protocol::SortPhase2Event::AllAnnounced {
            slot_count: 0,
            memory_chunk_count: 2,
            total_records: 1,
        },
        crate::sort::protocol::SortPhase2Event::AllAnnounced {
            slot_count: 0,
            memory_chunk_count: 2,
            total_records: 1,
        },
    ];
    let err = collect_merge_batches(events, byte_limit, 256)
        .expect_err("duplicate AllAnnounced must error");
    let msg = err.to_string();
    assert!(msg.contains("duplicate AllAnnounced"), "unexpected error: {msg}");
}

/// The `Arc::try_unwrap` guard in `absorb_phase2_event` has no coverage from the
/// other tests: `coordinate_memory_chunk_event` always mints a fresh `Arc`, so
/// only the success path runs. The protocol moves memory chunks rather than
/// cloning them, so a shared `Arc` at the merge consumer means a producer kept a
/// handle — deep-cloning the record vector instead would silently double memory.
#[test]
fn test_sort_merge_fails_closed_on_shared_memory_chunk_arc() {
    use crate::sort::protocol::{MemoryChunkErased, SortPhase2Event};

    let (_header, records) = synthesize_records(8, 3);
    let chunk =
        Arc::new(MemoryChunkErased::Coordinate(fgumi_sort::InMemoryChunk::from_owned_records(
            records
                .into_iter()
                .map(|r| (fgumi_sort::RawCoordinateKey::default(), r.into_inner()))
                .collect(),
        )));

    // Two events sharing ONE Arc — the violation. `Arc::strong_count` is 2 when
    // the merge tries to take ownership of the first.
    let events = vec![
        SortPhase2Event::MemoryChunk { chunk: Arc::clone(&chunk), records_ingested_so_far: 8 },
        SortPhase2Event::MemoryChunk { chunk, records_ingested_so_far: 8 },
        SortPhase2Event::AllAnnounced { slot_count: 0, memory_chunk_count: 2, total_records: 8 },
    ];

    let err = run_merge_over_events(events).expect_err("a shared chunk Arc must fail closed");
    let msg = err.to_string();
    assert!(
        msg.contains("Arc unexpectedly shared"),
        "expected the shared-Arc guard message, got: {msg}"
    );
}

/// The buffer-sizing logic is duplicated in `next_batch` (the k-way `Merging`
/// path) and `next_fast_batch` (the single-chunk fast path).
/// `test_sort_merge_does_not_over_reserve_output_buffers` has zero slots and one
/// memory chunk, so it only ever exercises the fast path — a regression that
/// reintroduced full-budget reservation in `next_batch` would pass it. Two
/// memory chunks force `build_driver` and the real k-way merge.
#[test]
fn test_sort_merge_does_not_over_reserve_on_the_kway_path() {
    use fgumi_pipeline_core::item::HeapSize;

    let (_header, first) = synthesize_records(300, 7);
    let (_header2, second) = synthesize_records(300, 11);
    let byte_limit: u64 = 1 << 20; // 1 MiB
    let target = 256;

    let events = vec![
        coordinate_memory_chunk_event(first),
        coordinate_memory_chunk_event(second),
        crate::sort::protocol::SortPhase2Event::AllAnnounced {
            slot_count: 0,
            memory_chunk_count: 2,
            total_records: 600,
        },
    ];
    let batches = collect_merge_batches(events, byte_limit, target).expect("merge should succeed");

    let emitted: usize = batches.iter().map(|b| b.iter_record_bytes().count()).sum();
    assert_eq!(emitted, 600, "all records from both chunks must be emitted");
    assert!(
        batches.len() >= 2,
        "workload must span multiple count-bound batches, got {}",
        batches.len()
    );
    let total_heap: usize = batches.iter().map(HeapSize::heap_size).sum();
    assert!(
        (total_heap as u64) < byte_limit,
        "k-way output buffers over-reserved: {total_heap} bytes across {} batches",
        batches.len(),
    );
}

// ── Arena block-input front (ReadBlocks → InflateToArena → the split front) ──

/// `Exclusive` source draining a `Vec<BgzfBlock>` one block per `try_run`.
struct BgzfBlockSource {
    blocks: Vec<crate::types::BgzfBlock>,
    held: HeldSlot<Unpushed<crate::types::BgzfBlock>>,
    output_byte_limit: u64,
    /// `Serial` with reader affinity (worker 0) instead of `Exclusive`.
    reader: bool,
}

impl BgzfBlockSource {
    fn new(mut blocks: Vec<crate::types::BgzfBlock>, output_byte_limit: u64) -> Self {
        blocks.reverse();
        Self { blocks, held: HeldSlot::new(), output_byte_limit, reader: false }
    }

    /// Run as the chain's reader: a `Serial` step pinned to worker 0, which a
    /// split front's extraction then excludes.
    fn reader(mut self, on: bool) -> Self {
        self.reader = on;
        self
    }
}

impl Step for BgzfBlockSource {
    type Input = ();
    type Outputs = OrderedBytesSingle<crate::types::BgzfBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "BgzfBlockSource",
            kind: if self.reader { StepKind::Serial } else { StepKind::Exclusive },
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn affinity(&self) -> fgumi_pipeline_core::step::Affinity {
        if self.reader {
            fgumi_pipeline_core::step::Affinity::Reader
        } else {
            fgumi_pipeline_core::step::Affinity::None
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
        let Some(block) = self.blocks.pop() else {
            return Ok(StepOutcome::Finished);
        };
        match ctx.outputs.push(block) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }
}

/// What the arena front produced: one entry per emitted chunk (in emit order)
/// plus the terminal `AllAnnounced` counts.
#[derive(Default)]
struct ArenaFrontOutput {
    chunks: Vec<RecordBytes>,
    /// Each chunk's key kind, parallel to `chunks`.
    key_kinds: Vec<fgumi_sort::SpillKeyKind>,
    /// Each `Spill`'s `records_ingested_so_far`, in order: the spill boundaries.
    spill_ingested: Vec<u64>,
    slot_count: u32,
    total_records: u64,
}

/// Sink that copies each chunk's record bodies out and DROPS the chunk in the
/// same `try_run`. Retaining the chunks instead would pin their arena `Arc`:
/// `ReadBlocks` owns a capacity-1 arena pool, so run *k+1* cannot start until
/// run *k*'s chunk is released, and a hoarding sink wedges the pipeline.
struct ChunkEventSink {
    received: Arc<Mutex<ArenaFrontOutput>>,
}

impl Step for ChunkEventSink {
    type Input = SortChunkEvent;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ChunkEventSink",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(event) => {
                let mut out = self.received.lock();
                if let SortChunkEvent::Spill { records_ingested_so_far, .. } = &event {
                    out.spill_ingested.push(*records_ingested_so_far);
                }
                match event {
                    SortChunkEvent::Spill { chunk, .. }
                    | SortChunkEvent::Residual { chunk, .. } => {
                        out.chunks.push(
                            (0..chunk.len()).map(|i| chunk.record_bytes(i).to_vec()).collect(),
                        );
                        out.key_kinds.push(chunk.key_kind());
                    }
                    SortChunkEvent::AllAnnounced { slot_count, total_records, .. } => {
                        out.slot_count = slot_count;
                        out.total_records = total_records;
                    }
                }
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// Minimal BAM binary header (`magic + l_text=0 + n_ref` + one entry per ref).
/// Only `n_ref` matters to the front's `bam_header_len` scan and to the
/// coordinate key, so the reference names/lengths are placeholders.
fn minimal_binary_bam_header(n_ref: u32) -> Vec<u8> {
    let mut header = Vec::new();
    header.extend_from_slice(b"BAM\x01");
    header.extend_from_slice(&0u32.to_le_bytes()); // l_text = 0
    header.extend_from_slice(&n_ref.to_le_bytes());
    for _ in 0..n_ref {
        header.extend_from_slice(&2u32.to_le_bytes()); // l_name = 2
        header.extend_from_slice(b"r\0");
        header.extend_from_slice(&1_000_000u32.to_le_bytes()); // l_ref
    }
    header
}

/// Serialize `[binary header][[block_size][body]...]` and cut it into BGZF
/// blocks of at most `payload_bytes` uncompressed each. Records straddle block
/// boundaries by construction, which is what exercises the front's carry path.
fn bgzf_blocks_for(
    records: &[RawRecord],
    n_ref: u32,
    payload_bytes: usize,
) -> Vec<crate::types::BgzfBlock> {
    let mut stream = minimal_binary_bam_header(n_ref);
    for record in records {
        let bytes = record.as_ref();
        let block_size = u32::try_from(bytes.len()).expect("record fits u32");
        stream.extend_from_slice(&block_size.to_le_bytes());
        stream.extend_from_slice(bytes);
    }
    blocks_from_stream(&stream, payload_bytes)
}

/// `bgzf_blocks_for` with the stream's last `cut` bytes dropped, so the final
/// record is incomplete.
fn bgzf_blocks_truncated(
    records: &[RawRecord],
    n_ref: u32,
    payload_bytes: usize,
    cut: usize,
) -> Vec<crate::types::BgzfBlock> {
    let mut stream = minimal_binary_bam_header(n_ref);
    for record in records {
        let bytes = record.as_ref();
        stream.extend_from_slice(&u32::try_from(bytes.len()).expect("fits").to_le_bytes());
        stream.extend_from_slice(bytes);
    }
    stream.truncate(stream.len() - cut);
    blocks_from_stream(&stream, payload_bytes)
}

/// One single-block BGZF payload per `payload_bytes` of `stream`.
fn blocks_from_stream(stream: &[u8], payload_bytes: usize) -> Vec<crate::types::BgzfBlock> {
    stream
        .chunks(payload_bytes)
        .enumerate()
        .map(|(i, payload)| {
            let mut compressor = fgumi_bgzf::writer::InlineBgzfCompressor::new(1);
            compressor.write_all(payload).expect("compress payload");
            compressor.flush().expect("flush compressor");
            let mut blocks = compressor.take_blocks();
            assert_eq!(blocks.len(), 1, "payload must fit one BGZF block");
            crate::types::BgzfBlock {
                batch_serial: i as u64,
                bytes: blocks.remove(0).data.into(),
                uncompressed_size: u32::try_from(payload.len()).expect("payload fits u32"),
                index: None,
            }
        })
        .collect()
}

/// Drive `BgzfBlockSource → ReadBlocks → InflateToArena →` the coordinate
/// split front (`FindCoordinateBoundaries → ExtractCoordinateKeys →
/// SortCoordinateRun`, one sort thread) and return every emitted chunk's record bodies, chunk by chunk, plus the
/// terminal `AllAnnounced` counts.
fn drive_arena_front(
    blocks: Vec<crate::types::BgzfBlock>,
    n_ref: u32,
    memory_limit: usize,
    output_byte_limit: u64,
    threads: usize,
) -> Result<ArenaFrontOutput> {
    let received: Arc<Mutex<ArenaFrontOutput>> = Arc::new(Mutex::new(ArenaFrontOutput::default()));
    let ledger = RunLedger::new(CoordinateKeys::known_ref_width(&n_ref));

    let builder = Pipeline::builder();
    builder
        .chain(BgzfBlockSource::new(blocks, output_byte_limit))
        .chain(ReadBlocks::new(memory_limit, output_byte_limit, Arc::clone(&ledger)))
        .chain(InflateToArena::new(output_byte_limit))
        .chain(
            FindCoordinateBoundaries::new(n_ref, output_byte_limit)
                .with_run_ledger(Arc::clone(&ledger)),
        )
        .chain(ExtractCoordinateKeys::new(output_byte_limit))
        .chain(
            SortCoordinateRun::new(1, output_byte_limit)
                .with_run_ledger(Arc::clone(&ledger), memory_limit),
        )
        .chain(ChunkEventSink { received: Arc::clone(&received) })
        .into_sink_marker();
    let pipeline = builder.build()?;
    pipeline.run(PipelineConfig { threads, ..Default::default() })?;

    Ok(std::mem::take(&mut *received.lock()))
}

/// End-to-end parity for the block-input arena front, driven through the real
/// runtime rather than by poking `ingest_block`/`finalize` directly. A
/// `memory_limit` above the whole input yields exactly ONE run, so the single
/// residual chunk must be byte-identical to the reference sort.
///
/// `payload_bytes` varies where records fall relative to block boundaries: the
/// small case guarantees records straddle blocks (the front's carry path), the
/// large case puts the whole stream in one block.
#[rstest]
#[case::straddling_blocks(4096, 20)]
#[case::one_block_per_run(60_000, 2)]
fn arena_front_chain_single_run_matches_reference(
    #[case] payload_bytes: usize,
    #[case] min_blocks: usize,
) {
    const N_RECORDS: usize = 2_000;
    let (header, records) = synthesize_records(N_RECORDS, 0xA2E4_A100);
    let n_ref = u32::try_from(header.reference_sequences().len()).expect("n_ref fits u32");

    let blocks = bgzf_blocks_for(&records, n_ref, payload_bytes);
    assert!(
        blocks.len() >= min_blocks,
        "expected at least {min_blocks} blocks at {payload_bytes} B/block, got {}",
        blocks.len()
    );

    let out = drive_arena_front(blocks, n_ref, 256 * 1024 * 1024, 4 * 1024 * 1024, 4)
        .expect("arena front drives to completion");

    assert_eq!(out.slot_count, 0, "a budget-sized run never spills");
    assert_eq!(out.total_records, N_RECORDS as u64);
    assert_eq!(out.chunks.len(), 1, "one run ⇒ one residual chunk");

    let expected = reference_sort(SortOrder::Coordinate, &header, &records, None);
    assert_eq!(out.chunks[0].len(), expected.len(), "record count mismatch");
    for (i, (got, want)) in out.chunks[0].iter().zip(expected.iter()).enumerate() {
        assert_eq!(got, want, "record {i} bytes differ from the oracle");
    }
}

/// The same front under a `memory_limit` far below the input: `ReadBlocks` seals
/// several runs, so the front emits `Spill` chunks ahead of the residual. Each
/// run is independently sorted (the merge is a later stage), so the claim here
/// is that every record survives exactly once and each chunk is itself sorted.
#[test]
fn arena_front_chain_seals_multiple_runs_without_losing_records() {
    const N_RECORDS: usize = 4_000;
    let (header, records) = synthesize_records(N_RECORDS, 0x5EA1_5EA1);
    let n_ref = u32::try_from(header.reference_sequences().len()).expect("n_ref fits u32");

    let out = drive_arena_front(
        bgzf_blocks_for(&records, n_ref, 8192),
        n_ref,
        64 * 1024, // far below the input ⇒ several runs
        4 * 1024 * 1024,
        4,
    )
    .expect("arena front drives to completion");

    assert!(out.chunks.len() > 1, "a tiny memory limit must seal several runs");
    assert_eq!(out.slot_count as usize, out.chunks.len() - 1, "every run but the last spills");
    assert_eq!(out.total_records, N_RECORDS as u64);

    let mut got: Vec<Vec<u8>> = out.chunks.into_iter().flatten().collect();
    let mut want: Vec<Vec<u8>> = records.iter().map(|r| r.as_ref().to_vec()).collect();
    got.sort_unstable();
    want.sort_unstable();
    assert_eq!(got, want, "the sealed runs must carry every input record exactly once");
}

/// Drive the full BAM sort path — the benchmark path — end to end for a
/// coordinate sort: `BgzfBlockSource → ReadBlocks → InflateToArena →` the
/// production arena front `→ SpillGather → SpillBlockCompress → SpillWrite →`
/// merge supply `→ SortMerge`. Returns `(merged record bytes, run count)`.
/// This is the head the 43 GB `1kg-wgs` benchmark actually uses
/// (SAM/`SortBuffer` is the other head), so the run-former must coalesce here,
/// not only on the `SortBuffer` path. `caps` puts the phase caps on the pool
/// steps; the run sort uses two threads.
fn drive_arena_split_pipeline(
    blocks: Vec<crate::types::BgzfBlock>,
    header: &Header,
    memory_limit: usize,
    output_byte_limit: u64,
    threads: usize,
    caps: &ChainCaps,
) -> Result<(Vec<Vec<u8>>, usize)> {
    let run = drive_arena_chain(
        production_arena_front(SortOrder::Coordinate),
        buffer_config(SortOrder::Coordinate, memory_limit, 2),
        header,
        blocks,
        output_byte_limit,
        threads,
        MergeSink::Batches(StepKind::Serial),
        caps,
    )?;
    Ok((run.records, run.runs))
}

/// A binding shared phase cap on a 4-worker pool — including a cap of 1 across
/// `InflateToArena` AND `SpillBlockCompress` — still runs to completion (a
/// refused clone keeps the drain hand-off), delivers every record and is
/// byte-identical to the uncapped run. The caps' high-water marks are a sanity
/// check of the counter, not proof that each step honours it; the per-step
/// `full_cap_refuses_without_popping_*` tests pin that.
#[rstest]
#[case::phase1_capped_to_one(1, 4)]
#[case::both_capped_to_two(2, 2)]
#[case::phase2_capped_to_one(4, 1)]
fn shared_phase_caps_bound_pool_concurrency_without_changing_output(
    #[case] phase1_max: usize,
    #[case] phase2_max: usize,
) {
    use fgumi_pipeline_core::PhaseCap;
    let (header, records) = synthesize_sized_records(8_000, 0xCA9_0001, 40);
    let n_ref = u32::try_from(header.reference_sequences().len()).expect("n_ref fits u32");
    // `BgzfBlock` is not `Clone`; `bgzf_blocks_for` is deterministic, so build
    // the same block stream once per arm.
    let blocks = || bgzf_blocks_for(&records, n_ref, 4096);

    let uncapped = drive_arena_split_pipeline(
        blocks(),
        &header,
        64 * 1024,
        4 * 1024 * 1024,
        4,
        &ChainCaps::default(),
    )
    .expect("uncapped pipeline");

    let phase1 = PhaseCap::new("sort-phase1", phase1_max);
    let phase2 = PhaseCap::new("sort-phase2", phase2_max);
    let capped = drive_arena_split_pipeline(
        blocks(),
        &header,
        64 * 1024,
        4 * 1024 * 1024,
        4,
        &ChainCaps {
            phase1: Some(Arc::clone(&phase1)),
            phase2: Some(Arc::clone(&phase2)),
            ..ChainCaps::default()
        },
    )
    .expect("capped pipeline must run to completion");

    assert_eq!(capped.0.len(), records.len(), "every record delivered under the cap");
    assert_eq!(capped, uncapped, "caps must not change the output or the run count");
    assert!(capped.1 > 1, "fixture must spill so both phases do real work (runs={})", capped.1);
    assert!(
        phase1.peak() >= 1 && phase1.peak() <= phase1_max,
        "phase-1 peak {} > cap {phase1_max}",
        phase1.peak()
    );
    assert!(
        phase2.peak() >= 1 && phase2.peak() <= phase2_max,
        "phase-2 peak {} > cap {phase2_max}",
        phase2.peak()
    );
    assert_eq!(phase1.active(), 0);
    assert_eq!(phase2.active(), 0);
}

/// `new_worker_copy` must carry the cap, or only the seed clone would be
/// capped. Pinned by identity for all three steps: with the shared counter
/// exhausted, every copy's cap refuses too (an equal-looking new `PhaseCap`
/// would admit).
#[test]
fn worker_copies_share_the_phase_cap() {
    use fgumi_pipeline_core::PhaseCap;
    let cap = PhaseCap::new("sort-phase1", 2);
    let inflate = InflateToArena::new(4096).with_phase_cap(Some(Arc::clone(&cap)));
    let compress = SpillBlockCompress::new(SpillCodec::Zstd, 1, 4096)
        .with_step_cap(fgumi_pipeline_core::StepCap::from(Arc::clone(&cap)));
    let decompress = SortSpillDecompress::new(4096, &crate::sort::SpillSupply::new())
        .with_step_cap(fgumi_pipeline_core::StepCap::from(Arc::clone(&cap)));
    let copies =
        (inflate.new_worker_copy(), compress.new_worker_copy(), decompress.new_worker_copy());
    let _p1 = cap.try_acquire().expect("1 of 2");
    let _p2 = cap.try_acquire().expect("2 of 2");
    for (name, copy_cap) in [
        ("InflateToArena", copies.0.step_cap().and_then(fgumi_pipeline_core::StepCapRef::active)),
        (
            "SpillBlockCompress",
            copies.1.step_cap().and_then(fgumi_pipeline_core::StepCapRef::active),
        ),
        (
            "SortSpillDecompress",
            copies.2.step_cap().and_then(fgumi_pipeline_core::StepCapRef::active),
        ),
    ] {
        let copy_cap = copy_cap.unwrap_or_else(|| panic!("{name} copy lost its cap"));
        assert!(copy_cap.try_acquire().is_none(), "{name} copy must share the exhausted counter");
    }
}

/// The benchmark case, in miniature: a pre-sorted BAM re-sorted through the full
/// arena/BAM path under a small memory budget must coalesce every contiguous
/// sealed run into ONE spill file, and stay byte-identical to the oracle. This is
/// the arena split front (not `SortBuffer`), so it proves the
/// run-former covers the path the 43 GB `1kg-wgs` sort exercises.
#[test]
fn presorted_bam_coalesces_to_one_run_through_the_arena_split() {
    let (header, records) = synthesize_sized_records(8_000, 0xB00C_0001, 40);
    let n_ref = u32::try_from(header.reference_sequences().len()).expect("n_ref fits u32");
    let sorted = presorted_records(SortOrder::Coordinate, &header, &records);

    // Small per-block payload ⇒ many BGZF blocks; small memory ⇒ ReadBlocks seals
    // several contiguous runs the run-former must coalesce.
    let blocks = bgzf_blocks_for(&sorted, n_ref, 4096);
    let (out, runs) = drive_arena_split_pipeline(
        blocks,
        &header,
        64 * 1024,
        4 * 1024 * 1024,
        4,
        &ChainCaps::default(),
    )
    .expect("arena-split pipeline drives to completion");

    let expected = reference_sort(SortOrder::Coordinate, &header, &sorted, None);
    assert_eq!(out, expected, "arena/BAM-path coalesced output diverged from the oracle");
    assert_eq!(runs, 1, "contiguous sealed runs must coalesce into a single run on the BAM path");
}

/// The same full BAM path on a shuffled input: partial coalescing, but the merged
/// output must still match the oracle byte-for-byte.
#[test]
fn shuffled_bam_matches_the_oracle_through_the_arena_split() {
    let (header, records) = synthesize_sized_records(8_000, 0xB00C_0002, 40);
    let n_ref = u32::try_from(header.reference_sequences().len()).expect("n_ref fits u32");

    let blocks = bgzf_blocks_for(&records, n_ref, 4096);
    let (out, runs) = drive_arena_split_pipeline(
        blocks,
        &header,
        64 * 1024,
        4 * 1024 * 1024,
        4,
        &ChainCaps::default(),
    )
    .expect("arena-split pipeline drives to completion");

    let expected = reference_sort(SortOrder::Coordinate, &header, &records, None);
    assert_eq!(out, expected, "arena/BAM-path shuffled output diverged from the oracle");
    // Shuffled BAM input under a tiny per-block budget must open more than one
    // run (partial coalescing), otherwise the case no longer exercises it.
    assert!(runs >= 2, "arena/BAM shuffled input must form >=2 runs (got {runs})");
}

// ── Domain counters (T-BW2): ReadBlocks / FindCoordinateBoundaries / SpillWrite /
// SortMerge, driven together through the real arena-split chain ──────────────

/// Serial sink that sleeps briefly per received `RecordBatch` before recording
/// it — mirroring `fgumi_pipeline_core::tests::CountingSink`'s per-item
/// `thread::sleep`. Every step under test here sits upstream of this sink and
/// can burst through its own work well within the first sampler tick (see
/// `source::read_bam::tests::try_run_bumps_blocks_and_bytes_read_counters` and
/// `sink::write_bgzf::tests::try_run_bumps_bytes_written_counter` for the two
/// distinct ways that bites a naive telemetry test); throttling the terminal
/// sink keeps the run open long enough for the sampler to observe each
/// upstream step's counter at its frozen final value.
struct ThrottledVecSink {
    received: Arc<Mutex<Vec<RecordBatch>>>,
}

impl Step for ThrottledVecSink {
    type Input = RecordBatch;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ThrottledVecSink",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(batch) => {
                std::thread::sleep(std::time::Duration::from_micros(500));
                self.received.lock().push(batch);
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// Drives the full arena-split chain (`BgzfBlockSource -> ReadBlocks ->
/// InflateToArena -> FindCoordinateBoundaries -> ExtractCoordinateKeys ->
/// SortCoordinateRun -> SpillGather -> SpillBlockCompress
/// -> SpillWrite -> SortSpillDecompress -> SortMerge -> ThrottledVecSink`) with
/// telemetry enabled and asserts four of its steps' T-BW2 domain counters land
/// in the telemetry files with sane, bounded values: `ReadBlocks`
/// (`blocks`/`bytes_read`), `FindCoordinateBoundaries` (`records`), `SpillWrite`
/// (`spill_bytes_written`), and `SortMerge` (`records`). A tiny `memory_limit`
/// forces several sealed runs (several spill files), giving every one of these
/// steps repeated, spread-out `try_run` calls rather than a single batch.
#[test]
#[allow(clippy::too_many_lines)] // one end-to-end chain and its four counter checks
fn arena_split_pipeline_bumps_read_blocks_boundaries_spill_and_merge_counters() {
    use fgumi_pipeline_core::builder::{InstrumentationLevel, Pipeline, PipelineConfig};
    use fgumi_pipeline_core::runtime::telemetry::TelemetryConfig;
    use fgumi_sort::TmpDirAllocator;
    use std::time::Duration;

    const N_RECORDS: usize = 4_000;
    let (header, records) = synthesize_sized_records(N_RECORDS, 0xB00C_0003, 40);
    let n_ref = u32::try_from(header.reference_sequences().len()).expect("n_ref fits u32");

    let blocks = bgzf_blocks_for(&records, n_ref, 4096);
    let total_blocks = blocks.len() as u64;
    let total_bytes: u64 = blocks.iter().map(|b| b.bytes.len() as u64).sum();

    let memory_limit = 64 * 1024; // tiny ⇒ several sealed runs / spill files
    let output_byte_limit = 4 * 1024 * 1024;

    let dir = tempfile::TempDir::new().unwrap();
    let alloc =
        TmpDirAllocator::with_probe(vec![dir.path().to_path_buf()], Box::new(|_| Ok(u64::MAX)), 0)
            .unwrap();
    let temp_dirs = Arc::new(vec![dir]);
    let codec = SpillCodec::Zstd;

    let merge = SortMerge::<RecordBatchOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        output_byte_limit,
        256,
    );
    let received: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink = ThrottledVecSink { received: Arc::clone(&received) };

    let telemetry_dir =
        std::env::temp_dir().join(format!("fgumi-tbw2-arena-split-{}", std::process::id()));
    std::fs::create_dir_all(&telemetry_dir).unwrap();
    let stem = telemetry_dir.join("run");

    let (planner, pread, decompress, demand) =
        spill_supply(output_byte_limit, 1, SupplyTuning::default());
    let ledger = RunLedger::new(CoordinateKeys::known_ref_width(&n_ref));
    let builder = Pipeline::builder();
    let fan = builder
        .chain(BgzfBlockSource::new(blocks, output_byte_limit))
        .chain(ReadBlocks::new(memory_limit, output_byte_limit, Arc::clone(&ledger)))
        .chain(InflateToArena::new(output_byte_limit))
        .chain(
            FindCoordinateBoundaries::new(n_ref, output_byte_limit)
                .with_run_ledger(Arc::clone(&ledger)),
        )
        .chain(ExtractCoordinateKeys::new(output_byte_limit))
        .chain(
            SortCoordinateRun::new(1, output_byte_limit)
                .with_run_ledger(Arc::clone(&ledger), memory_limit),
        )
        .chain(SpillGather::new(output_byte_limit))
        .chain(SpillBlockCompress::new(codec, 3, output_byte_limit))
        .chain(SpillWrite::new(Arc::new(Mutex::new(alloc)), codec, output_byte_limit, temp_dirs))
        .chain(planner)
        .into_multi();
    fan.b0.chain(pread).chain(decompress).into_sink_marker();
    fan.b1.chain(merge.with_merge_demand(demand)).chain(sink).into_sink_marker();
    fan.b2.chain(crate::sort::NoConsolidation).into_sink_marker();
    let pipeline = builder.build().expect("pipeline builds");
    pipeline
        .run(PipelineConfig {
            threads: 1,
            instrumentation: InstrumentationLevel::Summary,
            telemetry: Some(TelemetryConfig {
                stem: stem.clone(),
                interval: Duration::from_millis(1),
            }),
            ..Default::default()
        })
        .expect("pipeline runs to completion");

    // Ground truth, independent of the sampled telemetry file.
    let collected = std::mem::take(&mut *received.lock());
    let total_out_records: usize = collected.iter().map(|b| b.iter_record_bytes().count()).sum();
    assert_eq!(total_out_records, N_RECORDS, "every record reached the sink exactly once");

    // Step indices, in declaration order: source=0, ReadBlocks=1,
    // InflateToArena=2, FindCoordinateBoundaries=3, ExtractCoordinateKeys=4,
    // SortCoordinateRun=5, SpillGather=6, SpillBlockCompress=7, SpillWrite=8,
    // SpillReadPlanner=9, PreadSpillSlices=10, SortSpillDecompress=11,
    // SortMerge=12, sink=13.
    let counters = std::fs::read_to_string(telemetry_dir.join("run.ticks.counters.tsv")).unwrap();
    let last_value_for = |step: &str, counter: &str| -> Option<u64> {
        counters
            .lines()
            .skip(1)
            .filter(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                f[2] == step && f[3] == counter
            })
            .last()
            .map(|l| l.split('\t').nth(5).unwrap().parse::<u64>().unwrap())
    };

    let read_blocks_blocks =
        last_value_for("1", "0").expect("ReadBlocks blocks counter recorded at least once");
    let read_blocks_bytes =
        last_value_for("1", "1").expect("ReadBlocks bytes_read counter recorded at least once");
    assert!(
        read_blocks_blocks > 0 && read_blocks_blocks <= total_blocks,
        "read_blocks_blocks={read_blocks_blocks}"
    );
    assert!(
        read_blocks_bytes > 0 && read_blocks_bytes <= total_bytes,
        "read_blocks_bytes={read_blocks_bytes}"
    );

    let scan_records =
        last_value_for("3", "0").expect("FindCoordinateBoundaries records counter recorded");
    assert!(scan_records > 0 && scan_records <= N_RECORDS as u64, "scan_records={scan_records}");

    let spill_bytes_written =
        last_value_for("8", "0").expect("SpillWrite spill_bytes_written counter recorded");
    assert!(spill_bytes_written > 0, "spill_bytes_written={spill_bytes_written}");

    let merge_records = last_value_for("12", "0").expect("SortMerge records counter recorded");
    assert!(
        merge_records > 0 && merge_records <= N_RECORDS as u64,
        "merge_records={merge_records}"
    );

    std::fs::remove_dir_all(&telemetry_dir).ok();
}

/// `Exclusive` source of a fixed `SortPhase1Event` list (in order).
struct Phase1EventSource {
    events: std::collections::VecDeque<crate::sort::protocol::SortPhase1Event>,
    held: HeldSlot<Unpushed<crate::sort::protocol::SortPhase1Event>>,
    output_byte_limit: u64,
}

impl Step for Phase1EventSource {
    type Input = ();
    type Outputs = fgumi_pipeline_core::outputs::Single<crate::sort::protocol::SortPhase1Event>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "Phase1EventSource",
            kind: StepKind::Exclusive,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take()
            && let Err(again) = ctx.outputs.retry(unpushed)
        {
            self.held.put(again);
            return Ok(StepOutcome::Contention);
        }
        let Some(event) = self.events.pop_front() else {
            return Ok(StepOutcome::Finished);
        };
        if let Err(unpushed) = ctx.outputs.push(event) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

/// A BGZF spill file whose records are at positions `file_id + k × stride`
/// (`per_block` records per block, `blocks` blocks), keyed by position.
fn interleaved_spill(
    dir: &std::path::Path,
    file_id: u32,
    stride: usize,
    blocks: usize,
    per_block: usize,
) -> Arc<fgumi_sort::SortMergeSlot> {
    use std::io::Write;
    let path = dir.join(format!("interleaved{file_id}.spill"));
    let mut file = std::fs::File::create(&path).unwrap();
    file.write_all(fgumi_sort::spill_magic(SpillCodec::Bgzf)).unwrap();
    let mut compressor = fgumi_sort::SpillBlockCompressor::new(SpillCodec::Bgzf, 1).unwrap();
    for b in 0..blocks {
        let mut raw = Vec::new();
        for i in 0..per_block {
            let pos = file_id as usize + (b * per_block + i) * stride;
            let p = i32::try_from(pos).unwrap();
            let record =
                make_bam_bytes(0, p, 0, format!("r{pos:06}").as_bytes(), &[], 4, -1, -1, &[]);
            let key = fgumi_sort::RawCoordinateKey { sort_key: pos as u64 };
            fgumi_sort::frame_keyed_record_into(&mut raw, &key, &record).unwrap();
        }
        file.write_all(&compressor.compress_block(&raw).unwrap()).unwrap();
    }
    file.write_all(fgumi_sort::spill_trailer(SpillCodec::Bgzf)).unwrap();
    drop(file);
    fgumi_sort::open_spill_slot(&path, file_id).unwrap()
}

/// Liveness of the planner after its input drains: once every phase event
/// has been forwarded, the planner is revisited only when a worker passes over
/// it — the claims that lower a slot's stash below its allowance push nothing
/// to any queue. Four spill files of at least four cold allowances each (fills
/// shrunk so the files stay small), two threads, every worker and driver on
/// 10 s park timers, a Detached sink (the merge's output wakes no worker): the
/// run must finish, every record in key order, within the 10 s watchdog.
#[test]
fn planner_tops_up_after_input_drains_with_long_timers() {
    use fgumi_pipeline_core::runtime::{TestBackoff, TestBackoffTarget};
    const K: usize = 4;
    const PER_BLOCK: usize = 8;
    const BLOCKS: usize = 300;
    const COLD_FILL: u64 = 1 << 10;
    let tmp = tempfile::tempdir().unwrap();
    let slots: Vec<_> = (0..K)
        .map(|f| interleaved_spill(tmp.path(), u32::try_from(f).unwrap(), K, BLOCKS, PER_BLOCK))
        .collect();
    for s in &slots {
        assert!(s.len() - s.body_start() >= 4 * 2 * COLD_FILL, "four cold allowances");
    }
    let mut events: std::collections::VecDeque<_> = slots
        .iter()
        .map(|s| crate::sort::protocol::SortPhase1Event::SpillReady {
            slot: Arc::clone(s),
            path: std::path::PathBuf::from("unused"),
            records_ingested_so_far: 0,
        })
        .collect();
    let total = K * BLOCKS * PER_BLOCK;
    events.push_back(crate::sort::protocol::SortPhase1Event::AllAnnounced {
        slot_count: u32::try_from(K).unwrap(),
        memory_chunk_count: 0,
        total_records: total as u64,
    });
    let limit = 1 << 20;
    let source = Phase1EventSource { events, held: HeldSlot::new(), output_byte_limit: limit };
    let tuning =
        SupplyTuning { fills: Some((4 * COLD_FILL, COLD_FILL)), ..SupplyTuning::default() };
    let (planner, pread, decompress, demand) = spill_supply(limit, 2, tuning);
    let merge =
        SortMerge::<RecordBatchOutput>::with_target_batch_count(SortOrder::Coordinate, limit, 64)
            .with_merge_demand(demand);
    let received: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink =
        VecSink { received: Arc::clone(&received), kind: StepKind::Detached, progress: None };
    let builder = Pipeline::builder();
    let fan = builder.chain(source).chain(planner).into_multi();
    fan.b2.chain(crate::sort::NoConsolidation).into_sink_marker();
    fan.b0.chain(pread).chain(decompress).into_sink_marker();
    fan.b1.chain(merge).chain(sink).into_sink_marker();
    let pipeline = builder.build().unwrap();
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(pipeline.run(PipelineConfig {
            threads: 2,
            test_backoff: vec![
                TestBackoff { target: TestBackoffTarget::AllDrivers, us: 10_000_000 },
                TestBackoff { target: TestBackoffTarget::AllWorkers, us: 10_000_000 },
            ],
            ..Default::default()
        }));
    });
    rx.recv_timeout(std::time::Duration::from_secs(10))
        .expect("WEDGED: the planner was not revisited after its input drained")
        .unwrap();
    let positions: Vec<i32> = received
        .lock()
        .iter()
        .flat_map(|b| b.iter_record_bytes().map(fgumi_raw_bam::pos).collect::<Vec<_>>())
        .collect();
    assert_eq!(positions, (0..i32::try_from(total).unwrap()).collect::<Vec<_>>());
}

// ── Output batching: a stall never closes a batch ───────────────────────────

/// [`collect_merge_batches_with`] for the terminal `SortMerge<BlockOutput>`:
/// each emitted block's `(batch_serial, bytes)`, the bytes being the
/// `[u32 len][record]` framing the merge writes before `BgzfCompress`.
///
/// It honors only `demand`, `threads`, `test_backoff` and `sink_kind`.
///
/// # Panics
///
/// Panics if `opts` sets a supply ledger, a sink progress counter, a feeder
/// or `hold_until_parked`, rather than silently dropping them.
fn collect_merge_blocks_with(
    events: Vec<crate::sort::protocol::SortPhase2Event>,
    output_byte_limit: u64,
    target_batch_count: usize,
    opts: MergeRunOptions,
) -> Result<Vec<(u64, Vec<u8>)>> {
    let MergeRunOptions {
        demand,
        ledger,
        threads,
        test_backoff,
        progress,
        feeder,
        hold_until_parked,
        sink_kind,
    } = opts;
    assert!(
        ledger.is_none() && progress.is_none() && feeder.is_none() && !hold_until_parked,
        "collect_merge_blocks_with supports only demand, threads, test_backoff and sink_kind"
    );
    let received: Arc<Mutex<Vec<crate::types::DecompressedBlock>>> =
        Arc::new(Mutex::new(Vec::new()));
    let source = Phase2EventSource::new(events, output_byte_limit);
    let mut merge = SortMerge::<BlockOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        output_byte_limit,
        target_batch_count,
    );
    if let Some(d) = demand {
        merge = merge.with_merge_demand(d);
    }
    let sink =
        BlockSink { received: Arc::clone(&received), kind: sink_kind.unwrap_or(StepKind::Serial) };
    let builder = Pipeline::builder();
    builder.chain(source).chain(merge).chain(sink).into_sink_marker();
    let pipeline = builder.build()?;
    pipeline.run(PipelineConfig { threads: threads.max(1), test_backoff, ..Default::default() })?;
    let blocks = std::mem::take(&mut *received.lock());
    Ok(blocks.into_iter().map(|b| (b.batch_serial, b.bytes)).collect())
}

/// A valid BAM record body of exactly `len` bytes at (tid 0, `pos`), read name
/// `{tag}{i:07}`, so a test can tell the run from byte 32 of the body (the
/// first byte of the name).
fn test_record_body(tag: u8, i: usize, pos: i32, len: usize) -> Vec<u8> {
    let name = format!("{}{i:07}", tag as char);
    let fixed = 32 + name.len() + 1;
    assert!(len >= fixed + 4, "len too small");
    // Pad with one Z tag so the body is exactly `len`.
    let aux_len = len - fixed;
    let mut aux = b"XZZ".to_vec();
    aux.extend(std::iter::repeat_n(b'x', aux_len - 4));
    aux.push(0);
    let body = make_bam_bytes(0, pos, 0, name.as_bytes(), &[], 0, -1, -1, &aux);
    assert_eq!(body.len(), len);
    body
}

/// The run a [`test_record_body`] record came from (the name's first byte).
fn run_tag(body: &[u8]) -> u8 {
    body[32]
}

/// Two coordinate-sorted runs whose positions interleave with cross-run ties:
/// run A at 0,1,2,2,3,4,4,…; run B at 0,1,1,2,3,3,4,5,5,…. Bodies carry their
/// run in the name's first byte.
fn interleaved_runs(n_per_run: usize, body_len: usize) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let a = (0..n_per_run)
        .map(|i| test_record_body(b'A', i, i32::try_from(i - i / 3).unwrap(), body_len))
        .collect();
    let b = (0..n_per_run)
        .map(|i| test_record_body(b'B', i, i32::try_from(i - (i + 1) / 3).unwrap(), body_len))
        .collect();
    (a, b)
}

/// Oracle: a stable merge by `RawCoordinateKey` (A before B on ties, then
/// position within the run) — `a` then `b`, stable-sorted by key, computed
/// without `MergeDriver`.
fn oracle_merge(a: &[Vec<u8>], b: &[Vec<u8>]) -> Vec<Vec<u8>> {
    use fgumi_sort::{RawCoordinateKey, RawSortKey};
    let mut all: Vec<&Vec<u8>> = a.iter().chain(b).collect();
    all.sort_by_cached_key(|r| RawCoordinateKey::extract_from_record(r));
    all.into_iter().cloned().collect()
}

/// The `(records, bytes)` of every batch a merge with no stall flush must
/// produce, from the record sequence alone: after pushing a record the batch
/// closes when it holds `count` records or `byte_limit` bytes (as `measure`d);
/// the end of input closes the final partial.
fn expected_batches(
    records: &[Vec<u8>],
    count: usize,
    byte_limit: u64,
    measure: impl Fn(&[u8]) -> usize,
) -> Vec<(usize, usize)> {
    let (mut out, mut n, mut bytes) = (Vec::new(), 0usize, 0usize);
    for r in records {
        n += 1;
        bytes += measure(r);
        if n >= count || bytes as u64 >= byte_limit {
            out.push((n, bytes));
            (n, bytes) = (0, 0);
        }
    }
    if n > 0 {
        out.push((n, bytes));
    }
    out
}

/// One emitted batch, from either output type.
#[derive(Debug)]
struct BatchView {
    ordinal: u64,
    records: usize,
    bytes: usize,
    record_bytes: Vec<Vec<u8>>,
}

/// Embedded-key spill framing for coordinate slots: `[u32 len][body]` each.
fn frame_block(records: &[Vec<u8>]) -> Vec<u8> {
    let mut out = Vec::new();
    for r in records {
        fgumi_sort::frame_keyed_record_into(&mut out, &fgumi_sort::RawCoordinateKey::default(), r)
            .expect("frame record");
    }
    out
}

/// Run the merge over `events` and view every emitted batch, from the terminal
/// `BlockOutput` (`terminal`) or the intermediate `RecordBatchOutput`.
fn merge_batch_views(
    events: Vec<crate::sort::protocol::SortPhase2Event>,
    byte_limit: u64,
    count: usize,
    opts: MergeRunOptions,
    terminal: bool,
) -> Result<Vec<BatchView>> {
    use fgumi_pipeline_core::Ordered as _;
    if terminal {
        collect_merge_blocks_with(events, byte_limit, count, opts).map(|v| {
            v.into_iter()
                .map(|(ordinal, framed)| {
                    let recs = unframe_block_records(&framed);
                    BatchView {
                        ordinal,
                        records: recs.len(),
                        bytes: framed.len(),
                        record_bytes: recs,
                    }
                })
                .collect()
        })
    } else {
        collect_merge_batches_with(events, byte_limit, count, opts).map(|run| {
            run.batches
                .iter()
                .map(|b| {
                    let recs: Vec<Vec<u8>> = b.iter_record_bytes().map(<[u8]>::to_vec).collect();
                    BatchView {
                        ordinal: b.ordinal(),
                        records: recs.len(),
                        bytes: recs.iter().map(Vec::len).sum(),
                        record_bytes: recs,
                    }
                })
                .collect()
        })
    }
}

/// What [`merge_with_a_stall_per_block`] observed.
struct StallRun {
    batches: Vec<BatchView>,
    records: Vec<Vec<u8>>,
    demand: Arc<fgumi_sort::MergeDemand>,
}

/// Feed two runs through two slots, one block (`block_records` records) per
/// stall. The feeder loop:
///   1. wait (≤ 30 s) until `demand.awaited() == Some(id)` for a slot that
///      still has blocks;
///   2. push that slot's next block through the slot's own API
///      (`bp_commit_read` + `bp_insert_drain_finalize`);
///   3. `demand.notify_delivered(id)` — the only call that clears `awaited`;
///   4. count one feed.
///
/// Before the very last block it sleeps `hold_last` (the merge is stalled on
/// it, holding its partial batch). An empty run is finalized at EOF up front.
/// The merge runs on a spawned thread under a 30 s watchdog; a merge that
/// stops stalling makes step 1 time out, so the fixture cannot silently stop
/// stalling, and a merge that fails or finishes early is reported from step 1
/// at once, with its result.
fn merge_with_a_stall_per_block(
    runs: (&[Vec<u8>], &[Vec<u8>]),
    block_records: usize,
    terminal: bool,
    byte_limit: u64,
    count: usize,
    hold_last: std::time::Duration,
) -> StallRun {
    use crate::sort::protocol::SortPhase2Event;
    use std::time::{Duration, Instant};
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let slots = two_open_slots();
    let blocks: [Vec<Vec<u8>>; 2] = [
        runs.0.chunks(block_records).map(frame_block).collect(),
        runs.1.chunks(block_records).map(frame_block).collect(),
    ];
    for i in 0..2 {
        if blocks[i].is_empty() {
            slots[i].bp_commit_read(0, true);
            slots[i].bp_drain_and_finalize();
        }
    }
    let total = (runs.0.len() + runs.1.len()) as u64;
    let events = vec![
        SortPhase2Event::SpillReady {
            slot: Arc::clone(&slots[0]),
            path: std::path::PathBuf::new(),
            records_ingested_so_far: runs.0.len() as u64,
        },
        SortPhase2Event::SpillReady {
            slot: Arc::clone(&slots[1]),
            path: std::path::PathBuf::new(),
            records_ingested_so_far: total,
        },
        SortPhase2Event::AllAnnounced {
            slot_count: 2,
            memory_chunk_count: 0,
            total_records: total,
        },
    ];
    let opts =
        MergeRunOptions { demand: Some(Arc::clone(&demand)), threads: 1, ..Default::default() };
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(merge_batch_views(events, byte_limit, count, opts, terminal));
    });
    let (mut next, mut seq) = ([0usize; 2], [0u64; 2]);
    let last_total = blocks[0].len() + blocks[1].len();
    let mut fed = 0usize;
    while fed < last_total {
        let deadline = Instant::now() + Duration::from_secs(30);
        let id = loop {
            if let Some(id) = demand.awaited()
                && next[id as usize] < blocks[id as usize].len()
            {
                break id as usize;
            }
            if let Ok(result) = rx.try_recv() {
                panic!(
                    "the merge ended with {fed} of {last_total} blocks fed: {:?}",
                    result.map(|views| views.len())
                );
            }
            assert!(Instant::now() < deadline, "the merge stopped stalling per block");
            std::thread::yield_now();
        };
        let last = next[id] + 1 == blocks[id].len();
        if fed + 1 == last_total && !hold_last.is_zero() {
            std::thread::sleep(hold_last);
        }
        slots[id].bp_commit_read(1, last);
        slots[id].bp_insert_drain_finalize(seq[id], vec![blocks[id][next[id]].clone()], 1);
        seq[id] += 1;
        next[id] += 1;
        fed += 1;
        demand.notify_delivered(u32::try_from(id).unwrap());
    }
    let views = rx.recv_timeout(Duration::from_secs(30)).expect("merge wedged").unwrap();
    let records = views.iter().flat_map(|v| v.record_bytes.clone()).collect();
    // Not an equality: one stall episode can span two feeds (priming stalls on
    // slot 0 then slot 1 with no record produced between them), and the
    // merge's re-check after declaring a slot awaited can find a block the
    // feeder pushed just after `awaited` flipped (a feed with no episode). The
    // feeder only feeds on `awaited() == Some(id)` and times out otherwise, so
    // it is itself the proof of a stall per block; the episode count must
    // merely be non-zero and never exceed it.
    let s = demand.snapshot();
    let episodes = usize::try_from(s.stall_episodes).unwrap();
    assert!(episodes > 0 && episodes <= fed, "stall episodes {episodes} vs feeds {fed}");
    // The merge's own count of the same thing: each feed waited for
    // `awaited() == Some(id)`, which only the merge's registration sets, and
    // the feed's `notify_delivered` cleared it, so every feed consumed a
    // registration of its own.
    let registrations = usize::try_from(s.registrations).unwrap();
    assert!(registrations >= fed, "{registrations} registrations for {fed} feeds");
    StallRun { batches: views, records, demand }
}

/// With a stall forced after every block, the emitted batches are exactly the
/// full batches plus one final partial — the boundaries computed from the
/// record sequence alone — so no stall produced a batch. Records equal the
/// oracle; ordinals are dense; the merge switched sources and, as the fixture
/// checks, really stalled once per block.
#[rstest]
#[case::terminal_block_output(true)]
#[case::intermediate_record_batches(false)]
fn stalls_never_close_a_batch(#[case] terminal: bool) {
    let (a, b) = interleaved_runs(10_000, 300);
    let count = 1024;
    // The count cap binds: 1024 × 304 B ≈ 304 KiB.
    let byte_limit = 4 << 20;
    let run = merge_with_a_stall_per_block(
        (&a, &b),
        40,
        terminal,
        byte_limit,
        count,
        std::time::Duration::ZERO,
    );
    let want = oracle_merge(&a, &b);
    assert!(run.records == want, "records differ from the oracle");
    let switches = run.records.windows(2).filter(|w| run_tag(&w[0]) != run_tag(&w[1])).count();
    assert!(switches > 1, "the fixture must interleave sources");
    let measure = |r: &[u8]| if terminal { r.len() + 4 } else { r.len() };
    let expected = expected_batches(&want, count, byte_limit, measure);
    let got: Vec<(usize, usize)> = run.batches.iter().map(|v| (v.records, v.bytes)).collect();
    assert_eq!(got, expected, "batches must be full batches plus one final partial");
    assert_eq!(run.batches.len(), 20_000_usize.div_ceil(count), "⌈records / count⌉ items");
    assert!(run.batches.windows(2).all(|w| w[1].ordinal == w[0].ordinal + 1), "dense ordinals");
}

/// The byte cap binding instead of the count cap: the same identity.
#[test]
fn stalls_never_close_a_batch_under_the_byte_cap() {
    let (a, b) = interleaved_runs(4_000, 300);
    // ≈ 215 records per batch, well under the count cap of 1024.
    let byte_limit = 64 * 1024;
    let run = merge_with_a_stall_per_block(
        (&a, &b),
        40,
        true,
        byte_limit,
        1024,
        std::time::Duration::ZERO,
    );
    let want = oracle_merge(&a, &b);
    assert!(run.records == want, "records differ from the oracle");
    let expected = expected_batches(&want, 1024, byte_limit, |r| r.len() + 4);
    assert_eq!(run.batches.iter().map(|v| (v.records, v.bytes)).collect::<Vec<_>>(), expected);
}

/// A long consumer stall does not flush: the merge holds a partial batch while
/// the feeder withholds the last block for 2 s; nothing is emitted until EOF,
/// and then exactly one batch (the final partial) carries every record.
#[rstest]
#[case::terminal(true)]
#[case::intermediate(false)]
fn a_long_stall_does_not_flush(#[case] terminal: bool) {
    // 500 records < count: one final partial.
    let (a, _) = interleaved_runs(500, 300);
    let run = merge_with_a_stall_per_block(
        (&a, &[]),
        100,
        terminal,
        4 << 20,
        1024,
        std::time::Duration::from_secs(2),
    );
    assert!(run.records == a, "records differ from the run");
    let stall_ns = run.demand.snapshot().stall_ns;
    assert!(stall_ns >= 1_000_000_000, "the merge was stalled ≥ 1 s, saw {stall_ns} ns");
    assert_eq!(run.batches.len(), 1, "no flush during the stall; Done emits the only batch");
    assert_eq!(run.batches[0].records, 500);
}

/// A Serial sink that pops only after the test sends a ticket on `release`,
/// and freely once `release` is dropped.
struct GatedSink {
    received: Arc<Mutex<Vec<RecordBatch>>>,
    release: std::sync::mpsc::Receiver<()>,
    have_ticket: bool,
    drain_freely: bool,
}

impl Step for GatedSink {
    type Input = RecordBatch;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "GatedSink",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        use std::sync::mpsc::TryRecvError;
        if !self.have_ticket {
            match self.release.try_recv() {
                Ok(()) => self.have_ticket = true,
                Err(TryRecvError::Disconnected) => {
                    self.have_ticket = true;
                    self.drain_freely = true;
                }
                Err(TryRecvError::Empty) => return Ok(StepOutcome::NoProgress),
            }
        }
        match ctx.input.pop() {
            Some(b) => {
                self.received.lock().push(b);
                if !self.drain_freely {
                    self.have_ticket = false;
                }
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// Fully populated two-slot events (no stalls): each run framed into
/// `block_records`-record blocks, all inserted, both slots finalized.
fn populated_two_slot_events(
    a: &[Vec<u8>],
    b: &[Vec<u8>],
    block_records: usize,
) -> Vec<crate::sort::protocol::SortPhase2Event> {
    let slots = two_open_slots();
    for (slot, run) in slots.iter().zip([a, b]) {
        let blocks: Vec<Vec<u8>> = run.chunks(block_records).map(frame_block).collect();
        let n = blocks.len();
        // `bp_insert_drain_finalize` drains at most `PHASE2_DECOMP_CAP` blocks
        // into the FIFO; with no decompressor in this chain the rest would
        // never drain.
        assert!(n <= fgumi_sort::PHASE2_DECOMP_CAP, "{n} blocks exceed the FIFO cap");
        slot.bp_commit_read(n, true);
        slot.bp_insert_drain_finalize(0, blocks, n);
    }
    let mut events = announce_slots(&slots, a.len() + b.len());
    if let crate::sort::protocol::SortPhase2Event::SpillReady { records_ingested_so_far, .. } =
        &mut events[0]
    {
        *records_ingested_so_far = a.len() as u64;
    }
    events
}

/// Backpressure: with a one-batch output edge and a sink released one item at
/// a time, every batch is delivered exactly once, in order, full except the
/// last, and the edge itself (its hold clock, under pipeline stats) records
/// that it refused the merge's pushes, so the held-batch path really ran. A
/// regression guard: the slots are fully populated, so the merge never stalls;
/// it catches a batch lost or duplicated at a full output edge.
#[test]
fn full_batches_under_backpressure_are_delivered_once() {
    use fgumi_pipeline_core::Ordered as _;
    use std::time::Duration;
    let (a, b) = interleaved_runs(3_000, 300);
    let want = oracle_merge(&a, &b);
    let received = Arc::new(Mutex::new(Vec::<RecordBatch>::new()));
    let (release_tx, release_rx) = std::sync::mpsc::channel::<()>();
    let sink = GatedSink {
        received: Arc::clone(&received),
        release: release_rx,
        have_ticket: false,
        drain_freely: false,
    };
    // 30 blocks per slot, within the FIFO cap; the merge's output edge holds
    // one full batch's bytes.
    let events = populated_two_slot_events(&a, &b, 100);
    let source = Phase2EventSource::new(events, 1 << 20);
    let merge = SortMerge::<RecordBatchOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        1024 * 300,
        1024,
    );
    let builder = Pipeline::builder();
    builder.chain(source).chain(merge).chain(sink).into_sink_marker();
    let pipeline = builder.build().unwrap();
    let stats = pipeline.stats();
    let run_stats = Arc::clone(&stats);
    let runner = std::thread::spawn(move || {
        pipeline.run(PipelineConfig { threads: 2, stats: Some(run_stats), ..Default::default() })
    });
    // `held_retries` counts each rejected retry while a batch is held, so it
    // moves before the release (`holds` counts only once the push succeeds).
    let merge_held_retries = || {
        stats
            .snapshot()
            .steps
            .iter()
            .find(|(name, _)| *name == "SortMerge")
            .map_or(0, |s| s.1.held_retries)
    };
    for round in 0..6_000_usize.div_ceil(1024) {
        // The first release waits until the merge has filled the edge and a
        // held retry was refused, so the run is provably backpressured before
        // the sink drains anything. Later rounds only pace the sink: near the
        // end the merge may have nothing left to push, so a refusal is not
        // owed there (the final `holds > 0` assertion still covers the run).
        let before = merge_held_retries();
        let patience = if round == 0 { Duration::from_secs(10) } else { Duration::from_millis(50) };
        let deadline = std::time::Instant::now() + patience;
        while merge_held_retries() == before {
            if std::time::Instant::now() >= deadline {
                assert!(round > 0, "the output edge never refused a push");
                break;
            }
            std::thread::sleep(Duration::from_millis(1));
        }
        release_tx.send(()).unwrap();
    }
    drop(release_tx);
    let (done_tx, done_rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || done_tx.send(runner.join()).unwrap());
    done_rx
        .recv_timeout(Duration::from_secs(30))
        .expect("merge wedged under backpressure")
        .unwrap()
        .unwrap();
    let batches = std::mem::take(&mut *received.lock());
    let got: Vec<Vec<u8>> =
        batches.iter().flat_map(|b| b.iter_record_bytes().map(<[u8]>::to_vec)).collect();
    assert!(got == want, "records differ from the oracle");
    assert!(batches.windows(2).all(|w| w[1].ordinal() == w[0].ordinal() + 1), "dense ordinals");
    let expected = expected_batches(&want, 1024, 1024 * 300, <[u8]>::len);
    let shape: Vec<(usize, usize)> = batches
        .iter()
        .map(|b| (b.len(), b.iter_record_bytes().map(<[u8]>::len).sum::<usize>()))
        .collect();
    assert_eq!(shape, expected);
    let snapshot = stats.snapshot();
    let merge_stats =
        &snapshot.steps.iter().find(|(name, _)| *name == "SortMerge").expect("SortMerge").1;
    assert!(merge_stats.holds > 0, "the output edge never refused a push: {merge_stats:?}");
}

// ── The split arena front (template-coordinate) ─────────────────────────────

/// A template-record fixture's lane content.
#[derive(Clone, Copy, Debug)]
enum LaneCase {
    AutoPlain,
    AutoCb,
    AutoMi,
    AutoBoth,
    Full,
    NonePlain,
    ExplicitCb,
    ExplicitTertiary,
}

impl LaneCase {
    /// `(key types, cell tag, records carry CB, records carry MI, expected lane)`.
    fn params(
        self,
    ) -> (KeyTypesSpec, Option<fgumi_raw_bam::SamTag>, bool, bool, fgumi_sort::SpillKeyKind) {
        use fgumi_raw_bam::SamTag;
        use fgumi_sort::SpillKeyKind::{TemplateCb32, TemplateK24, TemplateK40, TemplateTert32};
        match self {
            Self::AutoPlain => (KeyTypesSpec::Auto, None, false, false, TemplateK24),
            Self::AutoCb => (KeyTypesSpec::Auto, Some(SamTag::CB), true, false, TemplateCb32),
            Self::AutoMi => (KeyTypesSpec::Auto, None, false, true, TemplateTert32),
            Self::AutoBoth => (KeyTypesSpec::Auto, Some(SamTag::CB), true, true, TemplateK40),
            Self::Full => (KeyTypesSpec::Full, Some(SamTag::CB), true, true, TemplateK40),
            Self::NonePlain => (KeyTypesSpec::None, None, false, false, TemplateK24),
            Self::ExplicitCb => (
                KeyTypesSpec::Explicit { cb: true, tertiary: false },
                Some(SamTag::CB),
                true,
                false,
                TemplateCb32,
            ),
            Self::ExplicitTertiary => (
                KeyTypesSpec::Explicit { cb: false, tertiary: true },
                None,
                false,
                true,
                TemplateTert32,
            ),
        }
    }
}

/// `n` shuffled paired template records, optionally carrying `CB:Z` and
/// `MI:Z` values that vary across records.
fn template_records(n: usize, seed: u64, cb: bool, mi: bool) -> (Header, Vec<RawRecord>) {
    let mut rng = fixture_rng(seed);
    let records = (0..n)
        .map(|i| {
            #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
            let (tid, pos) = (
                rng.random_range(0..N_TEST_REFS as u64) as i32,
                rng.random_range(0..500_000) as i32,
            );
            let mut aux = Vec::new();
            if cb {
                aux.extend_from_slice(b"CBZ");
                aux.extend_from_slice(["AAAC", "CCAG", "GGTA"][i % 3].as_bytes());
                aux.push(0);
            }
            if mi {
                aux.extend_from_slice(b"MIZ");
                aux.extend_from_slice(format!("{}", rng.random_range(0..40)).as_bytes());
                aux.push(0);
            }
            let flags = if i % 2 == 0 { 0x63 } else { 0x93 };
            let name = format!("tpl{}", rng.random_range(0..1 << 20));
            RawRecord::from(make_bam_bytes(
                tid,
                pos,
                flags,
                name.as_bytes(),
                &[],
                30,
                tid,
                pos + 200,
                &aux,
            ))
        })
        .collect();
    (synthetic_header(), records)
}

/// The split front at every lane width and thread count seals every run in the
/// lane the configuration selects, and the merged output is the reference
/// stable sort. Small batches (64 records) over several runs, with records
/// straddling block and run boundaries.
#[rstest]
fn template_split_matches_reference_across_runs(
    #[values(
        LaneCase::AutoPlain,
        LaneCase::AutoCb,
        LaneCase::AutoMi,
        LaneCase::AutoBoth,
        LaneCase::Full,
        LaneCase::NonePlain,
        LaneCase::ExplicitCb,
        LaneCase::ExplicitTertiary
    )]
    case: LaneCase,
    #[values(1usize, 2, 4)] threads: usize,
) {
    let (key_types, cell_tag, cb, mi, lane) = case.params();
    let (header, records) = template_records(3_000, 0x7E4A + threads as u64, cb, mi);
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let config = SortBufferConfig {
        sort_order: SortOrder::TemplateCoordinate,
        memory_limit: 96 * 1024,
        sort_threads: threads,
        cell_tag,
        key_types,
    };
    let caps = ChainCaps { batch_records: Some(64), ..ChainCaps::default() };
    let blocks = || bgzf_blocks_for(&records, n_ref, 4096);
    let split =
        drive_arena_chunks(ArenaFront::Split, config, &header, blocks(), threads, &caps).unwrap();
    assert!(split.chunks.len() >= 3, "{case:?}: several runs (got {})", split.chunks.len());
    let spilled = &split.chunks[..split.slot_count as usize];
    assert!(
        spilled.iter().all(|c| c.len() > 64),
        "{case:?}: every spilled run spans several batches"
    );
    assert!(split.key_kinds.iter().all(|&k| k == lane), "{case:?}: {:?}", split.key_kinds);
    assert_eq!(split.total_records, records.len() as u64, "{case:?}");
    let pieces = if split.slot_count > 0 && threads > 1 { threads } else { 1 };
    assert_eq!(residual_count(&split), pieces, "{case:?}: the final run's pieces");
    let merged = drive_arena_chain(
        ArenaFront::Split,
        config,
        &header,
        blocks(),
        4 * 1024 * 1024,
        threads,
        MergeSink::Batches(StepKind::Serial),
        &caps,
    )
    .unwrap();
    assert_record_parity(
        &format!("{case:?} t{threads}"),
        &merged.records,
        &reference_sort(SortOrder::TemplateCoordinate, &header, &records, cell_tag),
    );
}

/// A dropped-lane violation under `--key-types none` is reported with the same
/// text by both template fronts — the record-input front and the split front
/// at one and four threads — and names the first offending record (7), though
/// another (300) sits in a later batch.
#[test]
fn violation_report_is_identical_across_template_paths() {
    let records: Vec<RawRecord> = (0..600i32)
        .map(|i| {
            let mi = match i {
                7 => 2,
                300 => 3,
                _ => 1,
            };
            RawRecord::from(record_with_mi(1000 + i, format!("rd_{i}").as_bytes(), mi))
        })
        .collect();
    let header = synthetic_header();
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let config = SortBufferConfig {
        key_types: KeyTypesSpec::None,
        ..buffer_config(SortOrder::TemplateCoordinate, 64 * 1024, 2)
    };
    let caps = ChainCaps { batch_records: Some(64), ..ChainCaps::default() };
    let arena_err = |front, threads| {
        let err = drive_arena_chain(
            front,
            config,
            &header,
            bgzf_blocks_for(&records, n_ref, 4096),
            4 * 1024 * 1024,
            threads,
            MergeSink::Batches(StepKind::Serial),
            &caps,
        )
        .expect_err("a violation fails the sort");
        format!("{err:#}")
    };
    // The innermost violation: the record-input front wraps it in its own
    // push context, which is not part of the report.
    let violation = |msg: &str| {
        msg.rsplit_once("record rd_").map(|(_, rest)| format!("record rd_{rest}")).expect(msg)
    };
    let split = violation(&arena_err(ArenaFront::Split, 1));
    assert!(split.starts_with("record rd_7 "), "{split}");
    assert_eq!(violation(&arena_err(ArenaFront::Split, 4)), split, "t4");
    let buffer = drive_buffer_chain(
        config,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        2,
        MergeSink::Batches(StepKind::Serial),
        SupplyTuning::default(),
        ZSTD_1,
    )
    .expect_err("the record-input front fails too");
    assert_eq!(violation(&format!("{buffer:#}")), split);
}

/// Edge inputs on the split front: a header-only BAM seals no records, a
/// single in-memory run is the reference sort, and a final record cut short
/// is the truncation error.
#[rstest]
#[case::header_only(0, usize::MAX, false)]
#[case::residual_only(500, 64 << 20, false)]
#[case::truncated(500, 64 * 1024, true)]
fn split_front_edge_cases(
    #[values(
        SortOrder::Coordinate,
        SortOrder::TemplateCoordinate,
        SortOrder::Queryname(QuerynameComparator::Lexicographic),
        SortOrder::Queryname(QuerynameComparator::Natural)
    )]
    order: SortOrder,
    #[case] n: usize,
    #[case] memory_limit: usize,
    #[case] truncate: bool,
) {
    let (header, records) = template_records(n, 0xED6E, false, false);
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let config = buffer_config(order, memory_limit.min(64 << 20), 2);
    let blocks = || {
        if truncate {
            bgzf_blocks_truncated(&records, n_ref, 4096, 5)
        } else {
            bgzf_blocks_for(&records, n_ref, 4096)
        }
    };
    let caps = ChainCaps { batch_records: Some(64), ..ChainCaps::default() };
    match drive_arena_chunks(ArenaFront::Split, config, &header, blocks(), 2, &caps) {
        Ok(out) => {
            assert!(!truncate);
            assert_eq!(out.slot_count, 0, "{order:?}: nothing spills");
            assert!(!out.key_kinds.is_empty(), "{order:?}: the residual is announced");
            assert!(
                out.key_kinds.iter().all(|&k| key_kind_is_of(order, k)),
                "{order:?}: every chunk is the order's kind: {:?}",
                out.key_kinds
            );
            assert_eq!(out.total_records, n as u64, "{order:?}");
            let sealed: Vec<Vec<u8>> = out.chunks.into_iter().flatten().collect();
            assert!(
                sealed == reference_sort(order, &header, &records, None),
                "{order:?}: the in-memory run is the reference sort"
            );
        }
        Err(err) => {
            assert!(truncate, "{order:?}: {err:#}");
            let err = format!("{err:#}");
            assert!(err.contains("truncated BAM: final record incomplete"), "{order:?}: {err}");
        }
    }
}

/// Under a phase-1 cap of one with four workers — inflate, extraction and
/// spill compression share one permit and the run sort reserves it — the
/// split front still completes, with the reference output.
#[rstest]
fn a_one_permit_phase1_cap_completes_on_the_split_front(
    #[values(
        SortOrder::Coordinate,
        SortOrder::TemplateCoordinate,
        SortOrder::Queryname(QuerynameComparator::Lexicographic),
        SortOrder::Queryname(QuerynameComparator::Natural)
    )]
    order: SortOrder,
) {
    let (header, records) = tie_heavy_records(4_000, 21);
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let cap = fgumi_pipeline_core::PhaseCap::new("sort-phase1", 1);
    let caps = ChainCaps { phase1: Some(Arc::clone(&cap)), ..ChainCaps::default() };
    let run = drive_arena_chain(
        ArenaFront::Split,
        buffer_config(order, 96 * 1024, 4),
        &header,
        bgzf_blocks_for(&records, n_ref, 4096),
        4 * 1024 * 1024,
        4,
        MergeSink::Batches(StepKind::Serial),
        &caps,
    )
    .unwrap();
    assert!(run.merge_sources >= 2, "the fixture spills");
    assert_record_parity("cap 1", &run.records, &reference_sort(order, &header, &records, None));
    assert_eq!(cap.peak(), 1);
}

/// Extraction overlaps the scan: the boundary step holds after cutting the
/// sort's first batch until that batch is extracted, without blocking its
/// driver. Were extraction serialised behind the scan, the batch would never
/// be extracted and the 30 s watchdog would fire. The census reports
/// overlapped records.
#[rstest]
fn split_extracts_while_the_scan_continues(
    #[values(
        SortOrder::Coordinate,
        SortOrder::TemplateCoordinate,
        SortOrder::Queryname(QuerynameComparator::Lexicographic),
        SortOrder::Queryname(QuerynameComparator::Natural)
    )]
    order: SortOrder,
    #[values(2usize, 4)] threads: usize,
) {
    let (header, records) = tie_heavy_records(3_000, 31);
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let probe = Arc::new(crate::sort::arena_split::SplitProbe::holding_scan());
    let caps = ChainCaps {
        batch_records: Some(64),
        probe: Some(Arc::clone(&probe)),
        ..ChainCaps::default()
    };
    let (tx, rx) = std::sync::mpsc::channel();
    let worker = std::thread::spawn(move || {
        let out = drive_arena_chain(
            ArenaFront::Split,
            buffer_config(order, 96 * 1024, threads),
            &header,
            bgzf_blocks_for(&records, n_ref, 4096),
            4 * 1024 * 1024,
            threads,
            MergeSink::Batches(StepKind::Serial),
            &caps,
        )
        .map(|r| r.records.len());
        let _ = tx.send(out);
    });
    let n = rx
        .recv_timeout(std::time::Duration::from_secs(30))
        .expect("watchdog: extraction never ran while the scan was held")
        .unwrap();
    worker.join().unwrap();
    assert!(n > 0);
    let census = probe.census.lock().expect("the run sort published its census");
    assert!(census.overlapped_records > 0, "{census:?}");
}

/// At one thread the reader holds the only pool worker, so extraction runs on
/// the `sort-coord` driver — the thread that did the work before the split.
#[rstest]
fn t1_split_runs_extract_on_the_coord_driver(
    #[values(
        SortOrder::Coordinate,
        SortOrder::TemplateCoordinate,
        SortOrder::Queryname(QuerynameComparator::Lexicographic),
        SortOrder::Queryname(QuerynameComparator::Natural)
    )]
    order: SortOrder,
) {
    let (header, records) = tie_heavy_records(2_000, 41);
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let probe = Arc::new(crate::sort::arena_split::SplitProbe::default());
    let caps =
        ChainCaps { probe: Some(Arc::clone(&probe)), reader_source: true, ..ChainCaps::default() };
    let run = drive_arena_chain(
        ArenaFront::Split,
        buffer_config(order, 96 * 1024, 1),
        &header,
        bgzf_blocks_for(&records, n_ref, 4096),
        4 * 1024 * 1024,
        1,
        MergeSink::Batches(StepKind::Serial),
        &caps,
    )
    .unwrap();
    assert_eq!(run.records.len(), records.len());
    let threads = probe.extract_threads.lock().clone();
    assert_eq!(
        threads.into_iter().collect::<Vec<_>>(),
        vec![format!("fgumi-driver-{}", crate::sort::SORT_COORD_GROUP)]
    );
}

/// The split front completes with the reference output at one and two
/// threads, uncapped and under an explicit `--sort-threads 1` (whose run sort
/// reserves the whole one-permit cap on the same driver that may host the
/// extraction).
#[rstest]
fn split_front_at_one_and_two_threads(
    #[values(
        SortOrder::Coordinate,
        SortOrder::TemplateCoordinate,
        SortOrder::Queryname(QuerynameComparator::Lexicographic),
        SortOrder::Queryname(QuerynameComparator::Natural)
    )]
    order: SortOrder,
    #[values(1usize, 2)] threads: usize,
    #[values(false, true)] capped: bool,
) {
    let (header, records) = tie_heavy_records(3_000, 51);
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let caps = ChainCaps {
        phase1: capped.then(|| fgumi_pipeline_core::PhaseCap::new("sort-phase1", 1)),
        reader_source: true,
        ..ChainCaps::default()
    };
    let run = drive_arena_chain(
        ArenaFront::Split,
        buffer_config(order, 96 * 1024, if capped { 1 } else { threads }),
        &header,
        bgzf_blocks_for(&records, n_ref, 4096),
        4 * 1024 * 1024,
        threads,
        MergeSink::Batches(StepKind::Serial),
        &caps,
    )
    .unwrap();
    assert!(run.merge_sources >= 2, "the fixture spills");
    assert_record_parity(
        &format!("t{threads} capped={capped}"),
        &run.records,
        &reference_sort(order, &header, &records, None),
    );
}

/// A large spilling template sort through the split front completes under a
/// watchdog at one, two and four pipeline threads (soak).
#[cfg(feature = "stress-tests")]
#[rstest]
fn three_step_chain_large_spill_completes_on_the_arena_split(
    #[values(1usize, 2, 4)] pipeline_threads: usize,
) {
    let (header, records) = template_records(60_000, 0x50A4, true, true);
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let expected = reference_sort(SortOrder::TemplateCoordinate, &header, &records, None);
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let out = drive_arena_chain(
            ArenaFront::Split,
            buffer_config(SortOrder::TemplateCoordinate, 1 << 20, pipeline_threads),
            &header,
            bgzf_blocks_for(&records, n_ref, 16 * 1024),
            4 * 1024 * 1024,
            pipeline_threads,
            MergeSink::Batches(StepKind::Serial),
            &ChainCaps { reader_source: true, ..ChainCaps::default() },
        );
        let _ = tx.send(out.map(|r| (r.merge_sources, r.records)));
    });
    let (merge_sources, out) = rx
        .recv_timeout(std::time::Duration::from_secs(120))
        .expect("watchdog: the split front did not finish")
        .unwrap();
    assert!(merge_sources >= 4, "the soak spills several runs (got {merge_sources})");
    assert_record_parity(&format!("t{pipeline_threads}"), &out, &expected);
}

/// Queryname records whose instrument prefix changes every 700 records of
/// input (so key batches re-learn their class mid-run), shuffled within each
/// prefix block, with read pairs, secondary copies sharing a name and flags
/// (ties), and digit runs of varying width (natural order). Each record
/// carries a unique `XI` tag, so tie order is observable.
fn instrument_switching_records(n: usize, seed: u64) -> (Header, Vec<RawRecord>) {
    const PREFIXES: [&str; 3] = ["A00132:53:HFHJ", "M04567:12:000000000-ABCDE", "NB551:7:HXY"];
    let mut rng = fixture_rng(seed);
    let mut specs: Vec<(u16, Vec<u8>)> = Vec::with_capacity(n + 8);
    let mut block = Vec::new();
    while specs.len() < n {
        let prefix = PREFIXES[(specs.len() + block.len()) / 700 % PREFIXES.len()];
        let width = 1 + u32::try_from(rng.random_range(0..5)).unwrap();
        let name = format!(
            "{prefix}:{}:{}:{}:{}",
            1 + rng.random_range(0..4),
            1101 + rng.random_range(0..30),
            rng.random_range(0..10u64.pow(width)),
            rng.random_range(0..40_000)
        )
        .into_bytes();
        block.push((0x41u16, name.clone()));
        block.push((0x81u16, name.clone()));
        if chance(&mut rng, 4) {
            block.push((0x141u16, name));
            block.push((0x141u16, block.last().unwrap().1.clone()));
        }
        if block.len() >= 700 {
            for i in (1..block.len()).rev() {
                #[allow(clippy::cast_possible_truncation)]
                let j = rng.random_range(0..i as u64 + 1) as usize;
                block.swap(i, j);
            }
            specs.append(&mut block);
        }
    }
    specs.append(&mut block);
    let records = specs
        .into_iter()
        .enumerate()
        .map(|(i, (flags, name))| {
            let mut aux = b"XIi".to_vec();
            aux.extend_from_slice(&u32::try_from(i).expect("index fits u32").to_le_bytes());
            RawRecord::from(make_bam_bytes(-1, -1, flags | 0x4, &name, &[], 40, -1, -1, &aux))
        })
        .collect();
    (synthetic_header(), records)
}

/// The queryname split front, lexicographic and natural, at one, two and four
/// threads, over 64-record key batches that re-learn their prefix class as the
/// instrument changes: the merged output is the reference stable sort, ties in
/// input order.
#[rstest]
fn queryname_split_matches_reference_across_runs(
    #[values(QuerynameComparator::Lexicographic, QuerynameComparator::Natural)]
    cmp: QuerynameComparator,
    #[values(1usize, 2, 4)] threads: usize,
) {
    let order = SortOrder::Queryname(cmp);
    let (header, records) = instrument_switching_records(4_000, 0x9A7E + threads as u64);
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let config = buffer_config(order, 96 * 1024, threads);
    let caps = ChainCaps { batch_records: Some(64), ..ChainCaps::default() };
    let blocks = || bgzf_blocks_for(&records, n_ref, 4096);
    let split =
        drive_arena_chunks(ArenaFront::Split, config, &header, blocks(), threads, &caps).unwrap();
    assert!(split.chunks.len() >= 3, "{cmp:?}: several runs (got {})", split.chunks.len());
    let kind = match cmp {
        QuerynameComparator::Lexicographic => fgumi_sort::SpillKeyKind::QuerynameLex,
        QuerynameComparator::Natural => fgumi_sort::SpillKeyKind::QuerynameNatural,
    };
    assert!(split.key_kinds.iter().all(|&k| k == kind), "{cmp:?}: {:?}", split.key_kinds);
    assert_eq!(split.total_records, records.len() as u64);
    let merged = drive_arena_chain(
        ArenaFront::Split,
        config,
        &header,
        blocks(),
        4 * 1024 * 1024,
        threads,
        MergeSink::Batches(StepKind::Serial),
        &caps,
    )
    .unwrap();
    let expected = reference_sort(order, &header, &records, None);
    assert_record_parity(&format!("{cmp:?} t{threads}"), &merged.records, &expected);
    assert_input_order_within_equal_keys(order, &header, &merged.records);
}

/// The coordinate split front at one, two and four threads, over 64-record key
/// batches across several runs: every run seals a coordinate chunk, and the
/// merged output is the reference stable sort, ties in input order.
#[rstest]
fn coordinate_split_matches_reference(#[values(1usize, 2, 4)] threads: usize) {
    let order = SortOrder::Coordinate;
    let (header, records) = tie_heavy_records(4_000, 0xC0 + threads as u64);
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let config = buffer_config(order, 96 * 1024, threads);
    let caps = ChainCaps { batch_records: Some(64), ..ChainCaps::default() };
    let blocks = || bgzf_blocks_for(&records, n_ref, 4096);
    let split =
        drive_arena_chunks(ArenaFront::Split, config, &header, blocks(), threads, &caps).unwrap();
    assert!(split.chunks.len() >= 3, "several runs (got {})", split.chunks.len());
    assert!(split.key_kinds.iter().all(|&k| k == fgumi_sort::SpillKeyKind::Coordinate));
    assert_eq!(split.total_records, records.len() as u64);
    let merged = drive_arena_chain(
        ArenaFront::Split,
        config,
        &header,
        blocks(),
        4 * 1024 * 1024,
        threads,
        MergeSink::Batches(StepKind::Serial),
        &caps,
    )
    .unwrap();
    let expected = reference_sort(order, &header, &records, None);
    assert_record_parity(&format!("t{threads}"), &merged.records, &expected);
    assert_input_order_within_equal_keys(order, &header, &merged.records);
}

/// Whether a sealed chunk's key kind belongs to `order` (any lane for template
/// order).
fn key_kind_is_of(order: SortOrder, kind: fgumi_sort::SpillKeyKind) -> bool {
    use fgumi_sort::SpillKeyKind as K;
    match order {
        SortOrder::Coordinate => kind == K::Coordinate,
        SortOrder::Queryname(QuerynameComparator::Lexicographic) => kind == K::QuerynameLex,
        SortOrder::Queryname(QuerynameComparator::Natural) => kind == K::QuerynameNatural,
        SortOrder::TemplateCoordinate => {
            matches!(kind, K::TemplateK24 | K::TemplateCb32 | K::TemplateTert32 | K::TemplateK40)
        }
    }
}

/// One BGZF block per `per_block` whole records, after a block holding only the
/// BAM header: no record straddles a block, so each block's uncompressed size
/// and record count are known. Returns the blocks and `(isize, records)` per
/// block.
fn whole_record_blocks(
    records: &[RawRecord],
    n_ref: u32,
    per_block: usize,
) -> (Vec<crate::types::BgzfBlock>, Vec<(u64, u64)>) {
    let mut payloads = vec![(minimal_binary_bam_header(n_ref), 0u64)];
    for group in records.chunks(per_block) {
        let mut payload = Vec::new();
        for record in group {
            let bytes = record.as_ref();
            payload.extend_from_slice(&u32::try_from(bytes.len()).unwrap().to_le_bytes());
            payload.extend_from_slice(bytes);
        }
        payloads.push((payload, group.len() as u64));
    }
    let shape = payloads.iter().map(|(p, n)| (p.len() as u64, *n)).collect();
    let blocks = payloads
        .iter()
        .enumerate()
        .map(|(i, (payload, _))| {
            let mut compressor = fgumi_bgzf::writer::InlineBgzfCompressor::new(1);
            compressor.write_all(payload).unwrap();
            compressor.flush().unwrap();
            let mut out = compressor.take_blocks();
            assert_eq!(out.len(), 1, "a group of records fits one BGZF block");
            crate::types::BgzfBlock {
                batch_serial: i as u64,
                bytes: out.remove(0).data.into(),
                uncompressed_size: u32::try_from(payload.len()).unwrap(),
                index: None,
            }
        })
        .collect();
    (blocks, shape)
}

/// The `Spill.records_ingested_so_far` sequence the charged caps imply, worked
/// by hand: run 0's cap is `M·258/(258+w0)`; run k > 0's is `M·m/(m+w)` with
/// `m = ⌊bytes_{k−1} / records_{k−1}⌋`. Mirroring `ReadBlocks`, a run's seal is
/// armed by the first block whose append makes the run's bytes reach the cap
/// (`run_cumulative >= run_cap`, tested after appending), and fires — as a
/// spill — only when a following block arrives.
fn expected_boundaries(blocks: &[(u64, u64)], memory_limit: u64, w0: u64, w: u64) -> Vec<u64> {
    let mut cap = memory_limit * 258 / (258 + w0);
    let (mut bytes, mut records, mut ingested) = (0u64, 0u64, 0u64);
    let mut boundaries = Vec::new();
    for (i, &(isize, n)) in blocks.iter().enumerate() {
        bytes += isize;
        records += n;
        ingested += n;
        if bytes >= cap && i + 1 < blocks.len() {
            boundaries.push(ingested);
            let m = bytes / records;
            cap = memory_limit * m / (m + w);
            (bytes, records) = (0, 0);
        }
    }
    boundaries
}

/// Spill boundaries follow the charged caps exactly, at one and four threads:
/// template order under `--key-types auto` (run 0 charged at the widest 56-byte
/// ref, later runs at the provisioned 40), queryname (40) and coordinate (24),
/// each over at least four runs. Two runs at four threads split identically.
#[rstest]
#[case::template_auto(SortOrder::TemplateCoordinate, 56, 40)]
#[case::queryname(SortOrder::Queryname(QuerynameComparator::Lexicographic), 40, 40)]
#[case::coordinate(SortOrder::Coordinate, 24, 24)]
fn spill_boundaries_match_the_hand_computed_caps(
    #[case] order: SortOrder,
    #[case] w0: u64,
    #[case] w: u64,
) {
    let (header, records) = match order {
        SortOrder::Queryname(_) => instrument_switching_records(3_000, 0xB0B),
        _ => template_records(3_000, 0xB0B, false, false),
    };
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let (_, shape) = whole_record_blocks(&records, n_ref, 16);
    let memory_limit = 64 * 1024;
    let expected = expected_boundaries(&shape, memory_limit as u64, w0, w);
    assert!(expected.len() >= 4, "{order:?}: at least four runs spill ({expected:?})");
    let config = buffer_config(order, memory_limit, 2);
    let run = |threads| {
        let (blocks, _) = whole_record_blocks(&records, n_ref, 16);
        drive_arena_chunks(
            ArenaFront::Split,
            config,
            &header,
            blocks,
            threads,
            &ChainCaps::default(),
        )
        .unwrap()
        .spill_ingested
    };
    assert_eq!(run(1), expected, "{order:?} t1");
    let t4 = run(4);
    assert_eq!(t4, expected, "{order:?} t4");
    assert_eq!(run(4), t4, "{order:?}: two t4 runs split identically");
}

/// Residual chunks a front announced: chunks beyond the spills.
fn residual_count(out: &ArenaFrontOutput) -> usize {
    out.chunks.len() - out.slot_count as usize
}

/// When the sort spilled, the final in-memory run is sealed as `sort_threads`
/// pre-sorted pieces (contiguous input ranges), each its own `Residual`; the
/// merge of the spills and the pieces is the reference stable sort, equal keys
/// in input order. Every order, at two and four sort threads.
#[rstest]
fn residual_pieces_merge_to_the_reference(
    #[values(
        SortOrder::Coordinate,
        SortOrder::TemplateCoordinate,
        SortOrder::Queryname(QuerynameComparator::Lexicographic),
        SortOrder::Queryname(QuerynameComparator::Natural)
    )]
    order: SortOrder,
    #[values(2usize, 4)] sort_threads: usize,
) {
    // Queryname ties (equal name and flags) are dense in the instrument
    // fixture's duplicated secondaries, so pieces of the final run share them.
    let (header, records) = match order {
        SortOrder::Queryname(_) => instrument_switching_records(8_000, 62),
        _ => tie_heavy_records(6_000, 61),
    };
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let config = buffer_config(order, 96 * 1024, sort_threads);
    let blocks = || bgzf_blocks_for(&records, n_ref, 4096);
    let out = drive_arena_chunks(
        ArenaFront::Split,
        config,
        &header,
        blocks(),
        sort_threads,
        &ChainCaps::default(),
    )
    .unwrap();
    assert!(out.slot_count >= 1, "{order:?}: the fixture spills");
    assert_eq!(residual_count(&out), sort_threads, "{order:?}: one residual per piece");
    let run = drive_arena_chain(
        ArenaFront::Split,
        config,
        &header,
        blocks(),
        4 * 1024 * 1024,
        sort_threads,
        MergeSink::Batches(StepKind::Serial),
        &ChainCaps::default(),
    )
    .unwrap();
    let expected = reference_sort(order, &header, &records, None);
    // The final run is the input tail after the last spill, split into
    // `sort_threads` contiguous pieces; some equal-key group of the tail must
    // straddle a piece boundary, or a cross-piece reorder would go unseen.
    let tail_start = *run.spill_boundaries.last().expect("the fixture spills");
    let piece = (records.len() as u64 - tail_start).div_ceil(sort_threads as u64);
    let piece_starts: Vec<u64> = (1..sort_threads as u64).map(|k| tail_start + k * piece).collect();
    let tail: Vec<Vec<u8>> =
        expected.iter().filter(|r| u64::from(xi_of(r)) >= tail_start).cloned().collect();
    assert!(
        equal_key_group_spans_a_boundary(order, &header, &tail, &piece_starts),
        "{order:?} t{sort_threads}: no equal-key group straddles a piece boundary"
    );
    assert_record_parity(&format!("{order:?} t{sort_threads}"), &run.records, &expected);
    assert_input_order_within_equal_keys(order, &header, &run.records);
    let buffer = drive_front(Front::SortBuffer, config, &header, &records, sort_threads).unwrap();
    assert_record_parity(
        &format!("{order:?} t{sort_threads} SortBuffer"),
        &buffer.records,
        &expected,
    );
}

/// A sort that never spills keeps its single in-memory chunk, at any thread
/// count: the merge's single-chunk fast path stays available.
#[rstest]
fn in_memory_sort_keeps_one_residual_chunk(
    #[values(SortOrder::Coordinate, SortOrder::TemplateCoordinate)] order: SortOrder,
    #[values(1usize, 4)] sort_threads: usize,
) {
    let (header, records) = tie_heavy_records(2_000, 62);
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let config = buffer_config(order, 64 << 20, sort_threads);
    let out = drive_arena_chunks(
        ArenaFront::Split,
        config,
        &header,
        bgzf_blocks_for(&records, n_ref, 4096),
        sort_threads,
        &ChainCaps::default(),
    )
    .unwrap();
    assert_eq!((out.slot_count, residual_count(&out)), (0, 1), "{order:?} t{sort_threads}");
}

/// At one sort thread the final run of a spilling sort is never split.
#[rstest]
fn t1_never_splits(
    #[values(SortOrder::Coordinate, SortOrder::Queryname(QuerynameComparator::Natural))]
    order: SortOrder,
) {
    let (header, records) = tie_heavy_records(6_000, 63);
    let n_ref = u32::try_from(header.reference_sequences().len()).unwrap();
    let out = drive_arena_chunks(
        ArenaFront::Split,
        buffer_config(order, 96 * 1024, 1),
        &header,
        bgzf_blocks_for(&records, n_ref, 4096),
        1,
        &ChainCaps::default(),
    )
    .unwrap();
    assert!(out.slot_count >= 1);
    assert_eq!(residual_count(&out), 1, "{order:?}");
}
