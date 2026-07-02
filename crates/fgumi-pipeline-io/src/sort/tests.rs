//! Tests for the runall-sort three-step chain
//! (`SortAndSpill` → `SortSpillDecompress` → `SortMerge`).
//!
//! `RawExternalSorter::sort` (driven here via [`sort_via_legacy`]) is retained
//! as the parity oracle the streaming three-step chain is validated against.

use std::io;
use std::sync::Arc;

use anyhow::Result;
use fgumi_raw_bam::RawRecord;
use fgumi_raw_bam::testutil::make_bam_bytes;
use fgumi_sort::{QuerynameComparator, RawExternalSorter, SortOrder, SpillCodec};
use noodles::sam::Header;
use parking_lot::Mutex;
use rstest::rstest;

use super::*;
use crate::types::RecordBatch;
use fgumi_pipeline_core::{
    Unpushed,
    builder::{Pipeline, PipelineConfig},
    held::HeldSlot,
    outputs::OrderedBytesSingle,
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
                self.received.lock().push(batch);
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// One sort's output as a flat list of raw BAM record-byte payloads, in output
/// order. Both the streaming pipeline and the legacy oracle produce this shape.
type RecordBytes = Vec<Vec<u8>>;

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

fn drive_sort_pipeline(
    sorter: RawExternalSorter,
    header: &Header,
    batches: Vec<RecordBatch>,
    output_byte_limit: u64,
    threads: usize,
    sink_kind: StepKind,
) -> Result<Vec<Vec<u8>>> {
    drive_sort_pipeline_tuned(
        sorter,
        header,
        batches,
        output_byte_limit,
        threads,
        sink_kind,
        SortDecompressTuning::default(),
        SpillCodec::Zstd,
    )
}

/// Drive the production sort chain (`VecSource` → `SortBuffer` →
/// `CompressSpill` → `SortSpillDecompress` → `SortMerge` → `VecSink`). The
/// legacy `SortAndSpill` Phase-1 head was retired in P7, so the Phase-2 tests
/// (decompress-granularity, out-of-order, soak) run through the same production
/// chain `fgumi sort` / `runall` build.
#[allow(clippy::too_many_arguments)]
fn drive_sort_pipeline_tuned(
    sorter: RawExternalSorter,
    header: &Header,
    batches: Vec<RecordBatch>,
    output_byte_limit: u64,
    threads: usize,
    sink_kind: StepKind,
    decompress_tuning: SortDecompressTuning,
    spill_codec: SpillCodec,
) -> Result<Vec<Vec<u8>>> {
    drive_sort_buffer_pipeline(
        sorter,
        header,
        batches,
        output_byte_limit,
        threads,
        sink_kind,
        decompress_tuning,
        spill_codec,
    )
}

/// Drive the P6 four-step buffer chain
/// (`VecSource` → `SortBuffer` → `CompressSpill` → `SortSpillDecompress` →
/// `SortMerge` → `VecSink`) and collect the merged record bytes. Exercises
/// every sort order `SortBuffer` supports (see
/// `sort_buffer_chain_matches_legacy_all_orders`).
#[allow(clippy::too_many_arguments)]
fn drive_sort_buffer_pipeline(
    sorter: RawExternalSorter,
    header: &Header,
    batches: Vec<RecordBatch>,
    output_byte_limit: u64,
    threads: usize,
    sink_kind: StepKind,
    decompress_tuning: SortDecompressTuning,
    spill_codec: SpillCodec,
) -> Result<Vec<Vec<u8>>> {
    use fgumi_sort::TmpDirAllocator;

    let received: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sort_order = sorter.sort_order();

    // Temp dir + allocator for CompressSpill, held alive by the step. The
    // deterministic always-ample probe avoids any dependency on host free space.
    let dir = tempfile::TempDir::new()?;
    let alloc =
        TmpDirAllocator::with_probe(vec![dir.path().to_path_buf()], Box::new(|_| Ok(u64::MAX)), 0)?;
    let temp_dirs = Arc::new(vec![dir]);

    let source = VecSource::new(batches, output_byte_limit);
    let sort_buffer = SortBuffer::from_sorter(sorter, header, output_byte_limit)?;
    // Codec/compression affect only intermediate spill bytes, not the final
    // sorted records, so any codec yields output parity. The caller passes the
    // codec so codec-specific tests (e.g. the BGZF block-parallel parity test)
    // actually exercise their codec end-to-end.
    let compress = CompressSpill::new(
        Arc::new(Mutex::new(alloc)),
        spill_codec,
        3,
        output_byte_limit,
        temp_dirs,
    );
    let decompress = SortSpillDecompress::new(output_byte_limit, decompress_tuning);
    let merge =
        SortMerge::<RecordBatchOutput>::with_target_batch_count(sort_order, output_byte_limit, 256);
    let sink = VecSink { received: Arc::clone(&received), kind: sink_kind };

    let builder = Pipeline::builder();
    builder
        .chain(source)
        .chain(sort_buffer)
        .chain(compress)
        .chain(decompress)
        .chain(merge)
        .chain(sink)
        .into_sink_marker();
    let pipeline = builder.build()?;
    pipeline.run(PipelineConfig { threads, ..Default::default() })?;

    let collected = std::mem::take(&mut *received.lock());
    let mut out = Vec::new();
    for batch in collected {
        for bytes in batch.iter_record_bytes() {
            out.push(bytes.to_vec());
        }
    }
    Ok(out)
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
/// path / legacy oracle (lever 1 parity gate).
#[allow(clippy::too_many_arguments)]
fn drive_sort_block_pipeline(
    sorter: RawExternalSorter,
    header: &Header,
    batches: Vec<RecordBatch>,
    output_byte_limit: u64,
    threads: usize,
    sink_kind: StepKind,
    decompress_tuning: SortDecompressTuning,
    spill_codec: SpillCodec,
) -> Result<Vec<Vec<u8>>> {
    use fgumi_sort::TmpDirAllocator;

    let received: Arc<Mutex<Vec<crate::types::DecompressedBlock>>> =
        Arc::new(Mutex::new(Vec::new()));
    let sort_order = sorter.sort_order();

    let dir = tempfile::TempDir::new()?;
    let alloc =
        TmpDirAllocator::with_probe(vec![dir.path().to_path_buf()], Box::new(|_| Ok(u64::MAX)), 0)?;
    let temp_dirs = Arc::new(vec![dir]);

    let source = VecSource::new(batches, output_byte_limit);
    let sort_buffer = SortBuffer::from_sorter(sorter, header, output_byte_limit)?;
    let compress = CompressSpill::new(
        Arc::new(Mutex::new(alloc)),
        spill_codec,
        3,
        output_byte_limit,
        temp_dirs,
    );
    let decompress = SortSpillDecompress::new(output_byte_limit, decompress_tuning);
    let merge =
        SortMerge::<BlockOutput>::with_target_batch_count(sort_order, output_byte_limit, 256);
    let sink = BlockSink { received: Arc::clone(&received), kind: sink_kind };

    let builder = Pipeline::builder();
    builder
        .chain(source)
        .chain(sort_buffer)
        .chain(compress)
        .chain(decompress)
        .chain(merge)
        .chain(sink)
        .into_sink_marker();
    let pipeline = builder.build()?;
    pipeline.run(PipelineConfig { threads, ..Default::default() })?;

    let mut collected = std::mem::take(&mut *received.lock());
    // `Detached` collapses ordering to None, so blocks arrive in merge-emit
    // order on the single sink; sort by serial defensively in case a future
    // change makes the sink parallel.
    collected.sort_by_key(|b| b.batch_serial);
    let mut out = Vec::new();
    for block in &collected {
        out.extend(unframe_block_records(&block.bytes));
    }
    Ok(out)
}

// ── Reference: RawExternalSorter::sort to bytes ─────────────────────────────

fn sort_via_legacy(
    sort_order: SortOrder,
    header: &Header,
    records: &[RawRecord],
    memory_limit: usize,
    threads: usize,
) -> Result<Vec<Vec<u8>>> {
    let tmp_in = tempfile::NamedTempFile::new()?;
    {
        let mut writer = fgumi_bam_io::create_raw_bam_writer(tmp_in.path(), header, 1, 1)?;
        for r in records {
            writer.write_raw_record(r.as_ref())?;
        }
        writer.finish()?;
    }

    let tmp_out = tempfile::NamedTempFile::new()?;
    let sorter = RawExternalSorter::new(sort_order)
        .memory_limit(memory_limit)
        .threads(threads)
        .output_compression(1)
        .temp_compression(1);
    sorter.sort(tmp_in.path(), tmp_out.path())?;

    let (mut reader, _hdr) = fgumi_bam_io::create_raw_bam_reader_with_opts(
        tmp_out.path(),
        1,
        fgumi_bam_io::PipelineReaderOpts::default(),
    )?;
    let mut record = RawRecord::default();
    let mut out = Vec::new();
    loop {
        match reader.read_record(&mut record)? {
            0 => break,
            _ => out.push(record.as_ref().to_vec()),
        }
    }
    Ok(out)
}

// ── Tests ───────────────────────────────────────────────────────────────────
//
// (The former `three_step_chain_{in_memory,multi_spill}_path_matches_legacy`
// tests, which drove the retired `SortAndSpill` head, are subsumed by
// `sort_buffer_chain_matches_legacy_all_orders` below — same orders, same
// regimes, against the same `sort_via_legacy` oracle, but through the
// production `SortBuffer` → `CompressSpill` chain.)

// ── P6 buffer-chain parity (SortBuffer → CompressSpill → decompress → merge) ─

/// The P6 four-step chain must produce byte-identical output to the legacy
/// oracle for EVERY sort order, across the in-memory and multi-spill regimes at
/// 2 and 4 threads. This is the inc-4/5 parity gate: it exercises `SortBuffer`'s
/// streaming-spill emission, `CompressSpill`'s inline writes + slot opens
/// (`file_id == seq`), and the single-residual fast path all the way through the
/// real `SortSpillDecompress` + `SortMerge` — for coordinate, template, and both
/// queryname comparators.
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
fn sort_buffer_chain_matches_legacy_all_orders(
    #[case] sort_order: SortOrder,
    #[case] n: usize,
    #[case] memory_limit: usize,
    #[case] threads: usize,
) {
    let (header, records) = synthesize_records(n, 0x5A17_C0DE);
    let sorter = RawExternalSorter::new(sort_order)
        .memory_limit(memory_limit)
        .threads(threads)
        .output_compression(1)
        .temp_compression(1);
    let new_out = drive_sort_buffer_pipeline(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        // Run the pipeline with exactly the case's thread count — production
        // sort uses `num_threads` workers (no `max(3)` floor), so the parity
        // test must too.
        threads,
        StepKind::Exclusive,
        SortDecompressTuning::default(),
        SpillCodec::Zstd,
    )
    .expect("buffer pipeline drives to completion");

    let legacy_out =
        sort_via_legacy(sort_order, &header, &records, memory_limit, threads).expect("legacy");

    assert_eq!(new_out.len(), legacy_out.len(), "{sort_order:?} record count mismatch");
    if is_stable_for_equal_keys(sort_order) {
        // Coordinate / template-coordinate use a stable radix sort: equal keys
        // keep input order deterministically, so output is byte-for-byte equal.
        for (i, (got, want)) in new_out.iter().zip(legacy_out.iter()).enumerate() {
            assert_eq!(got, want, "{sort_order:?} record {i} bytes differ");
        }
    } else {
        // Queryname uses an UNSTABLE comparator sort: the order of equal keys
        // (same name + segment flags) is unspecified and differs between the
        // parallel buffer chain and the `.sort()` oracle (and run-to-run). The
        // sound parity claim is multiset equality — same records, possibly in a
        // different equal-key order. Sortedness is covered by the production
        // `queryname_*_sort_matrix` integration tests.
        let mut got = new_out.clone();
        let mut want = legacy_out.clone();
        got.sort_unstable();
        want.sort_unstable();
        assert_eq!(got, want, "{sort_order:?} record multiset differs from the oracle");
    }
}

/// Regression guard for the `SortBuffer` peak-memory invariant: when a single
/// input `RecordBatch` is larger than `memory_limit`, `ingest_one_batch` seals
/// several chunks in one call (staging multiple `Spill` events into `pending`
/// before `emit_pending` runs). That transiently exceeds the "~one spill chunk"
/// production bound, but it MUST stay correct — every record emitted in sorted
/// order, none stranded by the loop. Packs all records into ONE oversized batch
/// against a small `memory_limit` (the case the 256-record-per-batch parity
/// matrix never hits) and asserts byte-for-byte parity with the legacy oracle
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

    let sorter = RawExternalSorter::new(sort_order)
        .memory_limit(memory_limit)
        .threads(threads)
        .output_compression(1)
        .temp_compression(1);
    let new_out = drive_sort_buffer_pipeline(
        sorter,
        &header,
        batches,
        4 * 1024 * 1024,
        threads,
        StepKind::Exclusive,
        SortDecompressTuning::default(),
        SpillCodec::Zstd,
    )
    .expect("buffer pipeline drives to completion");

    let legacy_out =
        sort_via_legacy(sort_order, &header, &records, memory_limit, threads).expect("legacy");

    assert_eq!(new_out.len(), records.len(), "every record must survive — none dropped");
    assert_eq!(new_out.len(), legacy_out.len(), "record count matches oracle");
    for (i, (got, want)) in new_out.iter().zip(legacy_out.iter()).enumerate() {
        assert_eq!(got, want, "record {i} bytes differ from oracle");
    }
}

/// Lever 1 parity gate: the terminal `SortMerge<BlockOutput>` path (framed
/// `DecompressedBlock`s, wired straight to `BgzfCompress`) must carry
/// byte-identical records to the legacy oracle for every sort order, across the
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
fn sort_merge_block_output_matches_legacy(
    #[case] sort_order: SortOrder,
    #[case] n: usize,
    #[case] memory_limit: usize,
    #[case] threads: usize,
) {
    let (header, records) = synthesize_records(n, 0x5A17_C0DE);
    let make_sorter = || {
        RawExternalSorter::new(sort_order)
            .memory_limit(memory_limit)
            .threads(threads)
            .output_compression(1)
            .temp_compression(1)
    };

    // Terminal framed-block path.
    let block_out = drive_sort_block_pipeline(
        make_sorter(),
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        threads,
        StepKind::Exclusive,
        SortDecompressTuning::default(),
        SpillCodec::Zstd,
    )
    .expect("block pipeline drives to completion");

    // Intermediate RecordBatch path — the framed bytes must un-frame to exactly
    // the records this path emits (cross-check the two SortMerge framings agree).
    let batch_out = drive_sort_buffer_pipeline(
        make_sorter(),
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        threads,
        StepKind::Exclusive,
        SortDecompressTuning::default(),
        SpillCodec::Zstd,
    )
    .expect("buffer pipeline drives to completion");

    let legacy_out =
        sort_via_legacy(sort_order, &header, &records, memory_limit, threads).expect("legacy");

    assert_eq!(block_out.len(), legacy_out.len(), "{sort_order:?} record count vs legacy");
    assert_eq!(block_out.len(), batch_out.len(), "{sort_order:?} record count vs RecordBatch path");

    if is_stable_for_equal_keys(sort_order) {
        // Stable orders: byte-for-byte identical output, in order.
        for (i, (got, want)) in block_out.iter().zip(legacy_out.iter()).enumerate() {
            assert_eq!(got, want, "{sort_order:?} record {i} (block path vs legacy) differs");
        }
        assert_eq!(block_out, batch_out, "{sort_order:?} block vs RecordBatch path bytes differ");
    } else {
        // Queryname's comparator sort is unstable on equal keys; assert multiset
        // equality (same records, possibly different equal-key order).
        let mut got = block_out.clone();
        let mut want = legacy_out.clone();
        got.sort_unstable();
        want.sort_unstable();
        assert_eq!(got, want, "{sort_order:?} block-path record multiset differs from oracle");
    }
}

/// `true` for sort orders whose sort is stable on equal keys (so byte-for-byte
/// output parity is deterministic). Queryname's comparator sort is unstable.
fn is_stable_for_equal_keys(order: SortOrder) -> bool {
    matches!(order, SortOrder::Coordinate | SortOrder::TemplateCoordinate)
}

/// L2.6: a `StepKind::Detached` SINK — the `WriteBgzfFile` analogue, i.e. the
/// detached-thread runtime driving a pure consumer — yields byte-identical
/// merged output to the legacy oracle, exactly like the pool-scheduled sink.
/// The chain's `SortMerge` is already `Detached`, so this drives the full chain
/// through `pipeline.run` with TWO off-pool threads (merge + sink) over a
/// multi-spill coordinate workload, pinning that `run_detached_driver` preserves
/// the output bytes for both the producing (merge) and consuming (sink) roles.
#[test]
fn detached_sink_chain_matches_legacy_coordinate() {
    let memory_limit = 64 * 1024; // small → forces many real spill files
    let threads = 4;
    let (header, records) = synthesize_records(20_000, 0xD17A_C4ED);
    let sorter = RawExternalSorter::new(SortOrder::Coordinate)
        .memory_limit(memory_limit)
        .threads(threads)
        .output_compression(1)
        .temp_compression(1);
    let detached_out = drive_sort_buffer_pipeline(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        threads,
        StepKind::Detached,
        SortDecompressTuning::default(),
        SpillCodec::Zstd,
    )
    .expect("detached-sink buffer pipeline drives to completion");

    let legacy_out =
        sort_via_legacy(SortOrder::Coordinate, &header, &records, memory_limit, threads)
            .expect("legacy");

    assert_eq!(detached_out.len(), legacy_out.len(), "record count mismatch");
    for (i, (got, want)) in detached_out.iter().zip(legacy_out.iter()).enumerate() {
        assert_eq!(got, want, "detached-sink record {i} bytes differ from oracle");
    }
}

/// The buffer chain must hold coordinate parity across BOTH decompress
/// granularities × block batches, with a multi-spill workload that forces real
/// spill files (so `CompressSpill`'s written chunks feed the block-parallel
/// reorder path). Guards against any spill-format / slot-ordering drift between
/// `CompressSpill` and the proven `SortSpillDecompress` reader.
#[rstest]
#[case::file_b1(true, 1)]
#[case::file_b4(true, 4)]
#[case::block_b1(false, 1)]
#[case::block_b4(false, 4)]
fn sort_buffer_chain_coordinate_matches_legacy_across_decompress_tunings(
    #[case] file_granularity: bool,
    #[case] block_batch: usize,
) {
    let (header, records) = synthesize_records(20_000, 0xBADD_CAFE);
    let memory_limit = 256 * 1024;
    let threads = 2;
    let sorter = RawExternalSorter::new(SortOrder::Coordinate)
        .memory_limit(memory_limit)
        .threads(threads)
        .output_compression(1)
        .temp_compression(1);
    let new_out = drive_sort_buffer_pipeline(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        // Exactly the case's thread count (see the sibling test): production
        // runs the pipeline with `num_threads`, no `max(3)` floor.
        threads,
        StepKind::Exclusive,
        SortDecompressTuning { file_granularity, block_batch },
        SpillCodec::Zstd,
    )
    .expect("buffer pipeline drives to completion");

    let legacy_out =
        sort_via_legacy(SortOrder::Coordinate, &header, &records, memory_limit, threads)
            .expect("legacy");

    assert_eq!(new_out.len(), legacy_out.len(), "record count mismatch");
    for (i, (got, want)) in new_out.iter().zip(legacy_out.iter()).enumerate() {
        assert_eq!(got, want, "record {i} bytes differ");
    }
}

/// Equal-key stability across spill boundaries — the output-identity-critical
/// tie-break the single-residual design depends on. Every record shares one
/// coordinate (tid 0, pos 0), so a stable sort must emit them in input order;
/// the buffer chain must preserve that across multiple spill chunks (tie-broken
/// by `file_id == seq`) and the residual. Pinned both directly (output == input
/// order) and against the legacy oracle.
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

    let sorter = RawExternalSorter::new(SortOrder::Coordinate)
        .memory_limit(memory_limit)
        .threads(2)
        .output_compression(1)
        .temp_compression(1);
    let out = drive_sort_buffer_pipeline(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        2,
        StepKind::Exclusive,
        SortDecompressTuning::default(),
        SpillCodec::Zstd,
    )
    .expect("buffer pipeline drives to completion");

    assert_eq!(out.len(), input_bytes.len(), "record count mismatch");
    assert_eq!(out, input_bytes, "equal-key records must preserve input order across spills");

    let legacy =
        sort_via_legacy(SortOrder::Coordinate, &header, &records, memory_limit, 2).expect("legacy");
    assert_eq!(out, legacy, "equal-key order must also match the legacy oracle");
}

// ── Decompression-granularity parity (file-granularity × block-batch) ────────

/// The streaming sort must produce byte-identical output for BOTH decompression
/// granularities (`file_granularity ∈ {true, false}`) across `block_batch ∈
/// {1, 4}`, validated against the legacy reference. The multi-spill workload
/// forces real spill files so the block-parallel reorder path is exercised (the
/// in-memory-only path never opens a slot reader).
#[rstest]
#[case::file_b1(true, 1)]
#[case::file_b4(true, 4)]
#[case::block_b1(false, 1)]
#[case::block_b4(false, 4)]
// block_batch == 0 is clamped to 1 in `SortSpillDecompress::new`. Without the
// clamp, the inline path declares a phantom EOF after reading zero blocks
// (silent record loss) and the block-parallel path livelocks (queue_eof never
// finalizes). These cases assert the clamp holds: identical to legacy, no hang.
#[case::file_b0(true, 0)]
#[case::block_b0(false, 0)]
fn three_step_chain_granularity_matrix_matches_legacy(
    #[case] file_granularity: bool,
    #[case] block_batch: usize,
) {
    let sort_order = SortOrder::Coordinate;
    let threads = 4;
    let (header, records) = synthesize_sized_records(30_000, 0x5EED_1234, 120);
    // Small per-thread memory ⇒ many spill files ⇒ many slot blocks.
    let memory_limit = 256 * 1024;

    let sorter = RawExternalSorter::new(sort_order)
        .memory_limit(memory_limit)
        .threads(2)
        .output_compression(1)
        .temp_compression(1);
    let new_out = drive_sort_pipeline_tuned(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        threads,
        StepKind::Exclusive,
        SortDecompressTuning { file_granularity, block_batch },
        SpillCodec::Zstd,
    )
    .expect("pipeline drives to completion");

    let legacy_out =
        sort_via_legacy(sort_order, &header, &records, memory_limit, 2).expect("legacy");

    assert_eq!(
        new_out.len(),
        legacy_out.len(),
        "record count mismatch (file_granularity={file_granularity}, block_batch={block_batch})"
    );
    assert_eq!(
        new_out, legacy_out,
        "sorted bytes differ (file_granularity={file_granularity}, block_batch={block_batch})"
    );
}

/// Block-parallel decompression over BGZF spill files (the non-default codec)
/// must also match the legacy sorter — the matrix/soak/proptest exercise the
/// default zstd spills, so this confirms the block-parallel reorder path is
/// codec-agnostic. (Output records are codec-independent: the spill codec only
/// affects temp files, not the sorted output.)
#[test]
fn block_parallel_bgzf_spill_matches_legacy() {
    let sort_order = SortOrder::Coordinate;
    let (header, records) = synthesize_sized_records(30_000, 0x5EED_BEEF, 120);
    let memory_limit = 256 * 1024;
    let sorter = RawExternalSorter::new(sort_order)
        .memory_limit(memory_limit)
        .threads(2)
        .output_compression(1)
        .temp_compression(1)
        .spill_codec(fgumi_sort::SpillCodec::Bgzf);
    let new_out = drive_sort_pipeline_tuned(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        8,
        StepKind::Exclusive,
        SortDecompressTuning { file_granularity: false, block_batch: 2 },
        SpillCodec::Bgzf,
    )
    .expect("pipeline drives to completion");
    let legacy_out =
        sort_via_legacy(sort_order, &header, &records, memory_limit, 2).expect("legacy");
    assert_eq!(new_out.len(), legacy_out.len(), "record count mismatch (bgzf block-parallel)");
    assert_eq!(new_out, legacy_out, "sorted bytes differ (bgzf block-parallel)");
}

/// Block-parallel decompression completes out of order (workers decompress one
/// file's blocks concurrently), yet the reassembled output must be byte-
/// identical to the in-order (file-granularity) result. Property test over a
/// range of record counts and `block_batch` sizes and a high pipeline-thread
/// count (more concurrent decompressors ⇒ more out-of-order completion). A
/// straggler worker hitting reader-EOF while another holds an in-flight block
/// must not truncate the output (record count is asserted equal).
#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #![proptest_config(ProptestConfig { cases: 24, ..ProptestConfig::default() })]

        #[test]
        fn block_parallel_matches_file_granularity(
            n_records in 2_000usize..18_000,
            block_batch in 1usize..=6,
            seed in any::<u64>(),
        ) {
            let sort_order = SortOrder::Coordinate;
            let (header, records) = synthesize_sized_records(n_records, seed, 100);
            // Force spilling so slots (and the reorder path) are exercised.
            let memory_limit = 256 * 1024;
            let pipeline_threads = 6;

            let make_sorter = || RawExternalSorter::new(sort_order)
                .memory_limit(memory_limit)
                .threads(2)
                .output_compression(1)
                .temp_compression(1);

            let in_order = drive_sort_pipeline_tuned(
                make_sorter(),
                &header,
                pack_batches(&records, 256),
                4 * 1024 * 1024,
                pipeline_threads,
                StepKind::Exclusive,
                SortDecompressTuning { file_granularity: true, block_batch },
                SpillCodec::Zstd,
            ).expect("file-granularity pipeline");

            let out_of_order = drive_sort_pipeline_tuned(
                make_sorter(),
                &header,
                pack_batches(&records, 256),
                4 * 1024 * 1024,
                pipeline_threads,
                StepKind::Exclusive,
                SortDecompressTuning { file_granularity: false, block_batch },
                SpillCodec::Zstd,
            ).expect("block-parallel pipeline");

            prop_assert_eq!(out_of_order.len(), records.len(), "no truncation");
            prop_assert_eq!(out_of_order, in_order, "block-parallel diverges from in-order");
        }
    }
}

/// Maximum-contention soak for the block-parallel decompress path
/// (`file_granularity == false`). Drives the path repeatedly under the most
/// adversarial settings the knobs allow — many spill files, a tiny reorder
/// window (so stragglers continuously hit `bp_reorder_admits` backpressure and
/// the Phase-B drain-only path), `block_batch == 1` (maximum per-block churn and
/// the most frequent `reader_eof`/`in_flight` transitions), and far more
/// pipeline worker threads (12) than sorter threads (so many workers race to
/// decompress one file's blocks concurrently and finalize out of order).
///
/// Each iteration uses a fresh random seed and is checked for *byte-identity*
/// against the legacy `RawExternalSorter::sort` oracle — so a lost, duplicated,
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
fn block_parallel_high_contention_soak_matches_legacy() {
    use std::time::Duration;

    const ITERATIONS: usize = 40;
    const RECORDS_PER_ITER: usize = 15_000;
    const PIPELINE_THREADS: usize = 12;
    const SORTER_THREADS: usize = 2;
    // Tiny per-thread sort memory ⇒ many spill files ⇒ many slot readers.
    const MEMORY_LIMIT: usize = 128 * 1024;
    // Tiny output/reorder-window budget ⇒ the block-parallel reorder window is
    // ~1 block, so `bp_reorder_admits` backpressures aggressively and workers
    // are repeatedly forced through the Phase-B drain-only path.
    const OUTPUT_BYTE_LIMIT: u64 = 64 * 1024;
    const BLOCK_BATCH: usize = 1;
    // Per-iteration watchdog: a livelock in any single iteration fails fast.
    const WATCHDOG: Duration = Duration::from_secs(60);

    let sort_order = SortOrder::Coordinate;
    for iter in 0..ITERATIONS {
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
            SortDecompressTuning { file_granularity: false, block_batch: BLOCK_BATCH },
            WATCHDOG,
        );
    }
}

/// R1 (Campaign 3, N+2): the capacity-1 arena cycle now spans two thread
/// domains — the COORDINATION driver runs `ReadBlocks` (acquire+admit) and
/// `FindBoundariesAndSort` (seal+free) while the POOL runs `InflateToArena` in
/// between. The single arena segment can't be re-acquired for run *k+1* until
/// run *k* seals AND its chunk is framed AND every `Arc` clone drops. This is a
/// three-party cross-domain cycle the pre-N+2 `detached_two_sided_no_deadlock`
/// (single detached step) does not cover.
///
/// This soak drives that cycle at maximum churn: a *tiny* memory limit seals an
/// arena after only a handful of blocks, so run acquire/inflate/seal/free cycles
/// as fast as possible while the driver round-robins admit vs seal and the pool
/// inflates. A wedge (e.g. the driver parking on `FindBoundariesAndSort` while
/// `SpillGather` holds the chunk that frees the arena — the G1 failure mode) trips
/// the per-iteration watchdog and fails fast instead of hanging. Byte-parity vs
/// the legacy oracle also proves no arena slot is aliased or lost under the churn.
#[test]
fn arena_capacity1_driver_pool_soak_no_deadlock() {
    use std::time::Duration;

    const ITERATIONS: usize = 24;
    const RECORDS_PER_ITER: usize = 12_000;
    const PIPELINE_THREADS: usize = 8;
    const SORTER_THREADS: usize = 4;
    // Very tight sort memory ⇒ an arena seals after only a few blocks ⇒ maximal
    // capacity-1 acquire/reuse churn between the coordination driver and pool.
    const MEMORY_LIMIT: usize = 24 * 1024;
    const OUTPUT_BYTE_LIMIT: u64 = 48 * 1024;
    const WATCHDOG: Duration = Duration::from_secs(60);

    for iter in 0..ITERATIONS {
        let seed = 0xC0FF_EE00_u64 ^ (iter as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        let (header, records) = synthesize_sized_records(RECORDS_PER_ITER, seed, 130);
        run_watchdogged_parity(
            &format!("arena-cap1-soak-i{iter}"),
            SortOrder::Coordinate,
            header,
            records,
            MEMORY_LIMIT,
            SORTER_THREADS,
            PIPELINE_THREADS,
            OUTPUT_BYTE_LIMIT,
            SortDecompressTuning { file_granularity: false, block_batch: 1 },
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
            "{label}: output diverges from legacy at record {idx}: \
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
/// legacy oracle. Owns `header`/`records` so the worker thread can take them.
///
/// Both the legacy oracle AND the streaming pipeline run *inside* the worker, so
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
    tuning: SortDecompressTuning,
    watchdog: std::time::Duration,
) {
    use std::sync::mpsc;

    let (tx, rx) = mpsc::channel();
    let worker = std::thread::Builder::new()
        .name(label.to_string())
        .spawn(move || {
            // (pipeline_out, legacy_out) — both computed under the watchdog.
            let result = (|| -> Result<(RecordBytes, RecordBytes)> {
                let legacy_out =
                    sort_via_legacy(sort_order, &header, &records, memory_limit, sorter_threads)?;
                let batches = pack_batches(&records, 256);
                let sorter = RawExternalSorter::new(sort_order)
                    .memory_limit(memory_limit)
                    .threads(sorter_threads)
                    .output_compression(1)
                    .temp_compression(1);
                let out = drive_sort_pipeline_tuned(
                    sorter,
                    &header,
                    batches,
                    output_byte_limit,
                    pipeline_threads,
                    StepKind::Serial,
                    tuning,
                    SpillCodec::Zstd,
                )?;
                Ok((out, legacy_out))
            })();
            let _ = tx.send(result);
        })
        .expect("spawn soak worker");

    match rx.recv_timeout(watchdog) {
        Ok(Ok((out, legacy_out))) => {
            assert_record_parity(label, &out, &legacy_out);
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

struct SoakParams {
    records: usize,
    seq_len: usize,
    memory_limit: usize,
    output_byte_limit: u64,
    block_batch: usize,
}

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
            // 96 KiB budget, with a tiny reorder window and block_batch == 1
            // (max per-block churn and the most `reader_eof`/`in_flight` events).
            SoakRegime::ManySmallFiles => SoakParams {
                records: 40_000,
                seq_len: 150,
                memory_limit: 96 * 1024,
                output_byte_limit: 128 * 1024,
                block_batch: 1,
            },
            // ~12k × 150B ≈ 1.8 MB spilled into ~4 large files at a 512 KiB
            // budget, with a roomy window and block_batch == 4.
            SoakRegime::FewLargeFiles => SoakParams {
                records: 12_000,
                seq_len: 150,
                memory_limit: 512 * 1024,
                output_byte_limit: 4 * 1024 * 1024,
                block_batch: 4,
            },
        }
    }
}

/// External-watchdog soak MATRIX for the Phase-2 decompress path. Crosses
/// pipeline-thread count × decompress granularity × spill regime, so both the
/// block-parallel reorder/in-flight/EOF protocol and the file-granularity FIFO
/// are hammered across {1, 2, 8} workers, {many small, few large} spill-file
/// shapes, and both code paths. Each (case × iteration) runs under a wall-clock
/// watchdog and is checked byte-for-byte against the legacy oracle, so a
/// livelock fails fast and any lost / duplicated / reordered block is caught.
///
/// `rstest` generates the full cross product (3 × 2 × 2 = 12 cases); nextest runs
/// them as independent parallel tests. The mandatory larger-than-budget run is
/// `ManySmallFiles` (~100 spill files at a 96 KiB budget). This is the real P5
/// hardening gate the OFF-default flip rests on; it complements the loom model
/// (exhaustive but tiny), the proptest (random sizes), and the single-corner
/// high-contention soak (12 threads, 1-block window).
#[rstest]
fn block_parallel_soak_matrix_matches_legacy(
    #[values(1, 2, 8)] pipeline_threads: usize,
    #[values(true, false)] file_granularity: bool,
    #[values(SoakRegime::ManySmallFiles, SoakRegime::FewLargeFiles)] regime: SoakRegime,
) {
    use std::time::Duration;

    const ITERATIONS: usize = 4;
    const SORTER_THREADS: usize = 2;
    const WATCHDOG: Duration = Duration::from_secs(120);

    let sort_order = SortOrder::Coordinate;
    let SoakParams { records, seq_len, memory_limit, output_byte_limit, block_batch } =
        regime.params();

    for iter in 0..ITERATIONS {
        // Distinct seed per (regime, threads, granularity, iter).
        let seed = 0x50A4_0000_u64
            .wrapping_add((iter as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15))
            .wrapping_add((pipeline_threads as u64) << 40)
            .wrapping_add(u64::from(file_granularity) << 32)
            ^ regime.seed_salt();
        let (header, recs) = synthesize_sized_records(records, seed, seq_len);
        let label = format!(
            "soak-{regime:?}-t{pipeline_threads}-fg{file_granularity}-bb{block_batch}-i{iter}"
        );
        run_watchdogged_parity(
            &label,
            sort_order,
            header,
            recs,
            memory_limit,
            SORTER_THREADS,
            pipeline_threads,
            output_byte_limit,
            SortDecompressTuning { file_granularity, block_batch },
            WATCHDOG,
        );
    }
}

fn synthesize_sized_records(n: usize, seed: u64, seq_len: usize) -> (Header, Vec<RawRecord>) {
    let header = Header::default();
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(1);
    let mut next_u32 = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        #[allow(clippy::cast_possible_truncation)]
        let v = state as u32;
        v
    };

    let mut records = Vec::with_capacity(n);
    for i in 0..n {
        let name = format!("rd_{}", next_u32() % 100_000);
        let pos: i32 = (next_u32() % 1_000_000).cast_signed();
        let is_paired = i % 2 == 0;
        let flags: u16 = 0x4 | if is_paired { 0x1 | 0x8 } else { 0 };
        let bytes = make_bam_bytes(-1, pos, flags, name.as_bytes(), &[], seq_len, -1, -1, &[]);
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
    // Default decompress tuning (block-parallel, block_batch 4) under a
    // per-case watchdog: the regression this pins is a deadlock at high
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
        SortDecompressTuning::default(),
        Duration::from_secs(90),
    );
}

#[test]
fn three_step_chain_empty_input_drains_cleanly() {
    let header = Header::default();
    let sorter = RawExternalSorter::new(SortOrder::Coordinate).memory_limit(1024 * 1024);
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
    let output_byte_limit = 1 << 20;
    let received: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let source = Phase2EventSource::new(events, output_byte_limit);
    let merge = SortMerge::<RecordBatchOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        output_byte_limit,
        256,
    );
    let sink = VecSink { received: Arc::clone(&received), kind: StepKind::Serial };

    let builder = Pipeline::builder();
    builder.chain(source).chain(merge).chain(sink).into_sink_marker();
    let pipeline = builder.build()?;
    pipeline.run(PipelineConfig { threads: 1, ..Default::default() })?;

    let collected = std::mem::take(&mut *received.lock());
    let mut out = Vec::new();
    for batch in collected {
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
    let received: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let source = Phase2EventSource::new(events, output_byte_limit);
    let merge = SortMerge::<RecordBatchOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        output_byte_limit,
        target_batch_count,
    );
    let sink = VecSink { received: Arc::clone(&received), kind: StepKind::Serial };

    let builder = Pipeline::builder();
    builder.chain(source).chain(merge).chain(sink).into_sink_marker();
    let pipeline = builder.build()?;
    pipeline.run(PipelineConfig { threads: 1, ..Default::default() })?;

    Ok(std::mem::take(&mut *received.lock()))
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
