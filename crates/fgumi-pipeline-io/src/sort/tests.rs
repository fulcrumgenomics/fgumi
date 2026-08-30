//! Tests for the runall-sort chains.
//!
//! The record-input chain is `SortBuffer` → `CompressSpill` →
//! `SortSpillDecompress` → `SortMerge`; the legacy `SortAndSpill` Phase-1 head
//! it replaced was retired in P7. The block-input arena front (`ReadBlocks` →
//! `InflateToArena` → `FindBoundariesAndSort`) is covered at the end of this
//! module.
//!
//! `RawExternalSorter::sort` (driven here via [`sort_via_legacy`]) is retained
//! as the parity oracle both chains are validated against.

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
use crate::sort::protocol::SortChunkEvent;
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
    use fgumi_sort::KeyTypesSpec;

    let records: Vec<RawRecord> = (0..8i32)
        .map(|i| {
            // Record 0 fixes the variant with MI 1; record 4 violates it.
            let mi = if i == 4 { 2 } else { 1 };
            RawRecord::from(record_with_mi(100 - i, format!("rd_{i}").as_bytes(), mi))
        })
        .collect();

    let sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
        .memory_limit(256 * 1024 * 1024)
        .threads(1)
        .key_types(KeyTypesSpec::None);
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
    let err = drive_sort_buffer_pipeline(
        sorter,
        &Header::default(),
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        1,
        StepKind::Serial,
        SortDecompressTuning::default(),
        SpillCodec::Zstd,
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

        // The two framings must also agree with EACH OTHER. Without this the
        // unstable branch only length-checks `batch_out`, so a framing divergence
        // between `SortMerge<BlockOutput>` and `SortMerge<RecordBatchOutput>`
        // that preserved record count would pass the queryname cases silently.
        let mut batch_sorted = batch_out.clone();
        batch_sorted.sort_unstable();
        assert_eq!(
            got, batch_sorted,
            "{sort_order:?} block vs RecordBatch path record multiset differs"
        );
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

// ── Production spill split + run extension (SpillGather run-former) ────────────
//
// `drive_sort_buffer_pipeline` above drives the retired monolithic `CompressSpill`
// for legacy Phase-2 coverage. The tests below drive the PRODUCTION spill split
// (`SortBuffer` → `SpillGather` → `SpillBlockCompress` → `SpillWrite` →
// `SortSpillDecompress` → `SortMerge`) that `fgumi sort` actually builds, so they
// exercise the run-forming `SpillGather` end-to-end: contiguous chunks coalesce
// into one run, and the merge sees one slot per run.

/// Drive the production spill split and return `(merged record bytes, run count)`.
/// The run count is the `SortMerge` stats `runs_written`, which equals the number
/// of spill files (`SpillReady`s) the merge saw — i.e. the number of runs the
/// `SpillGather` run-former formed. Mirrors `ChainBuilder::add_sort`'s wiring.
fn drive_spill_split_pipeline(
    sorter: RawExternalSorter,
    header: &Header,
    batches: Vec<RecordBatch>,
    output_byte_limit: u64,
    threads: usize,
) -> Result<(Vec<Vec<u8>>, usize)> {
    use fgumi_sort::TmpDirAllocator;

    let received: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let stats_slot: Arc<Mutex<Option<fgumi_sort::SortStats>>> = Arc::new(Mutex::new(None));
    let sort_order = sorter.sort_order();

    let dir = tempfile::TempDir::new()?;
    let alloc =
        TmpDirAllocator::with_probe(vec![dir.path().to_path_buf()], Box::new(|_| Ok(u64::MAX)), 0)?;
    let temp_dirs = Arc::new(vec![dir]);
    let codec = SpillCodec::Zstd;

    let source = VecSource::new(batches, output_byte_limit);
    let sort_buffer = SortBuffer::from_sorter(sorter, header, output_byte_limit)?;
    let gather = SpillGather::new(output_byte_limit);
    let compress = SpillBlockCompress::new(codec, 3, output_byte_limit);
    let write = SpillWrite::new(Arc::new(Mutex::new(alloc)), codec, output_byte_limit, temp_dirs);
    let decompress = SortSpillDecompress::new(output_byte_limit, SortDecompressTuning::default());
    let merge =
        SortMerge::<RecordBatchOutput>::with_target_batch_count(sort_order, output_byte_limit, 256)
            .with_stats_slot(Arc::clone(&stats_slot));
    let sink = VecSink { received: Arc::clone(&received), kind: StepKind::Serial };

    let builder = Pipeline::builder();
    builder
        .chain(source)
        .chain(sort_buffer)
        .chain(gather)
        .chain(compress)
        .chain(write)
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
    let runs = stats_slot.lock().take().expect("SortMerge published its stats").runs_written;
    Ok((out, runs))
}

/// Records sorted into `order`, obtained by running them through the oracle and
/// parsing its (byte-identical) sorted output back into `RawRecord`s. Feeding
/// this back in is the already-sorted-input case the run-former must coalesce.
fn presorted_records(order: SortOrder, header: &Header, records: &[RawRecord]) -> Vec<RawRecord> {
    sort_via_legacy(order, header, records, 256 * 1024 * 1024, 1)
        .expect("oracle sort")
        .into_iter()
        .map(RawRecord::from)
        .collect()
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

    let sorter = RawExternalSorter::new(order)
        .memory_limit(256 * 1024) // small ⇒ many contiguous chunks to coalesce
        .threads(2)
        .output_compression(1)
        .temp_compression(1);
    let (out, runs) =
        drive_spill_split_pipeline(sorter, &header, pack_batches(&sorted, 256), 4 * 1024 * 1024, 4)
            .expect("spill-split pipeline drives to completion");

    let expected = sort_via_legacy(order, &header, &sorted, 256 * 1024 * 1024, 1).expect("oracle");
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

    let sorter = RawExternalSorter::new(order)
        .memory_limit(256 * 1024)
        .threads(2)
        .output_compression(1)
        .temp_compression(1);
    let (out, runs) = drive_spill_split_pipeline(
        sorter,
        &header,
        pack_batches(&reversed, 256),
        4 * 1024 * 1024,
        4,
    )
    .expect("spill-split pipeline drives to completion");

    let expected =
        sort_via_legacy(order, &header, &reversed, 256 * 1024 * 1024, 1).expect("oracle");
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

    let sorter = RawExternalSorter::new(order)
        .memory_limit(256 * 1024)
        .threads(2)
        .output_compression(1)
        .temp_compression(1);
    let (out, runs) = drive_spill_split_pipeline(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        4,
    )
    .expect("spill-split pipeline drives to completion");

    let expected = sort_via_legacy(order, &header, &records, 256 * 1024 * 1024, 1).expect("oracle");
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

    let sorter = RawExternalSorter::new(order)
        .memory_limit(256 * 1024 * 1024) // above the whole input ⇒ no spill
        .threads(2)
        .output_compression(1)
        .temp_compression(1);
    let (out, runs) = drive_spill_split_pipeline(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        4,
    )
    .expect("spill-split pipeline drives to completion");

    let expected = sort_via_legacy(order, &header, &records, 256 * 1024 * 1024, 1).expect("oracle");
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

    let sorter = RawExternalSorter::new(order)
        .memory_limit(256 * 1024) // small ⇒ spills, exercising the merge at t1
        .threads(1)
        .output_compression(1)
        .temp_compression(1);
    let (out, runs) = drive_spill_split_pipeline(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        1,
    )
    .expect("spill-split pipeline drives to completion at one thread");

    let expected = sort_via_legacy(order, &header, &records, 256 * 1024 * 1024, 1).expect("oracle");
    assert_eq!(out, expected, "{order:?}: single-thread scheduled sort diverged from the oracle");
    // The small budget must actually spill; otherwise this degrades to an
    // in-memory sort and no longer exercises the merge path it names.
    assert!(runs >= 1, "{order:?}: the small budget must spill (got {runs} run(s))");
}

/// Block-parallel decompression completes out of order (workers decompress one
/// file's blocks concurrently), yet the reassembled output must be byte-
/// identical to the in-order (file-granularity) result. Property test over a
/// range of record counts and `block_batch` sizes and a high pipeline-thread
/// count (more concurrent decompressors ⇒ more out-of-order completion). A
/// straggler worker hitting reader-EOF while another holds an in-flight block
/// must not truncate the output (record count is asserted equal).
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
            SortDecompressTuning { file_granularity: false, block_batch: BLOCK_BATCH },
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
/// legacy oracle proves nothing is lost, duplicated, or reordered under churn.
///
/// NOTE: this drives `drive_sort_pipeline_tuned` (the record-input `SortBuffer`
/// chain), NOT the block-input arena front. The front's capacity-1 arena cycle —
/// `ReadBlocks` acquire+admit and `FindBoundariesAndSort` seal+free on the
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
    block_batch: usize,
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

#[cfg(feature = "stress-tests")]
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

    const SORTER_THREADS: usize = 2;
    const WATCHDOG: Duration = Duration::from_secs(120);

    let iterations = soak_iterations(4, 1);
    let sort_order = SortOrder::Coordinate;
    let SoakParams { records, seq_len, memory_limit, output_byte_limit, block_batch } =
        regime.params();

    for iter in 0..iterations {
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

/// Number of reference sequences the synthetic header declares.
const N_TEST_REFS: usize = 4;

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
    let header = {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;
        let len = NonZeroUsize::new(1_000_000).expect("nonzero");
        let mut builder = Header::builder();
        for i in 0..N_TEST_REFS {
            builder = builder
                .add_reference_sequence(format!("chr{i}"), Map::<ReferenceSequence>::new(len));
        }
        builder.build()
    };
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

// ── Arena block-input front (ReadBlocks → InflateToArena → FindBoundariesAndSort) ──

/// `Exclusive` source draining a `Vec<BgzfBlock>` one block per `try_run`.
struct BgzfBlockSource {
    blocks: Vec<crate::types::BgzfBlock>,
    held: HeldSlot<Unpushed<crate::types::BgzfBlock>>,
    output_byte_limit: u64,
}

impl BgzfBlockSource {
    fn new(mut blocks: Vec<crate::types::BgzfBlock>, output_byte_limit: u64) -> Self {
        blocks.reverse();
        Self { blocks, held: HeldSlot::new(), output_byte_limit }
    }
}

impl Step for BgzfBlockSource {
    type Input = ();
    type Outputs = OrderedBytesSingle<crate::types::BgzfBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "BgzfBlockSource",
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
                match event {
                    SortChunkEvent::Spill { chunk, .. }
                    | SortChunkEvent::Residual { chunk, .. } => {
                        out.chunks.push(
                            (0..chunk.len()).map(|i| chunk.record_bytes(i).to_vec()).collect(),
                        );
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
                bytes: blocks.remove(0).data,
                uncompressed_size: u32::try_from(payload.len()).expect("payload fits u32"),
                index: None,
            }
        })
        .collect()
}

/// Drive `BgzfBlockSource → ReadBlocks → InflateToArena → FindBoundariesAndSort`
/// and return every emitted chunk's record bodies, chunk by chunk, plus the
/// terminal `AllAnnounced` counts.
fn drive_arena_front(
    blocks: Vec<crate::types::BgzfBlock>,
    n_ref: u32,
    memory_limit: usize,
    output_byte_limit: u64,
    threads: usize,
) -> Result<ArenaFrontOutput> {
    let received: Arc<Mutex<ArenaFrontOutput>> = Arc::new(Mutex::new(ArenaFrontOutput::default()));

    let builder = Pipeline::builder();
    builder
        .chain(BgzfBlockSource::new(blocks, output_byte_limit))
        .chain(ReadBlocks::new(memory_limit, output_byte_limit))
        .chain(InflateToArena::new(output_byte_limit))
        .chain(FindBoundariesAndSort::new(CoordinateStrategy::new(n_ref), 1, output_byte_limit))
        .chain(ChunkEventSink { received: Arc::clone(&received) })
        .into_sink_marker();
    let pipeline = builder.build()?;
    pipeline.run(PipelineConfig { threads, ..Default::default() })?;

    Ok(std::mem::take(&mut *received.lock()))
}

/// End-to-end parity for the block-input arena front, driven through the real
/// runtime rather than by poking `ingest_block`/`finalize` directly. A
/// `memory_limit` above the whole input yields exactly ONE run, so the single
/// residual chunk must be byte-identical to the legacy oracle's coordinate sort.
///
/// `payload_bytes` varies where records fall relative to block boundaries: the
/// small case guarantees records straddle blocks (the front's carry path), the
/// large case puts the whole stream in one block.
#[rstest]
#[case::straddling_blocks(4096, 20)]
#[case::one_block_per_run(60_000, 2)]
fn arena_front_chain_single_run_matches_legacy_oracle(
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

    let legacy_out =
        sort_via_legacy(SortOrder::Coordinate, &header, &records, 256 * 1024 * 1024, 1)
            .expect("legacy oracle");
    assert_eq!(out.chunks[0].len(), legacy_out.len(), "record count mismatch");
    for (i, (got, want)) in out.chunks[0].iter().zip(legacy_out.iter()).enumerate() {
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

/// Drive the full BAM sort path — the benchmark path — end to end:
/// `BgzfBlockSource → ReadBlocks → InflateToArena → FindBoundariesAndSort →
/// SpillGather → SpillBlockCompress → SpillWrite → SortSpillDecompress →
/// SortMerge`. Returns `(merged record bytes, run count)`. This is the head the
/// 43 GB `1kg-wgs` benchmark actually uses (SAM/`SortBuffer` is the other head),
/// so the run-former must coalesce here, not only on the `SortBuffer` path.
fn drive_arena_split_pipeline(
    blocks: Vec<crate::types::BgzfBlock>,
    n_ref: u32,
    memory_limit: usize,
    output_byte_limit: u64,
    threads: usize,
) -> Result<(Vec<Vec<u8>>, usize)> {
    use fgumi_sort::TmpDirAllocator;

    let received: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let stats_slot: Arc<Mutex<Option<fgumi_sort::SortStats>>> = Arc::new(Mutex::new(None));

    let dir = tempfile::TempDir::new()?;
    let alloc =
        TmpDirAllocator::with_probe(vec![dir.path().to_path_buf()], Box::new(|_| Ok(u64::MAX)), 0)?;
    let temp_dirs = Arc::new(vec![dir]);
    let codec = SpillCodec::Zstd;

    let merge = SortMerge::<RecordBatchOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        output_byte_limit,
        256,
    )
    .with_stats_slot(Arc::clone(&stats_slot));

    let builder = Pipeline::builder();
    builder
        .chain(BgzfBlockSource::new(blocks, output_byte_limit))
        .chain(ReadBlocks::new(memory_limit, output_byte_limit))
        .chain(InflateToArena::new(output_byte_limit))
        .chain(FindBoundariesAndSort::new(CoordinateStrategy::new(n_ref), 1, output_byte_limit))
        .chain(SpillGather::new(output_byte_limit))
        .chain(SpillBlockCompress::new(codec, 3, output_byte_limit))
        .chain(SpillWrite::new(Arc::new(Mutex::new(alloc)), codec, output_byte_limit, temp_dirs))
        .chain(SortSpillDecompress::new(output_byte_limit, SortDecompressTuning::default()))
        .chain(merge)
        .chain(VecSink { received: Arc::clone(&received), kind: StepKind::Serial })
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
    let runs = stats_slot.lock().take().expect("SortMerge published its stats").runs_written;
    Ok((out, runs))
}

/// The benchmark case, in miniature: a pre-sorted BAM re-sorted through the full
/// arena/BAM path under a small memory budget must coalesce every contiguous
/// sealed run into ONE spill file, and stay byte-identical to the oracle. This is
/// the `FindBoundariesAndSort` head (not `SortBuffer`), so it proves the
/// run-former covers the path the 43 GB `1kg-wgs` sort exercises.
#[test]
fn presorted_bam_coalesces_to_one_run_through_the_arena_split() {
    let (header, records) = synthesize_sized_records(8_000, 0xB00C_0001, 40);
    let n_ref = u32::try_from(header.reference_sequences().len()).expect("n_ref fits u32");
    let sorted = presorted_records(SortOrder::Coordinate, &header, &records);

    // Small per-block payload ⇒ many BGZF blocks; small memory ⇒ ReadBlocks seals
    // several contiguous runs the run-former must coalesce.
    let blocks = bgzf_blocks_for(&sorted, n_ref, 4096);
    let (out, runs) = drive_arena_split_pipeline(blocks, n_ref, 64 * 1024, 4 * 1024 * 1024, 4)
        .expect("arena-split pipeline drives to completion");

    let expected = sort_via_legacy(SortOrder::Coordinate, &header, &sorted, 256 * 1024 * 1024, 1)
        .expect("oracle");
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
    let (out, runs) = drive_arena_split_pipeline(blocks, n_ref, 64 * 1024, 4 * 1024 * 1024, 4)
        .expect("arena-split pipeline drives to completion");

    let expected = sort_via_legacy(SortOrder::Coordinate, &header, &records, 256 * 1024 * 1024, 1)
        .expect("oracle");
    assert_eq!(out, expected, "arena/BAM-path shuffled output diverged from the oracle");
    // Shuffled BAM input under a tiny per-block budget must open more than one
    // run (partial coalescing), otherwise the case no longer exercises it.
    assert!(runs >= 2, "arena/BAM shuffled input must form >=2 runs (got {runs})");
}
