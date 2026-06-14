//! Tests for the runall-sort three-step chain
//! (`SortAndSpill` → `SortSpillDecompress` → `SortMerge`) and for the
//! single-step `SortBamFile` wrapper.
//!
//! Strategy: build a 5-step pipeline
//! `VecSource → SortAndSpill → SortSpillDecompress → SortMerge → VecSink`
//! and compare emitted record bytes against `RawExternalSorter::sort()`
//! byte-for-byte. The sort engine itself is covered by 360+ tests in
//! `fgumi-sort`; these tests pin the typed-step adapter (push, drain,
//! batch, retry, ordinal-restamp, inter-step event flow).
//!
//! We exercise both the in-memory fast path (`memory_limit` > total) and
//! the multi-spill / k-way merge path (`memory_limit` << total).

use std::io;
use std::sync::Arc;

use anyhow::Result;
use fgumi_raw_bam::RawRecord;
use fgumi_raw_bam::testutil::make_bam_bytes;
use fgumi_sort::{QuerynameComparator, RawExternalSorter, SortOrder};
use noodles::sam::Header;
use parking_lot::Mutex;
use rstest::rstest;

use super::*;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepProfile};
use crate::pipeline::steps::types::RecordBatch;

// ── In-memory source / sink test steps ──────────────────────────────────────

/// `Exclusive` source that drains a `Vec<RecordBatch>` one batch per
/// `try_run` call. Returns `Finished` when empty. Output ordering is
/// `ByOrdinal` so the framework's reorder stage exercises the same path
/// real producers go through.
struct VecSource {
    batches: Vec<RecordBatch>,
    held: HeldSlot<Unpushed<RecordBatch>>,
    output_byte_limit: u64,
}

impl VecSource {
    fn new(mut batches: Vec<RecordBatch>, output_byte_limit: u64) -> Self {
        // Pop from the end, so reverse to preserve user-supplied order.
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

/// Sink that appends every received batch into a shared `Vec<RecordBatch>`,
/// for byte-level equivalence checks. `kind` selects the step kind: most tests
/// use `Exclusive`; the `threads == 1` deadlock repro uses `Serial` so the
/// chain has only ONE `Exclusive` step (`VecSource`) and is allowed to run at
/// a single worker.
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

// ── Synthetic-record helpers ────────────────────────────────────────────────

/// Build `n` synthetic unmapped records with deterministic, sort-relevant
/// variation (read names, tid/pos, paired flags). Thin wrapper over
/// [`synthesize_sized_records`] with a zero-length sequence — each record is
/// the minimal ~50-byte header, so spills are ≤ 1 BGZF block.
fn synthesize_records(n: usize, seed: u64) -> (Header, Vec<RawRecord>) {
    synthesize_sized_records(n, seed, 0)
}

/// Pack records into batches of `batch_size` records each.
fn pack_batches(records: &[RawRecord], batch_size: usize) -> Vec<RecordBatch> {
    use crate::pipeline::steps::types::RecordBatchBuilder;
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

/// Drive the pipeline:
/// `VecSource(batches)` → `SortAndSpill` → `SortSpillDecompress` →
/// `SortMerge` → `VecSink(sink_kind)` → collected.
/// Returns the per-record bytes in emitted order. `sink_kind` is
/// `Exclusive` for the parity tests and `Serial` for the `threads == 1`
/// repro (so the chain has a single `Exclusive` step and one worker is legal).
fn drive_sort_pipeline(
    sorter: RawExternalSorter,
    header: &Header,
    batches: Vec<RecordBatch>,
    output_byte_limit: u64,
    threads: usize,
    sink_kind: StepKind,
) -> Result<Vec<Vec<u8>>> {
    let received: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));

    let sort_order = sorter.sort_order();
    let source = VecSource::new(batches, output_byte_limit);
    let and_spill = SortAndSpill::from_sorter(sorter, header, 32)?;
    // 64 events × ~256 KiB/block ≈ 16 MiB queue ceiling for spill blocks.
    let decompress = SortSpillDecompress::new(64);
    let merge = SortMerge::with_target_batch_count(sort_order, output_byte_limit, 256);
    let sink = VecSink { received: Arc::clone(&received), kind: sink_kind };

    let builder = Pipeline::builder();
    builder
        .chain(source)
        .chain(and_spill)
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

// ── Reference: RawExternalSorter::sort to bytes ─────────────────────────────

/// Run `RawExternalSorter::sort()` on `records` and return per-record
/// bytes in the resulting BAM.
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

#[rstest]
#[case::coordinate(SortOrder::Coordinate, 2)]
#[case::coordinate_t4(SortOrder::Coordinate, 4)]
#[case::queryname_lex(SortOrder::Queryname(QuerynameComparator::Lexicographic), 2)]
#[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural), 2)]
#[case::template_coordinate(SortOrder::TemplateCoordinate, 2)]
fn three_step_chain_in_memory_path_matches_legacy(
    #[case] sort_order: SortOrder,
    #[case] threads: usize,
) {
    // 5K small records, memory_limit large enough to keep everything
    // in memory — exercises the in-memory residual chunk (no spill files)
    // path of the new chain.
    let (header, records) = synthesize_records(5_000, 0x00C0_FFEE);
    let memory_limit = 256 * 1024 * 1024;

    let sorter = RawExternalSorter::new(sort_order)
        .memory_limit(memory_limit)
        .threads(threads)
        .output_compression(1)
        .temp_compression(1);
    let new_out = drive_sort_pipeline(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        threads.max(3),
        StepKind::Exclusive,
    )
    .expect("pipeline drives to completion");

    let legacy_out =
        sort_via_legacy(sort_order, &header, &records, memory_limit, threads).expect("legacy");

    assert_eq!(new_out.len(), legacy_out.len(), "{sort_order:?} record count mismatch");
    for (i, (got, want)) in new_out.iter().zip(legacy_out.iter()).enumerate() {
        assert_eq!(got, want, "{sort_order:?} record {i} bytes differ");
    }
}

#[rstest]
#[case::coordinate(SortOrder::Coordinate, 2)]
#[case::queryname_lex(SortOrder::Queryname(QuerynameComparator::Lexicographic), 2)]
#[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural), 2)]
#[case::template_coordinate(SortOrder::TemplateCoordinate, 2)]
fn three_step_chain_multi_spill_path_matches_legacy(
    #[case] sort_order: SortOrder,
    #[case] threads: usize,
) {
    // 20K records under a tight 256 KB memory_limit → many spills,
    // exercising the spill + decompress + k-way merge path of the new
    // chain. This is the hot path that motivated the split.
    let (header, records) = synthesize_records(20_000, 0xFEED_FACE);
    let memory_limit = 256 * 1024;

    let sorter = RawExternalSorter::new(sort_order)
        .memory_limit(memory_limit)
        .threads(threads)
        .output_compression(1)
        .temp_compression(1);
    let new_out = drive_sort_pipeline(
        sorter,
        &header,
        pack_batches(&records, 256),
        4 * 1024 * 1024,
        threads.max(3),
        StepKind::Exclusive,
    )
    .expect("pipeline drives to completion");

    let legacy_out =
        sort_via_legacy(sort_order, &header, &records, memory_limit, threads).expect("legacy");

    assert_eq!(new_out.len(), legacy_out.len(), "{sort_order:?} record count mismatch");
    for (i, (got, want)) in new_out.iter().zip(legacy_out.iter()).enumerate() {
        assert_eq!(got, want, "{sort_order:?} record {i} bytes differ");
    }
}

/// Like [`synthesize_records`] but each record carries a `seq_len`-base
/// sequence payload so the record byte size is tunable. Large records make
/// each spill chunk span MANY BGZF blocks (more than `PHASE2_DECOMP_CAP`),
/// which is the production condition the small-record tests above miss: a
/// slot does not reach `queue_eof` on its first cap-fill, so the merge
/// consumer drains it to empty and re-parks on `block_ready` while the
/// producer is still feeding the same slot.
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

/// Regression test for issue #330 BUG #3 — the fused-sort-chain deadlock.
///
/// Production symptom (fgumi-benchmarks, c7g.4xlarge, `--threads 8`, real
/// CODEC BAM): `runall --start-from sort` reads/sorts/spills then wedges
/// forever at the sort→group handoff. Root cause: the merge consumer
/// blocked on `slot.block_ready.wait()` whenever a slot's decompressed queue
/// was momentarily empty and not yet at EOF, parking a framework worker (and,
/// being `Serial`, holding the step lock) with no guarantee a sibling worker
/// would refill the slot.
///
/// The condition the small-record `*_multi_spill_*` tests miss is **large
/// spill files** — each spill chunk spans far more than `PHASE2_DECOMP_CAP`
/// (8) BGZF blocks, so a slot stays `!queue_eof` while the merge consumer
/// drains it to empty mid-merge. The fix makes the merge non-blocking and
/// resumable (`MergeDriver::try_step` → `SortMerge` returns `Contention`
/// instead of parking); this test asserts the chain completes at every worker
/// count, deterministically (a hang ⇒ the bug is back).
///
/// `threads == 1` is the M1 case the fix targets (a sole worker that parks in
/// the merge can never round-robin to refill). It runs here with a large
/// output-queue limit on purpose: at `threads == 1` the round-robin's
/// priority-restart separately starves the single worker's downstream drain
/// under output backpressure — a pre-existing, `threads == 1`-only artifact of
/// this minimal harness. It does not affect production: the real fused chain
/// has two `Exclusive` steps (BAM reader + writer) and so always runs at
/// `threads >= 2`, where a separate worker drains the output. Sizing the queue
/// above the total output isolates the merge-liveness property this test
/// guards (the old blocking merge hangs at `threads == 1` regardless of queue
/// size). The non-blocking-merge unit coverage is in `fgumi-sort`
/// (`from_slots_try_step_*`).
#[rstest]
#[case::t1(1)]
#[case::t2(2)]
#[case::t4(4)]
#[case::t8(8)]
fn three_step_chain_large_spill_completes(#[case] pipeline_threads: usize) {
    use std::sync::mpsc;
    use std::time::Duration;

    // ~320-byte records (seq_len=200 → 200 + 100 seq/qual bytes + header)
    // under a 1 MiB memory_limit ⇒ ~3K records/spill ⇒ ~1 MiB ⇒ ~15 BGZF
    // blocks/spill (> CAP 8). 60K records ⇒ ~20 spill files. The
    // many-blocks-per-spill shape is what makes the merge consumer drain a
    // slot to empty while that slot is still `!queue_eof` — the production
    // condition that triggered the original `block_ready` park.
    let total_records = 60_000usize;
    let memory_limit = 1024 * 1024;
    let sort_order = SortOrder::Coordinate;
    // Sort-engine internal threads are independent of the pipeline worker
    // count; keep them low so the spill shape is stable across cases.
    let sorter_threads = 2;

    let (header, records) = synthesize_sized_records(total_records, 0xBADD_CAFE, 200);
    let batches = pack_batches(&records, 256);
    let sorter = RawExternalSorter::new(sort_order)
        .memory_limit(memory_limit)
        .threads(sorter_threads)
        .output_compression(1)
        .temp_compression(1);

    // Reference output (byte-for-byte) so a hang ISN'T the only thing this
    // test catches — a merge/refill ordering regression must also fail.
    let legacy_out = sort_via_legacy(sort_order, &header, &records, memory_limit, sorter_threads)
        .expect("legacy sort");

    // Output-queue ceiling sized above the ~19 MiB total output so the
    // `threads == 1` case isn't confounded by the unrelated single-worker
    // output-backpressure starvation (see the test doc). The old blocking
    // merge hangs at `threads == 1` regardless of this size.
    let output_queue_limit = 256 * 1024 * 1024;
    let (tx, rx) = mpsc::channel();
    let worker = std::thread::Builder::new()
        .name(format!("bug3-repro-t{pipeline_threads}"))
        .spawn(move || {
            // Drive the chain directly at the requested pipeline worker count
            // — NOT `threads.max(3)` — with a `Serial` sink so the chain has a
            // single `Exclusive` step and `threads == 1` is legal, exercising
            // single-worker merge liveness rather than hiding it.
            let result = drive_sort_pipeline(
                sorter,
                &header,
                batches,
                output_queue_limit,
                pipeline_threads,
                StepKind::Serial,
            );
            let _ = tx.send(result);
        })
        .expect("spawn watchdog worker");

    match rx.recv_timeout(Duration::from_secs(90)) {
        Ok(Ok(out)) => {
            assert_eq!(
                out.len(),
                total_records,
                "record count mismatch (threads={pipeline_threads})"
            );
            assert_eq!(
                out, legacy_out,
                "sorted output diverges from legacy (threads={pipeline_threads})"
            );
            worker.join().expect("worker thread panicked");
        }
        Ok(Err(e)) => panic!("pipeline errored (threads={pipeline_threads}): {e:#}"),
        Err(mpsc::RecvTimeoutError::Timeout) => panic!(
            "DEADLOCK: fused sort chain did not complete within 90s \
             (pipeline_threads={pipeline_threads}, {total_records} records, \
             memory_limit={memory_limit})"
        ),
        Err(mpsc::RecvTimeoutError::Disconnected) => {
            panic!("watchdog worker dropped its sender (panicked?) at threads={pipeline_threads}")
        }
    }
}

#[test]
fn three_step_chain_empty_input_drains_cleanly() {
    let header = Header::default();
    let sorter = RawExternalSorter::new(SortOrder::Coordinate).memory_limit(1024 * 1024);
    let out = drive_sort_pipeline(sorter, &header, Vec::new(), 64 * 1024, 3, StepKind::Exclusive)
        .expect("empty pipeline");
    assert!(out.is_empty());
}

// ── SortBamFile (Exclusive single-step) tests ──────────────────────────────

/// Smoke test for the `SortBamFile` step. Constructs a synthetic input
/// BAM, sorts it via `[SortBamFile]` through the pipeline, then sorts
/// the same input via `RawExternalSorter::sort()` directly, and
/// verifies byte-identical output. Since `SortBamFile`'s body is just
/// `sorter.sort(input, output)`, this test mostly verifies the
/// framework wrapping doesn't corrupt anything — the heavy testing is
/// in `fgumi-sort`'s own 360+ tests.
#[test]
fn sort_bam_file_matches_legacy_sort() {
    let (header, records) = synthesize_records(2_000, 0xCAFE_F00D);

    let dir = tempfile::tempdir().expect("tempdir");
    let input = dir.path().join("input.bam");
    let legacy_out = dir.path().join("legacy.bam");
    let pipeline_out = dir.path().join("pipeline.bam");

    {
        let mut writer =
            fgumi_bam_io::create_raw_bam_writer(&input, &header, 1, 1).expect("write input");
        for r in &records {
            writer.write_raw_record(r.as_ref()).expect("write rec");
        }
        writer.finish().expect("finish");
    }

    RawExternalSorter::new(SortOrder::Coordinate)
        .memory_limit(64 * 1024 * 1024)
        .threads(1)
        .output_compression(1)
        .temp_compression(1)
        .sort(&input, &legacy_out)
        .expect("legacy sort");

    let sorter = RawExternalSorter::new(SortOrder::Coordinate)
        .memory_limit(64 * 1024 * 1024)
        .threads(1)
        .output_compression(1)
        .temp_compression(1);
    let step = SortBamFile::new(
        sorter,
        input.clone(),
        pipeline_out.clone(),
        Arc::new(parking_lot::Mutex::new(None)),
    );

    let builder = Pipeline::builder();
    builder.chain(step).into_sink_marker();
    let pipeline = builder.build().expect("Pipeline::build");
    pipeline.run(PipelineConfig { threads: 1, ..Default::default() }).expect("Pipeline::run");

    let legacy_records = read_all_records(&legacy_out);
    let pipeline_records = read_all_records(&pipeline_out);
    assert_eq!(legacy_records.len(), pipeline_records.len(), "record count mismatch");
    for (i, (a, b)) in legacy_records.iter().zip(pipeline_records.iter()).enumerate() {
        assert_eq!(a, b, "record {i} bytes differ");
    }
}

#[test]
fn sort_bam_file_profile_is_exclusive() {
    let dir = tempfile::tempdir().expect("tempdir");
    let sorter = RawExternalSorter::new(SortOrder::Coordinate);
    let step = SortBamFile::new(
        sorter,
        dir.path().join("in.bam"),
        dir.path().join("out.bam"),
        Arc::new(parking_lot::Mutex::new(None)),
    );
    let p = step.profile();
    assert_eq!(p.name, "SortBamFile");
    assert_eq!(p.kind, StepKind::Exclusive);
    assert!(!p.sticky);
    assert!(p.output_queues.is_empty());
    assert!(p.branch_ordering.is_empty());
}

/// Read every record's bytes from a BAM file.
fn read_all_records(path: &std::path::Path) -> Vec<Vec<u8>> {
    let (mut reader, _hdr) = fgumi_bam_io::create_raw_bam_reader_with_opts(
        path,
        1,
        fgumi_bam_io::PipelineReaderOpts::default(),
    )
    .expect("create reader");
    let mut out = Vec::new();
    let mut rec = RawRecord::default();
    loop {
        match reader.read_record(&mut rec).expect("read") {
            0 => break,
            _ => out.push(rec.as_ref().to_vec()),
        }
    }
    out
}
