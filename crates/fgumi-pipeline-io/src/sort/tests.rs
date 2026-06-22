//! Tests for the runall-sort three-step chain
//! (`SortAndSpill` → `SortSpillDecompress` → `SortMerge`) and for the
//! single-step `SortBamFile` wrapper.

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
use crate::types::RecordBatch;
use fgumi_pipeline_core::{
    Unpushed,
    builder::{Pipeline, PipelineConfig},
    held::HeldSlot,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Step, StepCtx, StepProfile},
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
    let received: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));

    let sort_order = sorter.sort_order();
    let source = VecSource::new(batches, output_byte_limit);
    let and_spill = SortAndSpill::from_sorter(sorter, header, output_byte_limit)?;
    let decompress = SortSpillDecompress::new(output_byte_limit);
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
    use std::sync::mpsc;
    use std::time::Duration;

    let total_records = 60_000usize;
    let memory_limit = 1024 * 1024;
    let sort_order = SortOrder::Coordinate;
    let sorter_threads = 2;

    let (header, records) = synthesize_sized_records(total_records, 0xBADD_CAFE, 200);
    let batches = pack_batches(&records, 256);
    let sorter = RawExternalSorter::new(sort_order)
        .memory_limit(memory_limit)
        .threads(sorter_threads)
        .output_compression(1)
        .temp_compression(1);

    let legacy_out = sort_via_legacy(sort_order, &header, &records, memory_limit, sorter_threads)
        .expect("legacy sort");

    let output_queue_limit = 256 * 1024 * 1024;
    let (tx, rx) = mpsc::channel();
    let worker = std::thread::Builder::new()
        .name(format!("bug3-repro-t{pipeline_threads}"))
        .spawn(move || {
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
    let stats_slot = Arc::new(parking_lot::Mutex::new(None));
    let step =
        SortBamFile::new(sorter, input.clone(), pipeline_out.clone(), Arc::clone(&stats_slot));

    let builder = Pipeline::builder();
    builder.chain(step).into_sink_marker();
    let pipeline = builder.build().expect("Pipeline::build");
    pipeline.run(PipelineConfig { threads: 1, ..Default::default() }).expect("Pipeline::run");

    // `SortBamFile::try_run` must publish its stats into the shared slot; the
    // `SortFinalizeHook` reads this slot after `Pipeline::run` to log the summary.
    let stats = stats_slot.lock().take().expect("SortBamFile should publish SortStats");
    assert_eq!(stats.total_records, records.len() as u64);
    assert_eq!(stats.output_records, records.len() as u64);

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
    let merge = SortMerge::with_target_batch_count(SortOrder::Coordinate, output_byte_limit, 256);
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
