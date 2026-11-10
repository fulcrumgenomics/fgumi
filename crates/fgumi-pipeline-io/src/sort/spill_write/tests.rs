//! Unit tests for `SpillWrite`'s `StepCtx`-free core: per-file demux, codec
//! magic/trailer bracketing, runs held closed until `AllAnnounced`, the
//! consolidation-job lifecycle under `--max-temp-files`, and event mapping.
//! The jobs themselves run through the real tail in `consolidate_runs`' tests;
//! byte-exact readback of a consolidated sort is gated end-to-end by the sort
//! command's consolidation parity tests.

use super::*;
use crate::sort::merge_mode::RangeEligibility;
use crate::sort::protocol::MemoryChunkErased;
use crate::sort::spill_fixtures::{compress, keyed_run, write_run_file};
use fgumi_pipeline_core::testing::StepProbe;
use fgumi_sort::{InMemoryChunk, RawCoordinateKey};
use rstest::rstest;
use std::time::{Duration, Instant};
use tempfile::TempDir;

/// Build a `SpillWrite` writing into a fresh temp dir; returns it plus the dir
/// (kept alive by the caller) so written files can be inspected.
fn make_writer(codec: SpillCodec) -> (SpillWrite, TempDir) {
    let dir = TempDir::new().unwrap();
    let alloc = TmpDirAllocator::new(vec![dir.path().to_path_buf()]).unwrap();
    let writer = SpillWrite::new(Arc::new(Mutex::new(alloc)), codec, 1 << 20, Arc::new(Vec::new()));
    (writer, dir)
}

fn block(codec: SpillCodec, file_id: u32, is_last: bool, raw: &[u8]) -> SpillBlockEvent {
    SpillBlockEvent::Block {
        ordinal: 0,
        file_id,
        key_kind: SpillKeyKind::Coordinate,
        is_last_in_file: is_last,
        records_ingested_so_far: 42,
        bytes: compress(codec, raw),
        frames: None,
        index_stride: 1,
    }
}

fn announced(slot_count: u32, memory_chunk_count: u32, total_records: u64) -> SpillBlockEvent {
    SpillBlockEvent::AllAnnounced {
        ordinal: 0,
        slot_count,
        memory_chunk_count,
        total_records,
        index_stride: 1,
    }
}

#[rstest]
#[case(SpillCodec::Zstd)]
#[case(SpillCodec::Bgzf)]
fn a_closed_run_is_announced_only_at_all_announced(#[case] codec: SpillCodec) {
    let (mut w, dir) = make_writer(codec);
    // First (non-last) block: file opens, no event.
    w.process_event(block(codec, 5, false, &[1u8; 32])).unwrap();
    assert!(w.current.is_some(), "file must be open after first block ({codec:?})");

    // Last block: the file is finished and closed, but no slot is opened yet —
    // the merge cannot start before `AllAnnounced`, and an early slot would hold
    // a descriptor and read-ahead for the rest of the sort.
    w.process_event(block(codec, 5, true, &[2u8; 32])).unwrap();
    assert!(w.current.is_none(), "file closed after last block ({codec:?})");
    assert!(w.outbox.is_empty(), "no event until AllAnnounced ({codec:?})");

    w.process_event(announced(1, 0, 42)).unwrap();
    let Some(SortPhase1Event::SpillReady { slot, path, records_ingested_so_far }) =
        w.outbox.pop_front()
    else {
        panic!("expected SpillReady ({codec:?})");
    };
    assert_eq!(slot.file_id, 5, "slot file_id == logical seq ({codec:?})");
    assert_eq!(records_ingested_so_far, 42);
    assert!(path.starts_with(dir.path()), "spill file under temp dir ({codec:?})");
    assert_eq!(slot.codec, codec, "codec detected from written magic ({codec:?})");
    assert!(matches!(
        w.outbox.pop_front(),
        Some(SortPhase1Event::AllAnnounced {
            slot_count: 1,
            memory_chunk_count: 0,
            total_records: 42
        })
    ));
    assert!(w.outbox.is_empty());
}

#[test]
fn runs_are_announced_in_file_id_order_with_distinct_files() {
    let codec = SpillCodec::Zstd;
    let (mut w, _dir) = make_writer(codec);
    for file_id in [0, 1, 4] {
        w.process_event(block(codec, file_id, true, &[7u8; 16])).unwrap();
    }
    w.process_event(announced(3, 0, 3)).unwrap();
    let mut ids = Vec::new();
    let mut paths = Vec::new();
    while let Some(event) = w.outbox.pop_front() {
        match event {
            SortPhase1Event::SpillReady { slot, path, .. } => {
                ids.push(slot.file_id);
                paths.push(path);
            }
            SortPhase1Event::AllAnnounced { slot_count, .. } => assert_eq!(slot_count, 3),
            SortPhase1Event::MemoryChunk { .. }
            | SortPhase1Event::ConsolidationStart { .. }
            | SortPhase1Event::ConsolidationRetire { .. } => panic!("only runs were announced"),
        }
    }
    assert_eq!(ids, vec![0, 1, 4]);
    paths.dedup();
    assert_eq!(paths.len(), 3, "distinct file_ids must yield distinct paths");
}

#[test]
fn an_announced_run_count_that_disagrees_with_the_runs_written_fails_closed() {
    let codec = SpillCodec::Zstd;
    let (mut w, _dir) = make_writer(codec);
    w.process_event(block(codec, 0, true, &[7u8; 16])).unwrap();
    // Upstream claims two runs; only one was written — a run was lost.
    match w.process_event(announced(2, 0, 2)) {
        Err(err) => assert!(err.to_string().contains("announced 2"), "got: {err}"),
        Ok(()) => panic!("a run-count mismatch must fail closed"),
    }
}

#[test]
fn runs_with_different_key_types_fail_closed() {
    let codec = SpillCodec::Zstd;
    let (mut w, _dir) = make_writer(codec);
    w.process_event(block(codec, 0, true, &[7u8; 16])).unwrap();
    // A second run keyed differently could not be consolidated with the first.
    match w.process_event(keyed_run(codec, 1, 0..3)) {
        Err(err) => assert!(err.to_string().contains("keyed"), "got: {err}"),
        Ok(()) => panic!("mixed key types must fail closed"),
    }
}

/// Runs are held until `AllAnnounced`; an input that drains without it would
/// otherwise finish "successfully" having delivered none of them.
#[test]
fn draining_with_unannounced_runs_fails_closed() {
    let codec = SpillCodec::Zstd;
    let (mut w, _dir) = make_writer(codec);
    w.ensure_runs_announced().expect("a writer that closed no runs has nothing to strand");
    w.process_event(keyed_run(codec, 0, 0..3)).unwrap();
    let err = w.ensure_runs_announced().expect_err("a closed, unannounced run must fail");
    assert!(err.to_string().contains("never announced"), "got: {err}");

    w.process_event(announced(1, 0, 3)).unwrap();
    w.ensure_runs_announced().expect("announced runs are not stranded");
}

/// Everything after `AllAnnounced` would land in a run stack nobody reads.
#[test]
fn an_event_after_all_announced_fails_closed() {
    let codec = SpillCodec::Zstd;
    let (mut w, _dir) = make_writer(codec);
    w.process_event(announced(0, 0, 0)).unwrap();
    match w.process_event(keyed_run(codec, 0, 0..3)) {
        Err(err) => assert!(err.to_string().contains("after AllAnnounced"), "got: {err}"),
        Ok(()) => panic!("a block after AllAnnounced must fail closed"),
    }
}

#[test]
fn a_limit_below_two_never_consolidates() {
    let codec = SpillCodec::Zstd;
    let (w, _dir) = make_writer(codec);
    let registry = ConsolidationRegistry::new();
    let mut w = w.with_consolidation(1, registry);
    for file_id in 0..10 {
        w.process_event(keyed_run(codec, file_id, 0..3)).unwrap();
        assert!(w.active.is_none());
    }
    assert_eq!(w.runs.runs().len(), 10);
}

#[test]
fn block_for_wrong_file_id_while_open_errors() {
    let codec = SpillCodec::Zstd;
    let (mut w, _dir) = make_writer(codec);
    // Open file 0 with a non-last block, then feed a block for file 1 — a
    // contiguity violation that must fail loud, not silently corrupt file 0.
    w.process_event(block(codec, 0, false, &[1u8; 16])).unwrap();
    match w.process_event(block(codec, 1, false, &[2u8; 16])) {
        Err(err) => {
            assert!(
                err.to_string().contains("contiguous"),
                "expected contiguity error, got: {err}"
            );
        }
        Ok(()) => panic!("a block for a different open file_id must error"),
    }
}

#[test]
fn residual_while_file_open_errors() {
    let codec = SpillCodec::Zstd;
    let (mut w, _dir) = make_writer(codec);
    // Open a file with a non-last block, then feed a Residual — a missing
    // is_last_in_file terminator must fail loud, not drop the open spill.
    w.process_event(block(codec, 0, false, &[1u8; 16])).unwrap();
    let chunk = MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(vec![(
        RawCoordinateKey { sort_key: 1 },
        vec![9u8; 8],
    )]));
    match w.process_event(SpillBlockEvent::Residual {
        ordinal: 1,
        chunk,
        records_ingested_so_far: 1,
    }) {
        Err(err) => {
            assert!(err.to_string().contains("still open"), "expected open-file error, got: {err}");
        }
        Ok(()) => panic!("residual while a spill file is open must error"),
    }
}

#[test]
fn default_is_serial_writer_and_with_detached_flips_to_detached() {
    use fgumi_pipeline_core::step::{Affinity, Step, StepKind};
    let (w, _dir) = make_writer(SpillCodec::Zstd);
    // Default: pool-scheduled Serial + Affinity::Writer (pinned to worker N-1).
    assert_eq!(w.profile().kind, StepKind::Serial, "default spill writer is Serial");
    assert_eq!(w.affinity(), Affinity::Writer, "default spill writer pins to the writer worker");
    // `with_detached()` flips only the advertised kind to Detached (own thread,
    // off the pool); the write body — and hence the bytes it writes — is
    // unchanged, so full-sort parity still covers the on-disk format.
    let wd = w.with_detached();
    assert_eq!(wd.profile().kind, StepKind::Detached, "with_detached flips kind to Detached");
}

#[test]
fn residual_maps_to_memory_chunk_and_announced_passes_through() {
    let codec = SpillCodec::Zstd;
    let (mut w, _dir) = make_writer(codec);

    let chunk = MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(vec![(
        RawCoordinateKey { sort_key: 1 },
        vec![9u8; 8],
    )]));
    w.process_event(SpillBlockEvent::Residual { ordinal: 0, chunk, records_ingested_so_far: 3 })
        .unwrap();
    let Some(SortPhase1Event::MemoryChunk { chunk, records_ingested_so_far }) =
        w.outbox.pop_front()
    else {
        panic!("expected MemoryChunk");
    };
    assert_eq!(records_ingested_so_far, 3);
    assert_eq!(Arc::strong_count(&chunk), 1, "residual chunk wrapped in a fresh unique Arc");

    // No spill runs: `AllAnnounced` passes straight through with zero slots.
    w.process_event(announced(0, 1, 500)).unwrap();
    assert!(matches!(
        w.outbox.pop_front(),
        Some(SortPhase1Event::AllAnnounced {
            slot_count: 0,
            memory_chunk_count: 1,
            total_records: 500,
        })
    ));
}

// ── Step wiring: profile / affinity / detached group ─────────────────────────

#[test]
fn profile_defaults_to_a_pool_scheduled_serial_writer() {
    let (w, _dir) = make_writer(SpillCodec::Zstd);
    let profile = w.profile();
    assert_eq!(profile.name, "SpillWrite");
    assert_eq!(profile.kind, StepKind::Serial, "default is the pool-scheduled writer");
    assert!(profile.sticky);
    assert_eq!(profile.branch_ordering, vec![BranchOrdering::None]);
    match profile.output_queues.as_slice() {
        [QueueSpec::ByteBounded { limit_bytes }] => assert_eq!(*limit_bytes, 1 << 20),
        other => panic!("expected one byte-bounded queue, got {other:?}"),
    }
    // Affinity pins the pool-scheduled writer; ignored once detached.
    assert_eq!(w.affinity(), Affinity::Writer);
}

#[test]
fn with_detached_flips_only_the_step_kind() {
    let (w, _dir) = make_writer(SpillCodec::Zstd);
    let before = w.profile();
    let detached = w.with_detached();
    let after = detached.profile();

    assert_eq!(before.kind, StepKind::Serial);
    assert_eq!(after.kind, StepKind::Detached, "detached runs on its own thread");
    // Everything else about the step is unchanged — the doc promises the
    // `try_run` body and the bytes written are identical either way.
    assert_eq!(after.name, before.name);
    assert_eq!(after.sticky, before.sticky);
    assert_eq!(after.branch_ordering, before.branch_ordering);
    assert_eq!(detached.affinity(), Affinity::Writer);
}

#[test]
fn detached_writer_shares_the_sort_io_group() {
    // Phase-1 spill and phase-2 output writes are temporally disjoint, so both
    // ride the same driver thread rather than each taking one.
    let (w, _dir) = make_writer(SpillCodec::Zstd);
    assert_eq!(w.detached_group(), DetachedGroup::Shared(crate::sort::SORT_IO_GROUP));
    let (w2, _dir2) = make_writer(SpillCodec::Bgzf);
    assert_eq!(
        w2.with_detached().detached_group(),
        DetachedGroup::Shared(crate::sort::SORT_IO_GROUP)
    );
}

// ── Open-file bookkeeping ────────────────────────────────────────────────────

#[test]
fn ensure_no_open_file_passes_when_idle_and_fails_while_a_file_is_open() {
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    w.ensure_no_open_file("Residual").expect("idle writer has no open file");

    // Opening a file without its is_last block leaves it dangling.
    w.process_event(block(SpillCodec::Zstd, 3, false, &[7u8; 16])).unwrap();
    assert!(w.outbox.is_empty());
    assert!(w.current.is_some());

    let err = w.ensure_no_open_file("AllAnnounced").expect_err("dangling file must fail closed");
    let msg = err.to_string();
    assert!(msg.contains("AllAnnounced"), "error names the offending event: {msg}");
    assert!(msg.contains("file_id 3"), "error names the open file: {msg}");
}

#[test]
fn open_file_refuses_to_reuse_an_existing_path() {
    let (w, dir) = make_writer(SpillCodec::Zstd);
    // First open succeeds and creates the file on disk.
    let opened = w.open_file(9, SpillKeyKind::Coordinate).expect("first open succeeds");
    drop(opened);
    assert!(dir.path().join("chunk_0009.keyed").exists(), "spill file is created eagerly");

    // A reused file_id must fail closed rather than truncate the existing file:
    // silently overwriting a spill would drop records from the merge.
    // `OpenSpill` is not `Debug`, so match instead of using `expect_err`.
    match w.open_file(9, SpillKeyKind::Coordinate) {
        Ok(_) => panic!("reusing a file_id must fail"),
        Err(e) => assert_eq!(e.kind(), io::ErrorKind::AlreadyExists),
    }
}

/// `AllAnnounced` arriving while a spill file is still open must fail closed.
///
/// Without the guard, `AllAnnounced` reaches `SortMerge` before the matching
/// `SpillReady`, so the merge starts against an undercounted slot set and
/// silently drops a spill file's records.
#[test]
fn all_announced_while_a_file_is_open_fails_closed() {
    let codec = SpillCodec::Zstd;
    let (mut w, _dir) = make_writer(codec);
    // Open file 0 and never terminate it with an is_last_in_file block.
    w.process_event(block(codec, 0, false, &[1u8; 16])).unwrap();

    match w.process_event(SpillBlockEvent::AllAnnounced {
        ordinal: 1,
        slot_count: 1,
        memory_chunk_count: 0,
        total_records: 1,
        index_stride: 1,
    }) {
        Err(err) => {
            let msg = err.to_string();
            assert!(msg.contains("AllAnnounced"), "error names the event: {msg}");
            assert!(msg.contains("still open"), "error names the cause: {msg}");
            assert!(msg.contains("file_id 0"), "error names the open file: {msg}");
        }
        Ok(()) => panic!("AllAnnounced while a spill file is open must error"),
    }
}

// ── Domain counter (T-BW2): spill_bytes_written, driven through a real pipeline ──

/// `Exclusive` source draining a `Vec<SpillBlockEvent>`, one event per `try_run`.
struct SpillEventSource {
    events: Vec<SpillBlockEvent>,
    held: HeldSlot<Unpushed<SpillBlockEvent>>,
}

impl Step for SpillEventSource {
    type Input = ();
    type Outputs = Single<SpillBlockEvent>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SpillEventSource",
            kind: StepKind::Exclusive,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 20 }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take()
            && let Err(again) = ctx.outputs.retry(unpushed)
        {
            self.held.put(again);
            return Ok(StepOutcome::Progress);
        }
        let Some(event) = self.events.pop() else {
            return Ok(StepOutcome::Finished);
        };
        if let Err(unpushed) = ctx.outputs.push(event) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

/// Serial sink that sleeps briefly per received `SortPhase1Event` before
/// recording it — mirroring `fgumi_pipeline_core::tests::CountingSink`'s
/// per-item `thread::sleep`. `SpillWrite` (like `ReadBgzfBlocks` in
/// `source::read_bam::tests::try_run_bumps_blocks_and_bytes_read_counters`) is
/// the step under test and sits *upstream* of this sink, so it can burst
/// through every spill file it writes well within the first sampler tick;
/// throttling the terminal sink (rather than `SpillWrite` itself, which is
/// production code) keeps the run — and the sampler's window onto it — open
/// long enough to observe `SpillWrite`'s counter at its frozen final value.
struct ThrottledEventSink {
    received: Arc<Mutex<Vec<SortPhase1Event>>>,
}

impl Step for ThrottledEventSink {
    type Input = SortPhase1Event;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ThrottledEventSink",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(event) => {
                std::thread::sleep(std::time::Duration::from_micros(300));
                self.received.lock().push(event);
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// Drives `SpillEventSource -> SpillWrite -> ThrottledEventSink` with telemetry
/// enabled and asserts `SpillWrite`'s `spill_bytes_written` counter lands in the
/// telemetry files with a sane, bounded value: `N_FILES` single-block spill
/// files, each producing exactly one `SpillReady` (at `AllAnnounced`) for the
/// sink to (slowly) drain.
#[test]
fn try_run_bumps_spill_bytes_written_counter() {
    use fgumi_pipeline_core::builder::{InstrumentationLevel, Pipeline, PipelineConfig};
    use fgumi_pipeline_core::runtime::telemetry::TelemetryConfig;
    use std::time::Duration;

    const N_FILES: u32 = 30;
    const RAW_LEN: usize = 32;
    let codec = SpillCodec::Zstd;

    let dir = TempDir::new().unwrap();
    let alloc = TmpDirAllocator::new(vec![dir.path().to_path_buf()]).unwrap();
    let writer = SpillWrite::new(Arc::new(Mutex::new(alloc)), codec, 1 << 20, Arc::new(Vec::new()));

    let mut compressed_lens: Vec<u64> = Vec::with_capacity(N_FILES as usize);
    let mut events: Vec<SpillBlockEvent> = (0..N_FILES)
        .map(|file_id| {
            let fill = u8::try_from(file_id).expect("N_FILES fits in u8");
            let compressed = compress(codec, &[fill; RAW_LEN]);
            compressed_lens.push(compressed.len() as u64);
            SpillBlockEvent::Block {
                ordinal: u64::from(file_id),
                file_id,
                key_kind: SpillKeyKind::Coordinate,
                is_last_in_file: true,
                records_ingested_so_far: u64::from(file_id) + 1,
                bytes: compressed,
                frames: None,
                index_stride: 1,
            }
        })
        .collect();
    events.push(SpillBlockEvent::AllAnnounced {
        ordinal: u64::from(N_FILES),
        slot_count: N_FILES,
        memory_chunk_count: 0,
        total_records: u64::from(N_FILES),
        index_stride: 1,
    });
    events.reverse(); // `pop()` drains the tail first, so file_ids come out ascending.
    let expected_bytes: u64 = compressed_lens.iter().sum();

    let source = SpillEventSource { events, held: HeldSlot::new() };
    let received: Arc<Mutex<Vec<SortPhase1Event>>> = Arc::new(Mutex::new(Vec::new()));
    let sink = ThrottledEventSink { received: Arc::clone(&received) };

    let telemetry_dir =
        std::env::temp_dir().join(format!("fgumi-tbw2-spill-write-{}", std::process::id()));
    std::fs::create_dir_all(&telemetry_dir).unwrap();
    let stem = telemetry_dir.join("run");

    let builder = Pipeline::builder();
    builder.chain(source).chain(writer).chain(sink).into_sink_marker();
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

    // Ground truth, independent of the sampled telemetry file: one
    // `SpillReady` per file, then `AllAnnounced`, reached the sink.
    let collected = std::mem::take(&mut *received.lock());
    assert_eq!(collected.len(), N_FILES as usize + 1, "one SpillReady per spill file");

    // `SpillWrite` is step index 1 (source=0, writer=1, sink=2).
    let names = std::fs::read_to_string(telemetry_dir.join("run.ticks.counter_names.tsv")).unwrap();
    let name_rows: Vec<&str> = names.lines().skip(1).filter(|l| l.starts_with("1\t")).collect();
    assert_eq!(
        name_rows,
        vec!["1\t0\tspill_bytes_written\tbytes"],
        "SpillWrite declares its spill-bytes counter (consolidation records are ConsolidateRuns')"
    );

    let counters = std::fs::read_to_string(telemetry_dir.join("run.ticks.counters.tsv")).unwrap();
    let last_value = counters
        .lines()
        .skip(1)
        .filter(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f[2] == "1" && f[3] == "0"
        })
        .last()
        .map(|l| l.split('\t').nth(5).unwrap().parse::<u64>().unwrap())
        .expect("spill_bytes_written counter recorded at least once");
    assert!(last_value > 0 && last_value <= expected_bytes, "last_value={last_value}");

    std::fs::remove_dir_all(&telemetry_dir).ok();
}

/// `SpillBlockCompress` turns every range into a compressed block, so a range
/// reaching the writer is a wiring defect: it fails closed rather than write
/// unframed record indices as spill bytes.
#[test]
fn spill_write_rejects_an_unframed_range() {
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    let chunk = MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(vec![(
        RawCoordinateKey { sort_key: 1 },
        vec![7u8; 10],
    )]));
    let range = SpillBlockEvent::Range {
        ordinal: 0,
        file_id: 0,
        key_kind: SpillKeyKind::Coordinate,
        is_last_in_file: false,
        records_ingested_so_far: 1,
        chunk: Arc::new(chunk),
        start: 0,
        end: 1,
        framed_bytes: 14,
        frame_bytes: 0,
        index_stride: 1,
    };
    let err = w.process_event(range).expect_err("a range must not reach the writer");
    assert!(err.to_string().contains("unframed Range"), "got: {err}");
}

/// `n` mid-file blocks of file 0, each compressing `raw_len` bytes:
/// incompressible bytes when `noise`, so a block is about `raw_len` on disk.
fn mid_file_blocks(n: usize, raw_len: usize, noise: bool) -> VecDeque<SpillBlockEvent> {
    let codec = SpillCodec::Zstd;
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut byte = move |i: usize| {
        if !noise {
            return u8::try_from(i % 251).unwrap();
        }
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state.to_le_bytes()[0]
    };
    (0..n)
        .map(|i| {
            let raw: Vec<u8> = (0..raw_len).map(|_| byte(i)).collect();
            block(codec, 0, false, &raw)
        })
        .collect()
}

/// One `try_run` writes a batch: up to `WRITE_BLOCKS_PER_CALL` blocks, ending
/// early after a file's last block, once `WRITE_BYTES_PER_CALL` bytes are
/// written, or when the input runs dry. The bytes it reports are exactly those
/// of the blocks it consumed.
#[rstest]
#[case::count_cap(WRITE_BLOCKS_PER_CALL + 36, None, 32, false, WRITE_BLOCKS_PER_CALL)]
#[case::last_block_ends_the_batch(100, Some(10), 32, false, 11)]
#[case::input_runs_dry(5, None, 32, false, 5)]
#[case::byte_cap(12, None, 1 << 20, true, 8)]
fn a_write_batch_stops_at_its_caps_and_at_a_files_last_block(
    #[case] n: usize,
    #[case] last_at: Option<usize>,
    #[case] raw_len: usize,
    #[case] noise: bool,
    #[case] expected: usize,
) {
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    let mut input = mid_file_blocks(n, raw_len, noise);
    if let Some(i) = last_at
        && let SpillBlockEvent::Block { is_last_in_file, .. } = &mut input[i]
    {
        *is_last_in_file = true;
    }
    let total = input.len();
    let expected_bytes: u64 = input
        .iter()
        .take(expected)
        .map(|e| match e {
            SpillBlockEvent::Block { bytes, .. } => bytes.len() as u64,
            _ => 0,
        })
        .sum();
    let first = input.pop_front().unwrap();
    let written = w.write_batch(first, || input.pop_front()).unwrap();
    assert_eq!(total - input.len(), expected, "events consumed by one batch");
    assert_eq!(written, expected_bytes, "the batch reports the bytes it wrote");
}

/// The writer hands each written block's buffer back to the shared pool, and
/// the pool keeps at most its cap of cleared buffers.
#[test]
fn written_block_buffers_return_to_the_pool() {
    let codec = SpillCodec::Zstd;
    let pool = SpillBufPool::new();
    let (w, _dir) = make_writer(codec);
    let mut w = w.with_buf_pool(Arc::clone(&pool));
    let event = block(codec, 0, false, &[7u8; 300]);
    let SpillBlockEvent::Block { bytes, .. } = &event else { unreachable!() };
    let ptr = bytes.as_ptr();
    w.process_event(event).unwrap();
    let reused = pool.take();
    assert_eq!(reused.as_ptr(), ptr, "the written block's buffer is reused");
    assert!(reused.is_empty() && reused.capacity() > 0, "cleared, capacity kept");
}

fn job_writer(limit: usize) -> (SpillWrite, TempDir, Arc<ConsolidationRegistry>) {
    let (w, dir) = make_writer(SpillCodec::Zstd);
    let registry = ConsolidationRegistry::new();
    (w.with_consolidation(limit, Arc::clone(&registry)), dir, registry)
}

/// Review S2 (r1): `next_merge` fires at `n >= limit`, so with limit 2 the job
/// appears after the SECOND run. Catches: retiring before completion, losing the
/// merged run, or leaving absorbed inputs on disk.
#[test]
fn start_emits_a_job_and_retire_follows_completion() {
    let (mut w, dir, registry) = job_writer(2);
    w.process_event(keyed_run(SpillCodec::Zstd, 0, 0..10)).unwrap();
    assert!(w.outbox.is_empty(), "one run: no job");
    w.process_event(keyed_run(SpillCodec::Zstd, 1, 5..15)).unwrap();
    let Some(SortPhase1Event::ConsolidationStart { job_id, inputs, key_kind, index_stride }) =
        w.outbox.pop_front()
    else {
        panic!("expected a job")
    };
    assert_eq!(index_stride, 1, "the live runs are indexed densely");
    let output = registry.output_of(job_id).expect("the job's output is registered");
    assert_eq!(inputs.iter().map(|i| i.file_id).collect::<Vec<_>>(), [0, 1]);
    assert_eq!(key_kind, SpillKeyKind::TemplateK40);
    for i in &inputs {
        assert_eq!(i.len, std::fs::metadata(&i.path).unwrap().len(), "{:?}", i.path);
    }
    assert!(!output.exists(), "WriteConsolidatedRun creates it");

    assert!(!w.poll_registry().unwrap(), "no completion yet");
    assert!(w.outbox.is_empty());
    let bytes = write_run_file(&output, SpillCodec::Zstd, 0, 0..15);
    registry.complete(job_id, bytes, false);
    assert!(w.poll_registry().unwrap());
    assert!(matches!(
        w.outbox.pop_front(),
        Some(SortPhase1Event::ConsolidationRetire { job_id: j }) if j == job_id
    ));
    assert_eq!(
        w.runs.runs().iter().map(|r| (r.file_id, r.path.clone(), r.bytes)).collect::<Vec<_>>(),
        [(0, output.clone(), bytes)]
    );
    for i in &inputs {
        assert!(!i.path.exists(), "absorbed input {:?} removed", i.path);
    }
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 1);
    assert_eq!(w.stats.consolidations(), 1);
}

/// Backpressure: while a job is in flight, `try_run` takes no input. Catches: a
/// third run accepted mid-job (the phase switch could then precede a retire).
#[test]
fn no_input_is_taken_while_a_job_is_in_flight() {
    let (mut w, _dir, _registry) = job_writer(2);
    let probe = StepProbe::new(&w);
    probe.push_input(keyed_run(SpillCodec::Zstd, 0, 0..10));
    probe.push_input(keyed_run(SpillCodec::Zstd, 1, 0..10));
    probe.push_input(keyed_run(SpillCodec::Zstd, 2, 0..10));
    let mut outcomes = Vec::new();
    for _ in 0..8 {
        outcomes.push(probe.try_run(&mut w).unwrap());
    }
    assert!(w.active.is_some(), "the job is in flight: {outcomes:?}");
    assert!(!probe.input_is_empty(), "run 2 stays queued while the job runs");
    assert_eq!(*outcomes.last().unwrap(), StepOutcome::NoProgress);
}

#[test]
fn registry_complete_unparks_the_registered_waiter() {
    let registry = ConsolidationRegistry::new();
    let (tx, rx) = std::sync::mpsc::channel();
    let r2 = Arc::clone(&registry);
    let parked = std::thread::spawn(move || {
        let t0 = Instant::now();
        r2.register_waiter();
        let got = loop {
            if let Some(c) = r2.take_completed() {
                break c;
            }
            std::thread::park_timeout(Duration::from_secs(10));
        };
        tx.send((got, t0.elapsed())).unwrap();
    });
    std::thread::sleep(Duration::from_millis(50)); // let it park
    registry.complete(7, 4096, false);
    let (got, waited) =
        rx.recv_timeout(Duration::from_secs(5)).expect("the waiter must not sit out its 10 s park");
    assert_eq!(got, (7, 4096, false));
    assert!(waited < Duration::from_secs(1), "{waited:?}");
    parked.join().unwrap();
}

/// The lost-wake window, deterministically: the completion lands before the
/// waiter registers, and `poll_registry`'s re-check after `register_waiter`
/// returns it. Catches: register-then-park without the re-check.
#[test]
fn register_then_recheck_catches_a_completion_in_between() {
    let (mut w, _dir, registry) = job_writer(2);
    w.process_event(keyed_run(SpillCodec::Zstd, 0, 0..10)).unwrap();
    w.process_event(keyed_run(SpillCodec::Zstd, 1, 0..10)).unwrap();
    let Some(SortPhase1Event::ConsolidationStart { job_id, .. }) = w.outbox.pop_front() else {
        panic!()
    };
    let output = registry.output_of(job_id).expect("registered");
    let bytes = write_run_file(&output, SpillCodec::Zstd, 0, 0..20);
    registry.set_hook_before_register_for_test(Box::new(move |r: &ConsolidationRegistry| {
        r.complete(job_id, bytes, false);
    }));
    assert!(w.poll_registry().unwrap(), "the re-check after register_waiter finds the completion");
}

// ── Run index ────────────────────────────────────────────────────────────────

/// Records per frame in the run-index tests: every record frames to the same
/// size and the frame size is this many records, so each frame closes on it.
const RECORDS_PER_FRAME: usize = 3;

/// Body length of every record in the run-index tests.
const INDEX_BODY_LEN: usize = 60;

/// `key`'s serialized bytes.
fn key_bytes<K: fgumi_sort::RawSortKey>(key: &K) -> Vec<u8> {
    let mut out = Vec::new();
    key.write_to(&mut out).unwrap();
    out
}

/// A sorted chunk of `n` records of `kind`, keys ascending from `base`, plus
/// each record's serialized key (what an index stores for a frame it starts).
fn index_chunk(kind: SpillKeyKind, base: u64, n: u64) -> (MemoryChunkErased, Vec<Vec<u8>>) {
    use crate::sort::spill_fixtures::{coordinate_body_len, template_key};
    use fgumi_sort::{RawSortKey as _, TemplateKey24, TemplateMemChunk};
    let body = |i: u64| vec![u8::try_from(i % 251).unwrap(); INDEX_BODY_LEN];
    match kind {
        SpillKeyKind::Coordinate => {
            let recs: Vec<_> = (0..n)
                .map(|i| {
                    let pos = i32::try_from(base + i).unwrap();
                    let b = coordinate_body_len(0, pos, u32::try_from(i).unwrap(), INDEX_BODY_LEN);
                    (RawCoordinateKey::extract_from_record(&b), b)
                })
                .collect();
            let keys = recs.iter().map(|(k, _)| key_bytes(k)).collect();
            (MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(recs)), keys)
        }
        SpillKeyKind::TemplateK24 => {
            let recs: Vec<_> = (0..n)
                .map(|i| {
                    (TemplateKey24 { primary: base + i, secondary: 0, name_hash_upper: 0 }, body(i))
                })
                .collect();
            let keys = recs.iter().map(|(k, _)| key_bytes(k)).collect();
            let chunk = TemplateMemChunk::K24(InMemoryChunk::from_owned_records(recs));
            (MemoryChunkErased::TemplateCoordinate(chunk), keys)
        }
        SpillKeyKind::TemplateK40 => {
            let recs: Vec<_> = (0..n).map(|i| (template_key(base + i), body(i))).collect();
            let keys = recs.iter().map(|(k, _)| key_bytes(k)).collect();
            let chunk = TemplateMemChunk::K40(InMemoryChunk::from_owned_records(recs));
            (MemoryChunkErased::TemplateCoordinate(chunk), keys)
        }
        other => unreachable!("no index-test chunk for {other:?}"),
    }
}

/// Spill `chunks` (one run each when their keys do not extend the previous
/// run) through the production path — `SpillGather` at `block_size`,
/// `SpillBlockCompress` with the run index `enabled`, `SpillWrite` — and
/// return each announced run's slot and path, in announcement order.
fn spill_through_gather(
    codec: SpillCodec,
    enabled: bool,
    chunks: Vec<MemoryChunkErased>,
    block_size: usize,
) -> (TempDir, Vec<(Arc<fgumi_sort::SortMergeSlot>, PathBuf)>) {
    use crate::sort::protocol::SortChunkEvent;
    use crate::sort::spill_block_compress::SpillBlockCompress;
    use crate::sort::spill_gather::SpillGather;
    let (mut w, dir) = make_writer(codec);
    let mut gather = SpillGather::new(1 << 20).with_block_size_for_test(block_size);
    let mut compress = SpillBlockCompress::new(codec, 1, 1 << 20).with_run_index(enabled);
    let mut events: Vec<SortChunkEvent> = chunks
        .into_iter()
        .enumerate()
        .map(|(seq, chunk)| SortChunkEvent::Spill {
            seq: u32::try_from(seq).unwrap(),
            chunk,
            records_ingested_so_far: 0,
            index_stride: 1,
        })
        .collect();
    events.push(SortChunkEvent::AllAnnounced {
        slot_count: 0,
        memory_chunk_count: 0,
        total_records: 0,
        index_stride: 1,
    });
    for event in gather.drive_for_test(events) {
        w.process_event(compress.compress_event(event).unwrap()).unwrap();
    }
    let slots = w
        .outbox
        .drain(..)
        .filter_map(|e| match e {
            SortPhase1Event::SpillReady { slot, path, .. } => Some((slot, path)),
            _ => None,
        })
        .collect();
    (dir, slots)
}

/// Every frame a spill run's index records is the frame at that position in
/// the file: frame `i` starts at the absolute file offset `index.offset(i)`
/// (`comp_len` counting the `u32` length prefix), the frames tile the body up to the
/// trailer, and each frame's raw length, record count and first key are what
/// decompressing it shows. Two runs per case, so the index restarts per file;
/// 1 frame is a run of only the gather's raw final block, 7 and 40 add
/// multi-frame ranges.
#[rstest]
#[case::coordinate(SpillKeyKind::Coordinate)]
#[case::k24(SpillKeyKind::TemplateK24)]
#[case::k40(SpillKeyKind::TemplateK40)]
fn spilled_run_index_matches_its_file(
    #[case] kind: SpillKeyKind,
    #[values(1, 7, 40)] frames: usize,
) {
    let codec = SpillCodec::Zstd;
    let n = (frames * RECORDS_PER_FRAME) as u64;
    let framed_len = kind.serialized_key_bytes() + 4 + INDEX_BODY_LEN;
    // Run 1's keys precede run 0's, so it starts a new file.
    let (run0, keys0) = index_chunk(kind, 1_000, n);
    let (run1, keys1) = index_chunk(kind, 0, n);
    let (_dir, slots) =
        spill_through_gather(codec, true, vec![run0, run1], RECORDS_PER_FRAME * framed_len);
    assert_eq!(slots.len(), 2, "two runs announced");
    let mut dec = fgumi_sort::SpillBlockDecompressor::new();
    for ((slot, path), keys) in slots.iter().zip([keys0, keys1]) {
        let index = slot.run_index().expect("an indexed zstd run carries its index");
        assert_eq!(index.len(), frames, "{path:?}");
        assert_eq!(index.key_width(), fgumi_sort::spill_key_width(kind).unwrap());
        let bytes = std::fs::read(path).unwrap();
        // The footer follows the frames; the slot's frame bytes exclude it.
        assert_eq!(bytes.len() as u64, slot.len() + index.footer_bytes(), "{path:?}");
        let mut at = usize::try_from(slot.body_start()).unwrap();
        let mut offset = slot.body_start();
        for i in 0..frames {
            assert_eq!(at as u64, index.offset(i), "frame {i} starts at its absolute offset");
            assert_eq!(index.offset(i), offset, "frame {i} follows the previous frame");
            let len = u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap()) as usize;
            assert_eq!(index.comp_len(i), 4 + len as u64, "frame {i} comp_len");
            let raw = dec.decompress_one(codec, &bytes[at + 4..at + 4 + len]).unwrap();
            assert_eq!(raw.len(), index.raw_len(i) as usize, "frame {i} raw_len");
            let (mut pos, mut records) = (0usize, 0u32);
            while pos < raw.len() {
                let len_at = pos + kind.serialized_key_bytes();
                let body = u32::from_le_bytes(raw[len_at..len_at + 4].try_into().unwrap());
                pos = len_at + 4 + body as usize;
                records += 1;
            }
            assert_eq!(records as usize, RECORDS_PER_FRAME, "frame {i} holds a frame's records");
            assert_eq!(index.records(i), records, "frame {i} record count");
            assert_eq!(index.key_bytes(i), keys[i * RECORDS_PER_FRAME], "frame {i} first key");
            at += 4 + len;
            offset += index.comp_len(i);
        }
        assert_eq!(at as u64, slot.len(), "the frames tile the body");
        assert_eq!(
            index.end(index.len() - 1),
            slot.len(),
            "the index ends where the footer starts"
        );
        assert_eq!(index.total_records(), n);
    }
}

/// A run-final block of five `[u32 len][body]` coordinate records, keyed
/// `kind`, compressed for `codec` with the run index `enabled`.
fn raw_coordinate_block(kind: SpillKeyKind, codec: SpillCodec, enabled: bool) -> SpillBlockEvent {
    use crate::sort::spill_block_compress::SpillBlockCompress;
    use crate::sort::spill_fixtures::{coordinate_body, embedded_frame};
    let raw: Vec<u8> = (0..5)
        .flat_map(|i| embedded_frame(&coordinate_body(0, i, u32::try_from(i).unwrap())))
        .collect();
    SpillBlockCompress::new(codec, 1, 1 << 20)
        .with_run_index(enabled)
        .compress_event(SpillBlockEvent::Block {
            ordinal: 0,
            file_id: 0,
            key_kind: kind,
            is_last_in_file: true,
            records_ingested_so_far: 5,
            bytes: raw,
            frames: None,
            index_stride: 1,
        })
        .unwrap()
}

/// No index when it is disabled, the codec is BGZF, or the order is queryname
/// (whose keys an index does not hold): the block carries none and the
/// announced slot has none.
#[rstest]
#[case::disabled(false, SpillCodec::Zstd, SpillKeyKind::Coordinate)]
#[case::bgzf(true, SpillCodec::Bgzf, SpillKeyKind::Coordinate)]
#[case::queryname(true, SpillCodec::Zstd, SpillKeyKind::QuerynameLex)]
fn index_absent_when(#[case] enabled: bool, #[case] codec: SpillCodec, #[case] kind: SpillKeyKind) {
    let event = raw_coordinate_block(kind, codec, enabled);
    assert!(
        matches!(&event, SpillBlockEvent::Block { frames: None, .. }),
        "no frames on the block"
    );
    let (mut w, _dir) = make_writer(codec);
    w.process_event(event).unwrap();
    w.process_event(announced(1, 0, 5)).unwrap();
    let Some(SortPhase1Event::SpillReady { slot, .. }) = w.outbox.pop_front() else {
        panic!("expected SpillReady");
    };
    assert!(slot.run_index().is_none());
}

/// A file with any block that carries no frames has no index, whichever block
/// it is: an index must cover every frame of its file or none.
#[rstest]
#[case::first_block_unindexed(false, true)]
#[case::second_block_unindexed(true, false)]
fn a_file_with_an_unindexed_block_has_no_index(#[case] first: bool, #[case] second: bool) {
    let codec = SpillCodec::Zstd;
    let mut a = raw_coordinate_block(SpillKeyKind::Coordinate, codec, first);
    let mut b = raw_coordinate_block(SpillKeyKind::Coordinate, codec, second);
    if let SpillBlockEvent::Block { is_last_in_file, .. } = &mut a {
        *is_last_in_file = false;
    }
    if let SpillBlockEvent::Block { ordinal, .. } = &mut b {
        *ordinal = 1;
    }
    let (mut w, _dir) = make_writer(codec);
    w.process_event(a).unwrap();
    w.process_event(b).unwrap();
    w.process_event(announced(1, 0, 10)).unwrap();
    let Some(SortPhase1Event::SpillReady { slot, .. }) = w.outbox.pop_front() else {
        panic!("expected SpillReady");
    };
    assert!(slot.run_index().is_none());
}

/// A consolidated run is announced with an index covering its whole file, and
/// the runs consolidation left alone keep theirs. Driven through `KChain`, the
/// consolidation chain pumped by hand: limit 3 merges the two small newest runs
/// and leaves the large first one. With a BGZF temp codec the merged run has no
/// frames to index.
#[rstest]
#[case::zstd(SpillCodec::Zstd, true)]
#[case::bgzf(SpillCodec::Bgzf, false)]
fn consolidated_run_is_indexed(#[case] codec: SpillCodec, #[case] expect_index: bool) {
    use crate::sort::protocol::SortPhase2Event;
    use crate::sort::spill_block_compress::SpillBlockCompress;
    use crate::sort::spill_fixtures::{KChain, KOpts, framed_template_records_of};
    let mut k = KChain::new(KOpts { limit: 3, codec, ..KOpts::default() });
    let mut compress = SpillBlockCompress::new(codec, 1, 1 << 20).with_run_index(true);
    for (id, keys) in [(0u32, 0..400u64), (1, 0..5), (2, 3..8)] {
        let event = compress
            .compress_event(SpillBlockEvent::Block {
                ordinal: u64::from(id),
                file_id: id,
                key_kind: SpillKeyKind::TemplateK40,
                is_last_in_file: true,
                records_ingested_so_far: u64::from(id) + 1,
                bytes: framed_template_records_of(id, keys),
                frames: None,
                index_stride: 1,
            })
            .unwrap();
        assert_eq!(
            matches!(&event, SpillBlockEvent::Block { frames: Some(_), .. }),
            codec == SpillCodec::Zstd
        );
        k.feed(event);
        k.settle().unwrap();
    }
    assert_eq!(k.jobs.len(), 1, "one consolidation");
    assert_eq!(k.jobs[0].1.len(), 2, "of two runs");
    k.feed(announced(3, 0, 0));
    k.close();
    k.settle().unwrap();
    let slots: Vec<_> = k
        .phase2
        .iter()
        .filter_map(|e| match e {
            SortPhase2Event::SpillReady { slot, .. } => Some(slot.clone()),
            _ => None,
        })
        .collect();
    let indexed: Vec<(u32, bool)> =
        slots.iter().map(|s| (s.file_id, s.run_index().is_some())).collect();
    assert_eq!(indexed, [(0, codec == SpillCodec::Zstd), (1, expect_index)]);
    if expect_index {
        let merged = &slots[1];
        let index = merged.run_index().expect("indexed");
        assert_eq!(index.offset(0), merged.body_start(), "frames start after the magic");
        assert_eq!(index.end(index.len() - 1), merged.len(), "frames end at the footer");
        assert!(merged.has_index_footer(), "the merged run's index is its footer");
        for i in 1..index.len() {
            assert_eq!(index.offset(i), index.end(i - 1), "frame {i} abuts the previous");
        }
        assert_eq!(index.total_records(), 10, "the two merged runs' records");
    }
}

// ── Phase-1 index policy (followed, not decided) ─────────────────────────────

const MIB: u64 = 1 << 20;

/// An eligibility with indexed runs, not forced serial.
fn eligibility(phase2: usize, merge_memory: u64, serial_r: u64) -> RangeEligibility {
    RangeEligibility {
        indexable: true,
        phase2_threads: phase2,
        pool_threads: phase2,
        force_serial: false,
        merge_memory,
        serial_r,
    }
}

/// A one-block run `file_id` keyed `kind`, carrying a stand-in dense index of
/// `frames` frames (`key_width`-byte keys, 1000 raw bytes and one record
/// each) and stamped `index_stride` by the front: the writer stores the bytes
/// and the index without decoding either.
fn indexed_run_at(
    file_id: u32,
    kind: SpillKeyKind,
    frames: u32,
    index_stride: u32,
) -> SpillBlockEvent {
    let width = fgumi_sort::spill_key_width(kind).expect("an indexed kind");
    let mut index = RunIndex::new(width);
    for f in 0..frames {
        index.push_frame(10, 1000, 1, &vec![u8::try_from(f % 256).unwrap(); usize::from(width)]);
    }
    SpillBlockEvent::Block {
        ordinal: 0,
        file_id,
        key_kind: kind,
        is_last_in_file: true,
        records_ingested_so_far: u64::from(frames),
        bytes: vec![0u8; 10 * frames as usize],
        frames: Some(index),
        index_stride,
    }
}

/// [`indexed_run_at`] at stride 1 (dense).
fn indexed_run(file_id: u32, kind: SpillKeyKind, frames: u32) -> SpillBlockEvent {
    indexed_run_at(file_id, kind, frames, 1)
}

/// The live runs' index footers.
fn live_footers(w: &SpillWrite) -> Vec<bool> {
    w.runs.runs().iter().map(|r| r.index_footer).collect()
}

/// Announce `w`'s `n` closed runs (no residual) and return the announced
/// slots' index strides (`None`: no index attached).
fn announce_strides(w: &mut SpillWrite, n: u32) -> Vec<Option<u32>> {
    w.process_event(announced(n, 0, 0)).unwrap();
    w.outbox
        .drain(..)
        .filter_map(|e| match e {
            SortPhase1Event::SpillReady { slot, .. } => {
                Some(slot.run_index().map(|ix| ix.stride()))
            }
            _ => None,
        })
        .collect()
}

/// A closed run's index is written to its file's footer and dropped: the
/// writer holds no index between runs, the run's frame bytes end where the
/// footer starts, and the footer reads back as the index the blocks built
/// (offsets past the codec magic). Columns: frames in the run.
#[rstest]
#[case::one_frame(1)]
#[case::many_frames(1_000)]
fn close_run_writes_the_footer_and_the_slot_excludes_it(#[case] frames: u32) {
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    let SpillBlockEvent::Block { frames: Some(mut expected), .. } =
        indexed_run(0, SpillKeyKind::TemplateK24, frames)
    else {
        unreachable!()
    };
    let magic = spill_magic(SpillCodec::Zstd).len() as u64;
    expected.rebase(magic);
    w.process_event(indexed_run(0, SpillKeyKind::TemplateK24, frames)).unwrap();
    assert_eq!(live_footers(&w), [true]);
    let path = w.runs.runs()[0].path.clone();
    let body_end = magic + 10 * u64::from(frames);
    assert_eq!(w.runs.runs()[0].bytes, body_end);
    assert_eq!(std::fs::metadata(&path).unwrap().len(), body_end + expected.footer_bytes());
    let slot = fgumi_sort::open_spill_slot(&path, 0, fgumi_sort::FramesEnd::At(body_end)).unwrap();
    assert_eq!(slot.len(), body_end);
    assert_eq!(slot.read_run_index().unwrap(), Some(expected.clone()));
    w.process_event(announced(1, 0, 0)).unwrap();
    let Some(SortPhase1Event::SpillReady { slot, .. }) = w.outbox.pop_front() else {
        panic!("expected SpillReady");
    };
    assert_eq!(slot.len(), body_end);
    assert_eq!(slot.run_index().map(|ix| (**ix).clone()), Some(expected));
}

/// A run stamped stride 0 (the front found no range fit) drops the open
/// run's index and stops indexing for good: that run and later ones get no
/// footer, and `AllAnnounced` reads no footer back, so no slot is indexed.
#[test]
fn indexing_off_reads_no_footer() {
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    for id in 0..3 {
        w.process_event(indexed_run(id, SpillKeyKind::Coordinate, 4)).unwrap();
    }
    assert_eq!(live_footers(&w), vec![true; 3]);
    w.process_event(indexed_run_at(3, SpillKeyKind::Coordinate, 4, 0)).unwrap();
    assert_eq!(w.index_policy_for_test(), (false, 1));
    w.process_event(indexed_run(4, SpillKeyKind::Coordinate, 4)).unwrap();
    assert_eq!(live_footers(&w), [true, true, true, false, false], "off for good");
    assert_eq!(announce_strides(&mut w, 5), vec![None; 5]);
}

/// Runs closed at stride 1 and a run closed after the front raised it to 2
/// are all announced at stride 2: `AllAnnounced` sparsifies every footer to
/// the final stride, keeping each run's frames and raw bytes.
#[test]
fn announce_sparsifies_footers_to_the_final_stride() {
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    w.process_event(indexed_run(0, SpillKeyKind::Coordinate, 4)).unwrap();
    w.process_event(indexed_run(1, SpillKeyKind::Coordinate, 4)).unwrap();
    w.process_event(indexed_run_at(2, SpillKeyKind::Coordinate, 4, 2)).unwrap();
    assert_eq!(w.index_policy_for_test(), (true, 2));
    w.process_event(announced(3, 0, 0)).unwrap();
    let slots: Vec<_> = w
        .outbox
        .drain(..)
        .filter_map(|e| match e {
            SortPhase1Event::SpillReady { slot, .. } => Some(slot),
            _ => None,
        })
        .collect();
    assert_eq!(slots.len(), 3);
    for slot in &slots {
        let index = slot.run_index().expect("indexed");
        assert_eq!(
            (index.stride(), index.len(), index.covered_frames(), index.total_raw_bytes()),
            (2, 2, 4, 4_000),
            "four frames in two entries"
        );
    }
}

/// A stride lower than the live one changes nothing: an index cannot be made
/// denser again.
#[test]
fn spill_write_never_lowers_the_stride() {
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    w.process_event(indexed_run_at(0, SpillKeyKind::Coordinate, 4, 2)).unwrap();
    w.process_event(indexed_run_at(1, SpillKeyKind::Coordinate, 4, 1)).unwrap();
    assert_eq!(w.index_policy_for_test(), (true, 2));
    assert_eq!(announce_strides(&mut w, 2), vec![Some(2); 2]);
}

/// How a test damages a closed run's file before `AllAnnounced`.
#[derive(Clone, Copy, Debug)]
enum Damage {
    /// Drop the file's last byte: the trailer is gone.
    TruncateTrailer,
    /// Overwrite the index block's magic.
    BlockMagic,
    /// Make the first offset point past the frames' start.
    FirstOffset,
}

/// A run written with a footer whose file no longer reads as written fails
/// the announcement: it is never read as frames, nor silently left
/// unindexed. Columns: the damage, a needle of the expected error.
#[rstest]
#[case::truncated_trailer(Damage::TruncateTrailer, "do not end in a footer trailer")]
#[case::bad_block_magic(Damage::BlockMagic, "block magic")]
#[case::offsets_do_not_tile(Damage::FirstOffset, "the run's frames")]
fn a_damaged_footer_fails_the_announcement(#[case] damage: Damage, #[case] needle: &str) {
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    w.process_event(indexed_run(0, SpillKeyKind::Coordinate, 4)).unwrap();
    let run = w.runs.runs()[0].clone();
    let mut bytes = std::fs::read(&run.path).unwrap();
    let block = usize::try_from(run.bytes).unwrap();
    match damage {
        Damage::TruncateTrailer => {
            bytes.pop();
        }
        Damage::BlockMagic => bytes[block] ^= 0xFF,
        // The block's fixed header is 36 bytes; the first offset follows.
        Damage::FirstOffset => bytes[block + 36] ^= 0x01,
    }
    std::fs::write(&run.path, &bytes).unwrap();
    let err = w.process_event(announced(1, 0, 0)).unwrap_err();
    assert!(err.to_string().contains(needle), "{err}");
}

/// A block of run `file_id` keyed `kind` with `frames` stand-in frames of
/// the given raw lengths (10 bytes and one record each), stamped `stride`,
/// last in its run or not.
fn run_block(
    file_id: u32,
    kind: SpillKeyKind,
    raw: &[u32],
    is_last_in_file: bool,
    stride: u32,
) -> SpillBlockEvent {
    let width = fgumi_sort::spill_key_width(kind).expect("an indexed kind");
    let mut index = RunIndex::new(width);
    for &r in raw {
        index.push_frame(10, r, 1, &vec![0u8; usize::from(width)]);
    }
    SpillBlockEvent::Block {
        ordinal: 0,
        file_id,
        key_kind: kind,
        is_last_in_file,
        records_ingested_so_far: raw.len() as u64,
        bytes: vec![0u8; 10 * raw.len()],
        frames: Some(index),
        index_stride: stride,
    }
}

/// The open run's index is held at the live stride, block by block, never
/// dense: after a first block of 9 frames it has `9 / stride` entries
/// (rounded up), and the run's footer covers both blocks at that stride.
/// Columns: the stamped stride, entries after the first block, entries in
/// the footer.
#[rstest]
#[case::dense(1, 9, 12)]
#[case::stride_four(4, 3, 3)]
fn the_open_run_index_is_held_at_the_live_stride(
    #[case] stride: u32,
    #[case] open_entries: usize,
    #[case] footer_entries: usize,
) {
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    w.process_event(run_block(0, SpillKeyKind::Coordinate, &[1000; 9], false, stride)).unwrap();
    assert_eq!(w.open_index_shape_for_test(), Some((stride, open_entries)));
    w.process_event(run_block(0, SpillKeyKind::Coordinate, &[1000; 3], true, stride)).unwrap();
    assert_eq!(w.open_index_shape_for_test(), None, "the closed run's index is in its footer");
    let run = w.runs.runs()[0].clone();
    let slot =
        fgumi_sort::open_spill_slot(&run.path, 0, fgumi_sort::FramesEnd::At(run.bytes)).unwrap();
    let index = slot.read_run_index().unwrap().expect("a footer");
    assert_eq!((index.stride(), index.len(), index.covered_frames()), (stride, footer_entries, 12));
}

/// Phase 1 holds no closed run's index: 64 runs whose blocks carry
/// 16,384-frame indexes of 40-byte entries (640 KiB each, 40 MiB in all)
/// each go to their footer, and process RSS grows by well under the 40 MiB
/// the indexes would hold if they were kept. The RSS check runs only under
/// nextest (one process per test): in `cargo test`'s shared process,
/// concurrent tests move RSS by more than the margin.
#[test]
fn phase1_holds_no_index_memory() {
    const RUNS: u32 = 64;
    const FRAMES: u32 = 16_384;
    assert_eq!(RunIndex::entry_bytes(24), 40);
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    let mut rss_early = None;
    for id in 0..RUNS {
        w.process_event(indexed_run(id, SpillKeyKind::TemplateK24, FRAMES)).unwrap();
        if id == 3 {
            rss_early = fgumi_sort::process_rss_bytes();
        }
    }
    let rss_late = fgumi_sort::process_rss_bytes();
    assert_eq!(live_footers(&w), vec![true; RUNS as usize]);
    if std::env::var_os("NEXTEST").is_some()
        && let (Some(before), Some(after)) = (rss_early, rss_late)
    {
        let growth = after.saturating_sub(before);
        assert!(growth < 8 << 20, "RSS grew {growth} bytes over 60 runs");
    }
}

/// Spill `runs` through a writer and announce them as `slot_count` runs with
/// no residual chunk; return the announced slots.
fn announce(runs: Vec<SpillBlockEvent>) -> (TempDir, Vec<Arc<fgumi_sort::SortMergeSlot>>) {
    let (mut w, dir) = make_writer(SpillCodec::Zstd);
    let n = u32::try_from(runs.len()).unwrap();
    for run in runs {
        w.process_event(run).unwrap();
    }
    w.process_event(announced(n, 0, 0)).unwrap();
    let slots = w
        .outbox
        .drain(..)
        .filter_map(|e| match e {
            SortPhase1Event::SpillReady { slot, .. } => Some(slot),
            _ => None,
        })
        .collect();
    (dir, slots)
}

/// The merge decision fits at the indexes' own stride: after the writer
/// followed a stride-2 stamp, the published fit's stride is every slot
/// index's stride, and the published index bytes are the indexes' exact heap
/// bytes.
#[test]
fn the_merge_decision_fits_at_the_written_stride() {
    let elig = eligibility(8, 512 * MIB, 64 * MIB);
    let runs = (0..4).map(|id| indexed_run_at(id, SpillKeyKind::TemplateK40, 9, 2)).collect();
    let (_dir, slots) = announce(runs);
    let mode = crate::sort::decide_merge_mode(Some(&elig), &slots, 0);
    let crate::sort::MergeMode::Range { fit, index_bytes, .. } = mode else {
        panic!("range mode: {mode:?}");
    };
    let indexes: Vec<_> = slots.iter().map(|s| s.run_index().expect("indexed")).collect();
    assert!(indexes.iter().all(|ix| ix.stride() == fit.stride), "{fit:?}");
    assert_eq!(fit.stride, 2);
    let heap: u64 = indexes.iter().map(|ix| ix.heap_bytes() as u64).sum();
    assert_eq!(index_bytes, heap);
    // Shrunk to fit: five entries of a 40-byte key, six offsets, per run.
    assert_eq!(heap, 4 * (5 * (40 + 8) + 6 * 8));
}

/// A final run the front spilled keeps every index: the writer decides
/// nothing at the last close, so the merge runs by range over every run. At 4
/// workers a 10.6 MiB budget fits two of these runs but not three, so a
/// writer fitting `live + 1` runs at the last close (the old policy) dropped
/// both indexes after the front had already spilled the final run, and the
/// sort ended serial with no residual.
#[test]
fn final_spill_run_keeps_its_indexes() {
    let elig = eligibility(4, (10 << 20) + (600 << 10), MIB);
    let runs = (0..2).map(|id| indexed_run(id, SpillKeyKind::Coordinate, 1)).collect();
    let (_dir, slots) = announce(runs);
    assert_eq!(
        slots.iter().map(|s| s.run_index().map(|ix| ix.stride())).collect::<Vec<_>>(),
        vec![Some(1); 2]
    );
    let mode = crate::sort::decide_merge_mode(Some(&elig), &slots, 0);
    assert!(
        matches!(mode, crate::sort::MergeMode::Range { fit, .. } if fit.stride == 1),
        "{mode:?}"
    );
}

/// A consolidated run's footer is read back with the live runs' at
/// `AllAnnounced`, at their stride, and none is read once indexing is off.
/// Columns: the stride the front stamped on the two input runs, whether the
/// merged run carries a footer (its job indexes only while indexing is on),
/// the merged run's expected stride.
#[rstest]
#[case::sparsified(2, true, Some(2))]
#[case::dropped(0, false, None)]
fn a_consolidated_run_follows_the_policy(
    #[case] stride: u32,
    #[case] footer: bool,
    #[case] expected: Option<u32>,
) {
    let (mut w, _dir, registry) = job_writer(2);
    w.process_event(indexed_run_at(0, SpillKeyKind::TemplateK40, 9, stride)).unwrap();
    w.process_event(indexed_run_at(1, SpillKeyKind::TemplateK40, 9, stride)).unwrap();
    let Some(SortPhase1Event::ConsolidationStart { job_id, .. }) = w.outbox.pop_front() else {
        panic!("expected a job")
    };
    let output = registry.output_of(job_id).expect("registered");
    // The merged run: 18 stand-in frames of 10 bytes and, when indexed, its
    // dense index as the footer.
    let magic = spill_magic(SpillCodec::Zstd);
    let mut file = magic.to_vec();
    file.extend_from_slice(&[0u8; 180]);
    let body_end = file.len() as u64;
    if footer {
        let mut dense = RunIndex::new(40);
        for _ in 0..18 {
            dense.push_frame(10, 1000, 1, &[0u8; 40]);
        }
        dense.rebase(magic.len() as u64);
        dense.write_footer(&mut file).unwrap();
    }
    std::fs::write(&output, &file).unwrap();
    registry.complete(job_id, body_end, footer);
    assert!(w.poll_registry().unwrap());
    assert_eq!(live_footers(&w), [footer]);
    assert_eq!(w.index_policy_for_test(), (expected.is_some(), stride.max(1)));
    assert_eq!(announce_strides(&mut w, 2), [expected]);
}

/// `AllAnnounced` carries the front's last stride, which no block carries
/// for a resident final run. With the final run spilled, the footers are
/// read back at that stride (a raised one sparsifies them). With a resident
/// final run (a residual chunk) the merge is serial whatever the indexes
/// say, so no footer is read back, whether the front turned indexing off
/// (stride 0) or not, and the merge reports the residual. Columns: the
/// announced stride, residual chunks, expected slot strides, whether the
/// runs keep their footers, the expected serial reason (`None`: range mode).
#[rstest]
#[case::raised_with_the_final_run_spilled(2, 0, vec![Some(2), Some(2)], true, None)]
#[case::unchanged_with_the_final_run_spilled(1, 0, vec![Some(1), Some(1)], true, None)]
#[case::resident_final_run_reads_no_footer(
    1,
    1,
    vec![None, None],
    true,
    Some(crate::sort::SerialReason::ResidualInMemory { chunks: 1 })
)]
#[case::fit_failed_at_the_final_run(
    0,
    1,
    vec![None, None],
    true,
    Some(crate::sort::SerialReason::ResidualInMemory { chunks: 1 })
)]
fn the_announced_stride_is_followed(
    #[case] stride: u32,
    #[case] residual_chunks: u32,
    #[case] expected: Vec<Option<u32>>,
    #[case] footers: bool,
    #[case] reason: Option<crate::sort::SerialReason>,
) {
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    for id in 0..2 {
        w.process_event(indexed_run(id, SpillKeyKind::Coordinate, 4)).unwrap();
    }
    for _ in 0..residual_chunks {
        let residual = SpillBlockEvent::Residual {
            ordinal: 0,
            chunk: MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(vec![(
                RawCoordinateKey { sort_key: 9 },
                vec![1u8; 40],
            )])),
            records_ingested_so_far: 9,
        };
        w.process_event(residual).unwrap();
    }
    w.process_event(SpillBlockEvent::AllAnnounced {
        ordinal: 0,
        slot_count: 2,
        memory_chunk_count: residual_chunks,
        total_records: 9,
        index_stride: stride,
    })
    .unwrap();
    let slots: Vec<_> = w
        .outbox
        .drain(..)
        .filter_map(|e| match e {
            SortPhase1Event::SpillReady { slot, .. } => Some(slot),
            _ => None,
        })
        .collect();
    assert_eq!(
        slots.iter().map(|s| s.run_index().map(|ix| ix.stride())).collect::<Vec<_>>(),
        expected
    );
    assert!(slots.iter().all(|s| s.has_index_footer() == footers));
    let elig = eligibility(8, 512 * MIB, 64 * MIB);
    match (crate::sort::decide_merge_mode(Some(&elig), &slots, residual_chunks), reason) {
        (crate::sort::MergeMode::Range { .. }, None) => {}
        (crate::sort::MergeMode::Serial(got), Some(expected)) => assert_eq!(got, expected),
        (mode, expected) => panic!("mode {mode:?}, expected {expected:?}"),
    }
}

/// A run whose index overflows a `u32` sum when brought to a raised stride
/// goes unindexed instead of failing the sort: raised while the run is open,
/// `close_run` writes no footer; raised only at `AllAnnounced`, the footer
/// (written at the live stride) is not attached. Either way the
/// announcement succeeds and the merge is serial. Columns: the second
/// block's stamped stride, the announced stride, whether the run keeps a
/// footer.
#[rstest]
#[case::raised_while_open(2, 2, false)]
#[case::raised_at_the_announcement(1, 2, true)]
fn an_index_overflowing_the_raised_stride_leaves_the_run_unindexed(
    #[case] second_stride: u32,
    #[case] announced_stride: u32,
    #[case] footer: bool,
) {
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    let kind = SpillKeyKind::Coordinate;
    w.process_event(run_block(0, kind, &[u32::MAX - 10], false, 1)).unwrap();
    w.process_event(run_block(0, kind, &[11], true, second_stride)).unwrap();
    assert_eq!(live_footers(&w), [footer]);
    w.process_event(SpillBlockEvent::AllAnnounced {
        ordinal: 0,
        slot_count: 1,
        memory_chunk_count: 0,
        total_records: 2,
        index_stride: announced_stride,
    })
    .unwrap();
    let Some(SortPhase1Event::SpillReady { slot, .. }) = w.outbox.pop_front() else {
        panic!("expected SpillReady");
    };
    assert_eq!(slot.has_index_footer(), footer);
    assert!(slot.run_index().is_none(), "the run is unindexed");
    let elig = eligibility(8, 512 * MIB, 64 * MIB);
    assert!(matches!(
        crate::sort::decide_merge_mode(Some(&elig), &[slot], 0),
        crate::sort::MergeMode::Serial(_)
    ));
}

/// A footer whose key width is not the runs' key kind's fails the
/// announcement with `InvalidData` naming the run, rather than being read
/// as keys of the wrong width. Columns: the run's key kind, the footer's
/// key width, a needle of the error (empty: accepted).
#[rstest]
#[case::matching_width(SpillKeyKind::Coordinate, 8, "")]
#[case::template_width_on_a_coordinate_run(SpillKeyKind::Coordinate, 24, "24-byte keys")]
#[case::coordinate_width_on_a_template_run(SpillKeyKind::TemplateK40, 8, "8-byte keys")]
fn a_footer_of_the_wrong_key_width_fails_the_announcement(
    #[case] kind: SpillKeyKind,
    #[case] width: u16,
    #[case] needle: &str,
) {
    let (mut w, _dir) = make_writer(SpillCodec::Zstd);
    let mut index = RunIndex::new(width);
    index.push_frame(10, 1000, 1, &vec![0u8; usize::from(width)]);
    w.process_event(SpillBlockEvent::Block {
        ordinal: 0,
        file_id: 0,
        key_kind: kind,
        is_last_in_file: true,
        records_ingested_so_far: 1,
        bytes: vec![0u8; 10],
        frames: Some(index),
        index_stride: 1,
    })
    .unwrap();
    let path = w.runs.runs()[0].path.clone();
    match w.process_event(announced(1, 0, 0)) {
        Ok(()) => assert!(needle.is_empty(), "expected an error naming {needle:?}"),
        Err(e) => {
            assert!(!needle.is_empty(), "unexpected error: {e}");
            assert_eq!(e.kind(), io::ErrorKind::InvalidData, "{e}");
            let msg = e.to_string();
            assert!(msg.contains(needle) && msg.contains(&path.display().to_string()), "{msg}");
        }
    }
}
