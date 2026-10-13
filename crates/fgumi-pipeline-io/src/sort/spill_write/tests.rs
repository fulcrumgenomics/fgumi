//! Unit tests for `SpillWrite`'s `StepCtx`-free core: per-file demux, codec
//! magic/trailer bracketing, runs held closed until `AllAnnounced`, the
//! consolidation-job lifecycle under `--max-temp-files`, and event mapping.
//! The jobs themselves run through the real tail in `consolidate_runs`' tests;
//! byte-exact readback of a consolidated sort is gated end-to-end by the sort
//! command's consolidation parity tests.

use super::*;
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
    }
}

fn announced(slot_count: u32, memory_chunk_count: u32, total_records: u64) -> SpillBlockEvent {
    SpillBlockEvent::AllAnnounced { ordinal: 0, slot_count, memory_chunk_count, total_records }
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
            }
        })
        .collect();
    events.push(SpillBlockEvent::AllAnnounced {
        ordinal: u64::from(N_FILES),
        slot_count: N_FILES,
        memory_chunk_count: 0,
        total_records: u64::from(N_FILES),
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
    let Some(SortPhase1Event::ConsolidationStart { job_id, inputs, key_kind }) =
        w.outbox.pop_front()
    else {
        panic!("expected a job")
    };
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
    registry.complete(job_id, bytes);
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
    registry.complete(7, 4096);
    let (got, waited) =
        rx.recv_timeout(Duration::from_secs(5)).expect("the waiter must not sit out its 10 s park");
    assert_eq!(got, (7, 4096));
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
        r.complete(job_id, bytes);
    }));
    assert!(w.poll_registry().unwrap(), "the re-check after register_waiter finds the completion");
}
