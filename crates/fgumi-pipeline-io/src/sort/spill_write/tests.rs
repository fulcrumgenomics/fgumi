//! Unit tests for `SpillWrite::process_event` (the `StepCtx`-free core): per-file
//! demux, codec magic/trailer bracketing, finalization on `is_last_in_file`, and
//! event mapping. Byte-exact readback of the assembled file is gated end-to-end
//! by the full-sort parity test (the production reader is crate-private to
//! `fgumi-sort`).

use super::*;
use crate::sort::protocol::MemoryChunkErased;
use fgumi_sort::{InMemoryChunk, RawCoordinateKey, SpillBlockCompressor};
use rstest::rstest;
use tempfile::TempDir;

/// Build a `SpillWrite` writing into a fresh temp dir; returns it plus the dir
/// (kept alive by the caller) so written files can be inspected.
fn make_writer(codec: SpillCodec) -> (SpillWrite, TempDir) {
    let dir = TempDir::new().unwrap();
    let alloc = TmpDirAllocator::new(vec![dir.path().to_path_buf()]).unwrap();
    let writer = SpillWrite::new(Arc::new(Mutex::new(alloc)), codec, 1 << 20, Arc::new(Vec::new()));
    (writer, dir)
}

/// Kernel-compress one raw block for `codec`, mirroring `SpillBlockCompress`.
fn compress(codec: SpillCodec, raw: &[u8]) -> Vec<u8> {
    SpillBlockCompressor::new(codec, 1).unwrap().compress_block(raw).unwrap()
}

fn block(codec: SpillCodec, file_id: u32, is_last: bool, raw: &[u8]) -> SpillBlockEvent {
    SpillBlockEvent::Block {
        ordinal: 0,
        file_id,
        is_last_in_file: is_last,
        records_ingested_so_far: 42,
        bytes: compress(codec, raw),
    }
}

#[rstest]
#[case(SpillCodec::Zstd)]
#[case(SpillCodec::Bgzf)]
fn non_last_block_opens_file_emits_nothing_last_block_emits_spill_ready(#[case] codec: SpillCodec) {
    let (mut w, dir) = make_writer(codec);
    // First (non-last) block: file opens, no event.
    let out = w.process_event(block(codec, 5, false, &[1u8; 32])).unwrap();
    assert!(out.is_none(), "non-last block emits no event ({codec:?})");
    assert!(w.current.is_some(), "file must be open after first block ({codec:?})");

    // Last block: trailer written, slot opened, SpillReady emitted.
    let out = w.process_event(block(codec, 5, true, &[2u8; 32])).unwrap();
    let Some(SortPhase1Event::SpillReady { slot, path, records_ingested_so_far }) = out else {
        panic!("expected SpillReady ({codec:?})");
    };
    assert!(w.current.is_none(), "file closed after last block ({codec:?})");
    assert_eq!(slot.file_id, 5, "slot file_id == logical seq ({codec:?})");
    assert_eq!(records_ingested_so_far, 42);
    assert!(path.exists(), "spill file exists ({codec:?})");
    assert!(path.starts_with(dir.path()), "spill file under temp dir ({codec:?})");
    assert_eq!(slot.codec, codec, "codec detected from written magic ({codec:?})");
}

#[test]
fn distinct_file_ids_produce_distinct_files() {
    let codec = SpillCodec::Zstd;
    let (mut w, _dir) = make_writer(codec);
    // File 0 (single block), then file 1 (single block) — contiguous per file.
    let r0 = w.process_event(block(codec, 0, true, &[7u8; 16])).unwrap().unwrap();
    let r1 = w.process_event(block(codec, 1, true, &[8u8; 16])).unwrap().unwrap();
    let (
        SortPhase1Event::SpillReady { path: p0, slot: s0, .. },
        SortPhase1Event::SpillReady { path: p1, slot: s1, .. },
    ) = (r0, r1)
    else {
        panic!("expected two SpillReady events");
    };
    assert_ne!(p0, p1, "distinct file_ids must yield distinct paths");
    assert_eq!(s0.file_id, 0);
    assert_eq!(s1.file_id, 1);
}

#[test]
fn block_for_wrong_file_id_while_open_errors() {
    let codec = SpillCodec::Zstd;
    let (mut w, _dir) = make_writer(codec);
    // Open file 0 with a non-last block, then feed a block for file 1 — a
    // contiguity violation that must fail loud, not silently corrupt file 0.
    w.process_event(block(codec, 0, false, &[1u8; 16])).unwrap();
    // `SortPhase1Event` is not `Debug`, so match instead of `unwrap_err`.
    match w.process_event(block(codec, 1, false, &[2u8; 16])) {
        Err(err) => {
            assert!(
                err.to_string().contains("contiguous"),
                "expected contiguity error, got: {err}"
            );
        }
        Ok(_) => panic!("a block for a different open file_id must error"),
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
        Ok(_) => panic!("residual while a spill file is open must error"),
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
    let out = w
        .process_event(SpillBlockEvent::Residual { ordinal: 0, chunk, records_ingested_so_far: 3 })
        .unwrap();
    let Some(SortPhase1Event::MemoryChunk { chunk, records_ingested_so_far }) = out else {
        panic!("expected MemoryChunk");
    };
    assert_eq!(records_ingested_so_far, 3);
    assert_eq!(Arc::strong_count(&chunk), 1, "residual chunk wrapped in a fresh unique Arc");

    let out = w
        .process_event(SpillBlockEvent::AllAnnounced {
            ordinal: 1,
            slot_count: 4,
            memory_chunk_count: 1,
            total_records: 500,
        })
        .unwrap();
    assert!(matches!(
        out,
        Some(SortPhase1Event::AllAnnounced {
            slot_count: 4,
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
    let out = w.process_event(block(SpillCodec::Zstd, 3, false, &[7u8; 16])).unwrap();
    assert!(out.is_none());
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
    let opened = w.open_file(9).expect("first open succeeds");
    drop(opened);
    assert!(dir.path().join("chunk_0009.keyed").exists(), "spill file is created eagerly");

    // A reused file_id must fail closed rather than truncate the existing file:
    // silently overwriting a spill would drop records from the merge.
    // `OpenSpill` is not `Debug`, so match instead of using `expect_err`.
    match w.open_file(9) {
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
        Ok(_) => panic!("AllAnnounced while a spill file is open must error"),
    }
}
