//! Unit tests for `SpillWrite::process_event` (the `StepCtx`-free core): per-file
//! demux, codec magic/trailer bracketing, finalization on `is_last_in_file`, and
//! event mapping. Byte-exact readback of the assembled file is gated end-to-end
//! by the full-sort parity test (the production reader is crate-private to
//! `fgumi-sort`).

use super::*;
use crate::sort::protocol::MemoryChunkErased;
use fgumi_sort::{InMemoryChunk, RawCoordinateKey, SpillBlockCompressor};
use tempfile::TempDir;

/// Build a `SpillWrite` writing into a fresh temp dir; returns it plus the dir
/// (kept alive by the caller) so written files can be inspected.
fn make_writer(codec: SpillCodec) -> (SpillWrite, TempDir) {
    let dir = TempDir::new().unwrap();
    let alloc = TmpDirAllocator::new(vec![dir.path().to_path_buf()]).unwrap();
    let writer = SpillWrite::new(Arc::new(Mutex::new(alloc)), codec, 1 << 20, Arc::new(Vec::new()));
    (writer, dir)
}

/// Kernel-compress one raw block for `codec`, mirroring `SpillCompress`.
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

#[test]
fn non_last_block_opens_file_emits_nothing_last_block_emits_spill_ready() {
    for codec in [SpillCodec::Zstd, SpillCodec::Bgzf] {
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
