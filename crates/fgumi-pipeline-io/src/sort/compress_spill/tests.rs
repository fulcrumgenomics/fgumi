//! Unit tests for `CompressSpill::compress_event` (the `StepCtx`-free core).
//!
//! These drive the per-event compression directly on synthetic sorted chunks,
//! covering the three event arms: `Spill` (compress to disk, slot `file_id` =
//! logical `seq`), `Residual` (pass through as a uniquely-owned `MemoryChunk`),
//! and `AllAnnounced` (verbatim forward). The end-to-end chain parity vs the
//! legacy oracle is exercised later, once `SortBuffer` + `add_sort` are wired
//! (P6 increments 3–4).

use super::*;

use fgumi_raw_bam::RawRecord;
use fgumi_sort::{RawCoordinateKey, TemplateKey};
use rstest::rstest;
use tempfile::TempDir;

use crate::sort::protocol::MemoryChunkErased;

/// Which easily-constructed sort-key variant a parameterized case exercises.
/// (The two queryname variants are the same generic `write_sorted_chunk::<K>`
/// dispatch; they get real coverage from the inc-4 chain parity test, where
/// genuine queryname keys come from BAM data rather than fragile hand
/// construction.)
#[derive(Clone, Copy)]
enum SpillVariant {
    Coordinate,
    Template,
}

type CoordRecords = Vec<(RawCoordinateKey, RawRecord)>;

/// `n` coordinate records with distinct, sized payloads so the spill file has
/// real content (and multiple zstd frames for large `n`).
#[allow(clippy::cast_possible_truncation)] // payload byte is `% 251`, always fits u8
fn coord_records(n: usize) -> CoordRecords {
    (0..n)
        .map(|i| {
            // Distinct, ascending sort keys so the byte-identity checks exercise
            // key serialization (not just record-byte order).
            let key = RawCoordinateKey { sort_key: i as u64 };
            (key, RawRecord::from(vec![(i % 251) as u8; 100 + i % 64]))
        })
        .collect()
}

/// Deterministic coordinate records keyed off `seq`, so a concurrently-written
/// spill file can be checked against an independent reference write.
#[allow(clippy::cast_possible_truncation)] // payload byte is `% 251`, always fits u8
fn coord_records_for(seq: u32) -> CoordRecords {
    let n = 40 + (seq as usize % 24);
    (0..n)
        .map(|i| {
            let byte = (seq as usize + i) % 251;
            // Distinct keys per (seq, i) so key serialization is exercised.
            let key = RawCoordinateKey { sort_key: (u64::from(seq) << 32) | i as u64 };
            (key, RawRecord::from(vec![byte as u8; 80 + i % 32]))
        })
        .collect()
}

/// Pack owned coordinate records into the zero-copy arena-backed
/// [`InMemoryChunk`] the `Coordinate` protocol variant now carries. The spill
/// file this produces must be byte-identical to a `write_sorted_chunk` of the
/// same owned records (the tests assert exactly that).
fn coord_chunk(records: CoordRecords) -> fgumi_sort::InMemoryChunk<RawCoordinateKey> {
    fgumi_sort::InMemoryChunk::from_owned_records(
        records.into_iter().map(|(k, r)| (k, r.into_inner())).collect(),
    )
}

/// `n` template-coordinate records (the other easily-constructed key variant).
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
fn template_records(n: usize) -> Vec<(TemplateKey, RawRecord)> {
    (0..n)
        .map(|i| {
            let k = TemplateKey::new(
                i as i32,
                i as i32,
                false,
                i32::MAX,
                i32::MAX,
                false,
                0,
                0,
                (0, false),
                i as u64,
                false,
            );
            (k, RawRecord::from(vec![(i % 251) as u8; 120 + i % 48]))
        })
        .collect()
}

/// Wrap owned template records into the arena-backed `InMemoryChunk<TemplateKey>`
/// the `TemplateCoordinate` protocol variant now carries — the template analogue
/// of [`coord_chunk`]. The spill file this produces must be byte-identical to a
/// `write_sorted_chunk` of the same owned records (the test asserts exactly that).
fn template_chunk(
    records: Vec<(TemplateKey, RawRecord)>,
) -> fgumi_sort::InMemoryChunk<TemplateKey> {
    fgumi_sort::InMemoryChunk::from_owned_records(
        records.into_iter().map(|(k, r)| (k, r.as_ref().to_vec())).collect(),
    )
}

/// A `CompressSpill` over a single temp dir with an always-ample free-space
/// probe (deterministic, no dependency on the host's real free space), plus the
/// `TempDir` guard that must outlive the opened spill slots.
fn make_step(codec: SpillCodec, compression: u32) -> (CompressSpill, TempDir) {
    let dir = TempDir::new().expect("temp dir");
    let alloc =
        TmpDirAllocator::with_probe(vec![dir.path().to_path_buf()], Box::new(|_| Ok(u64::MAX)), 0)
            .expect("allocator builds");
    // The test holds `dir` alive itself, so the step's RAII set is empty.
    let step = CompressSpill::new(
        Arc::new(Mutex::new(alloc)),
        codec,
        compression,
        8 * 1024 * 1024,
        Arc::new(Vec::new()),
    );
    (step, dir)
}

/// A `Spill` event compresses the chunk to disk, opens a slot whose `file_id`
/// equals the logical `seq` (not write order), and writes bytes byte-identical
/// to a direct `write_sorted_chunk` of the same records — for every sort-key
/// variant the dispatch handles.
#[rstest]
#[case::coordinate(SpillVariant::Coordinate, 5, 1234)]
#[case::template(SpillVariant::Template, 2, 300)]
fn spill_event_writes_chunk_with_seq_file_id_byte_identical(
    #[case] variant: SpillVariant,
    #[case] seq: u32,
    #[case] records_ingested: u64,
) {
    let (step, dir) = make_step(SpillCodec::Zstd, 3);

    // Write an independent reference from the same (deterministic) records, then
    // move those records into the event chunk — borrow-then-move avoids cloning
    // the chunk (`MemoryChunkErased` is deliberately not `Clone`).
    let reference = dir.path().join("reference.keyed");
    let chunk = match variant {
        SpillVariant::Coordinate => {
            let recs = coord_records(500);
            fgumi_sort::write_sorted_chunk(&reference, SpillCodec::Zstd, 3, &recs).unwrap();
            MemoryChunkErased::Coordinate(coord_chunk(recs))
        }
        SpillVariant::Template => {
            let recs = template_records(300);
            fgumi_sort::write_sorted_chunk(&reference, SpillCodec::Zstd, 3, &recs).unwrap();
            MemoryChunkErased::TemplateCoordinate(fgumi_sort::TemplateMemChunk::K40(
                template_chunk(recs),
            ))
        }
    };

    let event = SortChunkEvent::Spill { seq, chunk, records_ingested_so_far: records_ingested };
    let forwarded = step.compress_event(event).expect("compress Spill event");

    let SortPhase1Event::SpillReady { slot, path, records_ingested_so_far } = forwarded else {
        panic!("Spill must forward as SpillReady");
    };
    assert_eq!(records_ingested_so_far, records_ingested, "records_ingested must be propagated");
    assert_eq!(slot.file_id, seq, "slot file_id must equal the logical spill seq");
    assert_eq!(slot.codec, SpillCodec::Zstd, "codec must be detected from the written magic");
    assert!(path.exists(), "spill file must exist");
    assert!(path.starts_with(dir.path()), "spill file must live under the allocated temp dir");
    assert!(
        path.file_name().unwrap().to_str().unwrap().contains(&format!("{seq:04}")),
        "spill file should be named by its seq"
    );

    // Byte-identical to a direct write of the same records (the step must add no
    // framing of its own — it delegates straight to write_sorted_chunk).
    assert_eq!(
        std::fs::read(&path).unwrap(),
        std::fs::read(&reference).unwrap(),
        "CompressSpill output must match write_sorted_chunk byte-for-byte"
    );
}

/// Distinct `seq` values yield distinct slots/paths, so concurrent workers never
/// collide and the merge can tie-break by `file_id`.
#[test]
fn distinct_seqs_yield_distinct_file_ids_and_paths() {
    let (step, _dir) = make_step(SpillCodec::Bgzf, 1);
    let mut paths = Vec::new();
    for seq in [0u32, 1, 7, 42] {
        let event = SortChunkEvent::Spill {
            seq,
            chunk: MemoryChunkErased::Coordinate(coord_chunk(coord_records(50))),
            records_ingested_so_far: u64::from(seq),
        };
        let SortPhase1Event::SpillReady { slot, path, .. } =
            step.compress_event(event).expect("compress")
        else {
            panic!("expected SpillReady");
        };
        assert_eq!(slot.file_id, seq);
        paths.push(path);
    }
    let unique: std::collections::HashSet<_> = paths.iter().collect();
    assert_eq!(unique.len(), paths.len(), "every spill path must be unique");
}

/// A `Residual` event passes through as a `MemoryChunk` wrapping a uniquely-owned
/// `Arc` (the invariant `SortMerge`'s `Arc::try_unwrap` relies on).
#[test]
fn residual_event_passes_through_as_unique_memory_chunk() {
    let (step, _dir) = make_step(SpillCodec::Zstd, 3);
    let records = coord_records(120);

    let event = SortChunkEvent::Residual {
        chunk: MemoryChunkErased::Coordinate(coord_chunk(records)),
        records_ingested_so_far: 99,
    };
    let forwarded = step.compress_event(event).expect("compress Residual event");

    let SortPhase1Event::MemoryChunk { chunk, records_ingested_so_far } = forwarded else {
        panic!("Residual must forward as MemoryChunk");
    };
    assert_eq!(records_ingested_so_far, 99);
    let inner = Arc::try_unwrap(chunk).unwrap_or_else(|_| {
        panic!("MemoryChunk Arc must be uniquely owned (SortMerge unwraps it)")
    });
    assert_eq!(inner.len(), 120, "residual chunk content must pass through intact");
}

/// `AllAnnounced` forwards verbatim — the counts `SortBuffer` computed are the
/// completion target `SortMerge` keys off of.
#[test]
fn all_announced_passes_through_verbatim() {
    let (step, _dir) = make_step(SpillCodec::Zstd, 3);
    let event =
        SortChunkEvent::AllAnnounced { slot_count: 4, memory_chunk_count: 1, total_records: 5000 };
    let forwarded = step.compress_event(event).expect("compress AllAnnounced");
    let SortPhase1Event::AllAnnounced { slot_count, memory_chunk_count, total_records } = forwarded
    else {
        panic!("AllAnnounced must forward as AllAnnounced");
    };
    assert_eq!((slot_count, memory_chunk_count, total_records), (4, 1, 5000));
}

/// Multiple cloned workers (the `new_worker_copy` fan-out) sharing one
/// allocator can compress spill chunks concurrently: every file is unique,
/// carries the right `file_id`, and is byte-identical to an independent
/// reference write — i.e. the shared `Mutex<TmpDirAllocator>` is the only
/// cross-worker state and it serializes cleanly.
#[test]
fn parallel_workers_share_allocator_without_collision() {
    let (base, dir) = make_step(SpillCodec::Zstd, 3);
    let total_seqs: u32 = 64;
    let num_workers = 8;

    // Each worker gets its own clone (fresh `held`), all sharing `base`'s
    // allocator Arc — exactly how the framework materializes Parallel workers.
    let results = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..num_workers)
            .map(|w| {
                let worker = base.clone();
                scope.spawn(move || {
                    let mut produced = Vec::new();
                    let mut seq = w;
                    while seq < total_seqs {
                        let records = coord_records_for(seq);
                        let event = SortChunkEvent::Spill {
                            seq,
                            chunk: MemoryChunkErased::Coordinate(coord_chunk(records)),
                            records_ingested_so_far: u64::from(seq),
                        };
                        let SortPhase1Event::SpillReady { slot, path, .. } =
                            worker.compress_event(event).expect("worker compress")
                        else {
                            panic!("expected SpillReady");
                        };
                        produced.push((seq, slot.file_id, path));
                        seq += num_workers;
                    }
                    produced
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().expect("worker thread")).collect::<Vec<_>>()
    });

    let mut all: Vec<(u32, u32, std::path::PathBuf)> = results.into_iter().flatten().collect();
    all.sort_by_key(|(seq, _, _)| *seq);

    assert_eq!(all.len(), total_seqs as usize, "every seq must produce exactly one spill file");

    let unique_paths: std::collections::HashSet<_> = all.iter().map(|(_, _, p)| p).collect();
    assert_eq!(unique_paths.len(), all.len(), "no two workers may collide on a spill path");

    for (seq, file_id, path) in &all {
        assert_eq!(file_id, seq, "file_id must equal the logical seq, not write order");
        // Independent reference write of this seq's deterministic records.
        let reference = dir.path().join(format!("ref_{seq:04}.keyed"));
        fgumi_sort::write_sorted_chunk(&reference, SpillCodec::Zstd, 3, &coord_records_for(*seq))
            .unwrap();
        assert_eq!(
            std::fs::read(path).unwrap(),
            std::fs::read(&reference).unwrap(),
            "concurrently-written chunk {seq} must match its reference byte-for-byte"
        );
    }
}

/// An empty residual chunk is still a valid passthrough (zero-record fast path).
#[test]
fn empty_residual_passes_through() {
    let (step, _dir) = make_step(SpillCodec::Zstd, 3);
    let event = SortChunkEvent::Residual {
        chunk: MemoryChunkErased::Coordinate(coord_chunk(Vec::new())),
        records_ingested_so_far: 0,
    };
    let SortPhase1Event::MemoryChunk { chunk, .. } =
        step.compress_event(event).expect("compress empty residual")
    else {
        panic!("expected MemoryChunk");
    };
    let inner =
        Arc::try_unwrap(chunk).unwrap_or_else(|_| panic!("MemoryChunk Arc must be uniquely owned"));
    assert!(inner.is_empty());
}
