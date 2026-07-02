//! Go/no-go microbenchmark for the parallel-inflate redesign (increment 3b.0).
//!
//! The redesign's whole wall-clock premise is that the *serial* Phase-1 ingest
//! step gets cheaper per record by NOT copying each record's bytes: today's
//! `SortBlockBuffer` does `BoundaryState::scan` → `CoordinateChunkSorter::push`,
//! and `push` memcpys every record body into the sorter's own arena (~the
//! 121 ns/rec the profile flagged). The redesign instead builds a lightweight
//! `(key, offset, len)` ref over bytes already in a shared arena — no copy.
//!
//! This bench measures EXACTLY that per-record delta with today's public API —
//! no new pipeline steps, no arena, no visibility changes:
//!   * `scan_copy_push`  — scan + `CoordinateChunkSorter::push` (copies body)
//!   * `scan_refbuild`   — scan + `extract_coordinate_key_inline` + `RecordRef::new` (no copy)
//!
//! It is a NECESSARY-condition gate, not the full wall proof: it shows whether
//! removing the copy cuts serial per-record cost. If `scan_refbuild` is not
//! materially faster than `scan_copy_push`, the redesign cannot beat legacy on
//! wall and we stop before building the pipeline. If it is faster, the full wall
//! win still has to be confirmed by the end-to-end gp3 measurement (the overlap
//! the freeze-per-run model trades away is not captured here).
//!
//!     cargo bench -p fgumi-pipeline-io --bench serial_ingest

use criterion::{BatchSize, Criterion, Throughput, criterion_group, criterion_main};
use fgumi_pipeline_io::boundaries::BoundaryState;
use fgumi_raw_bam::testutil::make_bam_bytes;
use fgumi_sort::{RawExternalSorter, RecordRef, SortOrder, extract_coordinate_key_inline};
use noodles::sam::Header;

/// Build a decompressed-BAM-style buffer of `n` framed records: each record is
/// `[block_size: u32 LE][body]`, exactly what `BoundaryState::scan` walks. Bodies
/// are ~100 bp aligned reads at scattered positions on tid 0 (realistic per-record
/// size; the sort would do real reordering work).
fn synth_framed(n: usize) -> Vec<u8> {
    let mut buf = Vec::new();
    for i in 0..n {
        let pos = (i as u64).wrapping_mul(2_654_435_761) % 5_000_000;
        let name = format!("r{i:08}");
        let body = make_bam_bytes(0, pos as i32, 0, name.as_bytes(), &[], 100, -1, -1, &[]);
        let block_size = u32::try_from(body.len()).expect("record fits u32");
        buf.extend_from_slice(&block_size.to_le_bytes());
        buf.extend_from_slice(&body);
    }
    buf
}

fn bench_serial_ingest(c: &mut Criterion) {
    const N: usize = 200_000;
    let framed = synth_framed(N);
    let n_ref = 1u32;

    // Sanity: both paths see the same record count (so the throughput is over the
    // same work). Computed once, untimed.
    {
        let mut bs = BoundaryState::new_no_header();
        let (offsets, _range) = bs.scan(&framed).expect("scan");
        assert_eq!(offsets.len().saturating_sub(1), N, "scan must yield N records");
    }

    let mut group = c.benchmark_group("serial_ingest");
    group.throughput(Throughput::Elements(N as u64));
    group.sample_size(20);

    // OLD path: scan + CoordinateChunkSorter::push — `push` copies each body into
    // the sorter's arena (plus key-extract + internal ref-push).
    group.bench_function("scan_copy_push", |b| {
        b.iter_batched(
            || {
                RawExternalSorter::new(SortOrder::Coordinate)
                    .threads(1)
                    .into_coordinate_chunk_sorter(&Header::default())
                    .expect("build coordinate chunk sorter")
            },
            |mut sorter| {
                let mut bs = BoundaryState::new_no_header();
                let (offsets, range) = bs.scan(&framed).expect("scan");
                let recs = bs.records_bytes(range);
                for w in offsets.windows(2) {
                    let body = &recs[w[0] + 4..w[1]];
                    sorter.push(body).expect("push");
                }
                std::hint::black_box(&mut sorter);
            },
            BatchSize::PerIteration,
        );
    });

    // NEW path: scan + key-extract + ref-build — NO body copy. This is the
    // per-record work the redesign's serial FindBoundariesAndSort step does.
    group.bench_function("scan_refbuild", |b| {
        b.iter_batched(
            || Vec::<RecordRef>::with_capacity(N),
            |mut refs| {
                refs.clear();
                let mut bs = BoundaryState::new_no_header();
                let (offsets, range) = bs.scan(&framed).expect("scan");
                let recs = bs.records_bytes(range);
                for w in offsets.windows(2) {
                    let body = &recs[w[0] + 4..w[1]];
                    let key = extract_coordinate_key_inline(body, n_ref);
                    // offset/len point at the body (prefix skipped) — the redesign's
                    // arena ref representation; here the offset is illustrative.
                    let offset = u64::try_from(w[0] + 4).expect("offset fits u64");
                    let len = u32::try_from(w[1] - w[0] - 4).expect("len fits u32");
                    refs.push(RecordRef::new(key, offset, len));
                }
                std::hint::black_box(&mut refs);
            },
            BatchSize::PerIteration,
        );
    });

    group.finish();
}

criterion_group!(benches, bench_serial_ingest);
criterion_main!(benches);
