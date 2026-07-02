//! Microbenchmark for the P6 `SortBuffer` chunk sorters' sort + materialize
//! step (`take_sorted_chunk`).
//!
//! This isolates the cost the A/B smoke campaign flagged: `take_sorted_chunk`
//! materialized a full owned `Vec<(K, RawRecord)>` (one heap copy + malloc per
//! record) while the source `RecordBuffer` arena was still alive — doubling
//! peak memory and paying a full per-record copy. The fix routes coordinate /
//! template through the zero-copy `InMemoryChunk` (arena-move). Run before and
//! after to confirm the materialize step gets cheaper:
//!
//!     cargo bench -p fgumi-sort --bench chunk_sorter

use criterion::{BatchSize, Criterion, Throughput, criterion_group, criterion_main};
use fgumi_raw_bam::testutil::make_bam_bytes;
use fgumi_sort::{RawExternalSorter, SortOrder};
use noodles::sam::Header;

/// Build `n` aligned records at scattered positions on tid 0 so the sort does
/// real reordering work; ~100 bp reads ≈ a realistic per-record size.
fn synth_records(n: usize) -> Vec<Vec<u8>> {
    (0..n)
        .map(|i| {
            let pos = (i as u64).wrapping_mul(2_654_435_761) % 5_000_000;
            let name = format!("r{i:08}");
            make_bam_bytes(0, pos as i32, 0, name.as_bytes(), &[], 100, -1, -1, &[])
        })
        .collect()
}

fn bench_materialize(c: &mut Criterion) {
    const N: usize = 500_000;
    let records = synth_records(N);

    let mut group = c.benchmark_group("chunk_sorter_materialize");
    group.throughput(Throughput::Elements(N as u64));
    group.sample_size(10);

    group.bench_function("coordinate_500k", |b| {
        b.iter_batched(
            || {
                let mut sorter = RawExternalSorter::new(SortOrder::Coordinate)
                    .threads(4)
                    .into_coordinate_chunk_sorter(&Header::default())
                    .expect("build coordinate chunk sorter");
                for r in &records {
                    sorter.push(r).expect("push");
                }
                sorter
            },
            |mut sorter| std::hint::black_box(sorter.take_sorted_chunk()),
            BatchSize::PerIteration,
        );
    });

    group.bench_function("template_coordinate_500k", |b| {
        b.iter_batched(
            || {
                let mut sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
                    .threads(4)
                    .into_template_chunk_sorter(&Header::default())
                    .expect("build template chunk sorter");
                for r in &records {
                    sorter.push(r).expect("push");
                }
                sorter
            },
            |mut sorter| std::hint::black_box(sorter.take_sorted_chunk()),
            BatchSize::PerIteration,
        );
    });

    group.finish();
}

criterion_group!(benches, bench_materialize);
criterion_main!(benches);
