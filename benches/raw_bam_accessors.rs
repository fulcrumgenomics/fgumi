#![deny(unsafe_code)]

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use fgumi_raw_bam::RawRecordView;
use fgumi_raw_bam::testutil::{encode_op, make_bam_bytes};

fn bench_view_accessors(c: &mut Criterion) {
    let bytes = make_bam_bytes(
        3,
        100_000,
        0x42,
        b"read_name_example",
        &[encode_op(0, 150)],
        150,
        5,
        100_300,
        b"NMc\x05MDZ150\0",
    );

    c.bench_function("view::flags", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes));
            black_box(v.flags())
        });
    });
    c.bench_function("view::pos", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes));
            black_box(v.pos())
        });
    });
    c.bench_function("view::read_name", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes));
            black_box(v.read_name())
        });
    });
    c.bench_function("view::cigar_ops_iter_sum", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes));
            let mut s = 0u32;
            for op in v.cigar_ops_iter() {
                s = s.wrapping_add(op);
            }
            black_box(s)
        });
    });
    c.bench_function("view::quality_scores_sum", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes));
            let s: u32 = v.quality_scores().iter().map(|&q| u32::from(q)).sum();
            black_box(s)
        });
    });
    c.bench_function("view::tags_find_int", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes));
            black_box(v.tags().find_int(b"NM"))
        });
    });

    // Shorter seq (below SIMD threshold)
    let bytes_short = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 16)], 16, -1, -1, &[]);
    c.bench_function("view::sequence_vec_16bp", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes_short));
            black_box(v.sequence_vec());
        });
    });

    // 150 bp read (typical Illumina short read; SIMD active)
    let bytes_150 = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 150)], 150, -1, -1, &[]);
    c.bench_function("view::sequence_vec_150bp", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes_150));
            black_box(v.sequence_vec());
        });
    });

    // 300 bp read (2× typical; SIMD loop iterations dominate)
    let bytes_300 = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 300)], 300, -1, -1, &[]);
    c.bench_function("view::sequence_vec_300bp", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes_300));
            black_box(v.sequence_vec());
        });
    });
}

criterion_group!(benches, bench_view_accessors);
criterion_main!(benches);
