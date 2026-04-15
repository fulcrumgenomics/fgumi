use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use fgumi_raw_bam::RawRecord;
use fgumi_raw_bam::testutil::{encode_op, make_bam_bytes};

fn bench_length_changing(c: &mut Criterion) {
    let baseline = || {
        let bytes = make_bam_bytes(
            0,
            0,
            0,
            b"my_long_read_name_typical_size",
            &[encode_op(0, 150)],
            150,
            -1,
            -1,
            b"NMi\x05\x00\x00\x00",
        );
        RawRecord::from(bytes)
    };

    c.bench_function("set_read_name_grow_then_shrink", |b| {
        b.iter(|| {
            let mut rec = baseline();
            rec.set_read_name(black_box(b"a_completely_different_longer_name_xxx"));
            rec.set_read_name(black_box(b"x"));
        });
    });

    c.bench_function("set_cigar_ops_replace_with_3_ops", |b| {
        b.iter(|| {
            let mut rec = baseline();
            rec.set_cigar_ops(black_box(&[encode_op(0, 50), encode_op(2, 1), encode_op(0, 100)]));
        });
    });

    c.bench_function("set_sequence_and_qualities_grow", |b| {
        let seq = vec![b'A'; 200];
        let qual = vec![30u8; 200];
        b.iter(|| {
            let mut rec = baseline();
            rec.set_sequence_and_qualities(black_box(&seq), black_box(&qual));
        });
    });

    c.bench_function("set_sequence_and_qualities_same_length", |b| {
        let seq = vec![b'C'; 150];
        let qual = vec![25u8; 150];
        b.iter(|| {
            let mut rec = baseline();
            rec.set_sequence_and_qualities(black_box(&seq), black_box(&qual));
        });
    });
}

criterion_group!(benches, bench_length_changing);
criterion_main!(benches);
