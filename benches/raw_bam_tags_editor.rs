use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use fgumi_raw_bam::RawRecord;
use fgumi_raw_bam::testutil::{encode_op, make_bam_bytes};

fn bench_editor(c: &mut Criterion) {
    c.bench_function("editor::update_int_same_type", |b| {
        let bytes =
            make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, b"NMi\x05\x00\x00\x00");
        let mut rec = RawRecord::from(bytes);
        b.iter(|| {
            let mut ed = rec.tags_editor();
            ed.update_int(b"NM", black_box(7));
        });
    });

    c.bench_function("editor::update_int_resize", |b| {
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, b"NMc\x05");
        let mut rec = RawRecord::from(bytes);
        b.iter(|| {
            let mut ed = rec.tags_editor();
            ed.update_int(b"NM", black_box(100_000));
            ed.update_int(b"NM", 5);
        });
    });

    c.bench_function("editor::append_then_remove", |b| {
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        let mut rec = RawRecord::from(bytes);
        b.iter(|| {
            let mut ed = rec.tags_editor();
            ed.append_string(b"XX", black_box(b"hello"));
            ed.remove(b"XX");
        });
    });

    c.bench_function("editor::update_float_in_place", |b| {
        let mut aux = b"ASf".to_vec();
        aux.extend_from_slice(&1.0f32.to_le_bytes());
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &aux);
        let mut rec = RawRecord::from(bytes);
        b.iter(|| {
            let mut ed = rec.tags_editor();
            ed.update_float(b"AS", black_box(99.25));
        });
    });
}

criterion_group!(benches, bench_editor);
criterion_main!(benches);
