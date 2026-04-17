//! Consolidated benchmarks for the `fgumi-raw-bam` crate.
//!
//! Groups:
//!   - `view::*` — zero-copy header/CIGAR/seq/qual accessors via `RawRecordView`
//!   - `tags::*` — tag iteration and batched string extraction via `RawTagsView`
//!   - `editor::*` — length-changing tag mutations via `RawTagsEditor`
//!   - `mut::*` — fixed-length tag mutations via `RawTagsMut`
//!   - `length::*` — length-changing record setters (read name, CIGAR, seq/qual)
//!   - `seq::*` — per-base / per-qual writes
//!   - `flags_read`, `update_int_tag`, `set_read_name`, `set_cigar_ops`,
//!     `set_sequence_and_qualities` — side-by-side raw-byte vs noodles `RecordBuf`
//!     groups (each containing `raw_bytes` and `recordbuf` members)
#![deny(unsafe_code)]

use std::hint::black_box;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use noodles::sam::alignment::record::cigar::op::{Kind, Op};
use noodles::sam::alignment::record_buf::Cigar as CigarBuf;

use fgumi_raw_bam::raw_records_to_record_bufs;
use fgumi_raw_bam::testutil::{encode_op, make_bam_bytes};
use fgumi_raw_bam::{RawRecord, RawRecordView};

// ============================================================================
// Shared helpers
// ============================================================================

fn baseline_bytes_mapped() -> Vec<u8> {
    make_bam_bytes(
        3,
        100_000,
        0x42,
        b"read_name_example",
        &[encode_op(0, 150)],
        150,
        5,
        100_300,
        b"NMc\x05MDZ150\0",
    )
}

/// Seven-tag aux section: NM, MD, RG, CR, MC, MI, AS.
fn multi_tag_aux() -> Vec<u8> {
    let mut aux = Vec::new();
    aux.extend_from_slice(b"NMi\x05\x00\x00\x00");
    aux.extend_from_slice(b"MDZ150\0");
    aux.extend_from_slice(b"RGZorg1\0");
    aux.extend_from_slice(b"CRZcellABC\0");
    aux.extend_from_slice(b"MCZ150M\0");
    aux.extend_from_slice(b"MIi\x05\x00\x00\x00");
    aux.extend_from_slice(b"ASi\x32\x00\x00\x00");
    aux
}

fn multi_tag_bytes() -> Vec<u8> {
    let aux = multi_tag_aux();
    make_bam_bytes(
        3,
        100_000,
        0x42,
        b"read_name_example",
        &[encode_op(0, 150)],
        150,
        5,
        100_300,
        &aux,
    )
}

fn baseline_bytes_unmapped() -> Vec<u8> {
    make_bam_bytes(
        -1,
        -1,
        0x42,
        b"read_name_example",
        &[encode_op(0, 150)],
        150,
        -1,
        -1,
        b"NMi\x05\x00\x00\x00",
    )
}

fn bi_array_aux(values: &[i32]) -> Vec<u8> {
    let mut aux = b"XYBi".to_vec();
    let count = u32::try_from(values.len()).expect("array length fits in u32");
    aux.extend_from_slice(&count.to_le_bytes());
    for v in values {
        aux.extend_from_slice(&v.to_le_bytes());
    }
    aux
}

fn decode_one(bytes: &[u8]) -> noodles::sam::alignment::RecordBuf {
    let mut bufs = raw_records_to_record_bufs(&[bytes.to_vec()]).expect("decode should succeed");
    bufs.pop().expect("should have one record")
}

fn encode_one(buf: &noodles::sam::alignment::RecordBuf) -> Vec<u8> {
    use noodles::sam::alignment::io::Write as AlignmentWrite;
    let header = noodles::sam::Header::default();
    let mut out: Vec<u8> = Vec::with_capacity(512);
    let mut writer = noodles::bam::io::Writer::from(&mut out);
    writer.write_header(&header).expect("write_header should succeed");
    writer.write_alignment_record(&header, buf).expect("write_alignment_record should succeed");
    // Strip 12-byte BAM file header prefix (magic + header text len + n_ref)
    // so the result matches the block_size-prefixed bytes from `make_bam_bytes`.
    out[12..].to_vec()
}

// ============================================================================
// view::* — zero-copy accessors
// ============================================================================

fn bench_view_accessors(c: &mut Criterion) {
    let bytes = baseline_bytes_mapped();

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
}

// ============================================================================
// tags::* — tag-section reads
// ============================================================================

fn bench_tags_read(c: &mut Criterion) {
    let bytes = multi_tag_bytes();

    c.bench_function("tags::iter_all", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes));
            let mut n = 0usize;
            let tags = v.tags();
            for entry in &tags {
                black_box(entry);
                n += 1;
            }
            black_box(n)
        });
    });

    c.bench_function("tags::extract_string_batch", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes));
            black_box(v.tags().extract_string_batch(b"CR"))
        });
    });
}

// ============================================================================
// editor::* — length-changing tag mutations
// ============================================================================

fn bench_tags_editor_scalar(c: &mut Criterion) {
    c.bench_function("editor::update_int_same_type", |b| {
        let bytes =
            make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, b"NMi\x05\x00\x00\x00");
        let mut rec = RawRecord::from(bytes);
        b.iter(|| {
            rec.tags_editor().update_int(b"NM", black_box(7));
        });
    });

    c.bench_function("editor::update_int_resize", |b| {
        b.iter_batched_ref(
            || {
                let bytes =
                    make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, b"NMc\x05");
                RawRecord::from(bytes)
            },
            |rec| {
                rec.tags_editor().update_int(b"NM", black_box(100_000));
            },
            BatchSize::SmallInput,
        );
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
            rec.tags_editor().update_float(b"AS", black_box(99.25));
        });
    });

    c.bench_function("editor::normalize_int_to_smallest_signed", |b| {
        b.iter_batched_ref(
            || {
                let bytes = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, b"NMi\x05\x00\x00\x00");
                RawRecord::from(bytes)
            },
            |rec| {
                rec.tags_editor().normalize_int_to_smallest_signed(b"NM");
            },
            BatchSize::SmallInput,
        );
    });
}

fn bench_tags_editor_variable(c: &mut Criterion) {
    c.bench_function("editor::update_string_same_length", |b| {
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, b"XXZabcde\0");
        let mut rec = RawRecord::from(bytes);
        b.iter(|| {
            rec.tags_editor().update_string(b"XX", black_box(b"ABCDE"));
        });
    });

    c.bench_function("editor::update_string_grow", |b| {
        b.iter_batched_ref(
            || {
                let bytes = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, b"XXZabc\0");
                RawRecord::from(bytes)
            },
            |rec| {
                rec.tags_editor().update_string(b"XX", black_box(b"a_much_longer_string"));
            },
            BatchSize::SmallInput,
        );
    });

    c.bench_function("editor::append_string", |b| {
        b.iter_batched_ref(
            || {
                let bytes = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
                RawRecord::from(bytes)
            },
            |rec| {
                rec.tags_editor().append_string(b"XX", black_box(b"hello"));
            },
            BatchSize::SmallInput,
        );
    });

    c.bench_function("editor::remove", |b| {
        b.iter_batched_ref(
            || {
                let bytes =
                    make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, b"NMi\x05\x00\x00\x00XXZhello\0");
                RawRecord::from(bytes)
            },
            |rec| {
                rec.tags_editor().remove(b"XX");
            },
            BatchSize::SmallInput,
        );
    });

    c.bench_function("editor::copy_from", |b| {
        let src_bytes = multi_tag_bytes();
        let src_tags = RawRecordView::new(&src_bytes).tags();
        b.iter_batched_ref(
            || {
                let bytes = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
                RawRecord::from(bytes)
            },
            |rec| {
                rec.tags_editor().copy_from(src_tags, black_box(&[b"NM"]));
            },
            BatchSize::SmallInput,
        );
    });

    c.bench_function("editor::update_array_i32_same_len", |b| {
        let aux = bi_array_aux(&[1, 2, 3, 4]);
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &aux);
        let mut rec = RawRecord::from(bytes);
        let replacement: [i32; 4] = [10, 20, 30, 40];
        b.iter(|| {
            rec.tags_editor().update_array_i32(b"XY", black_box(&replacement));
        });
    });
}

// ============================================================================
// mut::* — fixed-length tag mutations
// ============================================================================

fn bench_tags_mut(c: &mut Criterion) {
    c.bench_function("mut::set_array_element_i32", |b| {
        let aux = bi_array_aux(&[1, 2, 3, 4]);
        let bytes = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &aux);
        let mut rec = RawRecord::from(bytes);
        b.iter(|| {
            rec.tags_mut().set_array_element_i32(b"XY", 2, black_box(99));
        });
    });
}

// ============================================================================
// length::* — length-changing record setters
// ============================================================================

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

    c.bench_function("length::set_read_name_grow_then_shrink", |b| {
        b.iter(|| {
            let mut rec = baseline();
            rec.set_read_name(black_box(b"a_completely_different_longer_name_xxx"));
            rec.set_read_name(black_box(b"x"));
        });
    });

    c.bench_function("length::set_cigar_ops_replace_with_3_ops", |b| {
        b.iter(|| {
            let mut rec = baseline();
            rec.set_cigar_ops(black_box(&[encode_op(0, 50), encode_op(2, 1), encode_op(0, 100)]));
        });
    });

    c.bench_function("length::set_sequence_and_qualities_grow", |b| {
        let seq = vec![b'A'; 200];
        let qual = vec![30u8; 200];
        b.iter(|| {
            let mut rec = baseline();
            rec.set_sequence_and_qualities(black_box(&seq), black_box(&qual));
        });
    });

    c.bench_function("length::set_sequence_and_qualities_same_length", |b| {
        let seq = vec![b'C'; 150];
        let qual = vec![25u8; 150];
        b.iter(|| {
            let mut rec = baseline();
            rec.set_sequence_and_qualities(black_box(&seq), black_box(&qual));
        });
    });
}

// ============================================================================
// seq::* — per-base / per-qual writes
// ============================================================================

fn bench_per_base(c: &mut Criterion) {
    let seq_bytes = baseline_bytes_mapped();

    c.bench_function("seq::set_base", |b| {
        let mut rec = RawRecord::from(seq_bytes.clone());
        b.iter(|| {
            rec.set_base(black_box(75), black_box(b'C'));
        });
    });

    c.bench_function("seq::set_qual", |b| {
        let mut rec = RawRecord::from(seq_bytes.clone());
        b.iter(|| {
            rec.set_qual(black_box(75), black_box(40));
        });
    });
}

// ============================================================================
// vs_recordbuf: side-by-side raw-byte vs noodles `RecordBuf`
// ============================================================================

fn bench_flags_read(c: &mut Criterion) {
    let mut group = c.benchmark_group("flags_read");
    let bytes = baseline_bytes_unmapped();

    group.bench_function("raw_bytes", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes));
            black_box(v.flags())
        });
    });
    group.bench_function("recordbuf", |b| {
        b.iter(|| {
            let buf = decode_one(black_box(&bytes));
            black_box(buf.flags().bits())
        });
    });
    group.finish();
}

fn bench_update_int_tag(c: &mut Criterion) {
    let mut group = c.benchmark_group("update_int_tag");
    let bytes = baseline_bytes_unmapped();

    group.bench_function("raw_bytes", |b| {
        b.iter(|| {
            let mut rec = RawRecord::from(black_box(bytes.clone()));
            rec.tags_editor().update_int(b"NM", 7);
            black_box(rec)
        });
    });
    group.bench_function("recordbuf", |b| {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        let tag = Tag::from(*b"NM");
        b.iter(|| {
            let mut buf = decode_one(black_box(&bytes));
            buf.data_mut().insert(tag, Value::Int32(7));
            black_box(encode_one(&buf))
        });
    });
    group.finish();
}

fn bench_set_read_name_vs(c: &mut Criterion) {
    let mut group = c.benchmark_group("set_read_name");
    let bytes = baseline_bytes_unmapped();
    let new_name = b"a_completely_different_longer_name_for_this_read";

    group.bench_function("raw_bytes", |b| {
        b.iter(|| {
            let mut rec = RawRecord::from(black_box(bytes.clone()));
            rec.set_read_name(black_box(new_name));
            black_box(rec)
        });
    });
    group.bench_function("recordbuf", |b| {
        b.iter(|| {
            let mut buf = decode_one(black_box(&bytes));
            *buf.name_mut() = Some(bstr::BString::from(new_name.as_ref()));
            black_box(encode_one(&buf))
        });
    });
    group.finish();
}

fn bench_set_cigar_ops_vs(c: &mut Criterion) {
    let mut group = c.benchmark_group("set_cigar_ops");
    let bytes = baseline_bytes_unmapped();
    let raw_ops = [encode_op(0, 50), encode_op(2, 1), encode_op(0, 100)];
    let noodles_ops: CigarBuf =
        [Op::new(Kind::Match, 50), Op::new(Kind::Deletion, 1), Op::new(Kind::Match, 100)]
            .into_iter()
            .collect();

    group.bench_function("raw_bytes", |b| {
        b.iter(|| {
            let mut rec = RawRecord::from(black_box(bytes.clone()));
            rec.set_cigar_ops(black_box(&raw_ops));
            black_box(rec)
        });
    });
    group.bench_function("recordbuf", |b| {
        b.iter(|| {
            let mut buf = decode_one(black_box(&bytes));
            *buf.cigar_mut() = noodles_ops.clone();
            black_box(encode_one(&buf))
        });
    });
    group.finish();
}

fn bench_set_sequence_and_qualities_vs(c: &mut Criterion) {
    let mut group = c.benchmark_group("set_sequence_and_qualities");
    let bytes = baseline_bytes_unmapped();
    let seq = vec![b'A'; 200];
    let qual = vec![30u8; 200];
    let noodles_seq = noodles::sam::alignment::record_buf::Sequence::from(seq.clone());
    let noodles_qual = noodles::sam::alignment::record_buf::QualityScores::from(qual.clone());
    // CIGAR must match new seq length or noodles rejects re-encoding.
    let noodles_cigar_200m: CigarBuf = [Op::new(Kind::Match, 200)].into_iter().collect();

    group.bench_function("raw_bytes", |b| {
        b.iter(|| {
            let mut rec = RawRecord::from(black_box(bytes.clone()));
            rec.set_sequence_and_qualities(black_box(&seq), black_box(&qual));
            black_box(rec)
        });
    });
    group.bench_function("recordbuf", |b| {
        b.iter(|| {
            let mut buf = decode_one(black_box(&bytes));
            *buf.sequence_mut() = noodles_seq.clone();
            *buf.quality_scores_mut() = noodles_qual.clone();
            *buf.cigar_mut() = noodles_cigar_200m.clone();
            black_box(encode_one(&buf))
        });
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_view_accessors,
    bench_tags_read,
    bench_tags_editor_scalar,
    bench_tags_editor_variable,
    bench_tags_mut,
    bench_length_changing,
    bench_per_base,
    bench_flags_read,
    bench_update_int_tag,
    bench_set_read_name_vs,
    bench_set_cigar_ops_vs,
    bench_set_sequence_and_qualities_vs,
);
criterion_main!(benches);
