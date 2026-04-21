#![deny(unsafe_code)]

//! Side-by-side benchmarks: raw-byte API vs noodles `RecordBuf` decode-mutate-reencode.
//!
//! Each `BenchmarkGroup` contains two members:
//!   - `raw_bytes` — our new zero-copy / in-place API
//!   - `recordbuf` — noodles `RecordBuf` round-trip (decode → mutate → re-encode to bytes)
//!
//! Re-encoding uses `noodles::bam::io::Writer` against an in-memory `Vec<u8>` so both
//! paths start and end at raw bytes, giving an apples-to-apples comparison.

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use noodles::sam::alignment::record::cigar::op::{Kind, Op};
use noodles::sam::alignment::record_buf::Cigar as CigarBuf;

use fgumi_raw_bam::RawRecord;
use fgumi_raw_bam::RawRecordView;
use fgumi_raw_bam::raw_records_to_record_bufs;
use fgumi_raw_bam::testutil::{encode_op, make_bam_bytes};

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

/// Build the canonical benchmark record bytes.
///
/// 150 bp unmapped read (tid=-1, pos=-1), one 150M CIGAR op, NM:i:5 aux tag.
/// Using unmapped coordinates so the record can be re-encoded by noodles using
/// an empty `Header` (no reference sequences required).
fn baseline_bytes() -> Vec<u8> {
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

/// Decode raw bytes to a noodles `RecordBuf`.
fn decode_one(bytes: &[u8]) -> noodles::sam::alignment::RecordBuf {
    let mut bufs = raw_records_to_record_bufs(&[bytes.to_vec()]).expect("decode should succeed");
    bufs.pop().expect("should have one record")
}

/// Re-encode a mutated `RecordBuf` back to raw BAM record bytes (without the
/// 4-byte `block_size` prefix), matching production call sites (see
/// `crates/fgumi-sam/src/record_utils.rs` and `crates/fgumi-sam/src/clipper.rs`).
fn encode_one(buf: &noodles::sam::alignment::RecordBuf) -> Vec<u8> {
    let header = noodles::sam::Header::default();
    fgumi_raw_bam::encode_record_buf_to_raw(buf, &header)
        .expect("encode_record_buf_to_raw should succeed")
        .into_inner()
}

// ---------------------------------------------------------------------------
// Benchmark groups
// ---------------------------------------------------------------------------

/// 1. Header field read: `flags()`
///
/// Both paths start from raw bytes.  The raw path does two array reads
/// (construct view + 2-byte read).  The `RecordBuf` path decodes the full record.
fn bench_flags_read(c: &mut Criterion) {
    let mut group = c.benchmark_group("flags_read");
    let bytes = baseline_bytes();

    group.bench_function("raw_bytes", |b| {
        b.iter(|| {
            let v = RawRecordView::new(black_box(&bytes));
            black_box(v.flags())
        });
    });

    group.bench_function("recordbuf", |b| {
        b.iter(|| {
            let buf = decode_one(black_box(&bytes));
            // Read the flags and convert to bits so the value is actually used.
            black_box(buf.flags().bits())
        });
    });

    group.finish();
}

/// 2. Tag update (in-place, same size): `NM:i:7`
///
/// Both paths start from raw bytes and end with a mutated record in memory.
/// The raw path updates the integer tag in-place without any allocation.
/// The `RecordBuf` path decodes, inserts into the tag map, then re-encodes.
fn bench_update_int_tag(c: &mut Criterion) {
    let mut group = c.benchmark_group("update_int_tag");
    let bytes = baseline_bytes();

    group.bench_function("raw_bytes", |b| {
        b.iter(|| {
            let mut rec = RawRecord::from(black_box(bytes.clone()));
            rec.tags_editor().update_int(b"NM", 7);
            black_box(rec)
        });
    });

    group.bench_function("recordbuf", |b| {
        b.iter(|| {
            let mut buf = decode_one(black_box(&bytes));
            {
                use noodles::sam::alignment::record_buf::data::field::Value;
                let tag = noodles::sam::alignment::record::data::field::Tag::from(*b"NM");
                buf.data_mut().insert(tag, Value::Int32(7));
            }
            black_box(encode_one(&buf))
        });
    });

    group.finish();
}

/// 3. Length-changing edit: `set_read_name`
///
/// Both paths start from a fresh record and replace the read name with a
/// longer string, requiring a memmove of the variable-length data that follows.
fn bench_set_read_name(c: &mut Criterion) {
    let mut group = c.benchmark_group("set_read_name");
    let bytes = baseline_bytes();
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
            // name_mut() returns &mut Option<BString>; BString is bstr's growable byte string.
            *buf.name_mut() = Some(bstr::BString::from(new_name.as_ref()));
            black_box(encode_one(&buf))
        });
    });

    group.finish();
}

/// 4. Length-changing edit: `set_cigar_ops`
///
/// Replace a single 150M CIGAR with three operations (50M 1D 99M).
fn bench_set_cigar_ops(c: &mut Criterion) {
    let mut group = c.benchmark_group("set_cigar_ops");
    let bytes = baseline_bytes();

    // Raw: three raw u32 CIGAR ops — 50M + 1D + 100M = 150 query bases (matches l_seq=150).
    let raw_ops = [encode_op(0, 50), encode_op(2, 1), encode_op(0, 100)];

    // RecordBuf: noodles Op values — same 50M + 1D + 100M.
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

/// 5. Length-changing edit: `set_sequence_and_qualities`
///
/// Grow the sequence from 150 bp to 200 bp with fresh bases and quality scores.
/// Both paths also update the CIGAR to 200M so the final records are
/// equivalent (raw bytes vs. re-encoded `RecordBuf`); otherwise noodles'
/// re-encoder rejects a sequence-length / CIGAR-length mismatch.
fn bench_set_sequence_and_qualities(c: &mut Criterion) {
    let mut group = c.benchmark_group("set_sequence_and_qualities");
    let bytes = baseline_bytes();

    let seq = vec![b'A'; 200];
    let qual = vec![30u8; 200];

    // Raw: single 200M CIGAR op keeps the record self-consistent post-edit.
    let raw_cigar_200m = [encode_op(0, 200)];

    // RecordBuf: Sequence and QualityScores both wrap Vec<u8> and implement From<Vec<u8>>.
    let noodles_seq = noodles::sam::alignment::record_buf::Sequence::from(seq.clone());
    let noodles_qual = noodles::sam::alignment::record_buf::QualityScores::from(qual.clone());
    // Updated CIGAR to match new 200 bp read length so noodles re-encoder does
    // not reject the record.
    let noodles_cigar_200m: CigarBuf = [Op::new(Kind::Match, 200)].into_iter().collect();

    group.bench_function("raw_bytes", |b| {
        b.iter(|| {
            let mut rec = RawRecord::from(black_box(bytes.clone()));
            rec.set_sequence_and_qualities(black_box(&seq), black_box(&qual));
            rec.set_cigar_ops(black_box(&raw_cigar_200m));
            black_box(rec)
        });
    });

    group.bench_function("recordbuf", |b| {
        b.iter(|| {
            let mut buf = decode_one(black_box(&bytes));
            *buf.sequence_mut() = noodles_seq.clone();
            *buf.quality_scores_mut() = noodles_qual.clone();
            // CIGAR must be updated to match new sequence length or noodles
            // will reject re-encoding with "read length-sequence length mismatch".
            *buf.cigar_mut() = noodles_cigar_200m.clone();
            black_box(encode_one(&buf))
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_flags_read,
    bench_update_int_tag,
    bench_set_read_name,
    bench_set_cigar_ops,
    bench_set_sequence_and_qualities,
);
criterion_main!(benches);
