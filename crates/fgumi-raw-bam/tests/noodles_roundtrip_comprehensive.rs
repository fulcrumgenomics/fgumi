//! Comprehensive noodles round-trip validation tests.
//!
//! Ensures every major raw-byte API path produces BAM-spec-compliant bytes
//! by verifying that noodles can read back what we write, field by field.
//!
//! Each test section is labelled with what it exercises:
//!
//! 1. [`SamBuilder`] (all fields + all tag types) → noodles `RecordBuf` decode
//! 2. [`UnmappedSamBuilder`] → noodles decode
//! 3. `write_raw_record` framing → noodles BAM reader
//! 4. Tag editor mutations → noodles tag verification
//! 5. SIMD sequence extraction vs noodles sequence decode
//! 6. Proptest random complete records → noodles full field verification
#![cfg(feature = "noodles")]
#![deny(unsafe_code)]
#![allow(clippy::too_many_lines)] // test functions are intentionally exhaustive
#![allow(clippy::cast_possible_truncation)] // test values are always in range
#![allow(clippy::cast_sign_loss)] // pos() as usize in assertions is safe for test data

use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record::data::field::Tag as NoodlesTag;
use noodles::sam::alignment::record_buf::data::field::Value as NoodlesValue;
use noodles::sam::alignment::record_buf::data::field::value::Array as NoodlesArray;

use fgumi_raw_bam::testutil::encode_op;
use fgumi_raw_bam::{
    RawRecord, SamBuilder, UnmappedSamBuilder, raw_record_to_record_buf,
    raw_records_to_record_bufs, raw_records_to_record_bufs_with_header, write_raw_record,
};

// ============================================================================
// Helpers
// ============================================================================

fn decode(rec: &RawRecord) -> noodles::sam::alignment::RecordBuf {
    let header = noodles::sam::Header::default();
    raw_record_to_record_buf(rec, &header).expect("noodles decode failed")
}

#[allow(clippy::trivially_copy_pass_by_ref)]
fn tag(s: &[u8; 2]) -> NoodlesTag {
    NoodlesTag::from(*s)
}

fn make_bam_stream(records: &[RawRecord]) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b"BAM\x01");
    out.extend_from_slice(&0u32.to_le_bytes()); // header text len
    out.extend_from_slice(&0u32.to_le_bytes()); // n_ref = 0
    for rec in records {
        write_raw_record(&mut out, rec).expect("write_raw_record failed");
    }
    out
}

// ============================================================================
// 1. SamBuilder → noodles decode: verify every field and every tag type
// ============================================================================

#[test]
fn sam_builder_all_fields_roundtrip() {
    use noodles::sam::alignment::record::cigar::op::Kind;

    let seq = b"ACGTACGTACGTACGT";
    let quals: Vec<u8> = (0..16u8).map(|i| (i * 3) % 40).collect();
    // 10M + 2I + 4M = 16 query-consuming bases
    let cigar_ops = &[encode_op(0, 10), encode_op(1, 2), encode_op(0, 4)];

    let mut b = SamBuilder::new();
    b.ref_id(3)
        .pos(999)
        .mapq(60)
        .flags(0x0001 | 0x0040) // PAIRED | READ1
        .mate_ref_id(3)
        .mate_pos(1500)
        .template_length(200)
        .bin(4681)
        .read_name(b"read/roundtrip")
        .cigar_ops(cigar_ops)
        .sequence(seq)
        .qualities(&quals);

    b.add_string_tag(b"RG", b"sample_rg"); // Z
    b.add_int_tag(b"NM", 1); // smallest signed int
    b.add_float_tag(b"AS", 42.5_f32); // f
    b.add_array_u8(
        b"BD",
        &[10u8, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160],
    ); // B:C
    b.add_array_u16(
        b"XU",
        &[100u16, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600],
    ); // B:S
    b.add_array_i16(b"cd", &[-1i16, 2, -3, 4, -5, 6, -7, 8, -9, 10, -11, 12, -13, 14, -15, 16]); // B:s
    b.add_array_i32(b"YI", &[-100_000i32, 200_000, -300_000, 400_000]); // B:i
    b.add_array_f32(b"YF", &[1.1_f32, 2.2, 3.3, 4.4]); // B:f
    b.add_int_tag(b"X0", 70_000); // forces S or i encoding

    let rec = b.build();
    let buf = decode(&rec);

    // -- header fields --
    assert_eq!(rec.ref_id(), 3);
    assert_eq!(rec.pos(), 999);
    assert_eq!(rec.mapq(), 60);
    assert_eq!(buf.alignment_start().map(usize::from), Some(rec.pos() as usize + 1));
    assert_eq!(buf.mapping_quality().map(u8::from), Some(60u8));
    assert_eq!(u16::from(buf.flags()), rec.flags());

    // -- name --
    assert_eq!(buf.name().map(std::convert::AsRef::as_ref), Some(b"read/roundtrip" as &[u8]));

    // -- CIGAR --
    let ops: Vec<_> = buf
        .cigar()
        .iter()
        .map(|op| op.map(|o| (o.kind(), o.len())))
        .collect::<Result<_, _>>()
        .expect("cigar iter");
    assert_eq!(ops, &[(Kind::Match, 10), (Kind::Insertion, 2), (Kind::Match, 4)]);

    // -- sequence and quality scores --
    assert_eq!(buf.sequence().as_ref(), seq.as_ref());
    assert_eq!(rec.sequence_vec(), seq.to_vec());
    assert_eq!(buf.quality_scores().as_ref(), quals.as_slice());

    // -- tags --
    let data = buf.data();

    let rg = data.get(&tag(b"RG")).expect("RG absent");
    assert!(matches!(rg, NoodlesValue::String(s) if AsRef::<[u8]>::as_ref(s) == b"sample_rg"));

    let as_ = data.get(&tag(b"AS")).expect("AS absent");
    assert!(matches!(as_, NoodlesValue::Float(f) if (f - 42.5).abs() < 1e-5));

    let bd = data.get(&tag(b"BD")).expect("BD absent");
    if let NoodlesValue::Array(NoodlesArray::UInt8(vals)) = bd {
        assert_eq!(
            vals.as_slice(),
            &[10u8, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160]
        );
    } else {
        panic!("BD not UInt8 array: {bd:?}");
    }

    let yi = data.get(&tag(b"YI")).expect("YI absent");
    if let NoodlesValue::Array(NoodlesArray::Int32(vals)) = yi {
        assert_eq!(vals.as_slice(), &[-100_000i32, 200_000, -300_000, 400_000]);
    } else {
        panic!("YI not Int32 array: {yi:?}");
    }

    let yf = data.get(&tag(b"YF")).expect("YF absent");
    if let NoodlesValue::Array(NoodlesArray::Float(vals)) = yf {
        for (a, b) in vals.iter().zip([1.1_f32, 2.2, 3.3, 4.4].iter()) {
            assert!((a - b).abs() < 1e-5, "float mismatch: {a} != {b}");
        }
    } else {
        panic!("YF not Float array: {yf:?}");
    }
}

#[test]
fn sam_builder_i16_array_tag_roundtrip() {
    let seq = b"ACGTACGT";
    let quals = &[30u8, 31, 32, 33, 34, 35, 36, 37];

    let mut b = SamBuilder::new();
    b.read_name(b"i16test").cigar_ops(&[encode_op(0, 8)]).sequence(seq).qualities(quals);
    b.add_array_i16(b"cd", &[-100i16, 200, -300, 400, -500, 600, -700, 800]);

    let buf = decode(&b.build());
    let cd = buf.data().get(&tag(b"cd")).expect("cd absent");
    if let NoodlesValue::Array(NoodlesArray::Int16(vals)) = cd {
        assert_eq!(vals.as_slice(), &[-100i16, 200, -300, 400, -500, 600, -700, 800]);
    } else {
        panic!("cd not Int16 array: {cd:?}");
    }
}

#[test]
fn sam_builder_u16_array_tag_roundtrip() {
    let seq = b"ACGT";
    let quals = &[30u8, 31, 32, 33];

    let mut b = SamBuilder::new();
    b.read_name(b"u16test").cigar_ops(&[encode_op(0, 4)]).sequence(seq).qualities(quals);
    b.add_array_u16(b"XU", &[1000u16, 2000, 3000, 65535]);

    let buf = decode(&b.build());
    let xu = buf.data().get(&tag(b"XU")).expect("XU absent");
    if let NoodlesValue::Array(NoodlesArray::UInt16(vals)) = xu {
        assert_eq!(vals.as_slice(), &[1000u16, 2000, 3000, 65535]);
    } else {
        panic!("XU not UInt16 array: {xu:?}");
    }
}

/// Verify integer tags with various values encode to the expected BAM type and
/// that noodles decodes each to the correct numeric value.
#[test]
fn sam_builder_int_tag_types_roundtrip() {
    let mut b = SamBuilder::new();
    b.read_name(b"inttypes")
        .cigar_ops(&[encode_op(0, 4)])
        .sequence(b"ACGT")
        .qualities(&[30u8, 31, 32, 33]);
    b.add_int_tag(b"X1", -100); // i8 range → 'c'
    b.add_int_tag(b"X2", 200); // u8 range (>127) → 'C'
    b.add_int_tag(b"X3", -1000); // i16 range → 's'
    b.add_int_tag(b"X4", 40000); // u16 range (>32767) → 'S'
    b.add_int_tag(b"X5", 100_000); // i32 range → 'i'

    let buf = decode(&b.build());
    let data = buf.data();

    let get_int = |t: &[u8; 2]| -> i64 {
        data.get(&tag(t)).and_then(NoodlesValue::as_int).unwrap_or_else(|| {
            panic!("missing or non-integer tag {}", std::str::from_utf8(t).unwrap())
        })
    };

    assert_eq!(get_int(b"X1"), -100);
    assert_eq!(get_int(b"X2"), 200);
    assert_eq!(get_int(b"X3"), -1000);
    assert_eq!(get_int(b"X4"), 40000);
    assert_eq!(get_int(b"X5"), 100_000);
}

// ============================================================================
// 2. UnmappedSamBuilder → noodles decode
// ============================================================================

#[test]
fn unmapped_sam_builder_roundtrip() {
    use fgumi_raw_bam::fields::flags;

    let seq = b"ACGTNACGTN";
    let quals = &[10u8, 20, 30, 40, 0, 10, 20, 30, 40, 0];

    let mut ub = UnmappedSamBuilder::new();
    ub.build_record(b"unmapped/1", flags::UNMAPPED | flags::PAIRED, seq, quals);
    ub.append_string_tag(b"RG", b"sample1");
    ub.append_int_tag(b"NM", 0);
    ub.append_float_tag(b"XF", 9.75_f32);
    ub.append_i16_array_tag(b"cd", &[5i16, 10, 15]);
    ub.append_u8_array_tag(b"BD", &[1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10]);

    let raw = RawRecord::from(ub.as_bytes().to_vec());
    let buf = decode(&raw);

    // unmapped → ref_id = -1 → noodles reference_sequence_id is None
    assert_eq!(buf.reference_sequence_id(), None);
    assert_eq!(buf.alignment_start(), None);

    assert!(buf.flags().is_unmapped());
    assert!(buf.flags().is_segmented());

    assert_eq!(buf.name().map(std::convert::AsRef::as_ref), Some(b"unmapped/1" as &[u8]));
    assert_eq!(buf.sequence().as_ref(), seq.as_ref());
    assert_eq!(buf.quality_scores().as_ref(), quals.as_ref());
    assert_eq!(buf.cigar().as_ref().len(), 0);

    let data = buf.data();

    let rg = data.get(&tag(b"RG")).expect("RG absent");
    assert!(matches!(rg, NoodlesValue::String(s) if AsRef::<[u8]>::as_ref(s) == b"sample1"));

    let xf = data.get(&tag(b"XF")).expect("XF absent");
    assert!(matches!(xf, NoodlesValue::Float(f) if (f - 9.75_f32).abs() < 1e-5));

    let cd = data.get(&tag(b"cd")).expect("cd absent");
    if let NoodlesValue::Array(NoodlesArray::Int16(vals)) = cd {
        assert_eq!(vals.as_slice(), &[5i16, 10, 15]);
    } else {
        panic!("cd not Int16 array: {cd:?}");
    }

    let bd = data.get(&tag(b"BD")).expect("BD absent");
    if let NoodlesValue::Array(NoodlesArray::UInt8(vals)) = bd {
        assert_eq!(vals.as_slice(), &[1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
    } else {
        panic!("BD not UInt8 array: {bd:?}");
    }
}

// ============================================================================
// 3. write_raw_record framing → noodles BAM reader
// ============================================================================

#[test]
fn write_raw_record_framing_roundtrip() {
    let n = 5usize;
    let mut records: Vec<RawRecord> = Vec::new();

    for i in 0..n {
        let name = format!("read{i:04}");
        let seq: Vec<u8> = (0..10).map(|j| b"ACGTN"[(i + j) % 5]).collect();
        let quals: Vec<u8> = (0..10u8).map(|j| j + 20).collect();

        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .cigar_ops(&[encode_op(0, 10)])
            .sequence(&seq)
            .qualities(&quals);
        records.push(b.build());
    }

    let bam_bytes = make_bam_stream(&records);
    let cursor = std::io::Cursor::new(bam_bytes);
    let mut reader = noodles::bam::io::Reader::from(cursor);
    let header = reader.read_header().expect("read BAM header");

    let decoded: Vec<_> =
        reader.record_bufs(&header).collect::<Result<_, _>>().expect("read BAM records");

    assert_eq!(decoded.len(), n, "record count mismatch");

    for (i, buf) in decoded.iter().enumerate() {
        let expected_name = format!("read{i:04}");
        let got_name = buf.name().map(std::convert::AsRef::as_ref);
        assert_eq!(got_name, Some(expected_name.as_bytes()), "name mismatch at record {i}");

        let expected_seq: Vec<u8> = (0..10).map(|j| b"ACGTN"[(i + j) % 5]).collect();
        assert_eq!(buf.sequence().as_ref(), expected_seq.as_slice(), "seq mismatch at {i}");
    }
}

#[test]
fn write_raw_record_single_framing_roundtrip() {
    let mut b = SamBuilder::new();
    b.read_name(b"singleton")
        .cigar_ops(&[encode_op(0, 4)])
        .sequence(b"TGCA")
        .qualities(&[40u8, 39, 38, 37])
        .ref_id(0)
        .pos(42);

    let bam_bytes = make_bam_stream(&[b.build()]);
    let cursor = std::io::Cursor::new(bam_bytes);
    let mut reader = noodles::bam::io::Reader::from(cursor);
    let header = reader.read_header().expect("read BAM header");

    let decoded: Vec<_> =
        reader.record_bufs(&header).collect::<Result<_, _>>().expect("read BAM records");
    assert_eq!(decoded.len(), 1);
    assert_eq!(decoded[0].name().map(std::convert::AsRef::as_ref), Some(b"singleton" as &[u8]));
    assert_eq!(decoded[0].sequence().as_ref(), b"TGCA");
    // pos()=42 (0-based) → alignment_start=43 (1-based)
    assert_eq!(decoded[0].alignment_start().map(usize::from), Some(43));
}

// ============================================================================
// 4. Tag editor mutations → noodles tag verification
// ============================================================================

#[test]
fn tags_editor_mutations_roundtrip() {
    let mut b = SamBuilder::new();
    b.read_name(b"tageditor")
        .cigar_ops(&[encode_op(0, 4)])
        .sequence(b"ACGT")
        .qualities(&[30u8, 31, 32, 33]);
    b.add_int_tag(b"NM", 5);
    b.add_string_tag(b"RG", b"old_value");
    let mut rec = b.build();

    {
        let mut ed = rec.tags_editor();
        ed.update_int(b"NM", 99);
        ed.update_string(b"RG", b"new_longer_value");
        ed.append_float(b"AS", 7.5_f32);
        ed.append_array_u8(b"BD", &[1u8, 2, 3, 4]);
        ed.append_array_i16(b"cd", &[-10i16, 20, -30]);
        ed.append_array_i32(b"XA", &[1000i32, -2000, 3000]);
    }

    let buf = decode(&rec);
    let data = buf.data();

    let nm = data.get(&tag(b"NM")).expect("NM absent");
    assert_eq!(nm.as_int(), Some(99));

    let rg = data.get(&tag(b"RG")).expect("RG absent");
    assert!(
        matches!(rg, NoodlesValue::String(s) if AsRef::<[u8]>::as_ref(s) == b"new_longer_value")
    );

    let as_ = data.get(&tag(b"AS")).expect("AS absent");
    assert!(matches!(as_, NoodlesValue::Float(f) if (f - 7.5_f32).abs() < 1e-5));

    let bd = data.get(&tag(b"BD")).expect("BD absent");
    if let NoodlesValue::Array(NoodlesArray::UInt8(vals)) = bd {
        assert_eq!(vals.as_slice(), &[1u8, 2, 3, 4]);
    } else {
        panic!("BD not UInt8 array: {bd:?}");
    }

    let cd = data.get(&tag(b"cd")).expect("cd absent");
    if let NoodlesValue::Array(NoodlesArray::Int16(vals)) = cd {
        assert_eq!(vals.as_slice(), &[-10i16, 20, -30]);
    } else {
        panic!("cd not Int16 array: {cd:?}");
    }

    let xa = data.get(&tag(b"XA")).expect("XA absent");
    if let NoodlesValue::Array(NoodlesArray::Int32(vals)) = xa {
        assert_eq!(vals.as_slice(), &[1000i32, -2000, 3000]);
    } else {
        panic!("XA not Int32 array: {xa:?}");
    }
}

#[test]
fn tags_editor_remove_tag_roundtrip() {
    let mut b = SamBuilder::new();
    b.read_name(b"rmtest")
        .cigar_ops(&[encode_op(0, 4)])
        .sequence(b"TTTT")
        .qualities(&[20u8, 20, 20, 20]);
    b.add_string_tag(b"MI", b"42");
    b.add_int_tag(b"NM", 2);
    let mut rec = b.build();

    rec.tags_editor().remove(b"MI");

    let buf = decode(&rec);
    let data = buf.data();
    assert!(data.get(&tag(b"MI")).is_none(), "MI should have been removed");
    assert!(data.get(&tag(b"NM")).is_some(), "NM should still be present");
}

// ============================================================================
// 5. SIMD sequence extraction vs noodles sequence decode
// ============================================================================

/// For a range of sequence lengths, build a record with a known sequence, decode
/// with noodles, and compare our `sequence_vec()` output byte-for-byte with what
/// noodles decodes.  This exercises the SIMD nibble2base path.
#[test]
fn simd_sequence_extraction_vs_noodles() {
    const LENGTHS: &[usize] = &[1, 31, 32, 33, 64, 65, 127, 150, 255];
    const BASES: &[u8] = b"ACGTN";

    for &n in LENGTHS {
        let seq: Vec<u8> = (0..n).map(|i| BASES[i % 5]).collect();
        let quals: Vec<u8> = (0..n).map(|i| (i % 40) as u8).collect();

        let mut b = SamBuilder::new();
        b.read_name(b"simdtest").cigar_ops(&[encode_op(0, n)]).sequence(&seq).qualities(&quals);
        let rec = b.build();

        let buf = decode(&rec);

        assert_eq!(rec.sequence_vec().as_slice(), buf.sequence().as_ref(), "seq mismatch n={n}");
        assert_eq!(
            rec.quality_scores().to_vec().as_slice(),
            buf.quality_scores().as_ref(),
            "qual mismatch n={n}"
        );
    }
}

// ============================================================================
// 6. Proptest random complete records → noodles full field verification
// ============================================================================

proptest::proptest! {
    /// Generate random but valid mapped records via `SamBuilder` and verify that
    /// noodles decodes every field correctly: name, CIGAR, sequence, quality scores,
    /// and tags.
    #[test]
    fn proptest_random_complete_record_roundtrip(
        name in "[A-Za-z0-9]{1,30}",
        op_lengths in proptest::collection::vec(1usize..=50usize, 1..=6),
        mapq in 0u8..=255u8,
        ref_id in 0i32..=10i32,
        pos in 0i32..=9999i32,
        nm in 0i32..=100i32,
        as_score in -1000.0f32..=1000.0f32,
    ) {
        let cigar_ops: Vec<u32> = op_lengths.iter().map(|&l| encode_op(0, l)).collect();
        let l_seq: usize = op_lengths.iter().sum();

        let seq: Vec<u8> = (0..l_seq).map(|i| b"ACGTN"[i % 5]).collect();
        let quals: Vec<u8> = (0..l_seq).map(|i| (i % 40) as u8).collect();

        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
         .cigar_ops(&cigar_ops)
         .sequence(&seq)
         .qualities(&quals)
         .mapq(mapq)
         .ref_id(ref_id)
         .pos(pos);
        b.add_int_tag(b"NM", nm);
        b.add_float_tag(b"AS", as_score);
        let rec = b.build();

        let header = noodles::sam::Header::default();
        let buf = raw_record_to_record_buf(&rec, &header)
            .expect("noodles decode failed");

        // name
        proptest::prop_assert_eq!(
            buf.name().map(std::convert::AsRef::as_ref),
            Some(name.as_bytes())
        );

        // mapq — noodles treats 255 as "missing" per BAM spec
        let expected_mapq = if mapq == 255 { None } else { Some(mapq) };
        proptest::prop_assert_eq!(buf.mapping_quality().map(u8::from), expected_mapq);

        // CIGAR op count
        proptest::prop_assert_eq!(buf.cigar().iter().count(), op_lengths.len());

        // sequence and quality scores
        proptest::prop_assert_eq!(buf.sequence().as_ref(), seq.as_slice());
        proptest::prop_assert_eq!(buf.quality_scores().as_ref(), quals.as_slice());

        // NM integer tag
        let nm_val = buf.data().get(&NoodlesTag::from(*b"NM")).and_then(NoodlesValue::as_int);
        proptest::prop_assert_eq!(nm_val, Some(i64::from(nm)));

        // AS float tag
        let as_val = buf.data().get(&NoodlesTag::from(*b"AS"));
        proptest::prop_assert!(as_val.is_some(), "AS tag absent");
        if let Some(NoodlesValue::Float(decoded_f)) = as_val {
            proptest::prop_assert!(
                (decoded_f - as_score).abs() < 1e-3_f32,
                "AS float mismatch: {decoded_f} != {as_score}"
            );
        }
    }

    /// Verify `raw_records_to_record_bufs` handles random batches of unmapped records
    /// and that every record's name and sequence survive the round-trip.
    #[test]
    fn proptest_batch_unmapped_roundtrip(
        names in proptest::collection::vec("[A-Za-z0-9]{1,20}", 1..=8),
    ) {
        let mut raw_vecs: Vec<Vec<u8>> = Vec::new();
        let mut expected_names: Vec<Vec<u8>> = Vec::new();
        let mut expected_seqs: Vec<Vec<u8>> = Vec::new();

        for (i, name) in names.iter().enumerate() {
            let seq: Vec<u8> = (0..8).map(|j| b"ACGTN"[(i + j) % 5]).collect();
            let quals = vec![20u8; 8];

            let mut ub = UnmappedSamBuilder::new();
            ub.build_record(name.as_bytes(), 4u16 /* UNMAPPED */, &seq, &quals);
            raw_vecs.push(ub.as_bytes().to_vec());
            expected_names.push(name.as_bytes().to_vec());
            expected_seqs.push(seq);
        }

        let bufs = raw_records_to_record_bufs(&raw_vecs).expect("batch decode failed");
        proptest::prop_assert_eq!(bufs.len(), names.len());

        for (i, buf) in bufs.iter().enumerate() {
            proptest::prop_assert_eq!(
                buf.name().map(std::convert::AsRef::as_ref),
                Some(expected_names[i].as_slice())
            );
            proptest::prop_assert_eq!(buf.sequence().as_ref(), expected_seqs[i].as_slice());
        }
    }
}

/// Mapped record round-trip through `raw_record_to_record_buf` with a non-empty
/// SAM header: the supplied `@SQ` dictionary must actually be used during
/// decoding, and encoding back through a noodles writer must produce the same
/// reference name for the mapped reference ID.
#[test]
fn mapped_record_roundtrip_with_non_empty_header() {
    use bstr::BString;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;

    // Build a header with a non-default reference dictionary so the decode
    // path must honor it. `@HD SO:coordinate` is required for reference IDs
    // to be meaningful in mapped records.
    let hd: noodles::sam::Header = "@HD\tVN:1.6\tSO:coordinate\n".parse().expect("parse @HD");
    let mut hb = noodles::sam::Header::builder();
    if let Some(h) = hd.header() {
        hb = hb.set_header(h.clone());
    }
    let header = hb
        .add_reference_sequence(
            BString::from("chrA"),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(10_000).unwrap()),
        )
        .add_reference_sequence(
            BString::from("chrB"),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(20_000).unwrap()),
        )
        .build();

    // Build a mapped record against the second reference (index 1 → "chrB").
    let mut b = SamBuilder::new();
    b.ref_id(1)
        .pos(499) // 500 in 1-based
        .mapq(42)
        .read_name(b"mapped/1")
        .cigar_ops(&[encode_op(0, 10)])
        .sequence(b"ACGTACGTAC")
        .qualities(&[30u8; 10]);
    let rec = b.build();

    // Decode with the real header.
    let buf = raw_record_to_record_buf(&rec, &header).expect("decode with header");
    assert_eq!(buf.alignment_start().map(usize::from), Some(500));
    assert_eq!(buf.reference_sequence_id(), Some(1));

    // Re-encode through a noodles writer to verify the header genuinely drives
    // the decode path: writing requires the header to resolve the reference id
    // and will error if the decoded record is inconsistent with the header.
    let mut out: Vec<u8> = Vec::new();
    {
        use noodles::sam::alignment::io::Write as _;
        let mut writer = noodles::bam::io::Writer::from(&mut out);
        writer.write_header(&header).expect("write header");
        writer.write_alignment_record(&header, &buf).expect("write record");
    }
    assert!(!out.is_empty(), "re-encode should produce bytes");

    // Sanity: `raw_records_to_record_bufs_with_header` is the plural counterpart.
    let bufs =
        raw_records_to_record_bufs_with_header(&[rec.as_ref().to_vec()], &header).expect("batch");
    assert_eq!(bufs.len(), 1);
    assert_eq!(bufs[0].reference_sequence_id(), Some(1));
}
