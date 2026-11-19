//! Writer options and record shapes beyond the basic round trip: long CIGARs
//! (`CG:B,I`), reference-free output, the bzip2/lzma codecs, and the file ID.
mod common;

use common::{appended_md_nm, cram_layout, read_cram, record, reference, write_cram, write_fasta};
use fgumi_htslib::CramWriteOptions;
use fgumi_tag::SamTag;
use rstest::rstest;

const REF_LEN: usize = 60_000;
const HEADER: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrL\tLN:60000\n";

/// Require an unsigned `CG:B,I` array and return its payload.
fn cg_tag(rec: &[u8]) -> &[u8] {
    let tag = fgumi_raw_bam::find_array_tag(fgumi_raw_bam::aux_data_slice(rec), SamTag::CG)
        .expect("long CIGAR requires a CG array");
    assert_eq!(tag.elem_type, b'I', "CG array must use unsigned 32-bit elements");
    tag.data
}

#[test]
fn long_cigar_round_trips_through_cg_tag() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference(REF_LEN, 0xD1B5_4A32_D192_ED03);
    let fa = write_fasta(dir.path(), "chrL", &seq);
    // 70,000 ops alternating 1M/1I over 35,000 reference bases.
    let mut cigar = Vec::new();
    let mut query = Vec::new();
    for i in 0..35_000usize {
        cigar.push(1 << 4); // 1M
        cigar.push(1 << 4 | 1); // 1I
        query.push(seq[100 + i]);
        query.push(b'T');
    }
    let rec = record("long", 0, 0, 100, &cigar, &query, b"");
    let out = dir.path().join("long.cram");
    write_cram(
        &out,
        HEADER,
        &CramWriteOptions { reference: Some(&fa), ..Default::default() },
        std::slice::from_ref(&rec),
    );

    let (got, _) = read_cram(&out, Some(&fa)).unwrap();
    assert_eq!(got.len(), 1);
    // Same fixed fields, read name, placeholder CIGAR, sequence and qualities.
    let body = 32 + 5 + 8 + query.len().div_ceil(2) + query.len();
    assert_eq!(
        &got[0][..body],
        &rec[..body],
        "core fields, placeholder CIGAR, SEQ or QUAL changed"
    );
    assert_eq!(cg_tag(&got[0]), cg_tag(&rec), "the real CIGAR (CG:B,I) did not round trip");
}

fn mapped_records(seq: &[u8]) -> Vec<Vec<u8>> {
    let mut recs: Vec<Vec<u8>> = (0..3000usize)
        .map(|i| {
            let pos = (i * 17) % (REF_LEN - 200);
            record(
                &format!("r{i}"),
                0,
                0,
                i32::try_from(pos).unwrap(),
                &[100 << 4],
                &seq[pos..pos + 100],
                b"",
            )
        })
        .collect();
    recs.sort_by_key(|r| fgumi_raw_bam::pos(r));
    recs
}

/// A reference-free CRAM decodes with no reference at all, and (with nothing
/// to regenerate MD/NM from) byte-identically.
#[test]
fn reference_free_round_trips_exactly_without_a_reference() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference(REF_LEN, 0xD1B5_4A32_D192_ED03);
    let recs = mapped_records(&seq);
    let out = dir.path().join("o.cram");
    let opts =
        CramWriteOptions { no_ref: true, file_id: Some("final-name.cram"), ..Default::default() };
    write_cram(&out, HEADER, &opts, &recs);
    let bytes = std::fs::read(&out).unwrap();
    assert_eq!(&bytes[6..26], b"final-name.cram\0\0\0\0\0", "file-definition ID");
    let (got, _) = read_cram(&out, None).unwrap();
    assert_eq!(got, recs);
}

/// The requested codec appears in emitted blocks, and both paths preserve the
/// fixture records in order (apart from regenerated MD/NM tags, which must
/// hold the exact-match values `MD:Z:100` and `NM:0`).
#[rstest]
#[case::bzip2(true, false)]
#[case::lzma(false, true)]
fn optional_codecs_take_effect_and_round_trip(#[case] use_bzip2: bool, #[case] use_lzma: bool) {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference(REF_LEN, 0xD1B5_4A32_D192_ED03);
    let fa = write_fasta(dir.path(), "chrL", &seq);
    let recs = mapped_records(&seq);
    let plain = dir.path().join("plain.cram");
    let codec = dir.path().join("codec.cram");
    let base = CramWriteOptions { reference: Some(&fa), file_id: Some("x"), ..Default::default() };
    write_cram(&plain, HEADER, &base, &recs);
    write_cram(&codec, HEADER, &CramWriteOptions { use_bzip2, use_lzma, ..base }, &recs);
    let method = if use_bzip2 { 2 } else { 3 };
    let uses_method = |path: &std::path::Path| {
        cram_layout(&std::fs::read(path).unwrap())
            .iter()
            .flat_map(|container| &container.blocks)
            .any(|block| block.method == method)
    };
    assert!(!uses_method(&plain), "the baseline must not use the optional codec");
    assert!(uses_method(&codec), "the requested codec must appear in an emitted block");
    let (a, _) = read_cram(&plain, Some(&fa)).unwrap();
    let (b, _) = read_cram(&codec, Some(&fa)).unwrap();
    assert_eq!(a, b);
    for decoded in [&a, &b] {
        assert_eq!(decoded.len(), recs.len());
        for (actual, expected) in decoded.iter().zip(&recs) {
            // Every fixture read is a 100-base exact match.
            assert_eq!(
                appended_md_nm(actual, expected),
                ("100".to_string(), 0),
                "only regenerated MD/NM may differ from the input, with their true values"
            );
        }
    }
}
