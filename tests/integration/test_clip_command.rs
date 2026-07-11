//! End-to-end tests for the clip command.
//!
//! These tests invoke `Clip::execute()` in-process and validate:
//! 1. Basic read clipping
//! 2. Fixed clipping (5' and 3' ends)
//! 3. Metrics output

use clap::Parser;
use fgumi_lib::commands::clip::Clip;
use fgumi_lib::commands::command::Command;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use rstest::rstest;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_coordinate_sorted_header, create_minimal_header, create_test_reference, to_record_buf,
};

/// Create a BAM with paired reads using the shared minimal (query-grouped) header.
fn create_paired_bam(path: &Path, pairs: Vec<(RawRecord, RawRecord)>) {
    let header = create_minimal_header("chr1", 10000);
    write_paired_bam_with_header(path, &header, pairs);
}

/// Write paired reads under a caller-supplied header (used to exercise
/// header-less and coordinate-sorted input).
fn write_paired_bam_with_header(
    path: &Path,
    header: &noodles::sam::Header,
    pairs: Vec<(RawRecord, RawRecord)>,
) {
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(header).expect("Failed to write header");

    for (r1, r2) in pairs {
        writer.write_alignment_record(header, &to_record_buf(&r1)).expect("Failed to write R1");
        writer.write_alignment_record(header, &to_record_buf(&r2)).expect("Failed to write R2");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// CLIP3-05: clip must reject coordinate-sorted (non-query-grouped) input,
/// matching fgbio's `Bams.requireQueryGrouped`. On coordinate-sorted input mates
/// scatter, so pair clip / overlap / mate-fix silently no-op — hard-fail instead.
///
/// Exercised on both `Clip::execute` paths (the single-threaded fast path and the
/// `--threads` pipeline) since the guard is wired into each.
#[rstest]
#[case::single_threaded(&[])]
#[case::threaded(&["--threads", "2"])]
fn test_clip_rejects_coordinate_sorted_input(#[case] extra_args: &[&str]) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4])
            .mate_ref_id(0)
            .mate_pos(103)
            .template_length(12);
        b.build()
    };
    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(103)
            .mapq(60)
            .cigar_ops(&[8 << 4])
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-12);
        b.build()
    };
    let header = create_coordinate_sorted_header("chr1", 10000);
    write_paired_bam_with_header(&input_bam, &header, vec![(r1, r2)]);

    let mut args = vec![
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--ref",
        ref_path.to_str().unwrap(),
        "--clip-overlapping-reads",
    ];
    args.extend_from_slice(extra_args);
    let cmd = Clip::try_parse_from(args).expect("failed to parse clip args");

    let err = cmd.execute("fgumi clip").expect_err("must reject coordinate-sorted input");
    assert!(
        err.to_string().contains("queryname sorted or query grouped"),
        "unexpected error message: {err}"
    );
}

/// A header carrying only `@SQ` (no `@HD` line) -- mirrors header-less input.
fn headerless_header(ref_name: &str, ref_len: usize) -> noodles::sam::Header {
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReferenceSequence;

    let header = noodles::sam::Header::builder()
        .add_reference_sequence(
            ref_name,
            Map::<ReferenceSequence>::new(
                std::num::NonZero::new(ref_len).expect("non-zero reference length"),
            ),
        )
        .build();
    assert!(header.header().is_none(), "precondition: header must lack @HD");
    header
}

/// Test basic clip command with fixed-end clipping.
#[test]
fn test_clip_command_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Create a paired-end read
    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(103)
            .template_length(12);
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(103)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-12);
        b.build()
    };

    create_paired_bam(&input_bam, vec![(r1, r2)]);

    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--read-one-five-prime",
        "1",
        "--read-one-three-prime",
        "1",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify output has the reads
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert_eq!(count, 2, "Should have both reads in output");
}

/// CLIP3-05: clip must reject header-less input. A header-less BAM synthesizes
/// `@HD VN:1.6 SO:unsorted` (via `ensure_hd_record`), which is neither queryname
/// sorted nor query grouped, so `require_query_grouped` rejects it — matching
/// fgbio's `Bams.requireQueryGrouped`. The `@HD` synthesis itself is still exercised
/// end-to-end by `correct`/`review` (which have no query-grouped guard).
///
/// Exercised on both `Clip::execute` paths (the single-threaded fast path and the
/// `--threads` pipeline): both run `ensure_hd_record` before `require_query_grouped`,
/// so header-less input is rejected the same way regardless of execution mode.
#[rstest]
#[case::single_threaded(&[])]
#[case::threaded(&["--threads", "2"])]
fn test_clip_rejects_headerless_input(#[case] extra_args: &[&str]) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(103)
            .template_length(12);
        b.build()
    };
    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(103)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-12);
        b.build()
    };

    let header = headerless_header("chr1", 10000);
    write_paired_bam_with_header(&input_bam, &header, vec![(r1, r2)]);

    let mut args = vec![
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--read-one-five-prime",
        "1",
        "--read-one-three-prime",
        "1",
        "--compression-level",
        "1",
    ];
    args.extend_from_slice(extra_args);
    let cmd = Clip::try_parse_from(args).expect("failed to parse clip args");

    let err = cmd.execute("fgumi clip").expect_err("must reject header-less input");
    assert!(
        err.to_string().contains("queryname sorted or query grouped"),
        "unexpected error message: {err}"
    );
}

/// Test clip command with metrics output.
#[test]
fn test_clip_command_with_metrics() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");
    let ref_path = create_test_reference(temp_dir.path());

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(199)
            .template_length(108);
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(199)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-108);
        b.build()
    };

    create_paired_bam(&input_bam, vec![(r1, r2)]);

    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--read-one-five-prime",
        "2",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command with metrics failed");
    assert!(metrics_path.exists(), "Metrics file not created");
}

/// Write a BAM from an explicit list of records (any flags), using the shared minimal header.
fn create_bam_from_records(path: &Path, records: &[RawRecord]) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for r in records {
        writer.write_alignment_record(&header, &to_record_buf(r)).expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// One decoded output record: `(name, cigar-ops, raw-flags)`.
type OutputRecord = (String, Vec<(CigarKind, usize)>, u16);

/// Read output records so tests can assert the *actual* clip state (CIGAR changed /
/// unchanged), not merely the record count.
fn read_output_records(path: &Path) -> Vec<OutputRecord> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    reader
        .record_bufs(&header)
        .map(|r| {
            let rec = r.expect("Failed to read output record");
            let name = rec
                .name()
                .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                .unwrap_or_default();
            let ops: Vec<(CigarKind, usize)> =
                rec.cigar().as_ref().iter().map(|op| (op.kind(), op.len())).collect();
            (name, ops, u16::from(rec.flags()))
        })
        .collect()
}

/// Read output records as full noodles `RecordBuf`s so tests can assert mate metadata
/// (mate ref/pos/strand, MC, MQ, TLEN) — not just CIGAR — after clipping.
fn read_output_record_bufs(path: &Path) -> Vec<RecordBuf> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    reader.record_bufs(&header).map(|r| r.expect("Failed to read output record")).collect()
}

/// The record's CIGAR as `(kind, len)` ops (e.g. `[(Match, 98), (HardClip, 2)]`).
fn cigar_ops(rec: &RecordBuf) -> Vec<(CigarKind, usize)> {
    rec.cigar().as_ref().iter().map(|op| (op.kind(), op.len())).collect()
}

/// The record's mate-CIGAR (`MC`) tag as a string, if present.
fn mate_cigar(rec: &RecordBuf) -> Option<String> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    match rec.data().get(&Tag::MATE_CIGAR) {
        Some(Value::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

/// The record's mate-mapping-quality (`MQ`) tag as an integer, if present.
fn mate_mapq(rec: &RecordBuf) -> Option<i64> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    rec.data().get(&Tag::MATE_MAPPING_QUALITY).and_then(Value::as_int)
}

/// Build a mapped read with the given flags/position/CIGAR-length (all `M`).
fn mapped_read(name: &[u8], flags: u16, pos: i32, match_len: u32, mapq: u8) -> RawRecord {
    let len = match_len as usize;
    let mut b = SamBuilder::new();
    b.read_name(name)
        .sequence(&vec![b'A'; len])
        .qualities(&vec![30; len])
        .flags(flags)
        .ref_id(0)
        .pos(pos)
        .mapq(mapq)
        .cigar_ops(&[match_len << 4]) // <len>M
        .mate_ref_id(0)
        .mate_pos(pos);
    b.build()
}

/// Build a mapped read with an explicit CIGAR (raw BAM `(len << 4) | op` codes) and a
/// matching-length sequence/qualities. `read_len` must equal the CIGAR's read-consuming length.
fn read_with_cigar(
    name: &[u8],
    flags: u16,
    pos: i32,
    cigar_ops: &[u32],
    read_len: usize,
    mapq: u8,
) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name)
        .sequence(&vec![b'A'; read_len])
        .qualities(&vec![30; read_len])
        .flags(flags)
        .ref_id(0)
        .pos(pos)
        .mapq(mapq)
        .cigar_ops(cigar_ops)
        .mate_ref_id(0)
        .mate_pos(pos);
    b.build()
}

/// Runs `clip --upgrade-clipping` (Hard mode) over a template whose only clipped read is a
/// supplementary alignment (`10S40M`), and asserts the supplementary's leading soft clip is
/// upgraded to hard. fgbio `ClipBam` upgrades `template.allReads` (`ClipBam.scala:123`), not just
/// the primary R1/R2, so the supplementary must be upgraded too. `threads` selects the
/// single-threaded (`None`) vs `--threads` code path — the fix must hold on both.
fn run_supplementary_upgrade_case(threads: Option<&str>) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Primary pair, both 50M (no clipping to upgrade — they must stay 50M).
    let r1 = mapped_read(b"t", flags::PAIRED | flags::FIRST_SEGMENT, 99, 50, 60);
    let r2 = mapped_read(b"t", flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE, 99, 50, 60);
    // R1 supplementary with a leading soft clip: 10S40M. The only read carrying clipping.
    let supp = read_with_cigar(
        b"t",
        flags::PAIRED | flags::FIRST_SEGMENT | flags::SUPPLEMENTARY,
        199,
        &[(10u32 << 4) | 4, 40u32 << 4], // 10S40M
        50,
        60,
    );
    create_bam_from_records(&input_bam, &[r1, r2, supp]);

    let mut args = vec![
        "clip".to_string(),
        "--input".to_string(),
        input_bam.to_str().unwrap().to_string(),
        "--output".to_string(),
        output_bam.to_str().unwrap().to_string(),
        "--reference".to_string(),
        ref_path.to_str().unwrap().to_string(),
        "--clipping-mode".to_string(),
        "hard".to_string(),
        "--upgrade-clipping".to_string(),
        "true".to_string(),
        "--compression-level".to_string(),
        "1".to_string(),
    ];
    if let Some(t) = threads {
        args.push("--threads".to_string());
        args.push(t.to_string());
    }
    let cmd = Clip::try_parse_from(&args).expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");

    let recs = read_output_records(&output_bam);
    assert_eq!(recs.len(), 3, "all three reads retained");

    // The supplementary read's leading 10S must have been upgraded to 10H.
    let supp_out = recs
        .iter()
        .find(|(_, _, f)| f & flags::SUPPLEMENTARY != 0)
        .expect("supplementary read present in output");
    assert_eq!(
        supp_out.1,
        vec![(CigarKind::HardClip, 10), (CigarKind::Match, 40)],
        "supplementary soft clip must be upgraded to hard (fgbio ClipBam.scala:123 upgrades \
         template.allReads); got {:?}",
        supp_out.1
    );

    // The primaries carried no clipping, so upgrade is a no-op on them.
    for (_, ops, f) in &recs {
        if f & flags::SUPPLEMENTARY == 0 {
            assert_eq!(*ops, vec![(CigarKind::Match, 50)], "primary CIGAR must be unchanged");
        }
    }
}

/// Single-threaded path: `--upgrade-clipping` upgrades a supplementary alignment's clipping.
#[test]
fn test_upgrade_clipping_upgrades_supplementary_single_threaded() {
    run_supplementary_upgrade_case(None);
}

/// `--threads` path: `--upgrade-clipping` upgrades a supplementary alignment's clipping.
#[test]
fn test_upgrade_clipping_upgrades_supplementary_threads_mode() {
    run_supplementary_upgrade_case(Some("2"));
}

/// `--threads` mode routes through `execute_threads_mode`; exercise the `(Some, Some)` primary-pair
/// branch with fixed R1/R2 clipping, overlap clipping, mate-extension clipping and clip upgrading
/// all enabled. R1 (fwd, 150M) and R2 (rev, 100M) fully overlap, so overlap clipping engages.
#[test]
fn test_clip_command_threads_mode_primary_pair_all_options() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let mut r1 =
        mapped_read(b"p", flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE, 99, 150, 60);
    let mut r2 =
        mapped_read(b"p", flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE, 99, 100, 60);
    // A non-zero TLEN is required for `is_fr_pair_raw` to recognize the pair (htsjdk keys FR
    // detection off the insert size); without it overlap and mate-extension clipping would
    // silently no-op and only the fixed 5'/3' clips would apply.
    fgumi_raw_bam::set_template_length(r1.as_mut_vec(), 150);
    fgumi_raw_bam::set_template_length(r2.as_mut_vec(), -150);
    create_bam_from_records(&input_bam, &[r1, r2]);

    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--threads",
        "2",
        "--read-one-five-prime",
        "1",
        "--read-one-three-prime",
        "1",
        "--read-two-five-prime",
        "1",
        "--read-two-three-prime",
        "1",
        "--clip-overlapping-reads",
        "true",
        "--clip-bases-past-mate",
        "true",
        "--upgrade-clipping",
        "true",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");

    // Exact post-clip CIGARs, derived independently from the clip pipeline (Hard mode),
    // not just "some hard clip exists". Inputs are fixed: R1 fwd 150M@100 (1-based),
    // R2 rev 100M@100.
    //
    //   R1: fixed 5'/3' = 1 each  -> 1H148M1H, ref 101..248
    //       overlap clip trims R1's 3' end down to the pair midpoint
    //       (midpoint(r1_start=101, r2_end=198) = 149), clipping ref 150..248 (99 M bases)
    //       -> 1H49M100H, ref 101..149
    //   R2: fixed 5'/3' = 1 each  -> 1H98M1H, ref 101..198
    //       overlap clip trims R2's 5' (low-coord) end up to midpoint+1 = 150,
    //       clipping ref 101..149 (49 M bases) -> 50H49M1H, ref 150..198
    //
    // Mate-extension clipping then finds nothing past the mate (the reads now meet at the
    // midpoint), so it is a no-op here — but its code path still runs. Totals check out:
    // R1 = 1+49+100 = 150 read bases, R2 = 50+49+1 = 100 read bases.
    let recs = read_output_records(&output_bam);
    assert_eq!(recs.len(), 2, "both reads retained");
    let expected_r1 =
        vec![(CigarKind::HardClip, 1), (CigarKind::Match, 49), (CigarKind::HardClip, 100)];
    let expected_r2 =
        vec![(CigarKind::HardClip, 50), (CigarKind::Match, 49), (CigarKind::HardClip, 1)];
    for (name, ops, flag_bits) in &recs {
        assert_eq!(name, "p", "unexpected read name");
        let is_r1 = flag_bits & flags::FIRST_SEGMENT != 0;
        let (label, expected) = if is_r1 { ("R1", &expected_r1) } else { ("R2", &expected_r2) };
        assert_eq!(ops, expected, "{label} exact post-clip CIGAR mismatch; cigar={ops:?}");
    }
}

/// `--threads` mode with a lone unpaired fragment exercises the `(Some, None)` fragment branch.
#[test]
fn test_clip_command_threads_mode_fragment() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // No PAIRED flag: find_primary_pair_indices resolves R1 only.
    let frag = mapped_read(b"f", 0, 99, 100, 60);
    create_bam_from_records(&input_bam, &[frag]);

    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--threads",
        "2",
        "--read-one-five-prime",
        "2",
        "--read-one-three-prime",
        "2",
        "--upgrade-clipping",
        "true",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");

    // The lone forward-strand fragment is hard-clipped 2bp at each end (read-one 5'/3' = 2,
    // Hard mode) => 2H96M2H, proving the (Some, None) branch actually clipped rather than
    // just passing the record through.
    let recs = read_output_records(&output_bam);
    assert_eq!(recs.len(), 1, "fragment retained");
    assert_eq!(
        recs[0].1,
        vec![(CigarKind::HardClip, 2), (CigarKind::Match, 96), (CigarKind::HardClip, 2)],
        "fragment must be hard-clipped 2bp at each end (2H96M2H); cigar={:?}",
        recs[0].1
    );
}

/// `--threads` mode with a template that has no primary read (secondary only) exercises the
/// `(None, None)` no-op branch — the record passes through untouched.
#[test]
fn test_clip_command_threads_mode_secondary_only() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let secondary =
        mapped_read(b"s", flags::PAIRED | flags::FIRST_SEGMENT | flags::SECONDARY, 99, 100, 60);
    create_bam_from_records(&input_bam, &[secondary]);

    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--threads",
        "2",
        "--read-one-five-prime",
        "1",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");

    // The (None, None) branch is a true no-op: the secondary read passes through with its
    // CIGAR unchanged (still 100M, not clipped) and its SECONDARY flag preserved.
    let recs = read_output_records(&output_bam);
    assert_eq!(recs.len(), 1, "secondary read passed through");
    let (_name, ops, flag_bits) = &recs[0];
    assert_eq!(
        ops,
        &vec![(CigarKind::Match, 100)],
        "secondary read must be unchanged (100M no-op); cigar={ops:?}"
    );
    assert_ne!(flag_bits & flags::SECONDARY, 0, "SECONDARY flag must be preserved");
}

/// `--threads` mode with a chimeric/split template — primary R1 + primary R2 + a supplementary
/// R1 — exercises the PR's headline fix end-to-end: after the primary pair is clipped,
/// `fix_supplemental_mate_info` must repair the supplementary alignment's mate metadata to point
/// at the *post-clip* primary R2. The `clip.rs` unit test covers the helper in isolation; this
/// verifies the `execute_threads_mode` wiring (correct r1/r2 indices, and that the supplemental
/// mate snapshot is taken *after* the primary is clipped, not before). A wiring bug (wrong
/// index, or records mutated out of order) would be invisible to every other test in this file.
#[test]
fn test_clip_command_threads_mode_supplementary_mate_repair() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Chimeric template: a primary FR pair (R1 fwd 100M @100, R2 rev 100M @300) plus a
    // supplementary R1 mapped far away (50M @5000). `mapped_read` seeds each record's mate fields
    // to its *own* position with no MC/MQ, so the supplementary's input mate info is deliberately
    // stale (points at 5000, not the primary R2) — only a correct repair yields the expected
    // primary-R2-derived values.
    let r1 = mapped_read(
        b"chim",
        flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE,
        99,
        100,
        60,
    );
    let r2 =
        mapped_read(b"chim", flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE, 299, 100, 40);
    let supp = mapped_read(
        b"chim",
        flags::PAIRED | flags::FIRST_SEGMENT | flags::SUPPLEMENTARY,
        4999,
        50,
        30,
    );
    create_bam_from_records(&input_bam, &[r1, r2, supp]);

    // Clip only read-two's 5' end by 2bp. R2 is reverse, so its 5' end is the high-coordinate
    // (right) end: 100M @300 -> 98M2H, alignment start unchanged at 300 (1-based). The
    // supplementary R1 is never touched by the clip loop (only the primary pair is), so it keeps
    // 50M — but its mate fields must be rewritten to the post-clip primary R2.
    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--threads",
        "2",
        "--read-two-five-prime",
        "2",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");

    let recs = read_output_record_bufs(&output_bam);
    assert_eq!(recs.len(), 3, "all three records retained");

    let flag_bits = |r: &RecordBuf| u16::from(r.flags());
    let supp = recs
        .iter()
        .find(|r| flag_bits(r) & flags::SUPPLEMENTARY != 0)
        .expect("supplementary record present");
    let r2 = recs
        .iter()
        .find(|r| {
            let f = flag_bits(r);
            f & flags::LAST_SEGMENT != 0
                && f & flags::SUPPLEMENTARY == 0
                && f & flags::SECONDARY == 0
        })
        .expect("primary R2 present");

    // Independent oracle for the clip itself: primary R2 (reverse, 5' clip of 2) -> 98M2H @300,
    // start unchanged, MAPQ 40. Derived by hand from the fixed inputs, not from the code under test.
    assert_eq!(
        cigar_ops(r2),
        vec![(CigarKind::Match, 98), (CigarKind::HardClip, 2)],
        "primary R2 must be 98M2H after a 5' clip of 2"
    );
    assert_eq!(usize::from(r2.alignment_start().unwrap()), 300, "R2 alignment start unchanged");
    assert_eq!(r2.mapping_quality().map(u8::from), Some(40), "R2 MAPQ unchanged");
    assert_ne!(flag_bits(r2) & flags::REVERSE, 0, "primary R2 is reverse");

    // The supplementary alignment itself is not clipped (only the primary pair is).
    assert_eq!(cigar_ops(supp), vec![(CigarKind::Match, 50)], "supplementary unchanged (50M)");

    // Headline fix: the supplementary's mate metadata now points at the post-clip primary R2.
    // Absolute values are hand-derived from the inputs (mate ref 0, mate pos 300, MC 98M2H, MQ 40,
    // mate-reverse set) — so this is an independent oracle, not a self-consistency check.
    assert_eq!(supp.mate_reference_sequence_id(), Some(0), "mate ref = primary R2 ref");
    assert_eq!(
        usize::from(supp.mate_alignment_start().unwrap()),
        300,
        "mate pos = primary R2 start"
    );
    assert_ne!(
        flag_bits(supp) & flags::MATE_REVERSE,
        0,
        "mate-reverse must be set (primary R2 is reverse)"
    );
    assert_eq!(mate_cigar(supp).as_deref(), Some("98M2H"), "MC = primary R2 post-clip CIGAR");
    assert_eq!(mate_mapq(supp), Some(40), "MQ = primary R2 MAPQ");
    assert_eq!(
        supp.template_length(),
        -r2.template_length(),
        "supplementary TLEN = negation of primary R2 TLEN"
    );

    // Cross-check: the repaired mate fields agree with the actual primary R2 record, proving the
    // snapshot was taken from the correct (post-clip) primary and not a stale/wrong index.
    assert_eq!(supp.mate_reference_sequence_id(), r2.reference_sequence_id());
    assert_eq!(supp.mate_alignment_start(), r2.alignment_start());
}
