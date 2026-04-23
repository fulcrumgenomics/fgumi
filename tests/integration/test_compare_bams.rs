//! Integration tests for the `fgumi compare bams` command.
//!
//! These tests exercise all three compare modes (content, full, grouping)
//! using the raw byte comparison path.

use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

/// Writes a BAM file from the given header and records.
fn write_bam(path: &Path, header: &Header, records: &[RawRecord]) {
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(header).expect("Failed to write header");
    for record in records {
        writer
            .write_alignment_record(header, &to_record_buf(record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Runs `fgumi compare bams` and returns (success, stdout).
fn run_compare(bam1: &Path, bam2: &Path, mode: &str, extra_args: &[&str]) -> (bool, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["compare", "bams"])
        .arg(bam1)
        .arg(bam2)
        .args(["--mode", mode])
        .args(extra_args)
        .output()
        .expect("Failed to run fgumi");

    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    (output.status.success(), stdout)
}

/// Builds a simple mapped record with a given name and position.
fn mapped_record(name: &[u8], pos: i32) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name)
        .sequence(b"ACGTACGT")
        .qualities(&[30; 8])
        .ref_id(0)
        .pos(pos - 1) // pos is 1-based in tests, BAM uses 0-based
        .mapq(60);
    b.build()
}

/// Builds a simple mapped record with an MI tag.
fn mapped_record_with_mi(name: &[u8], pos: i32, mi: &str) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name)
        .sequence(b"ACGTACGT")
        .qualities(&[30; 8])
        .ref_id(0)
        .pos(pos - 1) // pos is 1-based in tests, BAM uses 0-based
        .mapq(60)
        .add_string_tag(SamTag::MI, mi.as_bytes());
    b.build()
}

// ---------------------------------------------------------------------------
// Content mode tests
// ---------------------------------------------------------------------------

#[test]
fn test_content_mode_identical_bams() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![mapped_record(b"read1", 100), mapped_record(b"read2", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, "content", &[]);
    assert!(success, "Expected success for identical BAMs, stdout:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
}

#[test]
fn test_content_mode_different_position() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![mapped_record(b"read1", 100)];
    let records2 = vec![mapped_record(b"read1", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, "content", &[]);
    assert!(!success, "Expected failure for different BAMs, stdout:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
}

#[test]
fn test_content_mode_different_record_count() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![mapped_record(b"read1", 100), mapped_record(b"read2", 200)];
    let records2 = vec![mapped_record(b"read1", 100)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, "content", &[]);
    assert!(!success, "Expected failure for different record counts, stdout:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
}

#[test]
fn test_content_mode_multithreaded() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records: Vec<RawRecord> =
        (0..20).map(|i| mapped_record(format!("read{i}").as_bytes(), 100 + i * 10)).collect();

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, "content", &["-t", "4"]);
    assert!(success, "Expected success for identical BAMs with threads, stdout:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
}

// ---------------------------------------------------------------------------
// Full mode tests
// ---------------------------------------------------------------------------

#[test]
fn test_full_mode_identical_bams() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![
        mapped_record_with_mi(b"read1", 100, "1"),
        mapped_record_with_mi(b"read2", 200, "1"),
        mapped_record_with_mi(b"read3", 300, "2"),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, "full", &[]);
    assert!(success, "Expected success for identical BAMs in full mode, stdout:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
}

#[test]
fn test_full_mode_different_mi_tags() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 =
        vec![mapped_record_with_mi(b"read1", 100, "1"), mapped_record_with_mi(b"read2", 200, "1")];
    let records2 =
        vec![mapped_record_with_mi(b"read1", 100, "1"), mapped_record_with_mi(b"read2", 200, "2")];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, "full", &[]);
    assert!(!success, "Expected failure for different MI tags in full mode, stdout:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
}

// ---------------------------------------------------------------------------
// Grouping mode tests
// ---------------------------------------------------------------------------

#[test]
fn test_grouping_mode_identical_bams() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // Build paired reads so that grouping mode can match R1/R2 flags.
    let records = vec![
        {
            let mut b = SamBuilder::new();
            b.read_name(b"read1")
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .add_string_tag(SamTag::MI, b"1")
                .add_string_tag(SamTag::RX, b"AAAA");
            b.build()
        },
        {
            let mut b = SamBuilder::new();
            b.read_name(b"read1")
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::LAST_SEGMENT)
                .ref_id(0)
                .pos(199)
                .mapq(60)
                .add_string_tag(SamTag::MI, b"1")
                .add_string_tag(SamTag::RX, b"AAAA");
            b.build()
        },
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, "grouping", &[]);
    assert!(success, "Expected success for identical grouped BAMs, stdout:\n{stdout}");
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT in output, got:\n{stdout}");
}

#[test]
fn test_grouping_mode_flag_mismatch_reported_separately() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // BAM1: read1 is R1, read2 is R2. BAM2: read1 is R2, read2 is R1 (swapped).
    // Names and order match, but R1/R2 flags do not — this is a flag mismatch,
    // not an order/name mismatch.
    let make = |name: &[u8], pos: i32, flags: u16, mi: &str| -> RawRecord {
        let mut b = SamBuilder::new();
        b.read_name(name)
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags)
            .ref_id(0)
            .pos(pos - 1)
            .mapq(60)
            .add_string_tag(SamTag::MI, mi.as_bytes());
        b.build()
    };

    let records1 = vec![
        make(b"read1", 100, flags::PAIRED | flags::FIRST_SEGMENT, "1"),
        make(b"read2", 200, flags::PAIRED | flags::LAST_SEGMENT, "1"),
    ];
    let records2 = vec![
        make(b"read1", 100, flags::PAIRED | flags::LAST_SEGMENT, "1"),
        make(b"read2", 200, flags::PAIRED | flags::FIRST_SEGMENT, "1"),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, "grouping", &[]);
    assert!(!success, "Expected failure for flag mismatches, stdout:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
    assert!(
        stdout.contains("Order/name mismatches: 0"),
        "Expected order/name count to remain zero, got:\n{stdout}"
    );
    assert!(
        stdout.contains("R1/R2 flag mismatches: 2"),
        "Expected flag mismatch count of 2, got:\n{stdout}"
    );
}

// ---------------------------------------------------------------------------
// Reordered tags tests (content and full modes)
// ---------------------------------------------------------------------------

#[test]
fn test_content_mode_reordered_tags_equivalent() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // Same record but tags in different order
    let records1 = vec![{
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .add_string_tag(SamTag::MI, b"1")
            .add_string_tag(SamTag::RX, b"AAAA");
        b.build()
    }];
    let records2 = vec![{
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .add_string_tag(SamTag::RX, b"AAAA")
            .add_string_tag(SamTag::MI, b"1");
        b.build()
    }];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, "content", &[]);
    assert!(success, "Expected success for reordered tags in content mode, stdout:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
    assert!(
        stdout.contains("tags in different order"),
        "Expected tag order diff note in output, got:\n{stdout}"
    );
}

#[test]
fn test_full_mode_reordered_tags_equivalent() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![{
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .add_string_tag(SamTag::MI, b"1")
            .add_string_tag(SamTag::RX, b"AAAA");
        b.build()
    }];
    let records2 = vec![{
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .add_string_tag(SamTag::RX, b"AAAA")
            .add_string_tag(SamTag::MI, b"1");
        b.build()
    }];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, "full", &[]);
    assert!(success, "Expected success for reordered tags in full mode, stdout:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
    assert!(
        stdout.contains("tags in different order"),
        "Expected tag order diff note in output, got:\n{stdout}"
    );
}

// ---------------------------------------------------------------------------
// Ignore-order grouping mode tests
// ---------------------------------------------------------------------------

#[test]
fn test_grouping_mode_ignore_order() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // BAM1: R1 then R2 for each pair
    let records1 = vec![
        {
            let mut b = SamBuilder::new();
            b.read_name(b"read1")
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .add_string_tag(SamTag::MI, b"1");
            b.build()
        },
        {
            let mut b = SamBuilder::new();
            b.read_name(b"read1")
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::LAST_SEGMENT)
                .ref_id(0)
                .pos(199)
                .mapq(60)
                .add_string_tag(SamTag::MI, b"1");
            b.build()
        },
    ];

    // BAM2: R2 then R1 (reversed order), same MI grouping
    let records2 = vec![
        {
            let mut b = SamBuilder::new();
            b.read_name(b"read1")
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::LAST_SEGMENT)
                .ref_id(0)
                .pos(199)
                .mapq(60)
                .add_string_tag(SamTag::MI, b"5");
            b.build()
        },
        {
            let mut b = SamBuilder::new();
            b.read_name(b"read1")
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .add_string_tag(SamTag::MI, b"5");
            b.build()
        },
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, "grouping", &["--ignore-order"]);
    assert!(success, "Expected success for ignore-order grouping mode, stdout:\n{stdout}");
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT in output, got:\n{stdout}");
}

// ---------------------------------------------------------------------------
// Missing MI tag regression tests
//
// When both BAMs lack MI tags, the compare should NOT report them as
// equivalent under grouping or full modes: no grouping has been verified.
// ---------------------------------------------------------------------------

#[test]
fn test_grouping_mode_fails_when_mi_missing_in_both_bams() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    // Paired reads with NO MI tag on either record.
    let records = vec![
        {
            let mut b = SamBuilder::new();
            b.read_name(b"read1")
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(99)
                .mapq(60);
            b.build()
        },
        {
            let mut b = SamBuilder::new();
            b.read_name(b"read1")
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::LAST_SEGMENT)
                .ref_id(0)
                .pos(199)
                .mapq(60);
            b.build()
        },
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, "grouping", &[]);
    assert!(!success, "Expected failure when MI tags are missing in both BAMs, stdout:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
    assert!(
        stdout.contains("Missing MI in BAM1: 2") && stdout.contains("Missing MI in BAM2: 2"),
        "Expected missing-MI counts of 2 each, got:\n{stdout}"
    );
}

#[test]
fn test_full_mode_fails_when_mi_missing_in_both_bams() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    // Two identical records, neither has an MI tag.
    let records = vec![mapped_record(b"read1", 100), mapped_record(b"read2", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, "full", &[]);
    assert!(
        !success,
        "Expected failure when MI tags are missing in both BAMs (full mode), stdout:\n{stdout}"
    );
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
}

#[test]
fn test_grouping_unordered_fails_when_mi_missing_in_both_bams() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![mapped_record(b"read1", 100), mapped_record(b"read2", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, "grouping", &["--ignore-order"]);
    assert!(
        !success,
        "Expected failure when MI tags are missing in both BAMs (ignore-order), stdout:\n{stdout}"
    );
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
}

// ---------------------------------------------------------------------------
// Paired-UMI MI strand-suffix regression tests
//
// Paired UMI grouping (fgumi and fgbio) emits MI as a Z-type string of the
// form `<id>/<A|B>`, where the suffix distinguishes the two strand
// orientations of the same double-stranded molecule. The comparator must:
//   1. Recognise the full encoding as a valid MI (not "missing").
//   2. Preserve the A/B distinction when checking grouping equivalence.
// ---------------------------------------------------------------------------

/// Build a paired record pair (R1 + R2) sharing a single MI tag value.
fn paired_record_pair(name: &[u8], mi: &str) -> [RawRecord; 2] {
    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(name)
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .add_string_tag(SamTag::MI, mi.as_bytes())
            .add_string_tag(SamTag::RX, b"AAAA");
        b.build()
    };
    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(name)
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT)
            .ref_id(0)
            .pos(199)
            .mapq(60)
            .add_string_tag(SamTag::MI, mi.as_bytes())
            .add_string_tag(SamTag::RX, b"AAAA");
        b.build()
    };
    [r1, r2]
}

#[test]
fn test_grouping_mode_paired_strand_suffix_equivalent() {
    // Two identical paired-UMI BAMs: one /A-strand pair and one /B-strand pair,
    // same base molecule id. Comparator must recognise MI:Z:0/A and 0/B as
    // present (not missing) and report EQUIVALENT.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let mut records: Vec<RawRecord> = Vec::new();
    records.extend(paired_record_pair(b"read1", "0/A"));
    records.extend(paired_record_pair(b"read2", "0/B"));

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, "grouping", &[]);
    assert!(success, "Expected EQUIVALENT for identical paired-strand BAMs, stdout:\n{stdout}");
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT, got:\n{stdout}");
    // Guards against regression to the old `parse::<i64>()` path where every
    // /A /B record was reported as missing MI.
    assert!(
        stdout.contains("Missing MI in BAM1: 0") && stdout.contains("Missing MI in BAM2: 0"),
        "Expected zero missing-MI counts for paired-strand MIs, got:\n{stdout}"
    );
}

#[test]
fn test_grouping_unordered_paired_strand_suffix_equivalent() {
    // Same as above but via the --ignore-order code path that builds
    // per-BAM MI maps in parallel.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let mut records: Vec<RawRecord> = Vec::new();
    records.extend(paired_record_pair(b"read1", "0/A"));
    records.extend(paired_record_pair(b"read2", "0/B"));

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, "grouping", &["--ignore-order"]);
    assert!(success, "Expected EQUIVALENT for ignore-order paired-strand BAMs, stdout:\n{stdout}");
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT, got:\n{stdout}");
    assert!(
        stdout.contains("Missing MI in BAM1: 0") && stdout.contains("Missing MI in BAM2: 0"),
        "Expected zero missing-MI counts for paired-strand MIs, got:\n{stdout}"
    );
}

#[test]
fn test_grouping_mode_paired_strand_a_flipped_to_b_is_mismatch() {
    // BAM1 groups the second read as /B; BAM2 groups it as /A. The /A vs /B
    // swap must be surfaced as a grouping mismatch (not silently accepted).
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let mut records1: Vec<RawRecord> = Vec::new();
    records1.extend(paired_record_pair(b"read1", "0/A"));
    records1.extend(paired_record_pair(b"read2", "0/B"));

    let mut records2: Vec<RawRecord> = Vec::new();
    records2.extend(paired_record_pair(b"read1", "0/A"));
    records2.extend(paired_record_pair(b"read2", "0/A"));

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    // Ordered grouping mode (no --ignore-order) is used deliberately so that
    // the mismatch-sample path (which prints the conflicting MI values) is
    // exercised, letting us assert they render via Display instead of Debug.
    let (success, stdout) = run_compare(&bam1, &bam2, "grouping", &[]);
    assert!(!success, "Expected DIFFER when /A and /B assignments disagree, stdout:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "Expected DIFFER, got:\n{stdout}");
    // Specifically assert the mismatch reason: the paired MIs must be parsed as
    // present (zero missing MI) and the failure must be a grouping mismatch, so
    // regressing to the old "every /A /B counted as missing MI" path still trips
    // this test.
    assert!(
        stdout.contains("Missing MI in BAM1: 0") && stdout.contains("Missing MI in BAM2: 0"),
        "Expected paired-strand MIs to be parsed as present, got:\n{stdout}"
    );
    assert!(
        stdout.contains("Grouping mismatches: 1"),
        "Expected the failure to be a grouping mismatch, got:\n{stdout}"
    );
    // The sampled MIs in the mismatch detail must render via Display
    // (`0/A`, `0/B`) not Debug (`Strand { base: 0, strand: 65 }`).
    assert!(
        stdout.contains("0/A") && stdout.contains("0/B"),
        "Expected Display-formatted MiKeys in mismatch detail, got:\n{stdout}"
    );
    assert!(
        !stdout.contains("Strand {"),
        "MiKey debug-formatting leaked into user-facing output:\n{stdout}"
    );
}
