//! Integration tests for the `fgumi compare bams` command.
//!
//! These tests exercise all three compare modes (content, full, grouping)
//! using the raw byte comparison path.

use fgumi_lib::sam::builder::RecordBuilder;
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::create_minimal_header;

/// Writes a BAM file from the given header and records.
fn write_bam(path: &Path, header: &Header, records: &[RecordBuf]) {
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(header).expect("Failed to write header");
    for record in records {
        writer.write_alignment_record(header, record).expect("Failed to write record");
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
fn mapped_record(name: &str, pos: usize) -> RecordBuf {
    RecordBuilder::new()
        .name(name)
        .sequence("ACGTACGT")
        .qualities(&[30; 8])
        .reference_sequence_id(0)
        .alignment_start(pos)
        .mapping_quality(60)
        .build()
}

/// Builds a simple mapped record with an MI tag.
fn mapped_record_with_mi(name: &str, pos: usize, mi: &str) -> RecordBuf {
    RecordBuilder::new()
        .name(name)
        .sequence("ACGTACGT")
        .qualities(&[30; 8])
        .reference_sequence_id(0)
        .alignment_start(pos)
        .mapping_quality(60)
        .tag("MI", mi)
        .build()
}

// ---------------------------------------------------------------------------
// Content mode tests
// ---------------------------------------------------------------------------

#[test]
fn test_content_mode_identical_bams() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![mapped_record("read1", 100), mapped_record("read2", 200)];

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

    let records1 = vec![mapped_record("read1", 100)];
    let records2 = vec![mapped_record("read1", 200)];

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

    let records1 = vec![mapped_record("read1", 100), mapped_record("read2", 200)];
    let records2 = vec![mapped_record("read1", 100)];

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
    let records: Vec<RecordBuf> =
        (0..20).map(|i| mapped_record(&format!("read{i}"), 100 + i * 10)).collect();

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
        mapped_record_with_mi("read1", 100, "1"),
        mapped_record_with_mi("read2", 200, "1"),
        mapped_record_with_mi("read3", 300, "2"),
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
        vec![mapped_record_with_mi("read1", 100, "1"), mapped_record_with_mi("read2", 200, "1")];
    let records2 =
        vec![mapped_record_with_mi("read1", 100, "1"), mapped_record_with_mi("read2", 200, "2")];

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
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .paired(true)
            .first_segment(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .tag("MI", "1")
            .tag("RX", "AAAA")
            .build(),
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .paired(true)
            .first_segment(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mapping_quality(60)
            .tag("MI", "1")
            .tag("RX", "AAAA")
            .build(),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, "grouping", &[]);
    assert!(success, "Expected success for identical grouped BAMs, stdout:\n{stdout}");
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT in output, got:\n{stdout}");
}

// ---------------------------------------------------------------------------
// Reordered tags tests (content and full modes)
// ---------------------------------------------------------------------------

#[test]
fn test_content_mode_reordered_tags_equivalent() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // Same record but tags in different order
    let records1 = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .tag("MI", "1")
            .tag("RX", "AAAA")
            .build(),
    ];
    let records2 = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .tag("RX", "AAAA")
            .tag("MI", "1")
            .build(),
    ];

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

    let records1 = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .tag("MI", "1")
            .tag("RX", "AAAA")
            .build(),
    ];
    let records2 = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .tag("RX", "AAAA")
            .tag("MI", "1")
            .build(),
    ];

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
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .paired(true)
            .first_segment(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .tag("MI", "1")
            .build(),
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .paired(true)
            .first_segment(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mapping_quality(60)
            .tag("MI", "1")
            .build(),
    ];

    // BAM2: R2 then R1 (reversed order), same MI grouping
    let records2 = vec![
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .paired(true)
            .first_segment(false)
            .reference_sequence_id(0)
            .alignment_start(200)
            .mapping_quality(60)
            .tag("MI", "5")
            .build(),
        RecordBuilder::new()
            .name("read1")
            .sequence("ACGTACGT")
            .qualities(&[30; 8])
            .paired(true)
            .first_segment(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .mapping_quality(60)
            .tag("MI", "5")
            .build(),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, "grouping", &["--ignore-order"]);
    assert!(success, "Expected success for ignore-order grouping mode, stdout:\n{stdout}");
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT in output, got:\n{stdout}");
}
