//! Integration tests for the `fgumi compare bams` command.
//!
//! These tests exercise the three CLI behaviors:
//! - default (no flags): strict content comparison (all fields + tags, including MI)
//! - `--check-grouping`: MI grouping equivalence + content comparison excluding MI
//! - `--check-grouping --ignore-order`: MI grouping equivalence only, unordered

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

/// Runs `fgumi compare bams` with the given args and returns (success, stdout).
fn run_compare(bam1: &Path, bam2: &Path, args: &[&str]) -> (bool, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["compare", "bams"])
        .arg(bam1)
        .arg(bam2)
        .args(args)
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
// Default mode tests (strict content comparison)
// ---------------------------------------------------------------------------

#[test]
fn test_default_identical_bams() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![mapped_record("read1", 100), mapped_record("read2", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, &[]);
    assert!(success, "Expected success for identical BAMs, stdout:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
}

#[test]
fn test_default_different_position() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![mapped_record("read1", 100)];
    let records2 = vec![mapped_record("read1", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, &[]);
    assert!(!success, "Expected failure for different BAMs, stdout:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
}

#[test]
fn test_default_different_record_count() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![mapped_record("read1", 100), mapped_record("read2", 200)];
    let records2 = vec![mapped_record("read1", 100)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, &[]);
    assert!(!success, "Expected failure for different record counts, stdout:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
}

#[test]
fn test_default_multithreaded() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records: Vec<RecordBuf> =
        (0..20).map(|i| mapped_record(&format!("read{i}"), 100 + i * 10)).collect();

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, &["-t", "4"]);
    assert!(success, "Expected success for identical BAMs with threads, stdout:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
}

#[test]
fn test_default_different_mi_tags_fails() {
    // In default mode, MI is part of the content comparison, so different MI values should fail.
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

    let (success, stdout) = run_compare(&bam1, &bam2, &[]);
    assert!(!success, "Expected failure for different MI in default mode, stdout:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
}

#[test]
fn test_default_reordered_tags_equivalent() {
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

    let (success, stdout) = run_compare(&bam1, &bam2, &[]);
    assert!(success, "Expected success for reordered tags, stdout:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
    assert!(
        stdout.contains("tags in different order"),
        "Expected tag order diff note in output, got:\n{stdout}"
    );
}

// ---------------------------------------------------------------------------
// --check-grouping tests (MI grouping + content excluding MI)
// ---------------------------------------------------------------------------

#[test]
fn test_check_grouping_identical_bams() {
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

    let (success, stdout) = run_compare(&bam1, &bam2, &["--check-grouping"]);
    assert!(
        success,
        "Expected success for identical BAMs with --check-grouping, stdout:\n{stdout}"
    );
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
}

#[test]
fn test_check_grouping_different_mi_values_same_grouping_ok() {
    // The key test for --check-grouping: different MI values with equivalent grouping
    // should succeed (MI is excluded from content comparison, grouping is verified separately).
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // BAM1: read1 and read2 share MI=1; read3 has MI=2
    let records1 = vec![
        mapped_record_with_mi("read1", 100, "1"),
        mapped_record_with_mi("read2", 200, "1"),
        mapped_record_with_mi("read3", 300, "2"),
    ];
    // BAM2: same grouping, but MI values are 7 and 9 (different from BAM1 but equivalent)
    let records2 = vec![
        mapped_record_with_mi("read1", 100, "7"),
        mapped_record_with_mi("read2", 200, "7"),
        mapped_record_with_mi("read3", 300, "9"),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, &["--check-grouping"]);
    assert!(
        success,
        "Expected success for equivalent grouping with different MI values, stdout:\n{stdout}"
    );
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
    assert!(
        !stdout.contains("tags in different order"),
        "MI-only drift should not report tag order differences, got:\n{stdout}"
    );
}

#[test]
fn test_check_grouping_different_content_fails() {
    // Different non-MI content (POS) should still fail with --check-grouping.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![mapped_record_with_mi("read1", 100, "1")];
    let records2 = vec![mapped_record_with_mi("read1", 500, "1")];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, &["--check-grouping"]);
    assert!(!success, "Expected failure for different POS, stdout:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
}

#[test]
fn test_check_grouping_inequivalent_grouping_fails() {
    // Groupings that are not equivalent should fail.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // BAM1: read1 and read2 share MI=1
    let records1 =
        vec![mapped_record_with_mi("read1", 100, "1"), mapped_record_with_mi("read2", 200, "1")];
    // BAM2: read1 and read2 have different MIs (not equivalent grouping)
    let records2 =
        vec![mapped_record_with_mi("read1", 100, "1"), mapped_record_with_mi("read2", 200, "2")];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, &["--check-grouping"]);
    assert!(!success, "Expected failure for inequivalent grouping, stdout:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
}

// ---------------------------------------------------------------------------
// --check-grouping --ignore-order tests (MI grouping equivalence only)
// ---------------------------------------------------------------------------

#[test]
fn test_check_grouping_ignore_order_reversed_records() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // BAM1: R1 then R2 for each pair, MI=1
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

    // BAM2: R2 then R1 (reversed order), MI=5 (different value, equivalent grouping)
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

    let (success, stdout) = run_compare(&bam1, &bam2, &["--check-grouping", "--ignore-order"]);
    assert!(success, "Expected success for --check-grouping --ignore-order, stdout:\n{stdout}");
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT in output, got:\n{stdout}");
}

#[test]
fn test_ignore_order_without_check_grouping_errors() {
    // --ignore-order requires --check-grouping
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![mapped_record("read1", 100)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, _stdout) = run_compare(&bam1, &bam2, &["--ignore-order"]);
    assert!(!success, "Expected failure when --ignore-order used without --check-grouping");
}

// ---------------------------------------------------------------------------
// --command preset tests
// ---------------------------------------------------------------------------

#[test]
fn test_command_preset_extract() {
    // extract preset: check_grouping=false, ignore_order=false. Records with no MI.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![mapped_record("read1", 100), mapped_record("read2", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare(&bam1, &bam2, &["--command", "extract"]);
    assert!(success, "Expected success for --command extract, stdout:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
}

#[test]
fn test_command_preset_simplex_allows_different_mi_values() {
    // simplex preset: check_grouping=true. Different MI values should pass if grouping matches.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records1 =
        vec![mapped_record_with_mi("read1", 100, "1"), mapped_record_with_mi("read2", 200, "1")];
    let records2 =
        vec![mapped_record_with_mi("read1", 100, "42"), mapped_record_with_mi("read2", 200, "42")];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare(&bam1, &bam2, &["--command", "simplex"]);
    assert!(success, "Expected success for --command simplex, stdout:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
}

#[test]
fn test_command_preset_group_different_mi_equivalent_grouping() {
    // group preset: check_grouping=true, ignore_order=true. Different MI values, reversed order.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

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

    let (success, stdout) = run_compare(&bam1, &bam2, &["--command", "group"]);
    assert!(success, "Expected success for --command group, stdout:\n{stdout}");
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT in output, got:\n{stdout}");
}

#[test]
fn test_command_preset_override_with_explicit_flag() {
    // Override: --command extract defaults to check_grouping=false, but --check-grouping
    // explicitly sets it to true. Records with different MI values should pass.
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records1 =
        vec![mapped_record_with_mi("read1", 100, "1"), mapped_record_with_mi("read2", 200, "1")];
    let records2 =
        vec![mapped_record_with_mi("read1", 100, "9"), mapped_record_with_mi("read2", 200, "9")];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    // With --command extract alone (check_grouping=false), different MI → fail
    let (success, _) = run_compare(&bam1, &bam2, &["--command", "extract"]);
    assert!(!success, "Expected failure with --command extract and different MI");

    // With --command extract --check-grouping (explicit override), different MI → pass
    let (success, stdout) =
        run_compare(&bam1, &bam2, &["--command", "extract", "--check-grouping"]);
    assert!(success, "Expected success when --check-grouping overrides preset, stdout:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
}
