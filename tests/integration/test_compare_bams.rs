//! Integration tests for the `fgumi compare bams` command.
//!
//! These tests exercise all three compare modes (content, full, grouping)
//! using the raw byte comparison path.

use clap::Parser;
use fgumi_lib::commands::command::Command as FgumiCommand;
use fgumi_lib::commands::compare::{
    CompareBams, CompareMismatch, ContentPredicate, KeyJoinConfig, KeyJoinOutcome,
    canonicalize_to_queryname, keyjoin_compare, positional_compare, sort_verify_compare,
};
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::sam::Header;
use rstest::rstest;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_coordinate_sorted_header, create_minimal_header, keyjoin_cfg, mi_record, write_bam,
};

/// Runs `fgumi compare bams` via subprocess and returns (success, stdout).
///
/// Used by tests that assert on specific substrings in the diagnostic report
/// (e.g. EQUIVALENT vs IDENTICAL, mismatch counts). For tests that only need
/// to know whether the BAMs matched, prefer `run_compare_in_process`.
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

/// Runs `CompareBams::execute()` in-process and returns true on match
/// (IDENTICAL or EQUIVALENT), false on `CompareMismatch` (DIFFER).
///
/// Panics if `execute` returns any other error.
fn run_compare_in_process(bam1: &Path, bam2: &Path, mode: &str, extra_args: &[&str]) -> bool {
    let mut args: Vec<&str> =
        vec!["bams", bam1.to_str().unwrap(), bam2.to_str().unwrap(), "--mode", mode];
    args.extend(extra_args);
    let cmd = CompareBams::try_parse_from(args).expect("failed to parse compare bams args");
    match cmd.execute("fgumi compare bams") {
        Ok(()) => true,
        Err(e) if e.is::<CompareMismatch>() => false,
        Err(e) => panic!("compare bams hit unexpected error: {e:#}"),
    }
}

/// Runs `fgumi compare bams --command <command>` via subprocess and returns (success,
/// stdout).
///
/// Unlike `run_compare` (which sets `--mode` directly), this drives the real
/// `--command <stage>` CLI surface — i.e. the path an actual `group` cross-tool
/// comparison would take — so tests using this helper exercise `CommandPreset`'s
/// defaults/`content_predicate` dispatch end to end, not just the underlying engine.
fn run_compare_command(
    bam1: &Path,
    bam2: &Path,
    command: &str,
    extra_args: &[&str],
) -> (bool, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["compare", "bams"])
        .arg(bam1)
        .arg(bam2)
        .args(["--command", command])
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

// `mapped_record_with_mi` is now shared as `mi_record` in `helpers::bam_generator`
// (see the top-of-file import) since `test_compare_mutation.rs` had an identical
// duplicate.

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

    assert!(
        run_compare_in_process(&bam1, &bam2, "content", &[]),
        "Expected match for identical BAMs in content mode"
    );
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

    assert!(
        !run_compare_in_process(&bam1, &bam2, "content", &[]),
        "Expected mismatch for different positions in content mode"
    );
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

    assert!(
        !run_compare_in_process(&bam1, &bam2, "content", &[]),
        "Expected mismatch for different record counts in content mode"
    );
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

    assert!(
        run_compare_in_process(&bam1, &bam2, "content", &["-t", "4"]),
        "Expected match for identical BAMs with --threads 4"
    );
}

// ---------------------------------------------------------------------------
// Full mode tests
// ---------------------------------------------------------------------------

#[test]
fn test_full_mode_identical_bams() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![
        mi_record(b"read1", 100, "1"),
        mi_record(b"read2", 200, "1"),
        mi_record(b"read3", 300, "2"),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    assert!(
        run_compare_in_process(&bam1, &bam2, "full", &[]),
        "Expected match for identical BAMs in full mode"
    );
}

#[test]
fn test_full_mode_different_mi_tags() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];
    let records2 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "2")];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    assert!(
        !run_compare_in_process(&bam1, &bam2, "full", &[]),
        "Expected mismatch for different MI tags in full mode"
    );
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

    // `content` mode now delegates to the positional engine (`positional_compare`),
    // whose `PositionalOutcome` does not carry a "tags in different order" counter —
    // only whether the pair is content-equal under the configured `ContentPredicate`
    // (which is itself already order-independent). The report-contract preserved
    // across that reimplementation is the `RESULT: ... IDENTICAL`/`DIFFER` line and
    // the exit code, not this now-retired note; see the compare-hardening Phase 2
    // plan's report-contract fix.
    let (success, stdout) = run_compare(&bam1, &bam2, "content", &[]);
    assert!(success, "Expected success for reordered tags in content mode, stdout:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
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

    assert!(
        !run_compare_in_process(&bam1, &bam2, "full", &[]),
        "Expected mismatch when MI tags are missing in both BAMs (full mode)"
    );
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

    assert!(
        !run_compare_in_process(&bam1, &bam2, "grouping", &["--ignore-order"]),
        "Expected mismatch when MI tags are missing in both BAMs (ignore-order)"
    );
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

// ---------------------------------------------------------------------------
// Positional engine tests (key-lockstep + presence, no resync)
//
// `positional_compare` is not yet wired into any CLI preset (that is Task
// 2.3); these tests call the engine directly to verify its sound core:
// paired records must match by `RecordKey` before content is ever compared,
// and a key mismatch stops pairing rather than attempting to resync.
// ---------------------------------------------------------------------------

/// Builds a mapped, paired-flag record with a given name/segment/position.
fn paired_record(name: &[u8], segment_flag: u16, pos: i32) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name)
        .sequence(b"ACGTACGT")
        .qualities(&[30; 8])
        .flags(flags::PAIRED | segment_flag)
        .ref_id(0)
        .pos(pos - 1)
        .mapq(60);
    b.build()
}

#[test]
fn test_positional_compare_identical_bams_match() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![
        paired_record(b"read1", flags::FIRST_SEGMENT, 100),
        paired_record(b"read2", flags::FIRST_SEGMENT, 200),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 100, ContentPredicate::Exact)
        .expect("positional_compare should succeed");

    assert_eq!(outcome.bam1_count, 2);
    assert_eq!(outcome.bam2_count, 2);
    assert_eq!(outcome.content_diffs, 0);
    assert_eq!(outcome.key_mismatch_at, None);
    assert!(outcome.is_match(), "identical BAMs must be a positional match: {outcome:?}");
}

#[test]
fn test_positional_compare_seq_difference_is_one_content_diff() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![paired_record(b"read1", flags::FIRST_SEGMENT, 100)];
    let records2 = vec![{
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGA") // last base differs
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(99)
            .mapq(60);
        b.build()
    }];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 100, ContentPredicate::Exact)
        .expect("positional_compare should succeed");

    assert_eq!(outcome.bam1_count, 1);
    assert_eq!(outcome.bam2_count, 1);
    assert_eq!(outcome.content_diffs, 1, "a single SEQ base difference must count as one diff");
    assert_eq!(outcome.key_mismatch_at, None, "same RecordKey on both sides, no key mismatch");
    assert!(!outcome.is_match(), "a content diff must not be reported as a match: {outcome:?}");
}

#[test]
fn test_positional_compare_swapped_records_is_key_mismatch() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // Two records with distinct RecordKeys (different QNAME AND different
    // segment), swapped between the two files. Index 0 pairs read1/FIRST
    // against read2/LAST — a genuine key mismatch that must be caught
    // immediately, not silently compared as content or resynced.
    let read1_r1 = paired_record(b"read1", flags::FIRST_SEGMENT, 100);
    let read2_r2 = paired_record(b"read2", flags::LAST_SEGMENT, 200);

    let records1 = vec![read1_r1.clone(), read2_r2.clone()];
    let records2 = vec![read2_r2, read1_r1];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 100, ContentPredicate::Exact)
        .expect("positional_compare should succeed");

    assert_eq!(outcome.bam1_count, 2);
    assert_eq!(outcome.bam2_count, 2);
    assert_eq!(
        outcome.key_mismatch_at,
        Some(0),
        "the swap must be caught as a key mismatch at index 0, not silently paired"
    );
    assert!(!outcome.is_match(), "a key mismatch must never be reported as a match: {outcome:?}");
}

#[test]
fn test_positional_compare_extra_trailing_record_is_presence_differ() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![
        paired_record(b"read1", flags::FIRST_SEGMENT, 100),
        paired_record(b"read2", flags::FIRST_SEGMENT, 200),
    ];
    let records2 = vec![paired_record(b"read1", flags::FIRST_SEGMENT, 100)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 100, ContentPredicate::Exact)
        .expect("positional_compare should succeed");

    assert_eq!(outcome.bam1_count, 2);
    assert_eq!(outcome.bam2_count, 1);
    assert!(
        !outcome.is_match(),
        "differing record counts must be a presence DIFFER, not a match: {outcome:?}"
    );
}

/// R2 header-comparison gap, wired into the content/positional engine: two BAMs whose
/// records are byte-identical must still DIFFER if their `@SQ` reference dictionaries
/// disagree (here, on length) — a record-level match says nothing about whether the two
/// files actually declare the same reference.
#[test]
fn test_positional_compare_header_sq_length_mismatch_is_a_diff() {
    let tmp = TempDir::new().unwrap();
    let header1 = create_minimal_header("chr1", 10000);
    let header2 = create_minimal_header("chr1", 20000);

    let records = vec![paired_record(b"read1", flags::FIRST_SEGMENT, 100)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header1, &records);
    write_bam(&bam2, &header2, &records);

    let outcome = positional_compare(&bam1, &bam2, 1, 64, 100, ContentPredicate::Exact)
        .expect("positional_compare should succeed");

    assert!(outcome.header_mismatch, "a @SQ length mismatch must be flagged: {outcome:?}");
    assert!(
        !outcome.is_match(),
        "identical records under mismatched @SQ headers must still DIFFER: {outcome:?}"
    );
}

/// Reads back every record's name from a BAM file, in on-disk order, via the raw
/// (non-noodles) reader that the sort/key-join engines use internally.
fn read_all_names(path: &Path) -> Vec<String> {
    let file = fs::File::open(path).expect("open canonicalized BAM");
    let mut reader = fgumi_sort::RawBamRecordReader::new(file).expect("open raw reader");
    reader.skip_header().expect("skip header");
    let mut names = Vec::new();
    while let Some(record) = reader.next_record().expect("read record") {
        names.push(String::from_utf8_lossy(record.read_name()).into_owned());
    }
    names
}

/// `--command group`'s key-join canonicalization step must (a) preserve the exact
/// set of input records and (b) leave the output in non-decreasing queryname order,
/// whether or not the configured memory budget forces the sort to spill to disk.
#[rstest]
#[case::ample_memory(1024 * 1024, false)]
#[case::tiny_memory_forces_spill(1, true)]
fn test_canonicalize_to_queryname_sorts_and_preserves_records(
    #[case] sort_memory: usize,
    #[case] expect_spill: bool,
) {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let input_names = ["read3", "read1", "read2", "read0", "read4"];
    let records: Vec<RawRecord> =
        input_names.iter().map(|n| mapped_record(n.as_bytes(), 100)).collect();
    let input = tmp.path().join("in.bam");
    write_bam(&input, &header, &records);

    let output = tmp.path().join("out.bam");
    let cfg = KeyJoinConfig { threads: 1, sort_memory, sort_tmp_dirs: vec![], max_diffs: 10 };
    let stats =
        canonicalize_to_queryname(&input, &output, &cfg).expect("canonicalize should succeed");

    assert_eq!(stats.total_records, input_names.len() as u64);
    assert_eq!(stats.output_records, input_names.len() as u64);
    assert_eq!(
        stats.chunks_written > 0,
        expect_spill,
        "chunks_written={} for sort_memory={sort_memory}",
        stats.chunks_written
    );

    let got_names = read_all_names(&output);
    let mut expected_order = got_names.clone();
    expected_order.sort();
    assert_eq!(got_names, expected_order, "canonicalized output must be queryname-sorted");

    let mut got_set = got_names;
    got_set.sort();
    let mut want_set: Vec<String> = input_names.iter().map(|s| (*s).to_string()).collect();
    want_set.sort();
    assert_eq!(got_set, want_set, "canonicalization must preserve the record set");
}

// ---------------------------------------------------------------------------
// keyjoin engine tests (merge-join + MI bijection for `group`'s key-join)
//
// `keyjoin_compare` is not yet wired into any CLI preset (that is Task 3.4); these
// tests call the engine directly to verify its sound core: the merge-join tolerates
// record reordering and MI renumbering (the whole point of `group`'s cross-tool
// comparison) while still catching a genuine content diff, a genuine MI-grouping
// divergence, a presence-only record, and (F1) a canonicalized stream that violates
// `RecordKey` ordering.
// ---------------------------------------------------------------------------

// `test_keyjoin_config` is now shared as `keyjoin_cfg` in `helpers::bam_generator`
// (see the top-of-file import) since `test_compare_mutation.rs` had an identical
// duplicate.

/// Builds a simple mapped record with a given name, position, MI tag, and SEQ, so
/// content-diff cases can mutate SEQ independently of the MI value under test.
fn mapped_record_with_mi_and_seq(name: &[u8], pos: i32, mi: &str, seq: &[u8]) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name)
        .sequence(seq)
        .qualities(&vec![30; seq.len()])
        .ref_id(0)
        .pos(pos - 1)
        .mapq(60)
        .add_string_tag(SamTag::MI, mi.as_bytes());
    b.build()
}

/// Builds a minimal record with an arbitrary raw flag word, for exercising
/// `RecordKey`'s segment classification directly (see
/// `test_keyjoin_compare_hard_errors_on_record_key_order_violation`).
fn record_with_flags(name: &[u8], raw_flags: u16) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name).sequence(b"ACGTACGT").qualities(&[30; 8]).flags(raw_flags);
    b.build()
}

/// (a) Identical grouping, records in a different physical order, same MI numbering on
/// both sides: canonicalization erases the ordering difference, so this must MATCH.
#[test]
fn test_keyjoin_compare_reordered_records_same_mi_numbering_match() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let r1 = mi_record(b"read1", 100, "1");
    let r2 = mi_record(b"read2", 200, "1");
    let r3 = mi_record(b"read3", 300, "2");

    let records1 = vec![r1.clone(), r2.clone(), r3.clone()];
    let records2 = vec![r3, r1, r2]; // same records, different physical order

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let cfg = keyjoin_cfg(&tmp);
    let outcome: KeyJoinOutcome =
        keyjoin_compare(&bam1, &bam2, &cfg).expect("keyjoin_compare should succeed");

    assert_eq!(outcome.bam1_count, 3);
    assert_eq!(outcome.bam2_count, 3);
    assert_eq!(outcome.matched, 3);
    assert_eq!(outcome.content_diffs, 0);
    assert_eq!(outcome.mi_bijection_mismatches, 0);
    assert!(outcome.is_match(), "reordered, identically-grouped BAMs must match: {outcome:?}");
}

/// (b) Same grouping partition, different MI numbering (fgumi assigns `1`/`2`; fgbio
/// assigns `5`/`9` for the same two groups), content otherwise identical: MI renumbering
/// alone must not be reported as a diff — this is the key case `group` comparison exists
/// to tolerate.
#[test]
fn test_keyjoin_compare_mi_renumbered_but_grouping_equivalent_match() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![
        mi_record(b"read1", 100, "1"),
        mi_record(b"read2", 200, "1"),
        mi_record(b"read3", 300, "2"),
    ];
    let records2 = vec![
        mi_record(b"read1", 100, "5"),
        mi_record(b"read2", 200, "5"),
        mi_record(b"read3", 300, "9"),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let cfg = keyjoin_cfg(&tmp);
    let outcome: KeyJoinOutcome =
        keyjoin_compare(&bam1, &bam2, &cfg).expect("keyjoin_compare should succeed");

    assert_eq!(outcome.content_diffs, 0, "MI renumbering must not be reported as a content diff");
    assert_eq!(outcome.mi_bijection_mismatches, 0, "the grouping partition is unchanged");
    assert!(
        outcome.is_match(),
        "MI-renumbered but equivalently-grouped BAMs must match: {outcome:?}"
    );
}

/// (c) BS1 regression proof: a non-MI content bug (SEQ mutated) on one matched pair,
/// with the MI partition completely untouched, must still DIFFER. Under the
/// pre-hardening `execute_grouping_unordered`, only MI equivalence was checked, so this
/// same input falsely reported EQUIVALENT — this is the load-bearing proof that the
/// key-join engine now actually checks content.
#[test]
fn test_keyjoin_compare_content_diff_with_intact_grouping_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![
        mapped_record_with_mi_and_seq(b"read1", 100, "1", b"ACGTACGT"),
        mapped_record_with_mi_and_seq(b"read2", 200, "1", b"ACGTACGT"),
    ];
    // Same grouping (both MI=1 on both sides), but read2's SEQ is mutated on bam2.
    let records2 = vec![
        mapped_record_with_mi_and_seq(b"read1", 100, "1", b"ACGTACGT"),
        mapped_record_with_mi_and_seq(b"read2", 200, "1", b"ACGTACGA"),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let cfg = keyjoin_cfg(&tmp);
    let outcome: KeyJoinOutcome =
        keyjoin_compare(&bam1, &bam2, &cfg).expect("keyjoin_compare should succeed");

    assert_eq!(
        outcome.mi_bijection_mismatches, 0,
        "the MI grouping itself is untouched by this mutation"
    );
    assert_eq!(outcome.content_diffs, 1, "the mutated SEQ must be reported as a content diff");
    assert!(
        !outcome.is_match(),
        "a content diff on a matched pair must DIFFER even with intact grouping: {outcome:?}"
    );
}

/// (d) A real grouping divergence: `read1` and `read2` share one molecule (MI=1) in
/// bam1, but are split into two distinct molecules (MI=1 and MI=2) in bam2. Content is
/// otherwise identical, so this must surface purely as an MI-bijection violation.
#[test]
fn test_keyjoin_compare_mi_bijection_violated_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];
    let records2 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "2")];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let cfg = keyjoin_cfg(&tmp);
    let outcome: KeyJoinOutcome =
        keyjoin_compare(&bam1, &bam2, &cfg).expect("keyjoin_compare should succeed");

    assert_eq!(outcome.content_diffs, 0, "MI is excluded from content comparison");
    assert!(
        outcome.mi_bijection_mismatches > 0,
        "a molecule split across two bam2 MIs must be flagged: {outcome:?}"
    );
    assert!(!outcome.is_match(), "a bijection violation must DIFFER: {outcome:?}");
}

/// (e) A record present in only one file: the merge-join must report a presence DIFFER
/// (not silently drop the record or resync past it).
#[test]
fn test_keyjoin_compare_record_dropped_in_bam2_presence_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];
    let records2 = vec![mi_record(b"read1", 100, "1")]; // read2 dropped

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let cfg = keyjoin_cfg(&tmp);
    let outcome: KeyJoinOutcome =
        keyjoin_compare(&bam1, &bam2, &cfg).expect("keyjoin_compare should succeed");

    assert_eq!(outcome.bam1_count, 2);
    assert_eq!(outcome.bam2_count, 1);
    assert_eq!(outcome.only_in_bam1, 1);
    assert_eq!(outcome.only_in_bam2, 0);
    assert!(!outcome.is_match(), "a record dropped on one side must DIFFER: {outcome:?}");
}

/// R2 header-comparison gap, wired into the key-join engine: `keyjoin_compare` internally
/// canonicalizes both inputs to queryname order via a rewritten `@HD`, so the header check
/// must compare the *original* inputs, not the canonicalized scratch copies — otherwise a
/// genuine `@RG` divergence (e.g. two different samples fed into `group`) would be masked
/// by the canonicalization always emitting the same `@HD`. Records and MI groupings here
/// are otherwise identical, isolating the header check as the sole source of the DIFFER.
#[test]
fn test_keyjoin_compare_header_rg_sample_mismatch_differs() {
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReadGroup;
    use noodles::sam::header::record::value::map::read_group::tag as rg_tag;

    let tmp = TempDir::new().unwrap();

    let mut rg1 = Map::<ReadGroup>::default();
    rg1.other_fields_mut().insert(rg_tag::SAMPLE, bstr::BString::from("sampleA"));
    let mut header1 = create_minimal_header("chr1", 10000);
    header1.read_groups_mut().insert(bstr::BString::from("rg1"), rg1);

    let mut rg2 = Map::<ReadGroup>::default();
    rg2.other_fields_mut().insert(rg_tag::SAMPLE, bstr::BString::from("sampleB"));
    let mut header2 = create_minimal_header("chr1", 10000);
    header2.read_groups_mut().insert(bstr::BString::from("rg1"), rg2);

    let records = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header1, &records);
    write_bam(&bam2, &header2, &records);

    let cfg = keyjoin_cfg(&tmp);
    let outcome: KeyJoinOutcome =
        keyjoin_compare(&bam1, &bam2, &cfg).expect("keyjoin_compare should succeed");

    assert_eq!(outcome.content_diffs, 0, "records themselves are identical");
    assert_eq!(outcome.mi_bijection_mismatches, 0, "MI grouping itself is untouched");
    assert!(outcome.header_mismatch, "an @RG SM divergence must be flagged: {outcome:?}");
    assert!(
        !outcome.is_match(),
        "an @RG divergence must DIFFER even with matching content/grouping: {outcome:?}"
    );
}

/// (f) F1 soundness guard: a single template with all three `RecordKey` segment
/// classifications — `FIRST_OF_PAIR` only, `LAST_OF_PAIR` only, and (pathologically)
/// *both* set — must hard-fail rather than silently desync.
///
/// `fgumi_sort`'s `RawQuerynameLexKey` (the key type behind
/// `QuerynameComparator::Lexicographic`, which the canonicalization sort uses) orders by
/// name and then by the *full* packed flag word `((flags & 0xc0) << 8) | ((flags & 0x100)
/// << 3) | ((flags & 0x800) >> 3)` (`queryname_flag_order`), so a record with both R1 and
/// R2 bits set (`0xc0`) sorts *after* a plain R2 record (`0x80`) in the canonicalized
/// physical order. But `record_key::record_key` maps `(true, true)` to
/// `Segment::Fragment`, which sorts *first* under `RecordKey`'s derived `Ord`. The two
/// orderings disagree for this record, so the canonicalized stream is not actually
/// non-decreasing under the key the merge-join advances on — exactly the case the F1
/// guard exists to catch as a loud failure instead of a silent miscompare.
#[test]
fn test_keyjoin_compare_hard_errors_on_record_key_order_violation() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records = vec![
        record_with_flags(b"read1", flags::FIRST_SEGMENT),
        record_with_flags(b"read1", flags::LAST_SEGMENT),
        record_with_flags(b"read1", flags::FIRST_SEGMENT | flags::LAST_SEGMENT),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let cfg = keyjoin_cfg(&tmp);
    let err = keyjoin_compare(&bam1, &bam2, &cfg)
        .expect_err("a RecordKey-order-violating stream must hard-fail, not silently desync");
    assert!(
        format!("{err:#}").contains("RecordKey ordering"),
        "error should name the violated invariant, got: {err:#}"
    );
}

// ---------------------------------------------------------------------------
// `--command group` end-to-end tests (Task 3.4: wiring the key-join engine onto the
// real CLI preset, not just calling `keyjoin_compare` directly as the tests above do).
//
// These prove: (1) the `group` preset now runs through the key-join engine and still
// tolerates a content-identical, MI-renumbered-and-reordered pair of BAMs; (2) BS1 is
// closed end to end via the actual `--command group` command path, not just at the
// engine level; and (3) the pre-hardening stdout tokens (`EQUIVALENT`, `DIFFER`,
// `Missing MI in BAM1: {n}`/`Missing MI in BAM2: {n}`) that other tooling greps still
// print unchanged.
// ---------------------------------------------------------------------------

#[test]
fn test_command_group_content_identical_mi_renumbered_and_reordered_matches() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let r1 = mi_record(b"read1", 100, "1");
    let r2 = mi_record(b"read2", 200, "1");
    let r3 = mi_record(b"read3", 300, "2");
    let records1 = vec![r1, r2, r3];

    // Same content, reordered, and MI values renumbered (fgumi's 1/1/2 relabeled as
    // fgbio's 9/5/5) — exactly what a cross-tool `group` comparison must tolerate.
    let records2 = vec![
        mi_record(b"read3", 300, "9"),
        mi_record(b"read1", 100, "5"),
        mi_record(b"read2", 200, "5"),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let scratch = tmp.path().join("scratch");
    fs::create_dir_all(&scratch).expect("create scratch dir");
    let (success, stdout) =
        run_compare_command(&bam1, &bam2, "group", &["--sort-tmp-dir", scratch.to_str().unwrap()]);

    assert!(
        success,
        "Expected EQUIVALENT for content-identical, MI-renumbered, reordered BAMs via \
         `--command group`, stdout:\n{stdout}"
    );
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT in output, got:\n{stdout}");
    assert!(
        stdout.contains("Missing MI in BAM1: 0") && stdout.contains("Missing MI in BAM2: 0"),
        "Expected zero missing-MI counts, got:\n{stdout}"
    );
}

/// BS1 regression proof, end to end: a non-MI content bug (SEQ mutated on one matched
/// pair) with the MI partition completely untouched must DIFFER via the real
/// `--command group` path, not just via a direct `keyjoin_compare` call. Under the
/// pre-hardening `execute_grouping_unordered` (MI-equivalence-only), this exact input
/// falsely reported EQUIVALENT.
#[test]
fn test_command_group_content_diff_with_intact_grouping_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![
        mapped_record_with_mi_and_seq(b"read1", 100, "1", b"ACGTACGT"),
        mapped_record_with_mi_and_seq(b"read2", 200, "1", b"ACGTACGT"),
    ];
    // Same grouping (both MI=1 on both sides), but read2's SEQ is mutated on bam2.
    let records2 = vec![
        mapped_record_with_mi_and_seq(b"read1", 100, "1", b"ACGTACGT"),
        mapped_record_with_mi_and_seq(b"read2", 200, "1", b"ACGTACGA"),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let scratch = tmp.path().join("scratch");
    fs::create_dir_all(&scratch).expect("create scratch dir");
    let (success, stdout) =
        run_compare_command(&bam1, &bam2, "group", &["--sort-tmp-dir", scratch.to_str().unwrap()]);

    assert!(
        !success,
        "Expected DIFFER for a content bug with intact MI grouping via `--command group` \
         (BS1 must be closed end-to-end), stdout:\n{stdout}"
    );
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
    assert!(
        stdout.contains("Content diffs (excluding MI): 1"),
        "Expected exactly one content diff, got:\n{stdout}"
    );
    assert!(
        stdout.contains("MI bijection mismatches: 0"),
        "the MI grouping itself is untouched by this mutation, got:\n{stdout}"
    );
}

/// Preserved-token proof: `--command group` still reports non-zero `Missing MI in BAM1:
/// {n}`/`Missing MI in BAM2: {n}` (and DIFFERs) when both BAMs are entirely ungrouped —
/// this is orthogonal to content and must not be silently accepted just because content
/// happens to match.
#[test]
fn test_command_group_missing_mi_in_both_bams_differs_with_preserved_tokens() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![mapped_record(b"read1", 100), mapped_record(b"read2", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let scratch = tmp.path().join("scratch");
    fs::create_dir_all(&scratch).expect("create scratch dir");
    let (success, stdout) =
        run_compare_command(&bam1, &bam2, "group", &["--sort-tmp-dir", scratch.to_str().unwrap()]);

    assert!(
        !success,
        "Expected DIFFER when MI tags are missing in both BAMs via `--command group`, \
         stdout:\n{stdout}"
    );
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
    assert!(
        stdout.contains("Missing MI in BAM1: 2") && stdout.contains("Missing MI in BAM2: 2"),
        "Expected missing-MI counts of 2 each, got:\n{stdout}"
    );
}

// ============================================================================
// `sort_verify_compare` / `--command sort` tests (Phase 4)
// ============================================================================
//
// `sort` order *is* the payload, but conforming sort implementations may legitimately
// break ties differently (coordinate: records equal on (tid, pos, reverse) are an
// unordered set; template-coordinate: fgumi tie-breaks the name lane with a hash where
// samtools/fgbio use lexical order — the documented SORT-01 residue). These tests prove
// the sort-verify engine tolerates intra-run tie reordering while still catching a
// genuine mis-sort, a missing/extra record within a tied run, or a content difference.

/// Builds an unpaired (fragment) record at `(ref_id, pos)`, optionally reverse strand.
/// Deliberately unpaired: `extract_template_key_inline`'s pairing logic (mate
/// unclipped-position lookup, upper/lower canonicalization) isn't what these tests
/// exercise, so fragment reads keep the fixtures focused on the tie-break behavior under
/// test — for both coordinate and template-coordinate order, a fragment read's core
/// template key reduces to `(tid, pos, reverse)`, mirroring the coordinate key exactly.
fn fragment_record(name: &[u8], ref_id: i32, pos: i32, reverse: bool) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name)
        .sequence(b"ACGTACGT")
        .qualities(&[30; 8])
        .flags(if reverse { flags::REVERSE } else { 0 })
        .ref_id(ref_id)
        .pos(pos - 1) // 1-based test convention -> 0-based BAM
        .mapq(60);
    b.build()
}

/// Builds a queryname-sorted SAM header (`@HD SO:queryname`), optionally declaring an
/// `SS` sub-sort tag to select `sort_verify_compare`'s `QuerynameComparator`. `ss: None`
/// matches `detect_sort_order`'s default when the tag is absent
/// (`QuerynameComparator::Lexicographic`, fgumi's current bare-SO writer output);
/// `ss: Some("natural")` selects `QuerynameComparator::Natural`.
fn create_queryname_sorted_header(ref_name: &str, ref_len: usize, ss: Option<&str>) -> Header {
    use bstr::BString;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::{
        Map, map::Header as HeaderRecord, map::ReferenceSequence,
    };
    use std::num::NonZeroUsize;

    let HeaderTag::Other(tag) = HeaderTag::from(*b"SO") else { unreachable!() };
    let mut builder = HeaderRecordMap::<HeaderRecord>::builder().insert(tag, "queryname");
    if let Some(ss_value) = ss {
        let HeaderTag::Other(tag) = HeaderTag::from(*b"SS") else { unreachable!() };
        builder = builder.insert(tag, ss_value);
    }
    let header_map = builder.build().expect("valid header map");

    let reference_sequence = Map::<ReferenceSequence>::new(
        NonZeroUsize::new(ref_len).expect("reference length must be non-zero"),
    );

    Header::builder()
        .set_header(header_map)
        .add_reference_sequence(BString::from(ref_name), reference_sequence)
        .build()
}

/// Coordinate, template-coordinate, and queryname sort orders each tolerate reordering
/// two records that tie on the full sort key (an "equal-key run" of size 2) — the
/// SORT-01 residue these three near-identical cases each guard against for a different
/// `SortOrder`. Not every original standalone test asserted `sort_order`/record counts,
/// so those checks are optional per case (`None` skips the assertion) to preserve each
/// case's exact original coverage.
#[rstest]
#[case::coordinate(
    create_coordinate_sorted_header("chr1", 10000),
    // readA and readB tie on (tid=0, pos=100, reverse=false): coordinate sort has no
    // tie-break lane beyond that triple, so either relative order is conforming output.
    fragment_record(b"readA", 0, 100, false),
    fragment_record(b"readB", 0, 100, false),
    None,
    Some((2, 2))
)]
#[case::template_coordinate(
    // SO:unsorted GO:query SS:template-coordinate (fgumi's current bare-form writer output).
    create_minimal_header("chr1", 10000),
    // Two fragment reads at the same (tid=0, pos=100, reverse=false) tie on every
    // TemplateKey core lane (tid1/pos1/neg1 equal; tid2/pos2/neg2 both MAX/false for
    // unpaired reads; no CB/library/MI). Only the name-hash tie-break lane (dropped by
    // `core_cmp`) would order readA vs readB — simulating the fgumi-hash-vs-samtools-
    // lexical SORT-01 residue by simply swapping their order between the two files.
    fragment_record(b"readA", 0, 100, false),
    fragment_record(b"readB", 0, 100, false),
    Some(fgumi_sort::SortOrder::TemplateCoordinate),
    None
)]
#[case::queryname(
    // No `SS` tag -> `detect_sort_order` defaults to `QuernameComparator::Lexicographic`
    // (fgumi's current bare-`SO` writer output), the comparator this engine picks by default.
    create_queryname_sorted_header("chr1", 10000, None),
    // r_a and r_b share the same queryname AND the same segment/record-type flag lane
    // (`queryname_flag_order` only encodes read1/read2 + secondary + supplementary, and
    // both records here are unpaired, non-secondary, non-supplementary alignments), so
    // they tie on the full queryname sort key even though they are genuinely different
    // records (different pos/strand) — an equal-key run of size 2, exactly like the
    // coordinate and template-coordinate ties above.
    fragment_record(b"read1", 0, 100, false),
    fragment_record(b"read1", 0, 50, true),
    Some(fgumi_sort::SortOrder::Queryname(fgumi_sort::QuerynameComparator::Lexicographic)),
    Some((2, 2))
)]
fn test_sort_verify_equal_key_run_tolerates_reorder(
    #[case] header: Header,
    #[case] r_a: RawRecord,
    #[case] r_b: RawRecord,
    #[case] expected_sort_order: Option<fgumi_sort::SortOrder>,
    #[case] expected_counts: Option<(u64, u64)>,
) {
    let tmp = TempDir::new().unwrap();
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[r_a.clone(), r_b.clone()]);
    write_bam(&bam2, &header, &[r_b, r_a]); // swapped order

    let outcome =
        sort_verify_compare(&bam1, &bam2, 100).expect("sort_verify_compare should succeed");

    if let Some(sort_order) = expected_sort_order {
        assert_eq!(outcome.sort_order, sort_order);
    }
    if let Some((bam1_count, bam2_count)) = expected_counts {
        assert_eq!(outcome.bam1_count, bam1_count);
        assert_eq!(outcome.bam2_count, bam2_count);
    }
    assert_eq!(outcome.bam1_violations, 0);
    assert_eq!(outcome.bam2_violations, 0);
    assert_eq!(outcome.run_mismatches, 0);
    assert!(
        outcome.is_match(),
        "tied records in a different order between files must match: {outcome:?}"
    );
}

/// Coordinate and queryname sort orders each catch a genuine mis-sort in `bam1` while
/// `bam2` remains correctly sorted.
#[rstest]
#[case::coordinate(
    create_coordinate_sorted_header("chr1", 10000),
    // bam1 is genuinely out of coordinate order (pos 200 appears before pos 100).
    vec![fragment_record(b"read1", 0, 200, false), fragment_record(b"read2", 0, 100, false)],
    // bam2 is correctly sorted.
    vec![fragment_record(b"read2", 0, 100, false), fragment_record(b"read1", 0, 200, false)]
)]
#[case::queryname(
    create_queryname_sorted_header("chr1", 10000, None),
    // Lexicographic order: "read1" < "read2" (byte-wise). bam1 is genuinely out of
    // queryname order (read2 appears before read1).
    vec![fragment_record(b"read2", 0, 100, false), fragment_record(b"read1", 0, 200, false)],
    // bam2 is correctly sorted.
    vec![fragment_record(b"read1", 0, 200, false), fragment_record(b"read2", 0, 100, false)]
)]
fn test_sort_verify_mis_sorted_bam1_differs(
    #[case] header: Header,
    #[case] records1: Vec<RawRecord>,
    #[case] records2: Vec<RawRecord>,
) {
    let tmp = TempDir::new().unwrap();
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        sort_verify_compare(&bam1, &bam2, 100).expect("sort_verify_compare should succeed");

    assert!(outcome.bam1_violations > 0, "bam1 is genuinely mis-sorted: {outcome:?}");
    assert_eq!(outcome.bam2_violations, 0, "bam2 is correctly sorted");
    assert!(!outcome.is_match(), "a genuine mis-sort must DIFFER: {outcome:?}");
}

#[test]
fn test_sort_verify_coordinate_missing_record_in_equal_key_run_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);

    // Both files are individually correctly coordinate-sorted, but the equal-key run at
    // (tid=0, pos=100) has 2 records in bam1 and only 1 in bam2 — a dropped record must
    // still DIFFER even though neither file violates its own sort order.
    let r_a = fragment_record(b"readA", 0, 100, false);
    let r_b = fragment_record(b"readB", 0, 100, false);
    let r_next = fragment_record(b"readC", 0, 200, false);

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[r_a.clone(), r_b, r_next.clone()]);
    write_bam(&bam2, &header, &[r_a, r_next]);

    let outcome =
        sort_verify_compare(&bam1, &bam2, 100).expect("sort_verify_compare should succeed");

    assert_eq!(outcome.bam1_violations, 0, "bam1 is itself correctly sorted");
    assert_eq!(outcome.bam2_violations, 0, "bam2 is itself correctly sorted");
    assert!(
        outcome.run_mismatches > 0,
        "a record dropped from an equal-key run must be caught: {outcome:?}"
    );
    assert!(!outcome.is_match(), "a missing record in a tied run must DIFFER: {outcome:?}");
}

/// Bonus coverage for the other `QuernameComparator` variant: `SS:natural` selects
/// `QuernameComparator::Natural`, whose numeric-run handling ("read2" < "read10") is
/// exactly where it disagrees with lexicographic byte order ("read10" < "read2"). A file
/// ordered `[read2, read10]` is natural-sorted but would violate lexicographic order,
/// so a clean (zero-violation) result here confirms `detect_sort_order` picked the
/// natural comparator from the `SS:natural` tag rather than silently defaulting to
/// lexicographic.
#[test]
fn test_sort_verify_queryname_natural_order_detected_from_ss_tag() {
    let tmp = TempDir::new().unwrap();
    let header = create_queryname_sorted_header("chr1", 10000, Some("natural"));

    let records =
        vec![fragment_record(b"read2", 0, 100, false), fragment_record(b"read10", 0, 200, false)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let outcome =
        sort_verify_compare(&bam1, &bam2, 100).expect("sort_verify_compare should succeed");

    assert_eq!(
        outcome.sort_order,
        fgumi_sort::SortOrder::Queryname(fgumi_sort::QuerynameComparator::Natural)
    );
    assert_eq!(outcome.bam1_violations, 0, "read2 < read10 under natural order: {outcome:?}");
    assert_eq!(outcome.bam2_violations, 0);
    assert!(outcome.is_match(), "{outcome:?}");
}

/// Coordinate and queryname sort orders each catch a genuine content difference (a single
/// base flipped) on the sole record in an equal-key run of size 1, even though neither
/// file violates its own sort order.
#[rstest]
#[case::coordinate(create_coordinate_sorted_header("chr1", 10000))]
#[case::queryname(create_queryname_sorted_header("chr1", 10000, None))]
fn test_sort_verify_content_diff_differs(#[case] header: Header) {
    let tmp = TempDir::new().unwrap();

    let r1_bam1 = fragment_record(b"read1", 0, 100, false);
    let r1_bam2 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGA") // last base differs from r1_bam1
            .qualities(&[30; 8])
            .ref_id(0)
            .pos(99)
            .mapq(60);
        b.build()
    };

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[r1_bam1]);
    write_bam(&bam2, &header, &[r1_bam2]);

    let outcome =
        sort_verify_compare(&bam1, &bam2, 100).expect("sort_verify_compare should succeed");

    assert_eq!(outcome.bam1_violations, 0);
    assert_eq!(outcome.bam2_violations, 0);
    assert!(
        outcome.run_mismatches > 0,
        "a content difference on the sole record in the run must be caught: {outcome:?}"
    );
    assert!(!outcome.is_match(), "a content difference must DIFFER: {outcome:?}");
}

/// R2 header-comparison gap, wired into the sort-verify engine: two files that agree on
/// sort order (so `detect_sort_order` succeeds and never bails) and on every record's
/// content can still DIFFER if their `@SQ` reference dictionaries disagree — a field
/// `detect_sort_order` never inspects.
#[test]
fn test_sort_verify_header_sq_length_mismatch_differs() {
    let tmp = TempDir::new().unwrap();
    let header1 = create_coordinate_sorted_header("chr1", 10000);
    let header2 = create_coordinate_sorted_header("chr1", 20000);

    let records = vec![fragment_record(b"read1", 0, 100, false)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header1, &records);
    write_bam(&bam2, &header2, &records);

    let outcome =
        sort_verify_compare(&bam1, &bam2, 100).expect("sort_verify_compare should succeed");

    assert_eq!(outcome.bam1_violations, 0);
    assert_eq!(outcome.bam2_violations, 0);
    assert_eq!(outcome.run_mismatches, 0, "record content itself is identical");
    assert!(outcome.header_mismatch, "a @SQ length mismatch must be flagged: {outcome:?}");
    assert!(
        !outcome.is_match(),
        "a @SQ divergence must DIFFER even with matching sort order/content: {outcome:?}"
    );
}

#[test]
fn test_sort_verify_mismatched_sort_orders_between_files_errors() {
    let tmp = TempDir::new().unwrap();
    let coord_header = create_coordinate_sorted_header("chr1", 10000);
    let template_header = create_minimal_header("chr1", 10000);

    let records = vec![fragment_record(b"read1", 0, 100, false)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &coord_header, &records);
    write_bam(&bam2, &template_header, &records);

    let err = sort_verify_compare(&bam1, &bam2, 100)
        .expect_err("mismatched declared sort orders must be a hard error, not a DIFFER");
    assert!(err.to_string().contains("different sort orders"), "unexpected error message: {err}");
}

/// End-to-end smoke test: `--command sort` on identical, correctly-sorted BAMs reports
/// IDENTICAL through the real CLI dispatch (`CommandPreset::Sort` → `execute_sort_verify`
/// → `sort_verify_compare`), not just the underlying engine function.
#[test]
fn test_command_sort_identical_bams_match() {
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let records =
        vec![fragment_record(b"read1", 0, 100, false), fragment_record(b"read2", 0, 200, false)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare_command(&bam1, &bam2, "sort", &[]);

    assert!(success, "identical sorted BAMs must be IDENTICAL via --command sort:\n{stdout}");
    assert!(stdout.contains("IDENTICAL"), "expected IDENTICAL in output, got:\n{stdout}");
}

/// End-to-end smoke test: `--command sort` reports DIFFER (and exits non-zero) when
/// bam1 is genuinely mis-sorted, through the real CLI dispatch.
#[test]
fn test_command_sort_mis_sorted_bam1_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let records1 =
        vec![fragment_record(b"read1", 0, 200, false), fragment_record(b"read2", 0, 100, false)];
    let records2 =
        vec![fragment_record(b"read2", 0, 100, false), fragment_record(b"read1", 0, 200, false)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (success, stdout) = run_compare_command(&bam1, &bam2, "sort", &[]);

    assert!(!success, "a genuine mis-sort must DIFFER via --command sort:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "expected DIFFER in output, got:\n{stdout}");
}

/// `--mode`/`--ignore-order` don't apply to the dedicated sort-verify engine and must be
/// rejected rather than silently ignored.
#[test]
fn test_command_sort_rejects_explicit_mode() {
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let records = vec![fragment_record(b"read1", 0, 100, false)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare_command(&bam1, &bam2, "sort", &["--mode", "content"]);
    assert!(!success, "--mode with --command sort should be rejected, got:\n{stdout}");
}

/// Mirrors `test_command_sort_rejects_explicit_mode` for `--ignore-order`: neither
/// concept applies to the dedicated sort-verify engine, so both must be rejected
/// explicitly rather than one being silently accepted while the other errors.
#[test]
fn test_command_sort_rejects_explicit_ignore_order() {
    let tmp = TempDir::new().unwrap();
    let header = create_coordinate_sorted_header("chr1", 10000);
    let records = vec![fragment_record(b"read1", 0, 100, false)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (success, stdout) = run_compare_command(&bam1, &bam2, "sort", &["--ignore-order"]);
    assert!(!success, "--ignore-order with --command sort should be rejected, got:\n{stdout}");
}
