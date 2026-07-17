//! Integration tests for the `fgumi compare bams` command.
//!
//! These tests exercise both compare modes (content, grouping)
//! using the raw byte comparison path.

use clap::Parser;
use fgumi_lib::commands::command::Command as FgumiCommand;
use fgumi_lib::commands::compare::{
    CompareBams, CompareMismatch, ContentPredicate, molecule_join_compare, positional_compare,
    sort_verify_compare,
};
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::sam::Header;
use rstest::rstest;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_coordinate_sorted_header, create_minimal_header, mi_record, write_bam,
};

/// Runs `fgumi compare bams` via subprocess and returns `(exit_code, stdout, stderr)`.
///
/// Returns the *exact* exit code (not just success/failure) so a negative test can tell a
/// clean DIFFER (`Some(1)`, set explicitly for `CompareMismatch` in `main`) apart from a
/// panic/abort (`Some(101)`/`Some(134)`), a usage error (`Some(2)`), or a signal (`None`) —
/// reducing all of those to `!success()` would let an unrelated crash satisfy a
/// "must-not-match" assertion. `stderr` is returned so an argument-rejection test can assert
/// the specific diagnostic (a generic `anyhow` runtime error also exits 1) rather than
/// merely "did not succeed". Used by tests that assert on specific substrings in the
/// diagnostic report (e.g. EQUIVALENT vs IDENTICAL, mismatch counts); for tests that only
/// need to know whether the BAMs matched, prefer `run_compare_in_process`.
fn run_compare(
    bam1: &Path,
    bam2: &Path,
    mode: &str,
    extra_args: &[&str],
) -> (Option<i32>, String, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["compare", "bams"])
        .arg(bam1)
        .arg(bam2)
        .args(["--mode", mode])
        .args(extra_args)
        .output()
        .expect("Failed to run fgumi");

    (
        output.status.code(),
        String::from_utf8_lossy(&output.stdout).to_string(),
        String::from_utf8_lossy(&output.stderr).to_string(),
    )
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

/// Runs `fgumi compare bams --command <command>` via subprocess and returns
/// `(exit_code, stdout, stderr)` (see [`run_compare`] for why the exact code and stderr are
/// returned rather than a bare `success()`).
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
) -> (Option<i32>, String, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["compare", "bams"])
        .arg(bam1)
        .arg(bam2)
        .args(["--command", command])
        .args(extra_args)
        .output()
        .expect("Failed to run fgumi");

    (
        output.status.code(),
        String::from_utf8_lossy(&output.stdout).to_string(),
        String::from_utf8_lossy(&output.stderr).to_string(),
    )
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
// Content mode: MI treated as an ordinary tag
// ---------------------------------------------------------------------------

#[test]
fn test_content_mode_different_mi_tags() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];
    let records2 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "2")];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    // Content mode compares MI as an ordinary tag value (exact), so read2's
    // MI:Z:1 vs MI:Z:2 is a tag difference and the BAMs DIFFER.
    assert!(
        !run_compare_in_process(&bam1, &bam2, "content", &[]),
        "Expected mismatch for different MI tags in content mode"
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

    let (code, stdout, _stderr) = run_compare(&bam1, &bam2, "grouping", &[]);
    assert_eq!(code, Some(0), "Expected success for identical grouped BAMs, stdout:\n{stdout}");
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT in output, got:\n{stdout}");
}

#[test]
fn test_grouping_mode_flag_mismatch_is_presence_differ() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // BAM1: read1 is R1, read2 is R2. BAM2: read1 is R2, read2 is R1 (swapped). Both
    // records share MI=1 on both sides, so each file has exactly one molecule (canon
    // "read1") that matches its counterpart by canonical id; `compare_molecule` keys
    // membership on RecordKey (which encodes the R1/R2 segment), so the flag swap makes
    // every member unmatchable within that molecule — the mismatch surfaces as membership
    // diffs, and the BAMs DIFFER.
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

    let (code, stdout, _stderr) = run_compare(&bam1, &bam2, "grouping", &[]);
    assert_eq!(code, Some(1), "Expected DIFFER (exit 1) for flag mismatches, stdout:\n{stdout}");
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
    assert!(
        stdout.contains("Molecules matched: 0"),
        "Expected the sole molecule to fail matching, got:\n{stdout}"
    );
}

// ---------------------------------------------------------------------------
// Reordered tags tests (content mode)
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
    let (code, stdout, _stderr) = run_compare(&bam1, &bam2, "content", &[]);
    assert_eq!(
        code,
        Some(0),
        "Expected success for reordered tags in content mode, stdout:\n{stdout}"
    );
    assert!(stdout.contains("IDENTICAL"), "Expected IDENTICAL in output, got:\n{stdout}");
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

    let (code, stdout, _stderr) = run_compare(&bam1, &bam2, "grouping", &["--ignore-order"]);
    assert_eq!(code, Some(0), "Expected success for ignore-order grouping mode, stdout:\n{stdout}");
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT in output, got:\n{stdout}");
}

/// Runs `CompareBams::execute()` in-process under `--mode grouping` (plus any `extra_args`,
/// e.g. `--ignore-order`) and returns the error's full display string. Panics if `execute`
/// succeeds, or if it fails with a [`CompareMismatch`] (a genuine record-level DIFFER) rather
/// than the fully-MI-less hard-fail this helper exists to capture (Finding A, Task 8) —
/// mirrors `run_compare_content_expect_err`'s success/DIFFER-vs-hard-error distinction for
/// the grouping-mode guard.
fn run_compare_grouping_expect_err(bam1: &Path, bam2: &Path, extra_args: &[&str]) -> String {
    let mut args: Vec<&str> =
        vec!["bams", bam1.to_str().unwrap(), bam2.to_str().unwrap(), "--mode", "grouping"];
    args.extend(extra_args);
    let cmd = CompareBams::try_parse_from(args).expect("failed to parse compare bams args");
    match cmd.execute("fgumi compare bams") {
        Ok(()) => panic!("expected compare bams to hard-fail on fully-MI-less input"),
        Err(e) if e.is::<CompareMismatch>() => {
            panic!("expected the fully-MI-less guard error, got a CompareMismatch (DIFFER): {e:#}")
        }
        Err(e) => format!("{e:#}"),
    }
}

// ---------------------------------------------------------------------------
// Missing MI tag behavior (streaming molecule-join engine)
//
// The old (retired) key-join engine explicitly counted "missing MI" records and always
// DIFFERed when both files entirely lacked MI tags, on the theory that "no grouping has
// been verified" must never read as EQUIVALENT. The streaming molecule-join engine
// (`molecule_runs`) has no such per-record check: consecutive records with no MI tag
// (`base == None`) are simply grouped into one molecule, like any other base-MI value —
// which, left unguarded, let two *content-identical* fully-MI-less BAMs MATCH (Task 7
// finding). `molecule_join_compare` closes that gap with an explicit whole-input guard
// (Finding A, Task 8): if **neither** side ever saw an MI tag on any record, the compare
// hard-fails rather than reporting MATCH or DIFFER — see the `saw_mi1`/`saw_mi2` tracking
// and the `bail!` in `molecule_join_compare`. This is a deliberate, user-facing
// oracle-semantics change: `--mode grouping`/`--command group` require already-grouped
// (same-MI-consecutive) input, and feeding it fully-ungrouped input is now a hard error
// instead of silently reporting MATCH. A *partial*-MI-less pair (one side grouped, the
// other not) is untouched by this guard and still DIFFERs via ordinary molecule-count/
// membership asymmetry (see `test_grouping_partial_mi_less_pair_still_differs` below).
// ---------------------------------------------------------------------------

#[test]
fn test_grouping_mode_content_identical_but_mi_missing_in_both_bams_errors() {
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

    // See this section's header comment / Finding A: fully-MI-less input on both sides is
    // now a hard error, not a MATCH.
    let err = run_compare_grouping_expect_err(&bam1, &bam2, &[]);
    assert!(err.contains("MI-tagged"), "got: {err}");
}

#[test]
fn test_grouping_unordered_content_identical_but_mi_missing_in_both_bams_errors() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![mapped_record(b"read1", 100), mapped_record(b"read2", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    // See this section's header comment / Finding A: `--ignore-order` doesn't change the
    // fully-MI-less guard outcome.
    let err = run_compare_grouping_expect_err(&bam1, &bam2, &["--ignore-order"]);
    assert!(err.contains("MI-tagged"), "got: {err}");
}

/// Finding A must NOT touch the partial-MI-less case: one side grouped, the other not, is
/// still an ordinary DIFFER (via molecule-count asymmetry — see `molecule_runs`' doc comment
/// on how consecutive no-MI records fold into one molecule), not the whole-input hard error.
#[test]
fn test_grouping_partial_mi_less_pair_still_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // bam1: two distinct MI-tagged molecules (2 runs).
    let grouped = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "2")];
    // bam2: same two records, but with NO MI tag at all -> molecule_runs folds them into
    // one run (consecutive base == None), so bam1_molecules (2) != bam2_molecules (1).
    let ungrouped = vec![mapped_record(b"read1", 100), mapped_record(b"read2", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &grouped);
    write_bam(&bam2, &header, &ungrouped);

    let outcome =
        molecule_join_compare(&bam1, &bam2, 10).expect("molecule_join_compare should succeed");
    assert_eq!(outcome.bam1_molecules, 2);
    assert_eq!(outcome.bam2_molecules, 1);
    assert!(!outcome.is_match(), "a partial-MI-less pair must DIFFER, not error: {outcome:?}");
}

/// Review fix (empty-vs-empty exemption): two EMPTY grouped BAMs (zero records on both
/// sides) must MATCH, not hit the fully-MI-less guard. Neither side ever sees an MI tag —
/// there are no records at all to see one on — but there is also nothing to compare, so this
/// is a vacuous MATCH, not the "misuse" case (non-empty, fully MI-less input) the guard
/// exists to catch. See `molecule_join_compare`'s guard condition, which only fires when at
/// least one side actually produced molecules.
#[test]
fn test_grouping_mode_both_empty_bams_match() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[]);
    write_bam(&bam2, &header, &[]);

    let outcome = molecule_join_compare(&bam1, &bam2, 10).expect("empty-vs-empty must not error");
    assert_eq!(outcome.bam1_molecules, 0);
    assert_eq!(outcome.bam2_molecules, 0);
    assert!(outcome.is_match(), "two empty BAMs should vacuously MATCH: {outcome:?}");

    // Also confirm via the real CLI entry point (`CompareBams::execute`), not just the
    // engine function directly.
    assert!(
        run_compare_in_process(&bam1, &bam2, "grouping", &[]),
        "two empty BAMs should MATCH via --mode grouping"
    );
}

/// Review fix (guard-vs-diff precedence): both existing MI-less-guard tests above use
/// *identical* content on both sides, so the fully-MI-less guard firing was never
/// distinguished from an ordinary content DIFFER also firing. This test uses genuinely
/// different content (same read names/keys so the records would otherwise pair up and
/// surface a POS content diff — see `RecordKey`'s doc comment on POS not being part of a
/// primary record's identity) on two fully MI-less BAMs, and locks that the fully-MI-less
/// guard still wins: the compare hard-fails with the guard's message rather than reporting
/// an ordinary DIFFER.
#[test]
fn test_grouping_mode_mi_missing_and_content_differs_still_errors_via_guard() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // No MI tag on either side, but genuinely different content (different positions) so
    // that, absent the guard, this would be an ordinary content DIFFER, not a MATCH.
    let records1 = vec![mapped_record(b"read1", 100), mapped_record(b"read2", 200)];
    let records2 = vec![mapped_record(b"read1", 150), mapped_record(b"read2", 250)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let err = run_compare_grouping_expect_err(&bam1, &bam2, &[]);
    assert!(err.contains("MI-tagged"), "got: {err}");
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

    let (code, stdout, _stderr) = run_compare(&bam1, &bam2, "grouping", &[]);
    assert_eq!(
        code,
        Some(0),
        "Expected EQUIVALENT for identical paired-strand BAMs, stdout:\n{stdout}"
    );
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT, got:\n{stdout}");
    // Note: unlike the retired key-join engine, this molecule-join engine has no
    // "missing MI" counter to assert against directly. Regression coverage for MI:Z:0/A
    // and 0/B being parsed as present (not missing) lives primarily in
    // `test_grouping_mode_paired_strand_a_flipped_to_b_is_mismatch` below: if
    // `get_mi_tag_raw` regressed to treat these as unparseable, the strand-partition
    // check there would degenerate to two empty (trivially "equal") sets and the
    // expected DIFFER would be silently lost.
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

    let (code, stdout, _stderr) = run_compare(&bam1, &bam2, "grouping", &["--ignore-order"]);
    assert_eq!(
        code,
        Some(0),
        "Expected EQUIVALENT for ignore-order paired-strand BAMs, stdout:\n{stdout}"
    );
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT, got:\n{stdout}");
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

    let (code, stdout, _stderr) = run_compare(&bam1, &bam2, "grouping", &[]);
    assert_eq!(
        code,
        Some(1),
        "Expected DIFFER (exit 1) when /A and /B assignments disagree, stdout:\n{stdout}"
    );
    assert!(stdout.contains("DIFFER"), "Expected DIFFER, got:\n{stdout}");
    // read1 and read2 share one base-MI-0 molecule on both sides (canonicalized by min
    // read name), so this surfaces as a duplex strand-partition diff on that single
    // matched molecule, not a presence/membership diff.
    assert!(
        stdout.contains("duplex strand partition differs"),
        "Expected the failure to be a duplex strand-partition diff, got:\n{stdout}"
    );
    assert!(
        stdout.contains("Molecules matched: 0"),
        "the sole molecule must fail to match: {stdout}"
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

// ---------------------------------------------------------------------------
// molecule-join engine tests (streaming per-molecule hash-join for `group`)
//
// These tests call `molecule_join_compare` directly to verify its sound core: the
// hash-join tolerates molecule reordering and MI renumbering (the whole point of
// `group`'s cross-tool comparison) while still catching a genuine content diff, a
// genuine grouping divergence (membership split across molecules), and a
// presence-only molecule. Unlike the retired key-join engine, `molecule_join_compare`
// never re-sorts either input: both inputs must already be grouped (same-MI reads
// consecutive within the file), so every fixture below keeps each MI's records
// contiguous per file even when reordering molecules relative to each other or
// between the two files.
// ---------------------------------------------------------------------------

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

/// (a) Identical grouping, molecules in a different physical order between the two files
/// (each MI's own records stay contiguous within each file), same MI numbering on both
/// sides: molecule matching is by canonical id, not physical order, so this must MATCH.
#[test]
fn molecule_join_compare_reordered_molecules_same_mi_numbering_match() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let r1 = mi_record(b"read1", 100, "1");
    let r2 = mi_record(b"read2", 200, "1");
    let r3 = mi_record(b"read3", 300, "2");

    let records1 = vec![r1.clone(), r2.clone(), r3.clone()]; // molecule MI1={read1,read2}, MI2={read3}
    let records2 = vec![r3, r1, r2]; // same molecules, different physical order

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        molecule_join_compare(&bam1, &bam2, 10).expect("molecule_join_compare should succeed");

    assert_eq!(outcome.bam1_molecules, 2);
    assert_eq!(outcome.bam2_molecules, 2);
    assert_eq!(outcome.matched, 2);
    assert!(
        outcome.diff_details.is_empty(),
        "reordered, identically-grouped BAMs must have no diffs: {outcome:?}"
    );
    assert!(outcome.is_match(), "reordered, identically-grouped BAMs must match: {outcome:?}");
}

/// (b) Same grouping partition, different MI numbering (fgumi assigns `1`/`2`; fgbio
/// assigns `5`/`9` for the same two groups), content otherwise identical: MI renumbering
/// alone must not be reported as a diff — this is the key case `group` comparison exists
/// to tolerate.
#[test]
fn molecule_join_compare_mi_renumbered_but_grouping_equivalent_match() {
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

    let outcome =
        molecule_join_compare(&bam1, &bam2, 10).expect("molecule_join_compare should succeed");

    assert!(
        outcome.diff_details.is_empty(),
        "MI renumbering must not be reported as a diff: {outcome:?}"
    );
    assert!(
        outcome.is_match(),
        "MI-renumbered but equivalently-grouped BAMs must match: {outcome:?}"
    );
}

/// (c) BS1 regression proof: a non-MI content bug (SEQ mutated) on one matched pair,
/// with the MI partition completely untouched, must still DIFFER. A pre-hardening
/// MI-equivalence-only grouping check reported EQUIVALENT for this same input — this is
/// the load-bearing proof that the molecule-join engine now actually checks content.
#[test]
fn molecule_join_compare_content_diff_with_intact_grouping_differs() {
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

    let outcome =
        molecule_join_compare(&bam1, &bam2, 10).expect("molecule_join_compare should succeed");

    assert!(
        !outcome.diff_details.is_empty(),
        "the mutated SEQ must be reported as a diff: {outcome:?}"
    );
    assert!(
        !outcome.is_match(),
        "a content diff on a matched pair must DIFFER even with intact grouping: {outcome:?}"
    );
}

/// (d) A real grouping divergence: `read1` and `read2` share one molecule (MI=1) in
/// bam1, but are split into two distinct molecules (MI=1 and MI=2) in bam2. Content is
/// otherwise identical, so this must surface as a membership diff on the matched
/// molecule plus a residual molecule present only in bam2.
#[test]
fn molecule_join_compare_grouping_split_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];
    let records2 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "2")];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        molecule_join_compare(&bam1, &bam2, 10).expect("molecule_join_compare should succeed");

    assert!(
        !outcome.diff_details.is_empty(),
        "a molecule split across two bam2 MIs must be flagged: {outcome:?}"
    );
    assert!(!outcome.is_match(), "a grouping split must DIFFER: {outcome:?}");
}

/// (e) A record present in only one file: the hash-join must report a presence DIFFER
/// (not silently drop the record).
#[test]
fn molecule_join_compare_record_dropped_in_bam2_presence_differs() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let records1 = vec![mi_record(b"read1", 100, "1"), mi_record(b"read2", 200, "1")];
    let records2 = vec![mi_record(b"read1", 100, "1")]; // read2 dropped

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let outcome =
        molecule_join_compare(&bam1, &bam2, 10).expect("molecule_join_compare should succeed");

    assert_eq!(outcome.bam1_molecules, 1);
    assert_eq!(outcome.bam2_molecules, 1);
    assert!(
        !outcome.diff_details.is_empty(),
        "a record dropped on one side must produce a diff: {outcome:?}"
    );
    assert!(!outcome.is_match(), "a record dropped on one side must DIFFER: {outcome:?}");
}

// ---------------------------------------------------------------------------
// `--command group` end-to-end tests: wiring the streaming molecule-join engine onto
// the real CLI preset, not just calling `molecule_join_compare` directly as the tests
// above do.
//
// These prove: (1) the `group` preset runs through the molecule-join engine and still
// tolerates a content-identical, MI-renumbered-and-reordered pair of grouped BAMs; (2)
// BS1 (a non-MI content bug must DIFFER even with intact grouping) is closed end to end
// via the actual `--command group` command path, not just at the engine level; and (3)
// the pre-hardening stdout tokens (`EQUIVALENT`/`DIFFER`) that other tooling greps still
// print unchanged. `--sort-tmp-dir`/`--sort-memory` no longer exist (there is no
// canonicalization sort to configure), so these tests no longer pass them.
// ---------------------------------------------------------------------------

#[test]
fn test_command_group_content_identical_mi_renumbered_and_reordered_matches() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    let r1 = mi_record(b"read1", 100, "1");
    let r2 = mi_record(b"read2", 200, "1");
    let r3 = mi_record(b"read3", 300, "2");
    let records1 = vec![r1, r2, r3]; // molecule MI1={read1,read2}, MI2={read3}, each contiguous

    // Same content, molecules reordered relative to each other, and MI values renumbered
    // (fgumi's 1/1/2 relabeled as fgbio's 9/5/5) — exactly what a cross-tool `group`
    // comparison must tolerate. Each MI's own records stay contiguous within the file.
    let records2 = vec![
        mi_record(b"read3", 300, "9"),
        mi_record(b"read1", 100, "5"),
        mi_record(b"read2", 200, "5"),
    ];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records1);
    write_bam(&bam2, &header, &records2);

    let (code, stdout, _stderr) = run_compare_command(&bam1, &bam2, "group", &[]);

    assert_eq!(
        code,
        Some(0),
        "Expected EQUIVALENT for content-identical, MI-renumbered, reordered BAMs via \
         `--command group`, stdout:\n{stdout}"
    );
    assert!(stdout.contains("EQUIVALENT"), "Expected EQUIVALENT in output, got:\n{stdout}");
}

/// Regression: an explicit `--mode content` overrides a preset's own predicate. Under
/// `--command group` alone, an MI-only difference is a tolerated grouping renumbering
/// (EQUIVALENT); but `--command group --mode content` asks for exact positional content
/// comparison, so the same MI-only difference must DIFFER (the preset's `ExactMinusMi` must
/// not silently weaken an explicitly-requested content compare).
#[test]
fn test_command_group_explicit_content_mode_flags_mi_only_difference() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // Two BAMs identical except for the MI tag value on the single record.
    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &[mi_record(b"read1", 100, "1")]);
    write_bam(&bam2, &header, &[mi_record(b"read1", 100, "2")]);

    // `--command group` alone: the MI renumbering is a consistent bijection -> EQUIVALENT.
    let (group_code, group_stdout, _) = run_compare_command(&bam1, &bam2, "group", &[]);
    assert_eq!(
        group_code,
        Some(0),
        "`--command group` alone must tolerate an MI-only renumbering, stdout:\n{group_stdout}"
    );

    // `--command group --mode content`: explicit content mode forces the exact predicate,
    // so the MI-only difference is now a real DIFFER.
    let (content_code, content_stdout, _) =
        run_compare_command(&bam1, &bam2, "group", &["--mode", "content"]);
    assert_eq!(
        content_code,
        Some(1),
        "explicit `--mode content` must DIFFER on an MI-only diff (exact predicate), \
         stdout:\n{content_stdout}"
    );
    assert!(content_stdout.contains("DIFFER"), "expected DIFFER in output, got:\n{content_stdout}");
}

/// BS1 regression proof, end to end: a non-MI content bug (SEQ mutated on one matched
/// pair) with the MI partition completely untouched must DIFFER via the real
/// `--command group` path, not just via a direct `molecule_join_compare` call. A
/// pre-hardening MI-equivalence-only grouping check reported EQUIVALENT for this exact
/// input.
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

    let (code, stdout, _stderr) = run_compare_command(&bam1, &bam2, "group", &[]);

    assert_eq!(
        code,
        Some(1),
        "Expected DIFFER (exit 1) for a content bug with intact MI grouping via `--command \
         group` (BS1 must be closed end-to-end), stdout:\n{stdout}"
    );
    assert!(stdout.contains("DIFFER"), "Expected DIFFER in output, got:\n{stdout}");
    assert!(
        stdout.contains("Molecules matched: 0"),
        "the sole (single-MI) molecule must fail to match on the SEQ diff, got:\n{stdout}"
    );
}

/// Finding A (Task 8), end to end: see the "Missing MI tag behavior" section above (around
/// `test_grouping_mode_content_identical_but_mi_missing_in_both_bams_errors`). Two
/// content-identical, fully MI-less BAMs must hard-fail via the real `--command group` path
/// too, not just at the `molecule_join_compare` engine level.
#[test]
fn test_command_group_content_identical_but_mi_missing_in_both_bams_errors() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);
    let records = vec![mapped_record(b"read1", 100), mapped_record(b"read2", 200)];

    let bam1 = tmp.path().join("a.bam");
    let bam2 = tmp.path().join("b.bam");
    write_bam(&bam1, &header, &records);
    write_bam(&bam2, &header, &records);

    let (code, stdout, stderr) = run_compare_command(&bam1, &bam2, "group", &[]);

    assert_ne!(
        code,
        Some(0),
        "content-identical, fully MI-less BAMs must hard-fail via `--command group`, \
         stdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(stderr.contains("MI-tagged"), "expected the MI-less guard message, got:\n{stderr}");
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

    let (code, stdout, _stderr) = run_compare_command(&bam1, &bam2, "sort", &[]);

    assert_eq!(
        code,
        Some(0),
        "identical sorted BAMs must be IDENTICAL (exit 0) via --command sort:\n{stdout}"
    );
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

    let (code, stdout, _stderr) = run_compare_command(&bam1, &bam2, "sort", &[]);

    assert_eq!(
        code,
        Some(1),
        "a genuine mis-sort must DIFFER (exit 1) via --command sort:\n{stdout}"
    );
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

    let (code, stdout, stderr) = run_compare_command(&bam1, &bam2, "sort", &["--mode", "content"]);
    // An argument rejection is a runtime `anyhow` error (exit 1 with an `Error:` diagnostic on
    // stderr), which must be distinguished from a clean DIFFER (also exit 1, but a `RESULT`
    // line on stdout and no `Error:`) and from a panic/usage error (exit 101/2). Assert both
    // the code and the specific diagnostic so an unrelated failure can't satisfy this test.
    assert_eq!(
        code,
        Some(1),
        "--mode with --command sort must be rejected (exit 1), got stdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(
        stderr.contains("--mode is not valid with --command sort"),
        "rejection must name the offending flag, got stderr:\n{stderr}"
    );
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

    let (code, stdout, stderr) = run_compare_command(&bam1, &bam2, "sort", &["--ignore-order"]);
    // Mirrors `test_command_sort_rejects_explicit_mode`: assert the exact exit code and the
    // specific `--ignore-order` diagnostic, so a panic or an unrelated error can't stand in
    // for the intended argument rejection.
    assert_eq!(
        code,
        Some(1),
        "--ignore-order with --command sort must be rejected (exit 1), got stdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(
        stderr.contains("--ignore-order is not valid with --command sort"),
        "rejection must name the offending flag, got stderr:\n{stderr}"
    );
}

// ============================================================================
// Header precondition gate (Task 2): `require_compatible_headers` is wired into
// `CompareBams::execute` as a hard exit before any engine/record comparison runs.
// ============================================================================

/// Builds a header declaring `SO:unsorted GO:query` with the given `SS` value and a single
/// `@SQ` (`chr1`, 10000) — used to exercise the bare-vs-prefixed template-coordinate `SS`
/// spelling through the real `@HD` a BAM file carries (mirrors
/// `create_queryname_sorted_header` above, but for the template-coordinate `SO`/`GO`
/// combination).
fn header_with_ss(ss: &str) -> Header {
    use bstr::BString;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::{
        Map, map::Header as HeaderRecord, map::ReferenceSequence,
    };
    use std::num::NonZeroUsize;

    let HeaderTag::Other(tag) = HeaderTag::from(*b"SO") else { unreachable!() };
    let mut builder = HeaderRecordMap::<HeaderRecord>::builder().insert(tag, "unsorted");
    let HeaderTag::Other(tag) = HeaderTag::from(*b"GO") else { unreachable!() };
    builder = builder.insert(tag, "query");
    let HeaderTag::Other(tag) = HeaderTag::from(*b"SS") else { unreachable!() };
    builder = builder.insert(tag, ss);
    let header_map = builder.build().expect("valid header map");

    let reference_sequence =
        Map::<ReferenceSequence>::new(NonZeroUsize::new(10000).expect("non-zero"));
    Header::builder()
        .set_header(header_map)
        .add_reference_sequence(BString::from("chr1"), reference_sequence)
        .build()
}

/// Writes two BAMs encoding the identical single paired molecule (one MI group), one with
/// fgumi's bare `SS:template-coordinate` `@HD` spelling and the other with fgbio's
/// SO-prefixed `SS:unsorted:template-coordinate` spelling. `@SQ` and every record are
/// otherwise identical — the only difference between the two files is the `SS` tag's
/// spelling, which `require_compatible_headers`/`sort_order_from_header` must normalize to
/// the same `SortOrder::TemplateCoordinate` (regression: the old byte-wise `@HD` compare in
/// `compare_hd` treated these as different `SS` values and false-failed every
/// fgumi-vs-fgbio grouping comparison).
fn write_same_molecule_two_ss_spellings(dir: &Path) -> (PathBuf, PathBuf) {
    let header_bare = header_with_ss("template-coordinate");
    let header_prefixed = header_with_ss("unsorted:template-coordinate");
    let records = paired_record_pair(b"read1", "1");

    let bam_bare = dir.join("bare.bam");
    let bam_prefixed = dir.join("prefixed.bam");
    write_bam(&bam_bare, &header_bare, &records);
    write_bam(&bam_prefixed, &header_prefixed, &records);
    (bam_bare, bam_prefixed)
}

/// Thin wrapper over [`run_compare_in_process`] naming the `grouping` mode explicitly at
/// call sites that only care about the header-gate behavior, not the other `run_compare_*`
/// knobs.
fn compare_bams_grouping(bam1: &Path, bam2: &Path) -> bool {
    run_compare_in_process(bam1, bam2, "grouping", &[])
}

/// Writes a minimal single-record BAM declaring `SO:coordinate` and one `@SQ` with the
/// given name/length, for exercising a genuine `@SQ` dictionary mismatch between two files.
fn write_bam_with_ref(dir: &Path, filename: &str, ref_name: &str, ref_len: usize) -> PathBuf {
    let header = create_coordinate_sorted_header(ref_name, ref_len);
    let path = dir.join(filename);
    write_bam(&path, &header, &[mapped_record(b"read1", 100)]);
    path
}

/// Runs `CompareBams::execute()` in-process under `--mode content` and returns the error's
/// full display string. Panics if `execute` succeeds, or if it fails with a
/// [`CompareMismatch`] (a genuine record-level DIFFER) rather than the hard header-gate
/// error this helper exists to capture — the two must never be confused (see
/// `run_compare_in_process`'s doc comment for the same distinction in the success path).
fn run_compare_content_expect_err(bam1: &Path, bam2: &Path) -> String {
    let cmd = CompareBams::try_parse_from([
        "bams",
        bam1.to_str().unwrap(),
        bam2.to_str().unwrap(),
        "--mode",
        "content",
    ])
    .expect("failed to parse compare bams args");
    match cmd.execute("fgumi compare bams") {
        Ok(()) => panic!("expected compare bams to hard-fail on incompatible headers"),
        Err(e) if e.is::<CompareMismatch>() => {
            panic!("expected a hard header-gate error, got a CompareMismatch (DIFFER): {e:#}")
        }
        Err(e) => format!("{e:#}"),
    }
}

/// Two grouped BAMs whose only `@HD SS` difference is the bare-vs-prefixed template-coordinate
/// spelling must NOT be reported as a header difference (regression: the old byte-wise SS
/// compare false-failed every fgumi-vs-fgbio grouping comparison).
#[test]
fn grouping_accepts_bare_vs_prefixed_template_coordinate_ss() {
    let tmp = TempDir::new().unwrap();
    // identical single molecule, one file bare-SS, one file prefixed-SS in the @HD.
    let (bam_bare, bam_prefixed) = write_same_molecule_two_ss_spellings(tmp.path());
    assert!(
        compare_bams_grouping(&bam_bare, &bam_prefixed),
        "same order, different SS spelling must MATCH, not report a header diff"
    );
}

/// A `@SQ` dictionary mismatch is a hard failure (nonzero), not a record cascade.
#[test]
fn compare_hard_fails_on_sq_dictionary_mismatch() {
    let tmp = TempDir::new().unwrap();
    let a = write_bam_with_ref(tmp.path(), "a.bam", "chr1", 1000);
    let b = write_bam_with_ref(tmp.path(), "b.bam", "chr1", 2000);
    let err = run_compare_content_expect_err(&a, &b);
    assert!(err.contains("@SQ"), "got: {err}");
}

// ============================================================================
// Universal sort-order verification precondition (Task 3): `require_compatible_headers`'
// declared order is a header claim, not proof — `Content` mode pairs records purely by
// position, so a file that DECLARES an order but whose records don't actually honor it must
// be rejected up front rather than silently corrupting the positional pairing. Orderless
// (undeterminable) pairs — e.g. `extract`/`fastq`/`zipper` output — have nothing to verify
// and must proceed untouched.
// ============================================================================

/// A BAM declaring `SO:coordinate` whose records are genuinely out of coordinate order
/// (pos 500 before pos 100) — the `@HD` claim must not be trusted at face value.
fn write_bam_declaring_coordinate_but_unsorted(dir: &Path) -> PathBuf {
    let header = create_coordinate_sorted_header("chr1", 10000);
    let records =
        vec![fragment_record(b"read1", 0, 500, false), fragment_record(b"read2", 0, 100, false)];
    let path = dir.join("bad_coordinate.bam");
    write_bam(&path, &header, &records);
    path
}

/// A companion BAM, genuinely coordinate-sorted, sharing the same header/`@SQ` as
/// [`write_bam_declaring_coordinate_but_unsorted`] so the pair is header-compatible.
fn write_sorted_coordinate_bam(dir: &Path) -> PathBuf {
    let header = create_coordinate_sorted_header("chr1", 10000);
    let records =
        vec![fragment_record(b"read1", 0, 100, false), fragment_record(b"read2", 0, 500, false)];
    let path = dir.join("good_coordinate.bam");
    write_bam(&path, &header, &records);
    path
}

/// Builds a header declaring `SO:unsorted GO:query` with no `SS` tag — the shape
/// `extract`/`fastq`/`zipper` output declare: genuinely unordered data with no `SortOrder`
/// this engine can determine (see `require_compatible_headers`'s `Ok(None)` case).
fn header_unsorted_query_no_ss(ref_name: &str, ref_len: usize) -> Header {
    use bstr::BString;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::{
        Map, map::Header as HeaderRecord, map::ReferenceSequence,
    };
    use std::num::NonZeroUsize;

    let HeaderTag::Other(tag) = HeaderTag::from(*b"SO") else { unreachable!() };
    let mut builder = HeaderRecordMap::<HeaderRecord>::builder().insert(tag, "unsorted");
    let HeaderTag::Other(tag) = HeaderTag::from(*b"GO") else { unreachable!() };
    builder = builder.insert(tag, "query");
    let header_map = builder.build().expect("valid header map");

    let reference_sequence =
        Map::<ReferenceSequence>::new(NonZeroUsize::new(ref_len).expect("non-zero"));
    Header::builder()
        .set_header(header_map)
        .add_reference_sequence(BString::from(ref_name), reference_sequence)
        .build()
}

/// Writes two identical, orderless (`extract`-style `SO:unsorted GO:query`, no `SS`) BAMs.
fn write_identical_orderless_pair(dir: &Path) -> (PathBuf, PathBuf) {
    let header = header_unsorted_query_no_ss("chr1", 10000);
    let records = vec![mapped_record(b"read1", 100), mapped_record(b"read2", 200)];
    let a = dir.join("orderless_a.bam");
    let b = dir.join("orderless_b.bam");
    write_bam(&a, &header, &records);
    write_bam(&b, &header, &records);
    (a, b)
}

/// A BAM that DECLARES coordinate order but whose records are out of coordinate order must be
/// rejected up-front (don't trust the @HD tag), independent of record content.
#[test]
fn content_mode_rejects_records_not_in_declared_coordinate_order() {
    let tmp = TempDir::new().unwrap();
    let bad = write_bam_declaring_coordinate_but_unsorted(tmp.path()); // pos 500 then pos 100
    let good = write_sorted_coordinate_bam(tmp.path());
    let err = run_compare_content_expect_err(&good, &bad);
    assert!(err.to_lowercase().contains("order"), "got: {err}");
}

/// An orderless (extract-style `SO:unsorted GO:query`, no SS) pair has no verifiable order:
/// content comparison must proceed (no spurious order-verification error).
#[test]
fn content_mode_accepts_orderless_extract_headers() {
    let tmp = TempDir::new().unwrap();
    let (a, b) = write_identical_orderless_pair(tmp.path()); // SO:unsorted GO:query, no SS
    assert!(run_compare_in_process(&a, &b, "content", &[]), "orderless identical pair must MATCH");
}

// ============================================================================
// Task 6: two-sided streaming molecule hash-join driver (`molecule_join_compare`).
// ============================================================================

/// Writes two grouped BAMs encoding the same two molecules (`M_a` = reads `mol_a1`/`mol_a2`,
/// `M_b` = reads `mol_b1`/`mol_b2`), but with bam2's molecule *order* reversed relative to
/// bam1's and each molecule's MI number changed between the files. Content
/// (position/sequence) is otherwise identical, so
/// [`compare_molecule`](fgumi_lib::commands::compare)'s MI-invariant checks should find both
/// molecules equivalent despite neither the file order nor the MI numbering lining up — this
/// is exactly the scenario the hash-join (rather than a simple merge-join assuming
/// synchronized order) exists to handle.
fn write_two_molecules_reordered_and_renumbered(dir: &Path) -> (PathBuf, PathBuf) {
    let header = create_minimal_header("chr1", 10000);

    let mut bam1_records = vec![mi_record(b"mol_a1", 100, "5"), mi_record(b"mol_a2", 100, "5")];
    bam1_records.extend([mi_record(b"mol_b1", 200, "6"), mi_record(b"mol_b2", 200, "6")]);

    // bam2: M_b before M_a (reversed file order), MI renumbered (5->7, 6->9).
    let mut bam2_records = vec![mi_record(b"mol_b1", 200, "9"), mi_record(b"mol_b2", 200, "9")];
    bam2_records.extend([mi_record(b"mol_a1", 100, "7"), mi_record(b"mol_a2", 100, "7")]);

    let bam1 = dir.join("molecule_join_reordered1.bam");
    let bam2 = dir.join("molecule_join_reordered2.bam");
    write_bam(&bam1, &header, &bam1_records);
    write_bam(&bam2, &header, &bam2_records);
    (bam1, bam2)
}

/// Writes two grouped BAMs sharing one molecule (`M_a` = reads `mol_a1`/`mol_a2`), but bam1
/// has a second molecule (`M_b` = reads `mol_b1`/`mol_b2`, canonical id `mol_b1`) that bam2
/// lacks entirely.
fn write_extra_molecule_in_bam1(dir: &Path) -> (PathBuf, PathBuf) {
    let header = create_minimal_header("chr1", 10000);

    let m_a = vec![mi_record(b"mol_a1", 100, "1"), mi_record(b"mol_a2", 100, "1")];
    let m_b = vec![mi_record(b"mol_b1", 200, "2"), mi_record(b"mol_b2", 200, "2")];

    let mut bam1_records = m_a.clone();
    bam1_records.extend(m_b);
    let bam2_records = m_a;

    let bam1 = dir.join("molecule_join_extra1.bam");
    let bam2 = dir.join("molecule_join_extra2.bam");
    write_bam(&bam1, &header, &bam1_records);
    write_bam(&bam2, &header, &bam2_records);
    (bam1, bam2)
}

/// Reordered molecules across files, MI renumbered, still MATCH — no external sort, no window.
#[test]
fn molecule_join_matches_reordered_renumbered_molecules() {
    let tmp = TempDir::new().unwrap();
    // file1 molecules in order [M_a, M_b]; file2 in order [M_b, M_a], different MI numbers.
    let (bam1, bam2) = write_two_molecules_reordered_and_renumbered(tmp.path());
    let out = molecule_join_compare(&bam1, &bam2, 10).unwrap();
    assert_eq!(out.matched, 2);
    assert!(out.is_match(), "{out:?}");
}

/// A molecule present in only one file is a DIFFER named by its canonical id.
#[test]
fn molecule_join_flags_molecule_only_in_one_file() {
    let tmp = TempDir::new().unwrap();
    let (bam1, bam2) = write_extra_molecule_in_bam1(tmp.path());
    let out = molecule_join_compare(&bam1, &bam2, 10).unwrap();
    assert!(!out.is_match());
    assert!(out.diff_details.iter().any(|d| d.contains("only in bam1")), "{:?}", out.diff_details);
}

/// Mirror of [`write_extra_molecule_in_bam1`] with bam1/bam2 swapped: bam2 has the extra
/// molecule this time. Exists solely so a bam1/bam2 swap bug in the residual-drain loop
/// (`pending1.drain()` vs `pending2.drain()` in `molecule_join_compare`) is caught — without
/// it, a bug that swapped those two drains would still pass `molecule_join_flags_molecule_only_in_one_file`
/// (which only ever produces a "bam1" residual).
fn write_extra_molecule_in_bam2(dir: &Path) -> (PathBuf, PathBuf) {
    let header = create_minimal_header("chr1", 10000);

    let m_a = vec![mi_record(b"mol_a1", 100, "1"), mi_record(b"mol_a2", 100, "1")];
    let m_b = vec![mi_record(b"mol_b1", 200, "2"), mi_record(b"mol_b2", 200, "2")];

    let bam1_records = m_a.clone();
    let mut bam2_records = m_a;
    bam2_records.extend(m_b);

    let bam1 = dir.join("molecule_join_extra_bam2_1.bam");
    let bam2 = dir.join("molecule_join_extra_bam2_2.bam");
    write_bam(&bam1, &header, &bam1_records);
    write_bam(&bam2, &header, &bam2_records);
    (bam1, bam2)
}

/// Task 6 Minor: the mirror of `molecule_join_flags_molecule_only_in_one_file` with the
/// extra molecule on bam2 instead of bam1 — catches a bam1/bam2 swap bug in the residual
/// drain that the bam1-only version can't.
#[test]
fn molecule_join_flags_molecule_only_in_bam2() {
    let tmp = TempDir::new().unwrap();
    let (bam1, bam2) = write_extra_molecule_in_bam2(tmp.path());
    let out = molecule_join_compare(&bam1, &bam2, 10).unwrap();
    assert!(!out.is_match());
    assert!(out.diff_details.iter().any(|d| d.contains("only in bam2")), "{:?}", out.diff_details);
}

/// Documents the accepted `RecordKey`-collision residual (`record_key.rs:20-24`): `RecordKey`
/// is collision-resistant, not collision-free — for a primary record it reduces to
/// `(name, segment)`, with no locus discriminator, so two records that share a name and
/// segment but differ in position/content collapse to one entry in `index_by_key`'s
/// `BTreeMap` (the later one, in file order, wins). This is a malformed molecule (`group`
/// should never actually emit two same-name/same-segment primaries in one molecule), but
/// when it happens, `compare_molecule` silently only ever "sees" the winner.
///
/// bam1's molecule holds THREE physical records: `r1`, plus two same-name/same-segment `dup`
/// records at different positions (200 and 300) that collide to one `RecordKey`. bam2's
/// molecule holds only TWO: `r1` and a `dup` matching the position (300) of whichever bam1
/// record wins the collapse (the later one in file order). The observed verdict is MATCH —
/// bam1's second physical `dup`@200 record vanishes from the comparison entirely, silently
/// hiding a real difference in record count between the two files. This test locks in that
/// *actual* observed behavior (not a should-be) so a future change to the collision
/// resolution order is a deliberate, visible decision rather than a silent behavior change.
#[test]
fn molecule_join_duplicate_record_key_within_a_molecule_collapses_silently() {
    let tmp = TempDir::new().unwrap();
    let header = create_minimal_header("chr1", 10000);

    // bam1: r1, plus two colliding same-name/same-segment "dup" records (pos 200 then 300).
    let bam1_records =
        vec![mi_record(b"r1", 100, "5"), mi_record(b"dup", 200, "5"), mi_record(b"dup", 300, "5")];
    // bam2: r1, plus only the "dup" record matching the collapse winner (pos 300).
    let bam2_records = vec![mi_record(b"r1", 100, "9"), mi_record(b"dup", 300, "9")];

    let bam1 = tmp.path().join("dup_key1.bam");
    let bam2 = tmp.path().join("dup_key2.bam");
    write_bam(&bam1, &header, &bam1_records);
    write_bam(&bam2, &header, &bam2_records);

    let out =
        molecule_join_compare(&bam1, &bam2, 10).expect("molecule_join_compare should succeed");
    // OBSERVED (not desired) verdict: MATCH. bam1's "dup"@200 record silently vanishes from
    // the comparison, so a genuine 3-vs-2 physical-record-count difference is not reported.
    assert!(
        out.is_match(),
        "documents the RecordKey-collision residual: bam1's colliding \"dup\"@200 record is \
         silently dropped by index_by_key's BTreeMap collapse, so this (malformed) case \
         currently reports MATCH despite bam1 having one more physical record than bam2: {out:?}"
    );
}
